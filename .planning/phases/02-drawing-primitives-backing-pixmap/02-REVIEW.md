---
phase: 02-drawing-primitives-backing-pixmap
reviewed: 2026-04-11T00:00:00Z
depth: standard
files_reviewed: 10
files_reviewed_list:
  - prpr/xvtest0.f
  - xvue/qt/CMakeLists.txt
  - xvue/qt/include/xvue_qt_api.h
  - xvue/qt/src/xvue_qt_api.cpp
  - xvue/qt/src/xvue_qt_canvas.cpp
  - xvue/qt/src/xvue_qt_canvas.h
  - xvue/qt/src/xvue_qt_state.cpp
  - xvue/qt/src/xvue_qt_state.h
  - xvue/qt/MEFISTO_POINT_AUDIT.md
  - xvue/qt/README_RESIZE.md
findings:
  critical: 1
  warning: 5
  info: 4
  total: 10
status: issues_found
---

# Phase 2: Code Review Report

**Reviewed:** 2026-04-11
**Depth:** standard
**Files Reviewed:** 10
**Status:** issues_found

## Summary

Phase 2 landed the backing-QPixmap + long-lived QPainter infrastructure and
the DRAW-01..09 primitive set. Architecture (single-backing model, DPR-aware
resize, long-lived painter ended/restarted across resize) is sound and the
ABI parity audit (57 symbols, MefistoPoint INTEGER*2(2,N) verification) is
thorough. The resize/DRAW-09 invariants in `xvue_qt_canvas.cpp` closely
match `README_RESIZE.md`, and the arc angle conversion (X11 1/64 deg ->
Qt 1/16 deg) is correctly documented and implemented.

The review found one Critical ABI-semantics drift: the four legacy rectangle
entry points (`xvbordrectangle_`, `xvfbordrectangle_`, `xvrectangle_`,
`xvfrectangle_`) are collapsed through a single helper that calls
`QPainter::drawRect`, which strokes AND fills. Legacy calls
`XDrawRectangle` (outline only) vs `XFillRectangle` (fill only) on separate
symbols — the fill/outline dichotomy is the whole reason there are four
symbols. This is invisible today only because `applyPen()` sets pen and
brush to the same `foreground_` color, but it will produce incorrect output
once Phase 3 unlocks independent pen/brush colors. Five Warnings cover
defensive null-guard gaps that are inconsistent with the rest of the file,
event-loop re-entrancy risks in per-primitive `processEvents`, and a subtle
`drawPolygon` outline drift in `xvface_` that mirrors the rectangle issue.

## Critical Issues

### CR-01: Rectangle ABI collapses fill and outline semantics

**File:** `xvue/qt/src/xvue_qt_api.cpp:42-50, 587-613`
**Issue:** All four legacy rectangle entry points route through
`xvue_qt_draw_rect_common`, which calls `st->painter_->drawRect(QRect(...))`.
`QPainter::drawRect` strokes the outline with the current pen **and** fills
with the current brush in a single call. The four legacy symbols are NOT
equivalent:

- `xvfbordrectangle_` (xvuelc.c:2440) -> `XDrawRectangle` — **outline only**
- `xvbordrectangle_`  (xvuelc.c:2457) -> `XDrawRectangle` — **outline only**
- `xvfrectangle_`     (xvuelc.c:2503) -> `XFillRectangle` — **fill only**
- `xvrectangle_`      (xvuelc.c:2521) -> `XFillRectangle` — **fill only**

The "bord" (border) vs "f/[none]" (fill) distinction is the only reason
four symbols exist. Under the current implementation `xvbordrectangle_`
will both fill and outline, and `xvrectangle_` will both fill and outline.
In Phase 2 this is invisible because `applyPen()` sets `pen_.setColor(foreground_)`
and `brush_ = QBrush(foreground_, ...)` — both `Qt::white` — so the stray
fill/outline has the same color as the intended primitive. As soon as
Phase 3 differentiates pen and brush colors (or once `xvfacetraits_` sets
separate edge/face colors), every caller of `xvbordrectangle_` will start
painting a filled interior and every caller of `xvrectangle_` will gain a
1-pixel outline that legacy does not draw.

**Fix:** Pass the intended operation into the helper and branch; never rely
on `drawRect` for single-mode output. Minimal patch:

```cpp
enum class RectMode { Outline, Fill };

inline void xvue_qt_draw_rect_common(int x, int y, int w, int h, RectMode m) {
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;
    const QRect r(x, y, w, h);
    if (m == RectMode::Outline) {
        // XDrawRectangle: stroke only, no fill
        QBrush saved = st->painter_->brush();
        st->painter_->setBrush(Qt::NoBrush);
        st->painter_->drawRect(r);
        st->painter_->setBrush(saved);
    } else {
        // XFillRectangle: fill only, no outline
        st->painter_->fillRect(r, st->painter_->brush());
    }
    if (win->canvas()) win->canvas()->update();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

// callers:
void proc(xvfbordrectangle)(...)  { xvue_qt_draw_rect_common(..., RectMode::Outline); }
void proc(xvbordrectangle) (...)  { xvue_qt_draw_rect_common(..., RectMode::Outline); }
void proc(xvfrectangle)    (...)  { xvue_qt_draw_rect_common(..., RectMode::Fill);    }
void proc(xvrectangle)     (...)  { xvue_qt_draw_rect_common(..., RectMode::Fill);    }
```

`xvtest0.f` exercises all four symbols in sequence at line 58-61 — add a
visual Check in the phase summary that the two "bord" rectangles show hollow
outlines and the two fill variants show solid filled blocks with no border.

## Warnings

### WR-01: xvface_ and xvfacetraits_ stroke an extra outline vs legacy XFillPolygon

**File:** `xvue/qt/src/xvue_qt_api.cpp:439, 534-535`
**Issue:** Legacy `xvface_` calls `XFillPolygon(..., Complex, CoordModeOrigin)`
which fills only (no edge stroke). Qt's `painter_->drawPolygon(poly, Qt::OddEvenFill)`
fills with the current brush **and** strokes the polygon boundary with the
current pen. Same latent-but-hidden drift as CR-01: pen and brush are both
`foreground_ = Qt::white` in Phase 2 so the extra stroke is invisible, but
it will appear as a 1-px (or wider, after `xvepaisseur_`) outline once
Phase 3 differentiates pen and brush colors. Also affects `xvfacetraits_`
at line 534 (fill step) — though line 535 then explicitly strokes the
outline with `drawPolygon(poly)` (no brush change), so the fill step's
implicit stroke is wasted work and the two calls may paint the same edge
pixels twice with different pens once colors diverge.

**Fix:** For `xvface_`, wrap the fill in a `Qt::NoPen` override:
```cpp
st->painter_->save();
st->painter_->setPen(Qt::NoPen);
st->painter_->drawPolygon(poly, Qt::OddEvenFill);
st->painter_->restore();
```
For `xvfacetraits_`, use the same `setPen(Qt::NoPen)` + `setBrush(face_color)`
for the fill step, then restore and draw the outline with `setBrush(Qt::NoBrush)` +
edge pen. This also lays the groundwork for honoring `ncf`/`nca` in Phase 3
without another semantic change.

### WR-02: Per-primitive processEvents enables re-entrant resize mid-batch

**File:** `xvue/qt/src/xvue_qt_api.cpp:49, 315, 359, 392, 442, 479, 515, 538, 576, 635, 659`
**Issue:** Nearly every draw primitive ends with
`QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents)`. A
pending `ResizeEvent` queued by the window manager during a Fortran drawing
loop will be delivered inside that `processEvents` call, which runs
`XvueCanvas::resizeEvent`, which **ends the painter, deletes the old
backing, allocates a new one, and re-begins the painter**. Any Fortran
caller that batches many draws (e.g. triacoul.f, fap3*d.f loops) will then
race the resize handler: primitives before the resize land on the
pre-resize backing (now deleted content-preserved) and primitives after
land on the new backing. The DRAW-09 "preserve top-left" invariant holds,
so the visible result is correct **if** the resize happens between calls.
The risk is that the Fortran main holds no pointer that outlives the call,
so this is probably safe today, but it costs effort to reason about and
it couples drawing correctness to Qt event delivery timing.

It is also inconsistent with the stated DRAW-01 "long-lived painter" model:
the point of a long-lived painter is to batch draws and flush at
`xvvoir_`/`xvpause_`. Flushing after every single primitive defeats the
model and risks a latent bug if someone adds a code path that dereferences
`st->backing_` or `st->painter_` around a call.

**Fix:** Remove `processEvents` from every primitive except `xvvoir_`
(explicit flush) and `xvpause_`. Replace `win->canvas()->update()` at the
tail of each primitive with a lightweight "mark dirty" flag that
`xvvoir_` honors, or simply rely on Qt's coalesced `update()` which
schedules a single paintEvent. The phase 2 primitives already call
`canvas()->update()` which is enough; the `processEvents` is what makes
the per-primitive sync visible. Document the new invariant in
`README_RESIZE.md`: "backing mutations are deferred until the next event
loop pump; Fortran callers that need a forced flush must call `xvvoir_`."

### WR-03: Primitive entry points dereference int*/float* without null checks

**File:** `xvue/qt/src/xvue_qt_api.cpp:446-467, 470-480, 588-613, 618-636, 642-660`
**Issue:** `xvtypetrait_`, `xvepaisseur_`, `xvtrait_`, `xvftrait_`,
`xv*rectangle_` (all four), `xvbordarcellipse_`, and `xvarcellipse_` all
dereference their pointer parameters without a null check:
```cpp
st->pen_style_ = *ptype;                          // line 453
st->painter_->drawLine(*x1, *y1, *x2, *y2);        // line 477
xvue_qt_draw_rect_common(*x, *y, *width, *height); // line 591, 598, 605, 612
const QRect bbox(*x - *width, *y - *height, ...);  // line 627, 651
```
This is inconsistent with every warn-once stub in the same file, which
carefully routes `(void)arg;` through to guard against stray null calls,
and with `xvpxfenetre_` (line 399) which has an explicit `if (!x || !y)`.
Fortran callers should never pass null, but defensive consistency matters
given the ABI surface.

**Fix:** Add symmetric null guards at the top of each primitive. Example:
```cpp
void proc(xvtrait)(int *x1, int *y1, int *x2, int *y2) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x1 || !y1 || !x2 || !y2) return;   // ADD
    auto& win = XvueApp::window_slot();
    // ... rest unchanged
}
```
Apply the same pattern to the other ten entry points listed above.

### WR-04: xvue_qt_draw_rect_common reads width/height before null-check

**File:** `xvue/qt/src/xvue_qt_api.cpp:588-613`
**Issue:** `xvue_qt_draw_rect_common(*x, *y, *width, *height)` at each
caller site dereferences all four pointers in the call expression itself.
If any pointer is null, the crash happens at the caller, not at the
helper's defensive checks. Even after WR-03 is fixed, the helper itself
takes raw ints and cannot protect itself. Low risk, but the helper's
signature is now unable to express "I need valid args" vs "I got nulls".

**Fix:** After fixing WR-03 at each caller, optionally change the helper
to take pointers and centralize the null guard, or keep the int signature
and rely on the callers being correct. Either works; pick one and apply
consistently across the file.

### WR-05: resizeEvent leaks painter_ if begin() fails on the new backing

**File:** `xvue/qt/src/xvue_qt_canvas.cpp:67-73`
**Issue:** If `state_->painter_->begin(new_backing)` fails (QPainter returns
false, painter left inactive), the canvas returns silently and every
subsequent primitive short-circuits on `!st->painter_->isActive()`. The
user sees a frozen display with no diagnostic — much worse than a warn-once
stderr line.

More subtly, if `new QPainter()` throws `std::bad_alloc` at line 68,
`new_backing` has already been assigned to `state_->backing_` at line 64
but `old_backing` has been deleted — state is fine. But if the allocation
at line 48 (`new QPixmap(device)`) throws **after** line 42 ended the old
painter, the old painter's device (the old backing) has been ended but
still exists, and the function propagates the exception without
re-begin()-ing the painter on any backing. `state_->painter_` ends up
pointing at a now-inactive painter with a dangling device reference.

**Fix:** Order the allocations before the teardown:
```cpp
QPixmap* new_backing = new QPixmap(device);   // allocate FIRST
new_backing->setDevicePixelRatio(dpr);
{ QPainter tmp(new_backing); tmp.fillRect(...); if (old_backing) tmp.drawPixmap(...); }

// Only now is it safe to tear down the old painter/backing
if (state_->painter_ && state_->painter_->isActive()) state_->painter_->end();
delete state_->backing_;
state_->backing_ = new_backing;

if (!state_->painter_) state_->painter_ = new QPainter();
if (!state_->painter_->begin(new_backing)) {
    std::fprintf(stderr,
        "xvue-qt: resizeEvent: QPainter::begin failed on %dx%d backing\n",
        device.width(), device.height());
    return;
}
state_->painter_->setRenderHint(QPainter::Antialiasing, true);
state_->applyPen();
```

## Info

### IN-01: xvfacetraits_ ignores ncf/nca face/edge colors (documented TODO)

**File:** `xvue/qt/src/xvue_qt_api.cpp:522`
**Issue:** Legacy `xvfacetraits_` (xvuelc.c:2055-2065) calls `xvcouleur(ncf)`,
fills with XFillPolygon, then `xvcouleur(nca)` and strokes with XDrawLines.
The Qt version comments `(void)ncf; (void)nca;  // TODO(phase 3)`. This is
correct for Phase 2 scope (no palette yet), but leaves the DRAW-03 fill/edge
dichotomy looking correct in the test driver only because foreground_ is
white and background_ is black. Worth an explicit check in the Phase 8 A/B
validation script.

**Fix:** No action this phase. Phase 3 should track this alongside CR-01
and WR-01 as a single "color divergence" test.

### IN-02: xvtest0.f cumulative sleep is ~28 s, no CI-friendly mode

**File:** `prpr/xvtest0.f:72, 91, 97`
**Issue:** `SLEEP(15) + SLEEP(10) + SLEEP(3)` = 28 s per run. Fine for
interactive visual verification but painful in any automated harness,
and `SLEEP` is a gfortran intrinsic extension (not F77 standard). Not a
bug; just be aware that re-running under `cbl_tout` will add ~30 s per
smoke run.

**Fix:** Optional — honor an env var like `MEFISTO_TEST_NOSLEEP` to skip
the holds when running non-interactively. Phase 3 test harness concern,
not a Phase 2 blocker.

### IN-03: XvueState destructor order is safe only if XvueWindow owns both

**File:** `xvue/qt/src/xvue_qt_state.cpp:34-44`, `xvue/qt/src/xvue_qt_canvas.h:23`
**Issue:** The destructor ends the painter and deletes the backing, but
`XvueCanvas::state_` is a raw pointer. If the canvas ever outlives its
state (e.g. a Phase 3 refactor that swaps owners), the `state_->backing_`
read in `paintEvent` dereferences freed memory. Today this is guarded by
"owned by XvueWindow (D-15)" contract in the header comment, but nothing
enforces it.

**Fix:** Optional — make `XvueState` non-copyable/non-movable explicitly
(`= delete`), and add an invariant comment in the destructor reminding
future editors that `XvueCanvas::state_` must be nulled or the canvas
must be destroyed first. Current hierarchy is fine; document the
constraint.

### IN-04: Phase 0 unused-macro residue in xvue_qt_api.h

**File:** `xvue/qt/include/xvue_qt_api.h:21-27`
**Issue:** The `#ifdef __GNUC__ / #  define proc(x) x##_ / #else / ... / #endif`
block is immediately followed by `#undef proc` + unconditional
`#define proc(x) x##_`. The first block is dead code — the `#undef` wipes
it, and only the unconditional definition survives. This was copied
verbatim from xvuelc.c:60-70 as stated in the header comment, but the
legacy comment is obscured by the double-definition dance and any reader
will wonder whether it's broken. Pure style.

**Fix:** Optional — collapse to a single definition with a comment:
```c
// Trailing-underscore Fortran name mangling (all supported compilers).
#define proc(x) x##_
```
Keep the "copied from xvuelc.c" reference to preserve the audit trail.

---

_Reviewed: 2026-04-11_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
