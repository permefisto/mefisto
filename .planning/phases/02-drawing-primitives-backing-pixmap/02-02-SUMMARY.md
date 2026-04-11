---
phase: 02-drawing-primitives-backing-pixmap
plan: 02
subsystem: xvue-qt-backing-pixmap-painter
tags: [qt6, backing-pixmap, single-painter, resize, antialiasing, wave-1]

requires:
  - phase: 02-drawing-primitives-backing-pixmap
    plan: 01
    provides: "ABI audit (29/29 INTEGER*2), README_RESIZE.md DRAW-09 convention, float* angle signatures, prpr/xvtest0.f DRAW-01..09 coverage"
provides:
  - "XvueState grown: backing_ (QPixmap*), painter_ (QPainter*), foreground_, pen_, brush_, pen_style_, pen_width_base_, applyPen(), ~XvueState()"
  - "XvueCanvas::resizeEvent — D-07 sequence: end/allocate/fill/blit/delete/swap/begin/hint/applyPen"
  - "XvueCanvas::paintEvent = single drawPixmap(0,0) blit (DRAW-01 invariant)"
  - "Real bodies: effacer_ (D-15), xvvoir_ (D-02), xvpxfenetre_ (D-23), xvfond_ D-24 extension"
  - "DRAW-01 single long-lived QPainter + DRAW-07 clear/flush/query + DRAW-08 antialiasing + DRAW-09 top-left preserve — all live"
affects:
  - "02-03 (Wave 2: primitive bodies dereference state_->painter_ and rely on isActive() == true when a window is open)"
  - "02-04 (Wave 3: resize gate can now be visually verified — drag-resize during SLEEP(1) preserves black background)"

tech-stack:
  added: []
  patterns:
    - "Single-long-lived-painter invariant (DRAW-01): painter_->begin() called ONLY from XvueCanvas::resizeEvent; end() called by resizeEvent and ~XvueState"
    - "Scoped temp QPainter for resize fill+blit: one begin/end pair for fillRect + old-content drawPixmap"
    - "Device-pixel backing tagged with setDevicePixelRatio(dpr) so painter commands take logical coordinates (Pitfall 8)"
    - "D-01 epilogue = canvas->update() + QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents); applied by every flushing entry point"

key-files:
  created:
    - xvue/qt/src/xvue_qt_state.cpp
    - .planning/phases/02-drawing-primitives-backing-pixmap/02-02-SUMMARY.md
  modified:
    - xvue/qt/src/xvue_qt_state.h
    - xvue/qt/src/xvue_qt_canvas.h
    - xvue/qt/src/xvue_qt_canvas.cpp
    - xvue/qt/CMakeLists.txt
    - xvue/qt/src/xvue_qt_api.cpp

key-decisions:
  - "D-04/D-05 applied: backing_ and painter_ live in XvueState, not in XvueCanvas. Owns the lifetime from first resizeEvent to ~XvueState."
  - "Scoped temp painter combines (d)(e) fillRect+drawPixmap into one begin/end pair on new_backing — avoids two paint sessions and Pitfall 7 device-delete-while-active risk."
  - "xvfond_ epilogue runs unconditionally (even when backing is absent) so the second open/close cycle in xvtest0.f stays symmetric with Phase 1."
  - "ExcludeUserInputEvents flag is non-negotiable in every D-01 epilogue site (Pitfall 6) — 4 call sites total: xvfermer_, effacer_, xvvoir_, xvfond_."

requirements-completed:
  - DRAW-01
  - DRAW-07
  - DRAW-08
  - DRAW-09

duration: ~20min
completed: 2026-04-11
---

# Phase 2 Plan 02: Wave 1 Backing-Pixmap + Single-Painter Lifecycle Summary

**DRAW-01 single long-lived QPainter + DRAW-09 top-left preserve + DRAW-08 antialiasing + DRAW-07 clear/flush/query landed. XvueState carries backing_/painter_/pen_/brush_/applyPen; XvueCanvas::resizeEvent enforces the D-07 sequence; effacer_/xvvoir_/xvpxfenetre_/xvfond_ have real bodies. `pp/ppxvtest0_qt` exits 0 with 11 warn-once lines remaining (down from 13) — the two that disappeared are DRAW-07 entry points and confirm Wave 1 bodies reached execution.**

## Performance

- **Duration:** ~20 min
- **Started:** 2026-04-11T18:10:00Z
- **Completed:** 2026-04-11T18:36:00Z
- **Tasks:** 3 / 3
- **Files modified:** 6 (2 created, 4 modified)

## Accomplishments

### Task 1 — XvueState grown (D-04..D-22)

`xvue/qt/src/xvue_qt_state.h` extended additively: `background_` remains the first field (D-04 invariant). New fields added in order: `foreground_ = Qt::white`, `QPixmap* backing_ = nullptr`, `QPainter* painter_ = nullptr`, `QPen pen_`, `QBrush brush_ = QBrush(Qt::white, SolidPattern)`, `int pen_style_ = 0`, `int pen_width_base_ = 0`. New method decls `void applyPen()` and `~XvueState()`.

`xvue/qt/src/xvue_qt_state.cpp` (new) implements `applyPen()` as the single pen+brush source of truth. Three cases:
- `case 0` → SolidLine, width = pen_width_base_
- `case 1` → DashLine, width = pen_width_base_
- `default` (type 2 or anything else) → DashLine, width = `std::max(1, pen_width_base_ * 2)` (D-17)

After rebuilding pen_ (color=foreground_, style, width) and brush_ = `QBrush(foreground_, SolidPattern)`, `applyPen()` pushes both to `painter_->setPen()/setBrush()` when `painter_->isActive()`.

`~XvueState()` enforces the DRAW-01 invariant: if `painter_` non-null and active, call `painter_->end()`; then `delete painter_; delete backing_`. Safe when either is already null.

`xvue/qt/CMakeLists.txt` — added `src/xvue_qt_state.cpp` to the `add_library(xvueqt STATIC …)` source list (alphabetical, between `xvue_qt_canvas.cpp` and `xvue_qt_window.cpp`).

### Task 2 — XvueCanvas paintEvent swap + resizeEvent lifecycle (D-04..D-09)

`xvue/qt/src/xvue_qt_canvas.h` — added `void resizeEvent(QResizeEvent* event) override;`, forward-declared `QResizeEvent`.

`xvue/qt/src/xvue_qt_canvas.cpp`:

- **paintEvent** body is now exactly one defensive-guarded line:
  ```cpp
  if (state_ && state_->backing_) {
      QPainter(this).drawPixmap(0, 0, *state_->backing_);
  }
  ```
  The Phase 1 `fillRect(rect(), background_)` is gone. No state mutation, no second operation (Pitfall 1).

- **resizeEvent** implements the D-07 sequence in exact order:
  1. **(a)** `state_->painter_->end()` (guarded by `isActive()`) — ends old painter BEFORE deleting its device (Pitfall 7)
  2. **(b)(c)** allocate `QPixmap* new_backing = new QPixmap(size() * devicePixelRatioF())`, then `new_backing->setDevicePixelRatio(dpr)` so painter commands on it take logical coordinates (Pitfall 8)
  3. **(d)(e)** scoped `QPainter tmp(new_backing)` — one begin/end pair — first `tmp.fillRect(new_backing->rect(), background_)` (README_RESIZE.md invariant 2, the uncovered region shows background), then `if (old_backing) tmp.drawPixmap(0, 0, *old_backing)` (DRAW-09 top-left preserve, README_RESIZE.md invariant 3)
  4. **(f)** `delete old_backing` — safe now that no painter is attached to it
  5. **(g)** `state_->backing_ = new_backing` — atomic swap
  6. **(h)(i)** if `painter_` is null allocate `new QPainter()`; then `state_->painter_->begin(new_backing)` — this is the **only** `painter_->begin()` call site in the entire codebase, enforcing DRAW-01 single-long-lived-painter
  7. **(j)** `state_->painter_->setRenderHint(QPainter::Antialiasing, true)` — DRAW-08, re-applied after every begin() because hints do not carry across begin/end (Pitfall 5)
  8. **(k)** `state_->applyPen()` — re-push pen+brush through the painter (D-22, Pitfall 5)

### Task 3 — Real bodies for effacer_/xvvoir_/xvpxfenetre_ + xvfond_ D-24 extension (DRAW-07)

`xvue/qt/src/xvue_qt_api.cpp` — added `#include <QPainter>` and `#include <QPixmap>` at the top.

- **effacer_ (D-15):** `ensure()` + `ASSERT_MAIN_THREAD()` + get window → if state/painter active/backing present: `st->painter_->fillRect(st->backing_->rect(), st->background_)` + D-01 epilogue. Warn-once dropped.

- **xvvoir_ (D-02):** `ensure()` + `ASSERT_MAIN_THREAD()` + canvas->update() + `QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents)`. Warn-once dropped.

- **xvpxfenetre_ (D-23):** `ensure()` + `ASSERT_MAIN_THREAD()` + null-param guard → if window/canvas missing set `*x = *y = 0`, else `*x = canvas->width(); *y = canvas->height()` (logical pixels, SHELL-06). Warn-once dropped.

- **xvfond_ (D-24):** Phase 1 body preserved — the icolor→QColor mapping with xvfond_range out-of-range warn-once stays intact. After `state->background_ = chosen` the Phase 2 extension adds: `painter_->fillRect(backing_->rect(), background_)` when painter active + backing present, then canvas->update() + `processEvents(ExcludeUserInputEvents)`.

All four bodies start with `XvueApp::ensure(); XVUE_QT_ASSERT_MAIN_THREAD();`. No body contains `lasopsc`, `courgb`, `counb`, `ypixels`, `concat`, or `fpo` — verified by grep (D-26, D-28 compliance).

Total `processEvents(QEventLoop::ExcludeUserInputEvents)` call sites in xvue_qt_api.cpp: **4** (xvfermer_ from Phase 1, effacer_, xvvoir_, xvfond_).

## Task Commits

1. **Task 1: Grow XvueState (backing_/painter_/pen_/brush_/applyPen/dtor)** — `774bbb6` (feat)
2. **Task 2: XvueCanvas paintEvent swap + resizeEvent backing lifecycle** — `b01acae` (feat)
3. **Task 3: Real bodies for effacer_/xvvoir_/xvpxfenetre_ + xvfond_ D-24 extension** — `1714a97` (feat)

## Files Created/Modified

- `xvue/qt/src/xvue_qt_state.h` — (modified) additive growth: foreground_, backing_, painter_, pen_, brush_, pen_style_, pen_width_base_, applyPen() decl, ~XvueState() decl
- `xvue/qt/src/xvue_qt_state.cpp` — (new) applyPen() body with 3-case switch; ~XvueState() ends painter, deletes painter + backing
- `xvue/qt/src/xvue_qt_canvas.h` — (modified) resizeEvent override declaration, forward-declares QResizeEvent
- `xvue/qt/src/xvue_qt_canvas.cpp` — (modified) paintEvent reduced to drawPixmap blit; resizeEvent implements D-07 sequence
- `xvue/qt/CMakeLists.txt` — (modified) xvue_qt_state.cpp registered in xvueqt target
- `xvue/qt/src/xvue_qt_api.cpp` — (modified) real bodies for effacer_, xvvoir_, xvpxfenetre_; xvfond_ D-24 extension; #include <QPainter>, <QPixmap> added

## Exact Lines Changed

### xvue/qt/src/xvue_qt_state.h — additive fields (after background_)
```cpp
QColor foreground_ = Qt::white;
QPixmap*  backing_ = nullptr;
QPainter* painter_ = nullptr;
QPen     pen_;
QBrush   brush_          = QBrush(Qt::white, Qt::SolidPattern);
int      pen_style_      = 0;
int      pen_width_base_ = 0;
void applyPen();
~XvueState();
```

### xvue/qt/src/xvue_qt_canvas.cpp — paintEvent (Phase 1 → Phase 2)
```cpp
// Before:
QPainter(this).fillRect(rect(), state_->background_);
// After:
if (state_ && state_->backing_) {
    QPainter(this).drawPixmap(0, 0, *state_->backing_);
}
```

### xvue/qt/src/xvue_qt_canvas.cpp — resizeEvent D-07 sequence
```cpp
if (state_->painter_ && state_->painter_->isActive()) state_->painter_->end();
QPixmap* old_backing = state_->backing_;
QPixmap* new_backing = new QPixmap(size() * devicePixelRatioF());
new_backing->setDevicePixelRatio(devicePixelRatioF());
{ QPainter tmp(new_backing);
  tmp.fillRect(new_backing->rect(), state_->background_);
  if (old_backing) tmp.drawPixmap(0, 0, *old_backing); }
delete old_backing;
state_->backing_ = new_backing;
if (!state_->painter_) state_->painter_ = new QPainter();
state_->painter_->begin(new_backing);
state_->painter_->setRenderHint(QPainter::Antialiasing, true);
state_->applyPen();
```

## pp/ppxvtest0_qt Run Proof (after Task 3)

```
Gtk-Message: 18:34:45.286: Failed to load module "colorreload-gtk-module"
Gtk-Message: 18:34:45.286: Failed to load module "window-decorations-gtk-module"
xvue-qt: stub xvtypetrait_ not implemented yet
xvue-qt: stub xvepaisseur_ not implemented yet
xvue-qt: stub xvtrait_ not implemented yet
xvue-qt: stub xvtraits_ not implemented yet
xvue-qt: stub xvface_ not implemented yet
xvue-qt: stub xvbordrectangle_ not implemented yet
xvue-qt: stub xvrectangle_ not implemented yet
xvue-qt: stub xvfrectangle_ not implemented yet
xvue-qt: stub xvfbordrectangle_ not implemented yet
xvue-qt: stub xvarcellipse_ not implemented yet
xvue-qt: stub xvbordarcellipse_ not implemented yet

 ===========================================
 Phase 1+2: cycle open/close + primitives
 ===========================================
 [xvtest0] premier appel XVINITGRAPHIQUE
 [xvtest0] premier appel XVFERMER
 [xvtest0] second appel XVINITGRAPHIQUE (reopen)
 [xvtest0] second appel XVFERMER

 [xvtest0] OK — cycle open/close/open/close + draws
EXIT: 0
```

**11 warn-once lines** remain (Wave 0 baseline was 13). The two that disappeared are `effacer_` and `xvvoir_` — confirming Wave 1 bodies executed. `xvpxfenetre_` and `xvfond_` are not called from xvtest0.f's current draw-coverage section so they do not appear as warn-once lines in either baseline (they only surface via other test drivers). The reopen cycle is clean; exit 0.

## Build Verification

- `cmake --build xvue/qt/build` — exit 0; verify_abi: nm count = 57, header count = 57; verify_no_exec: OK
- `bin/cbl_tout_qt` — exit 0 (full Qt linking, all pp/pp*_qt executables rebuilt)
- `bin/cbl_tout` — exit 0 (legacy X11 untouched; CLAUDE.md "compilation must never break" respected)
- `bin/cbxvtest0_qt` — exit 0; `pp/ppxvtest0_qt` ~94 KB
- `pp/ppxvtest0_qt` — exit 0; 11 remaining warn-once lines; reopen cycle clean

## Decisions Made

- **Scoped temp painter for (d)(e) fill+blit** — one begin/end pair on new_backing instead of two, as recommended in the plan's behavior notes. Avoids a second device-touch cycle and keeps the resize path single-pass.
- **xvfond_ D-01 epilogue runs unconditionally** — `processEvents` is called even when the window is absent or the backing hasn't been allocated yet. Keeps the flush semantics identical whether the palette change is a no-op or not, matching Phase 1 xvfermer_'s pattern.
- **Defensive null-guard in paintEvent** — the plan says "should not happen on X11/Wayland but costs nothing." Kept as-is. The guard checks both `state_` and `state_->backing_`; without it a pathological pre-first-resize paint would crash inside `*state_->backing_`.
- **CMake source list kept alphabetical** — `xvue_qt_state.cpp` inserted between `xvue_qt_canvas.cpp` and `xvue_qt_window.cpp` rather than appended at the end, matching the pattern other Phase 1 files used.

## Deviations from Plan

None — plan executed exactly as written. All three tasks landed atomically, each build and xvtest0_qt run passed on first try, both Qt and legacy X11 builds stayed green, verify_abi/verify_no_exec clean throughout.

## Issues Encountered

None. Environment variables (`MEFISTO`, `MEFISTOX`, `PATH`) exported at session start; no surprises from applyPen signature, no header-ordering issues, no moc re-run needed (QResizeEvent is a forward decl in the header, defined in .cpp via `<QResizeEvent>`).

## Authentication Gates

None — pure local C++/CMake/Fortran work, no network, no credentials.

## Next Phase Readiness

- **Wave 2 (02-03) unblocked**: primitive entry points (xvtrait_, xvtraits_, xvface_, 4×rectangle, xvarcellipse_/xvbordarcellipse_, xvtypetrait_/xvepaisseur_) can now assume `state_->painter_->isActive() == true` whenever a window is open. Each body is a 3-liner: dereference `state_->painter_`, call the Qt API with arguments from the Fortran integers, run the D-01 epilogue. MefistoPoint arrays dereference safely per Wave 0 audit.
- **Wave 3 (02-04) unblocked**: the resize-preserve path is live. Plan 02-04's manual checkpoint (drag-resize during SLEEP(1), expect black to survive) can now run once Wave 2 lands visible primitives.
- **DRAW-01 invariant grep-verifiable**: `painter_->begin` appears in exactly one place (`xvue_qt_canvas.cpp::resizeEvent`), and `painter_->end` appears in exactly two (same resizeEvent — for the previous backing — and `~XvueState`). Wave 2 reviewers can re-verify with a 2-line grep.
- **Progressive smoke gate advances**: warn-once baseline dropped from 13 to 11. Wave 2 should land 9 more drops, leaving only `xvpxfenetre_`/`xvfond_` which the xvtest0.f current section does not exercise (they are covered by other drivers).

## Self-Check: PASSED

Verified artifacts exist:
- FOUND: `xvue/qt/src/xvue_qt_state.h` (modified — backing_, painter_, applyPen, ~XvueState all present)
- FOUND: `xvue/qt/src/xvue_qt_state.cpp` (new — applyPen + dtor)
- FOUND: `xvue/qt/src/xvue_qt_canvas.h` (modified — resizeEvent override decl)
- FOUND: `xvue/qt/src/xvue_qt_canvas.cpp` (modified — drawPixmap blit, D-07 sequence)
- FOUND: `xvue/qt/CMakeLists.txt` (modified — xvue_qt_state.cpp in xvueqt target)
- FOUND: `xvue/qt/src/xvue_qt_api.cpp` (modified — 4 real bodies, #include <QPainter>/<QPixmap>)

Verified commits exist:
- FOUND: `774bbb6` — Task 1 (XvueState growth)
- FOUND: `b01acae` — Task 2 (XvueCanvas lifecycle)
- FOUND: `1714a97` — Task 3 (effacer/xvvoir/xvpxfenetre/xvfond)

Verified builds green:
- `cmake --build xvue/qt/build` — exit 0, verify_abi 57/57, verify_no_exec OK
- `bin/cbl_tout_qt` — exit 0
- `bin/cbl_tout` — exit 0 (legacy X11)
- `pp/ppxvtest0_qt` — exit 0, 11 warn-once lines (down from 13)

---
*Phase: 02-drawing-primitives-backing-pixmap*
*Plan: 02 (Wave 1)*
*Completed: 2026-04-11*
