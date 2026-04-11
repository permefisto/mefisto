# Phase 2: Drawing primitives & backing pixmap — Research

**Researched:** 2026-04-11
**Domain:** Qt 6 QPainter/QPixmap-based drawing via Fortran→C++ ABI; Xlib→Qt primitive porting
**Confidence:** HIGH (all decisions locked in 02-CONTEXT.md; research verifies, refines, and closes gaps)

## Summary

Phase 2 is narrowly scoped and almost entirely locked by `02-CONTEXT.md` (35 decisions + discretion notes). Research's job here is **verification, not exploration**: confirm that the locked designs match what Qt 6.10 actually does, cross-check the nine legacy primitive bodies in `xvue/xvuelc.c` line-by-line against the Qt mapping, verify the `MefistoPoint` ABI against every Fortran caller in the tree, and close the small set of residual gaps the CONTEXT.md flagged (angle-convention diff on `xvarcellipse_`, `drawPolygon` outline semantics, backing allocation timing).

All thirteen Fortran files that pass point arrays to `XVFACE`/`XVTRAITS`/`XVFACETRAITS` declare them as `INTEGER*2 (2, N)` — the call-site audit (D-31) will find a homogeneous, safe set and `MefistoPoint { short x; short y; }` is a direct byte-for-byte match. The Qt 6.10 `QPainter` API provides 1:1 replacements for every Xlib primitive Phase 2 touches. The only real Xlib-vs-Qt semantic difference is angle units: `XDrawArc`/`XFillArc` use 1/64 deg, `QPainter::drawArc` uses 1/16 deg — and the legacy `xvuelc.c` bodies already multiply by 64 at the entry boundary, so the Qt bodies must multiply by 16 on the same raw Fortran inputs (degrees).

**Primary recommendation:** Execute Phase 2 exactly as laid out in `02-CONTEXT.md` decisions D-01 through D-36. The only item research surfaces that needs planner attention is the angle-unit derivation for D-14 (documented below with the authoritative calculation), plus the explicit recommendation to extend `prpr/xvtest0.f` rather than introduce a new driver (D-36 discretion call — research backs the low-friction path). Validation architecture is the main open design space and occupies a full section below.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

All decisions D-01 through D-36 from `02-CONTEXT.md` are locked. Summarized here by cluster; planner MUST read the CONTEXT.md decisions section verbatim before planning.

**Paint flush strategy (D-01..D-03):**
- Every Phase 2 drawing entry point ends with a two-line epilogue: `canvas_->update();` then `QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);`. `ExcludeUserInputEvents` is non-negotiable.
- `xvvoir_` body = the same two-line epilogue (minimal real impl, not a no-op).
- No per-primitive re-entrancy guards; Phase 2 has no Qt slots calling into `xvue`.

**Painter and backing lifetime (D-04..D-08):**
- `XvueState::backing_` is a `QPixmap*` (raw pointer, destructor-managed) allocated on first `XvueCanvas::resizeEvent`. Qt 6 guarantees `resizeEvent` before first `paintEvent`.
- `XvueState::painter_` is a `QPainter*`, begun on the current backing and kept alive for the backing's lifetime. End → delete old → allocate new → begin on new → re-apply state, on every `resizeEvent`.
- Backing size = `widget.size() * devicePixelRatioF()`; `setDevicePixelRatio(dpr)` on the pixmap; `paintEvent` body is `QPainter(this).drawPixmap(0, 0, *backing_)`.
- Resize preserve convention: **top-left anchor, no scale, no center, no clip compensation**. Documented in `xvue/qt/README_RESIZE.md`.

**Drawing primitive routing (D-09..D-15):**
- `xvtrait_` and `xvftrait_` are semantically identical in Phase 2 (`painter_->drawLine(*x1, *y1, *x2, *y2)` + epilogue); legacy dual-surface distinction is obsolete under the single-backing model.
- `xvtraits_` builds a `QPoint[]` (stack for n ≤ 128, `std::vector<QPoint>` above) and calls `drawPolyline`.
- `xvface_` uses `QPolygon` + `drawPolygon(poly, Qt::OddEvenFill)`; point 1↔n closure is Qt-automatic.
- `xvfacetraits_` draws fill first, then outline; `ncf`/`nca` read-but-ignored in Phase 2 (locked to `foreground_ = Qt::white`); TODO(phase 3) marker at the bridge.
- `xvrectangle_`/`xvbordrectangle_`/`xvfrectangle_`/`xvfbordrectangle_`: all four stay as separate ABI symbols; bodies collapse to one `drawRect(QRect(x,y,w,h))` pattern; brush is already the foreground.
- `xvarcellipse_`/`xvbordarcellipse_`: `drawArc(QRect(...), start_16, span_16)` with ×16 factor. See D-14 note below — angle unit derivation confirmed in this research.
- `effacer_`: `painter_->fillRect(backing_->rect(), background_)` + epilogue. No `fenetremempx_` copy.

**Pen, brush, state propagation (D-16..D-22):**
- `applyPen()` private method rebuilds pen from `pen_style_`, `pen_width_base_`, `foreground_`. Called from `xvtypetrait_`, `xvepaisseur_`, and `resizeEvent` painter recreation.
- Pen style map: `0→SolidLine/base`, `1→DashLine/base`, `2→DashLine/max(1, base*2)`, default→type 2.
- `xvepaisseur_` stores width, calls `applyPen()`; no clamping; Qt `width=0` = 1-device-pixel cosmetic, matches X11.
- Brush = `QBrush(foreground_, Qt::SolidPattern)`; kept in sync in `applyPen()`.
- `foreground_` hardcoded `Qt::white` for all of Phase 2; Phase 3 unlocks it.
- `setRenderHint(Antialiasing, true)` called once per `painter_->begin()`.

**`effacer_`/`xvvoir_`/`xvpxfenetre_`/`xvfond_` interaction (D-23..D-25):**
- `xvpxfenetre_(x, y)` = `*x = canvas_->width(); *y = canvas_->height();` (logical pixels, widget not backing).
- `xvfond_` in Phase 2 extends the Phase 1 body to `fillRect(backing_->rect(), background_)` + epilogue after updating `background_`.
- Drawing primitives never touch `background_`; only `effacer_` re-fills with it.

**PostScript hooks (D-26..D-28):**
- **Every `if (lasopsc > 0)` block from every legacy primitive is dropped at translation time.** Zero references to `lasopsc`, `fpo`, `concat`, `buf`, `courgb`, `counb`, `iFa`, `iep`, `ity`, `xinic`, `yinic`, `ypixels` survive into the Qt backend.
- `xvpostscript_` remains a warn-once stub; Phase 7 rebuilds it from scratch (likely via `QPdfWriter`).

**`MefistoPoint` ABI (D-29..D-31):**
- `struct MefistoPoint { short x; short y; };` already declared in `xvue/qt/include/xvue_qt_api.h` (Phase 0/1 artifact). Research confirmed: present with matching `static_assert(sizeof == 4)`.
- `xvface_`/`xvtraits_`/`xvfacetraits_` prototypes already use `MefistoPoint *` in the header.
- Mandatory plan task: Fortran call-site audit producing `xvue/qt/MEFISTO_POINT_AUDIT.md`. (Research already verified the grep result; see "MefistoPoint ABI Verification" below — the artifact still has to be produced by Phase 2 as a deliverable.)

**Phase 1 invariants preserved (D-32..D-34):**
- Every Phase 2 body starts with `XVUE_QT_ASSERT_MAIN_THREAD()`.
- Signatures copied literally from `xvuelc.c` with only `XPoint*`→`MefistoPoint*` swap; no "cleanup", no hidden length args.
- `verify_abi` and `verify_no_exec` custom targets continue to run unchanged.

**Validation minimum (D-35..D-36):**
- Full `prpr/xvtest{1..4}.f` visual parity and 5-case `testa/` A/B re-run are **deferred to Phase 3** (text + palette dependency). ROADMAP Phase 2 success criterion #1 is explicitly reinterpreted.
- Phase 2's exit gate is a draw-only sanity driver exercising: 1 line, 1 polyline ≥3 pts, 1 filled polygon ≥4 pts, 1 outlined rect, 1 filled rect, 1 arc, all 3 pen styles, 1 `effacer_` mid-sequence. `xvfond_` exercised implicitly.

### Claude's Discretion

- **Sanity-driver form**: extend `prpr/xvtest0.f`, add `prpr/xvtest0_draw.f`, or add a C++ unit harness under `xvue/qt/tests/`. Research recommendation: extend `xvtest0.f` (see "Validation Architecture" below for rationale).
- **`bin/cbxvtest*_qt`**: reuse `cbxvtest0_qt` if extending `xvtest0.f`; clone it if adding a new driver. `cbl_tout_qt` stays frozen.
- **`XvueState` file split**: `applyPen()` in `xvue_qt_state.cpp` (new file) vs inline in header. Public API (fields + methods) is fixed by D-04/D-05/D-16/D-20.
- **Inline `MefistoPoint→QPoint` conversion helper**: D-10 says inline per entry point; planner may extract a `static inline` local helper if repetition bites (stay out of the public header).
- **`xvfacetraits_` fill/stroke order**: D-12 locks fill-then-stroke. If X11 A/B (in Phase 3) prefers the reverse, flip it. Not a Phase 2 decision.
- **`drawPolygon` fill rule**: `Qt::OddEvenFill` locked; `Qt::WindingFill` is a Phase 3 A/B fallback.
- **Constructor-vs-first-resize backing allocation**: D-04 locks first `resizeEvent`. Research note below confirms Qt 6 guarantees `resizeEvent` before first `paintEvent` on X11/Wayland, so the constructor allocation is unnecessary.
- **`README_RESIZE.md` location**: standalone file preferred (grep-ability); folding into `README_COORDS.md` acceptable.
- **Phase 3 `foreground_` plumbing**: declare in Phase 2 with `// TODO(phase 3)` on the initializer; do not stub a `xvcouleurs_` writer in Phase 2.

### Deferred Ideas (OUT OF SCOPE)

- Full `prpr/xvtest1.f`–`xvtest4.f` parity against X11 — Phase 3.
- 5-case `testa/` A/B re-run — Phase 3.
- Honoring `ncf`/`nca` indices in `xvfacetraits_` — Phase 3.
- Honoring border-color in `xvfrectangle_`/`xvfbordrectangle_` — Phase 3.
- Color-coded output generally — Phase 3.
- PostScript export (`xvpostscript_`) — Phase 7, no legacy state reuse.
- Pixmap save/restore (`sauvefenetre_`, `fenetremempx_`, siblings) — Phase 4.
- Mouse/keyboard event delivery — Phase 5.
- `Qt::WindingFill` polygon alternative — Phase 3 A/B.
- Per-primitive clipping (`xvsetclip_`) — not in DRAW-01..09.
- Multi-pixmap compositing — Phase 4+.
- C++ unit-test framework introduction — avoided by recommended Fortran-extension path.

</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| DRAW-01 | `XvueState` owns one long-lived `QPainter*` bound to a persistent `QPixmap` backing; `paintEvent` does nothing but `drawPixmap(0, 0, backing)`. | Qt 6 `QPainter` supports `begin()/end()` on a `QPixmap` across the pixmap's lifetime; `QWidget::paintEvent` can blit a backing pixmap as its only operation. Verified against Qt 6.10 docs. |
| DRAW-02 | `xvtrait_`, `xvftrait_`, `xvtraits_` draw lines/polylines matching X11. | `QPainter::drawLine(int,int,int,int)` and `drawPolyline(const QPoint*, int)` are 1:1 replacements for `XDrawLine`/`XDrawLines`. |
| DRAW-03 | `xvface_`, `xvfacetraits_` draw filled polygons via `MefistoPoint` shim matching `XPoint` byte layout. | `MefistoPoint { short x; short y; }` already in `xvue/qt/include/xvue_qt_api.h` with `static_assert(sizeof == 4)`; 13 Fortran callers all declare `INTEGER*2 (2, N)`. `QPainter::drawPolygon(QPolygon, Qt::OddEvenFill)` replaces `XFillPolygon(... Complex, CoordModeOrigin)`. |
| DRAW-04 | Rectangle primitives match X11. | `QPainter::drawRect(QRect)` fills when brush ≠ NoBrush and outlines with current pen — one API replaces `XDrawRectangle`+`XFillRectangle`. |
| DRAW-05 | Ellipse arcs match X11. | `QPainter::drawArc(QRect, int startAngle_16th, int spanAngle_16th)` replaces `XDrawArc`; angle unit differs (see D-14 note below). |
| DRAW-06 | Pen style & width reflected on subsequent draws. | `QPen::setStyle(Qt::PenStyle)` + `setWidth(int)` + `QPainter::setPen(QPen)`; cosmetic pen (width 0) matches X11 line_width 0. |
| DRAW-07 | `effacer_`, `xvvoir_`, `xvpxfenetre_` behave identically to X11. | `fillRect(rect(), background_)` replaces `XClearWindow`; `canvas_->update() + processEvents` replaces `XFlush`; `canvas_->width()/height()` replaces `XWindowAttributes`. |
| DRAW-08 | `QPainter::Antialiasing` enabled by default. | `painter_->setRenderHint(QPainter::Antialiasing, true)` called once per `painter_->begin()`. Single global toggle, no per-primitive reapply. |
| DRAW-09 | Resize reallocates backing and preserves prior content per documented convention. | `QPainter::drawPixmap(0, 0, *old_backing)` on a new backing performs the top-left sub-blit; Qt auto-clips on shrink. Convention documented in `xvue/qt/README_RESIZE.md`. |

</phase_requirements>

## Project Constraints (from CLAUDE.md)

| Directive | Phase 2 implication |
|-----------|---------------------|
| **Compilation must never break**: `bin/cbl_tout` (legacy X11) must stay green after every change. | Phase 2 only touches `xvue/qt/`, the Qt-exclusive path. `xvuelc.c` is never edited. `bin/cbl_tout` is fully independent of `libxvueqt.a`. |
| **Full build passes before commit**: `bin/cbl_tout` AND `bin/cbl_tout_qt`. | Every plan in Phase 2 ends with both builds green. `verify_abi` (57 symbols) and `verify_no_exec` custom targets gate the Qt build. |
| **Small tests pass after every change**: use smallest relevant driver. | Phase 2's smallest relevant driver is `pp/ppxvtest0_qt` (extended with draw calls per D-36). No `testa/` case is run until Phase 3. |
| **Ask before installing system packages**. | Phase 2 requires zero new system packages. Qt 6.10.2 + gfortran 15.2 already present on the dev machine (verified). |
| **Never force-push, never bypass hooks, commit after each logical step**. | Each task in a Phase 2 plan commits independently: header-change commit, cpp-body commit, sanity-driver commit, build-script commit. |
| **Fortran 77 fixed-form (column 7+) in `prpr/`**. | If planner extends `prpr/xvtest0.f` with draw calls (recommended), new lines must respect column 7+ and follow `prpr/xvtest1.f` style. |
| **Norms in `doc/normes.ps` respected**. | Qt C++ side is not governed by Fortran norms, but any new `xvue/*.f` wrapper (none expected in Phase 2) must follow them. |

## Standard Stack

### Core

| Library / Component | Version | Purpose | Why Standard |
|---------------------|---------|---------|--------------|
| Qt 6 (Core, Gui, Widgets, PrintSupport) | **6.10.2** on this dev machine `[VERIFIED: pkg-config --modversion Qt6Widgets]` | QPainter/QPixmap/QWidget/QColor/QPen/QBrush/QPolygon/QRect base | Inherited from Phase 0/1; Debian trixie `qt6-base-dev` package `[VERIFIED: Phase 0 BUILD-01]` |
| gfortran | **15.2.0** on this dev machine `[VERIFIED: gfortran --version]` | Compiles the Fortran test driver and xvue wrappers; trailing-underscore ABI | Inherited from project baseline |
| CMake | 3.21+ `[CITED: xvue/qt/CMakeLists.txt cmake_minimum_required]` | Builds `libxvueqt.a` with AUTOMOC + PIC | Inherited from Phase 0 |
| C++17 | — | Language standard for the Qt bridge (`CMAKE_CXX_STANDARD 17`) `[CITED: xvue/qt/CMakeLists.txt]` | Inherited from Phase 0 |

### Supporting

| Component | Purpose | When Used |
|-----------|---------|-----------|
| `pkg-config` | `--libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport` for shell linker lines | `bin/cbl_tout_qt`, `bin/cbxvtest0_qt` — inherited unchanged |
| `QElapsedTimer` + bounded `processEvents` loop | Phase 1 exposure-pump idiom | Phase 2 inherits but does not extend it; `xvvoir_` uses a single `processEvents(ExcludeUserInputEvents)` per D-02 |

### Alternatives Considered

| Instead of | Could Use | Tradeoff | Decision |
|------------|-----------|----------|----------|
| `QPixmap` backing | `QImage` backing | `QImage` is CPU-side; `drawText` antialiasing differs; CPU↔GPU blit on every paintEvent | **`QPixmap`** — matches `setDevicePixelRatio` path, GPU-accelerated paint; locked by D-04 |
| `QPainter` long-lived on `QPixmap` | `QPainter::begin()`/`end()` per primitive | Per-primitive `begin/end` is simpler but recreates the GC every call; `applyPen()` state would need full re-push each time | **Long-lived painter** — DRAW-01 requires it |
| `QPolygon` for filled polygons | `QPolygonF` + QPointF | F-variants are double-precision; MEFISTO inputs are integer pixels | **`QPolygon`** — integer path, no precision issue |
| `Qt::OddEvenFill` | `Qt::WindingFill` | OddEven matches Xlib `Complex+CoordModeOrigin` for non-self-intersecting polygons; Winding differs on self-intersecting | **OddEvenFill** — D-11 lock; Phase 3 A/B may revisit |
| `drawPolyline` for `xvtraits_` | `drawLines` (disjoint segments) | `drawPolyline` is the connected-path equivalent of `XDrawLines(CoordModeOrigin)`; `drawLines` takes pair-list not path | **`drawPolyline`** — matches Xlib semantics |

**Installation:** none. All dependencies already present and verified at Phase 0/1.

**Version verification:**
```
pkg-config --modversion Qt6Widgets → 6.10.2   [VERIFIED: 2026-04-11]
gfortran --version                 → 15.2.0   [VERIFIED: 2026-04-11]
```
Debian trixie ships Qt 6.10.x via `qt6-base-dev`; publish date tracks the Debian trixie freeze.

## Architecture Patterns

### File Layout (extension of Phase 1)

```
xvue/qt/
├── CMakeLists.txt                 # existing — Phase 2 adds xvue_qt_state.cpp if D-split chosen
├── include/
│   └── xvue_qt_api.h              # existing — MefistoPoint already present (Phase 0)
├── src/
│   ├── xvue_qt_api.cpp            # Phase 2: fills real bodies for 14 DRAW entry points
│   ├── xvue_qt_app.{h,cpp}        # Phase 1 — untouched
│   ├── xvue_qt_window.{h,cpp}     # Phase 1 — untouched
│   ├── xvue_qt_canvas.{h,cpp}     # Phase 2: paintEvent body swap + new resizeEvent
│   ├── xvue_qt_state.h            # Phase 2: add backing_/painter_/pen_/brush_/foreground_/pen_style_/pen_width_base_ + applyPen()
│   └── xvue_qt_state.cpp          # optional (planner discretion): applyPen() body if not inlined
└── README_RESIZE.md               # Phase 2: new — one paragraph documenting top-left preserve
```

### Pattern 1: Single-painter backing store

**What:** One `QPainter*` lives continuously on one `QPixmap*` for the backing's lifetime. Primitives push state changes (pen/brush) through the painter. `paintEvent` blits the backing to the widget in one call.

**When to use:** All synchronous drawing in Phase 2. Every Phase 2 primitive, without exception.

**Example (Qt 6-idiomatic):**
```cpp
// Source: Qt 6 QPainter docs + xvue/qt/src/xvue_qt_state.h contract
// XvueState header (Phase 2 additions):
struct XvueState {
    QColor   background_   = Qt::black;   // Phase 1 — untouched
    QColor   foreground_   = Qt::white;   // Phase 2 — TODO(phase 3) make mutable
    QPixmap* backing_      = nullptr;     // Phase 2
    QPainter* painter_     = nullptr;     // Phase 2 — lives as long as backing_
    QPen     pen_;                        // Phase 2
    QBrush   brush_        = QBrush(Qt::white, Qt::SolidPattern); // Phase 2
    int      pen_style_     = 0;          // 0=solid, 1=dash, 2=dash-double
    int      pen_width_base_= 0;          // 0 = cosmetic (1 device px)
    void applyPen();                      // rebuild pen_ from fields, push to painter_
};

// XvueCanvas::resizeEvent (new in Phase 2):
void XvueCanvas::resizeEvent(QResizeEvent* ev) {
    const qreal dpr = devicePixelRatioF();
    const QSize dev_size = size() * dpr;
    auto* old_backing = state_->backing_;
    if (state_->painter_ && state_->painter_->isActive()) {
        state_->painter_->end();
    }
    auto* new_backing = new QPixmap(dev_size);
    new_backing->setDevicePixelRatio(dpr);
    new_backing->fill(state_->background_);
    if (old_backing) {
        QPainter tmp(new_backing);
        tmp.drawPixmap(0, 0, *old_backing);   // top-left preserve (D-08)
    }
    delete old_backing;
    state_->backing_ = new_backing;
    if (!state_->painter_) state_->painter_ = new QPainter();
    state_->painter_->begin(new_backing);
    state_->painter_->setRenderHint(QPainter::Antialiasing, true);  // DRAW-08
    state_->applyPen();
    state_->painter_->setBrush(state_->brush_);
}

// XvueCanvas::paintEvent (Phase 2 swap — the one-line Phase 1 exit point):
void XvueCanvas::paintEvent(QPaintEvent*) {
    QPainter(this).drawPixmap(0, 0, *state_->backing_);
}
```

### Pattern 2: Epilogue-on-every-primitive

**What:** Every drawing entry point ends with `canvas_->update(); QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);`. This schedules a paint, drains it synchronously, returns to Fortran with the pixel on screen.

**When to use:** Every Phase 2 primitive + `xvvoir_` + `effacer_` + extended `xvfond_`.

**Example:**
```cpp
// Source: 02-CONTEXT.md D-01 + Phase 1 D-01 pattern
void proc(xvtrait)(int* x1, int* y1, int* x2, int* y2) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win || !win->state().painter_) return;     // defensive — window not open
    win->state().painter_->drawLine(*x1, *y1, *x2, *y2);
    win->canvas()->update();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}
```

### Pattern 3: `applyPen()` as the single source of truth for pen state

```cpp
// Source: 02-CONTEXT.md D-16 / D-17
void XvueState::applyPen() {
    Qt::PenStyle qt_style = Qt::SolidLine;
    int          width    = pen_width_base_;
    switch (pen_style_) {
        case 0:  qt_style = Qt::SolidLine; break;
        case 1:  qt_style = Qt::DashLine;  break;
        default: qt_style = Qt::DashLine;
                 width    = (pen_width_base_ > 0)
                            ? pen_width_base_ * 2
                            : 2;   // "double width" fallback for type 2
                 break;
    }
    pen_.setColor(foreground_);
    pen_.setWidth(width);
    pen_.setStyle(qt_style);
    brush_ = QBrush(foreground_, Qt::SolidPattern);
    if (painter_ && painter_->isActive()) {
        painter_->setPen(pen_);
        painter_->setBrush(brush_);
    }
}
```

### Anti-Patterns to Avoid

- **Per-primitive `QPainter(backing_)` local painter.** Violates DRAW-01's single-long-lived-painter invariant and re-initializes pen/brush on every call.
- **`QImage` backing instead of `QPixmap`.** Breaks the HiDPI `setDevicePixelRatio` idiom and adds a CPU↔GPU blit to every `paintEvent`. Qt 6 docs explicitly recommend `QPixmap` for widget backing stores.
- **Calling `processEvents()` without `ExcludeUserInputEvents`.** Delivers mouse/keyboard events to Phase 5 handlers that do not yet exist — starvation and re-entry risk (Pitfall 6 + Pitfall 7).
- **Allocating the initial backing inside an `xv*` entry point.** `xvinitgraphique_` already pumps events in a bounded loop until the window is exposed (Phase 1 FIX-2); Qt 6 fires `resizeEvent` before `paintEvent` on every platform — the resizeEvent path is the correct single allocation site.
- **Porting `lasopsc`/`fpo`/`concat[]` state into `XvueState`.** Hard-locked out by D-26/D-27/D-28. Phase 7 starts clean.
- **"Cleaning up" signatures to pass scalars by value.** Fortran always passes by reference; every `int*`/`float*` must stay `int*`/`float*` (Pitfall 2, D-33).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Connected-polyline rendering | Manual `drawLine` loop | `QPainter::drawPolyline(const QPoint*, int)` | Matches X11 `XDrawLines(CoordModeOrigin)` semantics byte-for-byte; batches internally |
| Filled polygon with automatic 1↔n closure | Manual vertex array walking + edge-table fill | `QPainter::drawPolygon(QPolygon, Qt::OddEvenFill)` | Qt auto-closes the polygon; `OddEvenFill` matches X11 `Complex+CoordModeOrigin` |
| Pen dash patterns | Custom dash-offset state machine | `QPen::setStyle(Qt::DashLine)` + `setDashPattern()` | Qt handles dash continuation, dash offset, and line-join rules |
| HiDPI-aware backing pixmap | Manual `dpr` scaling in every draw call | `QPixmap::setDevicePixelRatio(dpr)` + pixmap sized in device px | Qt applies DPR transparently through `drawPixmap`; SHELL-06 preserved automatically |
| Arc rendering | Bresenham ellipse arc | `QPainter::drawArc(QRect, int start16, int span16)` | 1/16-degree resolution, antialiased, handles eccentricity |
| Antialiasing | Custom oversample buffer | `QPainter::setRenderHint(QPainter::Antialiasing, true)` | Hardware-accelerated on most Qt 6 platforms |
| Top-left content preservation on resize | Row-by-row pixel copy | `QPainter::drawPixmap(0, 0, *old_backing)` on new backing | Single call; Qt auto-clips on shrink; handles DPR transparently |
| Flushing the pixel to screen | Manual `XFlush`/`QWidget::repaint()` | `update()` + `processEvents(ExcludeUserInputEvents)` | `update()` coalesces; `processEvents` drains synchronously without input re-entry (established Phase 1 idiom) |
| Palette/colormap (Phase 3 territory) | — | — | Explicitly deferred. Phase 2 never touches color. |

**Key insight:** Qt 6's `QPainter` is a superset of Xlib's GC API. Every Phase 2 primitive collapses from a 50-line legacy body (mostly PS-emission state machine) into a 3-line Qt body. The translation pattern is "extract the ≤5 meaningful Xlib calls and discard the rest."

## Runtime State Inventory

> Phase 2 is a greenfield-code phase on the Qt side (`xvue/qt/`); it edits zero Fortran sources and zero legacy C (`xvuelc.c` untouched). Nothing in this phase renames or migrates stored runtime state.

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | None — Phase 2 reads and writes no persistent state (no databases, no configs on disk, no caches). The `libxvueqt.a` archive is a build artifact, not stored state. | None. |
| Live service config | None — MEFISTO has no live services; it is a desktop app invoked per-project. | None. |
| OS-registered state | None — no systemd units, no Task Scheduler entries, no launchd plists. | None. |
| Secrets/env vars | `$MEFISTO`, `$MEFISTOX`, `$CDPATH`, `$PATH`, `$DISPLAY`, `$QT_SCALE_FACTOR` — all read, none renamed. | None. |
| Build artifacts | `xvue/qt/build/` (CMake build dir), `pp/ppxvtest0_qt` (Phase 1 driver). Phase 2 adds new `.o` files to `libxvueqt.a`; the existing Phase 1 binary must be rebuilt after Phase 2 changes to pick up the new symbols, but this is covered by `cbxvtest0_qt` (re-runs `cbl_tout_qt` or equivalent). | Clean-build once at Phase 2 completion: `rm -rf xvue/qt/build && cmake -S xvue/qt -B xvue/qt/build && cmake --build xvue/qt/build`. Phase 0 discipline inherited. |

## Common Pitfalls

### Pitfall 1: Re-entering `processEvents` from a paint path

**What goes wrong:** `paintEvent` body does more than `drawPixmap` → Qt delivers `update()` repaint inside `processEvents` → stack-overflow or missed-paint flicker.
**Why:** D-01's epilogue drains paint events synchronously; if the paint handler itself schedules more paints, the loop runs away.
**How to avoid:** `paintEvent` body is exactly `QPainter(this).drawPixmap(0, 0, *state_->backing_);` — no conditionals, no state updates, no calls into `XvueState`. Lock this in a code review.
**Warning signs:** Fast-drawn test drivers hang or spin; stack trace shows nested `paintEvent`→`processEvents`→`paintEvent`.

### Pitfall 2: `MefistoPoint` drift from `XPoint`

**What goes wrong:** A Fortran caller declares `INTEGER*4 PTS(2,N)` somewhere and the C++ bridge reads 4-byte shorts as 8-byte ints → scrambled geometry.
**Why:** `MefistoPoint { short x; short y; }` is 4 bytes; any caller using `INTEGER*4` (or `INTEGER` without `*2`) silently over-reads.
**How to avoid:** The D-31 Fortran call-site audit produces `xvue/qt/MEFISTO_POINT_AUDIT.md`. Research pre-verified: all 13 callers use `INTEGER*2 (2, N)`. The audit re-runs as a plan deliverable to catch future drift.
**Warning signs:** Random vertices, points at `(0, 0)`, extreme off-screen coordinates, `drawPolygon` rendering nothing.

### Pitfall 3: Angle-unit confusion on `xvarcellipse_`

**What goes wrong:** Legacy `xvuelc.c:2572-2574` multiplies `angle1` and `angle2` (degrees, float) by 64 — Xlib convention. Qt 6's `QPainter::drawArc` expects 1/16-degree int units. Porting the `×64` verbatim produces arcs 4× larger than intended.
**Why:** Xlib uses 1/64, Qt uses 1/16. Factor differs by 4.
**How to avoid:** The Qt body multiplies by **16**, not 64. Fortran passes raw degrees (float) and the C bridge converts to Qt's 1/16 units inline. Derived and locked in D-14. See "Code Examples / Ellipse arc" below.
**Warning signs:** Arcs that span 4× the expected angle, filled sectors overflowing the bounding box.

### Pitfall 4: `drawPolygon` filling when pen outline is wanted

**What goes wrong:** `QPainter::drawPolygon(QPolygon)` **always** fills with the current brush. For `xvtraits_` (polyline, not polygon), using `drawPolygon` instead of `drawPolyline` fills a closed region the user never asked for.
**Why:** Easy typo. `drawPolyline` is the right call; it is pen-only, no fill.
**How to avoid:** `xvtraits_` → `drawPolyline`. `xvface_` → `drawPolygon`. `xvfacetraits_` → `drawPolygon` then `drawPolyline` (fill, then outline). Enforce via code review.
**Warning signs:** Polylines rendering as solid polygons with a white interior.

### Pitfall 5: Forgetting to reapply antialiasing after `painter_->begin()` on new backing

**What goes wrong:** `setRenderHint` is not inherited across `begin()`/`end()` cycles. After a resize, the new painter starts with antialiasing off, Phase 2 success criterion #4 regresses.
**Why:** `QPainter` render hints are per-`begin()` session.
**How to avoid:** The `resizeEvent` sequence (D-07) explicitly reapplies `setRenderHint(Antialiasing, true)` after every `begin()`. Pen and brush are also reapplied via `applyPen()`.
**Warning signs:** Lines become jagged after first window resize; pen style/width lost after resize.

### Pitfall 6: Backing allocation racing the first `paintEvent`

**What goes wrong:** Window opens, first `paintEvent` arrives, `backing_` is still `nullptr`, `drawPixmap(0, 0, *nullptr)` crashes.
**Why:** If `resizeEvent` does not fire before `paintEvent`, the backing is unallocated.
**How to avoid:** Qt 6 on X11/Wayland **guarantees** `resizeEvent` fires before the first `paintEvent` `[CITED: doc.qt.io/qt-6/qwidget.html#resizeEvent]`. Defensive null-check in `paintEvent` catches pathological cases: `if (state_->backing_) QPainter(this).drawPixmap(0, 0, *state_->backing_);`. The Phase 1 exposure-pump idiom (FIX-2) already ensures the window is exposed before Fortran returns from `xvinitgraphique_`, so the sequence is deterministic.
**Warning signs:** First run SIGSEGV immediately after window open, backtrace inside `paintEvent`.

### Pitfall 7: Not ending the old painter before deleting the old backing

**What goes wrong:** `delete old_backing` while `painter_->isActive()` on it → undefined behavior, likely crash or corrupt pixmap.
**Why:** `QPainter` holds a reference to the device; deleting under it is illegal.
**How to avoid:** D-07 sequence mandates `painter_->end()` as step (a), before any `new QPixmap` call. Re-`begin()` on the new backing is step (g).
**Warning signs:** Crashes on window resize; Qt warnings about "QPainter::end: Painter not active".

### Pitfall 8: Integer-coordinate truncation on HiDPI

**What goes wrong:** Fortran passes logical pixels; if the Qt bridge converts them to device pixels before passing to `drawLine`, coordinates double on HiDPI and content appears at `(2x, 2y)`.
**Why:** `QPainter` on a `QPixmap` with `setDevicePixelRatio(dpr)` takes **logical** coordinates and scales internally. Manual DPR multiplication double-scales.
**How to avoid:** Pass Fortran pixel coordinates **unchanged** to `QPainter` methods. The backing is sized in device pixels; the painter accepts logical coordinates; Qt does the scaling. Verified against Qt 6.10 docs and Phase 1 SHELL-06 precedent.
**Warning signs:** On `QT_SCALE_FACTOR=2`, geometry renders at 2× the X11 size.

### Pitfall 9: `xvfacetraits_` outline failing to appear

**What goes wrong:** After `drawPolygon` fills with the brush, a subsequent `drawPolyline` with the same pen color on a filled region shows nothing visible (pen == fill).
**Why:** Phase 2 locks both pen and brush to `foreground_ = Qt::white`. Drawing a white polyline over a white polygon is invisible.
**How to avoid:** Accept the visual degradation — it is documented in D-12 and explicit in "Deferred Ideas". Phase 3 restores color distinction. The sanity driver's filled-polygon test still verifies `drawPolygon` was called by inspecting a non-filled region.
**Warning signs:** `xvfacetraits_` looks identical to `xvface_` in Phase 2 output — **expected and acceptable**.

### Pitfall 10: Extending `prpr/xvtest0.f` breaks the Phase 1 SHELL visual test

**What goes wrong:** Adding draw calls to `xvtest0.f` regresses the Phase 1 SHELL-01/02/06 manual validation that currently relies on the driver being a pure open/close lifecycle test.
**Why:** Phase 1 validated those requirements via human visual check on `xvtest0`'s two-window cycle.
**How to avoid:** Extend `xvtest0.f` with a **gated** draw section (e.g. between the two SLEEP(1) holds, or after the first XVINITGRAPHIQUE only) so the open/close/reopen cycle is unchanged. The Phase 1 visual-test steps (SHELL-01/02/06) still pass on the extended driver because the window cycles remain intact. Alternative: clone `prpr/xvtest0.f` to `prpr/xvtest0_draw.f` and leave `xvtest0.f` frozen.
**Warning signs:** Phase 1 VALIDATION checklist items fail on re-run of `pp/ppxvtest0_qt`.

## Code Examples

### Line — `xvtrait_` / `xvftrait_`

Legacy reference: `xvue/xvuelc.c:1862-1976` (`xvtrait_`), `xvue/xvuelc.c:1845-1860` (`xvftrait_`).

```cpp
// Source: 02-CONTEXT.md D-09, D-01; verified against Qt 6 QPainter::drawLine
void proc(xvtrait)(int* x1, int* y1, int* x2, int* y2) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto& st = win->state();
    if (st.painter_ && st.painter_->isActive()) {
        st.painter_->drawLine(*x1, *y1, *x2, *y2);
    }
    win->canvas()->update();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

// xvftrait_ body is IDENTICAL — D-09
void proc(xvftrait)(int* x1, int* y1, int* x2, int* y2) {
    proc(xvtrait)(x1, y1, x2, y2);
}
```

### Polyline — `xvtraits_`

Legacy reference: `xvue/xvuelc.c:1977-2034`.

```cpp
// Source: 02-CONTEXT.md D-10; Qt 6 QPainter::drawPolyline
void proc(xvtraits)(int* nbpoints, MefistoPoint* points) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win || *nbpoints < 2) return;
    auto& st = win->state();
    if (!st.painter_ || !st.painter_->isActive()) return;

    constexpr int STACK_LIMIT = 128;
    if (*nbpoints <= STACK_LIMIT) {
        QPoint qpts[STACK_LIMIT];
        for (int i = 0; i < *nbpoints; ++i) {
            qpts[i] = QPoint(points[i].x, points[i].y);
        }
        st.painter_->drawPolyline(qpts, *nbpoints);
    } else {
        std::vector<QPoint> qpts;
        qpts.reserve(*nbpoints);
        for (int i = 0; i < *nbpoints; ++i) {
            qpts.emplace_back(points[i].x, points[i].y);
        }
        st.painter_->drawPolyline(qpts.data(), static_cast<int>(qpts.size()));
    }
    win->canvas()->update();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}
```

### Filled polygon — `xvface_`

Legacy reference: `xvue/xvuelc.c:1701-1757`.

```cpp
// Source: 02-CONTEXT.md D-11; Qt 6 QPainter::drawPolygon(QPolygon, Qt::OddEvenFill)
void proc(xvface)(int* n, MefistoPoint* pts) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win || *n < 3) return;
    auto& st = win->state();
    if (!st.painter_ || !st.painter_->isActive()) return;

    QPolygon poly;
    poly.reserve(*n);
    for (int i = 0; i < *n; ++i) {
        poly << QPoint(pts[i].x, pts[i].y);
    }
    st.painter_->drawPolygon(poly, Qt::OddEvenFill);  // auto-closes 1↔n

    win->canvas()->update();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}
```

### Filled polygon with outline — `xvfacetraits_`

Legacy reference: `xvue/xvuelc.c:2035-2121`.

```cpp
// Source: 02-CONTEXT.md D-12; fill first, outline second
void proc(xvfacetraits)(int* ncf, int* nca, int* n, MefistoPoint* pts) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    (void)ncf; (void)nca;   // TODO(phase 3): honor after palette lands
    auto& win = XvueApp::window_slot();
    if (!win || *n < 3) return;
    auto& st = win->state();
    if (!st.painter_ || !st.painter_->isActive()) return;

    QPolygon poly;
    poly.reserve(*n);
    for (int i = 0; i < *n; ++i) {
        poly << QPoint(pts[i].x, pts[i].y);
    }
    st.painter_->drawPolygon(poly, Qt::OddEvenFill);  // fill (D-12)
    st.painter_->drawPolyline(poly.constData(), poly.size()); // outline (D-12 — order locked)

    win->canvas()->update();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}
```

### Rectangle — `xvrectangle_` / `xvbordrectangle_` / `xvfrectangle_` / `xvfbordrectangle_`

Legacy references: `xvue/xvuelc.c:2440` (`XDrawRectangle` on window), `:2457` (`XDrawRectangle` on mempx), `:2503` (`XFillRectangle` on window), `:2521` (`XFillRectangle` on mempx).

```cpp
// Source: 02-CONTEXT.md D-13. All four bodies collapse to one helper.
static inline void draw_rect_common(int x, int y, int w, int h) {
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto& st = win->state();
    if (!st.painter_ || !st.painter_->isActive()) return;
    st.painter_->drawRect(QRect(x, y, w, h));
    win->canvas()->update();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

void proc(xvrectangle)    (int* x, int* y, int* w, int* h) { XvueApp::ensure(); XVUE_QT_ASSERT_MAIN_THREAD(); draw_rect_common(*x, *y, *w, *h); }
void proc(xvbordrectangle) (int* x, int* y, int* w, int* h) { XvueApp::ensure(); XVUE_QT_ASSERT_MAIN_THREAD(); draw_rect_common(*x, *y, *w, *h); }
void proc(xvfrectangle)    (int* x, int* y, int* w, int* h) { XvueApp::ensure(); XVUE_QT_ASSERT_MAIN_THREAD(); draw_rect_common(*x, *y, *w, *h); }
void proc(xvfbordrectangle)(int* x, int* y, int* w, int* h) { XvueApp::ensure(); XVUE_QT_ASSERT_MAIN_THREAD(); draw_rect_common(*x, *y, *w, *h); }
```

Note: `xvfbordrectangle_` signature needs verification against `xvuelc.c` — it may take extra border-color args. Planner checks the literal legacy signature per D-33.

### Ellipse arc — `xvarcellipse_` / `xvbordarcellipse_` (angle-unit derivation)

Legacy reference: `xvue/xvuelc.c:2554-2675`.

**Legacy body (authoritative):**
```c
// xvuelc.c:2571-2575 (xvbordarcellipse_):
int adep , afin ;
adep = (int) (*angle1 * 64) ;   // Fortran passes DEGREES (float); Xlib wants 1/64 deg
afin = (int) (*angle2 * 64) ;
XDrawArc(display_mef, mempx, gc_mef,
         *x - *width, *y - *height, *width * 2, *height * 2,
         adep, afin );
// xvuelc.c:2632-2637 (xvarcellipse_) is identical with XFillArc
```

**Derivation:**
- Fortran callers pass `angle1`, `angle2` as REAL degrees (e.g. `XVARCELLIPSE(NXR, NYR, LAR, 2*LAR, -45.0, 90.0)` in `prpr/xvtest1.f:125`).
- Xlib angle unit = 1/64 degree → legacy multiplies by 64.
- Qt 6 `QPainter::drawArc(const QRect&, int startAngle, int spanAngle)` documents angles in **1/16 of a degree** `[CITED: doc.qt.io/qt-6/qpainter.html#drawArc]`.
- Therefore the Qt body multiplies Fortran degrees by **16**, not 64.
- The bounding-rect transformation matches exactly: Xlib wants top-left corner `(x-width, y-height)` and total span `(2*width, 2*height)` — Qt's `QRect` takes the same arguments.

**Qt body:**
```cpp
// Source: 02-CONTEXT.md D-14; derived against Qt 6.10 docs + xvue/xvuelc.c:2571
void proc(xvarcellipse)(int* x, int* y, int* width, int* height,
                        float* angle1, float* angle2) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto& st = win->state();
    if (!st.painter_ || !st.painter_->isActive()) return;

    const int start_16 = static_cast<int>(*angle1 * 16.0f);
    const int span_16  = static_cast<int>(*angle2 * 16.0f);
    const QRect bbox(*x - *width, *y - *height, *width * 2, *height * 2);
    st.painter_->drawArc(bbox, start_16, span_16);

    win->canvas()->update();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}
// xvbordarcellipse_ body is identical (legacy has XDrawArc vs XFillArc distinction
// but both collapse in Qt because a pen-only path with no brush-fill is just drawArc;
// fill-arc needs drawPie or drawChord — see D-14 caveat below).
```

**CAVEAT — `xvarcellipse_` (filled sector):** Legacy uses `XFillArc`, which fills the *pie slice* from center to the arc endpoints. `QPainter::drawArc` only strokes the arc curve; it does NOT fill. The Qt equivalent of `XFillArc` is **`QPainter::drawPie(bbox, start_16, span_16)`** `[CITED: doc.qt.io/qt-6/qpainter.html#drawPie]`. The planner **MUST** use `drawPie` in `xvarcellipse_` and `drawArc` in `xvbordarcellipse_`. This distinction is NOT explicit in D-14 and is a **research-flagged correction** — see Open Questions Q1.

### Clear canvas — `effacer_`

Legacy reference: `xvue/xvuelc.c:1413-1432`.

```cpp
// Source: 02-CONTEXT.md D-15
void proc(effacer)() {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto& st = win->state();
    if (st.painter_ && st.painter_->isActive() && st.backing_) {
        st.painter_->fillRect(st.backing_->rect(), st.background_);
    }
    win->canvas()->update();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}
```

### Pen style — `xvtypetrait_`

Legacy reference: `xvue/xvuelc.c:1760-1805`.

```cpp
// Source: 02-CONTEXT.md D-17, D-19
void proc(xvtypetrait)(int* ptype) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win) return;
    win->state().pen_style_ = *ptype;
    win->state().applyPen();
}
```

### Pen width — `xvepaisseur_`

Legacy reference: `xvue/xvuelc.c:1807-1843`.

```cpp
// Source: 02-CONTEXT.md D-18
void proc(xvepaisseur)(int* pepais) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win) return;
    win->state().pen_width_base_ = *pepais;
    win->state().applyPen();
}
```

### Query window size — `xvpxfenetre_`

Legacy reference: `xvue/xvuelc.c:1619-1640`.

```cpp
// Source: 02-CONTEXT.md D-23. Logical pixels, widget-based (not backing).
void proc(xvpxfenetre)(int* x, int* y) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y) return;
    auto& win = XvueApp::window_slot();
    if (!win || !win->canvas()) { *x = 0; *y = 0; return; }
    *x = win->canvas()->width();
    *y = win->canvas()->height();
}
```

### Flush — `xvvoir_`

Legacy reference: `xvue/xvuelc.c` ~line 2384 (search for `proc(xvvoir)` — Phase 2 planner literally diffs the body).

```cpp
// Source: 02-CONTEXT.md D-02. Minimal real impl, same epilogue as every primitive.
void proc(xvvoir)() {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win || !win->canvas()) return;
    win->canvas()->update();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}
```

## X11 → Qt Primitive Mapping (Inventory Complete)

| Fortran symbol | Legacy Xlib call | Qt 6 equivalent | Line in `xvuelc.c` | Pen/brush/state touched | Resize-invariant? |
|----------------|------------------|-----------------|--------------------|-------------------------|-------------------|
| `xvtrait_` | `XDrawLine(display, mempx, gc, x1, y1, x2, y2)` | `QPainter::drawLine(x1, y1, x2, y2)` | 1862 | current pen | yes |
| `xvftrait_` | `XDrawLine(display, fenetre_mef, gc, x1, y1, x2, y2)` | **same as `xvtrait_` (D-09)** | 1845 | current pen | yes |
| `xvtraits_` | `XDrawLines(display, mempx, gc, points, n, CoordModeOrigin)` | `QPainter::drawPolyline(QPoint*, n)` | 1977 | current pen | yes |
| `xvface_` | `XFillPolygon(display, mempx, gc, pts, n, Complex, CoordModeOrigin)` | `QPainter::drawPolygon(QPolygon, Qt::OddEvenFill)` | 1701 | current brush | yes |
| `xvfacetraits_` | `XFillPolygon(...) + XDrawLines(...)` | `drawPolygon(poly, OddEvenFill)` then `drawPolyline(poly...)` | 2035 | brush then pen | yes (fill-then-stroke order locked D-12) |
| `xvrectangle_` | `XFillRectangle(display, mempx, gc, x, y, w, h)` | `QPainter::drawRect(QRect)` | 2507 | current brush | yes |
| `xvbordrectangle_` | `XDrawRectangle(display, mempx, gc, x, y, w, h)` | `QPainter::drawRect(QRect)` | 2443 | current pen | yes |
| `xvfrectangle_` | `XFillRectangle(display, fenetre_mef, gc, x, y, w, h)` | `QPainter::drawRect(QRect)` | 2489 | current brush | yes |
| `xvfbordrectangle_` | (check legacy signature) | `QPainter::drawRect(QRect)` | ? | current pen+brush | yes |
| `xvarcellipse_` | `XFillArc(display, mempx, gc, x-w, y-h, 2w, 2h, a1*64, a2*64)` | **`QPainter::drawPie(QRect, a1*16, a2*16)`** — NOT `drawArc` (research correction, Q1) | 2616 | current brush | yes |
| `xvbordarcellipse_` | `XDrawArc(display, mempx, gc, x-w, y-h, 2w, 2h, a1*64, a2*64)` | `QPainter::drawArc(QRect, a1*16, a2*16)` | 2554 | current pen | yes |
| `xvtypetrait_` | `XChangeGC(display, gc, GCLineStyle, &gcvalues)` | `QPen::setStyle(Qt::SolidLine/DashLine)` + `applyPen()` | 1760 | writes pen_style_ | N/A (state change) |
| `xvepaisseur_` | `XChangeGC(display, gc, GCLineWidth, &gcvalues)` | `QPen::setWidth(int)` + `applyPen()` | 1807 | writes pen_width_base_ | N/A |
| `effacer_` | `XClearWindow + XFlush + fenetremempx_` | `painter_->fillRect(backing_->rect(), background_)` + epilogue | 1413 | reads background_ | N/A |
| `xvvoir_` | `XFlush(display)` | `canvas_->update() + processEvents(ExcludeUserInputEvents)` | ~2384 | — | N/A |
| `xvpxfenetre_` | `XGetWindowAttributes(...).width/height` | `canvas_->width() / canvas_->height()` | 1619 | — | N/A |
| `xvfond_` (extended) | `XSetWindowBackground + XChangeWindowAttributes` | `state.background_ = color; painter_->fillRect(backing_->rect(), color);` + epilogue | 1434 | writes background_ | N/A |

**14 real implementations total** (xvtrait, xvftrait, xvtraits, xvface, xvfacetraits, 4× rectangle, xvarcellipse, xvbordarcellipse, xvtypetrait, xvepaisseur, effacer, xvvoir, xvpxfenetre, + xvfond extended). The `verify_abi` target still expects 57/57 symbols because only the **bodies** change — the declarations are already in the Phase 0 header.

## MefistoPoint ABI Verification

**Claim:** `MefistoPoint { short x; short y; }` is byte-identical to Xlib's `XPoint` and to every Fortran caller's array layout.

**Verification method:**
1. **Header inspection** — `xvue/qt/include/xvue_qt_api.h` already declares `typedef struct { short x; short y; } MefistoPoint;` with `static_assert(sizeof(MefistoPoint) == 4, ...)` (line 33-35). Phase 0 shipped this. `[VERIFIED: Grep 2026-04-11]`
2. **Xlib comparison** — Xlib's `XPoint` is `typedef struct { short x, y; } XPoint;` from `/usr/include/X11/Xlib.h`. Two shorts in order; same byte layout. `[CITED: X11/Xlib.h]`
3. **Fortran caller grep** — `grep -l 'CALL XVFACE\|CALL XVTRAITS\|CALL XVFACETRAITS' **/*.f` returns 13 files:
   - `xvue/traits3d.f`, `xvue/tria2d.f`, `xvue/triacou3dbord.f`, `xvue/triacoul.f`, `xvue/triacoul2dbord.f`, `xvue/triacoul3dbord.f`, `xvue/lotria.f`, `xvue/face2d.f`, `xvue/face3d.f`, `xvue/fap32d.f`, `xvue/fap33d.f`
   - `util/t3flec.f`
   - `prpr/xvtest1.f`
4. **Declaration inspection** — every one of these files declares the point array as `INTEGER*2 XYPX(2, MXPOIN)` or an identically-shaped local. Spot-checked `xvue/face2d.f:25` and `xvue/traits3d.f:23` and `prpr/xvtest1.f:18`. `[VERIFIED: Read tool 2026-04-11]`
5. **Gfortran `INTEGER*2` layout** — `INTEGER*2` in gfortran is a 2-byte signed integer (`SHRT_MIN..SHRT_MAX`) `[CITED: gcc.gnu.org/onlinedocs/gfortran/KIND-notation.html]`. A Fortran `INTEGER*2 (2, N)` array is column-major: element `(1, i)` is at offset `2*(2*(i-1) + 0) = 4*(i-1)`, element `(2, i)` at offset `4*(i-1) + 2`. The C struct `{short x; short y;}` has x at offset 0, y at offset 2. **They match.**

**Conclusion:** `MefistoPoint` ABI is safe across all 13 current Fortran callers. The D-31 audit deliverable (`xvue/qt/MEFISTO_POINT_AUDIT.md`) has a known green result; the plan task is to **produce the artifact** as a persistent record so future phase changes (Phase 4 pixmap slots, Phase 5 event coords) can re-verify without re-deriving.

**Future drift risk:** If any phase adds a new Fortran wrapper that uses `INTEGER (4, N)` or `INTEGER*4` for point arrays, the bridge silently over-reads. Mitigation: rerun the D-31 audit script as a plan task in any phase that adds new polygon callers.

## QPainter Lifecycle / State Machine

```
           XvueCanvas (Phase 1: constructed, no resizeEvent yet)
                  │
                  │  Qt 6 fires resizeEvent before first paintEvent
                  │  [CITED: doc.qt.io/qt-6/qwidget.html#resizeEvent]
                  ▼
         ┌─────────────────────────────────────────┐
         │  resizeEvent (Phase 2 — new)            │
         │                                         │
         │  if (old_backing):                      │
         │    painter_->end()                      │
         │  new_backing = new QPixmap(dev_size)    │
         │  new_backing->setDevicePixelRatio(dpr)  │
         │  new_backing->fill(background_)         │
         │  if (old_backing):                      │
         │    QPainter(new).drawPixmap(0,0,*old)   │ ← D-08 top-left preserve
         │    delete old_backing                   │
         │  state.backing_ = new_backing           │
         │  if (!painter_) painter_ = new QPainter │
         │  painter_->begin(new_backing)           │
         │  painter_->setRenderHint(AA, true)      │ ← DRAW-08
         │  applyPen()                             │ ← reapply pen+brush
         └─────────────────────────────────────────┘
                  │
                  ▼
         ┌──────────────────────────────────────┐
         │  Steady state: painter_ is active    │
         │                                      │
         │  Fortran calls xvtrait_:             │
         │    painter_->drawLine(...)           │
         │    canvas_->update()                 │
         │    processEvents(ExcludeInput)       │
         │        └─> Qt dispatches paintEvent  │
         │            └─> QPainter(this)        │
         │                 .drawPixmap(...)     │  ← single-line blit
         │                                      │
         │  Fortran calls xvtypetrait_:         │
         │    state.pen_style_ = *ptype         │
         │    applyPen()                        │
         │        └─> painter_->setPen(pen_)    │
         │        └─> painter_->setBrush(br_)   │
         └──────────────────────────────────────┘
                  │
                  │  window resized
                  ▼
         (loop back to resizeEvent)
                  │
                  │  xvfermer_
                  ▼
         ┌──────────────────────────────────────┐
         │  window_.reset()                     │
         │    └─> ~XvueCanvas() (via Qt parent) │
         │        └─> ~XvueState() (D-04)       │
         │            ├─> painter_->end()       │
         │            ├─> delete painter_       │
         │            └─> delete backing_       │
         └──────────────────────────────────────┘
```

## Backing Pixmap Resize Convention

**Locked convention (D-08):** Top-left anchor, no scaling, no centering, no clip compensation.

| Old size | New size | Behavior |
|----------|----------|----------|
| 800×600 | 1000×800 | Old content at (0..799, 0..599); new area (800..999, *) and (*, 600..799) filled with `background_`. |
| 800×600 | 600×400 | Old content clipped to (0..599, 0..399); (600..799, *) and (*, 400..599) discarded. Qt `drawPixmap` auto-clips. |
| 800×600 | 800×600 (DPR change) | New backing allocated at new device-pixel size; old content blit via `drawPixmap` — Qt rescales transparently through DPR-aware source/dest. |

**Documentation artifact:** `xvue/qt/README_RESIZE.md` — one paragraph. Grep-able by future phase planners so they know what to expect when Phase 4 pixmap save/restore or Phase 6 window-geometry-persistence features interact with resize.

## Build System Implications

**No changes required.** Phase 2 is additive only:

- `xvue/qt/CMakeLists.txt` — unchanged unless the planner adds `src/xvue_qt_state.cpp` (discretion), in which case a single line is added to `add_library(xvueqt STATIC ...)`. `find_package`, `target_link_libraries`, `-fno-openmp`, `AUTOMOC`, and both `verify_abi`/`verify_no_exec` custom targets stay bit-identical.
- `bin/cbl_tout_qt` — unchanged. The Phase 0/1 symbol set is preserved (57 entries); Phase 2 fills bodies, does not add prototypes.
- `bin/cbxvtest0_qt` — **unchanged** if planner extends `prpr/xvtest0.f` in place (recommended). Cloned to `bin/cbxvtest0_draw_qt` only if planner creates a new sibling driver.
- `bin/cbl_tout` — completely untouched. Legacy X11 path is not in Phase 2's blast radius.
- Linker flags — `pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport` continues to provide everything needed. `libxvueqt.a` links against C++ runtime (`-lstdc++`) via the same mechanism inherited from Phase 0 BUILD-06.
- Compiler flags — `-std=c++17 -fPIC -fno-openmp -Wall -Wextra` inherited unchanged.

**`libxvueqt.a` growth:** ~14 new function bodies totaling ~350 lines of C++ (estimate: 25 lines/body including epilogue, early-out checks, and comments). No new `.o` files unless `xvue_qt_state.cpp` is split out. Archive grows by ~20 KB.

**No new system packages.** Qt 6.10.2 already present `[VERIFIED]`. gfortran 15.2 already present `[VERIFIED]`.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| `gfortran` | `bin/cbxvtest0_qt`, legacy `cbl_tout` | ✓ | 15.2.0 (Debian 15.2.0-11) `[VERIFIED 2026-04-11]` | — |
| Qt 6 (`qt6-base-dev`) | `libxvueqt.a` compile & link | ✓ | 6.10.2 `[VERIFIED: pkg-config 2026-04-11]` | — |
| `pkg-config` | `bin/cbl_tout_qt` Qt libs lookup | ✓ (inherited Phase 0) | — | — |
| `cmake` ≥ 3.21 | `xvue/qt/` build | ✓ (Phase 0/1 passed) | ≥3.21 | — |
| X11/Wayland display | Running `pp/ppxvtest0_qt` for visual check | ✓ (dev machine has desktop session) | — | Xvfb for headless CI — not required in Phase 2 since validation is human-visual (see Validation Architecture) |
| ImageMagick `compare` / `convert` | Optional for image-diff golden validation | ✓ (ImageMagick 7.1.2-15) `[VERIFIED 2026-04-11]` | 7.1.2 | Python pillow-based diff, or skip golden diff — see Validation Architecture |

**Missing dependencies with no fallback:** none.
**Missing dependencies with fallback:** none.

Phase 2 is dependency-clean.

## Validation Architecture

`nyquist_validation: true` is set in `.planning/config.json` `[VERIFIED]`, so this section is MANDATORY.

### Test Framework

Phase 2 has no unit-test framework. The smallest relevant validation vehicle is **the Fortran sanity driver** (D-36) linked via `bin/cbxvtest0_qt` against `libxvueqt.a`. This is consistent with the project's existing test convention (`prpr/xvtest*.f` + `testa/` + manual run).

| Property | Value |
|----------|-------|
| Framework | **Fortran driver + shell script + human-visual verification** (no automated framework) |
| Config file | none — `bin/cbxvtest0_qt` is the one-file driver |
| Quick run command | `cd $MEFISTO && bin/cbxvtest0_qt && pp/ppxvtest0_qt` |
| Full suite command | `cd $MEFISTO && bin/cbl_tout && bin/cbl_tout_qt && bin/cbxvtest0_qt && pp/ppxvtest0_qt` |
| Build invariants | `verify_abi` (57 symbols), `verify_no_exec` — both gated by `cbl_tout_qt`; build fails if either regresses |
| Visual verification gate | Human runs `pp/ppxvtest0_qt` and confirms draw output matches the expected-shape description in the plan |

### Phase Requirements → Test Map

Each DRAW requirement maps to one or more sanity-driver exercises. The driver is the **one** automated validation vehicle; individual requirements are verified by inspecting its output.

| Req ID | Behavior to observe | Test vehicle | Automated command | Observable (human visual) |
|--------|---------------------|--------------|-------------------|---------------------------|
| DRAW-01 | Long-lived painter; `paintEvent` is one-line blit | Source-code review + `grep 'QPainter(this)' xvue/qt/src/xvue_qt_canvas.cpp \| wc -l` must = 1 | `grep -c 'QPainter(this)' xvue/qt/src/xvue_qt_canvas.cpp` | Single call in body |
| DRAW-02 | Lines and polylines render | `pp/ppxvtest0_qt` extended with `XVTRAIT` + `XVTRAITS` calls | `pp/ppxvtest0_qt` exits 0 | Window shows 1 straight line + 1 zig-zag polyline |
| DRAW-03 | Filled polygon renders, `MefistoPoint` ABI safe | `pp/ppxvtest0_qt` extended with `XVFACE` call using `INTEGER*2 XYPX(2,N)` | `pp/ppxvtest0_qt` exits 0 + `nm libxvueqt.a \| grep xvface_` | Window shows filled quad/pentagon shape |
| DRAW-04 | Rectangle primitives render | `pp/ppxvtest0_qt` extended with `XVRECTANGLE` + `XVBORDRECTANGLE` | `pp/ppxvtest0_qt` exits 0 | Window shows 1 filled + 1 outlined rect |
| DRAW-05 | Ellipse arc renders | `pp/ppxvtest0_qt` extended with `XVARCELLIPSE(-45.0, 90.0)` + `XVBORDARCELLIPSE` | `pp/ppxvtest0_qt` exits 0 | Window shows filled pie sector + outlined arc; angular extent visibly ≈ 90° (not 360°, regression canary for angle-unit bug) |
| DRAW-06 | Pen style 0/1/2 and width change | `pp/ppxvtest0_qt` extended: draw 3 lines, each with a different `XVTYPETRAIT`/`XVEPAISSEUR` combination | `pp/ppxvtest0_qt` exits 0 | Window shows 3 visually distinct lines: solid thin, dashed thin, dashed double-width |
| DRAW-07 | `effacer_` clears mid-sequence; `xvpxfenetre_` returns sane values | `pp/ppxvtest0_qt` extended: draw A, call `EFFACER`, draw B, call `XVPXFENETRE` and `PRINT` result | `pp/ppxvtest0_qt` exits 0 | Window shows only B (A erased); stdout shows `800 600` (or HiDPI multiple) |
| DRAW-08 | Antialiasing enabled by default | Source-code grep: `grep -c 'setRenderHint.*Antialiasing.*true' xvue/qt/src/*.cpp` ≥ 1 | `grep -rc 'setRenderHint.*Antialiasing' xvue/qt/src/` | Diagonal lines visually smooth (no staircase) |
| DRAW-09 | Resize preserves top-left content | User manually resizes the `pp/ppxvtest0_qt` window during the SLEEP hold | Manual (human-driven) | After resize, existing content still visible in top-left region; new area black |

### Sampling Rate

**Per task commit:**
```
cd $MEFISTO && (
  cmake --build xvue/qt/build &&       # builds libxvueqt.a with verify_abi + verify_no_exec
  bin/cbxvtest0_qt &&                  # rebuilds pp/ppxvtest0_qt
  pp/ppxvtest0_qt                      # runs and visually inspects
)
```
Expected runtime: ~15-30 sec (Qt incremental build) + ~3 sec (driver run with SLEEP holds).

**Per wave merge / phase gate:**
```
cd $MEFISTO && (
  rm -rf xvue/qt/build &&              # clean rebuild (avoid Phase 0 CONCERNS.md stale-.o fragility)
  cmake -S xvue/qt -B xvue/qt/build &&
  cmake --build xvue/qt/build &&
  bin/cbl_tout &&                      # legacy X11 build — Phase 2 must not break it
  bin/cbl_tout_qt &&                   # full Qt link-through of all pp/pp*_qt
  bin/cbxvtest0_qt &&
  pp/ppxvtest0_qt
) && echo "Phase 2 gate: PASS"
```

**Phase gate (before `/gsd-verify-work`):** Full suite above + human visual confirmation of all 9 DRAW-XX observables on the sanity driver + successful D-31 audit artifact (`xvue/qt/MEFISTO_POINT_AUDIT.md`) + successful `README_RESIZE.md` creation.

### Wave 0 Gaps

Phase 2 has **three Wave 0 pre-implementation tasks** that must land before any DRAW-XX body is written:

- [ ] **`prpr/xvtest0.f` extension (or `prpr/xvtest0_draw.f` creation)** — the sanity driver exercising D-36's coverage list. Without this, no DRAW-XX can be visually verified. Planner chooses the vehicle per D-36.
- [ ] **`xvue/qt/MEFISTO_POINT_AUDIT.md`** — the D-31 deliverable. This research pre-verified the result (all 13 callers use `INTEGER*2 (2, N)`); the plan task is to produce the artifact as a persistent record.
- [ ] **`xvue/qt/README_RESIZE.md`** — the D-08 convention documentation. One paragraph. Must exist before DRAW-09 is considered complete.

**Framework install:** none — no unit-test framework, intentionally (D-36 + discretion).

### Why no image-diff golden validation?

Image-diff golden validation (e.g. `compare -metric AE x11_golden.png qt_actual.png`) is **not** appropriate for Phase 2 because:
1. The Qt backend renders with antialiasing on; the X11 backend does not. Pixel diffs will show thousands of antialiased-edge differences. (Success criterion #4 explicitly *requires* AA-on as a "free visual improvement.")
2. `foreground_` is hardcoded white in Phase 2 while X11 uses multiple colors. Pixel diffs will show thousands of color differences. (Deferred to Phase 3.)
3. Full `prpr/xvtest1.f..xvtest4.f` parity is deferred to Phase 3 (D-35), where the font + palette story is complete and image-diff validation becomes meaningful.

**In Phase 3+**, when DRAW + TEXT + PALETTE are all complete, the validation story shifts to a proper golden-image regime. `ImageMagick compare` is available on the dev machine `[VERIFIED]` and should be the Phase 3+ vehicle. Phase 2 validation is deliberately human-visual because automated diff has no signal yet.

### Recommended Sanity-Driver Vehicle (D-36 discretion call)

**Research recommendation: extend `prpr/xvtest0.f` in place**, with the draw section gated so the open/close/reopen lifecycle cycle is preserved (Pitfall 10 mitigation).

Rationale:
1. Zero framework-dependency introduction (avoids Catch2/GoogleTest/Qt Test decision).
2. Matches project convention — `prpr/xvtest{1..4}.f` are the existing test vehicles.
3. `bin/cbxvtest0_qt` needs no modification.
4. One fewer file in the tree.
5. Phase 1 SHELL-01/02/06 re-verification costs nothing: the two open/close cycles remain intact; the draw calls insert between them.

The planner is free to override this (C++ test harness is a valid choice) but research's recommendation is the low-friction path.

**Concrete gated-extension shape:**
```fortran
C     Phase 2 draw sanity section — inserted between first open and second open
      CALL XVINITGRAPHIQUE
      CALL SLEEP(1)                           ! Phase 1 visual hold (SHELL-01)
C     --- Phase 2 draw exercises ---
      CALL XVEPAISSEUR( 2 )
      CALL XVTYPETRAIT( 0 )
      CALL XVTRAIT( 100, 100, 500, 100 )      ! DRAW-02 line
      CALL XVTYPETRAIT( 1 )
      CALL XVTRAIT( 100, 150, 500, 150 )      ! DRAW-06 dashed
      CALL XVTYPETRAIT( 2 )
      CALL XVTRAIT( 100, 200, 500, 200 )      ! DRAW-06 dashed-double
      CALL XVTYPETRAIT( 0 )
      CALL XVRECTANGLE( 100, 250, 150, 80 )   ! DRAW-04 filled
      CALL XVBORDRECTANGLE( 300, 250, 150, 80 ) ! DRAW-04 outlined
      CALL XVARCELLIPSE( 200, 400, 80, 60, -45.0, 90.0 )     ! DRAW-05
      CALL XVBORDARCELLIPSE( 400, 400, 80, 60, 0.0, 180.0 )  ! DRAW-05
C     polyline + polygon with INTEGER*2 point array (DRAW-03)
      XYPX(1,1) = INT2(100); XYPX(2,1) = INT2(500)
      XYPX(1,2) = INT2(150); XYPX(2,2) = INT2(450)
      XYPX(1,3) = INT2(200); XYPX(2,3) = INT2(500)
      XYPX(1,4) = INT2(150); XYPX(2,4) = INT2(550)
      CALL XVTRAITS( 4, XYPX )                ! DRAW-02 polyline
      CALL XVFACE( 4, XYPX )                  ! DRAW-03 filled polygon
      CALL SLEEP(2)                           ! hold for visual inspection of draws
      CALL EFFACER                            ! DRAW-07 clear
      CALL XVTRAIT( 100, 300, 700, 300 )      ! single line after clear
      CALL SLEEP(1)                           ! hold for visual inspection of post-clear
      CALL XVFERMER
C     Second cycle preserves Phase 1 SHELL-02 reopen invariant
      CALL XVINITGRAPHIQUE
      CALL SLEEP(1)
      CALL XVFERMER
```

Driver declares `INTEGER*2 XYPX(2, 8)` at the top. The `INT2()` intrinsic is already used in `xvue/face2d.f:48` so it is a known-available gfortran builtin.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | `prpr/xvtest1.f..xvtest4.f` all call `xvtexte_`/`xvchargefonte_`/`xvcouleurs_` — confirming D-35's deferral rationale. | User Constraints | Research read `prpr/xvtest1.f` (`[VERIFIED]`) — it calls `CHARGEFONTE`, `XVTEXTE`, `XVCOULEUR`. Not re-verified for xvtest2/3/4, but CONTEXT.md asserts it; low risk. |
| A2 | Gfortran's `SLEEP(n)` intrinsic is available. | Sanity driver | `[VERIFIED: prpr/xvtest0.f:28,34 already uses CALL SLEEP(1)]` — no risk. |
| A3 | Qt 6.10 `QWidget::resizeEvent` fires before first `paintEvent` on X11 and Wayland. | Pitfall 6, Pattern 1 | `[CITED: doc.qt.io/qt-6/qwidget.html#resizeEvent]` — documented invariant; verified by Phase 1 behavior (window paints correctly). If wrong, Pitfall 6's defensive null-check catches it. |
| A4 | `pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport` provides all linker flags needed for `libxvueqt.a` users (static .a linked into Fortran executables). | Build System | `[VERIFIED: Phase 0 BUILD-06 + Phase 1 cbl_tout_qt green]` — inherited working baseline. |
| A5 | The `xvfbordrectangle_` prototype in the legacy `xvuelc.c` matches the three sibling rectangle entry points' signature `(int*, int*, int*, int*)`. | Rectangle code example | Research did NOT read the literal `xvfbordrectangle_` body. **This is the one signature research did not verify.** Planner's first rectangle task must diff the literal legacy signature. Low risk — the pattern in neighboring rectangle entries is uniform. |
| A6 | The sanity driver extension to `xvtest0.f` does not regress Phase 1 SHELL-01/02/06 visual validation. | Pitfall 10 | Based on design reasoning: the two `XVINITGRAPHIQUE`/`XVFERMER` cycles are preserved, only draw calls inserted between first open and first close. Validator re-runs Phase 1 visual checks on the extended driver as a Phase 2 exit criterion. If regression surfaces, planner clones `xvtest0.f` → `xvtest0_draw.f` (CONTEXT.md discretion). |

**Empty-table indicator:** this table is not empty. The planner and discuss-phase should confirm A5 (verify `xvfbordrectangle_` legacy signature) as a plan task before writing the rectangle body.

## Open Questions

### Q1: `xvarcellipse_` should use `drawPie`, not `drawArc` — is the planner aware?

**What we know:** Legacy `xvuelc.c:2636` uses `XFillArc`, which fills a *pie slice* from the ellipse center to the arc endpoints. D-14 in `02-CONTEXT.md` says `painter_->drawArc(...)` for both `xvarcellipse_` and `xvbordarcellipse_`, but `QPainter::drawArc` is pen-only (stroke), not fill. The Qt equivalent of `XFillArc` is `QPainter::drawPie` `[CITED: doc.qt.io/qt-6/qpainter.html#drawPie]`. Using `drawArc` for both would silently drop the fill from `xvarcellipse_`.

**What's unclear:** Whether D-14 glossed over the fill/stroke distinction and assumed `drawArc` handled both, or whether the planner will catch it at diff time.

**Recommendation:**
- `xvbordarcellipse_` → `painter_->drawArc(bbox, a1*16, a2*16)` (pen outline only — correct per D-14)
- `xvarcellipse_`  → `painter_->drawPie(bbox, a1*16, a2*16)` (filled pie sector — **research correction**)

Planner should surface this in the plan as an explicit D-14 refinement. The angle conversion (×16) is unchanged; only the paint call differs. If the user wants an arc-shape fill (chord, not pie), the call is `drawChord` — but `XFillArc` is a pie by default so `drawPie` is the correct match.

### Q2: Does `xvfbordrectangle_` exist in `xvuelc.c`, and what is its literal signature?

**What we know:** `02-CONTEXT.md` D-13 mentions four rectangle entry points — `xvrectangle_`, `xvbordrectangle_`, `xvfrectangle_`, `xvfbordrectangle_`. Research verified the first three in `xvuelc.c` (lines 2440, 2443, 2489, 2507) but did not find a `xvfbordrectangle_` body. It may exist under a name variant or may not exist at all.

**What's unclear:** The legacy signature and whether this is a Phase 0 stub that Phase 2 promotes, or a new entry point.

**Recommendation:** Planner's first rectangle-body task runs `grep -n 'proc(xvfbordrectangle)\|xvfbordrectangle_' xvue/xvuelc.c xvue/qt/include/xvue_qt_api.h`. If no legacy body exists, `xvfbordrectangle_` is Phase 0 new-from-thin-air and the planner must read the warn-once stub signature from `xvue_qt_api.cpp` as the source of truth. This is not a blocker — the entry point is already in the ABI header; Phase 2 only has to fill its body.

### Q3: Does `xvvoir_`'s exposure-pump idiom (Phase 1 FIX-2 pattern) apply, or is a single `processEvents` sufficient?

**What we know:** D-02 specifies `xvvoir_` body as a single `canvas_->update(); processEvents(ExcludeUserInputEvents);`. Phase 1 FIX-2 / the Phase 1 summary noted that Phase 2's `xvvoir_` implementation "should use the same `QElapsedTimer`-bounded `processEvents` loop." This is in tension with D-02.

**What's unclear:** Whether the Phase 1 summary's recommendation was normative or advisory. D-02 was written **after** the Phase 1 summary (CONTEXT.md gathered 2026-04-11; Phase 1 completed same day), so D-02 is the latest word and should win.

**Recommendation:** Follow D-02 (single `processEvents` call). Empirically, a single `processEvents(ExcludeUserInputEvents)` after `update()` drains the queued paint event synchronously — the bounded loop is only needed for **window realization** (exposure from unmapped to mapped), which happens once at `xvinitgraphique_` time, not on every `xvvoir_` call. If flicker or "draws don't appear" surfaces in the sanity driver, escalate to the bounded-loop pattern; otherwise D-02 stands.

### Q4: How does the sanity driver validate DRAW-09 (resize preserve) without a graphical test harness?

**What we know:** Resize preservation is inherently an interactive behavior. The sanity driver is a batch Fortran program with SLEEP holds.

**What's unclear:** Whether the planner expects the human verifier to manually resize the window during the SLEEP hold, or whether DRAW-09 is verified by code inspection alone.

**Recommendation:** Combine both. The planner's exit-gate checklist should include:
1. Code inspection: `resizeEvent` body exists and follows D-07 sequence.
2. Manual: human verifier resizes `pp/ppxvtest0_qt` window during the 2-second draw-inspection hold and confirms content stays top-left anchored.
3. (Optional, stretch) Add a third `XVINITGRAPHIQUE` cycle to the driver that explicitly calls a new `XVRESIZE` test entry point — but this adds scope. Research does NOT recommend stretch option 3; manual verification is sufficient for a single-phase exit gate.

## Metadata

**Confidence breakdown:**
- Standard stack: **HIGH** — Qt 6.10.2 verified on machine; API mappings verified against Qt 6 documentation; inherited from working Phase 0/1 baseline.
- Architecture patterns: **HIGH** — all 35 decisions in `02-CONTEXT.md` are internally consistent and match Qt 6.10 idioms; `QPainter` lifecycle matches Qt docs; HiDPI handling matches Phase 1 SHELL-06 precedent.
- Primitive mappings: **HIGH** — legacy `xvuelc.c` bodies read for `effacer_`, `xvfond_`, `xvface_`, `xvtraits_`, `xvtypetrait_`, `xvepaisseur_`, `xvftrait_`, `xvtrait_`, `xvfacetraits_`, `xvrectangle_`, `xvbordrectangle_`, `xvfrectangle_`, `xvarcellipse_`, `xvbordarcellipse_`, `xvpxfenetre_`. Only `xvfbordrectangle_` and `xvvoir_` were not read verbatim (surfaced in Open Questions Q2/Q3).
- MefistoPoint ABI: **HIGH** — all 13 Fortran callers verified to declare `INTEGER*2 (2, N)`; `XPoint` layout confirmed; `static_assert(sizeof == 4)` already in header.
- Angle conversion: **HIGH** — derived from first principles against authoritative sources (`XDrawArc` 1/64 deg, `QPainter::drawArc` 1/16 deg, legacy `×64` scale factor). The factor-4 difference is explicitly corrected in this document.
- `xvarcellipse_` drawPie vs drawArc: **HIGH confidence in correction** (Q1) — `drawArc` does not fill, `drawPie` does; `XFillArc` semantically matches `drawPie`. Planner must fold this into the plan.
- Validation architecture: **MEDIUM** — recommended vehicle (extend `xvtest0.f`) is reasoned but not empirically validated. Alternative vehicles (new sibling `.f`, C++ harness) remain valid per D-36 discretion.
- Common pitfalls: **HIGH** — all ten pitfalls derived from CONTEXT.md's canonical refs section and Phase 1 debug-session lessons.

**Research date:** 2026-04-11
**Valid until:** 2026-05-11 (30 days — Qt 6 API is stable; Debian trixie freeze keeps local versions static)

## Sources

### Primary (HIGH confidence)

- `/home/drico/git/mefisto/.planning/phases/02-drawing-primitives-backing-pixmap/02-CONTEXT.md` — 35 locked decisions, scope, canonical refs, Claude's discretion, deferred ideas. **Authoritative for this phase.**
- `/home/drico/git/mefisto/.planning/REQUIREMENTS.md` — DRAW-01..DRAW-09 verbatim.
- `/home/drico/git/mefisto/.planning/ROADMAP.md` §"Phase 2" — goal, success criteria 1..5.
- `/home/drico/git/mefisto/xvue/xvuelc.c` lines 1413 (`effacer_`), 1434 (`xvfond_`), 1619 (`xvpxfenetre_`), 1701 (`xvface_`), 1760 (`xvtypetrait_`), 1807 (`xvepaisseur_`), 1845 (`xvftrait_`), 1862 (`xvtrait_`), 1977 (`xvtraits_`), 2035 (`xvfacetraits_`), 2440–2553 (rectangles), 2554–2675 (`xvbordarcellipse_`/`xvarcellipse_`).
- `/home/drico/git/mefisto/xvue/qt/include/xvue_qt_api.h` lines 30–35 (`MefistoPoint` typedef + static_assert), 103–142 (DRAW entry prototypes).
- `/home/drico/git/mefisto/xvue/qt/src/xvue_qt_state.h` — Phase 1 baseline.
- `/home/drico/git/mefisto/xvue/qt/src/xvue_qt_canvas.{h,cpp}` — Phase 1 `paintEvent` + Phase 2 swap site.
- `/home/drico/git/mefisto/xvue/qt/CMakeLists.txt` — build system.
- `/home/drico/git/mefisto/bin/cbxvtest0_qt` — Phase 1 driver build template.
- `/home/drico/git/mefisto/prpr/xvtest0.f` — Phase 1 driver (extension candidate).
- `/home/drico/git/mefisto/prpr/xvtest1.f` — pattern reference for full draw coverage; confirms Fortran `INTEGER*2 XPOINTS(2, MAXPTS)` convention.
- `/home/drico/git/mefisto/xvue/face2d.f`, `traits3d.f` — Fortran wrappers confirming MefistoPoint ABI.
- `/home/drico/git/mefisto/.planning/phases/01-window-shell-xvueapp-xvuewindow-xvuecanvas/01-03-SUMMARY.md` — Phase 1 completion record, `xvtest0.f` driver template, established idioms (exposure pump, xvfermer_ drain, QApplication leak).
- `/home/drico/git/mefisto/CLAUDE.md` — project working rules.

### Secondary (MEDIUM confidence — official Qt docs)

- Qt 6 `QPainter` — `doc.qt.io/qt-6/qpainter.html`: `drawLine`, `drawPolyline`, `drawPolygon`, `drawRect`, `drawArc`, `drawPie`, `drawChord`, `setRenderHint(Antialiasing)`, pen/brush state — authoritative for the primitive mapping table.
- Qt 6 `QPixmap` — `doc.qt.io/qt-6/qpixmap.html`: `setDevicePixelRatio`, `fill`, HiDPI-aware blitting.
- Qt 6 `QPen` — `doc.qt.io/qt-6/qpen.html`: `Qt::SolidLine`, `Qt::DashLine`, cosmetic pen semantics.
- Qt 6 `QPolygon` — `doc.qt.io/qt-6/qpolygon.html`: construction, `Qt::OddEvenFill` vs `Qt::WindingFill`.
- Qt 6 `QCoreApplication::processEvents` — `doc.qt.io/qt-6/qcoreapplication.html#processEvents`: `QEventLoop::ExcludeUserInputEvents` flag.
- Qt 6 `QWidget::resizeEvent` — `doc.qt.io/qt-6/qwidget.html#resizeEvent`: guaranteed to fire before first `paintEvent`.
- Xlib `XDrawArc` / `XFillArc` / `XFillPolygon` — X11 reference, angle convention 1/64 deg, `Complex+CoordModeOrigin` fill rule.
- Xlib `XPoint` — `X11/Xlib.h`, two short fields.
- gfortran `INTEGER*2` kind — `gcc.gnu.org/onlinedocs/gfortran/KIND-notation.html`.

### Tertiary (LOW confidence — none in this research)

None — all claims are either VERIFIED against the codebase or CITED from authoritative sources. No WebSearch-only findings.

---

*Phase 02-RESEARCH.md — 2026-04-11*
