# Phase 2: Drawing primitives & backing pixmap - Context

**Gathered:** 2026-04-11
**Status:** Ready for planning

<domain>
## Phase Boundary

All pure synchronous drawing primitives render into a persistent off-screen `QPixmap` via one long-lived `QPainter` owned by `XvueState`, visually matching the X11 backend on draw-only test cases. Phase 2 delivers DRAW-01 through DRAW-09: the backing pixmap lifecycle, the single `QPainter` lifetime, line / polyline / polygon / rectangle / ellipse-arc primitives, pen style and width, clear-canvas + flush, default antialiasing, and a documented resize-preserve convention.

**In scope:**
- Grow `XvueState` additively from Phase 1 — adds `QPixmap* backing_`, `QPainter* painter_`, `QColor foreground_ = Qt::white`, `QPen pen_`, `QBrush brush_`, `int pen_style_ = 0`, `int pen_width_base_ = 0`, and a private `applyPen()` helper. `background_` from Phase 1 is untouched.
- `XvueCanvas::paintEvent` body becomes the one-line swap `QPainter(this).drawPixmap(0, 0, *state_->backing_)` (the Phase 1 D-05 exit point).
- `XvueCanvas::resizeEvent` reallocates `backing_`, preserves old content via top-left sub-blit, recreates the long-lived `QPainter`, and reapplies antialiasing + pen + brush.
- Real implementations of: `xvtrait_`, `xvftrait_` (line — DRAW-02), `xvtraits_` (polyline — DRAW-02), `xvface_`, `xvfacetraits_` (filled polygon — DRAW-03), `xvrectangle_`, `xvbordrectangle_`, `xvfrectangle_`, `xvfbordrectangle_` (rectangle — DRAW-04), `xvarcellipse_`, `xvbordarcellipse_` (ellipse arc — DRAW-05), `xvtypetrait_`, `xvepaisseur_` (pen style/width — DRAW-06), `effacer_`, `xvvoir_`, `xvpxfenetre_` (clear + flush + pixel dims — DRAW-07).
- `struct MefistoPoint { short x; short y; }` ABI shim exposed in `xvue/qt/include/xvue_qt_api.h`, consumed by the three polygon/polyline entry points that historically accepted `XPoint*` (DRAW-03, Pitfall 4).
- `QPainter::Antialiasing` render hint reasserted on every new painter instance (DRAW-08).
- Documented resize-preserve convention: **top-left anchor, no scaling, no centering** (DRAW-09).
- A Phase 2 minimum draw-only sanity test (planner decides the exact driver shape per discretion note below).

**Explicitly out of scope:**
- Text / fonts / font-metric queries — Phase 3. `xvchargefonte_`, `xvtexte_`, `xvftexte_`, `xvlongtextepx_` stay on the Phase 0 warn-once stub path.
- Colormap and palette plumbing — Phase 3. `xvcouleurs_`, `xvCouleursImposees_`, `xvStockeRGBtoColormap_`, `norgb[]` logic all stay warn-once stubs. Phase 2 uses a single hardcoded `XvueState::foreground_ = Qt::white` for every pen and brush. No caller can change the color in Phase 2; this bridge is removed in Phase 3.
- Pixmap save/restore slots — Phase 4. `sauvefenetre_`, `restaurefenetre_`, `fenetremempx_`, `mempxfenetre_` stay warn-once stubs. Phase 2's `effacer_` does **not** perform the legacy `fenetremempx_` copy step.
- Event loop / mouse / keyboard delivery — Phase 5. `xvsouris_`, `xvsouris2_`, `deplsouris_` stay warn-once stubs. Phase 2's `processEvents(ExcludeUserInputEvents)` pumps never deliver input events (Pitfall 8 rationale preserved from Phase 1 D-01).
- Menus / toolbars / dock widgets — Phase 6.
- Image / GIF / PostScript export — Phase 7. **Phase 2 primitives contain zero `lasopsc` / `fpo` / `concat[]` side effects**; every `if (lasopsc > 0)` block from the legacy `xvuelc.c` bodies is dropped at translation time. `xvpostscript_` remains a warn-once stub.
- Full `prpr/xvtest1.f`–`xvtest4.f` visual parity and 5-case `testa/` A/B re-run — deferred to Phase 3. Those tests call `xvtexte_` / `xvchargefonte_` / `xvcouleurs_` which only become real in Phase 3. Running them against Phase 2 would produce uninterpretable partial output. The "matching on `prpr/xvtest1.f`–`xvtest4.f`" language in ROADMAP Phase 2 success criterion #1 is reinterpreted: Phase 2 must not *regress* anything xvtest1–4 touch visually, but pixel-level A/B on these drivers is a Phase 3 gate.

</domain>

<decisions>
## Implementation Decisions

### Paint flush strategy

- **D-01:** Every Phase 2 public drawing entry point ends with a two-line epilogue: `canvas_->update();` then `QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);`. The `update()` call schedules a repaint of the `XvueCanvas` widget; the `processEvents` call drains the queued paint event synchronously so the user sees the change before control returns to Fortran. `ExcludeUserInputEvents` is non-negotiable (Phase 1 D-01, Pitfall 6, Pitfall 8 rationale).
- **D-02:** `xvvoir_` is a real implementation in Phase 2 but its body is minimal: `canvas_->update(); QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);` — the same epilogue as every drawing primitive. This preserves the Fortran-level contract "`xvvoir_` makes everything visible" even though every drawing primitive already flushes. Rationale: legacy `xvtest*.f` callers that *do* still call `xvvoir_` (the ones not commented out) must not break, and adding a no-op `xvvoir_` costs nothing.
- **D-03:** No per-primitive guard against the case "paint pump reenters a drawing primitive". `ExcludeUserInputEvents` already prevents user callbacks from reaching the Fortran layer, and Phase 2 has no Qt-side timers or slots that call `xvue` entry points. Phase 5 revisits this when real input events start flowing.

### Painter and backing lifetime

- **D-04:** `XvueState::backing_` is a `QPixmap*` owned by `XvueState` (raw pointer, manual delete in the destructor — matches the Phase 1 style where `XvueWindow` owns `XvueCanvas` via `setCentralWidget`). `backing_` is allocated on first `XvueCanvas::resizeEvent` after construction. Phase 2 does not attempt to allocate it lazily inside drawing primitives; the widget's first `resizeEvent` fires before any `xvtrait_` call can reach the canvas (Qt 6 guarantees a `resizeEvent` before the first `paintEvent`, and `xvinitgraphique_` pumps `processEvents` on open per Phase 1 D-01).
- **D-05:** `XvueState::painter_` is a `QPainter*` owned by `XvueState`. `painter_->begin(backing_)` is called immediately after `backing_` is allocated in `resizeEvent`; `painter_->end()` is called immediately before the old `backing_` is deleted during a subsequent `resizeEvent`. The painter is alive continuously between these two moments, matching DRAW-01's "one long-lived `QPainter`".
- **D-06:** `backing_` size in device pixels is `widget.size() * devicePixelRatioF()`; the pixmap's logical DPR is set via `backing_->setDevicePixelRatio(dpr)`. `paintEvent`'s `drawPixmap(0, 0, *backing_)` is the Qt-6-idiomatic way to blit a DPR-aware pixmap — no explicit scaling needed. This keeps HiDPI correct without any special case (SHELL-06 Phase 1 invariant preserved).
- **D-07:** On `resizeEvent`, the sequence is: (a) `painter_->end()` on the old backing, (b) allocate `new_backing = new QPixmap(new_size_device_px)` with `setDevicePixelRatio(dpr)`, (c) fill the new backing with `state_->background_` via a temporary `QPainter` on it, (d) blit the old backing at logical coordinate (0, 0) via `QPainter::drawPixmap(0, 0, *old_backing)`, (e) delete the old backing, (f) move the new backing pointer into `XvueState::backing_`, (g) `painter_->begin(backing_)` on the new one, (h) re-apply `QPainter::Antialiasing` + current pen + current brush via `applyPen()` and `painter_->setBrush(brush_)`. No caller sees the transition; `xvtrait_` called before and after a resize observes the same painter API.
- **D-08:** Resize-preserve convention is **top-left anchor, no scaling, no centering, no clipping compensation**. When the widget grows, the old content sits in the top-left and new area is filled with `background_`. When the widget shrinks, the old content is clipped against the new smaller rectangle (Qt's `drawPixmap` handles this automatically). Documented in `xvue/qt/README_RESIZE.md` created by Phase 2 (one-paragraph file).

### Drawing primitive routing

- **D-09:** `xvtrait_` and `xvftrait_` are **semantically identical** in Phase 2. Both bodies are the same two lines: `state_->painter_->drawLine(*x1, *y1, *x2, *y2);` followed by the D-01 epilogue. The legacy distinction (`xvftrait` → window, `xvtrait` → mempx) is obsolete under the single-backing model: there is no separate window surface to draw to. The two symbols are kept because Fortran wrappers in `xvue/*.f` call both; merging them would break the ABI. Document the equivalence in a header comment on both entry points referencing this D-09. Phase 4's pixmap save/restore work revisits whether `xvftrait_` ever needs distinct behaviour (expected: no).
- **D-10:** `xvtraits_` (polyline, `int *nbpoints, MefistoPoint *points`) builds a stack-allocated `QPoint` array (or `std::vector<QPoint>` when `nbpoints > 128`, the threshold picked to keep the hot path allocation-free) by copying `{short x, short y}` into `{int x, int y}`, then calls `painter_->drawPolyline(qpoints, *nbpoints)`. The `MefistoPoint → QPoint` conversion happens inline in the entry-point body; no helper function is introduced in Phase 2.
- **D-11:** `xvface_` (filled polygon, `int *n, MefistoPoint *pts`) converts to a `QPolygon` (Qt-idiomatic for this case; `QPolygon` is a `QList<QPoint>` and reserves well) and calls `painter_->drawPolygon(polygon, Qt::OddEvenFill)`. `Qt::OddEvenFill` is Qt's default and matches Xlib `XFillPolygon(display, …, Complex, CoordModeOrigin)` for the geometries MEFISTO draws (non-convex polygons with no self-intersection are the common case in mesher output). The polygon is closed automatically by Qt (point 1 ↔ point n), matching the legacy comment "le pt 1 est raccorde au pt n".
- **D-12:** `xvfacetraits_` (filled polygon with explicit edge, `int *ncf, int *nca, int *n, MefistoPoint *pts`) draws the fill with `painter_->drawPolygon(polygon, Qt::OddEvenFill)` using the current brush, then draws the outline with `painter_->drawPolygon(polygon)` (no fill path — Qt's `drawPolygon` without fill flag still outlines). The `ncf` (fill color index) and `nca` (edge color index) parameters are **read but not honored** in Phase 2: both the fill brush and the pen use `state_->foreground_` (white). A `// TODO(phase 3): honor ncf/nca after palette lands` comment marks the bridge point. This satisfies the ABI but produces visually-degraded output vs. X11; the deferred `testa/` A/B re-run is the catch-up gate.
- **D-13:** `xvrectangle_` / `xvbordrectangle_` draw an outlined rectangle via `painter_->drawRect(QRect(*x, *y, *width, *height))`. `xvfrectangle_` / `xvfbordrectangle_` draw a filled rectangle the same way (Qt's `drawRect` fills with the current brush when the brush is not `Qt::NoBrush`; Phase 2 `applyPen()` ensures the brush is `QBrush(foreground_)`). The four entry points collapse to two bodies internally but stay four separate exported symbols (same rationale as D-09: ABI preservation). For `xvfrectangle_` / `xvfbordrectangle_`, Phase 2 behaviour is "fill with foreground color, no separate border color" — the border-color parameter (when present in the legacy signature) is deferred to Phase 3 alongside the palette.
- **D-14:** `xvarcellipse_` / `xvbordarcellipse_` call `painter_->drawArc(QRect(*x, *y, *width, *height), start_16, span_16)` where `start_16 = *startAngle * 16` and `span_16 = *arcAngle * 16`. The ×16 factor matches Qt 6's angle convention (1/16th of a degree per unit). The legacy X11 `XDrawArc` uses the same 1/64th of a degree convention that `xvuelc.c` already accounts for — Phase 2 must check `xvuelc.c` at the entry-point site and keep the same angular domain at the Fortran ABI. The planner's first task on these two entry points is a literal side-by-side diff of the existing `xvuelc.c` bodies against the Qt body draft, confirming angle conversion is preserved.
- **D-15:** `effacer_` body is three lines: `state_->painter_->fillRect(state_->backing_->rect(), state_->background_);` then the D-01 epilogue (`canvas_->update()` + `processEvents`). The legacy `XClearWindow` + `XFlush` + `fenetremempx` sequence is replaced by a single fill because the Qt backing *is* the window surface under `drawPixmap`. No `lasopsc 101/102` PostScript side effect (dropped per D-26).

### Pen, brush, and state propagation

- **D-16:** `XvueState::applyPen()` is a private method that rebuilds `pen_` from `pen_style_`, `pen_width_base_`, and `foreground_`, then calls `painter_->setPen(pen_)`. It is called by `xvtypetrait_`, `xvepaisseur_`, and the `resizeEvent` painter recreation path. Rebuilding on every change (rather than mutating `pen_` in place) keeps the code straightforward and avoids subtle bugs from partial updates. The body is:
  ```cpp
  pen_.setColor(foreground_);
  pen_.setWidth(style_width_);   // see D-17
  pen_.setStyle(qt_style_);      // see D-17
  painter_->setPen(pen_);
  ```
- **D-17:** Pen style mapping (`xvtypetrait_`) preserves the X11 "type 2 = double width" semantic:
  - `type 0` (solid): `qt_style_ = Qt::SolidLine`, `style_width_ = pen_width_base_`
  - `type 1` (dashed): `qt_style_ = Qt::DashLine`, `style_width_ = pen_width_base_`
  - `type 2` (dashed, double width): `qt_style_ = Qt::DashLine`, `style_width_ = max(1, pen_width_base_ * 2)`
  - anything else: fall through to type 2 (matches the legacy `default:` branch at `xvuelc.c:1782`)
  Storing `pen_style_` as the raw Fortran int (0/1/2) in `XvueState` and deriving `qt_style_` / `style_width_` lazily in `applyPen()` is the single-source-of-truth choice; mutations to either `pen_style_` or `pen_width_base_` re-run the full derivation.
- **D-18:** `xvepaisseur_` body stores `*pepais` into `pen_width_base_` then calls `applyPen()`. Qt's `width=0` means "cosmetic pen, always 1 device pixel", which matches X11's `line_width=0` "thin line" semantic; no special case. No validation / clamping on the input — if Fortran sends a negative width, Qt's `setWidth` accepts it and clamps internally. Not worth guarding Phase 2.
- **D-19:** `xvtypetrait_` body stores `*ptype` into `pen_style_` then calls `applyPen()`. No side effects beyond that (the legacy `lasopsc` PS recording block is dropped).
- **D-20:** `XvueState::brush_` is `QBrush(foreground_, Qt::SolidPattern)`, installed once at construction and kept in sync with `foreground_` whenever Phase 3 eventually changes it. In Phase 2 `foreground_` is hardcoded to `Qt::white` so `brush_` is effectively constant. `painter_->setBrush(brush_)` is called from `applyPen()` for symmetry so that both pen and brush are always synchronized after any state-changing entry point or resize.
- **D-21:** `XvueState::foreground_` is introduced in Phase 2 with the hardcoded value `Qt::white` and is **not** mutable from any Phase 2 entry point. Phase 3's palette work introduces `xvcouleurs_` / `xvCouleursImposees_` as real implementations that write to `foreground_` and call `applyPen()`. Phase 2 callers that expect different colors (e.g. the legacy mesher colors mesh elements by region) produce white-only output. This is an explicit, documented visual degradation for one phase; the Phase 3 A/B gate is where it is reconciled.
- **D-22:** `QPainter::setRenderHint(QPainter::Antialiasing, true)` is called once immediately after every `painter_->begin(backing_)` call (once at initial `resizeEvent`, once per subsequent `resizeEvent`). This satisfies DRAW-08 without per-primitive reapplication. It does not need to be toggled off anywhere in Phase 2.

### `effacer_`, `xvvoir_`, `xvpxfenetre_`, `xvfond_` interaction

- **D-23:** `xvpxfenetre_(int *x, int *y)` returns the widget's logical-pixel size: `*x = canvas_->width(); *y = canvas_->height();`. It reads the widget (not the backing pixmap) so that Fortran callers see widget coordinates consistent with what they pass to drawing primitives. The legacy body at `xvuelc.c:1619` returns `XWindowAttributes` width/height — same semantic.
- **D-24:** `xvfond_` (Phase 1 D-14) is extended in Phase 2 to also `state_->painter_->fillRect(state_->backing_->rect(), state_->background_)` after updating `state_->background_`, so that changing the background color takes effect on already-drawn content by repainting from scratch. This matches the legacy semantic "setting the background color clears the window to that color" inherited from `XSetWindowBackground` + `XClearWindow`. The fill is followed by the D-01 epilogue. Phase 1's `xvfond_` only stored the color and called `update()`; Phase 2 replaces the body with the full fill-and-flush path. Phase 1 D-15 "canvas observes new value on next paint" is preserved as a fallback but no longer the primary mechanism.
- **D-25:** Drawing a primitive never touches `state_->background_`. `effacer_` is the only Phase 2 path that re-fills the backing with `background_`.

### PostScript hooks

- **D-26:** **Every `if (lasopsc > 0)` block from every legacy primitive is dropped at translation time.** Phase 2 `xvtrait_`, `xvface_`, `xvtypetrait_`, `xvepaisseur_`, `xvrectangle_` (and siblings), `xvarcellipse_` (and sibling), `effacer_`, `xvfond_` contain zero references to `lasopsc`, `fpo`, `concat[]`, `buf[]`, `courgb[]`, `counb`, `iFa`, `iep`, `ity`, or any of the PostScript-recording globals from `xvuelc.c`. The Qt bridge has no `FILE*` for PS output in Phase 2.
- **D-27:** `xvpostscript_` remains a warn-once stub in Phase 2 (unchanged from Phase 0). Phase 7 reintroduces it as a real implementation using a dedicated `QPainter` on a `QPdfWriter` or a direct-emit approach (Phase 7 decides); it will **not** share the drawing-primitive code path. The recorder-inside-primitives pattern that `xvuelc.c` uses is explicitly rejected for the Qt backend.
- **D-28:** The legacy variables `lasopsc`, `norgb[]`, `courgb[]`, `counb`, `fpo`, `chaine[]`, `concat[]`, `nbrcon`, `ypixels`, `xinic`, `yinic`, `iFa`, `iep`, `ity` are **not** ported to the Qt `XvueState` struct or anywhere in `xvue/qt/`. They are X11-backend artifacts that do not survive the migration. If Phase 7 needs state for PS recording, it introduces its own dedicated struct in its own header; it does not reuse these names.

### `MefistoPoint` ABI shim and Fortran call-site audit

- **D-29:** `struct MefistoPoint { short x; short y; };` is declared in `xvue/qt/include/xvue_qt_api.h` inside the existing `extern "C"` block, immediately after the `#define proc(x) x##_` macro and before the first Fortran-facing prototype. The declaration is a plain C struct (no C++ features) so the header stays compilable by both C (legacy Fortran-side consumers) and C++ (the Qt bridge implementation). A `static_assert(sizeof(struct MefistoPoint) == 4, …)` is added in the `xvue_qt_api.cpp` translation unit to catch layout drift at build time.
- **D-30:** The three entry points that historically took `XPoint*` are declared with `MefistoPoint*`: `void xvface_(int *n, struct MefistoPoint *pts);`, `void xvtraits_(int *nbpoints, struct MefistoPoint *points);`, `void xvfacetraits_(int *ncf, int *nca, int *n, struct MefistoPoint *pts);`. The C++ bodies convert inline to `QPoint`/`QPolygon` (D-10, D-11, D-12). No `#include <X11/Xlib.h>` appears anywhere under `xvue/qt/`.
- **D-31:** A Fortran call-site audit is a mandatory plan task: grep `xvue/*.f`, `prpr/*.f`, `mail/`, `elas/`, `flui/`, `ther/`, `nlse/`, `util/`, `reso/` for `CALL XVFACE`, `CALL XVTRAITS`, `CALL XVFACETRAITS` (case-insensitive), confirm every caller declares its point array as `INTEGER*2 … (2, N)` or compatible layout, and record the findings in a new `xvue/qt/MEFISTO_POINT_AUDIT.md` (one file, list of caller locations + declared type). Any caller that uses a non-short layout is a hard block — Phase 2 cannot ship until it is reconciled. This is the concrete form of Pitfall 4's "Audit every Fortran call site for these three functions once during Phase 2" recommendation.

### Unchanged from Phase 1

- **D-32:** Every Phase 2 entry point body starts with `XVUE_QT_ASSERT_MAIN_THREAD()` (Phase 1 D-12/D-13 pattern). No exceptions. The `xvue_qt_api.cpp` bulk-insert from Phase 1 already covers the stubs; Phase 2 adds the assertion to each new real body.
- **D-33:** Every Phase 2 entry point signature is copied **literally** from `xvue/xvuelc.c` with only the `XPoint*` → `struct MefistoPoint*` swap from D-30 and the `#define proc(x) x##_` macro invocation. No scalar-by-value "cleanup" (Pitfall 2). No hidden trailing length argument (Pitfall 3 — not relevant in Phase 2 anyway since no string arguments land yet).
- **D-34:** The `verify_abi` CMake custom target from Phase 0 D-12 continues to run — adding real bodies to seven dozen entry points does not change the trailing-underscore symbol set, so the target passes without modification. The `verify_no_exec` target from Phase 1 D-10 continues to run unchanged. Both are required for the Phase 2 build to succeed.

### Phase 2 validation minimum (planner discretion, see Claude's Discretion)

- **D-35:** Full `prpr/xvtest1.f`–`xvtest4.f` parity and the 5-case `testa/` A/B re-run are **deferred to Phase 3**. Rationale: each of those drivers calls `xvtexte_` / `xvchargefonte_` / `xvcouleurs_` at least once; Phase 2 leaves those as warn-once stubs; the output would be uninterpretable. Phase 3 is the natural gate because it delivers the last piece (fonts + palette) that unblocks visual A/B comparison.
- **D-36:** Phase 2 still needs a **draw-only** sanity driver that exercises the new primitives in isolation. The exact form — extending `prpr/xvtest0.f` with a handful of draw calls vs. adding a new `prpr/xvtest0_draw.f` sibling vs. a dedicated C++ unit-test binary — is left to the planner (Claude's Discretion below). The minimum coverage the driver must exercise: one line (`xvtrait_`), one polyline with ≥ 3 points (`xvtraits_`), one filled polygon with ≥ 4 points (`xvface_`), one outlined rectangle (`xvrectangle_`), one filled rectangle (`xvfrectangle_`), one arc (`xvarcellipse_`), both pen styles (`xvtypetrait_` with 0, 1, 2), and one `effacer_` between two draw sequences. `xvfond_` is exercised implicitly through the background color set at open.

### Claude's Discretion

- **Phase 2 sanity-driver form** — extend `prpr/xvtest0.f`, introduce `prpr/xvtest0_draw.f`, or add a C++ unit-test harness under `xvue/qt/tests/`. The D-36 coverage list is the contract; the vehicle is the planner's call. A Fortran extension of `xvtest0.f` is the lowest-friction option and matches the project's existing test convention; a C++ unit test would exercise the same code paths without needing the Fortran link step but introduces a test-framework dependency that is otherwise absent. No wrong answer, pick one and commit.
- **`bin/cbxvtest*_qt` script variant(s)** — Phase 1 D-20 added `bin/cbxvtest0_qt` for the shell test; Phase 2 either reuses that script (if extending `xvtest0.f`) or clones it for a new driver. Either way, `bin/cbl_tout_qt` is not modified (Phase 1 D-20 scope anchor preserved).
- **`XvueState` internal file split** — whether `XvueState::applyPen()` lives in `xvue_qt_state.cpp` (matching Phase 1's file layout) or is inlined in the header is a planner call. The public API contract (the set of methods and fields) is fixed by D-04/D-05/D-16/D-20.
- **Inline helper for `MefistoPoint → QPoint` conversion** — D-10 says "no helper function in Phase 2", but if the same six lines end up copied across three entry points the planner may extract a `static inline` helper in the same translation unit. Do **not** export it from a header; it is strictly a local convenience.
- **Exact fill/stroke order in `xvfacetraits_`** — D-12 specifies fill first then outline. If visual comparison against the legacy X11 output shows the opposite order produces a closer match, flip it. The 1-pixel antialiasing halo on the edge is where the two orders differ.
- **`drawPolygon` fill rule selection** — D-11 locks `Qt::OddEvenFill`. If the mesher ever hands the polygon entry points a self-intersecting polygon (unlikely but not impossible for user-drawn 2D regions), `Qt::WindingFill` may produce a closer X11 match. The planner may switch the default if an A/B case surfaces in Phase 3; Phase 2 ships with OddEvenFill and a one-line comment noting the alternative.
- **Allocating the initial `backing_` inside `XvueCanvas` constructor vs first `resizeEvent`** — D-04 specifies first `resizeEvent`. If the planner finds that `XvueCanvas` in Phase 1 already receives a `resizeEvent` before the first `paintEvent` (Qt 6 does this on every platform), the constructor-based allocation is unnecessary. If not, a constructor-based allocation with a sentinel 1×1 `QPixmap` is acceptable and gets replaced on the first real `resizeEvent`. Either path converges on the same steady state.
- **`xvue/qt/README_RESIZE.md` location and form** — D-08 mandates a one-paragraph file documenting the top-left sub-blit convention. Whether it lives at `xvue/qt/README_RESIZE.md`, gets folded into an existing `xvue/qt/README_COORDS.md` (Phase 0 D-03), or becomes a comment block in `xvue_qt_state.h` is at the planner's discretion. A separate file is preferred for grep-ability.
- **Thread of future Phase 3 color plumbing** — `applyPen()` reads `foreground_`. Phase 2 never writes `foreground_`. Phase 3 will. The planner should not add any Phase-3 stub that writes it; a clean Phase 3 diff produces `xvcouleurs_` at that time. Adding placeholder `// TODO(phase 3)` comments at the single entry point `XvueState::foreground_` declaration is sufficient.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing Phase 2.**

### Phase scope anchors
- `.planning/ROADMAP.md` §"Phase 2: Drawing primitives & backing pixmap" — phase boundary, goal, depends-on (Phase 1), success criteria 1–5
- `.planning/ROADMAP.md` §"Merge Phase 3 (text/colormap) into Phase 2 — rejected" and §"Merge Phase 4 (pixmap save/restore) into Phase 2 — rejected" — explicit scope-boundary precedents that inform D-35's deferral
- `.planning/REQUIREMENTS.md` §"Draw — drawing primitives and backing pixmap" DRAW-01 through DRAW-09 — the 9 requirements this phase delivers
- `.planning/phases/00-build-skeleton-abi-stubs/00-CONTEXT.md` — Phase 0 decisions Phase 2 inherits: D-01 (xvue/qt/ layout), D-02 (read-only scope), D-04 (single ABI header — Phase 2 grows it with `MefistoPoint` per D-29), D-05 (trailing-underscore macro), D-06 (pointer arg convention — non-negotiable per D-33), D-12 (verify_abi custom target), D-17 (warn-once stub pattern for every entry point Phase 2 does not implement)
- `.planning/phases/01-window-shell-xvueapp-xvuewindow-xvuecanvas/01-CONTEXT.md` — Phase 1 decisions Phase 2 extends: D-01 (`processEvents(ExcludeUserInputEvents)` policy, inherited by every Phase 2 primitive epilogue per D-01), D-04 (`XvueState` additive growth — Phase 2 adds `backing_`, `painter_`, `foreground_`, `pen_`, `brush_`, `pen_style_`, `pen_width_base_`), D-05 (`XvueCanvas::paintEvent` one-line swap — Phase 2 actually performs the swap per D-06 of this phase), D-12/D-13 (`XVUE_QT_ASSERT_MAIN_THREAD()` retrofit — every Phase 2 body starts with it per D-32), D-14/D-15 (`xvfond_` storage — Phase 2 extends the body per D-24), D-20 (`bin/cbxvtest0_qt` script template)

### Research synthesis
- `.planning/research/ARCHITECTURE.md` §"Components" — `XvueApp`/`XvueWindow`/`XvueCanvas`/`XvueState` split; Phase 2 grows `XvueState` additively within this split (D-04, D-05)
- `.planning/research/ARCHITECTURE.md` §"Recommended Project Structure" — per-component header/source file layout the planner inherits from Phase 1
- `.planning/research/PITFALLS.md` §"Pitfall 2: Fortran passes everything by reference" — drives D-33 "signatures copied literally from xvuelc.c"
- `.planning/research/PITFALLS.md` §"Pitfall 3: Fortran character strings" — not directly relevant in Phase 2 (no string entry points) but the general "copy signatures literally" principle applies
- `.planning/research/PITFALLS.md` §"Pitfall 4: `XPoint*` is a short/short struct" — drives D-29, D-30, D-31 (`MefistoPoint` shim, header placement, Fortran call-site audit)
- `.planning/research/PITFALLS.md` §"Pitfall 6: `QApplication::exec()` top-level" — drives D-01, D-02 (per-primitive `processEvents(ExcludeUserInputEvents)` pump, no `exec()`)
- `.planning/research/PITFALLS.md` §"Pitfall 7: `processEvents` starvation" — informs D-01's choice of `ExcludeUserInputEvents` flag
- `.planning/research/PITFALLS.md` §"Pitfall 8: Modal `QDialog::exec()` re-entrancy" — drives D-03's stance on re-entrancy guards (rely on `ExcludeUserInputEvents` + Phase 2 having no slots)
- `.planning/research/PITFALLS.md` §"Pitfall 9: PIE / `-fPIC`" — not directly in Phase 2 but the CMake flag discipline inherited from Phase 0 depends on it
- `.planning/research/PITFALLS.md` §"Pitfall 14: Y-axis direction" — relevant to D-14 (`xvarcellipse_` angle conversion): Qt and X11 both use Y-down, so no inversion is needed at the primitive level, but the angle sign convention still differs between `XDrawArc` (64ths of a degree) and `QPainter::drawArc` (16ths of a degree)
- `.planning/research/PITFALLS.md` §"Pitfall 15: HiDPI pixel reporting" — drives D-06 (`devicePixelRatioF()` backing sizing) and D-23 (`xvpxfenetre_` returns logical pixels)
- `.planning/research/PITFALLS.md` §"Pitfall 16: `processEvents` inside `resizeEvent`" — by default Phase 2 does **not** call `processEvents` inside `resizeEvent`; the D-07 sequence is purely local (end painter, allocate new, blit, begin new, reapply state). The per-primitive epilogue is the only place `processEvents` runs

### Codebase maps
- `.planning/codebase/ARCHITECTURE.md` — shared-data layer (read-only for Phase 2) and the launcher/entry-point/solver-module/utility/graphics layering that `xvue/qt/` plugs into
- `.planning/codebase/STRUCTURE.md` — `xvue/` and `prpr/` directory layouts (relevant to D-36 sanity driver location)
- `.planning/codebase/STACK.md` — gfortran/gcc versions, Qt 6.10 package names on Debian trixie (unchanged from Phase 0)
- `.planning/codebase/CONVENTIONS.md` — Fortran 77 fixed-form column discipline (relevant if Phase 2 adds a new `prpr/xvtest0_draw.f`), `bin/cb*` shell-script conventions, `#define proc(x) x##_` idiom
- `.planning/codebase/CONCERNS.md` — stale `.o` fragility (drives the clean-build discipline inherited from Phase 0)
- `.planning/codebase/INTEGRATIONS.md` — Qt 6 CMake integration, `find_package(Qt6 COMPONENTS Widgets Gui Core)` precedents
- `.planning/codebase/TESTING.md` — manual `testa/`/`testf/` workflow; Phase 2 defers the A/B resume per D-35

### Direct source — read these literally
- `/home/drico/git/mefisto/xvue/xvuelc.c:1413-1432` — legacy `effacer_` body. D-15 replaces the semantic with a single `fillRect` on the backing; the `fenetremempx_` copy and the `lasopsc` PS block are dropped.
- `/home/drico/git/mefisto/xvue/xvuelc.c:1434-1460` — legacy `xvfond_` body. D-24 extends the Phase 1 implementation to re-fill the backing after updating `background_`.
- `/home/drico/git/mefisto/xvue/xvuelc.c:1619-1650` (approx.) — legacy `xvpxfenetre_` body. D-23 replaces `XWindowAttributes`-based width/height with `canvas_->width()` / `canvas_->height()`.
- `/home/drico/git/mefisto/xvue/xvuelc.c:1701-1757` — legacy `xvface_` body. D-11 preserves the `Complex+CoordModeOrigin` semantic via `Qt::OddEvenFill`; the full `lasopsc` PS emission block (lines ~1720–1755) is dropped per D-26.
- `/home/drico/git/mefisto/xvue/xvuelc.c:1760-1805` — legacy `xvtypetrait_` body. D-17 preserves the three-style semantic including type-2 "double width" at the Qt level; the `lasopsc` PS block is dropped.
- `/home/drico/git/mefisto/xvue/xvuelc.c:1807-1843` — legacy `xvepaisseur_` body. D-18 preserves the semantic; the `lasopsc` PS block is dropped.
- `/home/drico/git/mefisto/xvue/xvuelc.c:1845-1860` — legacy `xvftrait_` body (draws to `fenetre_mef` directly). D-09 makes it semantically equivalent to `xvtrait_` in the single-backing Qt model.
- `/home/drico/git/mefisto/xvue/xvuelc.c:1862-1976` — legacy `xvtrait_` body and full `lasopsc` PS-emission tail. D-01 and D-26 reduce this to `painter_->drawLine(...)` plus the two-line epilogue.
- `/home/drico/git/mefisto/xvue/xvuelc.c:1977-2034` — legacy `xvtraits_` body (polyline). D-10 implementation reference.
- `/home/drico/git/mefisto/xvue/xvuelc.c:2035-2090` (approx.) — legacy `xvfacetraits_` body (filled polygon with edge). D-12 implementation reference.
- `/home/drico/git/mefisto/xvue/xvuelc.c:2384-2442` (approx.) — legacy `xvvoir_` body. D-02 reduces this to `canvas_->update()` + `processEvents`.
- `/home/drico/git/mefisto/xvue/xvuelc.c:2443-2488` — legacy `xvbordrectangle_` body. D-13 implementation reference.
- `/home/drico/git/mefisto/xvue/xvuelc.c:2489-2506` — legacy `xvfrectangle_` body. D-13 implementation reference.
- `/home/drico/git/mefisto/xvue/xvuelc.c:2507-2553` — legacy `xvrectangle_` body. D-13 implementation reference.
- `/home/drico/git/mefisto/xvue/xvuelc.c:2554-2615` — legacy `xvbordarcellipse_` body. D-14 implementation reference (angle convention check).
- `/home/drico/git/mefisto/xvue/xvuelc.c:2616-2680` (approx.) — legacy `xvarcellipse_` body. D-14 implementation reference (angle convention check).
- `/home/drico/git/mefisto/xvue/xvuelc.c:60-70` (approx.) — `#define proc(x) x##_` macro that Phase 2 reuses unchanged via `xvue/qt/include/xvue_qt_api.h`.
- `/home/drico/git/mefisto/xvue/qt/include/xvue_qt_api.h` — Phase 0 ABI header. Phase 2 adds `struct MefistoPoint` (D-29) and changes three prototypes (D-30).
- `/home/drico/git/mefisto/xvue/qt/src/xvue_qt_api.cpp` — Phase 0/1 stub file. Phase 2 fills in real bodies for the DRAW entry points and leaves all other stubs on the warn-once path.
- `/home/drico/git/mefisto/xvue/README_COORDS.md` (Phase 0 D-03) — Y-axis convention audit. Phase 2 is the first phase where primitives actually consume pixel coordinates; every Phase 2 drawing PR reads this file once and confirms the convention is honored (no double-inversion, no off-by-one).
- `/home/drico/git/mefisto/prpr/xvtest0.f` — Phase 1 D-19 sanity driver. Phase 2 either extends it or adds a sibling per D-36.
- `/home/drico/git/mefisto/prpr/xvtest1.f`, `xvtest2.f`, `xvtest3.f`, `xvtest4.f` — the drivers ROADMAP Phase 2 nominally targets. D-35 defers full parity to Phase 3; planner may still read them for the expected Fortran call patterns around drawing primitives.
- `/home/drico/git/mefisto/bin/cbxvtest0_qt` (Phase 1 D-20) — template for any new Phase 2 sanity-driver build script.
- `/home/drico/git/mefisto/CLAUDE.md` — working rules: never break the X11 build, run the smallest relevant test after every change, commit after every logical step, ask before installing system packages.

### External references
- Qt 6 `QPainter` documentation — authoritative source for `drawLine`, `drawPolyline`, `drawPolygon`, `drawRect`, `drawArc`, `setRenderHint(Antialiasing)`, pen/brush state (`doc.qt.io/qt-6/qpainter.html`)
- Qt 6 `QPixmap` documentation — authoritative source for `setDevicePixelRatio` and HiDPI-aware blitting (`doc.qt.io/qt-6/qpixmap.html`)
- Qt 6 `QPen` documentation — authoritative source for `Qt::SolidLine` / `Qt::DashLine` / width semantics, cosmetic pen (width 0) behavior (`doc.qt.io/qt-6/qpen.html`)
- Qt 6 `QPolygon` documentation — authoritative source for polygon construction and `Qt::OddEvenFill` vs `Qt::WindingFill` (`doc.qt.io/qt-6/qpolygon.html`)
- Qt 6 `QCoreApplication::processEvents` documentation — authoritative source for `QEventLoop::ExcludeUserInputEvents` flag semantics (`doc.qt.io/qt-6/qcoreapplication.html#processEvents`)
- Xlib `XDrawArc` / `XFillPolygon` documentation — legacy reference for angle convention (64ths of a degree) and fill rule interpretation, used to verify D-14 angular conversion is correct
- Xlib `XPoint` definition (X11/Xlib.h) — two `short` fields, the layout Phase 2's `MefistoPoint` must match byte-for-byte per D-29

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **`xvue/qt/include/xvue_qt_api.h`** (Phase 0 D-04) — the single ABI header. Phase 2 grows it with `struct MefistoPoint` (D-29) and modifies three prototypes to take `MefistoPoint*` (D-30). No new headers are added to the Fortran-facing ABI surface.
- **`xvue/qt/src/xvue_qt_api.cpp`** (Phase 0 D-17, Phase 1 retrofit) — the single stub file. Phase 2 replaces the warn-once bodies of the DRAW entry points with real implementations; every other stub stays on the warn-once path untouched.
- **`XvueState` struct** (Phase 1 D-04) — grows additively. Phase 2 adds `QPixmap* backing_`, `QPainter* painter_`, `QColor foreground_ = Qt::white`, `QPen pen_`, `QBrush brush_`, `int pen_style_ = 0`, `int pen_width_base_ = 0`, and a private `applyPen()` method. The existing `background_` field (Phase 1 D-04) is untouched. Phase 3 will add font/palette fields; Phase 4 will add pixmap-slot fields.
- **`XvueCanvas::paintEvent`** (Phase 1 D-05) — the single-line swap point. Phase 2 performs the swap: the body becomes `QPainter(this).drawPixmap(0, 0, *state_->backing_);`. The warn-once fallback (`fillRect(rect(), state_->background_)`) is removed in the same commit.
- **`XvueCanvas::resizeEvent`** (new in Phase 2) — the single surface where backing allocation, painter recreation, content preservation, and antialiasing reapplication converge. D-07 is its exact contract.
- **`XVUE_QT_ASSERT_MAIN_THREAD()` macro** (Phase 1 D-12/D-13) — retrofitted into every Phase 2 real entry point body as its first line. Reused verbatim.
- **Phase 1 SHELL-04 screen metric accessors** — `xvpxecran_`, `xvmmecran_`, `xvinfo_` are unchanged by Phase 2. `xvpxfenetre_` (Phase 2 D-23) follows the same pattern (logical pixels via a Qt widget method).
- **`xvfond_`** (Phase 1 D-14/D-15) — Phase 2 D-24 extends the body with a full-backing re-fill. The Phase 1 storage and `canvas_->update()` paths remain as the underlying mechanism.
- **`#define proc(x) x##_` macro** (Phase 0 D-05, reused from `xvuelc.c`) — every Phase 2 entry-point body is declared via this macro. No new macro is introduced.
- **`bin/cbxvtest0_qt`** (Phase 1 D-20) — reused or cloned for the Phase 2 sanity-driver build (per D-36 and Claude's Discretion).
- **`verify_abi` and `verify_no_exec` CMake custom targets** (Phase 0 D-12, Phase 1 D-10) — continue to run unchanged. D-34 documents that they pass without modification.

### Established Patterns
- **Warn-once stub discipline** (Phase 0 D-17) — every entry point Phase 2 does not implement (text, fonts, colormap, save/restore slots, mouse/keyboard, menus, export) continues to print `"xvue-qt: stub <fn> not implemented yet"` on first call. Phase 2 introduces zero regressions in this area.
- **Single-file ABI surface** (Phase 0 D-04, Phase 1 preserved) — Phase 2 keeps `xvue_qt_api.h` as the only Fortran-facing header. The `MefistoPoint` struct lives inside it, not in a separate header.
- **Additive `XvueState` growth** (Phase 1 D-04 pattern) — Phase 2 appends fields; it does not rename or remove anything from Phase 1. The `background_` default `Qt::black` and `foreground_` default `Qt::white` match the legacy X11 `BlackPixel`/`WhitePixel` defaults.
- **Signatures copied literally from `xvuelc.c`** (Phase 0 D-06, Pitfalls 2 & 3) — every Phase 2 body's signature is a literal copy of the legacy signature, with the single `XPoint*` → `struct MefistoPoint*` swap per D-30. No simplification, no hidden length arguments, no pass-by-value.
- **Legacy-side-effect stripping at translation time** (new in Phase 2) — D-26, D-27, D-28 codify the pattern: when porting a legacy body, the `lasopsc` / `fpo` / `concat` / PS-recording logic is deleted, not preserved. Each port becomes a ~3-line body instead of a 50-line body. Phase 7 reintroduces PS as a dedicated code path with no shared state.
- **CMake-enforced invariants** (Phase 0 D-12, Phase 1 D-10) — the `verify_abi` and `verify_no_exec` targets gate every Phase 2 build. No new build-time invariants are introduced.
- **Fortran call-site audit as a first-class task** (new in Phase 2, applying Pitfall 4) — D-31 turns the audit into a concrete deliverable (`xvue/qt/MEFISTO_POINT_AUDIT.md`) that can be reviewed and re-run.
- **`XVUE_QT_ASSERT_MAIN_THREAD()` at the top of every body** (Phase 1 D-12/D-13) — every new real body added in Phase 2 inherits this pattern without exception.

### Integration Points
- **`XvueCanvas::paintEvent`** — the Phase-1-declared swap point. Phase 2 flips it from `fillRect(rect(), background_)` to `drawPixmap(0, 0, *backing_)`. Phase 3+ paintEvent bodies remain this single line; state changes happen on the backing, not in the paint path.
- **`XvueCanvas::resizeEvent`** — new in Phase 2. Solely responsible for backing/painter lifecycle. Phase 4's pixmap slot feature must not bypass this path; slot restore operations go through `painter_->drawPixmap(...)` on the existing backing.
- **`XvueState` field growth** — Phase 2 adds the drawing-state fields; Phase 3 adds font + palette (`foreground_` becomes mutable via `xvcouleurs_`); Phase 4 adds the pixmap-slot vector; Phase 5 adds input-state fields (pending click, modifier mask, etc.). Every phase grows additively without renaming.
- **Fortran `xvue/*.f` wrappers** — untouched by Phase 2 per the project invariant. Phase 2 must not modify anything under `xvue/*.f`; the `MefistoPoint` struct layout must match whatever the wrappers already pass (D-31 audit confirms this before any code lands).
- **`xvfond_`** — the single point where `background_` is written; extended by Phase 2 D-24 to trigger a full re-fill. Every other code path reads `background_` and never writes it.
- **Phase 3 bridge point: `XvueState::foreground_`** — Phase 2 declares and initializes it; Phase 3 writes it from `xvcouleurs_` and friends. No other Phase 2 code mutates it.
- **Phase 7 bridge point: nothing** — Phase 2 ports no PS hooks into the Qt backend (D-26/D-27/D-28). Phase 7 starts from a blank slate; there is no `lasopsc` or recording state for it to consume. This is the explicit anti-coupling decision.
- **`bin/cb*_qt` linker lines** — unchanged from Phase 0. Phase 2 adds no new libraries and no new compiler flags. `libxvueqt.a` grows by the new `.o` files for the real drawing bodies.
- **CMake `xvue/qt/CMakeLists.txt`** — grows by the new source files the planner adds (if any new `.cpp` files are introduced beyond `xvue_qt_api.cpp` and the Phase 1 component files). No changes to `find_package`, `target_link_libraries`, or the existing custom targets.

</code_context>

<specifics>
## Specific Ideas

- **"Modern Qt-idiomatic, not a 1:1 port"** — the defining directive for Phase 2. Every gray area resolved toward the simpler Qt approach when the X11 semantic had legacy baggage (PS side effects, dual-surface drawing, `fenetremempx` coupling). Where the legacy semantic carried real visual intent (pen style #2 "double width", top-left resize anchor, antialiasing-as-free-improvement), it is preserved.
- **"Signatures copied literally, bodies rewritten from scratch"** — the translation pattern. Phase 2 never edits `xvue/xvuelc.c`. It reads the legacy body, extracts the ≤5 meaningful Xlib calls, and rewrites a clean Qt body. The rest (PS recording, postscript flush, dual-surface coordination) is discarded.
- **"Ship Phase 2 without palette, explicitly degraded"** — D-21 locks `foreground_ = Qt::white` as a deliberate short-lived bridge. Color-coded output looks wrong for one phase; this is documented, known, and reconciled by Phase 3's A/B gate. Attempting to shim color in Phase 2 would leak palette state across the phase boundary.
- **"The PS-hook stripping is load-bearing for Phase 7"** — D-26/D-27/D-28 aren't just simplification; they lock in that the Qt backend's PS path is independent of the drawing primitive path. Phase 7 gets to design PS output as a real Qt feature (`QPdfWriter`, `QPrinter::setOutputFormat(QPrinter::PdfFormat)`, or direct emit) without inheriting two decades of `lasopsc` state-machine coupling.
- **"Fortran call-site audit is a deliverable, not a checklist item"** — D-31 produces `xvue/qt/MEFISTO_POINT_AUDIT.md`. This is explicit because Pitfall 4 is the one place where a silent ABI mismatch can drift across the migration; a grep-and-document task with a persistent artifact makes re-running the check trivial when Phase 4+ touches the pixmap slots or Phase 5 touches mouse-event coordinates.
- **"Test strategy deferred"** — user-directed. D-35 pushes full parity to Phase 3; D-36 keeps a draw-only sanity minimum as the Phase 2 exit gate and hands the exact form to the planner.

</specifics>

<deferred>
## Deferred Ideas

- **Full `prpr/xvtest1.f`–`xvtest4.f` visual parity against X11** — deferred to Phase 3 (D-35). These drivers call text and color entry points that remain warn-once stubs in Phase 2.
- **5-case `testa/` A/B re-run** (Phase 1 D-21 "resumes at Phase 2") — deferred to Phase 3 (D-35). Same reason: text + colormap dependency.
- **Honoring `ncf` / `nca` color indices in `xvfacetraits_`** — D-12 reads and ignores them. Phase 3 delivers the palette and wires them to real colors.
- **Honoring the border-color parameter in `xvfrectangle_` / `xvfbordrectangle_`** — D-13 uses `foreground_` for both fill and outline. Phase 3 reintroduces the distinction once the palette is real.
- **Color-coded output generally** — `foreground_` is locked to `Qt::white` for all of Phase 2. Phase 3's `xvcouleurs_` / `xvCouleursImposees_` unlock it.
- **PostScript export (`xvpostscript_`)** — deferred to Phase 7. D-26/D-27/D-28 strip every legacy hook; Phase 7 starts fresh.
- **Pixmap save/restore slots (`sauvefenetre_`, `restaurefenetre_`, `fenetremempx_`, `mempxfenetre_`)** — deferred to Phase 4. Phase 2's `effacer_` intentionally does *not* do the legacy `fenetremempx_` copy; the single-backing model makes the copy unnecessary anyway.
- **Mouse / keyboard / event delivery** — deferred to Phase 5. Phase 2's `processEvents(ExcludeUserInputEvents)` pump explicitly excludes them.
- **`Qt::WindingFill` alternative for polygon entry points** — D-11 locks `Qt::OddEvenFill`. Phase 3's A/B gate is where a WindingFill switch would be evaluated if a case surfaces.
- **Alternate angle-convention handling in `xvarcellipse_`** — D-14 assumes the Fortran ABI passes whatever the legacy `xvuelc.c` already converted. If the planner finds the legacy body converts degrees→64ths at the entry boundary, Phase 2 mirrors that conversion at the Qt boundary (1/16th). If the Fortran caller passes raw degrees, the Qt body multiplies by 16. The planner's angle-conversion diff is the authoritative reference point.
- **Per-primitive clipping region** (`xvsetclip_` / equivalent) — not in DRAW-01..09 scope. If a future phase introduces one, it goes through `painter_->setClipRect` / `setClipPath`.
- **Multi-pixmap compositing** (overlays, underlays) — deferred to Phase 4+ if needed; Phase 2 has exactly one backing.
- **C++ unit-test framework introduction** — if the planner chooses the C++ harness option for the Phase 2 sanity driver (per Claude's Discretion), this introduces a framework dependency (Catch2, GoogleTest, or Qt Test). A Fortran extension of `xvtest0.f` avoids the question entirely and is the lowest-friction path.
- **`xvue/qt/README_RESIZE.md` vs inline documentation** — D-08 mandates the documentation exists; the location is at the planner's discretion.

### Reviewed Todos (not folded)
None — no pending todos matched Phase 2 scope at init time.

</deferred>

---

*Phase: 02-drawing-primitives-backing-pixmap*
*Context gathered: 2026-04-11*
