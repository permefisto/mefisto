# Requirements: xvue-qt

**Defined:** 2026-04-10
**Core Value:** Every MEFISTO workflow that works today through X11 — interactive meshing, elasticity/fluid/thermal/nonlinear solving, visualization, picking, plot export — keeps working through the new Qt 6 interface, with the Fortran solver code underneath completely unchanged.

## v1 Requirements

Requirements for the initial Qt 6 backend release. Each maps to exactly one roadmap phase.

### Build — CMake skeleton and ABI stubs

- [ ] **BUILD-01**: Developer can build `libxvueqt.a` from `xvue/CMakeLists.txt` on Debian trixie using `qt6-base-dev` / `qt6-base-dev-tools` from apt, without vendoring Qt or pulling non-apt packages.
- [ ] **BUILD-02**: `libxvueqt.a` is built with `CMAKE_POSITION_INDEPENDENT_CODE ON` so it links into PIE Fortran executables without relocation errors, and built **without** `-fopenmp` so it does not collide with OpenMP runtimes on `_OMP` executables.
- [ ] **BUILD-03**: `xvue/CMakeLists.txt` uses `CMAKE_AUTOMOC ON` set before `find_package(Qt6 ...)`, declares all `QObject` subclasses in headers (never in `.cpp` files), and regenerates `moc_*.cpp` files cleanly on rebuild.
- [ ] **BUILD-04**: `xvue/xvue_qt_api.h` declares every one of the ~60 `extern "C"` entry points currently in `xvue/xvuelc.c` with byte-identical names (trailing underscore preserved via `#define proc(x) x##_` or equivalent) and byte-identical pointer-based signatures (every scalar is `int*`/`float*`, every string is `char* + int*` length pair).
- [ ] **BUILD-05**: `xvue/xvue_qt_api.cpp` implements every declared entry point as a no-op stub that returns successfully, so the bridge is link-complete before any graphics logic lands.
- [ ] **BUILD-06**: A new `bin/cbl_tout_qt` variant of `bin/cbl_tout` builds all MEFISTO executables against `libxvueqt.a` (via `pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport` + `-lstdc++`), cleans `xvue/build/` and `pp/` before building, and produces working `pp/pp*_qt` executables that link cleanly.
- [ ] **BUILD-07**: The existing `bin/cbl_tout` + `xvue/xvuelc.c` + `libX11` build path continues to work unchanged, producing the same `pp/pp*` executables as before this project started.
- [ ] **BUILD-08**: A build-time sanity check (`nm libxvueqt.a | grep '_$'`) verifies every Fortran-facing symbol has the trailing underscore and fails the build if any is missing.
- [ ] **BUILD-09**: The Y-axis convention used by `xvuelc.c` is audited, documented in `xvue/README_COORDS.md` (Y-up vs Y-down, origin location, whether inversion happens in C or in Fortran), and the Qt bridge follows the same convention.
- [ ] **BUILD-10**: Five canonical `testa/` cases (one per module: mesher, elas, flui, ther, nlse) are selected and recorded as the validation baseline for every subsequent phase, with a `.planning/validation/BASELINE.md` listing them.

### Shell — QApplication and window lifecycle

- [x] **SHELL-01**: `xvinitgraphique_` opens a `QMainWindow` whose central widget is an `XvueCanvas` `QWidget`, creating the `QApplication` on first call via `std::call_once` with static fabricated `argc`/`argv`.
- [x] **SHELL-02**: `xvfermer_` closes the window but does **not** destroy the `QApplication`; the `QApplication` is torn down only at process exit via an `atexit`-registered handler. Reopening the window (second call to `xvinitgraphique_`) must not crash or log a "QApplication: there can only be one" assertion.
- [x] **SHELL-03**: `QApplication::exec()` is never called anywhere in `xvue/`, enforced by a pre-commit `grep` rule or CMake target that fails if the string appears in `xvue/*.cpp`.
- [x] **SHELL-04**: `xvpxecran_` and `xvmmecran_` return screen dimensions from `QScreen` in **logical** pixels (device-independent); the convention is documented and Fortran callers rely on it unchanged.
- [x] **SHELL-05**: `xvfond_` sets the background color of the `XvueCanvas` without corrupting the backing pixmap.
- [x] **SHELL-06**: The window correctly renders at `devicePixelRatioF > 1.0` on a HiDPI display (4K monitor or `QT_SCALE_FACTOR=2`) without doubling or halving the mesh size relative to the X11 backend.
- [x] **SHELL-07**: Every `extern "C"` entry point contains a debug-build assertion `Q_ASSERT(QThread::currentThread() == qApp->thread())` to catch accidental graphics calls from OpenMP worker threads.

### Draw — drawing primitives and backing pixmap

- [ ] **DRAW-01**: `XvueState` owns one long-lived `QPainter*` bound to a persistent off-screen `QPixmap` backing store; `QWidget::paintEvent` does nothing except `drawPixmap(0, 0, backing)`.
- [ ] **DRAW-02**: `xvtrait_`, `xvftrait_`, `xvtraits_` draw straight lines / polylines into the backing pixmap via `QPainter::drawLine` / `drawPolyline`, matching the X11 backend visually on `prpr/xvtest1.f`–`xvtest4.f`.
- [ ] **DRAW-03**: `xvface_`, `xvfacetraits_` draw filled polygons via `QPainter::drawPolygon` using a local shim `struct MefistoPoint { short x; short y; }` whose layout matches Xlib's `XPoint` byte-for-byte, so the Fortran wrappers in `xvue/*.f` remain untouched.
- [ ] **DRAW-04**: `xvrectangle_`, `xvbordrectangle_`, `xvfrectangle_`, `xvfbordrectangle_` draw filled and outlined rectangles matching the X11 backend.
- [ ] **DRAW-05**: `xvarcellipse_`, `xvbordarcellipse_` draw ellipse arcs matching the X11 backend.
- [ ] **DRAW-06**: `xvtypetrait_`, `xvepaisseur_` set pen style and width, reflected on subsequent draw calls.
- [ ] **DRAW-07**: `effacer_`, `xvvoir_`, `xvpxfenetre_` clear the canvas, flush pending draws to the screen, and query pixel dimensions respectively.
- [ ] **DRAW-08**: All drawing primitives use `QPainter::Antialiasing` rendering hint enabled by default (a free visual improvement over the X11 backend).
- [ ] **DRAW-09**: Resizing the canvas reallocates the backing pixmap and preserves the previous content via a sub-blit (optional preserve, documented convention).

### Text — text, fonts, and colormap

- [ ] **TEXT-01**: `xvchargefonte_` loads a bundled fixed Qt font (reproducible across installs) and sets it as the current font for subsequent text draws, returning its pixel width and height.
- [ ] **TEXT-02**: `xvnbpixeltexte_` returns the text extent of a given string using `QFontMetrics::horizontalAdvance` and `QFontMetrics::height`, and the returned values drive Fortran label positioning without clipping or overlap on `testa/nafems_le1` and `testa/pan2d`.
- [ ] **TEXT-03**: `xvtexte_` and `xvftexte_` render text at a given pixel position into the backing pixmap using `QPainter::drawText`.
- [ ] **TEXT-04**: `xvcouleur_` sets the current pen and brush color by index into an internal palette `std::array<QColor, MAX_PALETTE>` that mirrors the X11 indexed colormap semantics.
- [ ] **TEXT-05**: `xvCouleursImposees_`, `xvStockeRGBtoColormap_`, `xvColormapToRGB_`, `xvrecuprgbdec_`, `xvactivervb_` populate, query, and activate the internal color palette with RGB values matching the X11 backend on a 24-bit display within 1 bit per channel.
- [ ] **TEXT-06**: Scientific colormaps (stress, temperature, velocity) are frozen against system dark-mode: `QPalette` changes affect chrome only, never the backing-pixmap colors.

### Pixmap — double buffering

- [ ] **PIXMAP-01**: `fenetremempx_` and `mempxfenetre_` copy between the canvas backing pixmap and an off-screen `XvuePixmapStack` slot using `QPixmap::drawPixmap`.
- [ ] **PIXMAP-02**: `sauvefenetre_` and `restaurefenetre_` save and restore the full canvas content to/from a named off-screen slot.
- [ ] **PIXMAP-03**: `sauvemempx_`, `restauremempx_`, `effacemempx_` operate on secondary off-screen slots for intermediate buffering.
- [ ] **PIXMAP-04**: The mesher's interactive rubber-band-drag behavior works without flicker using the pixmap save/restore pattern, validated on `testa/cavity2d`.

### Event — event bridge and blocking reads

- [ ] **EVENT-01**: `XvueEventBridge` is a `QObject` installed as an event filter on the `XvueCanvas`, exposing a `waitForEvent()` method that runs a local `QEventLoop` until an awaited mouse/keyboard/pause event arrives.
- [ ] **EVENT-02**: `xvsouris_` blocks on `waitForEvent()` and returns the mouse button, coordinates, and event type — matching the X11 backend semantics of `XNextEvent`-based polling — without ever calling `QApplication::exec()` at top level.
- [ ] **EVENT-03**: `xvsouris2_` (the multi-item menu variant) blocks on `waitForEvent()` and returns item selections identically to the X11 backend.
- [ ] **EVENT-04**: `xvpause_` blocks until any event arrives, matching the X11 backend.
- [ ] **EVENT-05**: `deplsouris_` returns the current mouse position without blocking.
- [ ] **EVENT-06**: `initaccrochage_` initializes the snap/crosshair state on the canvas.
- [ ] **EVENT-07**: Mouse-motion events are coalesced via deferred `QTimer::singleShot(0, loop.quit)` so the nested loop quits after the last pending motion, mirroring Xlib's `XEventsQueued(QueuedAfterReading)` pattern — the mesher does not stutter on fast mouse drags in `testa/pan2d` or `testa/torus`.
- [ ] **EVENT-08**: A re-entrancy counter `XvueApp::blockingDepth()` tracks nested `waitForEvent()` calls so menu actions can refuse to open modal dialogs while an outer blocking read is active.

### UX — Level-3 Qt chrome and menu surface

- [ ] **UX-01**: The `QMainWindow` has a `QMenuBar` with module-appropriate top-level menus (File, Edit, View, Mesh, Solve, Help), a `QToolBar` sharing `QAction`s with the menu bar, keyboard shortcuts (`Ctrl+S`, `Ctrl+Q`, `Ctrl+O`, etc.), tooltips, and a `QStatusBar` showing current mouse coordinates and hover-pick feedback.
- [ ] **UX-02**: `XvueMenuBridge` exposes a `pendingCommands_` queue; `QAction::triggered` signals push lexicon command strings into the queue without ever calling Fortran or opening modal dialogs directly.
- [ ] **UX-03**: The next `xvsouris_` / `xvsouris2_` call drains the menu command queue and returns a synthetic keyboard event (`notypeevent = 2`) carrying the queued command — Fortran dispatches through its existing text-lexicon code path, with no changes to any solver driver.
- [ ] **UX-04**: Modal dialogs (file open/save, properties) refuse to open (and show a status-bar message instead) when `XvueApp::blockingDepth() > 0`, preventing re-entrant `QDialog::exec()` inside `xvsouris` nested loops.
- [ ] **UX-05**: The per-module lexicon audit of `mail/`, `elas/`, `flui/`, `ther/`, `nlse/` enumerates every interactive command currently typed by the user and produces a `QAction` entry for each, grouped into menus and toolbars. Deliverable: `.planning/phase-6/LEXICON-AUDIT.md` per module.
- [ ] **UX-06**: File operations (`Open Project…`, `Save Project`, `Save Project As…`, `Export…`) use `QFileDialog` with a project-directory filter respecting `$MEFISTOX/<project>/`, and a recent-projects list persisted via `QSettings`.
- [ ] **UX-07**: Window geometry, dock layout, recent-projects list, and user preferences are persisted via `QSettings` across sessions.
- [ ] **UX-08**: The existing bilingual flag (`$MEFISTO/td/m/anglais`) selects the UI language; all menu labels, tooltips, dialogs, status messages, and the About dialog are available in both French and English.
- [ ] **UX-09**: An About dialog credits Alain Perronnet (LJLL/UPMC Paris) and shows the MEFISTO version, Qt version, and build date. A Help menu entry launches the PDF documentation under `doc/` via `QDesktopServices::openUrl`.
- [ ] **UX-10**: A console/log `QDockWidget` displays the Fortran solver's stdout in real time via a `QProcess` stdout pipe-reader, exposing a "keystone" surface that also feeds progress reporting and error dialogs.
- [ ] **UX-11**: Error lines matching the existing Fortran diagnostic patterns (e.g. lines beginning `*** ERREUR`) are parsed from the stdout stream and surfaced as `QMessageBox` alerts.
- [ ] **UX-12**: Mouse wheel zooms, middle-drag pans, and right-click opens a context menu on the canvas. A coordinate readout in the status bar updates live on mouse move.
- [ ] **UX-13**: Chrome follows the system dark-mode palette via `QPalette`; scientific colormaps (REQ TEXT-06) are unaffected.

### Export — image, GIF, and PostScript

- [ ] **EXPORT-01**: At the start of Phase 7 work, a 10-line probe prints `QImageWriter::supportedImageFormats()`; the output (GIF writer present or absent) is recorded in `.planning/phase-7/PROBE.md` and determines whether animated GIF uses `QImageWriter` or a per-frame-PNG + `ffmpeg` fallback (never ImageMagick).
- [ ] **EXPORT-02**: `XvueExport` writes PNG and JPEG single-frame snapshots of the current canvas via `QImageWriter` / `QImage::save`.
- [ ] **EXPORT-03**: `XvueExport` writes animated GIFs from a sequence of frames via the strategy chosen by EXPORT-01, matching the visual output of the legacy `bin/convertepsgif` pipeline on `testa/wave` and `testa/cavity2d` post-processing.
- [ ] **EXPORT-04**: `xvpostscript_` is implemented in `xvue/xvue_qt_postscript.cpp` by moving the existing hand-rolled PostScript emitter (~120 lines of `fprintf`) **verbatim** from `xvuelc.c`; no switch to `QPrinter` PostScript output.
- [ ] **EXPORT-05**: PDF export is added as a bonus via a new entry point using `QPrinter::PdfFormat`, not by modifying `xvpostscript_`.
- [ ] **EXPORT-06**: After EXPORT-01..05 ship, the ImageMagick runtime dependency is no longer invoked by any Qt-backend code path.

### Validation — A/B against testa subset

- [ ] **VALID-01**: All 5 canonical `testa/` cases (BUILD-10 baseline) pass visually (side-by-side comparison of X11 and Qt windows) at the end of every phase from 0 through 7; results are logged to `.planning/phase-N/VALIDATION.md`.
- [ ] **VALID-02**: The X11 backend continues to pass the same 5 cases at the end of every phase, proving no accidental damage to the legacy path during the A/B window.
- [ ] **VALID-03**: At the end of Phase 7, all 5 cases are run through `ELASTICER_OMP` (the OpenMP executable) to verify the Qt layer's main-thread-only invariant holds under OpenMP solver parallelism.
- [ ] **VALID-04**: At the end of Phase 7, all 5 cases are run on a HiDPI display (4K or `QT_SCALE_FACTOR=2`) to verify no size / position drift vs the X11 backend.
- [ ] **VALID-05**: Color-accuracy spot checks (stress/temperature/velocity color bars) on `testa/nafems_le1`, `testa/heat1d`, and `testa/cavity2d` show no user-visible drift between backends.
- [ ] **VALID-06**: Font-metric spot checks (node-number labels) on `testa/pan2d` and `testa/hexahedron` show no clipping or overlap in the Qt backend.
- [ ] **VALID-07**: A validation checklist (`.planning/phase-8/CHECKLIST.md`) records pass/fail per case per backend and is the gate for declaring v1 shippable.

### Retire — X11 backend removal

- [ ] **RETIRE-01**: After the one-release-cycle A/B window closes, `xvue/xvuelc.c` and `bin/ccxvue` are deleted.
- [ ] **RETIRE-02**: `libX11` linker lines and the hardcoded `/usr/X11R6/lib64` path are removed from all `bin/cb*` scripts.
- [ ] **RETIRE-03**: `bin/convertepsgif` and all `convert` shell-outs are audited across `bin/`, `td/`, `testa/`, `testf/` (`grep -rn 'convert' bin/ td/ testa/ testf/`) and replaced or removed. ImageMagick is dropped from the install documentation.
- [ ] **RETIRE-04**: `README`, `LISEZMOI`, and install scripts are updated to reference only Qt 6 runtime dependencies (`qt6-base`, `libqt6imageformats6-plugins`); `libX11-dev` and ImageMagick are removed from the dependency list.

## v2 Requirements

Deferred differentiators worth revisiting post-v1 after the A/B window closes. Tracked but not in the current roadmap.

### v1.x Usability

- **V1X-01**: Crash-recovery lockfile and "Recover last project?" prompt on startup (mitigates the legacy Fortran `STOP`-loses-state fragility from CONCERNS.md).
- **V1X-02**: Command log dock that records the lexicon-equivalent of every GUI action, with save-and-replay to produce scripted sessions — the "scripting" answer without introducing an embedded scripting language.
- **V1X-03**: Scientific colormap presets (viridis, plasma, grayscale) selectable per-plot; the legacy rainbow remains the default so A/B comparison stays valid.
- **V1X-04**: Color bar widget with log-scale toggle for solver result plots.
- **V1X-05**: Mesh statistics dock (element count, quality histogram, edge lengths) fed from the existing `mail/` internal state.
- **V1X-06**: Per-module saved workspace layouts via `QSettings`.
- **V1X-07**: Progress reporting bar fed by parsing iteration/residual lines from the solver stdout pipe-reader (UX-10).
- **V1X-08**: Animation scrubber for time-stepping solvers (`flui`, `ther`, `nlse`) using the existing per-timestep frame output.
- **V1X-09**: Single-level snapshot-based undo for mesh edits using scratch-file snapshots (zero Fortran-side change).

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| Multiple viewports / split view | Depends on `xvuelc.c` state being fully instance-scoped; defer until post-v1 when the backend is proven stable. |
| Mesh-quality color overlay | Requires per-element quality metrics to be exposed from `mail/` via a new API — violates "no Fortran solver changes." |
| 3D clip planes / section cuts | Requires a `QOpenGLWidget` rewrite of the display path; out of scope for v1 which stays on `QPainter` / `QPixmap`. |
| Live-switch bilingual UI via `QTranslator` | The existing `$MEFISTO/td/m/anglais` flag + process restart is sufficient; live-switch is a nice-to-have not worth the complexity. |
| VTK / CSV export for ParaView | Separate initiative; not required by Core Value. |
| Cooperative solver cancellation (Stop button) | Requires cooperative Fortran-side changes; violates the Core Value invariant. |
| Windows / macOS ports | Out of scope per PROJECT.md — Linux x86_64 only. |
| Changes to numerical solver code | Core Value invariant: Fortran side must not notice the migration. |
| Full CMake migration of the non-`xvue/` build | Out of scope per PROJECT.md — shell scripts continue to drive the Fortran build. |
| Automated CI (GitHub Actions / GitLab CI / Jenkins) | Out of scope per PROJECT.md — quality is maintained by running `testa/` A/B manually. |
| Fortran modernization (`IMPLICIT NONE`, F90 modules) | Orthogonal to the graphics migration; separate initiative. |
| New physics features or new solvers | This is a UI migration, not a feature release. |
| Custom Qt theming engine | Unbounded scope; system dark/light mode + `QPalette` is the table-stakes answer. |
| Plugin system / embedded scripting language (Python/Lua/QtScript) | Massive complexity; the command log dock (V1X-02) delivers the "scripted workflow" use case without the runtime. |
| Cloud sync / real-time collaboration | Not aligned with a single-user desktop scientific tool; out of scope forever. |
| Auto-updater / crash-reporter upload / telemetry | Out of scope; scientific users expect no phone-home. |
| Integrated code editor | Wrong tool for the job; users have their own editor. |
| OpenCASCADE / built-in CAD kernel | Separate project; not in scope. |
| MDI (multi-document interface) | Conflicts with the single `QApplication` singleton discipline; one project per process is the model. |
| Qt re-skinning of the bash launcher scripts (`bin/INITIER`, etc.) | Launchers stay as bash — only the graphics windows they spawn become Qt. |
| Rewriting the text-lexicon parser | The menu-queue pattern dispatches through the existing lexicon parser unchanged. |
| Integrated test runner for `testa/` | Tests are run manually; automating them is a separate initiative (see Out of Scope: CI). |

## Traceability

Populated during roadmap creation — all v1 requirements mapped to exactly one phase.

| Requirement | Phase | Status |
|-------------|-------|--------|
| BUILD-01 | Phase 0 | Pending |
| BUILD-02 | Phase 0 | Pending |
| BUILD-03 | Phase 0 | Pending |
| BUILD-04 | Phase 0 | Pending |
| BUILD-05 | Phase 0 | Pending |
| BUILD-06 | Phase 0 | Pending |
| BUILD-07 | Phase 0 | Pending |
| BUILD-08 | Phase 0 | Pending |
| BUILD-09 | Phase 0 | Pending |
| BUILD-10 | Phase 0 | Pending |
| SHELL-01 | Phase 1 | Complete — Validated in Phase 01-03 (human visual: two 800x600 MEFISTO windows, exit 0) |
| SHELL-02 | Phase 1 | Complete — Validated in Phase 01-03 (human visual: second window, no singleton assertion) |
| SHELL-03 | Phase 1 | Complete — Validated in Phase 01-01 (verify_no_exec injection test) |
| SHELL-04 | Phase 1 | Complete — Validated in Phase 01-02 (xvpxecran_/xvmmecran_ compiled and linked) |
| SHELL-05 | Phase 1 | Complete — Validated in Phase 01-02 (xvfond_ compiled and linked) |
| SHELL-06 | Phase 1 | Complete — Validated in Phase 01-03 (human visual: QT_SCALE_FACTOR=2 visibly larger windows, exit 0) |
| SHELL-07 | Phase 1 | Complete — Validated in Phase 01-02 (57/57 stubs carry XVUE_QT_ASSERT_MAIN_THREAD) |
| DRAW-01 | Phase 2 | Pending |
| DRAW-02 | Phase 2 | Pending |
| DRAW-03 | Phase 2 | Pending |
| DRAW-04 | Phase 2 | Pending |
| DRAW-05 | Phase 2 | Pending |
| DRAW-06 | Phase 2 | Pending |
| DRAW-07 | Phase 2 | Pending |
| DRAW-08 | Phase 2 | Pending |
| DRAW-09 | Phase 2 | Pending |
| TEXT-01 | Phase 3 | Pending |
| TEXT-02 | Phase 3 | Pending |
| TEXT-03 | Phase 3 | Pending |
| TEXT-04 | Phase 3 | Pending |
| TEXT-05 | Phase 3 | Pending |
| TEXT-06 | Phase 3 | Pending |
| PIXMAP-01 | Phase 4 | Pending |
| PIXMAP-02 | Phase 4 | Pending |
| PIXMAP-03 | Phase 4 | Pending |
| PIXMAP-04 | Phase 4 | Pending |
| EVENT-01 | Phase 5 | Pending |
| EVENT-02 | Phase 5 | Pending |
| EVENT-03 | Phase 5 | Pending |
| EVENT-04 | Phase 5 | Pending |
| EVENT-05 | Phase 5 | Pending |
| EVENT-06 | Phase 5 | Pending |
| EVENT-07 | Phase 5 | Pending |
| EVENT-08 | Phase 5 | Pending |
| UX-01 | Phase 6.0 | Pending |
| UX-02 | Phase 6.0 | Pending |
| UX-03 | Phase 6.0 | Pending |
| UX-04 | Phase 6.0 | Pending |
| UX-05 | Phase 6.1-6.5 | Pending |
| UX-06 | Phase 6.0 | Pending |
| UX-07 | Phase 6.0 | Pending |
| UX-08 | Phase 6.0 | Pending |
| UX-09 | Phase 6.0 | Pending |
| UX-10 | Phase 6.0 | Pending |
| UX-11 | Phase 6.0 | Pending |
| UX-12 | Phase 6.0 | Pending |
| UX-13 | Phase 6.0 | Pending |
| EXPORT-01 | Phase 7 | Pending |
| EXPORT-02 | Phase 7 | Pending |
| EXPORT-03 | Phase 7 | Pending |
| EXPORT-04 | Phase 7 | Pending |
| EXPORT-05 | Phase 7 | Pending |
| EXPORT-06 | Phase 7 | Pending |
| VALID-01 | Phase 8 | Pending |
| VALID-02 | Phase 8 | Pending |
| VALID-03 | Phase 8 | Pending |
| VALID-04 | Phase 8 | Pending |
| VALID-05 | Phase 8 | Pending |
| VALID-06 | Phase 8 | Pending |
| VALID-07 | Phase 8 | Pending |
| RETIRE-01 | Phase 9 | Pending |
| RETIRE-02 | Phase 9 | Pending |
| RETIRE-03 | Phase 9 | Pending |
| RETIRE-04 | Phase 9 | Pending |

**Coverage:**
- v1 requirements: 72 total
- Mapped to phases: 72
- Unmapped: 0 ✓

---
*Requirements defined: 2026-04-10*
*Last updated: 2026-04-11 — SHELL-01/02/06 validated in Phase 01-03; SHELL-03/04/05/07 validated in Phase 01-01/01-02*
