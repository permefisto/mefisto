# Roadmap: xvue-qt

## Overview

xvue-qt replaces MEFISTO's ~3619-line `xvue/xvuelc.c` Xlib backend with a Qt 6 C++ implementation behind a byte-identical `extern "C"` ABI, then layers a Level-3 UX modernization (`QMainWindow`, `QMenuBar`, `QToolBar`, `QDialog`, dock widgets, HiDPI/dark-mode chrome, Qt-native image/GIF/PDF export) on top. The Fortran solver tree (`mail/`, `elas/`, `flui/`, `ther/`, `nlse/`, `reso/`, `util/`) stays untouched — the migration's invariant is "Fortran must not notice the change." Phases follow a near-linear dependency chain: build plumbing (0) gates everything; a window (1) must exist before drawing (2, 3, 4); drawing must work before blocking event reads (5); menus (6) build on the same event plumbing as `xvsouris`; export (7) is the one independent slot; A/B validation (8) is the ship gate; X11 retirement (9) is trailing cleanup after the one-release-cycle A/B window closes.

## Phases

**Phase Numbering:**
- Integer phases (0, 1, 2, ...): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [ ] **Phase 0: Build skeleton & ABI stubs** - CMake + `libxvueqt.a` + no-op `extern "C"` stubs + `bin/cbl_tout_qt` link-through
- [x] **Phase 1: Window shell** - `XvueApp`/`XvueWindow`/`XvueCanvas` — an empty Qt window opens via `xvinitgraphique_`
- [ ] **Phase 2: Drawing primitives & backing pixmap** - `XvueState` + persistent `QPainter` on off-screen `QPixmap`; lines, polygons, rectangles, ellipses
- [ ] **Phase 3: Text, fonts, colormap** - Font loading, text metrics, indexed palette mirroring X11 colormap semantics
- [ ] **Phase 4: Pixmap save/restore** - Double-buffering slots for `fenetremempx`/`sauvefenetre` rubber-band workflow
- [ ] **Phase 5: Event bridge & blocking reads** - Nested `QEventLoop` pattern for `xvsouris`/`xvpause`; the architectural pivot
- [ ] **Phase 6: Level-3 UX chrome** - `QMenuBar`/`QToolBar`/`QAction` menu surface + per-module lexicon audit + dock/console/dialogs
- [ ] **Phase 7: Image, GIF, PostScript export** - Qt-native PNG/JPEG/PDF, preserved verbatim PostScript emitter, ImageMagick drop
- [ ] **Phase 8: A/B validation on testa subset** - Side-by-side X11 vs Qt on 5 canonical cases; HiDPI and `_OMP` sweeps
- [ ] **Phase 9: Retire X11 backend** - Delete `xvuelc.c`, `bin/ccxvue`, `libX11` linker lines, ImageMagick deps (gated on A/B window close)

## Phase Details

### Phase 0: Build skeleton & ABI stubs
**Goal**: Prove the Qt/CMake build plumbing integrates cleanly with the legacy shell linker and produces link-complete `pp/pp*_qt` executables, before any graphics logic exists.
**Depends on**: Nothing (first phase)
**Requirements**: BUILD-01, BUILD-02, BUILD-03, BUILD-04, BUILD-05, BUILD-06, BUILD-07, BUILD-08, BUILD-09, BUILD-10
**Success Criteria** (what must be TRUE):
  1. Developer can run `bin/cbl_tout_qt` on Debian trixie and produce `pp/ppmail_qt`, `pp/ppelas_qt`, `pp/ppflui_qt`, `pp/ppther_qt`, `pp/ppnlse_qt` that link cleanly against `libxvueqt.a` using only apt-provided Qt 6 packages.
  2. `nm libxvueqt.a | grep '_$'` shows every one of the ~60 entry points from `xvuelc.c` present with trailing underscore; a CMake target fails the build if any is missing.
  3. Running any `pp/pp*_qt` executable proceeds past the link stage and invokes the no-op ABI stubs without crashing (no graphics output yet — the bridge is link-complete only).
  4. `bin/cbl_tout` (the legacy X11 path) still builds the original `pp/pp*` executables unchanged, and all 5 canonical `testa/` cases selected as baseline still run on X11.
  5. `xvue/README_COORDS.md` documents the audited Y-axis convention (origin, up/down, where inversion happens) and `.planning/validation/BASELINE.md` lists the 5 canonical `testa/` cases (one per module: mesher, elas, flui, ther, nlse).
**Plans**: 4 plans
  - [x] 00-01-PLAN.md — Prerequisites: apt install qt6-base-dev + verify legacy X11 baseline (BUILD-07 anchor)
  - [x] 00-02-PLAN.md — CMake scaffold + xvue_qt_api.h (57 entries) + warn-once stubs → libxvueqt.a (BUILD-01..05, BUILD-08)
  - [x] 00-03-PLAN.md — bin/cbl_tout_qt + 5 per-module cb*_qt clones → pp/pp*_qt smoke-tested (BUILD-06)
  - [x] 00-04-PLAN.md — xvue/README_COORDS.md + .planning/validation/BASELINE.md + clean-tree rebuild + 5-case legacy run (BUILD-09, BUILD-10)

### Phase 1: Window shell (`XvueApp` + `XvueWindow` + `XvueCanvas`)
**Goal**: A blank Qt window opens through `xvinitgraphique_` and closes cleanly through `xvfermer_`, proving the `QApplication` singleton discipline and HiDPI convention in isolation before any drawing logic.
**Depends on**: Phase 0
**Requirements**: SHELL-01, SHELL-02, SHELL-03, SHELL-04, SHELL-05, SHELL-06, SHELL-07
**Success Criteria** (what must be TRUE):
  1. Developer can run `pp/ppmail_qt` on a test driver that calls `xvinitgraphique_` and observe a blank `QMainWindow` with an `XvueCanvas` central widget appear on screen.
  2. Calling `xvfermer_` then `xvinitgraphique_` a second time reopens the window with no "QApplication: there can only be one" assertion and no crash on process exit.
  3. `grep -rn 'QApplication::exec' xvue/*.cpp` returns zero matches, enforced by a pre-commit/CMake check that fails the build otherwise.
  4. `xvpxecran_`/`xvmmecran_` return screen dimensions in logical pixels from `QScreen`, and the window renders identically-sized to the X11 backend on a `QT_SCALE_FACTOR=2` or 4K HiDPI display.
  5. Every `extern "C"` entry point contains a debug-build `Q_ASSERT(QThread::currentThread() == qApp->thread())` to catch off-main-thread graphics calls.
**Plans**: 3 plans
  - [x] 01-01-PLAN.md — Scaffold XvueApp/Window/Canvas/State + CMake sources + verify_no_exec guard + SHELL-07 macro body
  - [x] 01-02-PLAN.md — Rewrite 7 SHELL entry bodies in xvue_qt_api.cpp + bulk macro retrofit into all stubs
  - [x] 01-03-PLAN.md — prpr/xvtest0.f + bin/cbxvtest0_qt + full-suite validation wave (visual gate)
**UI hint**: yes

### Phase 2: Drawing primitives & backing pixmap
**Goal**: All pure synchronous drawing primitives render into a persistent off-screen `QPixmap` via one long-lived `QPainter`, matching the X11 backend visually on `prpr/xvtest1.f`–`xvtest4.f`.
**Depends on**: Phase 1
**Requirements**: DRAW-01, DRAW-02, DRAW-03, DRAW-04, DRAW-05, DRAW-06, DRAW-07, DRAW-08, DRAW-09
**Success Criteria** (what must be TRUE):
  1. Running `prpr/xvtest1.f` through `xvtest4.f` against the Qt backend produces lines, polylines, filled polygons, rectangles, and ellipse arcs visually indistinguishable from the X11 backend on the same inputs.
  2. `XvueState` holds exactly one `QPainter*` bound to the backing `QPixmap` for the widget's lifetime; `QWidget::paintEvent` body is a single `drawPixmap(0, 0, backing)` call.
  3. `xvface_` and `xvfacetraits_` work unchanged from their Fortran callers because the C bridge uses a `struct MefistoPoint { short x; short y; }` byte-identical to `XPoint`.
  4. Pen style (`xvtypetrait_`), pen width (`xvepaisseur_`), and `effacer_`/`xvvoir_`/`xvpxfenetre_` behave identically to X11; `QPainter::Antialiasing` is enabled by default.
  5. Resizing the canvas window reallocates the backing pixmap preserving prior content per a documented convention.
**Plans**: 4 plans
  - [x] 02-01-PLAN.md — Wave 0: MEFISTO_POINT_AUDIT + README_RESIZE + ABI ellipse signature fix (int*→float*) + extended xvtest0.f draw-coverage driver
  - [x] 02-02-PLAN.md — Wave 1: XvueState growth + XvueCanvas paintEvent/resizeEvent backing lifecycle + effacer_/xvvoir_/xvpxfenetre_/xvfond_ real bodies (DRAW-01, DRAW-07, DRAW-08, DRAW-09)
  - [x] 02-03-PLAN.md — Wave 2: line/polyline/polygon/rectangle/arc primitives with drawPie/drawArc x16 angle correction + pen state (DRAW-02..06)
  - [x] 02-04-PLAN.md — Wave 3: clean rebuild, human visual checkpoint, 02-VALIDATION.md closure
**UI hint**: yes

### Phase 3: Text, fonts, colormap
**Goal**: Text rendering and the indexed-palette colormap faithfully mirror the X11 backend, with scientific colormaps frozen against system dark-mode so color-encoded physical data stays accurate.
**Depends on**: Phase 2
**Requirements**: TEXT-01, TEXT-02, TEXT-03, TEXT-04, TEXT-05, TEXT-06
**Success Criteria** (what must be TRUE):
  1. `xvchargefonte_` loads a bundled fixed Qt font reproducibly across installs and `xvnbpixeltexte_` returns extents consistent with `QFontMetrics::horizontalAdvance`/`height`.
  2. Node-number labels rendered via `xvtexte_` on `testa/nafems_le1` and `testa/pan2d` show no clipping or overlap, matching the X11 layout byte-for-byte where possible.
  3. `xvcouleur_` plus the `xvCouleursImposees_`/`xvStockeRGBtoColormap_`/`xvColormapToRGB_`/`xvrecuprgbdec_`/`xvactivervb_` family populate and query a `std::array<QColor, MAX_PALETTE>` whose RGB values match X11 within 1 bit per channel on a 24-bit display.
  4. Enabling system dark-mode changes only the window chrome via `QPalette`; the backing-pixmap stress/temperature/velocity colormaps are unchanged.
**Plans**: 4 plans
  - [x] 03-01-PLAN.md — Wave 0: bundle DejaVuSansMono TTF + qt_add_resources + verify_no_exec palette grep + XvueState palette/font grow + XvueApp font load + XvueCanvas ctor/resize guards + xvtest0.f TEXT coverage extension
  - [x] 03-02-PLAN.md — Wave 1: 7 real bodies (xvchargefonte_, xvnbpixeltexte_, xvtexte_/xvftexte_ collapsed, xvcouleur_, xvactivervb_ bulk-load, xvrecuprgbdec_) + xvinfo_ palette/font fill + retire xvfond_ hack
  - [x] 03-03-PLAN.md — Wave 2: clean rebuild + xvtest0 visual checkpoint (fonts + 8 colored lines + xvactivervb + xvnbpixeltexte box + resize preserve)
  - [ ] 03-04-PLAN.md — Wave 3: A/B catch-up gate vs prpr/xvtest1..4 + 5 canonical testa cases (D-26 hard gate) + 03-VALIDATION.md closure
**UI hint**: yes

### Phase 03.1: Build xvtest1..4 driver infrastructure on both backends (Qt + legacy X11) to unblock Phase 3 A/B catch-up gate (INSERTED)

**Goal:** Produce pp/ppxvtest{0..4}(_qt) on both Qt and legacy X11 backends (9 new cb* scripts + 2 top-level edits), fix the Qt-side xtinit_ warn-once stub so XVOUVRIR opens a real window, and harden plan 03-04 Task 1 smoke (timeout 5s, enriched grep including xtinit_, xvtest1 STOP tolerance) so the Phase 3 A/B visual gate can execute.
**Requirements**: BUILD-07, VALID-02, TEXT-01..06 (inherited — no new IDs)
**Depends on:** Phase 3
**Plans:** 3 plans

Plans:
- [ ] 03.1-01-PLAN.md — Wave 0: A6 grep verification + Qt-side xtinit_ promotion (xvue/qt/src/xvue_qt_api.cpp), X11 untouched
- [ ] 03.1-02-PLAN.md — Wave 1: 9 new cb* scripts (5 legacy + 4 Qt) + top-level wiring into cbl_tout{,_qt}
- [ ] 03.1-03-PLAN.md — Wave 2: hardened headless smoke + ABI recheck + 03-04 Task 1 patch + NOTES.md + SUMMARY.md

### Phase 4: Pixmap save/restore (double-buffering)
**Goal**: Named off-screen pixmap slots support the mesher's rubber-band drag workflow without flicker, reproducing the X11 double-buffering semantics.
**Depends on**: Phase 2
**Requirements**: PIXMAP-01, PIXMAP-02, PIXMAP-03, PIXMAP-04
**Success Criteria** (what must be TRUE):
  1. `fenetremempx_`/`mempxfenetre_` copy bytes between the canvas backing pixmap and an `XvuePixmapStack` slot via `QPixmap::drawPixmap` and survive a round trip unchanged.
  2. `sauvefenetre_`/`restaurefenetre_` save/restore the full canvas to a named slot without color or geometry drift.
  3. `sauvemempx_`/`restauremempx_`/`effacemempx_` handle secondary intermediate slots, allowing nested save/restore patterns.
  4. Running `testa/cavity2d` through the mesher exhibits flicker-free rubber-band-drag behavior, matching the X11 backend.
**Plans**: TBD

### Phase 5: Event bridge & blocking reads
**Goal**: Blocking calls (`xvsouris`, `xvpause`, `deplsouris`) run on a nested local `QEventLoop` with proper re-entrancy and motion coalescing — the architectural pivot that makes the mesher interactive end-to-end.
**Depends on**: Phase 4
**Requirements**: EVENT-01, EVENT-02, EVENT-03, EVENT-04, EVENT-05, EVENT-06, EVENT-07, EVENT-08
**Success Criteria** (what must be TRUE):
  1. Developer can run `pp/ppmail_qt testa/pan2d` and interactively create a mesh with the mouse — button presses, mouse moves, and keyboard events all dispatch back into Fortran through `xvsouris_`/`xvsouris2_` on a nested `QEventLoop` without ever calling `QApplication::exec()` at top level.
  2. Fast mouse drags on `testa/pan2d` and `testa/torus` do not stutter; mouse-motion events are coalesced via deferred `QTimer::singleShot(0, quit)` mirroring `XEventsQueued(QueuedAfterReading)` semantics.
  3. `xvpause_` blocks until any event arrives and `deplsouris_` returns current mouse position without blocking, both matching X11 backend semantics.
  4. `XvueApp::blockingDepth()` counter correctly tracks nested `waitForEvent()` calls so subsequent phases can gate modal dialogs on it.
  5. End-to-end A/B run of all 5 baseline `testa/` cases on the Qt backend produces interactive behavior indistinguishable from X11 to the developer's eye.
**Plans**: TBD (research-flag: empirical mouse-motion coalescing validation needed at plan time)
**UI hint**: yes

### Phase 6: Level-3 UX chrome & menu surface
**Goal**: A full Qt `QMainWindow` shell (menu bar, toolbar, status bar, dock widgets, dialogs, HiDPI/dark-mode) dispatches user actions into the existing Fortran text-lexicon through a `QAction`-queues-command-strings pattern, without any changes to solver drivers.
**Depends on**: Phase 5
**Requirements**: UX-01, UX-02, UX-03, UX-04, UX-05, UX-06, UX-07, UX-08, UX-09, UX-10, UX-11, UX-12, UX-13
**Success Criteria** (what must be TRUE):
  1. User can click any menu item or toolbar button in `pp/ppmail_qt` and observe the same Fortran subroutine execute that the corresponding typed lexicon command would trigger — `QAction::triggered` pushes a command string into `XvueMenuBridge::pendingCommands_` which the next `xvsouris_` drains and returns as a synthetic `notypeevent=2` keyboard event.
  2. Modal `QFileDialog`/`QDialog` entry points refuse to open with a status-bar message when `XvueApp::blockingDepth() > 0`, preventing re-entrant `QDialog::exec()` inside a nested `xvsouris_` loop.
  3. A `.planning/phase-6/LEXICON-AUDIT.md` per module (`mail/`, `elas/`, `flui/`, `ther/`, `nlse/`) enumerates every interactive lexicon command and maps it to a `QAction` with keyboard shortcut and tooltip; every enumerated action is reachable from the menu bar or toolbar.
  4. Window geometry, dock layout, recent-projects list, and preferences persist via `QSettings` across process restarts; the bilingual flag `$MEFISTO/td/m/anglais` correctly selects FR/EN labels across all menus, dialogs, and the About dialog (which credits Alain Perronnet / LJLL / UPMC Paris).
  5. A console `QDockWidget` displays Fortran solver stdout in real time via a `QProcess` pipe-reader; lines matching `*** ERREUR` surface as `QMessageBox` alerts. Mouse wheel zoom, middle-drag pan, right-click context menu, live coordinate readout, and system dark-mode chrome all work on the canvas without affecting scientific colormaps.
**Plans**: TBD (research-flag: per-module lexicon audit may split into 5 sub-phases — one per solver module — at plan time)
**UI hint**: yes

### Phase 7: Image, GIF, and PostScript export
**Goal**: Qt-native PNG/JPEG/PDF plus a preserved-verbatim PostScript emitter and a runtime-probed animated GIF path fully replace the legacy `xvuelc.c` export code and `bin/convertepsgif` ImageMagick shell-out.
**Depends on**: Phase 6
**Requirements**: EXPORT-01, EXPORT-02, EXPORT-03, EXPORT-04, EXPORT-05, EXPORT-06
**Success Criteria** (what must be TRUE):
  1. `.planning/phase-7/PROBE.md` records the output of a `QImageWriter::supportedImageFormats()` probe run at phase kickoff, and the chosen animated-GIF strategy (`QImageWriter` loop or per-frame PNG + `ffmpeg` fallback — never ImageMagick) is documented and implemented.
  2. Developer can export the current canvas to PNG and JPEG via `XvueExport` and the files open correctly in a standard image viewer; running `testa/wave` and `testa/cavity2d` post-processing produces animated GIFs visually matching the legacy `bin/convertepsgif` pipeline.
  3. `xvpostscript_` is implemented in `xvue/xvue_qt_postscript.cpp` by moving the existing ~120-line hand-rolled `fprintf` emitter **verbatim** from `xvuelc.c`; diffing an exported `.ps` against the X11 backend shows byte-for-byte (or near byte-for-byte) parity on the same scene.
  4. PDF export via a new `QPrinter::PdfFormat` entry point works as an additive bonus without modifying `xvpostscript_`.
  5. No Qt-backend code path invokes ImageMagick's `convert` — verified by `grep -rn 'convert' xvue/` returning only the PDF keyword or similar.
**Plans**: TBD (research-flag: run GIF writer probe before committing to strategy)

### Phase 8: A/B validation on testa subset
**Goal**: The 5 canonical `testa/` cases run side-by-side on both backends — including HiDPI and `_OMP` sweeps — with no user-visible drift, producing a checklist that gates declaring v1 shippable.
**Depends on**: Phase 7
**Requirements**: VALID-01, VALID-02, VALID-03, VALID-04, VALID-05, VALID-06, VALID-07
**Success Criteria** (what must be TRUE):
  1. All 5 baseline `testa/` cases pass a visual side-by-side comparison between X11 and Qt backends; `.planning/phase-8/VALIDATION.md` logs pass/fail per case per backend.
  2. All 5 cases also pass when run through `ELASTICER_OMP` (or the equivalent `_OMP` executable per module) proving the Qt layer's main-thread-only invariant holds under OpenMP solver parallelism.
  3. All 5 cases render correctly on a HiDPI display (4K monitor or `QT_SCALE_FACTOR=2`) with no size or position drift vs X11.
  4. Color-bar spot checks on `testa/nafems_le1`, `testa/heat1d`, and `testa/cavity2d` show no user-visible drift; font-metric checks on `testa/pan2d` and `testa/hexahedron` show no clipping or overlap.
  5. `.planning/phase-8/CHECKLIST.md` is signed off as the explicit gate for declaring v1 shippable and opening the one-release-cycle A/B window.
**Plans**: TBD
**UI hint**: yes

### Phase 9: Retire X11 backend
**Goal**: After the one-release-cycle A/B window has closed, delete every X11 and ImageMagick code path, linker line, script, and documentation reference, leaving Qt 6 as the single graphics backend.
**Depends on**: Phase 8 **AND** one-release-cycle A/B window closed (process gate — not a date, not started until the maintainer explicitly confirms the window is closed)
**Requirements**: RETIRE-01, RETIRE-02, RETIRE-03, RETIRE-04
**Success Criteria** (what must be TRUE):
  1. `xvue/xvuelc.c`, `bin/ccxvue`, and every `libX11` / `/usr/X11R6/lib64` linker line across `bin/cb*` scripts are deleted; a clean build from the tree succeeds with no X11 references.
  2. `grep -rn 'convert' bin/ td/ testa/ testf/` returns no ImageMagick invocations; `bin/convertepsgif` is deleted or replaced.
  3. `README`, `LISEZMOI`, and install scripts list only Qt 6 runtime dependencies (`qt6-base`, `libqt6imageformats6-plugins`); `libX11-dev` and ImageMagick are removed from the dependency list.
  4. Running `bin/cbl_tout` produces a working MEFISTO build with Qt 6 as the only graphics backend and all 5 baseline `testa/` cases still pass.
**Plans**: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 0 → 1 → 2 → 3 → 4 → 5 → 6 → 7 → 8 → 9

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 0. Build skeleton & ABI stubs | 0/TBD | Not started | - |
| 1. Window shell | 3/3 | Complete | 2026-04-11 |
| 2. Drawing primitives & backing pixmap | 4/4 | Complete | 2026-04-11 |
| 3. Text, fonts, colormap | 0/TBD | Not started | - |
| 4. Pixmap save/restore | 0/TBD | Not started | - |
| 5. Event bridge & blocking reads | 0/TBD | Not started | - |
| 6. Level-3 UX chrome | 0/TBD | Not started | - |
| 7. Image, GIF, PostScript export | 0/TBD | Not started | - |
| 8. A/B validation on testa subset | 0/TBD | Not started | - |
| 9. Retire X11 backend | 0/TBD | Not started (gated) | - |

## Coverage

- v1 requirements: 72 total
- Mapped to phases: 72
- Unmapped: 0 ✓

All 72 v1 requirements from REQUIREMENTS.md are mapped to exactly one phase. The 9-phase structure proposed in `.planning/research/SUMMARY.md` is adopted as-is — its near-linear dependency chain (0 → 1 → {2, 3, 4} → 5 → 6 → 7 → 8 → 9) is forced by the architecture (backing pixmap requires a window; blocking reads require drawing; menus share the event bridge; export is the only independent slot; retirement is gated on A/B) and is reinforced by PITFALLS.md.

**Granularity note:** Configured granularity is `standard` (5–8 phases). At 9 phases this roadmap is at the upper edge of that band. Considered compressions were:

- **Merge Phase 3 (text/colormap) into Phase 2 (drawing)** — rejected. The colormap's dark-mode-freeze invariant and the bundled-font reproducibility work are a distinct concern from pure geometry primitives, and font-metric drift is a named pitfall. Keeping them separate isolates the debugging surface.
- **Merge Phase 4 (pixmap save/restore) into Phase 2 (drawing)** — rejected. Save/restore cannot be validated without interactive mesher sessions, which need the drawing primitives already stable. A separate phase gives a clean A/B checkpoint before the event-bridge pivot.
- **Merge Phase 8 (A/B validation) into Phase 7** — rejected. Phase 8 includes HiDPI and `_OMP` sweeps that are orthogonal to export work and is the explicit ship gate for the one-release-cycle window.

The 9-phase layout is accepted with explicit justification: each phase matches one ARCHITECTURE.md component boundary, each leaves the tree buildable and A/B-validatable, and each has a crisp definition of done that fits a solo developer's rollback granularity.

---
*Roadmap created: 2026-04-10*
