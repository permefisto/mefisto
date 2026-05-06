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
- [ ] **Phase 6.0: Shared shell, menu bridge, dialogs, persistence** - `QMainWindow` chrome + `XvueMenuBridge` + `XvueConsoleDock` + `QSettings` + bilingual i18n + `xvue_module_init_` ABI hook + shared `{File, View, Help}` menus
- [x] **Phase 6.1: Mesher (mail) menu wiring** - `registerMailActions()` + `LEXICON-AUDIT-mail.md` + mesh SVG icons + `CALL XVUE_MODULE_INIT('mail')` (completed 2026-04-22)
- [ ] **Phase 6.2: Elasticity (elas) menu wiring** - `registerElasActions()` + `LEXICON-AUDIT-elas.md` + elas SVG icons + `CALL XVUE_MODULE_INIT('elas')`
- [ ] **Phase 6.3: Fluid (flui) menu wiring** - `registerFluiActions()` + `LEXICON-AUDIT-flui.md` + flui SVG icons + `CALL XVUE_MODULE_INIT('flui')`
- [ ] **Phase 6.4: Thermal (ther) menu wiring** - `registerTherActions()` + `LEXICON-AUDIT-ther.md` + ther SVG icons + `CALL XVUE_MODULE_INIT('ther')`
- [ ] **Phase 6.5: Nonlinear (nlse) menu wiring** - `registerNlseActions()` + `LEXICON-AUDIT-nlse.md` + nlse SVG icons + `CALL XVUE_MODULE_INIT('nlse')`
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

### Phase 02.1: Qt drawing-primitive A/B gaps (xvface color state, multi-object 3D composition, dashed pen style) (INSERTED)

**Goal:** Fix the single Qt-backend root cause (`xvfacetraits_` dropping its `ncf`/`nca` color arguments via `(void)ncf; (void)nca;`) that blocks the Phase 3 A/B visual gate on xvtest2 and xvtest4. Research (02.1-RESEARCH.md) identified this as the only real bug; the "multi-object 3D" and "dashed pen" symptoms are downstream artifacts of the same color-drop and resolve automatically once `ncf`/`nca` are honored. Legacy X11 and Fortran wrappers stay bit-identical (BUILD-07 / VALID-02 / ABI freeze at 57).
**Requirements**: BUILD-07, VALID-02 (preservation guards — no new IDs)
**Depends on:** Phase 2
**Plans:** 1/1 plans complete

Plans:
- [x] 02.1-01-PLAN.md — Fix xvfacetraits_ ncf/nca color switch + A/B re-capture (single wave, single file: xvue/qt/src/xvue_qt_api.cpp)

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
  - [x] 03-04-PLAN.md — Wave 3: A/B catch-up gate vs prpr/xvtest1..4 + 5 canonical testa cases (D-26 hard gate) + 03-VALIDATION.md closure
**UI hint**: yes

### Phase 03.1: Build xvtest1..4 driver infrastructure on both backends (Qt + legacy X11) to unblock Phase 3 A/B catch-up gate (INSERTED)

**Goal:** Produce pp/ppxvtest{0..4}(_qt) on both Qt and legacy X11 backends (9 new cb* scripts + 2 top-level edits), fix the Qt-side xtinit_ warn-once stub so XVOUVRIR opens a real window, and harden plan 03-04 Task 1 smoke (timeout 5s, enriched grep including xtinit_, xvtest1 STOP tolerance) so the Phase 3 A/B visual gate can execute.
**Requirements**: BUILD-07, VALID-02, TEXT-01..06 (inherited — no new IDs)
**Depends on:** Phase 3
**Plans:** 4/4 plans complete

Plans:
- [x] 03.1-01-PLAN.md — Wave 0: A6 grep verification + Qt-side xtinit_ promotion (xvue/qt/src/xvue_qt_api.cpp), X11 untouched
- [x] 03.1-02-PLAN.md — Wave 1: 9 new cb* scripts (5 legacy + 4 Qt) + top-level wiring into cbl_tout{,_qt}
- [x] 03.1-03-PLAN.md — Wave 2: hardened headless smoke + ABI recheck + 03-04 Task 1 patch + NOTES.md + SUMMARY.md

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
**Plans**: 6 plans
  - [x] 05-01-PLAN.md — Wave 0 scaffolding: Qt test infrastructure + blockingDepth counter + XvueState mempxaccro fields + canvas StrongFocus (EVENT-01, EVENT-06, EVENT-08 infra)
  - [x] 05-02-PLAN.md — XvueEventBridge class + BlockingDepthGuard RAII + button/key dispatch + install on canvas (EVENT-01, EVENT-08)
  - [x] 05-03-PLAN.md — Motion coalescing via QTimer::singleShot deferred-quit + diagnostic counter (EVENT-07)
  - [x] 05-04-PLAN.md — Real bodies xvsouris_/xvpause_/deplsouris_ + AUTOEXIT extension to xvpause_ both backends (EVENT-02, EVENT-04, EVENT-05)
  - [x] 05-05-PLAN.md — initaccrochage_ + xvsouris2_ Strategy B accrochage save/restore (EVENT-03, EVENT-06)
  - [x] 05-06-PLAN.md — Clean rebuild + manual A/B drag test (pan2d, torus) + VALIDATION.md sign-off + README Wayland caveat (EVENT-07 human)
**UI hint**: yes

### Phase 6.0: Shared shell, menu bridge, dialogs, persistence
**Goal**: A full Qt `QMainWindow` shell (menu bar mechanism, toolbar, status bar, dock widgets, dialogs, `QSettings` persistence, bilingual i18n, system dark-mode chrome, wheel/middle-drag/right-click canvas interactions) is in place and ready for per-module menu wiring. `XvueMenuBridge` queues lexicon commands for `xvsouris_` to drain as synthetic `notypeevent=2` events. Modal dialogs respect `XvueApp::blockingDepth()`. One new `extern "C"` entry point `xvue_module_init_` takes ABI 57→58.
**Depends on**: Phase 5
**Requirements**: UX-01, UX-02, UX-03, UX-04, UX-06, UX-07, UX-08, UX-09, UX-10, UX-11, UX-12, UX-13
**Success Criteria** (what must be TRUE):
  1. `pp/ppmail_qt` launches with a `QMenuBar` showing `{File, View, Help}` (no module-specific menus yet; those ship in 6.1..6.5), a `QToolBar` with shared actions, a `QStatusBar` with live coordinate readout, and a `QDockWidget` console receiving solver stdout.
  2. `XvueMenuBridge::pushCommand(QString)` adds chars to a queue drained by `XvueEventBridge::waitForEvent` via a pre-`exec()` hook, returning each character as `notypeevent=2, nbc=<ascii>` across successive `xvsouris_` calls.
  3. Modal `QFileDialog`/`QDialog` entry points refuse to open with a 3-second status-bar message when `XvueApp::blockingDepth() > 0`.
  4. Window geometry, dock layout, recent-projects list (≤10 entries), and UI preferences persist via `QSettings` across process restarts; `xvueIsEnglish()` reads `$MEFISTO/td/m/anglais` file-existence and all shared chrome renders in FR or EN accordingly. About dialog credits Alain Perronnet / LJLL / UPMC Paris.
  5. Console `QDockWidget` shows Fortran stdout via `freopen(stdout)` + `QSocketNotifier` in-process pipe; `*** ERREUR` lines surface as `QMessageBox::warning`. Canvas wheel=zoom, middle-drag=pan, right-click=context menu, live coord readout in status bar. System dark-mode via `QPalette` affects chrome only; canvas scientific colormaps untouched.
  6. `nm xvue/qt/build/libxvueqt.a | grep '_$' | wc -l` returns 58 (was 57 through Phase 5). `xvue_module_init_` stub accepts a module-name string argument and records it in `XvueMenuBridge`. Build-time lint scripts (`verify_shortcut_modifiers.sh`, `verify_icon_source.sh`) exist and exit 0 against the 6.0 tree.
**Plans**: 7 plans
  - [x] 06.0-01-PLAN.md — Wave 0 scaffold: 9 class stub pairs + xvue_module_init_ ABI stub (57→58) + 2 lint scripts + xvue_qt_menu_tests target + empty xvue_icons.qrc
  - [x] 06.0-02-PLAN.md — Wave 2: fill i18n 48-row table + XvueMenuBridge queueLexicon/popChar + XvuePrefs QSettings wrapper (UX-02, UX-06, UX-07, UX-08)
  - [x] 06.0-03-PLAN.md — Wave 3: XvueEventBridge::waitForEvent menu-queue pre-drain + MouseMove coord-signal emit + xvuelc.c A6 audit (UX-02, UX-03, UX-12)
  - [x] 06.0-04-PLAN.md — Wave 2: XvueConsoleDock pipe+QSocketNotifier + XvueErrorBatcher deferred QMessageBox + About/Preferences/RecentProjects dialog bodies (UX-04, UX-06, UX-09, UX-10, UX-11)
  - [x] 06.0-05-PLAN.md — Wave 2: XvueCanvas wheel/middle-drag/contextMenu/mouseCoords signal/resetView + XvueState view_transform_/has_user_content_ + empty-state paintEvent (UX-12, UX-13)
  - [x] 06.0-06-PLAN.md — Wave 4: XvueApp menuBridge accessor + applyColorSchemePreference + XvueWindow full chrome composition + 13 shared QActions + xvue_module_init_ dispatch body (UX-01, UX-04, UX-06, UX-07, UX-09, UX-13)
  - [x] 06.0-07-PLAN.md — Wave 5: clean rebuild + 7-binary test sweep + human A/B verdict on pp/ppmail_qt + 06.0-VALIDATION.md sign-off + xvue/qt/README.md extension
**UI hint**: yes

### Phase 6.1: Mesher (mail) menu wiring
**Goal**: `ppmail_qt` shows a fully-populated menu bar with mesher-specific top-level menus (`{File, Mesh, View, Help}` — Edit dropped per 06.1 D-01), toolbar buttons for the top-5 most-used leaf lexicon commands, and `CALL XVUE_MODULE_INIT('mail')` at startup. The typed lexicon continues to work unchanged for every command, GUI-surfaced or not.
**Depends on**: Phase 6.0
**Requirements**: UX-05 (mail slice)
**Success Criteria** (what must be TRUE):
  1. `.planning/phases/06.1-*/LEXICON-AUDIT-mail.md` enumerates every interactive command in `mail/` drivers with FR+EN descriptions, frequency column, and QAction flag per D-03 cutoff.
  2. `registerMailActions(XvueMenuBridge*, QMenuBar*, QToolBar*)` wires 20–40 QActions, each connected via `QAction::triggered` → `XvueMenuBridge::pushCommand(<command-string>)`.
  3. `prpr/ppmail.f` (or equivalent entry) contains `CALL XVUE_MODULE_INIT('mail')` before the interactive loop.
  4. Custom mesh SVG icons live under `xvue/qt/resources/icons/mail/` and are registered in `xvue_icons.qrc`; none triggers `verify_shortcut_modifiers.sh` or `verify_icon_source.sh` failure.
  5. A mesher-specific QTest case simulates a menu click and verifies the synthetic event reaches `xvsouris_` with the correct ASCII sequence.
**Plans**: 3 plans
  - [x] 06.1-01-PLAN.md — Full LIMTCL tree walk LEXICON-AUDIT-mail.md + tools/validate_audit_md.py + user review checkpoint to freeze frequency bucketing and top-5 toolbar
  - [x] 06.1-02-PLAN.md — registerMailActions_stub_ strong-symbol body + bilingual menu-file parser (D-12) + [menu] echo (D-07) + 10 custom SVG icons + xvue_icons.qrc append + CMake wiring
  - [x] 06.1-03-PLAN.md — xvue/xvmodi.f X11 no-op stub + prpr/ppmail.f CALL XVUE_MODULE_INIT + D-09 closeEvent/onFileQuit rewrite + 4 D-13 QTest cases + manual A/B sign-off
**UI hint**: no (inherits 6.0 contract)

### Phase 6.2: Elasticity (elas) menu wiring
**Goal**: `ppelas_qt` shows an elasticity-specific top-level menu bar (`{File, Solve, View, Help}`), 80/20 QAction coverage, and `CALL XVUE_MODULE_INIT('elas')`.
**Depends on**: Phase 6.1
**Requirements**: UX-05 (elas slice)
**Success Criteria**: same pattern as 6.1, substituting `elas` for `mail`.
**Plans**: 5 plans (3 implementation + 2 gap-closure)
  - [x] 06.2-01-PLAN.md — LEXICON-AUDIT-elas.md full tree walk (debuelas + 12 elas-unique sub-menus + 20 shared util compressed) + user review checkpoint (reuses tools/validate_audit_md.py from 6.1)
  - [x] 06.2-02-PLAN.md — registerElasActions_stub_ strong body + xvue_qt_elas_actions_keepalive + 7 custom SVG icons + xvue_icons.qrc append + CMake wiring + tools/validate_audit_md.py generalized (module-aware)
  - [x] 06.2-03-PLAN.md — prpr/ppelas.f CALL XVUE_MODULE_INIT('elas', 4) insertion (xvue/xvmodi.f reused from 6.1) + 4+1 D-13 QTest cases + manual A/B sign-off
  - [x] 06.2-04-PLAN.md — gap-closure: insert Solve/Mesh menus before View (gap-1 menu order) + testMenuOrder QTest slots in elas+mail + IN-01/IN-03/IN-04 cleanups
  - [x] 06.2-05-PLAN.md — gap-closure: XvueMenuFileParser consults td/ma/<name> when anglais flag set (gap-2 bilingual labels) + testBilingualLabelsEnglish QTest slot
**UI hint**: no (inherits 6.0 contract)

### Phase 6.3: Fluid (flui) menu wiring
**Goal**: `ppflui_qt` shows a fluid-specific top-level menu bar (`{File, Fluid, View, Help}`), 80/20 QAction coverage, and `CALL XVUE_MODULE_INIT('flui')`.
**Depends on**: Phase 6.2
**Requirements**: UX-05 (flui slice)
**Success Criteria**: same pattern as 6.1, substituting `flui` for `mail`. testMenuOrder + testBilingualLabelsEnglish codify the gates that required manual sign-off in 6.2/6.1 — Phase 6.3 is fully autonomous.
**Plans**: 3 plans
  - [x] 06.3-01-PLAN.md — LEXICON-AUDIT-flui.md full tree walk (debuflui 15 codes + 20 flui-unique sub-menus fully expanded + 21 shared util compressed) + user review checkpoint
  - [x] 06.3-02-PLAN.md — registerFluiActions_stub_ strong body + xvue_qt_flui_actions_keepalive + 6 custom SVG icons + xvue_icons.qrc append + CMake wiring (keepalive alphabetical: elas/flui/mail)
  - [x] 06.3-03-PLAN.md — prpr/ppflui.f CALL XVUE_MODULE_INIT('flui', 4) insertion + 7 QTest slots (D-13 + testMenuOrder + testBilingualLabelsEnglish) — fully autonomous, no manual A/B checkpoint
**UI hint**: no (inherits 6.0 contract)

### Phase 6.4: Thermal (ther) menu wiring
**Goal**: `ppther_qt` shows a thermal-specific top-level menu bar (`{File, Thermal, View, Help}`), 80/20 QAction coverage, and `CALL XVUE_MODULE_INIT('ther')`.
**Depends on**: Phase 6.3
**Requirements**: UX-05 (ther slice)
**Success Criteria**: same pattern as 6.1, substituting `ther` for `mail`.
**Plans**: TBD (~2-3 plans — same shape as 6.1)
**UI hint**: no

### Phase 6.5: Nonlinear (nlse) menu wiring
**Goal**: `ppnlse_qt` shows a nonlinear-specific top-level menu bar (`{File, Nonlinear, View, Help}`), 80/20 QAction coverage, and `CALL XVUE_MODULE_INIT('nlse')`.
**Depends on**: Phase 6.4
**Requirements**: UX-05 (nlse slice)
**Success Criteria**: same pattern as 6.1, substituting `nlse` for `mail`.
**Plans**: 3 plans
  - [x] 06.5-01-PLAN.md — LEXICON-AUDIT-nlse.md full LIMTCL tree walk + domain-review frequency bucketing + Help-allowlist {97;} hand-off (no .nlse fixtures exist; pure domain-review unlike ther's gourd evidence)
  - [x] 06.5-02-PLAN.md — registerNlseActions_stub_ strong body + xvue_qt_nlse_actions_keepalive (alphabetical: between mail and ther) + 5 custom SVG icons + xvue_icons.qrc append + CMake wiring
  - [x] 06.5-03-PLAN.md — prpr/ppnlse.f CALL XVUE_MODULE_INIT('nlse', 4) insertion (after single-language LOGO, follows flui pattern NOT ther bilingual) + 7 QTest slots (4 D-13 + closeEvent + testMenuOrder + testBilingualLabelsEnglish) — fully autonomous, no manual A/B checkpoint
**UI hint**: no

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
**Plans**: 6 plans
  - [x] 07-01-PLAN.md — QImageWriter probe binary + bin/cb_probe_qt + PROBE.md kickoff (EXPORT-01)
  - [x] 07-02-PLAN.md — PsEmitter skeleton + xvpostscript_ body verbatim port from xvuelc.c:1187-1304 (EXPORT-04)
  - [x] 07-03-PLAN.md — PsEmitter per-primitive helpers + 15 emit-site wirings in xvue_qt_api.cpp (EXPORT-04)
  - [x] 07-04-PLAN.md — XvueExport class (PNG/JPEG/PDF) + File→Export submenu + bilingual messages (EXPORT-02, EXPORT-05)
  - [x] 07-05-PLAN.md — Animated GIF (probe-driven dispatch + auto-snapshot + ffmpeg fallback + frame caps) (EXPORT-03)
  - [x] 07-06-PLAN.md — Golden tests + EXPORT-06 grep gate + README Phase 7 section + manual A/B sign-off (EXPORT-03, EXPORT-06) — autonomous portion complete; A/B sign-off deferred (see VALIDATION-LOG.md)

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
**Plans**: 7 plans
  - [ ] 08-01-PLAN.md — Wave 1 bootstrap: bin/cbl_tout_qt freshness (D-08, Phase 7 Gap-A) + 3 Phase-7 deferred goldens (D-06: scene01.eps + wave_legacy.gif + cavity2d_legacy.gif) + ctest QSKIP→PASS flip (D-07) + 3 sweep-harness scripts (bin/ab_compare_pair.sh, ab_capture_x11.sh, ab_sweep_phase8.sh) + 5-case AUTOEXIT smoke probes (VALID-01, VALID-02, VALID-07)
  - [ ] 08-02-PLAN.md — Wave 2 X11 baseline column: 5 canonical cases captured under Xvfb + OMP_NUM_THREADS=1 deterministic baseline (VALID-02)
  - [ ] 08-03-PLAN.md — Wave 2 Qt 1x column: 5 captures under QT_QPA_PLATFORM=offscreen + AE compare vs X11 baseline at fuzz=5% (VALID-01)
  - [ ] 08-04-PLAN.md — Wave 2 Qt HiDPI 2x column: 5 captures under QT_SCALE_FACTOR=2 + downsample-then-AE-compare with dimension guard (VALID-04)
  - [ ] 08-05-PLAN.md — Wave 2 OMP column: 4 OMP-eligible cases captured on both backends under OMP_NUM_THREADS=8 (cavity2d N-A — no FLUIDER_OMP per D-05) + main-thread guard verification (VALID-03)
  - [ ] 08-06-PLAN.md — Wave 3 spot checks: 3 colorbar reports (nafems_le1/heat1d/cavity2d) per Pitfall 6 + 2 font reports (pan2d/hexahedron) per Pitfall 7; hexahedron captured separately as VALID-06 spot-only (NOT extending BUILD-10 baseline) (VALID-05, VALID-06)
  - [ ] 08-07-PLAN.md — Wave 4 ship-gate: compose 08-CHECKLIST.md per-cell verdict matrix (D-10) + 08-VALIDATION.md VALID-NN coverage + maintainer sign-off checkpoint that opens the one-release-cycle A/B window (VALID-07)
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
| 8. A/B validation on testa subset | 7/7 | Complete (2026-05-05) | dricoco — v1 ship gate signed; 5 overrides accepted; Phase 9 unblocked |
| 9. Retire X11 backend | 9/9 | Complete (2026-05-06) | dricoco — xvuelc.c + libX11 + ImageMagick + LVIDEO retired; v1.0-pre-retire tag for rollback; 3 P7 goldens deferred (P7 source defects) |

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
*Phase 8 plans created: 2026-05-05*
