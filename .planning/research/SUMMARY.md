# Project Research Summary

**Project:** xvue-qt
**Domain:** Qt 6 / C++ desktop GUI replacing an X11/Xlib graphics layer driven synchronously by Fortran 77/95
**Researched:** 2026-04-10
**Confidence:** HIGH overall (HIGH on stack, ABI, build integration, table-stakes UX; MEDIUM on event-loop nesting behaviour, image-export plugin availability, and per-module lexicon audit cost)

## Executive Summary

xvue-qt replaces the 3619-line `xvue/xvuelc.c` Xlib backend with a Qt 6 C++ implementation behind a **byte-identical `extern "C"` ABI** (~60 entry points, trailing-underscore gfortran mangling, all arguments passed by pointer). The Fortran solver tree (`mail/`, `elas/`, `flui/`, `ther/`, `nlse/`, `reso/`, `util/`) is not touched — the migration's core invariant is "the Fortran side must not notice the change." On top of that frozen ABI, the project also delivers a full Level-3 UX modernization: a real `QMainWindow` with `QMenuBar`, `QToolBar`, `QDialog`, status bar, docks, HiDPI/dark-mode chrome, PNG/PDF/GIF export via Qt, and mouse-driven pan/zoom/pick — replacing the current text-lexicon menu system while still dispatching to the same Fortran subroutines.

The recommended approach is a **retained-mode backing-pixmap architecture**: all drawing primitives paint into a persistent off-screen `QPixmap` owned by an `XvueState` object, and `QWidget::paintEvent` does nothing except blit that pixmap. Blocking calls (`xvsouris`, `xvpause`) are resolved with a **local nested `QEventLoop`** — `QApplication::exec()` is never called at top level, preserving Fortran's synchronous-subroutine calling model. The Qt layer is built as a single static library (`libxvueqt.a`) by a scoped `xvue/CMakeLists.txt` (CMake ≥ 3.21, C++17, `AUTOMOC` on, `-fPIC`, no `-fopenmp`), consumed by the existing shell-script Fortran linker via `pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport`. The legacy X11 backend stays alive in parallel for one release cycle for A/B validation.

The dominant risks are all ABI-and-event-loop issues concentrated at the start: name mangling, pointer-vs-value argument types, `XPoint` layout (`short x, y`), `QApplication` singleton lifetime, and the `QApplication::exec()` trap that would invert control flow. A single developer with no CI must also guard against regression drift (Pitfall 20) with a 5-canonical-`testa/`-case end-of-phase checklist. Secondary risks — font metric drift, indexed-colormap quantization differences, HiDPI double-scaling, GIF-writer plugin availability on Debian trixie — are handled by explicit probes and by preserving the hand-rolled PostScript emitter verbatim.

## Key Findings

### Recommended Stack

Debian trixie supplies everything from apt; no vendored builds, no source compilation of Qt. The Qt layer is strictly single-threaded on the main thread; OpenMP stays on the Fortran side only.

**Core technologies:**
- **Qt 6.10.2** (`qt6-base-dev`, `qt6-base-dev-tools`) — windowing, painting, events, image I/O, printing — distro-pinned, Qt 5 explicitly rejected as maintenance-mode
- **Qt6::Widgets + Qt6::Gui + Qt6::Core + Qt6::PrintSupport** — one-to-one mapping from `xvuelc.c` primitives onto `QPainter`/`QPixmap`/`QImage`/`QFontMetrics`, plus `QMainWindow`/`QMenuBar`/`QToolBar`/`QDialog` chrome; Qt Quick/QML rejected for this dense desktop scientific app
- **C++17** (`-std=c++17`) — Qt 6.10's official floor, no justification to move to C++20 for a single-developer project
- **CMake ≥ 3.21** (trixie ships 3.31+) — first release with robust Qt 6 `find_package` / `AUTOMOC` support, scoped to `xvue/` only (the rest of the Fortran build stays shell-based)
- **gcc/g++ 13 + gfortran 13** (trixie defaults) — must match majors so `libstdc++` and `libgfortran` are ABI-compatible on the final link line; do not mix compiler versions
- **`libqt6imageformats6-plugins`** — needed at runtime for GIF/TIFF/WebP; PNG/JPEG/BMP/PPM are built into `Qt6::Gui`

Full details: `.planning/research/STACK.md`.

### Expected Features

`.planning/research/FEATURES.md` — benchmarked against Gmsh 4.x, ParaView 5.x, Salome 9.x, and FreeCAD FEM.

**Must have (table stakes — v1 / P1):**
- Qt reimplementation of all ~60 `extern "C"` entry points with byte-identical names/signatures — non-negotiable
- `XvueApp`/`XvueWindow`/`XvueCanvas` shell with backing `QPixmap` and instance-scoped state
- `QMainWindow` + `QMenuBar` + `QToolBar` + `QAction` (shared between menu and toolbar), keyboard shortcuts, tooltips, status bar
- `QFileDialog` project open/save/export (respecting `$MEFISTOX/<project>/`), recent projects, `QSettings` geometry/layout persistence
- Context menus, About dialog, Help → `doc/` PDF launcher, bilingual (reuse existing `$MEFISTO/td/m/anglais` flag)
- Console/log dock with stdout pipe-reader from Fortran (keystone: also enables error dialogs and progress reporting)
- PNG/JPEG via `QImage::save()`, PDF via `QPrinter::PdfFormat`, animated GIF via `QImageWriter` loop — drops ImageMagick runtime dep
- Pan / zoom / rotate via mouse wheel + middle-drag (`QTransform`), mouse-coordinate readout in status bar, hover picking feedback on top of existing `xvsouris*`
- HiDPI via `devicePixelRatioF`; dark-mode chrome while keeping scientific colormaps frozen (colors encode physical data)
- Anti-aliased `QPainter` rendering (free win)
- Visual A/B on selected `testa/` canonical cases (process, not feature)

**Should have (v1.x differentiators — P2, post-A/B):**
- Crash-recovery lockfile + "Recover last project?" prompt on startup (mitigates CONCERNS.md STOP issue)
- Command log dock showing the lexicon-equivalent of each GUI action; save/replay as scripted session (the "scripting" answer without adding a scripting language)
- Scientific colormap presets (viridis, plasma, grayscale; rainbow stays as default for A/B)
- Color bar widget with log-scale toggle
- Mesh statistics dock (element count, quality histogram, edge lengths)
- Per-module saved workspace layouts
- Progress reporting parsed from solver stdout via the pipe-reader
- Animation scrubber for time-stepping solvers (`flui`, `ther`, `nlse`)
- Single-level snapshot-based undo for mesh edits (scratch-file snapshot, zero Fortran change)

**Defer (v2+ / P3):**
- Multiple viewports / split view (depends on instance-scoped `xvuelc` state proving rock-solid)
- Mesh-quality color overlay (depends on per-element metric exposure from `mail/`)
- 3D clip planes / section cuts (requires `QOpenGLWidget` rewrite of the display path)
- Live-switch bilingual UI via `QTranslator` (no-restart)
- VTK / CSV export for ParaView post-processing
- Cooperative solver cancellation (requires Fortran-side changes — never, unless Core Value is renegotiated)

**Explicit anti-features (kept out of scope forever):**
- Custom theming engine, plugin system, embedded scripting language (Python/Lua/QtScript), cloud sync, real-time collaboration, auto-updater, crash-reporter upload, telemetry, integrated code editor, OpenCASCADE / built-in CAD kernel, MDI multi-document interface, Qt re-skinning of bash launcher scripts, rewriting the text-lexicon parser, integrated test runner for `testa/`

### Architecture Approach

Details in `.planning/research/ARCHITECTURE.md`. The design reconciles Fortran's synchronous-blocking calling model with Qt's event-loop-owns-main-thread paradigm via two non-negotiable pattern choices: (1) a **persistent off-screen backing `QPixmap`** owned by the library with one long-lived `QPainter` bound to it, and (2) a **nested local `QEventLoop` per blocking entry point** (never a top-level `QApplication::exec()`).

**Major components (all in `xvue/`, one concern per file, ~12 files total):**

1. **`XvueAbiShim`** (`xvue_qt_api.cpp`) — ~60 thin `extern "C"` entry points (≤10 lines each) with byte-identical names and pointer-based signatures; no Qt types cross this boundary, every body is `try`-wrapped so exceptions cannot propagate into Fortran
2. **`XvueApp`** — process-lifetime singleton owning the single `QApplication`, created lazily via `std::call_once` with fabricated static `argc`/`argv`, torn down via `atexit` (never in `xvfermer`)
3. **`XvueWindow`** — `QMainWindow` with menu bar, toolbar, status bar, dock widgets; created in `xvinitgraphique`, destroyed in `xvfermer`
4. **`XvueCanvas`** — `QWidget` whose `paintEvent` is a single `drawPixmap(0, 0, backing)` blit; resize reallocates the backing pixmap preserving old content
5. **`XvueState`** — pen/brush/font/palette + the one persistent `QPainter*` bound to the backing pixmap; equivalent to Xlib's `GC`
6. **`XvuePixmapStack`** — named off-screen slots backing `fenetremempx`/`mempxfenetre`/`sauvefenetre`/`restaurefenetre`/`sauvemempx`/`restauremempx`/`effacemempx`
7. **`XvueEventBridge`** — `QObject::eventFilter` on the canvas; `waitForEvent()` runs a local `QEventLoop` with a `targetLoop_` guard for re-entrancy; mouse-move coalescing via deferred `QTimer::singleShot(0, quit)` to mirror `XEventsQueued` semantics
8. **`XvueMenuBridge`** — `QAction::triggered` **queues** lexicon command strings; the next `xvsouris` call drains the queue and synthesises a `notypeevent=2` keyboard return — menu clicks look to Fortran like the user typed the command (critical to avoid re-entrancy Pitfall 8)
9. **`XvueExport`** — PNG/JPEG via `QImageWriter`, animated GIF via per-frame write, PostScript via the **preserved verbatim hand-rolled emitter** (Pitfall 13 — not `QPrinter` PostScript); PDF is an additive bonus via `QPrinter::PdfFormat`

The shared header `xvue_qt_api.h` carries the `extern "C"` declarations; both `xvuelc.c` (X11) and `xvue_qt_api.cpp` (Qt) include it; a build-time switch (`bin/cbl_tout_qt` variant) picks which object ends up on the link line. No parent `CMakeLists.txt`.

### Critical Pitfalls

Top items from `.planning/research/PITFALLS.md` (20 pitfalls total, severity-weighted):

1. **Trailing-underscore name mangling mismatch (Phase 0)** — gfortran emits `xvtrait_`; preserve `#define proc(x) x##_` from `xvuelc.c` verbatim and add an `nm libxvueqt.a` build-time check. Forbid `-fno-underscoring`.
2. **Fortran pass-by-reference pointer signatures (Phase 0)** — every scalar is `int*`/`float*`, every string is `char* + int*` length pair. Copy `xvuelc.c` signatures literally.
3. **`XPoint` is two shorts (Phase 2)** — define `struct MefistoPoint { short x; short y; };` in the bridge, convert to `QPoint` only at `QPainter::drawPolyline` call sites. Affects `xvface_`, `xvtraits_`, `xvfacetraits_`.
4. **`QApplication::exec()` at top level inverts control flow (Phase 1)** — never call it; blocking entry points use a nested local `QEventLoop` whose filter calls `loop.quit()` on the awaited event. Pre-commit grep that fails on any `QApplication::exec` in `xvue/*.cpp`.
5. **`QApplication` double-initialization (Phase 1)** — guard construction with `std::call_once` and static fabricated `argc`/`argv`; **never destroy the `QApplication` in `xvfermer_`** (only in `atexit`).
6. **Modal `QDialog::exec()` nested inside `xvsouris` re-entrancy (Phase 5/6)** — use the "queue, do not execute" pattern: `QAction::triggered` pushes a lexicon command string into `XvueMenuBridge::pendingCommands_`; `xvsouris` drains it and synthesises a keyboard return. Enforce with an `XvueApp::blockingDepth()` counter.
7. **PostScript fidelity drift (Phase 7)** — keep the existing hand-rolled `xvpostscript` emitter verbatim in `xvue_qt_postscript.cpp`. Do **not** switch to `QPrinter` PostScript (Qt 6 reduced PS support, produces byte-different output). PDF is an additive bonus only.
8. **`-fPIC` missing on `libxvueqt.a` (Phase 0)** — Debian trixie defaults to PIE; `set(CMAKE_POSITION_INDEPENDENT_CODE ON)` at the top of `CMakeLists.txt`.
9. **`-fopenmp` collision on `_OMP` executables (Phase 0 + Phase 8)** — do **not** pass `-fopenmp` to the Qt library; add `target_compile_options(xvueqt PRIVATE -fno-openmp)`. Document: all graphics calls happen on the Fortran main thread only. Debug assert `QThread::currentThread() == qApp->thread()` at every entry point.
10. **Regression drift in a single-developer no-CI project (every phase)** — pick 5 canonical `testa/` cases at Phase 0 (one per module); run all 5 through **both** backends at the end of every phase; log to `.planning/phase-N/VALIDATION.md`. This is a process rule — the single biggest risk.

Additional medium-severity items handled in their phases: AUTOMOC missing vtables (Phase 0), GIF writer plugin availability (Phase 7 probe), indexed-colormap → RGBA drift (Phases 3+8), font metric drift (Phases 3+8), Y-up vs Y-down convention audit (Phase 0), HiDPI double-scaling (Phases 1+8), stale `.o` files after `incl/*.inc` edits (Phase 0), ImageMagick removal audit (Phase 9).

## Implications for Roadmap

STACK.md and ARCHITECTURE.md both propose a ~10-phase layout. Reconciled into **nine phases** with a strict dependency chain (0 → 1 → {2 ∥ 3} → 4 → 5 → 6 → 7 → 8 → 9). Each phase lands one well-defined layer of the backing-pixmap + nested-event-loop architecture and leaves the tree in a buildable, A/B-validatable state.

### Phase 0: CMake skeleton + empty ABI stubs
**Rationale:** Proves the build plumbing (CMake + Qt + `-fPIC` + `AUTOMOC` + pkg-config-based link into the existing shell linker) before any graphics logic exists. Disproportionately important because it catches nine of the twenty pitfalls by construction.
**Delivers:** `xvue/CMakeLists.txt`, `xvue_qt_api.{h,cpp}` with all ~60 symbols as no-op stubs, `libxvueqt.a`, `bin/cbl_tout_qt` variant that cleans and rebuilds, one test executable (e.g. `pp/ppmail_qt`) linking successfully. 5-case validation baseline captured. X11 build still works unchanged.
**Addresses:** Build-system-owns-only-`xvue/` constraint; byte-identical ABI invariant.
**Avoids:** Pitfalls 1 (mangling), 2 (pointer args), 5 (singleton plumbing), 9 (AUTOMOC), 10 (`-fPIC`), 11 (`-fopenmp`), 16 (Y-axis convention), 18 (clean-build discipline).

### Phase 1: `XvueApp` + `XvueWindow` + `XvueCanvas` shell
**Rationale:** A window must exist before drawing can be tested. Exercises the Qt singleton discipline (`std::call_once`, fabricated `argc`/`argv`, `atexit`) in isolation.
**Delivers:** `xvinitgraphique` opens an empty `QMainWindow` with a blank central `XvueCanvas`. `xvfermer` closes it without destroying the `QApplication`. `xvpxecran`/`xvmmecran` return screen metrics via `QScreen` (logical pixels). `xvfond` sets background. HiDPI convention documented in `xvue/README_COORDS.md`.
**Uses:** `Qt6::Widgets` + `Qt6::Core`.
**Avoids:** Pitfalls 5 (double init), 6 (never `exec()`), 17 (HiDPI convention up front).

### Phase 2: Drawing primitives + backing pixmap + `XvueState`
**Rationale:** Roughly half the ABI is pure synchronous drawing. Independently testable via the existing `prpr/xvtest[1-4].f` programs — no event handling needed yet. One persistent `QPainter` on the backing pixmap is the central invariant established here.
**Delivers:** `XvueState`, `XvuePixmapStack` (partial), all pure draw entry points: `xvtrait`, `xvftrait`, `xvtraits`, `xvface`, `xvfacetraits`, `xvrectangle`, `xvbordrectangle`, `xvfrectangle`, `xvbordarcellipse`, `xvarcellipse`, `xvcouleur`, `xvtypetrait`, `xvepaisseur`, `effacer`, `xvvoir`, `xvpxfenetre`. Anti-aliasing enabled.
**Uses:** `QPainter`, `QPixmap`, `QTransform`.
**Avoids:** Pitfalls 3 (string lengths), 4 (`MefistoPoint` shim), and the per-call `QPainter` anti-pattern.

### Phase 3: Text, fonts, colormap
**Rationale:** Independent of event handling; can be done in parallel with Phase 2 but the colormap is used by drawing primitives. Colors encode physical data — must be frozen against system dark-mode.
**Delivers:** `xvchargefonte`, `xvnbpixeltexte`, `xvtexte`, `xvftexte`, `xvCouleursImposees`, `xvStockeRGBtoColormap`, `xvColormapToRGB`, `xvrecuprgbdec`, `xvactivervb`. Internal `std::array<QColor, MAX_PALETTE>` mirrors the X11 colormap semantically. Bundled fixed font for reproducibility.
**Uses:** `QFontMetrics`, `QColor`, `QFont`.
**Avoids:** Pitfalls 14 (colormap drift), 15 (font metrics).

### Phase 4: Pixmap save/restore (double-buffering)
**Rationale:** Depends on Phase 2. Implements the double-buffering the mesher already uses heavily; needed before interactive mesher sessions can be validated visually.
**Delivers:** `fenetremempx`, `mempxfenetre`, `sauvefenetre`, `restaurefenetre`, `sauvemempx`, `restauremempx`, `effacemempx`.
**Uses:** `QPixmap::drawPixmap`, `XvuePixmapStack` named slots.

### Phase 5: Event bridge + blocking reads — THE architectural pivot
**Rationale:** Once this lands, the mesher becomes interactive. Depends on Phases 1–4. This is the phase where the nested `QEventLoop` pattern proves itself against real `testa/` cases; expect one round of empirical iteration on mouse-motion coalescing.
**Delivers:** `XvueEventBridge` with `QObject::eventFilter`, `waitForEvent()` wrapping a local `QEventLoop`, `targetLoop_` re-entrancy guard, deferred `singleShot(0, quit)` for mouse coalescing. Entry points: `xvsouris`, `xvsouris2`, `xvpause`, `deplsouris`, `initaccrochage`.
**Uses:** `QEventLoop`, `QObject::installEventFilter`, `QTimer::singleShot`.
**Avoids:** Pitfalls 6 (no top-level `exec()`), 7 (`processEvents` starvation), 8 (nested modal re-entrancy — bridge built here).

### Phase 6: Menu / toolbar surface (Level 3 UX)
**Rationale:** Depends on Phase 5 because menu clicks go through the same event-bridge plumbing. **Long pole is the per-module lexicon audit** (mail, elas, flui, ther, nlse) — itself research-grade work. Can be split into sub-phases per solver module.
**Delivers:** `XvueMenuBridge`, `QMenuBar`/`QToolBar`/`QAction` definitions per module, shared `QAction` between menu and toolbar, keyboard shortcuts, tooltips, context menus, About dialog with Perronnet/LJLL credit, status bar, `QFileDialog` with recent-projects list, `QSettings` window/layout persistence, bilingual UI using existing flag, `doc/` PDF launcher on F1, console/log dock with stdout pipe-reader, error dialogs parsed from stdout, dark-mode chrome. Pan/zoom/rotate via mouse wheel + middle-drag, coordinate readout, hover picking feedback.
**Uses:** `QMainWindow`, `QMenuBar`, `QToolBar`, `QAction`, `QDialog`, `QDockWidget`, `QSettings`, `QFileDialog`, `QMessageBox`, `QProcess` (pipe reader).
**Avoids:** Pitfall 8 (queue-do-not-execute `QAction` pattern; `blockingDepth()` gate on modals).

### Phase 7: Image / GIF / PostScript export
**Rationale:** Independent of event handling; defer to keep early phases focused on interactive correctness. **First task is a 10-line `QImageWriter::supportedImageFormats()` probe** to confirm GIF writer availability on trixie — if absent, fall back to per-frame PNG + `ffmpeg` (never ImageMagick).
**Delivers:** `XvueExport`: PNG/JPEG via `QImageWriter`, animated GIF via per-frame write loop (or ffmpeg fallback), PostScript via the **preserved hand-rolled emitter** moved verbatim from `xvuelc.c` into `xvue_qt_postscript.cpp`. Additive PDF bonus via `QPrinter::PdfFormat` under a new entry point (not by modifying `xvpostscript`). Drops the ImageMagick runtime dep.
**Uses:** `QImageWriter`, `QPrinter`, `QImage`.
**Avoids:** Pitfalls 12 (GIF plugin probe first), 13 (PostScript verbatim).

### Phase 8: A/B validation across the `testa/` subset
**Rationale:** Gate for declaring the migration done. Treat as "compare, tune, re-compare" with 1–2 iterations expected, mostly around colormap, font, HiDPI drift. Must validate on HiDPI (4K or `QT_SCALE_FACTOR=2`) and on `ELASTICER_OMP`.
**Delivers:** Run canonical cases per module (mesher, elas, flui, ther, nlse) through both backends. Visual side-by-side comparison. Validation log in `.planning/phase-8/VALIDATION.md`.
**Process, not code.**

### Phase 9: Retire X11 backend
**Rationale:** Can only run after the one-release-cycle A/B window closes. Trailing cleanup.
**Delivers:** Delete `xvuelc.c`, `bin/ccxvue`, `libX11` linker lines from `bin/cb*`, `bin/convertepsgif`. Audit all ImageMagick `convert` call sites (Pitfall 19) under `bin/`, `td/`, `testa/`, `testf/` before dropping the dependency. Update `README`, `LISEZMOI`, install scripts. Remove the hardcoded `/usr/X11R6/lib64` path.
**Avoids:** Pitfall 19 (audit before removal).

### Phase Ordering Rationale

- **Why this order:** The dependency graph is near-linear. The ABI bridge (0) gates everything. A window (1) must exist before drawing (2, 3, 4). Drawing must work before blocking event reads (5) can be tested meaningfully. Menus (6) depend on the event bridge because `QAction` triggers and `xvsouris` share the same plumbing. Export (7) is the only truly independent slot. A/B (8) and retirement (9) are trailing gates.
- **Why this grouping:** Each phase matches one component boundary from ARCHITECTURE.md — no phase crosses two components. Keeps commits small and each phase's "definition of done" crisp for a solo developer.
- **How this avoids pitfalls:** Phase 0 front-loads nine of the twenty pitfalls. Phase 1 isolates `QApplication` lifetime bugs before drawing obscures them. Phase 5 isolates the event-loop design at the moment it can be empirically validated against `testa/`. Phase 8 is the gate where color/font/HiDPI drift surfaces and can be tuned before retirement.

### Research Flags

Phases likely needing deeper research during planning (`/gsd-research-phase`):

- **Phase 5 (event bridge)** — MEDIUM confidence in the nested `QEventLoop` + mouse-motion-coalescing design. Empirical validation against a real mesher session required.
- **Phase 6 (menu surface)** — the per-module lexicon audit IS the research. Each interactive driver in `mail/`, `elas/`, `flui/`, `ther/`, `nlse/` must be enumerated. Expect this phase to split into 5 sub-phases.
- **Phase 7 (export)** — MEDIUM confidence on GIF writer plugin availability. Run the probe at phase kickoff before committing to `QImageWriter` vs `ffmpeg` fallback.

Phases with standard patterns (skip research, follow STACK.md / ARCHITECTURE.md directly):

- **Phases 0, 1, 2, 3, 4** — CMake + Qt 6 + `QPainter`/`QPixmap` are well-documented; STACK.md contains the copy-pasteable CMake template and exact ABI rules.
- **Phases 8, 9** — process/cleanup, checklist-driven.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | Direct reads of `xvue/xvuelc.c` and `bin/ccxvue` confirm mangling, flags, and entry-point pattern. Qt 6.10 / CMake 3.21 / C++17 pinned by distro and by official Qt 6 docs. Only long-term event-loop integration is MEDIUM. |
| Features | HIGH for table stakes and anti-features (well-established Qt/desktop conventions grounded in PROJECT.md constraints); MEDIUM for differentiators (Gmsh/ParaView/Salome are authoritative but trade-offs depend on MEFISTO workflows revealed only during the per-module lexicon audit). |
| Architecture | HIGH for Qt core patterns (`QApplication`/`QPainter`/`QPixmap` — stable since Qt 4); MEDIUM for nested-`QEventLoop` discipline around `QDialog::exec` re-entrancy (queue-don't-execute pattern mitigates but needs empirical validation); LOW on `xvpostscript` strategy (preserve verbatim is conservative). |
| Pitfalls | HIGH on ABI and build-integration traps (direct source reads); MEDIUM on event-loop re-entrancy and image-export parity (depend on empirical behaviour and on the trixie `qt6-imageformats` plugin bundle). |

**Overall confidence:** HIGH. The ABI surface, build integration, stack, and architectural patterns are all either directly verified against the source tree or pinned by authoritative Qt 6 documentation. The remaining MEDIUM areas are concentrated in Phases 5, 6, 7 and have explicit mitigation paths.

### Gaps to Address

- **`xvuelc.c` Y-axis convention (Pitfall 16).** One-hour grep + read during Phase 0 to document before the Qt bridge commits.
- **Debian trixie `qt6-imageformats` GIF writer plugin availability (Pitfall 12).** 10-line `QImageWriter::supportedImageFormats()` probe at Phase 7 kickoff. Fallback: per-frame PNG + `ffmpeg`.
- **Nested modal re-entrancy under the queue-do-not-execute pattern (Pitfall 8).** Empirical test during Phase 5 before Phase 6 builds on top.
- **`XPoint` call sites in `xvue/*.f` (Pitfall 4).** Audit the three wrappers (`XVFACE`, `XVTRAITS`, `XVFACETRAITS`) before Phase 2. Audit is small (three files).
- **Any `!$OMP PARALLEL` region calling `CALL XV*` in `_OMP` executables (Pitfall 11).** Grep `elas/`, `flui/`, `ther/` for `!$OMP` near `CALL XV*`. Expected none, verify before making the rule normative.
- **User preference on `.eps` vs `.pdf` as archival export format (Pitfall 13).** Recommendation: preserve EPS verbatim in v1, add PDF as additive bonus. Confirm before Phase 7.
- **`xvnbpixeltexte_` fallback behaviour for missing fonts.** Determines whether Qt version must carry a fallback.
- **Frame-selection mechanism for time-stepping solvers.** v1.x animation scrubber assumes `flui`, `ther`, `nlse` expose per-timestep frame selection through the text lexicon. Confirm in Phase 6 per-module audit.

## Sources

### Primary (HIGH confidence — direct source reads or pinned distro packages)
- `xvue/xvuelc.c` lines 1–70 and entry-point table — trailing-underscore ABI, `-fPIC -O -m64` flags, ~60 `extern "C"` entry points, `XPoint` short-short layout, existing hand-rolled PostScript emitter
- `bin/ccxvue` — current gcc flags and output location
- `/home/drico/git/mefisto/.planning/PROJECT.md` — Core Value, constraints, Out of Scope, Level 3 modernization target, Key Decisions
- `/home/drico/git/mefisto/.planning/codebase/` (`ARCHITECTURE.md`, `STRUCTURE.md`, `CONVENTIONS.md`, `TESTING.md`, `CONCERNS.md`) — module layout, text-lexicon dispatch pattern, manual A/B workflow, STOP-on-error and stale-`.o` concerns
- `/home/drico/git/mefisto/CLAUDE.md` — working rules, `99;` save-exit convention, Qt migration as active goal
- Debian trixie `qt6-base-dev` 6.10.2, `qt6-base-dev-tools`, `libqt6imageformats6-plugins`, CMake 3.31+ — version pins confirmed by apt
- gfortran Linux x86_64 mixed-language documentation — trailing-underscore, pass-by-reference, explicit length for `CHARACTER` strings
- Qt 6 upstream documentation (doc.qt.io/qt-6/) — `find_package(Qt6 COMPONENTS ...)`, `CMAKE_AUTOMOC`, `QPainter` lifetime on `QPixmap`, `QEventLoop::exec()` nesting, HiDPI / `devicePixelRatioF`

### Secondary (MEDIUM confidence)
- Qt forums and qt-project mailing list on embedding Qt in non-Qt main — fabricated `argc`/`argv`, `std::call_once`, `atexit` singleton pattern
- Gmsh 4.x, ParaView 5.x, Salome 9.x, FreeCAD FEM — reference points for scientific-FE GUI chrome
- `processEvents(WaitForMoreEvents)` and nested `QEventLoop` pattern in library-mode Qt — documented idiom with known re-entrancy gotchas

### Tertiary (LOW confidence — needs empirical validation)
- `QPrinter` PostScript output fidelity in Qt 6.10 vs Qt 5 — Qt 6 reduced PS support; hence "keep the hand-rolled emitter verbatim"
- Trixie `qt6-imageformats` GIF writer availability — needs Phase 7 probe
- Exact mouse-motion coalescing parity between `XEventsQueued` and `QTimer::singleShot(0, quit)` — needs Phase 5 empirical validation

---
*Research completed: 2026-04-10*
*Ready for roadmap: yes*
