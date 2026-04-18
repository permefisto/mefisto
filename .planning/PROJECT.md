# xvue-qt

## What This Is

`xvue-qt` is the initiative to replace MEFISTO's legacy X11/Motif graphics layer (`xvue/xvuelc.c` plus its Fortran wrappers and text-lexicon menu system) with a Qt 6 GUI that preserves every existing MEFISTO operation while delivering a modernized user experience (real menu bar, toolbars, and dialogs). The underlying numerical solvers in `mail/`, `elas/`, `flui/`, `ther/`, `reso/`, `nlse/`, and `util/` are explicitly out of scope and stay untouched — this is a graphics-and-UI migration, not a scientific rewrite.

It is targeted at a single developer on Debian trixie, x86_64 only, and will ship alongside the existing X11 build for one release cycle so users can A/B compare before the X11 layer is retired.

## Core Value

Every MEFISTO workflow that works today through X11 — interactive meshing, elasticity/fluid/thermal/nonlinear solving, visualization, picking, plot export — keeps working through the new Qt 6 interface, with the Fortran solver code underneath completely unchanged.

If everything else about the UX is debatable, this invariant is not: the Fortran side must not notice the migration.

## Requirements

### Validated

<!-- Existing MEFISTO capabilities inferred from codebase map (.planning/codebase/). -->
<!-- These are "validated" in the sense that they exist and must keep working. -->

- ✓ 2D/3D interactive mesh generation via `MAILLER` (`pp/ppmail`, `mail/`) — existing
- ✓ Linear elasticity solver via `ELASTICER` / `ELASTICER_OMP` (`pp/ppelas`, `elas/`) — existing
- ✓ Fluid mechanics solver via `FLUIDER` (`pp/ppflui`, `flui/`) — existing
- ✓ Thermal solver via `THERMICER` / `HEATER` (`pp/ppther`, `ther/`) — existing
- ✓ Nonlinear solver via `NLSER` (`pp/ppnlse`, `nlse/`) — existing
- ✓ Project initialization via `INITIER` (`pp/ppinit`) — existing
- ✓ Shared linear-system solvers (`reso/`) — existing
- ✓ Fortran subroutine calls into the ~60 `extern "C"` graphics entry points (`xvue/xvuelc.c`: drawing primitives, colors, fonts, pixmap buffers, mouse events, PostScript export) — existing
- ✓ Text-lexicon menu system driving solver commands (the `99;` save-exit convention etc.) — existing
- ✓ Bilingual French/English interface selected at runtime via `$MEFISTO/td/m/anglais` — existing
- ✓ Shell-based build via `bin/cbl_tout` and per-module `bin/cb*` scripts — existing
- ✓ Manual end-to-end tests in `testa/` and `testf/` (~60 cases each) — existing
- ✓ Tutorials and demo data under `td/` — existing

### Active

<!-- v1 scope for xvue-qt. Hypotheses until proven by running the real testa/ cases. -->

- [ ] **Qt 6 reimplementation of `xvue/xvuelc.c`** — every `extern "C"` entry point currently in `xvuelc.c` (~60 functions: drawing primitives, color/font, pixmap buffers, event loop, PostScript export, etc.) is reimplemented in C++ using Qt 6 (`QPainter`, `QPixmap`, `QWidget`, `QImage`, `QApplication`, `QTimer`) with byte-identical names and signatures so `xvue/*.f` Fortran wrappers are untouched.
- [ ] **CMake build owning only `xvue/`** — a new `xvue/CMakeLists.txt` compiles the Qt/C++ sources (running `moc` automatically), links against `Qt6::Widgets`/`Qt6::Gui`/`Qt6::Core`, and produces a static library consumed by the existing shell linker in `bin/cb*`. The rest of the Fortran build stays shell-based and unchanged.
- [ ] **Image and animated-GIF export via Qt** — `QImage` for single-frame PNG/JPEG/EPS, `QMovie` or `QImageWriter`-loop for animated GIF. The ImageMagick `convert` shell-out (`bin/convertepsgif`) becomes unnecessary and is dropped.
- [ ] **Real Qt UI chrome replacing the text lexicon** — `QMainWindow` with `QMenuBar`, `QToolBar`, and `QDialog` entry points wired to the same Fortran solver subroutines today reached via typed menu commands. Coverage pass across `mail/`, `elas/`, `flui/`, `ther/`, `nlse/` menus. Keyboard shortcuts and a modern About dialog included.
- [ ] **Parallel X11 build kept alive for one release cycle** — the legacy `bin/ccxvue` + `xvue/xvuelc.c` + `libX11` linker path stays functional during the transition. A build-time switch (env var or a `bin/cbl_tout_qt` variant) picks the backend. Both backends can be built from the same tree.
- [ ] **Visual A/B validation on a subset of `testa/`** — selected canonical cases per module (mesher, elasticity, fluid, thermal, nonlinear) are run through both backends and compared side-by-side by eye. The chosen subset, per-case checklist, and a short validation log are captured as part of the final phase.
- [ ] **Drop `libX11`/ImageMagick runtime dependencies after X11 retirement** — once the one-release-cycle A/B window closes, remove the X11 linker lines, delete `xvue/xvuelc.c`, remove `bin/convertepsgif` / ImageMagick invocations, and update `README` / `LISEZMOI` / install scripts.

### Out of Scope

<!-- Explicit boundaries with reasoning. Prevents re-adding later. -->

- **Windows and macOS ports** — Linux-x86_64 is the only target of this migration. Cross-platform is a separate future initiative and would multiply the test matrix. Kept out to keep this project finishable.
- **Changes to numerical solver code** — `mail/`, `elas/`, `flui/`, `ther/`, `reso/`, `nlse/`, and `util/` stay bit-identical. The invariant is: Fortran code must not notice the migration.
- **Rewriting the shell build system for non-xvue modules** — `bin/cbl_tout` and the per-module `bin/cb*` scripts continue to drive the Fortran build. CMake is introduced only for `xvue/`. A full CMake migration is deferred.
- **Automated CI (GitHub Actions / GitLab CI / Jenkins)** — nice to have, explicitly not part of this initiative. Quality is maintained by running the `testa/` A/B subset manually.
- **Fortran modernization** — `IMPLICIT NONE` rollout, F77→F90 module conversion, removing `GOTO`s. Out of scope; orthogonal to the graphics migration.
- **Qt-based re-skinning of the terminal launchers** — `bin/INITIER`, `bin/MAILLER`, etc., stay as bash scripts. Only the graphics windows they spawn become Qt.
- **New physics features / new solvers** — this is a UI migration, not a feature release.
- **Non-`xvue/` bug fixes** — unrelated bugs in `mail/` / `flui/` / etc. are tracked separately and not bundled into xvue-qt phases (except when they block A/B validation).

## Context

### Technical environment

- **Platform**: Debian trixie/sid, Linux x86_64
- **Compilers available**: `gfortran` (system), `gcc` (system), `g++` (system)
- **Qt available from apt**: `qt6-base-dev` 6.10.2, `qtbase5-dev` 5.15.17 (Qt 6 chosen)
- **CMake**: available from apt; will be used only inside `xvue/`
- **Current graphics**: X11 via `libX11-dev`, linked from `/usr/X11R6/lib64` (a compat symlink on modern Debian — known fragility)
- **Current image export**: EPS written by `xvuelc.c`, animated GIF via ImageMagick `convert` shell-out

### Prior work / codebase state

Full codebase map lives in `.planning/codebase/`:

- `STACK.md`, `INTEGRATIONS.md` — languages, compilers, external deps
- `ARCHITECTURE.md`, `STRUCTURE.md` — module layout, entry points, data flow
- `CONVENTIONS.md` — Fortran 77 fixed-form norms, naming, include-file discipline
- `TESTING.md` — manual `testa/`/`testf/` validation workflow, no CI
- `CONCERNS.md` — legacy F77 debt, tight `COMMON`-block coupling, planned Qt migration explicitly listed as future work, X11 obsolescence risk

`xvue/xvuelc.c` is **3619 lines** with approximately **60 `extern "C"` entry points**: drawing primitives (points, lines, faces, rectangles, ellipses, text), color and font management (`xvchargefonte`, `xvCouleursImposees`, `xvStockeRGBtoColormap`), pixmap double-buffering (`fenetremempx`, `mempxfenetre`, `sauvefenetre`, `restaurefenetre`), mouse event capture (`xvsouris`, `xvsouris2`, `deplsouris`), PostScript export (`xvpostscript`), and the interactive input loop that backs the text-lexicon menu. This is the entire API surface to reimplement.

The text-lexicon menu is not a single file — it is the pattern "read a command string from a graphics-window text area, dispatch to a Fortran subroutine" scattered across `mail/`, `elas/`, `flui/`, `ther/`, `nlse/` interactive drivers. Level 3 modernization means auditing each module and producing a QMenuBar/QToolBar/QDialog surface that dispatches to the same Fortran subroutines.

### User feedback / motivation

- `CLAUDE.md` already lists a future Qt migration as an active project goal.
- X11 is in long-term maintenance; Wayland compositors increasingly treat it as legacy. Replacing it removes a portability risk.
- Dropping ImageMagick as a hard runtime dependency simplifies installs.
- Real menus/toolbars are expected to be a noticeable UX improvement for end users over the current text-command lexicon.

## Constraints

- **Tech stack**: Qt 6 (6.10 from Debian apt) — Qt 5 explicitly rejected as maintenance-mode. Why: multi-year horizon, future-proofing.
- **Tech stack**: C++ (for the Qt side) and Fortran 77/95 (unchanged). Bridge is `extern "C"` with the same function names as today's `xvuelc.c`. Why: zero churn on the thousands of Fortran call sites.
- **Tech stack**: CMake owns only `xvue/`, producing a static library. Why: Qt projects need `moc`, which shell scripts handle poorly; a full CMake migration of the Fortran build is out of scope for this project.
- **Platform**: Linux-x86_64 only. Why: finite scope, and the Fortran build's `-mcmodel=large` + hard-coded lib paths constrain portability anyway.
- **Compatibility**: X11 build must remain functional throughout this project and for one release cycle after the first Qt release. Why: user-facing A/B validation and safe rollback.
- **Compatibility**: `extern "C"` entry-point names and signatures in the Qt implementation must be byte-identical to today's `xvuelc.c`. Why: Fortran wrappers in `xvue/*.f` stay untouched.
- **Performance**: Qt render path must handle the same mesh sizes that the X11 build handles today on the same hardware. No measured latency regression on the `testa/` subset.
- **Dependencies**: gfortran, gcc, g++, Qt 6 from Debian apt, CMake. No vendored builds, no source compilation of Qt, no non-apt packages required.
- **Team**: single developer, no CI, no code review partner. Implies: every phase must leave the tree in a working buildable state, and commits must be small enough to roll back.
- **Testing**: no automated test harness; validation is manual A/B on `testa/` subset, eyeballed.

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| **Qt 6 (not Qt 5)** | Qt 6.10 is in Debian trixie and is the current generation. Qt 5 is maintenance-only. Multi-year horizon favors Qt 6 despite slightly smaller tutorial corpus. | — Pending |
| **Replace `xvue/xvuelc.c` entirely** (not coexist forever) | End state is one graphics backend. Coexistence is a temporary bridge, not a permanent architecture. | — Pending |
| **Parallel X11 build during transition** | Safety net: A/B comparison, rollback, user-visible continuity. Both backends live in the tree for one release cycle. | — Pending |
| **Preserve `extern "C"` names and signatures** | Decouples the Fortran side from the migration. `xvue/*.f` wrappers untouched; thousands of call sites in `mail/`/`elas/`/`flui/`/`ther/` untouched. | — Pending |
| **CMake owns only `xvue/`** | Qt tooling (`moc`, `uic`, `rcc`) is hostile to shell-script builds. Rest of the Fortran build stays as-is to keep blast radius contained. | — Pending |
| **Full UX modernization (Level 3)** | This is a once-per-decade opportunity. Real `QMenuBar`/`QToolBar`/`QDialog` replace the text lexicon across all interactive modules. | — Pending |
| **Qt takes over image/GIF export** | Drops ImageMagick runtime dep, unifies the code path, `QImage`/`QImageWriter`/`QMovie` cover everything `xvuelc.c` + `convertepsgif` currently produce. | — Pending |
| **Linux x86_64 only** | Fortran build already constrained to this target by `-mcmodel=large` and `/usr/X11R6/lib64`. Cross-platform is a separate future initiative. | — Pending |
| **Big-bang Qt landing + one-release A/B window** | Incremental porting inside a single backend would leave the graphics layer inconsistent. Clean replacement, validated in parallel with the old backend. | — Pending |
| **Validation = manual A/B on `testa/` subset** | No CI, no reference outputs. Eyeballing selected canonical cases (mesher, elas, flui, ther, nlse) is the realistic bar. | — Pending |
| **Debian trixie + apt Qt 6** | Zero-friction install for the single developer. No vendored Qt. | — Pending |

## Evolution

This document evolves at phase transitions and milestone boundaries.

**After each phase transition** (via `/gsd-transition`):
1. Requirements invalidated? → Move to Out of Scope with reason
2. Requirements validated? → Move to Validated with phase reference
3. New requirements emerged? → Add to Active
4. Decisions to log? → Add to Key Decisions
5. "What This Is" still accurate? → Update if drifted

**After each milestone** (via `/gsd-complete-milestone`):
1. Full review of all sections
2. Core Value check — still the right priority?
3. Audit Out of Scope — reasons still valid?
4. Update Context with current state

---
*Last updated: 2026-04-18 after Phase 5 (event bridge & blocking reads) completion — `XvueEventBridge` QObject with `BlockingDepthGuard` RAII depth counter, motion coalescing via `QTimer::singleShot` deferred-quit (X11 `XEventsQueued(QueuedAfterFlush)` parity), and real Fortran-ABI bodies for `xvsouris_` / `xvsouris2_` / `xvpause_` / `deplsouris_` / `initaccrochage_`. Strategy B accrochage (save tile → blit sprite → restore on motion) reuses Phase 4's `saved_canvas_` pattern rather than emulating X11 XOR raster ops. `MEFISTO_XVSOURIS_AUTOEXIT` extended to `xvpause_` in both backends for headless-capture parity. 31 QTest cases pass (2 documented skips), ABI frozen at 57 symbols, SHELL-03 clean, human A/B drag verdict approved on pan2d + torus. Phase A kwin-mcp smoke captured 449 real event-filter round-trips with `depth=1` throughout — Assumption A2 (Qt `AA_CompressHighFrequencyEvents` vs bridge deferred-quit ordering) resolved empirically. D-09 Wayland caveat (`QCursor::setPos` no-op on pure Wayland, affects `mail/trfasevo.f:202`) documented in `xvue/qt/README.md`. Requirements EVENT-01..EVENT-08 closed. Phase 6 (Level-3 UX chrome — `QMenuBar`/`QToolBar`/`QDialog` + per-module lexicon audit) is next.*
