# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is MEFISTO

MEFISTO (MEsh and Finite element Software TO...) is a scientific computing package by Alain Perronnet (LJLL, UPMC Paris). It provides interactive 2D/3D mesh generation and finite element solvers for elasticity, fluid mechanics, thermal, nonlinear, and wave problems, with Qt 6 graphical visualization (the legacy X11 backend was retired in Phase 9 — git tag `v1.0-pre-retire` for rollback).

## Environment variables (required)

```bash
export MEFISTO=/path/to/mefistosource   # this repo root
export MEFISTOX=$HOME/mefistox          # user project working directory
export PATH=.:$PATH:$MEFISTO/bin
export CDPATH=.:$HOME:$MEFISTO:$MEFISTOX
```

## Build system

The build uses `gfortran` + `gcc` + `g++` with Qt 6 (CMake drives `xvue/qt/`; the rest of the build is shell scripts). All compilation scripts are shell scripts located in `bin/` (or `bin.lnx64/` for the 64-bit Linux prebuilt variants).

| Command | Effect |
|---|---|
| `bin/cbl_tout` | Compile **all** modules and link all executables |
| `bin/cbinit` | Compile the INITIER module only |
| `bin/cbmail` | Compile the mesher (MAILLER) module only |
| `bin/cbelas` | Compile the elasticity (ELASTICER) module only |
| `bin/cbflui` | Compile the fluid (FLUIDER) module only |
| `bin/cbther` | Compile the thermal (THERMICER) module only |
| `bin/cbnlse` | Compile the nonlinear solver (NLSER) module only |

Compiled executables land in `pp/` (e.g. `pp/ppelas`, `pp/ppflui`, `pp/ppinit`, …). The shell scripts in `bin/` (INITIER, MAILLER, ELASTICER, …) wrap these executables and set up the project context.

### Dependencies

- `gfortran` (Fortran 77/95 + OpenMP)
- `gcc`
- `g++` (C++ compiler for the Qt 6 backend in `xvue/qt/`)
- `qt6-base-dev` (Qt 6 graphics backend headers + libraries)
- `qt6-image-formats-plugins` (TIFF/WebP/etc. image-format plugins for `XvueExport`)
- `ffmpeg` (animated GIF export via `XvueExport::saveGifTo`)
- `cmake` (>= 3.16, drives `xvue/qt/` build)

## Repository structure

```
incl/       Fortran include files (.inc) — shared data structures across all modules
mail/       Mesher Fortran sources (.f)
elas/       Elasticity solver sources (.f)
flui/       Fluid solver sources (.f)
ther/       Thermal solver sources (.f)
reso/       Linear solver sources (.f)
util/       Shared utility Fortran routines
xvue/       Qt 6 graphical display layer (Fortran wrappers + C++ in xvue/qt/)
prpr/       Main program entry points (ppinit.f, ppmail.f, ppelas.f, …)
pp/         Compiled executables (output of build)
bin/        Shell scripts: launchers (INITIER, MAILLER, …) + compilation scripts (cb*)
bin.lnx64/  Prebuilt Linux 64-bit binaries and scripts
doc/        Active documentation (symlinked to doca/ or docf/ at install time)
doca/       English documentation
docf/       French documentation
td/         Tutorials and demo data (da/df=demos, ia/if=init, ma/mf=menus — a=English, f=French)
testa/      English test cases
testf/      French test cases
```

The `incl/homdir.inc` file is **generated at build time** by `cbl_tout` — it encodes `$MEFISTO` as a Fortran `PARAMETER` string. Do not edit it manually.

## Language and module conventions

All solver and utility code is **Fortran 77** (fixed-form, column 7+) with some Fortran 95 extensions for OpenMP parallel variants (`*_OMP` executables). The Qt 6 backend lives in `xvue/qt/` (C++); the legacy `xvue/xvuelc.c` C source was retired in Phase 9. Include files in `incl/` declare the shared common blocks and parameters that connect all modules.

## Programming norms

The project's coding norms are documented in `doc/normes.ps` (PostScript format — view with `evince` or `gs`). These norms **must be respected** in all modifications: naming conventions, comment style, fixed-form Fortran column layout, etc.

## Active project goals

- **Bug fixes**: identify and correct existing bugs without altering the overall behaviour.
- **Qt migration (completed in Phase 9)**: the X11/Motif graphical layer was replaced by Qt 6 (`xvue/qt/`). The legacy `xvue/xvuelc.c` is retired; `git tag v1.0-pre-retire` preserves the X11 path for rollback. The Fortran solver wrappers in `xvue/*.f` are unchanged.

## Working rules

### Compilation must never break
After every change, verify that the affected module still compiles with its `cb*` script. Before committing, the full build (`bin/cbl_tout`) must succeed.

### Tests
- Small tests in `testa/` or `testf/` **must continue to pass** after every change.
- For large or long-running tests, ask the user to run them before declaring a change complete.
- Always prefer the smallest relevant test case when checking a fix.

### Asking before acting
- If a C include or an Ubuntu system package is missing to compile, **ask the user** to install it — do not try to work around missing headers.
- Explain clearly what you are doing at each step before doing it.

### Git discipline
- Commit after every logical step where rolling back would be useful (working build, passing test, completed sub-fix…).
- Commit messages must describe *what* changed and *why*.
- Never force-push; never bypass hooks.

## Running a project

User projects live under `$MEFISTOX/ProjectName/`. The typical workflow is:

1. `INITIER` — initialise a new project
2. `MAILLER` — interactive mesh generation (opens a Qt 6 window)
3. `ELASTICER` / `FLUIDER` / `THERMICER` / `NLSER` / `WAVER` — run the chosen solver
4. To save and exit any interactive module: return to the main menu and type `99;` (**never** Ctrl-C)

### Quick test recipe (for maintainers verifying a build)

After `bin/cbl_tout` finishes, the smoke recipe is:

```bash
# 1. Make sure env is exported (each terminal session)
export MEFISTO=/home/mefisto/git/mefisto
export MEFISTOX=$HOME/mefistox
export PATH=.:$PATH:$MEFISTO/bin

# 2. Set up a project once (creates $MEFISTOX/<name>/ and seeds defaults)
INITIER
# Type a project name when prompted, e.g.: nafems_le1

# 3. Launch the mesher interactively
MAILLER
# Type the same project name; the Qt window opens with menus + canvas.
# In the canvas: click menu items, draw points/lines/surfaces, etc.
# To save and exit: return to the main menu and type 99;

# 4. Then run a solver (same project name flow)
ELASTICER     # elasticity
FLUIDER       # fluid
THERMICER     # thermal
NLSER         # nonlinear
```

### Direct binary launch (advanced — Phase 9 auto-batch)

The `pp/pp<module>` binaries can be invoked directly **without** the launcher
scripts. The Qt port detects when stdin is not a TTY (KDE launcher, Wayland
session, pipeline, kwin-mcp) and auto-sets `MEFISTO_BATCH_X11=1` +
`MEFISTO_XVSOURIS_AUTOEXIT=1` so Fortran takes the batch path:

```bash
cd $MEFISTOX/<project>
$MEFISTO/pp/ppmail <data-file>      # e.g. ppmail nafems_le1.mail
```

When stdin **is** a TTY (normal terminal launch), the binary stays
interactive and reads prompts from the user — same as the launcher scripts.

To force batch mode in an interactive terminal: `MEFISTO_BATCH_X11=1
MEFISTO_XVSOURIS_AUTOEXIT=1 $MEFISTO/pp/ppmail <data-file>`.

### Test restart procedure

If a launched module hangs or you need to restart cleanly:

```bash
# 1. Find any stale process
pgrep -a 'pp(mail|elas|flui|ther|nlse)' | head

# 2. Kill it
kill <PID>          # or pkill -f 'pp(mail|elas|flui|ther|nlse)'

# 3. Rebuild after source edits
cd $MEFISTO
bin/cbl_tout        # full rebuild; exit 0 expected

# 4. Re-launch via the launcher script (recommended)
MAILLER             # then type project name
```

### Verification gates (post-Phase-9)

```bash
cd $MEFISTO
bash bin/test_no_x11_in_build.sh         # no libX11 in active build
bash bin/test_no_imagemagick_in_qt.sh    # no convert under xvue/qt/
bash bin/test_no_lvideo.sh               # no LVIDEO Fortran callbacks
bash xvue/qt/cmake/verify_abi.sh \
     xvue/qt/build/libxvueqt.a \
     xvue/qt/include/xvue_qt_api.h        # nm count: 58 header count: 58
```

All four must exit 0 on a healthy tree.
