# Architecture

## Overview

MEFISTO (MEsh and FInite element Software TO...) is a **modular, monolithic scientific computing package** written primarily in Fortran 77 (with Fortran 95 + OpenMP extensions in a few places) and a thin C layer for X11 graphics. It provides interactive 2D/3D mesh generation and a family of finite-element solvers (elasticity, fluid mechanics, thermal, nonlinear, wave) with on-screen visualization.

The architecture is classic scientific-Fortran: **independent executables** — one per functional module — that **share data via Fortran `COMMON` blocks declared in `.inc` include files**, plus filesystem-based persistence of per-project state under `$MEFISTOX/<ProjectName>/`. The user drives the system through **shell launcher scripts** which set up the environment and invoke the right compiled executable from `pp/`.

## High-level pattern

```
   user
    │
    ▼
┌────────────────────────┐
│ bin/<LAUNCHER>         │   shell script (sets env, cd into project dir, runs pp/pp*)
└──────────┬─────────────┘
           │
           ▼
┌────────────────────────┐
│ pp/pp<module>          │   compiled Fortran executable (one per module)
└──────────┬─────────────┘
           │                 reads/writes common blocks via incl/*.inc
           ▼
┌────────────────────────┐
│ Solver module (.f)     │   mail/, elas/, flui/, ther/, reso/, nlse/, ...
│  + util/ (shared)      │
│  + xvue/ (X11 display) │
└──────────┬─────────────┘
           │
           ▼
    $MEFISTOX/<Project>/   project files (mesh, solution, plots)
```

There is no plugin system, no dynamic dispatch, and no runtime configuration framework — the coupling between modules is done **at compile time** through shared include files and at **runtime through on-disk project files**.

## Layers

1. **Launcher layer** (`bin/`, `bin.lnx64/`)
   Shell scripts that set environment variables, resolve the project directory, and invoke the compiled executable. Examples: `bin/INITIER`, `bin/MAILLER`, `bin/ELASTICER`, `bin/FLUIDER`, `bin/THERMICER`, `bin/NLSER`, `bin/HEATER`. They also handle bilingual (French/English) prompts by probing `$MEFISTO/td/m/anglais`.

2. **Entry-point layer** (`prpr/`)
   One short Fortran `PROGRAM` per module that initializes common blocks and dispatches to the module's driver. Examples:
   - `prpr/ppinit.f` — INITIER (project initialization)
   - `prpr/ppmail.f` — MAILLER (mesh generator)
   - `prpr/ppelas.f` — ELASTICER (elasticity solver)
   - `prpr/ppflui.f` — FLUIDER (fluid solver)
   - `prpr/ppther.f` — THERMICER (thermal solver)
   - `prpr/ppnlse.f` — NLSER (non-linear solver)
   - `prpr/pppoba.f`, `prpr/pppbo.f`, `prpr/pppara.f` — auxiliary entry points
   - `prpr/ppadap.f`, `prpr/ppdvsr.f`, `prpr/ppquat.f` — adaptive/support programs
   - `prpr/xvtest[1-4].f` — standalone X11 display test programs

3. **Solver-module layer** (one directory per scientific domain)
   Each module holds the domain-specific Fortran 77 subroutines, all fixed-form:
   - `mail/` — mesher (mesh generation, refinement, geometry — ~1100 files)
   - `elas/` — linear elasticity solver (~76 files)
   - `flui/` — fluid mechanics (~298 files)
   - `ther/` — thermal/heat transfer (~241 files)
   - `reso/` — linear system solvers used by all modules (~190 files)

4. **Utility layer** (`util/`)
   Large shared toolbox (~856 files) of Fortran subroutines used by every solver: I/O helpers, memory management, string/number conversions, sort/search, error reporting, time measurement, etc. No single module "owns" `util/` — it is effectively a global library linked into every executable.

5. **Graphics layer** (`xvue/`)
   X11-based display module. Mostly Fortran wrappers (~224 files) plus one large C file (`xvue/xvuelc.c`) that calls the Xlib API directly. This is the **only** C source in the project and is the layer slated for future Qt replacement. Any refactoring should keep X11-specific calls confined here.

6. **Shared-data layer** (`incl/`)
   ~174 Fortran include files (`a___*.inc`, `a_*__*.inc`, etc.) declaring `COMMON` blocks, `PARAMETER` constants, and derived-type-like structures. Every source file in the solver layers `INCLUDE`s whatever common blocks it needs. This is the **spine** of the architecture: changing a common block touches every module that includes it. `incl/homdir.inc` is **generated at build time** by `bin/cbl_tout` and encodes the `$MEFISTO` install path as a Fortran `PARAMETER` string — it must not be edited manually.

## Data flow

A typical user session:

1. **INITIER** creates a project directory under `$MEFISTOX/<ProjectName>/` and writes initial scratch files.
2. **MAILLER** opens an interactive X11 window, lets the user define geometry (points, lines, surfaces, volumes), and generates a mesh. The resulting mesh is written to the project directory.
3. **ELASTICER / FLUIDER / THERMICER / NLSER / WAVER** reads the mesh from the project directory, runs the finite-element solver, writes back the solution, and displays results via the X11 layer. Several of these have `_OMP` variants (built from `cblg_ompf` or `bin.lnx64` equivalents) for OpenMP-parallel execution.
4. All interactive modules exit cleanly via the menu command `99;` (typed by the user) — **never Ctrl-C**, which would leave project files in an inconsistent state.

Animated GIF output from post-processing is produced by shelling out to ImageMagick's `convert`.

## Key abstractions

- **Fortran `COMMON` blocks in `incl/`** — the universal data-sharing mechanism. Conceptually they play the role that structs/classes play in modern languages: each `.inc` file defines one shared structure (e.g. `a___arete.inc` for edges, `a___face.inc` for faces, `a___contrainte.inc` for stress tensors, `a___lestemps.inc` for time parameters, `a___materiaux.inc` for materials). Files named `*.inc95` are the Fortran 95 variants used by `_OMP` executables.
- **Lexicon/menu system** — most interactive modules share a text-command lexicon (ADAM-style) implemented in `util/`, which parses user commands like `99;` and dispatches to the matching Fortran subroutine.
- **Secondary memory / segment management** — large mesh and matrix data are stored in scratch files in the project directory and paged in/out through helper routines in `util/`. This lets the solvers handle problems larger than RAM, at the cost of tight coupling between solver code and the on-disk format.
- **Fixed-form Fortran 77 column layout** — all `.f` files in the solver layers use strict column-7-onwards source, column-1 comments (`C`), and 6-character labels. This is a hard convention enforced by the compiler flags in `bin/cb*` scripts.
- **Bilingual identifiers and prompts** — variable names and menu strings are a mix of French (the original LJLL/UPMC Paris language) and English. Launchers detect `$MEFISTO/td/m/anglais` to pick the interface language.

## Entry points

The **user-facing** entry points are the launcher scripts in `bin/`:

| Launcher | Purpose | Underlying executable | Main program |
|---|---|---|---|
| `bin/INITIER` | Create / initialize a project | `pp/ppinit` | `prpr/ppinit.f` |
| `bin/MAILLER` | Interactive 2D/3D mesher | `pp/ppmail` | `prpr/ppmail.f` |
| `bin/ELASTICER`, `bin/ELASTICER_OMP` | Linear elasticity solver | `pp/ppelas` (+ `_OMP`) | `prpr/ppelas.f` |
| `bin/FLUIDER` | Fluid mechanics solver | `pp/ppflui` | `prpr/ppflui.f` |
| `bin/THERMICER`, `bin/HEATER` | Thermal solver | `pp/ppther` | `prpr/ppther.f` |
| `bin/NLSER` | Nonlinear solver | `pp/ppnlse` | `prpr/ppnlse.f` |
| `bin/DELASTICER`, `bin/DTHERMICER`, `bin/dinitier`, `bin/dmailler` | Debug / development variants | corresponding `pp/pp*` | same as above |

The **developer-facing** entry points are the compilation scripts, also in `bin/`:

| Script | Effect |
|---|---|
| `bin/cbl_tout` | Rebuild every module and link every executable (the "build everything" command) |
| `bin/cbinit`, `bin/cbmail`, `bin/cbelas`, `bin/cbflui`, `bin/cbther`, `bin/cbnlse` | Per-module rebuilds |
| `bin/ccxvue` | Compile the C part of the X11 display layer (`xvue/xvuelc.c`) |
| `bin/cblg_tout`, `bin/cblg_ompf` | Variants that build the OpenMP executables |

Standalone test programs (not wrapped by launchers) live next to the main entry points: `prpr/xvtest1.f`..`prpr/xvtest4.f` exercise the X11 display primitives directly.

## Cross-cutting concerns

- **I/O and project state** — handled uniformly through `util/` helpers; every module reads and writes the same project directory, which is why launchers take care to `cd` there first.
- **Error handling** — traditional Fortran: status integers, printed diagnostics, and `STOP`. There is no structured exception mechanism.
- **Parallelism** — only through OpenMP in the `_OMP` executables, compiled with `gfortran -fopenmp` from the Fortran 95 variants of the include files and solver sources.
- **Graphics isolation** — every call into Xlib goes through `xvue/`, so the planned Qt migration can in principle replace `xvue/` without touching solver logic. In practice, the Fortran wrappers around X primitives are referenced from many solver files, so the migration surface is the `xvue/` Fortran API, not Xlib itself.
- **No build system beyond shell** — there is no Makefile, CMakeLists, or autotools. Dependency tracking is implicit in the order of file lists inside the `bin/cb*` scripts. Rebuilds are coarse-grained.
- **No automated test harness** — correctness is checked by running the demo cases in `testa/` (English) and `testf/` (French) by hand and comparing output.
