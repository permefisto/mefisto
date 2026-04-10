# Structure

## Top-level layout

```
mefisto/
├── CLAUDE.md              project instructions for Claude Code
├── README, LISEZMOI       English / French readme
├── install.bash           end-user install script
├── instalsource.bash      source install script
│
├── incl/                  ~174 Fortran include files (.inc, .inc95)
├── prpr/                  main PROGRAM entry points (one per executable)
├── mail/                  mesher sources                (~1097 files)
├── elas/                  elasticity solver              (~76 files)
├── flui/                  fluid mechanics solver        (~298 files)
├── ther/                  thermal solver                (~241 files)
├── reso/                  linear-system solvers shared by modules (~190 files)
├── util/                  shared Fortran utility toolbox (~856 files)
├── xvue/                  X11 graphics layer (Fortran + 1 C file) (~224 files)
│
├── prpr/                  main programs (ppinit.f, ppmail.f, ppelas.f, ...)
├── pp/                    compiled executables (build output)
├── bin/                   launchers + compilation scripts (shell)
├── bin.lnx64/             prebuilt / Linux-64 variant of bin/
│
├── doc/                   active doc (symlink into doca/ or docf/)
├── doca/                  English documentation (including normes.ps)
├── docf/                  French documentation
├── td/                    tutorials & demo data
│   ├── da/, df/           demos (a = English, f = French)
│   ├── ia/, if/           init data
│   └── ma/, mf/           menu data
├── testa/                 English test cases
└── testf/                 French test cases
```

## Directory purposes

### Source — shared declarations

- **`incl/`** — Fortran include files declaring all `COMMON` blocks, `PARAMETER`s, and struct-like shared state used across the solvers. Naming:
  - `a___<topic>.inc` — single-topic shared data (`a___arete.inc` edges, `a___face.inc` faces, `a___materiaux.inc` materials, `a___lestemps.inc` time parameters, `a___contrainte.inc` stress, `a___morse.inc` sparse storage…)
  - `a_<object>__<aspect>.inc` — object with sub-aspect (`a_ligne__bspline.inc`, `a_objet__definition.inc`)
  - `*.inc95` — Fortran 95 variants for OpenMP (`_OMP`) executables (e.g. `a___npef.inc95`)
  - `a_fonction__arbre.inc`, `a~_fonction__arbre.inc` — trailing-tilde file is a hand-kept backup/variant
  - **`incl/homdir.inc`** — GENERATED at build time by `bin/cbl_tout`; encodes `$MEFISTO` as a Fortran `PARAMETER` string. Never edit by hand.

### Source — executables and solvers

- **`prpr/`** — short Fortran `PROGRAM` units, one per binary:
  - Core modules: `ppinit.f`, `ppmail.f`, `ppelas.f`, `ppflui.f`, `ppther.f`, `ppnlse.f`
  - Auxiliary: `ppadap.f`, `ppdvsr.f`, `pppara.f`, `pppbo.f`, `pppoba.f`, `ppquat.f`
  - Specialized: `ppbrezfort2d.f`, `ppbrezfort3d.f`, `ppbrezfort3dmoy.f`, `ppst1c6c.f`, `ppst3c6c.f`, `ppst5c6c.f`, `prracpol3.f`, `picarre.f`, `q1.f`
  - X11 test drivers: `xvtest1.f`..`xvtest4.f`
- **`mail/`** — mesh generation and geometry (largest module, ~1100 files). Contains everything from edge/face routines to Delaunay triangulation and 3D tetrahedralization.
- **`elas/`** — linear elasticity assembly and solve.
- **`flui/`** — fluid mechanics (Navier–Stokes variants, Boussinesq coupling).
- **`ther/`** — thermal / heat-transfer finite elements.
- **`reso/`** — linear solvers (Morse sparse storage, iterative and direct methods) shared across all physics modules.
- **`util/`** — pervasive utility toolbox: I/O, memory management, string/number conversion, sorting, error printing, timers. Every executable links it in.
- **`xvue/`** — X11 graphical display layer. All Fortran files are wrappers around Xlib calls exposed from **`xvue/xvuelc.c`** (the only C source file in the project). This directory is the target of the planned Qt migration; keep Xlib calls confined here.

### Build, launch, and binaries

- **`bin/`** — shell scripts. Two groups:
  1. **Launcher scripts** (uppercase, e.g. `bin/INITIER`, `bin/MAILLER`, `bin/ELASTICER`, `bin/ELASTICER_OMP`, `bin/FLUIDER`, `bin/THERMICER`, `bin/HEATER`, `bin/NLSER`, `bin/DELASTICER`, `bin/DTHERMICER`) — end-user entry points. They detect the language via `$MEFISTO/td/m/anglais`, check memory/disk, `cd` into the project directory, and exec the matching `pp/pp*` executable.
  2. **Compilation scripts** (lowercase `cb*`): `bin/cbl_tout` (build everything), `bin/cbinit`, `bin/cbmail`, `bin/cbelas`, `bin/cbflui`, `bin/cbther`, `bin/cbnlse`, `bin/cbadap`, `bin/cbonde`, `bin/cbpara`, `bin/ccxvue` (C part of X11 layer), plus OpenMP/variant builders (`bin/cblg_tout`, `bin/cblg_ompf`, `bin/cbl_tous_f`, `bin/cblg_all`, …).
  3. Small color/window helpers (`bin/fenrouge`, `bin/fenvert`, …) and conversion helpers (`bin/convertepsgif`, `bin/chnomfich`, `bin/chgsufx`).
- **`bin.lnx64/`** — Linux-64 prebuilt variant of `bin/` with the same script layout, shipped for users who cannot recompile.
- **`pp/`** — output directory for compiled executables produced by the `bin/cb*` scripts. Git-ignored after the cleanup in commit `ac282f8`; only `pp/pxyz` (or similar small helper binaries) is kept in-tree.

### Documentation

- **`doc/`** — active documentation, symlinked at install time to either `doca/` (English) or `docf/` (French).
- **`doca/`, `docf/`** — parallel English and French documentation trees. The **coding norms document `doc/normes.ps`** lives here and is authoritative for style questions.

### Tests and tutorials

- **`td/`** — tutorials and demo data used by interactive sessions:
  - `td/da/`, `td/df/` — demo geometries / meshes (`a` = English, `f` = French)
  - `td/ia/`, `td/if/` — initial data for tutorials
  - `td/ma/`, `td/mf/` — menu scripts for guided walkthroughs
  - `td/m/anglais` — presence of this file switches launchers to English
- **`testa/`, `testf/`** — English and French test cases. Each test is a small project directory plus expected results, driven manually through launcher scripts. These must keep passing after every change; long-running tests are run on request.

## Where to add new code

| Kind of change | Add to |
|---|---|
| New shared data structure / common block | `incl/a___<name>.inc` (+ `a___<name>.inc95` if used by `_OMP`) |
| New subroutine in an existing solver | the module's directory (`mail/`, `elas/`, `flui/`, `ther/`, `reso/`) |
| New cross-module helper | `util/` |
| New main program / executable | new file in `prpr/`, new build script `bin/cb<name>`, new launcher `bin/<NAME>` |
| X11 display tweak | `xvue/` (Fortran wrapper) and, if needed, `xvue/xvuelc.c` (Xlib call) |
| New test case | `testa/<name>/` (English) and/or `testf/<nom>/` (French) |
| New tutorial / demo | `td/d{a,f}/<name>/` with matching `td/m{a,f}/` menu file |

## Naming and file conventions

- **Source files**: lowercase, `.f` for Fortran 77 fixed-form, `.inc` for include files, `.inc95` for Fortran 95 variants, `.c` and `.h` for the X11 C layer.
- **Fixed-form Fortran 77**: comments in column 1 (`C`), labels in columns 1–5, continuation character in column 6, code in columns 7–72. Enforced by compiler flags in `bin/cb*`.
- **Launcher scripts**: uppercase (`INITIER`, `MAILLER`, …) for interactive end-user entry points; lowercase `cb*` for build; lowercase `d*`/`D*` for debug or development variants.
- **Executables in `pp/`**: `pp<module>` (e.g. `pp/ppelas`, `pp/ppflui`), with optional `_OMP` suffix for OpenMP builds.
- **Include-file prefixes** in `incl/`: the leading `a___` / `a_` / `a~_` is a sort-key used to group related common blocks alphabetically; do not rename without updating every `INCLUDE` statement that references them.
- **Bilingual identifiers**: many variable and subroutine names mix French and English (`TPSINI`/`TPSFIN` for initial/final time, `aretefr` for "free edge"). Keep the existing language when editing a routine rather than translating wholesale.
