# Conventions

## Authoritative reference

The coding norms for MEFISTO are documented in **`doc/normes.ps`** (PostScript — view with `evince` or `gs`). That document is authoritative; everything below summarizes what is observable from the source tree and must stay consistent with it. Per `CLAUDE.md`, the norms **must be respected** in every modification.

## Language and source form

- **Fortran 77 fixed-form** for the vast majority of `.f` files in `mail/`, `elas/`, `flui/`, `ther/`, `reso/`, `util/`, `xvue/`, and `prpr/`. Columns 1–5 for labels, column 6 for the continuation character, columns 7–72 for code, column 1 for comments (`C` in the first column).
- **Fortran 95 variants** for OpenMP-parallel builds. These live alongside the F77 sources; include files used by parallel executables have the suffix **`.inc95`** (e.g. `incl/a___npef.inc95`).
- **One subroutine per file** is the strong convention in `mail/`, `util/`, etc. The filename matches the subroutine name (e.g. `mail/a1lgarti.f` contains `SUBROUTINE A1LGARTI(...)`).
- **Single C source**: only `xvue/xvuelc.c` is C. All other graphics code is Fortran wrappers around the routines exposed by that file.
- **Line continuation**: the character in column 6 is typically `%` (e.g. `% L1ARET, L2ARET, LARETE,`) rather than the more common `&`. Keep the existing style when editing.

## Comment style

- Comments live in column 1 with `C` (sometimes `C+++...++` banners for section separators). Inline `!` comments are not used in the solver layers.
- Every subroutine begins with a structured French (or bilingual) header describing:
  - `C BUT :` / `C -----` — purpose
  - `C ENTREES:` — inputs
  - `C MODIFIES :` — in-out arguments
  - `C SORTIES:` — outputs
  - Each argument documented one per line (e.g. `C MXSOM  : NOMBRE MAXIMAL DE SOMMETS DE LA TRIANGULATION`).
  - Common trailing fields: author, date, version.
- Banners use lines of `+` characters (`C+++++++++++...`) to delimit the header from the body.
- This header format is effectively a contract — do not remove it, and when adding a new argument, document it in the header at the same time as the signature change.

Example from `mail/a1lgarti.f`:

```
      SUBROUTINE A1LGARTI( MXSOM,  NBSOM,  XYZSOM,
     %                     L1ARET, L2ARET, LARETE,
     %                     MXTRIA, NBTRIA, NSTRIA,  IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AVEC FONCTION TAILLE_IDEALE(x,y,z) ou DARETE VALEUR PAR DEFAUT
C -----    de include"./incl/darete.inc"
...
C ENTREES:
C --------
C MXSOM  : NOMBRE MAXIMAL DE SOMMETS DE LA TRIANGULATION
...
```

## Naming

- **Subroutines and filenames**: short mnemonic names, usually ≤8 characters, uppercase in code and lowercase on disk (`A1LGARTI` ↔ `a1lgarti.f`). Many names encode a domain prefix: `ad*` for adaptation, `af*` for face, `ai*` for arête (edge), `a1*`/`a2*` for 1D/2D variants, `pp*` for main programs, `cb*` for compile scripts, `xv*` for X11/xvue.
- **Variables**: 1–8 characters, uppercase in declarations, mixing French and English. Typical patterns:
  - `NB<noun>` for counts (`NBSOM` = number of vertices, `NBTRIA` = number of triangles).
  - `MX<noun>` for maxima / dimensions (`MXSOM`, `MXTRIA`, `MXTETR`).
  - `N<object>` for index/numbering arrays (`NSTRIA`, `NOTETR`, `NVOLTE`).
  - `L<object>` for list/table arrays (`LARETE`, `LEFACO`).
  - `PT<object>` for coordinates/point arrays (`PTXYZD`).
  - `TPS<...>` for time variables (`TPSINI`, `TPSFIN`) — **these must be declared with a consistent type** across modules; recent commit `88884ec` homogenized them to `REAL` in `flui/`. Any new `TPS*` variable must follow the same type.
  - `IERR` is the universal error-status output integer.
- **Include files** in `incl/` follow a sort-key scheme:
  - `a___<topic>.inc` — single-topic shared data (`a___arete.inc`, `a___face.inc`, `a___materiaux.inc`).
  - `a_<object>__<aspect>.inc` — object/aspect split (`a_ligne__bspline.inc`, `a_objet__definition.inc`).
  - `*.inc95` — Fortran 95 variants consumed by `_OMP` executables.
  - `a~_*.inc` — backup / alternative variant files kept in-tree (e.g. `a~_fonction__arbre.inc`).
  - Do not rename the `a___` prefix — many `INCLUDE` statements reference these filenames literally.

## Include files and shared state

- Every solver source file that needs a shared structure uses `INCLUDE "./incl/<name>.inc"` (note the relative path — it relies on compilation running from `$MEFISTO`).
- Shared state is carried in `COMMON` blocks declared inside the include files. The agent-level exploration found roughly 2,466 `COMMON` declarations across `incl/` and the solver directories — this is the project's universal data-sharing mechanism.
- `incl/homdir.inc` is **generated at build time** by `bin/cbl_tout` and encodes `$MEFISTO` as a Fortran `PARAMETER` string. Do not edit it manually; if you need the install path in a new include file, reference this generated file.
- When changing a common block, **recompile every module that includes it**. There is no dependency tracker in the shell build, so stale object files can silently desync — in doubt, run `bin/cbl_tout`.

## Error handling

- Traditional Fortran: each public subroutine takes an integer output `IERR` and sets it to 0 on success, non-zero on failure. Callers check `IERR` and propagate upward. Hard failures use `STOP` with a printed diagnostic.
- There is no structured exception mechanism and no logging framework. Diagnostics go to standard output / standard error via `WRITE(*, ...)` or `PRINT *`.
- `IMPLICIT NONE` is used inconsistently across the tree. New or heavily-edited files should prefer explicit typing, but do not mass-rewrite existing code just to add it.

## Build-script conventions (shell)

The `bin/cb*` scripts are the project's build system — there is no Makefile, CMake, or autotools. Their conventions:

- **Bash**, starting with `#!/bin/bash`.
- Always `cd $MEFISTO` first so that relative include paths (`./incl/...`) resolve.
- Create the `pp/` output directory if missing.
- Detect the language via `test -f $MEFISTO/td/m/anglais` and set `LANGAGE=0|1` to pick French or English messages.
- Ensure `xvue/xvuelc.o` is built (via `bin/ccxvue`) before linking anything that uses the X11 layer.
- Compile flags are baked into the scripts, not centralized. The canonical flags for mesher-class executables are:
  ```
  gfortran -Wall -mcmodel=large -m64 -O -fopenmp -I. \
      prpr/ppmail.f mail/lib util/lib xvue/lib xvue/xvuelc.o \
      -lgfortran -L/usr/X11R6/lib64 -lX11 -o pp/ppmail
  ```
  Note: `-mcmodel=large` is required because of the size of the static common blocks; `-fopenmp` is always enabled, even for non-OMP executables. The `<module>/lib` tokens are **pre-archived `.a` libraries** (or ordered object lists) produced by the per-module compile scripts — they are not CMake targets.
- The launcher scripts (`bin/INITIER`, `bin/MAILLER`, …) follow the same `LANGAGE` detection pattern and delegate execution to `pp/pp<module>`.

When adding a new solver or utility file, add its object to the matching module's `cb*` script — there is no auto-discovery.

## Bilingual identifiers and messages

MEFISTO originates from Alain Perronnet (LJLL, UPMC Paris), so much of the code uses French identifiers (`arete`, `face`, `sommet`, `volume`, `maille`, `contrainte`, `pression`, `vitesse`). English equivalents appear in newer or translated files but the French forms are load-bearing — **do not wholesale translate them** when editing. Keep new code in whichever language the surrounding file uses.

User-facing messages in launcher scripts and interactive menus are always bilingual, selected at runtime via the `td/m/anglais` marker file. When adding a prompt, always add both the French and English versions in the same `if [ $LANGAGE = 0 ]` branch structure used everywhere else.

## Miscellaneous

- **Executable permissions**: every script in `bin/`, every `.f` source, and every `.inc` file has the executable bit set in the tree (a historical convention from the installer). `chmod 755` on new files matches this.
- **No linting / formatter**: there is no `fprettify`, no `clang-format`, no pre-commit hook. Consistency is maintained by hand, using the existing surrounding code as the style reference.
- **Never Ctrl-C an interactive module**: exit via the `99;` menu command so project files are flushed cleanly. This is a user-level convention but worth preserving in any new interactive code.
