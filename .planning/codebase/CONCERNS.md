# Concerns

This document captures technical debt, fragile areas, known issues, and modernization risks in MEFISTO. It is meant as a reference for planning: things to be careful about, and candidate targets for future work.

## Legacy language and paradigms

- **Fortran 77 fixed-form dominates.** The vast majority of source in `mail/`, `elas/`, `flui/`, `ther/`, `reso/`, `util/`, and `xvue/` is strict F77 with column-based layout. Editors must handle the 72-column limit and column-6 continuations; copy-paste from modern editors can silently introduce column-8 tab bugs.
- **Mix of F77 and F95/OpenMP.** `_OMP` executables are built from Fortran 95 variants of the same files, with parallel include files in `incl/*.inc95`. The two worlds share source files, so a change in a common block can require an edit to both `*.inc` and `*.inc95`. Forgetting the 95 side silently desyncs the OpenMP build.
- **Heavy use of `GOTO` and computed branches.** Rough count from grep across the solver modules is in the ~13,000 range — this is the house style of the mesher and solver code and cannot be refactored away in isolated patches. Read-before-edit carefully; do not reformat surrounding blocks.
- **Inconsistent `IMPLICIT NONE`.** Many older files still rely on implicit typing (`I-N` for integers). New edits should prefer explicit declarations where the surrounding file allows, but do not mass-convert existing files — that is a separate initiative.
- **Tight `COMMON`-block coupling.** Shared state lives in ~2,500 `COMMON` declarations spread across `incl/*.inc`. A change in one common block touches every module that includes it, and because the shell build has no dependency tracker, stale `.o` files can produce silent mismatches. When in doubt, `bin/cbl_tout`.

## Build system debt

- **No Makefile, CMake, or autotools.** The build is a pile of hand-written bash scripts in `bin/` (and `bin.lnx64/`). Compile flags are baked into each script, there is no incremental compilation beyond checking `.o` timestamps, and there is no way to query the build graph.
- **Hard-coded paths.** `bin/cbmail` (and most other `cb*` scripts) link against `/usr/X11R6/lib64` directly. On modern Linux distributions this path is usually a compat symlink; on systems without the symlink the build fails until the path is patched by hand.
- **Hard-coded `-mcmodel=large`.** Required because of the static size of the common blocks. It means the code cannot target non-x86_64 architectures without additional work.
- **`-fopenmp` on every build.** OpenMP is enabled for non-parallel executables too — a latent footgun if a subroutine is ever marked `!$OMP` without being intended for the parallel path.
- **Install path baked into the binary.** `incl/homdir.inc` is generated at build time encoding `$MEFISTO` as a Fortran `PARAMETER` string. Relocating the install directory requires a full rebuild.
- **`pp/` was recently cleaned of committed binaries** (commit `ac282f8`). That is the right direction, but the `.gitignore` story is still young — watch for accidental re-adds of compiled artifacts.

## Graphics layer — X11/Motif

- **`xvue/` wraps Xlib directly.** `xvue/xvuelc.c` is ~150 KB of C code calling Xlib, with ~220 Fortran wrappers around it. Xlib itself is in long-term maintenance mode and is increasingly poorly supported under Wayland compositors.
- **Planned Qt migration.** `CLAUDE.md` explicitly lists the Qt replacement as a future initiative. Keep all new graphics-related calls inside `xvue/` so the migration surface stays bounded to this directory — do not sprinkle Xlib or window-system assumptions into solver code.
- **Dependency on `/usr/X11R6/lib64`.** Same path concern as the build system (above).
- **`libX11-dev` and `convert` are required at runtime**, which means MEFISTO cannot run headless or inside minimal containers without extra work.

## Numerical / algorithmic fragility

- **Type inconsistencies fixed recently.** Commit `88884ec` homogenized the type of the time variables `TPSINI`/`TPSFIN` to `REAL` in `flui/`. Similar inconsistencies likely still lurk elsewhere — any new `TPS*`, `DT*`, or `TIME*` variable is a candidate for careful review. **When introducing a new time variable, match the `REAL` type established by the fix.**
- **Very large source files.** A handful of files in `prpr/`, `flui/`, and `ther/` are in the 3,000–4,500-line range (e.g. `prpr/pppoba.f`). They concentrate a lot of logic and are disproportionate contributors to build time and review difficulty. Avoid further growth; extract helpers into `util/` when editing.
- **Mesh-quality algorithms are heuristic and bug-prone.** Routines like `mail/a1teqm.f` ("AMELIORER LA QUALITE DES TETRAEDRES") are the usual suspects when a mesh generation step fails on a corner case. Their headers describe input/output semantics precisely — **read the header before changing any parameter list**.
- **Secondary memory / paging.** Large meshes and matrices are paged to scratch files in the project directory. The on-disk format is implicit in the code; there is no versioning. Changing the layout of a paged structure breaks backward compatibility with existing saved projects.

## Error handling and observability

- **No structured logging.** Diagnostics are `WRITE(*, ...)` or `PRINT *` statements, with no severity levels, no timestamps, and no way to redirect per subsystem. Post-mortem debugging relies on scrolling through the terminal.
- **`IERR` discipline is informal.** Every subroutine takes a status integer, but there is no central catalogue of error codes. A `1` in one routine and a `1` in another mean different things. New code should at least document its error codes in the subroutine header.
- **`STOP` on fatal errors.** Deep routines can abort the whole process, which in the context of an interactive mesher means losing unsaved project state.

## Security and portability

- **X11 authentication** is whatever the host provides; there is no additional auth layer in `xvue/`. Running MEFISTO over an untrusted X display is as risky as any X11 client.
- **Shell scripts do not quote paths carefully.** A project directory name with spaces or shell metacharacters may break some launchers. Worth keeping in mind when fixing or extending `bin/*`.
- **No input validation in the mesher interactive loop.** Buffer-size parameters come from include files; user typos at the menu can produce confusing errors rather than clean rejections.
- **Portability limited to x86_64 Linux with gfortran.** No attempt is made to support other compilers (ifort, Cray, NAG) or architectures in the shipped build scripts. The `bin.lnx64/` tree is Linux-specific.

## Missing modern tooling

- **No CI.** No GitHub Actions, GitLab CI, or Jenkinsfile. All quality signals depend on the committer remembering to run `bin/cbl_tout` and the relevant `testa/` cases.
- **No automated test framework.** See `TESTING.md` — everything is manual end-to-end runs.
- **No static analysis.** No `-Wall` warning budget, no `fcheck=all` in the default flags, no `flint`, no `gfortran -fsanitize=*` coverage.
- **No code-coverage pipeline.**
- **No dependency lockfile** — you get whatever `gfortran`, `libX11-dev`, and `ImageMagick` the distribution provides. Upstream breakage lands at build time.

## Bilingual identifiers

- **French and English mixed freely.** Subroutines, common blocks, variable names, and comments are written in the language the original author was thinking in — mostly French, with English pockets. Grep/search is harder because there is no single keyword to look for (e.g. "edge" is sometimes `arete`, sometimes `edge`, occasionally `ARETE`). When fixing a bug, `grep` in both languages.
- **Risk for non-French-speaking contributors.** Understanding the mesher and solver code without at least basic technical French is significantly harder. This is not a bug — it is the project's heritage — but it is a real contribution barrier.

## Known recent-history issues

From the git log since the initial import:
- `88884ec` — **`TPSINI`/`TPSFIN` type homogenization** in `flui/`. Indicates that time-variable typing used to be inconsistent across files. Treat any new time-related code as high-risk until audited.
- `dbba1a4`, `ac282f8` — cleanup of committed binaries and `.gitignore` additions. The repository was shipped with built artifacts for a long time; be alert to compiled objects sneaking back in via new build scripts.

There is no issue tracker in-tree; `CLAUDE.md` states that the current active goal is **bug fixes without altering overall behavior**. That framing means:

- **Prefer minimal, local, surgical fixes** over sweeping refactors.
- **Preserve numerical output** unless the change is an explicit correctness fix — scientific users care about reproducibility across runs.

## Modernization risks and opportunities (informational, not prescriptive)

If/when modernization is attempted, the most load-bearing moves in rough order of risk:

1. **Introduce a real build system** (CMake or Meson) that calls into the existing object files without rewriting source — lowest disruption, biggest debugging win.
2. **Wire up CI** that runs `bin/cbl_tout` and a handful of small `testa/` cases. Even a smoke CI catches most regressions.
3. **Qt migration of `xvue/`** — contained if (and only if) no X11 leakage exists outside the directory.
4. **Incremental `IMPLICIT NONE` rollout**, one file at a time, alongside any touch for bug-fix work.
5. **Extract the largest `prpr/` and `flui/` files** into smaller units — only when a bug fix already forces touching them.

Each of these is a multi-phase effort; none should be attempted as a drive-by commit inside a bug-fix PR.
