# Testing

## Overview

MEFISTO has **no automated test framework, no CI, and no unit tests**. Correctness is verified by running a collection of **end-to-end demo projects** by hand, through the same launcher scripts that end users use, and comparing the output (meshes, solver results, plots) against expected behavior documented in the tutorials.

Per `CLAUDE.md`:
- Small tests in `testa/` / `testf/` **must continue to pass** after every change.
- For large or long-running tests, ask the user to run them before declaring a change complete.
- Always prefer the smallest relevant test case when checking a fix.

## Test directories

### `testa/` — English test suite

~60 project directories, each a complete MEFISTO project that exercises one or more modules. Representative cases:

| Test | Exercises |
|---|---|
| `testa/demoa/` | General demo — `el`, `ln`, `ob`, `pt`, `sf`, `th`, `vl1`, `vl2` subcases covering elements, lines, objects, points, surfaces, thermal, volumes |
| `testa/cavity2d/`, `testa/cavity3d/` | 2D/3D cavity flow — fluid solver (`FLUIDER`) against the classic lid-driven cavity benchmark; files like `cavity2d.meshbf`, `cavity2d.meshth`, `cavity2d.stoke56cr` |
| `testa/heat1d/`, `testa/eigv3dheat/` | Thermal solver (`THERMICER`) — 1D heat conduction and 3D eigenvalue heat problems |
| `testa/nafems_le1/`, `testa/nafems_t4/`, `testa/nafems_te31/`, `testa/nafems_te51/` | NAFEMS benchmark cases — canonical elasticity validation problems |
| `testa/hexa/`, `testa/hexahedron/`, `testa/hexapyra/`, `testa/pyra/`, `testa/torus/`, `testa/torusstl/` | Mesher (`MAILLER`) — specific element types and geometries |
| `testa/tube2d/`, `testa/tube3d/`, `testa/pipe3d/`, `testa/pan2d/`, `testa/pan3d/` | Pipe / panel geometries — mesher + solver |
| `testa/nlsecu/`, `testa/nlsecuvv/`, `testa/nlsecylbreak/`, `testa/nlsesqri/`, `testa/nlsgpe/` | Nonlinear solver (`NLSER`) |
| `testa/wave/`, `testa/solitonsph20/` | Wave solver |
| `testa/step/`, `testa/marche/`, `testa/rectbou/`, `testa/recth/`, `testa/recyl/`, `testa/rol2dbfcg/`, `testa/rol2dthcg/`, `testa/rol3dbfcg/`, `testa/rol3dthcg/` | Coupled thermal/fluid/elasticity, rolling contact, etc. |
| `testa/1dwall/`, `testa/2cubFOAM/`, `testa/habc/`, `testa/cylball/`, `testa/gourd/`, `testa/lane3d/`, `testa/gpesigc/`, `testa/optionput/`, … | Specific bug reproductions and miscellaneous scenarios |

### `testf/` — French test suite

Same structure and purpose as `testa/`, with French project names (`bielle/`, `bm/`, `boule/`, `cardan/`, `carrecercle/`, …). Many cases are French-language duplicates of their English counterparts; others are unique.

The choice of suite is driven by the language flag (`$MEFISTO/td/m/anglais`); most fixes should be validated against both when a matching case exists.

## Tutorials and demo data — `td/`

`td/` is not strictly a test directory, but tutorials are often the easiest way to exercise a module end-to-end. Layout:

```
td/
├── da/   English demo data (mesh files, geometry files, solver inputs)
├── df/   French demo data
├── ia/   English initial data (bootstrap data for INITIER)
├── if/   French initial data
├── ma/   English menu scripts (guided walkthroughs)
├── mf/   French menu scripts
├── d/, i/, m/, g/, p/ — shared or language-neutral data
└── m/anglais   marker file whose presence flips launchers to English
```

The `td/m{a,f}/` menu files drive interactive sessions: they contain the same text commands a user would type, so they serve as scripted smoke tests for the interactive modules. Running a menu file against `MAILLER` is a common way to check that a mesher change has not regressed the tutorial flow.

## How tests are run

There is no `make test` equivalent. To run a test case:

1. Ensure the environment is set up (`$MEFISTO`, `$MEFISTOX`, `PATH`, `CDPATH`) — see `CLAUDE.md`.
2. `cd $MEFISTOX` (or wherever you want the project to live).
3. Copy or symlink the test project from `testa/<name>/` or `testf/<nom>/` into a working directory under `$MEFISTOX`.
4. Run the relevant launcher: `INITIER <project>`, then `MAILLER <project>`, then `ELASTICER <project>` / `FLUIDER <project>` / etc., in the order the test expects.
5. Exit each interactive module with `99;` — never Ctrl-C.
6. Inspect the resulting mesh / solution / plot files in the project directory, and compare against the tutorial output documented under `doc/` or `doca/` / `docf/`.

For interactive mesher tests, the menu files in `td/ma/` and `td/mf/` can be redirected into the launcher to automate the user input.

## Build-time "tests"

The coarsest but most important check is that **the build still succeeds**. Per `CLAUDE.md`:

- After every change, verify that the affected module compiles with its specific `bin/cb*` script (e.g. `bin/cbmail` after editing `mail/`).
- Before committing, the full build (`bin/cbl_tout`) must succeed.

Because the shell build has no dependency tracker, a changed include file or common block can produce silent stale-object-file mismatches. When in doubt, rebuild everything with `bin/cbl_tout` to validate.

The standalone graphics smoke tests `prpr/xvtest1.f`..`prpr/xvtest4.f` are the closest thing to unit tests in the repo: they compile into tiny executables that open an X11 window and exercise a subset of `xvue/` primitives. They are useful to confirm that the X11 layer still links and draws after changes there.

## Coverage gaps and limitations

- **No numerical accuracy regression harness** — no recorded "golden" output files to diff against automatically. Comparisons are visual or informal.
- **No unit tests** — subroutines cannot be exercised in isolation; every check goes through a full solver run.
- **No continuous integration** — there is no GitHub Actions config, no `.gitlab-ci.yml`, no Jenkinsfile. Quality is maintained by the committer running the relevant test cases by hand.
- **No code coverage** — `gcov` / `gfortran -fprofile-arcs` are not wired into the build scripts.
- **Long-running tests** — several of the 3D tests (`testa/torus/`, `testa/lane3d/`, `testa/rol3dthcg/`, `testa/eigv3dheat/`) take minutes to hours. Per `CLAUDE.md`, delegate those to the user.

## Adding a new test case

1. Pick a name that describes the scenario (`testa/<shortname>/`).
2. Add the minimum set of input files needed (mesh geometry, material data, solver input files).
3. If the test is a bug reproduction, include a `README` or short text file documenting the expected behavior and the commit / issue it guards against.
4. Add a matching French case under `testf/` only if the bug is language-sensitive.
5. Do **not** add expected-output binary dumps — they are fragile across platforms and solver versions. Describe expected qualitative behavior instead.
