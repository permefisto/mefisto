# Phase 0: Build skeleton & ABI stubs - Context

**Gathered:** 2026-04-10
**Status:** Ready for planning

<domain>
## Phase Boundary

Prove the Qt/CMake build plumbing integrates cleanly with the legacy shell linker and produces link-complete `pp/pp*_qt` executables **before any graphics logic exists**. This phase is load-bearing infrastructure for every later phase: it validates that `libxvueqt.a` (a Qt 6 static library built by CMake inside `xvue/qt/`) can be consumed by the existing shell-script Fortran linker, that the `extern "C"` ABI matches gfortran's trailing-underscore convention, that the X11 backend remains untouched, and that the validation baseline for every subsequent phase is captured. Not in scope for Phase 0: opening any window, drawing any pixel, reading any event, rendering any text, handling any color, exporting any image, or rewiring any menu.

Every `extern "C"` entry point in `xvue/xvuelc.c` gets a no-op stub in C++ with a byte-identical signature. Linking `pp/ppmail_qt`, `pp/ppelas_qt`, `pp/ppflui_qt`, `pp/ppther_qt`, and `pp/ppnlse_qt` must succeed. Running any of them must not crash when they call into a stub — they must proceed past the link stage and exercise the ABI.

</domain>

<decisions>
## Implementation Decisions

### Directory layout (A1)
- **D-01:** All new C++ sources, headers, and the scoped `CMakeLists.txt` live in a new `xvue/qt/` subfolder. The existing `xvue/*.f` wrappers (~224 files) and `xvue/xvuelc.c` stay untouched at the top level of `xvue/`.
- **D-02:** No file outside `xvue/qt/` is edited by Phase 0 *except* the new `bin/cbl_tout_qt` and `bin/cb*_qt` scripts. The rest of the repository — including `incl/`, `mail/`, `elas/`, `flui/`, `ther/`, `nlse/`, `reso/`, `util/`, `prpr/` — is read-only in this phase.
- **D-03:** `xvue/qt/README_COORDS.md` documents the audited Y-axis convention (origin location, Y-up vs Y-down, whether inversion happens in C or in Fortran) derived from reading `xvue/xvuelc.c`. This file is created in Phase 0 and referenced from every subsequent phase.

### ABI header organization (B1)
- **D-04:** One single header `xvue/qt/include/xvue_qt_api.h` declares every `extern "C"` entry point that exists in `xvue/xvuelc.c`. No category splits. This matches the one-big-blob nature of `xvuelc.c` itself and keeps grep/diff simple.
- **D-05:** The header uses `#define proc(x) x##_` (or an equivalent trailing-underscore macro) copied verbatim from `xvue/xvuelc.c` lines ~60–70 so that gfortran's trailing-underscore name mangling is honored identically across backends.
- **D-06:** Every scalar argument is declared as a pointer (`int*`, `float*`, `double*`). Every string argument is declared as `char* + int*` length pair in declared order (no hidden trailing length). Signatures are copied verbatim from `xvuelc.c` without simplification. A code-review checklist item enforces this: any `int` or `float` by value in the new header is a bug.
- **D-07:** The `XPoint` shim is declared inline in `xvue_qt_api.h` as `typedef struct { short x; short y; } MefistoPoint;` right next to the three entry points that use it (`xvface_`, `xvtraits_`, `xvfacetraits_`), with a `static_assert(sizeof(MefistoPoint) == 4, "MefistoPoint must match Xlib XPoint layout")` so byte drift fails at compile time.

### Build-system strategy (C1)
- **D-08:** `bin/cbl_tout_qt` is a literal clone of `bin/cbl_tout` with the `xvue/xvuelc.o` token replaced by `-Lxvue/qt/build -lxvueqt $(pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport) -lstdc++` on the final link line, and with a preceding step that runs `cmake -S xvue/qt -B xvue/qt/build && cmake --build xvue/qt/build`.
- **D-09:** Per-module clones `bin/cbmail_qt`, `bin/cbelas_qt`, `bin/cbflui_qt`, `bin/cbther_qt`, `bin/cbnlse_qt` exist alongside their `bin/cb*` originals. They follow the same language-detection (`$MEFISTO/td/m/anglais`) and directory-creation pattern as the existing scripts. No conditional logic is added to the original scripts.
- **D-10:** `bin/cbl_tout_qt` cleans `xvue/qt/build/` and `pp/` before every run so that stale `.o` files from the Fortran side (a known CONCERNS.md fragility) cannot mask a broken Qt build. This is cheaper than adding dependency tracking.
- **D-11:** `xvue/qt/CMakeLists.txt` sets `CMAKE_POSITION_INDEPENDENT_CODE ON`, `CMAKE_CXX_STANDARD 17`, `CMAKE_CXX_STANDARD_REQUIRED ON`, `CMAKE_AUTOMOC ON`, and explicitly `target_compile_options(xvueqt PRIVATE -fno-openmp)`. It sets `CMAKE_AUTOMOC` **before** `find_package(Qt6 ...)`. It requires `cmake_minimum_required(VERSION 3.21)`.
- **D-12:** A CMake custom target `verify_abi` runs `nm libxvueqt.a | grep -c '_$'` at the end of the build and compares the count against the number of `extern "C"` declarations grepped from `xvue_qt_api.h`. The build fails if any Fortran-facing symbol is missing the trailing underscore or if the count drifts.
- **D-13:** `xvue/qt/CMakeLists.txt` exposes an installable target: `install(TARGETS xvueqt DESTINATION xvue/qt/lib)` (or an `EXPORT` so the shell scripts can locate the `.a` by a stable relative path). The shell scripts assume the library is at `xvue/qt/build/libxvueqt.a`.

### Validation baseline (D-rec)
- **D-14:** Five canonical `testa/` cases form the Phase 0 validation baseline, recorded in `.planning/validation/BASELINE.md`:
  - **Mesher** — `testa/pan2d` (interactive 2D panel — exercises the mesher's pick-and-drag workflow)
  - **Elasticity** — `testa/nafems_le1` (NAFEMS LE1 benchmark — standard linear-elasticity reference case)
  - **Fluid** — `testa/cavity2d` (classic lid-driven cavity — widely-recognized visual output)
  - **Thermal** — `testa/heat1d` (1D heat conduction — simplest diagnostic case for the thermal solver)
  - **Nonlinear** — `testa/nlsecu` (nonlinear solver canonical case)
- **D-15:** For Phase 0, the validation gate is: all 5 baseline cases **still run successfully on the legacy X11 build** (`bin/cbl_tout` + `pp/pp*`). The Qt build does not yet produce any graphics output — running any `pp/pp*_qt` succeeds if it proceeds past the link stage and exercises the no-op ABI stubs without crashing. Every subsequent phase (1–7) runs the same 5 cases through **both** backends at the end of the phase and logs results to `.planning/phases/{padded}-*/VALIDATION.md`.
- **D-16:** The baseline document `.planning/validation/BASELINE.md` lists, for each of the 5 cases: the project directory path, the launcher script(s) that run it, the expected qualitative behavior (what a human eye looks for), and any known-flaky touchpoints. This is write-once in Phase 0 and read-only for every subsequent phase.

### Stub diagnostic behavior (E1)
- **D-17:** Every no-op stub in `xvue/qt/src/xvue_qt_api.cpp` prints `"xvue-qt: stub <function_name> not implemented yet\n"` to `stderr` on its **first** call, then stays silent for the remainder of the process. Implementation uses a per-function `static bool warned = false;` guard or a shared `std::once_flag`. This makes Phases 1–7 immediately spot accidental regressions ("I forgot to implement `sauvefenetre_`") without flooding the terminal during a long mesher session.
- **D-18:** Stubs return `void` (most entry points are void) or a sensible default (0, 0.0, or a null-equivalent) for non-void entry points. Stubs do **not** abort, do **not** set error flags, do **not** touch any Qt object. They are the thinnest possible ABI placeholder.

### Claude's Discretion

- **Debug thread-affinity assertion macro** — define `XVUE_QT_ASSERT_MAIN_THREAD()` as `Q_ASSERT(QThread::currentThread() == qApp->thread())` in debug builds (`#ifdef QT_DEBUG`) and empty in release. For Phase 0 stubs it's a no-op because there is no `qApp` yet, but the macro is declared in `xvue_qt_api.h` so Phase 1+ can drop it into every entry point body.
- **File naming and internal organization inside `xvue/qt/`** — the expected layout is `xvue/qt/CMakeLists.txt`, `xvue/qt/include/xvue_qt_api.h`, `xvue/qt/src/xvue_qt_api.cpp` (the stubs), `xvue/qt/README_COORDS.md`, `xvue/qt/build/` (gitignored). Additional files land in Phase 1+.
- **Stub categorization for future phases** — the ~60 stubs in `xvue_qt_api.cpp` may be grouped by category with `// === Phase 2: Drawing primitives ===` style comment banners to make later phases easy to navigate, but this is a readability choice, not a commitment.
- **pkg-config vs find_package for Qt detection in shell scripts** — the final shell linker line uses `$(pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport)` for portability; CMake itself uses `find_package(Qt6 ... COMPONENTS Widgets Gui Core PrintSupport)` internally. Both mechanisms are acceptable.
- **Formatting of `xvue/qt/README_COORDS.md`** — short free-form Markdown documenting what the Y-axis audit found. Content matters; layout doesn't.

</decisions>

<specifics>
## Specific Ideas

- **"Keep the shell-script convention"** — the user explicitly chose C1 clone-and-modify over C2 parametrize, signalling a preference for the existing one-script-per-purpose pattern used throughout `bin/`. Downstream phases should respect this: do not introduce `case $BACKEND in` constructs into the shell scripts.
- **Warn-once stub behavior (E1)** is the user's explicit choice over pure silence (E2) and over hard-abort (E3). The signal: Phase 0 must be instantly diagnostic without being noisy, and must not break the link-and-run sanity check.
- **`xvue/qt/` as a fresh subfolder (A1)** signals that the user wants the Qt port to be visually isolated from the legacy Fortran wrappers on disk, making git history, blame, and eventual retirement (Phase 9) as surgical as possible.
- **Single ABI header (B1)** signals optimization for grep/diff against `xvuelc.c` rather than conceptual cleanliness — the comparability with the existing C file is the load-bearing property.
- **`testa/` baseline choices (D-rec)** are anchored to benchmark cases with recognizable reference outputs (NAFEMS LE1, lid-driven cavity, 1D heat) so visual A/B drift is diagnosable by eye.

</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing Phase 0.**

### Project-level specs
- `.planning/PROJECT.md` — Core Value (Fortran-must-not-notice invariant), scope boundaries, Key Decisions table (Qt 6, CMake-scoped-to-`xvue/`, parallel X11 build, etc.)
- `.planning/REQUIREMENTS.md` §"Build — CMake skeleton and ABI stubs" — the 10 BUILD-* requirements this phase delivers (BUILD-01..10)
- `.planning/ROADMAP.md` §"Phase 0: Build skeleton & ABI stubs" — the canonical phase boundary, goal, requirements list, and success criteria

### Research synthesis
- `.planning/research/SUMMARY.md` §"Phase 0" — reconciled phase structure and rationale for front-loading 9 pitfalls into this phase
- `.planning/research/STACK.md` — the copy-pasteable CMake template, the gfortran trailing-underscore convention confirmation, the Fortran↔C++ ABI rules, the "what NOT to use" list (no Qt 5, no QML, no qmake, no vendored Qt)
- `.planning/research/ARCHITECTURE.md` §"Components" — the 12-file layout the Qt layer will eventually grow into (only the ABI shim and the singleton plumbing land in Phase 0)
- `.planning/research/PITFALLS.md` — specifically Pitfalls 1 (name mangling), 2 (pointer args), 3 (string lengths), 4 (`XPoint` shim), 5 (`QApplication` singleton — mostly Phase 1 but the plumbing starts here), 9 (AUTOMOC), 10 (`-fPIC`), 11 (`-fopenmp`), 16 (Y-axis convention), 18 (clean-build discipline), 20 (regression drift / validation baseline). Phase 0 front-loads nine of the twenty pitfalls.

### Codebase maps
- `.planning/codebase/ARCHITECTURE.md` — shared-data layer (`incl/` common blocks), launcher + entry-point + solver-module + utility + graphics layering, data flow, existing build scripts
- `.planning/codebase/STRUCTURE.md` — `xvue/` directory layout (where `xvuelc.c` lives, ~224 `.f` wrappers), `bin/` layout for launchers and `cb*` compile scripts
- `.planning/codebase/STACK.md` — existing compiler flags, gfortran version, gcc version, Qt package names on Debian trixie, X11 library locations
- `.planning/codebase/CONVENTIONS.md` — shell script conventions in `bin/` (language-detection with `$MEFISTO/td/m/anglais`, `cd $MEFISTO` first, etc.) and the Fortran fixed-form norms that Phase 0 must not disturb
- `.planning/codebase/CONCERNS.md` — stale `.o` fragility driving D-10, hard-coded `/usr/X11R6/lib64` path informing D-11, bilingual identifiers, single-developer no-CI risk driving D-14/15

### Direct source
- `/home/drico/git/mefisto/xvue/xvuelc.c` — **the** reference file. Every entry-point name, signature, and the `proc()` macro live here. Phase 0 must copy these byte-identically into `xvue/qt/include/xvue_qt_api.h`. Also carries the hand-rolled PostScript emitter that Phase 7 will preserve verbatim (not yet relevant here) and the Y-axis handling that Phase 0 must audit (D-03).
- `/home/drico/git/mefisto/bin/ccxvue` — current compile command for `xvuelc.c` (`gcc -Wall -O -m64 -c -o xvue/xvuelc.o xvue/xvuelc.c -I/usr/X11R6/include`). The `-I`/`-m64`/`-O` flags informing CMake target flags.
- `/home/drico/git/mefisto/bin/cbl_tout` and `/home/drico/git/mefisto/bin/cbmail` — the canonical shell-script patterns that `bin/cbl_tout_qt` and `bin/cb*_qt` clones must follow (language detection, `cd $MEFISTO`, creating `pp/` if missing, `gfortran -Wall -mcmodel=large -m64 -O -fopenmp -I. ... -o pp/pp*` invocation, error messages).
- `/home/drico/git/mefisto/CLAUDE.md` — working rules: never break the X11 build, run the smallest relevant test after every change, ask before installing system packages, commit after every logical step.

### External references (not in repo — read from upstream docs if needed during planning)
- Qt 6 upstream CMake integration guide (`doc.qt.io/qt-6/cmake-get-started.html`, `cmake-manual.html`) — the authoritative source for `find_package(Qt6 COMPONENTS ...)`, `CMAKE_AUTOMOC`, `target_link_libraries(Qt6::Widgets)`.
- Debian trixie package list: `qt6-base-dev` 6.10.2, `qt6-base-dev-tools` 6.10.2, `libqt6imageformats6-plugins`, `cmake` 3.31+ — version pins confirmed at init time, not re-checked in Phase 0 planning unless apt state drifts.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets

- **`xvue/xvuelc.c`** — the ground-truth source of every `extern "C"` entry-point name, signature, `XPoint` usage, and the `#define proc(x) x##_` trailing-underscore macro. Phase 0 copies these **literally** into `xvue/qt/include/xvue_qt_api.h`; any divergence is a Phase 0 bug.
- **`xvue/*.f`** (~224 Fortran wrapper files) — not touched by Phase 0, but they constrain the ABI surface. They call into the `extern "C"` names defined in `xvuelc.c`; the Qt backend must expose the same names with the same byte-level signatures, or the link fails.
- **`bin/cbl_tout`, `bin/cbmail`, `bin/cbelas`, `bin/cbflui`, `bin/cbther`, `bin/cbnlse`** — the existing per-module compile-and-link scripts. Phase 0 clones them into `bin/cbl_tout_qt`, `bin/cbmail_qt`, etc. (D-08, D-09). Each clone keeps the same language-detection and error-reporting idioms and differs only on the final linker line.
- **`bin/ccxvue`** — current `xvue/xvuelc.c` compile script. Carries the `-I/usr/X11R6/include` X11 path that the Qt clones do not need. Serves as the closest template for "what Phase 0 is replacing."
- **`.planning/codebase/` maps** — codebase-analysis documents written during brownfield init; read-only inputs to phase planning.
- **Existing `testa/` baseline cases** (`testa/pan2d`, `testa/nafems_le1`, `testa/cavity2d`, `testa/heat1d`, `testa/nlsecu`) — self-contained project directories with input files already staged; Phase 0 registers them in `.planning/validation/BASELINE.md` without modifying anything inside them.

### Established Patterns

- **One-script-per-purpose shell discipline in `bin/`** — every `bin/cb*` file is a self-contained bash script with `#!/bin/bash`, language detection, `cd $MEFISTO`, per-step echo output, and explicit error messages. The Phase 0 Qt variants inherit this exactly; no conditional branches for "which backend" inside shared scripts (per user decision C1).
- **Trailing-underscore ABI via `#define proc(x) x##_`** — the existing convention in `xvuelc.c`. The Qt header reuses the same macro so gfortran's name-mangling assumption holds byte-for-byte.
- **`static` file-scope for per-call-site state** — `xvuelc.c` uses file-scope `static` variables extensively. The Phase 0 stub file does not need them (stubs are stateless except for the per-function `warned` flag for E1), but Phase 1+ will restore them inside `XvueState`.
- **Bilingual language detection via `$MEFISTO/td/m/anglais`** — every `bin/` script checks for this file to choose French vs. English messages. The Phase 0 Qt variants preserve the same detection, even though they produce only English/French build-time diagnostics.
- **Fortran 77 fixed-form surrounding the C bridge** — Phase 0 does not touch any Fortran source; the wrappers in `xvue/*.f` are read-only. This is the Core Value invariant made concrete.

### Integration Points

- **Linker line in `bin/cb*_qt` scripts** — the only touch point between the Fortran build and the Qt layer. The Qt clones replace `xvue/xvuelc.o` in the existing `gfortran ... -lX11 -o pp/pp*` line with a Qt-linked static-library reference (D-08). Every phase after 0 only **grows** `libxvueqt.a` — the linker line does not change again until Phase 9 retirement.
- **CMake inside shell** — the Qt clones call `cmake -S xvue/qt -B xvue/qt/build && cmake --build xvue/qt/build` as a preliminary step before the gfortran link. This is the **only** place where CMake touches the build.
- **`xvue/qt/build/libxvueqt.a`** — the handoff artifact. Produced by CMake, consumed by shell scripts. Path is stable and referenced by relative path from `$MEFISTO`.
- **`xvue/qt/README_COORDS.md`** — produced by Phase 0, consumed by every phase that touches drawing (Phases 1–7). Encodes the Y-axis convention audit result as a read-only reference.
- **`.planning/validation/BASELINE.md`** — produced by Phase 0, consumed by every phase for end-of-phase validation. Stable list of 5 test cases.

</code_context>

<deferred>
## Deferred Ideas

- **Full CMake migration of the Fortran build** — explicitly out of scope per PROJECT.md. Not reconsidered in Phase 0.
- **Automated CI (GitHub Actions / GitLab CI)** — out of scope per PROJECT.md. Not reconsidered in Phase 0. Manual validation baseline (D-14, D-15) replaces CI for this project.
- **`pkg-config` wrapper or a generated `xvue/qt/qt_link_flags.txt`** — considered as a cleaner handoff than inlining `$(pkg-config --libs ...)` into each shell script, rejected for Phase 0 because it adds moving parts without commensurate benefit. Revisit if the shell scripts become hard to maintain during Phase 2–6.
- **Stub categorization banner comments in `xvue_qt_api.cpp`** — nice-to-have for readability, not a Phase 0 commitment. Will be added ad-hoc as later phases land.
- **`dpkg -l | grep qt6-base` runtime check in `bin/cbl_tout_qt`** — helpful for users on machines without Qt 6 installed, but out of scope for Phase 0. The script assumes `qt6-base-dev` is present and fails with the raw CMake error if not. Can be added later as a one-line preflight.
- **Windows/macOS CMake targets in `xvue/qt/CMakeLists.txt`** — out of scope per PROJECT.md (Linux x86_64 only). Not added to the CMakeLists even as guards.
- **Integration-testing the `verify_abi` CMake target against the actual `xvuelc.c` symbol count** — a nice safety net, but the user's baseline for "correct" is already "count matches the declaration count in `xvue_qt_api.h`"; cross-checking against `xvuelc.c` itself is redundant if the header is copied faithfully.
- **Per-executable link-time checks** (beyond the five baseline executables) — out of scope. Phase 0 only requires that `pp/ppmail_qt`, `pp/ppelas_qt`, `pp/ppflui_qt`, `pp/ppther_qt`, and `pp/ppnlse_qt` link and run. Other `pp/pp*` targets (e.g. `ppadap`, `pppoba`, `ppdvsr`) are out of scope until the user explicitly asks for them.

### Reviewed Todos (not folded)
None — no pending todos matched Phase 0 scope at init time.

</deferred>

---
*Phase: 00-build-skeleton-abi-stubs*
*Context gathered: 2026-04-10*
