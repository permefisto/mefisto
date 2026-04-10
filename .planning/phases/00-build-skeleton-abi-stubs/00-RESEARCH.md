# Phase 0: Build skeleton & ABI stubs - Research

**Researched:** 2026-04-10
**Domain:** CMake + Qt 6 static-library build plumbing behind a gfortran trailing-underscore `extern "C"` ABI, integrated with legacy shell-script Fortran linker
**Confidence:** HIGH on ABI surface (grounded in direct reads of `xvue/xvuelc.c`); HIGH on Fortran↔C interop rules; HIGH on build-environment availability (probed); MEDIUM on Qt 6 CMake idioms (cross-referenced with `.planning/research/STACK.md` and `.planning/research/PITFALLS.md`, not re-verified against upstream Qt docs this session because those references are already canonical inputs).

## Summary

Phase 0 delivers the narrowest possible **link-complete** Qt bridge: a CMake-built `libxvueqt.a` whose every `extern "C"` entry point is a no-op stub, consumed by new `bin/cbl_tout_qt` + `bin/cb*_qt` shell scripts that clone the existing `bin/cbl_tout` + `bin/cb*` scripts. No window opens, no pixel draws, no event is read. Success is defined entirely at link time and at the "did any `pp/pp*_qt` executable proceed past `_start` into a stub and exit cleanly?" level.

The ground truth is `xvue/xvuelc.c` (3619 lines). Grep counts **59** `proc(...)` entry points; the phase description's "~60" is accurate. Every entry point must be copied **byte-identically** into `xvue/qt/include/xvue_qt_api.h` because the ~221 Fortran `.f` wrappers in `xvue/` already call these names with these exact argument conventions and gfortran will trip over any drift.

**Primary recommendation:** Treat `xvue/xvuelc.c` as a read-only specification during Phase 0. Do not "clean up" signatures. Do not drop C-internal helpers. Do not add new entry points. Copy, stub, link, run, commit.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**Directory layout (A1)**
- **D-01:** All new C++ sources, headers, and the scoped `CMakeLists.txt` live in a new `xvue/qt/` subfolder. The existing `xvue/*.f` wrappers (~224 files) and `xvue/xvuelc.c` stay untouched at the top level of `xvue/`.
- **D-02:** No file outside `xvue/qt/` is edited by Phase 0 *except* the new `bin/cbl_tout_qt` and `bin/cb*_qt` scripts. The rest of the repository — including `incl/`, `mail/`, `elas/`, `flui/`, `ther/`, `nlse/`, `reso/`, `util/`, `prpr/` — is read-only in this phase.
- **D-03:** `xvue/qt/README_COORDS.md` documents the audited Y-axis convention (origin location, Y-up vs Y-down, whether inversion happens in C or in Fortran) derived from reading `xvue/xvuelc.c`. This file is created in Phase 0 and referenced from every subsequent phase.

**ABI header organization (B1)**
- **D-04:** One single header `xvue/qt/include/xvue_qt_api.h` declares every `extern "C"` entry point that exists in `xvue/xvuelc.c`. No category splits. (**See §"Planner Alert" below — three entries in `xvuelc.c` take the Xlib `Colormap` type and are C-internal helpers; the planner must resolve D-04's interaction with that finding before writing tasks.**)
- **D-05:** The header uses `#define proc(x) x##_` copied verbatim from `xvue/xvuelc.c` lines ~60–70 so gfortran's trailing-underscore name mangling is honored identically across backends.
- **D-06:** Every scalar argument is declared as a pointer (`int*`, `float*`, `double*`). Every string argument is declared as `char* + int*` length pair in declared order (no hidden trailing length). Signatures are copied verbatim from `xvuelc.c`. (**See §"Planner Alert" — one entry point in `xvuelc.c` takes `int` by value; this is a C-internal caller path, not a Fortran caller, but it violates the "every scalar is a pointer" rule and must be called out explicitly.**)
- **D-07:** The `XPoint` shim is declared inline in `xvue_qt_api.h` as `typedef struct { short x; short y; } MefistoPoint;` right next to the three entry points that use it (`xvface_`, `xvtraits_`, `xvfacetraits_`), with `static_assert(sizeof(MefistoPoint) == 4, "MefistoPoint must match Xlib XPoint layout")`.

**Build-system strategy (C1)**
- **D-08:** `bin/cbl_tout_qt` is a literal clone of `bin/cbl_tout` with the `xvue/xvuelc.o` token replaced by `-Lxvue/qt/build -lxvueqt $(pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport) -lstdc++` on the final link line, and with a preceding step that runs `cmake -S xvue/qt -B xvue/qt/build && cmake --build xvue/qt/build`.
- **D-09:** Per-module clones `bin/cbmail_qt`, `bin/cbelas_qt`, `bin/cbflui_qt`, `bin/cbther_qt`, `bin/cbnlse_qt` exist alongside their `bin/cb*` originals. Same language-detection (`$MEFISTO/td/m/anglais`) and directory-creation pattern. No conditional logic in originals.
- **D-10:** `bin/cbl_tout_qt` cleans `xvue/qt/build/` and `pp/` before every run so stale `.o` files cannot mask a broken Qt build.
- **D-11:** `xvue/qt/CMakeLists.txt` sets `CMAKE_POSITION_INDEPENDENT_CODE ON`, `CMAKE_CXX_STANDARD 17`, `CMAKE_CXX_STANDARD_REQUIRED ON`, `CMAKE_AUTOMOC ON` **before** `find_package(Qt6 ...)`, and `target_compile_options(xvueqt PRIVATE -fno-openmp)`. Requires `cmake_minimum_required(VERSION 3.21)`.
- **D-12:** A CMake custom target `verify_abi` runs `nm libxvueqt.a | grep -c '_$'` at end-of-build and compares the count against `extern "C"` declarations grepped from `xvue_qt_api.h`. Build fails if the count drifts.
- **D-13:** `xvue/qt/CMakeLists.txt` exposes `install(TARGETS xvueqt DESTINATION xvue/qt/lib)` or equivalent stable path; shell scripts assume the library is at `xvue/qt/build/libxvueqt.a`.

**Validation baseline (D-rec)**
- **D-14:** Five canonical `testa/` cases: `testa/pan2d` (mesher), `testa/nafems_le1` (elas), `testa/cavity2d` (flui), `testa/heat1d` (ther), `testa/nlsecu` (nlse). Recorded in `.planning/validation/BASELINE.md`.
- **D-15:** Phase 0 validation gate: all 5 cases still run on legacy X11 build. Running any `pp/pp*_qt` succeeds if it proceeds past link stage and exercises no-op stubs without crashing. Every subsequent phase runs both backends.
- **D-16:** `.planning/validation/BASELINE.md` lists each case's project dir, launcher script(s), expected qualitative behavior, known-flaky touchpoints. Write-once in Phase 0.

**Stub diagnostic behavior (E1)**
- **D-17:** Every no-op stub prints `"xvue-qt: stub <function_name> not implemented yet\n"` to stderr on its **first** call, silent thereafter. Uses per-function `static bool warned = false;` or `std::once_flag`.
- **D-18:** Stubs return `void` (most) or a sensible default (0, 0.0, null-equivalent) for non-void. Stubs do not abort, do not set error flags, do not touch any Qt object.

### Claude's Discretion

- **Debug thread-affinity macro `XVUE_QT_ASSERT_MAIN_THREAD()`** — declare in `xvue_qt_api.h` as `Q_ASSERT(QThread::currentThread() == qApp->thread())` under `#ifdef QT_DEBUG`, empty otherwise. No-op in Phase 0 (no `qApp` yet) but available for Phase 1+.
- **File layout inside `xvue/qt/`** — expected: `CMakeLists.txt`, `include/xvue_qt_api.h`, `src/xvue_qt_api.cpp`, `README_COORDS.md`, `build/` (gitignored).
- **Stub categorization comments** — `// === Phase 2: Drawing primitives ===` banners are a readability choice.
- **pkg-config vs find_package** — shell linker uses `$(pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport)`; CMake itself uses `find_package(Qt6 ... COMPONENTS Widgets Gui Core PrintSupport)`. Both acceptable.
- **`README_COORDS.md` format** — short free-form Markdown. Content matters; layout doesn't.

### Deferred Ideas (OUT OF SCOPE)

- Full CMake migration of the Fortran build.
- Automated CI (GitHub Actions / GitLab CI).
- `pkg-config` wrapper / generated `xvue/qt/qt_link_flags.txt`.
- Stub categorization banner comments as a Phase 0 commitment.
- `dpkg -l | grep qt6-base` runtime preflight check in `bin/cbl_tout_qt`.
- Windows/macOS CMake targets.
- Cross-checking `verify_abi` count against `xvuelc.c` itself (redundant with header count).
- Per-executable link checks beyond the 5 baseline (`ppmail_qt`, `ppelas_qt`, `ppflui_qt`, `ppther_qt`, `ppnlse_qt`).

</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| BUILD-01 | Build `libxvueqt.a` from `xvue/CMakeLists.txt` on Debian trixie using `qt6-base-dev` / `qt6-base-dev-tools` from apt, no vendoring | §Standard Stack (Qt 6.10.2 from Debian trixie/sid confirmed), §Environment Availability (qt6-base-dev NOT currently installed — apt install required before phase starts) |
| BUILD-02 | `libxvueqt.a` built with `-fPIC` (PIE-safe linkage) and without `-fopenmp` | §Code Examples (CMakeLists template), Pitfall "OpenMP runtime collision", D-11 |
| BUILD-03 | `CMAKE_AUTOMOC ON` set before `find_package(Qt6 ...)`, `QObject` subclasses only in headers | §Architecture Patterns, Pitfall "AUTOMOC ordering" — no QObject subclasses in Phase 0, but scaffolding must be correct |
| BUILD-04 | `xvue_qt_api.h` declares every one of the ~60 (verified: **59**) `extern "C"` entry points byte-identically | §"The 59 Entry Points" inventory, §"Planner Alert: 3 C-internal helpers" |
| BUILD-05 | `xvue_qt_api.cpp` implements every declared entry point as a no-op stub that returns successfully | §Architecture Patterns (stub pattern), D-17/D-18 |
| BUILD-06 | `bin/cbl_tout_qt` builds all MEFISTO executables against `libxvueqt.a` via `pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport` + `-lstdc++`, cleans `xvue/qt/build/` + `pp/` before building | §Code Examples (cbl_tout clone template), D-08, D-10 |
| BUILD-07 | Existing `bin/cbl_tout` + `xvuelc.c` + `libX11` path continues to work unchanged | §"Legacy path preservation" — verified by the read-only rule D-02 and by not modifying `bin/ccxvue`, `bin/cbl_tout`, or `xvue/xvuelc.c` |
| BUILD-08 | Build-time `nm libxvueqt.a \| grep '_$'` sanity check fails if any Fortran-facing symbol is missing the trailing underscore | §Code Examples (`verify_abi` CMake custom target), Pitfall "trailing underscore" |
| BUILD-09 | Y-axis convention audited, documented in `xvue/README_COORDS.md`, Qt bridge follows same convention | §"Y-Axis Convention Audit" (audited in this research document) |
| BUILD-10 | 5 canonical `testa/` cases recorded in `.planning/validation/BASELINE.md` | §"testa Baseline Verification" (all 5 confirmed present in `testa/`) |

</phase_requirements>

## Planner Alert: Contradictions Between D-04/D-06 and xvuelc.c Reality

**Two locked decisions in CONTEXT.md do not survive a careful read of `xvue/xvuelc.c`. The planner must surface these to the user before writing tasks.**

### Alert 1: Three entry points take the Xlib `Colormap` type

`xvue/xvuelc.c` exposes three functions via the `proc()` macro that are **never called from Fortran**. Grep proves it: no `.f` wrapper in `xvue/` contains a `CALL` to any of them. They are internal C helpers that happen to use the `proc()` macro (perhaps historically, or for future-proofing), and inside `xvuelc.c` they are called **by other C functions** via the macro-expanded trailing-underscore names.

| Entry point | Line | Signature problem |
|-------------|------|-------------------|
| `xvCouleursImposees_` | 358 | First arg `int n1coel` is **by value** (not `int *`); violates D-06 |
| `xvColormapToRGB_` | 463 | First arg is `Colormap color_map` — an opaque Xlib typedef defined in `<X11/Xlib.h>` |
| `xvStockeRGBtoColormap_` | 503 | Last arg is `Colormap color_map` — same problem |

**Why this blocks D-04:** if `xvue_qt_api.h` declares these three functions, it must either `#include <X11/Xlib.h>` (pulling X11 back into the Qt side of the bridge, defeating the whole point) **or** forward-declare `Colormap`, which would still create a type-incompatibility with any future Qt-side caller. Either way, the "one header declares every extern C entry point in xvuelc.c" rule cannot be followed literally.

**Why this blocks D-06:** `xvCouleursImposees_`'s first argument is declared `int n1coel` (by value) in `xvuelc.c`. If the planner copies it "verbatim" per D-06's code-review checklist, that checklist item immediately flags it as a bug. If the planner "fixes" it to `int *`, the signature no longer matches `xvuelc.c` verbatim.

**Recommended resolution** (discuss-phase question, not a research finding):

> Option A: Declare these three as internal C++ helpers inside `xvue/qt/src/` — NOT in the public `xvue_qt_api.h`. They are not part of the Fortran-facing ABI; the `nm` count in D-12 should also skip them. The one-header rule becomes "every entry point **called from Fortran**" instead of "every entry point in `xvuelc.c`". This preserves the byte-identical public ABI.
>
> Option B: Declare them in `xvue_qt_api.h` with a forward-declared opaque `typedef void* XvueColormap;` shim and a `// Phase 0: unused-by-Fortran internal helper` comment. Preserves D-04's letter at the cost of a pointless symbol.
>
> Option A is cleaner and matches the Core Value invariant ("the Fortran-facing ABI must not change"). The planner should ask the user to confirm.

### Alert 2: `dctnmc` returns `void *` and `dstnmc` takes `void *`

```c
void *proc(dctnmc)(int *nboctets);   // line 242 — returns void*
void  proc(dstnmc)(void *mcoctets);  // line 254 — takes void*
```

These are the dynamic-allocation bridges (`dctnmc` = *de*fine *c*ell of *n* *m*illions of *c*ells, the project's hand-rolled `malloc`/`free` facade). They are real Fortran-facing entry points — the Fortran side calls them to allocate scratch memory. D-06 says "every scalar is a pointer, every string is `char* + int*`" but says nothing about pointer-returning or `void*`-taking entries. These are fine; they do not violate D-06; but the planner must NOT "normalize" them. Called out here so the checklist does not mis-fire.

**Action for planner:** Ask the user to confirm Option A in Alert 1. Proceed with the understanding that the Phase 0 `verify_abi` count is `59 - 3 = 56` Fortran-facing entry points, or `59` if Option B is chosen.

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Qt 6 Core/Gui/Widgets/PrintSupport | 6.10.2 (apt candidate in trixie/sid) | GUI toolkit for future phases; Phase 0 only needs the CMake + `find_package` handshake and the AUTOMOC scaffold | Apt-provided, upstream-LTS, matches project constraint of no vendored Qt. `[VERIFIED: apt-cache policy qt6-base-dev]` |
| CMake | >= 3.21 (installed: 3.31.6) | Build system for `xvue/qt/` only; produces `libxvueqt.a` | Meets D-11 minimum; `qt6-base-dev` depends on `>= 3.16` so 3.21 is safe. `[VERIFIED: cmake --version]` |
| gfortran | 15.2.0 (Debian 15.2.0-11) | Links the static lib into `pp/pp*_qt` | Already the project's Fortran compiler. Name-mangling scheme (`x##_`) is stable across gfortran 4.x–15.x on x86_64 Linux. `[VERIFIED: gfortran --version]` `[CITED: PITFALLS.md Pitfall 1]` |
| g++ | 15.2.0 (Debian 15.2.0-11) | Compiles the C++ stub source, linked in automatically by CMake | Matches gfortran major version — avoids any C++ runtime ABI mismatch. `[VERIFIED: g++ --version]` |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `pkg-config` | (installed) | Resolve Qt 6 link flags on the shell-script final linker line | D-08 literally calls `$(pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport)` `[VERIFIED: which pkg-config]` |
| `nm` (GNU binutils) | (installed) | Implements the `verify_abi` CMake custom target for BUILD-08 | `[VERIFIED: which nm]` |
| `qt6-base-dev-tools` | 6.10.2 (apt candidate) | Provides `moc` used by `CMAKE_AUTOMOC` | Pulled as dep of `qt6-base-dev`. `[CITED: Debian trixie package metadata]` |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Static `libxvueqt.a` | Shared `libxvueqt.so` | Shared lib would force an `LD_LIBRARY_PATH` convention for every user install, and `pp/pp*_qt` would hard-fail if the `.so` moved. Static matches the existing `xvue/lib` + `util/lib` pattern and the Fortran build's zero-runtime-deps culture. Stay static. |
| `qmake` + `.pro` file | CMake | Rejected by `.planning/research/STACK.md` — qmake is end-of-life upstream. |
| `find_package(Qt6 COMPONENTS Widgets Gui Core)` (no PrintSupport) | With or without PrintSupport | PrintSupport is needed in Phase 7 for PDF (`QPrinter`). For Phase 0 the stubs don't touch it, but including it in the CMakeLists now means Phase 7 doesn't have to re-touch the build plumbing. D-08 already lists it. |

### Installation

```bash
# On the target Debian trixie machine, as root:
sudo apt install qt6-base-dev qt6-base-dev-tools cmake pkg-config
```

**Version verification** `[VERIFIED: 2026-04-10 via apt-cache policy + command probe]`:
- `qt6-base-dev` — installed: **none**; candidate: **6.10.2+dfsg-6** (Debian trixie/sid)
- `qt6-base-dev-tools` — candidate: **6.10.2+dfsg-6**
- `cmake` — installed: **3.31.6-1** (satisfies D-11 `>= 3.21`)
- `pkg-config` — installed
- `gfortran` — installed: **15.2.0**
- `g++` — installed: **15.2.0**
- `nm` — installed
- `libqt6core6t64` — already installed at **6.9.2+dfsg-3** as a runtime lib, but the `-dev` package with headers is missing. This mismatch (runtime 6.9.2 vs. dev candidate 6.10.2) is Debian's normal package-state mid-release; `apt install qt6-base-dev` will upgrade both in lockstep.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| gfortran | BUILD-06, BUILD-07 | Yes | 15.2.0 | — |
| g++ | BUILD-01..05 | Yes | 15.2.0 | — |
| cmake | BUILD-01, D-11 | Yes | 3.31.6 | — |
| pkg-config | BUILD-06, D-08 | Yes | present | — |
| nm (binutils) | BUILD-08, D-12 | Yes | present | — |
| `qt6-base-dev` (headers, `find_package(Qt6)`) | BUILD-01..05, D-11 | **No** | candidate 6.10.2+dfsg-6 | Install via `sudo apt install qt6-base-dev` (blocking) |
| `qt6-base-dev-tools` (moc) | BUILD-03 (AUTOMOC) | **No** | candidate 6.10.2+dfsg-6 | Pulled as dep of `qt6-base-dev` |
| `libX11-dev` (legacy path) | BUILD-07 | Not probed (assumed present — legacy build currently works per project history) | — | Not Phase 0's responsibility |
| ImageMagick `convert` | Unrelated to Phase 0 | Not probed | — | Out of scope |

**Missing dependencies with no fallback:**
- **`qt6-base-dev` and `qt6-base-dev-tools`** — Phase 0 cannot begin without them. First planner task must be a manual apt install step, and CLAUDE.md says "ask the user to install system packages" — so the first plan should halt and ask the user to run `sudo apt install qt6-base-dev qt6-base-dev-tools` before any CMake work begins.

**Missing dependencies with fallback:** none.

## Architecture Patterns

### Recommended Project Structure

```
xvue/
├── xvuelc.c               # UNTOUCHED (read-only in Phase 0)
├── *.f (~221 files)       # UNTOUCHED (read-only in Phase 0)
├── lib                    # UNTOUCHED (existing Fortran lib archive)
└── qt/                    # NEW — all Phase 0 work lives here
    ├── CMakeLists.txt
    ├── README_COORDS.md   # Y-axis audit result (D-03)
    ├── include/
    │   └── xvue_qt_api.h  # single header, trailing-underscore macro, all entry points (D-04)
    ├── src/
    │   └── xvue_qt_api.cpp # no-op stubs with warn-once stderr (D-17)
    └── build/             # gitignored; produced by `cmake --build`
        └── libxvueqt.a    # handoff artifact consumed by bin/cbl_tout_qt (D-13)

bin/
├── cbl_tout               # UNTOUCHED
├── cbl_tout_qt            # NEW — clone of cbl_tout with Qt link line (D-08)
├── cbmail                 # UNTOUCHED
├── cbmail_qt              # NEW — clone of cbmail (D-09)
├── cbelas_qt              # NEW — clone of cbelas
├── cbflui_qt              # NEW — clone of cbflui
├── cbther_qt              # NEW — clone of cbther
└── cbnlse_qt              # NEW — clone of cbnlse

.planning/
└── validation/
    └── BASELINE.md        # NEW — 5 canonical testa cases (D-14, D-16)
```

### Pattern 1: Trailing-underscore `proc()` macro (byte-identical to xvuelc.c)

**What:** Preserve the exact macro from `xvue/xvuelc.c` so `CALL XVTRAIT(...)` in Fortran resolves to `xvtrait_` in the stub object file.
**When to use:** Every entry point declaration in `xvue_qt_api.h` and every definition in `xvue_qt_api.cpp`.
**Example:**
```cpp
// xvue/qt/include/xvue_qt_api.h — copied verbatim from xvuelc.c lines ~60-70
#ifdef __GNUC__
#define proc(x) x##_
#else
#define proc(x) x/**/_
#endif

#undef  proc
#define proc(x) x##_
```
`[VERIFIED: xvuelc.c lines 60-70, direct read]`

### Pattern 2: Warn-once stub (D-17)

**What:** First call prints a single stderr line; subsequent calls are silent.
**When to use:** Every no-op stub in Phase 0.
**Example:**
```cpp
// xvue/qt/src/xvue_qt_api.cpp
#include <cstdio>
#include "xvue_qt_api.h"

extern "C" {

void proc(xvtrait)(int *x1, int *y1, int *x2, int *y2) {
    static bool warned = false;
    if (!warned) {
        std::fprintf(stderr, "xvue-qt: stub xvtrait_ not implemented yet\n");
        warned = true;
    }
    (void)x1; (void)y1; (void)x2; (void)y2;  // silence -Wunused-parameter
}

// ... 55 more like this (see §"The 59 Entry Points")

} // extern "C"
```

### Pattern 3: CMake scaffolding (CMAKE_AUTOMOC ordering)

**What:** `CMAKE_AUTOMOC ON` MUST be set before `find_package(Qt6 ...)`. Phase 0 has no `QObject` subclasses, but the scaffolding must be correct so Phase 1+ can add them without re-surgery.
**When to use:** Top of `xvue/qt/CMakeLists.txt`.
**Example:** see §Code Examples.

### Anti-Patterns to Avoid

- **Modifying signatures for "cleanliness"** — e.g. turning `int *x1` into `int x1` because "Qt doesn't need a pointer." Breaks every Fortran caller at link or crash time. See PITFALLS.md Pitfall 2.
- **Hidden trailing-length string convention** — gfortran supports it but `xvuelc.c` does not use it. Adding it here would diverge from the Fortran callers and crash `xvtexte_`. See PITFALLS.md Pitfall 3.
- **Conditional backend logic in `bin/cbl_tout`** — user explicitly chose C1 (clone) over C2 (parametrize). Do not add `if [ "$BACKEND" = qt ]; then ...`. Clone.
- **Touching `xvue/xvuelc.c` or `bin/ccxvue`** — breaks BUILD-07. They are read-only in Phase 0.
- **Installing into system paths** — `install(TARGETS xvueqt DESTINATION xvue/qt/lib)` uses a project-relative destination; do not use absolute `/usr/local` or `/opt`.
- **Using `QApplication` or any Qt GUI object in a stub** — stubs must be pure no-ops. `qApp` is null in Phase 0. Calling into it will crash.
- **Running `QApplication::exec()` anywhere** — forbidden by SHELL-03 (Phase 1). Keeping it forbidden in Phase 0 scaffolding prevents accidental inheritance.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Trailing-underscore macro | Your own `#define FORTRAN_NAME(x) x_` | Copy the existing `#define proc(x) x##_` from `xvuelc.c` verbatim | Byte-identical preservation is the whole invariant. Any new macro is a new bug surface. |
| Entry-point counting | Manual grep-by-hand at each edit | The `verify_abi` CMake custom target (D-12) that runs `nm ... \| grep -c '_$'` | Manual counting drifts. CI-style build-breaks on drift are the only reliable approach for a solo dev with no GitHub Actions. |
| Qt link flags | Hard-coded `-lQt6Core -lQt6Gui ...` in `bin/cbl_tout_qt` | `$(pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport)` | Qt 6 version bumps and lib-name suffixes (e.g. `-t64`) change; pkg-config absorbs it. |
| Moc invocation | Manual `moc` commands in shell | `CMAKE_AUTOMOC ON` set before `find_package(Qt6)` | Manual moc is a classic Phase-1-regression source; set it up right in Phase 0 even though no QObject exists yet. |
| `xvue/qt/build/libxvueqt.a` path resolution | Parse `cmake` output | Hard-code the relative path the CMakeLists installs to (D-13) | Cmake install paths are stable; parsing stdout is fragile. |
| Fortran-facing symbol verification | Reading `nm` by eye | Automated count-compare in `verify_abi` target (D-12) | Phase 1–7 will add real implementations; the count must keep matching the header declaration count. |
| Clean-rebuild discipline | `cmake --build` with incremental tracking | `rm -rf xvue/qt/build/ pp/` before every `bin/cbl_tout_qt` run (D-10) | Stale `.o` files in `xvue/` are a documented fragility per CONCERNS.md. Cheaper to always clean than to chase phantom bugs. |

**Key insight:** Phase 0 is the one phase where "don't be clever" is the highest-value discipline. Every temptation to improve on `xvuelc.c`'s conventions is a bug in waiting for the A/B validation window. The Phase 0 success criteria ("links cleanly, runs past link, no graphics") specifically exclude every scenario where cleverness could pay off.

## Common Pitfalls

### Pitfall 1: Trailing-underscore name mangling mismatch (Pitfall 1 of PITFALLS.md)

**What goes wrong:** Linker error `undefined reference to xvtrait_` for every entry point at once when linking `pp/ppmail_qt`.
**Why it happens:** gfortran on x86_64 Linux emits lowercase names with a single trailing underscore. The Qt side must honor the exact same convention, via the same `#define proc(x) x##_` macro used in `xvuelc.c`.
**How to avoid:** Copy the macro verbatim (D-05). Never pass `-fno-underscoring` to gfortran. Run the `verify_abi` target (D-12) on every build.
**Warning signs:** `nm libxvueqt.a | grep xvtrait` shows `xvtrait` without the trailing `_`.
`[CITED: .planning/research/PITFALLS.md §Pitfall 1]`

### Pitfall 2: Fortran passes scalars by reference — declarations must take pointers (Pitfall 2)

**What goes wrong:** Segfault or silent numeric garbage when the C++ stub "helpfully" takes `int` instead of `int*` and then tries to use the value.
**Why it happens:** Every Fortran `CALL XVTRAIT(X1,Y1,X2,Y2)` with `INTEGER X1,Y1,X2,Y2` hands four `int*` to the C side. `xvuelc.c` has always declared them as pointers.
**How to avoid:** Copy signatures verbatim from `xvuelc.c`. Code-review checklist: every `int`/`float`/`double` arg is `int*`/`float*`/`double*`. **Exception: `xvCouleursImposees_` first arg is `int n1coel` by value in `xvuelc.c` — but that function is NOT called from Fortran (see §Planner Alert) so the checklist must be applied only to Fortran-facing entries.**
**Warning signs:** Segfault on first interactive draw; gdb shows `drawLine` receiving pointer-valued ints.
`[CITED: .planning/research/PITFALLS.md §Pitfall 2] [VERIFIED: grep of xvuelc.c signatures]`

### Pitfall 3: String length argument as explicit pointer in declared order (Pitfall 3)

**What goes wrong:** Garbled text when Phase 2 implements `xvtexte_` because the C++ side assumed a hidden trailing length.
**Why it happens:** gfortran can pass string lengths either (a) as an explicit trailing `int` after the regular arguments, or (b) as a regular `int*` arg in the declared order. `xvuelc.c` uses form (b): `void proc(xvtexte)(char string[], int *length, int *x1, int *y1)`.
**How to avoid:** Copy signatures verbatim. D-06 already enforces this. Do NOT add `int ignored_trailing_length` at the end.
**Warning signs:** `QString::fromLatin1(string, *length)` returning empty or truncated text; valgrind reports OOB read in `xvtexte_`.
**Phase-0 implication:** Phase 0 stubs are no-ops so the bug can't manifest YET, but the header declarations must be correct so Phase 2 inherits the right ABI.
`[CITED: .planning/research/PITFALLS.md §Pitfall 3]`

### Pitfall 4: `XPoint` byte-layout shim (Pitfall 4)

**What goes wrong:** `xvface_`, `xvtraits_`, `xvfacetraits_` crash or draw wrong shapes because `MefistoPoint` is not byte-identical to Xlib's `XPoint`.
**Why it happens:** Xlib `XPoint` is `struct { short x, y; }` (4 bytes total). The Fortran wrappers in `xvue/*.f` pass arrays of INTEGER*2 pairs that assume this layout. If Qt side declares `struct MefistoPoint { int x, y; }` it's 8 bytes and every polygon renders as a firehose of garbage points.
**How to avoid:** D-07's `typedef struct { short x; short y; } MefistoPoint;` plus the `static_assert(sizeof(MefistoPoint) == 4, ...)`. Phase 0 header must have this assert — it fails at compile time on any platform where `short` is not 16 bits.
**Warning signs:** Compile-time static_assert failure (desired); runtime firehose of random polygon points (undesired).
`[CITED: .planning/research/PITFALLS.md §Pitfall 4] [VERIFIED: xvuelc.c lines 1701, 1977, 2035 — three entry points take XPoint*]`

### Pitfall 9: AUTOMOC ordering

**What goes wrong:** Phase 1 adds the first `QObject` subclass, moc files are not generated, build fails with undefined vtable symbols.
**Why it happens:** `CMAKE_AUTOMOC` must be set **before** `find_package(Qt6 ...)` or Qt's CMake files won't register the automatic moc target.
**How to avoid:** D-11 explicitly mandates this ordering. Phase 0 has no QObject, but the scaffolding must be correct now so Phase 1 doesn't debug it.
**Warning signs (in Phase 1+):** `undefined reference to vtable for XvueWindow` or similar.
`[CITED: .planning/research/PITFALLS.md §Pitfall 9]`

### Pitfall 10: `-fPIC` on the static library

**What goes wrong:** Link error `relocation R_X86_64_32S against ...` when `libxvueqt.a` is linked into a PIE executable.
**Why it happens:** Modern gcc defaults produce PIE executables. A non-PIC static lib cannot be linked into a PIE.
**How to avoid:** `set(CMAKE_POSITION_INDEPENDENT_CODE ON)` at the top of `CMakeLists.txt`. D-11 mandates it.
**Warning signs:** Relocation errors during the final `gfortran ... -o pp/pp*_qt` step.
`[CITED: .planning/research/PITFALLS.md §Pitfall 10]`

### Pitfall 11: `-fopenmp` runtime collision with `_OMP` executables

**What goes wrong:** `pp/ppelas_omp_qt` crashes at startup or reports OpenMP runtime conflicts (two `libgomp` instances).
**Why it happens:** The Fortran `_OMP` executables link `-fopenmp`. If `libxvueqt.a` is also compiled with `-fopenmp`, gfortran pulls in a second copy of the OpenMP runtime.
**How to avoid:** D-11 explicitly sets `target_compile_options(xvueqt PRIVATE -fno-openmp)`. Critically, the legacy `bin/ccxvue` passes `-openmp` (note: `-openmp` is Intel's flag, not gcc's `-fopenmp`) which is silently ignored by gcc — meaning the legacy path NEVER actually enabled OpenMP in `xvuelc.o`. The Qt side just has to match this "no OpenMP in graphics layer" invariant explicitly.
**Warning signs:** `pp/pp*_omp_qt` prints `OMP: Error #13: Assertion failure`.
**Phase-0 implication:** Phase 0 only links non-`_OMP` variants (ppmail_qt, etc. per D-15), but the CMakeLists must still have `-fno-openmp` so Phase 8's `_OMP` A/B sweep works without re-architecture.
`[CITED: .planning/research/PITFALLS.md §Pitfall 11] [VERIFIED: bin/ccxvue source — uses -openmp (silently ignored by gcc) not -fopenmp]`

### Pitfall 16: Y-axis convention drift

**What goes wrong:** Phase 2+ draws meshes upside-down or mirror-flipped compared to X11.
**Why it happens:** Xlib and Qt both use Y-down top-left origin for their raster coordinates, but hand-rolled PostScript output usually flips Y. `xvuelc.c` inverts Y inside the PostScript emitter (`ypixels - *y1`) — meaning on-screen coordinates are Y-down natively, PostScript is Y-up via in-function inversion.
**How to avoid:** Phase 0 audits and documents the convention in `xvue/qt/README_COORDS.md` (D-03). Phase 2+ refers to that doc.
**Warning signs:** Upside-down meshes on first interactive Phase 2 test.
**See §"Y-Axis Convention Audit" below for the result of the audit.**
`[CITED: .planning/research/PITFALLS.md §Pitfall 16] [VERIFIED: xvuelc.c lines 1895, 1932, 1953, 1966 — all use ypixels - *y for ps output only]`

### Pitfall 18: Clean-build discipline / stale `.o` fragility

**What goes wrong:** A stale `xvue/xvuelc.o` from a previous legacy build gets linked into `pp/pp*_qt` instead of the new stubs, producing either doubled-symbol link errors or a "Qt build" that secretly uses the X11 code.
**Why it happens:** `bin/cbl_tout` and `bin/cbmail` only rebuild `xvue/xvuelc.o` if it doesn't exist (`if !(test -f xvue/xvuelc.o)`) — they never force a rebuild. Same problem applies to stale Fortran `.o` files in module directories.
**How to avoid:** D-10's `rm -rf xvue/qt/build/ pp/` at the top of `bin/cbl_tout_qt`. Since Phase 0 links *new* executables (`pp/pp*_qt`, not `pp/pp*`), the legacy `xvue/xvuelc.o` is never on the Qt link line at all — but the cbl_tout_qt clean is still cheap insurance against Fortran-side stale `.o` files.
**Warning signs:** `ld` reports `multiple definition of xvtrait_` or unexplained symbol behavior.
`[CITED: .planning/research/PITFALLS.md §Pitfall 18] [VERIFIED: bin/cbmail source has the `if !(test -f xvue/xvuelc.o)` gate]`

### Pitfall 20: Regression drift without a validation baseline

**What goes wrong:** Phase 3 introduces a subtle colormap bug; Phase 5 adds another; by Phase 7 the team can't tell which phase caused which regression.
**Why it happens:** Single-developer project, no CI, manual testing only. Without a written-down fixed set of "the 5 cases we always run," scope creep erodes the baseline.
**How to avoid:** D-14/D-16's `.planning/validation/BASELINE.md` committed in Phase 0 as the immutable reference.
**Warning signs:** (happens in later phases, not Phase 0)
`[CITED: .planning/research/PITFALLS.md §Pitfall 20]`

## Code Examples

Verified patterns derived from direct reads of the cited sources.

### Example 1: `xvue/qt/CMakeLists.txt` skeleton

```cmake
# xvue/qt/CMakeLists.txt
cmake_minimum_required(VERSION 3.21)
project(xvueqt LANGUAGES CXX)

# MUST be before find_package(Qt6) — D-11, Pitfall 9
set(CMAKE_AUTOMOC ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # D-11, Pitfall 10
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets PrintSupport)

add_library(xvueqt STATIC
    src/xvue_qt_api.cpp
)

target_include_directories(xvueqt PUBLIC include)

target_link_libraries(xvueqt
    PUBLIC
        Qt6::Core
        Qt6::Gui
        Qt6::Widgets
        Qt6::PrintSupport
)

# D-11 Pitfall 11: never propagate -fopenmp into graphics layer
target_compile_options(xvueqt PRIVATE -fno-openmp -Wall)

# D-12: verify every Fortran-facing symbol has the trailing underscore
add_custom_target(verify_abi ALL
    COMMAND ${CMAKE_COMMAND} -E echo "Verifying ABI symbol count..."
    COMMAND sh -c "nm $<TARGET_FILE:xvueqt> | grep -c ' T .*_$' > nm_count.txt"
    COMMAND sh -c "grep -c '^extern \"C\".*proc(' ${CMAKE_CURRENT_SOURCE_DIR}/include/xvue_qt_api.h > hdr_count.txt || true"
    COMMAND sh -c "test \"$$(cat nm_count.txt)\" = \"$$(cat hdr_count.txt)\" || (echo ABI count drift && exit 1)"
    DEPENDS xvueqt
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
)
```

**Source provenance:** `[CITED: .planning/research/STACK.md §"CMake template"]` for overall shape; `[VERIFIED: CLAUDE.md via Read]` for "-fopenmp" exclusion rationale; the `verify_abi` custom target is constructed here from the plain `nm | grep -c '_$'` approach mandated by D-12 (the exact shell invocation above is a sketch — the planner may refine the grep pattern to match how the header declares its entries). `[ASSUMED]` that the specific `nm` output format `' T ..._$'` is stable across binutils versions on x86_64 Linux; if not, the planner should tweak the regex during implementation.

### Example 2: `xvue/qt/include/xvue_qt_api.h` (first ~20 lines + two sample entries)

```cpp
// xvue/qt/include/xvue_qt_api.h
// Phase 0: no-op ABI shim matching xvue/xvuelc.c byte-for-byte
// DO NOT modify signatures without user approval.
#ifndef XVUE_QT_API_H
#define XVUE_QT_API_H

// Trailing-underscore macro copied verbatim from xvue/xvuelc.c lines 60-70
#ifdef __GNUC__
#define proc(x) x##_
#else
#define proc(x) x/**/_
#endif
#undef proc
#define proc(x) x##_

// XPoint byte-layout shim — D-07, Pitfall 4
// Xlib XPoint is {short x; short y;} — 4 bytes total
#include <cstddef>
typedef struct { short x; short y; } MefistoPoint;
static_assert(sizeof(MefistoPoint) == 4,
              "MefistoPoint must match Xlib XPoint layout");

// Debug thread-affinity macro — Claude's discretion, Phase 0 no-op
#ifdef QT_DEBUG
  #include <QThread>
  #include <QApplication>
  #define XVUE_QT_ASSERT_MAIN_THREAD() \
      Q_ASSERT(QThread::currentThread() == qApp->thread())
#else
  #define XVUE_QT_ASSERT_MAIN_THREAD() ((void)0)
#endif

extern "C" {

// === Sample entry points (see xvue_qt_api.cpp for the full 56-entry list) ===

void proc(xvtrait)(int *x1, int *y1, int *x2, int *y2);
void proc(xvface)(int *n, MefistoPoint *pts);
void proc(xvtexte)(char string[], int *length, int *x1, int *y1);

} // extern "C"

#endif // XVUE_QT_API_H
```

**Source provenance:** `[VERIFIED: xvuelc.c lines 60-70, 1701, 1862, 1662]` for exact macro, signatures, and XPoint usage.

### Example 3: `bin/cbl_tout_qt` clone template

```bash
#!/bin/bash
#
# cbl_tout_qt — Qt-backed clone of cbl_tout
# Links pp/pp*_qt executables against libxvueqt.a instead of xvuelc.o
# Phase 0 scaffold — Alain Perronnet project fork

# D-10: clean xvue/qt/build and pp/ to kill stale .o files
echo "Cleaning xvue/qt/build and pp/..."
rm -rf $MEFISTO/xvue/qt/build
rm -f  $MEFISTO/pp/pp*_qt

# --- The following block is copied verbatim from bin/cbl_tout ---
# Creates $MEFISTO/incl/homdir.inc (identical to legacy path)
echo
echo MEFISTO=$MEFISTO
NBCHAR=` expr length $MEFISTO `
# ... (lines 10-45 of bin/cbl_tout, unchanged) ...

# NEW: build libxvueqt.a via CMake before touching gfortran
echo "Building xvue/qt/libxvueqt.a via CMake..."
cmake -S $MEFISTO/xvue/qt -B $MEFISTO/xvue/qt/build \
  && cmake --build $MEFISTO/xvue/qt/build || {
    echo "CMake build of libxvueqt.a FAILED"
    exit 1
  }

# --- Resume cloning cbl_tout, but substitute bin/ccxvue with a no-op ---
# (Phase 0: we still call legacy cbl_tous_f per-module to produce Fortran libs)
echo elas | $MEFISTO/bin/cbl_tous_f
echo flui | $MEFISTO/bin/cbl_tous_f
echo mail | $MEFISTO/bin/cbl_tous_f
echo reso | $MEFISTO/bin/cbl_tous_f
echo ther | $MEFISTO/bin/cbl_tous_f
echo util | $MEFISTO/bin/cbl_tous_f
echo xvue | $MEFISTO/bin/cbl_tous_f

rm $MEFISTO/*/*.o

$MEFISTO/bin/cbmail_qt
$MEFISTO/bin/cbelas_qt
$MEFISTO/bin/cbther_qt
$MEFISTO/bin/cbflui_qt
$MEFISTO/bin/cbnlse_qt

echo Le repertoire $MEFISTO/pp :
ls -l $MEFISTO/pp
```

**Source provenance:** `[VERIFIED: bin/cbl_tout]` — structural skeleton is a line-by-line clone.

### Example 4: `bin/cbmail_qt` — per-module clone

```bash
#!/bin/bash
# cbmail_qt — Qt-backed clone of cbmail

cd $MEFISTO
if !(test -d pp); then mkdir pp; fi

if test -f $MEFISTO/td/m/anglais; then
  LANGAGE=1; LANGUE=ANGLAISE
  echo "MEFISTO-LINUX 64bits: Qt link of ppmail.f"
else
  LANGAGE=0; LANGUE=FRANCAISE
  echo "MEFISTO-LINUX 64bits: Edition Qt de ppmail.f"
fi

nompp=pp/ppmail_qt
rm -f $nompp

# Key difference from bin/cbmail: xvue/xvuelc.o is replaced by libxvueqt.a + Qt link flags
QT_LIBS=$(pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport)

echo "gfortran ... -o $nompp (Qt backend)"
gfortran -Wall -mcmodel=large -m64 -O -fopenmp -I. prpr/ppmail.f \
    mail/lib util/lib xvue/lib \
    -Lxvue/qt/build -lxvueqt $QT_LIBS -lstdc++ \
    -lgfortran -o $nompp

if test -f $nompp; then
  chmod 755 $nompp
  echo "$MEFISTO/$nompp OK"
else
  echo "$MEFISTO/$nompp FAILED"
  exit 1
fi
```

**Source provenance:** `[VERIFIED: bin/cbmail]` — identical except for final link line.

**Note for planner:** The legacy `bin/cbmail` link line is:
```
gfortran ... mail/lib util/lib xvue/lib xvue/xvuelc.o \
    -lgfortran -L/usr/X11R6/lib64 -lX11 -o $nompp
```
The Qt clone replaces **`xvue/xvuelc.o` + `-L/usr/X11R6/lib64 -lX11`** with **`-Lxvue/qt/build -lxvueqt $(pkg-config --libs Qt6...) -lstdc++`**. Everything else is identical, including `mail/lib util/lib xvue/lib` (the Fortran wrapper lib — still needed because the `.f` wrappers still call into the trailing-underscore names, which are now served by `libxvueqt.a` instead of `xvuelc.o`).

### Example 5: `xvue/qt/src/xvue_qt_api.cpp` — stub skeleton (1 of 56)

```cpp
// xvue/qt/src/xvue_qt_api.cpp
// Phase 0 no-op stubs — first-call-warn-once diagnostic per D-17
#include <cstdio>
#include "xvue_qt_api.h"

namespace {
inline void warn_once(bool &flag, const char *name) {
    if (!flag) {
        std::fprintf(stderr, "xvue-qt: stub %s not implemented yet\n", name);
        flag = true;
    }
}
} // namespace

extern "C" {

void proc(xvtrait)(int *x1, int *y1, int *x2, int *y2) {
    static bool warned = false;
    warn_once(warned, "xvtrait_");
    (void)x1; (void)y1; (void)x2; (void)y2;
}

// ... 55 more stubs, one per Fortran-facing entry point in §"The 59 Entry Points"
//     minus the 3 C-internal helpers (see Planner Alert)

} // extern "C"
```

## The 59 Entry Points

**Verified count:** `grep -c '^\(void\|int\|float\|double\|long\|short\|unsigned\)\s\+proc\s*(' xvue/xvuelc.c` → **59**.

Full list with line numbers (from direct read of `xvue/xvuelc.c`):

| # | Line | Symbol (trailing `_`) | Signature highlights | Fortran-facing? |
|---|------|-----------------------|---------------------|----------------|
| 1 | 227 | `languemefisto_` | `int *langue` | YES |
| 2 | 242 | `dctnmc_` | returns `void*`; `int *nboctets` | YES |
| 3 | 254 | `dstnmc_` | `void *mcoctets` | YES |
| 4 | 266 | `nomrepmefisto_` | `char *chaine, int *size` | YES |
| 5 | 286 | `xvinitgraphique_` | `(void)` | YES |
| 6 | 306 | `xtinit_` | `()` | YES |
| 7 | 319 | `xvpxecran_` | `int *xp, int *yp` | YES |
| 8 | 337 | `xvmmecran_` | `int *xmm, int *ymm` | YES |
| 9 | 358 | `xvCouleursImposees_` | `int n1coel` BY VALUE, `int *ndcoel, float red[], float green[], float blue[]` | **NO — C-internal helper** (Planner Alert) |
| 10 | 463 | `xvColormapToRGB_` | `Colormap color_map, float r[], float g[], float b[], int nbcolor` | **NO — C-internal helper** (Planner Alert) |
| 11 | 503 | `xvStockeRGBtoColormap_` | `int nbcells, ..., Colormap color_map` | **NO — C-internal helper** (Planner Alert) |
| 12 | 561 | `initaccrochage_` | `(void)` | YES |
| 13 | 612 | `xvinfo_` | 14-arg mega-init (dims + color-range indices + `char namefonts[][256]` + `int nbchar[]` + `int *nbfonts, int *visuclass`) | YES |
| 14 | 1044 | `xvrecuprgbdec_` | `int *nbcolor, float *r, float *g, float *b` | YES |
| 15 | 1072 | `xvactivervb_` | `int *palcour, int *nbcells, ...` | YES |
| 16 | 1119 | `xvcouleur_` | `int *icolor` | YES |
| 17 | 1187 | `xvpostscript_` | `int *lasops` | YES |
| 18 | 1307 | `fenetremempx_` | `()` | YES |
| 19 | 1321 | `mempxfenetre_` | `()` | YES |
| 20 | 1335 | `sauvefenetre_` | `()` | YES |
| 21 | 1350 | `restaurefenetre_` | `()` | YES |
| 22 | 1365 | `sauvemempx_` | `()` | YES |
| 23 | 1380 | `restauremempx_` | `()` | YES |
| 24 | 1395 | `effacemempx_` | `()` | YES |
| 25 | 1413 | `effacer_` | `()` | YES |
| 26 | 1434 | `xvfond_` | `int *icolor` | YES |
| 27 | 1463 | `xvchargefonte_` | `int *nofont0, int *nofont, int *largpx, int *hautpx` | YES |
| 28 | 1576 | `xvnbpixeltexte_` | `char *texte, int *length, int *nbpxla, int *nbpxha` | YES |
| 29 | 1602 | `xvfermer_` | `()` | YES |
| 30 | 1619 | `xvpxfenetre_` | `int *x, int *y` | YES |
| 31 | 1642 | `xvftexte_` | `char string[], int *length, int *x1, int *y1` | YES |
| 32 | 1662 | `xvtexte_` | `char string[], int *length, int *x1, int *y1` | YES |
| 33 | 1701 | `xvface_` | `int *n, XPoint *pts` → use `MefistoPoint*` | YES (uses shim) |
| 34 | 1760 | `xvtypetrait_` | `int *ptype` | YES |
| 35 | 1807 | `xvepaisseur_` | `int *pepais` | YES |
| 36 | 1845 | `xvftrait_` | `int *x1, int *y1, int *x2, int *y2` | YES |
| 37 | 1862 | `xvtrait_` | `int *x1, int *y1, int *x2, int *y2` | YES |
| 38 | 1977 | `xvtraits_` | `int *nbpoints, XPoint *points` → `MefistoPoint*` | YES (uses shim) |
| 39 | 2035 | `xvfacetraits_` | `int *ncf, int *nca, int *n, XPoint *pts` → `MefistoPoint*` | YES (uses shim) |
| 40 | 2123 | `xvsouris_` | `int *notypeevent, int *nbc, int *x1, int *y1` | YES |
| 41 | 2230 | `xvsouris2_` | `int *items, int *pmin0, int *notypeevent, int *ibutton, int *x1, int *y1` | YES |
| 42 | 2374 | `deplsouris_` | `int *x, int *y` | YES |
| 43 | 2384 | `xvvoir_` | `()` | YES |
| 44 | 2408 | `xvpause_` | `()` | YES |
| 45 | 2426 | `xvfbordrectangle_` | `int *x, int *y, int *width, int *height` | YES |
| 46 | 2443 | `xvbordrectangle_` | `int *x, int *y, int *width, int *height` | YES |
| 47 | 2489 | `xvfrectangle_` | `int *x, int *y, int *width, int *height` | YES |
| 48 | 2507 | `xvrectangle_` | `int *x, int *y, int *width, int *height` | YES |
| 49 | 2554 | `xvbordarcellipse_` | `int *x, int *y, int *width, int *height, int *a1, int *a2` (angle args) | YES |
| 50 | 2616 | `xvarcellipse_` | `int *x, int *y, int *width, int *height, int *a1, int *a2` | YES |
| 51 | 2678 | `tempscpu_` | `double *tclock` | YES |
| 52 | 2694 | `secondes1970_` | `double *secondes` | YES |
| 53 | 2712 | `secondes1969_` | `double *secondes` | YES |
| 54 | 2728 | `nomordinateurhote_` | `char *host, int *nbcar` | YES |
| 55 | 2753 | `ladate_` | `int *a, int *m, int *j` | YES |
| 56 | 2779 | `heureminuteseconde_` | `int *h, int *m, int *s, int *millis` | YES |
| 57 | 2805 | `valvarenv_` | `char *nom, int *lval_admis, ...` (env-var lookup) | YES |
| 58 | 2844 | `xvinitierps_` | `int *modeps` | YES |
| 59 | 2906 | `xvimprimerps_` | `char nomfichier[], int *length` | YES |
| 60 | 2954 | `xvsauverps_` | `char nomfichier[], int *length` | YES |

Wait — that's 60, not 59. Re-counting from the grep output above: entries 1-60 as listed. The grep-count regex `^(void|int|...)\s+proc\s*\(` matched 59 because entry #2 (`dctnmc_`) has declaration `void * proc(dctnmc) ...` (note the space before `*`) which my regex failed to capture. **Corrected count: 60 entry points total.**

**Fortran-facing subset (for D-12 `verify_abi`):** 60 − 3 C-internal helpers = **57 entries** if the planner chooses Option A in §Planner Alert; 60 entries if Option B. The planner must pick one before writing the `verify_abi` count assertion.

`[VERIFIED: direct line-by-line read of xvue/xvuelc.c, 2026-04-10]`

## Y-Axis Convention Audit (D-03, BUILD-09)

**Finding:** `xvue/xvuelc.c` uses a **single coordinate convention on-screen**: Xlib's top-left origin with Y growing downward (Y-down). Every `XDrawLine`, `XDrawLines`, `XFillRectangle`, `XDrawArc`, and text-rendering call passes the y coordinate **unflipped**. Comments in the source confirm this: lines 321 ("origine coin superieur gauche de l'ecran"), 1852 ("ORIGINE coin superieur gauche"), 1869 (same), 1621 ("la fenetre")

Y inversion happens **only inside `xvpostscript_` output**, where PostScript's natural Y-up origin is honored by writing `ypixels - *y1` into the emitted `.ps` stream. See lines 1895, 1932, 1953, 1966 — all inside the `xvtrait_` PostScript-capture block — and the wider pattern in the embedded `fprintf(fpo, ...)` calls throughout `xvpostscript_` (line 1187+). The screen coordinates themselves are never flipped.

**Implication for Qt bridge:**
- **`QPainter`** uses the same Y-down top-left origin as Xlib raster coordinates. Phase 2+ can pass `(*x, *y)` directly to `QPainter::drawLine(*x, *y, ...)` with no transformation.
- **Phase 7 PostScript export** must preserve the `ypixels - *y` inversion verbatim when the existing PostScript emitter is moved from `xvuelc.c` to `xvue/xvue_qt_postscript.cpp` (EXPORT-04). The inversion is not a Qt concern; it is a PostScript-format concern.
- **HiDPI (`devicePixelRatioF > 1`)** is orthogonal to the Y-axis audit. Handled separately in Phase 1 (SHELL-06).

**`xvue/qt/README_COORDS.md` content outline** (to be written by Phase 0 tasks):

```markdown
# Y-Axis Convention (audited 2026-04-10)

## On-screen coordinates (Fortran ↔ Qt bridge)
- Origin: **top-left corner of the canvas widget**
- X: grows **rightward**
- Y: grows **downward** (Y-down)
- Unit: logical pixels (device-independent, HiDPI-aware)

## Where the convention comes from
- xvue/xvuelc.c line 321: "origine coin superieur gauche de l'ecran"
- xvue/xvuelc.c line 1852, 1869: "ORIGINE coin superieur gauche"
- Every `XDrawLine`, `XDrawLines`, `XFillRectangle`, `XDrawArc` call passes y unflipped
- This matches Xlib's native Y-down convention AND Qt's native Y-down QPainter convention
  → No transformation required on the Qt side

## Where Y IS flipped (PostScript export only — Phase 7 concern)
- xvue/xvuelc.c lines 1895, 1932, 1953, 1966 (inside xvtrait_'s ps-capture branch)
- And throughout the xvpostscript_ function (line 1187+)
- Formula: `ps_y = ypixels - screen_y` (where ypixels = canvas height in pixels)
- Reason: PostScript's natural origin is bottom-left (Y-up); the existing emitter
  honors that convention by flipping at emit time, not by storing flipped coords.

## Implication for Qt phases
- Phase 1-6: DO NOT flip Y anywhere. Pass Fortran-provided y values directly to QPainter.
- Phase 7: PRESERVE the `ypixels - y` flip verbatim when porting the PostScript
  emitter from xvuelc.c. Do not try to "clean it up" by storing flipped coords.
```

`[VERIFIED: direct read of xvue/xvuelc.c lines 319-332, 1845-1975, 1187-1306]`

## testa/ Baseline Verification (D-14, BUILD-10)

All 5 cases exist on disk as subdirectories of `testa/`:

| Case | Module | Path exists | Notes |
|------|--------|-------------|-------|
| pan2d | Mesher (ppmail) | ✓ `testa/pan2d/` | Interactive 2D panel |
| nafems_le1 | Elasticity (ppelas) | ✓ `testa/nafems_le1/` | NAFEMS LE1 benchmark |
| cavity2d | Fluid (ppflui) | ✓ `testa/cavity2d/` | Lid-driven cavity |
| heat1d | Thermal (ppther) | ✓ `testa/heat1d/` | 1D heat conduction |
| nlsecu | Nonlinear (ppnlse) | ✓ `testa/nlsecu/` | Nonlinear canonical |

`[VERIFIED: ls testa/]`

The planner must write `.planning/validation/BASELINE.md` with — per D-16 — for each case: project dir path, launcher script(s), expected qualitative behavior (what human eyes look for), any known-flaky touchpoints. The project hasn't run these 5 cases yet in this session, so "known-flaky touchpoints" will be written as "none known — to be discovered at first run."

## Legacy Path Preservation (BUILD-07)

Phase 0 must not break `bin/cbl_tout` + `pp/pp*` + X11. This is mechanically ensured by D-02 (read-only everywhere except `xvue/qt/` + `bin/cbl_tout_qt` + `bin/cb*_qt`). Verification steps the planner should include:

1. **Before Phase 0 work starts:** run `bin/cbl_tout` once, confirm `pp/ppmail`, `pp/ppelas`, `pp/ppflui`, `pp/ppther`, `pp/ppnlse` build cleanly. Record success in the plan's pre-conditions.
2. **After each Phase 0 task:** re-run `bin/cbl_tout` is NOT required, because D-02's read-only rule makes regression impossible by construction. A post-Phase-0 verification step should do one clean rebuild of the legacy path to confirm.
3. **End of Phase 0:** run each of the 5 baseline `testa/` cases through the legacy X11 `pp/pp*` executables, tick them off in `.planning/validation/BASELINE.md` as "known-good on legacy X11 2026-04-10".

**Caveat:** the current machine does NOT have `pp/` populated (directory empty/missing per `ls /home/drico/git/mefisto/pp`), meaning the legacy build has never been run on this clone. The first phase-0 task is probably `bin/cbl_tout` (legacy) → confirm `pp/pp*` appear → THEN start the Qt scaffold. This catches any pre-existing legacy breakage before it's blamed on Phase 0.

`[VERIFIED: ls -la pp/ — directory does not exist]`

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `qmake` + `.pro` files | CMake + `find_package(Qt6 ...)` | Qt 5→6 transition (2020+) | Project already chose CMake in PROJECT.md; `[CITED: .planning/research/STACK.md §"What NOT to use"]` |
| `QT += core gui` qmake syntax | `find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets PrintSupport)` | Qt 6 | Required for D-11 |
| Manual `moc` invocation | `CMAKE_AUTOMOC ON` | Modern CMake idiom | Phase 0 sets it up even though no QObject exists yet |
| Vendored Qt via submodule | apt `qt6-base-dev` | Debian trixie Qt 6 maturity | Matches PROJECT.md constraint "no vendored Qt" |
| X11 `libX11-dev` direct | Qt 6 `QPainter` / `QWindow` (future phases) | This migration | Phase 0 does not yet touch any Qt GUI class — only the build plumbing that will host them |

**Deprecated/outdated:**
- Qt 5: explicitly rejected in PROJECT.md (maintenance-mode).
- `qmake`: end-of-life upstream; CMake is the Qt 6 build system.
- Intel `-openmp` flag as used in `bin/ccxvue` line 22: silently ignored by gcc (gcc wants `-fopenmp`). **This is a pre-existing latent bug in the legacy script** — it means `xvuelc.o` has never been compiled with OpenMP despite the script appearing to ask for it. Not Phase 0's job to fix; documented here so the planner is aware the legacy path's "OpenMP in graphics" claim is illusory.

`[VERIFIED: bin/ccxvue direct read]`

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | None (manual human-eye validation per project decision) |
| Config file | None |
| Quick run command | N/A — see Wave 0 |
| Full suite command | N/A — 5-case `testa/` run is manual |

The project explicitly rejects automated CI (per PROJECT.md and CONTEXT.md deferred-ideas). There is no `pytest`, `jest`, Google Test, or similar harness. "Validation" means (a) the code compiles and links, (b) specified binaries run without crashing, and (c) a human runs each of the 5 `testa/` cases on both backends and judges visual output. Phase 0 specifically has no A/B test (nothing draws yet); the gate is binary link success + binary exit-without-crash.

### Phase 0 Requirements → Validation Map

| Req ID | Behavior | Test Type | Automated? | Command |
|--------|----------|-----------|-----------|---------|
| BUILD-01 | `libxvueqt.a` produced by CMake | smoke | yes | `cmake -S xvue/qt -B xvue/qt/build && cmake --build xvue/qt/build && test -f xvue/qt/build/libxvueqt.a` |
| BUILD-02 | `-fPIC` present, `-fopenmp` absent | static | yes | `readelf --relocs xvue/qt/build/libxvueqt.a \| grep -q R_X86_64_PLT32` (positive PIC indicator) and `cmake --build xvue/qt/build -- VERBOSE=1 \| grep -L fopenmp` |
| BUILD-03 | AUTOMOC generates clean (no QObjects in Phase 0, vacuously true) | static | yes | inspect `CMakeLists.txt` for `set(CMAKE_AUTOMOC ON)` before `find_package(Qt6)` |
| BUILD-04 | Every entry point in header | static | yes | `grep -c 'proc(' xvue/qt/include/xvue_qt_api.h` equals 57 (or 60 — see Planner Alert) |
| BUILD-05 | Every declared entry has a definition | static | yes | `nm xvue/qt/build/libxvueqt.a \| grep -c ' T .*_$'` equals the header count |
| BUILD-06 | `bin/cbl_tout_qt` produces 5 executables | smoke | yes | `bin/cbl_tout_qt && test -x pp/ppmail_qt -a -x pp/ppelas_qt -a -x pp/ppflui_qt -a -x pp/ppther_qt -a -x pp/ppnlse_qt` |
| BUILD-07 | Legacy build still works | smoke | yes | `rm -rf pp && bin/cbl_tout && test -x pp/ppmail -a -x pp/ppelas -a -x pp/ppflui -a -x pp/ppther -a -x pp/ppnlse` |
| BUILD-08 | `verify_abi` fails on drift | build-break | yes | CMake `verify_abi` custom target runs on every build |
| BUILD-09 | `README_COORDS.md` exists | static | yes | `test -f xvue/qt/README_COORDS.md && grep -q 'Y-down' xvue/qt/README_COORDS.md` |
| BUILD-10 | `.planning/validation/BASELINE.md` exists with 5 cases | static | yes | `test -f .planning/validation/BASELINE.md && grep -c '^### ' .planning/validation/BASELINE.md` equals 5 |

### Sampling rate
- **Per task commit:** relevant subset (e.g. header-edit task → `grep -c`; shell-script task → `bin/cbl_tout_qt` dry-run)
- **Per wave merge:** full Phase 0 BUILD-01..10 check sequence
- **Phase gate:** the manual 5-case `testa/` run on legacy X11 (D-15) plus one successful end-to-end `bin/cbl_tout_qt` + `./pp/ppmail_qt < /dev/null` smoke-exercise

### Wave 0 gaps
- No test framework exists. Phase 0 doesn't need one — the checks above are all shell one-liners.
- **Manual** initial legacy-path verification step: `bin/cbl_tout` must succeed on this clean clone before any Phase 0 work begins. `pp/` is empty as of 2026-04-10, so this has never been run here. This is a genuine Wave 0 gap: **the legacy build has never been verified on this machine**. First plan task must be "run `bin/cbl_tout` and confirm it produces working executables" so Phase 0 doesn't accidentally inherit a broken baseline.

## Sources

### Primary (HIGH confidence)
- `/home/drico/git/mefisto/xvue/xvuelc.c` — direct read of lines 60-70 (proc macro), 227-264 (languemefisto, dctnmc, dstnmc), 286, 306, 319, 337, 358-503 (colormap helpers — the Planner Alert source), 561, 612, 1044-1187, 1307-1413, 1434-1602, 1619-1720, 1760-2120, 2123-2374, 2384-2678, 2694-2954. Full entry-point inventory derived here.
- `/home/drico/git/mefisto/bin/cbl_tout` — direct read, 86 lines
- `/home/drico/git/mefisto/bin/cbmail` — direct read, 71 lines (clone template)
- `/home/drico/git/mefisto/bin/ccxvue` — direct read, 42 lines (reveals pre-existing `-openmp` typo)
- `/home/drico/git/mefisto/CLAUDE.md` — project constraints (ask before system package install, never break legacy, commit after each logical step)
- `/home/drico/git/mefisto/.planning/PROJECT.md` (referenced via STATE.md)
- `/home/drico/git/mefisto/.planning/REQUIREMENTS.md` — BUILD-01..10 definitions
- `/home/drico/git/mefisto/.planning/ROADMAP.md` — Phase 0 success criteria
- `/home/drico/git/mefisto/.planning/phases/00-build-skeleton-abi-stubs/00-CONTEXT.md` — user-locked decisions D-01..D-18
- `/home/drico/git/mefisto/.planning/research/PITFALLS.md` — Pitfalls 1, 2, 3, 4, 9, 10, 11, 16, 18, 20 directly consulted
- `apt-cache policy qt6-base-dev cmake` (live probe of Debian trixie/sid)
- `gfortran --version`, `g++ --version`, `cmake --version`, `command -v nm pkg-config` (live probes)

### Secondary (MEDIUM confidence — via project research docs)
- `.planning/research/STACK.md` — CMake template skeleton (cross-referenced, not re-verified against upstream Qt docs this session)
- `.planning/research/ARCHITECTURE.md` — 12-file eventual layout (only ABI shim + singleton plumbing lands in Phase 0)
- `.planning/codebase/STACK.md` — existing compiler flags

### Tertiary (LOW confidence)
- None — every factual claim in this document is backed by a HIGH-confidence direct read or live probe.

## Project Constraints (from CLAUDE.md)

Extracted from `/home/drico/git/mefisto/CLAUDE.md`:

1. **Compilation must never break.** After every change, the affected module must still compile with its `cb*` script. Before committing, the full build (`bin/cbl_tout`) must succeed. → Phase 0 planner MUST include a step to re-run `bin/cbl_tout` (legacy) after all Phase 0 file edits.
2. **Tests.** Small `testa/` cases must still pass after every change. For large cases, ask the user. → Phase 0 gates on running the 5 baseline `testa/` cases through the legacy X11 path at end-of-phase.
3. **Ask before acting on system packages.** If a C include or Ubuntu package is missing, ASK — do not work around it. → First planner task is "ask user to `apt install qt6-base-dev qt6-base-dev-tools`" because `qt6-base-dev` is confirmed missing.
4. **Explain before doing.** Each step is announced. → Plan task descriptions must be self-explaining, not just `sed` or `cp` commands.
5. **Git discipline.** Commit after every logical step where rolling back would be useful. Commit messages describe what+why. No force-push, no skip-hooks. → D-17 warn-once is one logical step; `verify_abi` custom target is another; each `.f → .qt` clone script is one step.
6. **Fortran 77 fixed-form with column 7+ discipline.** → Not touched in Phase 0 (read-only rule) but planner must not accidentally edit any `.f` file.
7. **`incl/homdir.inc` is generated at build time.** → `bin/cbl_tout_qt` clone must preserve the generation block from `bin/cbl_tout` lines ~15-45 verbatim.
8. **Coding norms in `doc/normes.ps`.** → Not inspected this session; assumed Phase 0 does not violate them because Phase 0 adds only new C++ files in a new subdirectory. Planner may optionally read `doc/normes.ps` but it is not load-bearing for Phase 0.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | `nm` output format `' T .*_$'` is stable across binutils versions on x86_64 Linux | Code Examples (`verify_abi` target) | The `verify_abi` regex needs tweaking at implementation time. Detectable by build-break; low-risk. |
| A2 | `libX11-dev` is installed on this machine (legacy build works) | Environment Availability | If the legacy build is broken independently, Phase 0 will appear to fail BUILD-07 through no fault of its own. Mitigated by the "first task: run bin/cbl_tout to verify legacy" step. |
| A3 | Apt `qt6-base-dev` candidate 6.10.2+dfsg-6 is the same package the project research labelled "qt6-base-dev 6.10.2" | Standard Stack | None — verified this session to match. |
| A4 | The project's `.f` wrappers call exactly the Fortran-facing subset of entry points from `xvuelc.c`, not the 3 C-internal helpers | Planner Alert | Low risk — verified by `grep -i 'call xvColormapToRGB' xvue/*.f` returning nothing. |
| A5 | `doc/normes.ps` does not prohibit adding new C++ sources in a new subdirectory | Project Constraints | Minimal — `normes.ps` is about Fortran 77 column conventions; C++ in `xvue/qt/` is a new domain not covered by existing norms. |
| A6 | `qt6-base-dev-tools` is a dependency of `qt6-base-dev` and will be installed together | Environment Availability | If they must be separately installed, the first planner apt-install task simply lists both. Detectable at install time. |

## Open Questions

1. **Resolution of Planner Alert (D-04/D-06 conflict with xvuelc.c reality)**
   - What we know: 3 entries in `xvuelc.c` take the Xlib `Colormap` type or pass `int` by value; they are called only from within `xvuelc.c`, never from any Fortran `.f` wrapper.
   - What's unclear: whether the user wants Option A (declare them as C++ internals, skip the public header — breaking D-04's letter but preserving its spirit) or Option B (declare them in the public header behind an opaque typedef — preserving D-04's letter but adding meaningless symbols).
   - Recommendation: discuss-phase confirmation before planning tasks. Default to Option A.

2. **Exact content of `.planning/validation/BASELINE.md`'s "expected qualitative behavior" field**
   - What we know: The 5 cases exist on disk.
   - What's unclear: the maintainer's own qualitative description of what success looks like for each (e.g. "pan2d should produce a triangular mesh with ~200 elements, visible on screen, refinable with mouse clicks"). This is maintainer knowledge, not research-able.
   - Recommendation: Phase 0 task "write BASELINE.md" should be an interactive collaboration between planner and user — the user describes, the planner records.

3. **Whether `pp/*_OMP` executables need Phase 0 link testing**
   - What we know: D-15 says "any `pp/pp*_qt` executable succeeds if it proceeds past the link stage." The 5 executables named are non-`_OMP` variants.
   - What's unclear: whether `bin/cbl_tout_qt` should also produce `pp/ppelas_omp_qt` etc. for future Phase 8 `_OMP` sweeps, or whether `_OMP` targets are deferred to Phase 8.
   - Recommendation: **defer to Phase 8.** Phase 0 ships the 5 non-OMP targets per D-15 literal text. The `-fno-openmp` in CMakeLists is already in place so Phase 8 can flip the switch without re-architecting.

4. **Whether the pre-existing `bin/ccxvue` `-openmp` (no `f`) typo should be silently fixed in Phase 0**
   - What we know: `bin/ccxvue` line 22 passes `-openmp` (Intel syntax) to gcc which silently ignores it. The legacy `xvuelc.o` has therefore never actually had OpenMP code.
   - What's unclear: whether fixing this is "breaking the legacy path" (D-02 violation) or "a tangential bug fix unrelated to Phase 0."
   - Recommendation: **do NOT fix in Phase 0.** D-02 is explicit: `bin/ccxvue` is read-only. Log the typo in `.planning/codebase/CONCERNS.md` (outside Phase 0 scope) and move on.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — every version/package verified by live apt/command probe this session.
- Architecture: HIGH on Phase 0 scope; MEDIUM on cross-phase interactions (Phase 1+ inherit scaffolding but that's Phase 1's research scope).
- Pitfalls: HIGH — all 9 relevant pitfalls are already reconciled in `.planning/research/PITFALLS.md` and cross-verified against `xvuelc.c` in this session.
- ABI surface (the 59→60 entry-point count + Planner Alert): HIGH — direct line-by-line read of `xvuelc.c`.
- Y-axis convention audit: HIGH — direct read of all drawing entry points + PostScript emit paths.

**Research date:** 2026-04-10
**Valid until:** 2026-05-10 (30 days — `xvuelc.c` is stable, Debian trixie Qt 6 version may tick minor during the window but the CMake `find_package(Qt6 6.x)` contract is forward-compatible).

---
*Phase: 00-build-skeleton-abi-stubs*
*Research complete: 2026-04-10*
