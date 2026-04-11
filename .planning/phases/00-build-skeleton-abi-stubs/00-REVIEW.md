---
phase: 00-build-skeleton-abi-stubs
reviewed: 2026-04-10T00:00:00Z
depth: standard
files_reviewed: 11
files_reviewed_list:
  - bin/cbl_tout_qt
  - bin/cbmail_qt
  - bin/cbelas_qt
  - bin/cbflui_qt
  - bin/cbther_qt
  - bin/cbnlse_qt
  - xvue/qt/CMakeLists.txt
  - xvue/qt/cmake/verify_abi.sh
  - xvue/qt/.gitignore
  - xvue/qt/include/xvue_qt_api.h
  - xvue/qt/src/xvue_qt_api.cpp
  - xvue/README_COORDS.md
findings:
  critical: 0
  warning: 3
  info: 6
  total: 9
status: issues_found
---

# Phase 0: Code Review Report

**Reviewed:** 2026-04-10
**Depth:** standard
**Files Reviewed:** 11 (12 listed in scope — xvue/qt/src/xvue_qt_api.cpp counted once)
**Status:** issues_found

## Summary

Phase 0 delivers a clean build skeleton for the Qt 6 replacement of xvuelc.o. The core deliverables are consistent and internally coherent:

- **ABI count is correct.** `xvue/qt/include/xvue_qt_api.h` declares exactly 57 Fortran-facing entry points (matches Planner Alert Option A). `xvue/qt/src/xvue_qt_api.cpp` provides exactly 57 stub definitions. Spot-checked signatures against `xvue/xvuelc.c` (xvinfo@612, xvsouris2@2230, valvarenv@2805, xvface@1701, xvtraits@1977, xvfacetraits@2035, xvtrait@1862) all match byte-for-byte including parameter types and layout.
- **Trailing-underscore mangling is preserved.** The `proc(x) x##_` macro is copied verbatim from `xvuelc.c` lines 60-70 and applies uniformly to every extern "C" entry.
- **MefistoPoint shim is sound.** `struct { short x; short y; }` with `static_assert(sizeof == 4)` correctly shadows Xlib's `XPoint` binary layout (D-07, Pitfall 4).
- **CMake configuration follows D-11/D-12.** `AUTOMOC` precedes `find_package(Qt6)`, `POSITION_INDEPENDENT_CODE` is ON, C++17 no extensions, Qt components Core/Gui/Widgets/PrintSupport. `verify_abi` runs on every build as an ALL target and fails on header/archive drift.
- **verify_abi.sh regex matches the header.** I ran the header-count regex against `xvue_qt_api.h` and it returns exactly 57 (verified independently).
- **Stub behavior follows D-17/D-18.** Each stub uses the warn-once `static bool warned` pattern, returns safe defaults (nullptr, no writes), and never touches Qt state. Parameters are fully `(void)`-cast to suppress `-Wunused-parameter`.
- **Y-axis convention document conclusion is correct.** Confirmed by reading `xvuelc.c` that every on-screen draw call passes Y unflipped and every PostScript emit applies `ypixels - y` — matching the document's Phase 1-7 action table.

The issues below are all in the Warning/Info bands. No Critical findings.

## Warnings

### WR-01: cbl_tout_qt does not validate $MEFISTO before rm -rf

**File:** `bin/cbl_tout_qt:49-52`
**Issue:** Lines 51-52 run `rm -rf $MEFISTO/xvue/qt/build` and `rm -f $MEFISTO/pp/pp*_qt` without verifying that `$MEFISTO` is set and non-empty. If a user runs `cbl_tout_qt` without sourcing their environment (or with a typo in an `export`), `$MEFISTO` expands to empty and the commands become `rm -rf /xvue/qt/build` and `rm -f /pp/pp*_qt`. On this particular system those paths do not exist so the commands are harmless, but the pattern is a foot-gun the moment anyone runs the script from an unexpected shell. The legacy `cbl_tout` has the same pattern, but Phase 0 is the right moment to fix it in the _qt variant since it is a fresh copy.
**Fix:**
```bash
# near the top of the script, right after `echo MEFISTO=$MEFISTO`
: "${MEFISTO:?MEFISTO is not set — source your environment first}"
if [ ! -d "$MEFISTO" ]; then
    echo "MEFISTO=$MEFISTO is not a directory — aborting" >&2
    exit 1
fi
```
Apply the same `: "${MEFISTO:?...}"` guard at the top of `cbmail_qt`, `cbelas_qt`, `cbflui_qt`, `cbther_qt`, `cbnlse_qt` since they all `cd $MEFISTO` on line 7-8.

### WR-02: cbl_tout_qt does not check exit status of the five sub-build scripts

**File:** `bin/cbl_tout_qt:91-95`
**Issue:** After successfully building `libxvueqt.a` via CMake (which does fail the script on error, lines 56-67), the script invokes `cbmail_qt`, `cbelas_qt`, `cbther_qt`, `cbflui_qt`, `cbnlse_qt` sequentially without any `|| exit 1` guard. If any one of those fails to produce its `pp/pp*_qt` binary, the script prints a localized "INCORRECT" message from that sub-script but `cbl_tout_qt` still exits 0 after the final `ls -l $MEFISTO/pp`. This breaks any CI or Makefile wrapper that relies on exit status to detect a bad build.
**Fix:**
```bash
$MEFISTO/bin/cbmail_qt || { echo "cbmail_qt FAILED" >&2; exit 1; }
$MEFISTO/bin/cbelas_qt || { echo "cbelas_qt FAILED" >&2; exit 1; }
$MEFISTO/bin/cbther_qt || { echo "cbther_qt FAILED" >&2; exit 1; }
$MEFISTO/bin/cbflui_qt || { echo "cbflui_qt FAILED" >&2; exit 1; }
$MEFISTO/bin/cbnlse_qt || { echo "cbnlse_qt FAILED" >&2; exit 1; }
```
The sub-scripts themselves should also `exit 1` on gfortran failure — currently they only print a message. Consider adding `set -e` at the top of each `cb*_qt`, or explicitly testing `$?` after the `gfortran` call and exiting non-zero when the target file is absent.

### WR-03: README_COORDS.md cites xvuelc.c line numbers that don't contain the quoted text

**File:** `xvue/README_COORDS.md:23-26`
**Issue:** The document's "Where the convention comes from" section attributes specific comment strings to specific line numbers, but two of the four citations are inaccurate when verified against `xvue/xvuelc.c`:

1. Line 321 is cited as containing "origine coin superieur gauche de l'ecran" inside `xvpxecran_`. Actual line 321 reads `BUT :   RECUPERER le nombre de pixels de la largeur et hauteur de l'ecran total`. The quoted phrase "origine coin superieur gauche de l'ecran" does not appear anywhere in the `xvpxecran_` block (verified via grep -i for "origine coin superieur gauche" across the whole file — only four hits, all in the line-drawing block at 1852, 1853, 1869, 1870).
2. Line 1621 is cited as containing "la fenetre" comment "inside `xvpxfenetre_`" with "origin implied top-left". Actual line 1621 reads `BUT :     RECUPERER le nombre de pixels de la largeur et hauteur de`. This is a size-query function; it does not document the origin at all.

Citations 2 and 3 (lines 1852, 1869) are accurate, and the PostScript-flip citations (1895, 1932, 1953, 1966) are all accurate. The conclusion of the document is still correct — Xlib is Y-down top-left by construction, and the PS emit does apply `ypixels - y` — but a phase 1 reader who trusts the direct-evidence table and follows the footnote back to `xvuelc.c:321` will be confused.
**Fix:** Drop citations 1 and 4 (xvpxecran, xvpxfenetre). They do not add evidence because those functions do not document origin semantics. Keep citations 2 and 3 (xvtrait block at 1852/1869) and add a stronger line: "Every `XDrawLine` / `XDrawLines` / `XFillRectangle` / `XDrawArc` call in `xvuelc.c` passes the Fortran-supplied `*y` / `*y1` / `*y2` as-is to Xlib — which is the whole proof. No call site flips y on the on-screen path." Or alternatively, replace citation 1 with the actual xvtrait comments and remove the xvpxecran/xvpxfenetre bullets.

## Info

### IN-01: Dead preprocessor block in xvue_qt_api.h

**File:** `xvue/qt/include/xvue_qt_api.h:21-27`
**Issue:** Lines 21-25 define `proc(x)` inside `#ifdef __GNUC__` / `#else` / `#endif`, then line 26 `#undef proc`, then line 27 unconditionally `#define proc(x) x##_`. The first block is entirely dead — its definitions are undefined before any use. The comment on line 20 says "copied verbatim from xvuelc.c lines 60-70. Do not change." which is the intent, but the dead branch adds zero semantic value to a C++-only header and is a small maintenance cost (it tempts future readers to reason about the non-GCC path that will never exist).
**Fix:** Since the phase's explicit goal is byte-identical copy of xvuelc.c's mangling block, leave it alone for Phase 0. Add a one-line comment after line 27: `// Lines 21-25 kept only for byte-identical parity with xvuelc.c:60-70; they are dead.` This preserves auditability while warning the reader.

### IN-02: xvue_qt_api.h uses C++ features but has a .h extension

**File:** `xvue/qt/include/xvue_qt_api.h:15,34,41-48`
**Issue:** The header includes `<cstddef>` (C++-only), uses `static_assert` with two arguments (C++11 form; C11 spells it `_Static_assert`), and includes `<QThread>` / `<QCoreApplication>` inside `#ifdef QT_DEBUG`. A `.h` extension conventionally implies "safe to include from C", but this file is not. The only current consumer is `xvue_qt_api.cpp` so there's no real breakage, but a future tool or IDE that treats `.h` as C will choke. Naming it `.hpp` would signal the intent.
**Fix:** Low-priority rename to `xvue_qt_api.hpp` when convenient. If the rename is deferred, add a `#ifndef __cplusplus` / `#error "xvue_qt_api.h is C++ only"` / `#endif` guard at the top of the file.

### IN-03: verify_abi.sh header regex is fragile against new return types

**File:** `xvue/qt/cmake/verify_abi.sh:18`
**Issue:** The BRE alternation `void\|int\|float\|double\|long\|short\|unsigned\|void[[:space:]]*\*` only counts declarations whose return type starts with one of those tokens. If a future entry returns e.g. `char *`, `size_t`, `MefistoPoint`, or `bool`, the grep will silently under-count the header, and the drift check will hard-fail _even though_ the header and implementation are consistent. For Phase 0 the 57 signatures happen to all fit the allowlist, so it works today.
**Fix:** Replace the header count with a simpler proxy that doesn't depend on return-type text. Candidates:
- Count lines matching `^[^#/].*\bproc\(` (exclude `#define` / `#undef` / `//` comments).
- Or count pragma-guarded entries by wrapping each declaration in a macro, e.g. `FORTRAN_ENTRY(void, languemefisto, (int *langue))`, then grep for `FORTRAN_ENTRY(`.
Not urgent — document the limitation in a comment for now so the next maintainer who adds a `char *` entry isn't surprised by the false drift.

### IN-04: verify_abi.sh nm pattern will break if a symbol has digits in the first character

**File:** `xvue/qt/cmake/verify_abi.sh:17`
**Issue:** The regex `' T [a-zA-Z_][a-zA-Z0-9_]*_$'` requires the symbol name to start with a letter or underscore and end with `_`. All 57 current Fortran entry points match, so this is purely future-proofing. The regex also excludes `t` (lowercase, local symbols), which is correct. Grade: fine for Phase 0, just noting the constraint exists.
**Fix:** No action needed. Comment mentioning "extern trailing-underscore convention only — not local symbols" would help the next reader.

### IN-05: bin/cb*_qt scripts duplicate QT_LIBS invocation and link flags

**File:** `bin/cbmail_qt:37`, `bin/cbelas_qt:38`, `bin/cbflui_qt:37`, `bin/cbther_qt:37`, `bin/cbnlse_qt:37`
**Issue:** Five scripts each re-run `pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport` and each spell out `-Lxvue/qt/build -lxvueqt $QT_LIBS -lstdc++` independently. If Phase 1 adds a Qt component (e.g. `Qt6OpenGLWidgets`) or drops one, all five scripts need to change in lockstep — easy to miss. Also: if `pkg-config` is not installed or Qt6 `.pc` files are not found, `QT_LIBS` silently becomes empty and gfortran will fail with a long cryptic undefined-reference error, not a helpful message.
**Fix:** Extract the shared bits into `bin/cb_qt_common.sh` (sourced by each `cb*_qt`), which sets `QT_LIBS` and `XVUEQT_LINK_FLAGS` once and validates that `pkg-config` succeeded:
```bash
# bin/cb_qt_common.sh — sourced by cb*_qt
QT_LIBS=$(pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport) || {
    echo "pkg-config failed to find Qt6 components" >&2
    exit 1
}
[ -n "$QT_LIBS" ] || { echo "QT_LIBS is empty — Qt6 .pc files missing?" >&2; exit 1; }
XVUEQT_LINK_FLAGS="-Lxvue/qt/build -lxvueqt $QT_LIBS -lstdc++"
```
Not urgent for Phase 0, but flagging before Phase 1 copies and multiplies the pattern.

### IN-06: Inconsistent gfortran flags between cb*_qt variants

**File:** `bin/cbmail_qt:39`, `bin/cbelas_qt:40`, `bin/cbflui_qt:39`, `bin/cbther_qt:39`, `bin/cbnlse_qt:39`
**Issue:** Each script uses a slightly different subset of gfortran flags:
- `cbmail_qt`: `-Wall -mcmodel=large -m64 -O -fopenmp`
- `cbelas_qt`: `-m64 -mcmodel=large -fPIC -O -fopenmp` (no -Wall)
- `cbflui_qt`: `-m64 -Wall -mcmodel=large -fPIC -O -fopenmp`
- `cbther_qt`: `-m64 -fPIC -O -fopenmp` (no -Wall, no -mcmodel=large)
- `cbnlse_qt`: `-m64 -fPIC -O -fopenmp` (no -Wall, no -mcmodel=large)

Notably `cbther_qt` and `cbnlse_qt` drop `-mcmodel=large`, which may cause relocation truncation errors on large static data arrays once linking against libxvueqt.a (pulling in Qt's large static data sections). `cbmail_qt` drops `-fPIC`, which may cause an issue because libxvueqt.a was built with `POSITION_INDEPENDENT_CODE ON` — link-compatible but slightly wasteful. These inconsistencies are inherited from the legacy `cb*` scripts and are not a Phase 0 regression, but since the five _qt scripts are fresh copies and pull in a common new dependency (libxvueqt.a + Qt libs), harmonizing them now would save debugging later.
**Fix:** Adopt a uniform set for all five: `-m64 -mcmodel=large -fPIC -Wall -O -fopenmp`. Verify the full build still succeeds afterward. Defer if `cbther_qt` / `cbnlse_qt` currently rely on the absence of `-mcmodel=large` for some compatibility reason, but the research doc doesn't note any such requirement.

---

_Reviewed: 2026-04-10_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
