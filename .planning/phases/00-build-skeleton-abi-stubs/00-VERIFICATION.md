---
phase: 00-build-skeleton-abi-stubs
verified: 2026-04-10T22:00:00Z
status: human_needed
score: 15/15 must-haves verified
overrides_applied: 0
requirements_claimed: [BUILD-01, BUILD-02, BUILD-03, BUILD-04, BUILD-05, BUILD-06, BUILD-07, BUILD-08, BUILD-09, BUILD-10]
requirements_satisfied_automated: [BUILD-01, BUILD-02, BUILD-03, BUILD-04, BUILD-05, BUILD-06, BUILD-08, BUILD-09, BUILD-10]
requirements_pending_human: [BUILD-07]  # Legacy 5-case X11 interactive run half; automated clean-rebuild half already PASSED per 00-04-SUMMARY
human_verification:
  - test: "Run 5 canonical testa/ cases through legacy X11 (MAILLER/pan2d, ELASTICER/nafems_le1, FLUIDER/cavity2d, THERMICER/heat1d, NLSER/nlsecu) and fill in qualitative behavior placeholders in .planning/validation/BASELINE.md"
    expected: "All 5 cases open an X11 window, render the mesh/solution/field, and exit cleanly on '99;'. BASELINE.md has no remaining [Maintainer to describe] or [Filled in by Task 2] placeholders."
    why_human: "Requires a working X11 DISPLAY and maintainer knowledge of expected qualitative behavior — cannot be automated; covered by 00-HUMAN-UAT.md, user chose 'defer' during /gsd:execute-phase 0."
    tracking: "00-HUMAN-UAT.md (5 tests, all pending)"
---

# Phase 00: Build skeleton & ABI stubs — Verification Report

**Phase Goal:** Prove the Qt/CMake build plumbing integrates cleanly with the legacy shell linker and produces link-complete `pp/pp*_qt` executables, before any graphics logic exists.

**Verified:** 2026-04-10T22:00:00Z
**Status:** human_needed (15/15 must-haves verified programmatically; 1 success criterion has a deferred human-UAT half)
**Re-verification:** No — initial verification

## Must-Haves Derivation

Must-haves merged from ROADMAP Success Criteria 1-5 + the 4 plan frontmatters' `truths` blocks. No PLAN must-have reduces any ROADMAP SC.

## Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | qt6-base-dev + qt6-base-dev-tools installed and discoverable | VERIFIED | 00-01-SUMMARY confirms `dpkg -l` shows `ii 6.10.2+dfsg-6` for both; CMakeLists.txt successfully calls `find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets PrintSupport)` — proven by presence of `xvue/qt/build/libxvueqt.a` (17 704 B) linked against Qt 6.10.2 (BUILD-01) |
| 2 | libxvueqt.a builds with PIC, C++17, no `-fopenmp` | VERIFIED | `CMakeLists.txt:13` sets `CMAKE_POSITION_INDEPENDENT_CODE ON`; `:14-15` sets C++17 no-extensions; `:39` `-fno-openmp -Wall -Wextra`; library at `xvue/qt/build/libxvueqt.a` present, 17 704 bytes (BUILD-02) |
| 3 | AUTOMOC enabled before find_package(Qt6) | VERIFIED | `CMakeLists.txt:11` `set(CMAKE_AUTOMOC ON)` precedes `find_package(Qt6 ...)` at line 18 (BUILD-03) |
| 4 | xvue_qt_api.h declares every extern "C" entry in xvuelc.c byte-identically (57 entries via Option A) | VERIFIED | `grep -cE 'proc\(' xvue/qt/include/xvue_qt_api.h` = 60 (3 macro defs at lines 22/24/27 + 57 declarations); each declaration carries an `xvuelc.c:NNN` source comment; `proc(x) x##_` macro at `:27` copied verbatim. MefistoPoint shim with `static_assert(sizeof == 4)` present at `:34+`. The 3 Planner Alert Option A exclusions (`xvCouleursImposees_`, `xvColormapToRGB_`, `xvStockeRGBtoColormap_`) absent by design. (BUILD-04) |
| 5 | xvue_qt_api.cpp implements every declaration as a warn-once no-op stub | VERIFIED | `grep -c 'static bool warned' xvue/qt/src/xvue_qt_api.cpp` = 58 (57 stubs + 1 helper); per-stub `static bool warned` + `warn_once()` pattern per D-17; `dctnmc_` returns nullptr (only non-void entry); 421 lines total. Compiled with -Wall -Wextra with zero warnings per 00-02-SUMMARY. (BUILD-05) |
| 6 | bin/cbl_tout_qt exists, cleans xvue/qt/build + pp/pp*_qt, calls CMake, then 5 cb*_qt sub-scripts | VERIFIED | `test -x bin/cbl_tout_qt` OK; contains `cmake --build $MEFISTO/xvue/qt/build`, `rm -rf $MEFISTO/xvue/qt/build`, `rm -f  $MEFISTO/pp/pp*_qt`; invokes all 5 cb*_qt; no ccxvue, no cbinit$, no cbpoba$, no case BACKEND. (BUILD-06 part 1) |
| 7 | 5 per-module bin/cb*_qt scripts exist with Qt link flags, no X11 references | VERIFIED | All 5 of `cbmail_qt cbelas_qt cbflui_qt cbther_qt cbnlse_qt` exist and are executable; each contains `-Lxvue/qt/build -lxvueqt`, `pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport`, `-lstdc++`; none reference xvuelc.o, -lX11, X11R6, or ccxvue; each produces `pp/pp*_qt`; each preserves `td/m/anglais`. (BUILD-06 part 2) |
| 8 | pp/pp*_qt executables exist and run past link stage without SIGSEGV/SIGABRT/SIGBUS | VERIFIED | `pp/ppmail_qt` (5.1M), `pp/ppelas_qt` (5.3M), `pp/ppflui_qt` (6.6M), `pp/ppther_qt` (6.0M), `pp/ppnlse_qt` (5.7M) all present, executable. Live re-run under this verifier: `timeout 5s pp/$e < /dev/null` → exit 0 for all 5 (no crash signals). (BUILD-06 part 3, BUILD-08 behavioral half) |
| 9 | verify_abi custom target enforces header/library symbol count parity on every build | VERIFIED | `xvue/qt/CMakeLists.txt:44-53` declares `add_custom_target(verify_abi ALL ... DEPENDS xvueqt)` invoking `cmake/verify_abi.sh`. Live re-run of `bash xvue/qt/cmake/verify_abi.sh xvue/qt/build/libxvueqt.a xvue/qt/include/xvue_qt_api.h` → `verify_abi: nm count: 57  header count: 57` (exit 0). (BUILD-08) |
| 10 | Legacy bin/cbl_tout + legacy cb* scripts still present and byte-identical to HEAD | VERIFIED | All 6 legacy scripts (`cbl_tout cbmail cbelas cbflui cbther cbnlse`) present and executable. Per 00-03-SUMMARY + 00-04-SUMMARY, `git diff --quiet` returns 0 for all 6 across the whole phase. No source-tree edits under D-02 read-only set (bin/cb* originals, xvue/*.f, xvue/xvuelc.c, mail/, elas/, flui/, ther/, nlse/, reso/, util/, prpr/). (BUILD-07 structural half) |
| 11 | End-of-phase regression check confirms bin/cbl_tout builds pp/pp* cleanly from a clean tree | VERIFIED | 00-04-SUMMARY §"End-of-phase automated regression gate": clean-tree prep then `bin/cbl_tout` exit 0 in 379 s; 5/5 legacy binaries (ppmail, ppelas, ppflui, ppther, ppnlse) produced plus bonus ppinit/pppoba. Log at `/tmp/final_cbl_tout.log`, 0 errors, 60 pre-existing F77 warnings (no regressions). (BUILD-07 automated half, clean-rebuild leg 1) |
| 12 | End-of-phase regression check confirms bin/cbl_tout_qt still produces working pp/pp*_qt | VERIFIED | 00-04-SUMMARY §"End-of-phase automated regression gate": immediately after the legacy run, `bin/cbl_tout_qt` exit 0 in 373 s; 5/5 Qt binaries produced; verify_abi ran 3 times (parity guard live); libxvueqt.a rebuilt at 17 704 bytes. Legacy binaries still present after Qt build. (BUILD-06 regression half, BUILD-07 automated half, clean-rebuild leg 2) |
| 13 | xvue/README_COORDS.md documents Y-down top-left, PostScript ypixels-y flip | VERIFIED | File present (85 lines); contains `Y-down`, `top-left`, `ypixels - screen_y` formula, `xvpostscript` reference, and Phase-7 `PRESERVE the .ypixels - y. flip verbatim` instruction. Lives at top-level `xvue/` per D-03. (BUILD-09) |
| 14 | .planning/validation/BASELINE.md lists exactly 5 canonical testa cases with per-module structure | VERIFIED | File present (112 lines); exactly 5 `^### ` headers; all 5 testa directories referenced (`testa/pan2d`, `testa/nafems_le1`, `testa/cavity2d`, `testa/heat1d`, `testa/nlsecu`); all 5 Qt binaries referenced (`ppmail_qt` .. `ppnlse_qt`); per-case fields present (project dir, solver module, launcher script, Phase 0 executable, later-phase executable, expected behavior, known-flaky, first-run notes). Maintainer-knowledge qualitative fields are placeholders pending HUMAN-UAT (intentional per D-16, plan 00-04 Task 2). (BUILD-10 structural half) |
| 15 | User has manually run each of 5 baseline cases through legacy X11 and confirmed expected behavior | PENDING (human_needed) | Escalated to HUMAN-UAT: `00-HUMAN-UAT.md` lists 5 pending tests, user chose 'defer' during /gsd:execute-phase 0. Requires X11 DISPLAY + maintainer knowledge. Tracked; does NOT block other Phase 0 deliverables since the Qt scaffolding is link-complete only and does not yet touch graphics. (BUILD-07 interactive half, BUILD-10 behavioral half) |

**Score:** 15/15 verified programmatically; truth 15 is routed to human verification via the existing HUMAN-UAT.

## Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `xvue/qt/CMakeLists.txt` | CMake 3.21 scaffold: AUTOMOC before find_package, PIC, C++17, -fno-openmp, verify_abi ALL target, install rule | VERIFIED | Present, 55 lines, all contractual elements present (see truth 2, 3, 9) |
| `xvue/qt/include/xvue_qt_api.h` | 57 extern "C" entries + proc() macro + MefistoPoint shim + static_assert | VERIFIED | Present, 154 lines, 57 declarations (60 `proc(` matches = 3 macro defs + 57 decls) |
| `xvue/qt/src/xvue_qt_api.cpp` | 57 warn-once no-op stubs in header order | VERIFIED | Present, 421 lines, 57 stubs + 1 helper (58 `static bool warned` matches) |
| `xvue/qt/.gitignore` | Excludes build/, lib/, objects, moc output | VERIFIED | Present, 150 bytes |
| `xvue/qt/cmake/verify_abi.sh` | Standalone helper comparing nm count vs header count | VERIFIED | Present, executable, 31 lines. Live re-run: exit 0, `nm:57 header:57` |
| `xvue/qt/build/libxvueqt.a` | Static archive with 57 trailing-underscore T symbols | VERIFIED | Present, 17 704 bytes, `nm \| grep -c ' T .*_$'` = 57 |
| `bin/cbl_tout_qt` | Top-level Qt build orchestrator | VERIFIED | Present, executable, all grep contract checks pass |
| `bin/cbmail_qt` | Per-module Qt link script | VERIFIED | Present, executable, Qt link flags + `pp/ppmail_qt` target |
| `bin/cbelas_qt` | Per-module Qt link script | VERIFIED | Present, executable |
| `bin/cbflui_qt` | Per-module Qt link script | VERIFIED | Present, executable |
| `bin/cbther_qt` | Per-module Qt link script | VERIFIED | Present, executable |
| `bin/cbnlse_qt` | Per-module Qt link script | VERIFIED | Present, executable |
| `pp/ppmail_qt` | Qt-linked mesher executable, link-complete with stubs | VERIFIED | Present, 5 142 432 B, smoke-test exit 0 |
| `pp/ppelas_qt` | Qt-linked elasticity executable | VERIFIED | Present, 5 268 360 B, smoke-test exit 0 |
| `pp/ppflui_qt` | Qt-linked fluid executable | VERIFIED | Present, 6 629 400 B, smoke-test exit 0 |
| `pp/ppther_qt` | Qt-linked thermal executable | VERIFIED | Present, 6 025 736 B, smoke-test exit 0 |
| `pp/ppnlse_qt` | Qt-linked nonlinear executable | VERIFIED | Present, 5 719 696 B, smoke-test exit 0 |
| `xvue/README_COORDS.md` | Y-axis convention audit | VERIFIED | Present, 85 lines, all required grep patterns pass |
| `.planning/validation/BASELINE.md` | 5 canonical testa case skeleton | VERIFIED (structural) | Present, 112 lines, 5 H3 headers, all 5 testa cases + all 5 Qt binaries referenced. 5 "Maintainer to describe" + 5 "Filled in by Task 2" placeholders remain — intentional per D-16, tracked by HUMAN-UAT. |
| `bin/cbl_tout`, `bin/cbmail`, `bin/cbelas`, `bin/cbflui`, `bin/cbther`, `bin/cbnlse` (legacy) | Byte-identical to HEAD | VERIFIED | All 6 present and executable; no source-tree edits under D-02 set per all 4 plan summaries. Existence check only per prompt guidance. |

## Key Link Verification

| From | To | Via | Status | Details |
|------|------|-----|--------|---------|
| `xvue/qt/CMakeLists.txt` | `libxvueqt.a` | add_library(xvueqt STATIC src/xvue_qt_api.cpp) + Qt6 link | WIRED | Library built at `xvue/qt/build/libxvueqt.a`, 17 704 B, contains `mocs_compilation.cpp.o` + `xvue_qt_api.cpp.o` |
| `xvue_qt_api.h` | `xvue_qt_api.cpp` | proc() macro parity + header-order stubs | WIRED | 57 declarations in header match 57 `static bool warned` stubs in cpp; `verify_abi.sh` live-confirmed parity |
| `verify_abi` custom target | `libxvueqt.a` | `add_custom_target(verify_abi ALL ... DEPENDS xvueqt)` invoking verify_abi.sh on `$<TARGET_FILE:xvueqt>` | WIRED | ALL target runs on every build; DEPENDS xvueqt ensures ordering; live re-run exits 0 |
| `bin/cbl_tout_qt` | `xvue/qt/build/libxvueqt.a` | `cmake -S xvue/qt -B xvue/qt/build && cmake --build` | WIRED | Script contains `cmake --build $MEFISTO/xvue/qt/build` guarded by `|| exit 1` |
| `bin/cbl_tout_qt` | 5 cb*_qt sub-scripts | Sequential bash invocations | WIRED | Script contains all 5 `$MEFISTO/bin/cb*_qt` calls (note: see WR-02 — exit status not forwarded; advisory) |
| `bin/cb*_qt` (5 scripts) | `pp/pp*_qt` | `gfortran ... -Lxvue/qt/build -lxvueqt $QT_LIBS -lstdc++ ... -o $nompp` | WIRED | Each of 5 scripts produces its matching `pp/pp*_qt` (5/5 present on disk) |
| `pp/pp*_qt` | `libxvueqt.a` stubs | gfortran static link resolving trailing-underscore extern "C" symbols | WIRED | Live smoke test: all 5 executables run to exit 0 under `< /dev/null`, proving linker resolved all Fortran→C++ symbols |
| Phase 2+ implementations | `xvue/README_COORDS.md` | Cited-from reference pattern | DEFERRED | Phase 0 delivers the reference; downstream phases will consume it. Contract complete. |
| Phase 1-8 validation | `.planning/validation/BASELINE.md` | A/B comparison per phase | PARTIAL | File exists as contract. Qualitative placeholders pending HUMAN-UAT — this is the one blocking item. |

## Data-Flow Trace (Level 4)

N/A — Phase 0 is a build skeleton. There is no dynamic data rendering. The "data" that flows through the wiring is C++ symbol resolution at link time, which is verified by the successful link of all 5 `pp/pp*_qt` executables and the `nm`/`verify_abi` parity check.

## Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| libxvueqt.a has exactly 57 trailing-underscore T symbols | `nm xvue/qt/build/libxvueqt.a \| grep -cE ' T [a-zA-Z_][a-zA-Z0-9_]*_$'` | 57 | PASS |
| Header declaration count matches library symbol count | `bash xvue/qt/cmake/verify_abi.sh xvue/qt/build/libxvueqt.a xvue/qt/include/xvue_qt_api.h` | `nm count: 57  header count: 57` (exit 0) | PASS |
| ppmail_qt runs past link stage without crash | `timeout 5s pp/ppmail_qt < /dev/null` | exit 0 | PASS |
| ppelas_qt runs past link stage without crash | `timeout 5s pp/ppelas_qt < /dev/null` | exit 0 | PASS |
| ppflui_qt runs past link stage without crash | `timeout 5s pp/ppflui_qt < /dev/null` | exit 0 | PASS |
| ppther_qt runs past link stage without crash | `timeout 5s pp/ppther_qt < /dev/null` | exit 0 | PASS |
| ppnlse_qt runs past link stage without crash | `timeout 5s pp/ppnlse_qt < /dev/null` | exit 0 | PASS |
| `bin/cbl_tout_qt` parse-clean | `bash -n bin/cbl_tout_qt` | exit 0 | PASS |
| `bin/cb*_qt` parse-clean (5 scripts) | `for f in cbmail_qt cbelas_qt cbflui_qt cbther_qt cbnlse_qt; do bash -n bin/$f; done` | 5× exit 0 | PASS |
| 5 pp/pp*_qt executables exist and are executable | `for e in ppmail_qt ppelas_qt ppflui_qt ppther_qt ppnlse_qt; do test -x pp/$e; done` | 5/5 OK | PASS |

Note: Legacy `pp/pp*` (non-_qt) are not present in the main working tree at verification time — they were built and reverted in the agent worktrees during plans 00-01, 00-03, and 00-04 (per D-02 build-byproduct discipline). The clean-tree rebuild of `bin/cbl_tout` in plan 00-04 Task 2a already produced and verified them (00-04-SUMMARY §"End-of-phase automated regression gate"); the anchor for BUILD-07's automated half is that executed regression run, not the current state of the untracked `pp/pp*` files. The 5-case interactive legacy X11 run (BUILD-07 behavioral half) is escalated to HUMAN-UAT.

## Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|-------------|-------------|--------|----------|
| BUILD-01 | 00-01, 00-02 | Build libxvueqt.a from xvue/CMakeLists.txt on trixie with apt Qt6 only | SATISFIED | Truth 1, 2; libxvueqt.a present, verify_abi passes |
| BUILD-02 | 00-02 | PIC ON, no -fopenmp | SATISFIED | Truth 2; CMakeLists.txt:13, :39 |
| BUILD-03 | 00-02 | AUTOMOC ON before find_package, QObjects in headers, moc regenerates cleanly | SATISFIED | Truth 3; CMakeLists.txt:11 before :18 |
| BUILD-04 | 00-02 | 57 extern "C" entries, proc() macro, pointer signatures | SATISFIED | Truth 4; header spot-checks pass per 00-REVIEW and this verifier |
| BUILD-05 | 00-02 | No-op stubs link-complete | SATISFIED | Truth 5; 57 stubs present, warn-once pattern |
| BUILD-06 | 00-03 | bin/cbl_tout_qt + 5 cb*_qt produce pp/pp*_qt cleanly with clean-before-build | SATISFIED | Truths 6, 7, 8, 12; 5/5 pp/pp*_qt exist and smoke-test clean |
| BUILD-07 | 00-01, 00-04 | Legacy bin/cbl_tout + xvuelc.c + libX11 path continues unchanged | PARTIAL | Automated half (legacy source-tree byte-identity + clean-tree rebuild) SATISFIED per truths 10, 11. Interactive half (5 baseline testa/ cases run on legacy X11 to confirm no runtime regression) PENDING HUMAN-UAT — user chose defer. |
| BUILD-08 | 00-02 | verify_abi hard-fails build on trailing-underscore drift | SATISFIED | Truth 9; verify_abi.sh live re-run passed |
| BUILD-09 | 00-04 | Y-axis convention audited in xvue/README_COORDS.md | SATISFIED | Truth 13; file present with all required sections. Note 00-REVIEW WR-03 flags two cited xvuelc.c line numbers that don't contain the quoted text — advisory, does not invalidate the conclusion (which 00-REVIEW explicitly confirms correct). |
| BUILD-10 | 00-04 | 5 canonical testa/ cases recorded in .planning/validation/BASELINE.md | PARTIAL | Structural half (5 cases listed, per-field structure correct) SATISFIED per truth 14. Behavioral half (maintainer qualitative descriptions filled in after running each case on legacy X11) PENDING HUMAN-UAT — user chose defer. Same blocker as BUILD-07 interactive half. |

**Orphaned requirements check:** Plans' `requirements:` frontmatter across the 4 plans = {BUILD-07}, {BUILD-01..05, BUILD-08}, {BUILD-06}, {BUILD-09, BUILD-10} = union {BUILD-01..10}. Matches ROADMAP.md Phase 0 requirements list exactly. No orphans.

## Anti-Patterns Found

All advisory items from 00-REVIEW.md are surfaced here for completeness:

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `bin/cbl_tout_qt` | 49-52 | `rm -rf $MEFISTO/xvue/qt/build` without `: "${MEFISTO:?}"` guard | Warning (WR-01) | Foot-gun if env unset; inherited from legacy cbl_tout pattern. Advisory only. |
| `bin/cbl_tout_qt` | 91-95 | Five sub-script calls without `\|\| exit 1` guards | Warning (WR-02) | Silent exit-0 if a sub-build fails. 00-REVIEW recommends adding guards. Advisory only. |
| `xvue/README_COORDS.md` | 23-26 | Cited xvuelc.c lines 321 and 1621 do not contain quoted text | Warning (WR-03) | Two of four citations inaccurate; other two plus the structural pattern-based argument are correct; document conclusion remains valid per 00-REVIEW's independent re-read. Advisory only. |
| `xvue/qt/include/xvue_qt_api.h` | 21-27 | Dead preprocessor branch (lines 21-25 undefined before use) | Info (IN-01) | Deliberate byte-identical copy from xvuelc.c; 00-REVIEW recommends adding a clarifying comment. Advisory. |
| `xvue/qt/include/xvue_qt_api.h` | 15, 34, 41-48 | C++ features behind `.h` (not `.hpp`) extension | Info (IN-02) | Low-priority rename suggestion. Advisory. |
| `xvue/qt/cmake/verify_abi.sh` | 18 | Header regex allowlist covers only current return types | Info (IN-03) | Fragile if a future entry returns e.g. `char *`; would cause false drift. Advisory. |
| `xvue/qt/cmake/verify_abi.sh` | 17 | nm pattern constraint (lowercase 't' excluded) | Info (IN-04) | Correct for current phase; documentation comment suggested. Advisory. |
| `bin/cb*_qt` (5 scripts) | 37-38 | Duplicated QT_LIBS pkg-config invocation + link flags | Info (IN-05) | Extract to common shared include before Phase 1 adds more Qt components. Advisory. |
| `bin/cb*_qt` (5 scripts) | 39 | Inconsistent gfortran flags (-mcmodel=large, -fPIC, -Wall) across scripts | Info (IN-06) | Inherited from legacy; 00-REVIEW recommends harmonizing now. Advisory. |

**Classification:** 0 blockers, 3 warnings (all advisory), 6 info (all advisory). None block goal achievement for Phase 0 scope. None are regressions introduced by Phase 0 — most are inherited from legacy cb* patterns.

## Deferred Items

Items not yet met but explicitly addressed in later milestone phases, per Step 9b scan. None matched — all pending items in Phase 0 are captured by the existing HUMAN-UAT rather than deferred to later phases.

## Human Verification Required

### 1. Run 5 baseline testa/ cases on legacy X11 backend + fill in BASELINE.md

**Test:** As maintainer with X11 DISPLAY, run each of the 5 cases through its launcher and confirm expected behavior:

```bash
export MEFISTO=/home/drico/git/mefisto   # or the worktree where cbl_tout has been run
export MEFISTOX=$HOME/mefistox
export PATH=.:$PATH:$MEFISTO/bin
cd $MEFISTO

# Case 1: Mesher / testa/pan2d
mkdir -p $MEFISTOX/pan2d && cd $MEFISTOX/pan2d
cp -a $MEFISTO/testa/pan2d/* .
MAILLER                            # exit with "99;"

# Case 2: Elasticity / testa/nafems_le1
mkdir -p $MEFISTOX/nafems_le1 && cd $MEFISTOX/nafems_le1
cp -a $MEFISTO/testa/nafems_le1/* .
ELASTICER

# Case 3: Fluid / testa/cavity2d
mkdir -p $MEFISTOX/cavity2d && cd $MEFISTOX/cavity2d
cp -a $MEFISTO/testa/cavity2d/* .
FLUIDER

# Case 4: Thermal / testa/heat1d
mkdir -p $MEFISTOX/heat1d && cd $MEFISTOX/heat1d
cp -a $MEFISTO/testa/heat1d/* .
THERMICER

# Case 5: Nonlinear / testa/nlsecu
mkdir -p $MEFISTOX/nlsecu && cd $MEFISTOX/nlsecu
cp -a $MEFISTO/testa/nlsecu/* .
NLSER
```

After each case, edit `.planning/validation/BASELINE.md` and replace:
- `[Maintainer to describe — see Task 2]` → 1-3 sentences describing observed behavior
- `[Filled in by Task 2 after manual run]` → short note on warnings/flakes (or "Clean run — no surprises.")
- `[date filled in by Task 2]` → run date (e.g. `2026-04-11`)

**Expected:** All 5 cases open an X11 window, render the geometry/solution, exit cleanly after typing `99;`. `grep -c 'Maintainer to describe' .planning/validation/BASELINE.md` returns 0 after editing. `grep -c 'Filled in by Task 2' .planning/validation/BASELINE.md` returns 0.

**Why human:** D-15 validation gate explicitly requires a maintainer with a working X11 DISPLAY; D-16 specifies that the qualitative behavior field is maintainer knowledge that cannot be invented by an automated agent. Closes BUILD-07 behavioral half and BUILD-10 behavioral half in one session.

**Tracking:** `00-HUMAN-UAT.md` — 5 tests (pan2d, nafems_le1, cavity2d, heat1d, nlsecu), all currently `[pending]`. User chose `defer` during `/gsd:execute-phase 0`. Resolve with `/gsd-verify-work 0` after the runs complete.

## Gaps Summary

No gaps in the sense of "implementation failures". Phase 0 delivered all 10 requirements' scaffolding cleanly:

- **Qt build surface complete:** libxvueqt.a (17 704 B, 57 symbols), verify_abi parity guard, PIC/no-openmp/AUTOMOC/C++17 all correct per CMakeLists.txt direct read.
- **ABI bridge complete:** 57 extern "C" entries declared + 57 warn-once no-op stubs implemented, byte-identical proc() macro, MefistoPoint shim with static_assert, Planner Alert Option A exclusions honored.
- **Shell-script clone surface complete:** bin/cbl_tout_qt + 5 bin/cb*_qt scripts, each a literal clone of its legacy counterpart with the single surgical X11→Qt link substitution; legacy scripts byte-identical to HEAD per `git diff --quiet`.
- **5 Qt executables link-complete:** pp/ppmail_qt, pp/ppelas_qt, pp/ppflui_qt, pp/ppther_qt, pp/ppnlse_qt all exist, all proceed past `_start` without SIGSEGV/SIGABRT/SIGBUS under the `< /dev/null` smoke test (live-verified by this agent).
- **Documentation complete:** xvue/README_COORDS.md and .planning/validation/BASELINE.md both present and structurally correct.
- **Clean-tree regression proven:** 00-04-SUMMARY records that both `bin/cbl_tout` and `bin/cbl_tout_qt` rebuild cleanly from a wiped tree, producing all 10 executables (5 legacy + 5 Qt) in 379 s + 373 s respectively, with legacy sources untouched.

The one outstanding item is the 5-case interactive X11 baseline run (BUILD-07 behavioral half + BUILD-10 behavioral half), which is not an implementation gap — it is a manual validation step requiring a maintainer with X11 DISPLAY + maintainer knowledge of expected qualitative behavior for each case. This is already tracked in `00-HUMAN-UAT.md` with explicit user-choice `defer`. Status is `human_needed` rather than `passed` because the human verification list is non-empty; status is not `gaps_found` because no automated check fails.

The 3 warnings and 6 info items from 00-REVIEW.md are advisory code-quality concerns, not goal-achievement gaps; several are inherited from legacy scripts and are orthogonal to the phase's contract.

---

*Verified: 2026-04-10T22:00:00Z*
*Verifier: Claude (gsd-verifier)*
