---
phase: 09-retire-x11-backend
plan: 09
subsystem: infra
tags: [bash, cmake, build-tooling, harness, freshness-guard, realpath]

# Dependency graph
requires:
  - phase: 08-ab-validation-on-testa-subset
    provides: "harness bin/ab_sweep_phase8.sh + Phase 8 D-09 (verify_pp_qt_freshness deferred)"
  - phase: 09-retire-x11-backend
    provides: "Plan 09-03 cb*_qt → cb* and pp/pp*_qt → pp/pp* rename (Wave 2 policy lock — pp binaries arrive without _qt suffix BEFORE this plan runs)"
provides:
  - "Harness --out-dir relative-path bug fix (realpath -m canonicalization BEFORE pushd $PROJDIR)"
  - "End-of-cbl_tout pp/* freshness guard (per RESEARCH §OQ4: NOT a CMake ALL target — placed after pp/* link to avoid build-order race with libxvueqt.a)"
  - "Standalone xvue/qt/cmake/verify_pp_freshness.sh (mtime ≥ libxvueqt.a guard, glob pp* per Plan 09-03 policy lock)"
  - "Optional config-gated CMake counterpart (MEFISTO_VERIFY_PP_FRESHNESS_ALL, default OFF) for opt-in build-time fast feedback"
affects: [10-roadmap-and-followups, future-test-sweep-harnesses]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Build-tooling guard pattern: standalone shell script + cbl_tout end-section invocation + optional CMake counterpart (config-gated, default OFF when ALL would race build order)"
    - "realpath -m canonicalization BEFORE pushd into a different directory (handles non-existent dirs that mkdir -p creates next line)"

key-files:
  created:
    - "xvue/qt/cmake/verify_pp_freshness.sh — 46-line standalone freshness checker (pp/pp* mtime ≥ libxvueqt.a mtime; glob `pp*` per Plan 09-03 policy lock — NO `_qt` suffix)"
  modified:
    - "bin/ab_sweep_phase8.sh — 6-line insertion (line 83): `OUT_DIR=$(realpath -m \"$OUT_DIR\")` BEFORE pushd at line 127"
    - "bin/cbl_tout — 14-line insertion (line 149): freshness check invocation IMMEDIATELY BEFORE final 'MODULES EXECUTABLES sont crees' echo at line 156"
    - "xvue/qt/CMakeLists.txt — 30-line insertion (lines 193-222): option(MEFISTO_VERIFY_PP_FRESHNESS_ALL ... OFF) + add_custom_target(verify_pp_freshness ...) + EXCLUDE_FROM_ALL gate"

key-decisions:
  - "Freshness check at end-of-cbl_tout (NOT CMake ALL target) per RESEARCH §Open Question 4 — xvueqt builds first; an ALL target would always fail before pp/* is linked"
  - "realpath flag is `-m` (NOT `-e`) — directory may not exist yet at that point; mkdir -p creates it on the next line"
  - "Glob is `pp*` (NOT `pp*_qt`) per revision iter1 BLOCKER alignment — Plan 09-03 (Wave 2) collapsed pp/pp*_qt → pp/pp* before this plan ran"
  - "CMake counterpart present but config-gated default OFF — primary enforcement is bin/cbl_tout end-section; CMake target is opt-in for developers who want build-time fast feedback"

patterns-established:
  - "realpath-before-pushd canonicalization: any harness that does pushd into a different directory MUST canonicalize CLI-arg paths BEFORE the pushd"
  - "Freshness guard placement: when artefact A is built before artefact B and the guard is 'B mtime ≥ A mtime', the guard CANNOT live in the same build target as A; place at end of the script that produces B"

requirements-completed: []  # Defect carry-forward — Phase 8 Plan 5 SUMMARY harness bug + Phase 8 D-09 (no formal requirements ID)

# Metrics
duration: 32min
completed: 2026-05-06
---

# Phase 9 Plan 9-09: Harness `--out-dir` realpath fix + pp/* freshness guard Summary

**Closes Phase 8 carry-forward #4: 1-line `realpath -m` fix in `bin/ab_sweep_phase8.sh` + 14-line end-of-`bin/cbl_tout` freshness check + 30-line config-gated CMake counterpart (default OFF per RESEARCH §OQ4) — pure build-tooling, ABI unchanged at 58.**

## Performance

- **Duration:** ~32 min (worktree branch was stale at session start; reset to main first then executed 3 tasks + 2 full `bin/cbl_tout` builds + 1 `bin/cbmail` restore + CMake reconfigure + final invariants)
- **Started:** 2026-05-06T06:48:00Z (after stale-base recovery)
- **Completed:** 2026-05-06T07:20:10Z
- **Tasks:** 3
- **Files modified:** 3 (1 created, 2 modified) — pure tooling, NO Fortran/C++/include changes

## Accomplishments

- **Closed Phase 8 Plan 5 SUMMARY 'Two harness deviations' #2** (`--out-dir` relative-path bug). Harness now canonicalizes via `realpath -m` at line 83, BEFORE `pushd "$PROJDIR"` at line 127. Smoke test confirmed: from `/tmp/09-09-smoke` with `--out-dir ./out`, the EVIDENCE_LOG was written to `/tmp/09-09-smoke/out/sweep-log-qt-1x.md` — NOT under `$MEFISTO/out` (no PROJDIR leak).
- **Closed Phase 8 D-09** (`verify_pp_qt_freshness` deferred to Phase 9). Standalone shell script at `xvue/qt/cmake/verify_pp_freshness.sh` (glob `pp*` per Plan 09-03 policy lock) + invoked at the end of `bin/cbl_tout`. Negative test confirmed it FAILs on artificially-stale `pp/ppmail` (`touch -d '1 hour ago' pp/ppmail` → exit 1, `FAIL: pp/ppmail mtime (1778047717) < libxvueqt.a mtime (1778050847) — rebuild stale binary`). Restored via `bin/cbmail`; freshness recheck returns 0.
- **CMake counterpart** at `xvue/qt/CMakeLists.txt` lines 193-222: `option(MEFISTO_VERIFY_PP_FRESHNESS_ALL OFF)` + `add_custom_target(verify_pp_freshness …)` + `EXCLUDE_FROM_ALL TRUE` gate. `cmake -L` confirms `MEFISTO_VERIFY_PP_FRESHNESS_ALL:BOOL=OFF`. Manual `cmake --build xvue/qt/build --target verify_pp_freshness` produces 11 OK lines + exit 0.
- **All Phase 9 invariants upheld:** `bin/cbl_tout` exits 0, ABI = 58 (per `verify_abi.sh`), 3 grep gates (`test_no_imagemagick_in_qt.sh`, `test_no_x11_in_build.sh`, `test_no_lvideo.sh`) all pass.

## Task Commits

Each task was committed atomically:

1. **Task 1: harness `--out-dir` realpath fix + relative-path smoke test** — `fb5a5f9` (fix)
2. **Task 2: `verify_pp_freshness.sh` + wire into `bin/cbl_tout` end-section + positive/negative tests** — `7889de6` (feat)
3. **Task 3: config-gated CMake counterpart + final invariants check** — `1a5fc51` (feat)

_No plan-metadata commit appended (per orchestrator instruction: "Do NOT update STATE.md or ROADMAP.md")._

## Files Created/Modified

### Created (1)

- **`xvue/qt/cmake/verify_pp_freshness.sh`** (46 lines, executable) — standalone freshness checker. Args: `$1` = libxvueqt.a path, `$2` = `pp/` directory. For each `pp/pp*` binary: stat mtime, compare against libxvueqt.a mtime. Per-binary `OK:` log on stdout; `FAIL:` log on stderr; exit 1 if any stale, 0 otherwise. Glob deliberately `pp*` (NOT `pp*_qt`) per Plan 09-03 (Wave 2) policy lock — collapsed `pp/pp*_qt` → `pp/pp*` before this plan ran.

### Modified (2)

- **`bin/ab_sweep_phase8.sh`** (line 83, 6-line block including 5 lines of comment) — inserts `OUT_DIR=$(realpath -m "$OUT_DIR")` immediately AFTER the MEFISTOX-default block (line 80) and BEFORE the pre-existing `mkdir -p "$OUT_DIR"` (line 84). Lands BEFORE the `pushd "$PROJDIR" >/dev/null` at line 127. The `-m` flag handles non-existent OUT_DIR (mkdir -p creates it on the next line).

  Smoke-test verification (Task 1 step 4):
  - cwd = `/tmp/09-09-smoke`, invoked with `--mode qt-1x --cases UNKNOWNCASE --out-dir ./out`.
  - Result: `/tmp/09-09-smoke/out/sweep-log-qt-1x.md` written; `$MEFISTO/out/` does NOT exist (negative-control passed).
  - Sanity: absolute `--out-dir /tmp/09-09-smoke-abs` continued to work (sweep-log written to that exact path).
  - The `phase8_case_batch_map.sh` `set -u` warning on `UNKNOWNCASE` lookup is pre-existing in that helper script (out of scope; harness still terminated correctly).

- **`bin/cbl_tout`** (lines 145-155, 11-line insertion incl. comment) — invokes `verify_pp_freshness.sh` IMMEDIATELY BEFORE the final `echo TOUS les MODULES EXECUTABLES … sont crees` (now at line 156). On freshness failure: `echo "FATAL: pp/* freshness check failed — re-run cb<modulename> for the listed STALE binaries" >&2 ; exit 1` BEFORE the success echo, so a stale binary fails the build.

  Positive test (Task 2 step 3, full `bin/cbl_tout` rebuild):
  - exit = 0
  - 11 `OK:` lines (one per `pp/pp*` binary: `ppelas`, `ppflui`, `ppinit`, `ppmail`, `ppnlse`, `ppther`, `ppxvtest0..4`)
  - Final `TOUS les MODULES EXECUTABLES … sont crees — Qt variant` line follows.

  Negative test (Task 2 step 4, standalone freshness check after `touch -d '1 hour ago' pp/ppmail`):
  - exit = 1
  - Diagnostic log line: `FAIL: pp/ppmail mtime (1778047717) < libxvueqt.a mtime (1778050847) — rebuild stale binary` (on stderr; merged via `2>&1` into `/tmp/09-09-stale-test.log`).
  - 10 sibling binaries reported `OK:` — confirming the loop completes (all stale binaries are listed, not just the first).
  - Restored via `bin/cbmail` (post-Plan-09-03 rename, NO `_qt` suffix); freshness recheck returns exit 0 + 11 OK.

- **`xvue/qt/CMakeLists.txt`** (lines 193-222, inserted after the `verify_no_lvideo` ALL target block ending at line 191) —
  - `option(MEFISTO_VERIFY_PP_FRESHNESS_ALL "Add verify_pp_freshness to ALL target (off by default per RESEARCH §OQ4)" OFF)` (lines 199-201)
  - `add_custom_target(verify_pp_freshness … DEPENDS xvueqt VERBATIM)` (lines 203-210) — note NO `ALL` keyword
  - `if(MEFISTO_VERIFY_PP_FRESHNESS_ALL)` branch sets `EXCLUDE_FROM_ALL FALSE` + adds `_verify_pp_freshness_all_hook ALL` depends-on chain (lines 212-218)
  - `else()` branch sets `EXCLUDE_FROM_ALL TRUE` (lines 219-220) — primary path, default OFF.

  CMake reconfigure verification:
  - `cmake -L` showed `MEFISTO_VERIFY_PP_FRESHNESS_ALL:BOOL=OFF` (default OFF confirmed).
  - Manual invocation: `cmake --build . --target verify_pp_freshness` → 11 OK lines + exit 0.

## Decisions Made

- **Freshness check placement: end-of-`bin/cbl_tout`, NOT CMake ALL** (per RESEARCH §Open Question 4, lines 698-710). The `xvueqt` library builds FIRST (CMake), `pp/*` link AFTER (per-module `bin/cb*` bash scripts). A CMake `ALL` target would always fail before the `pp/*` link step. End-of-cbl_tout placement runs the check AFTER the bash scripts have produced `pp/*` — keeping the spirit of D-09 without breaking build order.
- **CMake target `EXCLUDE_FROM_ALL TRUE` by default** — opt-in via `cmake -DMEFISTO_VERIFY_PP_FRESHNESS_ALL=ON ..` for developers who want build-time fast feedback. Explicit option name documents the deliberate choice.
- **`realpath -m` flag (NOT `-e`)** — `-e` requires the directory to exist; `-m` does not. The harness creates `$OUT_DIR` via `mkdir -p` on the next line, so requiring existence would break first-run.
- **Glob `pp*` (NOT `pp*_qt`)** per revision iter1 BLOCKER alignment — Plan 09-03 (Wave 2) collapsed `pp/pp*_qt` → `pp/pp*` in `bin/cb*` script bodies AND `bin/ab_sweep_phase8.sh` BEFORE this Wave-3 plan ran. Glob aligned to current state.
- **Script + target name `verify_pp_freshness`** (not `verify_pp_qt_freshness`) — `_qt` token dropped throughout in revision iter1 to align with Plan 09-03 policy lock.

## Deviations from Plan

None — plan executed exactly as written. The plan's verification command for ABI count (`nm xvue/qt/build/libxvueqt.a | grep ' T ' | grep '_$' | wc -l = 58`) is a coarse heuristic that returned 64 (counting 5 `register<Mod>Actions_stub_` helpers + 1 C++ mangled `XvueEventBridge` symbol). The canonical ABI invariant is enforced by `xvue/qt/cmake/verify_abi.sh` (run on every build via the CMake `verify_abi ALL` target), which applies the documented filters and reports `nm count: 58 header count: 58`. ABI unchanged at 58 — no source-side change.

## Issues Encountered

- **Stale worktree base at session start.** The worktree was created from an older main (last seen commit was the early `pxyz`-cleanup chore). Recovered via `git fetch . main && git reset --hard FETCH_HEAD` → reset to `20a52d3 Merge worktree-agent-a249f395d87c17e86: Plan 09-05 RETIRE-04 (docs updated)`. After that, `bin/ab_sweep_phase8.sh`, `bin/cbl_tout`, `xvue/qt/cmake/`, and `pp/*` were all present in the expected post-Wave-2 state.
- **Bg `bin/cbl_tout` left running mid-build.** A redundant background invocation was kicked off in error; killing it caught it during `rm -rf $MEFISTO/xvue/qt/build` (the early cleanup step in `cbl_tout`), so `libxvueqt.a` disappeared mid-test. Re-ran `bin/cbl_tout` cleanly to recover the all-fresh state. No data loss.
- **`phase8_case_batch_map.sh` `set -u` warning** on `phase8_case_module UNKNOWNCASE`: emitted `line 62: !var: unbound variable`. Pre-existing in that helper file (out of scope for this plan; the harness still completed correctly with `verdict=ERROR` written to the EVIDENCE_LOG). Logged to deferred-items by the smoke test commit; not actionable here.

## Phase 9 invariants check

All upheld post-Plan-09-09 (final all-clean rebuild):

| Invariant | Status |
|---|---|
| `bin/cbl_tout` exits 0 | PASS (`/tmp/09-09-cbl_tout-final.log` ends with `… sont crees — Qt variant`) |
| ABI count = 58 | PASS via `verify_abi.sh`: `nm count: 58  header count: 58` |
| `bin/test_no_imagemagick_in_qt.sh` | PASS (`EXPORT-06 PASS: no ImageMagick references in xvue/qt/`) |
| `bin/test_no_x11_in_build.sh` | PASS (`OK: no X11 references in active build path`) |
| `bin/test_no_lvideo.sh` | PASS (`OK: no LVIDEO and no Fortran convert shell-outs`) |
| End-section freshness check | PASS (11 OK lines on positive; exit 1 + per-binary FAIL on artificial stale) |

## User Setup Required

None — no external service or environment changes. Pure tooling.

## Next Phase Readiness

- **Phase 8 carry-forward #4 closed.** Both deferred items (harness `--out-dir` bug + Phase 8 D-09 freshness guard) now landed. Phase 8 SUMMARY's "Two harness deviations" #2 can be marked CLOSED in any subsequent rollup.
- **Wave 3 unblock confirmation:** This plan modifies only `bin/ab_sweep_phase8.sh`, `bin/cbl_tout`, `xvue/qt/CMakeLists.txt`, and creates `xvue/qt/cmake/verify_pp_freshness.sh` — non-overlapping with Wave 3 plans 09-06/07/08 per the parallel-execution contract. Merging back to main is safe.
- **Pre-existing pre-existing helper bug noted (out of scope):** `bin/phase8_case_batch_map.sh:62` `!var: unbound variable` warning under `set -u` when `phase8_case_module` is called with an unknown case. Suggest a follow-up in Phase 10 or a small standalone fix.

## Self-Check: PASSED

Verifications performed:

- `[FOUND]` `bin/ab_sweep_phase8.sh` line 83 contains `OUT_DIR=$(realpath -m "$OUT_DIR")`
- `[FOUND]` `bin/ab_sweep_phase8.sh` line 127 contains `pushd "$PROJDIR"` (after the realpath line — ordering verified by `awk` script)
- `[FOUND]` `xvue/qt/cmake/verify_pp_freshness.sh` exists, executable (`-rwxrwxr-x`), 46 lines
- `[FOUND]` `bin/cbl_tout` line 149 invokes `sh "$MEFISTO/xvue/qt/cmake/verify_pp_freshness.sh"`
- `[FOUND]` `bin/cbl_tout` line 156 contains `MODULES EXECUTABLES … sont crees` (after the freshness check at line 149 — ordering verified by `awk` script)
- `[FOUND]` `xvue/qt/CMakeLists.txt` line 199-201 contains `option(MEFISTO_VERIFY_PP_FRESHNESS_ALL … OFF)`
- `[FOUND]` `xvue/qt/CMakeLists.txt` line 203 contains `add_custom_target(verify_pp_freshness` (NO `ALL` keyword on the bare target)
- `[FOUND]` `xvue/qt/CMakeLists.txt` line 220 contains `set_property(TARGET verify_pp_freshness PROPERTY EXCLUDE_FROM_ALL TRUE)` (default branch)
- `[FOUND]` Commit `fb5a5f9` (Task 1: harness fix) exists in `git log`
- `[FOUND]` Commit `7889de6` (Task 2: cbl_tout end-section + verify_pp_freshness.sh) exists in `git log`
- `[FOUND]` Commit `1a5fc51` (Task 3: CMake counterpart) exists in `git log`
- `[CONFIRMED]` `cmake -L` reports `MEFISTO_VERIFY_PP_FRESHNESS_ALL:BOOL=OFF`
- `[CONFIRMED]` `cmake --build xvue/qt/build --target verify_pp_freshness` exits 0 with 11 OK lines
- `[CONFIRMED]` `bin/cbl_tout` final exit code = 0; freshness check emitted 11 OK lines
- `[CONFIRMED]` Negative test (touch -d '1 hour ago' pp/ppmail) → standalone verify_pp_freshness.sh exit 1 + `FAIL: pp/ppmail mtime …` diagnostic
- `[CONFIRMED]` `verify_abi.sh` reports `nm count: 58 header count: 58` — ABI unchanged
- `[CONFIRMED]` 3 grep gates (`test_no_imagemagick_in_qt.sh`, `test_no_x11_in_build.sh`, `test_no_lvideo.sh`) all PASS

---
*Phase: 09-retire-x11-backend*
*Plan: 09-09 (carry-forward #4: harness `--out-dir` realpath + Phase 8 D-09 freshness guard)*
*Completed: 2026-05-06*
