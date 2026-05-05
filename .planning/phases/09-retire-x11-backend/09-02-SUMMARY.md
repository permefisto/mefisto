---
phase: 09-retire-x11-backend
plan: 02
subsystem: build/graphics-backend
tags: [retire-x11, xlib, c-deletion, shell-script-edit, abi-invariant, qt-only-build]

# Dependency graph
requires:
  - phase: 09-retire-x11-backend
    provides: "v1.0-pre-retire tag (rollback artifact); 09-01-AUDIT-BASELINE.md (empirical reference counts); A/B window closed 2026-05-06"
provides:
  - "xvue/xvuelc.c (3749-line Xlib backend) deleted from working tree and git index"
  - "xvue/xvuelc.o (228 KB tracked binary) deleted from working tree and git index"
  - "bin/ccxvue (42-line C compile wrapper) deleted from working tree and git index"
  - "bin/cbl_tout no longer invokes ccxvue (3-line block removed)"
  - "Empirical proof that no Fortran INCLUDE/COMMON references xvuelc symbols (T-09-01 mitigated)"
  - "Empirical proof that bin/cbl_tout_qt + 4 testa cases (pan2d, nafems_le1, cavity2d, heat1d) survive deletion"
affects: [09-03-retire-x11-build-scripts, 09-04-retire-imagemagick-lvideo, 09-05-retire-docs-deps]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Mechanical-deletion-with-atomic-commit (RESEARCH §Pattern 1)"
    - "Pre-deletion ABI gate via verify_abi.sh (58 Fortran-facing symbols stable)"

key-files:
  created:
    - ".planning/phases/09-retire-x11-backend/09-02-SUMMARY.md"
  modified:
    - "bin/cbl_tout (removed lines 48-50: comment + ccxvue invocation)"
  deleted:
    - "xvue/xvuelc.c (3749 lines, 156 KB)"
    - "xvue/xvuelc.o (228 KB tracked build artifact)"
    - "bin/ccxvue (42 lines, 1.1 KB)"

key-decisions:
  - "Used `git rm xvue/xvuelc.o` (Rule 1 auto-fix: plan claimed .o was gitignored but it was tracked)"
  - "Re-ran ab_sweep_phase8.sh with --smoke-only (Rule 3 auto-fix: X11 baselines no longer exist post-deletion; same fix as 09-01)"
  - "Filtered ABI count (58, via verify_abi.sh logic) is the operative invariant; raw nm count (64) is a superset including 6 register*Actions_stub_ internal helpers (matches 09-01-AUDIT-BASELINE.md §3 row #7 with both numbers explained)"

patterns-established:
  - "Phase 9 wave 2 commit-message convention: feat(09-NN): RETIRE-MM — <short title>; body documents pre/post invariants and explicitly notes which legacy paths are temporarily broken pending downstream plans"
  - "ABI invariant assertion before AND after wholesale deletion (verify_abi.sh exit 0 both sides)"

requirements-completed: ["RETIRE-01"]

# Metrics
duration: ~10 min
completed: 2026-05-06
---

# Phase 9 Plan 02: RETIRE-01 X11 Backend Deletion Summary

**Deleted 3749-line `xvue/xvuelc.c` Xlib backend, 228 KB `xvuelc.o` tracked build artifact, and 42-line `bin/ccxvue` C compile wrapper; removed the corresponding 3-line `cbl_tout` invocation block; verified zero Fortran ghost references and stable 58-symbol Qt ABI on both sides of the deletion.**

## Performance

- **Duration:** ~10 min
- **Started:** 2026-05-06T00:21:00Z (approx)
- **Completed:** 2026-05-06T00:32:36Z (commit timestamp)
- **Tasks:** 3 (all autonomous)
- **Files modified:** 1 (bin/cbl_tout)
- **Files deleted:** 3 (xvuelc.c, xvuelc.o, ccxvue)

## Accomplishments

- **3749-line Xlib backend retired.** The single largest C file in the repository (`xvue/xvuelc.c`, 156 KB) is gone from the working tree and the git index.
- **C compile wrapper retired.** `bin/ccxvue` (the only `cb*` script that compiled `xvuelc.c`) is deleted; `bin/cbl_tout` no longer invokes it.
- **Tracked build artifact removed.** `xvue/xvuelc.o` (228 KB) was incorrectly checked into git in a prior era; this plan removes it from the index.
- **T-09-01 mitigated empirically.** `grep -rln 'xvuelc' xvue/*.f xvue/*.inc incl/*.inc` returns 0 — no Fortran source references the deleted C file by name in any active path.
- **ABI invariant held.** `verify_abi.sh` exit 0 both pre-deletion and post-deletion: 58 Fortran-facing T-symbols ending in `_` (excluding the 6 `register*Actions_stub_` internal dispatch helpers that the audit baseline counted unfiltered as 64).
- **Qt build green.** `bin/cbl_tout_qt` exit 0 post-deletion — all `pp*_qt` executables linked successfully.
- **4/5 testa baselines captured.** `pan2d`, `nafems_le1`, `cavity2d`, `heat1d` produced PNGs via the Phase 8 harness in smoke-only mode (`nlsecu` carries forward to Plan 9-06 per Phase 8 override #5 — ppnlse_qt offscreen deadlock).

## Task Commits

Each task contributed to a single atomic commit per the plan's "atomic deletion + edit (one logical commit)" instruction (Task 2 step 1-5 staged everything; Task 3 step 5 committed):

1. **Task 1: Pre-deletion ABI gate** — read-only verification, no commit
2. **Task 2: Deletion + cbl_tout edit (staged)** — folded into single commit
3. **Task 3: Post-deletion build/sweep + commit** — `769b54a` (feat)

**The single commit:**

| Hash | Type | Subject |
| --- | --- | --- |
| `769b54a` | feat | feat(09-02): RETIRE-01 — delete xvuelc.c + ccxvue + remove from cbl_tout |

A trailing metadata commit will be made for this SUMMARY itself.

## Files Created/Modified

- **Created:** `.planning/phases/09-retire-x11-backend/09-02-SUMMARY.md` (this file)
- **Modified:** `bin/cbl_tout` (3-line block removed: comment, blank, `$MEFISTO/bin/ccxvue` invocation)
- **Deleted:** `xvue/xvuelc.c` (3749 lines), `xvue/xvuelc.o` (228 KB), `bin/ccxvue` (42 lines)

## Empirical Counts vs. 09-01-AUDIT-BASELINE.md §3

| # | Reference | Pre (Plan 9-01) | Post (Plan 9-02) | Expected | Status |
|---|-----------|-----------------|------------------|----------|--------|
| 1 | `wc -l xvue/xvuelc.c` | 3749 | file does not exist | 0 (file gone) | ✓ |
| 7 | ABI raw count `nm \| grep ' T ' \| grep '_$' \| wc -l` | 64 | 64 | 64 unchanged | ✓ |
| — | ABI filtered count (verify_abi.sh) | 58 | 58 | 58 unchanged | ✓ |
| — | `grep -c ccxvue bin/cbl_tout` | 1 | 0 | 0 | ✓ |
| — | `grep -rln 'xvuelc' xvue/*.f xvue/*.inc incl/*.inc \| wc -l` | 0 | 0 | 0 | ✓ |

The remaining baseline rows (#2, #3, #4, #5, #6, #8, #9) are out of scope for this plan — they belong to Plans 9-03 (X11 build scripts), 9-04 (ImageMagick/LVIDEO), and 9-05 (docs).

## Decisions Made

- **Single atomic commit for Tasks 2 + 3** — the plan text explicitly mandates "atomic deletion + edit (one logical commit)" in Task 2's `<action>` block, and Task 3 step 5 is the commit step for the same logical change. Following the plan as written.
- **`xvue/xvuelc.o` removed via `git rm`** rather than `rm -f` (Rule 1 auto-fix; see Deviations).
- **Sweep run in `--smoke-only` mode** (Rule 3 auto-fix; see Deviations).
- **ABI count expressed as both raw (64) and filtered (58)** in the commit message and this summary, to bridge the apparent discrepancy between the plan text (which uses 58) and the 09-01 audit baseline (which uses 64) — they are the same underlying invariant measured with two different filters.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug in plan instructions] `xvue/xvuelc.o` was tracked, not gitignored**
- **Found during:** Task 2 (Delete xvuelc.c + xvuelc.o + ccxvue)
- **Issue:** The plan's Task 2 step 1 says: "rm -f xvue/xvuelc.o (.o is .gitignored — do not git-rm; just rm)". Empirically false: `git ls-files xvue/xvuelc.o` returned the file, and `git check-ignore -v xvue/xvuelc.o` returned no ignore rule. The .o was tracked in git.
- **Fix:** Used `git rm xvue/xvuelc.o` to stage the deletion in the index alongside the other two files.
- **Files modified:** Index entry for `xvue/xvuelc.o` removed; working tree file removed.
- **Verification:** `! test -f xvue/xvuelc.o`, `git ls-files xvue/xvuelc.o` returns empty, `git show --stat HEAD` shows the deletion in the commit.
- **Committed in:** `769b54a`

**2. [Rule 3 - Blocking issue] `ab_sweep_phase8.sh` requires `--smoke-only` post-deletion**
- **Found during:** Task 3 (Verify Qt build green + Qt-path testa baselines)
- **Issue:** The plan's Task 3 step 2 invokes `bin/ab_sweep_phase8.sh --mode qt-1x --cases ... --out-dir /tmp/09-02-post-retire01` (no `--smoke-only` flag). The harness requires X11 baselines to exist for diff comparison. Plan 9-01 established Qt-1x as the post-A/B-window baseline — there are no X11 PNGs to compare against. First sweep returned `verdict=ERROR reason="baseline missing: ...-x11.png"` for all 5 cases (exit 0 but no useful output).
- **Fix:** Re-ran with `--smoke-only` (the same auto-fix Plan 9-01 applied for the same root cause; documented in 09-01-AUDIT-BASELINE.md §5 paragraph "the plan's Task 3 step 3 used the non-smoke mode by mistake; auto-fix Rule 3 applied").
- **Files modified:** None (test-time invocation only)
- **Verification:** Second sweep produced `verdict=SMOKE` for all 5 cases with PNGs at `/tmp/09-02-post-retire01/{pan2d,nafems_le1,cavity2d,heat1d}-qt-1x.png` (4-of-5 capture; nlsecu at 0 bytes per Phase 8 override #5 carry-forward to Plan 9-06).
- **Committed in:** N/A (test-time only)

---

**Total deviations:** 2 auto-fixed (1 plan-text bug, 1 blocking-issue carry-forward).
**Impact on plan:** Both auto-fixes were essential to complete the plan as designed. No scope creep — the underlying intent (delete the 3 files, edit cbl_tout, verify Qt path green, capture testa baselines) was satisfied verbatim.

## Issues Encountered

- **Plan-text vs. baseline ABI count discrepancy.** Plan 09-02 Task 1 step 2 asserts `nm xvue/qt/build/libxvueqt.a | grep ' T ' | grep '_$' | wc -l == 58`. Empirically, that command returns **64** on commit `53bdb5b` (the v1.0-pre-retire baseline) and remains **64** post-deletion. The 09-01 audit baseline §3 row #7 records 64 with explicit drift acknowledgement (Phase 7 XvueExport added 6 entry points). The reconciliation: the project's own gate, `xvue/qt/cmake/verify_abi.sh`, applies two additional filters (`grep -v ' T _Z'` for C++ mangled names, `grep -v ' T register*Actions_stub_$'` for internal dispatch helpers) and yields **58**, which equals the header-declared count of 58. Both numbers describe the same stable ABI surface; plan text uses the filtered semantic, baseline uses the raw command. Recorded both in commit message and summary so future readers don't re-discover the apparent contradiction.

## Threat Flags

None. All file deletions stay within the plan's `<threat_model>` scope (T-09-01 mitigated, T-09-02 / T-09-03 accepted by design). No new network endpoints, auth paths, file access patterns, or schema changes introduced.

## Next Phase Readiness

- **Plan 09-03 (RETIRE-02) unblocked.** The 27-29 `bin/cb*` scripts that still link `-lX11` (and the `bin/Makefile*` files that mention X11R6) can now be deleted/edited freely; the `xvuelc.o` they reference no longer exists, so any latent attempt to use legacy `bin/cbl_tout` will surface immediately as a link error (loud-failure mode, intentional).
- **Plans 09-04, 09-05 unaffected.** This plan touched neither LVIDEO/ImageMagick callers nor docs.
- **Phase 8 carry-forward unchanged.** `nlsecu` ppnlse_qt offscreen deadlock (Phase 8 override #5) still pending Plan 9-06; not introduced or worsened by this plan.

### Known limitations introduced (intentional, by plan design)

- `bin/cbl_tout` (legacy X11 build entry) now FAILS at the `cb{init,mail,elas,flui,ther,nlse,...}` link stage because those `cb*` scripts still reference `-lX11` and (in some cases) the now-deleted `xvuelc.o`. Plan 9-03 deletes those legacy scripts and renames the `cb*_qt` variants to drop the suffix.
- This is documented in the commit message body and is the explicit Phase 9 wave 2 contract.

## Self-Check: PASSED

- File `.planning/phases/09-retire-x11-backend/09-02-SUMMARY.md`: FOUND (this file, post-write)
- Commit `769b54a`: FOUND (`git log --oneline | grep 769b54a` returns the feat commit)
- Anti-artifact `xvue/xvuelc.c`: confirmed absent (`! test -f xvue/xvuelc.c` ✓)
- Anti-artifact `xvue/xvuelc.o`: confirmed absent (`! test -f xvue/xvuelc.o` ✓)
- Anti-artifact `bin/ccxvue`: confirmed absent (`! test -f bin/ccxvue` ✓)
- `git ls-files xvue/xvuelc.c bin/ccxvue xvue/xvuelc.o`: empty (deletions in index ✓)
- `grep -c 'ccxvue' bin/cbl_tout`: 0 ✓
- `grep -rln 'xvuelc' xvue/*.f xvue/*.inc incl/*.inc`: 0 matches ✓
- `verify_abi.sh xvue/qt/build/libxvueqt.a xvue/qt/include/xvue_qt_api.h`: exit 0, 58/58 ✓
- `bin/cbl_tout_qt`: exit 0 ✓
- 4 PNGs at `/tmp/09-02-post-retire01/{pan2d,nafems_le1,cavity2d,heat1d}-qt-1x.png`: present ✓

---
*Phase: 09-retire-x11-backend*
*Completed: 2026-05-06*
