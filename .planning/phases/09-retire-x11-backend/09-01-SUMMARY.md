---
phase: 09-retire-x11-backend
plan: 01
subsystem: infra
tags: [git-tag, rollback-safety, audit-baseline, x11-retirement, lvideo, abi]

# Dependency graph
requires:
  - phase: 08-ab-validation-on-testa-subset
    provides: ship-gate sign-off (5 overrides accepted), one-release-cycle A/B window opened, override #5 carry-forward (ppnlse_qt deadlock), Phase 8 sweep harness
provides:
  - "git tag v1.0-pre-retire (annotated) at commit 53bdb5b — single-command rollback for the entire phase"
  - "git branch retire-restore-point at the same commit"
  - "09-01-AUDIT-BASELINE.md with empirical pre-retirement reference counts (9 metrics) — drives Plan 9-02..9-05 self-validation"
  - "Independent confirmation that bin/cbl_tout exit 0 on the retirement tag commit (dual-backend build co-existence proven)"
  - "4/5 testa BUILD-10 baselines captured under Qt-1x smoke mode at the tag commit"
affects: [09-02, 09-03, 09-04, 09-05, 09-06, 09-07, 09-08, 09-09]

# Tech tracking
tech-stack:
  added: []
  patterns: [
    "Annotated-tag rollback contract pattern: tag + branch at same commit, dereference via ^{commit} for SHA equality",
    "Smoke-only sweep mode for self-baseline capture (no upstream X11 reference exists when capturing the rollback tag baseline)",
    "Empirical-vs-RESEARCH delta documentation in audit baselines (drift acknowledgement table)"
  ]

key-files:
  created: [
    ".planning/phases/09-retire-x11-backend/09-01-AUDIT-BASELINE.md (15.5 KB, 6 sections + verification log)"
  ]
  modified: []

key-decisions:
  - "Tag + branch creation routed through parent repo path (`git -C $MEFISTO`) so refs are visible from main immediately, satisfying the worktree cross-branch visibility requirement"
  - "Empirical truth supersedes RESEARCH heuristics in §3 of the audit baseline — 3 measured deltas (#7 ABI +6, #8 cb*-X11 −1, #9 README −2) recorded with explanations"
  - "--smoke-only sweep is the correct invocation for Plan 9-01 testa baselines; the plan-text non-smoke command would have failed on missing X11 baseline (auto-fix Rule 3 applied)"
  - "Origin push deferred to maintainer (SSH host-key gate in agent env); D-08 'pushed to origin' is best-effort, local rollback works independently"

patterns-established:
  - "Annotated-tag SHA caveat documentation: `git rev-parse <tag>` returns the tag-object SHA; semantic equality requires `^{commit}` dereference"
  - "Empirical-count audit baseline pattern with drift table: each row carries (Expected from RESEARCH, Empirical, Command) so Plans 9-02..9-05 can re-run identical commands and assert post-deletion delta == expected"

requirements-completed: []  # PLAN frontmatter listed `requirements: []` — no requirement IDs to mark complete

# Metrics
duration: ~13min
completed: 2026-05-05
---

# Phase 9 Plan 01: Pre-Phase-9 rollback safety + audit baseline Summary

**Annotated tag `v1.0-pre-retire` + branch `retire-restore-point` created on main HEAD `53bdb5b`; 9-metric empirical audit baseline committed; bin/cbl_tout exit 0 + 4/5 testa BUILD-10 baselines captured. Wave 2 (Plans 9-02..9-05) unblocked.**

## Performance

- **Duration:** ~13 min
- **Started:** 2026-05-05T22:01:54Z (worktree branch check + plan read)
- **Completed:** 2026-05-05T22:15:18Z (Task 5 commit + summary)
- **Tasks:** 5/5 complete (Task 1 was a checkpoint — D-11 process gate already satisfied per orchestrator brief and STATE.md line 134)
- **Files modified:** 2 (1 created: 09-01-AUDIT-BASELINE.md; 1 created via git: v1.0-pre-retire tag + retire-restore-point branch refs in parent repo)

## Accomplishments

- **D-08/D-09 rollback contract established:** annotated git tag `v1.0-pre-retire` (object `73d32e3`, target commit `53bdb5b`) + branch `retire-restore-point` at the same commit. Single-command revert: `git reset --hard v1.0-pre-retire`.
- **D-10 build-health sanity check:** `bin/cbl_tout` exit 0 on the tag commit (~6 min build, all 8 legacy executables `pp{init,mail,elas,flui,ther,nlse,poba,xvtest{0..4}}` AND 8 Qt counterparts `pp*_qt` link cleanly — last commit on which both backends co-exist).
- **D-11 process gate confirmed:** A/B window closure recorded in STATE.md line 134 (maintainer dricoco, 2026-05-06). Audit baseline §1 cites the closure source.
- **9-metric empirical audit baseline** committed to `.planning/phases/09-retire-x11-backend/09-01-AUDIT-BASELINE.md` with detailed file lists for each retirement work surface — exact starting values for Plans 9-02..9-05 to drive to zero/target.
- **4/5 BUILD-10 testa baselines** captured under Qt-1x smoke mode (pan2d 71KB, nafems_le1 320KB, cavity2d 44KB, heat1d 21KB; nlsecu 0B per Phase 8 override #5 deadlock — explicit Plan 9-06 carry-forward).

## Task Commits

Per orchestrator brief: **only one git commit lands in this plan** — the audit baseline file (Task 5). The tag + branch (Task 2) are non-commit refs created directly in the parent repo. No per-task commits for Tasks 1-4 because Tasks 1, 2, 3, 4 produce no tracked-file changes (Task 1 = process gate confirmation, Task 2 = git refs, Task 3 = `/tmp/` build/sweep logs, Task 4 = audit doc that lands as Task 5's commit).

1. **Task 1: Maintainer A/B window closure confirmation** — checkpoint pre-cleared by orchestrator brief and `.planning/STATE.md` line 134.
2. **Task 2: Create v1.0-pre-retire tag + retire-restore-point branch** — refs in parent repo (no commit produced):
   - `v1.0-pre-retire` (annotated) → `53bdb5b59ecbb6a7af210d1bde3ded7857d376c5`
   - `retire-restore-point` → `53bdb5b59ecbb6a7af210d1bde3ded7857d376c5`
3. **Task 3: Build + 5 testa baselines** — `/tmp/09-01-cbl_tout-pre.log` (exit 0), `/tmp/09-01-pre-retire-baseline/{pan2d,nafems_le1,cavity2d,heat1d}-qt-1x.png`, `/tmp/09-01-sweep-pre.log`. No tracked-file changes; sweep verdict line: `case=$X mode=qt-1x verdict=SMOKE size=$N` for all 5.
4. **Task 4: Write 09-01-AUDIT-BASELINE.md** — file written, staged in worktree.
5. **Task 5: Commit Plan 9-01 artifacts** — `933bdbb` (`docs(09-01): commit pre-retirement audit baseline (v1.0-pre-retire gate)`).

## Files Created/Modified

**Created:**
- `.planning/phases/09-retire-x11-backend/09-01-AUDIT-BASELINE.md` (15.5 KB) — 6 sections + verification log:
  - §1 Process gate confirmation (D-11 source citation)
  - §2 Tag + branch SHAs + annotated-tag deref caveat
  - §3 9 empirical reference counts + detailed file lists + drift acknowledgement
  - §4 Build verification (cbl_tout exit 0, ~6 min, 8+8 executables linked)
  - §5 testa baseline capture log (4/5 PNGs + nlsecu disposition)
  - §6 Carry-forward acknowledgment (Wave 2 expectations + Plans 9-06..9-09 reuse)

**Created (parent-repo refs, not commits):**
- `v1.0-pre-retire` (annotated tag) at `73d32e3ad483bed391317d4525d96c5ae94d3b30`
- `retire-restore-point` (branch) at `53bdb5b59ecbb6a7af210d1bde3ded7857d376c5`

**Modified:** None (per orchestrator constraint: STATE.md and ROADMAP.md not touched).

## Decisions Made

- **Tag/branch routed through parent repo (`git -C $MEFISTO`)** — the worktree branch is a per-agent feature branch, so creating tag/branch on the worktree HEAD would not be visible from main. Routing creation through the parent repo path makes the refs visible from main immediately, satisfying the worktree cross-branch visibility requirement called out in the orchestrator's worktree_branch_check.
- **Worktree HEAD recovery** — initial worktree HEAD was at stale commit `ac282f8` (an old commit, not main). Per worktree_branch_check guidance, recovered via `git fetch . main && git reset --hard FETCH_HEAD` to sync to current main `53bdb5b`. After recovery, HEAD assertion passes (worktree-agent-* namespace, not protected ref).
- **`--smoke-only` sweep mode** — the plan-text Task 3 step 3 invoked the harness without `--smoke-only`, which fails on missing X11 baseline (`reason="baseline missing: ${case}-x11.png"`). Auto-fix Rule 3 (blocking issue): for a self-baseline capture (Qt-only at the rollback tag), the X11 left-operand of comparison does not exist by definition; `--smoke-only` is the correct mode. Re-ran with the flag, all 5 cases produced expected verdict=SMOKE output (4/5 with non-zero PNGs, nlsecu zero per documented override #5).
- **Empirical-truth supersedes RESEARCH** — 3 of 9 reference counts (ABI, cb*-X11 union, README/LISEZMOI) drifted from RESEARCH-predicted values. Empirical values recorded as canonical; drift acknowledged with explanation (Phase 6.5/7 ABI growth; minor scope corrections for the other two). Plans 9-02..9-05 anchor against the empirical values.
- **Origin push deferred** — SSH host-key gate in agent env blocked `git push origin v1.0-pre-retire retire-restore-point`. D-08 wording is "pushed to origin if remote configured" (best-effort); local rollback works independently. Maintainer follow-up documented in §2 of the audit baseline.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 — Blocking] Sweep harness `--smoke-only` flag required for self-baseline mode**
- **Found during:** Task 3 (Build the v1.0-pre-retire commit + run 5 testa baselines)
- **Issue:** Plan-text Task 3 step 3 invoked `bin/ab_sweep_phase8.sh --mode qt-1x ... --out-dir /tmp/09-01-pre-retire-baseline` without `--smoke-only`. The harness defaults to A/B-comparing the qt-1x capture against `${OUT_DIR}/${CASE}-x11.png`, which does not exist (Plan 9-01 IS the baseline; no upstream X11 reference is available). All 5 cases failed with `verdict=ERROR reason="baseline missing"`, sweep exit code 4.
- **Fix:** Re-ran with `--smoke-only`. Per harness header docs (lines 24-28): "Skip the X11-baseline-compare step entirely. Used by Plan 1 Task 3 BEFORE any X11 baseline exists. Verdict emitted is `verdict=SMOKE`." This is exactly the Plan-9-01 use case.
- **Files modified:** None (only changed the invocation flag).
- **Verification:** Re-run produced all 5 expected verdict=SMOKE lines (sweep exit 0); 4/5 PNGs land at expected paths with non-zero sizes; nlsecu produces 0-byte file matching Phase 8 override #5.
- **Committed in:** documented in audit baseline §5; no separate commit needed (sweep produces /tmp/ artifacts, not tracked files).

**2. [Rule 3 — Blocking] Worktree HEAD recovery from stale base**
- **Found during:** pre-Task-2 worktree branch check
- **Issue:** Worktree was created at commit `ac282f8` (stale, predates current main `53bdb5b` by many phase). Tag creation here would land on the wrong commit and fail the rollback contract.
- **Fix:** Per worktree_branch_check guidance: `git fetch . main && git reset --hard FETCH_HEAD`. HEAD now matches main `53bdb5b`.
- **Files modified:** None (only worktree state changed).
- **Verification:** `git rev-parse HEAD` = `53bdb5b59ecbb6a7af210d1bde3ded7857d376c5` (matches `git -C $MEFISTO rev-parse main`); HEAD assertion passes (`worktree-agent-adb43afc8c0cdfdc1` is in the per-agent namespace, not a protected ref).
- **Committed in:** N/A (worktree state recovery, no tracked-file changes).

**3. [Rule 1 — Bug, documentation] Annotated-tag SHA equality test wording**
- **Found during:** Task 2 verify step
- **Issue:** Plan-text automated check `[ "$(git rev-parse v1.0-pre-retire)" = "$(git rev-parse retire-restore-point)" ]` returns false for annotated tags because `git rev-parse <tag>` returns the tag-object SHA (`73d32e3`), not the target commit. Naive interpretation would suggest the verify step failed.
- **Fix:** Documented in audit baseline §2: the semantically correct comparison uses `git rev-parse v1.0-pre-retire^{commit}` and matches `retire-restore-point` SHA exactly. Acceptance-criterion intent (tag and branch point at the same commit) is satisfied. Tag, branch, and target commit all anchor on `53bdb5b`.
- **Files modified:** 09-01-AUDIT-BASELINE.md §2 (caveat section).
- **Verification:** `git rev-parse v1.0-pre-retire^{commit}` = `git rev-parse retire-restore-point` = `53bdb5b59ecbb6a7af210d1bde3ded7857d376c5`. Equal.

---

**Total deviations:** 3 auto-fixed (3 Rule 3/Rule 1 — all blocking or documentation-correctness items, no architectural changes).
**Impact on plan:** All three fixes are essential for correctness and reproducibility. No scope creep — the underlying retirement contract (tag + branch + audit) is unchanged.

## Authentication Gates

**1. Origin push (Task 2 step 4)**
- **Source:** sandboxed agent env lacks GitHub SSH host key trust
- **Command attempted:** `git push origin v1.0-pre-retire` and `git push origin retire-restore-point`
- **Error:** `Host key verification failed. fatal: Could not read from remote repository.`
- **Resolution:** D-08 wording explicitly makes origin push best-effort ("pushed to origin if remote configured" — and the plan's Task 2 step 4 itself wraps the push in `if git remote get-url origin >/dev/null 2>&1; then ... ; else echo "No origin remote — skipping push"; fi`). Local refs satisfy the rollback contract. Maintainer follow-up: `git push origin v1.0-pre-retire retire-restore-point` from the main worktree (where SSH agent is set up).
- **Documented:** audit baseline §2 "Origin push status: DEFERRED" with maintainer-follow-up command.

## Issues Encountered

- **Initial worktree HEAD on stale commit `ac282f8`:** noted under Deviation #2 above. Recovery was clean.
- **Path mismatch between Write tool and worktree:** the Write tool created `09-01-AUDIT-BASELINE.md` at the parent-repo path `/home/mefisto/git/mefisto/.planning/...` (not the worktree path `/home/mefisto/git/mefisto/.claude/worktrees/agent-adb43afc8c0cdfdc1/.planning/...`). Resolved by `cp` from parent to worktree, then `git add` from worktree (the worktree shares git history with parent, so the commit lands correctly on the per-agent branch).

## User Setup Required

**Maintainer follow-up (single command):**

```bash
cd $MEFISTO
git push origin v1.0-pre-retire retire-restore-point
```

This publishes the rollback artifacts to GitHub. Not a blocker for Phase 9 execution — local rollback (`git reset --hard v1.0-pre-retire`) works without origin.

## Next Phase Readiness

- **Wave 2 (Plans 9-02..9-05) UNBLOCKED.** Each retirement plan can re-run the §3 grep commands and assert post-deletion delta vs the recorded baselines.
- **Plans 9-06..9-09 carry-forward** can proceed independently. Plan 9-06 inherits the documented nlsecu deadlock disposition; Plan 9-07 inherits the v1.0-pre-retire tag for cross-tag worktree creation (per CONTEXT.md D-03 + plan 9-07 procedure).
- **Rollback safety:** any of the 4 retirement plans (9-02..9-05) can be backed out individually (per-plan revert) or wholesale (`git reset --hard v1.0-pre-retire`) without losing this audit-baseline commit (it lives at `933bdbb` on the per-agent branch and will land on main when the worktree merges back).

## Self-Check: PASSED (re-verified post-commit)

```
created files: 09-01-AUDIT-BASELINE.md (15888 B), 09-01-SUMMARY.md (15922 B) — FOUND
commits: 933bdbb, 6f50e4a — FOUND in git log --all
refs: v1.0-pre-retire (tag), retire-restore-point (branch) — FOUND in parent repo, target commit EQUAL on 53bdb5b
testa: 4/5 PNGs FOUND with non-zero size; nlsecu MISSING per documented override #5
```


**1. Created file exists:**
- `.planning/phases/09-retire-x11-backend/09-01-AUDIT-BASELINE.md`: FOUND (15,888 bytes, 6 sections + verification log, 17 references to v1.0-pre-retire)

**2. Commits exist:**
- `933bdbb docs(09-01): commit pre-retirement audit baseline (v1.0-pre-retire gate)`: FOUND on per-agent branch `worktree-agent-adb43afc8c0cdfdc1`

**3. Refs exist (in parent repo):**
- `v1.0-pre-retire` (annotated tag): FOUND at `73d32e3` → target `53bdb5b`
- `retire-restore-point` (branch): FOUND at `53bdb5b`
- Equality (commit-deref): SHA match ✓

**4. Build + sweep evidence:**
- `/tmp/09-01-cbl_tout-pre.log`: FOUND (last line "TOUS les MODULES EXECUTABLES ... sont crees", exit 0)
- `/tmp/09-01-sweep-pre.log`: FOUND (5 verdict=SMOKE lines, exit 0)
- `/tmp/09-01-pre-retire-baseline/{pan2d,nafems_le1,cavity2d,heat1d}-qt-1x.png`: FOUND (sizes 71/320/44/21 KB)
- `/tmp/09-01-pre-retire-baseline/nlsecu-qt-1x.png`: MISSING (per Phase 8 override #5 — documented disposition, NOT a failure)

---

*Phase: 09-retire-x11-backend*
*Plan: 01 (Pre-Phase-9 rollback safety + audit baseline)*
*Completed: 2026-05-05*
