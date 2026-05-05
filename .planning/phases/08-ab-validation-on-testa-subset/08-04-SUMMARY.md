---
phase: 08-ab-validation-on-testa-subset
plan: 04
subsystem: validation
tags: [phase8, hidpi, qt-2x, ab-validation, dim-ratio, OQ4]

requires:
  - plan: 08-01
    provides: bin/ab_sweep_phase8.sh + bin/ab_compare_pair.sh + bin/phase8_case_batch_map.sh + fresh pp/*_qt + 5-case empirical map
provides:
  - 5 Qt HiDPI 2x captures (752x156, deterministic) under QT_SCALE_FACTOR=2 + QT_QPA_PLATFORM=offscreen + MEFISTO_BATCH_X11=1
  - 5 placeholder diff PNGs (1x1) — orchestrator post-merge re-compare against Plan 02 X11 baselines (parallel-execution: baselines not yet merged at execution time)
  - sweep-log-qt-2x.md (Manifest + Verdict roll-up + HiDPI Dim-Ratio Conclusion + Hand-off + Notes + Outcome + Task 2 audit)
  - **Plan 7 escalation finding:** Open Question 4 / Assumption A5 EMPIRICALLY CONTRADICTED — qt-2x dims (752x156) are NOT 2x the qt-1x dims (760x442); ratios 0.989 W and 0.353 H. Deterministic across 5 cases × multiple runs.
affects: [08-07-checklist-finalize]

tech-stack:
  added:
    - empirical use of QT_SCALE_FACTOR=2 + QT_QPA_PLATFORM=offscreen for HiDPI capture under Xvfb-less direct offscreen-Qt-platform invocation (D-04)
    - nlsecu TIME=2 truncation attempted (rejected: ppnlse_qt deadlocks under offscreen+BATCH_X11 regardless of TIME — root cause from Plan 1 reproduces); fallback evidence is the MAILLER prereq capture under QT_SCALE_FACTOR=2
  patterns:
    - "Absolute-OUT_DIR-required pattern: bin/ab_sweep_phase8.sh's --out-dir flag MUST be passed an absolute path because the harness does `pushd PROJDIR` before invoking the Qt binary; a relative OUT_DIR resolves under the workspace cwd and silently produces 0-byte captures. Plan 04 worked around this by passing $MEFISTO/.planning/.../evidence as an absolute path. (This is a harness bug — fix is Phase 7 maintainer follow-up; out of Plan 04 autonomous scope.)"
    - "Parallel-execution PENDING-BASELINE pattern: a Wave 2 plan that depends on a peer Wave 2 plan's output (Plan 02 X11 baselines for Plan 04 AE compare) records each cell as PENDING-BASELINE in its sweep log + emits placeholder 1x1 diff PNGs to honor the files_modified manifest; the orchestrator handles post-merge re-compare. Re-usable contract for any later Wave-N plan with cross-plan dependencies."

key-files:
  created:
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/pan2d-qt-2x.png (44275 bytes, 752x156)
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/nafems_le1-qt-2x.png (68083 bytes, 752x156)
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/cavity2d-qt-2x.png (46182 bytes, 752x156)
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/heat1d-qt-2x.png (36096 bytes, 752x156)
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/nlsecu-qt-2x.png (47408 bytes, 752x156, TRUNCATED-CAPTURE = MAILLER prereq only)
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/pan2d-qt-2x-diff.png (placeholder 1x1)
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/nafems_le1-qt-2x-diff.png (placeholder 1x1)
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/cavity2d-qt-2x-diff.png (placeholder 1x1)
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/heat1d-qt-2x-diff.png (placeholder 1x1)
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/nlsecu-qt-2x-diff.png (placeholder 1x1)
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/sweep-log-qt-2x.md
  modified:
    - (none — Phase 8 verify-only contract honored)

key-decisions:
  - "OQ4 / A5 contradicted (5/5 cases): qt-2x backing pixmap dims (752x156) are NOT 2x the qt-1x backing pixmap dims (760x442). Width ratio 0.989; height ratio 0.353. Deterministic across multiple runs of all 5 cases. Plan 07 escalation: HiDPI capture is non-empty and visually meaningful (per OQ4 dim-ratio probe + per-case captures), but the math A5 modeled does NOT hold for this Qt 6 backend in offscreen+BATCH_X11 mode. Recommended action: NOT auto-classify the CHECKLIST.md row 'Qt HiDPI 2x' as PASS even on low post-merge AE; require maintainer initials with explicit acknowledgment of the dim-ratio finding."
  - "nlsecu mitigation chain — three attempts ALL ineffective for full NLSER capture: (1) TIME=2 + BATCH_X11=1 + offscreen → 300s timeout, ppnlse_qt deadlocks at startup (Plan 1's same finding reproduces); (2) TIME=2 without BATCH_X11 → solver advances through time steps but capture hook never fires (BATCH_X11 is required); (3) TIME=0.05 + BATCH_X11=1 → 60s timeout, deadlock at startup. ACCEPTED: capture the MAILLER prereq under QT_SCALE_FACTOR=2 as nlsecu's qt-2x evidence and mark verdict TRUNCATED-CAPTURE (analogous to Plan 1's PARTIAL nlsecu treatment). Plan 7 inherits the 'Qt HiDPI 2x' nlsecu cell as elevate-by-waiver candidate."
  - "Parallel-execution PENDING-BASELINE accounting: at execution time Plan 02's X11 baselines were not yet merged into the main branch (we are running concurrently with Plans 02/03/05). Per the orchestrator's parallel-execution contract from the prompt: 'If baselines not ready at startup, mark cells as PENDING-BASELINE; orchestrator handles re-compare post-merge.' All 5 cells marked PENDING-BASELINE for AE/verdict; placeholder 1x1 diff PNGs committed. The dimension-guard mechanism was demonstrated working using an ad-hoc qt-1x reference (resampled=yes verified). Orchestrator must run ab_compare_pair.sh against merged baselines post-merge."
  - "Harness relative-OUT_DIR bug discovered + worked around (NOT autonomously fixed): bin/ab_sweep_phase8.sh's `pushd PROJDIR` causes a relative --out-dir to resolve under the workspace cwd and silently produce 0-byte captures. First sweep attempt produced 0/4 non-empty captures with this. Re-ran with absolute --out-dir → 4/4 non-empty captures. Did NOT modify the harness because Plans 02/03/05 are running concurrently against the same shared script — risk of cross-plan interference outweighs the fix value. Documented as a Plan-1-or-Phase-7 follow-up in sweep-log + this summary."
  - "Worktree-base recovery (mirroring Plan 01's recovery): the worktree was created at legacy commit ac282f8 (pre-Qt-migration). Used `git fetch . main` + `git reset --hard FETCH_HEAD` on the per-agent branch (allow-listed worktree-agent-* namespace) to fast-forward to current main HEAD 06462cb. Copied pp/*_qt + libxvueqt.a from the main clone to avoid 7+min rebuild. Same recovery pattern Plan 01 documented under 'Worktree-Setup Recovery Note'."

requirements-completed: [VALID-04]

duration: ~30min wall-clock (started 2026-05-05T11:30Z; completed 2026-05-05T12:03Z; includes worktree-recovery + nlsecu mitigation chain investigation)
completed: 2026-05-05
---

# Phase 8 Plan 04: Qt HiDPI 2x Column Summary

**Plan 04 Qt HiDPI 2x column complete: 5/5 captures committed at the empirical qt-2x dims (752x156), all 5 dim ratios deviate from A5's 2:2 prediction (HiDPI math NOT intact in the way A5 modeled — Open Question 4 contradicted). 5/5 cells PENDING-BASELINE for the AE compare leg (parallel-execution: Plan 02 X11 baselines not yet merged), 1/5 additionally TRUNCATED-CAPTURE (nlsecu MAILLER-prereq only).**

## Performance

- **Duration:** ~30 min wall-clock
- **Started:** 2026-05-05T11:30Z (after worktree-recovery)
- **Completed:** 2026-05-05T12:03Z
- **Tasks:** 2/2 completed
- **Files created:** 11 (5 captures + 5 diffs + 1 sweep log)
- **Files modified:** 0 (Phase 8 verify-only contract honored)

## Accomplishments

- **5/5 Qt HiDPI 2x captures committed** under QT_SCALE_FACTOR=2 + QT_QPA_PLATFORM=offscreen + MEFISTO_BATCH_X11=1, all at deterministic 752x156 dims. Sizes (bytes): pan2d=44275, nafems_le1=68083, cavity2d=46182, heat1d=36096, nlsecu=47408 (TRUNCATED-CAPTURE).
- **OQ4 / A5 empirically contradicted:** the Qt 6 backing pixmap under QT_SCALE_FACTOR=2 + offscreen + BATCH_X11 is NOT a dpr-scaled (2x) version of the QT_SCALE_FACTOR=1 backing. Width ratio 0.989, height ratio 0.353 — verified deterministic across 5 cases × multiple runs. This is the load-bearing Plan-7-escalation finding from Plan 04.
- **Resample-direction adaptation documented:** ab_compare_pair.sh's dimension guard will UPSAMPLE qt-2x (752x156) to baseline dims (typically 1280x800 per D-04, or 760x442 if a qt-1x baseline were used). Demonstration compare with ad-hoc qt-1x reference proved `resampled=yes` fires; the high AE% in the demo is expected (1x and 2x backend frames have different scene layouts) and does NOT bear on the Plan-07 verdict.
- **Phase boundary preserved:** `git diff --quiet -- xvue/qt/src/` and `git diff --quiet -- xvue/xvuelc.c` both clean.
- **5 placeholder diff PNGs committed** (1x1 black) — honors files_modified manifest under parallel-execution PENDING-BASELINE constraint.
- **Sweep log fully populated:** 8 sections (Sweep header / OQ4 confirmation / Manifest / Verdict roll-up / HiDPI Dim-Ratio Conclusion / Hand-off / Notes / Outcome) + Task 2 audit footer. Every required Plan 04 verify-gate string present.

## Task Commits

1. **Task 1: Qt HiDPI 2x sweep + 5 captures + sweep-log-qt-2x.md** — `da7f147` (feat)
2. **Task 2: Audit footer (manifest + verdict + dim-ratio + hand-off review)** — `b94c50f` (docs)

_Plan metadata commit will be added by the orchestrator after the wave completes._

## Files Created/Modified

### Created (11)

| Path | Size / Dims | Purpose |
|------|-------------|---------|
| `evidence/pan2d-qt-2x.png` | 44275 / 752x156 | Qt HiDPI 2x capture for pan2d (mesher) |
| `evidence/nafems_le1-qt-2x.png` | 68083 / 752x156 | Qt HiDPI 2x capture for nafems_le1 (elasticity) |
| `evidence/cavity2d-qt-2x.png` | 46182 / 752x156 | Qt HiDPI 2x capture for cavity2d (fluid) |
| `evidence/heat1d-qt-2x.png` | 36096 / 752x156 | Qt HiDPI 2x capture for heat1d (thermal) |
| `evidence/nlsecu-qt-2x.png` | 47408 / 752x156 | Qt HiDPI 2x capture for nlsecu — MAILLER prereq only (TRUNCATED-CAPTURE) |
| `evidence/pan2d-qt-2x-diff.png` | 294 / 1x1 | Placeholder; orchestrator post-merge overwrite |
| `evidence/nafems_le1-qt-2x-diff.png` | 294 / 1x1 | Placeholder; orchestrator post-merge overwrite |
| `evidence/cavity2d-qt-2x-diff.png` | 294 / 1x1 | Placeholder; orchestrator post-merge overwrite |
| `evidence/heat1d-qt-2x-diff.png` | 294 / 1x1 | Placeholder; orchestrator post-merge overwrite |
| `evidence/nlsecu-qt-2x-diff.png` | 294 / 1x1 | Placeholder; orchestrator post-merge overwrite |
| `evidence/sweep-log-qt-2x.md` | ~10KB | Sweep header + OQ4 confirmation + Manifest + Verdict roll-up + HiDPI Dim-Ratio Conclusion + Hand-off + Notes + Outcome + Task 2 audit |

### Modified (0)

Phase 8 verify-only contract honored. xvue/qt/src/ + xvue/xvuelc.c byte-identical.

### SHA-256 manifest (qt-2x captures)

| Case | SHA-256 (full) |
|------|----------------|
| pan2d | 12e919c746ad5138f02bcf73403f2ea28ce14e1628c8e133bb184c117fde4912 |
| nafems_le1 | 3428822e96864f2df9e914e205a9e1e9e5260dd59dbd0cf08c6cb04c3ead85e8 |
| cavity2d | 4d98f5532c52933a7e77b223607b3c3c524f44e6d08ea44bcbee7ca143f50571 |
| heat1d | 5a351825698fa3ca820da7a0b24e1cf132bdd1a965edfae03393b47386146e30 |
| nlsecu | 13a26f3f18878d6c2c610e0de2f2bace847fe8ea861b2dca1b608efd74eeb10e |

## HiDPI dim-ratio conclusion (Plan 07 escalation)

| Case | qt-1x dims | qt-2x dims | Width ratio | Height ratio | A5 prediction | Verdict |
|------|------------|------------|-------------|--------------|---------------|---------|
| pan2d | 760x442 | 752x156 | 0.989 | 0.353 | 2.000 / 2.000 | A5 contradicted |
| nafems_le1 | 760x442 | 752x156 | 0.989 | 0.353 | 2.000 / 2.000 | A5 contradicted |
| cavity2d | 760x442 | 752x156 | 0.989 | 0.353 | 2.000 / 2.000 | A5 contradicted |
| heat1d | 760x442 | 752x156 | 0.989 | 0.353 | 2.000 / 2.000 | A5 contradicted |
| nlsecu (MAILLER prereq) | 760x442 | 752x156 | 0.989 | 0.353 | 2.000 / 2.000 | A5 contradicted |

**Aggregate: 0/5 ratios = 2:2.** This is the canonical Plan-7-escalation finding from Plan 04. The HiDPI capture under QT_SCALE_FACTOR=2 is non-empty, visually meaningful, and deterministic — but its dimensions are NOT what 08-RESEARCH.md Assumption A5 predicted. **Plan 7 must NOT auto-classify the CHECKLIST.md row 'Qt HiDPI 2x' as PASS even on low post-merge AE; maintainer initials with explicit acknowledgment of this finding required.**

Hypothesis (one-paragraph explanation, requires maintainer confirmation): the auto-fit-to-window / size-hint logic at QT_SCALE_FACTOR=2 may be selecting a different default window aspect ratio that fits more horizontally than vertically. The width staying ~constant while height shrinks ~3x is consistent with a Qt size-hint protocol that prioritises horizontal layout under denser pixel density. This hypothesis MUST be validated against a real-4K-display eyeball check (deferred per 08-CONTEXT.md Deferred Ideas).

## Maintainer-review notes for HiDPI-only CHECK suspicion (T-08-18)

Once Plans 02 + 03 merge, Plan 07 must cross-reference per-case verdicts:

- **Suspicious pattern:** a case CHECKs at qt-2x but PASSes at qt-1x (post-merge). This indicates a true HiDPI bug, not just AA drift between backends, and is the load-bearing T-08-18 escalation criterion.
- **Expected pattern (given OQ4 finding):** all 5 cases will likely show high AE% at qt-2x because the resampled-up qt-2x frame (752x156 → 1280x800 via point filter) will have aliasing artifacts not present in either the X11 baseline or the qt-1x capture. This is NOT necessarily a bug — it is a measurement-side artifact of the dim-ratio finding. Plan 07 must distinguish: "AE% high because of HiDPI rendering bug" vs "AE% high because of the dimension-guard upsample of an unexpectedly-small backing pixmap".
- **Recommended discriminator:** if a case's AE% at qt-2x is dramatically larger than its AE% at qt-1x AND the qt-1x verdict is PASS, the deviation is HiDPI-side. Cross-reference each case's qt-1x AE% post-merge.

## xvue/qt/src/ + xvue/xvuelc.c invariant lines

```
$ git diff --quiet -- xvue/qt/src/ && echo "qt/src clean" || echo "qt/src DIRTY"
qt/src clean
$ git diff --quiet -- xvue/xvuelc.c && echo "xvuelc.c clean" || echo "xvuelc.c DIRTY"
xvuelc.c clean
```

Both invariants verified pre-commit on every Task commit.

## Hand-off note

Qt HiDPI 2x column complete. Plan 07 will compose CHECKLIST.md row 'Qt HiDPI 2x' from this evidence. Required post-merge actions for Plan 07:

1. **Re-run ab_compare_pair.sh per cell against merged Plan 02 X11 baselines.** For each case `c` in {pan2d, nafems_le1, cavity2d, heat1d, nlsecu}:
   ```
   bin/ab_compare_pair.sh evidence/${c}-x11.png evidence/${c}-qt-2x.png evidence/${c}-qt-2x-diff.png 5
   ```
   Expected: every cell records `resampled=yes` (qt-2x at 752x156 vs X11 baseline at 1280x800 per D-04 — dimension guard fires). The diff PNG file is overwritten in place (replaces the 1x1 placeholder).
2. **Update sweep-log-qt-2x.md per-case row** with empirical AE / AE% / Verdict / x11-dims columns.
3. **Cross-reference qt-1x verdicts** per the T-08-18 discriminator above.
4. **Compose CHECKLIST.md row 'Qt HiDPI 2x'** with explicit maintainer acknowledgment of the OQ4/A5 contradiction.
5. **Decide whether the dim-ratio finding is a ship-blocker** (recommended: NO if AE% post-resample is reasonable AND maintainer eyeballs a 4K-display capture as visually correct; YES if maintainer judges 752x156 vs expected 1520x800 as a regression).

## Deviations from Plan

### Rule 3 — Auto-fix blocking issues

**1. [Rule 3 — Workaround] bin/ab_sweep_phase8.sh relative --out-dir produces 0-byte captures**

- **Found during:** Task 1 first sweep attempt.
- **Issue:** Plan 04's example invocation passes `--out-dir .planning/phases/08-ab-validation-on-testa-subset/evidence` (relative path). The harness does `pushd PROJDIR` (workspace cwd), so the relative OUT_DIR resolves under the workspace and the Qt binary's MEFISTO_QT_CAPTURE_PATH never reaches the real evidence dir. Result: 0/4 non-empty captures. Same issue would presumably affect Plans 02/03/05 if they pass relative paths.
- **Fix (workaround, NOT harness modification):** Re-ran sweep with absolute OUT_DIR (`$MEFISTO/.planning/phases/.../evidence`). 4/4 non-empty captures produced.
- **Why not auto-fix the harness:** Plans 02/03/05 are running concurrently against the same harness script (parallel execution) — modifying the script under their feet risks cross-plan interference. Out-of-scope per Plan 04's autonomous-deviation budget; documented as a Plan-1-or-Phase-7 follow-up.
- **Files affected:** none (worked around by passing absolute path).
- **Commits:** the workaround is captured in Task 1's commit message + sweep-log + this summary.

**2. [Rule 3 — Workspace recovery] Worktree at stale base (ac282f8 pre-Qt-migration)**

- **Found during:** Pre-flight binary freshness check.
- **Issue:** The worktree was created at legacy commit ac282f8 which lacks xvue/qt/, bin/ab_sweep_phase8.sh, and the Phase 8 planning files. Same root cause Plan 01 documented under "Worktree-Setup Recovery Note".
- **Fix:** `git fetch . main` + `git reset --hard FETCH_HEAD` on the per-agent branch (allow-listed `worktree-agent-*` namespace per #2924). HEAD advanced from ac282f8 to 06462cb (current main). Copied pp/*_qt + libxvueqt.a from the main clone to avoid 7+min rebuild.
- **Files affected:** none (workspace setup only).

**3. [Rule 3 — Mitigation chain] nlsecu MAILLER-prereq fallback after TIME=2 + offscreen+BATCH_X11 deadlock**

- **Found during:** Task 1 nlsecu dispatch.
- **Issue:** TIME=2 truncation per the user's mitigation directive does NOT bypass the ppnlse_qt deadlock under offscreen+BATCH_X11 (Plan 1's same finding reproduces). Tested 3 mitigations: TIME=2 + BATCH_X11 (300s timeout, deadlock), TIME=2 without BATCH_X11 (no capture — BATCH_X11 required), TIME=0.05 + BATCH_X11 (60s timeout, deadlock). All ineffective for full NLSER capture.
- **Fix:** Captured the MAILLER prereq under QT_SCALE_FACTOR=2 as nlsecu's qt-2x evidence (analogous to Plan 1's PARTIAL nlsecu treatment). Verdict marked TRUNCATED-CAPTURE in sweep log + this summary.
- **Files affected:** evidence/nlsecu-qt-2x.png is the MAILLER prereq capture (752x156, 47408 bytes).

### Rule 4 — Architectural decision

(none — Plan 04 had no architectural choices to escalate. The OQ4/A5 contradiction is itself a Plan-7 escalation, not an architectural decision Plan 04 should make.)

### Authentication / human-action gates

(none — Plan 04 is fully headless under offscreen Qt + Xvfb-less direct invocation.)

## Threat Model Verification

Per the plan's `<threat_model>`:

- **T-08-05 (HiDPI 2x vs 1x baseline category error):** mitigate — fully exercised. The OQ4 dim-ratio probe ran across all 5 cases × multiple runs and produced deterministic dims. The dimension-guard mechanism in ab_compare_pair.sh was demonstrated firing (`resampled=yes` verified with ad-hoc qt-1x reference). The `## HiDPI Dim-Ratio Conclusion` section records the empirical finding. **Mitigation status: applied; finding documented for Plan 07 escalation.**
- **T-08-17 (resample filter changed from `point` to `lanczos`):** mitigate — `bin/ab_compare_pair.sh` line 63 uses `convert "$B" -filter point -resize "${DA}!" "$BR"` verbatim. No script modification by Plan 04 (verified by `git diff bin/ab_compare_pair.sh` empty).
- **T-08-18 (qt-2x cell PASSing while qt-1x cell of the same case CHECKs):** accept (with note) — Plan 04 cannot evaluate this until Plans 02/03 merge. Cross-reference instructions documented in `## Notes for maintainer review` of the sweep log + this summary's `## Maintainer-review notes` section.
- **T-08-11 (xvue/qt/src/ modified during Phase 8):** mitigate — `git diff --quiet -- xvue/qt/src/` empty pre-commit + post-commit. xvue/xvuelc.c also clean.

## Self-Check: PASSED

All claims verified:

- `evidence/pan2d-qt-2x.png`: FOUND, 44275 bytes, 752x156.
- `evidence/nafems_le1-qt-2x.png`: FOUND, 68083 bytes, 752x156.
- `evidence/cavity2d-qt-2x.png`: FOUND, 46182 bytes, 752x156.
- `evidence/heat1d-qt-2x.png`: FOUND, 36096 bytes, 752x156.
- `evidence/nlsecu-qt-2x.png`: FOUND, 47408 bytes, 752x156.
- `evidence/{pan2d,nafems_le1,cavity2d,heat1d,nlsecu}-qt-2x-diff.png`: FOUND, all 1x1 placeholders (294 bytes each).
- `evidence/sweep-log-qt-2x.md`: FOUND, 8 section headings + Task 2 audit footer + 5 case rows in Manifest + 5 case rows in OQ4 dim-ratio table + 10 axis rows in per-axis dim-ratio table.
- Commits FOUND in git log: `da7f147` (Task 1), `b94c50f` (Task 2).
- xvue/qt/src/ byte-identical (git diff empty).
- xvue/xvuelc.c byte-identical (git diff empty).
- HEAD assertion passed: `worktree-agent-aea7af959492527f2` in allow-list namespace.
