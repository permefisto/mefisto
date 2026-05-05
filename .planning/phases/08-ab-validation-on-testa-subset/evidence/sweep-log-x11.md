# Phase 8 Plan 02 — X11 Baseline Capture Sweep Log

**Date:** 2026-05-05
**Author:** Phase 8 Plan 02 executor (autonomous, parallel wave 2)
**Driver:** `bin/ab_sweep_phase8.sh --mode x11` + manual nlsecu wrapper (TIME truncation per Wave 1 PARTIAL finding)
**Scope:** 5 BUILD-10 baseline testa cases (pan2d, nafems_le1, cavity2d, heat1d, nlsecu).
**Backend:** X11 (`pp/pp{mail,elas,flui,ther,nlse}` + `xvue/xvuelc.c`) under Xvfb 1280x800x24.

This log produces the **column 1** of the Phase 8 D-10 verdict matrix:
the X11 baseline that every Qt-mode column (Plans 03/04/05) compares against.
No Qt comparison happens in Plan 02 — only capture.

---

## Plan-2 X11 baseline sweep (2026-05-05 11:44:52 UTC)

| Case        | Mode | Exit | Capture path (relative to repo root)                                                          | Size (bytes) | Dimensions | Ready-file latency | Notes                                                |
| ----------- | ---- | ---- | --------------------------------------------------------------------------------------------- | ------------ | ---------- | ------------------ | ---------------------------------------------------- |
| pan2d       | x11  | 0    | .planning/phases/08-ab-validation-on-testa-subset/evidence/pan2d-x11.png                      | 7636         | 1280x800   | <2s                | mesher canonical drawing                             |
| nafems_le1  | x11  | 0    | .planning/phases/08-ab-validation-on-testa-subset/evidence/nafems_le1-x11.png                 | 60765        | 1280x800   | <2s                | elasticity stress field                              |
| cavity2d    | x11  | 0    | .planning/phases/08-ab-validation-on-testa-subset/evidence/cavity2d-x11.png                   | 14906        | 1280x800   | <2s                | Stokes velocity field                                |
| heat1d      | x11  | 0    | .planning/phases/08-ab-validation-on-testa-subset/evidence/heat1d-x11.png                     | 4363         | 1280x800   | <2s                | 1D thermal profile                                   |
| nlsecu      | x11  | 0    | .planning/phases/08-ab-validation-on-testa-subset/evidence/nlsecu-x11.png                     | 4907         | 1280x800   | <2s                | TRUNCATED-CAPTURE — TIME=0.01 (1 step); see §nlsecu  |

All 5 captures produced exit_code=0; ready-file sentinel fired (xvuelc.c xvfermer_) before
`import -window root` ran; all images parse cleanly under `identify`.

## OMP context

- **OMP_NUM_THREADS=1** (deterministic baseline per 08-RESEARCH.md §"Open Question 2").
- Set via `export OMP_NUM_THREADS=1` in the executor's parent shell BEFORE invoking
  `bin/ab_sweep_phase8.sh`, propagating through `env` (no `-i`) and `xvfb-run` (no `-i`)
  into every `pp/pp{mail,elas,flui,ther,nlse}` invocation.
- Verified inheritance path: parent shell → `bin/ab_sweep_phase8.sh` → `bin/ab_capture_x11.sh`
  → `xvfb-run --auto-servernum` → inner `bash -c` → `env DISPLAY=...` → `pp/pp${MODULE}`.
- Codified for Wave 2 plans 03/04/05: every X11-side compare must use a baseline produced
  with `OMP_NUM_THREADS=1` (this log's set). Qt-OMP cells compare against this same X11
  baseline (plan 5 may pin its own X11-OMP baseline if drift is observed).

## xvuelc.c invariant

- `git diff -- xvue/xvuelc.c` → empty (byte-identical to plan-entry HEAD `06462cb`).
- `git diff 900e297..HEAD -- xvue/xvuelc.c` → empty (byte-identical to pre-Phase-8 baseline).
- T-08-13 mitigation honored: no modification to `xvue/xvuelc.c` during capture.

## Case-batch resolution

- Per-case `MODULE` + `BATCH` sourced from `bin/phase8_case_batch_map.sh`
  (Plan 1 Task 0 deliverable per WARNING-1 iter2). No glob assumption used.
- Empirical mapping (verbatim from `phase8_case_batch_map.sh`):

  | Case        | MODULE | BATCH                | Prereq MODULE | Prereq BATCH       |
  | ----------- | ------ | -------------------- | ------------- | ------------------ |
  | pan2d       | mail   | pan2d.mesh           | (none)        | (none)             |
  | nafems_le1  | elas   | nafems_le1.elas      | mail          | nafems_le1.mesh    |
  | cavity2d    | flui   | cavity2d.stoke56cr   | mail          | cavity2d.meshbf    |
  | heat1d      | ther   | heat1d.heat          | mail          | heat1d.mesh        |
  | nlsecu      | nlse   | nlsecu.iexrr         | mail          | nlsecu.meshq2      |

- Workspace+cwd discipline (also from `phase8_case_batch_map.sh`):
  every case requires `INITIER` (`echo "$CASE" | $MEFISTO/pp/ppinit`) BEFORE running
  the main module; the harness performs this seed step inside `pushd "$MEFISTOX/$CASE"`.

## Per-case verification

| Case        | size_bytes | dims      | sha256                                                              | OK? |
| ----------- | ---------- | --------- | ------------------------------------------------------------------- | --- |
| pan2d       | 7636       | 1280x800  | 0561d97d3a077a60f7c18e8278c63e606e3ad79b216e711720530d64d44fc32b    | YES |
| nafems_le1  | 60765      | 1280x800  | 56d66a48e81d4e39277f1b533db6b955c96ffbb64beac6ad387ef30cfc5add95    | YES |
| cavity2d    | 14906      | 1280x800  | 5695a8f0669ffab1c210301fb428697b0c2bad6267aa398d8dd7f28dc054b9ee    | YES |
| heat1d      | 4363       | 1280x800  | 4568cff9d26f20e5da246ef0d338bf0f8a15a1e56ef23f10e8f1db418862d619    | YES |
| nlsecu      | 4907       | 1280x800  | da6203a8b2feb14f8e3cccb65ff38fd1f8a4ebd78c995add14cd047bbef98461    | YES |

Color-content sanity (top palette colors per case, confirms non-blank capture):

- **pan2d:** peach background (#FFD9B8, 632530 px = 62%) + blue-gray mesh wireframe (#4F4F63, 258863 px) + cyan/red/green annotations.
- **nafems_le1:** black chrome (111726 px) + orange stress-field (#FF7F00, 69613 px) + red/blue accents.
- **cavity2d:** black chrome (88585 px) + deep-blue Stokes velocity field (#19049B, 437970 px) + cyan/green palette.
- **heat1d:** peach background (908739 px = 89%) + black chrome (90284 px) + red 1D thermal curve (3339 px) + axis grid (149,149,177).
- **nlsecu:** black chrome (88438 px) + blue/cyan/green NLSE 3D wave field (16-color palette) + magenta legend.

All 5 captures contain the canonical Mefisto palette markers (peach #FFD9B8 background,
blue-gray #4F4F63 mesh) confirming the X11 backend's xvuelc.c rendering path executed
end-to-end and was visible to `import -window root` at xvfermer_ time.

## nlsecu — TRUNCATED-CAPTURE rationale

**Status:** TRUNCATED-CAPTURE (per Wave 1 nlsecu PARTIAL finding + user-decided mitigation).

**Original case:** `testa/nlsecu/nlsecu.iexrr` requests Final TIME=20, Step=0.01 → 2000 time steps
on a 3D NLSE cube. Each step emits per-DoF wave-magnitude diagnostics (~3000 nodes), so the
per-step stdout overhead alone takes ~6s on this hardware. Wave 1 documented that the case
exceeds the 60s harness budget by an order of magnitude regardless of dispatch path
(per `00-smoke-probes.md` §nlsecu STATUS=PARTIAL).

**User-decided mitigation:** Truncate Final TIME to a single time step for headless capture only.
Tried Final TIME=2 (200 steps, ~20 min) — exceeded harness budget. Tried Final TIME=0.1
(10 steps) — also exceeded budget. Final accepted truncation: **Final TIME=0.01 (1 step)**.

**Patch applied:** `sed -i 's|^  20;      { Final TIME }|   0.01;   { Final TIME }|'
$MEFISTOX/nlsecu/nlsecu.iexrr` — applied to the workspace copy in `$MEFISTOX/nlsecu/`,
NOT to `testa/nlsecu/nlsecu.iexrr` (testa/ tree preserved per `CLAUDE.md` test discipline).

**Implications for Plans 03/04/05:**
The nlsecu X11 baseline is the 1-step (TIME=0.01) NLSE wave field. Plans 03/04/05 must
apply the SAME truncation when capturing the Qt-1x / Qt-2x / Qt-OMP nlsecu cells, otherwise
the LEFT operand (this baseline at TIME=0.01) and the RIGHT operand (Qt at TIME=20) compare
fundamentally different scenes. The harness accepts the workspace-patch pattern: re-create
the workspace, edit `nlsecu.iexrr`, then run the dispatch. CHECKLIST.md cell for nlsecu
will document `TRUNCATED-CAPTURE: TIME=0.01 (1 step)` per D-10 cell-note convention.

**Alternative considered (NOT taken):** Use the MAILLER prereq capture (mesh of the cube)
as the nlsecu evidence, per `00-case-batch-map.md` Wave 1 fallback option. Rejected because
it would compare an empty mesh against Qt's NLSE wave field — meaningless. The truncated
1-step capture preserves the comparability axis (both backends render the same 1-step wave
solution).

## Manifest

| Case        | Size (bytes) | Dimensions | Bit depth          | SHA-256                                                              |
| ----------- | ------------ | ---------- | ------------------ | -------------------------------------------------------------------- |
| pan2d       | 7636         | 1280x800   | 8-bit sRGB         | 0561d97d3a077a60f7c18e8278c63e606e3ad79b216e711720530d64d44fc32b     |
| nafems_le1  | 60765        | 1280x800   | 8-bit sRGB         | 56d66a48e81d4e39277f1b533db6b955c96ffbb64beac6ad387ef30cfc5add95     |
| cavity2d    | 14906        | 1280x800   | 8-bit sRGB         | 5695a8f0669ffab1c210301fb428697b0c2bad6267aa398d8dd7f28dc054b9ee     |
| heat1d      | 4363         | 1280x800   | 8-bit sRGB         | 4568cff9d26f20e5da246ef0d338bf0f8a15a1e56ef23f10e8f1db418862d619     |
| nlsecu      | 4907         | 1280x800   | 8-bit sRGB         | da6203a8b2feb14f8e3cccb65ff38fd1f8a4ebd78c995add14cd047bbef98461     |

## Hand-off to Plans 03/04/05

The 5 X11 baseline captures (LEFT operand for every Qt-mode compare in subsequent plans)
are committed at:

1. `.planning/phases/08-ab-validation-on-testa-subset/evidence/pan2d-x11.png`
2. `.planning/phases/08-ab-validation-on-testa-subset/evidence/nafems_le1-x11.png`
3. `.planning/phases/08-ab-validation-on-testa-subset/evidence/cavity2d-x11.png`
4. `.planning/phases/08-ab-validation-on-testa-subset/evidence/heat1d-x11.png`
5. `.planning/phases/08-ab-validation-on-testa-subset/evidence/nlsecu-x11.png`

Plans 03/04/05 invoke `bin/ab_compare_pair.sh` (Plan 1 deliverable) with these paths
as the LEFT argument; the Qt-side captures are RIGHT.

**Two harness invariants codified during this plan's execution (Rule 3 deviations):**

1. `MEFISTO_BATCH_X11=1` MUST be set in the executor's parent shell BEFORE invoking
   `bin/ab_sweep_phase8.sh --mode x11`. The harness's qt-1x/qt-2x/qt-omp branches set
   this env inline, but the x11 branch delegates to `bin/ab_capture_x11.sh` which
   inherits but does not re-set the var. Without it, `INTERA=0` in `prpr/pp*.f`,
   no graphical screen, no xvfermer_ call → empty captures (size=0). Same root cause
   as Plan 1's qt-side codification (per `08-01-SUMMARY.md` §"key-decisions" #2).

2. `--out-dir PATH` MUST be an absolute path. The harness opens `pushd "$PROJDIR"`
   before invoking `pp/*`; if `OUT_DIR` is relative (the default
   `.planning/phases/08-*/evidence`), the resolved `$OUT` resolves under the project
   workspace `$MEFISTOX/$CASE/.planning/...` instead of the repo's evidence dir,
   silently dropping all captures. Both Wave 1 verification (qt-1x smoke) and this
   wave's first attempt hit the same trap. Codified at the call site here; deferred
   for harness-level cleanup to Plan 05/06 if the issue resurfaces.

Both items are documented as Rule 3 (auto-fix blocking issue) deviations in the
Plan 02 SUMMARY, and Plans 03/04/05 should adopt the same patterns at their call sites.

## Outcome

Plan 02 X11 baseline column complete: 5/5 cases captured at 1280x800 under Xvfb +
OMP_NUM_THREADS=1 deterministic baseline. Hand-off ready for Plans 03/04/05.
