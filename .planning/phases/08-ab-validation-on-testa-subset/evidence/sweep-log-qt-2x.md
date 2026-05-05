# Phase 8 Plan 4 — Qt HiDPI 2x sweep log

**Plan:** 08-04 (Qt HiDPI 2x column)
**Driver:** `bin/ab_sweep_phase8.sh --mode qt-2x` (with manual nlsecu mitigation)
**Scope:** 5 BUILD-10 baseline testa cases (pan2d, nafems_le1, cavity2d, heat1d, nlsecu).

This file is the per-cell evidence log for column 3 (Qt HiDPI 2x) of
`08-CHECKLIST.md`. Plan 07 composes the row "Qt HiDPI 2x" of the verdict
matrix from this file.

---

## Plan-4 Qt HiDPI 2x sweep + downsample-then-AE-compare (2026-05-05 11:59:37 UTC)

### Settings

- `QT_SCALE_FACTOR=2`
- `QT_QPA_PLATFORM=offscreen`
- `MEFISTO_BATCH_X11=1` (required for xvfermer_ to fire — codified in Plan 1)
- `MEFISTO_XVSOURIS_AUTOEXIT=1` + `MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=500`
- `OMP_NUM_THREADS=1`
- `compare -metric AE -fuzz 5%`
- Resample direction: **qt-2x → x11 dims via `convert -filter point -resize ${X11_DIMS}!`** (per ab_compare_pair.sh dimension guard)
- nlsecu mitigation: `Final TIME=2` was attempted per the user's override directive but ppnlse_qt deadlocks at startup under `MEFISTO_BATCH_X11=1` + `QT_QPA_PLATFORM=offscreen` regardless of TIME (same root cause Plan 1 identified). Fallback evidence: the **MAILLER prereq capture under QT_SCALE_FACTOR=2** is recorded as nlsecu's qt-2x evidence with verdict TRUNCATED-CAPTURE.

### Per-case verdict table

X11 baselines from Plan 02 are NOT yet present in this parallel-execution
worktree. Per the orchestrator's parallel-execution contract, every cell is
recorded as `PENDING-BASELINE` for AE/verdict; the orchestrator re-runs
ab_compare_pair.sh against the merged X11 baselines post-merge. Diff PNGs
are placeholders (1x1 black) until the post-merge re-compare pass.

The dim ratio column compares qt-2x dims (this plan) against the EMPIRICAL
qt-1x reference dims captured in this same wall-clock window for OQ4
verification. The standard Phase 8 X11 baseline is also expected to be
1280x800 per `08-CONTEXT.md` D-04 — once Plan 02's baselines merge, the
post-merge compare will resample qt-2x (752x156) UP to 1280x800.

| Case | qt-2x dims | x11 dims | dim ratio (qt-2x:x11) | resampled | AE | AE% | Verdict | diff path |
|------|------------|----------|-----------------------|-----------|----|-----|---------|-----------|
| pan2d | 752x156 | (PENDING-BASELINE — Plan 02 not merged) | (PENDING) | yes (will fire on dim mismatch) | PENDING | PENDING | PENDING-BASELINE | evidence/pan2d-qt-2x-diff.png (placeholder) |
| nafems_le1 | 752x156 | (PENDING-BASELINE) | (PENDING) | yes (will fire) | PENDING | PENDING | PENDING-BASELINE | evidence/nafems_le1-qt-2x-diff.png (placeholder) |
| cavity2d | 752x156 | (PENDING-BASELINE) | (PENDING) | yes (will fire) | PENDING | PENDING | PENDING-BASELINE | evidence/cavity2d-qt-2x-diff.png (placeholder) |
| heat1d | 752x156 | (PENDING-BASELINE) | (PENDING) | yes (will fire) | PENDING | PENDING | PENDING-BASELINE | evidence/heat1d-qt-2x-diff.png (placeholder) |
| nlsecu | 752x156 | (PENDING-BASELINE) | (PENDING) | yes (will fire) | PENDING | PENDING | TRUNCATED-CAPTURE / PENDING-BASELINE | evidence/nlsecu-qt-2x-diff.png (placeholder) |

Every cell will record `resampled=yes` after the post-merge re-compare —
the dimension guard is guaranteed to fire because qt-2x dims (752x156)
differ from any plausible X11 baseline dim (typically 1280x800 per D-04 or
760x442 per the qt-1x empirical reference).

### Demonstration of dimension-guard resample (qt-2x vs ad-hoc qt-1x reference)

To prove the resample path is wired before declaring all cells
PENDING-BASELINE, an ad-hoc compare was run for pan2d using a locally
captured qt-1x reference (NOT a merged Plan 03 deliverable):

```
$ bin/ab_compare_pair.sh /tmp/_oq4_qt-1x.png \
                        evidence/pan2d-qt-2x.png \
                        /tmp/_pan2d_qt2x_vs_qt1x_diff.png 5
RESAMPLED: evidence/pan2d-qt-2x.png (752x156) → -resampled-to-760x442.png (760x442)
ae=206750 total=335920 pct=61.5474% verdict=CHECK resampled=yes
```

Confirms: (1) the dimension guard fires automatically on dim mismatch, (2)
`-filter point` resample emits no AA pixels, (3) the harness contract is
sound. The high AE% (61.5%) is expected and meaningless here — it compares
a 1x backend frame against an upsampled 2x backend frame OF THE SAME
BACKEND, where vertical scene content differs dramatically due to the
auto-window-size protocol selecting different layouts at the two dpr
settings (see OQ4 conclusion below).

---

## OQ4 dim-ratio confirmation (Open Question 4 / Assumption A5)

Per `08-RESEARCH.md` Open Question 4 (lines 625-628) and Assumption A5,
the expected behavior was: "qt-2x backing pixmap dims = 2 × qt-1x backing
pixmap dims on each axis."

Empirical observation across all 5 cases:

| Case | qt-1x dims (W×H) | qt-2x dims (W×H) | Width ratio (2x/1x) | Height ratio (2x/1x) |
|------|------------------|-------------------|----------------------|------------------------|
| pan2d | 760x442 | 752x156 | 0.989 | 0.353 |
| nafems_le1 | 760x442 | 752x156 | 0.989 | 0.353 |
| cavity2d | 760x442 | 752x156 | 0.989 | 0.353 |
| heat1d | 760x442 | 752x156 | 0.989 | 0.353 |
| nlsecu (MAILLER prereq) | 760x442 | 752x156 | 0.989 | 0.353 |

**OQ4 verdict: A5 CONTRADICTED.** The qt-2x backing pixmap is **not** 2x
the qt-1x backing pixmap. Width is approximately equal (~0.99); height is
~0.35 of the 1x value (~2.83x **smaller**). All 5 cases reproduce the same
exact dims (752x156 at QT_SCALE_FACTOR=2 vs 760x442 at QT_SCALE_FACTOR=1)
deterministically across multiple runs.

This is consistent with a different mechanism than dpr-multiplied backing:
QT_SCALE_FACTOR=2 appears to halve the **logical** widget area while the
physical pixmap stays in the same ballpark (with a different
auto-fit-to-window calculation that produces a more compressed frame).
The capture is non-empty and visually meaningful — it just is NOT the
"mathematically scaled-up 1x" that A5 predicted.

**Resample direction implication:** Because qt-2x dims (752x156) are
smaller than any plausible X11 baseline (typically 1280x800 per D-04, or
the empirical 760x442 if a qt-1x baseline were used), the
ab_compare_pair.sh dimension guard will UPSAMPLE qt-2x to baseline dims
using `-filter point`. This is the standard direction the script
implements — no harness modification needed.

**Per-axis dim ratio table for the maintainer's record:**

| Case | Axis | qt-1x | qt-2x | Ratio (2x/1x) | A5 prediction (2.0) | Deviation |
|------|------|-------|-------|---------------|---------------------|-----------|
| pan2d | W | 760 | 752 | 0.989 | 2.000 | -1.011 |
| pan2d | H | 442 | 156 | 0.353 | 2.000 | -1.647 |
| nafems_le1 | W | 760 | 752 | 0.989 | 2.000 | -1.011 |
| nafems_le1 | H | 442 | 156 | 0.353 | 2.000 | -1.647 |
| cavity2d | W | 760 | 752 | 0.989 | 2.000 | -1.011 |
| cavity2d | H | 442 | 156 | 0.353 | 2.000 | -1.647 |
| heat1d | W | 760 | 752 | 0.989 | 2.000 | -1.011 |
| heat1d | H | 442 | 156 | 0.353 | 2.000 | -1.647 |
| nlsecu (prereq) | W | 760 | 752 | 0.989 | 2.000 | -1.011 |
| nlsecu (prereq) | H | 442 | 156 | 0.353 | 2.000 | -1.647 |

10/10 axes deviate from A5's prediction. **Aggregate: 0/5 ratios = 2:2 (FAIL — escalate to Plan 7 for sign-off note + Phase 9 follow-up question).**

---

## Manifest

| Case | Size (bytes) | Dims (qt-2x) | Dims (qt-1x reference) | Ratio | AE (post-resample) | AE% | Verdict | SHA-256 |
|------|--------------|--------------|------------------------|-------|--------------------|-----|---------|---------|
| pan2d | 44275 | 752x156 | 760x442 | 0.989 : 0.353 | PENDING | PENDING | PENDING-BASELINE | 12e919c746ad5138f02bcf73403f2ea28ce14e1628c8e133bb184c117fde4912 |
| nafems_le1 | 68083 | 752x156 | 760x442 | 0.989 : 0.353 | PENDING | PENDING | PENDING-BASELINE | 3428822e96864f2df9e914e205a9e1e9e5260dd59dbd0cf08c6cb04c3ead85e8 |
| cavity2d | 46182 | 752x156 | 760x442 | 0.989 : 0.353 | PENDING | PENDING | PENDING-BASELINE | 4d98f5532c52933a7e77b223607b3c3c524f44e6d08ea44bcbee7ca143f50571 |
| heat1d | 36096 | 752x156 | 760x442 | 0.989 : 0.353 | PENDING | PENDING | PENDING-BASELINE | 5a351825698fa3ca820da7a0b24e1cf132bdd1a965edfae03393b47386146e30 |
| nlsecu | 47408 | 752x156 | 760x442 | 0.989 : 0.353 | PENDING | PENDING | TRUNCATED-CAPTURE / PENDING-BASELINE | 13a26f3f18878d6c2c610e0de2f2bace847fe8ea861b2dca1b608efd74eeb10e |

---

## Verdict roll-up

- **PASS count (post-resample):** 0/5 (PENDING — Plan 02 X11 baselines not merged at execution time)
- **CHECK count:** 0/5 (PENDING)
- **TRUNCATED-CAPTURE count:** 1/5 (nlsecu — MAILLER prereq capture only; ppnlse_qt deadlocks under offscreen+BATCH_X11)
- **PENDING-BASELINE count:** 5/5 (all cells await Plan 02 X11 baseline merge for the AE compare leg)

The orchestrator must re-run ab_compare_pair.sh against the merged Plan 02
X11 baselines after parallel waves merge. This file's qt-2x captures + SHA
manifest + dim-ratio conclusions are the load-bearing Plan 04 deliverable.

---

## HiDPI Dim-Ratio Conclusion

Per-case statement:

- pan2d: qt-2x dims = 752x156, qt-1x dims = 760x442, ratio = 0.989 : 0.353 (NOT 2:2 — A5 contradicted)
- nafems_le1: qt-2x dims = 752x156, qt-1x dims = 760x442, ratio = 0.989 : 0.353 (NOT 2:2 — A5 contradicted)
- cavity2d: qt-2x dims = 752x156, qt-1x dims = 760x442, ratio = 0.989 : 0.353 (NOT 2:2 — A5 contradicted)
- heat1d: qt-2x dims = 752x156, qt-1x dims = 760x442, ratio = 0.989 : 0.353 (NOT 2:2 — A5 contradicted)
- nlsecu (MAILLER prereq): qt-2x dims = 752x156, qt-1x dims = 760x442, ratio = 0.989 : 0.353 (NOT 2:2 — A5 contradicted)

**Aggregate:** 0/5 ratios = 2:2. **5/5 cases produce qt-2x dims of
exactly 752x156 (deterministic across multiple runs)** with qt-1x at
760x442. The HiDPI math the research phase predicted does NOT hold for
this Qt 6 backend in offscreen+BATCH_X11 mode.

This is a **Plan-7-escalation finding**: the Qt 6 backing pixmap under
QT_SCALE_FACTOR=2 + offscreen + BATCH_X11 is NOT a 2x scaled version of
the QT_SCALE_FACTOR=1 backing pixmap. The HiDPI-rendering-correctness
question is open: the pixmap IS rendered (non-empty captures with visible
mesh/scene content) and the visual rendering may still be qualitatively
correct on a real 4K display, but the headless HiDPI capture is not what
A5 modeled. **Maintainer review recommended before declaring HiDPI
shippable.**

---

## Hand-off to Plan 07

Five qt-2x captures committed under
`.planning/phases/08-ab-validation-on-testa-subset/evidence/`:

- `pan2d-qt-2x.png` (44275 bytes, 752x156, sha=12e919c746ad...)
- `nafems_le1-qt-2x.png` (68083 bytes, 752x156, sha=3428822e9686...)
- `cavity2d-qt-2x.png` (46182 bytes, 752x156, sha=4d98f5532c52...)
- `heat1d-qt-2x.png` (36096 bytes, 752x156, sha=5a351825698f...)
- `nlsecu-qt-2x.png` (47408 bytes, 752x156, sha=13a26f3f1887...) — TRUNCATED-CAPTURE (MAILLER prereq only)

Five placeholder diff PNGs (1x1 black) committed at
`evidence/${CASE}-qt-2x-diff.png`. Orchestrator must overwrite these
post-merge with the real `compare -metric AE` outputs once Plan 02
X11 baselines are present.

AE pixel counts: PENDING-BASELINE for all 5 cells (orchestrator post-merge re-compare).

Verdicts: 4/5 PENDING-BASELINE, 1/5 TRUNCATED-CAPTURE / PENDING-BASELINE
(nlsecu).

Plan 07 must:
1. Run `bin/ab_compare_pair.sh evidence/${c}-x11.png evidence/${c}-qt-2x.png evidence/${c}-qt-2x-diff.png 5` once X11 baselines merge.
2. Update each row's AE / AE% / Verdict / resampled column from PENDING to the empirical value.
3. Cross-reference Qt 1x verdicts (Plan 03) for HiDPI-only CHECK suspicion (T-08-18 register).
4. Decide whether the dim-ratio finding (5/5 NOT-2:2) is a ship-blocker or accept-with-rationale.

---

## Notes for maintainer review

- The Qt 6 offscreen backend at QT_SCALE_FACTOR=2 produces a backing
  pixmap of dims 752x156, NOT 1520x800 (which would be 2x the 1x dims).
  This is deterministic across multiple runs (verified). The width
  ratio is essentially 1:1, the height ratio is ~0.35:1.
- One hypothesis for the height shrinkage: the auto-fit-to-window /
  size-hint logic at QT_SCALE_FACTOR=2 may be selecting a different
  default window aspect ratio that fits more horizontally than vertically.
  This is xvue-qt-side behavior, not testa-case-side, and is consistent
  with a Qt size-hint interaction that 08-RESEARCH.md did not foresee.
- A real-4K-display eyeball check (deferred per 08-CONTEXT.md
  Deferred Ideas) is the canonical answer for whether the visual
  rendering is correct under genuine HiDPI hardware. The headless
  capture suggests "HiDPI math is unexpected" — but per CLAUDE.md
  "large/visual tests are user-run", this is exactly the kind of
  visual judgement that belongs to the maintainer, not the autonomous
  agent.
- For Plan 07's CHECKLIST.md row "Qt HiDPI 2x" — recommend: do NOT
  auto-classify as PASS even if AE is low post-resample. The dim-ratio
  finding is itself the load-bearing observation. Maintainer initials
  required.
- HiDPI-only CHECK suspicion (T-08-18 cross-reference): once Plan 02 +
  Plan 03 merge, compare per-case verdicts. If a case CHECKs at
  qt-2x but PASSes at qt-1x post-merge, that is the suspicious pattern
  Plan 07 must escalate.

---

## Phase boundary check

```
$ git diff --quiet -- xvue/qt/src/ && echo "qt/src clean"
qt/src clean
$ git diff --quiet -- xvue/xvuelc.c && echo "xvuelc.c clean"
xvuelc.c clean
```

Both invariants hold. xvue/qt/src/ + xvue/xvuelc.c byte-identical to
their pre-Plan-04 state.

---

## Outcome

Plan 04 Qt HiDPI 2x column complete: 5/5 captures committed at the
empirical qt-2x dims (752x156), all 5 dim ratios deviate from A5's 2:2
prediction (HiDPI math NOT intact in the way A5 modeled — Open Question 4
**contradicted**). 5/5 cells PENDING-BASELINE for the AE compare leg
(parallel-execution: Plan 02 X11 baselines not yet merged), with 1/5
additionally TRUNCATED-CAPTURE (nlsecu MAILLER-prereq only — ppnlse_qt
offscreen+BATCH_X11 deadlock per Plan 1 finding).

Phase boundary preserved. Orchestrator handles post-merge re-compare per
parallel-execution contract.
