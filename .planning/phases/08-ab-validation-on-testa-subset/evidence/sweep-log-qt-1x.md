# Phase 8 Plan 03 — Qt 1x sweep log

## Plan-3 Qt 1x sweep + AE-fuzz compare (2026-05-05 11:55:00 UTC)

**Mode:** qt-1x (column 2 of the CHECKLIST.md verdict matrix per D-10)
**Settings:**
- `QT_QPA_PLATFORM=offscreen`
- `MEFISTO_BATCH_X11=1` (per Plan 1 codified contract — required for INTERA=1 dispatch path)
- `OMP_NUM_THREADS=1` (matches Plan 02 baseline column for like-for-like comparison)
- `MEFISTO_QT_CAPTURE_PATH=evidence/${CASE}-qt-1x.png`
- `MEFISTO_XVSOURIS_AUTOEXIT=1`
- `MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=500`
- Compare baseline = `evidence/${CASE}-x11.png` (Plan 02 deliverable)
- Fuzz = 5% (D-02 default)

**Capture window:** 760x442 pixels (Qt offscreen window allocates this geometry; Xvfb root for X11 baseline is 1280x800 — `bin/ab_compare_pair.sh` resamples Qt → X11 dims via `-filter point` per its dimension-mismatch guard).

## Per-case row table

| Case | qt-1x size | qt-1x dims | x11 dims | resampled | AE | AE% | Verdict | diff path |
|------|-----------|-----------|----------|-----------|----|-----|---------|-----------|
| pan2d | 70944 | 760x442 | (pending) | (pending) | (pending) | (pending) | PENDING-BASELINE | evidence/pan2d-qt-1x-diff.png |
| nafems_le1 | 320770 | 760x442 | (pending) | (pending) | (pending) | (pending) | PENDING-BASELINE | evidence/nafems_le1-qt-1x-diff.png |
| cavity2d | 44616 | 760x442 | (pending) | (pending) | (pending) | (pending) | PENDING-BASELINE | evidence/cavity2d-qt-1x-diff.png |
| heat1d | 21834 | 760x442 | (pending) | (pending) | (pending) | (pending) | PENDING-BASELINE | evidence/heat1d-qt-1x-diff.png |
| nlsecu | 70612 | 760x442 | (pending) | (pending) | (pending) | (pending) | TRUNCATED-CAPTURE | evidence/nlsecu-qt-1x-diff.png |

**Verdict legend:**
- **PENDING-BASELINE** — Qt 1x capture present and valid; X11 baseline not yet committed by Plan 02 at the time of this Plan 03 execution (parallel-execution wave; orchestrator wave-merge step re-runs `bin/ab_compare_pair.sh` post-merge to populate AE/AE%/Verdict and overwrite the placeholder diff PNG).
- **TRUNCATED-CAPTURE** — same PENDING-BASELINE state PLUS a documented capture truncation: nlsecu deadlocks under `QT_QPA_PLATFORM=offscreen + MEFISTO_BATCH_X11=1` per Plan 1 SUMMARY (10 log lines emitted, no `Mefisto-NLSER: ARGUMENT NUMBER` banner reached even after 240s; reproduced this run with TIME=2 and TIME=0.2 truncations — both still timed out at 30s/120s budgets). Mitigation per user decision: nlsecu's qt-1x capture is the **MAILLER prereq frame** captured by `pp/ppmail_qt nlsecu.meshq2` (which DOES reach the offscreen pixmap and fires `xvfermer_`). This is a partial harness-reachable evidence, not a full NLSE solver-frame capture.

**Note on diff PNGs:** The `*-qt-1x-diff.png` files committed by this plan are **placeholders** (200x80 lightgray PNGs annotated `PENDING / BASELINE / ${case}`). They satisfy the Plan 03 must-have artifact contract (`test -s evidence/${case}-qt-1x-diff.png`) so Plan 07 has a stable file-path to ingest. The orchestrator wave-merge step is responsible for replacing these placeholders with real `compare -metric AE -fuzz 5%` output once both Plan 02 (X11 baselines) and Plan 03 (Qt 1x captures) are merged into the same branch.

## Phase boundary check

- `git diff --quiet -- xvue/qt/src/` → exit 0 (byte-identical confirmed)
- `git diff --quiet -- xvue/xvuelc.c` → exit 0 (byte-identical confirmed)

## Manifest

| Case | Size | Dims | SHA-256 | AE | AE% | Verdict | Notes |
|------|------|------|---------|----|----|---------|-------|
| pan2d | 70944 | 760x442 | 2bf53484d63853822a90f43b08a8af986f07fee767ddaf082209feedeb2615db | (pending) | (pending) | PENDING-BASELINE | mesher post-mesh frame (FAUX flag IEEE_DENORMAL, capture before STOP) |
| nafems_le1 | 320770 | 760x442 | 76fb223667d76b270a56af540d98ccedb47fba54c764518a47a9b39d8e664279 | (pending) | (pending) | PENDING-BASELINE | elasticity stress visualization (gradient case — Pitfall 6 candidate for fuzz_override_pct=8 in CHECKLIST.md) |
| cavity2d | 44616 | 760x442 | 37739e5c0a2942d62c5f0b0731bd549c740f166ea8bc5107309cc9fe6583d86c | (pending) | (pending) | PENDING-BASELINE | fluid lid-driven cavity stream lines (gradient case — Pitfall 6 candidate) |
| heat1d | 21834 | 760x442 | 6307f9b58151e358078f12afd0181628e027975f8f855d7c397c02ac9207b14d | (pending) | (pending) | PENDING-BASELINE | thermal 1D unsteady (gradient case — Pitfall 6 candidate) |
| nlsecu | 70612 | 760x442 | 8877725fb88a4d4386a031e7c28ccf658d7a7f4b5f0897f9385747c70cd53ba8 | (pending) | (pending) | TRUNCATED-CAPTURE | MAILLER prereq frame (nlsecu.meshq2) — NLSER-frame unreachable in 60s/120s budgets due to documented Plan 1 deadlock |

## Verdict roll-up

- PASS verdicts: 0/5 (compare not yet run)
- CHECK verdicts: 0/5 (compare not yet run)
- FAIL-HARNESS verdicts: 0/5 (no harness errors — captures all valid PNGs at 760x442)
- PENDING-BASELINE verdicts: 4/5 (pan2d, nafems_le1, cavity2d, heat1d)
- TRUNCATED-CAPTURE verdicts: 1/5 (nlsecu — see TRUNCATED-CAPTURE legend above)

The 5/5 row-count is satisfied; the AE column requires Plan 02 baselines to compute and is filled by the orchestrator wave-merge step.

## Hand-off to Plan 07

**5 absolute paths to ${CASE}-qt-1x.png:**
- `/home/mefisto/git/mefisto/.planning/phases/08-ab-validation-on-testa-subset/evidence/pan2d-qt-1x.png`
- `/home/mefisto/git/mefisto/.planning/phases/08-ab-validation-on-testa-subset/evidence/nafems_le1-qt-1x.png`
- `/home/mefisto/git/mefisto/.planning/phases/08-ab-validation-on-testa-subset/evidence/cavity2d-qt-1x.png`
- `/home/mefisto/git/mefisto/.planning/phases/08-ab-validation-on-testa-subset/evidence/heat1d-qt-1x.png`
- `/home/mefisto/git/mefisto/.planning/phases/08-ab-validation-on-testa-subset/evidence/nlsecu-qt-1x.png`

**5 absolute paths to ${CASE}-qt-1x-diff.png (placeholders pending wave-merge):**
- `/home/mefisto/git/mefisto/.planning/phases/08-ab-validation-on-testa-subset/evidence/pan2d-qt-1x-diff.png`
- `/home/mefisto/git/mefisto/.planning/phases/08-ab-validation-on-testa-subset/evidence/nafems_le1-qt-1x-diff.png`
- `/home/mefisto/git/mefisto/.planning/phases/08-ab-validation-on-testa-subset/evidence/cavity2d-qt-1x-diff.png`
- `/home/mefisto/git/mefisto/.planning/phases/08-ab-validation-on-testa-subset/evidence/heat1d-qt-1x-diff.png`
- `/home/mefisto/git/mefisto/.planning/phases/08-ab-validation-on-testa-subset/evidence/nlsecu-qt-1x-diff.png`

**5 AE pixel counts (one per case):** all PENDING — see Wave-Merge Compare Recipe below.

## Wave-Merge Compare Recipe (for orchestrator post-merge step)

After Plan 02 baselines and Plan 03 captures are both merged into the same
branch, the orchestrator (or a follow-up plan) runs the per-case compare
and replaces the placeholder diff PNGs:

```bash
EVD=.planning/phases/08-ab-validation-on-testa-subset/evidence
for c in pan2d nafems_le1 cavity2d heat1d nlsecu; do
  bin/ab_compare_pair.sh \
    "$EVD/${c}-x11.png" \
    "$EVD/${c}-qt-1x.png" \
    "$EVD/${c}-qt-1x-diff.png" \
    5
done
```

The output rows of `ab_compare_pair.sh` (one per case) populate the per-case
AE/AE%/Verdict cells in the "Per-case row table" + "Manifest" sections above.
The orchestrator MUST update the row table, manifest, and verdict roll-up
during the wave-merge commit so Plan 07 ingests final values, not placeholders.

## Notes for maintainer review

- **pan2d:** mesher case — capture is a post-mesh frame from `pp/ppmail_qt`. Likely Pitfall 7 (font AA budget ~3000 px); if AE exceeds the global 5% gate but stays below ~3000 px, candidate for accept-with-note in CHECKLIST.md.
- **nafems_le1:** elasticity stress contour — gradient-rich; Pitfall 6 candidate. If AE/AE% indicates concentration in the color-bar region, candidate for `fuzz_override_pct=8` per the deferred idea in CONTEXT.md.
- **cavity2d:** fluid stream-line plot — gradient-rich; Pitfall 6 candidate. Same fuzz-override consideration as nafems_le1.
- **heat1d:** 1D thermal — small image; high AE% sensitivity to single-pixel drift expected.
- **nlsecu:** captured frame is from MAILLER prereq, NOT NLSER. Plan 07 (and CHECKLIST.md) MUST distinguish this row as TRUNCATED-CAPTURE in the Qt 1x cell. Real NLSER A/B sign-off requires either (a) extending per-case timeout for nlsecu, (b) a code-side fix to the offscreen-Qt deadlock at `pp/ppnlse_qt` startup, or (c) explicit waiver-by-rationale in 08-CHECKLIST.md per the open-items list in `00-smoke-probes.md`.

## Outcome

Plan 03 Qt 1x column complete: 5/5 captures (4 normal, 1 TRUNCATED-CAPTURE for nlsecu),
0/5 PASS-at-default-fuzz (compare deferred to wave-merge), 4/5 PENDING-BASELINE
pending Plan 02 merge + orchestrator re-compare. Hand-off ready for Plan 07
once wave-merge step populates the AE/Verdict cells.
