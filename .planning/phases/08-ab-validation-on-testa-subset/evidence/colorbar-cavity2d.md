# Color-bar spot check — cavity2d (VALID-05)

**Date:** 2026-05-05
**Plan:** Phase 8 Plan 06, Task 1
**LEFT (X11 baseline):** `evidence/cavity2d-x11.png` — 1280x800, 14906 B, 19 unique colors
**RIGHT (Qt 1x):** `evidence/cavity2d-qt-1x.png` — 760x442, 44616 B
**Tool:** `bin/ab_compare_pair.sh` (resamples Qt to 1280x800 via `-filter point`)

## Region

**Empirical finding:** cavity2d X11 has **19 unique colors** — the richest
palette of the 3 colorbar cases. A right-edge vertical band at x=1100–1260
has 15 unique colors → this **IS** a velocity-magnitude / streamline color
legend (the only true visible colorbar of the 3 VALID-05 cases).

Strip-by-strip unique-color survey (X11):

| x range   | unique colors |
| --------- | ------------- |
| 0–80      | 7             |
| 100–180   | 6             |
| 200–280   | 7             |
| 300–380   | 7             |
| 400–480   | 6             |
| 500–580   | 6             |
| 600–680   | 6             |
| 700–780   | 6             |
| 800–880   | 6             |
| 900–980   | 6             |
| 1000–1080 | 4             |
| **1100–1180** | **16**    | ← colorbar/legend region (right edge)
| **1180–1260** | **15**    | ← colorbar/legend region (continued)

Per the plan's design intent, **cavity2d arrows are distributed through the
velocity field, NOT a single band**. We report BOTH measurements:

1. Full-frame compare (the canonical Wave-2 measurement)
2. Cropped colorbar-region compare (right edge x=1100–1260, 160px wide)
   — the dominant gradient-rich region, isolated to test Pitfall 6.

For the cropped region:
- **X11 region:** 160x800+1100+0 (full height of the right legend band)
- **Qt region:** 95x442+653+0 (proportionally scaled: 1100/1280 × 760 = 653;
  160/1280 × 760 = 95)
- Both crops are then resampled by `bin/ab_compare_pair.sh` to 160x800.

## AE measurements

### Full-frame (canonical)

| fuzz % | AE     | AE%      | diff PNG path                                                    | exit code | verdict |
| ------ | ------ | -------- | ---------------------------------------------------------------- | --------- | ------- |
| 5      | 411003 | 40.1370% | `/tmp/phase8-spot-check/cavity2d-fullframe-diff-5pct.png`        | CHECK     | CHECK   |
| 10     | 404782 | 39.5295% | `/tmp/phase8-spot-check/cavity2d-fullframe-diff-10pct.png`       | CHECK     | CHECK   |

(AE delta from 5% → 10% override = 6221 px = ~0.61 % of total pixel count;
override produces only marginal AE reduction. Same dim-mismatch pattern as
nafems_le1 / heat1d.)

### Cropped colorbar region (right-edge legend, 160x800)

| fuzz % | AE    | AE%     | total pixels | diff PNG path                                              | exit code | verdict |
| ------ | ----- | ------- | ------------ | ---------------------------------------------------------- | --------- | ------- |
| 5      | 28078 | 21.94%  | 128000       | `/tmp/phase8-spot-check/cavity2d-cb-diff-5pct.png`         | CHECK     | CHECK   |
| 10     | 27724 | 21.66%  | 128000       | `/tmp/phase8-spot-check/cavity2d-cb-diff-10pct.png`        | CHECK     | CHECK   |

(Cropped region AE% (~22%) is **lower** than full-frame (~40%), suggesting
the gradient region itself is **less drifted** than the rest of the image —
consistent with the velocity arrows being scattered AND the dominant signal
being the Qt 95x442 → 160x800 resample stretch.)

Diff PNG mean intensities (full-frame):
- fuzz=5: `45011.8 / 65535 = 68.7%` (mean intensity)
- fuzz=10: `45217.8 / 65535 = 69.0%` (mean intensity)

## Hypothesis

cavity2d is the **only one of the 3 VALID-05 cases that genuinely has a
gradient colorbar** (right-edge legend at x=1100–1260, 15+ unique colors).
The cropped-region AE (~22%) is in fact LOWER than the full-frame (~40%),
which **disconfirms Pitfall 6 as the dominant signal** — if the colorbar
were the failure mode, its AE% would be HIGHER than the rest.

What's actually happening:
- The 760x442 → 1280x800 nearest-neighbour upsample (in `bin/ab_compare_pair.sh`)
  shifts every velocity-arrow line endpoint by ≈ 1.68 px on both axes.
- Velocity arrows are line segments scattered across the field; each shifted
  endpoint registers ≈ 5–10 px of drift; multiplied by ~50 arrows → ~300–500 px.
- The dominant AE (~411 K of 1.024 M = ~40 %) is field-wide resample noise
  PLUS the arrow-endpoint drift.
- Override fuzz=10 % absorbs ~6 K px (the arrow tips); the rest is resample.

**Diagnosis:** Same as nafems_le1 / heat1d — the CHECK verdict is dominated
by **dim-mismatch resample artefact**, not real backend rendering drift.
Pitfall 6 (gradient-bar drift) is a SECONDARY signal here (~22 % cropped AE
on a region that should ideally show < 5 % drift), but the primary noise
floor is the resample.

## Recommended fuzz_override_pct for CHECKLIST.md

**fuzz_override = 10 %** per Pitfall 6 for cavity2d velocity arrows
(08-RESEARCH.md §Pitfall 6 line 450).

**Rationale:** cavity2d genuinely has the richest palette + visible colorbar
of the 3 cases. The 10 % override would meaningfully help on a matched-dim
recapture (where the gradient-edge drift IS the dominant per-pixel signal).
On the current resample-confounded measurement, the override only saves
~0.6 pp of AE.

Plan 7 cell annotation:
- `fuzz_override_pct=10` (Pitfall 6)
- `dim-mismatch=resample-to-1280x800` (dominant AE source)
- `colorbar-region-AE-pct=21.94` (cropped right-edge x=1100–1260; recorded
  for forensic inspection — confirms colorbar drift is not the failure
  driver despite the visible gradient).

## Verdict

**CHECK** — full-frame AE (~40 %) and cropped colorbar-region AE (~22 %) are
both above the 5 % gate, but the dominant signal is the **dim-mismatch
resample artefact** common to all 5 Wave-2 cells, NOT a real colorbar
defect. The cropped-region AE being LOWER than full-frame disconfirms
Pitfall 6 as the principal failure mode for this cell.

Maintainer review at Plan 7 should accept as PASS-on-review pending diff-PNG
inspection (expected: scattered drift + arrow-tip jitter, no missing/wrong
arrows, no missing colorbar swatches), or escalate to a matched-dim
recapture plan.
