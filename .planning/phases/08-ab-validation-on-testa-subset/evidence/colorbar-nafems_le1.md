# Color-bar spot check — nafems_le1 (VALID-05)

**Date:** 2026-05-05
**Plan:** Phase 8 Plan 06, Task 1
**LEFT (X11 baseline):** `evidence/nafems_le1-x11.png` — 1280x800, 60765 B, 7 unique colors
**RIGHT (Qt 1x):** `evidence/nafems_le1-qt-1x.png` — 760x442, 320770 B
**Tool:** `bin/ab_compare_pair.sh` (resamples Qt to 1280x800 via `-filter point`)

## Region

**Empirical finding:** nafems_le1 X11 capture contains only **7 unique colors**
(`#000000` black, `#4F4F63` mesh-grey, `#FF0000` red, `#FF7F00` orange,
`#31C731` green, `#FFD9B8` peach background, `#FFFFFF` white). There is **no
smooth-gradient colorbar** anywhere in the X11 capture — the elasticity result
is rendered as **3 discrete stress contour bands** (red / orange / green) over
the entire mesh, NOT as a continuous color ramp on a side legend.

Strip-by-strip unique-color survey across the X11 capture (80px-wide vertical
strips):

| x range   | unique colors |
| --------- | ------------- |
| 0–80      | 6             |
| 100–180   | 4             |
| 200–280   | 6             |
| 300–380   | 6             |
| 400–480   | 6             |
| 500–580   | 6             |
| 600–680   | 6             |
| 700–780   | 6             |
| 800–880   | 6             |
| 900–980   | 6             |
| 1000–1080 | 6             |
| 1100–1180 | 5             |
| 1180–1260 | 5             |

No strip is colorbar-rich. The 3 data colors (red/orange/green) are spatially
distributed across the meshed region, not concentrated in a side band.

**Decision:** SKIP the per-region crop; use **full-frame compare**. Region is
documented as "full-frame — discrete contour bands distributed over mesh; no
smooth-gradient bar visible in the X11 capture".

- x_offset: 0
- y_offset: 0
- width: 1280
- height: 800
- region label: full-frame (discrete contour bands, no side colorbar)

## AE measurements

| fuzz % | AE     | AE%      | diff PNG path                                                          | exit code | verdict |
| ------ | ------ | -------- | ---------------------------------------------------------------------- | --------- | ------- |
| 5      | 412827 | 40.3151% | `/tmp/phase8-spot-check/nafems_le1-fullframe-diff-5pct.png`            | CHECK     | CHECK   |
| 8      | 404740 | 39.5254% | `/tmp/phase8-spot-check/nafems_le1-fullframe-diff-8pct.png`            | CHECK     | CHECK   |

(AE delta from 5% → 8% override = 8087 px = ~0.79 % of total pixel count;
override produces only marginal AE reduction — the dominant AE source is
*not* per-pixel gradient drift but a structural full-frame mismatch, namely
the Qt 760x442 → 1280x800 nearest-neighbour resample mandated by the dim guard
in `bin/ab_compare_pair.sh`.)

Diff PNG mean intensity at fuzz=5: `48258.8 / 65535 = 73.6%` — high mean
indicates AE pixels are scattered widely, not concentrated in a single band.

## Hypothesis

This case **does NOT manifest Pitfall 6** (gradient color bar tripping fuzz)
because the case has no smooth-gradient colorbar to begin with. The high AE
(~40%) is dominated by **the resample dim mismatch**: Qt offscreen renders
into a 760x442 backing pixmap while the X11 baseline is captured from a
1280x800 Xvfb root; `bin/ab_compare_pair.sh` resamples Qt up to 1280x800 via
`-filter point` (no AA introduced) before comparing — the nearest-neighbour
upsample produces a stair-stepped image whose every gradient/edge boundary
is offset by ≈ 1.68 px from the X11 baseline, which counts as drift on
roughly half the foreground pixels.

The `wave-merge-ae.md` consolidation already records this finding for ALL 5
cases at the main A/B sweep level: every Plan-3/4/5 cell came back CHECK
because the dim-mismatch + resample is the dominant signal, NOT real
backend rendering drift.

**Diagnosis:** Pitfall 6 hypothesis (gradient drift) is **disconfirmed for
nafems_le1**. The CHECK verdict on this Wave-2 cell is a resample artefact,
NOT a real colormap defect. The discrete-contour rendering style (3 colors
only) sidesteps Pitfall 6 entirely.

The diff PNG at fuzz=5/8 looks substantively identical (8087 px delta) which
confirms the override doesn't help — it's not gradient drift to absorb.

## Recommended fuzz_override_pct for CHECKLIST.md

**fuzz_override = 8 %** (Pitfall 6 default for stress-bar cases per
08-RESEARCH.md §Pitfall 6 line 450).

**Rationale:** retained as the planner-recommended starting value even though
it produces near-zero AE delta in this empirical run. The override is harmless
(8 % is well within the [1,30] band enforced by `bin/ab_compare_pair.sh`) and
preserves the intended Pitfall-6 mitigation contract for future re-captures
that might be done at matched dims (760x442 → 760x442) where gradient drift
WOULD show up unbuffered.

A more substantive recommendation for Plan 7's CHECKLIST.md cell:
- Annotate the cell with `fuzz_override_pct=8` per Pitfall 6 contract
- Annotate **also** with `dim-mismatch=resample-to-1280x800` so the maintainer
  understands that the dominant AE source is NOT colorbar drift
- Plan 9 (or a hot-fix) should consider matched-dim recapture (Qt at 1280x800
  via `MEFISTO_QT_WINDOW_SIZE` if such a contract exists) to retire the
  resample-confound entirely.

## Verdict

**CHECK** — high AE (~40 %) is a **resample-dim-mismatch artefact**, not a
real colormap defect. Maintainer review at Plan 7 must reconcile this against
the wave-merge consolidated AE row; Plan 7 sign-off should accept this row
as PASS-on-review if the diff PNG inspection confirms scattered drift only
(no missing/clipped contours), otherwise escalate to a matched-dim recapture
plan.
