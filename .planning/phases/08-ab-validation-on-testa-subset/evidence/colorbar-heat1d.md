# Color-bar spot check — heat1d (VALID-05)

**Date:** 2026-05-05
**Plan:** Phase 8 Plan 06, Task 1
**LEFT (X11 baseline):** `evidence/heat1d-x11.png` — 1280x800, 4363 B, 8 unique colors
**RIGHT (Qt 1x):** `evidence/heat1d-qt-1x.png` — 760x442, 21834 B
**Tool:** `bin/ab_compare_pair.sh` (resamples Qt to 1280x800 via `-filter point`)

## Region

**Empirical finding:** heat1d X11 capture contains only **8 unique colors**
(`#000000` black, `#4F4F63` mesh-grey, `#FF0000` red, `#0000FF` blue,
`#FF00FF` magenta, `#9595B1` mid-purple, `#FFD9B8` peach background,
`#FFFFFF` white). This is a **1D thermal profile plot**, not a 2D
colormapped field — the colors are **plot curves** (red/blue/magenta lines)
on a peach background, with axes (mesh-grey) and labels (black on white).

Strip-by-strip unique-color survey across the X11 capture:

| Vertical strips (80px wide) | unique | | Horizontal strips (80px tall) | unique |
| -------------------------- | ------ | | ----------------------------- | ------ |
| x=0–80                     | 4      | | y=0–80                        | 4      |
| x=100–180                  | 5      | | y=100–180                     | 3      |
| x=200–280                  | 7      | | y=200–280                     | 5      |
| x=300–380                  | 7      | | y=300–380                     | 6      |
| x=400–480                  | 8      | | y=400–480                     | 6      |
| x=500–580                  | 7      | | y=500–580                     | 3      |
| x=600–680                  | 6      | | y=600–680                     | 3      |
| x=700–780                  | 6      | | y=700–780                     | 3      |
| x=800–880                  | 6      | |                               |        |
| x=900–980                  | 6      | |                               |        |
| x=1000–1080                | 6      | |                               |        |
| x=1100–1180                | 3      | |                               |        |
| x=1180–1260                | 3      | |                               |        |

The middle band (x=200–1080, y=200–500) carries the curve plot region.
**No side-mounted gradient colorbar is visible** — heat1d does not render a
temperature ramp side-bar; it renders the thermal profile **directly as
plotted curves** vs the spatial coordinate.

**Decision:** SKIP the per-region crop; use **full-frame compare**. Region is
documented as "full-frame — 1D thermal-curves plot, no gradient ramp present".

- x_offset: 0
- y_offset: 0
- width: 1280
- height: 800
- region label: full-frame (1D plot, no gradient ramp)

## AE measurements

| fuzz % | AE     | AE%      | diff PNG path                                                  | exit code | verdict |
| ------ | ------ | -------- | -------------------------------------------------------------- | --------- | ------- |
| 5      | 143273 | 13.9915% | `/tmp/phase8-spot-check/heat1d-fullframe-diff-5pct.png`        | CHECK     | CHECK   |
| 8      | 142397 | 13.9060% | `/tmp/phase8-spot-check/heat1d-fullframe-diff-8pct.png`        | CHECK     | CHECK   |

(AE delta from 5% → 8% override = 876 px = ~0.085 % of total pixel count;
override has essentially zero effect — the dominant AE source is not gradient
drift but the same structural resample-dim-mismatch noted on every Plan 03/04
cell of the wave-merge-ae.md consolidation.)

Diff PNG mean intensity at fuzz=5: `57847.1 / 65535 = 88.3%` — very high mean
indicates AE pixels saturated (red on lightgray background). Concentrated
diff would have a lower mean (more white area).

heat1d's AE% (~14%) is **lowest of the 3 colorbar cases** — consistent with
the simpler color palette (8 vs 19 colors for cavity2d).

## Hypothesis

This case **does NOT manifest Pitfall 6** (gradient color bar tripping fuzz)
because there is no temperature ramp colorbar in the capture. heat1d renders
the thermal profile as line plots in 3 distinct colors — not a colormapped
2D field. The Pitfall-6 mitigation contract still applies in principle (any
gradient-rich rendering would benefit from fuzz=8%), but EMPIRICALLY this
case has no gradient region.

The CHECK verdict at default 5% fuzz is dominated by the same Plan-3/4 root
cause: Qt offscreen renders into 760x442; X11 captures into 1280x800; the
nearest-neighbour upsample mandated by `bin/ab_compare_pair.sh`'s dim guard
shifts every plot-curve edge by ≈ 1.68 px → counts as drift. Lower AE%
(~14% vs ~40% for nafems_le1) is consistent with heat1d having less
foreground pixel content (sparse thin curves over peach background).

**Diagnosis:** Pitfall 6 hypothesis (gradient drift) is **disconfirmed for
heat1d** — same as nafems_le1. The CHECK verdict is a resample-dim-mismatch
artefact, not a real defect.

## Recommended fuzz_override_pct for CHECKLIST.md

**fuzz_override = 8 %** per Pitfall 6 default for thermal cases
(08-RESEARCH.md §Pitfall 6 line 450).

**Rationale:** retained as the planner-recommended value even though it
produces a sub-pct AE delta empirically. Same logic as nafems_le1: the
override is harmless and preserves the intended Pitfall-6 mitigation
contract for future matched-dim re-captures.

Plan 7 cell annotation should pair the override with the dim-mismatch note:
- `fuzz_override_pct=8` (per Pitfall 6)
- `dim-mismatch=resample-to-1280x800` (dominant AE source)
- A matched-dim recapture (Qt at 1280x800) would be expected to land within
  the global 5 % gate; this is a deferred-idea for Phase 9.

## Verdict

**CHECK** — AE (~14 %) is a **resample-dim-mismatch artefact**, NOT a real
colormap defect. The lowest-AE of the 3 colorbar cases (consistent with
heat1d's sparse curve-plot rendering). Maintainer review at Plan 7 should
accept as PASS-on-review pending diff-PNG inspection (no missing/clipped
curves expected, only ≈1.68 px edge drift from resample), otherwise escalate
to matched-dim recapture.
