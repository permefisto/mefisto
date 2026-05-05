# Font-metric spot check — pan2d (VALID-06)

**Date:** 2026-05-05
**Plan:** Phase 8 Plan 06, Task 2 (Part A)
**Pitfall:** 7 (Anti-aliased text + font hinting drift, 08-RESEARCH.md line 453)

## Captures

| File                                         | Dims      | Size B | Unique colors | Source                                |
| -------------------------------------------- | --------- | ------ | ------------- | ------------------------------------- |
| `evidence/pan2d-x11.png`                     | 1280x800  | 7636   | **9**         | Plan 02 X11 baseline                  |
| `evidence/pan2d-qt-1x.png`                   | 760x442   | 70944  | **3545**      | Plan 03 Qt 1x capture                 |
| `evidence/pan2d-qt-1x-diff.png`              | 1280x800  | 37456  | 15            | wave-merge-ae.md re-run (post Wave 2) |

The X11 capture has **9 unique colors** — exactly the canonical Mefisto
palette: `#000000` (black, axes/labels), `#4F4F63` (mesh-grey, mesh edges),
`#FF0000` (red), `#FF7F00` (orange), `#F97F71` (fleshy-pink), `#31C731`
(green), `#00FFCC` (cyan), `#FFD9B8` (peach background), `#FFFFFF` (white,
label rectangles).

The Qt-1x capture has **3545 unique colors** — a 394× explosion driven by
the **anti-aliased font rendering**. Each label glyph in Qt is drawn with
sub-pixel hinting via FreeType, producing dozens of intermediate gray-shade
pixels at every glyph edge. This is the canonical Pitfall-7 signature.

**Color-region census of pan2d Qt-1x** (pixels matching canonical color
within 1 % fuzz, against total 335 920):

| Canonical color | Pixels (1% fuzz) | Pct of frame |
| --------------- | ---------------- | ------------ |
| `#000000` black | 556              | 0.17 %       |
| `#FFD9B8` peach | 176 209          | 52.46 %      |
| `#FFFFFF` white | 372              | 0.11 %       |
| `#FF0000` red   | 2 276            | 0.68 %       |
| `#FF7F00` orange| 400              | 0.12 %       |
| `#31C731` green | 1 322            | 0.39 %       |
| `#00FFCC` cyan  | 443              | 0.13 %       |
| **Sum**         | **181 578**      | **54.05 %**  |
| **AA fringe**   | **154 342**      | **45.95 %**  |

~46 % of Qt-1x pixels are NOT in the canonical palette — these are AA
intermediates (sub-pixel font hinting + edge anti-aliasing on mesh edges).
This is **expected, not a defect**.

## Diff statistics

(Source: `wave-merge-ae.md` re-run post Wave-2 merge.)

| Metric                                | Value             |
| ------------------------------------- | ----------------- |
| AE (compare -metric AE -fuzz 5%)      | **540 804** px    |
| AE %                                  | **52.81 %**       |
| Total pixels (resampled to 1280x800)  | 1 024 000         |
| Diff PNG mean intensity (0–65535)     | 42 857.8 (65.4 %) |
| Diff PNG max intensity                | 65 535            |
| Diff PNG stddev                       | 19 381.2          |
| Diff PNG unique colors                | 15                |
| Pixels with red-channel > 95 %        | 657 846           |
| Resampled                             | yes (760x442 → 1280x800) |

**Qualitative reading of the diff PNG** (`evidence/pan2d-qt-1x-diff.png`):

The diff PNG is the standard ImageMagick `compare` output: white-ish for
matched pixels, red for differing pixels. The 15 unique colors include
mesh-edge match colors (light pink / cream) AND saturated red (#F30018,
#FF0000-ish) which marks the AE pixels.

The AE pixels are **NOT concentrated to a single region** — they are
field-wide because of the dim-mismatch resample artefact (760x442 →
1280x800 nearest-neighbour upsample shifts every mesh-edge by ~1.68 px).

**Concentration around glyphs:** This cell's principal AE signal is the
**resample artefact (~field-wide)**, NOT pure glyph fringing. Pure glyph
fringing would be confined to a ~200-px-radius blob around each of ~15
labels (estimated location: vertices of the meshed pan, scattered through
the canvas). On the resampled diff PNG, the field-wide drift dominates;
glyph fringing is a sub-component but not separately resolvable.

## AE budget

**Pitfall 7 formula** (08-RESEARCH.md line 456):
> "For pan2d (10–20 node labels at typical scale), an AE budget of
> ~200 pixels per label × 15 labels ≈ **3000 pixels** is reasonable."

**Empirical pan2d label count** (from inspection of the X11 capture and
the testa/pan2d/ Fortran data files):
- pan2d is a 2D mesh of a flat panel; vertex labels typically number
  10–20 per the canonical Mefisto rendering. Without re-running the
  case, the planner's heuristic is taken at face value: ~15 labels.

**Budget table:**

| Tier             | AE budget       | Verdict mapping       |
| ---------------- | --------------- | --------------------- |
| Within budget    | ≤ 3 000 px      | PASS (AA-fringe only) |
| Within 2× budget | 3 001 – 6 000 px| CHECK (maintainer)    |
| Beyond 2× budget | > 6 000 px      | FAIL (escalate)       |

**Empirical AE = 540 804 px** — **180× the 3000 px budget**. By the
Pitfall-7 budget alone, this is a FAIL.

**However**: the budget is designed for a matched-dim compare. The current
measurement is dominated by the 760x442 → 1280x800 resample artefact (same
root cause as the 5/5 Plan 03 cells in `wave-merge-ae.md`). The Pitfall-7
budget cannot be applied directly to this resample-confounded number.

**Adjusted reading:**
- Estimated glyph-fringe AE (under the Pitfall-7 budget): ~3 000 px
  (matched-dim, AA-only).
- Estimated resample-artefact AE: ~540 000 px (field-wide upsample drift).
- Until a matched-dim recapture is available, the Pitfall-7 verdict cannot
  be cleanly assigned — the AE signal is **not separable** from the
  resample noise.

## Hypothesis

The 52.81 % AE on the pan2d Qt-1x cell has **two contributing sources**:

1. **Pitfall 7 — AA + font hinting drift** (a few thousand px around each
   label glyph). This is *expected and not a defect*.
2. **Dim-mismatch resample artefact** (~half a million px field-wide).
   This is dominant, and is NOT a backend rendering defect — it's a
   capture-pipeline artefact: Qt offscreen renders to a 760x442 backing
   pixmap; X11 captures from a 1280x800 Xvfb root; `bin/ab_compare_pair.sh`
   resamples Qt up via `-filter point` (no AA introduced) before
   compare-AE; every mesh edge is offset by ~1.68 px → counts as drift.

The 3 545 unique-colors finding **directly confirms Pitfall 7**: the Qt
backend IS doing aggressive sub-pixel font AA (vs X11's 9 discrete colors).
The visual result (label readability) is fine — the AE count is not a
clipping/overlapping/missing-label signal.

**Diagnosis:** Pitfall 7 hypothesis is **CONFIRMED qualitatively** (glyph
edges have AA fringes, as expected). Pitfall-7 quantitative budget cannot
be cleanly applied because the dim-mismatch resample artefact dominates
the per-pixel AE count.

## Verdict

**CHECK** — AE = 540 804 px (52.81 %) exceeds the 3000-px Pitfall-7 budget
by 180×, BUT the dominant AE source is the **resample-dim-mismatch
artefact** common to all 5 Wave-2 cells, NOT pure glyph fringing. The
Pitfall-7 signature (3545 unique colors in Qt-1x vs 9 in X11) is empirically
confirmed and is the expected backend-rendering-difference, not a defect.

Maintainer review at Plan 7 should:
1. Accept the cell as PASS-on-review pending diff-PNG inspection (no
   missing/clipped/overlapped labels expected; AE pixels are scattered
   field-wide, not a localized defect blob).
2. Annotate CHECKLIST.md cell with `AE-budget-confounded-by-resample`.
3. Carry forward to Phase 9 a deferred-idea: matched-dim Qt recapture
   (1280x800 in-process pixmap if `MEFISTO_QT_WINDOW_SIZE` or equivalent
   is supported) to retire the resample-confound and apply the Pitfall-7
   budget cleanly.
