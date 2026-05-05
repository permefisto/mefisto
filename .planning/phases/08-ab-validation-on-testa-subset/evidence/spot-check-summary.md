# Spot-Check Summary (VALID-05 + VALID-06)

**Date:** 2026-05-05
**Plan:** Phase 8 Plan 06, Task 3
**Source reports:** 5 spot-check markdown files + 1 hexahedron probe report
under `.planning/phases/08-ab-validation-on-testa-subset/evidence/`.

## Aggregated Verdicts

| Spot Check          | Type      | Cases       | Region                                | AE     | AE%      | Recommended override               | Verdict         | Report                                   |
|---------------------|-----------|-------------|---------------------------------------|--------|----------|------------------------------------|-----------------|------------------------------------------|
| colorbar-nafems_le1 | Pitfall 6 | nafems_le1  | full-frame (no smooth gradient bar)   | 412827 | 40.3151% | fuzz=8 + dim-mismatch annotation   | CHECK           | colorbar-nafems_le1.md                   |
| colorbar-heat1d     | Pitfall 6 | heat1d      | full-frame (1D plot, no ramp)         | 143273 | 13.9915% | fuzz=8 + dim-mismatch annotation   | CHECK           | colorbar-heat1d.md                       |
| colorbar-cavity2d   | Pitfall 6 | cavity2d    | full-frame + cropped right-edge x=1100-1260 (160x800) | 411003 | 40.1370% | fuzz=10 + dim-mismatch annotation  | CHECK           | colorbar-cavity2d.md                     |
| font-pan2d          | Pitfall 7 | pan2d       | full-frame                            | 540804 | 52.8129% | AE budget=3000 (confounded by dim-mismatch) | CHECK | font-pan2d.md                            |
| font-hexahedron     | Pitfall 7 | hexahedron  | N/A (probe = GAP)                     | N/A    | N/A      | AE budget=1600 (reserved, unapplied) | GAP-DOCUMENTED | font-hexahedron.md + hexahedron-probe.md |

**Roll-up counts:** PASS=0, CHECK=4, FAIL=0, GAP-DOCUMENTED=1.

## VALID-05 Roll-up

- **Coverage:** 3/3 color-bar cases verified (nafems_le1, heat1d, cavity2d).
- **Verdict counts:** PASS=0, CHECK=3, FAIL=0.
- **Common finding:** All 3 cells exhibit a dominant AE signal driven by the
  **760x442 → 1280x800 nearest-neighbour resample artefact** (Qt offscreen
  backing pixmap dim vs Xvfb root capture dim mismatch), NOT real backend
  rendering drift. The Pitfall 6 hypothesis (gradient color bar tripping
  the 5% fuzz gate) is **disconfirmed** for all 3 cases:
  - **nafems_le1** has only 7 unique colors → discrete contour bands, no
    gradient bar exists in the X11 capture.
  - **heat1d** has only 8 unique colors → 1D thermal-curves plot, no
    temperature-ramp side legend.
  - **cavity2d** has 19 unique colors with a real right-edge legend at
    x=1100–1260 (15 unique in that band), BUT the cropped colorbar-region
    AE (21.94%) is LOWER than the full-frame AE (40.14%), disconfirming
    Pitfall 6 as the principal failure mode.

- **Recommended fuzz overrides** (for Plan 7 to encode in CHECKLIST.md
  cell-Notes column):
  | Case        | fuzz_override_pct | Pair-with annotation                      |
  |-------------|-------------------|-------------------------------------------|
  | nafems_le1  | 8                 | `dim-mismatch=resample-to-1280x800`       |
  | heat1d      | 8                 | `dim-mismatch=resample-to-1280x800`       |
  | cavity2d    | 10                | `dim-mismatch=resample-to-1280x800` + `colorbar-region-AE-pct=21.94` |

- **Deferred idea (Phase 9 candidate):** matched-dim Qt recapture (1280x800
  in-process backing pixmap if `MEFISTO_QT_WINDOW_SIZE` or equivalent is
  supported) to retire the resample-confound and apply Pitfall 6 fuzz
  overrides cleanly to a like-for-like compare.

## VALID-06 Roll-up

- **Coverage:** 1/2 (pan2d only). hexahedron is GAP-DOCUMENTED — see
  escalation in `08-VALIDATION.md ## Open escalations`.
- **Verdict counts:** PASS=0, CHECK=1, FAIL=0, GAP-DOCUMENTED=1.
- **pan2d Pitfall 7 finding:** Qt-1x has 3545 unique colors vs X11's 9 →
  394× explosion → empirically **CONFIRMS** Pitfall 7 (anti-aliased text
  + font hinting drift). The pan2d cell's AE (540 804 px = 52.81%) is
  dominated by the same dim-mismatch resample artefact that drives the 3
  color-bar cells. Pitfall 7's 3000-px AE budget cannot be applied cleanly
  on the resample-confounded measurement.
- **hexahedron probe outcome:** GAP. All 5 alphabetical candidates
  ({ln, ob, pt, sf, vl}) probed; 4 crashed at unknown-name lookup with
  2045-byte all-peach background captures (xvfermer_ fired before any
  drawing primitive painted the mesh); 1 (`pt`) hung at the 60-s timeout
  (no rendering trigger; expects interactive menu pick to draw). The 5
  files are interlocking building blocks (pt → ln → sf → vl → ob), NOT
  independent single-file batch entry points. The legacy MAILLER protocol's
  single-batch invocation pattern (`pp/ppmail batchfile`) does not apply.

- **Recommended AE budgets** (for Plan 7 to encode in CHECKLIST.md):
  | Case        | Pitfall-7 AE budget | Status                                   |
  |-------------|---------------------|------------------------------------------|
  | pan2d       | 3000 px             | confounded by dim-mismatch resample      |
  | hexahedron  | 1600 px (reserved)  | unapplied — gap-documented per probe     |

- **Cross-reference:** `08-VALIDATION.md ## Open escalations` records the
  hexahedron escalation with the maintainer-decision options:
  (a) CONTEXT.md amendment naming the canonical driver (wrapper script,
  stdin-piped menu sequence, or new multi-file batch protocol), OR
  (b) waiver-by-rationale (e.g., "pan2d coverage is sufficient for AA
  drift gate"). Plan 7 carries this forward into the v1 ship-gate sign-off
  matrix as a maintainer-decision item.

## Hand-off to Plan 07

### Per-main-case CHECKLIST.md cell annotations

For each main-case canonical row (the 5 BUILD-10 baseline cases), Plan 7
encodes the following annotations in the matching CHECKLIST.md cell's
Notes column:

| Case        | Notes annotation                                                                          |
|-------------|-------------------------------------------------------------------------------------------|
| pan2d       | `pitfall-7-confirmed; ae-budget=3000-confounded-by-resample; fuzz_override=N/A`           |
| nafems_le1  | `pitfall-6-disconfirmed; fuzz_override_pct=8; dim-mismatch=resample-to-1280x800`          |
| cavity2d    | `pitfall-6-secondary; fuzz_override_pct=10; dim-mismatch=resample-to-1280x800; colorbar-region-AE-pct=21.94` |
| heat1d      | `pitfall-6-disconfirmed; fuzz_override_pct=8; dim-mismatch=resample-to-1280x800`          |
| nlsecu      | (no spot-check coverage in Plan 06; nlsecu's TRUNCATED-CAPTURE rationale handled by Plan 02) |

### hexahedron spot-check row (separate from main 5-case grid)

Per the design constraint that hexahedron is NOT in the BUILD-10 baseline
(T-08-23 mitigation), Plan 7 places the hexahedron verdict in a separate
"Spot-check rows" section of CHECKLIST.md (NOT in the canonical 5-case
verdict matrix). Suggested row content:

```
| Case       | Module | Verdict | Note                                                               |
|------------|--------|---------|--------------------------------------------------------------------|
| hexahedron | mail   | GAP     | escalation per evidence/hexahedron-probe.md and 08-VALIDATION.md   |
```

The maintainer's Plan 07 sign-off either:
- accepts the gap with rationale-(b) ("pan2d coverage is sufficient for
  the AA-drift gate"), OR
- refuses and demands gap-closure-(a) (CONTEXT.md amendment naming the
  canonical hexahedron driver).

## Outcome

Spot-check sweep complete: **5/5 reports filed**, **0 PASS / 4 CHECK / 0 FAIL /
1 GAP-DOCUMENTED**. Plan 07 ingests this summary verbatim. The 4 CHECK
verdicts are pending maintainer review at Plan 07 (PASS-on-review expected
on diff-PNG inspection given the resample-dim-mismatch dominant signal).
The 1 GAP-DOCUMENTED verdict (hexahedron) is escalated to the maintainer
in `08-VALIDATION.md ## Open escalations` for binary accept-or-gap-close
decision; Phase 7 sign-off is gated on that resolution.
