# Wave-2 Post-Merge AE Re-Run

**Date:** 2026-05-05 (post Wave-2 worktree merge)
**Trigger:** Plans 03/04 placeholder diff PNGs replaced with real AE compares now that Plan 02 X11 baselines are merged into main.
**Tool:** `bin/ab_compare_pair.sh` at fuzz=5% (D-02 tolerance band, BLOCKER-5 ORDER guarantee).

## Plan 03 — Qt 1x column (vs X11 baseline)

| Case       | AE     | AE%      | Verdict | Resampled | Diff PNG                                    |
|------------|--------|----------|---------|-----------|---------------------------------------------|
| pan2d      | 540804 | 52.8129% | CHECK   | yes       | evidence/pan2d-qt-1x-diff.png               |
| nafems_le1 | 412827 | 40.3151% | CHECK   | yes       | evidence/nafems_le1-qt-1x-diff.png          |
| cavity2d   | 411003 | 40.1370% | CHECK   | yes       | evidence/cavity2d-qt-1x-diff.png            |
| heat1d    | 143273 | 13.9915% | CHECK   | yes       | evidence/heat1d-qt-1x-diff.png              |
| nlsecu     | 728737 | 71.1657% | CHECK (TRUNCATED-CAPTURE) | yes | evidence/nlsecu-qt-1x-diff.png |

## Plan 04 — Qt HiDPI 2x column (vs X11 baseline, downsampled)

| Case       | AE     | AE%      | Verdict | Resampled | Diff PNG                                    |
|------------|--------|----------|---------|-----------|---------------------------------------------|
| pan2d      | 578116 | 56.4566% | CHECK   | yes       | evidence/pan2d-qt-2x-diff.png               |
| nafems_le1 | 413201 | 40.3517% | CHECK   | yes       | evidence/nafems_le1-qt-2x-diff.png          |
| cavity2d   | 661110 | 64.5615% | CHECK   | yes       | evidence/cavity2d-qt-2x-diff.png            |
| heat1d     | 277787 | 27.1276% | CHECK   | yes       | evidence/heat1d-qt-2x-diff.png              |
| nlsecu     | 794335 | 77.5718% | CHECK (TRUNCATED-CAPTURE) | yes | evidence/nlsecu-qt-2x-diff.png |

## Notes

- All cells CHECK at fuzz=5% — every Qt capture is 760x442 (1x) or 752x156 (HiDPI), while X11 baseline is 1280x800. Resample to common 1280x800 dims dominates the AE diff (kpx-scale visual diff from nearest-neighbor upsample, independent of any backend rendering drift).
- HiDPI dim ratio 0.989 × 0.353 (NOT 2x as Assumption A5 expected — Plan 04 OQ4 finding). Plan 7 must NOT auto-pass HiDPI rows without maintainer acknowledgment of the dim-ratio finding.
- nlsecu rows carry TRUNCATED-CAPTURE annotation: TIME=0.1 truncation per user mitigation; ppnlse_qt offscreen+BATCH_X11 deadlock in canonical run.
- All 10 cells require Plan 7 maintainer review per per-cell sign-off matrix (D-10) — CHECK is NOT FAIL; it means "subjective review required to decide PASS-on-review or escalation".

## Plan 7 ingest

CHECKLIST.md per-cell rows: read AE/AE%/Verdict from this file for `qt-1x` and `qt-2x` columns. The placeholder PNGs from Plans 03/04 have been overwritten with real AE diff output via `bin/ab_compare_pair.sh` post-merge. SHA-256 of the real diff PNGs:

```
$ sha256sum evidence/{cavity2d,heat1d,nafems_le1,nlsecu,pan2d}-qt-{1x,2x}-diff.png
```
(maintainer can compute on-demand from the committed real PNGs)
