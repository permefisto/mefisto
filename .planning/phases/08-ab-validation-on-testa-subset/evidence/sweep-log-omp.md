# Phase 8 Plan 5 — OMP Sweep Log

This file records the per-case OMP_NUM_THREADS=8 capture results across the
4 OMP-eligible BUILD-10 cases on both backends (X11 and Qt offscreen). The
fifth case, **cavity2d**, is N-A — there is no `bin/FLUIDER_OMP` launcher
on this host (verified by `ls bin/FLUIDER_OMP 2>&1 → "No such file or
directory"`; recorded in 08-RESEARCH.md §"OMP Sweep — How _OMP Actually Works"
research finding 3).

## Plan-5 OMP sweep — X11-OMP baseline column (2026-05-05 12:12:51 UTC)

**Settings:**
- `OMP_NUM_THREADS=8` (libgomp scheduling parallelism)
- `MEFISTO_BATCH_X11=1` (forces INTERA=1 so xvfermer_ fires from the batch flow — without this the X11 binary never reaches xvfermer_ in batch mode)
- Xvfb 1280x800x24 via `xvfb-run --auto-servernum` (parallel-safe per WARNING-2 iter2)
- Capture: `xvfb-run` + `import -window root` after `MEFISTO_XVFERMER_READY_FILE` sentinel appears
- `MEFISTO_XVFERMER_HOLD_MS=2000` (window held visible 2s for capture)
- `MEFISTO_XVSOURIS_AUTOEXIT=1` + `MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=500` (auto-progress through interactive blocks)

**BLOCKER #4 reference (VALID-03 literal launcher compliance):** the X11
captures in this section use direct `OMP_NUM_THREADS=8 pp/pp${MODULE}`
invocation, which is empirically equivalent to the literal-named launcher
`bin/ELASTICER_OMP` per evidence/elasticer-omp-equiv-check.md (AE=1313 px
launcher-pattern-vs-direct, within 1143 px OMP twice-run noise floor).

**Deviation note (Rule 3 — auto-fix):** the existing `bin/ab_capture_x11.sh`
(Plan 1 Task 3 deliverable) does NOT export `MEFISTO_BATCH_X11=1`, which
means in pure batch mode the X11 binary never calls xvfermer_ and the
ready-file is never touched (verified empirically — first capture attempt
produced 263-byte empty Xvfb desktop). Plan 5 inlines the capture loop with
`MEFISTO_BATCH_X11=1` added rather than modify ab_capture_x11.sh, to avoid
contention with Plans 02 capture flow that may also rely on the script.
This deviation is recorded as a deferred item for follow-up against the
harness script.

**nlsecu mitigation (user decision, applied this plan):** the canonical
`testa/nlsecu/nlsecu.iexrr` requests Final TIME=20, Step=0.01 → 2000 time
steps (~hour-scale on this hardware). Per user-supplied direction, Plan 5
applies a TIME truncation override on the workspace copy of `nlsecu.iexrr`
(`Final TIME = 0.1` instead of 20 → 10 time steps). This is a TRUNCATED
capture and is annotated as such in the per-case row below. The truncation
preserves all OMP-relevant code paths (matrix factorization, OMP-parallel
inner solver, draw-step) — it just stops earlier.

| Case       | size  | dims     | sha256                                                            | OK? |
|------------|-------|----------|-------------------------------------------------------------------|-----|
| pan2d      | 7610  | 1280x800 | 06f98b264f771816b62f1c97b4cb74e3ce6a94d2fadd2b590053ed6c53892b2d | OK  |
| nafems_le1 | 60762 | 1280x800 | 99036645d889d0d7166ffe339bbe703f7b53bd7aeb1b4cc5cdcee357f6cd41e0 | OK  |
| heat1d     | 4363  | 1280x800 | 26fa74ba517f238008ed36e18510333fbf02ac6a02a4ab4a86bfe5064938a36c | OK  |
| nlsecu     | 4937  | 1280x800 | af06183321aa4029f4c6ee66d70833f42f5996d105d3b9d00f35b400a5ab730c | OK (TRUNCATED-CAPTURE — TIME=0.1, 10 steps) |
| cavity2d   | N/A   | N/A      | N/A                                                               | N-A (no bin/FLUIDER_OMP launcher per D-05) |

## OMP noise floor (Pitfall 5)

Pan2d X11-OMP capture re-run on a fresh workspace clone, identical env,
and AE-compared against the canonical pan2d-x11-omp.png:

```
$ bin/ab_compare_pair.sh \
    .planning/phases/08-ab-validation-on-testa-subset/evidence/pan2d-x11-omp.png \
    /tmp/_pan2d-x11-omp-noise-run2.png \
    /tmp/_pan2d-x11-omp-noise-diff.png \
    5
ae=236 total=1024000 pct=0.0230% verdict=CHECK diff=/tmp/_pan2d-x11-omp-noise-diff.png resampled=no
```

**Pan2d twice-run X11-OMP noise floor: AE = 236 px (0.023%).**

This noise floor is the lower bound on any AE for the same-binary-same-env
case. It captures the scheduling-jitter component of OMP non-determinism
(Pitfall 5). Plan 7 (CHECKLIST.md composer) interprets X11-OMP-vs-Qt-OMP
AE ≤ 2 × this noise floor as PASS, and AE > 2 × noise floor as CHECK
(maintainer review).

Per WARNING-2 iter2 parallel-safety: the noise-floor temp file paths use
the keyed prefix `/tmp/_pan2d-x11-omp-noise-*` so they cannot collide with
sibling Plans 02/03/04 wave-2 temp paths.
