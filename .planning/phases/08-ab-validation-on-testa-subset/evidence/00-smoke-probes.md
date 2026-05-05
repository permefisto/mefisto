# Phase 8 Plan 1 Task 3 — Per-Case AUTOEXIT Smoke Probes

**Date:** 2026-05-05
**Author:** Phase 8 Plan 1 executor (autonomous)
**Driver:** `bin/ab_sweep_phase8.sh --mode qt-1x --smoke-only`
**Scope:** 5 BUILD-10 baseline testa cases (pan2d, nafems_le1, cavity2d,
heat1d, nlsecu).

The smoke probes drive each case through the same workspace-prep + INITIER
seed + (optional) MAILLER prereq + main-module dispatch chain that Plans 2-5
will use. With `--smoke-only`, no compare against an X11 baseline is
attempted — the verdict is purely "did the harness produce a non-empty
capture file?".

Probe environment per case (replicated by the harness, sourced from
`bin/phase8_case_batch_map.sh`):

```
QT_QPA_PLATFORM=offscreen
MEFISTO_BATCH_X11=1
MEFISTO_QT_CAPTURE_PATH=/tmp/phase8-smoke/${CASE}-qt-1x.png
MEFISTO_XVSOURIS_AUTOEXIT=1
MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=500
timeout 60 $MEFISTO/pp/pp${MODULE}_qt $BATCH
```

Sweep-log line emitted per case:

```
case=pan2d mode=qt-1x verdict=SMOKE out=/tmp/phase8-smoke/pan2d-qt-1x.png size=70318
case=nafems_le1 mode=qt-1x verdict=SMOKE out=/tmp/phase8-smoke/nafems_le1-qt-1x.png size=321164
case=cavity2d mode=qt-1x verdict=SMOKE out=/tmp/phase8-smoke/cavity2d-qt-1x.png size=44723
case=heat1d mode=qt-1x verdict=SMOKE out=/tmp/phase8-smoke/heat1d-qt-1x.png size=21621
case=nlsecu mode=qt-1x verdict=SMOKE out=/tmp/phase8-smoke/nlsecu-qt-1x.png size=0
```

---

## pan2d

- exit_code=0
- capture_size_bytes=70318
- capture_call_count=1 (xvfermer_ called once at end-of-mesher)
- STATUS=OK

The `IEEE_DENORMAL` STOP-after-error path documented in 00-case-batch-map.md
fires after the canonical post-mesh drawing finishes — the capture lands
before STOP and is the canonical mesh+drawing frame.

## nafems_le1

- exit_code=0
- capture_size_bytes=321164
- capture_call_count=1 (xvfermer_ called once at end of ELASTICER)
- STATUS=OK

MAILLER prereq with nafems_le1.mesh runs first to seed the workspace. The
ELASTICER pass then renders the stress-distribution scene before the
documented IEEE STOP fires.

## cavity2d

- exit_code=0
- capture_size_bytes=44723
- capture_call_count=1 (xvfermer_ called once at end of FLUIDER)
- STATUS=OK

MAILLER prereq with cavity2d.meshbf seeds the workspace. The FLUIDER
Stokes/Navier-Stokes pass renders the velocity field before the IEEE STOP
fires.

## heat1d

- exit_code=0
- capture_size_bytes=21621
- capture_call_count=1 (xvfermer_ called once at end of THERMICER)
- STATUS=OK

MAILLER prereq with heat1d.mesh seeds the workspace. The THERMICER
unsteady-1D thermal pass renders the temperature distribution at the final
time step.

## nlsecu

- exit_code=124 (timeout) at the harness's 60s budget
- capture_size_bytes=0
- capture_call_count=0 (no xvfermer_ reached within 60s)
- STATUS=PARTIAL (deferred per 00-case-batch-map.md analysis)

The nlsecu case requests Final TIME=20, Step=0.01 → 2000 time steps. The
timeline on this hardware exceeds 60s by an order of magnitude regardless of
dispatch path. The MAILLER prereq capture (size 136593 in Task 0 evidence)
proves the workspace+cwd discipline is sound; the NLSER main-step capture
remains unreachable within the Plan 1 60s harness budget.

---

## Open Items

- **nlsecu**: Plan 2 (NLSER A/B sweep) inherits the long-running constraint.
  Mitigation options for Plan 2 to choose from:
    1. Use a longer per-case timeout (300-600s) for nlsecu specifically.
    2. Substitute a smaller NLSE testa case from BUILD-10 scope (none
       currently exists; would require a BUILD-10 amendment).
    3. Accept the MAILLER prereq capture as the case's headless evidence
       and run the full NLSER comparison out-of-band per CLAUDE.md
       "long-running tests... ask the user to run".

- **Cross-reference with 00-case-batch-map.md Task 0 verdict:** consistent
  (both surfaces report nlsecu as long-running PARTIAL; 4/5 PASS for the
  others).
