# Phase 8 Plan 1 — Bootstrap Log

**Plan:** 08-01 — Bootstrap: pp/*_qt freshness + 3 deferred goldens + harness scripts + 5-case smoke
**Start:** 2026-05-05
**Author:** Phase 8 Plan 1 executor (autonomous)

This document records per-task bootstrap evidence for the four Phase 8 Plan 1
tasks (Task 0 — case-batch map; Task 1 — pp/*_qt freshness; Task 2 — 3 deferred
goldens + ctest re-verify; Task 3 — sweep-harness scripts + smoke probes).

---

## Task 1 — pp/*_qt freshness (D-08)

**Goal:** Refresh `pp/*_qt` against the current `libxvueqt.a` so the A/B
sweep cannot get a false-PASS from stale Qt binaries (Phase 7 Gap-A).

### bin/cbl_tout_qt result

```
cbl_tout_qt_exit=0
```

Build duration: ~7 min (CMake configure + libxvueqt.a build + 4 test-binary
targets + 7 Fortran library rebuilds + 5 cb*_qt + 5 cbxvtest*_qt linker
invocations). Last log line:

```
TOUS les MODULES EXECUTABLES sans debogueur de /home/mefisto/git/mefisto sous LINUX sont crees — Qt variant
```

### Freshness check (pp/*_qt mtime >= libxvueqt.a mtime)

```
OK: pp/ppmail_qt (mtime 1777977638 >= lib 1777977218)
OK: pp/ppelas_qt (mtime 1777977639 >= lib 1777977218)
OK: pp/ppflui_qt (mtime 1777977641 >= lib 1777977218)
OK: pp/ppther_qt (mtime 1777977640 >= lib 1777977218)
OK: pp/ppnlse_qt (mtime 1777977642 >= lib 1777977218)
```

(Mtimes recorded post-build in the main clone; identical artifacts staged in
this per-agent worktree by file copy after the build completed — both clones
share the same `.git`. After staging into the worktree the freshness ratio is
preserved: all 5 pp/*_qt mtime equal the libxvueqt.a mtime.)

### ABI count

```
verify_abi: nm count: 58  header count: 58
```

(Stable Phase 7 baseline; no ABI drift introduced by Phase 8 Plan 1.)

### ISO date

`2026-05-05`

---

## Task 2 — Phase-7 deferred goldens (D-06, D-07) — DEFERRED-TO-HUMAN

**Status:** NOT COMPLETED autonomously — defers to Phase 7 VERIFICATION.md §9
human-bootstrap procedure (the original designation pre-iter1 revision).

### Step A — scene01.eps (PostScript byte parity)

**Attempted:** compile `xvue/qt/tests/golden/scene01_driver.f` against the X11
backend per the BOOTSTRAP NOTE in the file header.

**Sub-step 1 (compile driver):** PASS.
```
gfortran -I$MEFISTO/incl -c scene01_driver.f
```
produces `scene01_driver.o` cleanly.

**Sub-step 2 (link to executable):** FAIL on this dev host.

The BOOTSTRAP NOTE link line (`xvue/*.o util/*.o`) cannot be honored because
`bin/cbl_tout_qt` in its current form (line 124, `rm $MEFISTO/*/*.o`) deletes
all loose `.o` files at the end of every Qt build. Three workarounds were
attempted, all failed:

1. **Substitute `.lib` archives** for `xvue/lib util/lib`: link still fails
   with segfault at startup (no symbol-resolution error, but the linked
   `scene01_x11` binary segfaults before reaching its first instruction —
   suggests an unresolved data-section initialization or missing OMP runtime
   reference that the archive form does not surface as a link-time error).

2. **Extract `.o` files from archives** (`ar x xvue/lib && ar x util/lib`)
   and link verbatim per the BOOTSTRAP NOTE: link fails on undefined
   references — `lxouvr_`, `pafcle_`, `pafcec_`, `coorle_`, `remass2_`,
   `x_`. None of these symbols are defined under `util/`, `xvue/`, or any
   `prpr/*.f` main program; they are not present in the live `pp/ppmail`
   binary either (verified via `nm`). Their definition site is unclear and
   resolution is non-trivial without invasive build-system investigation.

3. **Curated minimal link line** per the GOTCHA in 08-RESEARCH.md §Pattern 3
   (`xvue/xvuelc.o xvue/xvinit.o xvue/xvouvrir.o xvue/xvfermer.o`): fails
   immediately because `xvue/xvfermer.f` does NOT exist as a Fortran file
   (`xvfermer_` is a C function defined in `xvue/xvuelc.c`; the GOTCHA's
   curated set was speculative). And `xvue/xvinit.f` itself depends on the
   broader Fortran library so the minimal set keeps growing toward the
   already-failing full link.

**Conclusion:** the link-line for `scene01_driver.f` cannot be solved by an
autonomous executor on this dev host without invasive build-system work
(adding a permanent CMake/make target that owns the `scene01_x11` link, or
preserving `.o` files past `cbl_tout_qt` cleanup). Both options require
modifying the Fortran build infrastructure — out of Phase 8's verify-only
scope per the explicit Phase Boundary Discipline (08-RESEARCH.md).

This step IS the original Phase 7 VERIFICATION.md §9.1
`checkpoint:human-action` — which explicitly states
*"Why human: Requires X11 display (Xvfb or real) to materialise TEMPORAIRE.EPS
from the legacy backend. Autonomous executor has no attached terminal."*
The Phase 8 Plan 1 iter1 revision attempted to autonomize this checkpoint,
but the underlying build-infrastructure block remains. Recommended action:
re-establish the Phase 7 §9.1 human-bootstrap as the resolution path.

### Step B — wave_legacy.gif

**Attempted:** run `testa/wave` through the X11 backend + `bin/convertepsgif`.

`INITIER` + `MAILLER wave.mesh` execute under `xvfb-run --auto-servernum`
with `MEFISTO_XVSOURIS_AUTOEXIT=1`, but the `MAILLER` step does not produce
a saved object (the AUTOEXIT path on the legacy backend exits via the same
documented `IEEE_DENORMAL` STOP-after-error path before invoking `99;` save).
The subsequent `FLUIDER wave.wave` invocation reports `ERROR: NO OBJECT`.
Consequently, no `zfxy0*.eps` frames are produced and `bin/convertepsgif`
has no inputs to combine.

The legacy AUTOEXIT chain on multi-step driver test cases (testa/wave needs
INITIER → MAILLER → wave-solver, with the mesher saving via `99;` before
the next step takes over) has the same root-cause as the iter1
case-batch-map findings: testa/* `.mesh` files do not include `99;` save
commands, so the chain breaks unless the user issues them interactively.

**Conclusion:** Phase 7 VERIFICATION.md §9.2 designates this as a
`checkpoint:human-action` for the same root reason ("Visual GIF
comparison... requires X11+Qt side-by-side session"). The autonomous path
cannot reliably drive the multi-module legacy chain without the user
issuing `99;` between each step (or modifying the testa/* batch files,
which is out of scope per BUILD-10 immutability of the 5 cases). Recommended
action: human-bootstrap per Phase 7 VERIFICATION.md §9.2.

### Step C — cavity2d_legacy.gif

Same blocker as Step B (multi-module chain on the legacy backend cannot be
driven autonomously without `99;` saves between steps). Designation:
human-bootstrap per Phase 7 VERIFICATION.md §9.3.

### ctest re-verify (D-07)

Not executed (the 3 goldens are not in place; QSKIPs would still be present
and 0-SKIP gate would not flip). The Phase 7 ctest result published in
07/VALIDATION-LOG.md remains the latest authoritative state.

### VALIDATION-LOG.md edit

Not executed (the 3 DEFERRED rows are not yet eligible for elevation to
PASS — the goldens that would flip them are still un-bootstrapped).

### ISO date

`2026-05-05`

---

## Task 3 — Sweep harness + smoke probes

### Shipped scripts (4, including Task 0 deliverable)

| Script | Lines | +x | Purpose |
|--------|-------|----|---------|
| `bin/phase8_case_batch_map.sh` | 75 | yes | Sourceable case-batch map (Task 0 deliverable) |
| `bin/ab_compare_pair.sh` | 80 | yes | tolerance-band gate; ORDER: argc → fuzz → identify → compare |
| `bin/ab_capture_x11.sh` | 100 | yes | Xvfb-aware X11 capture wrapper using READY_FILE polling + import -window root |
| `bin/ab_sweep_phase8.sh` | 200 | yes | top-level (case, backend, mode) harness; supports --smoke-only and --baseline |

### Verify-gate semantic checks

```
$ bash bin/ab_compare_pair.sh /tmp/_pair_a.png /tmp/_pair_b.png /tmp/_diff.png 50
ab_compare_pair: fuzz must be in [1,30] (got: 50)
$ grep -q '^# ORDER: argc → fuzz → identify → compare' bin/ab_compare_pair.sh && echo OK
OK
$ grep -q '\-\-smoke-only' bin/ab_sweep_phase8.sh && echo OK
OK
$ grep -q '\-\-baseline' bin/ab_sweep_phase8.sh && echo OK
OK
```

### Smoke probe sweep result

Full harness invocation:

```
bin/ab_sweep_phase8.sh --mode qt-1x --smoke-only \
    --cases pan2d,nafems_le1,cavity2d,heat1d,nlsecu \
    --out-dir /tmp/phase8-smoke
```

Per-cell verdict (5 lines):

```
case=pan2d mode=qt-1x verdict=SMOKE out=/tmp/phase8-smoke/pan2d-qt-1x.png size=70318
case=nafems_le1 mode=qt-1x verdict=SMOKE out=/tmp/phase8-smoke/nafems_le1-qt-1x.png size=321164
case=cavity2d mode=qt-1x verdict=SMOKE out=/tmp/phase8-smoke/cavity2d-qt-1x.png size=44723
case=heat1d mode=qt-1x verdict=SMOKE out=/tmp/phase8-smoke/heat1d-qt-1x.png size=21621
case=nlsecu mode=qt-1x verdict=SMOKE out=/tmp/phase8-smoke/nlsecu-qt-1x.png size=0
```

4/5 OK; 1/5 PARTIAL (nlsecu — see 00-smoke-probes.md). Detail and per-case
status in `.planning/phases/08-ab-validation-on-testa-subset/evidence/
00-smoke-probes.md`.

### ISO date

`2026-05-05`
