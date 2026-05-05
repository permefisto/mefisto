---
phase: 08-ab-validation-on-testa-subset
plan: 01
subsystem: validation
tags: [phase8, bootstrap, harness, ab-validation, autoexit, headless, ImageMagick]

requires:
  - phase: 07-image-gif-and-postscript-export
    provides: Phase 7 PsEmitter + XvueExport landed; pp/*_qt link infrastructure stable; xvue/qt/tests/golden/scene01_driver.f committed; Phase 7 VERIFICATION.md §9 catalogues 3 deferred goldens
provides:
  - bin/phase8_case_batch_map.sh — empirical 5-case MODULE+BATCH+(prereq) map sourced by the harness; replaces fragile *.{mesh,elas,flui,heat,nlse} glob (BLOCKER #1 iter1)
  - bin/ab_compare_pair.sh — tolerance-band gate; ORDER `argc → fuzz → identify → compare` (BLOCKER #5 iter1); fuzz constrained to [1,30] per D-02
  - bin/ab_capture_x11.sh — Xvfb-aware X11 capture wrapper using READY_FILE polling + import -window root
  - bin/ab_sweep_phase8.sh — top-level (case, backend, mode) harness; --smoke-only (BLOCKER #2 iter1) + --baseline with literal `${CASE}` substitution (BLOCKER-B iter2)
  - .planning/phases/08-ab-validation-on-testa-subset/evidence/00-case-batch-map.md (Task 0 evidence with empirical ls + canonical-batch rationale + 5/5 smoke verdicts)
  - .planning/phases/08-ab-validation-on-testa-subset/evidence/00-bootstrap-log.md (Task 1 freshness log + Task 2 deferred-to-human note + Task 3 harness manifest)
  - .planning/phases/08-ab-validation-on-testa-subset/evidence/00-smoke-probes.md (5-case AUTOEXIT smoke probe per-case status)
  - Phase 7 Gap-A closed: pp/*_qt fresh against current libxvueqt.a (mtime ratio OK on all 5)
affects: [08-02-mail-elas-1x-sweep, 08-03-flui-ther-1x-sweep, 08-04-hidpi-2x-sweep, 08-05-omp-sweep, 08-06-checklist-finalize]

tech-stack:
  added:
    - bin/ab_*.sh sweep-harness shell scripts (Phase 8-only, kept under bin/ per EXPORT-06 scope split)
    - bin/phase8_case_batch_map.sh sourceable map (replaces glob assumption)
    - empirical use of MEFISTO_BATCH_X11=1 (existing infrastructure in prpr/pp{mail,elas,flui,ther,nlse}.f) to flip INTERA from 0 to 1 so xvfermer_ fires the capture hook
  patterns:
    - "Empirical-map-as-shell-source pattern: bin/phase8_case_batch_map.sh exports per-case PHASE8_CASE_${CASE}_{MODULE,BATCH,PREREQ_MODULE,PREREQ_BATCH} variables + 4 helper functions (phase8_case_module/batch/prereq_module/prereq_batch). Harness sources via `. \"$MEFISTO/bin/...\"`. Refutes the iter1 glob assumption verbatim."
    - "ORDER-as-load-bearing-comment pattern: ab_compare_pair.sh's `# ORDER: argc → fuzz → identify → compare` header documents an executable contract — the fuzz=50 reject path MUST fire before any file probe, so the verify gate's `bash bin/ab_compare_pair.sh ... 50 | grep 'fuzz must be in'` only succeeds when ORDER is preserved."
    - "Literal-token substitution pattern (BLOCKER-B iter2): --baseline accepts the literal string `${CASE}` (single-quoted by caller) and replaces via bash parameter expansion `${BASELINE_PATH//\\$\\{CASE\\}/$CURRENT_CASE}` — pure string substitution, NOT eval, NOT command substitution."
    - "Smoke-only-flag pattern (BLOCKER #2 iter1): the harness's --smoke-only flag is a Plan-1-bootstrap escape valve that skips the X11-baseline-compare step (no baseline exists yet at Plan 1). The verdict emitted is `verdict=SMOKE`. Plans 2+ omit the flag and the harness enforces baseline presence."

key-files:
  created:
    - bin/phase8_case_batch_map.sh (Task 0 — 77 lines, +x)
    - bin/ab_compare_pair.sh (Task 3 — 92 lines, +x)
    - bin/ab_capture_x11.sh (Task 3 — 111 lines, +x)
    - bin/ab_sweep_phase8.sh (Task 3 — 218 lines, +x)
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/00-case-batch-map.md
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/00-bootstrap-log.md
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/00-smoke-probes.md
  modified:
    - (none — Phase 8 verify-only contract honored; xvue/xvuelc.c byte-identical to 900e297; xvue/qt/src/ untouched in working tree)

key-decisions:
  - "Task 0 BEFORE Task 1 reordering (Rule 3 deviation): the iter1 plan order Task 0→1→2→3 implied probing case-batch behavior on the (initially STALE) pp/*_qt binaries. Empirically discovered that ppnlse_qt at the original mtime would deadlock at startup — the smoke probes only became meaningful after Task 1 refreshed the 4 stale binaries. Effective execution order this run: Task 0-partial → Task 1 → Task 0-completed → Task 3. Task 1 evidence (5 OK lines) was the gate enabling reliable Task 0 evidence."
  - "MEFISTO_BATCH_X11=1 codified into harness contract: empirically required to flip INTERA from 0 to 1 in batch mode so the Qt offscreen window is created and xvfermer_ fires the capture hook. Without it, MEFISTO_QT_CAPTURE_PATH is silently never honored. Existing prpr/pp{mail,elas,flui,ther,nlse}.f infrastructure — Phase 8 Plan 1 just codified its required-not-optional status."
  - "INITIER pre-step codified into harness contract: every case requires `echo \"$CASE\" | $MEFISTO/pp/ppinit` to seed ms10..ms15 storage files before the main module runs. Without it: `ERROR: EXECUTE INITIER BEFORE`. Implicit in the legacy bin/INITIER launcher, now codified in ab_sweep_phase8.sh."
  - "Task 2 (3 deferred goldens) re-deferred to human-bootstrap (Rule 4 deviation): the iter1 plan attempted to autonomize Phase 7 VERIFICATION.md §9's three `checkpoint:human-action` items. The scene01_driver.f link line cannot be solved on this dev host without invasive build-system modification (cbl_tout_qt at line 124 deletes all .o files; .lib archives segfault when used in scene01_x11; curated minimal link line GOTCHA references a non-existent xvue/xvfermer.f). The wave/cavity2d GIF chains require human-issued `99;` saves between MAILLER and the solver step (testa/* batch files have CLOSE; not 99;). All three items are restored to their original Phase 7 §9 human-bootstrap designation."
  - "nlsecu PARTIAL status accepted: NLSER on the 3D cube [-1,1]^3 with 2000 time steps is genuinely too long-running for the harness's 60s budget, on either dispatch path. The MAILLER prereq capture (136593 bytes) is the case's harness-reachable evidence. Plan 2 (NLSER A/B) inherits the long-running constraint with documented mitigation options."

patterns-established:
  - "Empirical-map-as-shell-source: future per-case discovery that exposes hidden prerequisites (e.g., when a new BUILD-10 case is added) extends bin/phase8_case_batch_map.sh in the same MODULE+BATCH+(PREREQ_MODULE+PREREQ_BATCH) shape — no harness changes needed."
  - "Verify-gate ORDER documentation: any future shell script with multiple validation steps documents the order as a load-bearing header comment so re-orderings are caught by the gate exercising the early-reject path."
  - "Literal-token substitution for cross-plan baseline overrides: any future Plan-N flag accepting a `${CASE}`-bearing path uses bash `${VAR//\\$\\{CASE\\}/$CURRENT_CASE}` parameter expansion."

requirements-completed: [VALID-01, VALID-02, VALID-07]

duration: 70min (estimated wall-clock incl. ~7 min cbl_tout_qt rebuild + 5-case smoke sweep)
completed: 2026-05-05
---

# Phase 8 Plan 01: Bootstrap (case-batch map + pp/*_qt freshness + harness scripts + 5-case smoke) Summary

**Phase 8 Wave 2 unblocked: empirical 5-case map + fresh pp/*_qt + 3 sweep-harness scripts shipped; 4/5 testa cases smoke-probed OK on offscreen-Qt; nlsecu PARTIAL (long-running 2000-step NLSE), 3 Phase-7 goldens DEFERRED-TO-HUMAN per original Phase 7 §9 designation.**

## Performance

- **Duration:** ~70 min wall-clock (start 2026-05-05T10:26:59Z; harness shipping completed by ~T11:12:59Z; remaining time was Task 0 nlsecu deep-investigation + worktree-setup recovery)
- **Started:** 2026-05-05T10:26:59Z
- **Completed:** 2026-05-05 (Task 3 commit at 42649d0)
- **Tasks:** 3/4 fully completed (Task 0, Task 1, Task 3); Task 2 deferred-to-human per Rule 4 deviation
- **Files created:** 7
- **Files modified:** 0 (Phase 8 verify-only contract honored)

## Accomplishments

- **Wave 2 unblocking deliverables shipped:** Plans 2-5 can now invoke `bin/ab_sweep_phase8.sh --mode {x11|qt-1x|qt-2x|qt-omp} --cases CSV` against fresh pp/*_qt binaries with a sound case-batch map.
- **5/5 cases empirically mapped (Task 0):** pan2d → mail/pan2d.mesh; nafems_le1 → elas/nafems_le1.elas (prereq mail/nafems_le1.mesh); cavity2d → flui/cavity2d.stoke56cr (prereq mail/cavity2d.meshbf); heat1d → ther/heat1d.heat (prereq mail/heat1d.mesh); nlsecu → nlse/nlsecu.iexrr (prereq mail/nlsecu.meshq2). Refutes iter1 glob (`*.{mesh,elas,flui,heat,nlse}` would silently miss cavity2d's `.meshbf/.stoke56cr` and nlsecu's `.iexrr/.meshq2`).
- **Phase 7 Gap-A closed (Task 1):** `bin/cbl_tout_qt` ran end-to-end exit 0; all 5 pp/*_qt binaries fresh against current libxvueqt.a; ABI count stable at 58.
- **5/5 testa cases smoke-probed (Task 3):** `bin/ab_sweep_phase8.sh --mode qt-1x --smoke-only` ran clean across pan2d/nafems_le1/cavity2d/heat1d/nlsecu. 4/5 captured non-empty PNGs (sizes 70318, 321164, 44723, 21621 bytes); 1/5 PARTIAL (nlsecu — see Deviations).
- **Phase 8 byte-identity invariants verified:** `git diff 900e297..HEAD -- xvue/xvuelc.c` empty (T-07-08 / Phase 8 invariant); `git diff HEAD -- xvue/qt/src/` empty in working tree.
- **Two existing-but-undocumented contracts codified:**
  1. `MEFISTO_BATCH_X11=1` is REQUIRED for headless Qt-side capture (otherwise xvfermer_ never fires). Existing infrastructure in `prpr/*.f`; now baked into the harness.
  2. `INITIER` pre-step (`echo "$CASE" | pp/ppinit`) is REQUIRED for every case to seed ms10..ms15. Implicit in the legacy `bin/INITIER` launcher; now baked into the harness.

## Task Commits

Each task was committed atomically:

1. **Task 0: Empirical case-batch map (5 BUILD-10 cases)** — `90dfd0c` (feat)
2. **Task 1: pp/*_qt freshness + bootstrap-log evidence (D-08)** — `a05ed50` (chore)
3. **Task 3: Sweep harness + smoke probes (3 scripts + 00-smoke-probes.md)** — `42649d0` (feat)

Task 2 was NOT committed — see Deviations below for the Rule 4 escalation.

_Plan metadata commit will be added by the orchestrator after the wave completes._

## Files Created/Modified

### Created (7)

| Path | Lines | Purpose |
|------|-------|---------|
| `bin/phase8_case_batch_map.sh` | 77 | Sourceable case-batch map (Task 0 deliverable). Defines per-case MODULE+BATCH+(PREREQ_MODULE+PREREQ_BATCH) for the 5 BUILD-10 cases plus 4 lookup helpers. |
| `bin/ab_compare_pair.sh` | 92 | Tolerance-band gate wrapping ImageMagick compare. Header `# ORDER: argc → fuzz → identify → compare` documents the load-bearing execution order. -fuzz constrained to [1,30] per D-02. Dimension-mismatch guard via point-resample. |
| `bin/ab_capture_x11.sh` | 111 | Xvfb-aware X11 capture wrapper. Polls MEFISTO_XVFERMER_READY_FILE before invoking `import -window root`; 60s overall timeout; reuses caller-provided $DISPLAY when set, else owns its own xvfb-run session. |
| `bin/ab_sweep_phase8.sh` | 218 | Top-level (case, backend, mode) harness. Sources phase8_case_batch_map.sh; --smoke-only (BLOCKER #2 iter1); --baseline PATH with literal `${CASE}` token substitution via bash parameter expansion (BLOCKER-B iter2). |
| `.planning/phases/08-ab-validation-on-testa-subset/evidence/00-case-batch-map.md` | 354 | Task 0 evidence: verbatim ls -la for 5 cases + canonical-batch selection rationale + per-case end-to-end smoke verdicts (4/5 PASS, 1/5 PARTIAL). |
| `.planning/phases/08-ab-validation-on-testa-subset/evidence/00-bootstrap-log.md` | 192 | Task 1 freshness log (5 OK lines + ABI count) + Task 2 deferred-to-human note (3 sub-sections detailing why) + Task 3 harness manifest + smoke-probe summary. |
| `.planning/phases/08-ab-validation-on-testa-subset/evidence/00-smoke-probes.md` | 95 | Per-case smoke-probe STATUS section (4× OK, 1× PARTIAL) with capture sizes + open-items list for Plan 2. |

### Modified (0)

Phase 8 verify-only contract honored — no source files in the project's
working tree (xvue/, prpr/, util/, mail/, elas/, flui/, ther/, nlse/) were
edited. Specifically:

- `xvue/xvuelc.c` byte-identical to `900e297` (verified by `git diff 900e297..HEAD -- xvue/xvuelc.c` → empty).
- `xvue/qt/src/` no working-tree changes (verified by `git diff HEAD -- xvue/qt/src/` → empty).

## Deviations from Plan

### Rule 3 — Auto-fix blocking issues

**1. [Rule 3 — Blocking issue] Task 0 ↔ Task 1 reorder (effective order: Task 0-partial → Task 1 → Task 0-completed → Task 3)**

- **Found during:** Task 0 — first attempts at the 5-case end-to-end smoke probe.
- **Issue:** The plan as written sequenced Task 0 → Task 1 → Task 2 → Task 3. Empirically discovered that `pp/ppnlse_qt`, `pp/ppelas_qt`, `pp/ppflui_qt`, `pp/ppther_qt` were ALL stale (pre-libxvueqt.a-rebuild May 4 23:50; libxvueqt.a was at May 5 12:25). The stale `pp/ppnlse_qt` deadlocked at startup under `MEFISTO_BATCH_X11=1` — no banner reached stdout even after 240s — invalidating any case-batch-map smoke verdict captured against it.
- **Fix:** Ran Task 1 (`bin/cbl_tout_qt`) before completing Task 0's end-to-end smoke probes. Re-ran nlsecu against fresh `pp/ppnlse_qt`. (NLSER is still long-running on its own merits — see deviation #4 — but the deadlock disappeared.)
- **Effective execution order this run:** Task 0 partial probe revealed staleness → Task 1 rebuilt → Task 0 completed end-to-end probes → Task 3 harness shipped.
- **Files affected:** none — purely an order-of-operations adjustment.
- **Commits:** Task 0 evidence at `90dfd0c` includes the post-rebuild capture sizes; Task 1 at `a05ed50` records the rebuild outcome.

**2. [Rule 3 — Blocking issue] MEFISTO_BATCH_X11=1 baked into harness contract**

- **Found during:** Task 0 — first attempts at pan2d smoke probe.
- **Issue:** The iter1 plan's smoke command did NOT include `MEFISTO_BATCH_X11=1`. Without it, `INTERA = 0` in `prpr/pp{mail,elas,flui,ther,nlse}.f`, the Qt offscreen window is never opened, `xvfermer_` is never called from the main flow, and `MEFISTO_QT_CAPTURE_PATH` is silently never honored. Result: every probe produced a 0-byte capture even when the case ran successfully.
- **Fix:** Added `MEFISTO_BATCH_X11=1` to every dispatch in `bin/ab_sweep_phase8.sh` (qt-1x, qt-2x, qt-omp branches) AND to the prereq-MAILLER branch. Codified in `00-case-batch-map.md` "Smoke-probe environment" section.
- **Files affected:** `bin/ab_sweep_phase8.sh` (built with the env baked in from the start).
- **Commit:** `42649d0`.

**3. [Rule 3 — Blocking issue] INITIER pre-step baked into harness contract**

- **Found during:** Task 0 — every case (not just one) reported `ERROR: EXECUTE INITIER BEFORE`.
- **Issue:** The iter1 plan's smoke command listed workspace prep as `cp -r $MEFISTO/testa/${CASE}/* $MEFISTOX/${CASE}/` followed directly by the main pp/* invocation. But MEFISTO requires `pp/ppinit` to seed `ms10`-`ms15` storage files in the project directory before any other module runs. Implicit in `bin/INITIER` (which is what users normally invoke first), but absent from the iter1 harness contract.
- **Fix:** Added `echo "$CURRENT_CASE" | "$MEFISTO/pp/ppinit"` to every per-case dispatch in `bin/ab_sweep_phase8.sh`. Documented in `00-case-batch-map.md` "Smoke-probe environment".
- **Files affected:** `bin/ab_sweep_phase8.sh`.
- **Commit:** `42649d0`.

### Rule 4 — Architectural decision (escalated, not autonomously fixed)

**4. [Rule 4 — Architectural] Task 2 (3 deferred goldens) re-deferred to human-bootstrap per Phase 7 VERIFICATION.md §9 original designation**

- **Found during:** Task 2 — first attempts at scene01_driver.f link + testa/wave/cavity2d X11+convertepsgif chain.
- **Issue:** The iter1 plan attempted to autonomize three items that Phase 7 explicitly flagged as `checkpoint:human-action`. Each step has a structural blocker that can only be resolved with build-infrastructure modification (out of Phase 8 scope per Phase Boundary Discipline) or human-in-the-loop.
  - **Step A (scene01.eps):** the scene01_driver.f BOOTSTRAP NOTE link line `xvue/*.o util/*.o` cannot be honored because `bin/cbl_tout_qt` line 124 deletes all loose `.o` files at end-of-build. Three workarounds tried — substitute `.lib` archives (segfault at startup), extract `.o` from archives (link fails on undefined `lxouvr_`/`pafcle_`/etc. — defined nowhere in the codebase or in `pp/ppmail` nm output), curated minimal link line per the GOTCHA (references non-existent `xvue/xvfermer.f` — `xvfermer_` is a C function in `xvue/xvuelc.c`, not a Fortran .f file).
  - **Step B (wave_legacy.gif):** `testa/wave` requires INITIER → MAILLER → solver chain, but `wave.mesh` ends with `CLOSE;` not `99;`. The legacy backend's AUTOEXIT path on multi-step chains exits via the documented IEEE_DENORMAL STOP-after-error before invoking `99;` save. So FLUIDER on `wave.wave` reports `ERROR: NO OBJECT`. No `zfxy0*.eps` frames produced; convertepsgif has nothing to combine.
  - **Step C (cavity2d_legacy.gif):** same multi-module legacy chain blocker as Step B.
- **Decision:** RESTORE the original Phase 7 VERIFICATION.md §9 designation — these are `checkpoint:human-action` items by their nature. Recommended action: a future plan or a Phase-7-Plan-06-Task-3 reopen with a human present runs the bootstrap procedure documented in `xvue/qt/tests/golden/scene01_driver.f` BOOTSTRAP NOTE + `bin/convertepsgif` + manual `99;` between MAILLER and solver steps.
- **Files affected:** documentation only — `00-bootstrap-log.md` records the autonomous-attempt failure modes verbatim.
- **Commit:** the Task 3 commit `42649d0` carries the updated bootstrap-log.

### Authentication / human-action gates

(none — Plan 1 had no auth requirements.)

## Worktree-Setup Recovery Note

A non-trivial side-event during execution: the orchestrator's worktree
creation set the per-agent branch `worktree-agent-a53694b17e416fa7e` to
the legacy commit `ac282f8` (pre-Qt-migration), which lacks `xvue/qt/`,
`bin/cbl_tout_qt`, and `.planning/phases/08-*`. All Bash commands
unintentionally ran against the main clone (`/home/mefisto/git/mefisto`)
because the shell tool's CWD reset between invocations bypassed the
worktree. After discovering this mid-Task-0:

1. Surfaced the blocker per #2924 prohibition on `git update-ref` for
   protected refs (HEAD on main is protected; HEAD on `worktree-agent-*`
   is in the allow-list namespace).
2. Used `git fetch . main` + `git reset --hard FETCH_HEAD` on the
   per-agent branch (allow-listed namespace) to fast-forward the worktree
   to current main HEAD `f1d9fe6`. This is recovery from a broken
   worktree-setup, not destruction of meaningful work — the per-agent
   branch had no commits beyond the stale `ac282f8` snapshot.
3. Copied the rebuilt `pp/*_qt` + `xvue/qt/build/libxvueqt.a` from the
   main clone (where Task 1's `bin/cbl_tout_qt` had run) into the
   worktree's gitignored paths to avoid a redundant 7+min rebuild.
4. Re-ran all verify gates against the worktree; redid the smoke probe
   sweep via `bin/ab_sweep_phase8.sh --mode qt-1x --smoke-only` to
   produce evidence sourced from the worktree binaries.
5. All commits landed on `worktree-agent-a53694b17e416fa7e` per the
   pre-commit assertion.

## Threat Model Verification

Per the plan's `<threat_model>`:

- **T-08-01 (stale pp/*_qt):** mitigated — Task 1 freshness gate enforced; all 5 binaries fresh; `00-bootstrap-log.md` "Task 1" carries 5 `OK: pp/pp` lines.
- **T-08-29 (glob misses cases):** mitigated — empirical map shipped; `00-case-batch-map.md` records the verbatim ls output that refutes the glob.
- **T-08-30 (workspace cwd discipline):** mitigated — `bin/ab_sweep_phase8.sh` does `pushd "$PROJDIR"` per case; both ab_capture_x11.sh and ab_sweep_phase8.sh document the cwd contract in their headers.
- **T-08-08 (process death before READY_FILE):** mitigated — `bin/ab_capture_x11.sh` polls both `READY_FILE` existence AND `/proc/$PID`; emits a stderr warning if the process died first.
- **T-08-09 (scene01_driver.f link line):** ACCEPTED at design time, FAILED in execution → restored to Phase 7 §9.1 human-bootstrap (Rule 4 deviation #4 above).
- **T-08-10 (testa case hangs under AUTOEXIT):** mitigated — `timeout 60` in every dispatch; smoke probe per case ran cleanly except nlsecu (long-running, see deviation).
- **T-08-11 (xvue/qt/ source modified during Phase 8):** mitigated — `git diff HEAD -- xvue/qt/src/` and `git diff 900e297..HEAD -- xvue/xvuelc.c` both empty.
- **T-08-12 (Phase 7 VALIDATION-LOG.md DEFERRED rows persist):** NOT YET MITIGATED — the rows remain DEFERRED because the Task 2 goldens were not bootstrapped autonomously (Rule 4 deviation). This is the explicit downstream impact of deferring Task 2 to a human.
- **T-08-31 (ab_compare_pair.sh order changed):** mitigated — `# ORDER: argc → fuzz → identify → compare` header comment shipped; verify gate exercises the fuzz=50 reject path which only succeeds when ORDER is preserved.
- **T-08-32 (--smoke-only flag missing):** mitigated — flag shipped, verify gate `grep -q '\-\-smoke-only|smoke_only' bin/ab_sweep_phase8.sh` passes.
- **T-08-36 (--baseline `${CASE}` substitution):** mitigated — script uses `${BASELINE_PATH//\$\{CASE\}/$CURRENT_CASE}` parameter expansion (literal-token replacement, NOT eval); concrete example demonstrated in the script header documentation.

## Hand-off

**Wave 2 (Plans 2-5) unblocked.** `bin/ab_sweep_phase8.sh` is the canonical
entry point. Plans 2-5 source `bin/phase8_case_batch_map.sh`; no glob
assumptions remain.

Before Plans 2-5 launch, the executor must remember to:

1. Export `MEFISTO=/path/to/mefistosource` and `MEFISTOX=/scratch/working/dir`.
2. Re-run `bin/cbl_tout_qt` if any commit since `f1d9fe6` touched
   `xvue/qt/` (the freshness check `for b in pp/pp*_qt; do [ "$(stat -c
   '%Y' "$b")" -ge "$(stat -c '%Y' xvue/qt/build/libxvueqt.a)" ] && echo
   OK || echo STALE; done` is the contract).
3. Plans 2-5 will need to address two carry-overs from Plan 1:
   - **nlsecu long-running** — pick from the 3 mitigations in
     `00-smoke-probes.md` "## Open Items".
   - **3 deferred goldens** — schedule the Phase 7 §9 human-bootstrap
     before Phase 8 closes; OR explicitly elevate-by-waiver in
     08-CHECKLIST.md.

## Self-Check: PASSED

All claims verified:

- `bin/phase8_case_batch_map.sh`: FOUND, +x, sources cleanly, 5/5 cases defined.
- `bin/ab_compare_pair.sh`: FOUND, +x, ORDER header present, fuzz=50 reject works.
- `bin/ab_capture_x11.sh`: FOUND, +x, syntax-OK.
- `bin/ab_sweep_phase8.sh`: FOUND, +x, --smoke-only AND --baseline flags grep clean.
- `00-case-batch-map.md`: FOUND, 5 `- CASE:` bullets present, 0 unboxed BLOCKED tokens.
- `00-bootstrap-log.md`: FOUND, "Task 1" section has 5 `OK: pp/pp` lines.
- `00-smoke-probes.md`: FOUND, 5+ `## ` headings, 0 `STATUS=BLOCKED` tokens.
- Commits FOUND in git log: `90dfd0c` (Task 0), `a05ed50` (Task 1), `42649d0` (Task 3).
- `xvuelc.c` byte-identical (git diff 900e297..HEAD empty).
- `xvue/qt/src/` clean in working tree (git diff HEAD empty).
