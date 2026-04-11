---
phase: 02-drawing-primitives-backing-pixmap
plan: 04
subsystem: xvue-qt-phase2-closure
tags: [qt6, validation, manual-checkpoint, phase-closure, wave-3]

requires:
  - phase: 02-drawing-primitives-backing-pixmap
    plan: 03
    provides: "All 13 DRAW-XX entry points shipped with real bodies; zero DRAW warn-once lines; pp/ppxvtest0_qt exits 0"
provides:
  - "Clean-rebuild regression proof: rm -rf xvue/qt/build && bin/cbl_tout_qt green on first try"
  - "Legacy X11 A/B invariant (BUILD-07) preserved: bin/cbl_tout still green"
  - "Human visual checkpoint approved — DRAW-02..06 geometry visible, DRAW-09 resize-preserve confirmed on a live Qt 6 display, DRAW-07 effacer verified, SHELL-02 reopen cycle clean"
  - "02-VALIDATION.md closed: status complete, nyquist_compliant true, wave_0_complete true, 13-row Per-Task Verification Map, dated Approval line"
  - "Extended prpr/xvtest0.f driver holds (SLEEP 15/10/3) — human-workable visual checkpoint without growing DRAW coverage"
affects:
  - "Phase 2 closed. Phase 3 (text, fonts, colormap) can be planned against a stable backing-pixmap + single-long-lived-painter substrate."

tech-stack:
  added: []
  patterns:
    - "Test-driver hold tuning as deviation Rule 3 (blocking issue): human checkpoint was unworkable with SLEEP(1); bumped to SLEEP(15/10/3) in a separate atomic commit — test coverage unchanged, only usability"
    - "Manual visual checkpoint wedged between two automated gates (02-04-01 clean rebuild, 02-04-03 docs close-out) to preserve Nyquist sampling continuity"

key-files:
  created:
    - .planning/phases/02-drawing-primitives-backing-pixmap/02-04-SUMMARY.md
  modified:
    - prpr/xvtest0.f
    - .planning/phases/02-drawing-primitives-backing-pixmap/02-VALIDATION.md

key-decisions:
  - "xvtest0.f driver hold extension is committed as a fix() not a test() because it is a usability correction to an existing test driver — the DRAW-01..09 coverage remains identical"
  - "Check 4 (effacer clears mid-sequence) accepted as implicitly verified by the earlier Wave 1 shorter-hold run rather than re-run with extended holds; the post-effacer-only state was already directly observed"
  - "STATE.md / ROADMAP.md intentionally NOT updated in this plan — the wave orchestrator owns those writes after the full Wave 3 completes"

requirements-completed:
  - DRAW-01
  - DRAW-02
  - DRAW-03
  - DRAW-04
  - DRAW-05
  - DRAW-06
  - DRAW-07
  - DRAW-08
  - DRAW-09

duration: ~40min (incl. human checkpoint wait)
completed: 2026-04-11
---

# Phase 2 Plan 04: Wave 3 Validation Checkpoint Summary

**Phase 2 closed. Clean rebuild green, legacy X11 A/B still green, every DRAW-01..09 requirement verified (headless + human visual), 02-VALIDATION.md sealed with `nyquist_compliant: true` and a dated 2026-04-11 approval line. One tactical test-driver usability fix (xvtest0.f SLEEP hold extensions) landed as a separate atomic commit so `pp/ppxvtest0_qt` became workable for interactive visual inspection without altering DRAW coverage.**

## Performance

- **Duration:** ~40 min wall-clock (including human checkpoint wait)
- **Tasks:** 3 / 3 (+ 1 bonus driver-usability fix)
- **Files modified:** 2 (1 driver, 1 validation doc) + 1 summary created

## Task Commits

| # | Task / Item                                                                        | Commit   | Type  |
|---|------------------------------------------------------------------------------------|----------|-------|
| 0 | Driver hold extension for human checkpoint (xvtest0.f SLEEP 15/10/3)               | `a82247b`| fix   |
| 1 | Clean rebuild + headless regression sweep (Qt + legacy X11 + pp/ppxvtest0_qt)       | (no code change — rebuild + smoke only; results captured in this SUMMARY) | — |
| 2 | Human visual checkpoint — 7 checks                                                  | (no commit — human gate) | — |
| 3 | Close 02-VALIDATION.md — Per-Task Verification Map + Sign-Off                       | `c84741b`| docs  |

The final metadata commit for this plan (containing the SUMMARY below) follows after the self-check. Plan 02-04 thus contributes **3 commits** to the repo: `a82247b` (driver fix), `c84741b` (validation close), and the SUMMARY commit.

## Task 0 — Driver hold extension (uncommitted fix discovered at checkpoint)

The plan's Task 2 `<how-to-verify>` block assumed the existing `CALL SLEEP(1)` in `prpr/xvtest0.f` would give the human operator enough time to:

1. visually audit all six DRAW-02..06 primitives on the pre-effacer scene,
2. perform a drag-resize to verify DRAW-09 top-left preserve, and
3. observe the post-effacer state.

One second is not enough wall-clock for a human to do any of this, let alone all three back-to-back. Under deviation Rule 3 (fix blocking issues that prevent completing the current task), the driver was amended to a three-hold structure:

1. **Pre-effacer scene:** `CALL XVVOIR` + `CALL SLEEP(15)` inserted AFTER `xvbordarcellipse` and BEFORE `CALL EFFACER`. The operator can inspect DRAW-02..06 geometry (Check 1), pen styles (Check 2), antialiasing (Check 3), AND drag-resize for DRAW-09 top-left preserve (Check 5) all within this 15-second window.
2. **Post-effacer dashed lines:** `SLEEP(1)` → `SLEEP(10)` on the existing post-effacer hold so the operator can confirm DRAW-07 effacer fully cleared the pre-effacer scene (Check 4).
3. **Reopen cycle:** `SLEEP(1)` → `SLEEP(3)` on the reopen hold so the blank reopened canvas is observable (Check 6 — SHELL-02 regression guard).

DRAW-01..09 coverage is unchanged — only SLEEP durations and one added `XVVOIR + SLEEP` pause between the draws and EFFACER. Committed as `a82247b` (`fix(02-04): extend xvtest0.f holds for manual visual checkpoint`) before proceeding to Task 3 so the driver improvement has its own audit trail.

## Task 1 — Clean rebuild + headless regression sweep

Executed before the handoff to the checkpoint agent:

```
rm -rf xvue/qt/build
bin/cbl_tout_qt          → exit 0   (libxvueqt.a + 5 pp*_qt executables rebuilt from scratch)
bin/cbl_tout             → exit 0   (legacy X11 path — BUILD-07 preserved)
pp/ppxvtest0_qt          → exit 0   (13/11/0 warn-once progression confirmed at terminal state: 0)
```

Grep proofs (all returning the expected results):

- `grep -E 'xvue-qt: (xvtrait|xvftrait|xvtraits|xvface|xvfacetraits|xvrectangle|xvbordrectangle|xvfrectangle|xvfbordrectangle|xvarcellipse|xvbordarcellipse|xvtypetrait|xvepaisseur|effacer|xvvoir|xvpxfenetre|xvfond)_' /tmp/02-04-xvtest0.out` → **empty** (zero DRAW warn-once lines)
- `grep -E 'lasopsc|courgb|counb|ypixels|iFa|iRe|ire|iEl|iel' xvue/qt/src/xvue_qt_api.cpp` → expected D-26/D-28 clean result (the sole `ire` substring hit on the word "fires" in a Phase 1 comment is a documented known baseline, not a real PS-recording state port)
- `verify_abi` — 57 symbols, header count = nm count
- `verify_no_exec` — OK

## Task 2 — Human visual checkpoint (resume-signal: approved)

The human operator ran `pp/ppxvtest0_qt` against the rebuilt binary (pp/ppxvtest0_qt, 136512 bytes, rebuilt 20:29 local) on a live Qt 6 display session. Results:

| Check | Requirement | Method | Result |
|-------|-------------|--------|--------|
| 1. Geometry present | DRAW-02..05 | First SLEEP(15) hold — visual audit | ✅ PASS — horizontal line, polyline, filled polygon, 4 rectangles, filled pie slice, outlined ellipse arc all visible on black background |
| 2. Pen style variety | DRAW-06 | First SLEEP(15) hold — visual audit | ✅ PASS — dashed patterns distinguishable, thicker stroke after XVEPAISSEUR(2) visible |
| 3. Antialiasing | DRAW-08 | First SLEEP(15) hold — visual audit | ✅ PASS — polyline and arc edges smooth, no stair-stepping (Pitfall 5 not triggered) |
| 4. Effacer mid-sequence clear | DRAW-07 | Earlier Wave 1 shorter-hold run + current post-effacer SLEEP(10) | ✅ PASS (implicit) — the earlier short-hold run showed ONLY the post-effacer dashed lines during the SLEEP hold, confirming the effacer fully cleared the pre-effacer scene |
| 5. Resize preserves top-left | DRAW-09 | First SLEEP(15) hold — drag-resize larger | ✅ PASS — drawn geometry stayed anchored at top-left, new area filled black, no stretch/recenter/disappearance |
| 6. Reopen cycle clean | SHELL-02 | Headless via Task 1 | ✅ PASS — no `QApplication: there can only be one` assertion, exit 0 |
| 7. Process exit clean | — | Headless via Task 1 | ✅ PASS — exit 0, no SIGSEGV, no `QPainter not active` warnings |

The resume signal from the operator was "approved" for the live-display checks (1, 2, 3, 5) with the note that Check 4 was already covered by the earlier Wave 1 observation and Checks 6/7 were already covered by the automated Task 1 sweep. This resolution is recorded in the 02-VALIDATION.md Approval line.

## Task 3 — Close 02-VALIDATION.md

Frontmatter flipped:

```yaml
status: draft          → status: complete
nyquist_compliant: false → nyquist_compliant: true
wave_0_complete: false → wave_0_complete: true
completed: 2026-04-11  (new)
```

The Per-Task Verification Map was populated with **13 rows** spanning plans 02-01 through 02-04 (4 + 3 + 3 + 3 = 13 tasks including the 02-04-03 docs close-out self-reference). Every row has a `DRAW-XX` requirement, a real automated command copied from the corresponding plan's `<verify><automated>` block (or "Human visual checkpoint" for the single manual gate), and `Status = ✅`.

Wave 0 requirements checklist fully ticked with commit references:

- xvtest0.f draw-coverage → `b556037` (02-01 task 4) + `a82247b` (02-04 hold extension)
- MEFISTO_POINT_AUDIT.md → `e8accd6` (02-01 task 1)
- README_RESIZE.md → `4ee34a6` (02-01 task 2)

Manual-Only Verifications table expanded from 2 to 3 rows (geometry present, resize preserve, effacer mid-sequence clear), each with an "Approved 2026-04-11" status. Automated fallback noted for Checks 6/7 (covered by Task 1 clean-rebuild sweep, no live display required).

Validation Sign-Off: all 6 checkboxes ticked. Approval line reads:

> **Approval:** approved 2026-04-11 by human visual checkpoint (plan 02-04, task 2) — all 7 visual checks passed. Checks 1/2/3/5 confirmed directly against `pp/ppxvtest0_qt` on a live display session (driver commit `a82247b`, binary rebuilt 20:29 local). Check 4 implicitly verified by the earlier shorter-hold run. Checks 6/7 verified headlessly by the plan 02-04 task 1 clean-rebuild regression sweep.

Committed as `c84741b` (`docs(02-04): close 02-VALIDATION.md — all Phase 2 DRAW reqs verified`).

## Phase 2 Requirements Closure

| Req ID | Name                              | Closed by         |
|--------|-----------------------------------|-------------------|
| DRAW-01 | Single long-lived QPainter        | Plan 02-02        |
| DRAW-02 | Lines / polylines                 | Plan 02-03 task 1 |
| DRAW-03 | Polygons (fill + fill+outline)    | Plan 02-03 task 1 |
| DRAW-04 | Rectangles (4 variants)           | Plan 02-03 task 2 |
| DRAW-05 | Ellipse arcs (drawPie + drawArc)  | Plan 02-03 task 3 |
| DRAW-06 | Pen style + width (3 cases)       | Plan 02-03 task 2 |
| DRAW-07 | effacer (clear + flush)           | Plan 02-02 task 3 |
| DRAW-08 | Antialiasing hint on painter      | Plan 02-02 task 2 |
| DRAW-09 | Resize preserves top-left         | Plan 02-02 task 2 |

All 9 DRAW-XX requirements are green, with both headless automated verification and the Wave 3 human visual checkpoint approving the live-display behavior.

## Build Verification

| Gate                                        | Result |
|---------------------------------------------|--------|
| `rm -rf xvue/qt/build && bin/cbl_tout_qt`   | exit 0 |
| `bin/cbl_tout` (legacy X11 regression)      | exit 0 |
| `pp/ppxvtest0_qt` headless                  | exit 0, zero DRAW warn-once lines |
| `pp/ppxvtest0_qt` on live display           | approved by human operator |
| `verify_abi`                                | 57 symbols (nm count = header count) |
| `verify_no_exec`                            | OK |

## Decisions Made

- **Driver usability fix committed as a standalone `fix()` commit** (`a82247b`) rather than squashed into Task 3's docs commit. Rationale: a future `git bisect` on DRAW coverage should find exactly one commit per concern; the hold extension is a test-driver usability correction, not a docs change.
- **Check 4 accepted as implicitly verified** rather than re-run with extended holds. The earlier Wave 1 shorter-hold run had already shown only post-effacer content during its single SLEEP, so re-running with longer holds for the same check would only add wall-clock without new evidence.
- **STATE.md / ROADMAP.md intentionally NOT touched in this plan.** The prompt directive from the orchestrator was explicit: the wave orchestrator owns state/roadmap writes after the full wave completes. This plan writes only the SUMMARY + VALIDATION docs.

## Deviations from Plan

- **[Rule 3 — Blocking issue] Extended xvtest0.f SLEEP holds for human checkpoint.** The plan's Task 2 assumed `SLEEP(1)` would suffice for the human visual checks; it did not. Extended to `SLEEP(15)` (pre-effacer) + `SLEEP(10)` (post-effacer) + `SLEEP(3)` (reopen) and added a `XVVOIR + SLEEP(15)` flush point between the draws and EFFACER. Committed as `a82247b`. DRAW coverage unchanged, only usability.
- **[Rule 3 — Blocking issue] Check 4 resolution path.** The plan implied Check 4 would be verified during the same run as Checks 1/2/3/5, but the extended driver re-runs that once DRAW-07 effacer runs mid-sequence, the pre-effacer scene is already gone — the operator cannot directly compare "before effacer" to "after effacer" within a single live-display run unless they quickly glance at the SLEEP(15) window first and the SLEEP(10) window next. Accepted the earlier Wave 1 shorter-hold run as the direct evidence for Check 4 and documented this explicitly in the Validation Sign-Off Approval line.

## Issues Encountered

None beyond the two deviations above.

## Authentication Gates

None — pure local C++/CMake/Fortran work, no network, no credentials, no CLI auth.

## Next Phase Readiness

- **Phase 2 closed.** All 9 DRAW-XX requirements green, 02-VALIDATION.md sealed with dated approval, `pp/ppxvtest0_qt` is the progressive smoke baseline (exit 0, zero DRAW-XX warn-once lines) that Phase 3 will inherit.
- **Substrate for Phase 3 (text, fonts, colormap):** the backing-pixmap + single-long-lived-painter invariants from Wave 1 are stable and visually verified. Phase 3 entry points (text rendering, colormap `xvcoul*`, `xvfond_` extensions for non-white-only palette) can assume `state_->painter_->isActive() == true` whenever a window is open, exactly as Wave 2 primitives did.
- **Legacy X11 A/B invariant (BUILD-07) still green.** `bin/cbl_tout` builds all 5 legacy executables unchanged. Phase 3 will continue to respect the A/B comparison window until the X11 backend is retired.
- **Test driver ready for extension.** `prpr/xvtest0.f` extended hold structure (SLEEP 15/10/3) is workable for future manual checkpoints. Phase 3 text/colormap tests can extend this driver with a new coverage section before the existing draw section without breaking the current DRAW-XX gate.
- **No blockers.** `libxvueqt.a` still reports 57 ABI symbols; Phase 3 planning can begin.

**Phase 2 complete — ready for Phase 3.**

## Self-Check: PASSED

Verified commits exist in `git log --oneline -5`:
- FOUND: `a82247b` — Task 0 (driver hold extension)
- FOUND: `c84741b` — Task 3 (02-VALIDATION.md closure)

Verified artifacts exist:
- FOUND: `.planning/phases/02-drawing-primitives-backing-pixmap/02-VALIDATION.md` (status: complete, nyquist_compliant: true, wave_0_complete: true)
- FOUND: `.planning/phases/02-drawing-primitives-backing-pixmap/02-04-SUMMARY.md` (this file)
- FOUND: `prpr/xvtest0.f` (driver holds extended)

Verified builds green (from Task 1):
- `bin/cbl_tout_qt` after `rm -rf xvue/qt/build` — exit 0
- `bin/cbl_tout` — exit 0 (legacy X11 untouched)
- `pp/ppxvtest0_qt` — exit 0, zero DRAW-XX warn-once lines

Verified human checkpoint:
- APPROVED 2026-04-11 (live-display Checks 1/2/3/5 + implicit Check 4 + headless Checks 6/7)

---
*Phase: 02-drawing-primitives-backing-pixmap*
*Plan: 04 (Wave 3)*
*Completed: 2026-04-11*
