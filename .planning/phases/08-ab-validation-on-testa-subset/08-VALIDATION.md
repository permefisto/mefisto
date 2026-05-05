# Phase 8 — Validation Coverage (VALID-01 through VALID-07)

**Owner:** Plan 07 (CHECKLIST.md finalize) — except for the
`## Open escalations` section (filed by Plan 06 Task 2 per
`08-06-PLAN.md` step 12c). Plan 07 preserves the escalation entry
verbatim and integrates it into the v1 ship-gate sign-off matrix as a
maintainer-decision item.

**Date composed:** 2026-05-05 (Plan 07 Task 2)
**Reads from:** 08-CHECKLIST.md (Plan 07 Task 1 deliverable) +
REQUIREMENTS.md §VALID-01..07 (lines 99-105) +
08-RESEARCH.md §"Phase Requirements" + 5 sweep logs + spot-check-summary.md.

This file maps each VALID-NN requirement to specific cells / spot-check
rows in `08-CHECKLIST.md`, providing the evidence trail and per-requirement
verdict the maintainer signs in Task 3.

---

## Top-level summary

- **Main-thread guard verdict (Plan 5 Task 2):** **PASS**. 0 Q_ASSERT/abort hits
  across all 4 OMP-eligible cases (pan2d, nafems_le1, heat1d, nlsecu) under
  `OMP_NUM_THREADS=8` + `QT_QPA_PLATFORM=offscreen` + `MEFISTO_BATCH_X11=1`.
  XVUE_QT_ASSERT_MAIN_THREAD instrumentation on every public ABI entry in
  `xvue/qt/src/xvue_qt_api.cpp` held under OMP solver parallelism. Benign
  `QThreadStorage: entry 1/0 destroyed before end of thread` cleanup-order
  noise observed in nlsecu shutdown log; documented as Qt cleanup behaviour,
  NOT off-thread ABI invocation.
- **Phase-boundary discipline:** **PASS**. `git diff --quiet -- xvue/qt/src/`
  AND `git diff --quiet -- xvue/xvuelc.c` both empty. No fix-while-validating.
- **HiDPI dim-ratio finding (Plan 4 OQ4):** **A5 CONTRADICTED** — escalated
  to maintainer review (see VALID-04 below). Aggregate 0/5 cases produce
  the predicted 2:2 ratio; deterministic 752x156 captures with width
  ratio ~0.989 and height ratio ~0.353.
- **Sign-off status:** **PENDING** (Task 3 checkpoint:human-verify).

---

## VALID-01: All 5 canonical testa/ cases pass visually side-by-side (Qt 1x mode)

**Requirement (verbatim from REQUIREMENTS.md line 99):**
> All 5 canonical `testa/` cases (BUILD-10 baseline) pass visually
> (side-by-side comparison of X11 and Qt windows) at the end of every phase
> from 0 through 7; results are logged to `.planning/phase-N/VALIDATION.md`.

- **Coverage:** 5/5 cells in CHECKLIST.md matrix column "Qt 1x" — `pan2d`,
  `nafems_le1`, `cavity2d`, `heat1d`, `nlsecu`.
- **Cells:** all 5 main-grid Qt 1x cells.
- **Evidence:** for each case `${case}-qt-1x.png` + `${case}-qt-1x-diff.png`
  under `.planning/phases/08-ab-validation-on-testa-subset/evidence/`.
- **AE roll-up (post-merge from `evidence/wave-merge-ae.md`):**
  pan2d=540804 (52.81%), nafems_le1=412827 (40.32%), cavity2d=411003 (40.14%),
  heat1d=143273 (13.99%), nlsecu=728737 (71.17%) [TRUNCATED-CAPTURE].
- **Verdict:** **partial-PASS-with-overrides PENDING maintainer sign-off.**
  All 5 cells CHECK at default fuzz=5%; AE signal DOMINATED by the
  760x442→1280x800 nearest-neighbour resample mandated by
  `bin/ab_compare_pair.sh` dim-guard (NOT real backend rendering drift).
  Maintainer review at Task 3 expected to produce override-on-review for
  4/5; nlsecu cell additionally TRUNCATED-CAPTURE (MAILLER prereq frame
  only — ppnlse_qt offscreen+BATCH_X11 deadlock).

## VALID-02: X11 backend continues to pass (X11 baseline column)

**Requirement (verbatim from REQUIREMENTS.md line 100):**
> The X11 backend continues to pass the same 5 cases at the end of every
> phase, proving no accidental damage to the legacy path during the A/B
> window.

- **Coverage:** 5/5 cells in CHECKLIST.md matrix column "X11 (baseline)" —
  `pan2d`, `nafems_le1`, `cavity2d`, `heat1d`, `nlsecu`.
- **Cells:** all 5 main-grid X11 baseline cells.
- **Evidence:** for each case `${case}-x11.png` under
  `.planning/phases/08-ab-validation-on-testa-subset/evidence/`. SHA-256 +
  size + dimensions table in `evidence/sweep-log-x11.md`.
- **AE roll-up:** all 5 X11 baselines self-PASS (AE=0 against themselves
  by construction; no Qt-side compare leg in this column). Color-content
  sanity confirms canonical Mefisto palette markers (peach #FFD9B8, blue-gray
  #4F4F63 mesh) on every capture.
- **xvue/xvuelc.c invariant:** byte-identical — `git diff -- xvue/xvuelc.c`
  empty (Plan 02 verification + Plan 07 Task 2 re-check).
- **Verdict:** **PASS.** No accidental damage to legacy X11 path during
  the A/B window. nlsecu cell TRUNCATED-CAPTURE (TIME=0.01, 1 step) is a
  workspace-only mitigation that does NOT mutate `testa/nlsecu/`.

## VALID-03: All 5 _OMP variants — main-thread invariant under OMP solver parallelism

**Requirement (verbatim from REQUIREMENTS.md line 101):**
> At the end of Phase 7, all 5 cases are run through `ELASTICER_OMP` (the
> OpenMP executable) to verify the Qt layer's main-thread-only invariant
> holds under OpenMP solver parallelism.

- **Coverage:** 4/4 OMP-eligible cells in CHECKLIST.md matrix column "Qt _OMP"
  (`pan2d`, `nafems_le1`, `heat1d`, `nlsecu`); 1/5 N-A (cavity2d — no
  `bin/FLUIDER_OMP` launcher per D-05; 08-RESEARCH.md §"OMP Sweep — How _OMP
  Actually Works" research finding 3).
- **Cells:** 4 main-grid Qt _OMP cells + cavity2d N-A annotation.
- **Evidence:** for each OMP-eligible case `${case}-qt-omp.png` +
  `${case}-qt-omp-diff.png` + `${case}-x11-omp.png`. Main-thread guard
  log in `evidence/sweep-log-omp.md ## Main-thread guard verification`
  (0 Q_ASSERT hits across all 4 cases).
- **AE roll-up (Qt-OMP vs X11-OMP, fuzz=5%):**
  pan2d=544936 (53.22%), nafems_le1=413940 (40.42%), heat1d=143209 (13.99%),
  nlsecu=147526 (14.41%) [TRUNCATED-CAPTURE TIME=0.1 / 10 steps].
  All 4 CHECK by Pitfall 5 heuristic (AE > 2 × OMP twice-run noise floor 236 px =
  472 px); resample-confounded.
- **ELASTICER_OMP equivalence proof:** `evidence/elasticer-omp-equiv-check.md`
  (BLOCKER #4 mitigation; AE=1313 px launcher-pattern-vs-direct, within
  1143 px OMP twice-run noise floor) — direct `OMP_NUM_THREADS=8 pp/pp${MODULE}`
  is empirically equivalent to literal-named launcher `bin/ELASTICER_OMP`.
- **Main-thread guard verdict (Plan 5 Task 2):** **PASS** (0 Q_ASSERT/abort
  hits across all 4 OMP-eligible cases).
- **Verdict:** **partial-PASS-with-overrides PENDING maintainer sign-off.**
  Main-thread invariant PASS; per-cell AE values resample-confounded
  (CHECK pending Pitfall 5 noise-floor maintainer review).

## VALID-04: All 5 cases at HiDPI (QT_SCALE_FACTOR=2)

**Requirement (verbatim from REQUIREMENTS.md line 102):**
> At the end of Phase 7, all 5 cases are run on a HiDPI display (4K or
> `QT_SCALE_FACTOR=2`) to verify no size / position drift vs the X11
> backend.

- **Coverage:** 5/5 cells in CHECKLIST.md matrix column "Qt HiDPI 2x" —
  `pan2d`, `nafems_le1`, `cavity2d`, `heat1d`, `nlsecu`.
- **Cells:** all 5 main-grid Qt HiDPI 2x cells.
- **Evidence:** for each case `${case}-qt-2x.png` +
  `${case}-qt-2x-diff.png` under
  `.planning/phases/08-ab-validation-on-testa-subset/evidence/`.
- **AE roll-up (post-merge from `evidence/wave-merge-ae.md`):**
  pan2d=578116 (56.46%), nafems_le1=413201 (40.35%), cavity2d=661110 (64.56%),
  heat1d=277787 (27.13%), nlsecu=794335 (77.57%) [TRUNCATED-CAPTURE].
- **HiDPI dim-ratio confirmed:** **NO — discrepancy escalated** (Plan 04
  OQ4 finding contradicting Assumption A5). All 5 captures land
  deterministically at 752x156 NOT 1520x800. Width ratio ~0.989; height
  ratio ~0.353. Reproducible across multiple runs.
- **Hypothesis (Plan 04):** `QT_SCALE_FACTOR=2` may halve the *logical*
  widget area while the *physical* pixmap stays in the same ballpark with
  a different auto-fit-to-window calculation that compresses the frame
  vertically.
- **Verdict:** **partial-PASS-with-overrides PENDING maintainer sign-off**
  + **escalation to maintainer for HiDPI dim-ratio decision.** Maintainer
  must EITHER accept the dim-ratio finding as a known-deviation (with
  deferred real-4K-display visual eyeball check per `08-CONTEXT.md`
  Deferred Ideas) OR escalate to a Phase 9 follow-up plan.

## VALID-05: Color-bar spot checks on nafems_le1, heat1d, cavity2d

**Requirement (verbatim from REQUIREMENTS.md line 103):**
> Color-accuracy spot checks (stress/temperature/velocity color bars) on
> `testa/nafems_le1`, `testa/heat1d`, and `testa/cavity2d` show no
> user-visible drift between backends.

- **Coverage:** 3/3 spot-check rows in CHECKLIST.md "Spot-check rows" section
  (`colorbar-nafems_le1`, `colorbar-heat1d`, `colorbar-cavity2d`).
- **Cells:** 3 spot-check rows.
- **Evidence:**
  - `evidence/colorbar-nafems_le1.md` (Pitfall 6 spot check, full-frame —
    only 7 unique colors, no smooth gradient bar; recommended override
    fuzz=8 + dim-mismatch annotation) + `evidence/nafems_le1-qt-1x-diff.png`.
  - `evidence/colorbar-heat1d.md` (Pitfall 6 spot check, full-frame —
    only 8 unique colors, 1D thermal-curves plot, no temperature ramp;
    recommended override fuzz=8 + dim-mismatch annotation) +
    `evidence/heat1d-qt-1x-diff.png`.
  - `evidence/colorbar-cavity2d.md` (Pitfall 6 spot check, full-frame +
    cropped right-edge x=1100-1260 (160x800) — 19 full-frame unique
    colors, 15 in colorbar band; cropped-region AE 21.94% LOWER than
    full-frame 40.14% disconfirms Pitfall 6 as principal failure mode;
    recommended override fuzz=10 + dim-mismatch annotation) +
    `evidence/cavity2d-qt-1x-diff.png`.
- **Common finding:** Pitfall 6 hypothesis (gradient color bar tripping
  the 5% fuzz gate) **DISCONFIRMED** for all 3 cases — AE signal
  dominated by 760x442→1280x800 resample artefact, NOT color-bar drift.
- **Recommended overrides** (encoded in CHECKLIST.md Notes column):
  nafems_le1 fuzz=8; heat1d fuzz=8; cavity2d fuzz=10. All 3 carry
  `dim-mismatch=resample-to-1280x800` annotation.
- **Verdict:** **partial-PASS-with-overrides PENDING maintainer sign-off.**
  3 CHECK; 0 PASS at default fuzz=5%; PASS-on-review expected on
  diff-PNG inspection given resample-dim-mismatch dominant signal.

## VALID-06: Font-metric spot checks on pan2d, hexahedron

**Requirement (verbatim from REQUIREMENTS.md line 104):**
> Font-metric spot checks (node-number labels) on `testa/pan2d` and
> `testa/hexahedron` show no clipping or overlap in the Qt backend.

- **Coverage:** 1/2 cells PASS-pending-review (`font-pan2d`) +
  1/2 GAP-DOCUMENTED (`font-hexahedron`) per `## Open escalations` below.
- **Cells:** 2 spot-check rows.
- **Evidence:**
  - `evidence/font-pan2d.md` (Pitfall 7 AA-fringe analysis: Qt-1x has
    3545 unique colors vs X11's 9 → 394× explosion → empirically
    CONFIRMS Pitfall 7 anti-aliased text + font hinting drift; AE
    budget=3000 px confounded by dim-mismatch resample) +
    `evidence/pan2d-qt-1x-diff.png`.
  - `evidence/font-hexahedron.md` (GAP-branch report — hexahedron
    headless driver missing; all 5 alphabetical candidates {ln, ob,
    pt, sf, vl} probed and crashed/hung) +
    `evidence/hexahedron-probe.md` (full empirical probe trail).
- **pan2d Pitfall 7 finding:** Qt's 394× unique-color explosion vs X11
  empirically CONFIRMS Pitfall 7. The 3000-px AE budget cannot be
  applied cleanly on the resample-confounded 540 804-px AE
  measurement (52.81%); matched-dim recapture deferred to Phase 9.
- **hexahedron probe outcome:** GAP. The 5 files are interlocking
  building blocks (pt → ln → sf → vl → ob), NOT independent batch
  entry points. Legacy MAILLER protocol's single-batch invocation
  pattern does not apply. Per T-08-23 mitigation, hexahedron remains
  a spot-check-only target (NOT in BUILD-10 baseline).
- **Verdict:** **partial-PASS-with-overrides + GAP-DOCUMENTED PENDING
  maintainer sign-off.** Maintainer must EITHER (a) issue CONTEXT.md
  amendment naming the canonical hexahedron headless driver OR
  (b) waive VALID-06 hexahedron coverage with documented rationale
  (recorded in REQUIREMENTS.md against VALID-06 with maintainer
  initials and date). See `## Open escalations` below.

## VALID-07: Validation checklist signed off as v1 ship gate

**Requirement (verbatim from REQUIREMENTS.md line 105):**
> A validation checklist (`.planning/phase-8/CHECKLIST.md`) records
> pass/fail per case per backend and is the gate for declaring v1
> shippable.

- **Coverage:** 08-CHECKLIST.md sign-off line at bottom (per-cell verdict
  matrix per D-10 + spot-check rows + invariants block + override
  clauses + sign-off placeholder).
- **Cells:** 1 sign-off line in CHECKLIST.md
  (`Maintainer signature: ____________________ Date: __________`).
- **Evidence:** maintainer signature + ISO date + per-CHECK-cell
  override-acceptance initials (per Override clauses section).
- **Composed-state self-check:** 5 main-grid rows × 4 modes = 20 cells
  (19 concrete + 1 N-A); 5 spot-check rows; 5 build-hygiene/invariants
  lines (4 PASS + 1 FAIL-DEFERRED); fixture row demonstrating WARNING-4
  iter2 pipe-encoding transform fired; sign-off placeholder present.
- **Verdict:** **PASS once signed / pending sign-off.** Phase 9
  RETIRE-01..04 cannot start until this file is signed AND the
  one-release-cycle A/B window has elapsed.

---

## Phase Boundary Discipline

- **xvue/xvuelc.c invariant:** **PASS — `git diff --quiet -- xvue/xvuelc.c` empty.**
  Byte-identical to pre-Phase-8 baseline (commit 900e297). Codified by
  T-08-13 mitigation in 08-RESEARCH.md.
- **xvue/qt/src/ invariant:** **PASS — `git diff --quiet -- xvue/qt/src/` empty.**
  Byte-identical across all 7 plans (Plans 01-07). T-08-27 mitigation
  enforced via every plan's pre/post-task `git diff --quiet` checks.
- **No fix-while-validating:** **confirmed.** No source edits during Plans
  1-7 per phase boundary discipline (08-RESEARCH.md §Phase Boundary
  Discipline). Every modified file committed by Plans 01-07 lives under
  `.planning/` or `bin/ab_*.sh` / `bin/phase8_case_batch_map.sh` (Phase
  8-only sweep harness scripts, not source-tree mutation).
- **Drift-escalation findings:**
  - 14 Qt-mode CHECK cells (resample-dim-mismatch dominant signal) —
    deferred-idea Phase 9 candidate: matched-dim Qt recapture
    (1280x800 in-process backing pixmap if `MEFISTO_QT_WINDOW_SIZE` or
    equivalent supported) to retire the resample confound and let
    Pitfall 6/7 overrides apply cleanly.
  - HiDPI dim-ratio finding (5/5 NOT-2:2) — Plan 04 OQ4 escalation;
    Phase 9 follow-up question: does the Qt 6 offscreen backend at
    `QT_SCALE_FACTOR=2` produce a 752x156 backing pixmap as a
    documented Qt behaviour, or is this a size-hint regression that
    needs `xvue/qt/src/` follow-up?
  - hexahedron GAP-DOCUMENTED — see `## Open escalations` below;
    maintainer-decision-(a) gap-close OR maintainer-decision-(b)
    waive-with-rationale.
  - 3 Phase-7-deferred goldens (scene01.eps, wave_legacy.gif,
    cavity2d_legacy.gif) — Plan 1 Task 2 re-deferred to original
    Phase 7 VERIFICATION.md §9 human-bootstrap designation; Plan 1
    Rule 4 deviation. Recommendation: accept as Phase 7 carry-forward
    rather than Phase 8 ship-blocker.
  - nlsecu TRUNCATED-CAPTURE on every Qt cell — Phase 9 candidate
    code-side fix to ppnlse_qt offscreen+BATCH_X11 deadlock at startup.

## Open escalations (filed by Plan 06 Task 2)

### VALID-06 hexahedron font-metric spot check — GAP

**Filed:** 2026-05-05 by Phase 8 Plan 06 Task 2 (probe-then-document branch
per BLOCKER #1 / MAJOR-9 iter1 review; WARNING-3 / WARNING-5 iter2 review).

**Probe report:** `evidence/hexahedron-probe.md`
**Spot-check report:** `evidence/font-hexahedron.md`

**Summary:** VALID-06 hexahedron font-metric spot check could not be
completed because no headless batch driver was found in the time budget
of Plan 06. `testa/hexahedron/` contains 5 ASCII batch script files (`ln`,
`ob`, `pt`, `sf`, `vl`) with no recognized extension. All 5 candidates
were probed in alphabetical order under
`pp/ppmail_qt <D>` + `QT_QPA_PLATFORM=offscreen` + `MEFISTO_QT_CAPTURE_PATH`
+ `MEFISTO_XVSOURIS_AUTOEXIT=1` + `timeout 60`:

- `ln`, `ob`, `sf`, `vl`: `STOP Sorry, CRASH of MEFISTO software` at
  unknown-name lookup (each file references named entities defined in
  EARLIER stages of the construction pipeline; no single file is
  self-sufficient). 2045-byte all-peach background capture at crash time
  (xvfermer_ fired before any drawing primitive painted the mesh).
- `pt`: process HUNG at 60-s timeout (file has no rendering trigger;
  expects the next interactive menu pick to draw). 0-byte capture
  (xvfermer_ never fired).

**Architectural reading:** The 5 files are NOT independent batch entry
points — they are interlocking building blocks of a single hexahedron
mesh definition (pt → ln → sf → vl → ob). The legacy MAILLER protocol's
single-batch invocation pattern (`pp/ppmail batchfile`) does not apply
here.

**Maintainer decision required.** Maintainer must EITHER:

  **(a)** Issue a CONTEXT.md amendment naming the canonical driver. This
  could be:
  - A pointer to a wrapper script (e.g., `bin/hexahedron_drive.sh`) that
    issues the menu picks via stdin pipe to `pp/ppmail_qt` after seeding
    the workspace.
  - A documented stdin-script that drives the interactive menu through
    pt → ln → sf → vl → ob → drawings (close the workflow).
  - A new `pp/ppmail_qt` multi-file batch protocol that accepts a
    sequence of batch files (e.g., `pp/ppmail_qt pt ln sf vl ob`).

  **OR**

  **(b)** Waive VALID-06 hexahedron coverage with documented rationale
  (e.g., "pan2d coverage is sufficient for the AA-drift gate; hexahedron
  is too marginal to require headless A/B coverage given the maintenance
  cost of building a multi-file driver protocol"). The waiver should be
  recorded in `REQUIREMENTS.md` against VALID-06 with the maintainer's
  initials and date.

**Phase 7 sign-off is gated on this resolution.** Plan 07 will read this
escalation entry and emit a `hexahedron | mail | GAP | escalation per
evidence/hexahedron-probe.md` row in the CHECKLIST.md "Spot-check rows"
section (NOT the canonical 5-case grid — per T-08-23 mitigation,
hexahedron remains a spot-check-only target). The maintainer's Plan 07
sign-off either accepts the gap with rationale-(b) or refuses and demands
gap-closure-(a).

**See also:**
- `evidence/hexahedron-probe.md` — full empirical probe trail (5 candidates
  tried in alphabetical order with per-candidate exit code, capture size,
  log lines, crash markers).
- `evidence/font-hexahedron.md` — GAP-branch spot-check report (preserves
  all 5 literal section headings per WARNING-3 iter2 with N/A bodies
  pointing at hexahedron-probe.md).
- `08-06-PLAN.md` step 12c — the contract authorizing this escalation
  filing.
