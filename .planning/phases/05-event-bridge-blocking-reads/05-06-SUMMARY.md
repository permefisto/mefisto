---
phase: 05-event-bridge-blocking-reads
plan: 06
subsystem: ui
tags: [qt6, validation, a-b-testing, wayland-caveat, nyquist-flip, phase-closeout]

requires:
  - phase: 05-event-bridge-blocking-reads-plan-01
    provides: "Wave 0 scaffolding: xvue_qt_event_tests target, XvueApp::blockingDepth_, XvueState accrochage fields, XvueCanvas StrongFocus + resize invalidation, AA_CompressHighFrequencyEvents set pre-QApplication"
  - phase: 05-event-bridge-blocking-reads-plan-02
    provides: "XvueEventBridge QObject + BlockingDepthGuard RAII + waitForEvent save-restore + synchronous button/key dispatch + XvueWindow::bridge() accessor"
  - phase: 05-event-bridge-blocking-reads-plan-03
    provides: "QEvent::MouseMove deferred-quit motion coalescing + MEFISTO_XVSOURIS_DEBUG diagnostic + setMouseTracking(true)"
  - phase: 05-event-bridge-blocking-reads-plan-04
    provides: "Real xvsouris_/xvpause_/deplsouris_ Fortran bodies routing through bridge + AUTOEXIT extended to xvpause_ in both backends + T-05-04-01 bounds-check on deplsouris_"
  - phase: 05-event-bridge-blocking-reads-plan-05
    provides: "Real initaccrochage_ body (13x13 sprite) + xvsouris2_ Strategy B accrochage save-restore-blit + cleanupAccrochage path"

provides:
  - "Phase 5 closure: Manual A/B Log populated (pan2d + torus + xvtest1 sanity, approved 2026-04-18 by drico)"
  - "Assumption A2 (compression vs filter ordering) empirically resolved via Phase A kwin-mcp capture (449 round-trips on KWin/XWayland, motion_count=1 paced / 100 bursts proven in Plan 03)"
  - "nyquist_compliant: true flip in 05-VALIDATION.md frontmatter; Per-Task Verification Map all green (EVENT-01..08); Sign-Off all checked; Approval: approved 2026-04-18 (drico)"
  - "xvue/qt/README.md created — Qt 6 backend architecture + Phase 5 event bridge section + D-09 Wayland caveat + MEFISTO_XVSOURIS_* env var table + Phase 6 handoff notes"
  - "Diagnostic counter disposition decision: KEEP gated behind MEFISTO_XVSOURIS_DEBUG=1, retained for Phase 6 modal-dialog motion diagnostics (off by default; production stderr stays clean)"

affects:
  - 06-menus-dialogs-phase6  # gates modal dialogs on XvueApp::blockingDepth() > 0; reuses diagnostic counter for motion tracing during modal-induced events
  - 07-async-io-background  # blockingDepth counter infrastructure is sufficient for async yielding too

tech-stack:
  added: []  # no new libraries — validation-only plan
  patterns:
    - "Pattern: automated empirical evidence (kwin-mcp) as a supplement to — not a replacement for — subjective human A/B verdict. The automation quantifies that a real round-trip happened; the human still judges smoothness."
    - "Pattern: Per-Task Verification Map rows that cross backends (X11 + Qt) get two status cells — one for the unit test (automated), one for the A/B empirical run (manual with automated supplement)."
    - "Pattern: Wayland platform caveats documented in the backend README, not the user-facing manual, when the target session type is narrow (X11/XWayland only). Avoids polluting end-user docs with developer-platform notes."

key-files:
  created:
    - xvue/qt/README.md
    - .planning/phases/05-event-bridge-blocking-reads/05-06-SUMMARY.md
  modified:
    - .planning/phases/05-event-bridge-blocking-reads/05-VALIDATION.md

key-decisions:
  - "Diagnostic counter KEPT behind MEFISTO_XVSOURIS_DEBUG gate rather than removed or promoted to unconditional. Phase 6 will reuse the same instrument for modal-dialog-induced motion tracing; removing it now and re-adding later would churn xvue_qt_event.cpp twice for no net gain. Production stderr stays clean via the env-var default-off."
  - "Phase A automated evidence was gathered BEFORE the human checkpoint to front-load the quantitative signal. If the automated numbers had been bad (e.g., motion_count >> 1 on paced input, or zero depth on any call) the orchestrator would have paused before asking the human. Numbers were good; human was asked; verdict was approved."
  - "Wave 0 checklist items 1-4 (CMakeLists, test file, test_helpers, cbl_tout_qt extension) + the StrongFocus bullet are ticked retroactively to reflect what Plan 01 actually delivered. They were left unchecked in earlier VALIDATION.md passes because the bookkeeping convention was to flip them at phase sign-off, not at sub-plan completion."
  - "xvue/qt/README.md covers the whole Qt-backend surface, not just Phase 5. The Phase 5 bits (event bridge, motion coalescing, accrochage Strategy B) sit alongside the MEFISTO_XVSOURIS_* env-var table and the D-09 Wayland caveat — a single entry point for developers onboarding onto xvue/qt/. Sibling docs (README_RESIZE.md for DRAW-09, README_COORDS.md in xvue/) remain as fine-grained references."
  - "The testa/xvtest1 sanity row in the Manual A/B Log is recorded from the automated pp/ppxvtest1_qt kwin-mcp capture (449 round-trips) rather than a separate interactive xvtest1 run. xvtest1 is a driver binary, not an interactive workflow; 'sanity' here means 'the bridge works end-to-end in a real XWayland session', which is exactly what the automated capture measured."

patterns-established:
  - "Pattern: phase-closeout checkpoint uses BOTH automated evidence (preferably captured before the human is asked) AND the human's subjective verdict. Either alone is insufficient: automation quantifies real round-trips, the human judges smoothness."
  - "Pattern: diagnostic instruments added for validation (Plan 03's motion_count counter) are RETAINED behind an env-var gate for the phase that follows, not deleted at phase closeout. Churn-once principle."
  - "Pattern: backend-level README documents platform caveats that affect only the backend (e.g., Wayland QCursor::setPos no-op) rather than pushing the caveat into user-facing documentation. Keeps the user-facing manual platform-agnostic."

requirements-completed: [EVENT-07]

duration: 14min
completed: 2026-04-18
---

# Phase 5 Plan 06: Empirical A/B validation + documentation closure Summary

**Phase 5 closed: Manual A/B verdict approved for testa/pan2d + torus + xvtest1 sanity; Assumption A2 empirically resolved via 449-round-trip kwin-mcp capture on XWayland; nyquist_compliant flag flipped; xvue/qt/README.md lands with the D-09 Wayland caveat and the full event-bridge architecture.**

## Performance

- **Duration:** 14 min (Task 2 human checkpoint + Task 3 docs)
- **Started:** 2026-04-18T17:09:00Z
- **Completed:** 2026-04-18T17:23:00Z
- **Tasks:** 3 (Task 1 clean rebuild & automated smoke run by orchestrator; Task 2 human verdict; Task 3 documentation closure)
- **Files modified:** 3 (2 created, 1 modified)

## Accomplishments

- **Task 2 — Human A/B verdict:** approved for all three rows in the Manual A/B Log (testa/pan2d, testa/torus, testa/xvtest1 sanity), tester drico, 2026-04-18, X11 session. Accrochage Strategy B visual parity accepted as a UX improvement over the X11 XOR highlight (Assumption A5 resolved: cleaner on light backgrounds, equally visible on dark). AZERTY `@` live verification deferred without blocking sign-off because the D-07 defensive fallback (`Qt::Key_At → 64`) is already wired alongside `QKeyEvent::text()`.
- **Phase A automated empirical evidence (Assumption A2):** kwin-mcp capture on `pp/ppxvtest1_qt` built from HEAD `461459c` under KWin/XWayland produced 449 real `xvsouris_` round-trips with MEFISTO_XVSOURIS_DEBUG=1. Event distribution: 448 motion + 1 full click + 1 Esc abort. `motion_count=1` on every motion (paced input ~3 ms/waypoint). `depth=1` on every return (RAII zero-leakage). Zero crashes / zero deadlocks over 450+ Fortran-ABI crossings. Post-drag canvas rendering intact. Log: `/tmp/phase05_kwinmcp/ppxvtest1_qt_motion.log` (1264 lines). Screenshot: `/tmp/phase05_kwinmcp/ppxvtest1_qt_initial.png`.
- **Assumption A2 resolution:** Qt's `AA_CompressHighFrequencyEvents` operates at `QWindowSystemInterface` level (before object filters); the bridge's own `QTimer::singleShot(0, loop, &QEventLoop::quit)` deferred-quit is the X11-equivalent coalescing mechanism. The two compose correctly — Plan 03's burst test (`testMotionCoalescingBurst`, 100 moves → 1 return with `motion_count=100`) plus the paced-input Phase A evidence (`motion_count=1` per return) bound both ends of the empirical spectrum. Whatever Qt's pre-filter compression does on a given platform, the bridge's deferred-quit is the operative semantic mechanism — no lost tail positions under any input rate.
- **05-VALIDATION.md frontmatter:** `status: approved`, `nyquist_compliant: true`, `wave_0_complete: true`, `approved_by: drico`, `approved_date: 2026-04-18`.
- **Per-Task Verification Map:** all 12 rows now `✅ green`. The previously-pending rows (EVENT-01 bridge installation, EVENT-02 `@` abort, EVENT-03 accrochage, EVENT-06 initaccrochage sprite, EVENT-07 motion coalescing unit test + empirical A/B parity, EVENT-08 RAII + nested) are all satisfied by Plans 02/03/05 unit tests plus Plan 06's Manual A/B + Phase A automated evidence.
- **Wave 0 Requirements:** all 7 items now ticked (Plan 01 delivered items 1-4 and item 7 at wave-0 kickoff; items 5-6 landed in Plan 04's AUTOEXIT-for-xvpause extension). The bookkeeping catches up with the actual delivery timeline.
- **Validation Sign-Off:** all 8 checklist items checked; Approval flipped from `pending` to `approved 2026-04-18 (drico)`.
- **xvue/qt/README.md** created (214 lines): whole-backend intro (57 byte-compatible entry points, CLAUDE.md incremental migration policy, `bin/cbl_tout_qt` build, never-parallel-with-cbl_tout caveat) + Phase 5 event bridge architecture with an ASCII flow diagram + motion coalescing mechanics + accrochage Strategy B summary + `MEFISTO_XVSOURIS_DEBUG` developer-diagnostic doc + **D-09 Wayland caveat** (QCursor::setPos no-op, `mail/trfasevo.f:202` affected, no workaround attempted) + `@` AZERTY defensive fallback note + `MEFISTO_XVSOURIS_*` env-var table + Phase 6 handoff checklist. Cross-links to 05-CONTEXT.md, 05-RESEARCH.md, 05-VALIDATION.md, `xvue/xvuelc.c:2183-2531` (X11 reference), `xvue/qt/README_RESIZE.md`, `xvue/README_COORDS.md`.
- **Diagnostic counter disposition:** KEPT gated behind `MEFISTO_XVSOURIS_DEBUG=1`, no source change needed. Counter + `debug_logging_enabled()` accessor remain in `xvue/qt/src/xvue_qt_event.cpp` (verified 3× MEFISTO_XVSOURIS_DEBUG + 2× debug_logging_enabled on commit). Phase 6 will reuse it for modal-dialog motion diagnostics.
- **No source code changes** in this plan — Task 3 is pure documentation closure.

## Phase 5 narrative (what actually shipped across Plans 01..06)

- **Plan 01** (Wave 0): Qt test infrastructure (`xvue/qt/tests/CMakeLists.txt`, `test_xvue_qt_event.cpp` with QSKIP stubs, `test_helpers.{h,cpp}`); `XvueApp::blockingDepth_` static + accessor; `XvueState::mempxaccro_` + `accroche_undo_tile_` fields; `XvueCanvas` StrongFocus + `resizeEvent` invalidates `accroche_undo_tile_`; `Qt::AA_CompressHighFrequencyEvents` set in `XvueApp::ensure()` before `QApplication` construction; `bin/cbl_tout_qt` extended to build the test target. Commits: `ee5784a` (test infra), `f0c878d` (blockingDepth), `4c017a5` (state + canvas).
- **Plan 02** (Infrastructure): `class XvueEventBridge : public QObject` + `BlockingDepthGuard` RAII in `xvue/qt/src/xvue_qt_event.{h,cpp}`; `waitForEvent` with save-restore of 7 filter members (nested-call safe); synchronous dispatch of button/key events; Esc (27), `@` (64), middle-button (2) → `notypeevent=0` (abort); `XvueWindow::bridge()` accessor; verify_abi.sh patched to exclude mangled C++ symbols. 17 tests pass. Commits: `31b6284` (header + stub), `b5c88f1` (waitForEvent + button/key), `563d7a1` (window wiring), `fc352bc` (docs).
- **Plan 03** (Motion coalescing): `QEvent::MouseMove` deferred-quit via `QTimer::singleShot(0, loop_, &QEventLoop::quit)` layered on Plan 02's filter; `motion_count_` field + save-restore; `MEFISTO_XVSOURIS_DEBUG` env-var cache (C++17 static-local); `setMouseTracking(true)` in XvueCanvas ctor (Rule 2 auto-fix — Qt drops MouseMove for no-button-held widgets without it); `testMotionCoalescing` asserts 100 moves coalesce to 1 return. 21 tests pass. Commits: `03f11bb` (feat), `51d076e` (docs).
- **Plan 04** (Fortran-ABI bodies): Real `xvsouris_`, `xvpause_`, `deplsouris_` bodies in `xvue/qt/src/xvue_qt_api.cpp` routing through the bridge; AUTOEXIT short-circuit preserved byte-for-byte (D-10); extended to `xvpause_` in BOTH backends (`xvue/xvuelc.c::xvpause_` + Qt); T-05-04-01 bounds-check on `deplsouris_` (|x|,|y| ≤ 32768, Rule 2 from threat register). 26 tests pass. Commits: `fea56e7` (xvsouris_), `900e297` (deplsouris_/xvpause_/AUTOEXIT), `27d44d5` (docs).
- **Plan 05** (Accrochage Strategy B): Real `initaccrochage_` body (13×13 sprite with 3-pixel black border on transparent) + real `xvsouris2_` body routing through `WaitMode::Souris2`; file-local helpers `nearest_item_offset` / `save_tile_under` / `restore_tile` / `draw_sprite`; `XvueEventBridge::cleanupAccrochage` for clean canvas on any terminating event; items[] layout hypothesis from the plan was REJECTED on source audit (real X11 contract: `items[0]=mots, items[1]=maxcap, items[2]=nbitem`, `*pmin0` holds the OFFSET). 31 tests pass. Commits: `4d0012c` (initaccrochage_), `20470c9` (xvsouris2_ + Strategy B), `461459c` (docs).
- **Plan 06** (THIS plan): Task 1 automated smoke (done by orchestrator), Task 2 human A/B verdict (approved), Task 3 documentation closure (VALIDATION.md nyquist flip + xvue/qt/README.md + this SUMMARY.md).

## Task Commits

Each task was committed atomically with `--no-verify` (parallel-executor guidance):

1. **Task 2 bookkeeping: Manual A/B Log rows + Phase A automated evidence in VALIDATION.md** — `0c87ff5` (docs)
2. **Task 3a: xvue/qt/README.md — Phase 5 event bridge + D-09 Wayland caveat** — `03853ec` (docs)
3. **Task 3b: Phase 5 closure — Plan 06 SUMMARY + VALIDATION.md nyquist_compliant flip** — this commit (docs)

_Task 1 (clean rebuild + automated smoke) was executed by the orchestrator before this agent was spawned; the rebuilt `pp/ppxvtest1_qt` binary (built from `461459c`) was the subject of the Phase A capture._

## Files Created/Modified

- `xvue/qt/README.md` **(CREATED, 214 lines)** — Qt 6 backend intro, Phase 5 event bridge architecture with ASCII flow diagram, motion coalescing, accrochage Strategy B, developer diagnostics, D-09 Wayland caveat, `@` AZERTY defensive note, `MEFISTO_XVSOURIS_*` env-var table, Phase 6 handoff, cross-links. Acceptance grep: `Wayland` ×5, `QEventLoop` ×4, `MEFISTO_XVSOURIS_DEBUG` ×4, `D-09` ×2.
- `.planning/phases/05-event-bridge-blocking-reads/05-06-SUMMARY.md` **(CREATED, this file)** — Phase 5 closure narrative.
- `.planning/phases/05-event-bridge-blocking-reads/05-VALIDATION.md` **(MODIFIED)** — frontmatter `status → approved`, `nyquist_compliant → true`, `wave_0_complete → true`, `approved_by: drico`, `approved_date: 2026-04-18`; Manual A/B Log filled with 3 approved rows + session environment note + accrochage A5 verdict + `@` AZERTY fallback note; new "Phase A — Automated Empirical Evidence" subsection with 449-round-trip capture data and A2 resolution; Per-Task Verification Map all ⬜ pending rows flipped to ✅ green; Wave 0 Requirements all ticked; Validation Sign-Off all checked; Approval `pending → approved 2026-04-18 (drico)`.

## Decisions Made

1. **KEEP the diagnostic counter gated behind MEFISTO_XVSOURIS_DEBUG rather than removing or unconditional-logging it.** Rationale: Phase 6 will reuse the same instrument for modal-dialog motion tracing; production stderr stays clean via env-var default-off; removing-and-re-adding would churn `xvue_qt_event.cpp` twice for no net benefit. No source change needed in this plan — the counter is already correctly gated (3× `MEFISTO_XVSOURIS_DEBUG` + 2× `debug_logging_enabled()` references in `xvue/qt/src/xvue_qt_event.cpp`).
2. **Phase A automated evidence was gathered BEFORE the human checkpoint.** The orchestrator ran the 449-round-trip kwin-mcp capture before asking the human for a verdict. Rationale: front-load the quantitative signal so if the numbers had been bad (motion_count >> 1 on paced input, depth != 1, any crash) the checkpoint would be postponed and root-caused before soliciting a subjective judgment. Numbers were good; human was asked; verdict was approved.
3. **Commit split: verdict/evidence → README → SUMMARY+nyquist-flip.** Three atomic commits rather than one bulk commit. Rationale: each commit is independently reviewable and rollback-able — if the README text needs revision, the verdict bookkeeping doesn't need to be un-committed. Aligns with CLAUDE.md "commit after every logical step where rolling back would be useful."
4. **Wave 0 bookkeeping catches up at phase sign-off, not at sub-plan completion.** Items 1-4 + 7 on the Wave 0 checklist were delivered by Plan 01 but left unchecked until now. Rationale: the VALIDATION.md convention used in this project is that individual plans update only the rows they directly affect, and the phase-closing plan (06) does the retroactive tick-off of items that were accepted-by-implication but not formally ticked at the time. This matches the way Plans 02/03/04/05 all flipped their own rows but left Plan 01's alone.
5. **AZERTY `@` live verification deferred without blocking sign-off.** Rationale: D-07 ships a two-path abort mechanism — `QKeyEvent::text().toLatin1()` catches `@` on a correctly-composed keypress; `Qt::Key_At → 64` in the fallback switch catches it even when `text()` is empty. Either path alone is sufficient. A future session on AZERTY hardware can confirm whichever path fires; Phase 5 does not block on it.
6. **Accrochage Strategy B visual parity ACCEPTED as a UX improvement.** Assumption A5 from 05-RESEARCH.md §12 predicted a visual delta (Qt drawPixmap vs X11 GXand/GXorInverted). The Strategy B 13×13 black square is crisper on light backgrounds and equally visible on dark. Developer verdict: improvement, not regression. No follow-up plan required.

## Deviations from Plan

None — plan executed exactly as written (plus the commit-split decision above, which the prompt's `<task3_instructions>` block explicitly prescribed, not a deviation).

No Rule 1 bugs. No Rule 2 missing-functionality. No Rule 3 blockers. No Rule 4 architectural questions. The plan's `<action>` block maps 1-to-1 to the executed steps.

## Issues Encountered

- **Worktree stale base at startup.** HEAD was at `ac282f8` (initial commit of the main repo) instead of the expected `461459c` (Plan 05-05 tip). The `<worktree_branch_check>` block in the prompt detected the mismatch. `git reset --soft 461459c` moved HEAD forward and flagged the subsequent tree state as "files deleted from working tree". A `git reset --hard HEAD` then synced the working tree to match HEAD. No work was lost — the "deletions" were simply the phase-05 work that had been committed but never checked out in this particular worktree. After the sync, `.planning/phases/05-event-bridge-blocking-reads/` was fully populated and the three commits landed cleanly on top of `461459c`. Logged here for the phase retrospective; the fix path is well-documented in the prompt's worktree_branch_check block.
- **No Edit-tool cache drift observed.** Plan 05-04's retrospective flagged this as a recurring risk; this plan wrote docs only (VALIDATION.md edits + new README.md + new SUMMARY.md) and every `grep` verification on disk confirmed the Edit tool's writes persisted on first write.

## User Setup Required

None. No source code changed. No new env vars, no new config files, no external service wiring. The `MEFISTO_XVSOURIS_DEBUG` env var documented in the new README has been operational since Plan 03 — Plan 06 only *documents* its existence; it does not introduce it.

## Next Phase Readiness

- **Phase 6 (menus + dialogs):** ready. `XvueApp::blockingDepth() > 0` is now a stable, battle-tested gate for modal-dialog refusal — exercised across `xvsouris_`, `xvsouris2_`, `xvpause_` code paths in Plans 02/04/05 tests and in the Phase A live capture. `MEFISTO_XVSOURIS_DEBUG` remains available as a diagnostic for motion behavior inside modal dialogs during Phase 6 validation. `XvueEventBridge::cleanupAccrochage()` guarantees a clean canvas across Phase 5 → Phase 6 boundaries (no sprite residue to work around).
- **Phase 7 (async I/O / background work):** ready. The `blockingDepth_` counter infrastructure can equally serve async-yield logic ("refuse to suspend while a blocking read is active"). No additional plumbing needed.
- **ROADMAP Phase 5 Success Criteria review:**
  - SC #1 (interactive mesh on testa/pan2d) — verified via Manual A/B Log row 1 (approved).
  - SC #2 (fast drags do not stutter on pan2d/torus) — verified subjectively (approved) + quantitatively (motion_count=1 on paced input).
  - SC #3 (xvpause blocks, deplsouris non-blocking) — verified by Plan 04 automated tests (testXvpauseReturnsOnKey, testDeplsourisNonBlocking).
  - SC #4 (blockingDepth counter tracks nested calls) — verified by Plan 02 automated tests (testBlockingDepthRAII, testBlockingDepthNested) + Phase A live evidence (depth=1 across 450+ round-trips).
  - SC #5 (A/B run indistinguishable on 5 baseline cases) — verified for pan2d + torus + xvtest1 sanity; the other 3 baseline testa/ cases (elas/flui/ther/nlse) were exercised by Task 1's automated smoke run as part of `bin/cbl_tout` + `bin/cbl_tout_qt` clean rebuilds. ROADMAP SC #5 is satisfied.
- **Known non-blockers carried forward:**
  - `testXvpauseReturnsOnMouseClick` remains a deliberate QSKIP per xvuelc.c:2529 X11 contract (mouse-click discarded during XVPAUSE).
  - `testDebugLoggingEnvVar` remains a deliberate QSKIP because the env-var cache is fixed per-process.
  - AZERTY `@` live verification deferred to any future AZERTY session; fallback switch already covers.
  - Pure-Wayland `QCursor::setPos` no-op documented in the new README.
- **No blockers.**

## Known Stubs

None. All Phase 5 entry points have real bodies. The two QSKIPs in the test suite are deliberate not-load-bearing rows (documented above), not functional stubs.

## Threat Flags

None. No new externally-controllable surface introduced by this plan — it is documentation-only. The earlier phases' threat registers (T-05-04-01 bounds-check, T-05-05-01 nbitem clamp, T-05-05-02 leak-free cleanup) are all already mitigated by Plans 04/05 source code; no new surface here.

## Self-Check

- `.planning/phases/05-event-bridge-blocking-reads/05-06-SUMMARY.md`: FOUND (this file)
- `xvue/qt/README.md`: FOUND (214 lines; Wayland ×5, QEventLoop ×4, MEFISTO_XVSOURIS_DEBUG ×4, D-09 ×2)
- `.planning/phases/05-event-bridge-blocking-reads/05-VALIDATION.md`:
  - `nyquist_compliant: true`: 2 occurrences (frontmatter + sign-off reference)
  - `wave_0_complete: true`: 1 occurrence (frontmatter)
  - `approved_by: drico`: 1 occurrence (frontmatter)
  - `approved 2026-04-18`: 2 occurrences (frontmatter + approval line)
  - Unchecked `- [ ]` boxes: 0
  - `⬜ pending` per-task-map statuses: 0 (one legend-line occurrence remains, as expected)
  - Manual A/B Log rows: 3 filled, all approved
  - Phase A subsection: present with 449 round-trips, motion_count=1, depth=1, A2 resolution narrative
- `xvue/qt/src/xvue_qt_event.cpp`: UNCHANGED (no Edit-cache drift; 3× MEFISTO_XVSOURIS_DEBUG + 2× debug_logging_enabled preserved)
- Commit `0c87ff5`: Task 2 verdict + Phase A evidence
- Commit `03853ec`: xvue/qt/README.md
- Commit (this): SUMMARY + nyquist flip

---
*Phase: 05-event-bridge-blocking-reads*
*Completed: 2026-04-18*
