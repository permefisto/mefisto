---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: milestone
status: completed
stopped_at: Phase 9 complete — X11 backend retired; Qt 6 single graphics backend
last_updated: "2026-05-06T08:00:00.000Z"
last_activity: 2026-05-06
progress:
  total_phases: 17
  completed_phases: 16
  total_plans: 67
  completed_plans: 68
  percent: 100
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-04-10)

**Core value:** Every MEFISTO workflow that works today through X11 keeps working through the new Qt 6 interface, with Fortran solver code unchanged.
**Current focus:** Milestone v1.0 X11→Qt6 migration COMPLETE. Phase 9 retirement done; xvuelc.c + libX11 + ImageMagick + LVIDEO all retired. Qt 6 is single graphics backend. Outstanding: 3 P7 goldens DEFERRED (Phase 7 design defects in scene01_driver.f + interactive testa pipelines empirically falsified by Plan 09-08 cross-tag attempt; carried forward).

## Current Position

Phase: 9 (complete)
Plan: 09-09 done; 9/9 plans (1 prereq + 4 RETIRE-NN + 4 carry-forward); v1.0-pre-retire tag exists for rollback
Status: Phase 9 closed 2026-05-06 — milestone v1.0 X11→Qt6 migration COMPLETE. xvuelc.c + ccxvue deleted (RETIRE-01); libX11 stripped from cb scripts + cb*_qt → cb* renames + pp suffix collapsed (RETIRE-02); convertepsgif + png2eps + png2jpg + LVIDEO selectively retired with non-LVIDEO drawing logic preserved (RETIRE-03); README + LISEZMOI + CLAUDE.md updated for Qt 6 single-backend reality (RETIRE-04). Carry-forwards: matched-dim Qt window-size knob shipped (Plan 09-06; Phase 8 override #1 closed); ppnlse_qt classified case-(c) long-runtime not deadlock (Plan 09-07; Phase 8 override #5 documentation-only closure); P7 goldens cross-tag attempted and DEFERRED with empirical evidence (Plan 09-08; Phase 7 source defects); harness realpath fix + verify_pp_freshness shipped (Plan 09-09). All 3 grep gates green; ABI 58/58; bin/cbl_tout exit 0; pp/* binaries fresh.
Last activity: 2026-05-06

Progress: [██████████] 100% (45/45 plans, Phase 6.3 complete)

## Performance Metrics

**Velocity:**

- Total plans completed: 43
- Average duration: —
- Total execution time: —

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 00 | 4 | - | - |
| 01 | 3 | - | - |
| 03.1 | 4 | - | - |
| 02.1 | 1 | - | - |
| 03 | 4 | - | - |
| 05 | 6 | - | - |
| 06.0 | 7 | - | - |
| 6.1 | 3 | - | - |
| 06.2 | 5 | - | - |
| 06.4 | 3 | - | - |
| 06.5 | 3 | - | - |

**Recent Trend:**

- Last 5 plans: —
- Trend: —

*Updated after each plan completion*
| Phase 01-window-shell-xvueapp-xvuewindow-xvuecanvas P01 | 19min | 3 tasks | 10 files |
| Phase 01 P02 | 15min | 2 tasks | 1 files |
| Phase 01 P03 | ~60min | 3 tasks | 5 files |
| Phase 06.2 P04 | 23m15s | 2 tasks | 5 files |
| Phase 06.2 P05 | 23m01s | 2 tasks | 3 files |
| Phase 06.3 P01 | 12min | 2 tasks | 1 files |
| Phase 06.3 P02 | ~30min | 3 tasks | 10 files |
| Phase 06.3 P03 | ~25min | 2 tasks | 3 files |

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Init: Qt 6 (not Qt 5) chosen; Qt 5 explicitly rejected as maintenance-mode
- Init: CMake owns only `xvue/`; Fortran shell build stays unchanged
- Init: `extern "C"` names and signatures must be byte-identical to `xvuelc.c`
- Init: Parallel X11 build kept alive for one release cycle for A/B validation
- Init: 9-phase roadmap adopted from research SUMMARY.md (dependency-forced)
- [Phase 01-window-shell-xvueapp-xvuewindow-xvuecanvas]: Phase 1 scaffolding landed: XvueApp/Window/Canvas/State quartet + SHELL-03 verify_no_exec guard + SHELL-07 real macro body
- [Phase 01]: xvfond_ with no open window is a documented no-op in Phase 1 (XvueState owned by XvueWindow)
- [Phase 01]: Macro retrofit via single replace_all on 'static bool warned = false;' — 51 stubs updated in one shot preserving D-18 ordering
- [Phase 01-03]: D-08 revised — QApplication deliberately leaked at atexit (qapp_.release() not reset()). Destroying QApplication in atexit alongside libgfortran races and crashes on Linux/Qt 6.
- [Phase 01-03]: D-01 revised — xvinitgraphique_ uses bounded-loop exposure pump (QElapsedTimer + isExposed()) not a single processEvents. X11 needs MapRequest/ConfigureNotify/Expose sequence.
- [Phase 01-03]: D-06 addendum — xvfermer_ drains deleteLater queue with processEvents after window_.reset(). Prevents stale DeferredDelete events at atexit.
- [Phase 01-03]: SHELL-01/02/06 empirically validated by human visual check — two 800x600 "MEFISTO" windows, exit 0, HiDPI ~2x scaling confirmed.
- [Phase 03.1-01]: xtinit_ promoted from warn-once stub to proc(xvinitgraphique)() forwarding body — unblocks XVOUVRIR → XVINIT chain on Qt backend.
- [Phase 03.1-02]: 9 new cb* build scripts (5 legacy + 4 Qt) + cbl_tout/cbl_tout_qt wiring; Pitfall 6 closed (legacy cbxvtest0 was missing).
- [Phase 03.1-03]: xvinfo_ no longer writes *maxfonts — the legacy contract is a pure input/capacity parameter. Writing through a Fortran PARAMETER-sourced actual crashed under gfortran -O in xvtest1..4_qt.
- [Phase 03.1-03]: 03-04 Task 1 smoke exit code tolerance expanded to {0, 124, 143} — timeout --preserve-status forwards SIGTERM (143), not 124.
- [Phase 03.1-03]: Legacy xvtest0/xvtest1 pre-existing xvue/xvuelc.c crashes deferred to a future legacy-hardening pass (see 03.1-NOTES.md + deferred-items.md).
- [Phase 03-04 reopen close 2026-04-14]: The 2 "MISMATCH" Qt solver traces from the 2026-04-13 first reopen were NOT real Qt rendering gaps. Root cause: Debian sid `apt upgrade` pulled `libgfortran5 = 16-20260322-1` (gcc-16 snapshot) which exposed latent UB in MEFISTO Fortran runtime. Fix: `sudo apt install libgfortran5_15.2.0-9_amd64.deb` (cached) + `apt-mark hold libgfortran5` + rebuild via `/tmp/gfortran-14-shim` PATH (gcc-15/gfortran-15 packages were also removed by the downgrade). Plus 2 long-standing batch-file bugs in `testa/nafems_le1/` (mesher missing `1;` between `5;90;`; elas missing `15;` between `1;90;`) — fixed in-place. Final D-27: 12/12 PASS, 1 DEFERRED (nlsecu).
- [Phase 03-04 reopen close 2026-04-14]: Build-environment pinning is now an open deferred item — `libgfortran5` must stay on `15.2.0-9` until either Debian sid stabilizes gcc-16 or the latent Fortran UB sites (uninitialized `TPSINI` in `ther/thed1t.f`, FPE traps in elasticity stress paths) are properly initialized in source.
- [Phase 06.2-04]: ensureTopLevelMenu generalized with optional `insertBefore` QMenu* anchor — when non-null, routes through QMenuBar::insertMenu(anchor->menuAction(), m). 6.3/6.4/6.5 must call ensureTopLevelMenu(..., viewMenuForAnchor) so per-module menus land between File and View.
- [Phase 06.2-04]: Helper duplication preserved across elas/mail (per 06.2-02 patterns-established). Updated in lockstep this plan; future refactor may hoist to a shared header but is out of scope.
- [Phase 06.2-04]: testMenuOrder regression-guard pattern codified — slot collects QMenuBar::actions() filtered on non-null QAction::menu(), takes first 4 objectNames, QCOMPARE against {File, <Module>, View, Help}. 6.3/6.4/6.5 should add the same slot to their menu test binaries.
- [Phase 06.2-04]: 6.1 mail co-fix folded into the same plan (same root cause, same hazard inherited from 6.1) rather than split out; lockstep helper update is the explicit ask.
- [Phase 06.2-05]: XvueMenuFileParser::loadFor() now language-aware — when xvueIsEnglish() returns true and td/ma/<name> exists, prefer the EN tree; fall back to td/m/<name> otherwise. Single-point fix benefits 6.1 mail and the upcoming 6.3/6.4/6.5 waves automatically.
- [Phase 06.2-05]: QFile::exists() pre-check gates the td/ma/ path so missing or unreadable EN file silently falls through to td/m/<name>; preserves existing-test compatibility without errno-based branching.
- [Phase 06.2-05]: testBilingualLabelsEnglish slot QSKIPs on missing preconditions (MEFISTO env / anglais flag / td/ma/debuelas absent) so it ports to constrained CI envs; on the live dev rig it PASSes (not SKIPs), gating real regressions.
- [Phase 06.2-05]: xvueIsEnglish() static-local cache caveat preserved: toggling the anglais flag mid-process is not supported; testBilingualLabelsEnglish must run before any slot that consumes xvueT() (verified for elas binary).
- [Phase 06.3]: Top-5 toolbar = {2;, 5;, 7;, 10;, 99;} (PHYSICAL-DATA / STEADY-STOKES / IMPLICIT-NS / DRAW-VEL+PRESS / SAVE-QUIT)
- [Phase 06.3]: Menu taxonomy {File, Fluid, View, Help} (Fluid replaces Solve relative to 6.2)
- [Phase 06.3]: 6 new flui SVGs queued + 1 reuse (mesh-draw.svg from 6.1); Stokes/NS icons consolidated (5;,6; share solve-stokes.svg; 7;,8;,9; share solve-navier-stokes.svg)
- [Phase 06.3-02]: Strong registerFluiActions_stub_ ships clean on first build attempt — 6.1 force-link Rule 3 lesson held; nm pp/ppflui_qt | grep ' T registerFluiActions_stub_' = 1
- [Phase 06.3-02]: Keepalive pair placed alphabetically in xvue_qt_app.cpp (elas/flui/mail); 6.4 ther will land after mail; 6.5 nlse after ther
- [Phase 06.3-02]: Three modules now pass AUDIT_STRICT_ICONS=1 with no WARN — validator's module-aware resolver from 6.2 needed zero further edits for flui (proves the 6.2 Plan 02 generalization was correctly future-proofed)
- [Phase 06.3]: [Phase 06.3-03]: testMenuOrder + testBilingualLabelsEnglish replace 6.2 manual A/B checkpoint -- Phase 6.3 fully autonomous via codified QTest gates
- [Phase 06.3]: [Phase 06.3-03]: testHelpNoQueue tightened from blanket 'no Help lexicon' to allowlist {97;} (Auto-fix Rule 1) -- per-module Help-lexicon set drawn from LEXICON-AUDIT, 6.4/6.5 inherit pattern
- [Phase 06.3]: [Phase 06.3-03]: Pattern proven idempotent across THREE modules (mail/elas/flui); 6.4 ther + 6.5 nlse are mechanical near-copies; xvmodi.f and xvue_qt_window.cpp continue to need ZERO edits
- [Phase 08]: Dimension 8e (VALIDATION.md pre-existence) WAIVED for Phase 8. Phase 8 deliverable IS 08-VALIDATION.md (Plan 7 output) — verification phases produce VALIDATION.md, do not consume it. Workflow rule incompatible with phase nature; user-explicit waiver 2026-05-05.

### Roadmap Evolution

- Phase 03.1 inserted after Phase 3: Build xvtest1..4 driver infrastructure on both backends (Qt + legacy X11) to unblock Phase 3 A/B catch-up gate (URGENT)
- Phase 03.1 CLOSED 2026-04-12: 3 plans, ABI stable at 57, Qt drivers green, 03-04 Task 1 hardened harness landed in-place.
- Phase 02.1 inserted after Phase 2: Qt drawing-primitive A/B gaps (xvface color state, multi-object 3D composition, dashed pen style) (URGENT — surfaced by 03-04 Task 2 A/B gate failure on xvtest2/xvtest4)

### Pending Todos

None yet.

### Blockers/Concerns

- **Phase 5** flagged for empirical validation of nested `QEventLoop` + mouse-motion coalescing during planning
- **Phase 6** per-module lexicon audit may split into 5 sub-phases (one per solver module) during planning
- **Phase 7** requires `QImageWriter::supportedImageFormats()` probe at phase kickoff to choose GIF strategy
- **Phase 9** A/B window closure: 2026-05-06 — maintainer dricoco. Window opened 2026-05-05 (Phase 8 sign-off), closed same dev-loop session. Phase 9 EXECUTE unblocked.
- **Qt menu UX carry-forward (root-caused 2026-05-06):** maintainer flagged "Open Project not working" + "some items missing text".
  - Open Project root cause: by-design `refuseIfBlocking()` modal guard (xvue_qt_window.cpp:133-139, Phase 6.1 D-09). When Fortran is in interactive mode (post-launch), Open Project click only flashes status bar 3s. UX gap: action not visibly disabled.
  - Missing-text symptom: not reproducible on ppmail_qt AT-SPI tree (all labels present). Possibly observed on ppinit (no window — tiny launcher) or transient state.
  - Phase-9.1 candidate fix: bidirectional state binding to setEnabled(false) when blockingDepth() > 0; or QMessageBox::information instead of status-bar timeout.
- **3 P7 deferred goldens (carry-forward; root-caused 2026-05-06 by Plan 09-08):** scene01_driver.f has 2 fundamental bugs (missing XTINIT/XVINFO init → SIGSEGV; XVCHARGEFONTE arity mismatch). testa/wave + cavity2d are interactive multi-module pipelines, not headless batch. Cross-tag attempt empirically falsified the Phase 7 §9 procedure. Phase-9.1 candidate fix: repair scene01_driver.f source + write proper headless batch wrappers for testa/wave + testa/cavity2d.

## Session Continuity

Last session: 2026-05-01T08:47:27.664Z
Stopped at: Phase 7 context gathered
Resume file: .planning/phases/07-image-gif-and-postscript-export/07-CONTEXT.md
