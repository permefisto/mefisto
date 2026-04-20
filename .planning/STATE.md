---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: milestone
status: executing
stopped_at: Phase 6.1 context gathered
last_updated: "2026-04-20T15:40:33.052Z"
last_activity: 2026-04-20
progress:
  total_phases: 17
  completed_phases: 9
  total_plans: 34
  completed_plans: 35
  percent: 100
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-04-10)

**Core value:** Every MEFISTO workflow that works today through X11 keeps working through the new Qt 6 interface, with Fortran solver code unchanged.
**Current focus:** Phase 06.0 — shared-shell-menu-bridge-dialogs

## Current Position

Phase: 6.1
Plan: Not started
Status: Executing Phase 06.0
Last activity: 2026-04-20

Progress: [██████████] 100% (Phase 03.1)

## Performance Metrics

**Velocity:**

- Total plans completed: 29
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

**Recent Trend:**

- Last 5 plans: —
- Trend: —

*Updated after each plan completion*
| Phase 01-window-shell-xvueapp-xvuewindow-xvuecanvas P01 | 19min | 3 tasks | 10 files |
| Phase 01 P02 | 15min | 2 tasks | 1 files |
| Phase 01 P03 | ~60min | 3 tasks | 5 files |

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
- **Phase 9** is gated on the one-release-cycle A/B window closing — process gate, not date-driven

## Session Continuity

Last session: 2026-04-20T15:40:33.046Z
Stopped at: Phase 6.1 context gathered
Resume file: .planning/phases/06.1-mesher-mail-menu-wiring/06.1-CONTEXT.md
