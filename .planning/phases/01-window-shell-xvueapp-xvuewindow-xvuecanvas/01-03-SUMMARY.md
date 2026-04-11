---
phase: 01-window-shell-xvueapp-xvuewindow-xvuecanvas
plan: 03
subsystem: xvue-qt
tags: [fortran, qt, qt6, shell, shell-01, shell-02, shell-06, validation, visual, hidpi]

requires:
  - phase: 01-window-shell-xvueapp-xvuewindow-xvuecanvas
    plan: 02
    provides: "Real xvinitgraphique_/xvfermer_/xvpxecran_/xvmmecran_/xvfond_/xvinfo_ bodies in libxvueqt.a; XVUE_QT_ASSERT_MAIN_THREAD macro retrofitted into all 57 stubs"
provides:
  - "prpr/xvtest0.f — Fortran 77 Phase 1 lifecycle test driver exercising open/close/reopen cycle"
  - "bin/cbxvtest0_qt — standalone build script producing pp/ppxvtest0_qt against libxvueqt.a"
  - "pp/ppxvtest0_qt — compiled Phase 1 test executable (build artifact)"
  - "Empirical proof of SHELL-01 (window opens), SHELL-02 (reopen without singleton crash), SHELL-06 (HiDPI scaling)"
  - "Three in-plan gap closures: FIX-1 (leak QApplication at atexit), FIX-2 (bounded-loop exposure pump), FIX-3 (xvfermer_ drains deleteLater queue)"
affects:
  - "02-drawing-primitives — xvtest0 remains in tree as permanent shell sanity driver; Phase 2+ can extend it"
  - "01-CONTEXT.md D-01/D-06/D-08 — revised with empirical findings; downstream phase planners must read revised decisions"
  - "all future phases — QApplication leak-at-atexit pattern is now the established precedent"

tech-stack:
  added: []
  patterns:
    - "Fortran test driver in prpr/ calling Phase 1 SHELL entry points via trailing-underscore ABI"
    - "Bounded exposure pump: QElapsedTimer + isExposed() loop replaces single processEvents call"
    - "QApplication deliberate leak at atexit: qapp_.release() not qapp_.reset()"
    - "xvfermer_ drains deleteLater queue with processEvents after window_.reset()"

key-files:
  created:
    - "prpr/xvtest0.f"
    - "bin/cbxvtest0_qt"
    - ".planning/phases/01-window-shell-xvueapp-xvuewindow-xvuecanvas/01-03-SUMMARY.md"
  modified:
    - "xvue/qt/src/xvue_qt_app.cpp"
    - "xvue/qt/src/xvue_qt_api.cpp"
    - "prpr/xvtest0.f (SLEEP hold added for human-visual check)"
    - ".planning/phases/01-window-shell-xvueapp-xvuewindow-xvuecanvas/01-CONTEXT.md"

key-decisions:
  - "QApplication is deliberately leaked at process exit — never destroyed in atexit (revised D-08). See debug session for backtrace evidence."
  - "xvinitgraphique_ pumps processEvents in bounded loop until isExposed(), not once (revised D-01). X11 MapRequest/ConfigureNotify/Expose requires multiple event-loop trips."
  - "xvfermer_ drains deleteLater queue after window_.reset() (D-06 addendum) to prevent stale events at atexit."
  - "prpr/xvtest0.f adds SLEEP(1) hold between open and close so each window is humanly observable (visual validation requirement)."

patterns-established:
  - "Embed Qt in Fortran main: construct QApplication once via call_once, leak at exit (never destroy)"
  - "Exposure pump idiom: QElapsedTimer + while(elapsed < 2000) { processEvents(ExcludeUserInputEvents, 20); if (isExposed()) break; }"
  - "xvfermer_ teardown sequence: window_.reset() -> processEvents(ExcludeUserInputEvents) -> return"

requirements-completed:
  - SHELL-01
  - SHELL-02
  - SHELL-06

duration: ~60min (Tasks 1+2 ~20min; debug investigation and fixes ~40min)
completed: 2026-04-11
---

# Phase 1 Plan 03: Fortran Test Driver + Full-Suite Validation Summary

**Fortran lifecycle driver prpr/xvtest0.f exercises open/close/reopen through the Qt 6 ABI; three in-plan gap closures fix bounded-loop window exposure, xvfermer_ deleteLater drain, and the QApplication deliberate-leak-at-atexit pattern; SHELL-01/02/06 empirically confirmed by human visual check.**

## Performance

- **Duration:** ~60 min
- **Started:** 2026-04-11 (after 01-02 checkpoint)
- **Completed:** 2026-04-11
- **Tasks:** 3 (Tasks 1+2 autonomous; Task 3 human-verify checkpoint — approved)
- **Files modified:** 5 (2 created, 3 modified for gap closures)

## Accomplishments

- Created `prpr/xvtest0.f` as a Fortran 77 fixed-form lifecycle driver: calls XVINITGRAPHIQUE/XVFERMER twice in sequence, with SLEEP(1) between each open/close to hold the window visible for human confirmation of SHELL-01/02/06.
- Created `bin/cbxvtest0_qt` cloned from `bin/cbmail_qt` with six substitutions (source file, exe name, dropped `mail/lib`, updated labels); produces `pp/ppxvtest0_qt` at 755.
- **In-plan gap closures (three commits):**
  - FIX-1 (`xvue_qt_app.cpp`): `teardown_atexit()` now calls `qapp_.release()` not `qapp_.reset()` — QApplication is leaked deliberately to prevent SIGSEGV in `__run_exit_handlers` caused by destroying QApplication alongside libgfortran's atexit chain on Linux/Qt 6.
  - FIX-2 (`xvue_qt_api.cpp`): `xvinitgraphique_` uses a `QElapsedTimer`-bounded loop pumping `processEvents(ExcludeUserInputEvents, 20)` until `win->windowHandle()->isExposed()` — necessary because X11 requires multiple event-loop trips (MapRequest → ConfigureNotify → Expose) before the window is visible.
  - FIX-3 (`xvue_qt_api.cpp`): `xvfermer_` calls `processEvents(ExcludeUserInputEvents)` after `window_.reset()` to drain the `deleteLater()` queue queued by widget teardown.
- **Human visual confirmation:** User ran `pp/ppxvtest0_qt` (two 800×600 black "MEFISTO" windows visible, exit 0) and `QT_SCALE_FACTOR=2 pp/ppxvtest0_qt` (visibly ~2× larger windows, exit 0, no core dump). SHELL-01, SHELL-02, SHELL-06 empirically validated.
- `verify_abi`: 57/57 symbols. `verify_no_exec`: OK. `cbl_tout` and `cbl_tout_qt`: both green.

## Task Commits

1. **Task 1: Create prpr/xvtest0.f Fortran lifecycle driver** — `673de2b` (feat)
2. **Task 2: Create bin/cbxvtest0_qt build script** — `4f23477` (feat)
3. **Task 3: Full-suite verification (human-verify checkpoint)** — no code commit (verification only)

**In-plan gap closures (debug session — between Task 2 commit and Task 3 approval):**
- `36f4fb6` — fix(xvue-qt): leak QApplication at exit to prevent teardown SIGSEGV (FIX-1)
- `a7b1c2c` — fix(xvue-qt): pump events until window exposed in xvinitgraphique_ (FIX-2 + FIX-3)
- `006ff3d` — test(phase-1): hold each xvtest0 open cycle with SLEEP(1) for visual check (FIX-4)
- `43551af` — docs(debug): capture phase-01 xvtest0 teardown-segfault investigation

## Files Created/Modified

### Created

- `prpr/xvtest0.f` — Fortran 77 fixed-form Phase 1 lifecycle driver (40 lines). Two XVINITGRAPHIQUE/XVFERMER cycles with SLEEP(1) holds. Column-7+ layout following project normes. Follows `prpr/xvtest{1..4}.f` convention.
- `bin/cbxvtest0_qt` — 0755 shell script, clone of `bin/cbmail_qt` with six substitutions: header comment, source file (`xvtest0.f`), exe name (`ppxvtest0_qt` ×3), label text (`PHASE 1 TEST`), and `mail/lib` dropped from linker line.

### Modified

- `xvue/qt/src/xvue_qt_app.cpp` — `teardown_atexit()` body revised: pumps processEvents before and after `window_.reset()`, then calls `qapp_.release()` instead of `qapp_.reset()`. QApplication is now a documented deliberate leak.
- `xvue/qt/src/xvue_qt_api.cpp` — `xvinitgraphique_` body revised with bounded exposure loop. `xvfermer_` body revised to drain deferred-delete queue after `window_.reset()`. Added `#include <QElapsedTimer>` and `#include <QWindow>`.
- `prpr/xvtest0.f` — `CALL SLEEP(1)` added after each `CALL XVINITGRAPHIQUE` so windows are visible for 1 second during human validation (FIX-4).

## Decisions Made

- **QApplication leak is the correct idiom for Fortran-hosted Qt on Linux.** The research plan prescribed `qapp_.reset()` at atexit (D-08), which matches the Qt embedding documentation at face value. In practice, Qt 6's internal static teardown on Linux interleaves destructively with libgfortran's `__run_exit_handlers`. The backtrace (frames 7–10 in program text inside `teardown_atexit` reaching into Qt internal state) confirms RC-2. The fix — `release()` not `reset()` — is the pattern used by PyQt, PySide, and Qt-plugin host applications. D-08 is revised in place in `01-CONTEXT.md`.
- **Single processEvents is insufficient for X11 window realization.** The X11 compositor requires: the kernel to dispatch the client's `MapRequest`, the server to generate `ConfigureNotify`, and then `Expose` to arrive and invoke `paintEvent`. One round of `processEvents` dispatches only the first of these steps. The bounded loop (≤2 s, break on `isExposed()`) is the minimal correct fix. D-01 is revised in place.
- **xvtest0.f keeps SLEEP(1) permanently.** The delay is acceptable for a human-visual validation driver; real interactive modules (Phase 2+) use `xvvoir_`'s pumped wait, not `SLEEP`. The file stays in the tree as a permanent shell sanity driver per the plan's deferred note.

## Deviations from Plan

### Auto-fixed Issues (in-plan gap closures)

**1. [Rule 1 - Bug] Fixed SIGSEGV in teardown_atexit: QApplication destruction race with libgfortran atexit chain**
- **Found during:** Task 3 (first run of pp/ppxvtest0_qt)
- **Issue:** `teardown_atexit()` called `qapp_.reset()`, destroying `QApplication` inside an atexit handler interleaved with libgfortran's own `__run_exit_handlers`. Qt's internal static teardown destructors crashed on partially-unwound state. SIGSEGV at frames 7–10 in program text.
- **Fix:** Changed `qapp_.reset()` to `qapp_.release()` — QApplication is leaked (documented deliberate leak). Added processEvents pumps around `window_.reset()` to drain pending deferred-delete events first.
- **Files modified:** `xvue/qt/src/xvue_qt_app.cpp`
- **Commit:** `36f4fb6`

**2. [Rule 1 - Bug] Fixed no-visible-window: single processEvents insufficient for X11 window realization**
- **Found during:** Task 3 (first run — user observed no window on screen)
- **Issue:** `xvinitgraphique_` called `processEvents` exactly once (per D-01 original). X11 compositor requires multiple event-loop trips for MapRequest → ConfigureNotify → Expose before `isExposed()` returns true. Window was mapped but not yet exposed by the time Fortran's CALL XVFERMER ran.
- **Fix:** Replaced single `processEvents` call with a `QElapsedTimer`-bounded loop (`while elapsed < 2000 ms: processEvents(ExcludeUserInputEvents, 20); if isExposed() break`).
- **Files modified:** `xvue/qt/src/xvue_qt_api.cpp`
- **Commit:** `a7b1c2c`

**3. [Rule 1 - Bug] Fixed latent stale-event accumulation in xvfermer_**
- **Found during:** Debug investigation of RC-2 (secondary contributing factor)
- **Issue:** `xvfermer_` did not drain the `deleteLater()` queue after `window_.reset()`. Stale `DeferredDelete` events accumulated in Qt's event queue and were only flushed when QApplication destructed — the worst possible time.
- **Fix:** Added `QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents)` after `window_.reset()` in `xvfermer_`.
- **Files modified:** `xvue/qt/src/xvue_qt_api.cpp`
- **Commit:** `a7b1c2c` (same commit as Fix 2)

**4. [Rule 2 - Missing Critical] Added SLEEP(1) hold to xvtest0.f for human-visual validation**
- **Found during:** Debug investigation of RC-1 (window imperceptibly brief)
- **Issue:** `xvtest0.f` called XVFERMER immediately after XVINITGRAPHIQUE with no hold. Even with the correct exposure pump, the window would be visible for microseconds — insufficient for human confirmation of SHELL-01/06.
- **Fix:** Added `CALL SLEEP(1)` after each `CALL XVINITGRAPHIQUE`.
- **Files modified:** `prpr/xvtest0.f`
- **Commit:** `006ff3d`

---

**Total deviations:** 4 auto-fixed (3 Rule 1 bugs in plan D-01/D-08 design; 1 Rule 2 missing visual hold)
**Design-level note:** Deviations 1, 2, and 3 are failures of the Phase 1 _design_ (D-01/D-08), not of the Phase 1 _executor_. The implementation was faithful to the decisions as written; the decisions themselves were underspecified or incorrect. D-01 and D-08 have been revised in `01-CONTEXT.md` with full history preserved.
**Impact on plan:** All four fixes are required for the plan's success criteria (two visible windows, exit 0, no core dump). No scope creep.

## Empirical Validation Results

| Check | Result |
|-------|--------|
| `pp/ppxvtest0_qt` exit code | 0 |
| Two visible 800×600 "MEFISTO" windows (SHELL-01 + SHELL-02) | Human confirmed |
| No "QApplication: there can only be one" in stderr | Confirmed |
| No core dump | Confirmed |
| `QT_SCALE_FACTOR=2 pp/ppxvtest0_qt` exit code | 0 |
| Visibly ~2× larger windows on non-HiDPI display (SHELL-06) | Human confirmed |
| `verify_abi` symbol count | 57/57 |
| `verify_no_exec` | OK |
| `cbl_tout` legacy X11 build | Green |
| `cbl_tout_qt` Qt build | Green |

## Issues Encountered

- **No visible window on first run + SIGSEGV** — Two coupled design bugs discovered in the first `pp/ppxvtest0_qt` run. Both were root-caused within the same debug session, fixed, and the fixes verified before presenting to the human verifier. See `.planning/debug/resolved/phase-01-xvtest0-teardown-segfault.md` for the full investigation.
- **MEFISTO environment variable unset in fresh shell** (recurring from Plans 01-01 and 01-02) — resolved by inlining `export MEFISTO=...` before build commands.

## Known Stubs

No new stubs introduced by this plan. All 51 warn-once stubs are the same Phase 0 stubs inherited from Plan 01-02, scheduled for Phases 2–8.

## User Setup Required

None — no external service configuration required. The test driver requires a running X11/Wayland display for the visual run; that is a runtime requirement, already satisfied by the user's desktop session.

## Next Phase Readiness

- **Phase 1 is complete.** All 7 SHELL requirements (SHELL-01..07) are validated:
  - SHELL-01: open window confirmed visually (xvtest0 first cycle)
  - SHELL-02: reopen confirmed visually (xvtest0 second cycle, no singleton assertion)
  - SHELL-03: verify_no_exec guard armed and tested in Plan 01-01
  - SHELL-04: xvpxecran_/xvmmecran_ compiled and linked (Plan 01-02)
  - SHELL-05: xvfond_ stores background color (Plan 01-02)
  - SHELL-06: QT_SCALE_FACTOR=2 run confirmed visually
  - SHELL-07: XVUE_QT_ASSERT_MAIN_THREAD in all 57 stubs (Plan 01-02)
- **Phase 2 (drawing primitives)** is unblocked. Key changes from this plan that Phase 2 inherits:
  - Exposure pump pattern: Phase 2's `xvvoir_` implementation should use the same `QElapsedTimer`-bounded `processEvents` loop.
  - `xvfermer_` drain pattern: established; Phase 2 must not bypass it.
  - QApplication leak: permanent — no phase may add `qapp_.reset()` or `qapp_.quit()` at atexit.
- **Legacy X11 build untouched**: `bin/cbl_tout` still green; `xvue/xvuelc.c` is clean.
- **No threat flags** — Plan 01-03 did not introduce new network endpoints, auth paths, or trust-boundary surface.

## Self-Check: PASSED

Verified:
- `prpr/xvtest0.f` exists: FOUND
- `bin/cbxvtest0_qt` exists: FOUND
- `xvue/qt/src/xvue_qt_app.cpp` exists: FOUND
- `xvue/qt/src/xvue_qt_api.cpp` exists: FOUND
- Task commits `673de2b` and `4f23477` in git log: FOUND
- Gap-closure commits `36f4fb6`, `a7b1c2c`, `006ff3d`, `43551af` in git log: FOUND
- Human verification result: confirmed fixed (both runs, exit 0, visible windows)

---
*Phase: 01-window-shell-xvueapp-xvuewindow-xvuecanvas*
*Completed: 2026-04-11*
