---
phase: 05-event-bridge-blocking-reads
plan: 03
subsystem: ui
tags: [qt6, motion-coalescing, qtimer-singleshot, deferred-quit, diagnostic-env-var, fortran-abi-preserve]

requires:
  - phase: 05-event-bridge-blocking-reads-plan-02
    provides: "XvueEventBridge QObject owning a stack-local QEventLoop, installed as event filter on XvueCanvas; waitForEvent with save-restore of six filter members; eventFilter dispatches button/key synchronously; MouseMove is pass-through awaiting Plan 03"

provides:
  - "QEvent::MouseMove branch in eventFilter uses QTimer::singleShot(0, loop_, &QEventLoop::quit) deferred-quit pattern (D-04)"
  - "100-move burst coalesces to ONE waitForEvent return with tail coordinate pair — X11 XEventsQueued(QueuedAfterFlush) parity"
  - "MEFISTO_XVSOURIS_DEBUG env var (cached on first getenv) enables [xvsouris] stderr diagnostic lines carrying motion_count= field"
  - "pending_ is reset to Result{} at every waitForEvent entry so x/y cannot be stale (Pitfall 8 unit-test portion closed)"
  - "quit_pending_ reset to false and motion_count_ reset to 0 at waitForEvent entry (Pitfall 9 unit-test portion closed)"
  - "save-restore of motion_count_ across nested waitForEvent calls alongside the other six filter members"
  - "XvueCanvas::setMouseTracking(true) in ctor — required for QEvent::MouseMove delivery without a held button (xvsouris_ nbc=0 contract)"

affects:
  - 05-event-bridge-blocking-reads-plan-04  # xvsouris_ Fortran body consumes bridge->waitForEvent() motion returns with notypeevent=-2
  - 05-event-bridge-blocking-reads-plan-05  # xvsouris2_ accrochage highlight redraw triggers on each motion tail
  - 05-event-bridge-blocking-reads-plan-06  # full-plan validation, resolves Assumption A2 using the motion_count= diagnostic

tech-stack:
  added: []  # no new libraries — QTimer/QEvent already in toolchain
  patterns:
    - "QTimer::singleShot(0, QObject*, PMF) deferred-quit — enqueue a zero-delay timer-event at the tail of the event queue so all already-queued events of the target type are processed first; the timer callback quits the QEventLoop with zero added latency"
    - "Env-var feature flag cached via C++17 static-local lambda initializer so runtime lookups are a single bool read and thread-safe init is guaranteed"
    - "Save-restore of per-call filter state (now six + motion_count_ = seven members) around loop.exec() so nested waitForEvent calls do not corrupt outer-call state"
    - "Capture-by-value in QTimer::singleShot lambdas: captures of pointers/values only, never by-reference to stack locals, because the timer may survive the test function's lifetime"

key-files:
  created: []
  modified:
    - xvue/qt/src/xvue_qt_event.h
    - xvue/qt/src/xvue_qt_event.cpp
    - xvue/qt/src/xvue_qt_canvas.cpp
    - xvue/qt/tests/test_xvue_qt_event.cpp

key-decisions:
  - "setMouseTracking(true) added to XvueCanvas ctor as a real code-path enabler, not a test hack: Qt 6 drops QEvent::MouseMove for a QWidget target when tracking is off AND no button is held. Without it the xvsouris_ motion path was dead in real mesher drags. X11 parity: the legacy fenetre_mef opens with PointerMotionMask always on (xvuelc.c). Tracked as Rule 2 deviation (missing critical functionality)."
  - "Esc terminator in testMotionCoalescing fires from a *separate, later* QTimer (200 ms) rather than the same callback as the 100 moves. Reason: if Esc queues in the same batch as the 100 moves, the filter sees Esc last and returns notypeevent=0 from the *first* waitForEvent — zero motion returns counted. Delaying the Esc lets the first waitForEvent's deferred-quit timer fire at the tail of the burst, then the drain loop consumes the remaining moves and the Esc terminates."
  - "QTimer::singleShot lambdas in tests capture canvas by VALUE (not reference). Observed regression during TDD: by-reference captures over stack locals survived the test function's return and segv'd the next test when `canvas` was read from freed stack memory."
  - "Debug env var cache uses C++17 static-local lambda. Trade-off accepted: the cache is set once per process and cannot be flipped mid-run, so a dedicated unit test cannot flip and retest — the testDebugLoggingEnvVar remains QSKIP, and the runtime behavior is validated by the outer-shell grep on `motion_count=` in stderr."
  - "fprintf uses `motion_count=N` as a literal token string so the verify block's grep pattern matches unambiguously (no number-formatting surprises)."

patterns-established:
  - "Pattern: zero-delay QTimer::singleShot(0, loop, &QEventLoop::quit) as the Qt equivalent of X11's XEventsQueued(QueuedAfterFlush) — coalesces bursts of same-type events to a single return, zero added latency, no hand-rolled de-duplication"
  - "Pattern: per-call reset of all stateful filter members at waitForEvent entry; save-restore around exec() so nested calls survive cleanly (now seven members)"
  - "Pattern: gated stderr diagnostic via MEFISTO_<SUBSYSTEM>_DEBUG env var with cache-on-first-access, off by default"

requirements-completed: [EVENT-07]

duration: 50min
completed: 2026-04-18
---

# Phase 5 Plan 03: Mouse-motion coalescing + diagnostic counter Summary

**QEvent::MouseMove deferred-quit coalescing (QTimer::singleShot(0, loop_, &QEventLoop::quit)) layered on Plan 02's bridge — 100-move burst returns in ONE waitForEvent call with tail position (X11 XEventsQueued(QueuedAfterFlush) parity) plus MEFISTO_XVSOURIS_DEBUG stderr diagnostic.**

## Performance

- **Duration:** 50 min
- **Started:** 2026-04-18T09:57:00Z
- **Completed:** 2026-04-18T10:47:19Z
- **Tasks:** 1
- **Files modified:** 4 (0 created, 4 modified)

## Accomplishments

- `QEvent::MouseMove` branch in `XvueEventBridge::eventFilter` replaces Plan 02's pass-through with the full D-04 deferred-quit pattern: stash `pending_.notypeevent=-2`, `pending_.x/y = me->pos().x()/y()`, increment `motion_count_`, and arm a zero-delay `QTimer::singleShot(0, loop_, &QEventLoop::quit)` if `!quit_pending_`. The timer enqueues loop.quit() at the tail of the event queue so any queued motion events are dispatched first (each overwriting pending_), and loop.exec() returns with the *last* coordinate pair in the burst.
- `XvueEventBridge::debug_logging_enabled()` static accessor caches the `getenv("MEFISTO_XVSOURIS_DEBUG")` result in a static-local bool on first call. Non-empty, non-"0" values enable logging; default-off.
- `waitForEvent()` emits one stderr line per return when the debug flag is set: `[xvsouris] mode=N notypeevent=N nbc=N x=N y=N motion_count=N depth=N`. The literal token `motion_count=` makes the outer-shell grep unambiguous.
- `motion_count_` reset to 0 at each `waitForEvent` entry and save-restored across nested calls (now seven filter members in the save-restore set).
- `pending_` was already reset to `Result{}` at waitForEvent entry in Plan 02 — the unit-test portion of Pitfall 8 (no stale x/y across calls) is verified by the new `testMotionFreshPerCall`.
- `quit_pending_` was already reset to false at waitForEvent entry in Plan 02 — the unit-test portion of Pitfall 9 (timer can re-arm on the next call) is verified by the new `testQuitPendingResetAcrossCalls`.
- `XvueCanvas::setMouseTracking(true)` added to the canvas constructor — required for QEvent::MouseMove delivery when no mouse button is held (xvsouris_ nbc=0 motion contract). Parity with X11 fenetre_mef opening with PointerMotionMask always on.
- Four new real motion tests PASS under xvfb-run (no QSKIPs left for EVENT-07):
  - `testXvsourisMotion` — single move → `(notypeevent=-2, nbc=0, x=42, y=73)`.
  - `testMotionCoalescing` — 100 rapid moves drain to **one** waitForEvent return at `(109, 119)` with `motion_count=100`; drain loop asserts `returns <= 20` (actual: 1).
  - `testMotionFreshPerCall` — two sequential calls with different coordinates; pending_ never leaks stale `(x, y)` across invocations.
  - `testQuitPendingResetAcrossCalls` — second pure-motion call succeeds without a fallback Esc; would time out (not deadlock silently) if `quit_pending_` were not reset.
- With `MEFISTO_XVSOURIS_DEBUG=1`: stderr contains **19** `motion_count=` lines across the suite. Without the env var: **0** such lines (production stderr stays clean).
- Full test suite: **21 passed, 0 failed, 7 skipped** (the seven skips are Plan 04/05/06 placeholders — xvpause_, xvsouris2_ accrochage, deplsouris_, initaccrochage_, debug env-var runtime flip).
- `bin/cbl_tout_qt` (Qt backend) and `bin/cbl_tout` (X11 baseline) both exit 0. SHELL-03 `verify_no_exec` guard still passes — zero matches for `QApplication::exec` or `qApp->exec()` anywhere under `xvue/qt/`.
- Acceptance-criteria grep counts all exceed threshold:
  - `QTimer::singleShot` in event.cpp: **2** (≥1 required)
  - `motion_count_` in event.cpp: **6** (≥3 required)
  - `MEFISTO_XVSOURIS_DEBUG` in event.cpp: **3** (≥1 required)
  - `pending_.notypeevent = -2` in event.cpp: **1** (≥1 required)

## Task Commits

Each task committed atomically with `--no-verify` (parallel-executor guidance — orchestrator runs hooks once post-merge):

1. **Task 1: MouseMove deferred-quit coalescing + MEFISTO_XVSOURIS_DEBUG diagnostic counter** — `03f11bb` (feat)

## Files Created/Modified

- `xvue/qt/src/xvue_qt_event.h` (MODIFIED) — added `int motion_count_ = 0;` field and `static bool debug_logging_enabled();` accessor declaration.
- `xvue/qt/src/xvue_qt_event.cpp` (MODIFIED, +40 / -3 lines) — added `#include <QTimer>`, `<cstdio>`, `<cstdlib>`; implemented `debug_logging_enabled()` with static-local cache; extended save-restore set to 7 members; reset `motion_count_ = 0` at waitForEvent entry; fprintf stderr diagnostic gated by the cache; replaced MouseMove pass-through with the full coalescing pattern (pending_ stash, ++motion_count_, QTimer::singleShot(0, loop_, &QEventLoop::quit) if !quit_pending_, return true).
- `xvue/qt/src/xvue_qt_canvas.cpp` (MODIFIED, +9 lines) — `setMouseTracking(true)` in ctor with comment explaining the Qt 6 "drop MouseMove without held button" behavior and X11 PointerMotionMask parity.
- `xvue/qt/tests/test_xvue_qt_event.cpp` (MODIFIED) — added `postMove()` static helper, replaced four QSKIPs with real bodies: testXvsourisMotion, testMotionCoalescing, testMotionFreshPerCall, testQuitPendingResetAcrossCalls; added testDebugLoggingEnvVar placeholder QSKIP explaining the cache-on-first-access limitation.

## Decisions Made

1. **setMouseTracking(true) in the canvas constructor (Rule 2 auto-fix).** The plan spec referenced Pitfall 4 (StrongFocus) but did not mention mouse tracking. Observed failure: testXvsourisMotion timed out at 300s because Qt's `QApplicationPrivate::notify_helper` drops `QEvent::MouseMove` for a QWidget target when `mouseTracking` is false AND no button is held. The xvsouris_ contract is motion-without-button (nbc=0), so the production code path was actually broken, not just the test. Fix: `setMouseTracking(true)` in the canvas ctor, mirroring X11's `PointerMotionMask`. Tracked as Rule 2 deviation (missing critical functionality — this is correctness, not a test-only fix).
2. **Esc terminator in testMotionCoalescing uses a SEPARATE, later timer (200 ms).** First attempt queued all 100 moves + Esc in one callback. Result: the filter dispatched the 100 moves (each arming the deferred-quit once, then updating pending_ with fresh x/y), then saw the Esc which reset pending_ to `(notypeevent=0, nbc=27)` and called loop.quit() directly — no motion return counted. Fix: schedule Esc at 200 ms so the first waitForEvent's deferred-quit fires on the burst's tail *before* any Esc arrives; subsequent drain waitForEvent calls eventually consume the Esc and break.
3. **Capture canvas by value in QTimer::singleShot lambdas.** Observed regression during TDD: by-reference captures to the test function's local `canvas*` variable survived the test function's return. When the timer fired later (in the next test), it dereferenced freed stack memory and segv'd in `QCoreApplicationPrivate::lockThreadPostEventList`. All motion-test lambdas now capture pointers/values only.
4. **Debug env var uses C++17 static-local cache.** Thread-safe initialization is guaranteed by the standard (and we are main-thread-only anyway per SHELL-07). The cache is fixed for the process — a dedicated unit test cannot flip and retest — so `testDebugLoggingEnvVar` remains QSKIP'd with a pointer to the acceptance grep that validates the runtime behavior from outside the test binary.
5. **Diagnostic format string includes the literal token `motion_count=`** so the acceptance grep pattern matches unambiguously regardless of integer formatting. Matches both the plan's verify block and Plan 06's future empirical-measurement harness.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing critical functionality] XvueCanvas was missing setMouseTracking(true)**
- **Found during:** Task 1 (first test run under xvfb — testXvsourisMotion timed out at 300s)
- **Issue:** Qt 6's `QApplicationPrivate::notify_helper` drops `QEvent::MouseMove` for a QWidget target when mouseTracking is false AND no mouse button is held. The xvsouris_ motion contract reports events with `nbc=0` (no button), so in production the motion path was actually broken for real interactive drags without a held button (rubber-band preview, hover snap feedback). The X11 fenetre_mef opens with `PointerMotionMask` always on — the Qt equivalent is setMouseTracking(true).
- **Fix:** Added `setMouseTracking(true)` to XvueCanvas ctor immediately after the `setFocusPolicy(Qt::StrongFocus)` call, with a comment documenting the Qt delivery rule and the X11 parity.
- **Files modified:** xvue/qt/src/xvue_qt_canvas.cpp
- **Verification:** testXvsourisMotion PASSes in <10ms after the fix; full coalescing test (100 moves → 1 return) also PASSes.
- **Committed in:** 03f11bb (Task 1 commit)

**2. [Rule 1 - Bug] Test singleShot lambdas captured stack pointers by reference**
- **Found during:** Task 1 (after the setMouseTracking fix unblocked testXvsourisMotion, the test suite exit-134'd with SIGSEGV in testMotionFreshPerCall)
- **Issue:** `QTimer::singleShot(200, [&]{ postKey(canvas, ...); });` inside testMotionCoalescing's drain loop captured `canvas` (a stack-local `QWidget*`) by reference. When testMotionCoalescing returned, its stack frame was reclaimed. The next test (testMotionFreshPerCall) ran — and then the leaked singleShot from testMotionCoalescing fired, reading `canvas` from freed stack memory and crashing in `QCoreApplicationPrivate::lockThreadPostEventList`.
- **Fix:** All motion-test singleShot lambdas now capture `canvas` by value (`[canvas]{...}`). Removed the 20 per-iteration fallback Esc timers in testMotionCoalescing — replaced by a single deferred Esc timer at 200 ms (captured by value).
- **Files modified:** xvue/qt/tests/test_xvue_qt_event.cpp
- **Verification:** Full test suite exits 0 consistently (3/3 runs). No more stack-frame use-after-free.
- **Committed in:** 03f11bb (Task 1 commit)

**3. [Rule 1 - Bug] testMotionCoalescing initially drained zero motion returns**
- **Found during:** Task 1 (after the first two fixes landed, testMotionCoalescing passed all invocations individually but failed in the full-suite run with `returns >= 1 returned FALSE`)
- **Issue:** When both the 100 moves and the Esc terminator are posted from the same callback, Qt's event dispatcher runs the 100 moves (each stashing pending_) and then the Esc, which synchronously overwrites pending_ to `(notypeevent=0, nbc=27)` and calls loop.quit() directly — *before* the deferred-quit timer fires on the motion tail. First waitForEvent returns Esc, not motion, and `returns` stays at 0.
- **Fix:** Split the posting into two timers: moves at t=10ms, Esc at t=200ms. The 190ms gap gives the first waitForEvent's deferred-quit time to fire on the motion tail and return with `(notypeevent=-2, motion_count=100)`; subsequent drain waitForEvent calls eventually see the Esc and break.
- **Files modified:** xvue/qt/tests/test_xvue_qt_event.cpp
- **Verification:** `testMotionCoalescing` consistently returns exactly 1 motion burst with `x=109, y=119, motion_count=100` (observed via MEFISTO_XVSOURIS_DEBUG=1 stderr).
- **Committed in:** 03f11bb (Task 1 commit)

---

**Total deviations:** 3 auto-fixed (1 missing-functionality, 2 bugs in the test design). All three stayed within the files Plan 03 already touches. No scope creep.

## Issues Encountered

- **xvfb MouseMove delivery without setMouseTracking:** the single biggest timesink of this plan. The plan's unit-test stubs assumed Qt would deliver synthetic MouseMove events through `QCoreApplication::postEvent` regardless of widget state. In reality Qt drops them at `notify_helper` when mouseTracking is off and no button is held. The fix is a one-liner in the canvas ctor but it took a 300s test-timeout to diagnose.
- **Test-function singleShot lifetime:** QTimer singleShots survive the test-function return. Always capture-by-value.
- **Event-dispatch order: same-callback posts are atomic from the filter's perspective.** A burst of 100 MouseMoves + 1 KeyPress in the same postEvent storm is dispatched as one run of 101 filter invocations, not 101 separate runs. This has design implications for coalescing tests — the terminator must be posted in a *later* timer to let the deferred-quit fire.

## User Setup Required

None — MEFISTO_XVSOURIS_DEBUG is documented in the commit message and the SUMMARY; developers wanting to see the diagnostic simply set the env var before running the test binary. No config files, no service setup.

## Next Phase Readiness

- **Plan 04 (Fortran ABI bodies — xvsouris_, xvpause_, deplsouris_):** Ready. `bridge->waitForEvent()` now returns correct motion codes alongside button/key codes. `xvsouris_` can map any (notypeevent, nbc, x, y) Result straight through to the four Fortran output integers.
- **Plan 05 (xvsouris2_ accrochage Strategy B):** Ready. `waitForEvent(WaitMode::Souris2, items, pmin0)` already plumbs the items/pmin0 arrays through to the filter via the saved members. Plan 05's handler code can run its nearest-item computation inside the MouseMove and ButtonPress branches and blit the accrochage highlight via the `mempxaccro_` + `accroche_undo_tile_` pair on XvueState. Motion coalescing is friendly to the blit pattern — the highlight redraws only once per burst tail, not once per raw move.
- **Plan 06 (empirical validation + Assumption A2 resolution):** Ready. The `motion_count=` stderr diagnostic is the exact instrument Plan 06 needs. Under `pp/ppmail_qt testa/pan2d` with MEFISTO_XVSOURIS_DEBUG=1, a developer can measure how many motion events Qt delivers to the filter during a 2-second rapid drag and compare to the X11 baseline.
- **No blockers or concerns.**

## Known Stubs

- `testDebugLoggingEnvVar` is a documented QSKIP'd test. The env-var cache cannot be flipped at runtime, so an in-process unit test cannot validate both branches. The runtime behavior is instead validated by the outer-shell grep on `motion_count=` from stderr — documented in the acceptance_criteria block of 05-03-PLAN.md and in the commit message. Not a functional stub; the debug counter itself is fully wired and has been observed emitting 19 lines in a full-suite run with the flag on, 0 with it off.

## Self-Check: PASSED

- `xvue/qt/src/xvue_qt_event.h`: contains `int motion_count_` and `static bool debug_logging_enabled()`
- `xvue/qt/src/xvue_qt_event.cpp`: contains `QTimer::singleShot`, `motion_count_`, `MEFISTO_XVSOURIS_DEBUG`, `pending_.notypeevent = -2`
- `xvue/qt/src/xvue_qt_canvas.cpp`: contains `setMouseTracking(true)` in ctor
- `xvue/qt/tests/test_xvue_qt_event.cpp`: contains testXvsourisMotion real body, testMotionCoalescing, testMotionFreshPerCall, testQuitPendingResetAcrossCalls
- Commit 03f11bb: FOUND
- `bin/cbl_tout_qt`: exits 0
- `bin/cbl_tout`: exits 0
- `xvfb-run -a xvue/qt/build/tests/xvue_qt_event_tests`: 21 passed, 0 failed, 7 skipped
- `MEFISTO_XVSOURIS_DEBUG=1` stderr `motion_count=` lines: 19
- Unset stderr `motion_count=` lines: 0
- `grep -c QTimer::singleShot xvue_qt_event.cpp`: 2 (>= 1 required)
- `grep -c motion_count_ xvue_qt_event.cpp`: 6 (>= 3 required)
- `grep -c MEFISTO_XVSOURIS_DEBUG xvue_qt_event.cpp`: 3 (>= 1 required)
- `grep -c 'pending_.notypeevent = -2' xvue_qt_event.cpp`: 1 (>= 1 required)
- `grep -r QApplication::exec xvue/qt/src/ xvue/qt/include/`: 0 matches (SHELL-03 clean)

---
*Phase: 05-event-bridge-blocking-reads*
*Completed: 2026-04-18*
