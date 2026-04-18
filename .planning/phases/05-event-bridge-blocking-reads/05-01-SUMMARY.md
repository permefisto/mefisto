---
phase: 05-event-bridge-blocking-reads
plan: 01
subsystem: testing
tags: [qt6, qtest, cmake, xvfb, event-bridge, wave-0, scaffolding]

# Dependency graph
requires:
  - phase: 01-window-shell-xvueapp-xvuewindow-xvuecanvas
    provides: XvueApp/XvueWindow/XvueCanvas/XvueState quartet + SHELL-03 verify_no_exec guard
  - phase: 02-drawing-primitives-backing-pixmap
    provides: XvueState::saved_canvas_ ownership pattern (model for mempxaccro_/accroche_undo_tile_)
  - phase: 04-pixmap-save-restore-double-buffering
    provides: Phase 4 save/restore pattern reused for accrochage tile
provides:
  - "QTest target xvue_qt_event_tests buildable via bin/cbl_tout_qt"
  - "XvueApp::blockingDepth() static accessor returning 0 in a fresh process"
  - "Qt::AA_CompressHighFrequencyEvents explicitly set before QApplication ctor (D-05)"
  - "XvueState::mempxaccro_ and XvueState::accroche_undo_tile_ QPixmap* slots (nullptr at ctor, deleted in dtor)"
  - "XvueCanvas Qt::StrongFocus policy (Pitfall 4 closed)"
  - "XvueCanvas::resizeEvent invalidates accroche_undo_tile_ (Pitfall 10 partially closed)"
  - "BlockingDepthGuard friend decl hook for Plan 02 RAII guard"
  - "QSKIP-stubbed test matrix for EVENT-01..EVENT-08 (Plans 02-06 drop in real bodies)"
affects: [05-02, 05-03, 05-04, 05-05, 05-06, phase-06-modal-dialogs]

# Tech tracking
tech-stack:
  added: [Qt6::Test (already present via qt6-base-dev), QTEST custom main, xvfb-run harness]
  patterns:
    - "Test gating via CMake option XVUE_QT_BUILD_TESTS (OFF by default, ON in cbl_tout_qt)"
    - "Custom main() wrapping QTest::qExec so XvueApp::ensure() owns QApplication construction"
    - "Static zero-init of blockingDepth_ + friend-scoped mutation"
    - "XvueState accrochage pixmaps mirror Phase 4 saved_canvas_ raw-pointer + manual-lifecycle pattern"
    - "Canvas resize-event invalidation of backing-coupled caches (accroche_undo_tile_)"

key-files:
  created:
    - xvue/qt/tests/CMakeLists.txt
    - xvue/qt/tests/test_xvue_qt_event.cpp
    - xvue/qt/tests/test_helpers.h
    - xvue/qt/tests/test_helpers.cpp
  modified:
    - xvue/qt/CMakeLists.txt
    - bin/cbl_tout_qt
    - xvue/qt/src/xvue_qt_app.h
    - xvue/qt/src/xvue_qt_app.cpp
    - xvue/qt/src/xvue_qt_state.h
    - xvue/qt/src/xvue_qt_state.cpp
    - xvue/qt/src/xvue_qt_canvas.cpp

key-decisions:
  - "Test binary uses a custom main() instead of QTEST_MAIN so XvueApp::ensure() owns QApplication construction and D-05's 'set AA_CompressHighFrequencyEvents BEFORE QApplication ctor' invariant actually holds. QTEST_MAIN would construct a second QApplication and trip Qt's single-instance assertion."
  - "Test subdirectory is opt-in via -DXVUE_QT_BUILD_TESTS=ON. Default OFF keeps the normal cmake configure clean; cbl_tout_qt passes the flag."
  - "BlockingDepthGuard is forward-friend-declared in xvue_qt_app.h so only Plan 02's RAII guard (living in xvue_qt_event.h) can touch the counter — encapsulation with zero runtime cost."
  - "mempxaccro_ is resolution-independent and NOT invalidated on resize; accroche_undo_tile_ IS invalidated because it references live backing pixels."

patterns-established:
  - "QTest targets link xvueqt static lib + Qt6::Test + Qt6::Widgets; include dirs reach into ../src and ../include"
  - "bin/cbl_tout_qt builds the test target as an explicit follow-on cmake --build target, with existence check"
  - "New static members on XvueApp added in the public section with friend-scoped mutator access"

requirements-completed: [EVENT-01, EVENT-06, EVENT-08]

# Metrics
duration: ~35min
completed: 2026-04-14
---

# Phase 5 Plan 01: Wave 0 Scaffolding — QTest infra + blockingDepth + accrochage state slots

**xvue_qt_event_tests QTest binary builds under bin/cbl_tout_qt, XvueApp::blockingDepth() + AA_CompressHighFrequencyEvents land, and XvueState/XvueCanvas grow the accrochage + StrongFocus hooks Plans 02-06 depend on — all without changing any existing behavior.**

## Performance

- **Duration:** ~35 min
- **Started:** 2026-04-14 (Phase 5 kickoff)
- **Completed:** 2026-04-14
- **Tasks:** 3
- **Files modified:** 7 (+ 4 new)

## Accomplishments

- **Wave 0 test harness** wired end-to-end: `cmake --build xvue/build --target xvue_qt_event_tests` produces a binary that runs under `xvfb-run -a` with 5 PASS + 17 SKIP + 0 failed; `bin/cbl_tout_qt` builds it as part of the normal Qt build.
- **QSKIP test matrix** for EVENT-01..EVENT-08 (17 tests) so every downstream plan has a concrete `<automated>` target to hang real assertions on. 7 tests mention `testXvsouris*` per the acceptance criterion.
- **XvueApp::blockingDepth()** static accessor + `blockingDepth_` counter (zero-init), ready for Plan 02's RAII guard. `BlockingDepthGuard` friend declaration ensures only the guard may touch the counter.
- **Qt::AA_CompressHighFrequencyEvents** explicitly set inside `XvueApp::ensure()`'s `call_once` lambda BEFORE `std::make_unique<QApplication>(...)` — D-05 invariant satisfied.
- **XvueState::mempxaccro_** and **XvueState::accroche_undo_tile_** QPixmap* slots added with nullptr default and matching dtor deletes (mirrors the Phase 4 saved_canvas_ pattern).
- **XvueCanvas::setFocusPolicy(Qt::StrongFocus)** called in ctor — Pitfall 4 closed, the Plan 02 event filter will now see QKeyEvents.
- **XvueCanvas::resizeEvent** invalidates `accroche_undo_tile_` after the painter/backing swap (Pitfall 10 partially closed; Plan 05 completes the full accrochage lifecycle).
- **Both backends still build:** `bin/cbl_tout` and `bin/cbl_tout_qt` exit 0. SHELL-03 grep guard clean (0 matches in `xvue/qt/src/` + `xvue/qt/include/`).

## Task Commits

1. **Task 1: QTest skeleton + CMake wiring + cbl_tout_qt extension** — `ee5784a` (test)
2. **Task 2: XvueApp::blockingDepth() + AA_CompressHighFrequencyEvents** — `f0c878d` (feat)
3. **Task 3: XvueState accrochage fields + canvas StrongFocus + resize invalidation** — `4c017a5` (feat)

## Files Created/Modified

### Created
- `xvue/qt/tests/CMakeLists.txt` — QTest target `xvue_qt_event_tests` linking `xvueqt` + `Qt6::Test` + `Qt6::Widgets`
- `xvue/qt/tests/test_xvue_qt_event.cpp` — QTest fixture with 17 EVENT-NN stubs + 5 real PASSes + custom main() owning `XvueApp::ensure()`
- `xvue/qt/tests/test_helpers.h` / `test_helpers.cpp` — `pumpEvents()` + canvas helper placeholders (Plan 02 wires real canvas creation)

### Modified
- `xvue/qt/CMakeLists.txt` — added `option(XVUE_QT_BUILD_TESTS)` + `add_subdirectory(tests)` gate
- `bin/cbl_tout_qt` — pass `-DXVUE_QT_BUILD_TESTS=ON` to cmake configure; add explicit `cmake --build --target xvue_qt_event_tests` with existence check
- `xvue/qt/src/xvue_qt_app.h` — `blockingDepth()` accessor + `blockingDepth_` private static + `BlockingDepthGuard` friend decl
- `xvue/qt/src/xvue_qt_app.cpp` — static zero-init of `blockingDepth_`, accessor body, `Qt::AA_CompressHighFrequencyEvents` call inside `ensure()` pre-ctor
- `xvue/qt/src/xvue_qt_state.h` — `mempxaccro_` + `accroche_undo_tile_` QPixmap* fields (nullptr default)
- `xvue/qt/src/xvue_qt_state.cpp` — dtor deletes both new pointers
- `xvue/qt/src/xvue_qt_canvas.cpp` — `setFocusPolicy(Qt::StrongFocus)` in ctor; `resizeEvent` invalidates `accroche_undo_tile_`

## Decisions Made

- **Custom main() vs QTEST_MAIN:** Required to satisfy D-05. `QTEST_MAIN` constructs `QApplication` itself before any test runs, which would mean `XvueApp::ensure()` either creates a second QApplication (Qt assertion) or gets bypassed entirely (defeating the D-05 "set BEFORE ctor" invariant). The custom main calls `XvueApp::ensure()` first, then `QTest::qExec`. This is the reusable pattern for every downstream Plan 02+ test too.
- **Test opt-in via CMake option:** `XVUE_QT_BUILD_TESTS` defaults OFF so manual `cmake -S xvue/qt -B xvue/qt/build` still works unchanged; `cbl_tout_qt` flips it ON. Keeps the normal developer loop clean.
- **BlockingDepthGuard friend decl is in Plan 01, not Plan 02:** lets Plan 02 add the guard struct without re-editing `xvue_qt_app.h` (Rule: later plans do not re-touch these files).
- **`mempxaccro_` NOT invalidated on resize, `accroche_undo_tile_` IS:** the sprite template is resolution-independent (13x13 logical pixels), but the undo tile is a copy of live backing pixels whose coordinates may be out of bounds after a shrink.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Task 1 forward-decl probe symbol instead of deferring the test**
- **Found during:** Task 1 (planning said "known failing test tolerated before Task 2 lands")
- **Issue:** The plan's original skeleton had `extern int xvue_qt_test_blocking_depth_probe();` declared but no local definition, so the Task 1 binary would fail to link entirely — not just "fail one test cleanly". That would have made Task 1's acceptance unverifiable until Task 2 landed, defeating the whole point of atomic task commits.
- **Fix:** Added a local `int xvue_qt_test_blocking_depth_probe() { return 0; }` in the test TU as a weak placeholder. Task 2 then replaced it with `#include "xvue_qt_app.h"` + direct `XvueApp::blockingDepth()` call (and removed the forward decl + local definition).
- **Files modified:** `xvue/qt/tests/test_xvue_qt_event.cpp`
- **Verification:** Task 1 `xvue_qt_event_tests` binary built and ran clean with 3 PASS + 16 SKIP + exit 0 (no tolerated failure); Task 2 rebuilt and the test still passes.
- **Committed in:** `ee5784a` (Task 1 commit) + `f0c878d` (Task 2 commit)

**2. [Rule 2 - Missing Critical] Custom main() in place of QTEST_MAIN**
- **Found during:** Task 2 (planning the compress-attribute test)
- **Issue:** Plan's Task 2 action block said "replace `QTEST_MAIN` with real XvueApp::blockingDepth() call". But QTEST_MAIN constructs its own QApplication; calling `XvueApp::ensure()` inside a test body would then try to construct a second QApplication and trip Qt's `qt_assert("Only one QApplication allowed")`. Either the D-05 attribute test would be meaningless (measuring Qt's default, not our explicit setAttribute) or the binary would abort.
- **Fix:** Replaced `QTEST_MAIN(TestXvueQtEvent)` with a custom `int main(int argc, char* argv[])` that calls `XvueApp::ensure()` first, then constructs the fixture and calls `QTest::qExec`. Net result: XvueApp owns the QApplication, D-05's "set BEFORE ctor" invariant is genuinely satisfied, and every downstream test inherits the pattern.
- **Files modified:** `xvue/qt/tests/test_xvue_qt_event.cpp`
- **Verification:** `testCompressHighFrequencyEventsSet` PASSes under xvfb-run, `testBlockingDepthAccessorZero` PASSes, no second-QApplication abort.
- **Committed in:** `f0c878d` (Task 2 commit)

**3. [Rule 3 - Blocking] cbl_tout_qt MEFISTO env + gfortran shim requirement**
- **Found during:** Task 1 first build
- **Issue:** `bin/cbl_tout_qt` requires `MEFISTO` env var set and `gfortran`/`g++`/`gcc` on PATH. Fresh shell had neither. Per STATE.md, Debian sid has libgfortran pinned to 15.2.0-9 and the project uses `/tmp/gfortran-14-shim` to surface the right compiler binaries.
- **Fix:** Every build invocation used `export MEFISTO=/home/drico/git/mefisto && export PATH=/tmp/gfortran-14-shim:$MEFISTO/bin:/usr/bin:/bin`. This is an executor-environment fix, not a source change — no files modified.
- **Files modified:** none
- **Verification:** Build succeeds.
- **Committed in:** n/a

**Total deviations:** 3 auto-fixed (2 Rule-2/3 source fixes, 1 Rule-3 environment)
**Impact on plan:** All three fixes were essential for correctness. The QTEST_MAIN fix is the most load-bearing — it's the reusable pattern every Plan 02+ test will adopt.

## Acceptance Criteria Note

The literal acceptance regex `grep -c "int XvueApp::blockingDepth_" xvue/qt/src/xvue_qt_app.cpp` returned 0 because the definition line uses alignment padding (`int[spaces]XvueApp::blockingDepth_ = 0;`) to match the other `XvueApp::` static definitions that precede it (`std::once_flag`, `std::unique_ptr<QApplication>`, etc.). The literal regex required a single space. The semantic criterion — "static definition of XvueApp::blockingDepth_ with zero-init" — is fully satisfied (`grep "XvueApp::blockingDepth_" xvue/qt/src/xvue_qt_app.cpp` finds the line at col 36). No code change made; flagged here for the verifier.

## Issues Encountered

- **QThreadStorage warnings at test exit:** `QThreadStorage: entry 1 destroyed before end of thread` printed after `Totals:` line. These are a documented side effect of the D-08 "leak QApplication at atexit" pattern (see `xvue/qt/src/xvue_qt_app.cpp` teardown comment). Exit code is still 0. Harmless.
- **First build failed with `MEFISTO=` empty:** executor environment did not inherit the project env vars. See deviation #3 above.

## User Setup Required

None — Qt6 Test module already available via the existing `qt6-base-dev` install. No new apt packages.

## Next Phase Readiness

**Ready for Plan 02 (XvueEventBridge + waitForEvent skeleton):**
- Test target `xvue_qt_event_tests` exists and PASSes cleanly → every Plan 02 verify block can drop real assertions into the QSKIP bodies.
- `XvueApp::blockingDepth()` + `BlockingDepthGuard` friend decl → Plan 02 drops in the RAII guard in `xvue_qt_event.h` without re-touching `xvue_qt_app.h`.
- `XvueState::mempxaccro_` / `accroche_undo_tile_` → Plan 05 fills these in during `initaccrochage_` + `xvsouris2_` without schema changes.
- `XvueCanvas::setFocusPolicy(Qt::StrongFocus)` → Plan 02's `installEventFilter` will see QKeyEvents on first try.

**Open items for Plan 02:**
- Create `xvue_qt_event.h` with `BlockingDepthGuard` struct (increments/decrements `XvueApp::blockingDepth_` via friendship).
- Create `XvueEventBridge` class (QObject) and wire `installEventFilter` in `XvueWindow` ctor (or `XvueCanvas` ctor after canvas_ is assigned).
- Replace QSKIP body of `testBridgeInstallation`, `testXvsourisButtonPress`, `testXvsourisButtonRelease`, `testXvsourisKeyPress`, `testBlockingDepthRAII`, `testBlockingDepthNested` with real assertions.

## Test Results (Task 3 final)

```
Totals: 5 passed, 0 failed, 17 skipped, 0 blacklisted, 0ms
```

PASSes: `initTestCase`, `testBlockingDepthAccessorZero`, `testXvueStateAccrochageFieldsNull`, `testCompressHighFrequencyEventsSet`, `cleanupTestCase`.

## Self-Check: PASSED

All created files verified on disk:
- `xvue/qt/tests/CMakeLists.txt` FOUND
- `xvue/qt/tests/test_xvue_qt_event.cpp` FOUND
- `xvue/qt/tests/test_helpers.h` FOUND
- `xvue/qt/tests/test_helpers.cpp` FOUND

All three task commits present in `git log --oneline`:
- `ee5784a` FOUND (Task 1)
- `f0c878d` FOUND (Task 2)
- `4c017a5` FOUND (Task 3)

Both backends build green: `bin/cbl_tout_qt` and `bin/cbl_tout` both exit 0.
SHELL-03 grep guard: 0 matches for `QApplication::exec` / `qApp->exec()` in `xvue/qt/src/` + `xvue/qt/include/`.
Test binary: `xvfb-run -a xvue/qt/build/tests/xvue_qt_event_tests` → exit 0, 5 PASS + 17 SKIP.

---
*Phase: 05-event-bridge-blocking-reads*
*Plan: 01*
*Completed: 2026-04-14*
