# Phase 06.0 Deferred Items

Items discovered during plan execution that are out of scope for the current
plan and tracked here for a future hotfix or follow-up plan.

## Resolution summary (2026-04-20)

All items from the 2026-04-19 log below are now resolved in the 2026-04-20
post-phase hotfix batch:

| Item | Status | Resolution commit |
|------|--------|-------------------|
| #1 `testWheelZoomIn` flake | ✅ RESOLVED | `c533ca7` — replaced `postEvent + processEvents(ExcludeUserInputEvents)` with `sendEvent` (synchronous, bypasses the user-input filter). 3/3 consecutive runs green. |
| #2 `testEmptyStateRendersText` flake | ✅ RESOLVED | `c533ca7` — replaced `repaint()` on hidden widget with `grab()` (forces paintEvent into offscreen pixmap regardless of visibility). |
| #3 A6 follow-up (Qt bridge MMB-abort diverges from xvuelc.c) | ✅ RESOLVED | `c1869a7` — removed both MMB-abort branches in `XvueEventBridge::eventFilter` (Souris + Souris2). Added two regression tests. |

**Full suite post-fix:** 109 PASS / 0 FAIL / 2 SKIP (was 106–107 P / 1–2 F
/ 2 S). ABI still 58. `bin/cbl_tout` and `bin/cbl_tout_qt` both green.

---

## From Plan 06.0-03 execution (2026-04-19)

### Pre-existing test failures in `xvue_qt_canvas_gestures_tests` (Plan 05)

**Discovered during:** Plan 06.0-03 verification (running full test suite for
non-regression check after Plan 03 source changes).

**Status:** Confirmed pre-existing on baseline commit a523f8a (verified by
reverting Plan 03 source changes and re-running the test binary). NOT caused
by Plan 03.

**Failing tests:**
1. `TestXvueQtCanvasGestures::testWheelZoomIn`
   - Failing assertion: `state_.view_transform_.m11() > 1.0` returned FALSE.
   - Likely cause: the test posts a wheel event but the canvas's
     `wheelEvent` may not have run (event-posting under offscreen QPA does
     not always deliver the wheel event the same way as a real input
     device). Plan 05's source-side verification was sandbox-blocked per
     06.0-05-SUMMARY.md so this regression slipped through.
2. `TestXvueQtCanvasGestures::testEmptyStateRendersText`
   - Failing assertion: `canvas_->lastPaintDrewEmptyState_` returned FALSE.
   - Likely cause: the canvas paint event did not run between the
     `update()` call and the assertion in the test, OR the empty-state
     branch was not entered because state's `has_user_content_` is true
     for some reason. Same root cause class as #1 — Plan 05 source-side
     verification was incomplete.

**Scope:** Plan 03 does NOT touch `xvue/qt/src/xvue_qt_canvas.{h,cpp}` and
its diff has zero overlap with the Plan 05 surface. The two failures are
Plan 05 follow-ups, not Plan 03 deviations.

**Suggested follow-up:** A small hotfix plan (or Plan 06's wave) should:
1. Wrap each test's event posting in `QCoreApplication::processEvents()` +
   `QTest::qWait(50)` to give Qt time to deliver the wheel/paint event.
2. Or replace the wheel-event post with a direct `canvas->wheelEvent(&ev)`
   call to bypass Qt's event-delivery quirks under offscreen QPA.
3. Or add `QSignalSpy` on `canvas->update`-triggered paintEvent and
   await the paint via `QTRY_VERIFY` instead of straight `QVERIFY`.

The two failures are deterministic across multiple runs, so they're real
test quality issues (not flake), but they do not block any Plan 03
acceptance criterion.
