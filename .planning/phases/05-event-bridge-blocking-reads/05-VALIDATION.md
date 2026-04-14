---
phase: 5
slug: event-bridge-blocking-reads
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-04-14
---

# Phase 5 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution. Derived from `05-RESEARCH.md` §10.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Qt `QTest` (ships with `qt6-base-dev`, no new apt packages). First C++ unit-test framework wired into the project — `xvue/qt/tests/` does not exist yet. |
| **Config file** | `xvue/qt/tests/CMakeLists.txt` — **Wave 0 creates this**. |
| **Quick run command** | `cmake --build xvue/build --target xvue_qt_event_tests && xvfb-run -a xvue/build/xvue_qt_event_tests` |
| **Full suite command** | `bin/cbl_tout_qt && xvfb-run -a xvue/build/xvue_qt_event_tests && xvfb-run -a bin/qt-capture.sh pp/ppmail_qt testa/pan2d /tmp/pan2d_qt.png` |
| **Estimated runtime** | ~5 seconds for unit tests; ~45 seconds for full suite |

---

## Sampling Rate

- **After every task commit:** Run `xvue_qt_event_tests` (<5s) + `bin/cbl_tout_qt` smoke build
- **After every plan wave:** Full test binary + `bin/xvtest-capture.sh pp/ppxvtest{1..4}_qt` + headless AUTOEXIT run of `pp/ppmail_qt testa/pan2d`
- **Before `/gsd-verify-work`:** Full suite green AND manual developer A/B drag test on `testa/pan2d` + `testa/torus` logged below
- **Max feedback latency:** 5 seconds (unit tests); 45 seconds (full wave)

---

## Per-Task Verification Map

Task IDs will be assigned by gsd-planner. This matrix maps every phase requirement to a concrete, automated check. Task IDs and plan numbers are filled in when PLAN.md files are written.

| Req | Behavior | Test Type | Automated Command / Assertion | File Exists | Status |
|-----|----------|-----------|-------------------------------|-------------|--------|
| EVENT-01 | `XvueEventBridge` exists, is a `QObject`, installed on `XvueCanvas` as event filter | unit | `testBridgeInstallation` — create XvueApp, open canvas, assert `bridge_ != nullptr` and installed via `canvas->installEventFilter(bridge_)` | ❌ W0 | ⬜ pending |
| EVENT-02 | `xvsouris_` returns correct `notypeevent / nbc / x1 / y1` for MotionNotify, ButtonPress, ButtonRelease, KeyPress | unit | `testXvsourisMotion`, `testXvsourisButtonPress`, `testXvsourisButtonRelease`, `testXvsourisKeyPress` — `QTest::mouseMove`, `QTest::mouseClick`, `QTest::keyClick` drive the filter; assert out-params | ❌ W0 | ⬜ pending |
| EVENT-02 | Esc (27) AND `@` (64) both map to `notypeevent = 0` (abort) — bit-for-bit X11 parity | unit | `testXvsourisEscapeAbort` with `QTest::keyClick(Qt::Key_Escape)`; `testXvsourisAtSignAbort` with `QTest::keyClick(Qt::Key_At)` and with literal `'@'` via `QTest::keyClicks` | ❌ W0 | ⬜ pending |
| EVENT-03 | `xvsouris2_` returns `notypeevent = 5` with `pmin0` updated during accrochage snap | unit | `testXvsouris2Accrochage` — populate `items[]` with 2 points, `QTest::mouseMove` near first point, assert `pmin0 == 0`; move near second, assert `pmin0 == 1` | ❌ W0 | ⬜ pending |
| EVENT-04 | `xvpause_` blocks until any event arrives, then returns | unit | `testXvpauseReturnsOnKey` with `QTest::keyClick(Qt::Key_Space)`; `testXvpauseReturnsOnMouseClick` with `QTest::mouseClick` | ❌ W0 | ⬜ pending |
| EVENT-04 | `xvpause_` honors `MEFISTO_XVSOURIS_AUTOEXIT` (headless short-circuit, NEW extension) | unit | `testXvpauseAutoexit` — `setenv("MEFISTO_XVSOURIS_AUTOEXIT", "1")`, call `xvpause_()`, assert returns immediately with no QEventLoop constructed | ❌ W0 | ⬜ pending |
| EVENT-05 | `deplsouris_` moves cursor without blocking, uses `QCursor::setPos(canvas->mapToGlobal(QPoint(x,y)))` | unit | `testDeplsourisNonBlocking` — call `deplsouris_(&x, &y)`, assert `QCursor::pos() == expected` under xvfb; assert call returns within 10ms | ❌ W0 | ⬜ pending |
| EVENT-06 | `initaccrochage_` allocates 13×13 `mempxaccro_` sprite pixmap on `XvueState` | unit | `testInitaccrochage` — call `initaccrochage_()`, assert `state->mempxaccro_ != nullptr && state->mempxaccro_->size() == QSize(13,13)` | ❌ W0 | ⬜ pending |
| EVENT-07 | Motion coalescing: fast drags produce a bounded number of returned motion events (≤ 20 for 100 input positions) | integration | `testMotionCoalescing` — `QTest::mouseMove` over 100 positions in a tight loop, drive `waitForEvent` repeatedly, assert the count of returned motion events is ≤ 20 (tunable; starts loose, tightens once A2 resolved) | ❌ W0 | ⬜ pending |
| EVENT-07 | Empirical A/B parity: `pp/ppmail_qt testa/pan2d` drag feels indistinguishable from `pp/ppmail testa/pan2d` | **human A/B** | Developer runs both binaries on the same machine, performs rapid rubber-band drags, records subjective verdict + event-count diagnostic (printf counter in eventFilter) in manual log below | manual | ⬜ pending |
| EVENT-08 | `XvueApp::blockingDepth()` increments on `waitForEvent` entry, decrements on exit, survives exceptions | unit | `testBlockingDepthRAII` — wrap `waitForEvent` in a try/catch, throw from within, assert `XvueApp::blockingDepth() == 0` after the catch | ❌ W0 | ⬜ pending |
| EVENT-08 | Nested `waitForEvent` call increments counter to 2, decrements back to 0 cleanly | unit | `testBlockingDepthNested` — enter `waitForEvent`, from event handler enter another `waitForEvent`, assert depths 1 and 2 observed via hook, then 0 after exit | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `xvue/qt/tests/CMakeLists.txt` — new QTest target `xvue_qt_event_tests`, linked against `xvueqt` static lib + `Qt6::Test` + `Qt6::Widgets`
- [ ] `xvue/qt/tests/test_xvue_qt_event.cpp` — full test file covering EVENT-01..EVENT-08 with the tests listed above (initial bodies may be stubs returning `QSKIP("not yet implemented")` to unblock the matrix)
- [ ] `xvue/qt/tests/test_helpers.{h,cpp}` — `runUntilBridgeReturns()` helper that drives the nested QEventLoop from the test side
- [ ] `bin/cbl_tout_qt` (or a sibling `bin/cbl_test_qt`) — extend to build the new test target so the quick-run command resolves
- [ ] `bin/xvtest-capture.sh` comment update — note that Phase 5 extends `MEFISTO_XVSOURIS_AUTOEXIT` to also short-circuit `xvpause_`
- [ ] Extend AUTOEXIT reader in both backends: `xvue/xvuelc.c::xvpause_` (X11) and `xvue/qt/src/xvue_qt_api.cpp::xvpause_` (Qt) — same env var, no new knob (§8 of RESEARCH.md)
- [ ] `xvue/qt/src/xvue_qt_canvas.cpp` — add `setFocusPolicy(Qt::StrongFocus)` in ctor (Pitfall 4) so keyboard events reach the canvas in real windows

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Rubber-band drag on `testa/pan2d` feels indistinguishable from X11 | EVENT-07 | Subjective "feels smooth" judgment cannot be automated; the Nyquist-flagged empirical validation the roadmap calls for | (1) `bin/cbl_tout && bin/cbl_tout_qt`. (2) Run `pp/ppmail testa/pan2d`, enter interactive zoom, perform 5 rapid rubber-band drags. (3) Run `pp/ppmail_qt testa/pan2d`, repeat. (4) Record subjective verdict in the log below. (5) Record motion-event counts per drag from the stderr printf counter in `eventFilter` (Task: add diagnostic counter). |
| Rubber-band drag on `testa/torus` feels indistinguishable from X11 | EVENT-07 | Same as above, second baseline case per ROADMAP success criterion #2 | As above, with `testa/torus`. |
| Accrochage highlight visual parity under mesh drawing | EVENT-03, EVENT-06 | Qt `QPainter::drawPixmap` vs X11 `GXand/GXorInverted` XOR raster-op may look subtly different — UX delta, not correctness (Assumption A5 in RESEARCH.md §12) | Run `pp/ppmail_qt testa/pan2d`, trigger accrochage via `xvsouris2_` code path (rubber-band pick). Compare side-by-side against X11 binary. Record whether the developer considers the difference acceptable. |
| `@` key aborts on French AZERTY keyboard | EVENT-02 | Assumption A3 — `QKeyEvent::text()` on AltGr+0 needs live-keyboard verification; defensively mitigated by `Qt::Key_At` in fallback switch but worth confirming on real hardware | Run `pp/ppmail_qt testa/pan2d` on a French AZERTY system, press AltGr+0, confirm the module returns to the parent menu the same way Esc does. |

### Manual A/B Log (to be filled during execution)

| Date | Case | X11 subjective feel | Qt subjective feel | X11 motion count / drag | Qt motion count / drag | Verdict | Tester |
|------|------|---------------------|--------------------|------------------------|------------------------|---------|--------|
| | testa/pan2d | | | | | ⬜ pending | |
| | testa/torus | | | | | ⬜ pending | |
| | testa/xvtest1 (sanity) | | | | | ⬜ pending | |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references (7 items above)
- [ ] No watch-mode flags
- [ ] Feedback latency < 5s for per-task sampling
- [ ] Manual A/B log has ≥ 3 entries before phase sign-off
- [ ] Assumption A2 (compression vs filter ordering) empirically resolved — diagnostic printf result recorded
- [ ] `nyquist_compliant: true` set in frontmatter once the above are complete

**Approval:** pending
