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
| EVENT-02 | `xvsouris_` returns correct `notypeevent / nbc / x1 / y1` for MotionNotify, ButtonPress, ButtonRelease, KeyPress | unit | `testXvsourisMotion`, `testXvsourisButtonPress`, `testXvsourisButtonRelease`, `testXvsourisKeyPress`, `testXvsourisFortranEntryPoint`, `testXvsourisNoWindow` — `QTest::mouseMove`, `QTest::mouseClick`, `QTest::keyClick` drive the filter; assert out-params | ✅ Plan 04 | ✅ green |
| EVENT-02 | Esc (27) AND `@` (64) both map to `notypeevent = 0` (abort) — bit-for-bit X11 parity | unit | `testXvsourisEscapeAbort` with `QTest::keyClick(Qt::Key_Escape)`; `testXvsourisAtSignAbort` with `QTest::keyClick(Qt::Key_At)` and with literal `'@'` via `QTest::keyClicks` | ❌ W0 | ⬜ pending |
| EVENT-03 | `xvsouris2_` returns `notypeevent = 5` with `pmin0` updated during accrochage snap | unit | `testXvsouris2Accrochage` — populate `items[]` with 2 points, `QTest::mouseMove` near first point, assert `pmin0 == 0`; move near second, assert `pmin0 == 1` | ❌ W0 | ⬜ pending |
| EVENT-04 | `xvpause_` blocks until any event arrives, then returns | unit | `testXvpauseReturnsOnKey` with `postKey(Qt::Key_Space)` (mouse-click variant skipped per xvuelc.c:2529 contract) | ✅ Plan 04 | ✅ green |
| EVENT-04 | `xvpause_` honors `MEFISTO_XVSOURIS_AUTOEXIT` (headless short-circuit, NEW extension) — extended to **both** Qt and X11 backends | unit | `testXvpauseAutoexit` — `qputenv("MEFISTO_XVSOURIS_AUTOEXIT", "1")`, call `xvpause_()`, assert elapsed < 500 ms; X11 side verified by grep on xvuelc.c | ✅ Plan 04 | ✅ green |
| EVENT-05 | `deplsouris_` moves cursor without blocking, uses `QCursor::setPos(canvas->mapToGlobal(QPoint(x,y)))` | unit | `testDeplsourisNonBlocking` — call `deplsouris_(&x, &y)`, assert elapsed < 100 ms; best-effort ±5 px tolerance on `QCursor::pos()` under xvfb (offscreen QPA may not move cursor at all) | ✅ Plan 04 | ✅ green |
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
- [x] `bin/xvtest-capture.sh` comment update — note that Phase 5 extends `MEFISTO_XVSOURIS_AUTOEXIT` to also short-circuit `xvpause_`
- [x] Extend AUTOEXIT reader in both backends: `xvue/xvuelc.c::xvpause_` (X11) and `xvue/qt/src/xvue_qt_api.cpp::xvpause_` (Qt) — same env var, no new knob (§8 of RESEARCH.md)
- [ ] `xvue/qt/src/xvue_qt_canvas.cpp` — add `setFocusPolicy(Qt::StrongFocus)` in ctor (Pitfall 4) so keyboard events reach the canvas in real windows

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Rubber-band drag on `testa/pan2d` feels indistinguishable from X11 | EVENT-07 | Subjective "feels smooth" judgment cannot be automated; the Nyquist-flagged empirical validation the roadmap calls for | (1) `bin/cbl_tout && bin/cbl_tout_qt`. (2) Run `pp/ppmail testa/pan2d`, enter interactive zoom, perform 5 rapid rubber-band drags. (3) Run `pp/ppmail_qt testa/pan2d`, repeat. (4) Record subjective verdict in the log below. (5) Record motion-event counts per drag from the stderr printf counter in `eventFilter` (Task: add diagnostic counter). |
| Rubber-band drag on `testa/torus` feels indistinguishable from X11 | EVENT-07 | Same as above, second baseline case per ROADMAP success criterion #2 | As above, with `testa/torus`. |
| Accrochage highlight visual parity under mesh drawing | EVENT-03, EVENT-06 | Qt `QPainter::drawPixmap` vs X11 `GXand/GXorInverted` XOR raster-op may look subtly different — UX delta, not correctness (Assumption A5 in RESEARCH.md §12) | Run `pp/ppmail_qt testa/pan2d`, trigger accrochage via `xvsouris2_` code path (rubber-band pick). Compare side-by-side against X11 binary. Record whether the developer considers the difference acceptable. |
| `@` key aborts on French AZERTY keyboard | EVENT-02 | Assumption A3 — `QKeyEvent::text()` on AltGr+0 needs live-keyboard verification; defensively mitigated by `Qt::Key_At` in fallback switch but worth confirming on real hardware | Run `pp/ppmail_qt testa/pan2d` on a French AZERTY system, press AltGr+0, confirm the module returns to the parent menu the same way Esc does. |

### Manual A/B Log (filled during Plan 06 Task 2)

| Date       | Case                   | X11 subjective feel                                          | Qt subjective feel                                                  | X11 motion count / drag                                                                     | Qt motion count / drag                                                                       | Verdict              | Tester |
|------------|------------------------|--------------------------------------------------------------|---------------------------------------------------------------------|---------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|----------------------|--------|
| 2026-04-18 | testa/pan2d            | Smooth rubber-band drags, no stutter on the X11 workstation. | Indistinguishable from X11 — feels identical on rapid drags.         | n/a (not measured via xev; Phase A automated smoke gives quantitative evidence, see below). | Paced input (3 ms/waypoint): `motion_count=1` per return across all 448 motion round-trips. | ✅ approved          | drico  |
| 2026-04-18 | testa/torus            | Smooth rubber-band drags; behaves identically to pan2d.      | Indistinguishable from X11 — same subjective feel as pan2d above.    | n/a (same rationale as pan2d).                                                              | Same distribution as pan2d — `motion_count=1` dominates on paced input.                      | ✅ approved          | drico  |
| 2026-04-18 | testa/xvtest1 (sanity) | Test driver exits cleanly, no graphical regression.          | Sanity test driver under `pp/ppxvtest1_qt` renders 32×32 grid; the 450-round-trip harness drives real `xvsouris_` returns through the bridge without crashing or dropping input. | n/a (driver is not an interactive drag target). | 449 real `xvsouris_` round-trips captured under MEFISTO_XVSOURIS_DEBUG=1; zero dropped events. | ✅ approved          | drico  |

**Session environment:** X11 session on drico's workstation (`XDG_SESSION_TYPE=x11`, `DISPLAY=:0`). Host compositor: stock X11. Wayland / XWayland variants not exercised in the subjective pass because D-09 documents the `QCursor::setPos` no-op caveat and the X11 baseline is the supported Phase 5 target.

**Accrochage visual parity (Assumption A5, EVENT-03 / EVENT-06):** Strategy B (Phase 4 save-restore blit) produces a crisp black 13×13 square highlight; developer verdict is that it reads cleaner than the X11 `GXand`/`GXorInverted` XOR highlight on light backgrounds and equally visible on dark backgrounds. Accepted as a UX improvement, not a regression. Strategy B verdict: ✅ approved.

**`@` key abort on French AZERTY (EVENT-02, Assumption A3):** Not exercised on a live AZERTY layout in this session. D-07 defensive fallback (`Qt::Key_At → 64`) is already wired in `translateKey()` alongside `QKeyEvent::text().toLatin1()`, so the `@` abort path has two independent triggers. Deferred to any future session that happens to run on AZERTY hardware; not a sign-off blocker.

### Phase A — Automated Empirical Evidence (Assumption A2 resolution)

Captured through kwin-mcp under KWin/XWayland on the live rebuilt `pp/ppxvtest1_qt` binary built from HEAD `461459c` (includes all Plans 05-01..05-05 code):

| Metric                                   | Value                                                                                                                                  |
|------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------|
| Binary                                   | `pp/ppxvtest1_qt` (built 2026-04-18 16:31 from `461459c`)                                                                              |
| Rebuild policy                           | `bin/cbl_tout` + `bin/cbl_tout_qt` run sequentially (see feedback_parallel_builds_share_pp).                                           |
| Total `xvsouris_` round-trips captured   | 449                                                                                                                                    |
| Event distribution                       | 448 motion (`notypeevent=-2`); 1 full click (`notypeevent=1 nbc=1`); 1 Esc abort (`notypeevent=0 nbc=27`)                              |
| `motion_count` distribution              | 448× `motion_count=1`, 2× `motion_count=0` (click + Esc — no motion preceded them in the window)                                       |
| Bridge depth (`XvueApp::blockingDepth`)  | `depth=1` across all 450+ returns — `BlockingDepthGuard` RAII zero-leakage verified                                                    |
| Crashes / deadlocks                      | 0 crashes, 0 deadlocks over 450+ Fortran-ABI round-trips                                                                               |
| Canvas visual state after drag           | Intact (32 colors × 32 fonts grid still visible in post-drag screenshot)                                                                |
| Log file                                 | `/tmp/phase05_kwinmcp/ppxvtest1_qt_motion.log` (1264 lines)                                                                            |
| Screenshot                               | `/tmp/phase05_kwinmcp/ppxvtest1_qt_initial.png`                                                                                        |

**Assumption A2 empirical resolution (RESEARCH §12 Q1/Q2 — compression vs filter ordering):**

- `Qt::AA_CompressHighFrequencyEvents` operates at `QWindowSystemInterface` level — **before** the object-specific event filter runs (confirmed by Qt source reading; matches the §4 [ASSUMED] note in RESEARCH.md).
- The bridge's own `QTimer::singleShot(0, loop, &QEventLoop::quit)` deferred-quit fires when the event queue drains, which is the X11-equivalent `XEventsQueued(QueuedAfterFlush)` coalescing.
- Paced input (~3 ms/waypoint in this harness): each motion reaches the filter individually, `motion_count=1` per `waitForEvent` return — matches X11 `XEventsQueued` semantics where the compositor delivered one motion per queue iteration.
- Burst input (QTest `0 ms` gap in Plan 03 `testMotionCoalescingBurst`): bridge collapses 100 moves into a single return (already proved in Plan 03 test suite — `motion_count=100` in one return).
- **Conclusion:** compression happens twice (Qt QWSI compression + bridge deferred-quit filter); they compose correctly. No lost final position. Whatever Qt's pre-filter behavior turns out to be on a given platform, the bridge's deferred-quit is the operative X11-parity mechanism.

This **resolves Assumption A2** (previously MEDIUM risk). The paired burst/paced empirical counts confirm that motion coalescing is load-bearing and working as designed.

**Note on the automated evidence vs the human verdict:** the Phase A numbers quantify *that* a real drag round-tripped end-to-end, they do not quantify "smoothness". The human verdict above is still the load-bearing subjective judgment for success criterion #5 ("A/B run indistinguishable"). The automated evidence supplements — it does not replace — the subjective verdict.

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
