---
phase: 05-event-bridge-blocking-reads
verified: 2026-04-18T18:00:00Z
status: passed
score: 8/8
overrides_applied: 0
human_verification_completed:
  - item: "EVENT-07 empirical A/B parity (testa/pan2d + testa/torus subjective feel)"
    verdict: approved
    tester: drico
    date: 2026-04-18
    evidence: "05-VALIDATION.md Manual A/B Log — 3 rows approved (pan2d, torus, xvtest1 sanity); Phase A kwin-mcp 449 round-trips, motion_count=1, depth=1 throughout"
---

# Phase 5: Event Bridge & Blocking Reads — Verification Report

**Phase Goal:** Blocking calls (`xvsouris`, `xvpause`, `deplsouris`) run on a nested local `QEventLoop` with proper re-entrancy and motion coalescing — the architectural pivot that makes the mesher interactive end-to-end.
**Verified:** 2026-04-18T18:00:00Z
**Status:** passed
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | `XvueEventBridge` QObject exists, is installed as event filter on `XvueCanvas` at `XvueWindow` construction | VERIFIED | `xvue/qt/src/xvue_qt_event.h` declares `class XvueEventBridge : public QObject`; `xvue/qt/src/xvue_qt_window.cpp:24-25` creates bridge with parent=this and calls `canvas_->installEventFilter(bridge_)` |
| 2 | `xvsouris_` blocks on `waitForEvent()` and returns the correct mouse/keyboard event codes, never calling `QApplication::exec()` | VERIFIED | `xvue/qt/src/xvue_qt_api.cpp:990` calls `win->bridge()->waitForEvent(WaitMode::Souris)`; SHELL-03 grep returns 0 matches for `QApplication::exec` in `xvue/qt/src/` + `xvue/qt/include/`; unit tests `testXvsourisFortranEntryPoint`, `testXvsourisButtonPress`, `testXvsourisButtonRelease`, `testXvsourisKeyPress*`, `testXvsourisEscapeAbort`, `testXvsourisAtSignAbort` all PASS |
| 3 | `xvsouris2_` dispatches through `WaitMode::Souris2` with Strategy B accrochage save-restore-blit; `pmin0` holds the correct X11-semantic offset | VERIFIED | `xvue/qt/src/xvue_qt_api.cpp:1047-1048` calls `waitForEvent(Souris2, items, pmin0)`; `xvue_qt_event.cpp` contains `nearest_item_offset`, `save_tile_under`, `restore_tile`, `draw_sprite`, `cleanupAccrochage`; `testXvsouris2Accrochage` PASSes with `pmin0=3` (offset) on click near (100,100) |
| 4 | `xvpause_` blocks until any event arrives; `MEFISTO_XVSOURIS_AUTOEXIT` short-circuits both Qt and X11 backends | VERIFIED | Qt: `xvue/qt/src/xvue_qt_api.cpp:1124` calls `waitForEvent(Pause)`; X11: `xvuelc.c:2533-2546` mirrors the AUTOEXIT pattern byte-for-byte; `testXvpauseReturnsOnKey` and `testXvpauseAutoexit` PASS |
| 5 | `deplsouris_` moves the cursor without blocking; bounds-checked against ±32768 | VERIFIED | `xvue/qt/src/xvue_qt_api.cpp` uses `QCursor::setPos(canvas->mapToGlobal(QPoint(x,y)))` with T-05-04-01 bounds-check; `testDeplsourisNonBlocking` PASSes (elapsed < 100 ms) |
| 6 | `initaccrochage_` allocates a 13x13 `mempxaccro_` sprite with a 3-pixel black border on transparent background | VERIFIED | `xvue/qt/src/xvue_qt_api.cpp:316-325` allocates `new QPixmap(13, 13)`, fills transparent, draws black border; `testInitaccrochage` PASSes (center alpha=0, border alpha=255, red=0) |
| 7 | Mouse-motion events are coalesced via `QTimer::singleShot(0, loop_, &QEventLoop::quit)`; 100-move burst collapses to one `waitForEvent` return with tail coordinates | VERIFIED | `xvue/qt/src/xvue_qt_event.cpp:462,474` arms deferred-quit timer; `testMotionCoalescing` asserts `returns <= 20` (actual: 1) with `motion_count=100`; Phase A kwin-mcp 449 round-trips showed `motion_count=1` per return on paced input; Assumption A2 empirically resolved in `05-VALIDATION.md` |
| 8 | `XvueApp::blockingDepth()` counter increments on `waitForEvent` entry and decrements on exit via `BlockingDepthGuard` RAII; survives nested calls (depth 2 → 0) | VERIFIED | `xvue/qt/src/xvue_qt_event.h:21-34` `BlockingDepthGuard` struct increments/decrements `blockingDepth_`; `testBlockingDepthRAII` and `testBlockingDepthNested` PASS; Phase A evidence shows `depth=1` throughout 450+ round-trips |

**Score:** 8/8 truths verified

---

## Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `xvue/qt/src/xvue_qt_event.h` | `XvueEventBridge` class + `BlockingDepthGuard` + `WaitMode` + `Result` | VERIFIED | File present; Q_ASSERT depth < 4; all three WaitMode enumerators; 7 private filter state members |
| `xvue/qt/src/xvue_qt_event.cpp` | `waitForEvent` + `eventFilter` + motion coalescing + Strategy B helpers | VERIFIED | `QTimer::singleShot` (2x), `motion_count_` (6x), `WaitMode::Souris2` (7x), `accroche_undo_tile_` (13x), `nearest_item` (3x) |
| `xvue/qt/src/xvue_qt_window.cpp` | Bridge instantiation + `installEventFilter` on canvas | VERIFIED | Line 24-25: `bridge_ = new XvueEventBridge(canvas_, this); canvas_->installEventFilter(bridge_)` |
| `xvue/qt/src/xvue_qt_api.cpp` | Real bodies for `xvsouris_`, `xvsouris2_`, `xvpause_`, `deplsouris_`, `initaccrochage_` | VERIFIED | All 5 entry points have real bodies dispatching through the bridge; AUTOEXIT short-circuits present; T-05-04-01 bounds-check in `deplsouris_` |
| `xvue/qt/src/xvue_qt_app.h` | `blockingDepth()` accessor + `blockingDepth_` private static + `BlockingDepthGuard` friend | VERIFIED | Lines 31-46 confirmed |
| `xvue/qt/src/xvue_qt_app.cpp` | Static zero-init of `blockingDepth_` + `AA_CompressHighFrequencyEvents` set pre-ctor | VERIFIED | Lines 32+34+62 confirmed |
| `xvue/qt/src/xvue_qt_state.h` | `mempxaccro_` and `accroche_undo_tile_` QPixmap* fields | VERIFIED | Lines 50-51 confirmed, both `= nullptr` default |
| `xvue/qt/src/xvue_qt_canvas.cpp` | `setFocusPolicy(Qt::StrongFocus)` + `setMouseTracking(true)` + resize-invalidates-tile | VERIFIED | Lines 26, 36, 134-136 confirmed |
| `xvue/xvuelc.c` | `xvpause_` AUTOEXIT short-circuit mirroring Qt backend | VERIFIED | Lines 2527-2546 confirmed |
| `xvue/qt/tests/test_xvue_qt_event.cpp` | 31 passing tests, 2 documented QSKIPs | VERIFIED | Final test count per 05-05-SUMMARY.md: 31 passed, 0 failed, 2 skipped |
| `xvue/qt/README.md` | Phase 5 event bridge architecture + D-09 Wayland caveat | VERIFIED | File present; `Wayland` x5, `QEventLoop` x4, `MEFISTO_XVSOURIS_DEBUG` x4, `D-09` x2 |

---

## Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `xvsouris_` / `xvpause_` / `xvsouris2_` (xvue_qt_api.cpp) | `XvueEventBridge::waitForEvent()` | `win->bridge()->waitForEvent(WaitMode::X)` | WIRED | xvue_qt_api.cpp:990, 1047, 1124 confirmed; all three WaitMode variants used |
| `XvueEventBridge::waitForEvent` | `XvueApp::blockingDepth_` | `BlockingDepthGuard` RAII increment/decrement | WIRED | xvue_qt_event.h:21-34 + xvue_qt_event.cpp:213 `BlockingDepthGuard depth_guard` confirmed |
| `XvueWindow` constructor | `XvueEventBridge` + `canvas_->installEventFilter` | ctor post-setCentralWidget | WIRED | xvue_qt_window.cpp:24-25 confirmed |
| `eventFilter MouseMove case` | `QEventLoop::quit` via `QTimer::singleShot(0)` | deferred-quit coalescing pattern | WIRED | xvue_qt_event.cpp:462 + 474 confirmed |
| `xvue_qt_event.cpp Souris2 branch` | `XvueState::accroche_undo_tile_` + `mempxaccro_` | Strategy B save-restore-blit helpers | WIRED | `save_tile_under`, `restore_tile`, `draw_sprite` at xvue_qt_event.cpp:107-145 confirmed |
| `xvuelc.c::xvpause_` | `MEFISTO_XVSOURIS_AUTOEXIT` | `getenv` per-call | WIRED | xvuelc.c:2533-2546 confirmed; matches Qt-backend contract |

---

## Data-Flow Trace (Level 4)

Not applicable to this phase. The event bridge is a control-flow pivot, not a data-rendering component. It receives Qt events and maps them to Fortran integer out-parameters — there is no separate "data source" to trace upstream.

---

## Behavioral Spot-Checks

| Behavior | Verification Basis | Result | Status |
|----------|--------------------|--------|--------|
| `xvsouris_` returns correct event codes for button/key | `testXvsourisFortranEntryPoint`, `testXvsourisButtonPress/Release`, `testXvsourisKeyPress*` (5 tests) — 31/31 PASS | All relevant tests PASS | PASS |
| Motion coalescing: 100 moves → 1 `waitForEvent` return | `testMotionCoalescing` asserts `returns <= 20` (actual: 1); `motion_count=100` confirmed | PASS | PASS |
| `BlockingDepthGuard` RAII keeps `blockingDepth()` at 0 after return | `testBlockingDepthRAII`, `testBlockingDepthNested`; Phase A live evidence: `depth=1` on all 450+ round-trips, returns to 0 | PASS | PASS |
| `initaccrochage_` 13x13 sprite correct | `testInitaccrochage` checks center alpha=0, border alpha=255+red=0 | PASS | PASS |
| `MEFISTO_XVSOURIS_AUTOEXIT` exits both backends headlessly | `testXvpauseAutoexit`; X11 smoke `MEFISTO_XVSOURIS_AUTOEXIT=1 timeout 15 pp/ppxvtest0` exits 0 per 05-04-SUMMARY.md | PASS | PASS |
| `QApplication::exec()` never called | `grep -r "QApplication::exec" xvue/qt/src/ xvue/qt/include/` = 0 matches (SHELL-03) | 0 matches confirmed | PASS |
| ABI symbol count frozen at 57 | `nm libxvueqt.a` pipeline in `verify_abi.sh` (patched in 05-02 to exclude mangled C++ names) | 57 confirmed per 05-04-SUMMARY.md | PASS |

---

## Requirements Coverage

| Requirement | Source Plans | Description | Status | Evidence |
|-------------|-------------|-------------|--------|---------|
| EVENT-01 | 05-01, 05-02 | `XvueEventBridge` QObject installed as event filter on `XvueCanvas` via `waitForEvent()` | SATISFIED | `xvue_qt_event.{h,cpp}` + `xvue_qt_window.cpp` wiring; `testBridgeInstallation` PASS |
| EVENT-02 | 05-02, 05-04 | `xvsouris_` blocks on `waitForEvent()` returning correct event types without `QApplication::exec()` | SATISFIED | `xvue_qt_api.cpp:990`; 7 xvsouris_ tests PASS; SHELL-03 clean |
| EVENT-03 | 05-05 | `xvsouris2_` returns `notypeevent=5` with `pmin0` updated on accrochage press | SATISFIED | `xvue_qt_event.cpp` Souris2 branch; `testXvsouris2Accrochage` PASS (pmin0=3) |
| EVENT-04 | 05-04 | `xvpause_` blocks until any event; AUTOEXIT parity in both backends | SATISFIED | `xvue_qt_api.cpp:1124`; `xvuelc.c:2533-2546`; `testXvpauseReturnsOnKey` + `testXvpauseAutoexit` PASS |
| EVENT-05 | 05-04 | `deplsouris_` returns current mouse position without blocking | SATISFIED | `xvue_qt_api.cpp` QCursor::setPos with bounds-check; `testDeplsourisNonBlocking` PASS (elapsed < 100 ms) |
| EVENT-06 | 05-01, 05-05 | `initaccrochage_` initializes 13x13 snap/crosshair sprite on `XvueState` | SATISFIED | `xvue_qt_api.cpp:310-332`; `testInitaccrochage` PASS |
| EVENT-07 | 05-03, 05-06 | Motion coalescing via `QTimer::singleShot(0, quit)` + empirical A/B parity on pan2d/torus | SATISFIED | Automated: `testMotionCoalescing` (100 moves → 1 return); Human: `05-VALIDATION.md` Manual A/B Log approved by drico 2026-04-18 (pan2d + torus + xvtest1 sanity); Phase A kwin-mcp 449 round-trips, motion_count=1, depth=1 throughout; Assumption A2 empirically resolved |
| EVENT-08 | 05-01, 05-02 | `XvueApp::blockingDepth()` tracks nested `waitForEvent()` calls via RAII | SATISFIED | `BlockingDepthGuard` in `xvue_qt_event.h`; `testBlockingDepthRAII` + `testBlockingDepthNested` PASS; Phase A depth=1 across 450+ live round-trips |

---

## Anti-Patterns Found

| File | Pattern | Severity | Impact |
|------|---------|----------|--------|
| `xvue/qt/tests/test_helpers.cpp` | `createTestCanvas()` always returns nullptr; `destroyTestCanvas()` is empty — documented dead code | Info (IN-02 from code review) | None — these helpers are unused; test file uses `XvueApp::window_slot()->canvas()` directly |
| `xvue/qt/src/xvue_qt_event.h:27` | `Q_ASSERT(blockingDepth_ < 4)` compiles out in release builds | Warning (WR-03 from code review) | Defense-in-depth guard absent in production; noted for Phase 6 hardening — not a current blocker since the depth invariant is proven by Phase A live evidence |
| `xvue/qt/src/xvue_qt_event.cpp:275-283` | Souris2 middle-button press short-circuits to `notypeevent=0, nbc=2` instead of falling through to the accrochage block (diverges from X11 semantics when saclav.f expects button=2 through the press path) | Warning (WR-01 from code review) | No current mesher call pattern exercises nested middle-button in Souris2 mode; saclav.f prints "PROBLEME AVEC LA SOURIS" on this path — not triggered in practice |

No blockers found. All anti-patterns are documented in `05-REVIEW.md`. Code-review warnings WR-01, WR-02, WR-03, WR-04 are all classified as non-critical by the review (0 critical findings); none prevent the phase goal from being achieved.

---

## Human Verification Status

EVENT-07 empirical A/B parity was the sole item requiring human judgment. It has been completed:

**Completed (not pending):**
- `testa/pan2d` — Qt subjective feel "indistinguishable from X11 on rapid drags" — approved 2026-04-18 by drico
- `testa/torus` — Same subjective verdict — approved 2026-04-18 by drico
- `testa/xvtest1` sanity — 449 automated `xvsouris_` round-trips captured under kwin-mcp/KWin+XWayland, zero dropped events, depth=1 throughout — approved 2026-04-18 by drico

Automated supplement (Assumption A2 resolution): `QTimer::singleShot(0, loop, &QEventLoop::quit)` is the operative coalescing mechanism regardless of `AA_CompressHighFrequencyEvents` platform behavior. Proven by the combination of Plan 03 burst test (100 moves → 1 return, `motion_count=100`) and Phase A paced-input evidence (`motion_count=1` per return at 3 ms/waypoint).

**No open human verification items remain.**

---

## Phase Narrative

Phase 5 delivered the complete event-bridge pivot across 6 plans and 15 commits:

- **Plan 01** (Wave 0 scaffold): QTest infrastructure, `blockingDepth_` counter, `XvueState` accrochage fields, `XvueCanvas` StrongFocus + resize-invalidation, `AA_CompressHighFrequencyEvents` set pre-ctor. 5 PASS + 17 SKIP.
- **Plan 02** (Infrastructure): `XvueEventBridge` QObject + `BlockingDepthGuard` RAII + `waitForEvent` with save-restore of 7 filter members for nested-call safety + synchronous button/key dispatch + `XvueWindow::bridge()` accessor. 17 PASS + 8 SKIP.
- **Plan 03** (Motion coalescing): `QTimer::singleShot(0, loop_, &QEventLoop::quit)` deferred-quit on every `QEvent::MouseMove`; `MEFISTO_XVSOURIS_DEBUG` diagnostic counter; `setMouseTracking(true)` in canvas ctor. 21 PASS + 7 SKIP.
- **Plan 04** (Fortran ABI bodies): Real `xvsouris_`, `xvpause_`, `deplsouris_` bodies routing through bridge; `AUTOEXIT` extended to `xvpause_` in **both** backends; T-05-04-01 bounds-check on `deplsouris_`. 26 PASS + 4 SKIP.
- **Plan 05** (Accrochage Strategy B): Real `initaccrochage_` (13×13 sprite) + `xvsouris2_` via `WaitMode::Souris2`; `nearest_item_offset` follows real X11 items[] layout (PLAN hypothesis was wrong, corrected on source audit); `cleanupAccrochage` on every terminating path. 31 PASS + 2 SKIP (both documented not-load-bearing).
- **Plan 06** (Closure): Manual A/B verdict approved; Assumption A2 resolved empirically; `nyquist_compliant: true` flipped; `xvue/qt/README.md` with D-09 Wayland caveat created.

All 8 EVENT-NN requirements are satisfied. Both `bin/cbl_tout` (X11) and `bin/cbl_tout_qt` (Qt) exit 0. ABI frozen at 57 symbols.

---

_Verified: 2026-04-18T18:00:00Z_
_Verifier: Claude (gsd-verifier)_
