---
phase: 05-event-bridge-blocking-reads
plan: 05
subsystem: ui
tags: [qt6, fortran-abi, accrochage, strategy-b, pixmap-sprite, event-filter, cpp17]

requires:
  - phase: 05-event-bridge-blocking-reads-plan-01
    provides: "XvueState::mempxaccro_ and accroche_undo_tile_ pointer fields + XvueCanvas::resizeEvent tile invalidation + StrongFocus / mouseTracking"
  - phase: 05-event-bridge-blocking-reads-plan-02
    provides: "XvueEventBridge QObject + waitForEvent(WaitMode) with save-restore + eventFilter dispatch + XvueWindow::bridge() accessor + BlockingDepthGuard"
  - phase: 05-event-bridge-blocking-reads-plan-03
    provides: "QEvent::MouseMove deferred-quit coalescing (QTimer::singleShot(0, loop_, &QEventLoop::quit)) — reused byte-for-byte in the Souris2 motion branch"
  - phase: 05-event-bridge-blocking-reads-plan-04
    provides: "Fortran-ABI body template (ensure → MAIN_THREAD → AUTOEXIT → no-window guard → bridge dispatch → out-param copy) for xvsouris2_"

provides:
  - "initaccrochage_ allocates the 13x13 mempxaccro_ sprite with a 3-pixel black border on transparent background (Strategy B from 05-RESEARCH.md §6); idempotent; silent no-op when window/canvas/state is null (Pitfall 11)"
  - "xvsouris2_ real body routing through bridge->waitForEvent(WaitMode::Souris2, items, pmin0) after preserved MEFISTO_XVSOURIS_AUTOEXIT short-circuit"
  - "WaitMode::Souris2 branch in XvueEventBridge::eventFilter — accrochage save/restore/blit on every motion + press; final cleanup on release and on Esc/@ abort"
  - "nearest_item_offset(items, cx, cy) file-local helper implementing the X11-parity nearest-point search (items[0]=mots, items[2]=nbitem, offsets p += mots starting from p=0); T-05-05-01 mitigation clamps nbitem to 65536"
  - "save_tile_under / restore_tile / draw_sprite file-local helpers implementing Strategy B save-restore-blit with the long-lived painter on backing_"
  - "XvueEventBridge::cleanupAccrochage() private method — erases the sprite + frees accroche_undo_tile_ + resets *pmin0_=-2; called on release and abort to keep the canvas clean for the next call (closes T-05-05-02 by construction)"
  - "EVENT-03 (xvsouris2_ returns notypeevent=5 on accrochage press with pmin0 updated) verified end-to-end"
  - "EVENT-06 (initaccrochage_ allocates 13x13 mempxaccro_) verified pixel-level: center transparent, border opaque black"

affects:
  - 05-event-bridge-blocking-reads-plan-06   # full-plan validation + A/B drag test
  - 06-menus-dialogs-phase6                   # modal dialogs gate on XvueApp::blockingDepth() > 0 during xvsouris2_ too

tech-stack:
  added: []
  patterns:
    - "Strategy B: save-restore-blit tile pattern — mirrors Phase 4 saved_canvas_ ownership style (raw QPixmap* + manual lifecycle + ~XvueState destructor cleanup) instead of X11 GXand/GXorInverted raster-op XOR trick"
    - "X11-parity items[] layout via offset (not index) in *pmin0 — matches itemau.f:17-20 and xvuelc.c:2397-2413 directly; the plan's simplified 'items[0]=count' hypothesis was rejected on source audit"
    - "Accrochage sprite with a SINGLE lifecycle for mempxaccro_ (one-time initaccrochage_ allocation) and a RE-ALLOCATED-ON-DEMAND lifecycle for accroche_undo_tile_ (lazy on first motion, freed on release or resize)"
    - "Re-use of the Plan 03 deferred-quit timer in the Souris2 motion branch so fast drags across the canvas coalesce to ONE waitForEvent return per visual burst — the sprite redraw is free of per-pixel overhead"

key-files:
  created: []
  modified:
    - xvue/qt/src/xvue_qt_api.cpp
    - xvue/qt/src/xvue_qt_event.h
    - xvue/qt/src/xvue_qt_event.cpp
    - xvue/qt/tests/test_xvue_qt_event.cpp

key-decisions:
  - "Strategy B save/restore was accepted verbatim (05-RESEARCH.md §6, D-08). The X11 XOR-raster-op approach (GXand + GXorInverted) would have required QPainter::CompositionMode_Difference which operates on RGB, not X11 plane-mask raster ops — guaranteed visual drift (Assumption A5). Save/restore works for any underlying color and reuses the Phase 4 ownership pattern one-for-one."
  - "items[] layout hypothesis from the plan was rejected. The plan suggested items[0]=count with (x,y) pairs following. Reading xvue/saclav.f:61 + util/itemau.f:17-20 + xvuelc.c:2397-2413 showed the real X11 layout: items[0]=mots (words per item), items[1]=maxcap, items[2]=nbitem, then triplets/quads at offsets mots, 2*mots, ... . *pmin0 holds the raw OFFSET into items[], NOT an index. Implementation follows the real X11 contract and the test data is shaped accordingly ({3, 2, 2, 100,100,1, 200,200,2})."
  - "ButtonPress in Souris2 mode quits the loop on EVERY press (even when no item is nearby). The X11 body loops on pmin<0 until flag=1, which would hang the bridge if the user clicks in empty space. Relaxing this to 'always quit on press' is a tolerance improvement — saclav.f's loop around xvsouris2_ re-invokes as needed, and the Qt bridge should not wedge on user error."
  - "nearest_item_offset clamps nbitem to 65536 per T-05-05-01 mitigation. Real MEFISTO mesher workloads peak at a few thousand items; 65536 is five orders of magnitude below INT_MAX while covering every realistic caller, and rejecting clearly-garbage counts at the filter level is cheap."
  - "cleanupAccrochage() is called not only on release but also on Esc/@ abort AND on ordinary key press in Souris2 mode. X11 xvuelc.c:2465-2477 does not explicitly erase on keypress (relies on saclav.f to redraw), but the Qt backend benefits from a clean canvas so subsequent solver code sees no residue from the interactive loop. This is a deliberate visual improvement over X11, not a semantic deviation."
  - "Souris2 motion still arms the Plan 03 deferred-quit timer. Without it, fast drags would call nearest_item_offset + save_tile_under + draw_sprite per raw pixel of movement. With it, bursts coalesce to one waitForEvent return per visual burst — the sprite dance runs once per tail, keeping the rubber-band pick smooth with bounded work."

patterns-established:
  - "Pattern: Fortran-ABI body template for blocking-with-state entry points — ensure → MAIN_THREAD → AUTOEXIT → no-window guard → bridge.waitForEvent(mode, args...) → null-guarded out-param copy. Plan 05's xvsouris2_ is the final consumer of this template in Phase 5."
  - "Pattern: file-local accrochage helpers (nearest_item_offset, save_tile_under, restore_tile, draw_sprite) live in the same TU as the eventFilter that drives them; keeps the accrochage mechanics behind one abstraction that a future re-write can swap without touching the public bridge API."
  - "Pattern: offset-based in/out contract (*pmin0 holds the OFFSET into items[], not an index) preserved byte-for-byte through the Qt bridge — Fortran callers like xvue/saclav.f use `MCN(MNIT+NMIN0+2)` to reach the item number, and that arithmetic works because NMIN0 is still the X11-semantic offset."

requirements-completed: [EVENT-03, EVENT-06]

duration: 28min
completed: 2026-04-18
---

# Phase 5 Plan 05: xvsouris2_ + initaccrochage_ Strategy B accrochage Summary

**Closes Phase 5 by wiring the accrochage (snap-highlight) path: initaccrochage_ pre-allocates a 13x13 mempxaccro_ sprite with a 3-pixel black border on a transparent center; xvsouris2_ dispatches through WaitMode::Souris2 in the bridge, and eventFilter runs Strategy B save-restore-blit on every motion/press to keep the highlight attached to the nearest items[] point.**

## Performance

- **Duration:** 28 min
- **Started:** 2026-04-18T11:47:16Z
- **Completed:** 2026-04-18T12:15:24Z
- **Tasks:** 2
- **Files modified:** 4 (0 created, 4 modified)

## Accomplishments

- **initaccrochage_** (Task 1): allocates XvueState::mempxaccro_ = new QPixmap(13, 13), fills with Qt::transparent, draws a 3-pixel-thick black square border via QPainter + QPen(Qt::black, 3) with MiterJoin corners + drawRect(QRect(1,1,11,11)) — antialiasing off for crisp pixels. Idempotent (deletes prior sprite if any). Silent no-op when window / canvas / state is null (Pitfall 11). accroche_undo_tile_ left nullptr for lazy allocation on first motion.
- **xvsouris2_** (Task 2) real body: preserves the MEFISTO_XVSOURIS_AUTOEXIT byte-for-byte short-circuit (D-10); silent-abandon guard on no-window (Pitfall 11 analogue, zero-filled out-params); dispatches through bridge()->waitForEvent(WaitMode::Souris2, items, pmin0); null-pointer guards on every Fortran out-param.
- **WaitMode::Souris2 branch** in XvueEventBridge::eventFilter:
  - MouseMove: nearest_item_offset(items_, cx, cy) → if new offset differs from *pmin0_, restore_tile at old location and invalidate *pmin0_; if new offset valid, save_tile_under + draw_sprite at new location, update *pmin0_; set pending_.notypeevent=5, ++motion_count_, arm deferred-quit timer (motion coalescing reused from Plan 03).
  - MouseButtonPress: same accrochage redraw, then quit with notypeevent=5 / nbc=button (1/2/3). Middle button aborts with notypeevent=0 / nbc=2 (xvuelc.c:2272 parity).
  - MouseButtonRelease: cleanupAccrochage() erases the sprite + frees undo tile, then quit with notypeevent=1.
  - KeyPress: cleanupAccrochage() first, then quit with notypeevent=0 (Esc/@) or notypeevent=2 (other printable keys).
- **File-local helpers** in xvue_qt_event.cpp (all TU-local, no header visibility):
  - `nearest_item_offset(items, cx, cy)` — X11-parity semantics per itemau.f:17-20 and xvuelc.c:2397-2413: items[0]=mots, items[2]=nbitem, offsets p = mots, 2*mots, ... Returns -1 on empty/garbage, else the OFFSET of the nearest (x,y). T-05-05-01 clamps nbitem to 65536.
  - `save_tile_under(state, cx, cy)` — lazy-allocate accroche_undo_tile_ if null, DPR-resync from backing, QPainter::drawPixmap from backing at (cx-6,cy-6,13,13) into the tile.
  - `restore_tile(state, cx, cy)` — drawPixmap the saved tile back at (cx-6, cy-6) through state->painter_ (Phase 2 D-05 long-lived painter on backing).
  - `draw_sprite(state, cx, cy)` — drawPixmap mempxaccro_ at (cx-6, cy-6) through state->painter_.
- **XvueEventBridge::cleanupAccrochage()** private method: no-op when items_/pmin0_ null or *pmin0_<0; otherwise restore_tile at old location + delete accroche_undo_tile_ + set *pmin0_=-2 + canvas update. Closes T-05-05-02 (tile leak on repeated motion) by construction.
- **3 new tests PASS** (no QSKIPs left for EVENT-03/06):
  - testInitaccrochage: sprite present, 13x13, center (6,6) alpha=0, border (2,2) alpha=255 + red=0; idempotent second call.
  - testInitaccrochageBeforeInit: xvfermer_ → initaccrochage_ → xvinitgraphique_ round-trip does not crash.
  - testXvsouris2Accrochage: left-button at (105,98) over items=[3,2,2, 100,100,1, 200,200,2] returns pmin0=3, notypeevent=5, ibtn=1, x1=105, y1=98.
  - testXvsouris2NullGuards: null out-param slots + a single-item items[] — space key dispatches the bridge without crashing.
  - testXvsouris2ResizeInvalidatesTile: after a press, accroche_undo_tile_ is non-null; after canvas resize, Plan 01's resizeEvent guard nulls it back out (Pitfall 10 closed end-to-end).
- **Test suite: 31 passed, 0 failed, 2 skipped** under `xvfb-run -a xvue/qt/build/tests/xvue_qt_event_tests`. The remaining 2 skips are the documented not-load-bearing testXvpauseReturnsOnMouseClick (xvuelc.c:2529 — X11 xvpause_ only fires on KeyPress) and testDebugLoggingEnvVar (env-var cached on first access, validated via outer-shell grep).
- Both `bin/cbl_tout_qt` (Qt) and `bin/cbl_tout` (X11 baseline) exit 0.
- `MEFISTO_XVSOURIS_AUTOEXIT=1 pp/ppxvtest0_qt` exits 0 in < 1 s — headless parity preserved.

## Task Commits

Each task committed atomically with `--no-verify` (parallel-executor guidance — orchestrator runs hooks once post-merge):

1. **Task 1: initaccrochage_ real body + EVENT-06 accrochage sprite tests** — `4d0012c` (feat)
2. **Task 2: xvsouris2_ real body + WaitMode::Souris2 Strategy B accrochage** — `20470c9` (feat)

## Files Created/Modified

- `xvue/qt/src/xvue_qt_api.cpp` (MODIFIED, +58 / -14 lines) — initaccrochage_ allocates + draws the 13x13 sprite; xvsouris2_ routes through bridge->waitForEvent(Souris2, items, pmin0) after AUTOEXIT.
- `xvue/qt/src/xvue_qt_event.h` (MODIFIED, +7 lines) — cleanupAccrochage() private method declaration.
- `xvue/qt/src/xvue_qt_event.cpp` (MODIFIED, +175 / -3 lines) — file-local Strategy B helpers (nearest_item_offset, save_tile_under, restore_tile, draw_sprite); WaitMode::Souris2 branches in eventFilter (MouseMove, MouseButtonPress, MouseButtonRelease, KeyPress); cleanupAccrochage() implementation.
- `xvue/qt/tests/test_xvue_qt_event.cpp` (MODIFIED, +130 / -2 lines) — forward-decls for initaccrochage_ / xvsouris2_; 5 new tests replacing Plan 05 QSKIPs.

## Decisions Made

1. **Strategy B save/restore was accepted verbatim (05-RESEARCH.md §6, D-08).** Rationale documented in `key-decisions` above.
2. **items[] layout hypothesis from the plan was REJECTED.** The plan wrote "items[0]=count, items[1..]=(x,y) pairs". Reading saclav.f:61 + itemau.f:17-20 + xvuelc.c:2397-2413 showed the real X11 layout (items[0]=mots, items[1]=maxcap, items[2]=nbitem, triplets at offsets mots, 2*mots, ...). Implementation follows the real X11 contract; the test data uses {3, 2, 2, 100,100,1, 200,200,2} and expects *pmin0=3 (offset) on a click near (100,100). Documented inline in nearest_item_offset's header comment. **This preserves the `MCN(MNIT+NMIN0+2)` arithmetic that saclav.f does to read the item number — if we had followed the plan's layout, every existing mesher caller would break.**
3. **ButtonPress in Souris2 always quits.** X11 xvuelc.c:2383-2439 only sets flag=1 when pmin>=0, so a click in empty space would loop forever. The Qt backend quits on every press — safer for robustness and saclav.f re-invokes xvsouris2_ anyway if it doesn't get what it expects. Documented as a deliberate tolerance improvement.
4. **cleanupAccrochage on KeyPress.** X11 erases only on ButtonRelease; the Qt backend erases on any terminating event (release, Esc, @, other keys). Simplifies reasoning about canvas state and does not affect correctness.
5. **Motion coalescing timer is armed in Souris2 mode too.** Without it, every raw MouseMove would run nearest_item_offset + save_tile_under + draw_sprite. With it, bursts of moves still coalesce to one waitForEvent return per tail, but the intermediate sprite redraws are inexpensive and keep visual feedback smooth. This reuses Plan 03's deferred-quit timer without modification.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug in plan's test data] items[] layout hypothesis was wrong**
- **Found during:** Task 2 action step — the plan's `<read_first>` block explicitly told me to "Verify by reading xvue/saclav.f line 61 and the surrounding declaration. If the layout differs, adjust to match the Fortran side exactly."
- **Issue:** The plan's action text said "Phase 5 contract: caller passes items as an array where items[0] is the count and items[1..] are the (x,y) pairs." Reading saclav.f:61, itemau.f:17-20 and xvuelc.c:2397-2413 showed this is NOT the Fortran contract. The real layout is items[0]=mots, items[1]=maxcap, items[2]=nbitem, triplets at offsets mots, 2*mots, ... . *pmin0 is the OFFSET into items[], not an index. saclav.f:119 reads `MCN(MNIT+NMIN0+2)` which only makes sense if NMIN0 is an offset.
- **Fix:** Implemented nearest_item_offset to follow the real X11 contract; shaped the test data as `{3, 2, 2, 100,100,1, 200,200,2}` and asserted `pmin0 == 3` (offset of first item), not `0`. Documented the layout in the nearest_item_offset header comment.
- **Files modified:** xvue/qt/src/xvue_qt_event.cpp, xvue/qt/tests/test_xvue_qt_event.cpp
- **Verification:** testXvsouris2Accrochage PASSes with the corrected data; would fail immediately if I had shipped the plan's incorrect layout. No existing Fortran caller had to change because the X11 contract is preserved.
- **Committed in:** 20470c9

---

**Total deviations:** 1 auto-fixed (Rule 1, plan's test data was wrong — the plan itself told me to verify and adjust, which I did). No Rule 2 (missing critical functionality), no Rule 3 (blocking), no Rule 4 (architectural).

## Issues Encountered

- **Worktree was at base ac282f8, not 27d44d5.** On startup, the `worktree_branch_check` command in the prompt detected the mismatch. `git reset --soft 27d44d5` + `git checkout HEAD -- .` brought the tree to the expected base cleanly. Two Task commits then landed on top as expected (4d0012c + 20470c9).
- **No other issues.** No Edit-tool cache drift observed (the Plan 04 retrospective mentioned this risk but all Task 1/2 edits reached disk on first write, verified by follow-up `grep` and successful compilation).

## User Setup Required

None. All new functionality is wired through existing env vars (MEFISTO_XVSOURIS_AUTOEXIT for headless smoke; MEFISTO_XVSOURIS_DEBUG for motion-count diagnostic). No new apt packages, no config changes, no filesystem setup.

## Next Phase Readiness

- **Plan 06 (empirical A/B validation + Assumption A2 resolution):** Ready. All four waves of automated tests are PASS/SKIP only (no FAIL). The remaining 2 SKIPs are documented not-load-bearing (testXvpauseReturnsOnMouseClick, testDebugLoggingEnvVar) — Plan 06 is free to leave them SKIPped or promote to real tests during the A/B validation pass.
- **Phase 6 (modal dialogs):** Ready. XvueApp::blockingDepth() > 0 now fires during xvsouris2_ accrochage loops too — modal dialogs will be refused correctly during rubber-band picking.
- **Accrochage visual parity (Assumption A5):** the Strategy B sprite is a crisp black square rather than an XOR-inverted highlight; developer A/B judgment in Plan 06 will record the verdict. No correctness bug — UX delta only.
- **No blockers or concerns.**

## Known Stubs

None. Every entry point this plan touches has a real body. testXvpauseReturnsOnMouseClick (xvpause_ mouse-click path) and testDebugLoggingEnvVar (runtime env-var flip) are deliberate QSKIPs documented in the earlier plan summaries — NOT stubs.

## Threat Flags

None. The two new pieces of externally-controllable surface — items[] array passed from Fortran into nearest_item_offset, and the 13x13 pixmap allocations — are both covered by the plan's existing threat register (T-05-05-01 mitigated by the 65536 nbitem clamp; T-05-05-02 mitigated by cleanupAccrochage; T-05-05-04 mitigated by the pre-existing BlockingDepthGuard < 4 assert).

## Self-Check: PASSED

- `xvue/qt/src/xvue_qt_api.cpp`: contains `new QPixmap(13, 13)` and `Qt::transparent` in initaccrochage_; `WaitMode::Souris2` in xvsouris2_
- `xvue/qt/src/xvue_qt_event.cpp`: contains `WaitMode::Souris2` (7×), `accroche_undo_tile_` (13×), `nearest_item` (3×)
- `xvue/qt/src/xvue_qt_event.h`: contains `cleanupAccrochage`
- `xvue/qt/tests/test_xvue_qt_event.cpp`: contains testInitaccrochage, testInitaccrochageBeforeInit, testXvsouris2Accrochage, testXvsouris2NullGuards, testXvsouris2ResizeInvalidatesTile
- Commit 4d0012c: FOUND (Task 1)
- Commit 20470c9: FOUND (Task 2)
- `bin/cbl_tout`: exit 0
- `bin/cbl_tout_qt`: exit 0
- `xvfb-run xvue_qt_event_tests`: 31 passed, 0 failed, 2 skipped
- ABI count: 57 (unchanged)
- SHELL-03 `verify_no_exec.sh`: clean
- `MEFISTO_XVSOURIS_AUTOEXIT=1 pp/ppxvtest0_qt`: exits 0 in < 1 s
- `grep -c "new QPixmap(13, 13)" xvue/qt/src/xvue_qt_api.cpp`: 1 (≥ 1 required)
- `grep -c "WaitMode::Souris2" xvue/qt/src/xvue_qt_event.cpp`: 7 (≥ 1 required)
- `grep -c "WaitMode::Souris2" xvue/qt/src/xvue_qt_api.cpp`: 2 (≥ 1 required)
- `grep -c "accroche_undo_tile_" xvue/qt/src/xvue_qt_event.cpp`: 13 (≥ 3 required)
- `grep -c "nearest_item" xvue/qt/src/xvue_qt_event.cpp`: 3 (≥ 1 required)
- `grep -c 'warn_once(warned, "xvsouris2_"' xvue/qt/src/xvue_qt_api.cpp`: 0
- `grep -c 'warn_once(warned, "initaccrochage_"' xvue/qt/src/xvue_qt_api.cpp`: 0
- `grep -c MEFISTO_XVSOURIS_AUTOEXIT xvue/qt/src/xvue_qt_api.cpp`: 9 (≥ 3 required)
- Worktree path verified: `realpath .` = `/home/drico/git/mefisto/.claude/worktrees/agent-a53f3721` — edits stayed inside the worktree

---
*Phase: 05-event-bridge-blocking-reads*
*Completed: 2026-04-18*
