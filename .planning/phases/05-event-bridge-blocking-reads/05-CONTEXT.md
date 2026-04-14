# Phase 5: Event bridge & blocking reads - Context

**Gathered:** 2026-04-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Replace the X11 blocking input primitives in `xvue/xvuelc.c` (`xvsouris_`, `xvsouris2_`, `xvpause_`, `deplsouris_` — lines 2183–2531) with a Qt implementation that:

- Blocks the Fortran caller on a nested local `QEventLoop` until a filtered mouse/keyboard event arrives
- Coalesces fast mouse-motion drags so `testa/pan2d` and `testa/torus` rubber-bands feel identical to X11
- Maintains a re-entrant `blockingDepth()` counter that Phase 6 will query to gate modal dialogs
- Preserves bit-for-bit Fortran ABI and menu key semantics (Esc, `@`, digits, letters)
- Preserves the `MEFISTO_XVSOURIS_AUTOEXIT` headless short-circuit already present in the current stubs

**Out of scope:** Modal dialogs, menus, toolbars (Phase 6). Async export and background I/O (Phase 7). Any new user-visible capabilities beyond parity with the X11 backend.

</domain>

<decisions>
## Implementation Decisions

### Event Loop & Bridge Structure
- **D-01:** `waitForEvent()` blocks via a **stack-local `QEventLoop`** (`loop.exec()` / `loop.quit()`). Each nested `xvsouris_` / `xvsouris2_` / `xvpause_` call constructs its own loop instance, giving natural re-entrancy without any global state.
- **D-02:** A dedicated `class XvueEventBridge : public QObject` owns the event filter. It is installed on `XvueCanvas` via `canvas->installEventFilter(bridge)` at canvas construction. The bridge owns: the active `QEventLoop*`, the last filtered event payload (type, coords, key code), the motion-coalescing deferred-quit timer, and the blocking-depth counter. Direct subclassing of `XvueCanvas` was explicitly rejected to keep Phase 6/7 event layering clean.
- **D-03:** `blockingDepth` counter location is **Claude's Discretion**. Default: `static int XvueApp::blockingDepth_` with static accessor `XvueApp::blockingDepth()`, incremented/decremented via a stack-scoped RAII guard at the top of `waitForEvent()`. Main-thread-only (SHELL-07) so no atomics required. Phase 6 gates modal dialogs with `if (XvueApp::blockingDepth() > 0) { /* refuse */ }`.

### Motion Coalescing
- **D-04:** Mouse-motion coalescing uses **`QTimer::singleShot(0, &loop, &QEventLoop::quit)`** — the "deferred quit" pattern. On each `QEvent::MouseMove` the bridge stashes `(x, y)` into its last-event slot and schedules a zero-delay single-shot quit. Any further motion events already queued get processed before the timer fires, each updating the stash. The timer ultimately quits the loop, which returns the *last* coordinate pair — faithful to X11's `XEventsQueued(QueuedAfterFlush)` semantics with zero added latency.
- **D-05:** Enable `Qt::AA_CompressHighFrequencyEvents` in `XvueApp::ensure()` before the `QApplication` is constructed, so Qt's native mouse-move compression is the primary deduplication mechanism. The deferred-quit pattern (D-04) provides the return-after-flush semantic on top of it — we do **not** hand-roll coalescing in parallel with Qt compression.

### Keyboard Mapping (KeySym → `nbc`)
- **D-06:** The Fortran menu system expects `nbc` as an integer ASCII code. **Esc (`27`) and `@` (`64`) must both map to the abort path** — parity with X11 (`xvuelc.c:~2260`), preserving French tutorial muscle memory and the "every workflow works unchanged" invariant from PROJECT.md.
- **D-07:** Keyboard translation mechanism is **Claude's Discretion**. Default: hybrid approach — try `QKeyEvent::text().toLatin1()` first (covers all printable ASCII plus Return/Tab/Esc when Qt populates `text()`), and fall through to a minimal `Qt::Key` switch for `Qt::Key_Escape → 27`, `Qt::Key_Return → 13`, `Qt::Key_Tab → 9` as a defensive net. Arrow keys, F-keys, modifiers return `nbc = 0` (unhandled / silently dropped) — matches X11 behavior where `XLookupString` returns length 0 for these.

### xvsouris2_ & deplsouris_ Scope
- **D-08:** **`xvsouris2_` is ported in full during Phase 5**, including `initaccrochage_` and the `mempxaccro` snap/rubber-band pixmap logic. `mempxaccro` reuses the `saved_canvas_` save/restore mechanism already delivered in Phase 4 — no new backing pixmap required, just a second managed `QPixmap` owned by the bridge (or by `XvueState`, planner's call). This delivers an interactive Qt mesher at end of Phase 5; no Phase 5.1 is needed.
- **D-09:** `deplsouris_` is a straightforward `QCursor::setPos(canvas->mapToGlobal(QPoint(*x, *y)))` port. Non-blocking, does not touch the event loop. **Wayland caveat:** `QCursor::setPos()` is a no-op on most Wayland compositors — document this in `xvue/qt/README.md` but do **not** attempt a workaround. X11 (directly or via XWayland) is the supported session type for Phase 5.

### Preserved Environment Hooks
- **D-10:** The `MEFISTO_XVSOURIS_AUTOEXIT` environment-variable short-circuit currently in the `xvsouris_` / `xvsouris2_` stubs (`xvue/qt/src/xvue_qt_api.cpp:~560`) must be retained byte-for-byte so the Phase 0 headless test harness and the Nyquist validation runners continue to work without opening any window.

### Claude's Discretion
- Exact `RAII` guard implementation for `blockingDepth` (helper struct, `std::scoped_exit`, or inline +/-).
- File layout inside `xvue/qt/src/` (`xvue_qt_event.{h,cpp}` vs merging into existing files).
- Whether `XvueEventBridge` is owned by `XvueApp`, `XvueWindow`, or the canvas itself (as long as its lifetime covers the canvas).
- Exact CMake target wiring for the new source file(s).

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### X11 Baseline (authoritative reference semantics)
- `xvue/xvuelc.c` §2183–2315 — `xvsouris_` body (motion coalescing, event filter, return codes)
- `xvue/xvuelc.c` §2317–2445 — `xvsouris2_` body (snap / `mempxaccro` / `initaccrochage_` / notypeevent=5)
- `xvue/xvuelc.c` §2482–2490 — `deplsouris_` body (`XWarpPointer`)
- `xvue/xvuelc.c` §2516–2531 — `xvpause_` body (blocking `XNextEvent` until KeyPress)

### Qt Backend State (to be extended)
- `xvue/qt/include/xvue_qt_api.h` — 57 `extern "C"` declarations; Phase 5 implements four of them
- `xvue/qt/src/xvue_qt_api.cpp` — current `warn_once` stubs + `MEFISTO_XVSOURIS_AUTOEXIT` short-circuit to preserve
- `xvue/qt/src/xvue_qt_app.h` — `XvueApp` singleton; add `blockingDepth_` static + accessor
- `xvue/qt/src/xvue_qt_canvas.h` — canvas widget the event filter is installed on
- `xvue/qt/src/xvue_qt_state.h` — runtime state struct; consider housing the snap `mempxaccro` `QPixmap` here alongside the existing `saved_canvas_`

### Cross-phase Constraints
- `.planning/phases/01-*/01-CONTEXT.md` — SHELL-03 (no `QApplication::exec()`), SHELL-07 (main-thread-only)
- `.planning/phases/02-*/02-CONTEXT.md` — single-backing-pixmap collapse, long-lived `QPainter`
- `.planning/phases/04-*/04-CONTEXT.md` — `saved_canvas_` save/restore (reusable for `mempxaccro`)

### Requirements
- `.planning/REQUIREMENTS.md` — EVENT-01 through EVENT-08
- `.planning/ROADMAP.md` — Phase 5 success criteria (5 items)
- `CLAUDE.md` — "Compilation must never break", "small `testa/` tests must continue to pass"

### Fortran Call Sites (read to verify ABI expectations)
- `mail/` — mesher interactive loops calling `xvsouris`, `xvsouris2`, `xvpause`
- `elas/`, `flui/`, `ther/` — solver interactive prompts calling `xvpause`
- `incl/*.inc` — any `COMMON` blocks touching event state

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **`QPixmap saved_canvas_` on `XvueState` (Phase 4):** Already implements the "snapshot backing, restore on release" pattern. The `mempxaccro` pixmap in `xvsouris2_` is structurally identical — planner should model it after `saved_canvas_`, not invent a new mechanism.
- **`XVUE_QT_ASSERT_MAIN_THREAD()` macro:** Install at the top of every new `extern "C"` entry point per SHELL-07.
- **`warn_once(warned, "name")` pattern:** Keep the pattern for any entry points intentionally left as stubs (there shouldn't be any in Phase 5 scope, but the helper is there if needed).
- **`MEFISTO_XVSOURIS_AUTOEXIT` short-circuit:** Lift the existing `getenv` guard out of the stubs into the new real implementation — don't rewrite it.

### Established Patterns
- **57 entry points stay byte-identical** between `xvue/xvuelc.c` and `xvue/qt/src/xvue_qt_api.cpp`. Phase 5 adds real bodies for four of them; the ABI signatures must not change.
- **No `QApplication::exec()` anywhere** (SHELL-03). CMake has a grep guard that fails the build if it appears — the nested `QEventLoop::exec()` is fine (different symbol, different semantics).
- **Fortran callers are imperative/blocking**, not event-driven. The whole purpose of the bridge is to preserve that imperative contract while the graphics toolkit underneath becomes event-driven.

### Integration Points
- **`XvueCanvas` construction site** — where the event filter is installed (in the canvas ctor, or in `XvueWindow` after creating the central widget).
- **`XvueApp::ensure()`** — where `Qt::AA_CompressHighFrequencyEvents` must be set *before* the `QApplication` constructor runs.
- **`xvue/qt/CMakeLists.txt`** — add new source file(s) for `xvue_qt_event.{h,cpp}` (or equivalent); ensure the grep guard for `QApplication::exec()` still passes.
- **`testa/pan2d`, `testa/torus`, and the other 3 baseline `testa/` cases** — end-to-end A/B validation targets per ROADMAP success criterion #5.

</code_context>

<specifics>
## Specific Ideas

- **"Feels indistinguishable from X11"** is the bar. Success criterion #5 is A/B parity on all 5 baseline `testa/` cases — the developer's eye is the final judge.
- **Mouse motion must not stutter on fast drags** — this is the specific failure mode the roadmap flags ("empirical mouse-motion coalescing validation needed at plan time"). The motion-coalescing decisions (D-04, D-05) exist exactly to address this.
- **French-language tutorial compatibility:** the `@` key alias for Esc is non-negotiable (D-06). Existing docs and muscle memory depend on it.
- **Planner should include an empirical validation step** where a developer runs `pp/ppmail_qt testa/pan2d`, does a rapid mouse drag, and confirms no stutter. This is not a unit test — it's a human check the roadmap explicitly calls for.

</specifics>

<deferred>
## Deferred Ideas

- **Wayland-native cursor warping** — deferred indefinitely. X11/XWayland is the supported session type. Documented as a caveat, not a task.
- **Modal-dialog gating via `blockingDepth()`** — Phase 5 delivers the counter; Phase 6 consumes it to refuse `QFileDialog::exec()` while a blocking read is active.
- **Async export / background I/O that must yield to events** — Phase 7 concern; Phase 5's `blockingDepth` counter will be sufficient infrastructure.
- **F-keys, arrow keys, modifier combinations** — currently return `nbc = 0` (matching X11 behavior). If the mesher ever needs to consume these, it's a future phase.
- **High-DPI / fractional scaling coordinate mapping** — assumed to be handled by Phase 2's canvas coordinate system; flag for verification during testa A/B but not a Phase 5 deliverable.

</deferred>

---

*Phase: 05-event-bridge-blocking-reads*
*Context gathered: 2026-04-14*
