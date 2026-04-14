# Phase 5: Event bridge & blocking reads — Research

**Researched:** 2026-04-14
**Domain:** Qt 6 nested event loops, installEventFilter, motion coalescing, Fortran-imperative ABI bridging
**Confidence:** HIGH on ABI / codebase; MEDIUM on Qt event-dispatch ordering (see Open Questions)

---

## User Constraints (from 05-CONTEXT.md)

### Locked Decisions (D-01..D-10)

- **D-01** — `waitForEvent()` blocks via a stack-local `QEventLoop`. Each nested call constructs its own loop instance; natural re-entrancy, no global state.
- **D-02** — dedicated `class XvueEventBridge : public QObject` installed via `canvas->installEventFilter(bridge)`. Owns: active `QEventLoop*`, last-event payload, coalescing timer, blocking-depth counter. Subclassing `XvueCanvas` was explicitly rejected.
- **D-03** — `blockingDepth` default: `static int XvueApp::blockingDepth_` + `XvueApp::blockingDepth()` accessor, RAII-guarded at top of `waitForEvent()`. Main-thread-only per SHELL-07. *(Claude's discretion — default accepted.)*
- **D-04** — Motion coalescing via `QTimer::singleShot(0, &loop, &QEventLoop::quit)`. On each `QEvent::MouseMove` stash `(x,y)` and schedule a zero-delay deferred quit. Loop returns last stashed coordinate.
- **D-05** — `Qt::AA_CompressHighFrequencyEvents` set in `XvueApp::ensure()` *before* `QApplication` ctor. Native compression is the primary dedup; D-04 is the return-after-flush semantic on top.
- **D-06** — Esc (27) **and** `@` (64) both map to abort. Parity with X11 `xvuelc.c:2308-2309` — non-negotiable for French tutorial muscle memory.
- **D-07** — Hybrid keymap: try `QKeyEvent::text().toLatin1()[0]` first; fall through to a minimal `Qt::Key` switch for `Escape→27`, `Return→13`, `Tab→9`. Arrow/F-keys/modifiers → `nbc=0`. *(Claude's discretion — default accepted.)*
- **D-08** — `xvsouris2_` ported **in full**, including `initaccrochage_` and `mempxaccro`. Reuses the Phase 4 `saved_canvas_` pattern (second managed `QPixmap`, same ownership style). No Phase 5.1.
- **D-09** — `deplsouris_` is `QCursor::setPos(canvas->mapToGlobal(QPoint(*x,*y)))`. Non-blocking. Wayland no-op caveat **documented**, no workaround.
- **D-10** — `MEFISTO_XVSOURIS_AUTOEXIT` environment short-circuit retained **byte-for-byte** from the current stubs.

### Claude's Discretion Areas

- Exact RAII guard implementation for `blockingDepth` (helper struct vs inline ±).
- File layout inside `xvue/qt/src/` (`xvue_qt_event.{h,cpp}` vs merge into existing files).
- Ownership of `XvueEventBridge` (`XvueApp`, `XvueWindow`, or canvas) — lifetime must cover the canvas.
- Exact CMake target wiring for new source file(s).

### Deferred (OUT OF SCOPE)

- Wayland-native cursor warping (deferred indefinitely; X11/XWayland only).
- Modal-dialog gating via `blockingDepth()` — Phase 5 delivers the counter; Phase 6 consumes it.
- Async export / background I/O — Phase 7.
- F-keys, arrow keys, modifier combos (currently return `nbc=0`; matches X11).
- High-DPI fractional scaling coord mapping — handled by Phase 2 canvas coord system.

---

## Phase Requirements

| ID | Behavior | Research Support |
|----|----------|------------------|
| **EVENT-01** | `XvueEventBridge` = QObject event filter on canvas, exposes `waitForEvent()` with local `QEventLoop` | §2 (Qt event loop), §4 (filter), §9 (file layout) |
| **EVENT-02** | `xvsouris_` returns `(notypeevent, nbc, x1, y1)` matching X11 `XNextEvent` semantics, no top-level `QApplication::exec()` | §5 (ABI contract), §2 (QEventLoop), Pitfall 1 |
| **EVENT-03** | `xvsouris2_` returns menu-item selections identically to X11 (`notypeevent=5` accrochage path) | §5, §6 (mempxaccro), CodeRef xvuelc.c:2317-2480 |
| **EVENT-04** | `xvpause_` blocks until any event arrives | §2, §8 (headless xvpause_ gap) |
| **EVENT-05** | `deplsouris_` returns/moves mouse position without blocking | D-09, Pitfall 5 (Wayland) |
| **EVENT-06** | `initaccrochage_` initializes snap/crosshair state on canvas | §6 (mempxaccro body at xvuelc.c:561-609) |
| **EVENT-07** | Mouse motion coalesced via deferred-quit — no stutter on fast drags in pan2d/torus | §3 (motion coalescing — flagged area) |
| **EVENT-08** | `XvueApp::blockingDepth()` re-entrancy counter for Phase 6 modal-dialog gating | D-03, §9 (XvueApp changes) |

---

## Project Constraints (from CLAUDE.md)

- **Compilation must never break.** Verify `bin/cbl_tout_qt` (and `bin/cbl_tout` for X11 parity) after every task. Full build green before phase closeout.
- **Small `testa/` tests must continue to pass.** Phase 5 end state must run `testa/pan2d`, `testa/torus`, plus the other 3 baseline cases on both backends.
- **Ask before installing system packages.** No new apt dependencies for Phase 5 (Qt 6 already covers it; no new headers).
- **Fortran ABI is frozen.** The four new entry points carry byte-identical names and signatures to `xvue/xvuelc.c` — verify with the existing `verify_abi.sh` guard.
- **Coding norms in `doc/normes.ps`** — C++ side is not covered there; the `xvue/qt/` convention is C++17, Qt style.
- **Git discipline:** commit after each working sub-step.

---

## 1. Summary

Phase 5 wires up a `QObject` event filter (`XvueEventBridge`) on `XvueCanvas` and four new bodies in `xvue_qt_api.cpp` (`xvsouris_`, `xvsouris2_`, `xvpause_`, `deplsouris_`) that call a single shared `waitForEvent()` helper. Each helper runs a stack-local `QEventLoop` and an RAII `blockingDepth` guard. Mouse-motion events are coalesced via `QTimer::singleShot(0, &loop, &QEventLoop::quit)` layered on top of Qt's built-in `AA_CompressHighFrequencyEvents` (default-true on X11). `xvsouris2_` is ported **in full** including `initaccrochage_` and a `mempxaccro` `QPixmap` that reuses the Phase 4 `saved_canvas_` ownership pattern. The `MEFISTO_XVSOURIS_AUTOEXIT` short-circuit in the existing stubs is lifted verbatim, and the Phase 5 plan extends it to `xvpause_` as well (§8 — gap identified).

**Primary recommendation:** implement in five tasks — (1) `XvueEventBridge` header + empty cpp + RAII guard; (2) install-filter wiring in `XvueWindow`/`XvueApp::ensure()` with `AA_CompressHighFrequencyEvents`; (3) `xvsouris_` + `xvpause_` bodies using `waitForEvent()`; (4) `xvsouris2_` + `initaccrochage_` + `mempxaccro` pixmap with GXand/GXorInverted emulated via `QPainter::CompositionMode_Difference` or two-pixmap blit; (5) `deplsouris_` + A/B empirical motion test on `testa/pan2d` + `testa/torus`.

---

## 2. Qt Event Loop Mechanics (deep dive)

### `QEventLoop::exec()` — authoritative API surface

From Qt 6 docs (https://doc.qt.io/qt-6/qeventloop.html) — `int QEventLoop::exec(ProcessEventsFlags flags = AllEvents)`:

> "Enters the main event loop and waits until exit() is called."

Supported `ProcessEventsFlags`:
- `AllEvents` — process everything
- `ExcludeUserInputEvents` — skip mouse/keyboard/tablet/touch
- `ExcludeSocketNotifiers`
- `WaitForMoreEvents` — block idle rather than return

**For Phase 5 we want `AllEvents`** — we *want* the user input that's the whole point. `ExcludeUserInputEvents` is wrong here. [CITED: doc.qt.io/qt-6/qeventloop.html]

### Re-entrancy

`QEventLoop` instances are cheap and independent. Nested `exec()` calls are the supported idiom for modal dialogs (`QDialog::exec()` is implemented on top of `QEventLoop::exec()`). Each nested loop processes events from the same thread-local queue in FIFO order; outer loops resume on inner quit. [CITED: doc.qt.io/qt-6/qeventloop.html — "Event loops are useful whenever you need to receive events … QDialog provides this functionality in a single-function call."]

### quit() before exec() race

Qt doc does not explicitly cover "quit() called before exec()." Empirical Qt behavior (Qt source `qeventloop.cpp::exec()`): `exec()` checks `d->exit.loadAcquire()` at the top of each iteration, so a pre-set `exit` flag causes `exec()` to return immediately without processing events. **This is actually desirable for Phase 5** — it means the deferred-quit timer cannot "fire too early and be lost" as long as we follow this invariant: construct the loop → arm the filter/timer → call `exec()`. [ASSUMED: based on Qt source-reading convention; needs verification if the planner wants HIGH confidence. Flag as empirical-risk.]

### Interaction with `QCoreApplication::processEvents`

`processEvents()` is fine to call from *outside* a nested `exec()`; the `xvvoir_` / `xvfermer_` teardown already does this (`xvue_qt_api.cpp:984`). Inside a nested `exec()` we should **not** call `processEvents()` from the filter — the loop is already pumping events for us.

### Interaction with SHELL-03 grep guard

The existing regex in `xvue/qt/cmake/verify_no_exec.sh:23`:
```
grep -R -n -E 'QApplication::exec|qApp->exec\(\)'
```
This matches **only** `QApplication::exec` and `qApp->exec()`. It does **not** match `loop.exec()`, `QEventLoop::exec()`, `QDialog::exec()`, or `bridge->waitForEvent()`. Phase 5 passes the guard unchanged. [VERIFIED: xvue/qt/cmake/verify_no_exec.sh:23]

---

## 3. Motion Coalescing Mechanics (empirical — this is the flagged area)

### X11 baseline (reference semantics to preserve)

From `xvue/xvuelc.c:2248-2263`:
```c
if (event.type == MotionNotify) {
    if (XEventsQueued(display_mef, QueuedAfterFlush) <= 0) {
        *notypeevent = -2; *x1 = event.xbutton.x; *y1 = event.xbutton.y;
        ... flag = 1;
    }
    // else: drop this motion event, loop back to XNextEvent
}
```
**Semantic:** "drop every motion event that has another event queued behind it; only return when we see the *last* motion in the burst." `QueuedAfterFlush` flushes pending output first so the server-side queue is a true snapshot. This is NOT "coalesce", it is "emit only the tail of the burst with zero added latency."

### Qt `AA_CompressHighFrequencyEvents`

From Qt 6 docs (https://doc.qt.io/qt-6/qt.html — `Qt::ApplicationAttribute`):

> "Enables compression of certain frequent events. On the X11 windowing system, the default value is true, which means that QEvent::MouseMove, QEvent::TouchUpdate, and changes in window size and position will be combined whenever they occur more frequently than the application handles them, so that they don't accumulate and overwhelm the application later."

[CITED: doc.qt.io/qt-6/qt.html]

**Critical nuances:**
1. **Default-true on X11 already.** Setting it explicitly in `XvueApp::ensure()` is defensive but not strictly required on our target platform. [VERIFIED against Qt 6 docs]
2. **Compression happens at dispatch time, not post time.** Qt merges consecutive `MouseMove` events queued for the same target into a single event. [CITED: doc.qt.io/qt-6/qt.html — "combined whenever they occur more frequently than the application handles them"]
3. **Must-be-set-before-QApplication claim in D-05 is NOT confirmed by Qt docs.** QCoreApplication docs say *"Some application attributes must be set before creating a QCoreApplication instance. Refer to the Qt::ApplicationAttribute documentation for more information."* The `Qt::ApplicationAttribute` page does NOT list `AA_CompressHighFrequencyEvents` in its before-construction list. [CITED: doc.qt.io/qt-6/qcoreapplication.html#setAttribute] **→ D-05's "BEFORE QApplication ctor" requirement is [ASSUMED]; safe to set, but setting it after is also likely fine. See Assumptions Log A1.**

### `QTimer::singleShot(0, &loop, &QEventLoop::quit)` — deferred quit pattern

**What we know from Qt docs** (not cited explicitly in the QTimer page for this idiom):
- A 0-ms `singleShot` enqueues a `QTimerEvent` (or invokes via `QMetaObject::invokeMethod` with `Qt::QueuedConnection` internally) at the **tail** of the thread's event queue.
- When `loop.exec()` next dispatches events, it processes the queue in FIFO order until the timer event fires `loop.quit()`.
- Any `MouseMove` events already queued **ahead** of the deferred-quit timer get dispatched first — the event filter stashes each new `(x,y)` without modifying the pending quit.
- Fresh `MouseMove` events that arrive *after* the timer is armed but *before* the queue has drained: these go to the tail too, *behind* the deferred-quit timer, so they get processed next iteration (if we re-arm).

**The correct implementation pattern** (to match X11 `QueuedAfterFlush`):

```cpp
bool XvueEventBridge::eventFilter(QObject* watched, QEvent* ev) {
    if (!loop_) return false;  // not in a waitForEvent — pass through

    switch (ev->type()) {
    case QEvent::MouseMove: {
        auto* me = static_cast<QMouseEvent*>(ev);
        pending_x_ = me->pos().x();
        pending_y_ = me->pos().y();
        pending_notypeevent_ = -2;
        pending_nbc_ = 0;  // X11 contract: motion carries no button
        if (!quit_pending_) {
            quit_pending_ = true;
            QTimer::singleShot(0, loop_, &QEventLoop::quit);
        }
        return true;  // eat the event — we've captured what we need
    }
    case QEvent::MouseButtonPress: /* decisive, quit immediately */ ...
    case QEvent::MouseButtonRelease: ...
    case QEvent::KeyPress: ...
    }
    return false;
}
```

**Key insight:** because `AA_CompressHighFrequencyEvents` is already coalescing motion events at Qt's level on X11, the deferred-quit timer is a **safety net** that handles the edge case where a burst of moves arrives *before* compression kicks in. In most real sessions the bridge will see one motion event per ~16ms vsync tick (already compressed), and the timer fires after that single move is stashed, and `exec()` returns.

**Empirical validation required** (see §10 Validation Architecture): record motion-event counts during a fast drag and compare X11 vs Qt. Target: within 2× of each other (X11 backend typically emits ~30–60 motions/sec on a pan drag; Qt should be in the same ballpark, not ~500/sec).

---

## 4. Event Filter Ordering & Return Semantics

From Qt 6 docs (https://doc.qt.io/qt-6/eventsandfilters.html):

> "An event filter gets to process events before the target object does."
> "If one of them stops processing (by returning true), the target and any later event filters do not get to see the event at all."
> "global event filters are called before the object-specific filters."

[CITED: doc.qt.io/qt-6/eventsandfilters.html]

### Phase 5 return-value strategy

| Event the bridge sees | Return value | Rationale |
|-----------------------|--------------|-----------|
| `MouseMove` / `MouseButtonPress` / `MouseButtonRelease` / `KeyPress` **while `loop_ != nullptr`** | **`true`** (eat) | The Fortran caller owns the blocking read; we do not want Qt's default handler or the canvas widget to see the event. Eating prevents double-handling and prevents focus-change side effects. |
| Any event **while `loop_ == nullptr`** (not inside `waitForEvent()`) | **`false`** (pass-through) | Phase 6 menu actions, canvas repaint, resize all flow unchanged. |
| `Close` / `Resize` / `Paint` | **`false`** always | These are never consumed by the bridge. |

### Filter install order

`bridge` is the **only** object-specific filter on `XvueCanvas` in Phase 5. Phase 6 may add a global filter (for menu-hotkey interception). The ordering doc guarantees that global filters run first — if a Phase 6 global filter eats an event before reaching the canvas filter, the Phase 5 bridge simply never sees it. This is acceptable and designed-in.

### Compression-vs-filter interaction (critical)

Does `AA_CompressHighFrequencyEvents` compress events **before** or **after** the event filter? Qt docs are silent on this specific ordering.

**[ASSUMED based on reading Qt source convention]:** Compression happens inside `QApplication::compressEvent()` (Qt source `qguiapplication.cpp`), which is called from the platform plugin's `postEvent()` path — i.e., *before* the event reaches `QCoreApplication::notify()`, which is *before* the object-specific event filters run. **Conclusion:** the bridge sees *already-compressed* motion events. This is the desired behavior. **Flag as A2 in Assumptions Log** — a 5-minute runtime printf in Task 5 can confirm (count `MouseMove` events observed by the filter during a 1-second drag and compare to X11 `XEventsQueued` counts on the same drag).

---

## 5. Fortran ABI Contract Verification

### `xvsouris_` — verified signature and event codes

**Header declaration** (`xvue/qt/include/xvue_qt_api.h`, byte-identical to `xvuelc.c:2183`):
```c
void proc(xvsouris)(int *notypeevent, int *nbc, int *x1, int *y1);
```

**Call sites verified** (example: `xvue/zoom2d3.f:28,65`, `util/searclic.f:68`, `prpr/xvtest1.f:156+`):
```fortran
CALL XVSOURIS( NOTYEV, NBC, NOPX, NOPY )
IF( NOTYEV .EQ. 0 .OR. NOTYEV .EQ. 2 ) ... abandon
IF( NOTYEV .EQ. -1 ) ... button pressed (NBC = button 1/2/3)
IF( NOTYEV .EQ. 1 ) ... button released
```

**Event code table (from `xvue/xvuelc.c:2183-2315` — the canonical source — and the Fortran doc comments in `util/searclic.f:70-75`):**

| `*notypeevent` | Meaning | `*nbc` semantic |
|---|---|---|
| **0** | ABANDON — Esc (27) or `@` (64) key, OR (historically) mouse button 2; see D-06 | ASCII of the abort key (27 or 64) |
| **1** | Mouse button pressed AND released (full click) | button number 1 / 2 / 3 |
| **-1** | Mouse button pressed only, not yet released | button number 1 / 2 / 3 |
| **-2** | Mouse motion, no button pressed | **0** (contract: "pas de bouton designe", `xvuelc.c:2198`) |
| **2** | Keyboard character typed | ASCII of the character; if 27 or 64, caller also sets notypeevent=0 |

*Wait* — the X11 code at `xvuelc.c:2257-2260` actually sets `*nbc` to the button number on a motion event, even though the doc comment at `xvuelc.c:2198` says `nbc=0` for motion. **The doc-comment contract is `nbc=0` for motion; the X11 code does something slightly different but in practice the motion path rarely carries button state the Fortran callers rely on.** The Qt bridge should follow the **doc-comment contract** (`nbc=0` on motion) for clean semantics, since no verified Fortran caller actually reads `nbc` after `notypeevent=-2`. [VERIFIED against doc comment and all grep-matched call sites; if planner finds a counter-example, revisit.]

### `xvsouris2_` — verified signature

**Header** (`xvue/qt/include/xvue_qt_api.h`):
```c
void proc(xvsouris2)(int *items, int *pmin0, int *notypeevent, int *ibutton, int *x1, int *y1);
```

**Call site** (`xvue/saclav.f:61`):
```fortran
CALL XVSOURIS2( MCN(MNIT), NMIN0, NOCODE, NASCII, NX, NY )
```

**Event codes** (from `xvue/xvuelc.c:2336-2342`, Fortran doc `saclav.f:64-72`):

| `*notypeevent` | Meaning |
|---|---|
| **0** | No event (error path — the Fortran caller retries after `XVPAUSE`) |
| **1** | Full mouse click (button pressed and released) — `ibutton` = 1/2/3 |
| **2** | Keyboard character — `ibutton` = ASCII |
| **5** | Button pressed **with or without motion** during accrochage (the snap path) — `ibutton` = 1/2/3, `pmin0` updated to nearest item index |

**`pmin0` contract:** in/out — `-2` means "not initialized, no item currently highlighted"; `>= 0` means "index of currently-highlighted item in `items[]`." The bridge must blit-erase the old highlight at `items[*pmin0..*pmin0+1]` before drawing the new one.

### `xvpause_` — verified signature

```c
void proc(xvpause)(void);  // no args
```

Call sites: 30+ instances in `mail/*.f`, `ther/*.f`, `elas/*.f`, `flui/*.f`. All call `CALL XVPAUSE` without any parenthesized argument list (Fortran 77 allows this). **Contract: block until any keypress, then return.** (`xvuelc.c:2516-2531`.)

### `deplsouris_` — verified signature

```c
void proc(deplsouris)(int *x, int *y);
```

Single call site found: `mail/trfasevo.f:202` — `CALL DEPLSOURIS( NOPX0, NOPY0 )`. Contract: warp the pointer to `(x,y)` in canvas coordinates. Non-blocking.

### `initaccrochage_` — verified signature

```c
void proc(initaccrochage)(void);
```

Body at `xvue/xvuelc.c:561-609`. Contract: **create the `mempxaccro` pixmap** (13×13 pixels, `lmempxaccro = hmempxaccro = 13` — `xvuelc.c:142-143`), fill it white, draw a 3-pixel-thick black square border inside, reset line width to 1. This is called **once** at the start of the interactive mesher, before any `xvsouris2_` call.

[VERIFIED: all four entry points against xvuelc.c source and Fortran grep.]

---

## 6. `mempxaccro` / `initaccrochage_` Semantics

### What `mempxaccro` actually is

A **13×13 pixmap** (`lmempxaccro = hmempxaccro = 13`) that contains a pre-drawn black square outline on white background. It's **not** a rubber-band line and **not** a snap indicator stored at a nearest-vertex location — it's a **reusable highlight sprite** that gets blitted onto the canvas wherever `xvsouris2_` decides an item needs to be visually marked as "currently snapped to." [VERIFIED: `xvuelc.c:561-609,1039,2415-2450`]

### How `xvsouris2_` uses it

From `xvuelc.c:2415-2450`:

1. On every `MotionNotify` or `ButtonPress` inside `xvsouris2_`:
2. Compute the nearest item to `(x1, y1)` → index `pmin`.
3. **If the previously-highlighted item is no longer the nearest:**
   - `XSetFunction(GXorInverted)` + `XCopyArea(mempxaccro → fenetre_mef)` at the **old** item's coordinates → this XOR-blits the sprite back over the old spot, **erasing** the highlight (because GXorInverted is self-inverse against GXand).
4. **If a new nearest item exists:**
   - `XSetFunction(GXand)` + `XCopyArea(mempxaccro → fenetre_mef)` at the **new** item's coordinates → this AND-blits the sprite onto the new spot, drawing a dark square highlight on top of the existing canvas content.
5. On `ButtonRelease`: erase one last time (step 3) and return `notypeevent=1`.

**The GXand / GXorInverted pair is an XOR-pair trick**: AND-blit to draw, OR-inverted-blit to erase. The canvas underneath is preserved because the sprite's white pixels AND with anything = the underlying pixel (no change outside the black border), and the black border pixels AND to zero (drawing black). The erase reverses the operation.

### Porting to Qt — two viable strategies

**Strategy A: composition-mode emulation (RISKY — raster ops ≠ Qt composition modes exactly).**
`QPainter::setCompositionMode(QPainter::CompositionMode_Difference)` or `CompositionMode_Xor` approximates XOR but operates on RGB values, not X11 plane-mask raster ops. **Visual result will differ from X11** — flag as validation risk.

**Strategy B (RECOMMENDED): save-restore two-pixmap blit, mirror the Phase 4 `saved_canvas_` pattern.**
Add a `QPixmap* mempxaccro_` field on `XvueState` next to `saved_canvas_`. `initaccrochage_` populates it with a 13×13 black square on transparent background (using `QPixmap::fill(Qt::transparent)` + `QPainter::drawRect(1,1,10,10)`). `xvsouris2_` uses the same "save tile under sprite → draw sprite → on next move, restore saved tile → save new tile → draw new sprite" pattern the Phase 4 rubber-band save/restore already uses.

This is strictly better than the X11 XOR trick because:
- It works regardless of Qt composition mode subtleties.
- It handles arbitrary underlying colors (X11's XOR only worked because the palette was indexed).
- It reuses the exact `saved_canvas_` lifecycle — no new moving parts.

**The `XvueState` addition** is a sibling field: `QPixmap* mempxaccro_ = nullptr;` initialized by `initaccrochage_`, plus `QPixmap* accroche_undo_tile_ = nullptr;` to hold the 13×13 bytes of canvas content currently under the highlight. Lifetime identical to `saved_canvas_` (`xvue_qt_state.{h,cpp}` `~XvueState` already deletes the saved pixmap — just mirror those two lines).

[VERIFIED against `xvuelc.c:561-609,2415-2450`; VERIFIED against `xvue/qt/src/xvue_qt_state.h:38`.]

---

## 7. Keyboard Mapping Edge Cases

### `QKeyEvent::text()` behavior

Qt 6 docs (https://doc.qt.io/qt-6/qkeyevent.html#text):
> "Returns the Unicode text that this key generated."

For a French AZERTY keyboard pressing AltGr+0 to produce `@`: Qt's platform plugin (xcb on X11) receives the already-composed Unicode `@` character from the input method layer, and `text()` returns `"@"`. `text().toLatin1()[0]` returns `0x40 = 64`. [ASSUMED based on Qt xcb plugin behavior; needs empirical verification on a live AZERTY layout — see A3.] **Mitigation:** the D-07 hybrid also checks `event->key() == Qt::Key_At` in the fallback switch, so even if `text()` is empty on a weird layout, `@` still maps correctly. **Defensive add to D-07 switch:** `Qt::Key_At → 64`.

### Edge cases for D-07 fallback switch

| Qt::Key | Fallback `nbc` | Rationale |
|---------|---------------|-----------|
| `Qt::Key_Escape` | 27 | X11 parity |
| `Qt::Key_Return` / `Qt::Key_Enter` | 13 | X11 `XLookupString` produces CR |
| `Qt::Key_Tab` | 9 | X11 parity |
| `Qt::Key_At` | 64 | D-06 abort, belt-and-braces |
| `Qt::Key_Backspace` | 8 | Used in `saclav.f:286` (`NASCII .EQ. 8`) |
| Arrow keys, F-keys, modifiers alone | 0 | `notypeevent` stays at "no event", loop continues |

### Testing without a live QApplication

`QTest::keyClick(widget, Qt::Key_Escape)` requires a `QApplication` and a real widget, but it does **not** require a visible window — `offscreen` QPA platform supports it. Unit tests for the bridge can use `QTest::keyClick(canvas, Qt::Key_At)` + assert `waitForEvent()` returns `notypeevent=0, nbc=64`. This is the recommended test pattern for EVENT-02, EVENT-03, EVENT-04.

[CITED: doc.qt.io/qt-6/qtest.html#keyClick]

---

## 8. Headless Test Harness (existing + extensions needed)

### Existing infrastructure (Phase 3 delivery)

`bin/xvtest-capture.sh:93-96`:
```bash
export MEFISTO_XVSOURIS_AUTOEXIT=1
export MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=100
export MEFISTO_XVFERMER_READY_FILE="$READY_FILE"
export MEFISTO_XVFERMER_HOLD_MS="$DELAY_MS"
```
Plus `bin/qt-capture.sh`, `bin/testa-capture.sh` all set the same env vars.

Current stub behavior (`xvue/qt/src/xvue_qt_api.cpp:904-933`): on AUTOEXIT, flush a pump loop for `delay_ms`, then synthesize `notypeevent=2, nbc=' '` and return.

### Gap for `xvpause_`

**Neither `xvpause_` (stub at `xvue_qt_api.cpp:988-993`) nor the X11 `xvpause_` (xvuelc.c:2516-2531) honors `MEFISTO_XVSOURIS_AUTOEXIT`.** Headless runs currently hit `CALL XVPAUSE` and hang forever (or rather: the stub is a no-op so it returns immediately in the Qt backend today; the X11 backend hangs — see `bin/xvtest-capture.sh` comment flow which relies on drivers that don't use `XVPAUSE`).

**Phase 5 must preserve or extend this:** the new real `xvpause_` body WILL block on a nested `QEventLoop` → a headless run hits it and hangs forever. Plan required:

1. **Extend the AUTOEXIT semantics** to cover `xvpause_`. Same env var, same delay, same return (immediate no-op short-circuit after `delay_ms` flush).
2. Also extend the **X11 side** `xvuelc.c::xvpause_` with the identical `getenv` guard, so both backends stay in parity (follow the Phase 3 pattern that added AUTOEXIT to xvuelc.c xvsouris lines 2218-2243).
3. Update `bin/xvtest-capture.sh` header comments to note `xvpause_` is now short-circuited too.

**Alternatively:** add a new env var `MEFISTO_XVPAUSE_AUTOEXIT` for finer control. **Recommendation: reuse `MEFISTO_XVSOURIS_AUTOEXIT`** — same test harness, same semantics, one fewer knob. The user's intent ("headless") is identical for both.

### Shell harness

No `ctest` / `make`-driven test runner. The pattern is:
```
bin/cbl_tout_qt && xvfb-run -a bin/qt-capture.sh pp/ppxvtest1_qt /tmp/out.png
```
Qt tests should follow the same pattern: `xvfb-run -a pp/ppmail_qt testa/pan2d` under `MEFISTO_XVSOURIS_AUTOEXIT=1` produces a deterministic exit and a capturable screenshot.

---

## 9. File Layout Proposal

### New files

| File | Purpose | LOC estimate |
|------|---------|--------------|
| `xvue/qt/src/xvue_qt_event.h` | `class XvueEventBridge : public QObject` declaration + `struct BlockingDepthGuard` RAII helper | ~60 |
| `xvue/qt/src/xvue_qt_event.cpp` | `eventFilter()` implementation, `waitForEvent()` helper, key→nbc translation | ~180 |

**Rationale for a separate TU:** keeps `xvue_qt_api.cpp` at ~1200 lines (already large), and isolates the `QObject` / `Q_OBJECT` macro behind CMake AUTOMOC (one new `moc_xvue_qt_event.cpp` generated). Follows the same separation as `xvue_qt_canvas.{h,cpp}`.

### Modifications to existing files

| File | Change |
|------|--------|
| `xvue/qt/src/xvue_qt_app.h` | Add `static int blockingDepth_;` + `static int blockingDepth();` + friend `BlockingDepthGuard` |
| `xvue/qt/src/xvue_qt_app.cpp` | Add `int XvueApp::blockingDepth_ = 0;` static init, accessor body |
| `xvue/qt/src/xvue_qt_app.cpp::ensure()` | Add `QCoreApplication::setAttribute(Qt::AA_CompressHighFrequencyEvents);` **before** `qapp_ = std::make_unique<QApplication>(...)` (defensive; Qt docs say it's default-true on X11 but setting explicitly documents intent) |
| `xvue/qt/src/xvue_qt_window.cpp` (or wherever `XvueCanvas` is constructed) | Instantiate `XvueEventBridge* bridge_` owned by `XvueWindow`, call `canvas->installEventFilter(bridge_)` after canvas construction |
| `xvue/qt/src/xvue_qt_state.h` | Add `QPixmap* mempxaccro_ = nullptr;` and `QPixmap* accroche_undo_tile_ = nullptr;` (§6 Strategy B); update `~XvueState` to delete both |
| `xvue/qt/src/xvue_qt_api.cpp::xvsouris_` | Replace lines 904-933 stub with real body: AUTOEXIT short-circuit (keep verbatim) → `bridge->waitForEvent(WaitMode::Souris, out_params)` |
| `xvue/qt/src/xvue_qt_api.cpp::xvsouris2_` | Replace lines 936-965 stub with real body including accrochage pixmap blit |
| `xvue/qt/src/xvue_qt_api.cpp::deplsouris_` | Replace lines 968-974 stub with `QCursor::setPos(canvas->mapToGlobal(QPoint(*x, *y)))` + `XVUE_QT_ASSERT_MAIN_THREAD()` |
| `xvue/qt/src/xvue_qt_api.cpp::xvpause_` | Replace lines 988-993 stub with AUTOEXIT short-circuit + `bridge->waitForEvent(WaitMode::Pause, out_params)` |
| `xvue/qt/src/xvue_qt_api.cpp::initaccrochage_` | Replace lines 281-286 stub with real body: allocate `state->mempxaccro_` 13×13 QPixmap, draw the square sprite |
| `xvue/xvuelc.c::xvpause_` | Add AUTOEXIT short-circuit for X11-side headless parity (mirror §8 gap fix) |
| `xvue/qt/CMakeLists.txt` | Add `src/xvue_qt_event.cpp` to the `add_library(xvueqt STATIC ...)` list (line 19-25); AUTOMOC picks up the new Q_OBJECT automatically |

### `XvueEventBridge` proposed header (sketch)

```cpp
// xvue/qt/src/xvue_qt_event.h
#pragma once
#include <QObject>

class QEventLoop;
class XvueCanvas;

class XvueEventBridge : public QObject {
    Q_OBJECT  // BUILD-03 / Phase 1 Pitfall 8
public:
    enum class WaitMode { Souris, Souris2, Pause };

    struct Result {
        int notypeevent = 0;
        int nbc         = 0;  // or ibutton, same slot
        int x           = 0;
        int y           = 0;
    };

    explicit XvueEventBridge(XvueCanvas* canvas, QObject* parent = nullptr);
    ~XvueEventBridge() override;

    // Blocks on a stack-local QEventLoop until a matching event arrives.
    // Increments/decrements XvueApp::blockingDepth_ via RAII.
    // Thread-affinity: main thread only (SHELL-07 asserted at entry).
    Result waitForEvent(WaitMode mode,
                        int* items = nullptr,
                        int* pmin0 = nullptr);

protected:
    bool eventFilter(QObject* watched, QEvent* event) override;

private:
    XvueCanvas* canvas_;
    QEventLoop* loop_ = nullptr;        // non-null only while waitForEvent is active
    WaitMode    mode_ = WaitMode::Souris;
    Result      pending_;
    bool        quit_pending_ = false;  // deferred-quit timer already armed?

    // For Souris2 mode only
    int*        items_ = nullptr;
    int*        pmin0_ = nullptr;

    static int  translateKey(class QKeyEvent* ev);
};
```

---

## 10. Validation Architecture (Nyquist)

### Test Framework

| Property | Value |
|----------|-------|
| Framework | No C++ unit-test framework wired today (`xvue/qt/tests/` does not exist). Phase 5 introduces the first. |
| Config file | **Wave 0 task: create** `xvue/qt/tests/CMakeLists.txt` using Qt's `QTest` (ships with `qt6-base-dev`, no extra apt). |
| Quick run command | `cmake --build xvue/build --target xvue_qt_event_tests && xvue/build/xvue_qt_event_tests` |
| Full suite command | `bin/cbl_tout_qt && xvfb-run -a xvue/build/xvue_qt_event_tests && xvfb-run -a bin/qt-capture.sh pp/ppmail_qt testa/pan2d /tmp/pan2d_qt.png` |

### Phase Requirements → Test Map

| Req | Behavior | Test Type | Automated Command | File Exists? |
|-----|----------|-----------|-------------------|--------------|
| EVENT-01 | Bridge exists, is QObject, installed on canvas | unit | `xvue_qt_event_tests::testBridgeInstallation` — create XvueApp, open canvas, assert `bridge_ != nullptr` and installed as filter | ❌ Wave 0 |
| EVENT-02 | `xvsouris_` returns correct `notypeevent/nbc/x/y` for each event type | unit | `xvue_qt_event_tests::testXvsouris_*` — `QTest::mouseClick`, `QTest::mouseMove`, `QTest::keyClick` drive the filter and assert return values | ❌ Wave 0 |
| EVENT-02 abort-key parity | Esc and @ both map to `notypeevent=0` | unit | `testXvsourisEscapeAbort`, `testXvsourisAtSignAbort` with `QTest::keyClick(Qt::Key_Escape)` and `QTest::keyClick(Qt::Key_At)` | ❌ Wave 0 |
| EVENT-03 | `xvsouris2_` returns `notypeevent=5` with pmin0 updated during accrochage | unit | `testXvsouris2Accrochage` — populate `items[]` with 2 points, `QTest::mouseMove` over the area, assert `pmin0` flips | ❌ Wave 0 |
| EVENT-04 | `xvpause_` returns on any key | unit | `testXvpauseReturnsOnKey` with `QTest::keyClick(Qt::Key_Space)` | ❌ Wave 0 |
| EVENT-05 | `deplsouris_` moves cursor without blocking | unit | `testDeplsourisNonBlocking` — call `deplsouris_(&x, &y)`, assert `QCursor::pos()` matches (X11/xvfb only) | ❌ Wave 0 |
| EVENT-06 | `initaccrochage_` allocates 13×13 `mempxaccro_` | unit | `testInitaccrochage` — assert `state->mempxaccro_->size() == QSize(13,13)` | ❌ Wave 0 |
| EVENT-07 | Motion coalescing produces ≤ N events per real drag | **integration** | `testMotionCoalescing` — `QTest::mouseMove` in a loop to 100 positions, assert filter sees ≤ 20 (tunable) motion events returned from `waitForEvent` | ❌ Wave 0 |
| EVENT-07 empirical | Real pp/ppmail_qt testa/pan2d drag does not stutter (subjective) | **human A/B** | Run `pp/ppmail testa/pan2d` and `pp/ppmail_qt testa/pan2d`, drag rapidly in rubber-band zoom, developer eye judges | **manual, logged to** `.planning/phases/05-*/05-VALIDATION.md` |
| EVENT-08 | `blockingDepth()` increments/decrements correctly, survives exception | unit | `testBlockingDepthRAII` — nested `waitForEvent`, throw from within, assert depth back to 0 | ❌ Wave 0 |

### Sampling Rate

- **Per task commit:** `xvue_qt_event_tests` (<2 seconds) + `bin/cbl_tout_qt` smoke build
- **Per wave merge:** full test binary + `bin/xvtest-capture.sh pp/ppxvtest{1..4}_qt` + `bin/qt-capture.sh pp/ppmail_qt testa/pan2d` headless AUTOEXIT run
- **Phase gate:** full suite green + **manual A/B drag test on `testa/pan2d` and `testa/torus`** logged to VALIDATION.md

### Wave 0 Gaps

- [ ] `xvue/qt/tests/CMakeLists.txt` — QTest target, linked against xvueqt + Qt6::Test
- [ ] `xvue/qt/tests/test_xvue_qt_event.cpp` — covers EVENT-01..EVENT-08
- [ ] `bin/cbl_tout_qt` extension to build the test target (or a sibling `bin/cbl_test_qt`)
- [ ] `.planning/phases/05-*/05-VALIDATION.md` template with A/B drag-test checklist
- [ ] `bin/xvtest-capture.sh` comment update — note `xvpause_` is now covered by AUTOEXIT
- [ ] AUTOEXIT extension to `xvuelc.c::xvpause_` and `xvue_qt_api.cpp::xvpause_`

---

## 11. Pitfalls

1. **SHELL-03 regex scope — safe as written.** `verify_no_exec.sh:23` matches `QApplication::exec|qApp->exec\(\)` **only**. `loop.exec()`, `QEventLoop::exec()`, `bridge->waitForEvent()` all pass. **DO NOT** loosen the regex to a bare `exec\(\)` — that would break Phase 5. [VERIFIED xvue/qt/cmake/verify_no_exec.sh:23]

2. **Re-entrancy stack overflow if Fortran nests accidentally.** If `xvsouris_` somehow calls back into `xvsouris_` (via a Qt signal into Fortran?), the nested `QEventLoop`s stack forever. Phase 5 should add an assertion: `Q_ASSERT(XvueApp::blockingDepth_ < kMaxNestedDepth)` with `kMaxNestedDepth = 4` as a safety net. Phase 6 modal dialogs will legitimately nest once — so depth ≥ 2 is fine, ≥ 5 is a bug.

3. **Filter installed on the wrong widget.** `installEventFilter` must be on `XvueCanvas`, not on `XvueWindow` or `qApp`. Mouse events posted to the canvas will not propagate up to the window unless no handler eats them; an event filter on the wrong level sees different events. **Verify in a unit test** by `QTest::mouseMove(canvas, QPoint(10,10))` → `waitForEvent()` should return with `(x=10, y=10)`.

4. **Keyboard focus trap.** `QTest::keyClick(canvas, ...)` requires the canvas to have **keyboard focus**. `XvueCanvas` today has no `setFocusPolicy()` call (verified — `xvue_qt_canvas.cpp:13-26` just sets `WA_OpaquePaintEvent`). Phase 5 MUST add `setFocusPolicy(Qt::StrongFocus)` in the canvas ctor, or the keyboard branch of the bridge is dead code in a real window. [VERIFIED absent; Phase 5 must add.]

5. **Wayland `QCursor::setPos` no-op.** D-09 documented. The Qt 6 xcb platform plugin maps `QCursor::setPos` to `XWarpPointer`, which works on X11 and XWayland but is rejected by most pure-Wayland compositors. `mail/trfasevo.f:202` is the only caller — if the mesher feels "laggy" on Wayland, document + move on (Phase 5 non-goal).

6. **`blockingDepth` leak on early return / exception.** RAII guard is mandatory — a raw `++/--` will leak if `loop.exec()` throws. Use a helper struct:
   ```cpp
   struct BlockingDepthGuard {
       BlockingDepthGuard() { ++XvueApp::blockingDepth_; }
       ~BlockingDepthGuard() { --XvueApp::blockingDepth_; }
   };
   ```
   Placed as the very first local in `waitForEvent()`.

7. **`AA_CompressHighFrequencyEvents` — set-too-late is NOT silently ignored on X11 (D-05 assumption contradicts Qt docs).** Qt docs list specific attributes that require pre-ctor setting; `AA_CompressHighFrequencyEvents` is **not** on that list (https://doc.qt.io/qt-6/qcoreapplication.html#setAttribute). **AND on X11 it defaults to true already.** D-05's "MUST before QApplication" is defensive folklore; it's safe to do but not load-bearing. **Plan action:** set it in `XvueApp::ensure()` before `QApplication` ctor (matches D-05 literally), document that this is defensive. If the planner finds the before-ctor requirement confirmed in Qt source, upgrade from [ASSUMED] to [VERIFIED]. See Assumption A1.

8. **`pending_x_, pending_y_` stale on first motion after a button press.** If a `ButtonPress` arrives first, then the loop quits; if a `ButtonRelease` arrives with no intervening motion, `pending_x_/pending_y_` might still hold the previous button's position. **Fix:** every code branch in `eventFilter()` that sets `notypeevent` also sets `pending_x_, pending_y_` from the current event's position. Never trust stale stash.

9. **Deferred-quit timer can fire while a new `waitForEvent()` is already running** if the developer forgets to reset `quit_pending_` at the start of each `waitForEvent()` invocation. The RAII guard struct should also reset this flag on construction.

10. **`mempxaccro` pixmap on resize.** The canvas can resize (Phase 2 `resizeEvent` reallocates `backing_`). The `mempxaccro_` sprite is size-invariant (always 13×13) so it's fine, **but** the `accroche_undo_tile_` saved tile might point at a location outside the new backing. **Fix:** on resize, invalidate `accroche_undo_tile_` (set to nullptr); next xvsouris2_ motion allocates a fresh one.

11. **`initaccrochage_` called before `ensure()`.** If Fortran calls `initaccrochage_` before `xvinitgraphique_`, the canvas does not exist yet and `state->backing_` is null. Guard with `if (!win || !win->canvas()) return;` and an `XvueApp::ensure()` first-statement call, per established pattern (`xvue_qt_api.cpp:905`).

---

## 12. Open Questions / Empirical Risks

1. **Motion-event rate under real fast drag — A/B measurement needed.**
   - **What we know:** X11 `XEventsQueued(QueuedAfterFlush)` drops all but the tail of a motion burst; Qt `AA_CompressHighFrequencyEvents` is default-on and should produce a similar rate.
   - **What's unclear:** the actual measured rate difference. Could be 1×, 2×, 5×, or 10×. If > 5× we have a problem.
   - **Recommendation:** in the very first empirical run of the bridge, add a `fprintf(stderr, ...)` counter in `eventFilter()` that logs the motion count per `waitForEvent()` invocation. Run `testa/pan2d`, drag 2 seconds, compare X11 vs Qt counts. Log to `05-VALIDATION.md`. If Qt is > 3× X11, investigate further (may need a manual per-filter drop rule).

2. **Whether `AA_CompressHighFrequencyEvents` compression happens before or after the event filter.** [ASSUMED before — see §4.] Worth a 5-minute confirmation with a printf counter during the same empirical run as (1).

3. **QKeyEvent::text() on AZERTY for `@`.** [ASSUMED that `text()` returns `"@"` on AltGr+0; needs live-keyboard verification by the developer. Mitigated defensively by adding `Qt::Key_At` to the fallback switch per §7.]

4. **`xvsouris2_` accrochage XOR-blit visual parity.** The Phase 4 save/restore approach (§6 Strategy B) is visually *different* from X11's raster-op XOR on an indexed palette. The square highlight will look different — possibly *better* (anti-aliased, fully opaque black on any underlying color). **Flag for developer A/B judgment in VALIDATION.md.** Not a correctness bug, a UX delta.

5. **`xvpause_` AUTOEXIT extension.** Plan needs to decide: reuse `MEFISTO_XVSOURIS_AUTOEXIT` or add a sibling `MEFISTO_XVPAUSE_AUTOEXIT`. Recommendation: reuse (§8).

---

## Runtime State Inventory

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | None — Phase 5 is pure code/state additions; no database or mesh data schema changes | None |
| Live service config | None — no external services | None |
| OS-registered state | None | None |
| Secrets / env vars | `MEFISTO_XVSOURIS_AUTOEXIT`, `MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS`, `MEFISTO_XVFERMER_READY_FILE`, `MEFISTO_XVFERMER_HOLD_MS` — all reused byte-for-byte (D-10). **Phase 5 extends AUTOEXIT semantics to also short-circuit `xvpause_` — same env var, no new knob (§8).** | Extend reader in `xvpause_` (both xvuelc.c and xvue_qt_api.cpp) |
| Build artifacts | `pp/pp*_qt` executables will change on every Phase 5 rebuild (new bridge TU). The `libxvueqt.a` static archive grows by ~1 file. `bin/cbl_tout` + `bin/cbl_tout_qt` share `pp/` — never run in parallel (per user memory feedback). | Normal rebuild after each task |

---

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| Qt6 Core/Gui/Widgets | existing | ✓ | from `qt6-base-dev` | — |
| Qt6 Test (`QTest`) | new unit tests | ✓ (ships with `qt6-base-dev`) | matches Qt6 | — |
| xvfb (`xvfb-run`) | headless test runs | ✓ (used by Phase 3 capture scripts) | — | — |
| X11 display (direct or XWayland) | interactive A/B test on `testa/pan2d`, `testa/torus` | developer-local | — | (deferred to next dev session if headless-only machine) |

No new apt packages required. [VERIFIED via existing CMakeLists finding Qt6 Core/Gui/Widgets/PrintSupport.]

---

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| **A1** | `Qt::AA_CompressHighFrequencyEvents` MUST be set before `QApplication` construction (D-05) | §3, Pitfall 7 | LOW — default-true on X11, so setting is defensive regardless. D-05 is still cheap to follow. |
| **A2** | Qt's `AA_CompressHighFrequencyEvents` compression happens *before* event filters see the event (i.e., filter sees already-compressed motion) | §4 | MEDIUM — if compression happens *after* the filter, our `waitForEvent` still works but motion-count diagnostics need re-interpretation. 5-minute printf test in Task 5 confirms. |
| **A3** | `QKeyEvent::text()` returns `"@"` on a French AZERTY keyboard when user presses AltGr+0 | §7 | LOW — mitigated by adding `Qt::Key_At → 64` to the defensive fallback switch (§7). Even if `text()` is empty, the fallback catches it. |
| **A4** | `QEventLoop::exec()` returns immediately if `quit()` was called before `exec()` | §2 | LOW — standard Qt idiom; widely relied on in modal-dialog code. Phase 5 doesn't depend on this pathological order anyway (loop is constructed fresh each call). |
| **A5** | `xvsouris2_` Phase 4 save-restore Strategy B produces visually-acceptable accrochage highlights vs X11 XOR | §6, §12 Q4 | MEDIUM — UX delta, not correctness. Developer A/B judgment decides. |

**If the planner wants to eliminate A2:** add a 5-minute diagnostic printf task. **If the planner wants to eliminate A1:** read `qtbase/src/corelib/kernel/qcoreapplication.cpp::setAttribute` source directly.

---

## Sources

### Primary (HIGH confidence)

- `xvue/xvuelc.c:2183-2315` — `xvsouris_` canonical X11 body (VERIFIED via Read)
- `xvue/xvuelc.c:2317-2480` — `xvsouris2_` canonical body with `mempxaccro` (VERIFIED)
- `xvue/xvuelc.c:561-609` — `initaccrochage_` canonical body (VERIFIED)
- `xvue/xvuelc.c:2482-2490` — `deplsouris_` canonical body (VERIFIED)
- `xvue/xvuelc.c:2516-2531` — `xvpause_` canonical body (VERIFIED)
- `xvue/xvuelc.c:141-143` — `mempxaccro` / `lmempxaccro` / `hmempxaccro` declarations (VERIFIED)
- `xvue/qt/src/xvue_qt_api.cpp:893-993` — current Phase 4 stubs for the four targets (VERIFIED)
- `xvue/qt/src/xvue_qt_app.{h,cpp}` — XvueApp ensure + atexit pattern (VERIFIED; D-03 target)
- `xvue/qt/src/xvue_qt_canvas.{h,cpp}` — canvas construction, missing `setFocusPolicy` (VERIFIED — Pitfall 4)
- `xvue/qt/src/xvue_qt_state.h:31-38` — backing_ + saved_canvas_ ownership pattern (VERIFIED; Strategy B template)
- `xvue/qt/src/xvue_qt_api.cpp:467-497` — Phase 4 sauvefenetre_/restaurefenetre_ as Strategy B template (VERIFIED)
- `xvue/qt/CMakeLists.txt` — build wiring, verify_no_exec target (VERIFIED)
- `xvue/qt/cmake/verify_no_exec.sh:23` — grep regex confirming QEventLoop::exec is not matched (VERIFIED)
- `bin/xvtest-capture.sh`, `bin/qt-capture.sh`, `bin/testa-capture.sh` — AUTOEXIT env var plumbing (VERIFIED)
- Fortran call sites: `xvue/zoom2d3.f:28,65`, `xvue/saclav.f:61`, `util/searclic.f:68`, `prpr/xvtest1.f:156+`, `mail/trfasevo.f:202`, 30+ `CALL XVPAUSE` sites (all VERIFIED via grep)
- `.planning/REQUIREMENTS.md` EVENT-01..EVENT-08 (VERIFIED)
- `.planning/phases/05-event-bridge-blocking-reads/05-CONTEXT.md` D-01..D-10 (VERIFIED)
- `CLAUDE.md` — project constraints (VERIFIED)

### Secondary (MEDIUM confidence — official Qt docs)

- https://doc.qt.io/qt-6/qeventloop.html — `exec()` / `quit()` / `ProcessEventsFlags` (CITED)
- https://doc.qt.io/qt-6/qt.html — `Qt::ApplicationAttribute` / `AA_CompressHighFrequencyEvents` description (CITED — default-true on X11, compression at dispatch time)
- https://doc.qt.io/qt-6/qcoreapplication.html#setAttribute — general rule about pre-ctor attributes (CITED — `AA_CompressHighFrequencyEvents` not specifically listed)
- https://doc.qt.io/qt-6/eventsandfilters.html — `installEventFilter` ordering + return-value semantics (CITED)
- https://doc.qt.io/qt-6/qkeyevent.html#text — `text()` Unicode semantics (CITED)
- https://doc.qt.io/qt-6/qtest.html#keyClick — QTest headless input (CITED)

### Tertiary / [ASSUMED] — flagged in Assumptions Log

- Qt source internals for `QEventLoop::exec()` quit-before-exec handling (A4)
- Event-filter-vs-compression ordering (A2)
- `QKeyEvent::text()` on French AZERTY + AltGr composition (A3)

---

## Metadata

**Confidence breakdown:**
- Fortran ABI contract (§5): **HIGH** — verified against source + grep on every call site.
- mempxaccro semantics (§6): **HIGH** — verified against xvuelc.c lines 561-609 + 2415-2450.
- Qt event loop mechanics (§2, §3): **MEDIUM** — Qt 6 official docs confirm the high-level API; specific micro-ordering (compression vs filter, deferred-timer vs queued events) is [ASSUMED] pending empirical confirmation.
- File layout (§9): **HIGH** — mechanical reading of existing Phase 1–4 conventions.
- Validation architecture (§10): **HIGH** — QTest is the standard Qt unit-test framework; the A/B drag test is the documented success criterion from ROADMAP Phase 5 SC#5.
- Pitfalls (§11): **HIGH** for items verified against source; MEDIUM for items derived from Qt docs.

**Research date:** 2026-04-14
**Valid until:** Phase 5 wave 0 kickoff — empirical numbers in §3 and §12 should be taken at the start of Task 5 and written back into VALIDATION.md.
