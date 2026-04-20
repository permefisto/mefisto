# xvue-qt — Qt 6 reimplementation of the MEFISTO X11 graphics layer

`xvue/qt/` is the Qt 6 reimplementation of the Fortran-facing graphics API
historically provided by `xvue/xvuelc.c` (Xlib + Motif). The Qt backend is
an in-progress replacement. The 57 `extern "C"` entry points listed in
`xvue/qt/include/xvue_qt_api.h` are byte-compatible with the X11 backend
so Fortran callers (`mail/`, `elas/`, `flui/`, `ther/`, `nlse/`, `prpr/`)
link identically against either backend — the choice of backend is made
at build time by `bin/cbl_tout` (X11) vs `bin/cbl_tout_qt` (Qt).

Per `CLAUDE.md`, the Qt port is an **incremental** migration goal: new
capabilities stay isolated inside `xvue/qt/`, no Fortran call site ever
has to change, and both backends must keep compiling until the X11
backend is retired.

## Build

```sh
bin/cbl_tout_qt         # full Qt build; outputs pp/pp*_qt binaries
```

Runtime requires `qt6-base-dev` (Core, Gui, Widgets, Test, PrintSupport).
No new apt packages beyond the X11 baseline — tests link against Qt6::Test
which ships with `qt6-base-dev`.

Never run `bin/cbl_tout` and `bin/cbl_tout_qt` in parallel; they share
`pp/` and will clobber each other's artifacts mid-link
(see `MEMORY.md :: feedback_parallel_builds_share_pp`).

## Phase 5 — Event bridge

The Qt backend uses a **nested `QEventLoop`** inside
`XvueEventBridge::waitForEvent` to implement the blocking Fortran
primitives (`xvsouris_`, `xvsouris2_`, `xvpause_`). Each blocking call
constructs its own stack-local `QEventLoop`, installs the bridge as an
event filter on `XvueCanvas`, and returns as soon as a matching event is
captured. `deplsouris_` is non-blocking and bypasses the bridge
(`QCursor::setPos(canvas->mapToGlobal(...))`).

### Architecture

```
Fortran       xvue_qt_api.cpp     xvue_qt_event.cpp     Qt 6
 CALL XVSOURIS -> xvsouris_ --+-> XvueEventBridge::waitForEvent
                              |       |
                              |       +-- BlockingDepthGuard (RAII on XvueApp::blockingDepth_)
                              |       +-- save-restore 7 filter members (nested-call safe)
                              |       +-- QEventLoop loop; loop.exec()
                              |                   |
                              |                   V
                              |            eventFilter on XvueCanvas
                              |              - MouseMove:    arm QTimer::singleShot(0, &loop, quit);  return true
                              |              - MouseButton:  loop.quit();                             return true
                              |              - KeyPress:     loop.quit();                             return true
                              |              - (anything else when loop_ null): pass-through          return false
                              |
                              +-- AUTOEXIT short-circuit (MEFISTO_XVSOURIS_AUTOEXIT=1 -> no event loop)
```

- `XvueEventBridge` is a `QObject`, parented to `XvueWindow`, installed
  as an event filter on `XvueCanvas` at window construction. Qt
  parent-child destruction cleans it up automatically when the window
  is destroyed. A fresh bridge is built on every
  `xvfermer_` → `xvinitgraphique_` reopen cycle.
- `BlockingDepthGuard` (RAII friend of `XvueApp`) increments / decrements
  `XvueApp::blockingDepth_` around every `waitForEvent` call.
  `XvueApp::blockingDepth()` is the Phase 6 gate for refusing modal
  dialogs while a Fortran blocking read is active (`blockingDepth() > 0`
  ⇒ refuse `QFileDialog::exec()` and friends).
- The 7 filter members (`loop_`, `mode_`, `pending_`, `quit_pending_`,
  `motion_count_`, `items_`, `pmin0_`) are saved on entry and restored on
  exit, so a single bridge can dispatch arbitrarily nested `waitForEvent`
  calls without corrupting outer-call state.

### Motion coalescing (X11 XEventsQueued parity)

Fast mouse drags could deliver dozens of `QEvent::MouseMove` per tick.
Returning each one to the Fortran caller would make `testa/pan2d` and
`testa/torus` rubber-band picks stutter. The bridge therefore **coalesces
a motion burst into a single `waitForEvent` return with the tail
position**, matching X11's `XEventsQueued(QueuedAfterFlush)` semantics
from `xvue/xvuelc.c:2248-2263`.

The mechanism is a zero-delay `QTimer::singleShot(0, loop_, &QEventLoop::quit)`
deferred quit, armed on the first `MouseMove` of a burst. The timer
enqueues `loop.quit()` at the tail of the event queue, so any motion
events already queued ahead of it are dispatched first (each overwriting
`pending_.x/y`) before the timer fires. Result: `loop.exec()` returns
with the last `(x, y)` in the burst, zero added latency, no hand-rolled
de-duplication. `Qt::AA_CompressHighFrequencyEvents` is also set in
`XvueApp::ensure()` before `QApplication` construction; it defaults to
true on X11 and provides a second layer of compression at the
`QWindowSystemInterface` level (before the filter sees the event).

Plan 06 Phase A empirically validated the composite: on paced input
(~3 ms / waypoint) the bridge returns once per motion with
`motion_count=1`; on a 100-move burst with 0 ms gap (Plan 03 unit test
`testMotionCoalescingBurst`) it returns once with `motion_count=100`.

### Accrochage (snap-highlight) — Strategy B

`xvsouris2_` drives the interactive mesh-picking loop in
`xvue/saclav.f`. On each mouse motion or press while inside a
`WaitMode::Souris2` bridge call, the filter runs a nearest-item search
over the `items[]` array, redraws a 13×13 highlight sprite
(`XvueState::mempxaccro_`) centered on the nearest item, and preserves
the canvas pixels under the sprite via the Phase 4 save-restore pattern
(`XvueState::accroche_undo_tile_`). The previous sprite's undo tile is
blitted back on motion-to-new-item and on release/abort, so the canvas
stays clean for the next call.

This is **Strategy B** in `05-RESEARCH.md §6` — a save-restore-blit
pattern that works on any underlying canvas color, unlike the X11
`GXand`/`GXorInverted` XOR raster-op trick that depended on an indexed
palette. Phase 4's `saved_canvas_` mechanism is reused byte-for-byte;
no new moving parts.

## Developer diagnostics — `MEFISTO_XVSOURIS_DEBUG`

Set this environment variable to see per-`waitForEvent` event counts on
stderr:

```sh
MEFISTO_XVSOURIS_DEBUG=1 pp/ppmail_qt testa/pan2d 2> /tmp/qt_motion.log
```

Each return emits one line:

```
[xvsouris] mode=<0..2> notypeevent=<code> nbc=<code> x=<px> y=<px> motion_count=<N> depth=<N>
```

The env var is cached on first access (C++17 static-local thread-safe
init), so it cannot be flipped mid-process. Off by default — production
stderr stays clean. The counter is retained behind the debug flag for
Phase 6 use: modal-dialog-induced motion behavior will benefit from the
same diagnostic during that phase's validation.

## Known caveats

### Wayland `QCursor::setPos` no-op (D-09)

`deplsouris_` (called from `mail/trfasevo.f:202` for programmatic
cursor repositioning during vertex editing) maps to
`QCursor::setPos(canvas->mapToGlobal(QPoint(*x, *y)))`. The Qt 6 xcb
platform plugin forwards this to `XWarpPointer`, which **works on X11
and XWayland** but is **rejected silently by most pure-Wayland
compositors**.

Phase 5 explicitly supports X11 / XWayland sessions only. On a
pure-Wayland compositor:

- Interactive vertex-repositioning workflows that depend on
  programmatic cursor warp will feel degraded (the mouse no longer
  snaps to the computed coordinate). Functional correctness is
  preserved — every event still flows through the bridge — but the
  visual cue is lost.
- Headless capture scripts (`MEFISTO_XVSOURIS_AUTOEXIT=1`) and
  non-warp-dependent interactive sessions are unaffected.

No workaround is attempted. The caveat is documented here rather than
in user-facing docs because MEFISTO's target session for this release
cycle is X11/XWayland. See `.planning/phases/05-event-bridge-blocking-reads/05-CONTEXT.md :: D-09`.

### `@` abort on French AZERTY (defensive mitigation)

Esc (27) and `@` (64) both map to the `notypeevent=0` abort path for
parity with `xvue/xvuelc.c:2308-2309` and French tutorial muscle
memory (see `05-CONTEXT.md :: D-06`). `translateKey()` tries
`QKeyEvent::text().toLatin1()` first and falls back to a
`Qt::Key_At → 64` switch case, so both paths cover the AZERTY-with-AltGr
layout. Live AZERTY verification was deferred from Plan 06; the
defensive fallback means either path alone is sufficient.

## Headless short-circuit — `MEFISTO_XVSOURIS_AUTOEXIT`

The bridge honors `MEFISTO_XVSOURIS_AUTOEXIT=1` across **all four**
blocking entry points (`xvsouris_`, `xvsouris2_`, `xvpause_`, plus
`deplsouris_` which is non-blocking and unaffected). Phase 5 Plan 04
extended this env var from Phase 3's `xvsouris_`-only coverage to
`xvpause_` in both backends (`xvue/xvuelc.c::xvpause_` and
`xvue/qt/src/xvue_qt_api.cpp::xvpause_`) so that headless capture
scripts (`bin/xvtest-capture.sh`, `bin/qt-capture.sh`,
`bin/testa-capture.sh`) never hang on `CALL XVPAUSE`.

Supported auxiliary env vars:

| Env var                          | Effect                                                                      |
|----------------------------------|-----------------------------------------------------------------------------|
| `MEFISTO_XVSOURIS_AUTOEXIT`      | `=1` short-circuits all blocking reads (synthesize a key event and return). |
| `MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS` | Delay in ms before the synthesized return (clamped 0..60000). Default 100. |
| `MEFISTO_XVFERMER_READY_FILE`    | File to `touch` when the window is fully up (used by capture scripts).      |
| `MEFISTO_XVFERMER_HOLD_MS`       | How long to hold the window open before `xvfermer_` returns.                |
| `MEFISTO_XVSOURIS_DEBUG`         | `=1` enables the `[xvsouris] motion_count=...` stderr diagnostic.           |

## Phase 6 handoff

- Menu / dialog code will gate modal spawn on `XvueApp::blockingDepth() > 0`
  (refuse `QFileDialog::exec()` and friends while a Fortran blocking
  read is active).
- The bridge's `XvueEventBridge::cleanupAccrochage()` is already invoked
  on any terminating event (release / Esc / `@` / other key) so canvas
  state is guaranteed clean across Phase 5 → Phase 6 boundaries.
- `MEFISTO_XVSOURIS_DEBUG` will be reused for Phase 6 modal-dialog
  motion diagnostics — no new env var needed.

## Phase 6.0 — Shared shell, menu bridge, dialogs, persistence

Phase 6.0 adds the full Qt `QMainWindow` chrome that every `pp*_qt`
executable shares: menu bar, toolbar, status bar, console dock, About
dialog, Preferences dialog, Recent Projects submenu, QSettings-backed
geometry/state persistence, and the `XvueMenuBridge` queue that routes
`QAction::triggered` events into the Fortran text-lexicon via
`XvueEventBridge::waitForEvent`'s pre-drain step.

The phase deliberately stops short of wiring any module-specific menus
(Mesh, Solve, etc.) — those land in 6.1..6.5, one module at a time.

### Components added by 6.0

| Component | File(s) | Plan |
|-----------|---------|------|
| Bilingual string table (UX-08) | `xvue_qt_i18n.{h,cpp}` | 02 |
| Menu → Fortran char queue (UX-02, UX-03) | `xvue_qt_menu_bridge.{h,cpp}` | 01 (stub) → 02 (body) |
| QSettings persistence (UX-06, UX-07) | `xvue_qt_prefs.{h,cpp}` | 01 (stub) → 02 (body) |
| Stdout-capture console dock (UX-10, UX-11) | `xvue_qt_console_dock.{h,cpp}`, `xvue_qt_error_batcher.{h,cpp}` | 01 (stub) → 04 (body) |
| About / Preferences / Recent Projects (UX-04, UX-06, UX-09) | `xvue_qt_about_dialog.{h,cpp}`, `xvue_qt_preferences.{h,cpp}`, `xvue_qt_recent_projects.{h,cpp}` | 01 (stub) → 04 (body) |
| Canvas gestures + empty-state (UX-12, UI-SPEC Flag #3) | `xvue_qt_canvas.{h,cpp}`, `xvue_qt_state.{h,cpp}` additions | 05 |
| Event-bridge pre-drain (UX-02, UX-03) | `xvue_qt_event.cpp` extension | 03 |
| Window/App shell + module init dispatch (UX-01, UX-04, UX-13, UX-06, UX-07) | `xvue_qt_window.{h,cpp}`, `xvue_qt_app.{h,cpp}`, `xvue_qt_api.cpp` | 06 |

### Menu-bridge drain flow (UX-02 / UX-03)

```
 user clicks a QAction
         │
         ▼
 XvueMenuBridge::queueLexicon(cmd)
         │   (for each QChar: pendingChars_.enqueue(ch.toLatin1()); enqueue(13))
         ▼
 XvueEventBridge::waitForEvent         ← pre-drain runs BEFORE QEventLoop::exec()
         │
         ▼
 xvsouris_ returns one char per call   (Fortran SACLAV accumulates into KLG(LHKLG))
         │
         ▼
 CR flushes the command → Fortran DONNMF parser
```

Typed-lexicon fallback (`99;`, `5;90;1;`, …) is preserved unchanged.
D-04 modifier rule: plain alphanumeric + digits + `;` + Esc + Return
flow to the canvas via `XvueEventBridge`; `Ctrl/Alt/Cmd+X` and F-keys
route to QActions.

### D-08 modal re-entrancy guard (UX-04)

When `XvueApp::blockingDepth() > 0` (inside a nested `xvsouris_`), any
`QDialog::exec` path refuses silently and shows a 3-second status-bar
message:

- EN: `"Finish current operation first (type 99;)"`
- FR: `"Terminez l'opération en cours (tapez 99;)"`

The QAction itself is NOT disabled — the guard lives inside
`XvueWindow::refuseIfBlocking()` which each File/Preferences/About slot
calls. The right-click `contextMenuEvent` on `XvueCanvas` also
short-circuits on `blockingDepth() > 0`.

### ABI contract (`xvue_module_init_`)

Phase 6.0 adds one new Fortran-ABI entry: `xvue_module_init_(char* name,
int* name_len)` — entry 58 (Phase 5 ended at 57). Its body:

1. Initialises `XvuePrefs` with a per-module `QSettings` group.
2. Applies the persisted color-scheme preference (UX-13).
3. Wires the console-dock stdout redirect.
4. Dispatches to a module-specific `registerXxxActions` stub (6.1..6.5
   each ship a stronger symbol that replaces the corresponding stub).
5. Marks the menu bridge with `markModuleRegistered()`.
6. Flips `state->has_user_content_` true so the empty-state hint stops
   rendering.

### 6.1..6.5 integration contract

Each per-module phase adds:

1. A single `CALL XVUE_MODULE_INIT('<mod>')` at the top of the module's
   Fortran main (e.g., `prpr/ppmail.f` for 6.1).
2. A stronger symbol `registerXxxActions` that replaces the warn-once
   stub in `xvue_qt_api.cpp`. Each implementation creates module-specific
   `QAction`s, assigns their lexicon via `menuBridge->addMenuItem(...)`,
   and plugs them into `XvueMenuBridge`.
3. A `LEXICON-AUDIT.md` catalog of every typed lexicon command, with
   the 80/20 subset wired to `QAction`s. Long-tail commands remain
   available via the typed lexicon.
4. A per-module manual A/B sweep against the module's `testa/` fixtures,
   covering UX-05 (the module-specific menus) plus the 8 manual-only
   items from `06.0-VALIDATION.md` that need real solver output (UX-10
   stdout capture end-to-end, UX-11 `*** ERREUR` surfacing, etc.).

### Pure-6.0 behavior and the `xvue_qt_mark_user_content()` hook

Pure 6.0 builds run without any `CALL XVUE_MODULE_INIT` (that line is
added by 6.1..6.5). Without a hook, the Plan 05 empty-state hint would
sit over the live mesher canvas chrome because `has_user_content_` would
stay `false` forever.

`xvue_qt_mark_user_content()` (anonymous-namespace helper in
`xvue_qt_api.cpp`) flips the flag on the first Fortran-side drawing
primitive and is wired into every draw entry: `xvue_qt_draw_rect_common`,
`xvue_qt_draw_text_common`, `xvue_qt_restore_from_slot`, `effacer`,
`effacemempx`, `xvfond`, `xvface`, `xvtrait`, `xvtraits`, `xvfacetraits`,
`xvbordarcellipse`, `xvarcellipse`. Idempotent: subsequent calls are a
single branch + return. Once `xvue_module_init_` lands in 6.1..6.5, this
hook remains harmless (the flag is already set).

## References

- `.planning/phases/05-event-bridge-blocking-reads/05-CONTEXT.md` — decisions D-01..D-10
- `.planning/phases/05-event-bridge-blocking-reads/05-RESEARCH.md` — Qt event-loop mechanics, assumptions log
- `.planning/phases/05-event-bridge-blocking-reads/05-VALIDATION.md` — EVENT-01..08 status, Manual A/B log, Phase A evidence
- `xvue/xvuelc.c:2183-2531` — X11 reference semantics for the four ported entry points
- `xvue/qt/README_RESIZE.md` — sibling doc on canvas resize-preserve convention (DRAW-09)
- `xvue/README_COORDS.md` — sibling doc on the coordinate-system convention
