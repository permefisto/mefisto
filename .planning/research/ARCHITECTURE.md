# Architecture Research — xvue-qt

**Domain:** Qt 6 reimplementation of a Fortran-driven X11/Xlib graphics library
**Researched:** 2026-04-10
**Confidence:** HIGH for Qt core patterns (QApplication / QPainter / QPixmap / QEventLoop — stable since Qt 4). MEDIUM for the exact `QEventLoop` nesting discipline around the blocking `xvsouris` call. LOW for one-shot edge cases flagged inline.

## The Core Problem

MEFISTO's Fortran solvers drive graphics through ~60 `extern "C"` subroutine calls in `xvue/xvuelc.c`. The calling contract is **strictly synchronous and blocking**:

- `xvtrait(x1,y1,x2,y2)` draws a line and returns. Fortran expects the pixel to be on screen (or at least in a backing buffer that the next call can read).
- `xvsouris(&kind, &nb, &x, &y)` blocks the Fortran caller inside a `while (!flag) XNextEvent(...)` loop until the user clicks or types a key.
- `xvpause()` blocks until one keystroke.
- `xvfermer()` tears the window down and returns.

There is **no Fortran event loop**. The Fortran code is the event loop — it calls the graphics library, reads an event, branches on it, calls more drawing primitives, and loops. The X11 implementation works because Xlib has no notion of "owning" the main loop: it only emits events when the client calls `XNextEvent`.

**Qt inverts this.** `QApplication::exec()` owns the main thread forever, and `QWidget::paintEvent()` is the only sanctioned place to create a `QPainter` on a widget. A naive port that calls `QApplication::exec()` from the first `xvinitgraphique()` would never return to the Fortran caller. A naive port that calls `QPainter painter(widget); painter.drawLine(...)` from inside `xvtrait()` would either assert (outside `paintEvent`) or paint onto a surface Qt is about to overwrite on the next repaint.

The architecture below resolves both problems with **one central design choice**: we paint into an **offscreen `QPixmap` backing store** that is owned by the library, and the widget's `paintEvent` does nothing but blit that pixmap to the screen. This matches the pre-existing MEFISTO double-buffering semantics (`fenetremempx` / `mempxfenetre` / `sauvefenetre` / `restaurefenetre`) exactly — those functions already assume a pixmap-to-window model.

The blocking event call (`xvsouris`) is resolved with a **local `QEventLoop` nested inside the main `QApplication::exec()` loop**, which is the Qt-blessed idiom for "block this function until a signal fires" without stopping the GUI. See the Event Loop Strategy section.

## Standard Architecture

### System Overview

```
┌───────────────────────────────────────────────────────────────────────┐
│                     FORTRAN SOLVERS  (unchanged)                       │
│    mail/  elas/  flui/  ther/  nlse/  reso/  util/                     │
│                                                                        │
│    CALL XVTRAIT(x1,y1,x2,y2)   CALL XVSOURIS(...)                      │
│    CALL XVCOULEUR(icolor)      CALL XVPAUSE()                          │
└──────────────┬────────────────────────────────────────────────────────┘
               │  Fortran → C calling convention (name-mangled `proc()`)
               ▼
┌───────────────────────────────────────────────────────────────────────┐
│               extern "C" ABI SHIM  (xvue/xvue_qt_api.cpp)              │
│   byte-identical names/signatures to today's xvuelc.c                  │
│   each function:  (1) ensures XvueApp is alive                         │
│                   (2) translates args                                  │
│                   (3) delegates to a C++ method on XvueState / XvueCanvas │
└──────────────┬────────────────────────────────────────────────────────┘
               │  no Qt types cross this line
               ▼
┌───────────────────────────────────────────────────────────────────────┐
│                      C++ / Qt 6 IMPLEMENTATION                         │
│                                                                        │
│   ┌─────────────┐   ┌─────────────┐   ┌────────────────────────────┐  │
│   │  XvueApp    │──▶│ XvueWindow  │──▶│      XvueCanvas            │  │
│   │ (singleton) │   │ QMainWindow │   │      QWidget                │  │
│   │ owns        │   │ menus,      │   │      holds backing QPixmap  │  │
│   │ QApplication│   │ toolbars    │   │      paintEvent = blit      │  │
│   └──────┬──────┘   └──────┬──────┘   └─────────────┬──────────────┘  │
│          │                  │                        │                  │
│          ▼                  ▼                        ▼                  │
│   ┌─────────────┐   ┌─────────────┐   ┌────────────────────────────┐  │
│   │ XvueState   │   │  XvueMenu   │   │      XvueEventBridge       │  │
│   │ pen, brush, │   │ Bridge      │   │  installs QObject filter,  │  │
│   │ font, color │   │ Fortran     │   │  wakes a nested QEventLoop │  │
│   │ map,        │   │ callback    │   │  on mouse/key events       │  │
│   │ current     │   │ registry    │   │  → fills xvsouris out-args │  │
│   │ QPainter    │   │             │   │                            │  │
│   └─────────────┘   └─────────────┘   └────────────────────────────┘  │
│                                                                        │
│   ┌────────────────────────────┐   ┌────────────────────────────────┐ │
│   │      XvuePixmapStack        │   │        XvueExport              │ │
│   │  sauvefenetre/restaurefenetre│  │  PNG/JPEG via QImageWriter     │ │
│   │  named pixmap slots          │  │  GIF animation (frame accum)   │ │
│   │  + hidden "mempx" scratch    │  │  PostScript via QPrinter /     │ │
│   │                              │  │     legacy hand-rolled writer  │ │
│   └────────────────────────────┘   └────────────────────────────────┘ │
└───────────────────────────────────────────────────────────────────────┘
```

### Component Responsibilities

| Component | Responsibility | One-line implementation |
|-----------|----------------|------------------------|
| **`XvueApp`** | Owns the single `QApplication` instance for the process lifetime; ensures it is created lazily on the first `extern "C"` call and destroyed via `atexit`. | `static std::unique_ptr<QApplication>` guarded by `std::call_once`, with fabricated `argc`/`argv` so Qt is happy. |
| **`XvueWindow`** | `QMainWindow` subclass holding the menu bar, toolbar, status bar, and the central `XvueCanvas`. Receives close events and forwards them through `XvueEventBridge`. | Constructed the first time `xvinitgraphique()` is called. |
| **`XvueCanvas`** | `QWidget` subclass. Its `paintEvent` does nothing but `QPainter(this).drawPixmap(0,0, *backing)`. No drawing logic lives here. Resize events resize the backing pixmap (preserving the old content as a sub-blit). | A thin widget; the real painting happens outside `paintEvent`. |
| **`XvueState`** | Current pen color, line width, dash style, brush, font, palette colormap, background, and the single long-lived `QPainter*` bound to the backing pixmap. Equivalent to Xlib's `GC` (graphics context). | One struct, pointer held by `XvueCanvas`. |
| **`XvuePixmapStack`** | Named off-screen pixmap slots backing `fenetremempx` / `mempxfenetre` / `sauvefenetre` / `restaurefenetre` / `sauvemempx` / `restauremempx` / `effacemempx`. Owns a primary backing pixmap and an auxiliary "mempx" pixmap. | `std::array<QPixmap, N>` + enum. |
| **`XvueEventBridge`** | Converts Qt events (`QMouseEvent`, `QKeyEvent`, `QCloseEvent`) into the integer quadruple `(notypeevent, nbc, x1, y1)` expected by `xvsouris`. Provides a single method `waitForEvent(out…)` that runs a local `QEventLoop` until the next matching event, then returns. | See Event Loop Strategy below. |
| **`XvueMenuBridge`** | Maps `QAction::triggered` signals to the Fortran-visible command strings that today come from the text lexicon. Each `QAction` carries a string payload that is pushed into the lexicon input queue, so the rest of the Fortran code — which reads its commands through `util/` lexicon routines — is untouched. | `QSignalMapper`-style dispatch, or direct lambda captures. |
| **`XvueExport`** | PNG/JPEG single-frame save (`QImageWriter`), animated GIF build-up (accumulate `QImage` frames, write with the GIF encoder), PostScript export. PostScript: keep the legacy hand-rolled writer initially (it is already present in `xvuelc.c`) to preserve byte-for-byte output, then migrate to `QPrinter` after A/B validation. | Pure functions; no widget interaction. |
| **`XvueAbiShim`** (`xvue/xvue_qt_api.cpp`) | ~60 `extern "C"` entry points, each a ≤10-line function that ensures `XvueApp` is up, looks up the `XvueCanvas`/`XvueState`, performs one Qt call, and returns. Nothing more. | Source of truth: `xvue/xvuelc.c` function table. |

### Singleton discipline for `QApplication`

**Problem.** Qt requires exactly one `QApplication` per process. It must be constructed on the main thread *before* any `QWidget`, and destroyed *after* the last `QWidget`. It cannot be re-created. The Fortran main program does not know any of this.

**Approach.**

```cpp
// xvue_qt_app.cpp
class XvueApp {
public:
  static XvueApp& instance();          // lazy, thread-safe via std::call_once
  QApplication& qapp()  { return *qapp_; }
  XvueWindow&   window(){ return *window_; }
private:
  XvueApp();
  static int    fake_argc_;
  static char*  fake_argv_[];
  std::unique_ptr<QApplication> qapp_;
  std::unique_ptr<XvueWindow>   window_;   // created lazily on first xvinitgraphique
};
```

- **Construction** is triggered by the first `extern "C"` call that needs it (typically `xvinitgraphique`). `std::call_once` prevents double-init.
- `argc`/`argv` are synthesized (Qt retains references, so they must have static storage).
- **Destruction** is hung off an `atexit` handler installed the first time `XvueApp::instance()` runs. That handler posts `QApplication::quit()` and joins. This is safer than relying on `~XvueApp()` at static-destruction time, because the order of static destructors relative to `QApplication` is hostile (Qt internals assume the application is still alive).
- `xvfermer()` **does not destroy** `XvueApp` — it only hides/closes the `XvueWindow` (and resets `XvueState`). A subsequent `xvinitgraphique()` from the same process re-creates a new window inside the still-alive `QApplication`. This matches what `ppmail` does: close the mesh display, pop back to the lexicon, then open a new display for the solution.

**Confidence:** HIGH. This is the documented Qt idiom for embedding Qt inside a non-Qt main (Qt docs "Qt in namespace" / "Embedding Qt in another application").

## Recommended Project Structure

```
xvue/
├── CMakeLists.txt                 # new — the only CMake file in the tree
├── xvuelc.c                       # PRESERVED during transition (legacy X11 backend)
├── xvue_qt_api.cpp                # extern "C" shim; implements all ~60 entry points
├── xvue_qt_api.h                  # declarations of the extern "C" ABI (shared by both backends)
├── xvue_qt_app.{h,cpp}            # XvueApp singleton, argc/argv fabrication, atexit handler
├── xvue_qt_window.{h,cpp}         # XvueWindow (QMainWindow), menu/toolbar shell
├── xvue_qt_canvas.{h,cpp}         # XvueCanvas (QWidget), paintEvent, resizeEvent
├── xvue_qt_state.{h,cpp}          # XvueState — pen/brush/font/palette + current QPainter
├── xvue_qt_pixmaps.{h,cpp}        # XvuePixmapStack — mempx/sauvefenetre/restaurefenetre
├── xvue_qt_events.{h,cpp}         # XvueEventBridge — waitForEvent / nested QEventLoop
├── xvue_qt_menu.{h,cpp}           # XvueMenuBridge — QAction → lexicon dispatch
├── xvue_qt_export.{h,cpp}         # XvueExport — PNG/JPEG/GIF/PS output
├── xvue_qt_colors.{h,cpp}         # Palette + xvCouleursImposees + xvStockeRGBtoColormap
├── xvue_qt_fonts.{h,cpp}          # xvchargefonte, xvnbpixeltexte
└── *.f                            # existing Fortran wrappers around the extern "C" API (UNTOUCHED)
```

### Structure Rationale

- **One Qt concern per file.** The Qt build must stay understandable for one developer. Each file owns exactly one of the responsibilities listed above. ~12 source files total is small enough to hold in your head and large enough that no single file is a 3000-line monolith like `xvuelc.c` became.
- **Shared header `xvue_qt_api.h`** carries the `extern "C"` declarations. Both `xvuelc.c` (X11) and `xvue_qt_api.cpp` (Qt) include it; whichever one is linked decides the backend. This makes the build-time switch (environment variable or `bin/cbl_tout_qt`) as simple as "link libxvue_qt.a or link xvuelc.o".
- **CMake confined to `xvue/`.** `xvue/CMakeLists.txt` calls `find_package(Qt6 COMPONENTS Widgets Gui Core REQUIRED)`, turns on `CMAKE_AUTOMOC`, compiles the C++ files into a **static archive** (`libxvue_qt.a`), and stops. The existing Fortran shell linker picks up this archive exactly the way it picks up `xvuelc.o` today. No parent `CMakeLists.txt`, no Fortran build in CMake.
- **Legacy `xvuelc.c` stays in-tree** for the one-release A/B window, as PROJECT.md mandates. The two backends never coexist in one binary; they are two parallel build targets.

## Data Flow: Fortran call → pixel on screen

### Drawing primitive: `CALL XVTRAIT(x1,y1,x2,y2)`

```
 Fortran: CALL XVTRAIT(ix1, iy1, ix2, iy2)
    │
    ▼  F77 → C name mangling (xvtrait_ or xvtrait depending on dialect)
 xvue_qt_api.cpp::xvtrait_(int* x1, int* y1, int* x2, int* y2)
    │
    │  1. XvueApp::instance();                        // no-op after first call
    │  2. auto& st = XvueApp::instance().window().canvas().state();
    │  3. st.painter().drawLine(*x1, *y1, *x2, *y2);  // paints INTO the backing QPixmap
    │  4. canvas().scheduleUpdate();                  // sets a dirty flag, posts QWidget::update()
    │
    ▼  returns immediately to Fortran
 Fortran: (next instruction)

 ...later, when control returns to the Qt event loop:

 XvueCanvas::paintEvent(QPaintEvent*)
    │
    │  QPainter p(this);
    │  p.drawPixmap(0, 0, backingPixmap_);
    │
    ▼
 Pixels on screen.
```

**Key invariants.**

- The backing `QPixmap` is allocated once (at `xvinitgraphique`) and re-allocated only on window resize. All drawing primitives (`xvtrait`, `xvface`, `xvtexte`, `xvrectangle`, `xvarcellipse`, `xvfrectangle`, `xvbordrectangle`, `xvtraits`, `xvfacetraits`, `xvbordarcellipse`) go through one `QPainter` that is `begin()`'d on the backing pixmap at window creation and `end()`'d only at `xvfermer` / window destruction. This is legal: `QPainter` can be kept open on a `QPixmap` across many unrelated function calls — the restriction ("only paint inside `paintEvent`") applies to painting on `QWidget`, not on `QPixmap`. **Confidence: HIGH** (Qt documentation on `QPainter` lifetime and paint devices).
- **`scheduleUpdate()`** is crucial. Without it, the screen does not refresh until the Fortran code reaches a point where Qt gets CPU (i.e. until `xvsouris` / `xvpause` blocks). We need an explicit hint — call `XvueCanvas::update()` (which posts a paint event) and optionally `QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents)` every N drawing calls or every T milliseconds to pump the paint queue without processing mouse/keyboard (we do not want reentrancy).
- `xvvoir()` — the existing "please refresh the window now" primitive — is the perfect explicit pump point. Implement it as `canvas().repaint()` (synchronous) plus a `QApplication::processEvents()` with `ExcludeUserInputEvents`. Fortran code that currently calls `xvvoir` after a batch of drawing will automatically get the refresh.

### Colour change: `CALL XVCOULEUR(icolor)`

Pure state update on `XvueState`: `painter_->setPen(state_.palette[icolor])`. No paint event. Returns immediately.

### Double buffer: `CALL FENETREMEMPX` / `CALL MEMPXFENETRE`

`fenetremempx` blits current window → mempx pixmap; `mempxfenetre` blits mempx pixmap → window. In the Qt version, both become `QPainter::drawPixmap` calls between two `XvuePixmapStack` slots, and `mempxfenetre` also triggers a canvas update. `sauvefenetre`/`restaurefenetre` work the same way on a second slot. `effacemempx` is `mempx.fill(backgroundColor)`.

### Blocking mouse read: `CALL XVSOURIS(type, nb, x1, y1)` — the hard one

See Event Loop Strategy.

## Event Loop Strategy

This is the single most important design decision of the migration. Four options exist; I recommend **Option B: nested QEventLoop**.

### Option A — Run `QApplication::exec()` once, drive everything from Qt signals

**What:** Convert the Fortran interactive loops into Qt event handlers. Fortran becomes callback-driven.

**Verdict: REJECTED.** This would require rewriting every interactive Fortran driver in `mail/`, `elas/`, `flui/`, `ther/`, `nlse/` to stop blocking on `XVSOURIS`. PROJECT.md is explicit: Fortran code must not notice the migration. This option violates the core invariant.

### Option B — Nested `QEventLoop` per blocking call  ← RECOMMENDED

**What:** `XvueApp` never runs `QApplication::exec()` at top level. Instead, each blocking `extern "C"` entry point (`xvsouris`, `xvsouris2`, `xvpause`, and implicitly `xvfermer`) runs a **local `QEventLoop`** that exits as soon as the event it is waiting for arrives.

```cpp
// xvue_qt_events.cpp
int XvueEventBridge::waitForEvent(int& kind, int& nb, int& x, int& y) {
    QEventLoop loop;
    pendingResult_ = {};                         // reset accumulator
    targetLoop_   = &loop;                       // event filter will quit this loop
    int rc = loop.exec(QEventLoop::AllEvents);   // blocks here; Qt still processes
                                                 // paint events, timers, other widgets
    targetLoop_ = nullptr;
    kind = pendingResult_.kind;
    nb   = pendingResult_.nb;
    x    = pendingResult_.x;
    y    = pendingResult_.y;
    return rc;
}

// Installed as a QObject::eventFilter on XvueCanvas:
bool XvueEventBridge::eventFilter(QObject*, QEvent* ev) {
    if (!targetLoop_) return false;              // not currently blocking
    switch (ev->type()) {
        case QEvent::MouseButtonRelease: /* fill pendingResult_, loop->quit(); */ return true;
        case QEvent::MouseButtonPress:   /* ditto for notypeevent=-1 */           return true;
        case QEvent::MouseMove:          /* ditto for notypeevent=-2 (with XEventsQueued-style
                                            coalescing: only quit if no more mouse events
                                            are pending — use QCoreApplication::hasPendingEvents
                                            or a single-shot QTimer(0) to defer) */
        case QEvent::KeyPress:           /* ditto */                              return true;
        default: return false;
    }
}
```

**Why this works:**
- `QEventLoop::exec()` runs the Qt event pump on the current thread. **Paint events still fire**, so the window stays responsive. **Timers still fire.** Other widgets (menu bar, toolbars) still respond — meaning the user can click a menu item *instead of* the canvas, and that click is routed through `XvueMenuBridge` which pushes a lexicon command.
- The outer Fortran call stack is undisturbed: `xvsouris` blocks, then returns, exactly like the X11 version.
- **Reentrancy is bounded.** `waitForEvent` sets `targetLoop_` before `exec()` and clears it after, so only the innermost call consumes the event. A menu `QAction` that happens to trigger another Fortran path while inside `xvsouris` is an edge case worth thinking about (see Anti-Patterns).
- **Mouse-motion coalescing** — the X11 version uses `XEventsQueued(..., QueuedAfterFlush) <= 0` to drop intermediate `MotionNotify` events and only return the last one. Qt already coalesces mouse moves at the platform layer, but to preserve exact semantics, defer the `loop.quit()` call one tick via `QTimer::singleShot(0, &loop, &QEventLoop::quit)` — this lets any further queued mouse move events update `pendingResult_` before quitting. **Confidence: MEDIUM** — needs validation against the existing behaviour during A/B testing.

**Trade-offs:**
- Nested event loops can reenter in subtle ways. Mitigate with `targetLoop_` guarding and by making `XvueMenuBridge` push commands into a queue instead of calling Fortran back directly.
- Qt warns against nested loops in some contexts (e.g. `QDialog::exec()` is itself a nested loop), but they are an officially-supported pattern. **Confidence: HIGH** (Qt docs on `QEventLoop` explicitly describe this use case: "The event loop can be started again by calling exec() from within a slot.").

### Option C — Separate GUI thread, Fortran thread blocks on a condition variable

**What:** Spawn a thread that runs `QApplication::exec()`. Fortran runs on the main thread and communicates via `QMetaObject::invokeMethod` + `std::condition_variable`.

**Verdict: REJECTED** for v1.
- Qt requires that `QApplication` and all `QWidget`s live on the thread that called `QApplication` constructor. Whether that is "main thread" or "some thread" is flexible, but on many platforms (Linux included) spawning QApplication off the main thread is hostile — Wayland, GLX, and some input plugins assume the main thread.
- Cross-thread marshalling of ~60 entry points (all of them) is an extra design burden with no user-visible benefit for a single developer.
- Race conditions between drawing and paint events become real.
- Only worth considering if Option B proves insufficient during validation. Reserve as a fallback.

### Option D — Poll via `QCoreApplication::processEvents()` in a busy loop

**What:** `xvsouris` loops `while (!event_received) { QApplication::processEvents(); usleep(1000); }`.

**Verdict: REJECTED.** Burns CPU, hurts responsiveness to low-priority events, and the timing is fragile. `QEventLoop::exec()` (Option B) does exactly this correctly without the busy-wait.

### Recommendation

**Option B, with these rules:**
1. Never call `QApplication::exec()` at top level.
2. Every blocking `extern "C"` entry point creates its own local `QEventLoop`.
3. Inside the loop, `XvueEventBridge` is the sole consumer of user events via an event filter on `XvueCanvas` (and on `XvueWindow` for close events).
4. Menu/toolbar actions are queued, not executed reentrantly, while a `waitForEvent` is in flight.
5. Mouse-move coalescing uses a deferred `singleShot(0, quit)` to mirror X11 `XEventsQueued` semantics.

## Thread Model

**Single-threaded. Main thread only.**

- `QApplication`, `XvueWindow`, `XvueCanvas`, all `QPixmap`s, all `QPainter`s, `XvueEventBridge`, everything — main thread.
- Fortran code runs on the main thread and calls into C++ on the main thread.
- `_OMP` executables already use OpenMP inside the solver layers. OpenMP parallel regions **must not** call graphics primitives. This is the current rule already (the X11 code is not thread-safe either); preserve it.
- `XvueExport` GIF encoding *could* be moved to a worker thread in a later phase, but v1 keeps it synchronous on the main thread for simplicity.

**Rationale.** Qt's thread affinity rules (`QObject::moveToThread`, etc.) are a well-known source of bugs. For a single-developer migration under time pressure, one thread is vastly simpler and matches the existing Xlib-era assumption that graphics is single-threaded.

## Phase Order Implications

The component graph above implies a dependency order that the roadmap should follow. Items earlier in the list are prerequisites for items later.

| # | Phase | What lands | Why this order |
|---|-------|-----------|----------------|
| **0** | **CMake + empty backend** | `xvue/CMakeLists.txt`, `xvue_qt_api.{h,cpp}` skeleton with all ~60 symbols present as no-op stubs, built into `libxvue_qt.a`, linked into one test executable (e.g. `pp/ppmail`) via a `bin/cbl_tout_qt` variant. Existing X11 build still works. | Proves the build plumbing before any graphics logic. Unblocks every subsequent phase. |
| **1** | **XvueApp + XvueWindow + XvueCanvas skeleton** | `xvinitgraphique` opens an empty `QMainWindow` with a blank central widget. `xvfermer` closes it. `xvpxecran`/`xvmmecran` return screen metrics via `QScreen`. `xvfond` sets the background. No drawing yet. | You need a window before you can draw into it. Exercises singleton/`atexit`/arg fabrication before the painting code piles on. |
| **2** | **Static drawing primitives + state + backing pixmap** | `XvueState`, `XvuePixmapStack`, and all the pure-drawing entry points: `xvtrait`, `xvftrait`, `xvtraits`, `xvface`, `xvfacetraits`, `xvrectangle`, `xvbordrectangle`, `xvfrectangle`, `xvbordarcellipse`, `xvarcellipse`, `xvcouleur`, `xvtypetrait`, `xvepaisseur`, `xvfond`, `effacer`, `xvvoir`, `xvpxfenetre`. | Half of the ABI. Can be validated by a standalone test (`prpr/xvtest1.f`..`xvtest4.f` already exist) without any event handling. |
| **3** | **Text, fonts, colormap** | `xvchargefonte`, `xvnbpixeltexte`, `xvtexte`, `xvftexte`, `xvCouleursImposees`, `xvStockeRGBtoColormap`, `xvColormapToRGB`, `xvrecuprgbdec`, `xvactivervb`. | Independent of event handling; can be built in parallel with Phase 2 if desired, but the colormap in particular is used by drawing primitives, so sequence it after or alongside. |
| **4** | **Pixmap save/restore** | `fenetremempx`, `mempxfenetre`, `sauvefenetre`, `restaurefenetre`, `sauvemempx`, `restauremempx`, `effacemempx`. | Depends on Phase 2 (backing pixmap exists). This is the double-buffering MEFISTO already uses heavily in the mesher. Needed before you can validate mesher interactions visually. |
| **5** | **Event bridge + blocking reads** | `XvueEventBridge`, nested `QEventLoop`, `xvsouris`, `xvsouris2`, `xvpause`, `deplsouris`, `initaccrochage`. | The architectural pivot. Once this lands, the mesher becomes interactive through the Qt backend. Depends on Phases 1-4 because the loop needs a real canvas with real drawing to be testable. |
| **6** | **Menu/toolbar surface (Level 3 UX)** | `XvueMenuBridge`, `QMenuBar`/`QToolBar`/`QAction` definitions per module (mail/elas/flui/ther/nlse), lexicon-command dispatch. `QDialog`-based About box. Keyboard shortcuts. | Depends on Phase 5 because menu clicks go through the same event bridge plumbing. Can be done per-module, module-at-a-time, using the A/B test cases. |
| **7** | **Image + GIF + PostScript export** | `XvueExport`, `xvinitierps`, `xvimprimerps`, `xvsauverps`, `xvpostscript`. Replace `bin/convertepsgif` (ImageMagick shell-out). | Independent of event handling. Could in principle happen in parallel with Phase 5, but defer until last to keep early phases focused on interactive correctness. |
| **8** | **A/B validation across `testa/` subset** | Run canonical cases per module through both backends, visual comparison, validation log. | Gate for declaring the migration done. |
| **9** | **Retire X11 backend** | After one release cycle: delete `xvuelc.c`, `bin/ccxvue`, `libX11` linker lines, `bin/convertepsgif`, update `README`/`LISEZMOI`/install scripts. | Post-milestone cleanup. |

**Strictly-required ordering:** 0 → 1 → (2 ∥ 3) → 4 → 5 → 6 → 8. Phase 7 is independent and slots in wherever convenient. Phase 9 is trailing cleanup.

**Phases most likely to need deeper research during roadmap execution:**
- **Phase 5** (event bridge) — the mouse-motion coalescing semantics and nested-loop reentrancy under `QAction` triggers are the two places where the design most needs empirical validation. Expect one round of "this doesn't behave exactly like X11" iteration.
- **Phase 7** (PostScript export) — keeping the existing hand-rolled PS output is safest for bit-for-bit parity; migrating to `QPrinter` is the right long-term move but is a separate research question (QPrinter PostScript support was reduced in Qt 6 vs Qt 5). **Confidence: LOW** — verify QPrinter EPS support in Qt 6.10 before committing.
- **Phase 6** (menu surface) — the text lexicon is scattered across modules; the per-module audit itself is the research.

## Architectural Patterns

### Pattern 1: Backing-pixmap retained-mode rendering

**What:** All drawing happens on a persistent offscreen `QPixmap`. The `QWidget::paintEvent` does not draw anything itself — it blits the pixmap.

**When to use:** When you are porting a synchronous/imperative graphics API (Xlib, GDI, OpenGL immediate-mode) to Qt without rewriting callers to become event-driven.

**Trade-offs:**
- **Pro:** Fortran callers see exactly the semantics they had — "draw a line" means "a line is now in the buffer".
- **Pro:** Free double-buffering.
- **Pro:** Matches the existing `fenetremempx` model.
- **Con:** Does not benefit from GPU-accelerated scene graphs (`QGraphicsScene`, `QQuickView`). For MEFISTO's 2D finite-element visualization, this does not matter — the X11 code was already CPU-drawn.
- **Con:** Window resize requires reallocating the backing pixmap. Decide: preserve old content (blit sub-rect) or redraw from scratch (force the Fortran caller to redraw). MEFISTO already has an explicit "redraw" convention in its mesher, so either choice is acceptable; default to preserving content.

### Pattern 2: Nested QEventLoop for synchronous-looking blocking calls

**What:** Inside a blocking C function, create a `QEventLoop loop;` on the stack and call `loop.exec()`. Install an event filter that calls `loop.quit()` on the awaited event.

**When to use:** Exactly this case — when a non-Qt main drives Qt through synchronous C entry points and needs to block until a user event.

**Trade-offs:**
- **Pro:** Keeps the Qt event loop pumping (paints, timers, menu updates).
- **Pro:** Officially supported Qt idiom (it is how `QDialog::exec` works internally).
- **Con:** Deep nesting can cause surprising reentrancy. Guard with a flag so only the innermost `waitForEvent` consumes the event.
- **Con:** Slightly harder to debug than a flat main loop.

### Pattern 3: Frozen ABI shim over a plastic C++ core

**What:** Freeze the `extern "C"` function names and signatures. Let the C++ implementation beneath them be refactored freely.

**When to use:** Gradual migrations where the caller base (here: Fortran) cannot be touched.

**Trade-offs:**
- **Pro:** Callers are decoupled from migration progress.
- **Pro:** Enables side-by-side A/B backends (X11 vs Qt) linked from the same ABI header.
- **Con:** The shim layer can become a dumping ground for awkward translations; discipline required to keep each shim function tiny (≤10 lines) and delegate immediately.

## State Management

```
XvueApp  (singleton, process lifetime)
  └─ QApplication
  └─ XvueWindow               (lifetime: one xvinitgraphique → xvfermer cycle)
       └─ QMenuBar, QToolBar, QStatusBar
       └─ XvueMenuBridge
       └─ XvueCanvas          (central widget)
            └─ XvueState      (pen, brush, font, palette, current QPainter*)
            └─ XvuePixmapStack (backing + mempx + named slots)
            └─ XvueEventBridge (installed as event filter)
```

**Lifetime rules.**
- `XvueApp` lives from the first `extern "C"` call until `atexit`.
- `XvueWindow` (and everything under it) is created in `xvinitgraphique` and destroyed in `xvfermer`. Subsequent `xvinitgraphique` calls create a fresh `XvueWindow` — do not try to keep stale state across a close/reopen cycle.
- `XvueState::painter_` is `begin()`'d on the backing pixmap when the canvas is constructed and `end()`'d in the canvas destructor. Never re-created per-call.
- `XvuePixmapStack` slots are allocated lazily on first `sauve…` call and freed with the canvas.

**Global state** is accessed only through `XvueApp::instance().window().canvas().state()` — never via free globals. This makes unit testing possible (construct a test `XvueApp`, run a drawing script, dump the pixmap, compare).

## Anti-Patterns

### Anti-Pattern 1: Creating a `QPainter` per drawing call

**What people do:** `void xvtrait(…) { QPainter p(canvas); p.drawLine(…); }`

**Why it's wrong:**
1. Painting on a `QWidget` outside `paintEvent` is illegal and asserts in debug builds.
2. Even on a `QPixmap`, constructing/destroying a `QPainter` per call is wasteful and discards pen/brush/font state between calls, so `xvcouleur` followed by `xvtrait` would not work.

**Do this instead:** One long-lived `QPainter` bound to the backing `QPixmap`, owned by `XvueState`, reused across all drawing calls. Recreate it only on window resize (when the backing pixmap is reallocated).

### Anti-Pattern 2: Calling `QApplication::exec()` from `xvinitgraphique`

**What people do:** Assume Qt apps always run `exec()` and call it from the first init function.

**Why it's wrong:** It never returns. The Fortran main is gone forever.

**Do this instead:** Never call `exec()` at the top level. Use `QEventLoop::exec()` locally inside each blocking entry point (Option B).

### Anti-Pattern 3: Letting `QAction::triggered` call Fortran directly from inside `xvsouris`

**What people do:** Connect a `QAction` to a slot that immediately calls a Fortran subroutine.

**Why it's wrong:** The slot fires *while* `xvsouris` is still running inside its nested `QEventLoop`. The Fortran call reenters the lexicon dispatcher *before* the current Fortran frame has returned, breaking the save/exit state machine.

**Do this instead:** `XvueMenuBridge` **queues** the lexicon command string into a `std::deque<std::string> pendingCommands_`. The `xvsouris` bridge checks that queue first — if a command is pending, synthesize a `notypeevent=2` (keyboard) return that injects the command string into the lexicon input buffer. Menu clicks thus look like "the user typed the command" to the rest of the Fortran code.

### Anti-Pattern 4: Destroying `XvueApp` in `xvfermer`

**What people do:** Treat `xvfermer` as "shut down graphics, free everything".

**Why it's wrong:** The user may call `xvinitgraphique` again in the same process (e.g. mesher opens display, closes it, solver opens a new one). Once `QApplication` is destroyed, recreating it in the same process is Qt-hostile — many platform plugins assume one-per-process.

**Do this instead:** `xvfermer` closes the window and resets state. `QApplication` is destroyed only at `atexit`.

### Anti-Pattern 5: Calling graphics primitives from OpenMP parallel regions in `_OMP` executables

**What people do:** Assume Qt is thread-safe because "modern framework".

**Why it's wrong:** `QPainter`, `QPixmap`, and `QWidget` are strictly not thread-safe. Calling them off the main thread is undefined behaviour.

**Do this instead:** Preserve the existing Fortran convention (already enforced in the X11 code): no graphics inside `!$OMP PARALLEL` regions. Document this in `xvue_qt_api.h`.

## Integration Points

### External

| Dependency | Integration | Notes |
|---|---|---|
| Qt 6 (Widgets/Gui/Core) | `find_package(Qt6 ... REQUIRED)` in `xvue/CMakeLists.txt`; static link into `libxvue_qt.a` | apt package `qt6-base-dev` 6.10.2 on Debian trixie |
| gfortran Fortran ABI | `extern "C"` with same name-mangling as today's `xvuelc.c` (trailing underscore on Linux) | Do not change the mangling convention — `xvue/*.f` wrappers depend on it |
| ImageMagick (legacy) | **Dropped** after Phase 7 | `XvueExport` replaces it |
| libX11 (legacy) | **Kept** in parallel during transition, **removed** in Phase 9 | Two build targets share one `xvue_qt_api.h` |

### Internal Boundaries

| Boundary | Communication | Notes |
|---|---|---|
| Fortran solver ↔ ABI shim | `extern "C"` subroutine call (by-reference ints) | Frozen forever |
| ABI shim ↔ C++ core | Direct method calls on `XvueApp::instance()...` | Plastic — refactor freely |
| `XvueEventBridge` ↔ Qt events | `QObject::eventFilter` installed on `XvueCanvas` and `XvueWindow` | Only one filter active at a time |
| `XvueMenuBridge` ↔ lexicon | `std::deque<std::string>` of pending commands, drained by `xvsouris` | Ensures no reentrancy into Fortran from slot context |
| Qt C++ core ↔ legacy `xvuelc.c` | None — they never coexist in one binary | Build-time switch selects one |

## Confidence Assessment

| Area | Confidence | Reason |
|---|---|---|
| `QApplication` singleton lifecycle, argc/argv fabrication, `atexit` shutdown | HIGH | Documented Qt pattern, stable since Qt 4 |
| Backing-pixmap retained-mode rendering with long-lived `QPainter` on `QPixmap` | HIGH | Core Qt paint system; `QPainter` on `QPixmap` is explicitly allowed outside `paintEvent` |
| Nested `QEventLoop` for blocking-style calls | HIGH | Officially supported; how `QDialog::exec` works internally |
| Mouse-move coalescing semantics exactly matching X11 `XEventsQueued` | MEDIUM | The `singleShot(0, quit)` deferral approximates the behaviour but needs empirical validation in A/B tests |
| Menu-click reentrancy via command queue | MEDIUM | Design is sound in principle; will surface edge cases when exercised against real mesher sessions |
| `QPrinter` PostScript support in Qt 6.10 for EPS export | LOW | Qt 6 reduced printing support vs Qt 5; recommend keeping legacy hand-rolled PS writer in Phase 7 and deferring `QPrinter` migration |
| Animated GIF via Qt `QImageWriter` | MEDIUM | `QImageWriter` supports GIF if the plugin is present; Debian ships it but verify at Phase 7 start |
| Exact Fortran name mangling match between `xvuelc.c` and new C++ shim | HIGH | Same compiler (gfortran), same convention; the existing `proc()` macro approach can be reused verbatim |

## Sources

- Qt 6 documentation — `QApplication`, `QEventLoop`, `QPainter`, `QPixmap`, `QWidget::paintEvent` (training data, stable across Qt 4/5/6)
- `/home/drico/git/mefisto/xvue/xvuelc.c` — ~60 `extern "C"` entry points, including the `xvsouris` `while (!flag) XNextEvent(...)` blocking pattern at line 2157, `xvpause` at line 2408, and the existing pixmap double-buffer primitives at lines 1307–1395
- `/home/drico/git/mefisto/.planning/PROJECT.md` — scope constraints: ABI frozen, Linux x86_64 only, CMake confined to `xvue/`, parallel X11 build during transition, Qt 6.10 from apt
- `/home/drico/git/mefisto/.planning/codebase/ARCHITECTURE.md` — confirms `xvue/` is the only C source and the sole X11 containment layer
- `/home/drico/git/mefisto/CLAUDE.md` — "keep graphics calls isolated in `xvue/` so the migration is incremental" is already a project rule

---
*Architecture research for: Qt 6 migration of MEFISTO's X11 graphics library*
*Researched: 2026-04-10*
