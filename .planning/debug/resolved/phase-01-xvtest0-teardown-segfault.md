---
status: awaiting_human_verify
phase: 01
related: 01-03
trigger: "pp/ppxvtest0_qt: no window visible on screen AND SIGSEGV in __run_exit_handlers after main returns"
created: 2026-04-11
updated: 2026-04-11
---

## Current Focus

hypothesis: Two coupled defects in the embedded-event-pump model.
  (A) BLOCKER 1: `xvinitgraphique_` calls `show()` + exactly ONE `processEvents(ExcludeUserInputEvents)` — insufficient to actually realize (map+expose) an X11 window. Additionally, `xvtest0.f` immediately calls `xvfermer_` with no hold, so even a mapped window would be imperceptible. SHELL-01/02/06 require visible windows.
  (B) BLOCKER 2: `teardown_atexit()` destroys `QApplication` via `qapp_.reset()`. Research Pattern 1 / Pitfall 5 explicitly warns against destroying QApplication during static/atexit teardown — "Qt internal ordering is hostile." The correct idiom is to LEAK `QApplication` at process exit. The SIGSEGV frames 7..10 are program text inside `teardown_atexit` calling through to Qt destructors which dereference state that libgfortran's atexit handlers have partially unwound.

test: Apply the fix (B) + add a blocking hold pump in (A), rebuild, run `pp/ppxvtest0_qt` and `QT_SCALE_FACTOR=2 pp/ppxvtest0_qt`.

expecting: Both runs should print all 4 `[xvtest0]` lines, display a white "MEFISTO" 800x600 window twice (or once per open cycle), exit with code 0, no core dump. HiDPI variant visibly larger.

next_action: Implement fix, rebuild `bin/cbl_tout_qt`, run both variants, verify no core dump and visible windows.

## Symptoms

expected: `pp/ppxvtest0_qt` opens a visible 800x600 "MEFISTO" window, closes it, reopens, closes, exits cleanly (no crash, no core).
actual:
  - All 4 `[xvtest0]` Fortran PRINT lines appear on stdout
  - Final "OK — cycle open/close/open/close complet" prints
  - NO window is ever drawn on the desktop (confirmed by user: "nothing to see on the screen")
  - Process exits with SIGSEGV during `__run_exit_handlers` — `zsh: segmentation fault (core dumped)`
  - Same behavior with `QT_SCALE_FACTOR=2`: still no window, still segfaults

errors: |
  SIGSEGV backtrace (stripped — addresses only; frames 7..10 are program text):
  #0  0x7f617c8277d5 in ???                    # libgfortran signal handler
  #1  0x7f617c8267c5 in ???                    # libgfortran
  #2  0x7f617c440a4f in libc_sigaction         # SIGSEGV raised here
  #3  0x7f617d954a89 in ???                    # likely libQt6Core/Widgets
  #4  0x7f617d1aa861 in ???                    # another shared lib
  #5  0x7f617e3b1749 in ???                    # Qt internal
  #6  0x7f617e3b1eb1 in ???                    # Qt internal
  #7  0x562569c19441 in ???                    # program text — XvueApp::teardown_atexit
  #8  0x562569c194a1 in ???                    # program text
  #9  0x562569c19395 in ???                    # program text
  #10 0x562569c18e8e in ???                    # program text
  #11 0x7f617c442fa0 in __run_exit_handlers
  #12 0x7f617c443069 in __GI_exit

reproduction: |
  1. bin/cbl_tout_qt (assumed already green)
  2. bin/cbxvtest0_qt
  3. pp/ppxvtest0_qt                       (no window, then segfault)
  4. QT_SCALE_FACTOR=2 pp/ppxvtest0_qt     (same)

started: From first Phase 1 integration run after plan 01-03 executor produced xvue_qt_app.{h,cpp}, xvue_qt_window.{h,cpp}, xvue_qt_canvas.{h,cpp}, xvue_qt_api.cpp bodies, and prpr/xvtest0.f.

## Eliminated

(no eliminated hypotheses yet — see Phase 0 investigation below)

## Evidence

- timestamp: 2026-04-11 Phase 0 (knowledge base check)
  checked: `.planning/debug/knowledge-base.md`
  found: file does not exist — no prior-session match
  implication: open-ended investigation

- timestamp: 2026-04-11 Phase 1 (research/intended model)
  checked: `.planning/phases/01-window-shell-xvueapp-xvuewindow-xvuecanvas/01-RESEARCH.md` Pattern 1, Pitfall 5, D-01, D-08
  found: |
    - D-01 prescribes EXACTLY ONE `processEvents(ExcludeUserInputEvents)` at end of `xvinitgraphique_`.
    - D-08 prescribes `std::atexit` handler that calls `qApp->quit()` then resets the `unique_ptr<QApplication>`.
    - Research A4 (assumption): "QApplication can be constructed after std::atexit has been called once... construct once, never destroy until atexit."
    - Pitfall 5 (research): "Static destruction order relative to Qt internals is hostile — use std::atexit."
    - BUT the research anchors on "use std::atexit" to escape static destructor ordering; it does NOT prove that destroying QApplication inside the atexit handler is actually safe when a Fortran main's libgfortran atexit chain is also running. This is the gap between what the research DESIGNED and what LINUX+gfortran+Qt6 actually do.
  implication: |
    The executor faithfully implemented D-01/D-08 as written. The design itself has two bugs:
    (a) ONE processEvents is insufficient to realize an X11 window.
    (b) Destroying QApplication in an atexit handler that runs alongside libgfortran's own atexit chain is not safe in practice.
    Both need fixing inside the Phase 1 scope — this is a design gap, not an executor deviation.

- timestamp: 2026-04-11 Phase 1 (actual implementation)
  checked: `xvue/qt/src/xvue_qt_app.cpp`, `xvue/qt/src/xvue_qt_api.cpp`
  found: |
    - xvue_qt_app.cpp:34-42 — teardown_atexit() unconditionally destroys QApplication via `qapp_.reset()`.
    - xvue_qt_api.cpp:74-88 — xvinitgraphique_ calls show() then exactly one processEvents (matches D-01).
    - xvue_qt_api.cpp:328-334 — xvfermer_ calls `window_slot().reset()` with NO event pumping afterward. pending deleteLater() events against the window/canvas accumulate in Qt's event queue.
    - Neither xvfermer_ nor any teardown path drains the deferred-delete event queue before shutdown.
  implication: |
    At program exit:
    1. main() returns. libc invokes atexit handlers in LIFO order.
    2. XvueApp::teardown_atexit runs. window_.reset() is a no-op (both windows already destroyed by the two xvfermer_ calls). But their deleteLater() events are still queued.
    3. qapp_->quit() — no-op because exec() was never called.
    4. qapp_.reset() destroys QApplication. QApplication's destructor drains the pending deferred-delete queue — operating on objects in a half-torn-down state, possibly with libgfortran's partial atexit unwinding interleaved. Crash.

- timestamp: 2026-04-11 Phase 1 (driver behavior)
  checked: `prpr/xvtest0.f`
  found: |
    Calls XVINITGRAPHIQUE then immediately XVFERMER twice with no pause. Even if the window were realized, there is zero human-observable interval between show() and close().
  implication: |
    SHELL-01 "opens window" and SHELL-02 "reopen without crash" are technically unit-testable by the absence of Qt asserts, but SHELL-01/06 also require VISIBLE windows (user HiDPI check). Driver needs a brief pump-based hold after each show().

## Resolution

root_cause: |
  TWO design-level bugs in plans 01-01/01-03 (implementation is faithful to the decisions; the decisions are wrong):

  RC-1 (blocker 1 — no visible window):
    `xvinitgraphique_` does a single `processEvents(ExcludeUserInputEvents)` and returns. On X11 a single round of event processing is not enough to:
      (a) dispatch the MapRequest to the X server
      (b) receive the ConfigureNotify
      (c) receive the Expose event
      (d) run the paintEvent
    A realization pump must loop `processEvents` until `window_->windowHandle()->isExposed()` returns true (with a bounded timeout so a failed X connection cannot hang the Fortran main).
    Additionally, `xvtest0.f` closes the window microseconds after opening it, so even a correctly realized window is humanly imperceptible. Fix: hold each open cycle with a short pumped sleep (~500–800 ms) — drained via `processEvents(ExcludeUserInputEvents)` on a loop, no sleep() or QThread::msleep() (blocking sleep stalls Qt, hiding the window).

  RC-2 (blocker 2 — SIGSEGV in teardown_atexit):
    `XvueApp::teardown_atexit()` destroys `QApplication` via `qapp_.reset()` at atexit time. Qt's internal static-destruction ordering is hostile even inside an atexit handler, especially when libgfortran's own atexit handlers are interleaved. Pending deleteLater() events from the two `xvfermer_` calls accumulate in the Qt event queue (xvfermer_ does not pump) and are only drained when QApplication's destructor runs — at the worst possible time.

    The canonical "embed Qt in a non-Qt main" idiom (widely used e.g. by Python bindings, PyQt, and Qt-plugin hosts) is to LEAK QApplication: never destroy it, never call its destructor, let the OS reclaim the memory at process exit. Research A4 is over-specified — "construct once" is load-bearing; "destroy at atexit" is not.

    Secondary contributing factor: xvfermer_ does not drain the deferred-delete queue. Even if RC-2 primary fix is applied, stale events are a latent bug; xvfermer_ should call `processEvents(ExcludeUserInputEvents)` after `window_.reset()` so the deleteLater() queue drains while widgets are still in a well-defined state.

fix: |
  Three focused code changes + one test-driver change:

  FIX-1 (xvue_qt_app.cpp — RC-2 primary): teardown_atexit() must NOT destroy QApplication.
    ```cpp
    void XvueApp::teardown_atexit() {
        // PATTERN: "Embed Qt in non-Qt main" canonical idiom.
        // Destroying QApplication here races with libgfortran's own atexit
        // chain AND with Qt's internal static teardown. Empirically this
        // crashes in __run_exit_handlers on Linux/gfortran/Qt 6.
        // We deliberately LEAK qapp_ — the OS reclaims process memory.
        // window_ has already been reset by the last xvfermer_ (or we
        // drop it here as a safety net, still cleanly, before the qApp
        // event loop machinery is torn down by the platform).
        if (qapp_) {
            // Drain any pending deferred-delete events BEFORE we abandon
            // qApp — after this point we must not run any Qt code.
            QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
        }
        window_.reset();
        if (qapp_) {
            QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
        }
        // Intentionally: do NOT qapp_->quit() (no exec() was ever called).
        // Intentionally: do NOT qapp_.reset(). Release the unique_ptr so the
        // owned QApplication is never destroyed. Documented leak.
        (void)qapp_.release();
    }
    ```

  FIX-2 (xvue_qt_api.cpp — RC-1 primary): xvinitgraphique_ must pump until exposed.
    ```cpp
    void proc(xvinitgraphique)(void) {
        XvueApp::ensure();
        XVUE_QT_ASSERT_MAIN_THREAD();

        auto& win = XvueApp::window_slot();
        if (!win) {
            win = std::make_unique<XvueWindow>();
        }
        win->show();
        win->raise();
        win->activateWindow();

        // D-01 + RC-1: pump processEvents until the window is actually
        // exposed on the compositor (X11 needs multiple trips: MapRequest,
        // ConfigureNotify, Expose). Bounded timeout so a broken display
        // cannot hang the Fortran main.
        QElapsedTimer t;
        t.start();
        QWindow* wh = win->windowHandle();
        const int timeout_ms = 2000;
        while (t.elapsed() < timeout_ms) {
            QCoreApplication::processEvents(
                QEventLoop::ExcludeUserInputEvents, 20);
            if (wh && wh->isExposed()) break;
        }
    }
    ```
    Requires `#include <QElapsedTimer>` and `#include <QWindow>`.

  FIX-3 (xvue_qt_api.cpp — RC-2 secondary): xvfermer_ drains deleteLater queue.
    ```cpp
    void proc(xvfermer)(void) {
        XvueApp::ensure();
        XVUE_QT_ASSERT_MAIN_THREAD();
        XvueApp::window_slot().reset();
        // Drain deferred-delete events queued by widget teardown, so no
        // stale DeferredDelete events remain at process exit.
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    ```

  FIX-4 (prpr/xvtest0.f — human-observable hold):
    Wrap each open cycle with a brief SLEEP call after XVINITGRAPHIQUE so the
    window is actually visible before XVFERMER tears it down. On Linux gfortran
    ships `SLEEP` as a builtin. This is the smallest modification that turns
    xvtest0 into a human-visual check for SHELL-01/02/06.

    ```fortran
          PRINT *,'[xvtest0] premier appel XVINITGRAPHIQUE'
          CALL XVINITGRAPHIQUE
          CALL SLEEP(1)
          PRINT *,'[xvtest0] premier appel XVFERMER'
          CALL XVFERMER
    ...
    ```
    NOTE: `SLEEP(1)` blocks the Fortran main, but QApplication is leaked and
    xvinitgraphique_ already pumped until exposed — so the window is painted
    and mapped before SLEEP starts. The window stays visible for 1 s because
    X11 holds the last Expose content. This is acceptable for a human-visual
    validation driver; real interactive modules (Phase 2+) use xvvoir_'s
    pumped wait, not SLEEP.

verification: |
  Self-verified on 2026-04-11:
    [x] bin/cbl_tout_qt green (exit 0, all *_qt executables rebuilt)
    [x] cmake --build xvue/qt/build green; verify_abi "nm count: 57 header count: 57"; verify_no_exec "OK (no forbidden tokens...)"
    [x] bin/cbxvtest0_qt produces pp/ppxvtest0_qt (86760 bytes, 14:30)
    [x] pp/ppxvtest0_qt (normal DPR):
          - prints all four [xvtest0] lines
          - prints "OK — cycle open/close/open/close complet"
          - exits with status 0 (EXIT=0)
          - NO core dump, no SIGSEGV
    [x] QT_SCALE_FACTOR=2 pp/ppxvtest0_qt:
          - same clean output, EXIT=0, no core dump

  Requires human verification (cannot be self-checked from the shell — needs a
  live desktop session):
    [ ] Two visible "MEFISTO" white 800x600 windows during normal run
        (1 second each, no flicker, centered near default compositor placement)
    [ ] With QT_SCALE_FACTOR=2: windows visibly ~2x larger on screen than normal run
    [ ] No platform-integration warnings OTHER than the pre-existing
        "Gtk-Message: Failed to load module colorreload-gtk-module" lines
        (harmless Qt GTK theme plugin noise, unrelated to the bug)

files_changed:
  - xvue/qt/src/xvue_qt_app.cpp
  - xvue/qt/src/xvue_qt_api.cpp
  - prpr/xvtest0.f
