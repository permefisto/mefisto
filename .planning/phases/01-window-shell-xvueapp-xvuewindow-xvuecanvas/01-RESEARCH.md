# Phase 1: Window shell (`XvueApp` + `XvueWindow` + `XvueCanvas`) - Research

**Researched:** 2026-04-11
**Domain:** Qt 6 Widgets — application singleton, top-level window, central widget, screen metrics, embedded (non-`exec()`) event pumping, Fortran↔C++ FFI boundary
**Confidence:** HIGH

## Summary

Phase 1 replaces the Phase 0 warn-once stubs for seven entry points (`xvinitgraphique_`, `xvfermer_`, `xvpxecran_`, `xvmmecran_`, `xvfond_`, `xvinfo_`-partial, plus the SHELL-07 macro retrofit across all stubs) with real Qt 6 code. The work is almost entirely prescribed by `01-CONTEXT.md` — this research document maps each locked decision to the exact Qt 6 API surface, the existing Phase 0 artifacts, the legacy `xvuelc.c` semantics the Qt version must match, and the pitfalls that apply.

The technical core is four Qt idioms, all stable and well-documented: (1) `std::call_once` + fabricated `argc/argv` + `std::atexit` to manage a process-lifetime `QApplication` embedded in a non-Qt (Fortran) main; (2) `QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents)` as the one-shot realization pump in place of a forbidden top-level `exec()`; (3) `QGuiApplication::primaryScreen()->size()` (logical pixels) and `->physicalSize()` (millimetres) for HiDPI-correct screen metrics; and (4) a minimal `QMainWindow` + central `QWidget` subclass whose `paintEvent` does exactly `fillRect(rect(), state_->background_)`. None of these require Qt 6 features beyond the base `qt6-base-dev` package already installed (6.10.2 verified).

**Primary recommendation:** Follow the 21 decisions in CONTEXT.md literally. Resist any temptation to add scope (no drawing, no backing `QPixmap`, no menu bar content, no `xvinfo_` palette work). The success criterion is the `prpr/xvtest0.f` reopen cycle — keep every file edit aimed at making that single driver pass.

## User Constraints (from CONTEXT.md)

### Locked Decisions

**Window open timing and sizing (D-01 to D-05):**
- D-01: `xvinitgraphique_` is the single lifecycle entry point: `std::call_once` → ensure `QApplication` → if `window_` null allocate fresh `XvueWindow` → `window_->show()` → one call to `QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents)` to realize on the compositor. `ExcludeUserInputEvents` is **mandatory** (Pitfall 6/8).
- D-02: Initial window size **800×600 logical pixels**, title **`"MEFISTO"`**. `XvueWindow` constructor sets both and sizes the central `XvueCanvas`.
- D-03: `xvinfo_` is a **partial real implementation**: calls `window_->resize(*ix, *iy)`, zeroes/defaults all other output parameters, prints one-shot `"xvue-qt: stub xvinfo_ palette outputs not implemented yet"`, returns. Palette/font outputs deferred to Phase 3.
- D-04: `XvueState` struct holds exactly **one field**: `QColor background_ = Qt::black;` (matches legacy `BlackPixel` at `xvue/xvuelc.c:935`). No other state in Phase 1.
- D-05: `XvueCanvas::paintEvent` body is exactly `QPainter(this).fillRect(rect(), state_->background_)`. One line. Phase 2 swaps it to `drawPixmap(0, 0, *backing_)`.

**`xvfermer_` and reopen semantics (D-06 to D-09):**
- D-06: `xvfermer_` does `window_.reset()` on the `std::unique_ptr<XvueWindow>` owned by `XvueApp`. Does **not** touch `QApplication`, does **not** call `qApp->quit()`, does **not** emit events. `XvueState` and (Phase 2+) backing pixmap die with the window — each reopen is a fresh session.
- D-07: Second `xvinitgraphique_` after `xvfermer_` observes `window_ == nullptr`, allocates fresh `XvueWindow`. `call_once` flag is already set — `QApplication` is not recreated.
- D-08: `QApplication` torn down by `std::atexit` handler installed inside `XvueApp::ensure()` on first call: `qApp->quit()` then reset the unique_ptr. No destructor at static-destruction time.
- D-09: `std::call_once` guards **only** the `QApplication` construction, not the window. Window is guarded by plain `if (!window_)` check. Load-bearing for the reopen story.

**SHELL-03 `exec()` ban (D-10, D-11):**
- D-10: New CMake custom target `verify_no_exec` runs as `POST_BUILD` on the `xvueqt` target, greps `xvue/qt/src` and `xvue/qt/include` for `QApplication::exec` or `qApp->exec()`, fails build on any match. Mirrors Phase 0 `verify_abi` pattern.
- D-11: **No** git pre-commit hook. CMake guard is sole enforcement.

**SHELL-07 thread assertion (D-12, D-13):**
- D-12: `XVUE_QT_ASSERT_MAIN_THREAD()` macro body fleshed out to real `Q_ASSERT(QThread::currentThread() == qApp->thread())` in debug, empty in release. **Bulk-inserted into every one of the ~57 existing stubs** as first statement of each function body, ahead of the `warn_once` line. Single-file bulk edit.
- D-13: Macro is null-guarded: `if (qApp) Q_ASSERT(...)` — the null guard is required because `xvpxecran_`/`xvmmecran_` are callable before `XvueApp::ensure()` runs on the very first call. Real Phase 1 implementations start with `XvueApp::ensure();` as the first statement, **before** the macro.

**`xvfond_` storage (D-14, D-15):**
- D-14: `xvfond_(int *icolor)` stores into `XvueState::background_` via a minimal 2-entry table: `0 → Qt::black`, `1 → Qt::white` (matches legacy `BlackPixel`/`WhitePixel`). Other values fall through to `Qt::black` with one-shot stderr warn. Disappears cleanly when Phase 3 wires the real palette.
- D-15: `XvueCanvas` holds a raw pointer to `XvueState` owned by `XvueWindow`. `xvfond_` ends with `canvas_->update()` to schedule a repaint.

**SHELL-04 screen metrics (D-16 to D-18):**
- D-16: `xvpxecran_(int *xp, int *yp)` → `QGuiApplication::primaryScreen()->size().width()/.height()` (logical pixels). Multi-monitor awareness explicitly deferred.
- D-17: `xvmmecran_(int *xmm, int *ymm)` → `QGuiApplication::primaryScreen()->physicalSize().width()/.height()` rounded to int (Qt 6 guarantees millimetres regardless of DPR).
- D-18: Both callable before `xvinitgraphique_` creates the window. Every real Phase 1 entry starts with `XvueApp::ensure();` **before** `XVUE_QT_ASSERT_MAIN_THREAD()` to guarantee `qApp` is non-null.

**Phase 1 test driver (D-19 to D-21):**
- D-19: New `prpr/xvtest0.f`, body: `CALL XVINITGRAPHIQUE; CALL XVFERMER; CALL XVINITGRAPHIQUE; CALL XVFERMER; STOP END`. Exercises SHELL-01, SHELL-02, and the `QApplication`-stays-alive invariant in the smallest unit. Existing `prpr/xvtest{1..4}.f` untouched.
- D-20: New `bin/cbxvtest0_qt` cloned from `bin/cbmail_qt` (thinnest Phase 0 clone), adjusting source file / exe name / linker line only. Produces `pp/ppxvtest0_qt`. **`bin/cbl_tout_qt` is NOT modified** — keeps the "five canonical modules" contract.
- D-21: Validation gate is: `bin/cbxvtest0_qt` builds cleanly, `pp/ppxvtest0_qt` runs without crash, human observes the blank 800×600 `"MEFISTO"` window **twice** (once per `XVINITGRAPHIQUE`). 5 baseline `testa/` cases **not** re-run in Phase 1 — resumes at Phase 2.

### Claude's Discretion

- Exact implementation of `XvueApp::ensure()` (function vs. static method, header/source split, member ordering)
- `XvueWindow` parent class specifics — simplest form is bare `QMainWindow` with only central widget set
- `xvtest0.f` exact line count — must do 2×(init, fermer) at minimum, extra `PRINT *,'...'` lines are OK
- `verify_no_exec` exact invocation — inlined CMake command, shell script helper, or `execute_process` — planner picks
- Header/source file fanout inside `xvue/qt/` — may follow the 4-file split (`xvue_qt_app`, `_window`, `_canvas`, `_state`) literally, or collapse; hard rule is no code leaks outside `xvue/qt/include/` + `xvue/qt/src/`
- Pass `QT_SCALE_FACTOR` from shell vs. set `Qt::AA_EnableHighDpiScaling` in code — Qt 6 default is on; defensive attribute call is optional

### Deferred Ideas (OUT OF SCOPE)

- Multi-monitor awareness for `xvpxecran_`/`xvmmecran_` (primaryScreen only)
- `QPalette`-based background handling (rejected — Phase 2 pixmap would cover it)
- Git pre-commit hook for `exec` guard (rejected — CMake only)
- Adding `pp/ppxvtest0_qt` to `bin/cbl_tout_qt`
- Defensive `Qt::AA_EnableHighDpiScaling` attribute (Qt 6 default)
- Full `XvueState` (pen/brush/font/palette/painter/pixmaps) — Phases 2/3/4
- `WA_DeleteOnClose` + `close()` for `xvfermer_` (rejected — async, non-deterministic)
- Running 5 baseline `testa/` cases in Phase 1
- `xvtest0.f` retirement post-Phase 1 (stays as permanent shell sanity driver)

## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| SHELL-01 | `xvinitgraphique_` opens `QMainWindow` with `XvueCanvas` central widget via `std::call_once` + static fake argc/argv | Standard Stack §Core, Architecture §Pattern 1 (Singleton), Code Examples §1 |
| SHELL-02 | `xvfermer_` closes window, does NOT destroy `QApplication`; reopen must not trigger "there can only be one" | Architecture §Pattern 2 (Reopen), Pitfall 1, Code Examples §2 |
| SHELL-03 | `QApplication::exec()` never called; enforced by CMake grep | Architecture §Pattern 3 (Embedded event pump), Code Examples §6 |
| SHELL-04 | `xvpxecran_`/`xvmmecran_` return logical px and physical mm from `QScreen` | Standard Stack §Core, Code Examples §3/§4, State of the Art (Qt 6 HiDPI defaults) |
| SHELL-05 | `xvfond_` sets background without corrupting backing pixmap | D-14/D-15; no backing pixmap exists yet in Phase 1, so trivially satisfied |
| SHELL-06 | Renders correctly at DPR > 1 (HiDPI) without size drift vs X11 | State of the Art §HiDPI, Pitfall 2 |
| SHELL-07 | Every `extern "C"` entry has debug-build `Q_ASSERT(thread == qApp->thread())` | Code Examples §7 (macro retrofit), Pitfall 3 |

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Qt6 Core | 6.10.2 [VERIFIED: `pkg-config --modversion Qt6Core`] | `QCoreApplication`, `QThread`, `QEventLoop`, `QObject` base | Already found via Phase 0 `find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets PrintSupport)` — nothing new to add |
| Qt6 Gui | 6.10.2 [VERIFIED] | `QGuiApplication`, `QScreen`, `QPainter`, `QColor`, `QPaintEvent` | Supplies the `primaryScreen()` API used by `xvpxecran_`/`xvmmecran_` |
| Qt6 Widgets | 6.10.2 [VERIFIED] | `QApplication`, `QMainWindow`, `QWidget` | Supplies the `QApplication` singleton (a `QGuiApplication` subclass) and `QMainWindow`/`QWidget` used by `XvueWindow`/`XvueCanvas` |
| Qt6 PrintSupport | 6.10.2 [VERIFIED] | Unused in Phase 1 — already linked by Phase 0 for future `xvpostscript_` (PDF export, Phase 7) | Keep in link line (matches `bin/cbmail_qt`) |
| C++ standard library (libstdc++) | GCC 15.2.0 [VERIFIED: `gfortran --version` → `GNU Fortran (Debian 15.2.0-11)`] | `std::call_once`, `std::once_flag`, `std::unique_ptr`, `std::atexit` (via `<cstdlib>`) | Standard C++17 idioms — `CMAKE_CXX_STANDARD 17` already set in Phase 0 |

**Installation:**
```bash
# Already installed — verified during research:
# qt6-base-dev 6.10.2+dfsg-6
# qt6-base-dev-tools 6.10.2+dfsg-6
# No new packages required for Phase 1.
```

**Version verification (2026-04-11):**
- `pkg-config --modversion Qt6Core` → `6.10.2` [VERIFIED]
- `dpkg -l | grep qt6-base-dev` → `qt6-base-dev 6.10.2+dfsg-6 amd64` [VERIFIED]
- `gfortran --version` → `GNU Fortran (Debian 15.2.0-11) 15.2.0` [VERIFIED]

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `<mutex>` | C++17 | `std::call_once` + `std::once_flag` for `XvueApp::ensure()` | First and only invocation site |
| `<memory>` | C++17 | `std::unique_ptr<QApplication>`, `std::unique_ptr<XvueWindow>` | RAII ownership inside `XvueApp` |
| `<cstdlib>` | C++17 | `std::atexit` for process-exit teardown | Registered once on first `XvueApp::ensure()` call |

### Alternatives Considered

| Instead of | Could Use | Tradeoff | Decision |
|------------|-----------|----------|----------|
| `std::call_once` | Static local `QApplication` inside a function | Simpler but less explicit; harder to add the `atexit` registration on the same code path | **Locked to `call_once`** by ARCHITECTURE.md §Singleton discipline and PITFALLS.md §Pitfall 5 |
| `std::unique_ptr<XvueWindow>` | Raw pointer + manual `delete` | Raw pointer loses RAII, leaks on exception | unique_ptr (standard C++17) |
| `WA_DeleteOnClose` + `close()` | Qt-idiomatic window-close pattern | Async deletion (posted to event queue), non-deterministic relative to `xvfermer_` return | **Rejected in D-06 deferred list** — `reset()` is synchronous and simpler |
| `QCoreApplication::processEvents()` (no flag) | Default-flag processEvents | Would run user input events → re-entrancy risk per Pitfall 6/8 | **Locked to `ExcludeUserInputEvents`** by D-01 |
| Singleton free function | Class with static method | Equivalent; class groups ownership + accessors | Planner discretion |

## Architecture Patterns

### Recommended Project Structure

```
xvue/qt/
├── CMakeLists.txt                    # Phase 0 file; grows by new sources + verify_no_exec target
├── .gitignore
├── cmake/
│   ├── verify_abi.sh                  # Phase 0 — unchanged
│   └── verify_no_exec.sh              # NEW Phase 1 — optional helper for the grep guard
├── include/
│   └── xvue_qt_api.h                   # Phase 0 — Phase 1 fleshes out the XVUE_QT_ASSERT_MAIN_THREAD body
└── src/
    ├── xvue_qt_api.cpp                 # Phase 0 — Phase 1 rewrites 7 SHELL entry bodies + retrofits all stubs
    ├── xvue_qt_app.{h,cpp}             # NEW — XvueApp singleton, argc/argv fab, atexit, ensure()
    ├── xvue_qt_window.{h,cpp}          # NEW — XvueWindow : QMainWindow
    ├── xvue_qt_canvas.{h,cpp}          # NEW — XvueCanvas : QWidget
    └── xvue_qt_state.{h,cpp}           # NEW — XvueState (single QColor field in Phase 1)
```

The 4-new-file split mirrors ARCHITECTURE.md §"Recommended Project Structure." The planner may collapse any pair(s) at discretion; the only hard rule is that no new file is added under `xvue/qt/include/` (that directory remains the single ABI surface — `xvue_qt_api.h` only). Internal `*.h` files for the new components live alongside their `.cpp` under `xvue/qt/src/`.

### Pattern 1: `QApplication` singleton via `std::call_once` + fabricated argc/argv

**What:** Ensure exactly one `QApplication` exists for process lifetime. Construct lazily on first need. Never re-create. Tear down only at `std::atexit`.

**When to use:** Every C++ entry point that needs Qt — must call `XvueApp::ensure()` as its first statement, before `XVUE_QT_ASSERT_MAIN_THREAD()`.

**Example:**
```cpp
// xvue/qt/src/xvue_qt_app.h
// Source: ARCHITECTURE.md §"Singleton discipline for QApplication"
//         PITFALLS.md §"Pitfall 5: QApplication double-initialization"
//         01-CONTEXT.md D-01, D-07, D-08, D-09
#pragma once
#include <memory>
#include <mutex>
#include <QApplication>

class XvueWindow;

class XvueApp {
public:
    static void  ensure();                 // idempotent; safe to call from any entry point
    static QApplication* qapp();           // non-null once ensure() has run
    static XvueWindow*&  window_slot();    // reference to the unique_ptr slot for lazy realloc
private:
    static std::once_flag             once_flag_;
    static std::unique_ptr<QApplication> qapp_;
    static std::unique_ptr<XvueWindow>   window_;
    static void teardown_atexit();
};
```

```cpp
// xvue/qt/src/xvue_qt_app.cpp
// Source: ARCHITECTURE.md §"Singleton discipline" + PITFALLS.md §Pitfall 5
#include "xvue_qt_app.h"
#include <cstdlib>

std::once_flag                       XvueApp::once_flag_;
std::unique_ptr<QApplication>        XvueApp::qapp_;
std::unique_ptr<XvueWindow>          XvueApp::window_;

void XvueApp::ensure() {
    std::call_once(once_flag_, []{
        // Static storage so Qt can safely hold references for app lifetime.
        static int   fake_argc = 1;
        static char  arg0[] = "mefisto";
        static char* fake_argv[] = { arg0, nullptr };
        qapp_ = std::make_unique<QApplication>(fake_argc, fake_argv);
        std::atexit(&XvueApp::teardown_atexit);
    });
}

QApplication* XvueApp::qapp() { return qapp_.get(); }
std::unique_ptr<XvueWindow>& XvueApp_window_slot_ref();  // (helper; or implement via friend / inline)

void XvueApp::teardown_atexit() {
    if (qapp_) qapp_->quit();
    window_.reset();
    qapp_.reset();
}
```

**Key points:**
- `std::call_once` guards **only** the `QApplication` construction (D-09). The `XvueWindow` is guarded by a plain `if (!window_)` inside `xvinitgraphique_`.
- `fake_argc`/`fake_argv` are `static` — their storage lives for the whole process, which Qt expects.
- `std::atexit` handler is the single teardown site. Do **not** rely on `~XvueApp()` at static-destruction time (Qt internal ordering is hostile — PITFALLS.md Pitfall 5).

### Pattern 2: Lazy `XvueWindow` reopen

**What:** A fresh `XvueWindow` is allocated every time `xvinitgraphique_` is called with a null window slot. `xvfermer_` resets the slot to null. The `QApplication` is orthogonal.

**Example:**
```cpp
// xvue/qt/src/xvue_qt_api.cpp (rewritten bodies for SHELL-01 + SHELL-02)
// Source: 01-CONTEXT.md D-01, D-06, D-07
#include "xvue_qt_api.h"
#include "xvue_qt_app.h"
#include "xvue_qt_window.h"
#include <QCoreApplication>
#include <QEventLoop>

extern "C" void proc(xvinitgraphique)(void) {
    XvueApp::ensure();                    // D-18: MUST be first, before the assertion
    XVUE_QT_ASSERT_MAIN_THREAD();

    auto& win = XvueApp::window_slot();   // reference to the unique_ptr<XvueWindow> slot
    if (!win) {
        win = std::make_unique<XvueWindow>();   // D-02: 800x600 "MEFISTO" set in ctor
    }
    win->show();

    // D-01: single realization pump, ExcludeUserInputEvents is mandatory (Pitfall 6/8)
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

extern "C" void proc(xvfermer)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    // D-06: destroy the window only; do NOT touch qApp, do NOT call qApp->quit()
    XvueApp::window_slot().reset();
}
```

### Pattern 3: Embedded event pump without `exec()`

**What:** Qt's default "app owns the event loop" is inverted: Fortran owns the main loop, Qt services events only when an entry point voluntarily pumps them. Phase 1 has exactly one pump call (at the end of `xvinitgraphique_`, D-01). Future phases add a pump in `xvvoir_` (Phase 2) and nested `QEventLoop` in `xvsouris_` (Phase 5).

**Rule:** No code in `xvue/qt/src/*.cpp` or `xvue/qt/include/*.h` ever writes the literal strings `QApplication::exec` or `qApp->exec()`. Enforced by `verify_no_exec` post-build CMake target (D-10).

### Pattern 4: Screen metrics from `QScreen`

**What:** Qt 6's `QGuiApplication::primaryScreen()` returns a `QScreen*` whose `size()` is logical pixels and whose `physicalSize()` is millimetres (as `QSizeF`, regardless of device pixel ratio).

**Example:**
```cpp
// Source: Qt 6 docs (doc.qt.io/qt-6/qscreen.html) + 01-CONTEXT.md D-16, D-17
#include <QGuiApplication>
#include <QScreen>

extern "C" void proc(xvpxecran)(int *xp, int *yp) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    QScreen* s = QGuiApplication::primaryScreen();
    *xp = s->size().width();    // logical pixels (DPR-aware)
    *yp = s->size().height();
}

extern "C" void proc(xvmmecran)(int *xmm, int *ymm) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    QScreen* s = QGuiApplication::primaryScreen();
    QSizeF mm = s->physicalSize();
    *xmm = static_cast<int>(mm.width() + 0.5);  // round to int
    *ymm = static_cast<int>(mm.height() + 0.5);
}
```

### Anti-Patterns to Avoid

- **Destroying `QApplication` inside `xvfermer_`**: Qt-hostile on all platforms — recreating after destroy crashes platform plugins (PITFALLS.md Pitfall 5, ARCHITECTURE.md §Anti-Pattern 4).
- **Calling `QApplication::exec()` from `xvinitgraphique_`**: Inverts control flow; Fortran never regains control (PITFALLS.md Pitfall 6, ARCHITECTURE.md §Anti-Pattern 2).
- **`processEvents()` without `ExcludeUserInputEvents`**: Opens the door to re-entrant entry-point calls before Phase 5's nested `QEventLoop` pattern is in place (PITFALLS.md Pitfalls 6, 7, 8).
- **Adding `QPalette` background plumbing instead of `XvueState::background_`**: Phase 2's backing pixmap will cover the palette, wasting the work (Deferred).
- **Destructor-time teardown of `QApplication`**: Static-destruction order relative to Qt internals is hostile — use `std::atexit` (PITFALLS.md Pitfall 5, ARCHITECTURE.md).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Process-wide app singleton | Custom mutex + static pointer + manual double-check | `std::call_once` + `std::once_flag` | C++17 standard idiom, thread-safe, zero corner cases |
| Fake argc/argv for embedded Qt | Heap-allocated `char**` per call | `static int argc; static char* argv[]` at file scope | Qt caches argv pointers; static storage outlives the app trivially |
| Background color storage | Parallel `QPalette` management | `QColor` field inside `XvueState` | D-04 locks this; palette would be stomped by Phase 2 backing pixmap anyway |
| Process-exit cleanup of Qt | `~XvueApp()` at static destruction | `std::atexit` registered on first `ensure()` | Static destruction order vs. Qt internals is unsafe (Pitfall 5) |
| Window-close detection for `xvfermer_` | Event filter on `QCloseEvent` | `unique_ptr::reset()` directly | D-06 makes teardown synchronous; no event round-trip needed |
| Thread affinity check on every entry | Hand-rolled `pthread_self()` tracking | `Q_ASSERT(QThread::currentThread() == qApp->thread())` (macro-wrapped) | Qt's canonical check; zero runtime cost in release (D-12) |
| HiDPI scaling toggle | `Qt::AA_EnableHighDpiScaling` attribute calls | Nothing — Qt 6 enables it by default | Qt 6 removed the attribute (it's a no-op); default behavior covers SHELL-06 |
| `exec()` grep guard | Git pre-commit hook with install dance | CMake `add_custom_command(TARGET xvueqt POST_BUILD ...)` | Single enforcement point, cannot be bypassed with `--no-verify` (D-10/D-11) |
| Screen dimensions | `XRandR`/`/sys/class/drm` parsing | `QGuiApplication::primaryScreen()->size()` + `physicalSize()` | Qt abstracts X11/Wayland/platform differences; HiDPI-correct out of the box |

**Key insight:** Phase 1 is the cheapest possible code that satisfies seven requirements. Every decision has been made to pick the Qt-canonical idiom and reject any hand-rolled alternative. The planner's job is to translate CONTEXT.md decisions into ordered file edits, not to re-litigate architecture.

## Common Pitfalls

### Pitfall 1: `QApplication` double-init on reopen
**What goes wrong:** Second call to `xvinitgraphique_` triggers `QApplication: there can only be one QApplication instance` and aborts.
**Why it happens:** Naïve implementation creates a fresh `QApplication` every call, or destroys it inside `xvfermer_` and tries to recreate later.
**How to avoid:** `std::call_once` (D-09) guards `QApplication` construction; `xvfermer_` only destroys the window (D-06). The `prpr/xvtest0.f` driver (D-19) exercises the exact failure mode with two init/fermer cycles — if Phase 1 passes that driver, Pitfall 1 is defeated.
**Warning signs:** Any code path in `xvue_qt_api.cpp` that touches `qapp_` outside `XvueApp::ensure()` or `teardown_atexit()`.
**Source:** PITFALLS.md §Pitfall 5 [CITED]; CONTEXT.md D-06, D-07, D-08, D-09.

### Pitfall 2: HiDPI size drift between X11 and Qt backends
**What goes wrong:** On a 4K monitor or with `QT_SCALE_FACTOR=2`, the Qt window renders at 2× the X11 size (or ½), violating SHELL-06.
**Why it happens:** Mixing physical pixels (what X11 returned) with logical pixels (what Qt returns by default) when `xvpxecran_`/`xvmmecran_` feed into Fortran sizing code.
**How to avoid:** Both `QScreen::size()` and `QWidget::resize(w, h)` use **logical** pixels in Qt 6 by default. Do not set `Qt::AA_DisableHighDpiScaling`. Do not apply any `devicePixelRatio` correction in Phase 1. Validate by running `pp/ppxvtest0_qt` under `QT_SCALE_FACTOR=2` and confirming the window is visibly 800×600 "logical."
**Warning signs:** Any `devicePixelRatio()` or `devicePixelRatioF()` call anywhere in Phase 1 code.
**Source:** Qt 6 HiDPI docs (doc.qt.io/qt-6/highdpi.html) [CITED]; SHELL-06.

### Pitfall 3: `qApp` null on first call to `xvpxecran_`/`xvmmecran_`
**What goes wrong:** A Fortran driver that calls `xvpxecran_` before `xvinitgraphique_` (legal in the X11 backend — `display_mef` was opened lazily) hits a null `qApp` when `XVUE_QT_ASSERT_MAIN_THREAD()` dereferences it.
**Why it happens:** The legacy X11 path handled this implicitly via `XOpenDisplay` being called in every metric function's prelude. Qt has no such idiom.
**How to avoid:** Every real Phase 1 implementation starts with `XvueApp::ensure();` as its first statement, **before** `XVUE_QT_ASSERT_MAIN_THREAD()` (D-18). The macro itself is null-guarded (`if (qApp)`) per D-13.
**Warning signs:** Any entry-point body that calls the macro before `ensure()`.
**Source:** CONTEXT.md D-13, D-18.

### Pitfall 4: `processEvents()` re-entering graphics entry points
**What goes wrong:** The realization pump at the end of `xvinitgraphique_` dispatches a queued mouse click back into `xvsouris_`, which nests a `QEventLoop`, which re-enters `processEvents`… stack corruption or deadlock.
**Why it happens:** Default `processEvents()` processes all pending events including user input.
**How to avoid:** **Always** use `QEventLoop::ExcludeUserInputEvents` in the pump call (D-01). Phase 5 will add separate nested-loop plumbing for the blocking reads; Phase 1 must not open that door prematurely.
**Warning signs:** Any `processEvents()` call without the `ExcludeUserInputEvents` flag.
**Source:** PITFALLS.md §Pitfalls 6, 7, 8 [CITED]; Qt 6 docs (doc.qt.io/qt-6/qcoreapplication.html#processEvents) [CITED].

### Pitfall 5: Legacy stub shadowing real implementation
**What goes wrong:** After rewriting `xvinitgraphique_` to a real body, the old warn-once line is left in place or the `static bool warned = false;` is left dangling; developer sees stderr spam and concludes the implementation didn't take.
**Why it happens:** Phase 0's stub pattern has a consistent 3-line template; rewriting requires touching all three lines and the function body.
**How to avoid:** For the 7 real entries, fully delete the `warn_once` block and replace with the real body. For the `xvinfo_` partial (D-03), keep **one** `warn_once` line with the adjusted text. For all other ~50 stubs, leave them untouched **except** for the SHELL-07 macro retrofit (first line of each body, ahead of existing `warn_once`).
**Warning signs:** Grepping `xvue/qt/src/xvue_qt_api.cpp` for `warn_once` after Phase 1 should yield ~50 matches (50 stubs minus 7 real implementations plus the `xvinfo_` partial; final count ~51).
**Source:** Phase 0 stub convention; CONTEXT.md D-12, D-17.

### Pitfall 6: Adding a menu bar, status bar, or toolbar in `XvueWindow`
**What goes wrong:** Code drift into Phase 6 territory; the "Phase 1 window shell" becomes entangled with chrome decisions that the per-module lexicon audit hasn't been done for yet.
**Why it happens:** `QMainWindow` offers `setMenuBar()`, `statusBar()`, `addToolBar()` one-liner conveniences. Tempting to sprinkle them in.
**How to avoid:** `XvueWindow` in Phase 1 is a bare `QMainWindow` with `setCentralWidget(canvas_)` and nothing else. No `setMenuBar`, no `statusBar()`, no `addToolBar`, no dock widgets.
**Warning signs:** Any Qt API call inside `xvue_qt_window.cpp` other than: ctor (`QMainWindow(parent)`), `setWindowTitle`, `resize`, `setCentralWidget`.
**Source:** CONTEXT.md Domain §"Explicitly out of scope"; ROADMAP.md Phase 6.

### Pitfall 7: CMake AUTOMOC skipping `xvue_qt_window.h` / `xvue_qt_canvas.h`
**What goes wrong:** `XvueCanvas : public QWidget` with `Q_OBJECT` macro in the header doesn't get moc-ed; link fails with `undefined reference to XvueCanvas::staticMetaObject`.
**Why it happens:** AUTOMOC only scans sources listed in the `add_library` / `add_executable` call. If the new `.cpp` files aren't added to the `xvueqt` target, AUTOMOC ignores their headers.
**How to avoid:** Every new `.cpp` under `xvue/qt/src/` must be added to `add_library(xvueqt STATIC ...)` in `xvue/qt/CMakeLists.txt`. Phase 0 D-11 already set `CMAKE_AUTOMOC ON` before `find_package(Qt6)` (Phase 0 Pitfall 9 guards this).
**Warning signs:** Link errors mentioning `staticMetaObject`, `qt_metacall`, or `moc_xvue_qt_*.cpp.o: No such file`.
**Source:** PITFALLS.md §Pitfall 9 [CITED]; Qt 6 CMake docs.

### Pitfall 8: Forgetting `Q_OBJECT` macro when declaring `XvueWindow`/`XvueCanvas`
**What goes wrong:** Qt meta-object system features (signals/slots, `QThread::currentThread() == qApp->thread()` checks via `QObject::thread()`) silently misbehave; Phase 5's event filter won't install.
**Why it happens:** Any class deriving from `QObject` (directly or via `QMainWindow`/`QWidget`) needs `Q_OBJECT` to get its `staticMetaObject` slot.
**How to avoid:** Both `XvueWindow` (derives `QMainWindow`) and `XvueCanvas` (derives `QWidget`) get `Q_OBJECT` in their class body, even though Phase 1 declares no signals or slots. This primes the hook for Phase 5 without retrofits.
**Warning signs:** Missing `Q_OBJECT` in a header that inherits any Qt class.
**Source:** Qt 6 QObject docs; PITFALLS.md §Pitfall 9.

### Pitfall 9: OpenMP-flag propagation into `libxvueqt.a`
**What goes wrong:** Phase 1 `.cpp` files pick up `-fopenmp` from somewhere, ELF-linking collides with the `_OMP` Fortran executables.
**Why it happens:** GCC propagates flags across directories if CMake inherits globals.
**How to avoid:** Phase 0 D-11 already set `target_compile_options(xvueqt PRIVATE -fno-openmp -Wall -Wextra)`. Phase 1 does not change this line. All new `.cpp` files inherit the same block because they're compiled as part of the `xvueqt` target.
**Warning signs:** Any `add_compile_options(-fopenmp)` or `set(CMAKE_CXX_FLAGS ... -fopenmp ...)` appearing in `xvue/qt/CMakeLists.txt`.
**Source:** PITFALLS.md §Pitfall 11 [CITED]; Phase 0 CMakeLists.txt:36.

### Pitfall 10: Accidentally breaking the X11 build path
**What goes wrong:** Editing `xvue/xvuelc.c` or touching a shared `bin/cb*` script inadvertently breaks `bin/cbl_tout`, violating CLAUDE.md's build-never-breaks rule.
**Why it happens:** Phase 1 has zero legitimate reason to edit `xvue/xvuelc.c`, but a misdirected refactor might.
**How to avoid:** Touch **nothing** under `xvue/` except inside `xvue/qt/`. Phase 0 D-02 established `xvue/qt/` as the read-only-outside-scope boundary; Phase 1 inherits it unchanged.
**Warning signs:** `git status` after Phase 1 work should show only files under `xvue/qt/`, new `prpr/xvtest0.f`, new `bin/cbxvtest0_qt`, and the one planning file.
**Source:** CLAUDE.md §"Compilation must never break"; Phase 0 D-02.

### Pitfall 11: `QT_QPA_PLATFORM` missing on the test bench
**What goes wrong:** Running `pp/ppxvtest0_qt` in a headless or misconfigured X environment throws `qt.qpa.plugin: Could not load the Qt platform plugin "xcb"`.
**Why it happens:** Phase 1 requires a real display — the validation gate (D-21) is **human observation** of a visible blank window. No `offscreen` fallback is attempted.
**How to avoid:** Run Phase 1 validation on a desktop session, with `DISPLAY` set (X11) or `WAYLAND_DISPLAY` set (Wayland). Document the requirement in any Phase 1 README note. Do **not** attempt `QT_QPA_PLATFORM=offscreen` — it bypasses the very rendering path SHELL-06 validates.
**Warning signs:** `cannot connect to X server` or `Could not load platform plugin` errors.
**Source:** Qt 6 platform plugin docs.

## Code Examples

Verified patterns. All match Qt 6 (6.10.2) stable API.

### 1. `XvueApp::ensure()` with `std::call_once`, fake argc/argv, `std::atexit`

```cpp
// xvue/qt/src/xvue_qt_app.cpp — FULL reference shape (planner may restructure)
// Source: ARCHITECTURE.md §"Singleton discipline", PITFALLS.md §Pitfall 5, CONTEXT.md D-08/D-09
#include "xvue_qt_app.h"
#include "xvue_qt_window.h"
#include <QApplication>
#include <cstdlib>
#include <memory>
#include <mutex>

namespace {
    std::once_flag                       g_once;
    std::unique_ptr<QApplication>        g_qapp;
    std::unique_ptr<XvueWindow>          g_window;

    void teardown_at_exit() {
        // D-08: qApp->quit() then reset unique_ptr. No static destructor reliance.
        if (g_qapp) g_qapp->quit();
        g_window.reset();
        g_qapp.reset();
    }
}

void XvueApp::ensure() {
    std::call_once(g_once, []{
        static int   fake_argc = 1;
        static char  arg0[] = "mefisto";
        static char* fake_argv[] = { arg0, nullptr };
        g_qapp = std::make_unique<QApplication>(fake_argc, fake_argv);
        std::atexit(&teardown_at_exit);
    });
}

QApplication* XvueApp::qapp() { return g_qapp.get(); }
std::unique_ptr<XvueWindow>& XvueApp::window_slot() { return g_window; }
```

### 2. `XvueWindow` (minimal `QMainWindow`)

```cpp
// xvue/qt/src/xvue_qt_window.h
// Source: CONTEXT.md D-02, D-04, D-15
#pragma once
#include <QMainWindow>
#include "xvue_qt_state.h"

class XvueCanvas;

class XvueWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit XvueWindow(QWidget* parent = nullptr);
    ~XvueWindow() override;
    XvueState* state() { return &state_; }
    XvueCanvas* canvas() { return canvas_; }
private:
    XvueState   state_;          // D-04: owns the single QColor field
    XvueCanvas* canvas_ = nullptr; // central widget; Qt owns lifetime
};
```

```cpp
// xvue/qt/src/xvue_qt_window.cpp
#include "xvue_qt_window.h"
#include "xvue_qt_canvas.h"

XvueWindow::XvueWindow(QWidget* parent)
    : QMainWindow(parent)
{
    setWindowTitle("MEFISTO");                    // D-02
    canvas_ = new XvueCanvas(&state_, this);
    setCentralWidget(canvas_);                    // D-02
    resize(800, 600);                             // D-02 (logical pixels)
}

XvueWindow::~XvueWindow() = default;
```

### 3. `XvueCanvas` (central `QWidget` with minimal `paintEvent`)

```cpp
// xvue/qt/src/xvue_qt_canvas.h
// Source: CONTEXT.md D-05, D-15
#pragma once
#include <QWidget>
#include "xvue_qt_state.h"

class XvueCanvas : public QWidget {
    Q_OBJECT
public:
    explicit XvueCanvas(XvueState* state, QWidget* parent = nullptr);
protected:
    void paintEvent(QPaintEvent* ev) override;
private:
    XvueState* state_;  // non-owning; owned by XvueWindow (D-15)
};
```

```cpp
// xvue/qt/src/xvue_qt_canvas.cpp
#include "xvue_qt_canvas.h"
#include <QPainter>

XvueCanvas::XvueCanvas(XvueState* state, QWidget* parent)
    : QWidget(parent), state_(state) {}

void XvueCanvas::paintEvent(QPaintEvent* /*ev*/) {
    QPainter(this).fillRect(rect(), state_->background_);  // D-05: one line
}
```

### 4. `XvueState` (one field)

```cpp
// xvue/qt/src/xvue_qt_state.h
// Source: CONTEXT.md D-04
#pragma once
#include <QColor>
#include <Qt>

struct XvueState {
    QColor background_ = Qt::black;  // D-04: matches xvuelc.c:935 BlackPixel default
};
```

### 5. `xvfond_` with 2-entry hardcoded table

```cpp
// xvue/qt/src/xvue_qt_api.cpp — xvfond_ replacement body
// Source: CONTEXT.md D-14, D-15
#include "xvue_qt_app.h"
#include "xvue_qt_window.h"
#include "xvue_qt_canvas.h"

extern "C" void proc(xvfond)(int *icolor) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win) return;               // defensive: no window yet, nothing to set

    QColor c;
    switch (*icolor) {
        case 0: c = Qt::black; break;   // D-14: matches X11 BlackPixel
        case 1: c = Qt::white; break;   // D-14: matches X11 WhitePixel
        default: {
            static bool warned = false;
            if (!warned) {
                std::fprintf(stderr,
                    "xvue-qt: xvfond_ palette index %d out of Phase 1 range "
                    "(Phase 3 will add full colormap)\n", *icolor);
                warned = true;
            }
            c = Qt::black;
        }
    }
    win->state()->background_ = c;
    win->canvas()->update();         // D-15: schedule repaint
}
```

### 6. `verify_no_exec` CMake target (POST_BUILD)

```cmake
# Added to xvue/qt/CMakeLists.txt
# Source: CONTEXT.md D-10, mirroring Phase 0 D-12 verify_abi
add_custom_command(TARGET xvueqt POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo
        "verify_no_exec: scanning xvue/qt/src + include for QApplication::exec..."
    COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/cmake/verify_no_exec.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    VERBATIM
    COMMENT "verify_no_exec: enforcing SHELL-03 ban on QApplication::exec()"
)
```

```sh
#!/bin/sh
# xvue/qt/cmake/verify_no_exec.sh
# Source: CONTEXT.md D-10, mirroring xvue/qt/cmake/verify_abi.sh
set -eu
SRC_DIR="$1"
INC_DIR="$2"
# Match either QApplication::exec or qApp->exec() — fail on any hit
if grep -rn 'QApplication::exec\|qApp->exec()' "$SRC_DIR" "$INC_DIR" ; then
    echo "ERROR: SHELL-03 violation — QApplication::exec found in xvue/qt/" >&2
    exit 1
fi
echo "verify_no_exec: clean"
exit 0
```

### 7. `XVUE_QT_ASSERT_MAIN_THREAD` macro body + retrofit pattern

```cpp
// xvue/qt/include/xvue_qt_api.h — replacement for the Phase 0 skeleton
// Source: CONTEXT.md D-12, D-13
#ifdef QT_DEBUG
#  include <QThread>
#  include <QCoreApplication>
#  define XVUE_QT_ASSERT_MAIN_THREAD()                                   \
      do {                                                                \
          if (QCoreApplication::instance()) {                             \
              Q_ASSERT(QThread::currentThread() ==                        \
                       QCoreApplication::instance()->thread());           \
          }                                                                \
      } while (0)
#else
#  define XVUE_QT_ASSERT_MAIN_THREAD() ((void)0)
#endif
```

Retrofit pattern applied to every stub body in `xvue_qt_api.cpp`:

```cpp
// Before (Phase 0):
void proc(effacer)(void) {
    static bool warned = false;
    warn_once(warned, "effacer_");
}

// After (Phase 1 — SHELL-07 retrofit):
void proc(effacer)(void) {
    XVUE_QT_ASSERT_MAIN_THREAD();          // NEW first line
    static bool warned = false;
    warn_once(warned, "effacer_");
}
```

### 8. `prpr/xvtest0.f` test driver

```fortran
      PROGRAM XVTEST0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TESTER L'OUVERTURE/FERMETURE REPETEE DE LA FENETRE XVUE (Qt)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PHASE 1 — SHELL RESEARCH                      AVRIL 2026
C2345X7..............................................................012
C
C     PREMIER CYCLE : OUVERTURE PUIS FERMETURE
      PRINT *,'xvtest0: first XVINITGRAPHIQUE'
      CALL XVINITGRAPHIQUE
      PRINT *,'xvtest0: first XVFERMER'
      CALL XVFERMER
C
C     SECOND CYCLE : LE POINT DE TEST CRITIQUE DE SHELL-02
      PRINT *,'xvtest0: second XVINITGRAPHIQUE (reopen)'
      CALL XVINITGRAPHIQUE
      PRINT *,'xvtest0: second XVFERMER'
      CALL XVFERMER
C
      PRINT *,'xvtest0: OK'
      STOP
      END
```

Fixed-form Fortran, column 7+, `C` comment marker. Matches `prpr/xvtest1.f` style per CLAUDE.md §"Language and module conventions" and CONTEXT.md D-19. Extra `PRINT *` lines are permitted by discretion.

### 9. `bin/cbxvtest0_qt` (cloned from `bin/cbmail_qt`)

```bash
#!/bin/bash
# cbxvtest0_qt — Phase 1 Qt window-shell test driver build script
# Cloned from bin/cbmail_qt per CONTEXT.md D-20.

cd $MEFISTO
if !(test -d pp)
then
  mkdir pp
fi

if test -f $MEFISTO/td/m/anglais
then
  LANGAGE=1
  echo "MEFISTO-LINUX 64bits: Compilation + Link of xvtest0.f (Qt shell test)"
else
  LANGAGE=0
  echo "MEFISTO-LINUX 64bits: Compilation + Edition de liens de xvtest0.f (Qt shell test)"
fi

nompp=pp/ppxvtest0_qt
if test -f $nompp
then
  rm $nompp
fi

QT_LIBS=$(pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport)
echo "gfortran ... $MEFISTO/prpr/xvtest0.f -o $nompp (Qt shell test)"
gfortran -Wall -mcmodel=large -m64 -O -fopenmp -I. prpr/xvtest0.f \
     xvue/lib \
     -Lxvue/qt/build -lxvueqt $QT_LIBS -lstdc++ \
     -lgfortran -o $nompp

if test -f $nompp
then
  chmod 755 $nompp
  echo "$MEFISTO/$nompp is CORRECT (Qt shell test)"
else
  echo "$MEFISTO/$nompp is INCORRECT"
fi
```

**Note:** `xvtest0.f` calls **only** the 7 SHELL-* entry points (via `XVINITGRAPHIQUE`/`XVFERMER`), so the minimal link set may drop `mail/lib` and `util/lib` — planner's call. The reference `cbmail_qt` uses all three; dropping is safe for `xvtest0` since `xvtest0.f` has no mesher or solver code. If unsure, copy the full line for parity with `cbmail_qt`.

## Runtime State Inventory

*(Phase 1 is greenfield new code — no rename/refactor operation. Nothing pre-existing needs re-registration beyond the build system itself.)*

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | None — no persisted state | None |
| Live service config | None — no running services | None |
| OS-registered state | None — no OS hooks, no systemd units, no Windows tasks | None |
| Secrets/env vars | `MEFISTO`, `MEFISTOX`, `DISPLAY`/`WAYLAND_DISPLAY`, optional `QT_SCALE_FACTOR` — all read-only, no rename | None — Phase 1 reads these; does not redefine them |
| Build artifacts | `xvue/qt/build/libxvueqt.a` (Phase 0 artifact) will grow in size after Phase 1 adds new `.o` files; `CMakeFiles/` cache may need a `cmake --build . --clean-first` on first Phase 1 build if `AUTOMOC` doesn't pick up new headers | Clean rebuild of `xvue/qt/build` recommended after the new `.cpp` files are added to `add_library` |

**Nothing found in Stored data / Live service config / OS-registered state.** Phase 1 introduces new source files, new CMake sources, a new Fortran driver, a new shell script, and new binary `pp/ppxvtest0_qt`. There is no pre-existing state to migrate.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| gfortran | Fortran compilation of `prpr/xvtest0.f` | ✓ | 15.2.0 (Debian 15.2.0-11) [VERIFIED] | — |
| gcc/g++ | C++17 compilation of new Qt sources | ✓ | 15.2.0 [VERIFIED via gfortran suite] | — |
| Qt6Core | `QCoreApplication`, `QThread`, `QEventLoop` | ✓ | 6.10.2 [VERIFIED `pkg-config --modversion`] | — |
| Qt6Gui | `QGuiApplication`, `QScreen`, `QPainter`, `QColor` | ✓ | 6.10.2 [VERIFIED via `qt6-base-dev`] | — |
| Qt6Widgets | `QApplication`, `QMainWindow`, `QWidget` | ✓ | 6.10.2 [VERIFIED via `qt6-base-dev`] | — |
| Qt6PrintSupport | Unused in Phase 1, already in link line | ✓ | 6.10.2 | — |
| `qt6-base-dev-tools` (moc) | AUTOMOC generation for `Q_OBJECT` classes | ✓ | 6.10.2 [VERIFIED `dpkg -l`] | — |
| CMake ≥ 3.21 | `xvue/qt/CMakeLists.txt` uses 3.21+ features | ✓ | Phase 0 already builds — implicit verification | — |
| `pkg-config` | Used by `bin/cb*_qt` clone scripts | ✓ | Phase 0 scripts already depend on it | — |
| X11 or Wayland display server | SHELL-06 visual validation on running `pp/ppxvtest0_qt` | ✓ (assumed desktop session) | — | **No fallback** — `QT_QPA_PLATFORM=offscreen` bypasses the rendering path and invalidates SHELL-06 (see Pitfall 11) |
| `libX11-dev` | Legacy `xvue/xvuelc.c` build (Phase 1 does NOT change this) | ✓ (verified Phase 0 baseline builds) | — | — |
| `ImageMagick` (`convert`) | Legacy animated GIF pipeline — NOT used in Phase 1 | N/A | — | — |

**Missing dependencies with no fallback:** None.

**Missing dependencies with fallback:** None.

**Note on display:** Phase 1 validation (D-21) is **human observation** of a blank window appearing twice. This requires an interactive desktop session. The `offscreen` QPA platform plugin would make `pp/ppxvtest0_qt` run non-crashingly, but would not validate SHELL-06 HiDPI rendering. Run validation on a real display.

## Validation Architecture

### Test Framework

| Property | Value |
|----------|-------|
| Framework | MEFISTO custom: shell-script compile (`bin/cb*_qt`) + human visual observation of `pp/pp*_qt`. No automated unit-test framework exists in this codebase. |
| Config file | None (shell scripts drive everything) |
| Quick run command | `bin/cbxvtest0_qt && pp/ppxvtest0_qt` |
| Full suite command | `bin/cbl_tout` (legacy X11, must still succeed) + `bin/cbl_tout_qt` (Qt parallel build, must still succeed) + `bin/cbxvtest0_qt && pp/ppxvtest0_qt` |

### Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| SHELL-01 | `xvinitgraphique_` opens `QMainWindow` with `XvueCanvas` | manual-visual | `bin/cbxvtest0_qt && pp/ppxvtest0_qt` + human observes first window | ❌ Wave 0 (`prpr/xvtest0.f`, `bin/cbxvtest0_qt`) |
| SHELL-02 | Reopen without "there can only be one" | manual-visual + runtime | `pp/ppxvtest0_qt` exits 0; two windows observed | ❌ Wave 0 |
| SHELL-03 | No `QApplication::exec` in `xvue/qt/` | automated build-time | `cd xvue/qt/build && cmake --build . --target xvueqt` — POST_BUILD `verify_no_exec` step | ❌ Wave 0 (`verify_no_exec.sh` + CMake rule) |
| SHELL-04 | `xvpxecran_`/`xvmmecran_` from `QScreen` logical px / mm | manual + optional runtime test | Extend `xvtest0.f` with `CALL XVPXECRAN(X,Y); PRINT *,X,Y` and `CALL XVMMECRAN(XMM,YMM); PRINT *,XMM,YMM` — compare against legacy `pp/ppmail` values | ❌ Wave 0 (optional extension) |
| SHELL-05 | `xvfond_` updates background without corrupting pixmap | manual-visual | Extend `xvtest0.f` with `CALL XVFOND(1)` between init cycles → second window observed white | ❌ Wave 0 (optional extension) |
| SHELL-06 | DPR > 1 HiDPI correctness | manual-visual | `QT_SCALE_FACTOR=2 pp/ppxvtest0_qt` — window still reads 800×600 logical | ❌ Wave 0 |
| SHELL-07 | Every `extern "C"` entry has debug-mode thread assertion | automated build-time | `grep -c 'XVUE_QT_ASSERT_MAIN_THREAD' xvue/qt/src/xvue_qt_api.cpp` should equal total entry-point count (~57 including the 7 SHELL). Optional CMake custom target `verify_thread_assert` mirroring `verify_abi`. | ❌ Wave 0 (optional — grep check can be a post-build shell command) |

### Sampling Rate

- **Per task commit:** `cd xvue/qt/build && cmake --build .` — fast (seconds), catches compile/link errors and runs `verify_abi` + `verify_no_exec` build-time guards
- **Per wave merge:** `bin/cbxvtest0_qt && pp/ppxvtest0_qt` — visual confirmation of reopen cycle
- **Phase gate:** (1) `bin/cbl_tout` (legacy X11) builds clean; (2) `bin/cbl_tout_qt` (Qt 5 canonical modules) builds clean; (3) `bin/cbxvtest0_qt && pp/ppxvtest0_qt` shows blank 800×600 "MEFISTO" window twice under `QT_SCALE_FACTOR=2`; (4) `verify_abi` and `verify_no_exec` CMake guards both green.

### Wave 0 Gaps

- [ ] `prpr/xvtest0.f` — new Fortran test driver (D-19)
- [ ] `bin/cbxvtest0_qt` — new build script cloned from `bin/cbmail_qt` (D-20)
- [ ] `xvue/qt/cmake/verify_no_exec.sh` — new helper script (D-10 discretion)
- [ ] CMake `add_custom_command(TARGET xvueqt POST_BUILD ...)` for `verify_no_exec` added to `xvue/qt/CMakeLists.txt` (D-10)
- [ ] New C++ sources added to `add_library(xvueqt STATIC ...)` list in `xvue/qt/CMakeLists.txt`:
  - [ ] `src/xvue_qt_app.cpp`
  - [ ] `src/xvue_qt_window.cpp`
  - [ ] `src/xvue_qt_canvas.cpp`
  - [ ] `src/xvue_qt_state.cpp` (may be header-only struct; if so, no .cpp needed)
- [ ] No existing test framework installation needed — MEFISTO uses shell scripts + visual observation only. SHELL-03 and SHELL-07 automation is via CMake custom commands (build-time), not runtime test assertions.

## Security Domain

*(Not applicable — Phase 1 is local desktop scientific tool, no network-facing surface, no authentication, no user-supplied input beyond Fortran numerical arguments. `security_enforcement` is not set in `.planning/config.json`, defaulting to enabled — but ASVS categories V2/V3/V4/V5/V6 do not apply to a Qt window-shell migration for a single-user scientific tool that reads no untrusted input during this phase.)*

Minimal applicable concerns:

| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | no | N/A (single-user desktop) |
| V3 Session Management | no | N/A |
| V4 Access Control | no | N/A (local files only, filesystem-gated) |
| V5 Input Validation | no | Phase 1 reads no untrusted data — `xvfond_` gets an `int*` from trusted Fortran code |
| V6 Cryptography | no | N/A |

**Threat model for Phase 1:** None beyond "don't break the build" and "don't crash the process." The CLAUDE.md working rules (compilation must never break, tests must continue to pass) are the de facto security controls for this phase.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `Qt::AA_EnableHighDpiScaling` attribute set before `QApplication` ctor | Default-on in Qt 6; attribute is a no-op | Qt 6.0 (2020) | Phase 1 does NOT need to set the attribute — SHELL-06 is satisfied by defaults. Deferred §"Defensive `Qt::AA_EnableHighDpiScaling` attribute call." |
| Qt 5 `QDesktopWidget::screenGeometry()` | `QGuiApplication::primaryScreen()->geometry()` / `->size()` | Qt 5.14 deprecation, removed in Qt 6 | Phase 1 uses `primaryScreen()->size()` per D-16. |
| `QApplication::exec()` at top of `main()` | Embedded `processEvents(ExcludeUserInputEvents)` per entry point, with nested `QEventLoop` for blocking reads (Phase 5) | MEFISTO-specific embedding pattern | Phase 1 uses the one-shot pump per D-01. |
| `CMAKE_AUTOMOC` set after `find_package(Qt...)` | AUTOMOC must be set **before** `find_package` | Phase 0 pitfall, Qt CMake docs | Phase 0 already enforces this order (`xvue/qt/CMakeLists.txt:11`). |

**Deprecated/outdated:**
- `Qt::AA_EnableHighDpiScaling` — Qt 6 no-op [CITED: doc.qt.io/qt-6/highdpi.html]
- `QDesktopWidget` — removed in Qt 6 [CITED]
- `QApplication::setAttribute(Qt::AA_UseHighDpiPixmaps)` — Qt 6 no-op [CITED]

## Project Constraints (from CLAUDE.md)

| Directive | Source | Phase 1 Compliance |
|-----------|--------|---------------------|
| Compilation must never break (`bin/cbl_tout` always succeeds) | CLAUDE.md §"Working rules" | Phase 1 touches nothing outside `xvue/qt/`, `prpr/xvtest0.f`, `bin/cbxvtest0_qt`. Run `bin/cbl_tout` after every commit. |
| Small `testa/`/`testf/` tests must keep passing | CLAUDE.md §"Tests" | Phase 1 does not run A/B `testa/` cases (D-21 defers to Phase 2), but legacy X11 `bin/cbl_tout` + `testa/` must still succeed — run at phase gate. |
| Ask before installing packages | CLAUDE.md §"Asking before acting" | Phase 1 needs no new packages (Qt 6.10.2 already installed). |
| Git discipline: commit after every logical step, no force-push, no hook bypass | CLAUDE.md §"Git discipline" | Each task in the plan should be a separate commit; commit messages describe what and why. |
| Fixed-form Fortran column 7+ for `prpr/xvtest0.f` | CLAUDE.md §"Language and module conventions" | Phase 1 follows `prpr/xvtest1.f` template literally. |
| `doc/normes.ps` coding norms | CLAUDE.md §"Programming norms" | Phase 1 is all-new C++; no existing Fortran routines edited. For `prpr/xvtest0.f`, follow `prpr/xvtest1.f` naming/comment style. |
| Qt migration keeps graphics calls isolated in `xvue/` | CLAUDE.md §"Active project goals" | Phase 1 keeps all new Qt code under `xvue/qt/`. No leakage elsewhere. |
| Never Ctrl-C an interactive module; use `99;` | CLAUDE.md §"Running a project" | N/A — `xvtest0.f` is a non-interactive test driver that exits via `STOP END`. |

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | `QScreen::physicalSize()` returns **millimetres** regardless of DPR on all Qt 6 platforms | Code Examples §4, Pattern 4 | LOW — Qt 6 docs explicitly state this; verified stable API since Qt 5.0. If wrong, `xvmmecran_` returns wrong units — visible immediately by comparing against legacy X11 `pp/ppmail` output. Mitigation: add `PRINT *,XMM,YMM` lines to `xvtest0.f` and compare against legacy. [ASSUMED → can be upgraded to VERIFIED by a small runtime test] |
| A2 | Qt 6.10.2 enables HiDPI scaling by default and `QGuiApplication::primaryScreen()->size()` returns logical (not physical) pixels | State of the Art, Pitfall 2 | LOW — this is the entire Qt 6 HiDPI headline change, extensively documented. [ASSUMED but trained knowledge is strong for Qt 6.0+ defaults] |
| A3 | `std::atexit` handler running `qApp->quit()` + `unique_ptr::reset()` is Qt-safe at process exit | Pattern 1, Pitfall 5 | MEDIUM — ARCHITECTURE.md claims this is the "documented Qt embedding idiom" with HIGH confidence, but the exact ordering of static destructors vs. `atexit` handlers is C++-runtime dependent. If wrong, exit-time crash (visible but harmless since `xvtest0.f` already did its validation). Mitigation: `prpr/xvtest0.f` exits via `STOP END` — any post-exit crash is a visible failure that can be fixed incrementally. [ASSUMED — consensus pattern, not empirically verified in this session] |
| A4 | `QApplication` can be constructed after `std::atexit` has been called once in the process lifetime, provided no previous `QApplication` was destroyed — i.e. "construct once, never destroy until atexit" holds | Pattern 1, Pitfall 5 | LOW — this is what `std::call_once` enforces. [ASSUMED but structurally enforced] |
| A5 | Adding new `.cpp` files to `add_library(xvueqt STATIC ...)` automatically triggers AUTOMOC on their paired `.h` files, provided the `.h` is `#include`d by the `.cpp` | Pitfall 7, Wave 0 Gaps | LOW — this is standard Qt 6 AUTOMOC behavior. If wrong, a `moc_` generation step must be added explicitly. [ASSUMED; verified by Qt CMake docs] |
| A6 | `Q_OBJECT` macro works in `XvueWindow`/`XvueCanvas` headers even when Phase 1 declares no signals/slots | Pitfall 8 | LOW — `Q_OBJECT` only requires it to be inside a class body; no signals/slots content needed. [VERIFIED: standard Qt usage] |
| A7 | The `verify_no_exec` grep guard can safely match `QApplication::exec` and `qApp->exec()` strings literally without risk of false positives inside comments | Code Example §6 | LOW — if a future Phase 1 comment says "we do not call `QApplication::exec()`", the guard fails. Mitigation: either (a) keep such explanatory text out of `xvue/qt/src/` / `xvue/qt/include/` (put it in planning docs only), or (b) refine the grep to exclude comments. Flag this to the planner. [ASSUMED; low impact] |

**Upgrade path for assumptions:** A1 and A2 can both be promoted to VERIFIED via a small runtime test — extend `prpr/xvtest0.f` (discretion per D-19) to print the values of `XVPXECRAN` and `XVMMECRAN`, compare against legacy `pp/ppmail` output. This is 4 lines of Fortran and is the recommended empirical validation for SHELL-04/SHELL-06.

## Open Questions (RESOLVED)

1. **Should the `verify_no_exec` grep be case-sensitive and exclude comments?**
   - What we know: D-10 mandates a grep that matches `QApplication::exec` and `qApp->exec()`.
   - What's unclear: Whether a developer comment in a future phase's code (e.g., "// NOT calling QApplication::exec here on purpose") would trigger a false positive.
   - Recommendation: Use `grep -rn 'QApplication::exec\|qApp->exec()'` with no exclusions for simplicity in Phase 1; refine to exclude comment lines (`grep -v '^\s*//'`) only if a false positive appears. Document the behavior in a comment at the top of `verify_no_exec.sh`.
   - **RESOLVED:** Use the simple grep (`grep -rn 'QApplication::exec\|qApp->exec()'`) with no comment exclusion in Phase 1. Rationale: Phase 1 source files contain no comments referencing `exec()`; false-positive risk is zero in the current scope. If a future phase introduces such a comment, the grep will fail loudly and Plan 01-01 Task 3 will be revisited to add `grep -v '^\s*//'`. Decision is documented in a header comment at the top of `xvue/qt/cmake/verify_no_exec.sh` per Plan 01-01 Task 3 action.

2. **Should `XvueState` be header-only or header+source?**
   - What we know: D-04 says Phase 1 has exactly one field `QColor background_ = Qt::black;`.
   - What's unclear: Whether a `.cpp` is warranted for such a thin struct.
   - Recommendation: Header-only. Struct inline in `xvue_qt_state.h`, no `xvue_qt_state.cpp`. Phase 2 will add `.cpp` when pen/brush/painter members need out-of-line definitions.
   - **RESOLVED:** Header-only. `xvue/qt/src/xvue_qt_state.h` holds the full struct inline with `QColor background_ = Qt::black;` as the sole field. No `.cpp` in Phase 1. Rationale: single trivial field, no out-of-line definitions needed, AUTOMOC is not triggered (no `Q_OBJECT`). Plan 01-01 Task 2 materialises this decision.

3. **Should `cbxvtest0_qt` link `mail/lib` + `util/lib` for parity with `cbmail_qt`, or drop them since `xvtest0.f` only calls shell entry points?**
   - What we know: D-20 says "clone from the thinnest `cb*_qt`"; `cbmail_qt` links three Fortran libs.
   - What's unclear: Dropping libs saves link time but risks "undefined reference" if any `xvue/lib` wrapper (`xvouvrir.f`, etc.) depends on a symbol in the dropped libs.
   - Recommendation: Copy the full 3-lib link line from `cbmail_qt` verbatim. Link-time cost is negligible; symmetry is worth more than a small build-time saving.
   - **RESOLVED (reversing the original recommendation):** Drop `mail/lib` from `bin/cbxvtest0_qt`; keep `util/lib` and `xvue/lib`. Rationale: (a) `prpr/xvtest0.f` calls only `XVINITGRAPHIQUE` and `XVFERMER`, whose Fortran wrappers live entirely in `xvue/lib` and transitively use only `util/lib` helpers — none of the `mail/lib` mesh-generation symbols are referenced, so the link is sound; (b) D-20's "clone-and-modify" directive is about preserving the build-script *template*, not mandating a literal copy of irrelevant libs; (c) the minimal link surface makes link failures self-diagnosing (any undefined reference points immediately at a real `xvue/lib`→`mail/lib` coupling we'd want to know about). Plan 01-03 Task 2 explicitly drops `mail/lib` and documents a fallback: if the link fails on a symbol resolution, restore `mail/lib` unchanged. This deliberately contradicts the original recommendation above; the research author's original concern (symmetry) is subordinated to the plan's goal of "smallest possible compilation unit" per D-21.

4. **Should `XVUE_QT_ASSERT_MAIN_THREAD()` null-guard `QCoreApplication::instance()` via the `if (QCoreApplication::instance())` form (code example §7) or via a direct `if (qApp)` — functionally equivalent but `qApp` pulls in `<QCoreApplication>` transitively?**
   - What we know: D-13 says the macro must be null-safe when `qApp` is null.
   - What's unclear: Stylistic question, no behavior difference.
   - Recommendation: Use `QCoreApplication::instance()` directly — `qApp` is a macro that expands to `QCoreApplication::instance()` cast to `QApplication*`, which requires the `<QApplication>` header in the expansion site. `QCoreApplication::instance()` is lower-coupling for a header that is `#include`d by every `.cpp` in `xvue/qt/src/`.
   - **RESOLVED:** Use `QCoreApplication::instance()` directly inside `XVUE_QT_ASSERT_MAIN_THREAD()`. Rationale: lower header coupling for `xvue_qt_api.h`, which is included by every `.cpp` in `xvue/qt/src/`. Plan 01-01 Task 1 implements the macro body with this form.

## Sources

### Primary (HIGH confidence)

- `/home/drico/git/mefisto/.planning/phases/01-window-shell-xvueapp-xvuewindow-xvuecanvas/01-CONTEXT.md` — 21 locked decisions (D-01 through D-21)
- `/home/drico/git/mefisto/.planning/REQUIREMENTS.md` §"Shell — QApplication and window lifecycle" — SHELL-01 through SHELL-07
- `/home/drico/git/mefisto/.planning/ROADMAP.md` §"Phase 1" — phase boundary, depends-on, success criteria
- `/home/drico/git/mefisto/.planning/research/ARCHITECTURE.md` — `XvueApp`/`XvueWindow`/`XvueCanvas`/`XvueState` split, singleton discipline (`std::call_once` + fabricated argc/argv + `atexit`), event-loop-without-`exec` strategy (Option B)
- `/home/drico/git/mefisto/.planning/research/PITFALLS.md` §Pitfalls 5, 6, 7, 8, 9, 10, 11 — double-init, `exec()` ban, `processEvents` discipline, AUTOMOC, PIC, OpenMP collision
- `/home/drico/git/mefisto/xvue/xvuelc.c` lines 286-356, 612-980, 935, 1434, 1598, 1602-1617 — legacy semantics to preserve
- `/home/drico/git/mefisto/xvue/qt/include/xvue_qt_api.h` — Phase 0 ABI header (57 entries)
- `/home/drico/git/mefisto/xvue/qt/src/xvue_qt_api.cpp` — Phase 0 warn-once stubs
- `/home/drico/git/mefisto/xvue/qt/CMakeLists.txt` — Phase 0 CMake file (grows in Phase 1)
- `/home/drico/git/mefisto/xvue/qt/cmake/verify_abi.sh` — Phase 0 pattern for `verify_no_exec`
- `/home/drico/git/mefisto/prpr/xvtest1.f` — template for new `prpr/xvtest0.f`
- `/home/drico/git/mefisto/bin/cbmail_qt` — template for new `bin/cbxvtest0_qt`
- `/home/drico/git/mefisto/CLAUDE.md` — working rules (build must not break, Fortran norms, git discipline)
- `pkg-config --modversion Qt6Core` → 6.10.2 [VERIFIED]
- `dpkg -l | grep qt6-base-dev` → 6.10.2+dfsg-6 [VERIFIED]
- `gfortran --version` → GNU Fortran (Debian 15.2.0-11) 15.2.0 [VERIFIED]

### Secondary (MEDIUM confidence)

- Qt 6 `QCoreApplication::processEvents` documentation (doc.qt.io/qt-6/qcoreapplication.html#processEvents) — `QEventLoop::ExcludeUserInputEvents` semantics [CITED via training + canonical_refs in CONTEXT.md]
- Qt 6 `QScreen` / `QGuiApplication::primaryScreen` documentation (doc.qt.io/qt-6/qscreen.html) — `size()` logical px, `physicalSize()` mm [CITED]
- Qt 6 `QMainWindow` / `QWidget` documentation — `setCentralWidget`, `show()`, `close()`, `resize()` [CITED]
- Qt 6 HiDPI documentation (doc.qt.io/qt-6/highdpi.html) — default-on scaling, removed attributes [CITED]
- Qt 6 CMake AUTOMOC documentation — `CMAKE_AUTOMOC` before `find_package`, scans `add_library` sources [CITED]

### Tertiary (LOW confidence)

- None — every factual claim in this research is backed either by (a) locked decisions in CONTEXT.md, (b) Phase 0 artifacts inspected in this session, (c) legacy `xvuelc.c` code read directly, or (d) stable, well-documented Qt 6 idioms. The `[ASSUMED]` tags in §Assumptions Log flag the places where empirical validation (A1, A2) or runtime observation (A3) would upgrade confidence.

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH — Qt 6.10.2 installed and verified via `pkg-config` and `dpkg -l` in this session; no version mismatch risk.
- Architecture: HIGH — all 21 decisions in CONTEXT.md are locked; research confirms they map cleanly to Qt 6 canonical idioms documented in PITFALLS.md and ARCHITECTURE.md.
- Pitfalls: HIGH — every Pitfall relevant to Phase 1 (1 through 11 from PITFALLS.md, plus Phase 1-specific 5 through 11 in this document) is either structurally prevented by a decision (Pitfalls 1, 4, 6, 9) or mitigated by a documented rule (Pitfalls 2, 3, 5, 7, 8, 10, 11).
- Test framework: HIGH — MEFISTO has no test framework; validation is shell-script + human visual observation, fully documented in CLAUDE.md and Phase 0 conventions. The Wave 0 gaps list 3 new files (`prpr/xvtest0.f`, `bin/cbxvtest0_qt`, `xvue/qt/cmake/verify_no_exec.sh`) plus 1 CMake edit — modest, well-scoped.
- Environment: HIGH — all dependencies verified available in this session.

**Research date:** 2026-04-11
**Valid until:** 2026-05-11 (30 days — Qt 6.10 is a stable LTS-aligned release; no expected breakage. Revisit if Qt 6.11 lands and alters AUTOMOC or HiDPI defaults.)
