# Phase 1: Window shell (`XvueApp` + `XvueWindow` + `XvueCanvas`) - Context

**Gathered:** 2026-04-11
**Status:** Ready for planning

<domain>
## Phase Boundary

A blank Qt `QMainWindow` opens through `xvinitgraphique_` and closes cleanly through `xvfermer_`, proving the `QApplication` singleton discipline, lazy lifetime, reopen semantics, and HiDPI convention **in isolation before any drawing logic exists**. Phase 1 delivers SHELL-01 through SHELL-07: the `XvueApp` singleton, the `XvueWindow`/`XvueCanvas` split, screen-metric accessors, the `xvfond_` background entry point, the `QApplication::exec()` ban with build-time enforcement, and the universal debug-build thread-affinity assertion.

**In scope:** `XvueApp::ensure()` (call_once + atexit), lazy `XvueWindow` construction, `XvueCanvas` widget with a minimal `paintEvent` that fills the widget background, a minimal `XvueState` holding only the background color, `xvinitgraphique_`/`xvfermer_`/`xvpxecran_`/`xvmmecran_`/`xvfond_` as real implementations, `xvinfo_` as a partial implementation that only resizes the window (palette/font outputs remain warn-once stubs until Phase 3), `SHELL-07` assertion macro body and retrofit into every existing `xvue_qt_api.cpp` stub, the CMake build-time `verify_no_exec` guard, and a minimal Fortran test driver `prpr/xvtest0.f` → `pp/ppxvtest0_qt` that exercises the reopen cycle.

**Explicitly out of scope:** drawing primitives (Phase 2), backing `QPixmap` (Phase 2), text/fonts (Phase 3), colormap and palette plumbing (Phase 3), pixmap save/restore (Phase 4), event loop and mouse/keyboard delivery (Phase 5), menus/toolbar/dock widgets (Phase 6), image/GIF/PS export (Phase 7). None of these land in Phase 1 even partially; every entry point not listed above stays on the Phase 0 warn-once stub path.

</domain>

<decisions>
## Implementation Decisions

### Window-open timing and sizing
- **D-01:** `xvinitgraphique_` is the single entry point that (a) ensures `XvueApp`'s `QApplication` via `std::call_once`, (b) allocates a fresh `XvueWindow` if `window_` is null, (c) calls `window_->show()`, and (d) calls `QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents)` exactly once to realize the window on the compositor before returning. The `ExcludeUserInputEvents` flag is mandatory — a bare `processEvents()` would re-enter graphics entry points and violate the "no input processing outside nested `QEventLoop`" invariant set up for Phase 5 (Pitfall 6, Pitfall 8).

  *(original, superseded — see revised text below)*

  **D-01 revised (after empirical verification in 01-03 — see .planning/debug/resolved/phase-01-xvtest0-teardown-segfault.md):** `xvinitgraphique_` pumps `processEvents(ExcludeUserInputEvents, 20)` in a **bounded loop** until `QWindow::isExposed()` is true or 2 seconds elapse. A single `processEvents` call is insufficient to realize an X11 window: the compositor requires multiple event-loop trips to dispatch MapRequest → ConfigureNotify → Expose before `QWindow::isExposed()` returns true. The `ExcludeUserInputEvents` flag and the Pitfall 6/8 rationale are unchanged. The implementation uses `QElapsedTimer` for the 2-second timeout and `win->windowHandle()->isExposed()` as the break condition. See fix FIX-2 in the debug session for the exact code pattern.
- **D-02:** Initial window size is **800×600 logical pixels** with title `"MEFISTO"`. This is a deliberate sentinel size small enough to fit any monitor and distinct enough from common defaults to be visibly "Phase 1" until `xvinfo_` resizes it. `XvueWindow` constructor sets the title and resizes the central `XvueCanvas` to match.
- **D-03:** `xvinfo_` is a **partial real implementation** in Phase 1: it calls `window_->resize(*ix, *iy)` (logical pixels) and leaves every other output parameter (palette indices, font tables, visual class) on the Phase 0 warn-once path — the function body writes zeroed/default values into the output pointers, prints the one-shot `"xvue-qt: stub xvinfo_ palette outputs not implemented yet"` line, and returns. Phase 3 replaces the warn-once block with the real colormap plumbing. Rationale: Fortran callers still drive window sizing through `xvinfo_` as the legacy X11 code does, so this hook cannot stay a pure stub without breaking the mesher's startup path; the non-sizing output parameters are not consumed by any Phase 1 code path.
- **D-04:** The `XvueState` struct introduced in Phase 1 holds exactly one field: `QColor background_ = Qt::black;` (matching the legacy X11 default `attributes.background_pixel = BlackPixel(display_mef, screen_mef)` at `xvue/xvuelc.c:935`). Phase 2 extends the same struct with pen, brush, line width, and the long-lived `QPainter*`. No other state is added in Phase 1.
- **D-05:** `XvueCanvas::paintEvent` in Phase 1 contains one line: `QPainter(this).fillRect(rect(), state_->background_)`. No backing `QPixmap` exists yet. Phase 2 replaces this body with `QPainter(this).drawPixmap(0, 0, *backing_)` as a single-line swap.

### `xvfermer_` and reopen semantics
- **D-06:** `xvfermer_` destroys the `XvueWindow` via `window_.reset()` on the `std::unique_ptr<XvueWindow>` owned by `XvueApp`. It does **not** touch the `QApplication`, does **not** call `qApp->quit()`, and does **not** emit any event back into Fortran. The `XvueState` and (in Phase 2+) the backing `QPixmap` die with the window — this is deliberate: each reopen is treated as a fresh display session, matching the real-world mesher → solver → post-processing flow where losing pixmap state between sessions is the expected behavior, not a bug.

  **D-06 addendum (after empirical verification in 01-03 — see .planning/debug/resolved/phase-01-xvtest0-teardown-segfault.md):** After `window_.reset()`, `xvfermer_` must call `QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents)` once to drain the deferred-delete events queued by widget teardown (`deleteLater()` posted by `QObject` internals during `QWidget` destruction). Without this drain, stale DeferredDelete events remain in Qt's event queue and are only flushed when `QApplication` destructs — at the worst possible time (atexit). See fix FIX-3 in the debug session.
- **D-07:** A second call to `xvinitgraphique_` after `xvfermer_` observes `window_ == nullptr` (the `QApplication` `call_once` flag is already set so the app is not re-created) and allocates a fresh `XvueWindow` through the normal lazy path. This satisfies SHELL-02's "no 'QApplication: there can only be one' assertion" requirement and is exactly the pattern the Phase 1 test driver `prpr/xvtest0.f` exercises.
- **D-08:** The `QApplication` is torn down by an `std::atexit`-registered handler installed inside `XvueApp::ensure()` on first call. The handler calls `qApp->quit()` then resets the `std::unique_ptr<QApplication>`. No destructor runs at static-destruction time — this is the documented Qt idiom for embedding Qt inside a non-Qt main (Pitfall 5, ARCHITECTURE.md §"Singleton discipline for `QApplication`").

  *(original, superseded — see revised text below)*

  **D-08 revised (after empirical verification in 01-03 — see .planning/debug/resolved/phase-01-xvtest0-teardown-segfault.md):** `teardown_atexit()` pumps `processEvents(ExcludeUserInputEvents)` around `window_.reset()` to drain any remaining widget deferred-delete events, then calls `qapp_.release()` — explicitly **not** `qapp_.reset()`. `QApplication` is **deliberately leaked** for the entire process lifetime. Destroying `QApplication` inside an `atexit` handler that runs alongside libgfortran's own atexit chain is empirically unsafe on Linux/Qt 6: Qt's internal static teardown interleaves with libgfortran's `__run_exit_handlers` and crashes. The canonical "embed Qt in a non-Qt main" idiom is: construct once, never destroy — let the OS reclaim memory at process exit. The original D-08 design was over-specified: "construct once" is the load-bearing invariant; "destroy at atexit" is not. See fix FIX-1 in the debug session for the exact implementation.
- **D-09:** `std::call_once` guards **only** the `QApplication` construction, not the `XvueWindow`. The window is guarded by a plain `if (!window_) …` check inside `xvinitgraphique_`. This separation is load-bearing for the reopen story: the `QApplication` must be constructed exactly once per process, but `XvueWindow` must be reconstructable 0-to-N times.

### SHELL-03 `exec()` ban enforcement
- **D-10:** A new CMake custom target `verify_no_exec` runs at the end of the `xvue/qt/` build (`add_custom_command(TARGET xvueqt POST_BUILD …)`). It invokes `grep -rn 'QApplication::exec\|qApp->exec()' xvue/qt/src xvue/qt/include` and fails the build with a clear error if any match is found. This follows the **same pattern** established by the Phase 0 `verify_abi` target (Phase 0 D-12) — a post-build grep-and-fail — so the Qt layer's CMake file grows by exactly one symmetric block, not a new enforcement mechanism.
- **D-11:** No git pre-commit hook is added. The CMake guard is the sole enforcement point. Rationale: single-developer project with no CI; every build already runs CMake, so the guard cannot be bypassed without noticing; adding a git hook would require install coordination and can be skipped with `--no-verify`, which is strictly worse than a build-time failure.

### SHELL-07 thread assertion scope
- **D-12:** `XVUE_QT_ASSERT_MAIN_THREAD()` macro body is fleshed out from its Phase 0 skeleton to the real `Q_ASSERT(QThread::currentThread() == qApp->thread())` in debug builds and empty in release, and is **retrofitted into every one of the ~60 existing stubs** in `xvue/qt/src/xvue_qt_api.cpp` as the first statement of each function body, ahead of the warn-once `stderr` print. This is a single-file bulk edit, cheap to apply and review, and satisfies the literal "every `extern "C"` entry point" requirement of SHELL-07 without leaving a trail of TODOs for Phase 2–7 to maintain.
- **D-13:** The `XVUE_QT_ASSERT_MAIN_THREAD()` call is placed **before** any access to `qApp`. In stubs where `qApp` might legitimately be null (only the very first call to `xvinitgraphique_` before `XvueApp::ensure()` has run), the macro is a no-op under `QT_NO_DEBUG` and under debug only asserts if `qApp != nullptr`. The macro definition therefore includes a null-guard: `if (qApp) Q_ASSERT(QThread::currentThread() == qApp->thread())`.

### `xvfond_` storage
- **D-14:** `xvfond_(int *icolor)` in Phase 1 stores the requested background color in `XvueState::background_`. Because Phase 1 has no palette (deferred to Phase 3), the `icolor` integer is mapped via a minimal hardcoded 2-entry table: `0 → Qt::black`, `1 → Qt::white` (matching the `BlackPixel`/`WhitePixel` convention of legacy X11). Any other value falls through to `Qt::black` with a one-shot `stderr` warning `"xvue-qt: xvfond_ palette index N out of Phase 1 range (Phase 3 will add full colormap)"`. The minimal mapping is sufficient for Phase 1 validation and disappears cleanly when Phase 3 wires in the real palette.
- **D-15:** `XvueCanvas` holds a raw pointer to the `XvueState` owned by `XvueWindow`, not a copy. When `xvfond_` updates the field, the canvas automatically observes the new value on its next paint. Phase 1 explicitly calls `canvas_->update()` at the end of `xvfond_` to schedule a repaint; the next `QCoreApplication::processEvents` invocation (either from `xvinitgraphique_` or from a future Phase 2 drawing call) flushes it.

### SHELL-04 screen metrics source
- **D-16:** `xvpxecran_(int *xp, int *yp)` returns `QGuiApplication::primaryScreen()->size().width()` and `->size().height()` — logical pixels from the primary screen, deterministic regardless of whether the window has been created yet. Multi-monitor awareness (switching to `XvueWindow::screen()` after the window is shown) is **explicitly deferred** until a later phase if it becomes needed.
- **D-17:** `xvmmecran_(int *xmm, int *ymm)` returns `QGuiApplication::primaryScreen()->physicalSize().width()` and `->physicalSize().height()` rounded to int (millimetres). Qt 6 guarantees `physicalSize()` is in millimetres regardless of DPR.
- **D-18:** Both entry points are callable **before** `xvinitgraphique_` creates the window, provided `QGuiApplication` exists. Because `XvueApp::ensure()` constructs `QApplication` (which is-a `QGuiApplication`) on the first call to **any** entry point, and because the Phase 0 stubs are retrofitted in D-12 to call `XVUE_QT_ASSERT_MAIN_THREAD()` which dereferences `qApp`, a simple fix applies: every real Phase 1 implementation starts with `XvueApp::ensure();` as its first statement, **before** `XVUE_QT_ASSERT_MAIN_THREAD()`. This guarantees `qApp` is non-null for the assertion and that `QGuiApplication::primaryScreen()` is callable.

### Phase 1 test driver
- **D-19:** A new minimal Fortran driver `prpr/xvtest0.f` is added following the existing `prpr/xvtest{1,2,3,4}.f` convention. Its body is: `CALL XVINITGRAPHIQUE; CALL XVFERMER; CALL XVINITGRAPHIQUE; CALL XVFERMER; STOP END`. This exercises SHELL-01 (first open), SHELL-02 (close-then-reopen without crash), and the `QApplication`-stays-alive invariant, in the smallest possible compilation unit. The existing `prpr/xvtest{1..4}.f` drivers (unused by Phase 0, reserved for Phase 2's drawing primitives per the ROADMAP research-flag) are untouched.
- **D-20:** A new `bin/cbxvtest0_qt` script is added, cloned from an existing `bin/cb*_qt` Phase 0 variant (whichever is the thinnest — likely `bin/cbmail_qt`) with the source file, executable name, and linker line adjusted. It produces `pp/ppxvtest0_qt`. `bin/cbl_tout_qt` is **not** modified to build `xvtest0_qt` — the new script is standalone and invoked directly for Phase 1 validation. Rationale: keeps `bin/cbl_tout_qt` scoped to the five canonical module executables (matching Phase 0 D-08's scope anchor); adding a sixth target would break that symmetry.
- **D-21:** Phase 1 validation gate is: `bin/cbxvtest0_qt` builds cleanly, `pp/ppxvtest0_qt` runs end-to-end without crash, and a human observes the blank 800×600 `"MEFISTO"` `QMainWindow` appear **twice** (once per `XVINITGRAPHIQUE` call). The five baseline `testa/` cases from `.planning/validation/BASELINE.md` are **not** re-run in Phase 1 — they are still on the warn-once Phase 0 path for everything except the five SHELL entry points, and nothing in Phase 1 changes legacy X11 behaviour. The 5-case A/B re-run resumes at Phase 2.

### Claude's Discretion

- **Exact implementation of `XvueApp::ensure()`** — function vs. static-method-on-singleton class, header/source split, member ordering. The only constraints are: `std::call_once` guards the `QApplication` construction, static `argc`/`argv` arrays survive the application, `atexit` registers a teardown on first call, and the function is safe to call from any Phase 1 entry point. Internal organization is up to the planner.
- **`XvueWindow` parent class choice** — `QMainWindow` is required by SHELL-01, but whether Phase 1 sets a menu bar, status bar, or dock frame (all empty) is a planner call. The simplest Phase 1 form is a bare `QMainWindow` with only the central widget set. Phase 6 brings the menu/toolbar surface.
- **`xvtest0.f` exact line count** — must call the two entry points twice in sequence per D-19, but whether additional logging calls (e.g., `PRINT *, 'after xvinitgraphique'`) are added for human diagnostic purposes is discretionary.
- **CMake `verify_no_exec` target exact invocation** — whether the grep is inlined as a CMake command, dispatched through a small shell script, or implemented via `execute_process` is a planner call. The contract is: it runs on every build, fails hard if `QApplication::exec` or `qApp->exec()` appears in `xvue/qt/src` or `xvue/qt/include`, and emits a clear error message naming the offending file and line.
- **Header/source file fanout inside `xvue/qt/`** — ARCHITECTURE.md proposes one file per component (`xvue_qt_app.{h,cpp}`, `xvue_qt_window.{h,cpp}`, `xvue_qt_canvas.{h,cpp}`, `xvue_qt_state.{h,cpp}`). Phase 1 may follow this split literally or collapse two or three of them into a single pair if the code volume is small. The only hard rule: no Phase 1 code leaks outside `xvue/qt/include/` and `xvue/qt/src/` (Phase 0 D-02).
- **Pass `QT_SCALE_FACTOR` from the shell vs. set `Qt::AA_EnableHighDpiScaling` in code** — Qt 6 enables high-DPI scaling by default; SHELL-06 can be validated by setting `QT_SCALE_FACTOR=2` at the shell level when running the test driver on a non-HiDPI display. No code change is required for default HiDPI support. The planner decides whether to add a defensive `QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling)` call ahead of `QApplication` construction; the default behavior already covers the requirement.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing Phase 1.**

### Phase scope anchors
- `.planning/ROADMAP.md` §"Phase 1: Window shell (`XvueApp` + `XvueWindow` + `XvueCanvas`)" — phase boundary, goal, depends-on, success criteria 1–5
- `.planning/REQUIREMENTS.md` §"Shell — QApplication singleton and window" SHELL-01 through SHELL-07 — the 7 requirements this phase delivers
- `.planning/phases/00-build-skeleton-abi-stubs/00-CONTEXT.md` — Phase 0 decisions that Phase 1 inherits, specifically D-01 (xvue/qt/ layout), D-02 (read-only scope), D-03 (README_COORDS.md location), D-04 (single ABI header), D-05 (trailing-underscore macro), D-06 (pointer arg convention), D-08/D-09 (shell-script clone pattern), D-12 (verify_abi custom target — template for D-10), D-17 (warn-once stub pattern), and the Claude-discretion note declaring the `XVUE_QT_ASSERT_MAIN_THREAD()` macro skeleton

### Research synthesis
- `.planning/research/ARCHITECTURE.md` §"Components" — the `XvueApp` / `XvueWindow` / `XvueCanvas` / `XvueState` split that Phase 1 instantiates for the first time
- `.planning/research/ARCHITECTURE.md` §"Singleton discipline for `QApplication`" — the `std::call_once` + fabricated `argc`/`argv` + `atexit` pattern, documented as the Qt embedding idiom (HIGH confidence)
- `.planning/research/ARCHITECTURE.md` §"Recommended Project Structure" — the per-component header/source file layout the planner may follow
- `.planning/research/PITFALLS.md` §"Pitfall 5: `QApplication` double-initialization / lifetime leak" — the exact failure Phase 1 must prevent and the code pattern to use
- `.planning/research/PITFALLS.md` §"Pitfall 6: `QApplication::exec()` top-level call inverts control flow" — drives D-10/D-11 enforcement, and the `processEvents(ExcludeUserInputEvents)` rule in D-01
- `.planning/research/PITFALLS.md` §"Pitfall 7: `processEvents` starvation" — informs D-01's choice of `ExcludeUserInputEvents` to avoid re-entrancy surprises
- `.planning/research/PITFALLS.md` §"Pitfall 8: Modal `QDialog::exec()` re-entrancy" — not directly in Phase 1 scope, but informs why D-01 uses `ExcludeUserInputEvents`
- `.planning/research/STACK.md` — Qt 6 CMake integration snippets, `find_package(Qt6 COMPONENTS …)`, the "what NOT to use" list (no Qt 5, no QML, no qmake) that Phase 1 inherits unchanged from Phase 0

### Codebase maps
- `.planning/codebase/ARCHITECTURE.md` — the shared-data layer (`incl/` common blocks, read-only for Phase 1) and the launcher/entry-point/solver-module/utility/graphics layering
- `.planning/codebase/STRUCTURE.md` — `xvue/` directory layout, `prpr/` entry-point conventions (relevant for D-19 `xvtest0.f` placement)
- `.planning/codebase/STACK.md` — gfortran/gcc versions, Qt package names on Debian trixie
- `.planning/codebase/CONVENTIONS.md` — `bin/cb*` shell script conventions that `bin/cbxvtest0_qt` (D-20) must follow (`#!/bin/bash`, language detection via `$MEFISTO/td/m/anglais`, `cd $MEFISTO` first, per-step echo output)
- `.planning/codebase/CONCERNS.md` — stale `.o` fragility (drives the clean-build discipline Phase 1 inherits from Phase 0 D-10)

### Direct source — read these literally
- `/home/drico/git/mefisto/xvue/xvuelc.c:286-303` — legacy `xvinitgraphique_` body (opens X11 display, not a window). Phase 1 replaces this entry point semantically: the Qt version opens a window in addition to the app-level setup.
- `/home/drico/git/mefisto/xvue/xvuelc.c:319-334` — legacy `xvpxecran_` body (pixel dimensions via `XDisplayHeight`/`XDisplayWidth`). Phase 1 D-16 swaps to `QGuiApplication::primaryScreen()->size()`.
- `/home/drico/git/mefisto/xvue/xvuelc.c:337-356` — legacy `xvmmecran_` body (millimetre dimensions). Phase 1 D-17 swaps to `QGuiApplication::primaryScreen()->physicalSize()`.
- `/home/drico/git/mefisto/xvue/xvuelc.c:612-980` (approx.) — legacy `xvinfo_` body, including the `XCreateWindow`/`XMapWindow` block at lines ~943–969. Phase 1 D-03 cherry-picks only the window-sizing semantics (`*ix`, `*iy`) and leaves the palette/font output pointers on the warn-once path.
- `/home/drico/git/mefisto/xvue/xvuelc.c:1434` — legacy `xvfond_` body. Phase 1 D-14 replicates only the background-color semantics with a 2-entry table; full palette support is Phase 3.
- `/home/drico/git/mefisto/xvue/xvuelc.c:1602` — legacy `xvfermer_` body. Phase 1 D-06 replaces the X11 close path with `window_.reset()`.
- `/home/drico/git/mefisto/xvue/xvuelc.c:935` — `attributes.background_pixel = BlackPixel(display_mef, screen_mef)` — the legacy default background that justifies `QColor = Qt::black` in D-04.
- `/home/drico/git/mefisto/xvue/README_COORDS.md` — Y-axis convention audit from Phase 0 D-03. Phase 1 is not yet drawing anything so the convention is not consumed here, but the file is referenced for completeness and every Phase 2+ drawing PR must read it.
- `/home/drico/git/mefisto/xvue/qt/include/xvue_qt_api.h` — the Phase 0 ABI header. Phase 1 adds the fleshed-out `XVUE_QT_ASSERT_MAIN_THREAD()` macro body here.
- `/home/drico/git/mefisto/xvue/qt/src/xvue_qt_api.cpp` — the Phase 0 warn-once stub file. Phase 1 rewrites the seven SHELL-* entry point bodies to real implementations (D-01, D-03 partial, D-06, D-14, D-16, D-17, plus SHELL-07 retrofit) and bulk-inserts `XVUE_QT_ASSERT_MAIN_THREAD()` into every other stub.
- `/home/drico/git/mefisto/prpr/xvtest1.f` — reference for the Fortran test-driver convention used by D-19's new `prpr/xvtest0.f`.
- `/home/drico/git/mefisto/bin/cbmail_qt` (and sibling `bin/cb*_qt` Phase 0 clones) — template for D-20's new `bin/cbxvtest0_qt` build script.
- `/home/drico/git/mefisto/CLAUDE.md` — working rules: never break the X11 build, run the smallest relevant test after every change, commit after every logical step.

### External references
- Qt 6 `QApplication` / `QCoreApplication::processEvents` docs — authoritative source for the `QEventLoop::ExcludeUserInputEvents` flag semantics (`doc.qt.io/qt-6/qcoreapplication.html#processEvents`) and the `std::call_once` embedding idiom
- Qt 6 `QScreen` / `QGuiApplication::primaryScreen` docs — authoritative source for `size()` (logical px) vs. `physicalSize()` (mm) and HiDPI scaling behaviour (`doc.qt.io/qt-6/qscreen.html`)
- Qt 6 `QMainWindow` / `QWidget` documentation — authoritative source for `show()`, `close()`, `WA_DeleteOnClose` semantics that inform D-06

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **`xvue/qt/include/xvue_qt_api.h`** (Phase 0 D-04) — the single ABI header. Phase 1 adds the fleshed-out `XVUE_QT_ASSERT_MAIN_THREAD()` macro body here; no new headers are added to the ABI surface. Any new component headers (`xvue_qt_app.h`, `xvue_qt_window.h`, `xvue_qt_canvas.h`, `xvue_qt_state.h`) are internal to `xvue/qt/` and not part of the Fortran-facing ABI.
- **`xvue/qt/src/xvue_qt_api.cpp`** (Phase 0 D-17) — the single stub file. Phase 1 replaces the bodies of the seven SHELL entry points and bulk-inserts the main-thread assertion into every other stub. The warn-once `stderr` pattern from Phase 0 is preserved for all stubs that remain unimplemented (most of them).
- **`xvue/qt/build/libxvueqt.a`** (Phase 0 D-13) — the handoff artifact. Phase 1 grows its size (new `.o` files for app/window/canvas/state sources) but does not change the CMake export path or the `bin/cb*_qt` linker lines. Phase 0 D-12's `verify_abi` custom target keeps checking that the trailing-underscore symbol count matches `xvue_qt_api.h`.
- **`prpr/xvtest{1,2,3,4}.f`** — existing Fortran test-driver convention reserved for Phase 2's drawing primitives (per ROADMAP Phase 2 goal mentioning them). Phase 1 adds a sibling `prpr/xvtest0.f` slot without touching the existing four files.
- **`bin/cbmail_qt`, `bin/cbelas_qt`, `bin/cbflui_qt`, `bin/cbther_qt`, `bin/cbnlse_qt`** (Phase 0 D-09) — the five per-module Qt-variant compile scripts. Phase 1 clones the thinnest one (likely `cbmail_qt`) into `bin/cbxvtest0_qt`, adjusting source file name, executable name, and the invoked CMake step (unchanged), but nothing else.
- **`xvue/xvuelc.c`** — the legacy X11 backend, untouched by Phase 1 and still built by the legacy `bin/cbl_tout` path. Serves as the ground-truth reference for the semantic intent of every SHELL entry point.
- **`xvue/README_COORDS.md`** (Phase 0 D-03) — Y-axis convention audit. Not consumed by Phase 1 code but referenced for completeness.
- **`.planning/validation/BASELINE.md`** (Phase 0 D-14) — the 5 canonical `testa/` cases. Not re-run in Phase 1 (see D-21); resumes from Phase 2.

### Established Patterns
- **Warn-once stub discipline** (Phase 0 D-17) — Phase 1's `xvinfo_` partial implementation keeps the warn-once `stderr` print for the non-sizing output parameters, following the exact pattern of every Phase 0 stub. Downstream phases will replace these messages one by one.
- **Single-file ABI surface** (Phase 0 D-04) — Phase 1 adds no new headers to `xvue/qt/include/`; all new Qt-internal headers live alongside their `.cpp` in `xvue/qt/src/` or are introduced in an unexposed internal include directory at the planner's discretion.
- **CMake custom targets for build-time guards** (Phase 0 D-12) — Phase 1's `verify_no_exec` target (D-10) follows the `verify_abi` precedent: a post-build grep-and-fail, invoked as part of the standard `xvueqt` target, with a clear error message naming the offending file and line. No new enforcement philosophy.
- **Shell-script clone-and-modify** (Phase 0 D-09, Phase 0 Specifics "Keep the shell-script convention") — `bin/cbxvtest0_qt` (D-20) is a literal clone of an existing `bin/cb*_qt` script, not a parametric rewrite. No `case $BACKEND in` branches anywhere.
- **Trailing-underscore Fortran ABI** (Phase 0 D-05) — every Phase 1 entry-point body is named via the same `#define proc(x) x##_` macro inherited from `xvuelc.c`. Byte-identical to the legacy backend.
- **Fortran `.f` fixed-form files in `prpr/`** — `prpr/xvtest0.f` follows the column-7+ fixed-form convention (CLAUDE.md §"Language and module conventions") and the project's normes from `doc/normes.ps`.

### Integration Points
- **`xvinitgraphique_` as the Fortran→Qt entry point for lifecycle** — every Fortran caller that wants to open a display goes through this symbol. Phase 1 makes it a real implementation that creates the `QApplication` and the `XvueWindow`; Phase 2+ entry points assume `XvueApp::ensure()` has already run by the time they are called.
- **`XvueCanvas::paintEvent`** — the single surface where Phase 1's "fill the background" and Phase 2's "blit the backing pixmap" meet. Phase 2 replaces the body as a one-line swap (`fillRect` → `drawPixmap`); the rest of Phase 1's `XvueCanvas` plumbing survives the transition.
- **`XvueState` struct** — Phase 1 introduces it with one field; Phase 2 adds pen/brush/font/painter; Phase 3 adds the indexed palette; Phase 4 adds the `XvuePixmapStack` slot ownership. Every phase grows this struct additively without renaming fields introduced earlier.
- **Shell-script linker line in `bin/cb*_qt`** — unchanged from Phase 0. The only new consumer of the linker line in Phase 1 is `bin/cbxvtest0_qt` (D-20), which uses the same `-Lxvue/qt/build -lxvueqt $(pkg-config --libs Qt6Widgets Qt6Gui Qt6Core Qt6PrintSupport) -lstdc++` block as the five canonical module scripts.
- **CMake `xvue/qt/CMakeLists.txt`** — grows by: (a) new source files for `xvue_qt_app.cpp`, `xvue_qt_window.cpp`, `xvue_qt_canvas.cpp`, `xvue_qt_state.cpp` (planner-chosen file split); (b) a new `verify_no_exec` custom target mirroring `verify_abi` from Phase 0. No change to `cmake_minimum_required`, `CMAKE_CXX_STANDARD`, or `find_package(Qt6 …)` calls.
- **`atexit`-registered QApplication teardown** — installed on first `XvueApp::ensure()` call. This is the single point where the `QApplication` is destroyed; no code path inside any `extern "C"` entry point (including `xvfermer_`) touches the `QApplication` lifetime.

</code_context>

<specifics>
## Specific Ideas

- **"Match the legacy X11 default"** — D-04 chooses `Qt::black` as the default background specifically because `xvue/xvuelc.c:935` uses `BlackPixel` as `attributes.background_pixel`. This preserves visual parity for the eventual Phase 8 A/B validation and avoids an avoidable drift source.
- **"Keep the warn-once trail alive through Phase 1"** — every stub that is not one of the seven SHELL entry points (or `xvinfo_`'s partial implementation) continues to print `"xvue-qt: stub <fn> not implemented yet"` on first call. The Phase 0 diagnostic philosophy (visible without being noisy) is unchanged. Phase 1 does not silence any existing warning.
- **"The reopen cycle is a first-class test"** — `prpr/xvtest0.f` calls `XVINITGRAPHIQUE` and `XVFERMER` twice by design (not once). The second call is where `QApplication: there can only be one` would crash, and where the `std::call_once` guard earns its keep. The test exists to make that failure mode reproducible from the Fortran side, not just under C++ unit tests.
- **"One CMake file, one ABI header, one stub file — still"** — Phase 1 adds new `.cpp` and `.h` files under `xvue/qt/src` / `xvue/qt/include` for the `XvueApp`/`XvueWindow`/`XvueCanvas`/`XvueState` components, but the ABI surface remains exactly one header (`xvue_qt_api.h`) and the stub file for unimplemented entry points remains exactly one `.cpp` (`xvue_qt_api.cpp`). The "grep/diff against `xvuelc.c`" property from Phase 0 B1 is preserved.
- **"`ExcludeUserInputEvents` is non-negotiable in D-01"** — this comes directly from Pitfall 6 and Pitfall 8, and is the mechanism that lets Phase 1 call `processEvents` once without opening a door that Phase 5's nested `QEventLoop` pattern would have to close again.

</specifics>

<deferred>
## Deferred Ideas

- **Multi-monitor awareness for `xvpxecran_`/`xvmmecran_`** — D-16/D-17 explicitly pick `QGuiApplication::primaryScreen()` and defer a `window->screen()` path. Revisit only if a concrete multi-monitor scenario surfaces during Phase 8 A/B validation.
- **`QPalette`-based background handling** — considered as an alternative to `XvueState::background_`, rejected for Phase 1 because Phase 2's backing pixmap `drawPixmap(0,0)` would cover the palette. `XvueState` field is cleaner and survives the Phase 2 transition.
- **Git pre-commit hook for `QApplication::exec` guard** — rejected (D-11) in favor of a CMake-only enforcement. Reconsider if the project ever grows CI, but not before.
- **Adding `pp/ppxvtest0_qt` to `bin/cbl_tout_qt`** — rejected (D-20) to keep the "big script builds the five canonical modules" contract intact. The Phase 1 test driver is built via its own standalone `bin/cbxvtest0_qt`.
- **Defensive `Qt::AA_EnableHighDpiScaling` attribute call** — Qt 6 enables HiDPI scaling by default; the attribute is a no-op on Qt 6 and merely signals intent. Planner may add it for documentation purposes but it is not required for SHELL-06.
- **Full `XvueState` struct (pen, brush, font, painter, palette)** — Phase 2 (pen/brush/painter), Phase 3 (font/palette), Phase 4 (pixmap slots). Phase 1 holds only `background_`.
- **`WA_DeleteOnClose` + `close()` idiom for `xvfermer_`** — considered as a Qt-idiomatic alternative to `window_.reset()`. Rejected because the timing of the actual deletion becomes async (posted to the event queue) and Phase 1 wants a synchronous teardown that is observable immediately after `xvfermer_` returns. Straight `reset()` is simpler and deterministic.
- **Running the 5 baseline `testa/` cases through `pp/pp*_qt` at Phase 1 end** — deferred (D-21). Nothing in Phase 1 changes legacy X11 behaviour and the Qt backend has no drawing yet, so the A/B re-run would produce the same blank-window result for all five cases. Resumes at Phase 2.
- **Retire `prpr/xvtest0.f` after Phase 1** — not deferred, just noted: `xvtest0.f` stays in the tree as a permanent shell sanity driver. Phase 2+ can extend it with drawing calls or leave it as a minimal lifecycle test.

### Reviewed Todos (not folded)
None — no pending todos matched Phase 1 scope at init time.

</deferred>

---

*Phase: 01-window-shell-xvueapp-xvuewindow-xvuecanvas*
*Context gathered: 2026-04-11*
