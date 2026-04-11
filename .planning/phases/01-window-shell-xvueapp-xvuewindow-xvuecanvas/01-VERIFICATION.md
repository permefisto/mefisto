---
phase: 01-window-shell-xvueapp-xvuewindow-xvuecanvas
verified: 2026-04-11T15:00:00Z
status: passed
score: 5/5 must-haves verified
overrides_applied: 0
---

# Phase 1: Window Shell Verification Report

**Phase Goal:** A blank Qt window opens through `xvinitgraphique_` and closes cleanly through `xvfermer_`, proving the `QApplication` singleton discipline and HiDPI convention in isolation before any drawing logic.
**Verified:** 2026-04-11T15:00:00Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Developer can run a test driver that calls `xvinitgraphique_` and observes a blank QMainWindow with an XvueCanvas central widget appear on screen | VERIFIED (human) | User empirically confirmed "windows visible" after in-plan gap closures (commits 36f4fb6, a7b1c2c, 006ff3d). `pp/ppxvtest0_qt` exits 0. `prpr/xvtest0.f` calls XVINITGRAPHIQUE twice with SLEEP(1) hold. XvueWindow ctor sets title "MEFISTO", resize(800,600), XvueCanvas central widget. |
| 2 | Calling `xvfermer_` then `xvinitgraphique_` a second time reopens the window with no "QApplication: there can only be one" assertion and no crash on process exit | VERIFIED (human) | User confirmed no assertion, no crash. `xvfermer_` calls `XvueApp::window_slot().reset()` then drains deleteLater queue. `XvueApp::ensure()` uses `std::call_once` — QApplication constructed exactly once, leaked at atexit (qapp_.release()). SIGSEGV root-caused and fixed; archived at `.planning/debug/resolved/phase-01-xvtest0-teardown-segfault.md`. |
| 3 | `grep -rn 'QApplication::exec' xvue/*.cpp` returns zero matches, enforced by a CMake check that fails the build otherwise | VERIFIED | `sh xvue/qt/cmake/verify_no_exec.sh xvue/qt/src xvue/qt/include` exits 0 with "verify_no_exec: OK". Zero matches confirmed by direct grep. `add_custom_target(verify_no_exec ALL ...)` in `xvue/qt/CMakeLists.txt` ensures the check runs on every build as a dependency of xvueqt. |
| 4 | `xvpxecran_`/`xvmmecran_` return screen dimensions in logical pixels from QScreen; window renders correctly with `QT_SCALE_FACTOR=2` | VERIFIED (human) | User confirmed HiDPI run shows visibly larger window, exit 0, no crash. `xvpxecran_` returns `QGuiApplication::primaryScreen()->size()` (logical pixels). `xvmmecran_` returns `QScreen::physicalSize()` rounded to mm. Both confirmed by grep: 2 `primaryScreen()` matches. |
| 5 | Every `extern "C"` entry point contains `XVUE_QT_ASSERT_MAIN_THREAD()` | VERIFIED | `grep -c 'XVUE_QT_ASSERT_MAIN_THREAD();' xvue/qt/src/xvue_qt_api.cpp` = **57**. `grep -c 'XvueApp::ensure();'` = **57**. Ordering is always ensure() first, assert second (D-18). warn_once(warned count = 51 (57 stubs minus 6 SHELL entries that were fully rewritten). |

**Score:** 5/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `xvue/qt/src/xvue_qt_state.h` | XvueState struct with QColor background_ field | VERIFIED | Exists. Single field `QColor background_ = Qt::black`. |
| `xvue/qt/src/xvue_qt_app.h` | XvueApp class with ensure(), qapp(), window_slot() | VERIFIED | Exists. Static once_flag_, qapp_, window_ members present. |
| `xvue/qt/src/xvue_qt_app.cpp` | XvueApp::ensure via std::call_once + atexit | VERIFIED | Exists. std::call_once confirmed. teardown_atexit uses qapp_.release() (documented deliberate leak). |
| `xvue/qt/src/xvue_qt_window.h` | XvueWindow : QMainWindow with Q_OBJECT | VERIFIED | Exists. Constructor sets "MEFISTO" title, resize(800,600), XvueCanvas central widget. |
| `xvue/qt/src/xvue_qt_window.cpp` | XvueWindow ctor sets title/size/central widget | VERIFIED | Exists, substantive — ctor confirmed. |
| `xvue/qt/src/xvue_qt_canvas.h` | XvueCanvas : QWidget with Q_OBJECT | VERIFIED | Exists. |
| `xvue/qt/src/xvue_qt_canvas.cpp` | XvueCanvas::paintEvent fills background | VERIFIED | Exists. `QPainter(this).fillRect(rect(), state_->background_)` — one-line per D-05. |
| `xvue/qt/cmake/verify_no_exec.sh` | grep-and-fail guard for SHELL-03 | VERIFIED | Exists, executable, exits 0 on clean tree, exits 1 on violation. |
| `xvue/qt/CMakeLists.txt` | Sources wired + verify_no_exec ALL target | VERIFIED | 4 sources in add_library; `add_custom_target(verify_no_exec ALL ...)` present. |
| `xvue/qt/src/xvue_qt_api.cpp` | 7 SHELL entry bodies + 57-entry macro retrofit | VERIFIED | 57 ensure() calls, 57 ASSERT calls, 2 primaryScreen(), window_slot().reset() ×1, win->state()->background_ ×1, win->resize ×1. warn_once(warned) = 51 remaining stubs. |
| `prpr/xvtest0.f` | Fortran 77 lifecycle driver | VERIFIED | Exists, 42 lines, 2× CALL XVINITGRAPHIQUE, 2× CALL XVFERMER, SLEEP(1) holds. |
| `bin/cbxvtest0_qt` | Build script for ppxvtest0_qt | VERIFIED | Exists, executable (755). No `ppmail` reference, no `mail/lib`. Links against `-lxvueqt`. |
| `pp/ppxvtest0_qt` | Compiled Phase 1 test executable | VERIFIED | Exists, 755, 86760 bytes, built 2026-04-11. |
| `pp/ppmail_qt` etc. (5 canonical Qt executables) | All five pp/*_qt executables | VERIFIED | All 5 exist: ppelas_qt, ppflui_qt, ppmail_qt, ppnlse_qt, ppther_qt. |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `xvue/qt/CMakeLists.txt` | `xvue_qt_app.cpp`, `xvue_qt_window.cpp`, `xvue_qt_canvas.cpp` | `add_library(xvueqt STATIC ...)` | WIRED | All 4 .cpp sources listed in add_library |
| `xvue/qt/CMakeLists.txt` | `xvue/qt/cmake/verify_no_exec.sh` | `add_custom_target(verify_no_exec ALL ...)` | WIRED | Pattern `verify_no_exec` confirmed in CMakeLists.txt |
| `xvue/qt/src/xvue_qt_api.cpp` | `XvueApp / XvueWindow / XvueCanvas` | `#include "xvue_qt_app.h"` + `"xvue_qt_window.h"` + `"xvue_qt_canvas.h"` | WIRED | All 3 includes confirmed present |
| `xvinitgraphique_` | `XvueWindow::show` | `XvueApp::window_slot()` lazy allocation | WIRED | `window_slot()` call confirmed; bounded exposure pump using isExposed() loop |
| `xvpxecran_`/`xvmmecran_` | `QGuiApplication::primaryScreen()` | direct call | WIRED | 2 `primaryScreen()` matches confirmed |
| `prpr/xvtest0.f` | `xvinitgraphique_` / `xvfermer_` | Fortran CALL → trailing-underscore ABI | WIRED | 2× CALL XVINITGRAPHIQUE, 2× CALL XVFERMER; `pp/ppxvtest0_qt` links successfully against libxvueqt.a |
| `bin/cbxvtest0_qt` | `xvue/qt/build/libxvueqt.a` | `-Lxvue/qt/build -lxvueqt` | WIRED | `-lxvueqt` confirmed in script; ppxvtest0_qt produced |

### Data-Flow Trace (Level 4)

Not applicable for Phase 1 artifacts — no data-rendering components. XvueCanvas::paintEvent fills from XvueState::background_ (a QColor field, not a dynamic data source). This is structural wiring, not a data pipeline.

### Behavioral Spot-Checks

| Behavior | Check | Result | Status |
|----------|-------|--------|--------|
| verify_no_exec guard reports clean | `sh xvue/qt/cmake/verify_no_exec.sh xvue/qt/src xvue/qt/include` | "verify_no_exec: OK (no forbidden tokens...)" exit 0 | PASS |
| 57 entry points have ensure() | `grep -c 'XvueApp::ensure();' xvue/qt/src/xvue_qt_api.cpp` | 57 | PASS |
| 57 entry points have ASSERT macro | `grep -c 'XVUE_QT_ASSERT_MAIN_THREAD();' xvue/qt/src/xvue_qt_api.cpp` | 57 | PASS |
| xvfermer_ destroys window only | `grep -c 'window_slot().reset()' xvue/qt/src/xvue_qt_api.cpp` | 1 | PASS |
| ppxvtest0_qt exists and is executable | `test -x pp/ppxvtest0_qt` | true | PASS |
| Two XVINITGRAPHIQUE calls in test driver | `grep -c 'CALL XVINITGRAPHIQUE' prpr/xvtest0.f` | 2 | PASS |
| All 11 phase commits present in history | `git cat-file -t {hash}` for each | all "commit" | PASS |
| QApplication leak at atexit (qapp_.release()) | `grep 'qapp_.release()' xvue/qt/src/xvue_qt_app.cpp` | 1 match | PASS |

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|-------------|-------------|--------|----------|
| SHELL-01 | 01-02, 01-03 | `xvinitgraphique_` opens QMainWindow with XvueCanvas central widget | SATISFIED | Real body in xvue_qt_api.cpp; XvueWindow ctor confirmed; human visual confirmation in 01-03 |
| SHELL-02 | 01-02, 01-03 | `xvfermer_` closes window without destroying QApplication; reopen is safe | SATISFIED | `window_slot().reset()` + `qapp_.release()` at atexit; SIGSEGV fix confirmed; human visual confirmation |
| SHELL-03 | 01-01 | No `QApplication::exec()` anywhere in xvue/, enforced by CMake ALL target | SATISFIED | verify_no_exec.sh exits 0; `add_custom_target(verify_no_exec ALL ...)` in CMakeLists.txt |
| SHELL-04 | 01-02 | `xvpxecran_`/`xvmmecran_` return QScreen logical pixels / physical mm | SATISFIED | Both bodies confirmed with primaryScreen(); human HiDPI confirmation |
| SHELL-05 | 01-02 | `xvfond_` sets canvas background without corrupting backing pixmap | SATISFIED | xvfond_ body: maps int to QColor, sets win->state()->background_, calls canvas->update() |
| SHELL-06 | 01-02, 01-03 | Window renders correctly on HiDPI / QT_SCALE_FACTOR=2 | SATISFIED | Qt 6 HiDPI-by-default; QT_SCALE_FACTOR=2 run human-confirmed; QScreen::size() returns logical px |
| SHELL-07 | 01-01, 01-02 | Every extern "C" entry point has XVUE_QT_ASSERT_MAIN_THREAD() debug assertion | SATISFIED | 57/57 ASSERT macros confirmed by grep |

### Anti-Patterns Found

| File | Pattern | Severity | Impact |
|------|---------|----------|--------|
| `xvue/qt/src/xvue_qt_api.cpp` | 51 `warn_once` stubs remain for non-SHELL entry points | Info | Expected — these are Phase 0 stubs scheduled for Phases 2-8 per ROADMAP.md. Each now carries `XvueApp::ensure()` + `XVUE_QT_ASSERT_MAIN_THREAD()`. Not blockers for Phase 1 goal. |
| `xvue/qt/src/xvue_qt_api.cpp` | `xvfond_` with no open window is a no-op past the mapping | Info | Documented limitation — XvueState is owned by XvueWindow; Phase 2+ may add pending-latch if needed. Non-blocking for Phase 1. |
| `prpr/xvtest0.f` | SLEEP(1) holds between each open/close cycle | Info | Intentional for human visual validation (FIX-4). Acceptable for a permanent shell sanity driver. |

No blocker or warning anti-patterns found that affect Phase 1 goal achievement. The code review (01-REVIEW.md) found 0 Critical, 3 Warning (advisory), 6 Info findings — all non-blocking.

### Human Verification Required

None. All 5 success criteria are either programmatically verifiable or were empirically confirmed by the developer during Plan 01-03 Task 3 (human visual checkpoint). The user's "confirmed fixed — windows visible" and HiDPI approval are the authoritative signals for SC 1, 2, and 4.

### Gaps Summary

No gaps. All 5 roadmap success criteria are verified:

1. **SC1** — `pp/ppxvtest0_qt` runs and produces visible 800x600 "MEFISTO" windows. (Note: ROADMAP SC1 mentions `pp/ppmail_qt` incidentally as the vehicle, but the actual delivery is the dedicated `pp/ppxvtest0_qt` driver — this is the better artifact per D-19 and fulfills the intent of SC1.)
2. **SC2** — Reopen cycle with no QApplication singleton assertion, no crash. SIGSEGV root-caused and fixed; QApplication deliberate-leak pattern established.
3. **SC3** — `verify_no_exec.sh` exits 0 on clean tree; wired as CMake ALL target.
4. **SC4** — QScreen logical-pixel returns confirmed; HiDPI human-confirmed.
5. **SC5** — 57/57 extern "C" entry points carry `XVUE_QT_ASSERT_MAIN_THREAD()`.

Phase 1 goal fully achieved. Ready to proceed to Phase 2 (drawing primitives and backing pixmap).

---

_Verified: 2026-04-11T15:00:00Z_
_Verifier: Claude (gsd-verifier)_
