---
phase: 01-window-shell-xvueapp-xvuewindow-xvuecanvas
reviewed: 2026-04-11T00:00:00Z
depth: standard
files_reviewed: 13
files_reviewed_list:
  - bin/cbxvtest0_qt
  - prpr/xvtest0.f
  - xvue/qt/CMakeLists.txt
  - xvue/qt/cmake/verify_no_exec.sh
  - xvue/qt/include/xvue_qt_api.h
  - xvue/qt/src/xvue_qt_api.cpp
  - xvue/qt/src/xvue_qt_app.cpp
  - xvue/qt/src/xvue_qt_app.h
  - xvue/qt/src/xvue_qt_canvas.cpp
  - xvue/qt/src/xvue_qt_canvas.h
  - xvue/qt/src/xvue_qt_state.h
  - xvue/qt/src/xvue_qt_window.cpp
  - xvue/qt/src/xvue_qt_window.h
findings:
  critical: 0
  warning: 3
  info: 6
  total: 9
status: issues_found
---

# Phase 1: Code Review Report

**Reviewed:** 2026-04-11
**Depth:** standard
**Files Reviewed:** 13
**Status:** issues_found

## Summary

Phase 1 delivers a coherent Qt 6 window shell: XvueApp singleton with process-lifetime QApplication (deliberately leaked, well documented in xvue_qt_app.cpp:1-14,47-69), XvueWindow/XvueCanvas/XvueState trio wired through a raw back-pointer (D-15), and 57 extern "C" ABI stubs each preceded by `XvueApp::ensure()` + `XVUE_QT_ASSERT_MAIN_THREAD()`. The no-exec policy is enforced in CMake via verify_no_exec.sh, and xvtest0.f exercises the reopen path twice as intended.

No Critical issues found. The main correctness gap is that C++ exceptions can cross the extern "C" boundary into Fortran in two places (xvinitgraphique_ allocation and the XvueApp::ensure() first-call path) — that is undefined behaviour and should be wrapped. A handful of Info items concern dead code in the proc-macro header block, the bash driver script's handling of an unset `$MEFISTO`, a missing explicit `<QColor>` include, and a suspected non-ASCII em-dash in the Fortran source comments that may conflict with `doc/normes.ps`.

## Warnings

### WR-01: C++ exceptions may propagate through extern "C" boundary

**File:** `xvue/qt/src/xvue_qt_api.cpp:76-107` (xvinitgraphique_), and implicitly every stub via `XvueApp::ensure()`
**Issue:** `std::make_unique<XvueWindow>()` (line 84) can throw `std::bad_alloc`; `QMainWindow`/`XvueCanvas` constructors can throw; and inside `XvueApp::ensure()` the `std::make_unique<QApplication>(...)` call (xvue_qt_app.cpp:34) can also throw. These entry points are `extern "C"` and are called directly from Fortran. Throwing a C++ exception through a `extern "C"` frame into a Fortran caller is undefined behaviour on every toolchain (gfortran has no C++ EH personality in its frames). On libstdc++/Linux the usual outcome is `std::terminate` with a confusing backtrace, but it can also corrupt the unwinder state.
**Fix:** Wrap every extern "C" body (or at minimum the ones that allocate — `xvinitgraphique_`, and `XvueApp::ensure()`) in a `try { ... } catch (const std::exception& e) { std::fprintf(stderr, "xvue-qt: fatal in <name>: %s\n", e.what()); std::abort(); } catch (...) { std::fprintf(stderr, "xvue-qt: fatal unknown exception in <name>\n"); std::abort(); }` shim. A small `XVUE_QT_EXTERN_C_BODY(name) { ... }` macro or a helper lambda keeps the 57 stubs readable. Calling `std::abort()` is preferable to letting the exception unwind across the C boundary — the Fortran caller will see a clean crash rather than UB.

### WR-02: `cd $MEFISTO` unguarded in bin/cbxvtest0_qt

**File:** `bin/cbxvtest0_qt:7`
**Issue:** `cd $MEFISTO` is unquoted and unchecked. If `MEFISTO` is unset (the CLAUDE.md environment variable requirement is not enforced at runtime), `cd` with no argument changes to `$HOME`, then line 11 runs `mkdir pp` in the user's home directory, and line 39 tries to compile `prpr/xvtest0.f` from `$HOME` — polluting the user's home and producing misleading "file not found" errors. The surrounding cb* scripts in `bin/` have the same shape, but new scripts should fail loud rather than propagate the legacy pattern.
**Fix:**
```bash
if [ -z "${MEFISTO:-}" ]; then
  echo "cbxvtest0_qt: MEFISTO environment variable is not set" >&2
  exit 1
fi
cd "$MEFISTO" || { echo "cbxvtest0_qt: cannot cd to \$MEFISTO=$MEFISTO" >&2; exit 1; }
```
Also quote the `$nompp` variable in the `rm` and `test -f` calls, and consider `set -eu` at the top now that this script is new (the legacy scripts predate that idiom and are out of scope).

### WR-03: xvinitgraphique_ silently returns on exposure timeout

**File:** `xvue/qt/src/xvue_qt_api.cpp:97-106`
**Issue:** The bounded-loop waits up to 2 s for `QWindow::isExposed()`, but on timeout the function returns with no diagnostic. The Fortran caller cannot distinguish "window is on screen" from "X server was wedged and we gave up" — both look like a successful init, and the next drawing call will land on an invisible (or not-yet-mapped) window. In Phase 1 the consequence is limited to the visual sleep test, but the same code will be reused by every solver.
**Fix:** After the loop, check the condition and emit a warn-once stderr line if `!(wh && wh->isExposed())`. Do NOT abort — the test explicitly requires the function to be non-blocking — but a single `xvue-qt: xvinitgraphique_ timed out waiting for expose (2000 ms)` line makes a broken display station debuggable. Consider the same for a null `windowHandle()` at loop exit.

## Info

### IN-01: Dead `#ifdef __GNUC__` block in xvue_qt_api.h proc macro

**File:** `xvue/qt/include/xvue_qt_api.h:21-27`
**Issue:** The header defines `proc(x)` conditionally under `__GNUC__`, then unconditionally `#undef`s and redefines it to the same GCC form. The `#ifdef __GNUC__` / `#else` branch is dead code the moment lines 26-27 execute.
**Fix:** Drop lines 21-25 and keep only `#define proc(x) x##_`. Add a one-line comment explaining that trailing-underscore mangling is gfortran's default ABI.

### IN-02: Missing explicit `#include <QColor>` in xvue_qt_api.cpp

**File:** `xvue/qt/src/xvue_qt_api.cpp:300` (use of `QColor`)
**Issue:** `xvfond_` constructs a `QColor`, but the file only transitively gets the type through `<QApplication>`/`<QGuiApplication>`. Any future Qt reorganisation that stops making `QColor` a transitive include would silently break this translation unit. A local rule elsewhere in the codebase (and Qt style-guide guidance) is to include what you use.
**Fix:** Add `#include <QColor>` near the other Qt includes (alphabetised with the existing list).

### IN-03: Unused pre-loop `windowHandle()` assignment

**File:** `xvue/qt/src/xvue_qt_api.cpp:100`
**Issue:** `QWindow* wh = win->windowHandle();` is read back on line 104 inside the loop before it is ever observed. The pre-loop assignment is dead and only serves to give `wh` a declared type.
**Fix:** Replace with `QWindow* wh = nullptr;` to make the intent obvious, or declare and assign inside the loop and hoist only if the post-loop diagnostic from WR-03 needs it.

### IN-04: Non-ASCII em dash in prpr/xvtest0.f comments

**File:** `prpr/xvtest0.f:4,39`
**Issue:** Lines 4 ("Phase 1 — SHELL-01, SHELL-02") and 39 ("OK — cycle open/close/open/close") contain the Unicode em-dash (U+2014). Fixed-form Fortran 77 as specified in `doc/normes.ps` is traditionally ASCII-only; gfortran tolerates UTF-8 in comments in practice but this is a departure from the project's programming norm. The same concern does NOT apply to the C++/CMake files where UTF-8 comments are standard.
**Fix:** Replace `—` with `--` or `-` in the Fortran source. Verify `doc/normes.ps` explicitly if the maintainer wants to confirm, but the CLAUDE.md entry "The project's coding norms are documented in doc/normes.ps ... must be respected" makes the conservative choice ASCII.

### IN-05: verify_no_exec.sh will false-positive on documentation comments

**File:** `xvue/qt/cmake/verify_no_exec.sh:21`
**Issue:** The grep pattern matches the literal tokens `QApplication::exec` and `qApp->exec()` anywhere, including inside comments that deliberately explain "we must not call QApplication::exec here". The current source happens to avoid that phrasing, but any future doc comment that references the forbidden API (e.g. "Unlike the legacy path that used QApplication::exec, ...") will break the build. This creates a pressure to avoid precise documentation, which is the opposite of what the no-exec guard is for.
**Fix:** Either (a) strip C/C++ line and block comments before grepping (e.g. `sed 's|//.*||; s|/\*.*\*/||'` pipeline or a short awk), or (b) restrict the match to lines that look like actual calls, e.g. require the pattern NOT to be preceded by `//` on the same line: `grep -R -n -E '^[^/]*\b(QApplication::exec|qApp->exec)\(' ...`. Option (b) is simpler and usually sufficient.

### IN-06: `xvfermer_` destroys window even when slot is already empty

**File:** `xvue/qt/src/xvue_qt_api.cpp:347-359`
**Issue:** Minor: `XvueApp::window_slot().reset()` is a no-op when the slot is already null (double close, or close before open), but the subsequent `processEvents` still fires. Not a bug — it's cheap and idempotent — but a one-line early return after the warn-once pattern used elsewhere would keep the "close with no open" case symmetric with the rest of the stubs and avoid pumping events on a cold path.
**Fix:** Optional. `if (!XvueApp::window_slot()) return;` before the reset/processEvents pair. Leaving as-is is also defensible since the code is trivially correct.

---

_Reviewed: 2026-04-11_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
