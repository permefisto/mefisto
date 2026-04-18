---
phase: 05-event-bridge-blocking-reads
reviewed: 2026-04-18T00:00:00Z
depth: standard
files_reviewed: 20
files_reviewed_list:
  - bin/cbl_tout_qt
  - bin/xvtest-capture.sh
  - xvue/qt/CMakeLists.txt
  - xvue/qt/cmake/verify_abi.sh
  - xvue/qt/README.md
  - xvue/qt/src/xvue_qt_api.cpp
  - xvue/qt/src/xvue_qt_app.cpp
  - xvue/qt/src/xvue_qt_app.h
  - xvue/qt/src/xvue_qt_canvas.cpp
  - xvue/qt/src/xvue_qt_event.cpp
  - xvue/qt/src/xvue_qt_event.h
  - xvue/qt/src/xvue_qt_state.cpp
  - xvue/qt/src/xvue_qt_state.h
  - xvue/qt/src/xvue_qt_window.cpp
  - xvue/qt/src/xvue_qt_window.h
  - xvue/qt/tests/CMakeLists.txt
  - xvue/qt/tests/test_helpers.cpp
  - xvue/qt/tests/test_helpers.h
  - xvue/qt/tests/test_xvue_qt_event.cpp
  - xvue/xvuelc.c
findings:
  critical: 0
  warning: 4
  info: 6
  total: 10
status: issues_found
---

# Phase 5: Code Review Report

**Reviewed:** 2026-04-18T00:00:00Z
**Depth:** standard
**Files Reviewed:** 20
**Status:** issues_found

## Summary

Phase 5 lands a well-designed Qt 6 event-bridge replacing the X11 XVSOURIS
event loop. Architecture is clean: the RAII BlockingDepthGuard, the nested
QEventLoop pattern, the save-restore of 7 filter members for nested calls,
the QTimer::singleShot deferred-quit coalescing, and the Strategy B
save/restore accrochage all line up with the design documents. Memory
ownership is correctly expressed through Qt parent-child (bridge parented
to window, canvas parented to window via setCentralWidget) and through the
XvueState destructor ordering that tears saved_canvas_ / mempxaccro_ /
accroche_undo_tile_ before painter_ / backing_. ABI byte-compatibility is
confirmed by spot-check against `xvue/qt/include/xvue_qt_api.h` — all five
Phase 5 entries (xvsouris_, xvsouris2_, xvpause_, deplsouris_,
initaccrochage_) match the public header signatures, and AUTOEXIT parity
between the X11 xvpause_ body (xvuelc.c:2516-2553) and the Qt xvpause_
body (xvue_qt_api.cpp:1101-1125) is clean.

No critical issues were found. Four warnings document real
semantic-divergence or state-management concerns that are unlikely to fire
in current drivers but will surface under nested re-entrancy or
non-mainline X11 parity scenarios. Six info items capture
smaller-impact findings: dead-code test stubs, comment-vs-code mismatches,
and robustness hardening opportunities.

## Warnings

### WR-01: xvsouris2_ Souris2 middle-button press diverges from X11 semantics

**File:** `xvue/qt/src/xvue_qt_event.cpp:275-283`
**Issue:** In Souris2 mode, a middle-button `MouseButtonPress` is
short-circuited to `pending_.notypeevent = 0; pending_.nbc = 2;` followed
by `loop_->quit()`. The X11 reference body at `xvue/xvuelc.c:2383-2440`
treats ALL three buttons identically inside the accrochage block — middle
button simply sets `*ibutton = 2` and continues through the nearest-item
search, returning `*notypeevent = 5` on press (and the `@`-abort remap at
`xvue/saclav.f:83-86` only fires on the RELEASE path where `NOCODE == 1`).
With the Qt early-return, saclav.f sees `NOCODE = 0` (line 74), prints
"SACLAV: PROBLEME AVEC LA SOURIS" and loops back through XVPAUSE — a
visible behavior change for middle-button-initiated drags on the mesher.
**Fix:** Remove the middle-button special-case in the Souris2 press path;
let the normal `accrochage redraw + return notypeevent=5` code handle all
three buttons uniformly. The @-abort remap belongs to saclav.f on release,
matching X11 parity.
```cpp
// Remove lines 275-283. Fall through to the nearest-item search at 284-315
// unconditionally; let saclav.f's NOCODE==1 (release) branch handle the
// btn=2 → @ abort remap as it does on the X11 backend.
```

### WR-02: nested xvsouris2_ calls share accroche_undo_tile_ across frames

**File:** `xvue/qt/src/xvue_qt_event.cpp:107-123, 236-251, 425-477`
**Issue:** `state->accroche_undo_tile_` is a process-global on `XvueState`
(not per-call), while `waitForEvent` save/restores the 7 bridge members
including `items_` and `pmin0_`. If a Fortran path ever calls xvsouris2_
while another xvsouris2_ is already blocking (nested re-entrancy — rare
but legal per the bridge design), the inner call's `save_tile_under` at
line 107-123 reuses the existing `accroche_undo_tile_` without reallocating
(line 109 `if (!state->accroche_undo_tile_)` — already non-null) and
OVERWRITES its contents with the inner's 13x13 backing snapshot. The
outer's saved tile is gone; when the outer resumes, its `*pmin0_` still
points to a drawn sprite on the canvas but `restore_tile` would blit the
wrong content. The inner's `cleanupAccrochage()` on release/abort also
deletes `accroche_undo_tile_` outright (line 248-249), leaving the outer's
drawn sprite with no undo path.
**Fix:** Either (a) push `accroche_undo_tile_` ownership down to the
bridge (stack-local to waitForEvent), or (b) save-restore it through the
7-member save/restore list alongside `items_` and `pmin0_`:
```cpp
// In waitForEvent, before loop.exec():
QPixmap* saved_tile = state ? state->accroche_undo_tile_ : nullptr;
int      saved_pmin = pmin0 ? *pmin0 : -2;
if (state) state->accroche_undo_tile_ = nullptr; // inner starts fresh
// ... loop.exec() ...
// On return, restore (leak-safe: inner cleanup already freed its own tile)
if (state) state->accroche_undo_tile_ = saved_tile;
if (pmin0) *pmin0 = saved_pmin;
```
This keeps the re-entrant contract the bridge already promises for the
other 7 scalars.

### WR-03: BlockingDepthGuard Q_ASSERT ceiling compiles out in release

**File:** `xvue/qt/src/xvue_qt_event.h:27`
**Issue:** The depth-DoS guard is `Q_ASSERT(XvueApp::blockingDepth_ < 4)`,
which is `#ifdef QT_DEBUG` — in a release/production build (which is the
default for `bin/cbl_tout_qt`) the assertion compiles to nothing. A
runaway Fortran caller that re-enters the bridge without returning can
drive `blockingDepth_` to arbitrary values, trivially defeating
T-05-02-01 depth-DoS mitigation in production. The `blockingDepth_`
counter itself is unbounded (plain `int`).
**Fix:** Replace the Q_ASSERT with a hard cap + stderr diagnostic that
compiles in all build modes, or increment with an explicit max-guard:
```cpp
BlockingDepthGuard() noexcept {
    if (XvueApp::blockingDepth_ >= 4) {
        std::fprintf(stderr,
            "xvue-qt: BlockingDepthGuard refusing entry at depth %d; "
            "nested waitForEvent limit exceeded (T-05-02-01)\n",
            XvueApp::blockingDepth_);
        std::abort();  // fail loud — a bug, not a workload
    }
    ++XvueApp::blockingDepth_;
}
```
Alternatively, if abort is too strong, keep the counter but refuse to
start a new loop via a status flag that `waitForEvent` honors to silently
return with `pending_ = Result{}`.

### WR-04: xvtest-capture.sh cleanup() traps lose exit code on SIGINT/SIGTERM

**File:** `bin/xvtest-capture.sh:61-74`
**Issue:** `cleanup()` captures `$?` into `rc` at entry, then `exit $rc`.
On a TERM or INT signal the running command's exit code is captured (often
0 if import was mid-run), not the signal-derived 130/143 a POSIX shell
would normally propagate. Capture harnesses running under a CI timeout
wrapper will see 0 and falsely report success when the wrapper killed
them. `set -u` is also active; the `DRIVER_PID=""` initialization on line
73 occurs AFTER `trap cleanup EXIT INT TERM` is registered on line 74 —
if an exit happens in the 1-instruction window between trap install and
DRIVER_PID init, `cleanup` dereferences an unset variable and itself
errors under `set -u`.
**Fix:** Swap the order so DRIVER_PID is initialized BEFORE the trap is
installed, and propagate signals:
```bash
DRIVER_PID=""
cleanup() {
  local rc=$?
  if [ -n "${DRIVER_PID:-}" ] && kill -0 "$DRIVER_PID" 2>/dev/null; then
    kill "$DRIVER_PID" 2>/dev/null; wait "$DRIVER_PID" 2>/dev/null
  fi
  if [ -n "${XVFB_PID:-}" ] && kill -0 "$XVFB_PID" 2>/dev/null; then
    kill "$XVFB_PID" 2>/dev/null; wait "$XVFB_PID" 2>/dev/null
  fi
  exit "$rc"
}
trap cleanup EXIT
trap 'rc=$?; cleanup; exit $((128+rc))' INT TERM
```

## Info

### IN-01: Souris2 middle-button comment contradicts the code (xvsouris2_)

**File:** `xvue/qt/src/xvue_qt_event.cpp:273-282`
**Issue:** Comment says "saclav.f:83-86 remaps btn=2 to notypeevent=-1 /
nbc=64 (@ abort path)" — implying the Qt code at line 278 should set
`pending_.nbc = 64`. The code actually sets `pending_.nbc = 2`, which is
correct IF the code path is preserved (saclav.f does the remap itself).
With the WR-01 fix applied this whole block disappears; otherwise the
comment should be rewritten to explain that the Qt layer returns the raw
button number and saclav.f handles the remap.
**Fix:** Either remove the block per WR-01, or update the comment to:
`// Middle button in Souris2: X11 parity would fall through to the
// accrochage block; we short-circuit to ntev=0 here for X. Keep nbc=2
// so saclav.f:83-86 can do the @-abort remap on the next release.`

### IN-02: test_helpers.cpp createTestCanvas/destroyTestCanvas are dead stubs

**File:** `xvue/qt/tests/test_helpers.cpp:20-29`
**Issue:** `createTestCanvas()` always returns nullptr and
`destroyTestCanvas()` is empty. The phase 5 test file
`test_xvue_qt_event.cpp` never calls either — it uses
`XvueApp::window_slot()->canvas()` directly. These helpers are unused
dead code declared in the header for Plan 02+ work that was then
redirected.
**Fix:** Either delete the two helpers and their header declarations, or
leave a `[[deprecated]]` tag with a comment pointing at the replacement
pattern. Low priority — does not affect behavior.

### IN-03: `verify_abi.sh` `grep -c` exit-code handling relies on non-pipefail

**File:** `xvue/qt/cmake/verify_abi.sh:22-23`
**Issue:** Under `set -eu`, `nm "$LIB" | grep ' T ...' | grep -vc ' T _Z'`
could theoretically abort the shell on a match-count of 0 if `pipefail`
were ever enabled. The `|| true` catches only the last stage. If someone
adds `set -o pipefail` for debugging and nm happens to produce no matches
(catastrophic build but possible), the script fails in a confusing place
before printing the drift diagnostic.
**Fix:** Wrap each grep with `|| true`:
```sh
NM_COUNT=$( { nm "$LIB" | grep ' T [a-zA-Z_][a-zA-Z0-9_]*_$' || true; } \
           | { grep -vc ' T _Z' || true; } )
```
Or compute count via `wc -l` at the end of the pipeline, which always
exits 0.

### IN-04: `cbl_tout_qt` unquoted `rm $MEFISTO/*/*.o` under set-less shell

**File:** `bin/cbl_tout_qt:103`
**Issue:** `rm $MEFISTO/*/*.o` runs under `#!/bin/bash` WITHOUT `set -e` or
`set -u`. If `$MEFISTO` is unset or empty, `rm //*.o` expands to `rm
/foo.o /bar.o ...` at root. If `$MEFISTO` contains whitespace, word-
splitting breaks the glob. The script earlier prints `MEFISTO=$MEFISTO`
but never validates that the variable is non-empty.
**Fix:** Add a guard at the top of the script:
```bash
: "${MEFISTO:?MEFISTO must be set to the source tree root}"
```
and quote the glob expansion in a way that still allows globbing via `set
--`, or simply drop the broad `*/*.o` cleanup in favour of per-module
scripts' own clean steps.

### IN-05: `initaccrochage_` does NOT free a prior `accroche_undo_tile_`

**File:** `xvue/qt/src/xvue_qt_api.cpp:298-335`
**Issue:** `initaccrochage_` correctly deletes and reallocates
`state->mempxaccro_` on a repeat call (lines 310-313), but leaves
`state->accroche_undo_tile_` untouched. If a previous xvsouris2_ session
terminated abnormally (e.g., the window was destroyed mid-drag before
`cleanupAccrochage()` ran, though the current code paths do call it),
`accroche_undo_tile_` could still hold a 13x13 tile from the prior run.
The resize-invalidation in xvue_qt_canvas.cpp:134-137 and the destructor
in xvue_qt_state.cpp:158-159 handle most termination paths; this is
belt-and-braces.
**Fix:** Free `accroche_undo_tile_` at the top of initaccrochage_ for
symmetry with mempxaccro_:
```cpp
if (state->accroche_undo_tile_) {
    delete state->accroche_undo_tile_;
    state->accroche_undo_tile_ = nullptr;
}
```

### IN-06: xvsouris_ no-window path returns `notypeevent = 0` without Esc semantics

**File:** `xvue/qt/src/xvue_qt_api.cpp:982-989`
**Issue:** When no window is open, xvsouris_ silently returns
`notypeevent = 0, nbc = 0, x = 0, y = 0`. saclav.f and similar callers
branch on `NOCODE .EQ. 0` → "PROBLEME AVEC LA SOURIS" → call XVPAUSE →
retry. Since there's also no window for XVPAUSE, the Fortran caller
loops forever. Consider returning `notypeevent = 0, nbc = 27` (Esc abort)
instead, so callers exit deterministically on the abort path.
**Fix:**
```cpp
if (!win || !win->bridge()) {
    if (notypeevent) *notypeevent = 0;
    if (nbc)         *nbc         = 27;  // Esc — callers treat as abort
    if (x1)          *x1          = 0;
    if (y1)          *y1          = 0;
    return;
}
```
The same pattern applies to `xvsouris2_` no-window guard (line 1038-1045)
and `xvpause_` (line 1123 — but xvpause has no out-params, just
`return;`, which is fine). Low priority because in practice no driver
calls these without a window open.

---

_Reviewed: 2026-04-18T00:00:00Z_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
