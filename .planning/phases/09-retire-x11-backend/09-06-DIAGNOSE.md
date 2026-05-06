# 09-06 — Matched-Dim Qt Recapture: Diagnose Log

Working scratch document for Plan 9-06 Task 1.

## Insertion points (xvue/qt/src/xvue_qt_api.cpp)

- `xvinitgraphique_` env-probe block: lines **350–381** (block lives AFTER
  `win->show()/raise()/activateWindow()` and BEFORE the bounded
  exposure-pump loop). Reads `MEFISTO_QT_WINDOW_SIZE`, gates on
  `MEFISTO_BATCH_X11=1` OR `QT_QPA_PLATFORM=offscreen`, calls
  `setMinimumSize/setMaximumSize` on the canvas + `win->adjustSize()`.
- `xvfermer_` clear block: lines **943–953** (just before
  `XvueApp::window_slot().reset()`). Resets canvas constraints to
  `(0,0)` / `(QWIDGETSIZE_MAX, QWIDGETSIZE_MAX)` so a subsequent reopen
  starts clean.

Insertion strategy is per RESEARCH §Pattern 4 + §Pitfall 5:
- env-probe in `xvinitgraphique_` (the canonical "open window" entry —
  fires once per reopen cycle), NOT in `xvinfo_` (where the original
  `win->resize(*ix, *iy)` call lives). Constraints applied here clip
  any subsequent `win->resize` from `xvinfo_` to the env-specified dim.
- `setMinimumSize` + `setMaximumSize` on the **canvas** (NOT
  `setFixedSize` on the **window**) so `xvfermer_` can zero them out and
  the X11-style chrome computation does not interfere.

## Empirical PNG dim (Task 1 step 6)

Command (from worktree root, `MEFISTO=$(pwd)`):

```sh
MEFISTO_QT_WINDOW_SIZE=1280x800 \
  bin/ab_sweep_phase8.sh --mode qt-1x --cases pan2d \
      --out-dir /tmp/09-06-verify-pan2d
identify -format '%wx%h' /tmp/09-06-verify-pan2d/pan2d-qt-1x.png
```

Output: **`1280x800`** (exact match — no chrome-accounting tolerance
needed).

Iteration count: **1** (first attempt produced the matched dim).
Open Question 5 (RESEARCH) is empirically resolved: `setMinimumSize` +
`setMaximumSize` on the canvas DOES propagate to
`XvueCanvas::resizeEvent`'s backing-pixmap allocation once
`win->adjustSize()` runs and `win->show()` triggers the natural
exposure-driven resize chain. No additional `processEvents` pump or
`QTimer::singleShot(0)` deferral was needed.

## ABI count (T-09-06 mitigation)

Pre-edit count: **64** (per 09-01-AUDIT-BASELINE §3 row 7 — the
empirical baseline; the *plan-quoted* "58" is the Phase-6.5-frozen value
that drifted during the Phase 7 export plans, +6 entries — see
09-01-AUDIT-BASELINE §"#7 ABI count drift").

Post-edit count: **64** (`nm xvue/qt/build/libxvueqt.a | grep ' T ' |
grep '_$' | wc -l = 64`). Delta: 0 — no new `extern "C"` entry was
added. The implementation is env-var-only per RESEARCH §Anti-Pattern
"Bumping ABI count from 58 to 59": adding a Fortran-callable setter
would have been the wrong fix.

The audit baseline's #7 row notes "Phase 9 retirement does not touch Qt
extern \"C\" entry points" — Plan 9-06 is consistent with that
expectation.

## Build status

`bin/cbl_tout` exit: **0**. All `pp/pp{init,mail,elas,flui,ther,nlse,
xvtest0..4}` rebuilt cleanly with the matched-dim code path. No new
compiler warnings introduced (`std::strcmp` and `std::sscanf` were
already in the file's reach via `<cstring>` and `<cstdio>` includes).

## Pitfall 5 mitigation evidence

The env probe is wrapped in
```cpp
if ((headless_batch || offscreen) && qt_window_size && qt_window_size[0]) {
```
Interactive sessions where neither `MEFISTO_BATCH_X11=1` nor
`QT_QPA_PLATFORM=offscreen` is set will short-circuit at the first
`&&` test even if `MEFISTO_QT_WINDOW_SIZE` is somehow set — the
constraints never fire and the existing `XvueWindow::resize(1024, 768)`
default from the constructor stays in effect.

T-09-05-A (input validation): the `sscanf("%dx%d")` failure path returns
`!=2` and the bounds check `0 < {w,h} < 8192` rejects negatives,
overflow-prone values, and zero. Malformed env values are silently
ignored — no error path is exposed to the headless harness.
