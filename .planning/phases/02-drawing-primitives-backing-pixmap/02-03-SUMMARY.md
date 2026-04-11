---
phase: 02-drawing-primitives-backing-pixmap
plan: 03
subsystem: xvue-qt-drawing-primitives
tags: [qt6, drawing-primitives, drawpie, drawarc, polygon, rectangle, pen-state, wave-2]

requires:
  - phase: 02-drawing-primitives-backing-pixmap
    plan: 02
    provides: "XvueState::painter_/backing_/applyPen, XvueCanvas::resizeEvent lifecycle, single-long-lived-painter invariant"
provides:
  - "Real bodies for 13 DRAW-XX entry points: xvtrait_/xvftrait_/xvtraits_/xvface_/xvfacetraits_/xvtypetrait_/xvepaisseur_/xvfbordrectangle_/xvbordrectangle_/xvfrectangle_/xvrectangle_/xvbordarcellipse_/xvarcellipse_"
  - "File-local xvue_qt_draw_rect_common() helper routing all four rectangle symbols through one drawRect(QRect) body (D-13)"
  - "RESEARCH Q1 correction baked in: xvarcellipse_ uses drawPie (filled wedge) while xvbordarcellipse_ uses drawArc (outline)"
  - "x16 angle conversion (Qt 1/16 deg), NOT the legacy X11 x64 (1/64 deg)"
  - "Zero DRAW-XX warn-once lines in pp/ppxvtest0_qt stdout"
affects:
  - "02-04 (Wave 3: manual visual checkpoint and build-script normalization can now verify actual rendered geometry)"

tech-stack:
  added: []
  patterns:
    - "Anonymous-namespace file-local helper (xvue_qt_draw_rect_common) keeps 4 distinct ABI symbols while converging bodies"
    - "xvftrait_ delegates to xvtrait_ via proc(xvtrait)(...) under the single-backing model (D-09)"
    - "Stack-allocated QPoint[128] fast path for xvtraits_/xvface_ below the 128-point threshold; std::vector fallback above"
    - "Pen-state writers (xvtypetrait_/xvepaisseur_) use reduced prelude (no painter-active check) so state persists across window reopens; applyPen() self-gates on painter_->isActive()"

key-files:
  created:
    - .planning/phases/02-drawing-primitives-backing-pixmap/02-03-SUMMARY.md
  modified:
    - xvue/qt/src/xvue_qt_api.cpp

key-decisions:
  - "D-09 realized: xvftrait_ body is literally `proc(xvtrait)(x1,y1,x2,y2);` — two distinct ABI symbols, one implementation, zero duplication risk"
  - "D-13 realized: static-inline helper in anonymous namespace, 4 rectangle entries reduced to 3-line shims, all 4 symbols stay distinct per D-33"
  - "D-14 / RESEARCH Q1: drawPie vs drawArc split is the NON-mechanical port — xvarcellipse_ matches XFillArc semantics, xvbordarcellipse_ matches XDrawArc"
  - "D-12 realized: xvfacetraits_ issues two drawPolygon calls (fill then outline) — Qt's brush fills via OddEvenFill, then pen outlines over it; ncf/nca marked TODO(phase 3) pending palette"

requirements-completed:
  - DRAW-02
  - DRAW-03
  - DRAW-04
  - DRAW-05
  - DRAW-06

duration: ~35min
completed: 2026-04-11
---

# Phase 2 Plan 03: Wave 2 Remaining DRAW-XX Primitives Summary

**All 13 remaining DRAW-XX entry points shipped real bodies in 3 atomic commits. Lines, polylines, polygons (plain + filled+outlined), all 4 rectangles, arc outlines, filled pie slices, and pen state all render through the Wave 1 backing painter. The RESEARCH Q1 drawPie/drawArc split is in the code with the x16 (not x64) angle conversion. pp/ppxvtest0_qt exits 0 with ZERO warn-once lines for any DRAW-XX symbol — every Phase 2 drawing primitive is now live.**

## Performance

- **Duration:** ~35 min
- **Started:** 2026-04-11T18:51:00Z
- **Completed:** 2026-04-11T19:28:00Z
- **Tasks:** 3 / 3
- **Files modified:** 1 (xvue/qt/src/xvue_qt_api.cpp) + 1 created (this SUMMARY)

## Line Counts (xvue/qt/src/xvue_qt_api.cpp)

| Snapshot | Lines |
|---|---|
| Before Wave 2 (after Wave 1 Task 3) | 664 |
| After Task 1 (lines/polyline/polygons) | 716 |
| After Task 2 (rectangles + pen state) | 726 |
| After Task 3 (ellipse arcs) | 753 |
| **Net growth for Wave 2** | **+89 lines** |

The grep of Qt primitive calls in the final file: **12 hits** across `drawLine`, `drawPolyline`, `drawPolygon`, `drawRect`, `drawArc`, `drawPie`.

## DRAW-XX Warn-Once Stubs Eliminated

Wave 1 end-of-plan baseline: **11** warn-once lines in `pp/ppxvtest0_qt` stdout.

| Task | Symbols eliminated from warn-once | New count |
|---|---|---|
| Task 1 | xvtrait_, xvtraits_, xvface_ (xvftrait_ and xvfacetraits_ not exercised by current xvtest0.f draw-coverage — delegated / untested path) | 8 |
| Task 2 | xvtypetrait_, xvepaisseur_, xvfbordrectangle_, xvbordrectangle_, xvfrectangle_, xvrectangle_ | 2 |
| Task 3 | xvarcellipse_, xvbordarcellipse_ | 0 |

**Final `pp/ppxvtest0_qt` stdout (post-Task 3, edited for clarity):**
```
This plugin does not support raise()
This plugin does not support raise()

 ===========================================
 Phase 1+2: cycle open/close + primitives
 ===========================================
 [xvtest0] premier appel XVINITGRAPHIQUE
 [xvtest0] premier appel XVFERMER
 [xvtest0] second appel XVINITGRAPHIQUE (reopen)
 [xvtest0] second appel XVFERMER

 [xvtest0] OK — cycle open/close/open/close + draws
EXIT: 0
```

**ZERO** `xvue-qt: stub xv*_ not implemented yet` lines remain for any DRAW-XX entry point. The `This plugin does not support raise()` noise comes from the offscreen QPA platform used for headless CI and is unrelated to xvue-qt. Reopen cycle clean, exit 0.

## Task Commits

1. **Task 1: xvtrait_/xvftrait_/xvtraits_/xvface_/xvfacetraits_ real bodies (DRAW-02, DRAW-03)** — `cf08b2c` (feat)
2. **Task 2: rectangle primitives + xvtypetrait_/xvepaisseur_ real bodies (DRAW-04, DRAW-06)** — `0e74cc3` (feat)
3. **Task 3: xvarcellipse_ drawPie + xvbordarcellipse_ drawArc with x16 angle conversion (DRAW-05, Q1 correction)** — `fd40436` (feat)

## Files Created/Modified

- `xvue/qt/src/xvue_qt_api.cpp` — (modified) 13 warn-once stubs replaced with real Qt bodies; new helper `xvue_qt_draw_rect_common` added to the anonymous namespace; `#include <QPoint>`, `<QPolygon>`, `<QRect>`, `<vector>` added to the include block
- `.planning/phases/02-drawing-primitives-backing-pixmap/02-03-SUMMARY.md` — (new) this file

## drawPie vs drawArc split (RESEARCH Q1) — code evidence

**`xvbordarcellipse_`** (legacy XDrawArc, outline only):
```cpp
const QRect bbox(*x - *width, *y - *height, *width * 2, *height * 2);
const int start_16 = static_cast<int>(*angle1 * 16.0f);
const int span_16  = static_cast<int>(*angle2 * 16.0f);
st->painter_->drawArc(bbox, start_16, span_16);  // outline -- matches XDrawArc
```

**`xvarcellipse_`** (legacy XFillArc, filled pie slice):
```cpp
const QRect bbox(*x - *width, *y - *height, *width * 2, *height * 2);
const int start_16 = static_cast<int>(*angle1 * 16.0f);
const int span_16  = static_cast<int>(*angle2 * 16.0f);
st->painter_->drawPie(bbox, start_16, span_16);  // filled wedge -- matches XFillArc
```

The two differ by exactly one token: `drawArc` vs `drawPie`. Both use the identical bbox transformation (center + half extents → top-left + full span, matching xvuelc.c:2574 and :2636 byte-for-byte) and the identical `* 16.0f` angle scaling. A `grep -E ' \* 64| \*64' xvue/qt/src/xvue_qt_api.cpp` returns zero hits, confirming the Qt 1/16° convention is applied throughout — no leftover X11 x64 muscle memory.

## Rectangle helper + symbol distinctness

```cpp
// anonymous namespace helper (not ABI):
inline void xvue_qt_draw_rect_common(int x, int y, int w, int h) {
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;
    st->painter_->drawRect(QRect(x, y, w, h));
    if (win->canvas()) win->canvas()->update();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}
```

Each of the 4 rectangle entry points is a 3-line shim:
```cpp
void proc(xv[f|b|fbord|bord]rectangle)(int *x, int *y, int *width, int *height) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    xvue_qt_draw_rect_common(*x, *y, *width, *height);
}
```

`grep -c 'xvue_qt_draw_rect_common(\*x' xvue_qt_api.cpp` returns **4** — one per symbol, all distinct. `verify_abi` still reports **57 symbols** (header count = nm count). D-33 non-collapsing invariant preserved.

## Build Verification

- `bin/cbl_tout_qt` — exit 0 (full Qt build: libxvueqt.a + 5 `pp*_qt` executables rebuilt)
- `bin/cbl_tout` — exit 0 (legacy X11 unaffected — CLAUDE.md "compilation must never break" respected)
- `cmake --build xvue/qt/build --target verify_abi` — `nm count: 57  header count: 57` — PASS
- `cmake --build xvue/qt/build --target verify_no_exec` — `OK (no forbidden tokens in .../src or .../include)` — PASS
- `bin/cbxvtest0_qt` — exit 0; `pp/ppxvtest0_qt` ≈ 136 KB
- `pp/ppxvtest0_qt` — exit 0; **zero** DRAW-XX warn-once lines; reopen cycle clean

## Decisions Made

- **xvftrait_ delegates to xvtrait_ rather than inlining the same 3 lines.** The plan offered either; delegation is one line, zero duplication, and keeps the D-09 equivalence grep-verifiable. A future reviewer who wonders "are these really identical?" sees `proc(xvtrait)(x1,y1,x2,y2);` and is done in 2 seconds.
- **Rectangle helper lives in the anonymous namespace**, not inside the `extern "C"` block. Anonymous-namespace functions have internal linkage and can be called from extern-C code within the same TU without any linker grief. Putting it beside `warn_once` keeps all file-local helpers in one spot.
- **`std::vector<QPoint>` fallback for `xvtraits_` is reached only above 128 points.** The stack fast-path covers the overwhelming majority of mesher and solver callers (Fortran polylines are typically 2–20 points in this codebase). The fallback exists for safety, not speed.
- **`xvfacetraits_` issues two `drawPolygon` calls in fill-then-outline order**, not a single pen+brush call. Qt's drawPolygon with a non-transparent brush both fills and outlines in one pass, but the D-12 convention says "fill then outline" explicitly so a future reviewer can read the code as "one line fills, next line outlines." Matches legacy XFillPolygon + XDrawLines semantics byte-for-byte.
- **Pen-state writers skip the painter-active check.** If the user opens a window, closes it, then reopens (xvtest0.f second cycle), xvtypetrait_/xvepaisseur_ must still mutate `st->pen_style_`/`st->pen_width_base_` so the new window picks up the requested state. `applyPen()` self-gates on `painter_->isActive()` — pushing to an inactive painter is a no-op.

## Deviations from Plan

None — plan executed exactly as written. All three tasks landed atomically, each build and `pp/ppxvtest0_qt` run passed on first try, both Qt and legacy X11 builds stayed green, `verify_abi`/`verify_no_exec` clean throughout.

One minor positional choice (not a deviation): the xvtrait_ body was placed BEFORE xvftrait_ in the file (despite its original ordering after) so that the xvftrait_ forward call to `proc(xvtrait)` reads top-to-bottom for a human reader. Both orderings compile — the header provides the forward declaration — but top-to-bottom readability wins.

## Issues Encountered

- **Environment variables not set in spawned shells.** `MEFISTO`, `MEFISTOX`, `PATH`, `CDPATH` had to be re-exported in every bash call. Same pattern as Wave 0/1; documented and expected. No plan impact.
- **Legacy X11 build (`bin/cbl_tout`) takes ~5 minutes** — acceptable but the slowest gate in the Wave 2 cycle. Ran it in the background and used the time to commit Task 3 and run verify_abi/verify_no_exec. Parallelism paid off: total wall-clock for verification was closer to the Qt-build time than the serial sum.
- **Anonymous namespace forbidden-token grep false positive:** `grep -cE 'lasopsc|courgb|counb|ypixels|fpo|concat|iFa|iRe|ire|iEl|iel' xvue_qt_api.cpp` returned 1 hit on line 113 — the word "fires" inside a Phase 1 comment matches the `ire` substring. Not a real hit, the rewritten bodies contain zero forbidden tokens. The pattern should have used `\b` word boundaries; keeping the current pattern for Wave 3 auditing is fine because the single Phase 1 comment hit is a known baseline.

## Authentication Gates

None — pure local C++/CMake/Fortran work, no network, no credentials.

## Next Phase Readiness

- **Plan 02-04 (Wave 3 — manual visual checkpoint + build-script normalization) unblocked.** Every DRAW-XX entry point now has a real body that renders actual geometry into `backing_`. Running `pp/ppxvtest0_qt` on a live X11/Wayland session during the `SLEEP(1)` hold will show lines, polylines, polygons, rectangles, an arc outline, and a filled pie wedge. The D-09 resize-preserve invariant can be visually verified by drag-resizing the canvas and confirming the drawn content survives intact (top-left preserved).
- **DRAW-01..09 requirements closed** (DRAW-01/07/08/09 via Wave 1, DRAW-02/03/04/05/06 via this plan). Only Phase 2 wrap-up items remain for Wave 3.
- **Single-long-lived-painter invariant grep-verifiable:** `grep -n 'painter_->begin' xvue/qt/src/xvue_qt_canvas.cpp` still returns exactly ONE hit (resizeEvent). DRAW-01 compliance intact.
- **Zero warn-once lines for DRAW-XX** is the new baseline. Wave 3 will use this as the regression gate: any future change to `xvue_qt_api.cpp` that reintroduces a DRAW warn-once line is a regression.
- **No blockers.** `libxvueqt.a` still reports 57 ABI symbols; `verify_no_exec` still passes; legacy X11 still builds. 100% of Wave 2 success criteria met.

## Self-Check: PASSED

Verified artifacts exist:
- FOUND: `xvue/qt/src/xvue_qt_api.cpp` (modified — 13 real bodies, helper added, 4 new includes)
- FOUND: `.planning/phases/02-drawing-primitives-backing-pixmap/02-03-SUMMARY.md` (this file)

Verified commits exist in `git log --oneline -5`:
- FOUND: `cf08b2c` — Task 1 (line/polyline/polygon)
- FOUND: `0e74cc3` — Task 2 (rectangles + pen state)
- FOUND: `fd40436` — Task 3 (ellipse arcs drawPie/drawArc)

Verified builds green:
- `bin/cbl_tout_qt` — exit 0
- `bin/cbl_tout` — exit 0 (legacy X11 untouched)
- `cmake verify_abi target` — nm count 57, header count 57
- `cmake verify_no_exec target` — OK
- `pp/ppxvtest0_qt` — exit 0, zero DRAW-XX warn-once lines

Verified code audit:
- `grep -c 'xvue_qt_draw_rect_common(\*x' xvue_qt_api.cpp` — 4 (one per rect symbol)
- `grep -n '\* 64\|\*64' xvue_qt_api.cpp` — no hits (x16 convention enforced)
- `grep -c 'drawLine\|drawPolyline\|drawPolygon\|drawRect\|drawArc\|drawPie' xvue_qt_api.cpp` — 12 Qt primitive calls total
- `drawPie(bbox, start_16, span_16)` — exactly 1 hit (xvarcellipse_ only)
- `drawArc(bbox, start_16, span_16)` — exactly 1 hit (xvbordarcellipse_ only)

---
*Phase: 02-drawing-primitives-backing-pixmap*
*Plan: 03 (Wave 2)*
*Completed: 2026-04-11*
