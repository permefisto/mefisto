---
phase: 02-drawing-primitives-backing-pixmap
fixed_at: 2026-04-11T21:40:00Z
review_path: .planning/phases/02-drawing-primitives-backing-pixmap/02-REVIEW.md
iteration: 1
findings_in_scope: 6
fixed: 6
skipped: 0
status: all_fixed
---

# Phase 2: Code Review Fix Report

**Fixed at:** 2026-04-11
**Source review:** `.planning/phases/02-drawing-primitives-backing-pixmap/02-REVIEW.md`
**Iteration:** 1
**Scope:** critical_warning (CR-01 + WR-01..05) — Info findings deferred.

**Summary:**
- Findings in scope: 6
- Fixed: 6
- Skipped: 0
- Qt build (`bin/cbl_tout_qt`): green after every fix
- Legacy X11 build (`bin/cbl_tout`): green after final fix
- `pp/ppxvtest0_qt` smoke (offscreen): exit 0, zero `xvue-qt: stub` warn-once lines, clean reopen cycle
- ABI symbol count (libxvueqt.a trailing-underscore `T` symbols): 57 (unchanged)

## Fixed Issues

### CR-01: Rectangle ABI collapses fill and outline semantics

**Files modified:** `xvue/qt/src/xvue_qt_api.cpp`
**Commit:** 6741c68
**Applied fix:** Added a `RectMode { Outline, Fill }` enum and refactored
`xvue_qt_draw_rect_common` to branch on the mode. `Outline` saves the current
brush, switches to `Qt::NoBrush`, calls `drawRect`, then restores the brush
(matches `XDrawRectangle`). `Fill` uses `fillRect(r, painter->brush())`
(matches `XFillRectangle`). The four legacy entry points now dispatch:

- `xvfbordrectangle_` → `RectMode::Outline`
- `xvbordrectangle_`  → `RectMode::Outline`
- `xvfrectangle_`     → `RectMode::Fill`
- `xvrectangle_`      → `RectMode::Fill`

This bundle also carries the WR-03 null-guards and WR-02 processEvents
removal for the rectangle helper, since they sit in the same code path
and must land together to keep the build green.

**Note:** Logic-correctness verification is by inspection + `xvtest0.f`
exercising all four symbols at lines 58-61. Human visual verification
recommended once Phase 3 differentiates pen/brush colors, per review
suggestion — flag as "fixed: requires human verification under Phase 3
palette".

### WR-01: xvface_ and xvfacetraits_ stroke extra outline vs legacy XFillPolygon

**Files modified:** `xvue/qt/src/xvue_qt_api.cpp`
**Commit:** fe52db4
**Applied fix:** `xvface_` now wraps `drawPolygon(poly, Qt::OddEvenFill)` in
a `setPen(Qt::NoPen)` / restore-pen pair, so only the brush paints the
interior (matches `XFillPolygon`). `xvfacetraits_` splits fill and outline
explicitly: fill with `Qt::NoPen` + current brush, then restore the pen and
draw the outline with `Qt::NoBrush` + current pen, then restore the brush.
This prevents double-painting of edge pixels and lays the groundwork for
Phase 3 ncf/nca color handling. Per-primitive `processEvents` also removed
from both (WR-02 overlap).

### WR-02: Per-primitive processEvents enables re-entrant resize mid-batch

**Files modified:** `xvue/qt/src/xvue_qt_api.cpp`
**Commit:** b9b37d6 (plus earlier bundles in 6741c68 and fe52db4)
**Applied fix:** Removed `QCoreApplication::processEvents(ExcludeUserInputEvents)`
from every draw primitive that carried it: `effacer_`, `xvfond_`, `xvtrait_`,
`xvtraits_`, `xvbordarcellipse_`, `xvarcellipse_`, plus the rectangle and
polygon entry points already cleaned up in CR-01 / WR-01. Flush points are
now `xvvoir_` (explicit), `xvpause_`, and the teardown drain in `xvfermer_`.
Primitives still call `canvas()->update()` which coalesces into a single
paint event. Each removed call site carries a `// WR-02: deferred flush`
comment for future editors.

### WR-03: Primitive entry points dereference int*/float* without null checks

**Files modified:** `xvue/qt/src/xvue_qt_api.cpp`
**Commit:** b9b37d6 (plus rectangles in 6741c68)
**Applied fix:** Added symmetric `if (!arg1 || !arg2 || ...) return;` null
guards to `xvtypetrait_`, `xvepaisseur_`, `xvtrait_`, `xvbordarcellipse_`,
`xvarcellipse_`, and the four `xv*rectangle_` entry points. Pattern matches
the pre-existing guard in `xvpxfenetre_` (line 399). `xvface_`, `xvtraits_`,
and `xvfacetraits_` already had full null-checks.

### WR-04: xvue_qt_draw_rect_common can't self-protect against null args

**Files modified:** `xvue/qt/src/xvue_qt_api.cpp`
**Commit:** 6741c68
**Applied fix:** Kept the helper's `int`-by-value signature and centralized
the null protection at every caller (xv*rectangle_). Review explicitly
offered "either works, pick one and apply consistently" — caller-side
guards are consistent with the rest of the file. No helper signature
change needed; this finding is fully addressed by the WR-03 rectangle
guards landed in the CR-01 commit.

### WR-05: resizeEvent leaks painter_ if begin() fails on new backing

**Files modified:** `xvue/qt/src/xvue_qt_canvas.cpp`
**Commit:** ccddb1d
**Applied fix:** Reordered `resizeEvent` so every throwing allocation
(`new QPixmap`, `new QPainter`) happens BEFORE the old painter/backing are
torn down. If `new QPixmap(device)` throws, the function returns early
with the old painter+backing still intact and logs a stderr diagnostic.
Added explicit check of `QPainter::begin()` return value with a loud
stderr diagnostic so a frozen canvas is no longer silent. Wrapped the
QPainter `new` in try/catch to avoid leaving the state with a dangling
inactive painter.

## Skipped Issues

None — all in-scope findings were fixed.

## Verification Matrix

| Check | Result |
|---|---|
| `bin/cbl_tout_qt` (Qt build) | green after every fix |
| `bin/cbl_tout` (legacy X11 regression) | green after final fix |
| `bin/cbxvtest0_qt` (rebuild smoke driver) | green |
| `pp/ppxvtest0_qt` exit code (offscreen platform) | 0 |
| `pp/ppxvtest0_qt` `xvue-qt: stub` warn-once lines | 0 |
| `pp/ppxvtest0_qt` reopen cycle | clean (two open/close passes) |
| `libxvueqt.a` trailing-underscore `T` symbol count | 57 (unchanged) |

## Info Findings Deferred

Per scope (`critical_warning`), IN-01..04 are not addressed in this fix
pass. IN-01 (xvfacetraits_ ncf/nca) pairs naturally with Phase 3 palette
work. IN-02..04 are optional polish.

---

_Fixed: 2026-04-11_
_Fixer: Claude (gsd-code-fixer)_
_Iteration: 1_
