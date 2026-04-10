# xvue — Y-Axis Coordinate Convention (audited 2026-04-10)

> This document records the audited Y-axis convention used by `xvue/xvuelc.c`
> and mirrored byte-for-byte by `xvue/qt/` (the Qt 6 backend). Every phase of
> the xvue-qt migration from Phase 1 onward MUST read this document before
> writing any QPainter-based drawing code.

## On-screen coordinates (Fortran ↔ C bridge ↔ graphics backend)

| Property | Value |
|----------|-------|
| Origin | Top-left corner of the canvas widget |
| X axis | Grows rightward |
| Y axis | Grows **downward** (Y-down) |
| Unit (X11) | Device pixels |
| Unit (Qt) | Logical pixels (device-independent, HiDPI-aware) |

## Where the convention comes from

Direct evidence from `xvue/xvuelc.c` (read 2026-04-10):

- Line 321 (comment inside `xvpxecran_`): "origine coin superieur gauche de l'ecran"
  — "origin = upper-left corner of the screen"
- Line 1852 (comment inside line-drawing block): "ORIGINE coin superieur gauche"
- Line 1869 (comment, same context): "ORIGINE coin superieur gauche"
- Line 1621 (comment inside `xvpxfenetre_`): "la fenetre" (window) — origin implied top-left

Every `XDrawLine`, `XDrawLines`, `XFillRectangle`, and `XDrawArc` call in
`xvue/xvuelc.c` passes the Y coordinate **unflipped**. This matches:

- **Xlib's native convention** — Y-down top-left origin, device pixels
- **Qt's native `QPainter` convention** — Y-down top-left origin, logical pixels

Therefore the Qt backend **does not transform Y coordinates** anywhere on the
on-screen path. Fortran-provided `(x, y)` values are passed directly to
`QPainter::drawLine(*x, *y, ...)` etc. in Phase 2+.

## Where Y IS flipped (PostScript export only — Phase 7 concern)

PostScript's natural origin is **bottom-left** (Y-up). The hand-rolled PostScript
emitter in `xvuelc.c` honors this by flipping Y at emit time, not by storing
flipped coordinates in any state.

Direct evidence:

- `xvue/xvuelc.c` lines 1895, 1932, 1953, 1966 — all inside the `xvtrait_`
  PostScript-capture branch — write `ypixels - *y1` (and analogous for the
  second endpoint) into the emitted `.ps` stream, where `ypixels` is the
  canvas height in pixels.
- Same pattern throughout the `xvpostscript_` function body (line 1187+).

Formula: `ps_y = ypixels - screen_y`.

## Implications for the Qt migration

| Phase | Action |
|-------|--------|
| Phase 1 (Window shell) | No drawing yet — no Y handling needed. |
| Phase 2 (Drawing primitives) | **DO NOT flip Y anywhere.** Pass Fortran-provided y values directly to `QPainter`. |
| Phase 3 (Text, fonts, colormap) | Same — text baseline y is top-left-origin, unflipped. |
| Phase 4 (Pixmap save/restore) | Same — pixmap coordinates are top-left-origin. |
| Phase 5 (Event bridge) | Mouse y from `QMouseEvent::y()` is already top-left-origin — no transform. |
| Phase 6 (UX chrome) | Widget coordinates are top-left-origin — no transform. |
| Phase 7 (Export) | **PRESERVE the `ypixels - y` flip verbatim** when porting the PostScript emitter from `xvuelc.c` to `xvue/xvue_qt_postscript.cpp`. Do NOT try to "clean it up" by storing flipped coordinates — the flip is a PostScript-format concern, not an internal-state concern. |
| HiDPI (`devicePixelRatioF > 1`) | Orthogonal to the Y-axis audit. Handled separately in Phase 1 (SHELL-06). |

## Anti-patterns (common mistakes to avoid)

- **"Store coordinates with Y-up because PostScript wants that"** — NO. Store
  coordinates as the Fortran side provides them (Y-down). Flip at emit time
  inside `xvpostscript_` only.
- **"Use `QTransform` to flip Y on the QPainter"** — NO. No y-flip transform
  anywhere in the on-screen path.
- **"Invert y at the Fortran boundary to 'normalize' to Y-up"** — NO. The
  Fortran `.f` wrappers and the C/C++ stubs are on the same side of the
  boundary, sharing the Y-down convention.

## Audit metadata

- **Audited:** 2026-04-10
- **Auditor:** Phase 0 planning (direct read of `xvue/xvuelc.c`)
- **Source lines inspected:** 319-332 (xvpxecran), 1187-1306 (xvpostscript),
  1619-1720 (window query), 1845-1975 (line drawing) — see `.planning/phases/00-build-skeleton-abi-stubs/00-RESEARCH.md §"Y-Axis Convention Audit"` for full details.
- **Read-only status:** This document is write-once in Phase 0 and read-only
  for every subsequent phase until Phase 9 retirement.
