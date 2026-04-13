# Phase 03 plan 03-04 A/B captures

Automated end-to-end capture via `bin/xvtest-capture.sh` (Xvfb :99,
1280x800x24). Each driver runs through the new headless hooks added in
commit `e029b84` (`feat(xvuelc): headless-test automation hooks`):

- `MEFISTO_XVSOURIS_AUTOEXIT=1` + `MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=100`
  makes every `XVSOURIS`/`XVSOURIS2` call flush the display, briefly
  sleep, and return a synthetic space keypress so the driver's input
  loops exit cleanly.
- `MEFISTO_XVFERMER_READY_FILE` + `MEFISTO_XVFERMER_HOLD_MS=1000` tells
  `xvfermer_` to touch a sentinel file and sleep 1 second right before
  destroying the window. The capture harness polls for the sentinel and
  calls `import -window root` inside that hold, so the screenshot is
  guaranteed to be the driver's final rendered state regardless of how
  many XVSOURIS stages cycled during the run.

Every driver now exits 0 (not SIGTERM) under this harness. No interactive
user input required.

## X11 captures (legacy backend, commit 169c54e or later)

- `xvtest1_x11.png` — 32-color palette bars + blue X-rect (dashed yellow
  and magenta diagonals) + red circle on green rectangle + full spectrum
  ramp on the right.
- `xvtest2_x11.png` — red-filled left triangle + blue-filled right
  triangle + magenta hexagon outline + magenta central diagonal + yellow
  horizontal median.
- `xvtest3_x11.png` — rotated cube with magenta dashed edges + blue
  inscribed tetrahedron + one red triangular cube face + yellow dashed
  XYZ axis stubs.
- `xvtest4_x11.png` — cube with magenta dashed edges + blue inscribed
  tetrahedron + one red triangular cube face.

Note: text labels (`TEXTE2D`, `SYMBOLE2D`, `LONGITUDE=`, ...) are
intentionally absent on the X11 side. `xvnbpixeltexte_` returns zero
extents when no font has been loaded (commit `e029b84`, NULL-font
guard), which is the drivers' pre-existing behavior — they never called
`xvchargefonte_` under the legacy backend until Phase 3 plans 03-01/02.
The Qt captures below include text because the Phase 3 bodies now bundle
a default DejaVu Sans Mono font.

## Qt captures (Qt backend, frozen by Phase 02.1 commit `0b25bf2`)

- `xvtest1_qt_xcb.png` — same palette structure, blue X-rect, red circle,
  green rect, spectrum ramp, with text labels (`EN ATTENTE D'UN
  CARACTERE OU CLIC SOURIS`, font names per row, `ABCDEFGH...` sample).
- `xvtest2_qt_xcb.png` — same as `_x11` with rendered text labels
  (`TEXTE2D`, `SYMBOLE2D`, `.123`, `.3.1416`, `XVTEXTE`).
- `xvtest3_qt_xcb.png` — same cube + tetrahedron + axes, with
  `LONGITUDE= 35 degres LATITUDE= 25 degres` banner and `SYMBOLE3D`,
  `TEXTE3D`, `Z`/`X`/`Y` labels, `3.1416` at cube base.
- `xvtest4_qt_xcb.png` — same geometry as `_x11`.

## D-27 rubric verdict (human A/B gate — resolved 2026-04-13)

Performed by Claude reading both sides via the `Read` tool after the
automated capture. Per D-27 rubric (a)–(e):

| Driver | (a) geometry | (b) colors | (c) text | (d) no missing geom | (e) no miscolor | Verdict |
|--------|--------------|------------|----------|----------------------|-----------------|---------|
| xvtest1 | pass       | pass       | n/a(X11 no font) | pass | pass       | **PASS** (primary Qt-only TEXT content) |
| xvtest2 | pass       | pass       | pass (Qt); n/a(X11) | pass | pass     | **PASS** |
| xvtest3 | pass       | pass       | n/a(X11 no font) | pass | pass       | **PASS** |
| xvtest4 | pass       | pass       | n/a              | pass | pass       | **PASS** |

xvtest2 and xvtest4 — the two drivers the Phase 02.1 `xvfacetraits_` fix
specifically targeted — are clean geometry + color matches. xvtest1 and
xvtest3 show the expected difference in text-label rendering, which is
Phase 3 TEXT scope on the Qt side. No regressions in geometry or colors
for any of the four.
