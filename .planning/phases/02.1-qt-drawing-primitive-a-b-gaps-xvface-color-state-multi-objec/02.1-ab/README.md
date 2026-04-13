# Phase 02.1 post-fix A/B captures

Captured after `fix(phase-02.1): xvfacetraits_ honors ncf/nca`.

- `xvtest2_qt_xcb.png` — Qt backend, QT_QPA_PLATFORM=xcb on Xvfb :99, post-fix.
  Must show red + blue filled triangles and magenta outline (see VALIDATION row 02.1-01-05).
- `xvtest4_qt_xcb.png` — Qt backend, xcb on Xvfb :99, post-fix.
  Must show colored cube faces and magenta dashed edges (see VALIDATION row 02.1-01-06).

Reference (legacy X11, frozen, DO NOT overwrite):
- `../../03-text-fonts-colormap/03-04-ab/xvtest2_x11.png`
- `../../03-text-fonts-colormap/03-04-ab/xvtest4_x11.png`

## Histogram results (recorded at capture time)

- xvtest2_qt_xcb.png: 682 unique colors (threshold ≥7 ✅), histogram contains #FF0000 (red), #0000FF (blue), #FF00FF (magenta).
- xvtest4_qt_xcb.png: 1105 unique colors (threshold ≥5 ✅), histogram contains #FF0000, #0000FF, #FF00FF (dashed edges).
