---
gate: phase-3-ab-visual
captured: 2026-04-13T01:00:00Z
captured_by: orchestrator (kwin-mcp + Xvfb)
qt_session: kwin_wayland --virtual 1280x800 (kwin-mcp session)
x11_session: Xvfb :99 1280x800x24
status: executed — divergences found, out of phase 03.1 scope
---

# Phase 3 A/B Visual Gate — Execution Record

Phase 03.1's deliverable was **make the A/B gate executable**. This directory records the first full execution of that gate after infrastructure was in place.

## Capture environments

- **Qt side:** kwin-mcp isolated session (virtual kwin_wayland 1280x800). Native Qt rendering through `libxvueqt.a`.
- **X11 side:** Xvfb :99 at 1280x800x24 with default font path. Captured via `import -window root`. No mouse input injected, so drivers that block on `XVSOURIS` before first draw show empty screens.

## Captures

| Driver | Qt | X11 (Xvfb) | Notes |
|--------|----|------------|-------|
| xvtest1 | [ppxvtest1_qt.png](ppxvtest1_qt.png) | — | Full font/color catalog — Qt shows fonts numbered 14..79, color bars, blue X-box, red circle, green triangle. X11 side not captured (blocks on `XVSOURIS`). |
| xvtest2 | [ppxvtest2_qt.png](ppxvtest2_qt.png) | [ppxvtest2_x11.png](ppxvtest2_x11.png) | **Divergence.** X11 renders filled red/blue triangles + magenta hex outline + inner diagonal; Qt only renders the white hex fill + yellow median line. Qt DOES render text labels (TEXTE2D, SYMBOLE2D, XVTEXTE, 123, 3.1416) where X11 falls back to placeholder lines (Xvfb lacks glyph cache). |
| xvtest3 | [ppxvtest3_qt.png](ppxvtest3_qt.png) | [ppxvtest3_x11.png](ppxvtest3_x11.png) | Qt shows 3D octahedron with coord axes + LATITUDE/LONGITUDE banner. X11 capture is blank because the driver blocks on an early `XVSOURIS` call and no mouse event was injected. Not a backend difference. |
| xvtest4 | [ppxvtest4_qt.png](ppxvtest4_qt.png) | [ppxvtest4_x11.png](ppxvtest4_x11.png) | **Divergence.** X11 renders a wireframe cube with an inscribed blue tetrahedron, red face, and dashed magenta edges. Qt renders a single solid white cube with no inscribed solid, no multi-color faces, and no dashed stroke style. |

## Divergences summary (out of phase 03.1 scope)

These findings belong to future drawing-primitive gap plans, **not** to phase 03.1:

1. **Filled polygon color fidelity (xvtest2)** — Qt skips the red/blue triangular fills. Suspected gap in the Qt implementation of `xvface_` or the foreground-color setter for filled primitives.
2. **Multi-object 3D composition (xvtest4)** — Qt draws only the outer cube, not the inscribed tetrahedron. Either a second-object draw call is being dropped, or the polygon stack is being cleared between objects.
3. **Dashed stroke style (xvtest4)** — Qt appears to use solid strokes only; legacy X11 uses dashed lines for the tetrahedron edges. Gap in the Qt pen-style mapping.

## What phase 03.1 confirmed

- ✓ All 5 Qt drivers run without crashing and without the `xtinit_` warn-once stub
- ✓ All 5 legacy X11 drivers build and initialize cleanly (512 fonts, 256 colors, TrueColor) under Xvfb
- ✓ The hardened 03-04 Task 1 smoke harness is sufficient to catch future regressions
- ✓ Visual A/B comparison is now physically possible — the gate is **executable**, which was the phase goal

## What phase 03.1 did NOT promise

Parity between Qt and legacy X11 rendering. That work belongs to subsequent gap plans under phase 2 (drawing primitives) or a new dedicated phase.
