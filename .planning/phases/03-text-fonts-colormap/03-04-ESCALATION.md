---
plan: 03-04
task: 2
status: blocking_fail
created: 2026-04-13
ab_capture_env: Xvfb :99 1280x800x24 + QT_QPA_PLATFORM=xcb (primary) and kwin-mcp Qt native Wayland (secondary — same results)
---

# Plan 03-04 Task 2 — Escalation

The xvtest1..4 human A/B visual gate failed per the D-27 rubric. Per the plan's `resume-signal` contract, this escalates and STOPS Task 2. Phase 3 cannot complete until the failures are resolved in a reopened fix plan.

## Capture environment

- **Qt side:** `QT_QPA_PLATFORM=xcb pp/ppxvtest{1..4}_qt` under Xvfb :99 (1280x800x24). Identical results reproduced under native Wayland Qt via kwin-mcp — ruling out the platform layer.
- **X11 side:** `pp/ppxvtest{1..4}` under Xvfb :99 (same display).
- **Capture tool:** `import -window root -display :99`.
- **Artifacts:** `.planning/phases/03-text-fonts-colormap/03-04-ab/{xvtest{1..4}_qt_xcb,xvtest{1..4}_x11}.png`

## Per-driver results

### xvtest1
- **Qt xcb:** full font/color catalog renders — fonts numbered 14..79, color bars, blue X-diagonaled square, red circle, green triangle. Title `EN ATTENTE D'UN CARACTERE OU CLIC SOURIS`.
- **X11:** pre-existing crash (`xvue/xvuelc.c` legacy bug — libc-start backtrace in log). Out of Phase 3 scope; documented in `.planning/phases/03.1-.../deferred-items.md`.
- **Rubric verdict:** deferred — cannot compare because legacy backend crashes. Not a Phase 3 regression.

### xvtest2 — **RUBRIC FAIL**
- **X11:** magenta hexagon outline + **red triangle (left)** + **blue triangle (right)** + white inner rectangle with magenta diagonal + yellow median line + color-bar placeholder lines.
- **Qt xcb:** white hexagon fill only + yellow median + rendered text labels (`TEXTE2D`, `SYMBOLE2D`, `123`, `3.1416`, `XVTEXTE`).
- **Missing from Qt:**
  - Red triangular fill (left tip of hexagon)
  - Blue triangular fill (right tip of hexagon)
  - Magenta hexagon outline
  - Inner magenta diagonal
- **Rubric failures:**
  - (a) Geometry present — Qt drops two filled triangles and the outline.
  - (b) Colors match — red and blue fills are absent, not just miscolored.
  - (e) No miscolored regions — hexagon appears all-white in Qt instead of outline-only-plus-filled-triangles.
- **Suspected root cause:** `xvue/qt/src/xvue_qt_api.cpp::xvface_` or whatever handles `CALL XVFACE` in xvtest2.f is ignoring or overwriting per-face brush color. Could also be an `xvcouleur_` state-flush problem before the fill, or a painter state that isn't saved/restored between successive face draws.
- **Files to inspect first:** `xvue/qt/src/xvue_qt_api.cpp::xvface_`, `xvue/qt/src/xvue_qt_api.cpp::xvcouleur_`, `xvue/qt/src/xvue_qt_canvas.cpp::applyPen`, `prpr/xvtest2.f:40-90` (polygon draw sequence).

### xvtest3
- **Qt xcb:** 3D octahedron with magenta/white faces, yellow X/Y/Z coordinate axes, `LONGITUDE= 35 degres  LATITUDE= 25 degres` banner, `SYMBOLE3D`, `TEXTE3D`, `234567`, `3.1416` labels.
- **X11:** blank — driver blocks on an early `XVSOURIS` call and no mouse event was injected. This is not a real A/B — the legacy backend never got past its first input wait.
- **Rubric verdict:** deferred — cannot compare without injecting mouse input to the X11 session. Not a Phase 3 regression per se, but xvtest3 needs to be re-run with an event injector (xdotool or equivalent) before it can be approved.

### xvtest4 — **RUBRIC FAIL**
- **X11:** wireframe cube with **inscribed blue tetrahedron** + red triangular face + dashed magenta cube edges + white triangular faces.
- **Qt xcb:** rotated solid white cube only.
- **Missing from Qt:**
  - Inscribed tetrahedron (the second 3D object in the plotlist)
  - Red triangular face
  - Blue tetrahedron faces
  - Dashed-stroke pen style on the cube edges (Qt uses solid strokes)
- **Rubric failures:**
  - (a) Geometry present — entire inscribed tetrahedron is absent.
  - (b) Colors match — red face and blue tetrahedron fills missing.
  - (e) No miscolored regions — cube appears solid-white instead of wireframe-with-inscribed-shape.
- **Suspected root cause:**
  1. Multi-object composition: the Qt backend may be clearing the polygon stack or the painter between object draws in a way that wipes the inscribed shape.
  2. Pen-style mapping: the legacy `CALL XVTRAITS` (or whatever sets dashed strokes) is likely a no-op on the Qt side.
  3. Color-per-face: same suspected root cause as xvtest2 — face color state not being preserved between successive `xvface_` calls.
- **Files to inspect first:** `xvue/qt/src/xvue_qt_api.cpp::xvface_`, line-attribute setters (`xvtraits`/`xvtypetrait`), `xvue/qt/src/xvue_qt_canvas.cpp::applyPen` (pen style), `prpr/xvtest4.f:20-80`.

## Aggregate verdict

| Driver | Qt xcb | X11 | D-27 Rubric | Action |
|--------|--------|-----|------------|--------|
| xvtest1 | renders | pre-existing crash | deferred | out of scope (legacy bug) |
| xvtest2 | partial | full | **FAIL (a, b, e)** | reopen fix plan |
| xvtest3 | renders | blocked on input | deferred | needs event injector to retest |
| xvtest4 | partial | full | **FAIL (a, b, e)** | reopen fix plan |

**Two hard rubric failures.** Phase 3 is NOT complete. The fix belongs in reopened plan 03-01 (imposed_defaults_fill / color-state wiring) or 03-02 (body implementations) — most likely 03-02, since the divergences are in the actual Qt-side drawing bodies (`xvface_`, `xvcouleur_`, `xvtraits_`).

## Recommended next steps

1. **Write a gap-closure plan for phase 03** via `/gsd-plan-phase 03 --gaps` once `03-VERIFICATION.md` carries the failing gaps. The gap plan would reopen the color-fidelity and multi-object-composition bodies in the Qt backend.

2. **Alternative:** these divergences are actually drawing-primitive gaps (phase 2 territory). A cleaner routing is `/gsd-insert-phase` to create a `02.1` polish phase dedicated to closing `xvface_` color state + multi-object 3D composition + dashed pen style, leaving phase 03 untouched.

3. **xvtest3 re-capture:** install `xdotool` (or equivalent) and rerun the xvtest3 A/B with a mouse-click injected after the first 2 seconds, so the legacy backend advances past its initial `XVSOURIS` and actually draws.

4. **xvtest1 legacy crash:** already deferred in `.planning/phases/03.1-.../deferred-items.md`. Not blocked by phase 03.

## Files referenced

- `.planning/phases/03-text-fonts-colormap/03-04-PLAN.md` (this task's source)
- `.planning/phases/03-text-fonts-colormap/03-04-ab/` (Xvfb xcb A/B captures)
- `.planning/phases/03.1-build-xvtest1-4-driver-infrastructure-on-both-backends-qt-le/visual-gate/` (earlier Wayland Qt captures — same results)
- `.planning/phases/03.1-.../deferred-items.md` (pre-existing legacy crashes)
