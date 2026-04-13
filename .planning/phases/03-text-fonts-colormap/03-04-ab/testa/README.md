# Phase 03 plan 03-04 Task 3 — testa/ 5-case A/B captures (REDONE)

**Task 3 was reopened 2026-04-13** because the initial pass committed
only X11-side captures and declared Task 3 green without ever performing
a Qt-vs-X11 A/B comparison. This README documents the honest Qt+X11
capture set and the resulting D-27 rubric application.

## Infrastructure (final, after Task 3 reopen)

Legacy X11 side — captured via `bin/testa-capture.sh`:
- Xvfb :99 1280x800x24
- `MEFISTO_BATCH_X11=1` upgrades the solver's `INTERA=0` batch mode
  to `INTERA=1` so the legacy X11 window opens while the batch file
  still drives the workflow
- `MEFISTO_XVSOURIS_AUTOEXIT=1` makes `xvsouris_` / `xvsouris2_` flush
  and return a synthetic keypress
- `MEFISTO_XVFERMER_READY_FILE` + `MEFISTO_XVFERMER_HOLD_MS` tell
  `xvfermer_` to touch a sentinel file + hold before destroying
- External capture via `import -display :99 -window root`

Qt side — captured via `bin/qt-capture.sh`:
- `QT_QPA_PLATFORM=offscreen` — no X server required
- Same `MEFISTO_BATCH_X11` + `MEFISTO_XVSOURIS_AUTOEXIT` vars (the same
  env-var hooks were added to `xvue/qt/src/xvue_qt_api.cpp`)
- NEW: `MEFISTO_QT_CAPTURE_PATH` env var tells the Qt `xvfermer_` hook
  to save `XvueState::backing_` (the canvas's authoritative backing
  pixmap) directly to that PNG path. This is a pure in-process grab
  with no external dependencies — works on CI without X or xcb-cursor0.

Both harnesses drive the exact same batch file (`.mesh` / `.heat` /
`.stoke56cr` / `.elas` / `.iexrr`) so the two sides are running
identical workflows, differing only in the graphics backend.

## Captures

All captures are in `03-04-ab/testa/` with filename convention
`<case>-<solver>_<backend>.png`. Qt captures are reference PNGs for the
backing pixmap at size 760x760 (from `XvueState::backing_`). X11 captures
are root-window grabs from Xvfb at 1280x800.

| Case       | Solver  | Qt PNG | X11 PNG |
|------------|---------|-------:|--------:|
| pan2d      | ppmail  | 592 KB | 7.6 KB  |
| nafems_le1 | ppmail  | 153 KB | 24 KB   |
| cavity2d   | ppmail  | 420 KB | 8.9 KB  |
| heat1d     | ppmail  | 48 KB  | 4.7 KB  |
| nlsecu     | ppmail  | 137 KB | 57 KB   |
| nafems_le1 | ppelas  | 174 KB | 21 KB   |
| cavity2d   | ppflui  | 306 KB | 15 KB   |
| heat1d     | ppther  | 21 KB  | 4.4 KB  |
| nlsecu     | ppnlse  | —      | —       |

`nlsecu-ppnlse` remains deferred — the batch file runs a 2000-step
complex-wave simulation that needs ~1h50 of wall time. Needs a shrunk
`.iexrr` variant or an offline cron run.

## D-27 rubric verdict (HONEST)

Performed by reading every pair through the `Read` tool and comparing
on (a) geometry, (b) colors, (c) text, (d) no missing geometry, (e) no
miscolored regions.

### xvtest1..4 (stored in the parent 03-04-ab/ directory)

| Driver    | Verdict | Notes |
|-----------|---------|-------|
| xvtest1   | **PASS** | Both sides show the 32-color palette bars, blue X-rect, red-on-green, spectrum ramp. Qt additionally renders the font catalog with crisp text; X11 only has the color bars. Qt gain is a TEXT/FONT bonus, not a regression. |
| xvtest2   | **PASS** | Both sides show red+blue triangles, magenta hex outline, yellow diagonal. Qt adds text labels (TEXTE2D/SYMBOLE2D/.123/.3.1416/XVTEXTE). |
| xvtest3   | **PASS** | Both sides show cube + inscribed blue tetrahedron + red face + dashed magenta edges + dashed yellow XYZ axes. Qt adds LONGITUDE/LATITUDE banner + SYMBOLE3D/TEXTE3D + 3.1416 labels. |
| xvtest4   | **PASS** | Both sides show cube + blue tetrahedron + red face + dashed magenta edges. Near-pixel-level match. |

**xvtest result: 4/4 PASS.**

### testa mesher (MAILLER/ppmail)

| Case       | Verdict | Notes |
|------------|---------|-------|
| pan2d      | **PASS** | Both show the 2D panel mesh with quality stats block, quad elements in dark blue. Geometry and color palette match. |
| nafems_le1 | **PASS** | Both show the quarter annulus with 200 quads colored by quality (deep blue high Q, cyan/green medium), quality histogram. Geometry + palette match. |
| cavity2d   | **PASS** | Both show the 2D square cavity with triangle elements, quality stats. Qt uses turquoise fills, X11 uses dark blue — both valid quality-coloring palettes for a single-bin histogram. Geometry matches. |
| heat1d     | **PASS** | Both show the 1D line from A to B with 10 segments, quality stats, axes, node labels. 1D test case — no color variation to compare. Geometry matches. |
| nlsecu     | **PASS** | Both show the 3D cube of cube elements with quality coloring, XYZ axes, SURFACE BOUNDARY label. X11 uses grey faces + blue borders; Qt uses filled cube with magenta surface. Same geometry, different choice of palette indices. |

**Mesher result: 5/5 PASS.** Differences in palette index choice exist
on cavity2d and nlsecu but these are aesthetic (both use values from
the shared palette); neither side shows missing geometry or miscolored
regions. Qt text rendering is consistently cleaner than X11's (Qt uses
the bundled DejaVu Sans Mono; X11 uses system X fonts that occasionally
clip the title bar).

### testa solver (non-mesher module)

| Case                | Verdict | Notes |
|---------------------|---------|-------|
| nafems_le1 ppelas   | **PASS (partial content)** | Both sides show the quarter annulus mesh with quality coloring. Neither side shows the deformed-shape / stress visualization implied by the batch command `8; 1; 90; 90;` (DESSIN Deformées Contraintes, sub-option 1 = "mesh quality"). So both are rendering the same sub-option content and the A/B matches. |
| cavity2d ppflui     | **MISMATCH**     | X11 shows the triangulated mesh in **dark blue** (uniform low pressure) with the "PRESSURE(t,X)" title and color bar legend 0..414. Qt shows the same mesh in **uniform turquoise/cyan** without the color bar legend. Qt's color choice maps to a different palette index and the color scale rendering is missing. Geometry matches, but colors differ and Qt drops the color bar. |
| heat1d ppther       | **MISMATCH**     | X11 shows the 1D mesh + flux arrows along the line + a red diagonal eigenvalue trace + "EIGENVALUE" and "NORMAL FLUX of TEMPERATURE" titles. Qt shows only the 1D mesh + AB labels — **the flux arrows and eigenvalue trace are entirely absent on the Qt side**. The batch file `8; 4; 15; 90;` dispatches (via trther.f) to `TR1DTER` / `TRFLUX` / `TRERTH` which draw via primitives. On X11 the arrows + trace are present; on Qt they are not. This is a real Phase 3 Qt-side rendering gap. |
| nlsecu ppnlse       | **DEFERRED**     | Both sides time out at ~1h50 compute cost. |

**Solver result: 1 PASS, 2 MISMATCH, 1 DEFERRED.**

## Honest summary

- 4/4 xvtest xvdrivers PASS
- 5/5 testa mesher PASS
- 1/4 testa solver PASS, 2/4 MISMATCH, 1/4 DEFERRED
- **Total A/B comparisons performed: 12 (excluding the deferred nlsecu)**
- **Total D-27 PASS: 10/12**
- **Total D-27 MISMATCH: 2/12**

## What the 2 solver mismatches reveal

**`cavity2d-ppflui` — color scale legend missing on Qt side.**
The batch `1; 2; 90;` sub-menu in `cavity2d.stoke56cr` triggers the
"ISO-PRESSURE COLOR ZONES" drawing. X11 renders a legend bar on the
right showing the pressure range 0..414. Qt does not. The ISO-pressure
sub-tracer in `flui/` probably calls a drawing primitive that either
is a warn-once stub on Qt or uses coordinates outside the Qt canvas
viewport. Needs investigation in `flui/trco2d.f` / `flui/trvi2d.f`
and the Qt `xvue_qt_api.cpp` stub list.

**`heat1d-ppther` — flux arrows and eigenvalue trace missing on Qt.**
The batch `8; 4; 15; 90;` in `heat1d.heat` triggers "DRAWING of ERRORS"
sub-option 4. On X11 this produces flux arrows (triangle markers) and
a diagonal red eigenvalue line plus a title. On Qt none of these
secondary visualisation elements appear, only the 1D mesh line. The
thermal trace sub-routine (`ther/tr1dter.f` or `ther/trflux.f` or
`ther/trerth.f`) must use a primitive that is not fully wired on Qt.

Neither mismatch is a `trelas.f`-style missing-`MEMPXFENETRE` issue
(the Qt capture reads `XvueState::backing_` directly, so the drawings
ARE missing from the backing itself — they never got painted).

Both mismatches are **Phase 3 Qt-side rendering gaps** that were never
caught by the xvtest0 coverage driver because xvtest0's Phase 3 section
only exercises `xvtexte_` / `xvcouleur_` / `xvactivervb_` / `xvnbpixel
texte_` / `xvtrait_` primitives, not the higher-level solver trace
sub-routines.

## Recommended follow-up

- **NOT** a silent-defer situation. Phase 3 should not be marked
  `nyquist_compliant: true` until either:
  1. The two Qt-side solver trace mismatches are investigated and
     either fixed (if a Qt primitive is missing) or explicitly
     documented as "same geometry, palette-index difference only,
     accepted deviation" with a commit message reasoning it through.
  2. OR a dedicated gap-closure phase (03.2?) is opened to own the
     trother/trcoefse/trco2d Qt trace investigation.
- `nlsecu-ppnlse` deferred status is separate and acceptable as-is
  (computational cost, needs a test-only shrunk variant).

## Reproducing (full capture set)

```bash
export MEFISTO=/path/to/mefisto
export MEFISTOX=/tmp/mefistox-testa
export PATH=$MEFISTO/bin:$PATH

# one-time: INITIER each project dir
for case in pan2d nafems_le1 cavity2d heat1d nlsecu; do
  rm -rf $MEFISTOX/$case
  mkdir -p $MEFISTOX/$case
  (cd $MEFISTOX/$case && for f in $MEFISTO/testa/$case/*; do ln -sf "$f" .; done
   echo $case | $MEFISTO/pp/ppinit > /dev/null)
done

AB=$MEFISTO/.planning/phases/03-text-fonts-colormap/03-04-ab/testa

# X11 mesher captures
pkill -9 Xvfb 2>/dev/null
for spec in "pan2d:pan2d.mesh" "nafems_le1:nafems_le1.mesh" \
            "cavity2d:cavity2d.meshbf" "heat1d:heat1d.mesh" \
            "nlsecu:nlsecu.meshq2"; do
  c=${spec%:*}; f=${spec#*:}
  bin/testa-capture.sh $MEFISTOX/$c ppmail $f $AB/${c}-mail_x11.png 2000 99
done

# X11 solver captures (skip nlsecu)
bin/testa-capture.sh $MEFISTOX/nafems_le1 ppelas nafems_le1.elas    $AB/nafems_le1-ppelas_x11.png 2000 99
bin/testa-capture.sh $MEFISTOX/cavity2d  ppflui cavity2d.stoke56cr $AB/cavity2d-ppflui_x11.png 2000 99
bin/testa-capture.sh $MEFISTOX/heat1d    ppther heat1d.heat        $AB/heat1d-ppther_x11.png   2000 99

# Qt mesher captures (no Xvfb needed; offscreen plugin)
for spec in "pan2d:pan2d.mesh" "nafems_le1:nafems_le1.mesh" \
            "cavity2d:cavity2d.meshbf" "heat1d:heat1d.mesh" \
            "nlsecu:nlsecu.meshq2"; do
  c=${spec%:*}; f=${spec#*:}
  bin/qt-capture.sh pp/ppmail_qt $AB/${c}-mail_qt.png $f 500 $MEFISTOX/$c
done

# Qt solver captures (skip nlsecu)
bin/qt-capture.sh pp/ppelas_qt $AB/nafems_le1-ppelas_qt.png nafems_le1.elas    500 $MEFISTOX/nafems_le1
bin/qt-capture.sh pp/ppflui_qt $AB/cavity2d-ppflui_qt.png  cavity2d.stoke56cr 500 $MEFISTOX/cavity2d
bin/qt-capture.sh pp/ppther_qt $AB/heat1d-ppther_qt.png    heat1d.heat        500 $MEFISTOX/heat1d
```
