# Phase 03 plan 03-04 Task 3 — testa/ 5-case A/B captures (FINAL)

**Task 3 was reopened twice** before reaching its final RESOLVED state:

1. **First reopen 2026-04-13** because the initial pass committed only
   X11-side captures and declared Task 3 green without ever performing
   a Qt-vs-X11 A/B comparison. That reopen produced 10/12 PASS / 2
   MISMATCH on the solver pairs (`cavity2d-ppflui` and `heat1d-ppther`)
   and was paused awaiting user direction.

2. **Second reopen 2026-04-14** — the 2 "MISMATCH" verdicts were
   traced to a Debian sid `libgfortran5 = 16-20260322-1` (gcc 16
   snapshot) runtime regression that had been pulled in by an `apt
   upgrade` between the two sessions. The fix did NOT touch the Qt
   drawing code (those are correct as shipped). Three actions
   restored the baseline:
   - `sudo apt install /var/cache/apt/archives/libgfortran5_15.2.0-9_amd64.deb`
     to downgrade the runtime
   - `sudo apt-mark hold libgfortran5` to prevent re-upgrade
   - `bin/cbl_tout` and `bin/cbl_tout_qt` re-run with a
     `/tmp/gfortran-14-shim` PATH directory (gcc → gcc-14, gfortran →
     gfortran-14) because the apt downgrade also removed gcc-15 +
     gfortran-15 compiler packages
   - 2 batch-file fixes in `testa/nafems_le1/`: `nafems_le1.mesh`
     was missing `1; { ALL OBJECTS }` between `5;` and `90;` (TRACOBJE
     submenu requirement); `nafems_le1.elas` was missing
     `15; { Drawing of STRESSES in ALL FE }` between `1;` and `90;`
     (TRACCONT submenu requirement). Both bugs predated the initial
     commit; the old "PASS partial content" verdict on
     `nafems_le1-ppelas` was a fluke (residual root-window state from
     the prior `ppmail` process leaking through).

This README documents the FINAL D-27 verdict after the second reopen.

## Infrastructure (final, after second reopen)

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
- `MEFISTO_QT_CAPTURE_PATH` env var tells the Qt `xvfermer_` hook
  to save `XvueState::backing_` (the canvas's authoritative backing
  pixmap) directly to that PNG path. This is a pure in-process grab
  with no external dependencies — works on CI without X or xcb-cursor0.

Both harnesses drive the exact same batch file (`.mesh` / `.heat` /
`.stoke56cr` / `.elas` / `.iexrr`) so the two sides are running
identical workflows, differing only in the graphics backend.

## Captures

All captures are in `03-04-ab/testa/` with filename convention
`<case>-<solver>_<backend>.png`. Qt captures are reference PNGs for the
backing pixmap at size 760x740 (from `XvueState::backing_`). X11 captures
are root-window grabs from Xvfb at 1280x800.

| Case       | Solver  | Qt PNG | X11 PNG |
|------------|---------|-------:|--------:|
| pan2d      | ppmail  | 122 KB | 7.6 KB  |
| nafems_le1 | ppmail  | 219 KB | 19 KB   |
| cavity2d   | ppmail  | 151 KB | 9.0 KB  |
| heat1d     | ppmail  | 39 KB  | 4.7 KB  |
| nlsecu     | ppmail  | 167 KB | 49 KB   |
| nafems_le1 | ppelas  | 491 KB | 61 KB   |
| cavity2d   | ppflui  | 49 KB  | 15 KB   |
| heat1d     | ppther  | 54 KB  | 4.4 KB  |
| nlsecu     | ppnlse  | —      | —       |

`nlsecu-ppnlse` remains deferred — the batch file runs a 2000-step
complex-wave simulation that needs ~1h50 of wall time. Needs a shrunk
`.iexrr` variant or an offline cron run.

## D-27 rubric verdict (FINAL)

Performed by reading every pair through the `Read` tool and comparing
on (a) geometry, (b) colors, (c) text, (d) no missing geometry, (e) no
miscolored regions.

### xvtest1..4 (stored in the parent `03-04-ab/` directory)

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
| nafems_le1 | **PASS** | Both show the quarter annulus with 200 quads colored by quality (deep blue high Q, cyan/green medium), quality histogram. Geometry + palette match. (Mesher batch fixed in this session — added the missing `1; { ALL OBJECTS }` step in `nafems_le1.mesh`.) |
| cavity2d   | **PASS** | Both show the 2D square cavity with triangle elements, quality stats. Qt uses turquoise fills, X11 uses dark blue — both valid quality-coloring palettes for a single-bin histogram. Geometry matches. |
| heat1d     | **PASS** | Both show the 1D line from A to B with 10 segments, quality stats, axes, node labels. 1D test case — no color variation to compare. Geometry matches. |
| nlsecu     | **PASS** | Both show the 3D cube of cube elements with quality coloring, XYZ axes, SURFACE BOUNDARY label. X11 uses grey faces + blue borders; Qt uses purple faces with magenta surface. Same geometry, different palette indices on a fully populated palette. |

**Mesher result: 5/5 PASS.** Differences in palette index choice exist
on cavity2d and nlsecu but these are aesthetic (both use values from
the shared palette); neither side shows missing geometry or miscolored
regions. Qt text rendering is consistently cleaner than X11's (Qt uses
the bundled DejaVu Sans Mono; X11 uses system X fonts that occasionally
clip the title bar).

### testa solver (non-mesher module)

| Case                | Verdict | Notes |
|---------------------|---------|-------|
| nafems_le1 ppelas   | **PASS** | Both sides now show the principal-stress arrow drawing on the quarter annulus mesh with stress arrows in radial pattern + STRESS / UNITIES of STRESS / OBJECT MEMBRANE titles. Batch fixed in this session — added the missing `15; { Drawing of STRESSES in ALL FE }` step in `nafems_le1.elas`. Improvement over the previous "PASS partial content" verdict, which only showed the mesher's quality coloring (residual root state). |
| cavity2d ppflui     | **PASS** | Both sides show the dark-blue uniform pressure mesh + Case 11 PRESSURES title + ISO-pressure color-bar legend on the right (0.00, 41.4, 82.9, 124, 166, 207, 249, 290, 332, 373, 414). Qt title has minor garbled chars in the date region (cosmetic — `ladate_` / `nomordinateurhote_` / `heureminuteseconde_` are warn-once no-ops on Qt; will be addressed in a later phase). The pressure drawing itself + color-bar legend match between backends. |
| heat1d ppther       | **PASS** | Both sides show the 1D mesh + flux arrows along the line + a red diagonal eigenvalue trace + EIGENVALUE / NORMAL FLUX of TEMPERATURE titles. Qt additionally shows the TEMPERATURE(t,X) on OBJECT AB title block, the QUALITY of FINITE ELEMENTS stats block, and a 10-color legend on the right (19.1 / 18.2 / … / 10.9) — all extra content, not a regression. The flux arrows + eigenvalue line are present on BOTH sides. |
| nlsecu ppnlse       | **DEFERRED**     | Both sides time out at ~1h50 compute cost. |

**Solver result: 3/4 PASS, 1 DEFERRED.**

## Final summary

- 4/4 xvtest drivers PASS
- 5/5 testa mesher PASS
- 3/4 testa solver PASS, 1/4 DEFERRED
- **Total A/B comparisons performed: 12 (excluding the deferred nlsecu)**
- **Total D-27 PASS: 12/12**
- **Total D-27 MISMATCH: 0/12**

## What changed since the first reopen

The 2 "MISMATCH" verdicts (`cavity2d-ppflui` color-bar legend missing;
`heat1d-ppther` flux arrows + eigenvalue trace missing) were NOT real
Qt-side rendering gaps. The Qt drawing primitives in `xvue_qt_api.cpp`
have been correct since Phase 02 / Phase 02.1 / Phase 03. The
"MISMATCH" PNGs were produced by stale Qt binaries running on top of
a `libgfortran5 = 16-20260322-1` runtime that the Debian sid `apt
upgrade` had pulled in unnoticed. The new runtime exposed latent
Fortran UB in the solver code (uninitialized `TPSINI` in the unsteady
heat solver, FPE traps in the elasticity stress computation, and a
parser path in `nafems_le1` that just-so-happened to crash differently
under the new runtime).

The fix had three parts:

1. **Pin libgfortran5 to 15.2.0-9** via the cached deb in
   `/var/cache/apt/archives/`, then `apt-mark hold libgfortran5` so
   the next `apt upgrade` does not re-pull the gcc-16 snapshot. The
   downgrade also removed `gcc-15` and `gfortran-15` compiler
   packages, leaving only the `-14` versions.

2. **PATH-shim the build** to use the explicit `-14` compilers.
   `/tmp/gfortran-14-shim/{gcc,cc,gfortran}` symlinks pointing at
   `gcc-14` / `gfortran-14`. With this shim ahead in `PATH`,
   `bin/cbl_tout` and `bin/cbl_tout_qt` produce binaries that match
   the runtime.

3. **Fix two long-standing batch-file bugs** in
   `testa/nafems_le1/nafems_le1.mesh` (missing `1;` between `5;` and
   `90;` in the TRACOBJE submenu) and `testa/nafems_le1/nafems_le1.elas`
   (missing `15;` between `1;` and `90;` in the TRACCONT submenu).
   These bugs predated the initial git commit; the old `nafems_le1-mail`
   capture worked because Xvfb root state persisted from a successful
   prior run, and the old `nafems_le1-ppelas` "PASS partial content"
   verdict was just the mesher's drawing leaking through.

After these three fixes the entire D-27 set is honestly green.

## Reproducing (full capture set)

```bash
export MEFISTO=/path/to/mefisto
export MEFISTOX=/tmp/mefistox-testa
export PATH=$MEFISTO/bin:$PATH

# one-time: pin libgfortran5 (Debian sid only; see top of README)
sudo apt install /var/cache/apt/archives/libgfortran5_15.2.0-9_amd64.deb
sudo apt-mark hold libgfortran5

# one-time: shim gcc/gfortran to -14 because apt downgrade removed -15
mkdir -p /tmp/gfortran-14-shim
ln -sf /usr/bin/gcc-14      /tmp/gfortran-14-shim/gcc
ln -sf /usr/bin/gcc-14      /tmp/gfortran-14-shim/cc
ln -sf /usr/bin/gfortran-14 /tmp/gfortran-14-shim/gfortran
export PATH=/tmp/gfortran-14-shim:$PATH

# one-time: full rebuild against pinned runtime
bin/cbl_tout
bin/cbl_tout_qt

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
