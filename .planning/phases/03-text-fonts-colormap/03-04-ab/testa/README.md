# Phase 03 plan 03-04 Task 3 — testa/ 5-case A/B captures

Automated captures produced via the new hybrid batch+X11 infrastructure
(commits `e029b84`, `3149e3f`, `fix(xvuelc)` series, `feat(prpr):
MEFISTO_BATCH_X11 env-var override`, `fix(flui/lifiviprte): format
descriptors`, `feat(bin): testa-capture.sh`).

## Infrastructure

For each case the harness:

1. Exports `MEFISTO_BATCH_X11=1` which tells `prpr/pp<solver>.f` to keep
   batch-file-driven flow (`INTERA=0` semantics: read `.mesh`/`.elas`/
   `.heat`/... from disk) but upgrade to `INTERA=1` so `XTINIT + XVINIT`
   open a legacy X11 window and the embedded drawing commands produce
   visible output. INTERA=1 (not 3) ensures `xvue/lereur.f:64` returns
   from errors without waiting for a click.
2. Exports `MEFISTO_XVSOURIS_AUTOEXIT=1` so any stray `xvsouris_` call
   returns a synthetic keypress.
3. Exports `MEFISTO_XVFERMER_READY_FILE` + `MEFISTO_XVFERMER_HOLD_MS`
   so `xvfermer_` touches a sentinel file + holds before destroying
   the window. The capture harness polls for the sentinel and runs
   `import -display :99 -window root`.
4. `xvfermer_` force-copies `mempx → fenetre_mef` right before the
   sentinel so whatever is in the off-screen pixmap lands on the
   visible window at capture time.

Xvfb is started fresh at 1280x800x24.

## Coverage

| Case             | Mesher (`ppmail <data>`)          | Solver (`pp<solver> <data>`)            |
|------------------|-----------------------------------|------------------------------------------|
| pan2d            | `pan2d-mail_x11.png` ✓ 7 643 B    | — (mesher-only case)                     |
| nafems_le1       | `nafems_le1-mail_x11.png` ✓ 24 277 B | `ppelas` — deferred (see below)       |
| cavity2d         | `cavity2d-mail_x11.png` ✓ 8 958 B | `cavity2d-ppflui_x11.png` ✓ 14 841 B    |
| heat1d           | `heat1d-mail_x11.png` ✓ 4 681 B   | `heat1d-ppther_x11.png` ✓ 4 391 B       |
| nlsecu           | `nlsecu-mail_x11.png` ✓ 56 737 B  | `ppnlse` — deferred (see below)         |

**Mesher captures: 5/5 PASS** on all canonical testa cases — proves the
hybrid batch+X11 infrastructure across the entire `pp/ppmail` path.

**Solver captures: 2/4 PASS** (`ppther` on heat1d, `ppflui` on cavity2d).

## What each PASS shows

- `pan2d-mail_x11.png`: 2D mesh of panel with multi-line object outline,
  node numbers, quality stats, coordinate axes.
- `nafems_le1-mail_x11.png`: quarter-annulus mesh with 200 quad elements,
  quality color-coded (deep blue = high Q, cyan = medium), quality
  histogram legend, XY axes, point labels A/B/C/D/Q.
- `cavity2d-mail_x11.png`: 2D square cavity mesh with triangular elements.
- `heat1d-mail_x11.png`: 1D mesh from point A to point B with 10 segments,
  node indices, quality statistics.
- `nlsecu-mail_x11.png`: 2D non-linear elasticity mesh.
- `heat1d-ppther_x11.png`: transient thermal result — `EIGENVALUE` title,
  `NORMAL FLUX of TEMPERATURE` arrows along the mesh, red line showing
  flux profile.
- `cavity2d-ppflui_x11.png`: Stokes pressure field — `PRESSURE(t,X) on the
  OBJ`, `Case 11 The PRESSURES at TIME 1.00000 MIN=0.00000 MAX=2.15786`,
  full color scale legend 0–414.

## Deferred solver cases

### `nafems_le1-ppelas` — mempx drawing-path semantics

`trelas.f:249` calls `EFFACEMEMPX` then delegates to sub-tracers
(`TRCONT`, `TRDEPL`, `TRVMTR`, ...) that draw to `mempx` without
explicit `MEMPXFENETRE` calls. With my `xvfermer_` copy hook active,
`mempx` should still contain the last drawings at teardown — but the
captured window is uniformly empty (background color). Hypothesis:
something in the batch+X11 teardown path issues another `EFFACEMEMPX`
after the trace, wiping `mempx` before `xvfermer_` fires.

Fixing this requires tracing exactly which sub-path of `trelas.f`
receives the `8; 1; 90; 90;` command tokens and whether the final
`FERMER;` token triggers an additional clear. Out of scope for this
session; tracked as a follow-up.

### `nlsecu-ppnlse` — computation time

`testa/nlsecu/nlsecu.iexrr` runs a non-linear wave simulation:
- Final time 20 s, step 0.01 s → 2 000 time steps
- Each step involves a complex-matrix solve on ~4 961 nodes
- Observed: step 88/2 000 reached in 300 s → ~1h50 total at this rate

Too long for on-line capture. Either the test case parameters need
shrinking (a test-only `.iexrr` variant) or the capture has to be a
cron-scheduled offline run.

## D-27 rubric verdict — partial

| Case                       | (a) geometry | (b) colors | (c) text | (d) no missing geom | (e) no miscolor | Verdict        |
|----------------------------|--------------|------------|----------|---------------------|-----------------|----------------|
| pan2d mesher               | pass         | pass       | pass     | pass                | pass            | **PASS**       |
| nafems_le1 mesher          | pass         | pass       | pass     | pass                | pass            | **PASS**       |
| cavity2d mesher            | pass         | pass       | pass     | pass                | pass            | **PASS**       |
| heat1d mesher              | pass         | pass       | pass     | pass                | pass            | **PASS**       |
| nlsecu mesher              | pass         | pass       | pass     | pass                | pass            | **PASS**       |
| heat1d thermal             | pass         | pass       | pass     | pass                | pass            | **PASS**       |
| cavity2d fluid             | pass         | pass       | pass     | pass                | pass            | **PASS**       |
| **nafems_le1 elasticity**  | —            | —          | —        | —                   | —               | **DEFERRED**   |
| **nlsecu non-linear**      | —            | —          | —        | —                   | —               | **DEFERRED**   |

**7/9 configurations PASS.** The two deferred cases surface pre-existing
drawing-path / compute-time issues unrelated to the Phase 3 TEXT/FONTS/
COLORMAP scope. Task 3 delivers the automation infrastructure and the
first working capture set; the two remaining cases can be completed in
a follow-up without touching Phase 3 code.

## Reproducing

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

# mesher captures (all 5)
AB_DIR=$MEFISTO/.planning/phases/03-text-fonts-colormap/03-04-ab/testa
for spec in \
  "pan2d:ppmail:pan2d.mesh" \
  "nafems_le1:ppmail:nafems_le1.mesh" \
  "cavity2d:ppmail:cavity2d.meshbf" \
  "heat1d:ppmail:heat1d.mesh" \
  "nlsecu:ppmail:nlsecu.meshq2" ; do
  c=${spec%%:*}; rest=${spec#*:}; s=${rest%%:*}; f=${rest#*:}
  bin/testa-capture.sh $MEFISTOX/$c $s $f $AB_DIR/${c}-mail_x11.png 2000 99
done

# solver captures (2 that work)
bin/testa-capture.sh $MEFISTOX/heat1d    ppther  heat1d.heat           $AB_DIR/heat1d-ppther_x11.png   3000 99
bin/testa-capture.sh $MEFISTOX/cavity2d  ppflui  cavity2d.stoke56cr    $AB_DIR/cavity2d-ppflui_x11.png 3000 99
```
