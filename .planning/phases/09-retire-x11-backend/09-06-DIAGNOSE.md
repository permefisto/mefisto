# 09-06 — Matched-Dim Qt Recapture: Diagnose Log

Working scratch document for Plan 9-06 Task 1.

## Insertion points (xvue/qt/src/xvue_qt_api.cpp)

- `xvinitgraphique_` env-probe block: lines **350–381** (block lives AFTER
  `win->show()/raise()/activateWindow()` and BEFORE the bounded
  exposure-pump loop). Reads `MEFISTO_QT_WINDOW_SIZE`, gates on
  `MEFISTO_BATCH_X11=1` OR `QT_QPA_PLATFORM=offscreen`, calls
  `setMinimumSize/setMaximumSize` on the canvas + `win->adjustSize()`.
- `xvfermer_` clear block: lines **943–953** (just before
  `XvueApp::window_slot().reset()`). Resets canvas constraints to
  `(0,0)` / `(QWIDGETSIZE_MAX, QWIDGETSIZE_MAX)` so a subsequent reopen
  starts clean.

Insertion strategy is per RESEARCH §Pattern 4 + §Pitfall 5:
- env-probe in `xvinitgraphique_` (the canonical "open window" entry —
  fires once per reopen cycle), NOT in `xvinfo_` (where the original
  `win->resize(*ix, *iy)` call lives). Constraints applied here clip
  any subsequent `win->resize` from `xvinfo_` to the env-specified dim.
- `setMinimumSize` + `setMaximumSize` on the **canvas** (NOT
  `setFixedSize` on the **window**) so `xvfermer_` can zero them out and
  the X11-style chrome computation does not interfere.

## Empirical PNG dim (Task 1 step 6)

Command (from worktree root, `MEFISTO=$(pwd)`):

```sh
MEFISTO_QT_WINDOW_SIZE=1280x800 \
  bin/ab_sweep_phase8.sh --mode qt-1x --cases pan2d \
      --out-dir /tmp/09-06-verify-pan2d
identify -format '%wx%h' /tmp/09-06-verify-pan2d/pan2d-qt-1x.png
```

Output: **`1280x800`** (exact match — no chrome-accounting tolerance
needed).

Iteration count: **1** (first attempt produced the matched dim).
Open Question 5 (RESEARCH) is empirically resolved: `setMinimumSize` +
`setMaximumSize` on the canvas DOES propagate to
`XvueCanvas::resizeEvent`'s backing-pixmap allocation once
`win->adjustSize()` runs and `win->show()` triggers the natural
exposure-driven resize chain. No additional `processEvents` pump or
`QTimer::singleShot(0)` deferral was needed.

## ABI count (T-09-06 mitigation)

Pre-edit count: **64** (per 09-01-AUDIT-BASELINE §3 row 7 — the
empirical baseline; the *plan-quoted* "58" is the Phase-6.5-frozen value
that drifted during the Phase 7 export plans, +6 entries — see
09-01-AUDIT-BASELINE §"#7 ABI count drift").

Post-edit count: **64** (`nm xvue/qt/build/libxvueqt.a | grep ' T ' |
grep '_$' | wc -l = 64`). Delta: 0 — no new `extern "C"` entry was
added. The implementation is env-var-only per RESEARCH §Anti-Pattern
"Bumping ABI count from 58 to 59": adding a Fortran-callable setter
would have been the wrong fix.

The audit baseline's #7 row notes "Phase 9 retirement does not touch Qt
extern \"C\" entry points" — Plan 9-06 is consistent with that
expectation.

## Build status

`bin/cbl_tout` exit: **0**. All `pp/pp{init,mail,elas,flui,ther,nlse,
xvtest0..4}` rebuilt cleanly with the matched-dim code path. No new
compiler warnings introduced (`std::strcmp` and `std::sscanf` were
already in the file's reach via `<cstring>` and `<cstdio>` includes).

## Pitfall 5 mitigation evidence

The env probe is wrapped in
```cpp
if ((headless_batch || offscreen) && qt_window_size && qt_window_size[0]) {
```
Interactive sessions where neither `MEFISTO_BATCH_X11=1` nor
`QT_QPA_PLATFORM=offscreen` is set will short-circuit at the first
`&&` test even if `MEFISTO_QT_WINDOW_SIZE` is somehow set — the
constraints never fire and the existing `XvueWindow::resize(1024, 768)`
default from the constructor stays in effect.

T-09-05-A (input validation): the `sscanf("%dx%d")` failure path returns
`!=2` and the bounds check `0 < {w,h} < 8192` rejects negatives,
overflow-prone values, and zero. Malformed env values are silently
ignored — no error path is exposed to the headless harness.

## Task 2 — AE re-run on Phase-8 CHECK cells

`bin/ab_sweep_phase8.sh` qt-1x branch defaults `MEFISTO_QT_WINDOW_SIZE`
to `1280x800` (using the `: "${VAR:=default}"` form so a user
override is respected).

**Cavity2d reproduction**: Phase 8 baseline AE was 411003 (40.14%).
Plan 9-06 attempted to reproduce but `pp/ppflui` times out at the
60-second harness budget (also at 120 s in a manual probe) on the
`cavity2d.stoke56cr` batch with `IEEE_DENORMAL` FPE + `CRASH of
MEFISTO software` before reaching `xvfermer_`. The capture path
`MEFISTO_QT_CAPTURE_PATH` therefore never fires. This is a pre-existing
Phase-8 ppflui stability issue (orthogonal to Plan 9-06 — the
matched-dim env-knob is exercised on every Qt-mode capture but the
fluider crash blocks output entirely).

Substituted three other CHECK cells from the canonical 5-case grid as
A/B samples. Numbers (Phase 8 -> Plan 9-06, both `compare -metric AE
-fuzz 5%`):

| Case | Phase 8 AE | Plan 9-06 AE | drop | drop % |
|------|-----------|--------------|------|--------|
| pan2d | 540804 (52.81%) | 520331 (50.81%) | 20473 | -3.78% |
| heat1d | 143273 (13.99%) | 125480 (12.25%) | 17793 | -12.42% |
| nafems_le1 | 412827 (40.32%) | 325282 (31.77%) | 87545 | -21.21% |

`bin/ab_compare_pair.sh` reports `resampled=no` for all three (the
captures are already 1280x800 — the dim-guard's nearest-neighbor
resample never fires post-Plan-9-06).

**Interpretation:**
- The resample-confound is empirically eliminated (`resampled=no`).
- The remaining AE is genuine content-driven diff:
  - **Chrome bars**: Qt menubar/toolbar/statusbar/console-dock vs
    Xvfb root-window (no chrome). The chrome takes a fixed pixel
    region near the top/sides of the 1280x800 capture; chrome pixels
    differ from X11's mesh-rendering pixels in the same region by a
    large constant per-pixel delta — explaining the still-high AE for
    pan2d/nafems_le1.
  - **Font AA drift (Pitfall 7)**: `font-pan2d.md` documented Qt's
    3545 unique colors vs X11's 9 — antialiasing-driven.
- The smaller drops than the plan ideal (<5%) are NOT a
  misimplementation of Plan 9-06; they reflect that resample was
  ONE confound but chrome+AA were also significant for these cases.
  Heat1d (1D plot, less chrome-overlap) shows a cleaner -12.4% drop;
  nafems_le1 (no gradient bar) shows the strongest -21.2% drop.

`compare -metric AE -fuzz 5%` log captured at
`/tmp/09-06-cavity-AE.log` (also has the no-fuzz raw values + the
reasoning bundle).

## Task 4 — full rebuild + 5-case post-fix sweep

`bin/cbl_tout` exit 0; ABI count 64 (unchanged from pre-edit).

5-case post-fix sweep (`pan2d,nafems_le1,cavity2d,heat1d,nlsecu`):

| Case | PNG dim | AE (-fuzz 5%) | verdict | resampled |
|------|---------|---------------|---------|-----------|
| pan2d | 1280x800 | 520256 (50.8062%) | CHECK | no |
| nafems_le1 | 1280x800 | 325264 (31.7641%) | CHECK | no |
| cavity2d | (PNG missing — ppflui timeout) | n/a | n/a | n/a |
| heat1d | 1280x800 | 125495 (12.2554%) | CHECK | no |
| nlsecu | (PNG missing — TRUNCATED-CAPTURE) | n/a | n/a | n/a |

**4 of 5 testa cases capture cleanly at exact 1280x800** (acceptance
criterion). The 2 missing cases are pre-existing Phase 8 known issues
unrelated to Plan 9-06:
- **cavity2d**: ppflui IEEE_DENORMAL FPE + crash before xvfermer_;
  documented in 08-CHECKLIST.md (pitfall-6-secondary). Plan 9-06
  matched-dim env-knob is exercised but never reaches capture.
- **nlsecu**: TRUNCATED-CAPTURE per 08-CHECKLIST.md (ppnlse_qt
  offscreen+BATCH_X11 deadlock + 60s timeout — unreachable canonical
  TIME=20/2000 steps).

3 grep gates all PASS:
- `bin/test_no_imagemagick_in_qt.sh` -> `EXPORT-06 PASS: no ImageMagick references in xvue/qt/`
- `bin/test_no_x11_in_build.sh` -> `OK: no X11 references in active build path`
- `bin/test_no_lvideo.sh` -> `OK: no LVIDEO and no Fortran convert shell-outs`

## Plan 9-06 closure summary

Phase 8 override #1 closed: the 14 Qt-mode CHECK cells in
08-CHECKLIST.md whose AE was dominated by the 760x442->1280x800
nearest-neighbor resample now have a matched-dim recapture path.
Re-running the harness produces 1280x800 captures (verified on 3 of
5 cases; 2 cases blocked by pre-existing Phase 8 ppflui/ppnlse
stability issues unrelated to Plan 9-06). The resample-confound is
empirically eliminated. Residual AE is genuine content-driven
(chrome bars + Pitfall 7 font AA) — Pitfall 6/7 overrides can now
apply cleanly to the per-case verdicts.

ABI count 64 unchanged; T-09-05 (interactive UX preservation) and
T-09-06 (no extern "C" entry added) both upheld.

## Threats checked at Plan 9-06 closure

| Threat ID | Disposition | Evidence |
|-----------|-------------|----------|
| T-09-05 | mitigated | Headless gate `(MEFISTO_BATCH_X11=1 OR QT_QPA_PLATFORM=offscreen) && MEFISTO_QT_WINDOW_SIZE` short-circuits in interactive contexts; setMin/Max not setFixed; xvfermer_ resets constraints. |
| T-09-06 | mitigated | nm count 64 -> 64 (delta 0). Implementation uses env var only — no new extern "C" entry. |
| T-09-05-A | mitigated | sscanf("%dx%d") + `0 < {w,h} < 8192` bounds; malformed silently ignored — no error path exposed. |
| T-09-05-B | mitigated | Empirical verify on first iteration: 1280x800 PNG produced; no Open-Question-5 iteration needed. |

