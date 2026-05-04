# Phase 7 — Validation Log

This file records the manual A/B sign-off verdict for the Phase 7 export
deliverables. One line per case + ISO date + verdict (PASS / FAIL <reason>
/ DEFERRED <reason>).

## Manual A/B sign-off (Plan 06 Task 3)

| Date (UTC) | Case | Verdict | Notes |
|------------|------|---------|-------|
| 2026-05-04 | testa/wave A/B | DEFERRED | Plan 06 was executed by an autonomous agent without an attached human terminal. The bootstrap procedure (compile scene01_driver.f against the X11 backend, run testa/wave under both X11 and Qt+XVUE_ANIM=1, eyeball-compare the two GIFs) requires a human verdict and cannot be performed headlessly. The byte-level golden + GIF A/B test slots in xvue/qt/tests/test_xvue_qt_postscript.cpp + test_xvue_qt_export.cpp QSKIP cleanly until the goldens land; once the human runs Task 3 and commits the goldens, those slots flip from QSKIP to PASS-required. Surfaced as a phase-verification gap for the orchestrator. |
| 2026-05-04 | testa/cavity2d A/B | DEFERRED | Same rationale as testa/wave above. |
| 2026-05-04 | scene01.eps byte-parity | DEFERRED | The byte-level golden requires running scene01_driver.f against the X11+xvuelc.c backend to materialize TEMPORAIRE.EPS as the reference. The Fortran driver source is committed at xvue/qt/tests/golden/scene01_driver.f with a header BOOTSTRAP NOTE documenting the procedure; the .eps golden file itself materializes when the human runs Task 3. |

## Automated verifications performed during Plan 06 (autonomous portion)

| Date (UTC) | Gate | Result |
|------------|------|--------|
| 2026-05-04 | EXPORT-06 grep gate (`bash bin/test_no_imagemagick_in_qt.sh`) | PASS — exit 0, "EXPORT-06 PASS: no ImageMagick references in xvue/qt/" |
| 2026-05-04 | EXPORT-06 CMake target (`verify_no_imagemagick_in_qt`) | PASS — built as part of `cmake --build xvue/qt/build`, target green |
| 2026-05-04 | ABI symbol count (`xvue/qt/cmake/verify_abi.sh`) | PASS — `nm count: 58 header count: 58` (unchanged) |
| 2026-05-04 | xvue/xvuelc.c byte-identical (T-07-08) | PASS — `git diff HEAD~5 HEAD -- xvue/xvuelc.c` is empty |
| 2026-05-04 | LVIDEO untouched (T-07-06, D-17) | PASS — `git diff --name-only HEAD~5 HEAD -- xvue/video1.f xvue/videofin.f xvue/videonm.f` is empty |
| 2026-05-04 | X11 full build (`bin/cbl_tout`) | PASS — exit 0; all 13 `pp/pp*` X11 binaries produced (ppelas / ppflui / ppinit / ppmail / ppnlse / pppoba / ppther / ppxvtest0..4 / pxyz) |
| 2026-05-04 | Qt full build (`cmake --build xvue/qt/build`) | PASS — exit 0; libxvueqt.a + 14 test binaries linked; all 4 verify_* gates green |
| 2026-05-04 | xvue_qt_postscript_tests (ctest) | PASS — 18 PASS + 1 SKIP (PsEmitter_postscriptVerbatim_golden — flips to PASS post-Task 3) |
| 2026-05-04 | xvue_qt_export_tests (ctest) | PASS — 16 PASS + 2 SKIP (XvueExport_gif_AB_compare_{wave,cavity2d} — flips to PASS post-Task 3) |
| 2026-05-04 | All other Qt test binaries | PASS at pre-Plan-06 rate; one pre-existing failure (`testPerModuleGroupIsolation` in `xvue_qt_i18n_menu_prefs_tests`) documented in Plan 02/03/05 SUMMARYs as out of scope per SCOPE BOUNDARY rule |

## Procedure to elevate DEFERRED → PASS

When a human is available to run the manual A/B sign-off, follow Plan 06
Task 3's `<how-to-verify>` block (steps 1-9). The headless executor has
already shipped the test harness (Task 1) and the EXPORT-06 grep gate
(Task 2); the remaining steps are:

1. Compile `scene01_driver.f` against the X11 backend (gfortran +
   `xvue/xvuelc.o` + `-lX11 -lXt`).
2. Run the resulting binary on Xvfb to produce `TEMPORAIRE.EPS`; commit
   it as `xvue/qt/tests/golden/scene01.eps`.
3. Run `testa/wave` and `testa/cavity2d` through the X11 backend +
   `bin/convertepsgif` to produce baseline GIFs; commit them as
   `xvue/qt/tests/golden/wave_legacy.gif` and `cavity2d_legacy.gif`.
4. Run the same testa cases through the Qt backend with `XVUE_ANIM=1`,
   eyeball-compare the resulting `animation.gif` against the legacy
   baselines.
5. Replace the DEFERRED rows above with PASS/FAIL verdicts.
6. Re-run `ctest -R '^xvue_qt_(postscript|export)_tests$'` and confirm
   the golden-dependent slots flip from QSKIP to PASS.

## Notes on test SKIP semantics

The QSKIP messages in the Phase 7 unit tests are explicit-by-design.
They do NOT silently mask gaps — every SKIP carries a pointer to the
exact bootstrap procedure (Plan 06 Task 3) that elevates the slot to
PASS-required. A future executor running `ctest` will see the SKIP
banner and know the gate is incomplete.
