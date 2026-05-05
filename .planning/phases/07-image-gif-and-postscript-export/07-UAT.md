---
status: partial
phase: 07-image-gif-and-postscript-export
source: [07-01-SUMMARY.md, 07-02-SUMMARY.md, 07-03-SUMMARY.md, 07-04-SUMMARY.md, 07-05-SUMMARY.md, 07-06-SUMMARY.md]
started: 2026-05-05T00:00:00Z
updated: 2026-05-05T00:30:00Z
---

## Current Test

[testing paused — 3 items outstanding (human-bootstrap goldens)]

## Tests

### 1. Cold Start Smoke Test
expected: From clean tree, `bin/cbl_tout` exits 0 and produces 13 X11 executables in pp/. `cmake --build xvue/qt/build` exits 0 and produces libxvueqt.a + pp*_qt binaries.
result: pass

### 2. ABI Symbol Count Invariant
expected: `bash xvue/qt/cmake/verify_abi.sh xvue/qt/build/libxvueqt.a xvue/qt/include/xvue_qt_api.h` exits 0 reporting "nm count: 58 header count: 58".
result: pass

### 3. EXPORT-06 Grep Gate
expected: `bash bin/test_no_imagemagick_in_qt.sh` exits 0 with message "EXPORT-06 PASS: no ImageMagick references in xvue/qt/". CMake build also runs verify_no_imagemagick_in_qt ALL target green.
result: pass

### 4. Automated Qt Test Suite
expected: `cd xvue/qt/build && xvfb-run --auto-servernum ctest -R '^xvue_qt_(postscript|export)_tests$' --output-on-failure` reports xvue_qt_postscript_tests 18 PASS + 1 SKIP and xvue_qt_export_tests 16 PASS + 2 SKIP. The 3 SKIPs are by-design golden-file-gated.
result: pass

### 5. X11 Backend Non-Regression
expected: One small testa case (e.g. testa/poutre or any quick mesher case) runs end-to-end through the X11 backend (INITIER → MAILLER → save 99;) without crashes or visual regressions versus pre-Phase-7 behaviour. xvuelc.c is byte-identical (`git diff 900e297..HEAD -- xvue/xvuelc.c` empty).
result: pass

### 6. File → Export Submenu Layout
expected: Launching any pp*_qt binary (e.g. ppmail_qt or ppinit_qt) shows a `File → Export` submenu containing four children in order: PNG…, JPEG…, PDF…, GIF…. Bilingual labels (FR/EN) follow xvueIsEnglish().
result: pass

### 7. File → Capture Animation Toggle
expected: Same Qt window's File menu shows a `Capture Animation` checkable entry directly under File (NOT inside Export submenu). Toggling it on starts an animation capture; toggling off ends it.
result: pass

### 8. PNG Export Round-Trip
expected: File → Export → PNG…, choose path, save. Resulting .png opens in any image viewer (eog, feh, gimp) and shows the current canvas content. Re-opening QFileDialog defaults to the last-used directory (QSettings persistence).
result: pass

### 9. JPEG Export Round-Trip
expected: File → Export → JPEG…, save to .jpg. Opens in viewer, shows canvas content (lossy compression acceptable).
result: pass

### 10. PDF Export at Canvas-Native Geometry
expected: File → Export → PDF…, save to .pdf. Opens in evince/okular at 72 dpi with page aspect matching the canvas (no fit-to-A4 distortion). PDF is raster (drawPixmap of backing_).
result: pass

### 11. Animated GIF End-to-End
expected: With `Capture Animation` toggled on, run a Fortran case that emits multiple xvpostscript_(0) closes (e.g. testa/wave or any animated solver). Toggle off, then File → Export → GIF…, save. ffmpeg produces a playable animated .gif viewable in any browser/viewer; frame count > 1.
result: pass

### 12. XVUE_ANIM=1 Env-Var Auto-Capture
expected: `XVUE_ANIM=1 ppflui_qt …` (or any pp*_qt) auto-starts capture at process boot. Running an animated case writes frames; on graceful shutdown an animation.gif is produced (or the GIF menu item dispatches with frames already buffered).
result: pass

### 13. TEMPORAIRE.EPS via xvpostscript_
expected: From any pp*_qt running session, triggering the legacy "save PostScript" Fortran call (xvpostscript_ entry) writes a TEMPORAIRE.EPS in cwd. File opens in evince/okular and renders the canvas content. Format strings match the legacy xvuelc.c per-primitive output.
result: pass

### 14. scene01.eps Golden Bootstrap (human)
expected: Compile xvue/qt/tests/golden/scene01_driver.f against X11 + xvuelc.o, run under Xvfb, copy TEMPORAIRE.EPS → xvue/qt/tests/golden/scene01.eps, commit. Re-run ctest -R xvue_qt_postscript_tests: Test 15 (PsEmitter_postscriptVerbatim_golden) flips QSKIP → PASS.
result: blocked
blocked_by: prior-phase
reason: "Pre-flagged in VERIFICATION.md §9 as human-bootstrap. Requires X11/Xvfb + manual gfortran link + scene01.eps commit. Deferred to Phase 8 A/B checkpoint per VERIFICATION.md verdict."

### 15. wave_legacy.gif Visual A/B (human)
expected: Run testa/wave through X11 + bin/convertepsgif → commit xvue/qt/tests/golden/wave_legacy.gif. Re-run with XVUE_ANIM=1 + Qt backend, eyeball-compare animation.gif vs baseline. Same mesh-evolution sequence; frame count matches ±1; no color/geometry drift.
result: blocked
blocked_by: prior-phase
reason: "Pre-flagged in VERIFICATION.md §9 as human-bootstrap. Visual GIF A/B inherently subjective; requires X11+Qt side-by-side session. Deferred to Phase 8 A/B checkpoint."

### 16. cavity2d_legacy.gif Visual A/B (human)
expected: Same procedure as wave for testa/cavity2d → commit xvue/qt/tests/golden/cavity2d_legacy.gif. Eyeball A/B Qt vs X11. Same parity criteria.
result: blocked
blocked_by: prior-phase
reason: "Pre-flagged in VERIFICATION.md §9 as human-bootstrap. Same rationale as wave. Deferred to Phase 8 A/B checkpoint."

## Summary

total: 16
passed: 13
issues: 0
pending: 0
skipped: 0
blocked: 3

## Gaps

[none — 13 pass, 3 blocked on pre-flagged human-bootstrap goldens (VERIFICATION.md §9). Not code issues; explicit prerequisite gates merged into Phase 8.]
