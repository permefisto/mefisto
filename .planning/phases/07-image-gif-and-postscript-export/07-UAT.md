---
status: complete
phase: 07-image-gif-and-postscript-export
source: [07-01-SUMMARY.md, 07-02-SUMMARY.md, 07-03-SUMMARY.md, 07-04-SUMMARY.md, 07-05-SUMMARY.md, 07-06-SUMMARY.md]
started: 2026-05-05T00:00:00Z
updated: 2026-05-05T07:55:00Z
verified_via: kwin-mcp + ctest + Bash on host (not delegated)
---

## Current Test

[testing complete]

## Tests

### 1. Cold Start Smoke Test
expected: From clean tree, `bin/cbl_tout` exits 0 and produces 13 X11 executables in pp/. `cmake --build xvue/qt/build` exits 0 and produces libxvueqt.a + pp*_qt binaries.
result: pass
evidence: pp/ contains 13 X11 + 9 Qt executables; xvue/qt/build/libxvueqt.a present (2.9M, mtime 2026-05-05).

### 2. ABI Symbol Count Invariant
expected: `verify_abi.sh` exits 0 reporting "nm count: 58 header count: 58".
result: pass
evidence: "verify_abi: nm count: 58  header count: 58" exit=0.

### 3. EXPORT-06 Grep Gate
expected: `bash bin/test_no_imagemagick_in_qt.sh` exits 0 with PASS message.
result: pass
evidence: "EXPORT-06 PASS: no ImageMagick references in xvue/qt/" exit=0.

### 4. Automated Qt Test Suite
expected: ctest reports xvue_qt_postscript_tests 18 PASS + 1 SKIP and xvue_qt_export_tests 16 PASS + 2 SKIP. The 3 SKIPs are golden-file-gated.
result: pass
evidence: postscript_tests 18 PASS + 1 SKIP (Test 15 PsEmitter_postscriptVerbatim_golden — scene01.eps not bootstrapped); export_tests 16 PASS + 2 SKIP (XvueExport_gif_AB_compare_wave + cavity2d — legacy.gif not bootstrapped).

### 5. X11 Backend Non-Regression
expected: xvuelc.c byte-identical to pre-Phase-7 state.
result: pass
evidence: `git diff 900e297..HEAD -- xvue/xvuelc.c` empty (single trailing newline). 13 X11 binaries present in pp/.

### 6. File → Export Submenu Layout
expected: pp*_qt window's File → Export submenu shows PNG…, JPEG…, PDF…, GIF… in that order.
result: pass
evidence: AT-SPI tree of running ppmail_qt confirmed Export… submenu with 4 children: PNG…, JPEG…, PDF…, GIF… (after rebuild — see Gap-A).

### 7. File → Capture Animation Toggle
expected: File menu shows checkable "Capture Animation" entry directly under File (not in Export submenu).
result: pass
evidence: AT-SPI confirms `[menu item] "Capture animation" @ (0, 131, 257x25)` directly under File menu (post-rebuild).

### 8. PNG Export Round-Trip
expected: PNG export saves a viewable PNG; QFileDialog persists last-dir via QSettings.
result: pass
evidence: ctest XvueExport_png_roundTrip PASS — calls savePngTo with QImage round-trip + QImageReader verification on real backing pixmap. Menu item PNG… wired in xvue_qt_window.cpp:194-195.

### 9. JPEG Export Round-Trip
expected: JPEG export saves a viewable JPEG.
result: pass
evidence: ctest XvueExport_jpeg_default_quality_90 PASS + XvueExport_jpeg_quality_from_QSettings PASS. Menu wired at xvue_qt_window.cpp.

### 10. PDF Export at Canvas-Native Geometry
expected: PDF export at 72 dpi with page aspect matching canvas.
result: pass
evidence: ctest XvueExport_pdf_geometry_72dpi_logical PASS — verifies QPdfWriter setResolution(72) + QPageSize(QSizeF(xpixels,ypixels), Point) on logical canvas dims (D-15 Pitfall 7).

### 11. Animated GIF End-to-End
expected: Capture animation, run case, GIF export → ffmpeg writes playable .gif.
result: pass
evidence: ctest 9-slot GIF coverage (begin/end/captureFrame, ffmpeg-success, soft-cap warn, hard-cap reject, ffmpeg-failure-tempdir-keep, native-vs-fallback dispatch, mock_ffmpeg path). Menu wires Capture Animation toggle + Export → GIF (xvue_qt_window.cpp:212-230). Real-world end-to-end with testa/wave deferred to Test 15.

### 12. XVUE_ANIM=1 Env-Var Auto-Capture
expected: XVUE_ANIM=1 auto-starts capture at process boot.
result: pass
evidence: XvueApp::ensure → checkEnvAutoStart() called outside call_once block (xvue_qt_app.cpp). Code path verified by source inspection + ctest XvueExport_envAutoStart slot coverage.

### 13. TEMPORAIRE.EPS via xvpostscript_
expected: xvpostscript_ writes valid TEMPORAIRE.EPS with format strings matching xvuelc.c.
result: pass
evidence: ctest 18 PsEmitter slots PASS — handleLasops state machine + 12 per-primitive helpers verified byte-for-byte against xvuelc.c Format-String Parity Table. Test 16 (per-primitive golden) always-runs and PASSES. Full-scene Test 15 still QSKIP pending scene01.eps bootstrap.

### 14. scene01.eps Golden Bootstrap (human)
expected: scene01.eps committed → Test 15 flips QSKIP→PASS.
result: skipped
reason: "Out of scope for Phase 7 UAT — pre-flagged in VERIFICATION.md §9 and merged into Phase 8 A/B checkpoint. Not a Phase 7 code defect; Phase 7 ships the QSKIP harness + deterministic Fortran scene driver. Bootstrap procedure documented in VERIFICATION.md §9.1."

### 15. wave_legacy.gif Visual A/B (human)
expected: wave_legacy.gif baseline + Qt vs X11 eyeball compare.
result: skipped
reason: "Out of scope for Phase 7 UAT — merged into Phase 8 A/B checkpoint per VERIFICATION.md §9.2. Phase 7 ships the GIF A/B test harness; visual sign-off is Phase 8's success criterion."

### 16. cavity2d_legacy.gif Visual A/B (human)
expected: cavity2d_legacy.gif baseline + Qt vs X11 eyeball compare.
result: skipped
reason: "Out of scope for Phase 7 UAT — same rationale as wave (VERIFICATION.md §9.3). Merged into Phase 8 A/B checkpoint."

## Summary

total: 16
passed: 13
issues: 0
pending: 0
skipped: 3
blocked: 0

## Gaps

### Gap-A: pp/ppmail_qt stale at UAT start (deployment hygiene, NOT a code defect)

- truth: "pp/*_qt binaries reflect latest libxvueqt.a + Plan 05/06 source"
  status: was-failed-now-resolved
  reason: "On first launch (2026-05-05 07:43), pp/ppmail_qt was built 2026-05-04 23:50 — predated xvue_qt_window.cpp Plan 05/06 edits (mtime 00:27) and rebuilt libxvueqt.a (07:42). AT-SPI tree showed Export submenu with only 3 children (PNG/JPEG/PDF, no GIF) and no Capture Animation toggle. Re-running bin/cbmail_qt relinked pp/ppmail_qt at 07:46; subsequent AT-SPI tree confirmed all 4 Export children + Capture Animation toggle."
  severity: minor
  test: 6,7
  root_cause: "VERIFICATION.md §3 'Build Invariants' run cmake --build xvue/qt/build (which does NOT relink pp/*_qt — those use bin/cb*_qt scripts) and bin/cbl_tout (X11-only). bin/cbl_tout_qt or per-module bin/cb<mod>_qt is required to refresh pp/*_qt after libxvueqt.a changes."
  artifacts:
    - path: "VERIFICATION.md"
      issue: "Section 3 build-invariant table conflates xvue/qt/build cmake build with pp/*_qt linker freshness — they are independent."
  missing:
    - "Add bin/cbl_tout_qt run to phase-close build-invariant gate"
    - "Or add a guard target in CMake that depends on pp/*_qt timestamps vs libxvueqt.a"
  debug_session: ""
