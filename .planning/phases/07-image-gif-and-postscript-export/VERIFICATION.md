---
phase: 07-image-gif-and-postscript-export
verified: 2026-05-04T00:00:00Z
status: human_needed
score: 4/5 must-haves verified (SC-2 GIF visual A/B is human-gated)
overrides_applied: 0
gaps: []
deferred: []
human_verification:
  - test: "Compile xvue/qt/tests/golden/scene01_driver.f against the X11 backend on Xvfb, run to produce TEMPORAIRE.EPS, commit as xvue/qt/tests/golden/scene01.eps; then re-run ctest -R xvue_qt_postscript_tests and confirm PsEmitter_postscriptVerbatim_golden (Test 15) flips from QSKIP to PASS."
    expected: "diff scene01.eps (Qt) vs X11 TEMPORAIRE.EPS shows byte-for-byte (or near byte-for-byte) parity on the deterministic scene01_driver.f primitive sequence."
    why_human: "Requires X11 display (Xvfb or real) to materialise TEMPORAIRE.EPS from the legacy backend. Autonomous executor has no attached terminal. CLAUDE.md §Tests explicitly delegates large/visual tests to the user."
  - test: "Run testa/wave through X11 backend + bin/convertepsgif to produce wave_legacy.gif; commit as xvue/qt/tests/golden/wave_legacy.gif. Then run with XVUE_ANIM=1 and Qt backend and eyeball-compare the resulting animation.gif against the legacy baseline."
    expected: "Both GIFs show the same mesh-evolution sequence; frame count matches (or is within ±1 due to open/close semantics)."
    why_human: "Visual GIF comparison is inherently subjective. Requires X11+Qt side-by-side session. Autonomous executor cannot render graphics."
  - test: "Run testa/cavity2d through X11 backend + bin/convertepsgif to produce cavity2d_legacy.gif; commit as xvue/qt/tests/golden/cavity2d_legacy.gif. Then run with XVUE_ANIM=1 and Qt backend and eyeball-compare."
    expected: "Both GIFs show the same cavity mesh-evolution sequence."
    why_human: "Same rationale as the wave case above."
---

# Phase 7: Image, GIF, and PostScript Export — Verification Report

**Phase Goal:** Qt-native PNG/JPEG/PDF plus a preserved-verbatim PostScript emitter and a runtime-probed animated GIF path fully replace the legacy `xvuelc.c` export code and `bin/convertepsgif` ImageMagick shell-out.
**Verified:** 2026-05-04
**Status:** PASS-WITH-GAPS (human verification required)
**Re-verification:** No — initial verification.

---

## Verdict: PASS-WITH-GAPS

All automated gates pass. Three human-required items remain open (the manual A/B visual sign-off and the three golden binary artifacts). These are explicitly pre-planned as `checkpoint:human-verify` work items documented in VALIDATION-LOG.md. No automated must-have is FAILED. No BLOCKER gaps exist. The phase may proceed to Phase 8 provided the human A/B sign-off and golden bootstrap are completed before Phase 8's own A/B gate.

---

## 1. Success Criteria Coverage (ROADMAP.md)

| # | ROADMAP Success Criterion | Status | Evidence |
|---|--------------------------|--------|----------|
| SC-1 | PROBE.md records `QImageWriter::supportedImageFormats()` output; animated-GIF strategy documented and implemented (QImageWriter loop or per-frame PNG + ffmpeg — never ImageMagick) | COVERED | `.planning/phases/07-image-gif-and-postscript-export/PROBE.md` exists, 32 lines, records `qt_version=6.10.2 gif_write_supported=0 ffmpeg version 8.1-3+b1`; D-11 ffmpeg path implemented in `xvue/qt/src/xvue_qt_export.cpp` via `QProcess::execute("ffmpeg", args)`. Path note: ROADMAP text says `.planning/phase-7/PROBE.md` — that directory does not exist; the actual path `.planning/phases/07-image-gif-and-postscript-export/PROBE.md` is a naming drift in the ROADMAP prose, not an implementation gap. |
| SC-2 | PNG and JPEG export via XvueExport opens correctly in image viewer; testa/wave + testa/cavity2d animated GIFs visually match legacy bin/convertepsgif pipeline | PARTIAL | PNG/JPEG: `XvueExport::savePngTo` / `saveJpegTo` implemented (xvue_qt_export.cpp:127-180), 5/5 unit tests pass including round-trip. GIF A/B: automated test harness shipped (test_xvue_qt_export.cpp:507-545); 2 slots QSKIP cleanly because `wave_legacy.gif` + `cavity2d_legacy.gif` golden baselines not yet committed — requires human (see Section 7). |
| SC-3 | `xvpostscript_` implemented in `xvue/xvue_qt_postscript.cpp` by moving the ~120-line fprintf emitter verbatim from xvuelc.c; .ps export byte-for-byte parity | COVERED (automated half) | `xvue/qt/src/xvue_qt_postscript.cpp` (1111 lines) contains `PsEmitter::handleLasops` — verbatim port of xvuelc.c:1187-1304 with three documented divergences (fopen-failure QMessageBox, chaine[-1] UB guard, post-free nullptr). Path note: ROADMAP says `xvue/xvue_qt_postscript.cpp`; actual path is `xvue/qt/src/xvue_qt_postscript.cpp` — naming drift in ROADMAP, not an implementation gap. Byte-parity Test 16 (per-primitive, always runs) passes; Test 15 (full-scene byte compare vs scene01.eps) QSKIPs until human commits the .eps golden. |
| SC-4 | PDF via QPrinter::PdfFormat as additive bonus without modifying xvpostscript_ | COVERED | `XvueExport::savePdfTo` uses `QPdfWriter` (the Qt 6 successor to `QPrinter::PdfFormat`) at 72 dpi canvas-native geometry (xvue_qt_export.cpp:190-230). `xvpostscript_` / `PsEmitter::handleLasops` is unmodified. ROADMAP used `QPrinter::PdfFormat` as planning shorthand; `QPdfWriter` is the correct Qt 6 API and satisfies the same intent. |
| SC-5 | No Qt-backend code path invokes ImageMagick `convert` — verified by grep | COVERED | `bash bin/test_no_imagemagick_in_qt.sh` exits 0: "EXPORT-06 PASS: no ImageMagick references in xvue/qt/". CMake `verify_no_imagemagick_in_qt ALL` target also green. Note: ROADMAP SC-5 says `grep -rn 'convert' xvue/` — running that broader scope does find `xvue/video1.f` and `xvue/videofin.f` (the LVIDEO pipeline). This is intentional and explicitly documented as Phase 9 RETIRE-03 scope (CONTEXT.md D-16/D-17). The Phase 7 grep gate is correctly scoped to `xvue/qt/` only. |

---

## 2. EXPORT-NN Requirement Coverage

| Requirement | Description | Status | Evidence |
|-------------|-------------|--------|----------|
| EXPORT-01 | QImageWriter probe run at phase kickoff | SATISFIED | PROBE.md committed at `ed1349c`; probe binary at `xvue/qt/probes/qimagewriter_probe.cpp`; `bin/cb_probe_qt` wrapper script exists and is executable |
| EXPORT-02 | PNG and JPEG export via XvueExport | SATISFIED | `savePngTo` + `saveJpegTo` bodies implemented; 5/5 unit tests pass (xvue_qt_export_tests) |
| EXPORT-03 | Animated GIF — auto-snapshot path + visual A/B | PARTIAL | Automated: `beginAnimation` / `captureFrame` / `saveGifTo` (ffmpeg fallback) implemented; 9/9 GIF unit tests pass. Manual A/B: DEFERRED — requires human (VALIDATION-LOG.md, 2 QSKIP slots) |
| EXPORT-04 | xvpostscript_ verbatim port from xvuelc.c | SATISFIED (automated) | `PsEmitter::handleLasops` at xvue_qt_postscript.cpp:80 is byte-for-byte port; 12 per-primitive helpers with verbatim xvuelc.c format strings; 17 unit tests pass + Test 15 QSKIP (golden-dependent) + Test 16 always-runs per-primitive parity passes |
| EXPORT-05 | PDF export as additive bonus via QPdfWriter | SATISFIED | `savePdfTo` uses `QPdfWriter::setResolution(72)` + canvas-native page geometry; 1/1 PDF geometry test passes |
| EXPORT-06 | No ImageMagick in xvue/qt/ | SATISFIED | `bash bin/test_no_imagemagick_in_qt.sh` exit 0; CMake ALL target green; word-boundary grep with Qt-API allowlist confirmed |

---

## 3. Threat-Model Adherence

| Threat | Description | Status | Evidence |
|--------|-------------|--------|----------|
| T-07-04 | No shell invocation for ffmpeg | ADHERED | `grep -n 'system(\|popen(\|QProcess::start.*sh' xvue_qt_export.cpp` returns 0 matches; ffmpeg invoked via `QProcess::execute("ffmpeg", QStringList{...})` with constants-only argv at export.cpp:490-502 |
| T-07-05 | Temp dir kept on failure, RAII-removed on success | ADHERED | `QTemporaryDir tempDir(...)` at export.cpp:440; `setAutoRemove(false)` on failure at export.cpp:505; Test 7 (`XvueExport_tempdir_kept_on_failure`) verifies the failure path |
| T-07-06 | LVIDEO scope: xvue/video1.f / videofin.f / videonm.f untouched | ADHERED | `git diff HEAD~10 HEAD -- xvue/video1.f xvue/videofin.f xvue/videonm.f` is empty; grep gate scope is `xvue/qt/` only (D-16/D-17); README Phase 7 section explicitly documents LVIDEO deferral to Phase 9 RETIRE-03 |
| T-07-07 | PsEmitter byte-output regression detection | ADHERED | Two-layer protection: Test 15 (full-scene, QSKIP until human commits scene01.eps) + Test 16 (per-primitive always-runs, 17 format-string assertions inline from Plan 03 Format-String Parity Table). `pyFlip(y)` called 35x inside xvue_qt_postscript.cpp; `grep 'ypixels.*-' xvue_qt_api.cpp` returns 0 |
| T-07-08 | xvue/xvuelc.c byte-identical (X11 backend untouched) | ADHERED | `git diff 900e297..HEAD -- xvue/xvuelc.c` is empty (0-line diff); last commit to xvuelc.c is `900e297` dated pre-Phase-7; X11 build `bin/cbl_tout` exit 0, 13 pp/* executables present |

---

## 4. Build Invariants

| Invariant | Expected | Actual | Status |
|-----------|----------|--------|--------|
| xvuelc.c byte-identical | No changes since Phase 7 start | `git diff 900e297..HEAD -- xvue/xvuelc.c` is empty | PASS |
| ABI symbol count | 58 | `verify_abi.sh` → "nm count: 58 header count: 58" exit 0 | PASS |
| X11 full build (`bin/cbl_tout`) | exit 0, 13 pp/* executables | VALIDATION-LOG.md records PASS on 2026-05-04; 13 executables in pp/ confirmed (ppelas ppflui ppinit ppmail ppnlse pppoba ppther ppxvtest0..4 pxyz) | PASS |
| Qt full build (`cmake --build xvue/qt/build`) | exit 0; libxvueqt.a linked | VALIDATION-LOG.md records PASS; libxvueqt.a present in xvue/qt/build/ | PASS |
| LVIDEO Fortran files untouched | 0 diff lines | `git diff HEAD~10 HEAD -- xvue/video*.f` empty | PASS |

---

## 5. Required Artifacts

| Artifact | Status | Details |
|----------|--------|---------|
| `xvue/qt/src/xvue_qt_postscript.h` | VERIFIED | Exists, 134 lines, declares PsEmitter class with all planned helper signatures |
| `xvue/qt/src/xvue_qt_postscript.cpp` | VERIFIED | Exists, 1111 lines — `handleLasops` verbatim port (lines 80-234) + 12 per-primitive helpers with documented format-string provenance |
| `xvue/qt/src/xvue_qt_export.h` | VERIFIED | Exists, 97 lines, declares all PNG/JPEG/PDF/GIF/animation surface methods |
| `xvue/qt/src/xvue_qt_export.cpp` | VERIFIED | Exists, 581 lines — `savePngTo` + `saveJpegTo` + `savePdfTo` + `saveGifTo` + animation state machine |
| `xvue/qt/src/xvue_qt_window.cpp` | VERIFIED | File→Export submenu with PNG/JPEG/PDF/GIF children wired at buildMenuBar():181-229; Capture Animation toggle at 219-229 |
| `xvue/qt/src/xvue_qt_api.cpp` | VERIFIED | `xvpostscript_` body at line 630-634 dispatches to `XvueApp::psEmitter().handleLasops(*lasops)`; 14 per-primitive psEmitter wiring sites confirmed |
| `.planning/phases/07-image-gif-and-postscript-export/PROBE.md` | VERIFIED | Exists, 32 lines, records empirical probe output (gif_write_supported=0, ffmpeg 8.1-3+b1) |
| `bin/test_no_imagemagick_in_qt.sh` | VERIFIED | Exists, 77 lines, executable; exit 0 confirmed by manual run |
| `xvue/qt/tests/golden/scene01_driver.f` | VERIFIED | Exists, 117 lines — deterministic Fortran scene driver with BOOTSTRAP NOTE |
| `xvue/qt/tests/golden/scene01.eps` | MISSING | Not yet committed — requires human Task 3 bootstrap (blocks Test 15 QSKIP flip) |
| `xvue/qt/tests/golden/wave_legacy.gif` | MISSING | Not yet committed — requires human Task 3 bootstrap (blocks GIF A/B Test slot flip) |
| `xvue/qt/tests/golden/cavity2d_legacy.gif` | MISSING | Not yet committed — requires human Task 3 bootstrap (blocks GIF A/B Test slot flip) |

---

## 6. Key Link Verification

| From | To | Via | Status |
|------|----|-----|--------|
| Fortran `CALL XVPOSTSCRIPT(lasops)` | `PsEmitter::handleLasops` | `proc(xvpostscript)` in xvue_qt_api.cpp:630-634 | WIRED |
| `PsEmitter::handleLasops(0)` | `XvueExport::captureFrame()` | Conditional hook at xvue_qt_postscript.cpp:136-138 | WIRED |
| `File → Export → PNG` menu action | `XvueExport::onMenuExportPng()` | `connect(actExportPng_, triggered, this, onFileExportPng)` in xvue_qt_window.cpp:194-195; slot calls `XvueExport::onMenuExportPng()` | WIRED |
| `File → Export → GIF` menu action | `XvueExport::onMenuExportGif()` | `connect(actExportGif_, triggered, this, onFileExportGif)` at xvue_qt_window.cpp:214-215 | WIRED |
| `XvueExport::saveGifTo` | ffmpeg subprocess | `QProcess::execute("ffmpeg", QStringList{...})` at export.cpp:499-502 | WIRED |
| EXPORT-06 grep gate | CMake build failure | `add_custom_target(verify_no_imagemagick_in_qt ALL ...)` at CMakeLists.txt:165-170 | WIRED |
| XvueCanvas resize | `PsEmitter::setCanvasDims` | `resizeEvent` in xvue_qt_canvas.cpp calls `psEmitter().setCanvasDims()` | WIRED (confirmed in 07-03-SUMMARY) |

---

## 7. Anti-Patterns Found

| File | Pattern | Severity | Assessment |
|------|---------|----------|------------|
| `xvue/qt/tests/test_xvue_qt_postscript.cpp` line 441 | `QSKIP("scene01.eps not yet bootstrapped...")` | INFO | By-design — QSKIP carries bootstrap-procedure pointer; flips to PASS-required when golden lands. Not a stub. |
| `xvue/qt/tests/test_xvue_qt_export.cpp` lines 517, 539 | `QSKIP("wave_legacy.gif / cavity2d_legacy.gif not yet bootstrapped...")` | INFO | By-design — same as above. Both carry inline Plan 06 Task 3 pointer. |

No blocker anti-patterns found. No placeholder returns, no empty handlers, no stub-only implementations in the production code paths.

---

## 8. Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| EXPORT-06: no ImageMagick in xvue/qt/ | `bash bin/test_no_imagemagick_in_qt.sh` | "EXPORT-06 PASS: no ImageMagick references in xvue/qt/" exit 0 | PASS |
| ABI symbol count = 58 | `bash xvue/qt/cmake/verify_abi.sh xvue/qt/build/libxvueqt.a xvue/qt/include/xvue_qt_api.h` | "nm count: 58 header count: 58" exit 0 | PASS |
| xvuelc.c byte-identical | `git diff 900e297..HEAD -- xvue/xvuelc.c` | empty diff (0 lines) | PASS |
| T-07-04: no shell invocation | `grep -nE '(system\(|popen\(|QProcess::start.*sh)' xvue/qt/src/xvue_qt_export.cpp` | 0 matches | PASS |
| LVIDEO untouched | `git diff HEAD~10 HEAD -- xvue/video1.f xvue/videofin.f xvue/videonm.f` | empty diff | PASS |
| Qt ctest (automated portion) | `ctest -R '^xvue_qt_(postscript\|export)_tests$'` | 18 PASS + 1 SKIP (postscript); 16 PASS + 2 SKIP (export) — SKIP slots are by-design golden-file-gated | PASS (with expected SKIPs) |

---

## 9. Human Verification Required

### 1. PostScript byte-parity golden — scene01.eps

**Test:** Compile `xvue/qt/tests/golden/scene01_driver.f` against the X11 backend (see BOOTSTRAP NOTE in the file header):

```bash
gfortran scene01_driver.f xvue/xvuelc.o -lX11 -lXt -o scene01_driver
Xvfb :99 &
DISPLAY=:99 ./scene01_driver
cp TEMPORAIRE.EPS xvue/qt/tests/golden/scene01.eps
git add xvue/qt/tests/golden/scene01.eps && git commit -m "test(07-06): commit scene01.eps PostScript golden"
cd xvue/qt/build && xvfb-run --auto-servernum ctest -R '^xvue_qt_postscript_tests$'
```

**Expected:** Test 15 (`PsEmitter_postscriptVerbatim_golden`) flips from QSKIP to PASS. The Qt-emitted TEMPORAIRE.EPS byte-sequence matches the X11 baseline within the three documented divergences (fopen-failure path, chaine[-1] guard, post-free null).

**Why human:** Requires an X11 display (Xvfb or real monitor). Autonomous executor has no attached terminal. Visual/byte-level judgment needed for any divergence outside the three documented ones.

### 2. GIF visual A/B — testa/wave

**Test:** Run `testa/wave` through the X11 backend + `bin/convertepsgif` to produce baseline GIF; commit as `xvue/qt/tests/golden/wave_legacy.gif`. Then run with `XVUE_ANIM=1` and the Qt backend; eyeball-compare `animation.gif` against the baseline. Record verdict in `VALIDATION-LOG.md`.

**Expected:** Both GIFs show the same mesh-evolution sequence. Frame count identical or within ±1. No color-band inversion, no missing frames, no geometry drift.

**Why human:** Animated GIF visual comparison cannot be automated. Requires side-by-side display and subjective judgment.

### 3. GIF visual A/B — testa/cavity2d

**Test:** Same procedure as testa/wave; commit as `xvue/qt/tests/golden/cavity2d_legacy.gif`.

**Expected:** Same parity criteria as the wave case.

**Why human:** Same rationale.

---

## 10. Verification Gaps Summary

### Open Human Items (3)

1. **scene01.eps PostScript golden** — blocks Test 15 (PsEmitter_postscriptVerbatim_golden). The per-primitive byte-parity (Test 16) already passes. The full-scene golden is the final verification of the ~120-line handleLasops verbatim port against a real X11 session.

2. **wave_legacy.gif** — blocks XvueExport_gif_AB_compare_wave. The GIF pipeline (ffmpeg dispatch, captureFrame hook, frame caps) is unit-tested; the visual A/B is the end-to-end integration check.

3. **cavity2d_legacy.gif** — same as above for the second canonical testa case.

### Do These Items Block Phase Close?

No — with one condition. The ROADMAP success criteria for Phase 7 are met in automated form for SC-1, SC-3 (automated half), SC-4, and SC-5. SC-2 (GIF visual A/B) and SC-3 (full-scene byte-parity) have QSKIP-gated test slots that are the explicit mechanism for tracking these open items. The test harness, bootstrap driver, and VALIDATION-LOG.md ledger were all committed as part of Plan 06.

**However:** Phase 8 ("A/B validation on testa subset") includes full visual A/B validation of testa/wave and testa/cavity2d as part of its own success criteria. The Phase 7 human items can therefore be merged into the Phase 8 human checkpoint — but they MUST be completed before Phase 8 can close.

**Recommended next action:** Before declaring Phase 7 complete and Phase 8 open, run Plan 06 Task 3 (the human checkpoint) using the procedure in VALIDATION-LOG.md §Procedure. This takes approximately 30 minutes and produces:
- `xvue/qt/tests/golden/scene01.eps` (flips Test 15 from QSKIP to PASS)
- `xvue/qt/tests/golden/wave_legacy.gif` (flips GIF A/B wave slot from QSKIP to PASS)
- `xvue/qt/tests/golden/cavity2d_legacy.gif` (flips GIF A/B cavity2d slot from QSKIP to PASS)
- Three PASS verdicts in VALIDATION-LOG.md replacing the current DEFERRED rows.

Once those three files are committed and ctest confirms 0 SKIP, Phase 7 can be declared complete and this VERIFICATION.md status updated to `passed`.

---

_Verified: 2026-05-04_
_Verifier: Claude (gsd-verifier), initial verification_
