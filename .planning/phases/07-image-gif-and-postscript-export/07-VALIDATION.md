---
phase: 7
slug: image-gif-and-postscript-export
status: ready
nyquist_compliant: true
wave_0_complete: false
created: 2026-05-04
last_updated: 2026-05-04
---

# Phase 7 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> Filled by gsd-planner after Plans 01–06 generation.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | GoogleTest + QTest (existing `xvue/qt/tests/`) + manual `bin/cbl_tout` build sweep |
| **Config file** | `xvue/qt/tests/CMakeLists.txt` (existing); Plans 02 / 04 add new test executables |
| **Quick run command** | `cd xvue/qt/build && xvfb-run --auto-servernum ctest -R '^xvue_qt_(postscript\|export)_tests$' --output-on-failure` |
| **Full suite command** | `bin/cbxvueqt && cd xvue/qt/build && xvfb-run --auto-servernum ctest -j4 --output-on-failure` |
| **Estimated runtime** | ~30 seconds (Qt unit tests); ~5 minutes including full xvueqt rebuild + bin/cbl_tout sweep |

---

## Sampling Rate

- **After every task commit:** Run quick test for the touched component (PsEmitter / Export / Probe / GIF).
- **After every plan wave:** Run full Qt test suite (`ctest -j4`) plus `bin/cbl_tout_qt` build sweep.
- **Before `/gsd-verify-work`:** Full suite green + `grep -rn 'convert' xvue/qt/` returns empty + PROBE.md committed + byte-level PS golden compare passes against an X11-emitted reference + `bin/cbl_tout` exits 0.
- **Max feedback latency:** ~30 seconds for unit-level feedback; ~5 minutes for the full A/B build sweep.

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|
| 07-01-01 | 01 | 1 | EXPORT-01 | T-07-01 | Probe writes only to PROBE.md inside .planning/phases/07-…/ | smoke | `cmake -S xvue/qt -B /tmp/xvueqt-probe-build -DXVUE_QT_BUILD_PROBES=ON && cmake --build /tmp/xvueqt-probe-build --target qimagewriter_probe` |
| 07-01-02 | 01 | 1 | EXPORT-01 | T-07-01 | bin/cb_probe_qt redirect target is hardcoded; PROBE.md committed | smoke | `bin/cb_probe_qt && test -s .planning/phases/07-image-gif-and-postscript-export/PROBE.md && grep -q 'gif_write_supported=' .planning/phases/07-image-gif-and-postscript-export/PROBE.md` |
| 07-02-01 | 02 | 2 | EXPORT-04 | T-07-02, T-07-09 | TEMPORAIRE.EPS hardcoded literal; fopen failure → QMessageBox + qApp->quit (NOT exit(1)) | unit | `cd xvue/qt/build && xvfb-run --auto-servernum ctest -R '^xvue_qt_postscript_tests$' --output-on-failure` |
| 07-03-01 | 03 | 3 | EXPORT-04 | T-07-07 | Per-primitive helpers use verbatim xvuelc.c format strings; Y-flip inside helpers only | unit | `cd xvue/qt/build && xvfb-run --auto-servernum ctest -R '^xvue_qt_postscript_tests$' --output-on-failure` |
| 07-03-02 | 03 | 3 | EXPORT-04 | T-07-08 | xvuelc.c untouched; bin/cbl_tout still green | smoke | `bin/cbl_tout && bin/cbl_tout_qt` |
| 07-04-01 | 04 | 2 | EXPORT-02 | T-07-03 | PNG/JPEG round-trip via QImageWriter; failure logged not silenced | integration | `cd xvue/qt/build && xvfb-run --auto-servernum ctest -R '^xvue_qt_export_tests$' --output-on-failure` |
| 07-04-02 | 04 | 2 | EXPORT-05 | (Pitfall 7) | PDF page geometry = canvas-native at 72 dpi (logical xpixels/ypixels, NOT backing) | integration | `cd xvue/qt/build && xvfb-run --auto-servernum ctest -R '^xvue_qt_export_tests$' --output-on-failure` |
| 07-05-01 | 05 | 3 | EXPORT-03 | T-07-04 | ffmpeg invoked via QProcess::execute(QStringList) — no shell, no user-typed args | integration | `cd xvue/qt/build && xvfb-run --auto-servernum ctest -R '^xvue_qt_export_tests$' --output-on-failure` |
| 07-05-02 | 05 | 3 | EXPORT-03 | T-07-05 | xvue-gif-XXXXXX QTemporaryDir cleaned on success, kept-with-log on failure | integration | `cd xvue/qt/build && xvfb-run --auto-servernum ctest -R '^xvue_qt_export_tests$' --output-on-failure` |
| 07-05-03 | 05 | 3 | EXPORT-03 | T-07-03 | Frame-count caps: warn at 100, hard-reject at 10000 | unit | `cd xvue/qt/build && xvfb-run --auto-servernum ctest -R '^xvue_qt_export_tests$' --output-on-failure` |
| 07-05-04 | 05 | 3 | EXPORT-03 | T-07-06 | LVIDEO files untouched: xvue/video1.f, videofin.f, videonm.f | smoke | `git diff --name-only HEAD~ -- xvue/video1.f xvue/videofin.f xvue/videonm.f` returns empty |
| 07-06-01 | 06 | 4 | EXPORT-04 | T-07-07 | Byte-level golden PsEmitter compare against scene01.eps | golden | `cd xvue/qt/build && xvfb-run --auto-servernum ctest -R 'PsEmitter_postscriptVerbatim_golden' --output-on-failure` |
| 07-06-02 | 06 | 4 | EXPORT-04 | T-07-07 | Per-primitive byte-level compare | golden | `cd xvue/qt/build && xvfb-run --auto-servernum ctest -R 'PsEmitter_perPrimitive_golden' --output-on-failure` |
| 07-06-03 | 06 | 4 | EXPORT-06 | T-07-06 | grep gate scoped to xvue/qt/ ONLY | smoke | `bash bin/test_no_imagemagick_in_qt.sh` |
| 07-06-04 | 06 | 4 | EXPORT-03 | (manual) | testa/wave + testa/cavity2d A/B sign-off vs X11+convertepsgif baseline | A/B (human) | `checkpoint:human-verify` recorded in `.planning/phases/07-image-gif-and-postscript-export/VALIDATION-LOG.md` |
| 07-06-05 | 06 | 4 | EXPORT-01..06 | T-07-08 | Final non-regression: full build green | smoke | `bin/cbl_tout && bin/cbl_tout_qt` |

*Status: each row marked ⬜ pending until Plan 06 close-out.*

---

## Wave 0 Requirements

All test files / fixtures must exist before Wave 1 starts so that ❌ markers above flip to ✅ as tasks land.

- [ ] `xvue/qt/probes/CMakeLists.txt` + `xvue/qt/probes/qimagewriter_probe.cpp` — probe binary (Plan 01 Task 1).
- [ ] `bin/cb_probe_qt` — shell wrapper for probe build + run + PROBE.md write (Plan 01 Task 2).
- [ ] `.planning/phases/07-image-gif-and-postscript-export/PROBE.md` — committed at Plan 01 Task 2.
- [ ] `xvue/qt/tests/test_xvue_qt_postscript.cpp` — Plan 02 Task 1 creates with 6 state-machine + Y-flip slots.
- [ ] `xvue/qt/tests/test_xvue_qt_export.cpp` — Plan 04 Task 1 creates with 5 PNG/JPEG/PDF slots; Plan 05 extends with 9 GIF slots.
- [ ] `xvue/qt/tests/golden/scene01.eps` — Plan 06 Task 1 bootstraps from X11 backend.
- [ ] `xvue/qt/tests/golden/scene01_driver.f` — Plan 06 Task 1 (reproducibility companion).
- [ ] `xvue/qt/tests/golden/wave_legacy.gif` — Plan 06 Task 3 (human checkpoint commits binary).
- [ ] `xvue/qt/tests/golden/cavity2d_legacy.gif` — Plan 06 Task 3.
- [ ] `bin/test_no_imagemagick_in_qt.sh` — Plan 06 Task 2 (EXPORT-06 grep gate).
- [ ] `xvue/qt/CMakeLists.txt` — Plans 02, 04, 05 each append source lines; Plan 01 adds `add_subdirectory(probes)` guarded; Plan 06 adds `verify_no_imagemagick_in_qt` ALL target.

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Animated-GIF visual match against legacy `bin/convertepsgif` on `testa/wave` and `testa/cavity2d` | EXPORT-03 | No CI; image-similarity diff is eyeball-only per project convention. Frame-count + first/last hash automated; full visual parity is human eye. | Plan 06 Task 3 — see plan body for full procedure. Result logged in `.planning/phases/07-image-gif-and-postscript-export/VALIDATION-LOG.md`. |
| `_OMP` invariant under Phase 7 export | VALID-03 (Phase 8 cross-ref) | OpenMP solver parallelism stress; not part of Phase 7's gating but flagged by validation pyramid. | Defer to Phase 8. Phase 7 only ensures the export code path is main-thread-only via `XVUE_QT_ASSERT_MAIN_THREAD()`. |
| HiDPI export resolution check | (Claude's Discretion in CONTEXT.md) | Requires 4K monitor or `QT_SCALE_FACTOR=2` env. | Defer to Phase 8 VALID-04. Phase 7 documents the math in `xvue/qt/README.md` Plan 06 Task 2. |

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify (filled by planner; `checkpoint:human-verify` for the testa/ A/B is the only manual case).
- [x] Sampling continuity: no 3 consecutive tasks without automated verify.
- [x] Wave 0 covers all MISSING references (probe binary, two new test files, goldens dir, grep-gate script).
- [x] No watch-mode flags (`ctest` is one-shot).
- [x] Feedback latency < 30s for unit; < 5min for full sweep.
- [x] `nyquist_compliant: true` set in frontmatter.

**Approval:** ready for execution.
