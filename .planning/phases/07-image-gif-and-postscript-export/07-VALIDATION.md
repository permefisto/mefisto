---
phase: 7
slug: image-gif-and-postscript-export
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-05-04
---

# Phase 7 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> Filled by gsd-planner from 07-RESEARCH.md `## Validation Architecture`.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | GoogleTest + QTest (existing `xvue/qt/tests/`) + manual `bin/cbl_tout` build sweep |
| **Config file** | `xvue/qt/tests/CMakeLists.txt` (existing); add new test executables for Phase 7 |
| **Quick run command** | `cd xvue/qt/build && ctest -R xvue_qt_postscript_tests -j4 --output-on-failure` |
| **Full suite command** | `bin/cbxvueqt && cd xvue/qt/build && ctest -j4 --output-on-failure` |
| **Estimated runtime** | ~30 seconds (Qt unit tests); ~5 minutes including full xvueqt rebuild |

---

## Sampling Rate

- **After every task commit:** Run quick test for the touched component (PsEmitter / Export / Probe / GIF)
- **After every plan wave:** Run full Qt test suite (`ctest -j4`) plus `bin/cbl_tout_qt` build sweep
- **Before `/gsd-verify-work`:** Full suite green + `grep -rn 'convert' xvue/qt/` returns empty + PROBE.md committed + at least one byte-level PS golden compare passes against an X11-emitted reference
- **Max feedback latency:** ~30 seconds for unit-level feedback; ~5 minutes for the full A/B build sweep

---

## Per-Task Verification Map

> Filled in by gsd-planner during Plan 01–06 generation. Skeleton below establishes the contract.

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 07-01-01 | 01 | 1 | EXPORT-01 | T-07-01 (file-write to .planning/) | Probe writes only PROBE.md inside .planning/phases/07-…/ | smoke | `xvue/qt/build/probes/qimagewriter_probe \| tee .planning/phases/07-image-gif-and-postscript-export/PROBE.md` | ❌ W0 (probe binary build needed) | ⬜ pending |
| 07-02-01 | 02 | 2 | EXPORT-04 | T-07-02 (path-traversal via Fortran filename) | TEMPORAIRE.EPS opened via QFile in cwd, no path concat | unit | `ctest -R PsEmitter_handleLasops` | ❌ W0 | ⬜ pending |
| 07-02-02 | 02 | 2 | EXPORT-04 | — | xvpostscript_ body byte-for-byte port | golden | `ctest -R PsEmitter_postscriptVerbatim_golden` | ❌ W0 | ⬜ pending |
| 07-03-01 | 03 | 3 | EXPORT-04 | — | Each PsEmitter primitive emits identical bytes to xvuelc.c counterpart | golden | `ctest -R PsEmitter_perPrimitive_golden` | ❌ W0 | ⬜ pending |
| 07-04-01 | 04 | 2 | EXPORT-02 | T-07-03 (disk-fill, write-failure) | PNG/JPEG export uses QFileDialog-validated path; failure logged not silenced | integration | `ctest -R XvueExport_pngJpeg_roundTrip` | ❌ W0 | ⬜ pending |
| 07-04-02 | 04 | 2 | EXPORT-05 | — | PDF page geometry = canvas-native at 72 dpi | integration | `ctest -R XvueExport_pdf_geometry` | ❌ W0 | ⬜ pending |
| 07-05-01 | 05 | 3 | EXPORT-03 | T-07-04 (QProcess argv hygiene) | ffmpeg invoked via QProcess with constant argv list, no shell | integration | `ctest -R Animation_ffmpegFallback` | ❌ W0 | ⬜ pending |
| 07-05-02 | 05 | 3 | EXPORT-03 | T-07-05 (temp-dir cleanup) | xvue-gif-XXXXXX cleaned on success, kept-with-log on failure | integration | `ctest -R Animation_tempDirCleanup` | ❌ W0 | ⬜ pending |
| 07-06-01 | 06 | 4 | EXPORT-06 | T-07-06 (D-17 LVIDEO scope creep) | grep gate scoped to xvue/qt/ NOT xvue/ | smoke | `! grep -rn 'convert' xvue/qt/ \|\| grep -c 'convert' xvue/qt/ == 0` | ❌ W0 | ⬜ pending |
| 07-06-02 | 06 | 4 | EXPORT-03 | — | testa/wave + testa/cavity2d animation visually matches X11 reference (frame count + first/last frame hash) | A/B | manual side-by-side recorded in `.planning/phases/07-…/VALIDATION-LOG.md` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

> The `File Exists ❌ W0` rows are addressed by the Wave 0 requirements below.

---

## Wave 0 Requirements

> All test files / fixtures must exist before Wave 1 starts so that ❌ markers above flip to ✅ as tasks land.

- [ ] `xvue/qt/tests/xvue_qt_postscript_tests.cpp` — GTest stubs for `PsEmitter_handleLasops`, `PsEmitter_postscriptVerbatim_golden`, `PsEmitter_perPrimitive_golden` (one TEST per emit gate group)
- [ ] `xvue/qt/tests/xvue_qt_export_tests.cpp` — GTest + QTest stubs for `XvueExport_pngJpeg_roundTrip`, `XvueExport_pdf_geometry`, `Animation_ffmpegFallback`, `Animation_tempDirCleanup`
- [ ] `xvue/qt/tests/golden/` — directory for byte-level PS golden references; bootstrapped by running the X11 backend on a small canonical scene before Phase 7 ships
- [ ] `xvue/qt/probes/CMakeLists.txt` + `xvue/qt/probes/qimagewriter_probe.cpp` — probe binary skeleton (test infrastructure dependency for 07-01-01)
- [ ] `bin/cb_probe_qt` — shell wrapper for probe build + run + PROBE.md write
- [ ] Test data: small canonical scene driver (Fortran or test harness) that produces a deterministic PS output for the golden compare. Reuse `xvue/qt/tests/` existing test fixtures where possible (Phase 1–6 patterns: `xvue_qt_window_chrome_tests`, `xvue_qt_event_tests`).

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Animated-GIF visual match against legacy `bin/convertepsgif` output on `testa/wave` and `testa/cavity2d` | EXPORT-03 | No CI; image-similarity diff is eyeball-only per project convention (see `.planning/codebase/TESTING.md`). Frame-count + first/last hash automated; full visual parity is human eye. | 1. Run `testa/wave` on X11 backend → record convertepsgif `.gif`. 2. Run `testa/wave` on Qt backend → record Qt `.gif`. 3. Open both in image viewer; compare frame-by-frame. 4. Log result in `.planning/phases/07-…/VALIDATION-LOG.md`. |
| `_OMP` invariant under Phase 7 export | VALID-03 (Phase 8 cross-ref) | OpenMP solver parallelism stress; not part of Phase 7's gating but flagged by validation pyramid. | Defer to Phase 8. Phase 7 only ensures the export code path is main-thread-only via `XVUE_QT_ASSERT_MAIN_THREAD()`. |
| HiDPI export resolution check | (Claude's Discretion in CONTEXT.md) | Requires 4K monitor or `QT_SCALE_FACTOR=2` env. | Defer to Phase 8 VALID-04. Phase 7 documents the math in `xvue/qt/README.md`. |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies (filled by planner)
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references (probe binary, two new test files, goldens dir)
- [ ] No watch-mode flags (`ctest` is one-shot)
- [ ] Feedback latency < 30s for unit; < 5min for full sweep
- [ ] `nyquist_compliant: true` set in frontmatter (set by planner after task verify-map is filled)

**Approval:** pending
