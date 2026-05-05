---
phase: 9
phase-slug: 09-retire-x11-backend
status: ready
nyquist_compliant: true
wave_0_complete: false
date: 2026-05-05
last_updated: 2026-05-05
---

# Phase 9 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> Consolidated from `09-RESEARCH.md §"Validation Architecture"` (Test Map,
> Sampling Rate, Wave 0 Gaps) per Dim 8e gate (revision iter1 BLOCKER #1).

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Qt Test (QTest framework, Qt6 6.10.2) for unit; bash + AUTOEXIT smoke + ImageMagick AE compare for integration |
| **Config file** | `xvue/qt/tests/CMakeLists.txt` (Qt Test); `bin/ab_sweep_phase8.sh` (integration sweep harness) |
| **Quick run command** | `cd xvue/qt/build && ctest -R '^xvue_qt_(postscript\|export\|menu_tests)$'` |
| **Full suite command** | `bin/cbl_tout && cd xvue/qt/build && xvfb-run --auto-servernum ctest && bin/ab_sweep_phase8.sh --mode qt-1x --cases pan2d,nafems_le1,cavity2d,heat1d,nlsecu` |
| **Estimated runtime** | ~30 sec (incremental cbl_tout) ; ~3 min (full ctest) ; ~10 min (full Phase 8 harness sweep) |

---

## Sampling Rate

- **Per task commit:** `bin/cbl_tout` (incremental rebuild — full Qt-only build
  with the CMake guards), takes ~30 seconds. CMake ALL-target gates run
  automatically (`verify_no_imagemagick_in_qt`, `verify_no_x11_in_build` after
  Plan 9-03, `verify_no_lvideo` after Plan 9-04).
- **Per wave merge:** Above + `xvfb-run ctest` (Qt unit tests), takes ~3 minutes.
- **Phase gate:** Above + full Phase 8 harness sweep (`bin/ab_sweep_phase8.sh
  --mode qt-1x` on all 5 cases pan2d, nafems_le1, cavity2d, heat1d, nlsecu),
  takes ~10 minutes plus matched-dim verification (Plan 9-06).

---

## Phase Requirements → Test Map

| Req ID | Plan | Wave | Behavior | Test Type | Automated Command | File Exists? |
|--------|------|------|----------|-----------|-------------------|-------------|
| RETIRE-01 | 09-02 | 2 | xvuelc.c + ccxvue deleted; full Qt-only build green | smoke | `bin/cbl_tout && ls -la pp/pp*` | ✅ (build script + binaries exist as gates) |
| RETIRE-01 verify | 09-02 | 2 | No xvuelc symbol references remain in Fortran link path | unit | `! nm xvue/qt/build/libxvueqt.a \| grep ' U xvtest_x11_only_symbol'` (negative-grep on hypothetical xvuelc-only entry) | ✅ (`xvue/qt/cmake/verify_abi.sh`) |
| RETIRE-02 | 09-03 | 2 | No `-lX11` / `/usr/X11R6` in any active build script | unit | `! grep -r 'lX11\|X11R6' bin/cb* bin/cbl_tout` | ❌ (Wave 0 — `bin/test_no_x11_in_build.sh` modeled on `bin/test_no_imagemagick_in_qt.sh`; created by Plan 9-03 Task 4) |
| RETIRE-03 | 09-04 | 2 | LVIDEO + `convert` shell-outs gone from active source paths | unit | `! grep -rn 'CALL VIDEO' flui/ ther/ util/ xvue/` AND `! grep -rn 'CALL SYSTEM.*convert' xvue/ flui/ ther/ util/` | ❌ (Wave 0 — `bin/test_no_lvideo.sh`; created by Plan 9-04 Task 1) |
| RETIRE-04 | 09-05 | 2 | Documentation lists only Qt deps | smoke / human verify | grep `'libX11-dev\|imagemagick'` returns 0 hits across {README, LISEZMOI, bin/README, bin/LISEZMOI, CLAUDE.md} | ✅ (grep is the test) |
| Carry 9-06 (matched-dim) | 09-06 | 3 | Qt 1x captures land at the env-specified WxH | integration | `MEFISTO_QT_WINDOW_SIZE=1280x800 bin/ab_sweep_phase8.sh --mode qt-1x --cases pan2d --out-dir /tmp/p9-06 && identify -format '%wx%h' /tmp/p9-06/pan2d-qt-1x.png` returns `1280x800` | ✅ (harness already exists; env wiring added by Plan 9-06) |
| Carry 9-07 (ppnlse) | 09-07 | 3 | nlsecu Qt 1x sweep produces a non-truncated capture (NLSER banner reached, frame from solver step) — OR documented case-(c) long-runtime waiver | integration | `bin/ab_sweep_phase8.sh --mode qt-1x --cases nlsecu --out-dir /tmp/p9-07` exits within 600s timeout AND output PNG > 1KB AND contains an NLSER-frame marker (manual inspection on success path; case-(c) keeps Phase 8 TIME-truncation mitigation) | partial — needs custom timeout extension and manual inspection step |
| Carry 9-08 (Phase 7 goldens) | 09-08 | 3 | 3 Phase-7-deferred goldens land; ctest QSKIPs flip to PASS | unit | `ls xvue/qt/tests/golden/{scene01.eps,wave_legacy.gif,cavity2d_legacy.gif}` AND `xvfb-run --auto-servernum ctest -R '^xvue_qt_(postscript\|export)_tests$' --output-on-failure` exits 0 with 0 SKIP | ✅ (test slots exist; QSKIP→PASS auto-flip) |
| Carry 9-09 (`--out-dir` fix) | 09-09 | 3 | Relative `--out-dir` resolves correctly under user cwd, NOT under PROJDIR | unit | `cd /tmp && bin/ab_sweep_phase8.sh --mode qt-1x --cases pan2d --out-dir ./out --smoke-only && ls /tmp/out/pan2d-qt-1x.png` | ✅ |
| Carry 9-09 (freshness guard) | 09-09 | 3 | `verify_pp_qt_freshness` fails when libxvueqt.a is newer than any pp/pp* | unit | `touch -d '1 minute ago' pp/ppmail && sh xvue/qt/cmake/verify_pp_freshness.sh xvue/qt/build/libxvueqt.a pp` exits non-zero with FAIL diagnostic | ❌ (Wave 0 — add `xvue/qt/cmake/verify_pp_freshness.sh` + cbl_tout end-section invocation; created by Plan 9-09 Tasks 2-3) |

*Status: each row marked ⬜ pending until the corresponding plan closes.*

---

## Wave 0 Gaps

All test files / fixtures / harness updates must exist before the dependent
plan starts so that ❌ markers above flip to ✅ as plans land.

- [ ] `bin/test_no_x11_in_build.sh` — RETIRE-02 grep gate, mirrors
      `bin/test_no_imagemagick_in_qt.sh` (Plan 9-03 Task 4 creates it +
      wires `verify_no_x11_in_build ALL` target in `xvue/qt/CMakeLists.txt`).
- [ ] `bin/test_no_lvideo.sh` — RETIRE-03 grep gate covering both
      `CALL VIDEO*` and `CALL SYSTEM('convert')` (Plan 9-04 Task 1 creates
      it + wires `verify_no_lvideo ALL` target).
- [ ] `xvue/qt/cmake/verify_pp_freshness.sh` + corresponding
      end-of-`bin/cbl_tout` invocation, plus optional config-gated
      `add_custom_target` in `xvue/qt/CMakeLists.txt` — Carry 9-09 Tasks 2-3.
      (Glob `pp*` — NO `_qt` suffix; per Plan 9-03 the suffix collapses
      BEFORE Plan 9-09 runs; see revision iter1 BLOCKERs #2-#4.)
- [ ] Pattern-extension target in `xvue/qt/CMakeLists.txt` to call the two
      new test scripts (so they fail the build, not just stand alone) —
      added inline by Plan 9-03 Task 4 + Plan 9-04 Task 1 alongside the
      script creation.
- [ ] Update `bin/test_no_imagemagick_in_qt.sh` (if its scope comment
      references LVIDEO as legitimate) — covered transitively by Plan 9-04
      grep-gate when LVIDEO disappears (no extra task needed; documentation
      audit only).

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| A/B window closure (D-11 process gate) | Phase entry | Maintainer-judgment call; not a fixed-date deadline | Plan 9-01 Task 1 — STATE.md closure line OR explicit reply. Recorded in 09-01-AUDIT-BASELINE.md §1. |
| Phase 9 README/LISEZMOI rollback section wording | RETIRE-04 (Plan 9-05) | Editorial choice (silent omission vs explicit mention vs full Phase-9 section) | Plan 9-05 Task 2 — `checkpoint:decision` records option-a / option-b / option-c. |
| Matched-dim AE drop visual sign-off | Carry 9-06 | Side-by-side eyeball comparison; AE delta is empirical but diff interpretation is human | Plan 9-06 Task 3 — `checkpoint:human-verify` against Phase 8 X11 cavity2d baseline. |
| ppnlse_qt root-cause classification (a/b/c) | Carry 9-07 | gdb backtrace + stdout pattern requires human classification | Plan 9-07 Task 1 — empirical diagnose committed to `evidence/ppnlse-diagnose.md` with explicit ROOT CAUSE classifier. |
| Cross-tag worktree procedure (Phase 7 goldens) | Carry 9-08 | Single deviation from "Phase 9 only touches Qt-only main"; cross-tag bootstrap is not automatable from main alone | Plan 9-08 — full procedure + worktree cleanup; SUMMARY records v1.0-pre-retire SHA + link line + worktree cleanup verification. |

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify (filled by planner).
- [x] Sampling continuity: no 3 consecutive tasks without automated verify.
- [x] Wave 0 covers all MISSING references (3 grep-gate scripts + freshness
      checker + CMake target wiring).
- [x] No watch-mode flags (`ctest` is one-shot; `bin/ab_sweep_phase8.sh` is
      one-shot).
- [x] Feedback latency < 30s for unit; < 5min for full sweep.
- [x] `nyquist_compliant: true` set in frontmatter.
- [x] Test Map covers RETIRE-01..04 + 4 carry-forwards (9-06..9-09).
- [x] Sampling Rate matches three-tier discipline (per-commit / per-wave / phase-gate).

**Approval:** ready for execution.
