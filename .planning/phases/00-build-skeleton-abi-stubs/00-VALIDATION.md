---
phase: 0
slug: build-skeleton-abi-stubs
status: approved
nyquist_compliant: true
wave_0_complete: true
created: 2026-04-10
updated: 2026-04-10
---

# Phase 0 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Shell scripts + CMake targets (no unit-test framework; link-time verification is the test) |
| **Config file** | `xvue/qt/CMakeLists.txt` (verify_abi target), `bin/cbl_tout_qt` |
| **Quick run command** | `cmake --build build-qt --target verify_abi` |
| **Full suite command** | `bin/cbl_tout_qt && for e in ppmail ppelas ppflui ppther ppnlse; do pp/${e}_qt --version-stub || echo NOOP-OK; done && bin/cbl_tout` |
| **Estimated runtime** | ~60–120 seconds (both Qt and legacy builds) |

---

## Sampling Rate

- **After every task commit:** Run `cmake --build build-qt --target verify_abi` (when build-qt exists) OR the task's automated check
- **After every plan wave:** Run the full Qt + legacy build suite
- **Before `/gsd-verify-work`:** Both `bin/cbl_tout_qt` AND `bin/cbl_tout` must succeed; `nm libxvueqt.a | grep -c '_$'` must match the documented ABI count (57 after D-04/D-06 resolution: C-internal helpers excluded from public header)
- **Max feedback latency:** 120 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | Status |
|---------|------|------|-------------|-----------|-------------------|--------|
| 00-01-01 | 01 | 1 | BUILD-01 (env) | human-checkpoint | `dpkg -s qt6-base-dev qt6-base-dev-tools` and `cmake --find-package -DNAME=Qt6 -DCOMPILER_ID=GNU -DLANGUAGE=CXX -DMODE=EXIST` | ⬜ pending |
| 00-01-02 | 01 | 1 | BUILD-07 (legacy anchor) | build | `cd $MEFISTO && bin/cbl_tout && test -x pp/ppmail && test -x pp/ppelas && test -x pp/ppflui && test -x pp/ppther && test -x pp/ppnlse` | ⬜ pending |
| 00-02-01 | 02 | 2 | BUILD-04 (header) | grep | `test -f xvue/qt/include/xvue_qt_api.h && grep -c 'proc(xv' xvue/qt/include/xvue_qt_api.h` (expect 57) | ⬜ pending |
| 00-02-02 | 02 | 2 | BUILD-05 (stubs) | compile | `g++ -c xvue/qt/src/xvue_qt_api.cpp -Ixvue/qt/include -fPIC -o /tmp/stubs.o` | ⬜ pending |
| 00-02-03 | 02 | 2 | BUILD-01, BUILD-02, BUILD-03, BUILD-08 | cmake+nm | `cmake -S xvue/qt -B xvue/qt/build && cmake --build xvue/qt/build && nm xvue/qt/build/libxvueqt.a \| grep -c ' T [a-zA-Z_][a-zA-Z0-9_]*_$'` (expect 57) | ⬜ pending |
| 00-03-01 | 03 | 3 | BUILD-06 (cbl_tout_qt) | script | `test -x bin/cbl_tout_qt && bash -n bin/cbl_tout_qt` | ⬜ pending |
| 00-03-02 | 03 | 3 | BUILD-06 (cb*_qt clones) | script | `for m in cbmail cbelas cbflui cbther cbnlse; do test -x bin/${m}_qt && bash -n bin/${m}_qt; done` | ⬜ pending |
| 00-03-03 | 03 | 3 | BUILD-06 (smoke test) | run | `bin/cbl_tout_qt && for e in ppmail ppelas ppflui ppther ppnlse; do test -x pp/${e}_qt; done` | ⬜ pending |
| 00-04-01 | 04 | 4 | BUILD-09, BUILD-10 | files | `test -f xvue/README_COORDS.md && test -f .planning/validation/BASELINE.md` | ⬜ pending |
| 00-04-02 | 04 | 4 | BUILD-07 (regression) | human-verify | Clean tree rebuild: `git clean -fdx pp xvue/qt/build && bin/cbl_tout_qt && bin/cbl_tout && run 5 canonical testa cases` | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `xvue/qt/CMakeLists.txt` — skeleton with `verify_abi` custom target
- [ ] `xvue/qt/xvue_qt_api.h` — 57 entry-point declarations (Option A: D-04/D-06 treated as C++ internals)
- [ ] `xvue/qt/xvue_qt_stubs.cpp` — no-op implementations with trailing-underscore symbol names
- [ ] `bin/cbl_tout_qt` — top-level Qt build script (clone of `bin/cbl_tout`)
- [ ] `bin/cbmail_qt`, `bin/cbelas_qt`, `bin/cbflui_qt`, `bin/cbther_qt`, `bin/cbnlse_qt` — per-module Qt link scripts
- [ ] `.planning/validation/BASELINE.md` — 5 canonical testa/ cases (one per module)
- [ ] `xvue/README_COORDS.md` — Y-axis convention (Y-down top-left, PS emitter inverts via `ypixels - *y`)

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Legacy X11 `bin/cbl_tout` still builds unchanged pp/pp* | BUILD-09 | Requires X11 DISPLAY or headless xvfb setup; baseline compare by hash | After Phase 0, run `bin/cbl_tout`, confirm `pp/ppmail pp/ppelas pp/ppflui pp/ppther pp/ppnlse` exist and re-run the 5 canonical `testa/` cases listed in BASELINE.md |
| `pp/pp*_qt` executables proceed past link stage without crash | BUILD-08 | Depends on apt-installed Qt 6 present | Run each `pp/pp*_qt` binary; expect clean exit from no-op stubs (no graphics output) |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify command or Wave 0 dependency
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references (CMake scaffold, cb*_qt scripts, header + stubs)
- [ ] No watch-mode flags
- [ ] Feedback latency < 120s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
