---
phase: 0
slug: build-skeleton-abi-stubs
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-04-10
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

*Populated by the planner. Each task must map to a command that can be run from the repo root and exits non-zero on failure.*

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| TBD | TBD | TBD | BUILD-01..10 | — | N/A (no runtime behavior — link-complete stubs only) | build/link | TBD | ❌ W0 | ⬜ pending |

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
