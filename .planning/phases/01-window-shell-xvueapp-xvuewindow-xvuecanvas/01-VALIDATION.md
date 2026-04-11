---
phase: 1
slug: window-shell-xvueapp-xvuewindow-xvuecanvas
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-04-11
---

# Phase 1 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Shell-script build verification + visual driver (`pp/ppxvtest0_qt`) + CMake build-time guards |
| **Config file** | `xvue/qt/CMakeLists.txt` (custom commands) + `bin/cbxvtest0_qt` |
| **Quick run command** | `bin/cbl_tout_qt` (legacy build still works + Qt build still works) |
| **Full suite command** | `bin/cbl_tout && bin/cbl_tout_qt && bin/cbxvtest0_qt && pp/ppxvtest0_qt` |
| **Estimated runtime** | ~60 seconds (builds) + manual visual confirmation |

---

## Sampling Rate

- **After every task commit:** Run `bin/cbl_tout_qt` (affected module rebuild)
- **After every plan wave:** Run full suite command
- **Before `/gsd-verify-work`:** Full suite green + visual driver confirmed
- **Max feedback latency:** ~60 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| TBD | — | — | SHELL-01..07 | — | Qt main-thread discipline enforced | build + visual | `bin/cbl_tout_qt` | ❌ W0 | ⬜ pending |

*Populated by gsd-planner in Step 8 — one row per task across all PLAN.md files.*

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `prpr/xvtest0.f` — minimal Fortran driver calling `xvinitgraphique_` / `xvfermer_` twice (D-19)
- [ ] `bin/cbxvtest0_qt` — build script for the test driver (D-20)
- [ ] `xvue/qt/cmake/verify_no_exec.sh` — CMake build-time guard for SHELL-03 (`QApplication::exec` ban)
- [ ] `xvue/qt/src/xvue_qt_app.{h,cpp}` — `XvueApp` singleton owner
- [ ] `xvue/qt/src/xvue_qt_window.{h,cpp}` — `XvueWindow` (QMainWindow subclass)
- [ ] `xvue/qt/src/xvue_qt_canvas.{h,cpp}` — `XvueCanvas` central widget
- [ ] `xvue/qt/src/xvue_qt_state.h` — header-only runtime state holder
- [ ] `XVUE_QT_ASSERT_MAIN_THREAD()` macro body in `xvue_qt_api.h`
- [ ] Macro retrofit across all ~57 `extern "C"` stubs in `xvue_qt_api.cpp`

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Blank QMainWindow appears on screen with XvueCanvas central widget | SHELL-01 | Requires display server | Run `pp/ppxvtest0_qt`, observe 800×600 window titled "MEFISTO" |
| Re-open after `xvfermer_` works without QApplication singleton assertion | SHELL-02 | Process-lifecycle behavior | Driver calls init→fermer→init→fermer; no assertion, no crash |
| HiDPI window renders identically at `QT_SCALE_FACTOR=2` | SHELL-06 | Requires display | Re-run driver with `QT_SCALE_FACTOR=2`, compare dimensions |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 90s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
