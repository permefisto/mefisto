---
phase: 1
slug: window-shell-xvueapp-xvuewindow-xvuecanvas
status: draft
nyquist_compliant: true
wave_0_complete: true
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
| **Quick run command** | `bin/cbl_tout_qt` (Qt build path — ~60s) |
| **Full suite command** | `bin/cbl_tout && bin/cbl_tout_qt && bin/cbxvtest0_qt && pp/ppxvtest0_qt` (phase-gate only) |
| **Quick run runtime** | ~60 seconds |
| **Full suite runtime** | unbounded (legacy X11 build + Qt build + driver build + visual run) — explicitly excluded from per-task latency budget |

---

## Sampling Rate (Tiered)

Phase 1 uses a **two-tier Nyquist sampling model** because the full legacy X11 build (`bin/cbl_tout`) alone typically exceeds 60 seconds, violating the per-task latency budget if run after every commit.

### Tier 1 — Fast per-task sampling (per commit)

- **Command:** `bin/cbl_tout_qt`
- **Runtime:** ~60 seconds
- **Use after:** Every task commit in Plans 01-01 and 01-02 (C++ / CMake / macro retrofit work)
- **Purpose:** Catches Qt build regressions immediately; preserves per-task feedback latency ≤ 60s

### Tier 2 — Slow phase-gate sampling (per plan completion)

- **Command:** `bin/cbl_tout && bin/cbl_tout_qt && bin/cbxvtest0_qt && pp/ppxvtest0_qt`
- **Runtime:** unbounded (declared out of per-task latency budget)
- **Use after:** Plan 01-03 completion (and once at end of Plan 01-02 as a pre-gate sanity check)
- **Purpose:** End-to-end validation of legacy X11 regression guard + Qt build + driver build + visual SHELL-01/02/06 confirmation
- **Runs:** once per plan wave minimum; once before `/gsd-verify-work` mandatory

**Per-task latency budget:** ≤ 60 seconds (Tier 1 only; Tier 2 is excluded from this budget by design)

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 01-01 T1 | 01-01 | 1 | SHELL-07 | T-01-03 | Thread-affinity macro refuses non-main-thread entry | grep + build | `grep -q 'XVUE_QT_ASSERT_MAIN_THREAD' xvue/qt/include/xvue_qt_api.h && bin/cbl_tout_qt` | ✅ W0 (Plan 01-01 T2 materialises headers) | ⬜ pending |
| 01-01 T2 | 01-01 | 1 | SHELL-01, SHELL-02 | T-01-01, T-01-02 | XvueApp singleton + lazy XvueWindow scaffolding compile under AUTOMOC | build | `bin/cbl_tout_qt` | ❌ W0 (creates xvue_qt_app/window/canvas/state sources) | ⬜ pending |
| 01-01 T3 | 01-01 | 1 | SHELL-03 | T-01-04 | Build-time grep-and-fail guard bans `QApplication::exec` and `qApp->exec()` | CMake custom target | `bin/cbl_tout_qt` (target runs verify_no_exec; fails red on violation) | ❌ W0 (creates xvue/qt/cmake/verify_no_exec.sh) | ⬜ pending |
| 01-02 T1 | 01-02 | 2 | SHELL-01, SHELL-02, SHELL-04, SHELL-05, SHELL-06 | T-02-01, T-02-02 | 7 SHELL entry points implement real bodies without regressing ABI | build | `bin/cbl_tout_qt` | ✅ (extends xvue_qt_api.cpp from Phase 0) | ⬜ pending |
| 01-02 T2 | 01-02 | 2 | SHELL-07 | T-02-03 | Every non-SHELL stub gains `XVUE_QT_ASSERT_MAIN_THREAD()` at entry | grep count + build | `grep -c 'XVUE_QT_ASSERT_MAIN_THREAD' xvue/qt/src/xvue_qt_api.cpp && bin/cbl_tout_qt` (expect ≥ 57) | ✅ (retrofits Phase 0 stubs) | ⬜ pending |
| 01-03 T1 | 01-03 | 3 | SHELL-01, SHELL-02, SHELL-06 | — | Fortran 77 lifecycle driver is well-formed fixed-form and calls XVINITGRAPHIQUE/XVFERMER twice each | grep + column lint | `test -f prpr/xvtest0.f && grep -c 'CALL XVINITGRAPHIQUE' prpr/xvtest0.f && grep -c 'CALL XVFERMER' prpr/xvtest0.f && awk 'length($0)>72 {exit 1}' prpr/xvtest0.f` | ❌ W0 (creates prpr/xvtest0.f) | ⬜ pending |
| 01-03 T2 | 01-03 | 3 | SHELL-01, SHELL-02 | T-03-01 | Build script clones cbmail_qt template, produces pp/ppxvtest0_qt, drops mail/lib | build | `bin/cbxvtest0_qt && test -x pp/ppxvtest0_qt` | ❌ W0 (creates bin/cbxvtest0_qt) | ⬜ pending |
| 01-03 T3 | 01-03 | 3 | SHELL-01, SHELL-02, SHELL-06 (+ regression SHELL-03, SHELL-04, SHELL-05, SHELL-07) | T-03-01, T-03-02, T-03-03 | Full-suite phase-gate: legacy X11 + Qt build + driver build + visual run succeed; two blank MEFISTO windows observed; no QApplication singleton assertion; HiDPI scale works | full-suite build + visual checkpoint | `cd $MEFISTO && bin/cbl_tout && bin/cbl_tout_qt && bin/cbxvtest0_qt && test -x pp/ppmail && test -x pp/ppmail_qt && test -x pp/ppxvtest0_qt` (Tier 2 — excluded from per-task latency budget) | ✅ (all deps materialised by earlier tasks) | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

**Sampling continuity check:** Every task (01-01 T1 through 01-03 T3) has at least one automated verify command. No three consecutive tasks go without a Tier 1 (`bin/cbl_tout_qt`) sample — Tasks 01-01 T1/T2/T3 and 01-02 T1/T2 all run the Qt build directly, and Tasks 01-03 T1/T2/T3 form the Tier 2 phase-gate where 01-03 T2 (`bin/cbxvtest0_qt`) also exercises the compile path.

---

## Wave 0 Requirements

- [x] `prpr/xvtest0.f` — minimal Fortran driver calling `xvinitgraphique_` / `xvfermer_` twice (D-19) — scheduled in Plan 01-03 Task 1
- [x] `bin/cbxvtest0_qt` — build script for the test driver (D-20) — scheduled in Plan 01-03 Task 2
- [x] `xvue/qt/cmake/verify_no_exec.sh` — CMake build-time guard for SHELL-03 (`QApplication::exec` ban) — scheduled in Plan 01-01 Task 3
- [x] `xvue/qt/src/xvue_qt_app.{h,cpp}` — `XvueApp` singleton owner — scheduled in Plan 01-01 Task 2
- [x] `xvue/qt/src/xvue_qt_window.{h,cpp}` — `XvueWindow` (QMainWindow subclass) — scheduled in Plan 01-01 Task 2
- [x] `xvue/qt/src/xvue_qt_canvas.{h,cpp}` — `XvueCanvas` central widget — scheduled in Plan 01-01 Task 2
- [x] `xvue/qt/src/xvue_qt_state.h` — header-only runtime state holder (per RESEARCH.md Q2 RESOLVED) — scheduled in Plan 01-01 Task 2
- [x] `XVUE_QT_ASSERT_MAIN_THREAD()` macro body in `xvue_qt_api.h` — scheduled in Plan 01-01 Task 1
- [x] Macro retrofit across all ~57 `extern "C"` stubs in `xvue_qt_api.cpp` — scheduled in Plan 01-02 Task 2

All Wave 0 gaps are fully scheduled across Plans 01-01 (Wave 1), 01-02 (Wave 2), and 01-03 (Wave 3). `wave_0_complete: true` reflects that every MISSING-reference in the per-task map is claimed by an explicit task before its downstream consumer runs.

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Blank QMainWindow appears on screen with XvueCanvas central widget | SHELL-01 | Requires display server | Run `pp/ppxvtest0_qt`, observe 800×600 window titled "MEFISTO" |
| Re-open after `xvfermer_` works without QApplication singleton assertion | SHELL-02 | Process-lifecycle behavior | Driver calls init→fermer→init→fermer; no assertion, no crash |
| HiDPI window renders identically at `QT_SCALE_FACTOR=2` | SHELL-06 | Requires display | Re-run driver with `QT_SCALE_FACTOR=2`, compare dimensions |

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies (per per-task map above)
- [x] Sampling continuity: no 3 consecutive tasks without automated verify (every task samples `bin/cbl_tout_qt` or an equivalent build-path check)
- [x] Wave 0 covers all MISSING references (9/9 items scheduled)
- [x] No watch-mode flags (all commands are one-shot)
- [x] Feedback latency tiered: Tier 1 (per-task) ≤ 60s via `bin/cbl_tout_qt`; Tier 2 (phase-gate) unbounded and explicitly excluded from per-task budget
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
