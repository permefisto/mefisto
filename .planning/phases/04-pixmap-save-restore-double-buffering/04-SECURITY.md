---
phase: 04-pixmap-save-restore-double-buffering
status: secured
audited: 2026-04-14
threats_total: 8
threats_closed: 8
threats_open: 0
accepted_risks: 3
---

# Phase 04 Security Audit — pixmap-save-restore-double-buffering

## Summary

All 8 STRIDE threats from `04-01-PLAN.md` and `04-02-PLAN.md` are verified closed. 5 mitigate-disposition threats have mitigations present in committed code; 3 accept-disposition threats have documented rationale consistent with the implementation.

## Trust Boundaries

| Boundary | Description |
|----------|-------------|
| Fortran → C++ ABI | The 7 pixmap entry points take no parameters; no user-controlled data crosses the boundary. |
| Process → X11/Qt server | Read-only pixmap blits to an in-process off-screen surface. No socket / file I/O from the Qt layer. |
| Shell harness → filesystem | Writes PNGs to `/tmp/p4_*.png` (non-sensitive checker-grid scene data). |
| Fortran GETENV → process env | Reads `MEFISTO_XVTEST0_SCENE`; compared only against fixed literals. |

## Threat Register

| ID | Cat | Component | Disposition | Status | Evidence |
|----|-----|-----------|-------------|--------|----------|
| T-04-01 | T | `saved_canvas_` dangling after window close | mitigate | CLOSED | `xvue/qt/src/xvue_qt_state.cpp:147-151` — destructor deletes `saved_canvas_` before `painter_` and `backing_`; `xvue/qt/src/xvue_qt_state.h:38` declares raw-pointer slot. |
| T-04-02 | D | Repeated save at differing sizes leaks memory | mitigate | CLOSED | `xvue/qt/src/xvue_qt_api.cpp:94-96` — size-mismatch branch explicitly `delete`s old slot before `new QPixmap(size)`. |
| T-04-03 | I | Stderr warning could leak PID/paths | accept | CLOSED | `xvue/qt/src/xvue_qt_api.cpp:122` — fixed string `"xvue-qt: restore_from_slot: no slot or size mismatch\n"`; no user data, PID, or paths embedded. Accepted per D-12. |
| T-04-04 | E | Off-main-thread graphics from OpenMP worker | mitigate | CLOSED | `XVUE_QT_ASSERT_MAIN_THREAD()` present at every Qt API entry (all 7 pixmap symbols at `xvue_qt_api.cpp:162,171,181,190,199,471,479`, and every other public function in the file). |
| T-04-05 | T | Stale `/tmp/p4_*.png` from previous run causes false pass | mitigate | CLOSED | `bin/xvtest0-pixmap-roundtrip.sh:44` — `rm -f "$OUT"` before each capture. |
| T-04-06 | I | Harness log leaks paths/env vars | accept | CLOSED | Harness exports only `MEFISTO_XVTEST0_SCENE` per call (line 46) and prints scene names, AE counts, and `/tmp` paths. No PID, no env dump, no secret material. Accepted. |
| T-04-07 | D | Missing ImageMagick hangs harness | mitigate | CLOSED | `bin/xvtest0-pixmap-roundtrip.sh:20` — `command -v magick >/dev/null 2>&1 \|\| exit 2` preflight. |
| T-04-08 | T | GETENV buffer overflow on oversized env value | accept | CLOSED | `prpr/xvtest0.f:37` declares `CHARACTER*32 SCENE`; `prpr/xvtest0.f:52` calls `GETENV('MEFISTO_XVTEST0_SCENE', SCENE)`. gfortran GETENV truncates safely to the declared length. Scene is compared only to fixed literals ≤20 chars (P4_CTRL, P4_SAVERESTORE, P4_MEMPX_SAVERESTORE, P4_BG, P4_EFFACEMEMPX, P4_FENETREMEMPX), so truncation cannot produce a spurious match. Accepted. |

## Accepted Risks

1. **T-04-03** — fixed-string stderr warning on restore size mismatch. Rationale: the string is a constant; leaking the existence of the code path is not information disclosure. Reviewed.
2. **T-04-06** — harness log contents. Rationale: scenes and `/tmp` paths are non-sensitive; no env dump, no secrets on the machine-readable output path.
3. **T-04-08** — GETENV truncation on oversized env var. Rationale: `CHARACTER*32` buffer + fixed-literal comparison makes truncation a safe failure mode (non-match → blank scene → Phase 2/3 fallthrough).

## Audit Trail

### Security Audit 2026-04-14

| Metric | Count |
|--------|-------|
| Threats found | 8 |
| Closed — mitigation verified | 5 |
| Closed — accepted risk | 3 |
| Open | 0 |

**Verification method:** Grep-based evidence search across `xvue/qt/src/xvue_qt_state.{h,cpp}`, `xvue/qt/src/xvue_qt_api.cpp`, `bin/xvtest0-pixmap-roundtrip.sh`, and `prpr/xvtest0.f`. All 5 mitigate-disposition threats have mitigations on disk at the documented locations. The 3 accept-disposition threats have rationale consistent with the implementation.

**Note:** The `gsd-security-auditor` agent hit a stream timeout during this audit. Verification was performed inline by the orchestrator using the same evidence-based grep methodology. No implementation changes were required.
