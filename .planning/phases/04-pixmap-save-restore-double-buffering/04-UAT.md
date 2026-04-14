---
status: complete
phase: 04-pixmap-save-restore-double-buffering
source:
  - 04-01-SUMMARY.md
  - 04-02-SUMMARY.md
started: 2026-04-14T14:10:00Z
updated: 2026-04-14T15:10:00Z
---

## Current Test

[testing complete]

## Tests

### 1. Pixmap round-trip harness — all 4 pairs
expected: |
  `bin/xvtest0-pixmap-roundtrip.sh` completes with all 4 PIXMAP pairs AE=0 and
  prints "ALL 4 round-trip pairs PASS — Phase 4 green". Exit code 0.
result: pass
notes: |
  Verified orchestrator-run. All 4 pairs PASS (PIXMAP-01/02/03a/03b, AE=0),
  exit 0. 6 scene PNGs captured under /tmp.

### 2. Qt backend rebuild
expected: |
  `bin/cbl_tout_qt` completes with exit code 0. Produces/refreshes
  `pp/ppxvtest0_qt` (~500 KB) and `xvue/qt/build/libxvueqt.a`.
result: pass
notes: |
  Serial `bin/cbl_tout_qt` exit 0. All 5 ppxvtest*_qt binaries present
  (529–541 KB) + solvers (ppelas_qt, ppflui_qt, ppmail_qt, ppnlse_qt,
  ppther_qt). `xvue/qt/build/libxvueqt.a` = 727 KB.

### 3. Legacy X11 backend still builds
expected: |
  `bin/cbl_tout` completes with exit code 0. Fortran 77 modules compile and
  link against the unchanged X11 xvue layer. No regressions from Phase 4
  changes (which are Qt-side only).
result: pass
notes: |
  Serial `bin/cbl_tout` exit 0. `ppxvtest4 is CORRECT`. All 5 legacy
  xvtest*_x11 binaries built plus full legacy solver set. Earlier
  "INCORRECT" signal traced to parallel build collision on shared `pp/`
  directory, not a Phase 4 regression. Lesson: do NOT run `cbl_tout`
  and `cbl_tout_qt` in parallel — they share `pp/`.

### 4. Legacy xvtest0 driver path (no Phase 4 scene)
expected: |
  Running `pp/ppxvtest0_qt` without `MEFISTO_XVTEST0_SCENE` set reaches the
  Phase 2 DRAW + Phase 3 TEXT blocks and prints
  `[xvtest0] OK — cycle open/close/open/close + draws`.
result: pass
notes: |
  `QT_QPA_PLATFORM=offscreen pp/ppxvtest0_qt` (unset env var) completes
  the open/close/open/close + draws cycle. Phase 4 coverage block falls
  through cleanly when no scene is selected. rc=0.

### 5. No regression in xvtest1..4 drivers
expected: |
  `pp/ppxvtest{1,2,3,4}_qt` run to completion headlessly without new errors
  or stub warnings. The 7 former-stub symbols no longer emit warn-once.
result: pass
notes: |
  All 4 drivers exit rc=0 under QT_QPA_PLATFORM=offscreen. No stub
  warnings for the 7 pixmap symbols. Normal `ENTRER UN CARACTERE POUR
  FINIR` prompts only.

### 6. ABI & exec-surface gates
expected: |
  `verify_abi` reports `nm count: 57  header count: 57`; `verify_no_exec`
  finds no forbidden tokens and passes the palette scan.
result: pass
notes: |
  `xvue/qt/cmake/verify_abi.sh build/libxvueqt.a include/xvue_qt_api.h`
  → `nm count: 57  header count: 57`, rc=0.
  `xvue/qt/cmake/verify_no_exec.sh src src` → no forbidden tokens +
  palette-leak scan clean, rc=0.

## Summary

total: 6
passed: 6
issues: 0
pending: 0
skipped: 0
blocked: 0

## Gaps

[none]
