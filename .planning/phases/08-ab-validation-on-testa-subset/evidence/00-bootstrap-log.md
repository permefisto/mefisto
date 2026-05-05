# Phase 8 Plan 1 — Bootstrap Log

**Plan:** 08-01 — Bootstrap: pp/*_qt freshness + 3 deferred goldens + harness scripts + 5-case smoke
**Start:** 2026-05-05
**Author:** Phase 8 Plan 1 executor (autonomous)

This document records per-task bootstrap evidence for the four Phase 8 Plan 1
tasks (Task 0 — case-batch map; Task 1 — pp/*_qt freshness; Task 2 — 3 deferred
goldens + ctest re-verify; Task 3 — sweep-harness scripts + smoke probes).

---

## Task 1 — pp/*_qt freshness (D-08)

**Goal:** Refresh `pp/*_qt` against the current `libxvueqt.a` so the A/B
sweep cannot get a false-PASS from stale Qt binaries (Phase 7 Gap-A).

### bin/cbl_tout_qt result

```
cbl_tout_qt_exit=0
```

Build duration: ~7 min (CMake configure + libxvueqt.a build + 4 test-binary
targets + 7 Fortran library rebuilds + 5 cb*_qt + 5 cbxvtest*_qt linker
invocations). Last log line:

```
TOUS les MODULES EXECUTABLES sans debogueur de /home/mefisto/git/mefisto sous LINUX sont crees — Qt variant
```

### Freshness check (pp/*_qt mtime >= libxvueqt.a mtime)

```
OK: pp/ppmail_qt (mtime 1777977638 >= lib 1777977218)
OK: pp/ppelas_qt (mtime 1777977639 >= lib 1777977218)
OK: pp/ppflui_qt (mtime 1777977641 >= lib 1777977218)
OK: pp/ppther_qt (mtime 1777977640 >= lib 1777977218)
OK: pp/ppnlse_qt (mtime 1777977642 >= lib 1777977218)
```

(Mtimes recorded post-build in the main clone; identical artifacts staged in
this per-agent worktree by file copy after the build completed — both clones
share the same `.git`. After staging into the worktree the freshness ratio is
preserved: all 5 pp/*_qt mtime equal the libxvueqt.a mtime.)

### ABI count

```
verify_abi: nm count: 58  header count: 58
```

(Stable Phase 7 baseline; no ABI drift introduced by Phase 8 Plan 1.)

### ISO date

`2026-05-05`

---

## Task 2 — Phase-7 deferred goldens (D-06, D-07)

(Filled in by Task 2 below.)

---

## Task 3 — Sweep harness + smoke probes

(Filled in by Task 3 below.)
