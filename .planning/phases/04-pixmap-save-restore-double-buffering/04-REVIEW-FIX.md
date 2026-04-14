---
phase: 04-pixmap-save-restore-double-buffering
fixed_at: 2026-04-14T00:00:00Z
review_path: .planning/phases/04-pixmap-save-restore-double-buffering/04-REVIEW.md
iteration: 1
findings_in_scope: 2
fixed: 2
skipped: 0
status: all_fixed
---

# Phase 4: Code Review Fix Report

**Fixed at:** 2026-04-14
**Source review:** `.planning/phases/04-pixmap-save-restore-double-buffering/04-REVIEW.md`
**Iteration:** 1

**Summary:**
- Findings in scope (Critical + Warning): 2
- Fixed: 2
- Skipped: 0
- Info findings (IN-01..IN-05): out of scope for this iteration

## Fixed Issues

### WR-01: Shell AE parser collapses multi-line stderr into a broken scalar

**Files modified:** `bin/xvtest0-pixmap-roundtrip.sh`
**Commit:** `23dd4df`
**Applied fix:** Replaced `magick compare ... 2>&1 || true` with a
stdout-dropping redirect (`2>&1 >/dev/null`) followed by an explicit `rc=$?`
check that distinguishes "images differ" (rc=1, handled by the COUNT check)
from "magick error" (rc>=2, reported as FAIL with the magick stderr). Changed
the count extractor from `echo "$AE" | awk '{print $1}'` (which prints $1 for
every stderr line) to `printf '%s\n' "$AE" | awk 'NR==1 {print $1; exit}'`,
which takes only the first line's leading token and so is immune to IM7
policy/delegate warning lines prepended before the AE count. Added a comment
documenting that `FAILED` is an intentional global accumulator (also addresses
the readability concern flagged by IN-04).

Verified:
- Tier 1: re-read the helper, fix text present, surrounding code intact.
- Tier 2: `bash -n` syntax check passed.
- Simulated-input check: old parser on `"@warning: ...\n0"` yields a
  multi-line string that never matches `"0"` (reproduces the bug); new parser
  yields `@warning:` on the polluted case and `0` on both the clean `"0"`
  case and the IM7-normalized `"0 (0)"` case.

Pending end-to-end verification: the user task asked for a re-run of
`bin/xvtest0-pixmap-roundtrip.sh` to confirm all 4 pairs still AE=0, but
`pp/ppxvtest0_qt` is not built in this sandbox (no `g++`/`gfortran` available;
`bin/cbl_tout_qt` fails at CMake `No CMAKE_CXX_COMPILER could be found`). The
harness change is pure bash, self-contained, and structurally/semantically
verified, but the end-to-end round-trip run remains to be performed by the
user (or the verifier phase) before Phase 4 is declared green.

### WR-02: `saved_canvas_` devicePixelRatio not refreshed on slot reuse

**Files modified:** `xvue/qt/src/xvue_qt_api.cpp`
**Commit:** `a04c26b`
**Applied fix:** In `xvue_qt_save_to_slot` (anonymous-namespace helper), moved
`setDevicePixelRatio(backing_->devicePixelRatio())` out of the lazy-realloc
`if` block and placed it unconditionally after the allocation check. The
allocation branch still creates a fresh `QPixmap(backing_->size())`, and the
reuse branch now re-syncs DPR every call. Added an inline comment explaining
why: Qt 6 can change `backing_->devicePixelRatio()` when the window moves
between monitors without changing `backing_->size()` (size is in device
pixels), so the reuse path must refresh DPR or `restore_from_slot`'s
`drawPixmap(0,0,...)` would misscale. The operation is cheap and idempotent,
and the 57-entry ABI surface is untouched.

Verified:
- Tier 1: re-read lines 86-112 of `xvue_qt_api.cpp`, fix text present,
  braces balanced, surrounding save/restore logic intact.
- Tier 2: not available — no C++ compiler installed in this sandbox
  (`g++`/`c++` not found; `bin/cbl_tout_qt` fails at CMake CXX-compiler
  detection). Fell back to Tier 1 only per verification strategy.

Pending end-to-end verification: the user task asked for a `bin/cbl_tout_qt`
build check after `xvue_qt_api.cpp` changes. This could not be performed in
the sandbox due to the missing toolchain and should be run by the user (or
the verifier phase) before Phase 4 is declared green. The change is a
one-line move of an already-correct call plus a comment, so risk is minimal.

## Skipped Issues

None — both in-scope findings were fixed.

---

_Fixed: 2026-04-14_
_Fixer: Claude (gsd-code-fixer)_
_Iteration: 1_
