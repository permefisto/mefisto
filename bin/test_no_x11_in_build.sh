#!/bin/bash
# bin/test_no_x11_in_build.sh — Phase 9 Plan 9-03 (RETIRE-02)
#
# CI-style assertion: the post-Phase-9 build must NEVER reference X11.
# Scope: bin/cb* + bin/Makefile* + xvue/qt/ (the active build path).
# Out of scope: .planning/ (historical docs may reference X11 contextually),
# xvue/qt/build/ (CMake-generated, gets clobbered), worktree-*/ +
# .claude/worktrees/ (cross-tag work for Plan 9-08 by design contains
# X11 from v1.0-pre-retire).
#
# This gate clones the bin/test_no_imagemagick_in_qt.sh structure with
# the regex 'lX11|/usr/X11R6'.
#
# Exit codes:
#   0  — gate passes
#   1  — gate fails (X11 references found)
#   2  — environment misconfigured

set -e
if [ -z "$MEFISTO" ]; then
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
    REPO_ROOT="$(dirname "$SCRIPT_DIR")"
    if [ -d "$REPO_ROOT/xvue/qt" ]; then
        MEFISTO="$REPO_ROOT"
    else
        echo "test_no_x11_in_build.sh: \$MEFISTO unset and could not auto-detect from $SCRIPT_DIR" >&2
        exit 2
    fi
fi
cd "$MEFISTO"

# Run the grep on bin/cb* + xvue/qt/. Exclude .planning/ (historical),
# xvue/qt/build/ (generated), worktree-*/ + .claude/worktrees/
# (cross-tag artifact dirs).
#
# Allowlist (legitimate references to the legacy X11 bootstrap procedure):
#   - xvue/qt/tests/golden/scene01_driver.f — Fortran driver whose HEADER
#     COMMENT documents the cross-tag (v1.0-pre-retire) bootstrap recipe
#     for materializing the EPS golden (Plan 06 Task 3 procedure). The
#     procedure runs from a separate v1.0-pre-retire worktree per
#     Phase 9 D-04; the comment is reference documentation, not an
#     active build dependency.
MATCHES=$(grep -rln '\-lX11\|\-lXt\|/usr/X11R6\|/usr/X11R5' \
    bin/cb* bin/Makefile* xvue/qt/ 2>/dev/null \
    | grep -v 'xvue/qt/build/' \
    | grep -v '\.planning/' \
    | grep -v 'worktree-' \
    | grep -v '\.claude/worktrees' \
    | grep -v 'xvue/qt/tests/golden/scene01_driver\.f$' \
    || true)
if [ -n "$MATCHES" ]; then
    echo "FAIL: post-Phase-9 build path still references X11:" >&2
    echo "$MATCHES" >&2
    exit 1
fi
echo "OK: no X11 references in active build path"
exit 0
