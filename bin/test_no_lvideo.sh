#!/bin/bash
# bin/test_no_lvideo.sh -- Phase 9 Plan 9-04 (RETIRE-03)
#
# CI-style assertion: the post-Phase-9 active source tree must NEVER
# reference the legacy LVIDEO pipeline (CALL VIDEO1/FIN/NM in Fortran)
# or the ImageMagick shell-out (CALL SYSTEM('convert ...')).
#
# Scope: Fortran sources under xvue/, flui/, ther/, util/, prpr/.
# Out of scope: .planning/, xvue/qt/, build/, worktree-*/. Phase 7
# XvueExport::saveGifTo (ffmpeg) is the supported animation export
# (CONTEXT.md D-07).
#
# Exit codes:
#   0  -- gate passes
#   1  -- gate fails (one or more references found; details printed)
#   2  -- environment misconfigured ($MEFISTO unset and auto-detect failed)

set -e

if [ -z "$MEFISTO" ]; then
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
    REPO_ROOT="$(dirname "$SCRIPT_DIR")"
    if [ -d "$REPO_ROOT/xvue/qt" ]; then
        MEFISTO="$REPO_ROOT"
    else
        echo "test_no_lvideo.sh: \$MEFISTO unset and auto-detect failed" >&2
        exit 2
    fi
fi

cd "$MEFISTO"

# Layer 1: any CALL VIDEO* (uppercase, fixed-form column 7+)
LVIDEO_HITS=$(grep -rln 'CALL[[:space:]]*VIDEO\(1\|FIN\|NM\)' \
    xvue/*.f flui/ ther/ util/ prpr/ 2>/dev/null \
    | grep -v '\.planning/' | grep -v 'xvue/qt/' | grep -v 'worktree-' \
    || true)
if [ -n "$LVIDEO_HITS" ]; then
    echo "FAIL: residual CALL VIDEO* in active tree:" >&2
    echo "$LVIDEO_HITS" >&2
    exit 1
fi

# Layer 2: any CALL SYSTEM('convert' -- Fortran shell-out
CONVERT_HITS=$(grep -rln "CALL[[:space:]]*SYSTEM[[:space:]]*([[:space:]]*'convert" \
    xvue/*.f flui/ ther/ util/ prpr/ 2>/dev/null \
    | grep -v '\.planning/' | grep -v 'xvue/qt/' | grep -v 'worktree-' \
    || true)
if [ -n "$CONVERT_HITS" ]; then
    echo "FAIL: residual CALL SYSTEM('convert ...') in active tree:" >&2
    echo "$CONVERT_HITS" >&2
    exit 1
fi

echo "OK: no LVIDEO and no Fortran convert shell-outs"
exit 0
