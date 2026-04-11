#!/bin/sh
# xvue/qt/cmake/verify_no_exec.sh
# Helper for the verify_no_exec CMake custom_target (SHELL-03, D-10).
# Fails the build if QApplication::exec or qApp->exec() appears anywhere
# under the xvue/qt Qt layer source tree. Single enforcement point — no git
# pre-commit hook is installed (D-11: CMake guard is the sole enforcement).
#
# Usage: verify_no_exec.sh <src-dir> <include-dir>
#
# Exit 0 — no forbidden tokens found.
# Exit 1 — at least one match (build fails).

set -eu

SRC_DIR="$1"
INC_DIR="$2"

# Pattern 1: QApplication::exec (static call)
# Pattern 2: qApp->exec() (pointer call)
# Use grep -R -n so the error message includes file and line.
MATCHES=$(grep -R -n -E 'QApplication::exec|qApp->exec\(\)' \
    "$SRC_DIR" "$INC_DIR" 2>/dev/null || true)

if [ -n "$MATCHES" ]; then
    echo "ERROR: SHELL-03 violation — QApplication::exec / qApp->exec() forbidden in xvue/qt/" >&2
    echo "Offending matches:" >&2
    echo "$MATCHES" >&2
    exit 1
fi

echo "verify_no_exec: OK (no forbidden tokens in $SRC_DIR or $INC_DIR)"
exit 0
