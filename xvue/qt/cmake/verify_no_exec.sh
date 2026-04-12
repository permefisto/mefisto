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

# Phase 3 D-19: palette-leak guard. XvueCanvas + xvue_qt_api.cpp must
# never read qApp->palette() or this->palette() — those would leak the
# system dark-mode QPalette into the backing pixmap (TEXT-06). Scoped
# grep (Pitfall 8) — window/app TUs may legitimately touch palette in
# Phase 6, so they are excluded.
PAL_MATCHES=$(grep -n -E 'qApp->palette|->palette\(\)' \
    "$SRC_DIR"/xvue_qt_canvas.cpp \
    "$SRC_DIR"/xvue_qt_canvas.h \
    "$SRC_DIR"/xvue_qt_api.cpp \
    2>/dev/null || true)

if [ -n "$PAL_MATCHES" ]; then
    echo "ERROR: TEXT-06/D-19 violation — qApp->palette / ->palette() forbidden in xvue_qt_canvas.* and xvue_qt_api.cpp" >&2
    echo "Offending matches:" >&2
    echo "$PAL_MATCHES" >&2
    exit 1
fi

echo "verify_no_exec: OK (palette-leak scan clean in canvas + api)"
exit 0
