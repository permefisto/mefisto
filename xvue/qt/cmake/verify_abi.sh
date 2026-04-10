#!/bin/sh
# xvue/qt/cmake/verify_abi.sh
# Helper for the verify_abi CMake custom_target. Kept outside of the
# CMakeLists.txt COMMAND string so we don't have to fight the
# VERBATIM + GNU make escape rules for `$`.
#
# Usage: verify_abi.sh <libxvueqt.a path> <xvue_qt_api.h path>
#
# Exit 0  — nm symbol count matches header declaration count.
# Exit 1  — drift detected (CMake build fails).

set -eu

LIB="$1"
HDR="$2"

NM_COUNT=$(nm "$LIB" | grep -c ' T [a-zA-Z_][a-zA-Z0-9_]*_$' || true)
HDR_COUNT=$(grep -c '^[[:space:]]*\(void\|int\|float\|double\|long\|short\|unsigned\|void[[:space:]]*\*\)[[:space:]]*proc(' "$HDR" || true)

echo "verify_abi: nm count: $NM_COUNT  header count: $HDR_COUNT"

if [ "$NM_COUNT" != "$HDR_COUNT" ]; then
    echo "ERROR: ABI symbol count drift between $LIB and $HDR" >&2
    exit 1
fi

if [ "$NM_COUNT" != "57" ]; then
    echo "WARNING: expected 57 Fortran-facing entries (Planner Alert Option A), got $NM_COUNT" >&2
fi

exit 0
