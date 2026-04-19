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

# Count Fortran-facing public text symbols ending in underscore, excluding
# C++ mangled names (which start with _Z and can end with _ via suffixes like
# S1_ that the name-mangler emits for repeated type references). Phase 5 added
# XvueEventBridge::waitForEvent(WaitMode, int*, int*) whose mangled name
# matches the bare [a-zA-Z_][a-zA-Z0-9_]*_$ pattern accidentally.
NM_COUNT=$(nm "$LIB" | grep ' T [a-zA-Z_][a-zA-Z0-9_]*_$' | grep -vc ' T _Z' || true)
HDR_COUNT=$(grep -c '^[[:space:]]*\(void\|int\|float\|double\|long\|short\|unsigned\|void[[:space:]]*\*\)[[:space:]]*proc(' "$HDR" || true)

echo "verify_abi: nm count: $NM_COUNT  header count: $HDR_COUNT"

if [ "$NM_COUNT" != "$HDR_COUNT" ]; then
    echo "ERROR: ABI symbol count drift between $LIB and $HDR" >&2
    exit 1
fi

if [ "$NM_COUNT" != "58" ]; then
    echo "WARNING: expected 58 Fortran-facing entries (Phase 6.0 adds xvue_module_init_), got $NM_COUNT" >&2
fi

exit 0
