#!/usr/bin/env sh
# Phase 9 Plan 9-09 (Phase 8 D-09 carry-forward): Fail if libxvueqt.a is
# newer than ANY pp/pp* binary. Indicates the maintainer rebuilt the
# library but forgot to re-link one or more pp/* executables — the next
# test sweep would silently use a stale binary.
#
# Per RESEARCH §Open Question 4: this script is invoked at the END of
# bin/cbl_tout (NOT as a CMake ALL target — would chicken-and-egg with
# libxvueqt.a vs pp/* link order).
#
# NOTE: Plan 09-03 (Wave 2) renames cb*_qt -> cb* AND collapses pp/pp*_qt ->
# pp/pp* in cb script bodies AND in bin/ab_sweep_phase8.sh BEFORE this
# Wave-3 plan runs. Glob aligned to pp* (no _qt suffix).
#
# Args:
#   $1 = path to libxvueqt.a
#   $2 = path to pp/ directory (containing the linked pp/pp* binaries)
#
# Exit codes:
#   0 = all pp/pp* binaries are fresh (mtime >= libxvueqt.a mtime)
#   1 = at least one pp/pp* binary is stale (mtime < libxvueqt.a mtime),
#       OR libxvueqt.a is missing
set -eu
LIB=$1
PPDIR=$2

if [ ! -f "$LIB" ]; then
    echo "verify_pp_freshness: libxvueqt.a not found at $LIB" >&2
    exit 1
fi

LIB_MTIME=$(stat -c '%Y' "$LIB")
EXIT=0
# Orphaned legacy binaries not built by any current bin/cb* script.
# Listed here so the freshness check skips them — they're frozen artifacts
# from older builds, not Qt-port consumers of libxvueqt.a.
SKIPLIST=" pppoba "
for binary in "$PPDIR"/pp*; do
    [ ! -f "$binary" ] && continue
    [ -d "$binary" ] && continue
    base=$(basename "$binary")
    case "$SKIPLIST" in
        *" $base "*)
            echo "SKIP: $binary (orphaned legacy; not in current cb* build)"
            continue ;;
    esac
    BIN_MTIME=$(stat -c '%Y' "$binary")
    if [ "$BIN_MTIME" -lt "$LIB_MTIME" ]; then
        echo "FAIL: $binary mtime ($BIN_MTIME) < libxvueqt.a mtime ($LIB_MTIME) — rebuild stale binary" >&2
        EXIT=1
    else
        echo "OK: $binary mtime ($BIN_MTIME) >= libxvueqt.a mtime ($LIB_MTIME)"
    fi
done
exit $EXIT
