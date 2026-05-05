#!/usr/bin/env bash
# bin/ab_capture_x11.sh — Phase 8 X11-side capture wrapper.
# Phase 8 Plan 1 Task 3 deliverable.
#
# Wraps the X11 capture pattern in 08-RESEARCH.md §"Pattern 1 / canonical
# X11-side":
#   - Set MEFISTO_XVFERMER_READY_FILE so xvfermer_ touches a sentinel file
#     just before destroying the display.
#   - Set MEFISTO_XVFERMER_HOLD_MS so xvfermer_ usleep()s long enough to
#     keep the rendered scene on-screen for `import -window root`.
#   - Set MEFISTO_XVSOURIS_AUTOEXIT=1 so blocking xvsouris_ reads return
#     a synthetic SPACE keypress on first iteration.
#
# Usage:  ab_capture_x11.sh --case CASE --module MODULE --batch BATCH
#                           --out OUT.png [--display :99] [--hold-ms 1500]
#
# WORKSPACE+CWD DISCIPLINE: caller MUST `cd $MEFISTOX/$CASE` BEFORE invoking
# this script. The legacy launcher pattern (bin/MAILLER, bin/ELASTICER, ...)
# enforces this with `direxec=$MEFISTOX/$nomproj; cd $direxec`. The script
# does NOT chdir on the caller's behalf — it inherits whatever cwd the
# caller has set up.
#
# This script wraps its body inside `xvfb-run --auto-servernum -s ...` only
# if no $DISPLAY is exported by the caller. Otherwise it reuses the caller's
# display.

set -euo pipefail

CASE=""
MODULE=""
BATCH=""
OUT=""
DISP=""
HOLD_MS=1500

while [ "$#" -gt 0 ]; do
    case "$1" in
        --case)     CASE=$2; shift 2;;
        --module)   MODULE=$2; shift 2;;
        --batch)    BATCH=$2; shift 2;;
        --out)      OUT=$2; shift 2;;
        --display)  DISP=$2; shift 2;;
        --hold-ms)  HOLD_MS=$2; shift 2;;
        *) echo "ab_capture_x11: unknown arg $1" >&2; exit 2;;
    esac
done

if [ -z "$CASE" ] || [ -z "$MODULE" ] || [ -z "$BATCH" ] || [ -z "$OUT" ]; then
    echo "usage: $0 --case CASE --module MODULE --batch BATCH --out OUT.png [--display :99] [--hold-ms 1500]" >&2
    exit 2
fi

if [ -z "${MEFISTO:-}" ]; then
    echo "ab_capture_x11: MEFISTO env var must be set" >&2
    exit 2
fi

# Decide whether we own the Xvfb session (no DISPLAY exported by caller) or
# we reuse the caller's session.
do_capture() {
    local READY_FILE="/tmp/ab_ready_${CASE}_$$"
    rm -f "$READY_FILE"

    # Background launch of pp/pp${MODULE}.
    env DISPLAY="${DISPLAY:-:99}" \
        MEFISTO_XVFERMER_READY_FILE="$READY_FILE" \
        MEFISTO_XVFERMER_HOLD_MS="$HOLD_MS" \
        MEFISTO_XVSOURIS_AUTOEXIT=1 \
        MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=500 \
      "$MEFISTO/pp/pp${MODULE}" "$BATCH" >/dev/null 2>&1 &
    local PID=$!

    # Overall 60s safety net so a hung case can't run away with the harness.
    ( sleep 60 && kill -TERM "$PID" 2>/dev/null ) &
    local KILLER=$!

    # Poll for the sentinel; bail if process dies first.
    while [ ! -f "$READY_FILE" ]; do
        sleep 0.1
        if [ ! -d "/proc/$PID" ]; then
            echo "ab_capture_x11: pid $PID exited before READY_FILE appeared" >&2
            break
        fi
    done

    # Capture root window — covers all decorations + the Mefisto canvas.
    DISPLAY="${DISPLAY:-:99}" import -window root "$OUT" 2>/dev/null || true

    # Wait for solver to wrap up; record its exit code.
    wait "$PID" 2>/dev/null || true
    local EXIT_CODE=$?
    kill "$KILLER" 2>/dev/null || true
    rm -f "$READY_FILE"

    local SIZE
    SIZE=$(stat -c '%s' "$OUT" 2>/dev/null || echo 0)
    echo "case=$CASE module=$MODULE backend=x11 out=$OUT exit=$EXIT_CODE size=$SIZE"
}

if [ -n "$DISP" ]; then
    DISPLAY="$DISP" do_capture
elif [ -z "${DISPLAY:-}" ]; then
    # Caller has no DISPLAY; we own the Xvfb session.
    xvfb-run --auto-servernum -s "-screen 0 1280x800x24" -- bash -c "
        export MEFISTO='$MEFISTO'
        $(declare -f do_capture)
        CASE='$CASE' MODULE='$MODULE' BATCH='$BATCH' OUT='$OUT' HOLD_MS='$HOLD_MS' do_capture
    "
else
    do_capture
fi
