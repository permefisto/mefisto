#!/bin/bash
# xvtest-capture.sh — run a MEFISTO X11 test driver under Xvfb and capture a PNG
#
# Usage:
#   bin/xvtest-capture.sh <driver> <output.png> [delay_ms] [display_num]
#
# Arguments:
#   driver       absolute or relative path to pp/ppxvtestN (or any MEFISTO
#                binary that calls XVSOURIS before its final XVFERMER)
#   output.png   path where the captured PNG should be written
#   delay_ms     how long xvsouris_ should wait in the autoexit short-circuit
#                before returning (default: 1200). Must be long enough for the
#                final render + the ImageMagick import call below.
#   display_num  Xvfb display number (default: 99)
#
# The script:
#   1. Starts Xvfb on :$display_num (1280x800x24) if not already running
#   2. Sets MEFISTO_XVSOURIS_AUTOEXIT so the driver does not block waiting for
#      keyboard input
#   3. Launches the driver in the background
#   4. Waits ~60% of delay_ms, then grabs the root window via `import`
#   5. Waits for the driver to exit and propagates its exit code
#
# Requires: Xvfb, ImageMagick's `import` or `magick import`, `pgrep`.

set -u

if [ $# -lt 2 ]; then
  echo "usage: $0 <driver> <output.png> [delay_ms=1200] [display_num=99]" >&2
  exit 2
fi

DRIVER="$1"
OUTPUT="$2"
DELAY_MS="${3:-1200}"
DISP_NUM="${4:-99}"

if [ ! -x "$DRIVER" ]; then
  echo "xvtest-capture: driver not executable: $DRIVER" >&2
  exit 2
fi

command -v Xvfb  >/dev/null 2>&1 || { echo "xvtest-capture: Xvfb not found" >&2; exit 2; }
command -v import >/dev/null 2>&1 || { echo "xvtest-capture: ImageMagick 'import' not found" >&2; exit 2; }

OUTDIR=$(dirname "$OUTPUT")
mkdir -p "$OUTDIR" 2>/dev/null || true

# Start Xvfb only if nothing is already listening on :$DISP_NUM
XVFB_PID=""
if ! xdpyinfo -display ":$DISP_NUM" >/dev/null 2>&1; then
  Xvfb ":$DISP_NUM" -screen 0 1280x800x24 >/tmp/xvtest-capture-xvfb.log 2>&1 &
  XVFB_PID=$!
  # wait up to 2s for Xvfb to be ready
  for i in 1 2 3 4 5 6 7 8 9 10; do
    if xdpyinfo -display ":$DISP_NUM" >/dev/null 2>&1; then break; fi
    sleep 0.2
  done
fi

cleanup() {
  local rc=$?
  if [ -n "$DRIVER_PID" ] && kill -0 "$DRIVER_PID" 2>/dev/null; then
    kill "$DRIVER_PID" 2>/dev/null
    wait "$DRIVER_PID" 2>/dev/null
  fi
  if [ -n "$XVFB_PID" ] && kill -0 "$XVFB_PID" 2>/dev/null; then
    kill "$XVFB_PID" 2>/dev/null
    wait "$XVFB_PID" 2>/dev/null
  fi
  exit $rc
}
DRIVER_PID=""
trap cleanup EXIT INT TERM

# Launch the driver headless. Two env-var mechanisms cooperate:
#
#   MEFISTO_XVSOURIS_AUTOEXIT + MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS (kept
#     small here, default 100ms) make every XVSOURIS call return a
#     synthetic keypress quickly so the driver streams through all its
#     drawing stages without waiting for user input.
#
#   MEFISTO_XVFERMER_READY_FILE + MEFISTO_XVFERMER_HOLD_MS create a
#     deterministic capture window at the very end: xvfermer_ touches the
#     sentinel file and sleeps delay_ms before destroying the X window.
#     We poll for the sentinel and call `import` inside that hold — the
#     screenshot is guaranteed to be the final rendered state regardless
#     of how many XVSOURIS stages the driver has.
READY_FILE=$(mktemp /tmp/xvtest-capture-ready.XXXXXX)
rm -f "$READY_FILE"

export DISPLAY=":$DISP_NUM"
export MEFISTO_XVSOURIS_AUTOEXIT=1
export MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=100
export MEFISTO_XVFERMER_READY_FILE="$READY_FILE"
export MEFISTO_XVFERMER_HOLD_MS="$DELAY_MS"

"$DRIVER" </dev/null >/tmp/xvtest-capture-driver.log 2>&1 &
DRIVER_PID=$!

# Poll for the sentinel file the Fortran driver's xvfermer_ touches just
# before it holds. Timeout after ~30s of total driver runtime.
WAITED_MS=0
POLL_MS=50
MAX_WAIT_MS=30000
while [ ! -f "$READY_FILE" ]; do
  if [ "$WAITED_MS" -ge "$MAX_WAIT_MS" ]; then
    echo "xvtest-capture: timed out waiting for xvfermer_ sentinel ($READY_FILE)" >&2
    rm -f "$READY_FILE"
    exit 3
  fi
  sleep 0.05
  WAITED_MS=$(( WAITED_MS + POLL_MS ))
  if ! kill -0 "$DRIVER_PID" 2>/dev/null; then
    # driver already exited — capture post-mortem may still show root if
    # the WM kept the last contents, but usually not. Try anyway.
    break
  fi
done

import -display ":$DISP_NUM" -window root "$OUTPUT"
IMPORT_RC=$?
rm -f "$READY_FILE"

wait "$DRIVER_PID"
DRIVER_RC=$?
DRIVER_PID=""

if [ "$IMPORT_RC" -ne 0 ]; then
  echo "xvtest-capture: import failed ($IMPORT_RC) capturing $OUTPUT" >&2
  exit "$IMPORT_RC"
fi
if [ "$DRIVER_RC" -ne 0 ]; then
  echo "xvtest-capture: driver $DRIVER exited $DRIVER_RC (see /tmp/xvtest-capture-driver.log)" >&2
  exit "$DRIVER_RC"
fi

echo "xvtest-capture: wrote $OUTPUT (driver exit 0)"
exit 0
