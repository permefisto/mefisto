#!/bin/bash
# testa-capture.sh — run a MEFISTO solver in hybrid batch+X11 mode under
# Xvfb and capture a PNG.
#
# Usage:
#   bin/testa-capture.sh <project_dir> <solver> <data_file> <output.png>
#                        [hold_ms=1500] [display_num=99]
#
# The solver is launched with `MEFISTO_BATCH_X11=1` so prpr/pp<solver>.f
# keeps its batch-mode data-file-driven flow but opens the legacy X11
# window (INTERA=1). The xvfermer_ hold sentinel + import grab the final
# rendered state.
#
# The script does NOT call ppinit — the caller must have already
# INITIAlised the project directory (and chdir'ed into it).
#
# Requires: Xvfb, ImageMagick `import`, xdpyinfo.

set -u

if [ $# -lt 4 ]; then
  echo "usage: $0 <project_dir> <solver> <data_file> <output.png> [hold_ms=1500] [display=99]" >&2
  exit 2
fi

PROJECT_DIR="$1"
SOLVER="$2"
DATA_FILE="$3"
OUTPUT="$4"
HOLD_MS="${5:-1500}"
DISP_NUM="${6:-99}"

: "${MEFISTO:?MEFISTO env var not set}"

if [ ! -d "$PROJECT_DIR" ]; then
  echo "testa-capture: project dir not found: $PROJECT_DIR" >&2
  exit 2
fi

SOLVER_BIN="$MEFISTO/pp/$SOLVER"
if [ ! -x "$SOLVER_BIN" ]; then
  echo "testa-capture: solver not executable: $SOLVER_BIN" >&2
  exit 2
fi

if [ ! -f "$PROJECT_DIR/$DATA_FILE" ] && [ ! -L "$PROJECT_DIR/$DATA_FILE" ]; then
  echo "testa-capture: data file not found: $PROJECT_DIR/$DATA_FILE" >&2
  exit 2
fi

OUTDIR=$(dirname "$OUTPUT")
mkdir -p "$OUTDIR" 2>/dev/null || true

XVFB_PID=""
if ! xdpyinfo -display ":$DISP_NUM" >/dev/null 2>&1; then
  Xvfb ":$DISP_NUM" -screen 0 1280x800x24 >/tmp/testa-capture-xvfb.log 2>&1 &
  XVFB_PID=$!
  for i in 1 2 3 4 5 6 7 8 9 10; do
    xdpyinfo -display ":$DISP_NUM" >/dev/null 2>&1 && break
    sleep 0.2
  done
fi

READY_FILE=$(mktemp /tmp/testa-capture-ready.XXXXXX)
rm -f "$READY_FILE"

cleanup() {
  local rc=$?
  if [ -n "${SOLVER_PID:-}" ] && kill -0 "$SOLVER_PID" 2>/dev/null; then
    kill "$SOLVER_PID" 2>/dev/null
    wait "$SOLVER_PID" 2>/dev/null
  fi
  if [ -n "$XVFB_PID" ] && kill -0 "$XVFB_PID" 2>/dev/null; then
    kill "$XVFB_PID" 2>/dev/null
    wait "$XVFB_PID" 2>/dev/null
  fi
  rm -f "$READY_FILE"
  exit $rc
}
SOLVER_PID=""
trap cleanup EXIT INT TERM

export DISPLAY=":$DISP_NUM"
export MEFISTO_BATCH_X11=1
export MEFISTO_XVSOURIS_AUTOEXIT=1
export MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=50
export MEFISTO_XVFERMER_READY_FILE="$READY_FILE"
export MEFISTO_XVFERMER_HOLD_MS="$HOLD_MS"

: "${TESTA_CAPTURE_TIMEOUT:=300}"
(cd "$PROJECT_DIR" && timeout "$TESTA_CAPTURE_TIMEOUT" "$SOLVER_BIN" "$DATA_FILE" </dev/null \
   >/tmp/testa-capture-driver.log 2>&1) &
SOLVER_PID=$!

SNAP_DONE=0
while kill -0 "$SOLVER_PID" 2>/dev/null; do
  if [ -f "$READY_FILE" ] && [ $SNAP_DONE -eq 0 ]; then
    import -display ":$DISP_NUM" -window root "$OUTPUT"
    SNAP_DONE=1
  fi
  sleep 0.1
done

wait "$SOLVER_PID"
SOLVER_RC=$?
SOLVER_PID=""

if [ $SNAP_DONE -eq 0 ]; then
  echo "testa-capture: solver finished without firing xvfermer_ sentinel — no PNG" >&2
  echo "testa-capture: driver log at /tmp/testa-capture-driver.log" >&2
  exit 3
fi

if [ "$SOLVER_RC" -ne 0 ]; then
  echo "testa-capture: $SOLVER exited $SOLVER_RC (capture still written)" >&2
fi

echo "testa-capture: wrote $OUTPUT (solver $SOLVER exit $SOLVER_RC)"
exit 0
