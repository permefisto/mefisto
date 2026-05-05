#!/bin/bash
# qt-capture.sh — run a MEFISTO Qt binary under QT_QPA_PLATFORM=offscreen
# and capture a PNG via the in-process MEFISTO_QT_CAPTURE_PATH hook.
#
# Usage:
#   bin/qt-capture.sh <driver_or_solver> <output.png> [data_file] [delay_ms=500] [project_dir=.]
#
# The capture uses the in-process QWidget::grab() hook inside xvfermer_,
# NOT an external screenshot tool. No X server is required — Qt's
# offscreen platform plugin handles the entire rendering pipeline. This
# is the reliable headless path when xcb-cursor0 or Xvfb is not available.
#
# For xvtest drivers pass the binary and the PNG path only:
#   bin/qt-capture.sh pp/ppxvtest2 /tmp/xvtest2_qt.png
#
# For solver runs pass the project_dir and the data file too:
#   bin/qt-capture.sh pp/ppmail /tmp/pan2d_qt.png pan2d.mesh /tmp/mefistox/pan2d
#
# (Pre-Phase-9 these examples named pp/ppxvtest2_qt / pp/ppmail_qt — RETIRE-02
# dropped the _qt suffix from pp/* binaries.)
#
# Requires: a built Qt target; no Xvfb or ImageMagick.

set -u

if [ $# -lt 2 ]; then
  echo "usage: $0 <binary> <output.png> [data_file] [delay_ms=500] [project_dir=.]" >&2
  exit 2
fi

BIN="$1"
OUTPUT="$2"
DATA_FILE="${3:-}"
DELAY_MS="${4:-500}"
PROJECT_DIR="${5:-.}"

: "${MEFISTO:?MEFISTO env var not set}"

if [ ! -x "$BIN" ]; then
  echo "qt-capture: binary not executable: $BIN" >&2
  exit 2
fi

# Absolute path for MEFISTO_QT_CAPTURE_PATH because the solver cd's into
# the project dir.
case "$OUTPUT" in
  /*) ABS_OUTPUT="$OUTPUT" ;;
  *)  ABS_OUTPUT="$PWD/$OUTPUT" ;;
esac

mkdir -p "$(dirname "$ABS_OUTPUT")" 2>/dev/null || true
rm -f "$ABS_OUTPUT"

case "$BIN" in
  /*) ABS_BIN="$BIN" ;;
  *)  ABS_BIN="$PWD/$BIN" ;;
esac

export QT_QPA_PLATFORM=offscreen
export MEFISTO_XVSOURIS_AUTOEXIT=1
export MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS="$DELAY_MS"
export MEFISTO_QT_CAPTURE_PATH="$ABS_OUTPUT"
export MEFISTO_BATCH_X11=1

: "${QT_CAPTURE_TIMEOUT:=300}"

if [ -n "$DATA_FILE" ]; then
  (cd "$PROJECT_DIR" && timeout "$QT_CAPTURE_TIMEOUT" "$ABS_BIN" "$DATA_FILE" \
     </dev/null >/tmp/qt-capture-driver.log 2>&1)
else
  timeout "$QT_CAPTURE_TIMEOUT" "$ABS_BIN" </dev/null \
    >/tmp/qt-capture-driver.log 2>&1
fi
RC=$?

if [ ! -s "$ABS_OUTPUT" ]; then
  echo "qt-capture: no PNG written to $ABS_OUTPUT (binary exit $RC — see /tmp/qt-capture-driver.log)" >&2
  exit 3
fi

echo "qt-capture: wrote $ABS_OUTPUT (binary $BIN exit $RC)"
exit 0
