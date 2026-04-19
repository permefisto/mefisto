#!/bin/sh
# xvue/qt/cmake/verify_icon_source.sh
# Phase 6.0 UI-SPEC §Iconography HiDPI: SVG only under resources/icons/.
# Plans 6.1..6.5 append SVG entries; this lint guards against accidental
# raster commits (PNG/JPG/BMP/GIF).
#
# Usage: verify_icon_source.sh <icons-dir>
#
# Exit 0 — clean (or directory does not exist yet — Plan 01 case).
# Exit 1 — at least one raster file found.

set -eu
ICONS_DIR="$1"

# Plan 01 ships an empty resources/xvue_icons.qrc; the icons/ subdir does
# not exist yet. Treat absence as success.
if [ ! -d "$ICONS_DIR" ]; then
    echo "verify_icon_source: $ICONS_DIR not present (empty Plan 01 state) — OK"
    exit 0
fi

HITS=$(find "$ICONS_DIR" -type f \
    \( -iname '*.png' -o -iname '*.jpg' -o -iname '*.jpeg' \
       -o -iname '*.bmp' -o -iname '*.gif' \) 2>/dev/null || true)

if [ -n "$HITS" ]; then
    echo "ERROR: UI-SPEC §Iconography violation — only SVG icons allowed under $ICONS_DIR:" >&2
    echo "$HITS" >&2
    exit 1
fi

echo "verify_icon_source: OK ($ICONS_DIR is SVG-clean)"
exit 0
