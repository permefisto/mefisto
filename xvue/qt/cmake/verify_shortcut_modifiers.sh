#!/bin/sh
# xvue/qt/cmake/verify_shortcut_modifiers.sh
# Phase 6.0 D-04 modifier-rule lint. Fails the build if any plain-character
# QShortcut/QKeySequence is found under the source tree (no Ctrl/Alt/Shift/F-key
# qualifier). Single-character bareword shortcuts collide with menu mnemonics
# and break the right-click canvas context lexicon.
#
# Usage: verify_shortcut_modifiers.sh <src-dir>
#
# Accept patterns:
#   QKeySequence("Ctrl+O"), "Ctrl+Shift+S", "F1", "F9", "Ctrl++", "Ctrl+-",
#   "Ctrl+0", "Ctrl+,", and any QKeySequence::StandardKey enum.
# Reject patterns:
#   QKeySequence("a"), "5", ";", "o" (bare single character).
#
# Exit 0 — clean.
# Exit 1 — at least one bareword found (build fails).

set -eu
SRC_DIR="$1"

# Pattern: QKeySequence(" + exactly one [a-zA-Z0-9;] + ").
# Then strip lines that contain a modifier token (Ctrl/Alt/Shift/Meta).
HITS=$(grep -rn 'QKeySequence("[a-zA-Z0-9;]")' "$SRC_DIR" 2>/dev/null \
    | grep -v 'Ctrl' | grep -v 'Alt' | grep -v 'Shift' | grep -v 'Meta' \
    || true)

if [ -n "$HITS" ]; then
    echo "ERROR: D-04 violation — bare-char QKeySequence found:" >&2
    echo "$HITS" >&2
    exit 1
fi

echo "verify_shortcut_modifiers: OK (no bareword shortcuts in $SRC_DIR)"
exit 0
