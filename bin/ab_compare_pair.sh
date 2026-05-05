#!/usr/bin/env bash
# ORDER: argc → fuzz → identify → compare
# This contract is load-bearing — Plan 1 Task 3 verify gate `bash bin/ab_compare_pair.sh
# /tmp/_pair_a.png /tmp/_pair_b.png /tmp/_diff.png 50 2>&1 | grep -q 'fuzz must be in'`
# depends on FUZZ validation happening BEFORE identify is called on (non-existent) PNGs.
# Re-ordering steps WILL break the verify gate.
#
# bin/ab_compare_pair.sh — Phase 8 sweep tolerance-band gate wrapper.
# Phase 8 Plan 1 Task 3 deliverable.
#
# Usage:  ab_compare_pair.sh A.png B.png DIFF.png [FUZZ%]
#
# Arguments:
#   A.png    LEFT operand  (typically the X11 baseline capture)
#   B.png    RIGHT operand (typically the Qt capture being validated)
#   DIFF.png OUT path for the visual diff image
#   FUZZ%    integer percent in [1,30] (default 5 per D-02). NOT 100% — see
#            08-RESEARCH.md Pitfall 3 (fuzz=100% silently passes everything).
#
# Behavior:
#   1. Validate argc >= 3.
#   2. Validate FUZZ in [1,30] (BEFORE any file probe).
#   3. identify dimensions of A and B.
#   4. If dims differ, resample B to A's dims via point sampling (no AA
#      introduced). Emit `RESAMPLED: B (DB) → BR (DA)` to stderr.
#   5. Run `compare -metric AE -fuzz FUZZ%`. Capture exit code.
#   6. Emit single machine-parseable status line on stdout:
#      `ae=N total=M pct=P% verdict={PASS|CHECK} diff=DIFF resampled={yes|no}`

set -euo pipefail

# 3) argc check
if [ "$#" -lt 3 ]; then
    echo "usage: $0 A.png B.png DIFF.png [FUZZ%]" >&2
    exit 2
fi

A=$1
B=$2
DIFF=$3
FUZZ=${4:-5}

# 5) FUZZ validation BEFORE any file probe.
#    Verify gate exercises this path with FUZZ=50 on non-existent inputs;
#    FUZZ check MUST fire before identify is called.
if ! [[ "$FUZZ" =~ ^[0-9]+$ ]]; then
    echo "ab_compare_pair: fuzz must be in [1,30] (got: $FUZZ)" >&2
    exit 2
fi
if [ "$FUZZ" -lt 1 ] || [ "$FUZZ" -gt 30 ]; then
    echo "ab_compare_pair: fuzz must be in [1,30] (got: $FUZZ)" >&2
    exit 2
fi

# 6) identify dim probe (only after fuzz check passed).
DA=$(identify -format "%wx%h" "$A")
DB=$(identify -format "%wx%h" "$B")

# 7) Dimension guard. -filter point => no AA pixels introduced by resample.
RESAMPLED="no"
if [ "$DA" != "$DB" ]; then
    BR="${B%.png}-resampled-to-${DA}.png"
    convert "$B" -filter point -resize "${DA}!" "$BR"
    echo "RESAMPLED: $B ($DB) → $BR ($DA)" >&2
    B="$BR"
    RESAMPLED="yes"
fi

# 8) compare invocation (writes metric to stderr).
EXIT=0
AE_LINE=$(compare -metric AE -fuzz "${FUZZ}%" "$A" "$B" "$DIFF" 2>&1) || EXIT=$?
AE=$(echo "$AE_LINE" | awk '{print $1}')

# Sanity-fix AE if compare reported anything non-numeric.
if ! [[ "$AE" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    AE=0
fi

TOTAL_PX=$(identify -format "%[fx:w*h]" "$A")
AE_PCT=$(awk "BEGIN { printf \"%.4f\", $AE / $TOTAL_PX * 100 }")

# `compare` exit codes: 0 = within tolerance (PASS), 1 = differ (CHECK),
# 2 = error.  We translate 0 → PASS and anything else → CHECK so the
# harness can emit a deterministic verdict.
if [ "$EXIT" -eq 0 ]; then
    VERDICT="PASS"
else
    VERDICT="CHECK"
fi

# 9) Final stdout line, machine-parseable for ab_sweep_phase8.sh.
echo "ae=$AE total=$TOTAL_PX pct=$AE_PCT% verdict=$VERDICT diff=$DIFF resampled=$RESAMPLED"
