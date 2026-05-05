#!/bin/bash
# xvtest0-pixmap-roundtrip.sh — Phase 4 pixmap save/restore round-trip test.
#
# Invokes pp/ppxvtest0 6 times under bin/qt-capture.sh, one per scene
# (P4_CTRL, P4_SAVERESTORE, P4_MEMPX_SAVERESTORE, P4_BG, P4_EFFACEMEMPX,
#  P4_FENETREMEMPX), then runs 4 pairwise `magick compare -metric AE`
# comparisons. Exit 0 iff all 4 pairs are pixel-identical (AE == 0).
#
# Covers: PIXMAP-01 (no-op flip ops), PIXMAP-02 (sauvefenetre/restaurefenetre),
# PIXMAP-03a (sauvemempx/restauremempx), PIXMAP-03b (effacemempx == background).
# PIXMAP-04 (interactive cavity2d rubber-band) is DEFERRED to Phase 5 per
# CONTEXT.md D-18.

set -u
set -o pipefail

SELF="$(basename "$0")"

# -- Preflight probes --------------------------------------------------
if ! command -v magick >/dev/null 2>&1; then
    echo "${SELF}: magick (ImageMagick) not found — install imagemagick" >&2
    exit 2
fi

if [ ! -x pp/ppxvtest0 ]; then
    echo "${SELF}: pp/ppxvtest0 not built — run bin/cbxvtest0 first" >&2
    exit 2
fi

if [ ! -x bin/qt-capture.sh ]; then
    echo "${SELF}: bin/qt-capture.sh missing or not executable" >&2
    exit 2
fi

: "${MEFISTO:?${SELF}: MEFISTO env var not set}"

# -- Capture 6 scenes --------------------------------------------------
SCENES="P4_CTRL P4_SAVERESTORE P4_MEMPX_SAVERESTORE P4_BG P4_EFFACEMEMPX P4_FENETREMEMPX"

for SCENE in $SCENES; do
    # /tmp/p4_ctrl.png, /tmp/p4_saverestore.png, etc.
    LC=$(echo "$SCENE" | tr 'A-Z' 'a-z')
    OUT="/tmp/${LC}.png"
    rm -f "$OUT"
    echo "== ${SELF}: capturing scene=${SCENE} -> ${OUT}"
    if ! MEFISTO_XVTEST0_SCENE="$SCENE" bin/qt-capture.sh pp/ppxvtest0 "$OUT"; then
        echo "${SELF}: FAILED to capture scene ${SCENE}" >&2
        exit 3
    fi
    if [ ! -s "$OUT" ]; then
        echo "${SELF}: empty capture for scene ${SCENE}" >&2
        exit 3
    fi
done

# -- Pairwise compare --------------------------------------------------
# Helper: compare_ae <label> <req> <A.png> <B.png>
# NOTE: FAILED is intentionally a global accumulator mutated from inside
# the function (not declared `local`). See IN-04 in 04-REVIEW.md.
FAILED=0
compare_ae() {
    local LABEL="$1" REQ="$2" A="$3" B="$4"
    local AE rc COUNT
    # WR-01: magick compare -metric AE writes the AE count to stderr. We
    # redirect stderr into $AE but DROP stdout (null: diff image). Using
    # `|| true` previously hid the IM7 exit code — rc=1 means "images
    # differ", rc>=2 means "magick error" (missing file, bad format, …);
    # we must distinguish them. Also, IM7 may prepend policy/delegate
    # warning lines, so parse only the FIRST line's leading integer.
    AE=$(magick compare -metric AE "$A" "$B" null: 2>&1 >/dev/null)
    rc=$?
    if [ "$rc" -ge 2 ]; then
        echo "FAIL  [${REQ}]  ${LABEL}  (magick error rc=${rc}: ${AE}, A=${A}, B=${B})"
        FAILED=$((FAILED + 1))
        return
    fi
    COUNT=$(printf '%s\n' "$AE" | awk 'NR==1 {print $1; exit}')
    if [ "$COUNT" = "0" ]; then
        echo "PASS  [${REQ}]  ${LABEL}  (AE=0)"
    else
        echo "FAIL  [${REQ}]  ${LABEL}  (AE=${AE}, A=${A}, B=${B})"
        FAILED=$((FAILED + 1))
    fi
}

echo
echo "== ${SELF}: pairwise magick compare -metric AE"
compare_ae "sauvefenetre/restaurefenetre round-trip == ctrl"  "PIXMAP-02"  "/tmp/p4_ctrl.png"  "/tmp/p4_saverestore.png"
compare_ae "sauvemempx/restauremempx round-trip == ctrl"      "PIXMAP-03a" "/tmp/p4_ctrl.png"  "/tmp/p4_mempx_saverestore.png"
compare_ae "effacemempx == background"                        "PIXMAP-03b" "/tmp/p4_bg.png"    "/tmp/p4_effacemempx.png"
compare_ae "fenetremempx/mempxfenetre no-op == ctrl"          "PIXMAP-01"  "/tmp/p4_ctrl.png"  "/tmp/p4_fenetremempx.png"

echo
if [ "$FAILED" -eq 0 ]; then
    echo "${SELF}: ALL 4 round-trip pairs PASS — Phase 4 green (PIXMAP-01, -02, -03; PIXMAP-04 deferred to Phase 5)"
    exit 0
else
    echo "${SELF}: ${FAILED} round-trip pair(s) FAILED — see /tmp/p4_*.png for post-mortem"
    exit 1
fi
