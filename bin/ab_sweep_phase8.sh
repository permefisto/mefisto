#!/usr/bin/env bash
# bin/ab_sweep_phase8.sh — Phase 8 top-level (case, backend, mode) sweep harness.
# Phase 8 Plan 1 Task 3 deliverable.
#
# Sources bin/phase8_case_batch_map.sh for per-case MODULE+BATCH+(prereq).
#
# Usage:
#   ab_sweep_phase8.sh --mode {x11|qt-1x|qt-2x|qt-omp} \
#                      --cases CSV \
#                      [--fuzz N] \
#                      [--out-dir DIR] \
#                      [--evidence-log PATH] \
#                      [--smoke-only] \
#                      [--baseline PATH]
#
# Per-cell behavior:
#   - x11    : invoke bin/ab_capture_x11.sh; OUT=${CASE}-x11.png
#   - qt-1x  : in-process Qt with QT_QPA_PLATFORM=offscreen + MEFISTO_BATCH_X11=1
#              + MEFISTO_QT_CAPTURE_PATH=${OUT_DIR}/${CASE}-qt-1x.png
#   - qt-2x  : qt-1x + QT_SCALE_FACTOR=2; OUT=${CASE}-qt-2x.png
#   - qt-omp : qt-1x + OMP_NUM_THREADS=8; OUT=${CASE}-qt-omp.png
#              (cavity2d skipped — no FLUIDER_OMP launcher per D-05)
#
# Flags:
#   --smoke-only   Skip the X11-baseline-compare step entirely. Used by
#                  Plan 1 Task 3 BEFORE any X11 baseline exists. Verdict
#                  emitted is `verdict=SMOKE`.
#   --baseline PATH  Override the default ${OUT_DIR}/${CASE}-x11.png LEFT
#                  operand for compare. PATH may contain the LITERAL token
#                  `${CASE}` (single-quoted by the caller so the shell
#                  does NOT interpolate). The harness performs PURE STRING
#                  SUBSTITUTION using bash parameter expansion
#                  `${BASELINE_PATH//\$\{CASE\}/$CURRENT_CASE}` — NOT eval.
#                  Plan 5 OMP cells use this to compare Qt-OMP captures
#                  against OMP-context-matched X11-OMP baselines.

set -euo pipefail

MODE=""
CASES_CSV=""
FUZZ=5
OUT_DIR=".planning/phases/08-ab-validation-on-testa-subset/evidence"
EVIDENCE_LOG=""
SMOKE_ONLY=0
BASELINE_PATH=""

while [ "$#" -gt 0 ]; do
    case "$1" in
        --mode)         MODE=$2; shift 2;;
        --cases)        CASES_CSV=$2; shift 2;;
        --fuzz)         FUZZ=$2; shift 2;;
        --out-dir)      OUT_DIR=$2; shift 2;;
        --evidence-log) EVIDENCE_LOG=$2; shift 2;;
        --smoke-only)   SMOKE_ONLY=1; shift;;
        --baseline)     BASELINE_PATH=$2; shift 2;;
        *) echo "ab_sweep_phase8: unknown arg $1" >&2; exit 2;;
    esac
done

case "$MODE" in
    x11|qt-1x|qt-2x|qt-omp) ;;
    *) echo "ab_sweep_phase8: --mode must be one of {x11,qt-1x,qt-2x,qt-omp} (got: $MODE)" >&2; exit 2;;
esac

if [ -z "$CASES_CSV" ]; then
    echo "ab_sweep_phase8: --cases CSV is required" >&2
    exit 2
fi

if [ -z "${MEFISTO:-}" ]; then
    echo "ab_sweep_phase8: MEFISTO env var must be set" >&2
    exit 2
fi
if [ -z "${MEFISTOX:-}" ]; then
    export MEFISTOX=/tmp/mefistox-phase8
fi

# Phase 9 Plan 9-09 (carry-forward #4): canonicalize OUT_DIR BEFORE pushd
# into PROJDIR. Without this, a relative --out-dir resolves under PROJDIR
# after pushd, so captures silently land in the wrong directory. Phase 8
# Plan 5 SUMMARY documented this as a deferred bug; Phase 9 cleans it.
# `-m` flag: do NOT require existence (mkdir -p creates dir next line).
OUT_DIR=$(realpath -m "$OUT_DIR")
mkdir -p "$OUT_DIR"
mkdir -p "$MEFISTOX"
[ -z "$EVIDENCE_LOG" ] && EVIDENCE_LOG="${OUT_DIR}/sweep-log-${MODE}.md"

# Source the empirical case-batch map (replaces fragile glob).
. "$MEFISTO/bin/phase8_case_batch_map.sh"

NOW_UTC=$(date -u +"%Y-%m-%d %H:%M:%S UTC")
{
    echo ""
    echo "## Sweep ${MODE} (${NOW_UTC})"
    echo ""
} >> "$EVIDENCE_LOG"

# Convert CSV to space-separated.
IFS=',' read -r -a CASES <<< "$CASES_CSV"

EXIT_CODE=0
for CURRENT_CASE in "${CASES[@]}"; do
    MODULE=$(phase8_case_module "$CURRENT_CASE")
    BATCH=$(phase8_case_batch "$CURRENT_CASE")
    PREREQ_MODULE=$(phase8_case_prereq_module "$CURRENT_CASE")
    PREREQ_BATCH=$(phase8_case_prereq_batch "$CURRENT_CASE")

    if [ -z "$MODULE" ] || [ -z "$BATCH" ]; then
        echo "case=$CURRENT_CASE mode=$MODE verdict=ERROR reason=\"case not in phase8_case_batch_map.sh\"" \
            | tee -a "$EVIDENCE_LOG"
        EXIT_CODE=3
        continue
    fi

    # qt-omp special-case: cavity2d has no FLUIDER_OMP launcher.
    if [ "$MODE" = "qt-omp" ] && [ "$CURRENT_CASE" = "cavity2d" ]; then
        echo "case=cavity2d mode=qt-omp verdict=N-A reason=\"no FLUIDER_OMP launcher (D-05)\"" \
            | tee -a "$EVIDENCE_LOG"
        continue
    fi

    # Workspace prep — fresh project tree + INITIER seed.
    PROJDIR="$MEFISTOX/$CURRENT_CASE"
    rm -rf "$PROJDIR"
    mkdir -p "$PROJDIR"
    cp -r "$MEFISTO/testa/$CURRENT_CASE/." "$PROJDIR/"
    pushd "$PROJDIR" >/dev/null
    # Phase 9 RETIRE-02: pp/ppinit is now Qt-linked (legacy X11 backend
    # retired); needs QT_QPA_PLATFORM=offscreen for headless CLI use.
    echo "$CURRENT_CASE" | env QT_QPA_PLATFORM=offscreen "$MEFISTO/pp/ppinit" >/dev/null 2>&1 || true

    # Optional MAILLER prereq (must run BEFORE the main module).
    if [ -n "$PREREQ_MODULE" ] && [ -n "$PREREQ_BATCH" ]; then
        env QT_QPA_PLATFORM=offscreen \
            MEFISTO_BATCH_X11=1 \
            MEFISTO_XVSOURIS_AUTOEXIT=1 \
            MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=500 \
            timeout 60 "$MEFISTO/pp/pp${PREREQ_MODULE}" "$PREREQ_BATCH" \
            >/dev/null 2>&1 || true
    fi

    OUT="${OUT_DIR}/${CURRENT_CASE}-${MODE}.png"
    rm -f "$OUT"

    case "$MODE" in
        x11)
            "$MEFISTO/bin/ab_capture_x11.sh" \
                --case "$CURRENT_CASE" \
                --module "$MODULE" \
                --batch "$BATCH" \
                --out "$OUT" >/dev/null 2>&1 || true
            ;;
        qt-1x)
            env QT_QPA_PLATFORM=offscreen \
                MEFISTO_BATCH_X11=1 \
                MEFISTO_QT_CAPTURE_PATH="$OUT" \
                MEFISTO_XVSOURIS_AUTOEXIT=1 \
                MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=500 \
                timeout 60 "$MEFISTO/pp/pp${MODULE}" "$BATCH" \
                >/dev/null 2>&1 || true
            ;;
        qt-2x)
            env QT_QPA_PLATFORM=offscreen \
                QT_SCALE_FACTOR=2 \
                MEFISTO_BATCH_X11=1 \
                MEFISTO_QT_CAPTURE_PATH="$OUT" \
                MEFISTO_XVSOURIS_AUTOEXIT=1 \
                MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=500 \
                timeout 60 "$MEFISTO/pp/pp${MODULE}" "$BATCH" \
                >/dev/null 2>&1 || true
            ;;
        qt-omp)
            env QT_QPA_PLATFORM=offscreen \
                OMP_NUM_THREADS=8 \
                MEFISTO_BATCH_X11=1 \
                MEFISTO_QT_CAPTURE_PATH="$OUT" \
                MEFISTO_XVSOURIS_AUTOEXIT=1 \
                MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=500 \
                timeout 60 "$MEFISTO/pp/pp${MODULE}" "$BATCH" \
                >/dev/null 2>&1 || true
            ;;
    esac

    popd >/dev/null

    SIZE=$(stat -c '%s' "$OUT" 2>/dev/null || echo 0)

    if [ "$SMOKE_ONLY" -eq 1 ]; then
        # No baseline compare in smoke mode.
        echo "case=$CURRENT_CASE mode=$MODE verdict=SMOKE out=$OUT size=$SIZE" \
            | tee -a "$EVIDENCE_LOG"
        continue
    fi

    if [ "$MODE" = "x11" ]; then
        # x11 mode IS the baseline; record the capture and continue.
        echo "case=$CURRENT_CASE mode=$MODE verdict=BASELINE out=$OUT size=$SIZE" \
            | tee -a "$EVIDENCE_LOG"
        continue
    fi

    # Resolve baseline path. ${CASE} literal-token substitution per BLOCKER-B
    # iter2: replace the literal characters `${CASE}` with the loop's
    # CURRENT_CASE — pure parameter expansion, NOT eval, NOT command sub.
    if [ -n "$BASELINE_PATH" ]; then
        RESOLVED_BASELINE="${BASELINE_PATH//\$\{CASE\}/$CURRENT_CASE}"
    else
        RESOLVED_BASELINE="${OUT_DIR}/${CURRENT_CASE}-x11.png"
    fi

    if [ ! -f "$RESOLVED_BASELINE" ]; then
        echo "case=$CURRENT_CASE mode=$MODE verdict=ERROR reason=\"baseline missing: $RESOLVED_BASELINE — run --mode x11 first or pass --baseline PATH\"" \
            | tee -a "$EVIDENCE_LOG"
        EXIT_CODE=4
        continue
    fi

    DIFF="${OUT_DIR}/${CURRENT_CASE}-${MODE}-diff.png"
    COMPARE_LINE=$("$MEFISTO/bin/ab_compare_pair.sh" \
        "$RESOLVED_BASELINE" "$OUT" "$DIFF" "$FUZZ" 2>/dev/null || true)

    echo "case=$CURRENT_CASE mode=$MODE $COMPARE_LINE" \
        | tee -a "$EVIDENCE_LOG"
done

exit $EXIT_CODE
