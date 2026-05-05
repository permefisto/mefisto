#!/usr/bin/env bash
# Phase 8 Plan 1 Task 0 — empirical batch-file map (per BLOCKER #1, #3 iter1 revision).
#
# Sourced by bin/ab_sweep_phase8.sh:
#   . "$MEFISTO/bin/phase8_case_batch_map.sh"
#
# Each entry exports per-case:
#   PHASE8_CASE_${CASE}_MODULE        — pp/pp${MODULE}_qt module suffix
#   PHASE8_CASE_${CASE}_BATCH         — argv[1] passed to pp/pp${MODULE}_qt
#   PHASE8_CASE_${CASE}_PREREQ_MODULE — optional; module to run BEFORE main module
#   PHASE8_CASE_${CASE}_PREREQ_BATCH  — optional; batch file for prereq run
#
# WORKSPACE DISCIPLINE: caller MUST `cd $MEFISTOX/$CASE` BEFORE running pp/*.
# Every case requires `INITIER` (echo "$CASE" | $MEFISTO/pp/ppinit) to seed
# the workspace MS files first; this is enforced by the harness.
#
# Empirically derived 2026-05-05; evidence under
# .planning/phases/08-ab-validation-on-testa-subset/evidence/00-case-batch-map.md.

# pan2d — mesher (no MAILLER prereq because pan2d IS the mesher case).
PHASE8_CASE_pan2d_MODULE=mail
PHASE8_CASE_pan2d_BATCH=pan2d.mesh

# nafems_le1 — elasticity NAFEMS LE1; prereq mail with nafems_le1.mesh.
PHASE8_CASE_nafems_le1_MODULE=elas
PHASE8_CASE_nafems_le1_BATCH=nafems_le1.elas
PHASE8_CASE_nafems_le1_PREREQ_MODULE=mail
PHASE8_CASE_nafems_le1_PREREQ_BATCH=nafems_le1.mesh

# cavity2d — fluid (Stokes/Navier-Stokes); prereq mail with cavity2d.meshbf.
PHASE8_CASE_cavity2d_MODULE=flui
PHASE8_CASE_cavity2d_BATCH=cavity2d.stoke56cr
PHASE8_CASE_cavity2d_PREREQ_MODULE=mail
PHASE8_CASE_cavity2d_PREREQ_BATCH=cavity2d.meshbf

# heat1d — thermal; prereq mail with heat1d.mesh.
PHASE8_CASE_heat1d_MODULE=ther
PHASE8_CASE_heat1d_BATCH=heat1d.heat
PHASE8_CASE_heat1d_PREREQ_MODULE=mail
PHASE8_CASE_heat1d_PREREQ_BATCH=heat1d.mesh

# nlsecu — nonlinear (NLSE on a cube); prereq mail with nlsecu.meshq2.
# NB: nlsecu.iexrr requests Final TIME=20, Step=0.01 → 2000 time steps;
# the case is genuinely too long-running to complete in 60s on this hardware.
# Plan 2 will need a longer per-case timeout for nlsecu (or substitute).
PHASE8_CASE_nlsecu_MODULE=nlse
PHASE8_CASE_nlsecu_BATCH=nlsecu.iexrr
PHASE8_CASE_nlsecu_PREREQ_MODULE=mail
PHASE8_CASE_nlsecu_PREREQ_BATCH=nlsecu.meshq2

# Helper: emit the batch filename for a case (used by harness).
phase8_case_batch() {
    local case=$1
    local var="PHASE8_CASE_${case}_BATCH"
    echo "${!var}"
}

# Helper: emit the module suffix for a case.
phase8_case_module() {
    local case=$1
    local var="PHASE8_CASE_${case}_MODULE"
    echo "${!var}"
}

# Helper: emit the prereq module suffix for a case (empty if none).
phase8_case_prereq_module() {
    local case=$1
    local var="PHASE8_CASE_${case}_PREREQ_MODULE"
    echo "${!var:-}"
}

# Helper: emit the prereq batch filename for a case (empty if none).
phase8_case_prereq_batch() {
    local case=$1
    local var="PHASE8_CASE_${case}_PREREQ_BATCH"
    echo "${!var:-}"
}
