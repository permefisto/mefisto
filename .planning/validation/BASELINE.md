# Validation Baseline — Phase 0 through Phase 8

> Immutable list of 5 canonical `testa/` cases (one per solver module) used
> as the A/B validation baseline for every phase of the xvue-qt migration.
> Write-once in Phase 0 (2026-04-10). Read-only for every subsequent phase
> through Phase 8.
>
> Per D-14, D-15, D-16 and REQUIREMENTS.md §BUILD-10.

## How this baseline is used

- **Phase 0:** Legacy X11 (`bin/cbl_tout` + `pp/pp*`) runs all 5 cases as the
  known-good starting reference (Task 2 below). Qt backend (`pp/pp*_qt`) is
  link-complete only — no graphics yet.
- **Phases 1-7:** End of each phase, all 5 cases are run through **both**
  backends. Results logged to `.planning/phases/{padded}-*/VALIDATION.md`.
- **Phase 8:** All 5 cases also run through `_OMP` variants and on a HiDPI
  display. The checklist in `.planning/phase-8/CHECKLIST.md` is the gate for
  declaring v1 shippable.
- **Phase 9:** After the one-release-cycle A/B window closes, legacy backend
  is retired; this baseline continues to anchor Qt-only regression testing.

## The 5 canonical cases

### 1. testa/pan2d — Mesher (ppmail)

| Field | Value |
|-------|-------|
| Project directory | `testa/pan2d/` |
| Solver module | `mail/` (mesher) |
| Launcher script | `MAILLER` (spawns `pp/ppmail` under legacy / `pp/ppmail_qt` under Qt) |
| Phase 0 executable | `pp/ppmail` (legacy, X11) |
| Later-phase executable | `pp/ppmail_qt` (Qt, Phase 1+) |
| Expected qualitative behavior | [Maintainer to describe — see Task 2] |
| Known-flaky touchpoints | [To be discovered — first-run notes below] |
| First-run notes (Phase 0, legacy) | [Filled in by Task 2 after manual run] |

### 2. testa/nafems_le1 — Elasticity (ppelas)

| Field | Value |
|-------|-------|
| Project directory | `testa/nafems_le1/` |
| Solver module | `elas/` (linear elasticity) |
| Launcher script | `ELASTICER` |
| Phase 0 executable | `pp/ppelas` |
| Later-phase executable | `pp/ppelas_qt` |
| Expected qualitative behavior | [Maintainer to describe — see Task 2] |
| Known-flaky touchpoints | [To be discovered] |
| First-run notes (Phase 0, legacy) | [Filled in by Task 2] |
| Reference | NAFEMS Linear Elastic benchmark LE1 (standard reference case) |

### 3. testa/cavity2d — Fluid (ppflui)

| Field | Value |
|-------|-------|
| Project directory | `testa/cavity2d/` |
| Solver module | `flui/` (fluid mechanics) |
| Launcher script | `FLUIDER` |
| Phase 0 executable | `pp/ppflui` |
| Later-phase executable | `pp/ppflui_qt` |
| Expected qualitative behavior | [Maintainer to describe — see Task 2] |
| Known-flaky touchpoints | [To be discovered] |
| First-run notes (Phase 0, legacy) | [Filled in by Task 2] |
| Reference | Classic lid-driven cavity flow |

### 4. testa/heat1d — Thermal (ppther)

| Field | Value |
|-------|-------|
| Project directory | `testa/heat1d/` |
| Solver module | `ther/` (thermal) |
| Launcher script | `THERMICER` |
| Phase 0 executable | `pp/ppther` |
| Later-phase executable | `pp/ppther_qt` |
| Expected qualitative behavior | [Maintainer to describe — see Task 2] |
| Known-flaky touchpoints | [To be discovered] |
| First-run notes (Phase 0, legacy) | [Filled in by Task 2] |
| Reference | 1D heat conduction (simplest diagnostic for thermal solver) |

### 5. testa/nlsecu — Nonlinear (ppnlse)

| Field | Value |
|-------|-------|
| Project directory | `testa/nlsecu/` |
| Solver module | `nlse/` (nonlinear solver) |
| Launcher script | `NLSER` |
| Phase 0 executable | `pp/ppnlse` |
| Later-phase executable | `pp/ppnlse_qt` |
| Expected qualitative behavior | [Maintainer to describe — see Task 2] |
| Known-flaky touchpoints | [To be discovered] |
| First-run notes (Phase 0, legacy) | [Filled in by Task 2] |
| Reference | Nonlinear solver canonical case |

## Amendment policy

This file is **write-once in Phase 0**. Amendments are only permitted if:

1. A case is discovered to be genuinely unrunnable on the legacy backend on a
   clean machine (environment issue, not a Phase 0 bug). The case is then
   replaced with another `testa/` case from the same solver module, and the
   replacement is logged here with the reason.
2. A new solver module is added to MEFISTO (out of scope for xvue-qt v1).

Beyond those two cases, the 5 baseline cases are immutable until Phase 9
retirement.

## Audit metadata

- **Created:** 2026-04-10 (Phase 0, Plan 04, Task 1)
- **First run through legacy X11:** [date filled in by Task 2]
- **Source of truth for existence check:** `ls testa/` on 2026-04-10
  confirmed all 5 project directories exist (see `.planning/phases/00-build-skeleton-abi-stubs/00-RESEARCH.md §"testa/ Baseline Verification"`)
