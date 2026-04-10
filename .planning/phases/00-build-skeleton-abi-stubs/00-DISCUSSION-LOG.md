# Phase 0: Build skeleton & ABI stubs - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions captured in CONTEXT.md — this log preserves the discussion.

**Date:** 2026-04-10
**Phase:** 00-build-skeleton-abi-stubs
**Mode:** discuss (standard, no advisor mode — no USER-PROFILE.md present)
**Areas discussed:** Directory layout, ABI header organization, Build-system strategy, Validation baseline, Stub diagnostic behavior

## Gray Areas Presented

### A. Directory layout inside `xvue/`
| Option | Description | Recommendation |
|--------|-------------|----------------|
| A1 | `xvue/qt/` subfolder — all new C++ isolated | ✓ Recommended |
| A2 | Flat inside `xvue/` — mix `.cpp`/`.h` with `.f`/`.c` | Rejected as messy for grep/blame |
| A3 | Top-level `xvue_qt/` sibling with `xvue/` | Rejected as diverging from tree convention |

### B. ABI header organization
| Option | Description | Recommendation |
|--------|-------------|----------------|
| B1 | Single header `xvue_qt_api.h` with all ~60 entry points | ✓ Recommended |
| B2 | Split by category (drawing, events, colors, pixmap, export, state) | Rejected as churn during bring-up |

### C. `bin/cbl_tout_qt` strategy
| Option | Description | Recommendation |
|--------|-------------|----------------|
| C1 | Clone-and-modify `bin/cbl_tout` → `bin/cbl_tout_qt` | ✓ Recommended |
| C2 | Parametrize with `BACKEND=qt` env var | Rejected — adds conditionals; harder to retire |

### D. Validation baseline — which 5 canonical `testa/` cases
| Option | Description | Recommendation |
|--------|-------------|----------------|
| D-rec | `pan2d`, `nafems_le1`, `cavity2d`, `heat1d`, `nlsecu` | ✓ Recommended (benchmark coverage per module) |
| D-alt | User picks alternative cases | Available |

### E. Stub diagnostic behavior
| Option | Description | Recommendation |
|--------|-------------|----------------|
| E1 | Print-once warning, then silent | ✓ Recommended |
| E2 | Completely silent no-ops | Rejected — silent regressions invisible |
| E3 | `assert(false)` — fail loudly | Rejected — breaks Phase 0 link-and-run sanity check |

## Claude's Discretion (presented, not discussed)

- ABI trailing-underscore macro copied verbatim from `xvuelc.c`
- `MefistoPoint` shim inline in the header with `static_assert` on size
- CMake custom target `verify_abi` running `nm` + symbol-count check
- `xvue/qt/README_COORDS.md` Y-axis convention audit output
- `XVUE_QT_ASSERT_MAIN_THREAD()` debug-only thread-affinity macro

## User Decisions

**Selected:** A1, B1, C1, D-rec, E1 — all recommended defaults.

No corrections, no alternative wording, no scope-creep suggestions, no deferred ideas raised interactively beyond the ones Claude surfaced in the "Deferred Ideas" discretionary list.

## Corrections Made

None — user accepted all five recommended defaults as-is.

## External Research

None — Phase 0 is infrastructure/plumbing; all evidence lives in `xvue/xvuelc.c`, `bin/cb*`, `.planning/codebase/`, and `.planning/research/`, which were already loaded before this discussion.

## Prior Phase Decisions Applied

None — Phase 0 is the first phase. Project-level decisions from PROJECT.md (Qt 6, CMake-scoped, parallel X11 build, byte-identical ABI) were already locked before discussion and are recorded in the CONTEXT.md `<decisions>` section as context, not revisited.

## Observations for Future Phases

- **User strongly prefers existing conventions.** Every choice (A1 isolation, B1 single header, C1 clone-and-modify) optimizes for minimal deviation from the current tree's one-purpose-per-file discipline. Future phases should default to "do what the existing code does" rather than "introduce a cleaner abstraction."
- **User accepted warn-once stubs (E1) over silent (E2) and hard-abort (E3).** Signal: wants diagnostics to be visible but non-disruptive. Apply the same calibration to future warning/error surfaces (e.g. Phase 6 dialog refusals).
- **User did not push back on the baseline `testa/` choices.** They are now locked for every subsequent phase's validation log. Changing the baseline later will invalidate prior validation records.
