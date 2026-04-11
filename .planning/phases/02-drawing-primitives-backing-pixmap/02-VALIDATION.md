---
phase: 02
slug: drawing-primitives-backing-pixmap
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-04-11
---

# Phase 02 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> See `02-RESEARCH.md` → "## Validation Architecture" for the authoritative design.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | MEFISTO legacy drivers (`prpr/xvtest*.f`) + shell-level build checks |
| **Config file** | none — uses `bin/cbl_tout` / `bin/cbl_tout_qt` build chain |
| **Quick run command** | `bin/cbxvue_qt && bin/verify_abi` |
| **Full suite command** | `bin/cbl_tout_qt && pp/ppxvtest0 && pp/ppxvtest1..4` (Qt backend) |
| **Estimated runtime** | ~60 seconds build + <10 seconds per driver |

The planner MUST refine these commands to whatever `02-RESEARCH.md` prescribes;
this file is a scaffold — the Validation Architecture section of RESEARCH.md is
the source of truth.

---

## Sampling Rate

- **After every task commit:** Run `bin/cbxvue_qt` (compile-only smoke)
- **After every plan wave:** Run `bin/verify_abi` + `pp/ppxvtest0` lifecycle driver
- **Before `/gsd-verify-work`:** Full `xvtest0..4` suite must be green under Qt backend
- **Max feedback latency:** ~60 seconds

---

## Per-Task Verification Map

*To be filled by the planner during plan generation. Each task in every
`02-NN-PLAN.md` must have a row here cross-referencing its `requirements`
and its automated verification command.*

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 02-NN-NN | NN   | W    | DRAW-XX     | —          | N/A (local mesh tool) | build/driver | `{command}` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

Per RESEARCH.md, three Wave 0 pre-implementation artifacts are required before
any DRAW-XX task can claim green status:

- [ ] `prpr/xvtest0.f` — extend with draw-coverage section between lifecycle cycles
      (or create `prpr/xvtest0_draw.f`) exercising every DRAW-01..09 primitive
- [ ] `xvue/qt/MEFISTO_POINT_AUDIT.md` — byte-layout audit document (D-31)
      enumerating all 13 Fortran call sites and confirming `INTEGER*2 (2,N)`
      matches `struct MefistoPoint { short x; short y; }`
- [ ] `xvue/qt/README_RESIZE.md` — documented convention for DRAW-09
      (backing pixmap reallocation on resize, content preservation rule)

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Visual parity xvtest1..4 vs X11 baseline | DRAW-01..08 | Antialiasing + white-only rendering make pixel diffs signal-free; image-diff golden validation deferred to Phase 3 | Run `pp/ppxvtest1` through `pp/ppxvtest4` under both `XVUE_BACKEND=x11` and Qt; observe side-by-side; confirm lines, polylines, fills, rectangles and ellipse arcs are visually indistinguishable |
| Resize preserves backing content | DRAW-09 | Requires interactive window manipulation during `SLEEP` hold | While `pp/ppxvtest0` (extended) is in its draw-hold SLEEP, drag the canvas window to a new size; confirm previously drawn shapes remain in place per the documented convention |

---

## Validation Sign-Off

- [ ] All tasks have automated build/driver verification or explicit Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers xvtest0 extension, MEFISTO_POINT_AUDIT.md, README_RESIZE.md
- [ ] Visual parity + resize listed as manual-only with clear instructions
- [ ] `bin/verify_abi` continues to pass (57 symbols, unchanged from Phase 1)
- [ ] `nyquist_compliant: true` set in frontmatter after planner fills per-task map

**Approval:** pending
