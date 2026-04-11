---
phase: 02
slug: drawing-primitives-backing-pixmap
status: complete
nyquist_compliant: true
wave_0_complete: true
created: 2026-04-11
completed: 2026-04-11
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
| **Full suite command** | `bin/cbl_tout_qt && pp/ppxvtest0_qt` |
| **Estimated runtime** | ~60 seconds build + <30 seconds driver (with extended holds) |

Phase 2 uses `pp/ppxvtest0_qt` as the single progressive smoke driver: every
Wave extends or flips its expected output (warn-once line count, rendered
geometry) and the human visual checkpoint in Wave 3 is the only non-headless
gate. Other `xvtest1..4` drivers are deferred to Phase 3 (they exercise text
and colormap paths that Phase 2 does not provide).

---

## Sampling Rate

- **After every task commit:** `bin/cbl_tout_qt` (or scoped `cmake --build xvue/qt/build`)
- **After every plan wave:** `pp/ppxvtest0_qt` (exit 0 + warn-once line delta check)
- **Before `/gsd-verify-work`:** clean rebuild (`rm -rf xvue/qt/build && bin/cbl_tout_qt`) + `bin/cbl_tout` regression guard + `pp/ppxvtest0_qt`
- **Max feedback latency:** ~60 seconds headless, ~30 seconds interactive (visual checkpoint)

---

## Per-Task Verification Map

| Task ID   | Plan | Wave | Requirement      | Threat Ref | Test Type     | Automated Command                                                                                                                                                                                                                      | Status |
|-----------|------|------|------------------|------------|---------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------|
| 02-01-01  | 01   | 0    | DRAW-03 (D-31)   | —          | audit         | `test -f xvue/qt/MEFISTO_POINT_AUDIT.md && grep -c 'INTEGER\*2' xvue/qt/MEFISTO_POINT_AUDIT.md` (≥29 live callers)                                                                                                                      | ✅     |
| 02-01-02  | 01   | 0    | DRAW-09 (D-08)   | —          | doc           | `test -f xvue/qt/README_RESIZE.md && grep -q 'top-left' xvue/qt/README_RESIZE.md`                                                                                                                                                      | ✅     |
| 02-01-03  | 01   | 0    | DRAW-05 (D-33)   | —          | build         | `grep -q 'float \*angle1' xvue/qt/include/xvue_qt_api.h && bin/cbl_tout_qt`                                                                                                                                                            | ✅     |
| 02-01-04  | 01   | 0    | DRAW-01..09      | —          | driver        | `bin/cbxvtest0_qt && pp/ppxvtest0_qt` (exit 0, 13 warn-once lines baseline)                                                                                                                                                            | ✅     |
| 02-02-01  | 02   | 1    | DRAW-01          | —          | build         | `grep -q 'backing_' xvue/qt/src/xvue_qt_state.h && grep -q 'applyPen' xvue/qt/src/xvue_qt_state.cpp && bin/cbl_tout_qt`                                                                                                                | ✅     |
| 02-02-02  | 02   | 1    | DRAW-01/08/09    | —          | build+driver  | `grep -q 'drawPixmap(0, 0' xvue/qt/src/xvue_qt_canvas.cpp && grep -q 'resizeEvent' xvue/qt/src/xvue_qt_canvas.cpp && bin/cbl_tout_qt && pp/ppxvtest0_qt`                                                                               | ✅     |
| 02-02-03  | 02   | 1    | DRAW-07          | —          | driver        | `grep -q 'fillRect' xvue/qt/src/xvue_qt_api.cpp && pp/ppxvtest0_qt` (effacer_/xvvoir_ warn-once gone → 11 lines)                                                                                                                        | ✅     |
| 02-03-01  | 03   | 2    | DRAW-02/03       | —          | driver        | `grep -q 'drawLine\|drawPolyline\|drawPolygon' xvue/qt/src/xvue_qt_api.cpp && pp/ppxvtest0_qt`                                                                                                                                          | ✅     |
| 02-03-02  | 03   | 2    | DRAW-04/06       | —          | driver        | `grep -q 'xvue_qt_draw_rect_common' xvue/qt/src/xvue_qt_api.cpp && pp/ppxvtest0_qt` (6 warn-once gone → 2 lines)                                                                                                                        | ✅     |
| 02-03-03  | 03   | 2    | DRAW-05          | —          | driver        | `grep -q 'drawPie(bbox' xvue/qt/src/xvue_qt_api.cpp && grep -q 'drawArc(bbox' xvue/qt/src/xvue_qt_api.cpp && ! grep -E ' \* 64\| \*64' xvue/qt/src/xvue_qt_api.cpp && pp/ppxvtest0_qt` (0 DRAW-XX warn-once)                            | ✅     |
| 02-04-01  | 04   | 3    | DRAW-01..09      | —          | build+driver  | `rm -rf xvue/qt/build && bin/cbl_tout_qt && bin/cbl_tout && pp/ppxvtest0_qt > /tmp/02-04-xvtest0.out && ! grep -qE 'xvue-qt: (xvtrait\|xvftrait\|xvtraits\|xvface\|xvfacetraits\|xv[fb]*rectangle\|xv[b]*arcellipse\|xvtypetrait\|xvepaisseur\|effacer\|xvvoir\|xvpxfenetre\|xvfond)_' /tmp/02-04-xvtest0.out` | ✅     |
| 02-04-02  | 04   | 3    | DRAW-02..09      | —          | manual        | Human visual checkpoint — see Manual-Only Verifications table (all 7 checks approved 2026-04-11)                                                                                                                                       | ✅     |
| 02-04-03  | 04   | 3    | DRAW-01..09      | —          | docs          | `grep -q 'status: complete' .planning/phases/02-drawing-primitives-backing-pixmap/02-VALIDATION.md && grep -q 'nyquist_compliant: true' .planning/phases/02-drawing-primitives-backing-pixmap/02-VALIDATION.md`                         | ✅     |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

**Sampling continuity check:** 13 tasks total, 11 with a headless automated
command (builds + greps + driver smoke), 1 manual visual checkpoint, 1
self-referential docs check. No 3 consecutive manual tasks; the single manual
gate is sandwiched between the automated clean-rebuild sweep (02-04-01) and
the automated docs close-out (02-04-03). Nyquist rate respected.

---

## Wave 0 Requirements

Per RESEARCH.md, three Wave 0 pre-implementation artifacts were required before
any DRAW-XX task could claim green status:

- [x] `prpr/xvtest0.f` — extended with draw-coverage section between lifecycle
      cycles, exercising every DRAW-01..09 primitive (plan 02-01 task 4,
      commit `b556037`; plan 02-04 extended holds, commit `a82247b`)
- [x] `xvue/qt/MEFISTO_POINT_AUDIT.md` — byte-layout audit document (D-31)
      enumerating all 29 live Fortran call sites and confirming
      `INTEGER*2 (2,N)` matches `struct MefistoPoint { short x; short y; }`
      (plan 02-01 task 1, commit `e8accd6`)
- [x] `xvue/qt/README_RESIZE.md` — documented convention for DRAW-09
      (backing pixmap reallocation on resize, top-left content preservation
      rule) (plan 02-01 task 2, commit `4ee34a6`)

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions | Status |
|----------|-------------|------------|-------------------|--------|
| Visual geometry present on canvas | DRAW-02..06 | Requires a live display session; pixel-diff golden deferred to Phase 3 | Run `pp/ppxvtest0_qt` on an X11/Wayland session. During the first SLEEP(15) hold, confirm the canvas shows: horizontal white line, 3-point polyline, filled 4-point polygon, four rectangles side-by-side, filled ellipse pie slice, outlined ellipse arc. | ✅ Approved 2026-04-11 |
| Resize preserves backing content | DRAW-09 | Requires interactive window manipulation during `SLEEP` hold | While `pp/ppxvtest0_qt` is in the first SLEEP(15) hold, drag the window edge to resize larger. Previously drawn geometry MUST remain anchored at top-left, new area filled background black — per `xvue/qt/README_RESIZE.md`. | ✅ Approved 2026-04-11 |
| Effacer clears mid-sequence | DRAW-07 | Visual confirmation the pre-effacer scene is gone after `CALL EFFACER` | During the second SLEEP(10) hold, confirm ONLY the post-effacer dashed lines are visible (no rectangles, no polygons, no pie slice). Implicitly verified by the earlier Wave 1 shorter-hold run where only post-effacer content was visible during the single SLEEP hold. | ✅ Approved 2026-04-11 (implicit) |

Checks 6 (clean reopen cycle, no `QApplication: there can only be one`
assertion) and 7 (exit 0, zero warn-once lines for DRAW-01..09) are covered
headlessly by 02-04-01 (clean-rebuild regression sweep) and do not require a
display — they were verified automatically.

---

## Validation Sign-Off

- [x] All tasks have automated build/driver verification or explicit Wave 0 dependencies
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 covers xvtest0 extension, MEFISTO_POINT_AUDIT.md, README_RESIZE.md
- [x] Visual parity + resize listed as manual-only with clear instructions
- [x] `bin/verify_abi` continues to pass (57 symbols, unchanged from Phase 1)
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** approved 2026-04-11 by human visual checkpoint (plan 02-04,
task 2) — all 7 visual checks passed. Checks 1 (geometry present), 2 (pen
style variety), 3 (antialiasing smoothness), and 5 (resize preserves
top-left) confirmed directly against `pp/ppxvtest0_qt` on a live display
session (driver commit `a82247b`, binary rebuilt 20:29 local). Check 4
(effacer clears mid-sequence) implicitly verified by the earlier shorter-hold
run where only the post-effacer dashed lines were visible during the single
SLEEP hold. Checks 6 (clean reopen cycle, no `QApplication` assertion) and 7
(exit 0, zero DRAW-XX warn-once lines) verified headlessly by the plan 02-04
task 1 clean-rebuild regression sweep.
