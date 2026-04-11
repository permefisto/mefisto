---
phase: 02-drawing-primitives-backing-pixmap
plan: 01
subsystem: xvue-qt-abi-scaffolding
tags: [qt6, fortran-abi, drawing-primitives, wave-0, audit, docs]

requires:
  - phase: 01-window-shell-xvueapp-xvuewindow-xvuecanvas
    provides: "XvueApp/XvueWindow/XvueCanvas quartet + 57 warn-once extern C stubs + verify_abi/verify_no_exec guards"
provides:
  - "MEFISTO_POINT_AUDIT.md — all 29 live XVFACE/XVTRAITS/XVFACETRAITS Fortran callers verified INTEGER*2 (2,N) ABI-safe"
  - "README_RESIZE.md — top-left anchor resize-preserve convention locked in (DRAW-09)"
  - "Corrected xvarcellipse_/xvbordarcellipse_ ABI signatures (float* angles, not int*) — literal parity with xvuelc.c:2554/2616 per D-33"
  - "Extended prpr/xvtest0.f DRAW-01..09 coverage driver — every primitive + 3 pen styles exercised, Phase 1 reopen cycle preserved"
affects:
  - 02-02 (Wave 1: line/polyline/polygon/rect real bodies — reads ABI header float* angles)
  - 02-03 (Wave 2: arc/ellipse + pen styles real bodies — dereferences float* angles)
  - 02-04 (Wave 3: resizeEvent/backing — reads README_RESIZE.md convention)

tech-stack:
  added: []
  patterns:
    - "Wave 0 pre-implementation convention: audit docs + ABI corrections land separately from behavior changes so git bisect isolates regressions"
    - "MefistoPoint { short x; short y; } validated as Xlib XPoint byte-identical across every live Fortran caller"

key-files:
  created:
    - xvue/qt/MEFISTO_POINT_AUDIT.md
    - xvue/qt/README_RESIZE.md
    - .planning/phases/02-drawing-primitives-backing-pixmap/02-01-SUMMARY.md
  modified:
    - xvue/qt/include/xvue_qt_api.h
    - xvue/qt/src/xvue_qt_api.cpp
    - prpr/xvtest0.f

key-decisions:
  - "D-31 empirical validation: 29/29 live callers use INTEGER*2 (2,N); 3 CCC-commented sites excluded from live ABI"
  - "D-33 literal-parity corrected: xvarcellipse_/xvbordarcellipse_ now declared float* angle1/angle2 matching xvuelc.c byte-for-byte; previous int* was a latent Phase 0 bug"
  - "D-36 draw-coverage section placed BEFORE SLEEP(1) so visual gate sees primitives once Wave 1/2 bodies land"

patterns-established:
  - "Wave 0 artifacts: grep-able docs (MEFISTO_POINT_AUDIT.md, README_RESIZE.md) + ABI corrections, no behavior change"
  - "xvtest0.f is the progressive smoke driver: run after each Wave, watch warn-once lines disappear as real bodies land"

requirements-completed:
  - DRAW-03
  - DRAW-09

duration: 12min
completed: 2026-04-11
---

# Phase 2 Plan 01: Wave 0 ABI Scaffolding Summary

**Wave 0 pre-implementation landed: full Fortran call-site audit (29/29 ABI-safe), DRAW-09 resize convention doc, xvarcellipse_/xvbordarcellipse_ float* ABI correction, and extended prpr/xvtest0.f draw-coverage driver — all while libxvueqt.a stays warn-once and builds green at 57 symbols.**

## Performance

- **Duration:** ~12 min
- **Started:** 2026-04-11T17:47:00Z
- **Completed:** 2026-04-11T17:59:00Z
- **Tasks:** 4 / 4
- **Files modified:** 5 (3 created, 2 modified)

## Accomplishments

- **D-31 Fortran call-site audit produced and committed.** Grepped `xvue/ prpr/ mail/ elas/ flui/ ther/ nlse/ util/ reso/` for `CALL XVFACE | XVTRAITS | XVFACETRAITS`. Found 32 hits total: 29 live calls + 3 `CCC`-commented-out `XVFACETRAITS` references in `xvue/fap32d.f`. Every live caller verified with surrounding-context reads: 100% use `INTEGER*2 (2, N)` layouts (`XYPX`, `XYF`, `XYSF`, `T`, `XPOINTS`). Zero `INTEGER*4` / `INTEGER` / 8-byte offenders — MefistoPoint `{ short x; short y; }` (4 bytes, matches Xlib XPoint) is confirmed ABI-safe.
- **D-08 resize-preserve convention locked in `xvue/qt/README_RESIZE.md`.** The five invariants (DPR-aware reallocation, background fill, top-left `drawPixmap(0,0)` blit, grow shows background, shrink clipped by Qt) are now grep-able alongside the Phase 0 `xvue/README_COORDS.md`.
- **ABI signature bug (D-33) corrected in lockstep.** `xvue/qt/include/xvue_qt_api.h` and `xvue/qt/src/xvue_qt_api.cpp` both now declare `xvbordarcellipse_` / `xvarcellipse_` with `float *angle1, float *angle2` — byte-identical to `xvue/xvuelc.c:2554` and `:2616`. Stub bodies remain warn-once; real bodies land in Wave 2 (plan 02-03).
- **`prpr/xvtest0.f` extended with DRAW-01..09 coverage section (D-36).** Between the first `XVINITGRAPHIQUE` and the first `SLEEP(1)`, every primitive is exercised: `XVTRAIT`, `XVTRAITS` (3-pt polyline), `XVFACE` (4-pt polygon), all four rect variants, `XVARCELLIPSE` / `XVBORDARCELLIPSE` with `REAL*4` angles, `EFFACER`, `XVVOIR`, and pen styles 0/1/2 via `XVTYPETRAIT` + `XVEPAISSEUR`. Phase 1 reopen cycle (`SHELL-01/02/06`) preserved intact — `pp/ppxvtest0_qt` runs exit 0 and prints 13 warn-once lines confirming every new stub is reachable and link-resolved.

## Task Commits

1. **Task 1: Produce MEFISTO_POINT_AUDIT.md Fortran call-site audit (D-31)** — `e8accd6` (docs)
2. **Task 2: Create README_RESIZE.md documenting D-08 resize-preserve convention** — `4ee34a6` (docs)
3. **Task 3: Fix xvarcellipse_/xvbordarcellipse_ ABI signature (int*→float*) per D-33** — `1f9acf0` (fix)
4. **Task 4: Extend prpr/xvtest0.f with DRAW-01..09 coverage section (D-36)** — `b556037` (test)

## Files Created/Modified

- `xvue/qt/MEFISTO_POINT_AUDIT.md` — (new) D-31 Fortran call-site audit: 29/29 live callers verified INTEGER*2
- `xvue/qt/README_RESIZE.md` — (new) D-08 resize-preserve convention (top-left, no scale, no center)
- `xvue/qt/include/xvue_qt_api.h` — (modified) xvarcellipse_/xvbordarcellipse_ declarations flipped to `float *angle1, float *angle2` per D-33 literal parity
- `xvue/qt/src/xvue_qt_api.cpp` — (modified) stub bodies for xvarcellipse_/xvbordarcellipse_ updated in lockstep with header; warn-once convention preserved
- `prpr/xvtest0.f` — (modified) draw-coverage section inserted between first XVINITGRAPHIQUE and SLEEP(1); reopen cycle preserved

## ABI header edit (Task 3) — before / after

**Before** (`xvue/qt/include/xvue_qt_api.h` lines 141–142):
```c
void  proc(xvbordarcellipse)(int *x, int *y, int *width, int *height, int *a1, int *a2); // xvuelc.c:2554
void  proc(xvarcellipse)(int *x, int *y, int *width, int *height, int *a1, int *a2);     // xvuelc.c:2616
```

**After:**
```c
void  proc(xvbordarcellipse)(int *x, int *y, int *width, int *height,
                             float *angle1, float *angle2);              // xvuelc.c:2554 (D-33: literal)
void  proc(xvarcellipse)(int *x, int *y, int *width, int *height,
                         float *angle1, float *angle2);                  // xvuelc.c:2616 (D-33: literal)
```

Stub bodies in `xvue/qt/src/xvue_qt_api.cpp` flipped in the same commit; `(void)angle1; (void)angle2;` casts added to match the existing warn-once idiom.

## Audit result (Task 1)

- **Live call sites:** 29 (20 × XVFACE, 7 × XVTRAITS, 2 × XVFACETRAITS)
- **All declared:** `INTEGER*2 <name>(2, N)` — zero exceptions
- **Commented-out (CCC):** 3 × XVFACETRAITS in `xvue/fap32d.f` lines 134/153/238 — not in live ABI, listed separately for completeness
- **Hard blockers:** None
- **Verdict:** `MefistoPoint { short x; short y; }` is ABI-safe; Wave 1 and Wave 2 may dereference `MefistoPoint*` in the real bodies without risking a short↔int misalignment

## pp/ppxvtest0_qt run proof (Task 4)

Build via `bin/cbxvtest0_qt` exits 0. Run stdout excerpt:

```
 ===========================================
 Phase 1+2: cycle open/close + primitives
 ===========================================
 [xvtest0] premier appel XVINITGRAPHIQUE
Gtk-Message: ... Failed to load module "colorreload-gtk-module"
Gtk-Message: ... Failed to load module "window-decorations-gtk-module"
xvue-qt: stub xvtypetrait_ not implemented yet
xvue-qt: stub xvepaisseur_ not implemented yet
xvue-qt: stub xvtrait_ not implemented yet
xvue-qt: stub xvtraits_ not implemented yet
xvue-qt: stub xvface_ not implemented yet
xvue-qt: stub xvbordrectangle_ not implemented yet
xvue-qt: stub xvrectangle_ not implemented yet
xvue-qt: stub xvfrectangle_ not implemented yet
xvue-qt: stub xvfbordrectangle_ not implemented yet
xvue-qt: stub xvarcellipse_ not implemented yet
xvue-qt: stub xvbordarcellipse_ not implemented yet
xvue-qt: stub effacer_ not implemented yet
xvue-qt: stub xvvoir_ not implemented yet
 [xvtest0] premier appel XVFERMER
 [xvtest0] second appel XVINITGRAPHIQUE (reopen)
 [xvtest0] second appel XVFERMER

 [xvtest0] OK — cycle open/close/open/close + draws
EXIT: 0
```

13 warn-once lines (one per DRAW entry point), all 13 link-resolved, reopen cycle clean. Gtk-Message warnings are environment noise (GTK theme plugin probe via Qt's xdg-portal integration) — unrelated to xvue-qt and present since Phase 1.

## Build verification

- `bin/cbl_tout_qt` exit 0 — full Qt build green, including `verify_abi` (57 symbols, header count = nm count) and `verify_no_exec`
- `bin/cbl_tout` exit 0 — legacy X11 build unaffected (CLAUDE.md invariant "Compilation must never break" respected)
- `bin/cbxvtest0_qt` exit 0 — `pp/ppxvtest0_qt` created, ~87 KB, executable
- `pp/ppxvtest0_qt` exit 0 — warn-once proof of link resolution

## Decisions Made

- **D-31 audit table includes line numbers of both the `CALL` and the `INTEGER*2` declaration** — makes audit re-verification trivial for future auditors (they don't have to re-grep the file to find the declaration).
- **Commented-out `CCC XVFACETRAITS` sites tracked in a separate table** rather than omitted — preserves the audit's claim of "ALL CALL XVFACE/XVTRAITS/XVFACETRAITS sites examined" without inflating the live count.
- **Draw-coverage section placed BEFORE `SLEEP(1)`**, not after, so once Wave 1/2 bodies land the SLEEP(1) hold becomes a visual gate for the drawings (human eyeballs the shapes while Fortran blocks).
- **Pen-style test uses three `XVTRAIT` calls wrapped with `XVTYPETRAIT(0/1/2)`** — minimal footprint that still exercises every pen branch in Wave 2.

## Deviations from Plan

None — plan executed exactly as written. All four tasks ran atomically, each verification command passed on first try, both the Qt and legacy X11 builds stayed green.

## Issues Encountered

- `bin/cbl_tout_qt` initially failed with "MEFISTO=" / `expr: syntax error` — cause: `MEFISTO` env var was not set in the spawned shell. Fixed by exporting `MEFISTO=/home/drico/git/mefisto` (and the companion `MEFISTOX` / `PATH`) before invoking the build script. This is environmental, not a plan issue; CLAUDE.md already documents the required env vars.

## Next Phase Readiness

- **Wave 1 (02-02-PLAN.md) unblocked**: real bodies for `xvtrait_`, `xvtraits_`, `xvface_`, `xvfacetraits_`, rectangle symbols can now dereference `MefistoPoint*` with ABI confidence (audit), and the header signatures match `xvuelc.c`.
- **Wave 2 (02-03-PLAN.md) unblocked**: `xvarcellipse_` / `xvbordarcellipse_` bodies can now dereference `float *angle1` / `float *angle2` directly — no further header surgery.
- **Wave 3 (02-04-PLAN.md) unblocked**: `resizeEvent` implementation has a grep-able convention doc (`README_RESIZE.md`) it must conform to; reviewers have a single-source-of-truth for DRAW-09.
- **Progressive smoke gate in place**: `pp/ppxvtest0_qt` run after each future Wave will show warn-once lines disappearing as bodies land. The 13-line baseline is the Wave 0 signature.
- No blockers. `libxvueqt.a` still reports 57 ABI symbols. 100% of the Wave 0 plan's success criteria met.

## Self-Check: PASSED

Verified artifacts exist:
- FOUND: `xvue/qt/MEFISTO_POINT_AUDIT.md`
- FOUND: `xvue/qt/README_RESIZE.md`
- FOUND: `xvue/qt/include/xvue_qt_api.h` (edited)
- FOUND: `xvue/qt/src/xvue_qt_api.cpp` (edited)
- FOUND: `prpr/xvtest0.f` (extended)

Verified commits exist in `git log --oneline -5`:
- FOUND: `e8accd6` — Task 1 (audit)
- FOUND: `4ee34a6` — Task 2 (README_RESIZE)
- FOUND: `1f9acf0` — Task 3 (ABI fix)
- FOUND: `b556037` — Task 4 (xvtest0.f extended)

Verified builds green:
- `bin/cbl_tout_qt` — exit 0 (57 symbols, verify_abi PASS, verify_no_exec PASS)
- `bin/cbl_tout` — exit 0 (legacy X11 untouched)
- `pp/ppxvtest0_qt` — exit 0, 13 warn-once lines printed

---
*Phase: 02-drawing-primitives-backing-pixmap*
*Plan: 01 (Wave 0)*
*Completed: 2026-04-11*
