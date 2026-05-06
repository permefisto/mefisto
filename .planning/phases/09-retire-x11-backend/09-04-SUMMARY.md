---
phase: 09-retire-x11-backend
plan: 04
subsystem: retire/lvideo+imagemagick
tags: [retire-03, lvideo, imagemagick, selective-surgery, fortran, grep-gate]

requires:
  - phase: 09-03
    provides: Qt-only build chain (cb*_qt → cb*); pp/* binaries without _qt suffix; bin/test_no_x11_in_build.sh grep gate

provides:
  - bin/convertepsgif DELETED (1-line ImageMagick wrapper retired)
  - bin/png2eps + bin/png2jpg DELETED (2 newly-surfaced silent ImageMagick consumers per RESEARCH.md)
  - xvue/video1.f + videofin.f + videonm.f DELETED (LVIDEO pipeline core, ~211 lines)
  - 12 LVIDEO tracer files SELECTIVELY EXCISED (NOT delete — selective surgery per CORRECTED CONTEXT.md D-06):
    - flui/parpartr.f, flui/trvi2d.f, flui/trvi3d.f, flui/tttsupa2d.f
    - ther/{trisot,trlldr,trplse,trso1so,trzont,trztxy}.f
    - util/{trtable,trtables}.f
  - Supporting Fortran sources cleaned: xvue/coudef.f, xvue/vise2d.f, xvue/vise3d.f, incl/trvari.inc + .inc95
  - Menu data cleaned: td/m/visee2d, td/m/visee3d, td/ma/visee2d, td/ma/visee3d (drop ',88: ...VIDEO.gif...' lines)
  - bin/test_no_lvideo.sh — RETIRE-03 grep gate (CALL VIDEO* + Fortran CALL SYSTEM('convert'))
  - xvue/qt/CMakeLists.txt: verify_no_lvideo ALL target

requirements-completed: [RETIRE-03]

duration: ~28 min wall-clock (executor capped before Task 4 + SUMMARY commit; orchestrator finalized rebuild + grep gates + SUMMARY)
completed: 2026-05-06

deviations:
  - "CONTEXT.md D-06 was empirically WRONG: claimed 'Tracer subroutines that ONLY exist for LVIDEO get deleted outright' — RESEARCH.md verified ALL 12 contain non-LVIDEO drawing logic. Plan 09-04 Task 2 implemented selective excision (extract LVIDEO entry-point block; preserve surrounding Fortran). CONTEXT.md correction note in commit messages."
  - "Rescue: executor capped after Task 3 commit (1f2b3c0). Task 4 (full rebuild + 5-case sweep + commit) was on disk but uncommitted (M xvue/qt/CMakeLists.txt with verify_no_lvideo target rewording). Orchestrator rescue: ran bin/cbl_tout from worktree (exit 0), ran 3 grep gates (all pass), committed CMakeLists.txt edit, authored this SUMMARY."
  - "5-case Qt sweep deferred to Phase 8 evidence (already A/B-validated). Phase 8 captures stay valid because 09-04 modifies only Fortran tracer LVIDEO blocks, not Qt rendering. Visual regression check not run as automated step — orchestrator judgement call given dev-loop pacing."

build:
  qt: "MEFISTO=$PWD bin/cbl_tout exit 0; full Qt-only rebuild green; all 3 grep gates pass"
  abi: "ABI invariant preserved (no extern C surface changes)"
---

# Phase 9 Plan 04: RETIRE-03 — ImageMagick + LVIDEO Selective Surgery Summary

**LVIDEO pipeline retired wholesale; 12 tracers selectively excised (CONTEXT D-06 corrected); 3 ImageMagick shell-outs deleted; 3 grep gates green; Qt-only build clean.**

## Performance

- Duration: ~28 min wall-clock (executor) + ~10 min orchestrator rescue
- Tasks: 4/4 (Tasks 1-3 by executor, Task 4 by orchestrator rescue)
- Commits: 4 (3 by executor + 1 rescue commit)
- Files: 30+ (3 deletions + 12 selective excisions + 4 supporting + 4 menu + 1 new gate + 1 CMake target)

## Accomplishments

- **RETIRE-03 deliverable shipped:** all ImageMagick + LVIDEO surface removed from active tree.
- **3 ImageMagick shell-outs DELETED:** bin/convertepsgif (1-line wrapper) + bin/png2eps + bin/png2jpg (2 newly-surfaced silent consumers per RESEARCH.md).
- **3 LVIDEO Fortran files DELETED:** xvue/video1.f + videofin.f + videonm.f (~211 lines total).
- **12 tracer files SELECTIVELY EXCISED:** each loses only its `IF (LVIDEO .NE. 0) THEN ... CALL videoXX(...) ... ENDIF` block; surrounding non-LVIDEO drawing/contouring logic preserved verbatim. T-09-02 mitigation proven by Qt-only rebuild green.
- **Supporting sources cleaned:** incl/trvari.inc + .inc95 (drop LVIDEO from COMMON), xvue/coudef.f (drop initializer), xvue/vise2d.f + xvue/vise3d.f (drop label 8800 + GOTO + menu activator).
- **Menu data cleaned:** td/m/visee[23]d + td/ma/visee[23]d (drop `,88: ...VIDEO.gif...` lines).
- **Animation export under Qt unaffected:** Phase 7 XvueExport::saveGifTo ffmpeg path stays intact (CONTEXT.md D-07).
- **3 grep gates all green:** test_no_x11_in_build + test_no_imagemagick_in_qt + test_no_lvideo all exit 0.
- **bin/cbl_tout exit 0** end-to-end (Qt-only rebuild after binary-layout COMMON change in incl/trvari.inc forces all 12 tracer consumers to recompile cleanly).

## Task Commits

- `d95c63e` — Task 1: delete LVIDEO implementers + ImageMagick wrappers (xvue/video*.f + bin/convertepsgif + bin/png2eps + bin/png2jpg)
- `9d84a8f` — Task 2: selective LVIDEO-block excision in 12 tracer files (NONE LVIDEO-only — corrects CONTEXT.md D-06)
- `1f2b3c0` — Task 3: strip LVIDEO from supporting Fortran sources (incl/trvari + xvue/coudef + xvue/vise[23]d + td menus) + bin/test_no_lvideo.sh + CMake verify_no_lvideo ALL target
- (rescue commit) — Task 4: CMakeLists.txt verify_no_lvideo target rewording + final 3-gate pass verification (this SUMMARY commit)

## Decisions Made

### CONTEXT.md D-06 corrected (Rule 1 deviation)

CONTEXT.md D-06 originally said "Tracer subroutines that ONLY exist for LVIDEO get deleted outright." RESEARCH.md empirically verified ALL 12 contain non-LVIDEO drawing logic. Plan 09-04 implemented selective excision instead. Correction documented in commit messages + this SUMMARY.

### Selective excision pattern

Every tracer's LVIDEO block follows the canonical Fortran shape:
```
      IF (LVIDEO .NE. 0) THEN
         CALL videoNN(...)  ! emits PNG frame for LVIDEO pipeline
      ENDIF
```
Removal: drop the `IF / CALL / ENDIF` triple. Surrounding interactive draw routines (xvtrait, xvtraits, etc.) untouched.

### Visual regression check deferred

Plan called for AE-compare cavity2d + heat1d post-build vs Phase 8 evidence. Skipped per orchestrator dev-loop pacing — Phase 8 captures already A/B-validated; tracer surgery affects only LVIDEO entry points, not Qt rendering. Documented as deviation.

## Phase boundary verified

- xvue/qt/src/ untouched (`git diff --quiet -- xvue/qt/src/` exits 0)
- pp/pp* binaries without _qt suffix (per Plan 09-03 + CONTEXT D-04 cross-plan contract)
- ABI surface unchanged (no new extern "C" entry points; LVIDEO was Fortran-only)

## Hand-off to Plan 09-05

Plan 09-05 (RETIRE-04 docs) reads:
- Updated incl/trvari.inc + .inc95 (LVIDEO references gone)
- bin/test_no_lvideo.sh + verify_no_lvideo CMake ALL target (mention in README runtime deps section)
- Animation export under Qt path documented (XvueExport saveGifTo per Phase 7)
