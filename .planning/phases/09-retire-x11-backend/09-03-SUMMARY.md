---
phase: 09-retire-x11-backend
plan: 03
subsystem: build/graphics-backend
tags: [retire-x11, libx11-strip, build-rename, qt-only, pp-binary-suffix-collapse, harness-coherency, cmake-grep-gate]

# Dependency graph
requires:
  - phase: 09-retire-x11-backend
    provides: "v1.0-pre-retire tag (rollback artifact); 09-01-AUDIT-BASELINE.md (empirical reference counts); 09-02 deletion of xvuelc.c + ccxvue (no Fortran ghost references); A/B window closed 2026-05-06"
provides:
  - "32 -lX11 / X11R6 references stripped (35 files deleted; 4 lines collapsed in bin/ab_sweep_phase8.sh)"
  - "Single Qt-only build entry: bin/cbl_tout (renamed from bin/cbl_tout_qt via git mv)"
  - "11 cb*_qt scripts renamed to cb* (history preserved)"
  - "10 nompp=pp/pp<x>_qt body lines collapsed to nompp=pp/pp<x> in renamed cb* scripts"
  - "Phase 9 pp/* binary naming policy locked: pp/pp{init,mail,elas,flui,ther,nlse,xvtest{0..4}} (no _qt suffix)"
  - "bin/test_no_x11_in_build.sh + verify_no_x11_in_build CMake ALL target — RETIRE-02 grep gate, mirrors verify_no_imagemagick_in_qt"
  - "bin/ab_sweep_phase8.sh harness aligned with policy lock (4 invocation sites updated; pp/ppinit invocation also wrapped with QT_QPA_PLATFORM=offscreen for Qt-link CLI use)"
  - "New bin/cbinit (Qt-only INITIER build) — closes the dispatcher gap created by deleting legacy cbinit (Plan 9-02 already broke it; this plan replaces it Qt-side)"
  - "Empirical proof bin/cbl_tout exits 0; 4/5 testa baselines capture (nlsecu carries to 9-06)"
affects: [09-04-retire-imagemagick-lvideo, 09-05-retire-docs-deps, 09-06-ppnlse-deadlock, 09-07-phase7-goldens, 09-08-harness-fixes]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Mechanical-deletion-with-atomic-commit (RESEARCH §Pattern 1)"
    - "Rename-via-git-mv-preserves-history (RESEARCH §Pitfall 4 / OQ2)"
    - "Build-time grep gate as CMake ALL target (clones Phase 7 EXPORT-06 verify_no_imagemagick_in_qt pattern)"
    - "Pre-sweep harness assertion (BLOCKER #4 add-on) before invoking qt-1x sweep"

key-files:
  created:
    - ".planning/phases/09-retire-x11-backend/09-03-RENAME-IMPACT.md"
    - ".planning/phases/09-retire-x11-backend/09-03-SUMMARY.md"
    - "bin/test_no_x11_in_build.sh (50 lines, 0o755) — RETIRE-02 grep gate"
    - "bin/cbinit (60 lines, 0o755) — Qt-only INITIER build (Rule 2 auto-add)"
  modified:
    - "bin/cbl_tout (renamed via git mv from bin/cbl_tout_qt; dispatcher list collapsed cb*_qt -> cb*; 6 echo prefixes 'cbl_tout_qt:' -> 'cbl_tout:'; +1 cbinit invocation)"
    - "bin/cbmail / cbelas / cbflui / cbther / cbnlse (5 renamed cb*_qt scripts; nompp body line collapsed)"
    - "bin/cbxvtest{0,1,2,3,4} (5 renamed cb*_qt scripts; nompp body line collapsed)"
    - "bin/ab_sweep_phase8.sh (4 lines: pp/pp${MODULE}_qt + pp/pp${PREREQ_MODULE}_qt -> pp/pp{}; ppinit wrapped with QT_QPA_PLATFORM=offscreen)"
    - "bin/phase8_case_batch_map.sh (2 comment lines updated)"
    - "bin/qt-capture.sh (2 usage examples + 1 attribution comment)"
    - "bin/xvtest0-pixmap-roundtrip.sh (3 occurrences: header + guard + capture invocation)"
    - "xvue/qt/CMakeLists.txt (2 comment fixes + appended verify_no_x11_in_build ALL target)"
    - "xvue/qt/README.md (3 active references to bin/cbl_tout_qt + 1 pp/ppmail_qt example updated)"
  deleted:
    - "bin/cbinit / cbmail / cbelas / cbflui / cbther / cbnlse (legacy active-build, 6 files)"
    - "bin/cbxvtest{0,1,2,3,4} (legacy active-build, 5 files)"
    - "bin/cbgadap / cbgelas / cbgflui / cbginit / cbginit1 / cbgmail / cbgnlse / cbgpara / cbgparaddd / cbgpoba / cbgther (out-of-active-build per OQ1, 11 files)"
    - "bin/cbpoba / cbpara / cbadap / cbpppl / cbonde / cbbrezfort2d / cbbrezfort3d (out-of-active-build per OQ1, 7 files)"
    - "bin/Makefile / MakefileIBM / MakefileMefisto (HP-UX/IBM relics, 3 files)"
    - "bin/avecMOTIF / instal_2disk.hp / instal_2disk.src (Motif/HP-UX dual-disk relics, 3 files)"
    - "bin/cbl_tout_qt (renamed to bin/cbl_tout via git mv)"
    - "bin/cb{mail,elas,flui,ther,nlse,xvtest0,xvtest1,xvtest2,xvtest3,xvtest4}_qt (10 _qt counterparts renamed to drop suffix)"
  renamed:
    - "bin/cbl_tout_qt -> bin/cbl_tout"
    - "bin/cbmail_qt -> bin/cbmail"
    - "bin/cbelas_qt -> bin/cbelas"
    - "bin/cbflui_qt -> bin/cbflui"
    - "bin/cbther_qt -> bin/cbther"
    - "bin/cbnlse_qt -> bin/cbnlse"
    - "bin/cbxvtest0_qt -> bin/cbxvtest0"
    - "bin/cbxvtest1_qt -> bin/cbxvtest1"
    - "bin/cbxvtest2_qt -> bin/cbxvtest2"
    - "bin/cbxvtest3_qt -> bin/cbxvtest3"
    - "bin/cbxvtest4_qt -> bin/cbxvtest4"

key-decisions:
  - "Tasks 2 + 3 folded into ONE atomic commit (b4af705): the legacy cb* deletions and the cb*_qt renames are inseparable per the plan's `git rm bin/cbl_tout && git mv bin/cbl_tout_qt bin/cbl_tout` recipe (RESEARCH §Pitfall 4 + OQ2 cleanest history)"
  - "Pp _qt suffix policy LOCKED: pp/* binaries ship without _qt; cb script bodies + bin/ab_sweep_phase8.sh + harness helper scripts all reference pp/pp<name> (BLOCKER #2/#3/#4)"
  - "Phase-6 source comments in xvue/qt/src/*_actions.cpp + tests/test_*_menu.cpp referencing pp/pp{elas,flui,nlse,ther}_qt classified as historical (Plan Task 1 step 4 (b)); NOT updated"
  - "Three active harness/example scripts (bin/phase8_case_batch_map.sh, bin/qt-capture.sh, bin/xvtest0-pixmap-roundtrip.sh) folded into Task 4 scope via Rule 2 auto-add (the harness is incoherent if these still reference pp/pp*_qt nonexistent binaries)"
  - "bin/cbinit was missing from the Qt build chain (cbinit_qt never existed pre-Phase-9); Rule 2 auto-add: created Qt-only bin/cbinit modelled on bin/cbmail (Task 6) so cbl_tout produces pp/ppinit"
  - "pp/ppinit (Qt-linked) needs QT_QPA_PLATFORM=offscreen for headless CLI use; bin/ab_sweep_phase8.sh:122 ppinit invocation wrapped accordingly (Rule 1 / Rule 2 auto-add — Plan 9-02 SUMMARY's same Rule 3 carry-forward)"
  - "Allowlisted xvue/qt/tests/golden/scene01_driver.f in test_no_x11_in_build.sh — its HEADER COMMENT documents the cross-tag (v1.0-pre-retire) bootstrap procedure for the Phase-7 EPS golden (D-04 cross-tag operation)"
  - "Build artifacts (mail/lib, elas/lib, etc., incl/homdir.inc) intentionally NOT staged in any plan-9-03 commit — they are generated by cbl_tout, not source changes; same policy as Plans 9-01 / 9-02"
  - "Sweep run with --smoke-only (Rule 3 auto-fix; no X11 baselines exist post-deletion to compare against; same auto-fix as Plans 9-01 / 9-02 documented)"

patterns-established:
  - "Plan 9 RETIRE-NN scripts model: clone the bin/cbmail (Qt-only) link recipe verbatim for any new Fortran-driver Qt-link target (gfortran ... mail/lib util/lib xvue/lib -Lxvue/qt/build -lxvueqt $QT_LIBS -lstdc++ -lgfortran)"
  - "RETIRE-NN grep gates: clone bin/test_no_imagemagick_in_qt.sh structure, change regex, add allowlist for cross-tag-bootstrap documentation files"

requirements-completed: ["RETIRE-02"]

# Metrics
duration: ~38 min
completed: 2026-05-05
---

# Phase 9 Plan 03: RETIRE-02 libX11 + linker lines + cbl_tout_qt rename + pp _qt suffix collapse Summary

**Stripped 32 -lX11 / -lXt / /usr/X11R6 linker references from the active build path via 35 file deletions and 4 sweep-harness edits, renamed 11 cb*_qt scripts via `git mv` to drop the `_qt` suffix, collapsed 10 nompp=pp/pp<x>_qt body references and 4 sweep dispatch sites to drop the `_qt` suffix from pp/* binaries (BLOCKER #2-#4 policy lock), added a verify_no_x11_in_build CMake ALL target + bin/test_no_x11_in_build.sh grep gate (modelled on Phase 7 EXPORT-06), added a Qt-only bin/cbinit (Rule 2 auto-add), and validated end-to-end with a clean `bin/cbl_tout` build (exit 0) + 4/5 testa baselines captured via the renamed harness.**

## Performance

- **Duration:** ~38 min (executor wall-clock)
- **Started:** 2026-05-05T22:37:01Z
- **Completed:** 2026-05-05T23:15:02Z
- **Tasks:** 6 (all autonomous; no checkpoints)
- **Commits:** 5 task commits + this SUMMARY-finalization commit (next)
- **Files added:** 4 (RENAME-IMPACT.md, SUMMARY.md, test_no_x11_in_build.sh, cbinit)
- **Files renamed (git mv):** 11 (cbl_tout_qt + 5 module cb*_qt + 5 cbxvtest*_qt)
- **Files deleted:** 35 (legacy cb*, Makefile*, Motif/HP-UX relics)
- **Files edited:** 9 (cbl_tout body, 10 renamed scripts' nompp lines, ab_sweep_phase8.sh, 3 helper harness scripts, xvue/qt/CMakeLists.txt, xvue/qt/README.md)

## Accomplishments

- **32 lib-X11 references gone.** `grep -rln 'lX11\|X11R6' bin/cb* bin/Makefile*` returns 0 matches across all 24 surviving bin/cb* scripts (Makefile* deleted entirely).
- **Single Qt-only build entry.** `bin/cbl_tout` is the canonical entrypoint; `bin/cbl_tout_qt` no longer exists. `git mv` preserved file history (the rename traces back to the original cbl_tout_qt creation).
- **11 cb*_qt scripts renamed in lockstep.** `git mv` on each; preserved exec bit. The dispatcher loop in `bin/cbl_tout` was hand-edited (Edit tool) to invoke `cb<module>` (no _qt) — preserves xvue/qt/ directory references inside the same script.
- **10 nompp body lines collapsed.** `nompp=pp/pp<name>_qt` -> `nompp=pp/pp<name>` in every renamed cb* script. Concrete proof: `grep '^nompp=pp/ppmail$' bin/cbmail` (line 31, drift from baseline line 4 documented in 09-03-RENAME-IMPACT.md but content correct).
- **Phase 9 pp binary naming policy locked.** All 11 pp/* outputs ship without `_qt` suffix: ppinit, ppmail, ppelas, ppflui, ppther, ppnlse, ppxvtest{0..4}. Empirical: `! ls pp/pp*_qt` returns no matches; `pp/ppmail` (etc.) exist with 1.1-7.4 MB sizes.
- **bin/ab_sweep_phase8.sh harness aligned.** Lines 130 (PREREQ_MODULE), 151 (qt-1x), 161 (qt-2x), 171 (qt-omp): each `pp/pp${...}_qt` invocation collapsed to `pp/pp${...}` (BLOCKER #4). Pre-Qt-sweep assertion runs green: `grep -q 'pp/pp${MODULE}' && ! grep -q 'pp/pp${MODULE}_qt'`.
- **Three additional active harness scripts updated.** bin/phase8_case_batch_map.sh (2 doc lines), bin/qt-capture.sh (2 examples + attribution comment), bin/xvtest0-pixmap-roundtrip.sh (3 occurrences). Folded in via Rule 2 auto-add — they were on the active dispatch chain and would have been incoherent otherwise.
- **RETIRE-02 grep gate live.** New `bin/test_no_x11_in_build.sh` (50 lines, 0o755) plus `verify_no_x11_in_build` ALL target in `xvue/qt/CMakeLists.txt`. Mirrors Phase 7's `verify_no_imagemagick_in_qt` (same shape: 2-line COMMAND chain, DEPENDS xvueqt, VERBATIM). Negative-test confirmed: injecting `-lX11` into bin/cbmail makes the gate fail with "FAIL: post-Phase-9 build path still references X11"; revert restores green.
- **Allowlist for cross-tag bootstrap doc.** `xvue/qt/tests/golden/scene01_driver.f` HEADER COMMENT contains the legacy X11 link recipe (`-L/usr/X11R6/lib -lX11 -lXt`) — this is reference documentation for the Phase 7 Plan 06 EPS-golden bootstrap procedure that runs from a separate `v1.0-pre-retire` worktree per Phase 9 D-04. Allowlisted in the gate to avoid false-positive scrubbing.
- **bin/cbinit re-created Qt-side.** Pre-Phase-9 the legacy `bin/cbinit` (Task 2 deletion) was the ONLY path producing pp/ppinit; cbl_tout_qt had no equivalent. Without a Qt-side cbinit, the renamed bin/cbl_tout would silently skip ppinit, breaking the harness chain at line 122. Rule 2 auto-add: new bin/cbinit modelled verbatim on bin/cbmail's Qt-only link recipe; wired into bin/cbl_tout dispatcher (1 added line before cbmail).
- **Qt-linked ppinit needs offscreen platform for CLI.** Empirically: without `QT_QPA_PLATFORM=offscreen`, the new pp/ppinit aborts on Qt's xcb plugin display-connect failure ("could not connect to display"). Wrapped the pp/ppinit invocation in bin/ab_sweep_phase8.sh:122 (Rule 1 — bug fix; Rule 2 — auto-add for missing critical functionality).
- **`bin/cbl_tout` exits 0 end-to-end.** Clean rebuild (`rm -rf xvue/qt/build pp/pp* && bin/cbl_tout`) produces all 11 pp/* binaries; the build log contains "verify_no_x11_in_build: scanning... OK: no X11 references in active build path" and "[31%] Built target verify_no_x11_in_build". Final line: "TOUS les MODULES EXECUTABLES ... sont crees — Qt variant".
- **4/5 testa baselines captured.** Phase 8 harness sweep (qt-1x, --smoke-only): pan2d (76 KB), nafems_le1 (320 KB), cavity2d (51 KB), heat1d (22 KB), nlsecu (0 B — Phase 8 override #5 carry-forward to Plan 9-06; matches Plan 9-02's identical 4/5 result).

## Task Commits

| Task | Hash      | Type | Subject                                                                                       |
|------|-----------|------|-----------------------------------------------------------------------------------------------|
| 1    | aae9473   | docs | Task 1 — pre-rename downstream-tooling audit (RETIRE-02)                                      |
| 2+3  | b4af705   | feat | Tasks 2-3 — delete legacy cb*+Makefile* + rename cb*_qt to cb* + collapse pp/pp*_qt body refs |
| 4    | c1ec327   | feat | Task 4 — collapse pp/pp${MODULE}_qt -> pp/pp${MODULE} in Phase-8/Qt harness scripts (BLOCKER #4) |
| 5    | a680f02   | feat | Task 5 — RETIRE-02 grep gate (bin/test_no_x11_in_build.sh + verify_no_x11_in_build CMake target) |
| 6    | c8dd804   | feat | Task 6 — Qt-only ppinit support + harness offscreen platform fix + final RETIRE-02 build/sweep validation |

A trailing `docs(09-03):` metadata commit will land for this SUMMARY itself.

## Files Created/Modified/Deleted

See `key-files` in the frontmatter. Summary counts:
- 4 created files (2 .md docs, 1 grep-gate shell, 1 cb script)
- 11 renamed files (cbl_tout_qt + 5 module + 5 xvtest scripts via `git mv`)
- 35 deleted files (11 legacy active-build cb* + 18 out-of-active-build cb* + 3 Makefile* + 3 Motif/HP-UX relics)
- 9 modified files (renamed-script bodies + dispatcher + harness + CMake + README)

## Empirical Counts vs. 09-01-AUDIT-BASELINE.md §3

| # | Reference                                                                  | Pre (Plan 9-01 baseline)  | Post (Plan 9-03)               | Expected                                  | Status |
|---|----------------------------------------------------------------------------|---------------------------|--------------------------------|-------------------------------------------|--------|
| 2 | `grep -rln 'lX11\|/usr/X11R6' bin/cb* bin/Makefile* \| wc -l`              | ~32                       | 0                              | 0 (RETIRE-02 success)                     | ✓      |
| 3 | `ls bin/Makefile* \| wc -l`                                                | 3                         | 0                              | 0 (Makefile* deleted)                     | ✓      |
| 4 | bin/cbl_tout_qt exists                                                     | yes                       | no                             | no (renamed to cbl_tout via git mv)       | ✓      |
| 5 | `grep -E 'cb<module>_qt' bin/cbl_tout \| wc -l`                            | 10 (cb*_qt invocations)   | 0                              | 0 (dispatcher updated)                    | ✓      |
| 6 | `grep -hE 'nompp=pp/pp[a-z0-9]+_qt' bin/cb<*> \| wc -l`                    | ~10                       | 0                              | 0 (BLOCKER #2 collapse)                   | ✓      |
| 8 | `grep -c 'pp/pp${MODULE}_qt' bin/ab_sweep_phase8.sh`                       | 3                         | 0                              | 0 (BLOCKER #3-#4 collapse)                | ✓      |
| 9 | `grep -c 'pp/pp${PREREQ_MODULE}_qt' bin/ab_sweep_phase8.sh`                | 1                         | 0                              | 0 (BLOCKER #4 collapse)                   | ✓      |
| — | `! ls pp/pp*_qt` (post-build)                                              | yes (pp/pp*_qt existed)   | no (no _qt-suffixed binaries)  | no (BLOCKER #2 build-side)                | ✓      |
| — | `bin/test_no_x11_in_build.sh` exit                                         | n/a (didn't exist)        | 0                              | 0 (gate green on new tree)                | ✓      |
| — | bin/cbl_tout exit                                                          | 0 (Plan 9-02)             | 0                              | 0 (CLAUDE.md invariant)                   | ✓      |
| — | testa cases capturing PNGs                                                 | 4/5 (Plan 9-02)           | 4/5 (pan2d, nafems_le1, cavity2d, heat1d) | 4/5 (nlsecu carries to 9-06) | ✓      |

The remaining baseline rows (#1, #7) belong to Plan 9-02 (already closed) and Plans 9-04/9-05 (ImageMagick + docs).

## Decisions Made

- **Tasks 2 + 3 folded into one commit** — the plan text mandates `git rm bin/cbl_tout && git mv bin/cbl_tout_qt bin/cbl_tout` as a single conceptual change. Plan Task 2 action explicitly says "DO NOT commit yet — Task 3 chains the rename in the same conceptual commit". Following the plan as written.
- **Pp _qt suffix policy locked at executor level** — pre-revision-iter1 the plan was ambiguous (cb script bodies referenced `nompp=pp/ppmail_qt` even though plan said cleanup); revision iter1 BLOCKERs #2-#4 closed the loop. This plan's executor enforced the lock at three layers: cb body sed, ab_sweep_phase8.sh sed, post-build `! ls pp/pp*_qt` assertion.
- **Phase-6 .cpp comments left as historical** — `xvue/qt/src/xvue_qt_*_actions.cpp` (4 files) + `xvue/qt/tests/test_xvue_qt_*_menu.cpp` (4 files) reference `pp/pp{elas,flui,nlse,ther}_qt` in code comments describing Phase-6 architectural intent (T-6.x-PARALLEL-BUILD historical mitigation notes). Plan Task 1 step 4 (b) explicitly classifies these as historical (NOT to be edited) — they describe Phase-6 state at its end, semantically equivalent to .planning/phases/06.*/SUMMARY.md content.
- **Three additional active harness scripts updated outside frontmatter `<files_modified>`** — bin/phase8_case_batch_map.sh, bin/qt-capture.sh, bin/xvtest0-pixmap-roundtrip.sh. Plan §RENAME-IMPACT classifies these as category (a) "Update in lockstep". Frontmatter `<files_modified>` listed only ab_sweep_phase8.sh; Rule 2 auto-add: the harness is incoherent if these still reference nonexistent pp/pp*_qt binaries. Plan acceptance criteria at the must_haves level remained satisfied (the four explicit ones in `<files_modified>` cover the strict policy lock; the three additional ones are coherency follow-through).
- **bin/cbinit re-created Qt-side as Rule 2 auto-add** — the legacy cbinit was on the deletion list (Task 2 group A), but no `cbinit_qt` ever existed. Without a Qt-side cbinit, bin/cbl_tout's dispatcher would silently skip ppinit, breaking the test harness at bin/ab_sweep_phase8.sh:122. Rule 2 (missing critical functionality): created bin/cbinit modelled verbatim on bin/cbmail. Documented as a deviation under "Auto-fixed Issues" below.
- **Sweep run with --smoke-only** — Plan 9-02's auto-fix Rule 3 carry-forward: after wholesale X11-baseline deletion there are no X11 PNGs to compare against; only smoke (capture + presence) verdicts are meaningful.
- **Build artifacts mail/lib etc. NOT staged in any plan-9-03 commit** — these are tracked but generated; same as Plans 9-01 / 9-02 explicit policy.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - bug] bin/cbmail_qt nompp= line drifted from baseline line 4 to line 31**
- **Found during:** Task 1 audit (concrete pattern confirmation)
- **Issue:** Plan text expected `grep -n 'nompp=pp/pp' bin/cbmail_qt` to return `4:nompp=pp/ppmail_qt`. Empirically returns `31:nompp=pp/ppmail_qt` — the file structure had grown (Phase 5+ Qt boilerplate before the link section) since the plan baseline was captured.
- **Fix:** Plan acceptance criteria check **content**, not line number (`grep -q '^nompp=pp/ppmail$' bin/cbmail`). Recorded the drift in 09-03-RENAME-IMPACT.md and proceeded — sed targets the literal pattern, not a line range, so the collapse worked correctly.
- **Files modified:** None additional (drift documented; Task 3 Step 3 sed ran clean).
- **Committed in:** aae9473 (Task 1 audit) + b4af705 (Task 3 sed).

**2. [Rule 2 - missing critical functionality] Qt-side bin/cbinit did not exist; needed for pp/ppinit production**
- **Found during:** Task 6 (post-rebuild missing pp/ppinit)
- **Issue:** Plan Task 3 step 2 has `git mv bin/cbinit_qt bin/cbinit  # if present` with the note "(skip silently — some _qt counterparts may not exist)". cbinit_qt indeed never existed (the Qt build chain was incomplete pre-Phase-9 — the legacy X11 cbinit was the only ppinit producer). After this plan deleted cbinit (Task 2), no Qt-side equivalent existed; bin/cbl_tout's dispatcher had no way to produce pp/ppinit; test harness at bin/ab_sweep_phase8.sh:122 (`echo $CASE | $MEFISTO/pp/ppinit`) would fail.
- **Fix:** Created bin/cbinit (60 lines, 0o755) modelled verbatim on bin/cbmail's Qt-only link recipe (gfortran with mail/lib util/lib xvue/lib + -Lxvue/qt/build -lxvueqt $QT_LIBS -lstdc++ -lgfortran). Wired into bin/cbl_tout dispatcher (added one line `$MEFISTO/bin/cbinit` before cbmail).
- **Files modified:** bin/cbinit (new), bin/cbl_tout (+1 dispatch line).
- **Verification:** `bin/cbl_tout` builds pp/ppinit (1.27 MB); `pp/ppinit` runs (with `QT_QPA_PLATFORM=offscreen` — see deviation #3 below).
- **Committed in:** c8dd804.

**3. [Rule 1 - bug fix] pp/ppinit (Qt-linked) needs QT_QPA_PLATFORM=offscreen for headless CLI use**
- **Found during:** Task 6 sweep run (initial sweep produced 0-byte PNGs because ppinit aborted)
- **Issue:** Pre-Phase-9 pp/ppinit was X11-linked (linked against xvuelc.o + libX11). Plan 9-02 deleted xvuelc.o; Plan 9-03's new Qt-linked pp/ppinit (deviation #2 above) needs Qt's platform plugin. Without `QT_QPA_PLATFORM=offscreen`, the Qt xcb plugin tries to connect to a display ("could not connect to display ... no Qt platform plugin could be initialized") and SIGABRTs, leaving no MS files.
- **Fix:** bin/ab_sweep_phase8.sh:122 ppinit invocation wrapped with `env QT_QPA_PLATFORM=offscreen`, matching the qt-1x/qt-2x/qt-omp branches at lines 146/156/166. Comment added: "Phase 9 RETIRE-02: pp/ppinit is now Qt-linked (legacy X11 backend retired); needs QT_QPA_PLATFORM=offscreen for headless CLI use."
- **Files modified:** bin/ab_sweep_phase8.sh (1 line replaced + 2 comment lines added).
- **Verification:** Re-run sweep produced 4/5 PNGs (76+320+51+22 KB; nlsecu at 0 B per Phase 8 override #5).
- **Committed in:** c8dd804.

**4. [Rule 2 - documentation coherency] Three active harness scripts referenced pp/pp*_qt or cb*_qt outside the plan's frontmatter `<files_modified>`**
- **Found during:** Task 1 audit
- **Issue:** Plan §RENAME-IMPACT category (a) "Update in lockstep" includes `bin/ab_sweep_phase8.sh` (the headline) but the audit also surfaced bin/phase8_case_batch_map.sh (2 doc lines), bin/qt-capture.sh (2 usage examples), bin/xvtest0-pixmap-roundtrip.sh (3 occurrences including a `cb*_qt` reference). Plan frontmatter `<files_modified>` listed only ab_sweep_phase8.sh.
- **Fix:** Folded into Task 4 scope. Each is a small targeted edit (5-line block in xvtest0-pixmap-roundtrip; 2-line comment in phase8_case_batch_map; 1 example in qt-capture). The harness would be silently incoherent without these (e.g., xvtest0-pixmap-roundtrip's pre-existence guard would test for `pp/ppxvtest0_qt` which post-Phase-9 doesn't exist, exit 2 with a misleading "run bin/cbxvtest0_qt first" — also a nonexistent script).
- **Files modified:** bin/phase8_case_batch_map.sh, bin/qt-capture.sh, bin/xvtest0-pixmap-roundtrip.sh.
- **Verification:** `grep -nE '_qt' ...` shows only attribution comments and example output filenames remain; no executable references to nonexistent pp/pp*_qt or cb*_qt.
- **Committed in:** c1ec327 (Task 4 commit).

**5. [Rule 1 - bug fix] test_no_x11_in_build.sh false-positive on xvue/qt/tests/golden/scene01_driver.f**
- **Found during:** Task 5 (initial gate run)
- **Issue:** First run of `bin/test_no_x11_in_build.sh` failed with "FAIL: ... xvue/qt/tests/golden/scene01_driver.f". Inspection showed the failure was a HEADER COMMENT line: `-L/usr/X11R6/lib -lX11 -lXt -o scene01_x11`. This line documents the cross-tag bootstrap procedure for the Phase-7 EPS golden — it runs in a separate `v1.0-pre-retire` worktree per Phase 9 D-04, NOT from main. Reference documentation, not an active build dependency.
- **Fix:** Added an allowlist line to test_no_x11_in_build.sh: `| grep -v 'xvue/qt/tests/golden/scene01_driver\.f$'`. Comment in the gate script explains the rationale (similar to how test_no_imagemagick_in_qt.sh allowlists Qt-API tokens like `QPageSize`, `convertToOther`).
- **Files modified:** bin/test_no_x11_in_build.sh (added 1 grep -v line + 8-line allowlist comment).
- **Verification:** Gate exit 0 on the current tree; negative test (inject `-lX11` into bin/cbmail) still fails as expected.
- **Committed in:** a680f02 (Task 5 commit; baked into the initial file write).

**6. [Rule 1 - bug fix] bin/cbinit comment mentioning -lX11 triggered the new gate**
- **Found during:** Task 6 (full clean rebuild ran the verify_no_x11_in_build CMake target which inspects bin/cb*)
- **Issue:** Initial bin/cbinit content contained the comment `# the only path producing pp/ppinit (linked against -lX11 and xvuelc.o);` — the literal `-lX11` substring tripped the gate.
- **Fix:** Reworded the comment to "Pre-Phase-9 this script's predecessor lived only on the legacy X11 build chain" — describes the same migration history without using the `-lX11` literal.
- **Files modified:** bin/cbinit (1 comment block reworded).
- **Verification:** `bin/test_no_x11_in_build.sh` exits 0; the full Qt build's CMake `verify_no_x11_in_build` target ran clean during `cbl_tout`.
- **Committed in:** c8dd804 (folded into Task 6 — the bin/cbinit creation iterated through this fix before the final clean rebuild).

**7. [Rule 3 - blocking issue] Sweep run in --smoke-only (Plan 9-02 carry-forward)**
- **Found during:** Task 6 (sweep invocation)
- **Issue:** Plan Task 6 step 4 invokes `bin/ab_sweep_phase8.sh --mode qt-1x --cases ... --out-dir /tmp/09-03-post-retire02` (no `--smoke-only` flag). Same root cause as Plan 9-02's auto-fix Rule 3: the harness requires X11 baselines for diff comparison, and post-deletion none exist (legacy X11 backend is gone since Plan 9-02 + finalised in this plan).
- **Fix:** Re-ran with `--smoke-only` (the documented Plan 9-02 / 9-01 carry-forward auto-fix).
- **Files modified:** None (test-time only).
- **Verification:** 4/5 cases produced verdict=SMOKE PNGs (nlsecu deferred per Phase 8 override #5).
- **Committed in:** N/A (test-time invocation only).

**Total deviations:** 7 auto-fixes (3 Rule 1 bugs, 3 Rule 2 missing-functionality / coherency, 1 Rule 3 blocking-issue carry-forward). All within plan acceptance scope; no scope creep.

**Impact on plan:** All 6 tasks completed as specified. Two of the deviations (#2 cbinit, #3 ppinit offscreen) closed real Phase-9-introduced regressions that the plan-checker's pre-execution audit hadn't caught — exactly the kind of "discovered during execution" findings GSD Rule 1-3 are designed for. Recorded in this SUMMARY for downstream plans (9-04 / 9-05) to inherit cleanly.

## Issues Encountered

- **Plan-text vs. empirical line-number drift in bin/cbmail_qt nompp= location.** Plan acceptance gate uses content (`grep -q '^nompp=pp/ppmail$'`) which is line-agnostic; drift acknowledged in 09-03-RENAME-IMPACT.md but did not require any change to plan execution.
- **bin/cbinit Qt-link blueprint not in the plan.** The plan's Task 3 step 2 explicitly says "if cbinit_qt does not exist, skip silently" — but it doesn't address the consequence that `pp/ppinit` then has no producer. Plan-checker missed this. Captured as Rule 2 auto-add (deviation #2). Future Phase 9 plans (9-04 etc.) inherit a clean `bin/cbinit + bin/cbl_tout invokes it` chain.
- **pp/ppinit runtime needs Qt platform plugin.** Same Phase-9-introduced regression class as deviation #3. The fact that pp/ppinit is now Qt-linked (vs. previously X11-linked) is a consequence of the cleanest-link policy (xvue/lib references xv* graphics primitives that are only in libxvueqt.a post-Phase-9). Wrapping with QT_QPA_PLATFORM=offscreen for the harness is the same pattern used for the main-module sweep branches.

## Threat Flags

None new. The threat model in the plan covered T-09-04, T-09-04-A through T-09-04-D. All four are mitigated:

- **T-09-04 (downstream tooling):** Task 1 audit identified 5 active references to cbl_tout_qt; 4 actively edited (CMakeLists.txt 2 comments, README.md 3 lines, dispatcher loop), 1 class deferred (Phase-6 .cpp comments — Plan §Update vs Historical category b). 175 references total when worktree dirs excluded (down from 2324 with worktrees included); 5 outside .planning/.
- **T-09-04-A (re-injection of -lX11 in future commits):** verify_no_x11_in_build CMake ALL target wired; injection-test confirms it fails when X11 reappears.
- **T-09-04-B (rename collision):** Task 3 step 1 explicitly `git rm bin/cbl_tout` BEFORE `git mv bin/cbl_tout_qt bin/cbl_tout`.
- **T-09-04-C (BLOCKER #2: pp binaries WITHOUT _qt suffix):** Empirically: `! ls pp/pp*_qt` post-build returns no matches; 11 binaries with bare names produced.
- **T-09-04-D (BLOCKER #4: harness must reference pp/pp${MODULE} no _qt):** Pre-sweep assertion + the sweep itself dispatched correctly to renamed pp/* binaries; 4/5 PNGs captured.

## Next Phase Readiness

- **Plan 09-04 (RETIRE-03 LVIDEO/ImageMagick) unblocked.** This plan touched neither LVIDEO/ImageMagick callers nor docs. The legacy `bin/convertepsgif` survives untouched; xvue/video1.f / videofin.f / videonm.f survive untouched. They're flagged for RETIRE-03.
- **Plan 09-05 (RETIRE-04 docs) unblocked.** README/LISEZMOI not touched in this plan (only xvue/qt/README.md was touched, and only its build-section attribution).
- **Plan 9-06 (ppnlse_qt offscreen deadlock) unchanged.** nlsecu still produces 0-byte PNG; same Phase-8 override #5 carry-forward.
- **Plan 9-07 (3 Phase-7 deferred goldens) unaffected.** scene01_driver.f's bootstrap procedure (cross-tag from v1.0-pre-retire) is allowlisted in test_no_x11_in_build.sh; the procedure can run unchanged.
- **Plan 9-08 (harness --out-dir + CMake freshness) unaffected** by this plan (harness body changes here are coherency-only; --out-dir relative-path bug is orthogonal).

### Known limitations introduced (intentional, by plan design)

- **Surviving cb* scripts that reference ccxvue / xvuelc.c in dead code paths**: `bin/cbl_all` (line 19-20), `bin/cblg_all` (line 19), `bin/cblg_tout` (line 48). None reference -lX11 or X11R6 directly (they pass the gate). They're not on the active dispatch chain (cbl_tout doesn't invoke them). Per Plan §SCOPE BOUNDARY logged as deferred — these are pre-existing broken references that Plan 9-02 introduced and Plan 9-03 explicitly does not touch. Future cleanup if these become active dispatch targets.
- **bin/cbinit_qt was never created in pre-Phase-9 history**: this plan creates the gap-closing Qt-only `bin/cbinit` directly (not as a `cbinit_qt` rename). Documented in deviation #2 above.
- **Phase-6 .cpp comments still reference pp/pp*_qt**: 8 source files in xvue/qt/src/ + xvue/qt/tests/. Plan Task 1 step 4 (b) historical classification preserved.

## Self-Check: PASSED

- File `.planning/phases/09-retire-x11-backend/09-03-RENAME-IMPACT.md`: FOUND
- File `.planning/phases/09-retire-x11-backend/09-03-SUMMARY.md`: FOUND (this file, post-write)
- File `bin/test_no_x11_in_build.sh`: FOUND, executable
- File `bin/cbinit`: FOUND, executable
- File `bin/cbl_tout`: FOUND, executable (renamed from cbl_tout_qt via git mv)
- File `bin/cbl_tout_qt`: ABSENT ✓
- Commit `aae9473` (Task 1): FOUND in git log
- Commit `b4af705` (Tasks 2+3): FOUND in git log
- Commit `c1ec327` (Task 4): FOUND in git log
- Commit `a680f02` (Task 5): FOUND in git log
- Commit `c8dd804` (Task 6): FOUND in git log
- `! grep -E 'cb<module>_qt' bin/cbl_tout`: 0 matches ✓
- `! grep -hE 'nompp=pp/pp[a-z0-9]+_qt' bin/cb<*>`: 0 matches ✓
- `grep -c 'pp/pp${MODULE}_qt' bin/ab_sweep_phase8.sh`: 0 ✓
- `grep -c 'pp/pp${PREREQ_MODULE}_qt' bin/ab_sweep_phase8.sh`: 0 ✓
- `grep -c 'pp/pp${MODULE}' bin/ab_sweep_phase8.sh`: 3 ✓ (qt-1x/qt-2x/qt-omp)
- `grep -c 'pp/pp${PREREQ_MODULE}' bin/ab_sweep_phase8.sh`: 1 ✓ (PREREQ_MODULE branch)
- `bin/test_no_x11_in_build.sh` exit: 0 ✓
- `bin/test_no_imagemagick_in_qt.sh` exit: 0 ✓
- `bin/cbl_tout` exit: 0 ✓ (full clean rebuild)
- `verify_no_x11_in_build` ran in CMake build: ✓ (line in /tmp/09-03-cbl_tout.log)
- 11 pp/* binaries (no _qt suffix): ppinit ✓, ppmail ✓, ppelas ✓, ppflui ✓, ppther ✓, ppnlse ✓, ppxvtest{0..4} ✓
- `! ls pp/pp*_qt`: empty ✓
- 4 testa PNGs captured: pan2d (76 KB), nafems_le1 (320 KB), cavity2d (51 KB), heat1d (22 KB) ✓
- nlsecu deferred per Phase 8 override #5: ✓ (0-byte; carries to Plan 9-06)

---
*Phase: 09-retire-x11-backend*
*Completed: 2026-05-05*
