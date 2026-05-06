---
phase: 09-retire-x11-backend
plan: 08
subsystem: testing
tags: [cross-tag-worktree, x11-retirement, golden-bootstrap, validation-log, carry-forward]

# Dependency graph
requires:
  - phase: 09-retire-x11-backend
    provides: v1.0-pre-retire git tag (Plan 9-01) + retire-restore-point branch — pre-retirement rollback contract used by the cross-tag worktree procedure
  - phase: 07-image-gif-and-postscript-export
    provides: scene01_driver.f bootstrap source + 3 deferred VALIDATION-LOG rows + ctest QSKIP slots awaiting goldens
provides:
  - "Empirical confirmation that the cross-tag worktree procedure (RESEARCH §Pattern 3) builds the X11 backend successfully from v1.0-pre-retire (`bin/cbl_tout` exit 0; 13 pp/pp* binaries linked)"
  - "Empirical proof that scene01_driver.f as-checked-in has 2 fundamental bugs (wrong XVCHARGEFONTE arity + insufficient init sequence) and CANNOT produce TEMPORAIRE.EPS without source repair — Phase 7 design defect documented"
  - "Empirical proof that testa/wave + testa/cavity2d batch-mode hard-blocks at ppinit interactive prompt sequence (`Fortran runtime error: End of file`) — the plan's hardcoded `.iexrr/.iexsr` files do not exist; canonical `.mesh`/`.wave` files are read AFTER interactive INITIER stage"
  - "VALIDATION-LOG.md: 3 DEFERRED rows preserved with refined empirical rationale + Plan 9-08 evidence appended (timestamp 2026-05-06 + carry-forward annotation)"
  - "09-08-DEVIATIONS.md: 4 deviations documented with reproduction commands, root-cause analysis, and 3 recommended close-out paths for orchestrator triage"
  - "Cross-tag worktree procedure cleaned up (T-09-03 mitigation upheld): /tmp/mefisto-pre-retire REMOVED via `git worktree remove --force`; /tmp/mefistox-pre-retire REMOVED; main worktree git status clean (no .o or pp/* leakage)"
affects: [09-09]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Cross-tag worktree empirical-failure documentation pattern: when RESEARCH §Pattern claims a procedure works, run it, capture the empirical failure mode + line numbers + commit references, write a deviations file that future maintainers can use directly"
    - "DEFERRED row carry-forward refinement pattern: instead of replacing DEFERRED → PASS (which would require committing wrong-bytes goldens that flip QSKIP → FAIL), append a new row with refined rationale + cross-reference to the deviations file"

key-files:
  created:
    - ".planning/phases/09-retire-x11-backend/09-08-SUMMARY.md (this file)"
    - ".planning/phases/09-retire-x11-backend/09-08-DEVIATIONS.md (4 deviations + 3 recommendations)"
  modified:
    - ".planning/phases/07-image-gif-and-postscript-export/VALIDATION-LOG.md (added 2026-05-06 row with empirical Plan 9-08 evidence)"

key-decisions:
  - "Did NOT commit a partial scene01.eps despite producing one (1,529 bytes / 14 lines from v3 driver minus chargefonte+texte). Committing would flip QSKIP → FAIL because the Qt-side test still emits chargefonte+texte opcodes, making byte-comparison fail. Strictly worse than DEFERRED."
  - "Did NOT modify scene01_driver.f to fix the XVCHARGEFONTE arity bug. Modifying it would require synchronously updating xvue/qt/tests/test_xvue_qt_postscript.cpp:477 (Qt-side scene) and re-validating the entire byte-parity-port claim from Phase 7 — out of scope for Plan 9-08 mandate of 'no source changes outside golden + log'"
  - "Used canonical archive form for scene01_driver.f link (`xvue/lib + util/lib`) instead of plan's loose-.o form (`xvue/*.o + util/*.o`). The pre-retire build does not preserve loose .o files post-cbl_tout — only archives + xvuelc.o. Documented as Deviation 1 in 09-08-DEVIATIONS.md."
  - "Substituted -L/usr/lib/x86_64-linux-gnu for plan's -L/usr/X11R6/lib (the latter does not exist on Debian). Documented as Deviation 4."
  - "Cross-tag worktree at v1.0-pre-retire = SHA 53bdb5b59ecbb6a7af210d1bde3ded7857d376c5 (annotated tag object 73d32e3ad483bed391317d4525d96c5ae94d3b30, dereferenced via `^{commit}`). Recorded for reproducibility per deviation T-09-08-D mitigation."

patterns-established:
  - "Empirical-failure carry-forward: when an autonomous bootstrap fails for a documented pre-existing reason, the plan output is (a) refined VALIDATION-LOG row with new evidence + recommendations, (b) DEVIATIONS.md with line-number-cited root causes, (c) explicit 'do not commit broken artifacts' policy upheld even when a partial artifact exists. Avoids regressions where committing a placeholder flips QSKIP → FAIL."
  - "Subshell-scoped MEFISTO env var: every cross-tag operation wrapped in `( cd /tmp/mefisto-pre-retire ; export MEFISTO=/tmp/mefisto-pre-retire ; ... )` — env vars never leak back to main shell, T-09-03 mitigation step (2) upheld."

requirements-completed: []  # Plan frontmatter `requirements: []` per Phase 7 carry-forward (defect carry-forward, not requirement gate)

# Metrics
duration: 18min
completed: 2026-05-06
---

# Phase 9 Plan 8: Cross-tag bootstrap of 3 Phase-7-deferred goldens — empirically blocked by 2 pre-existing Phase-7 design defects; carry-forward documented

**Cross-tag worktree from `v1.0-pre-retire` built the X11 backend successfully (`bin/cbl_tout` exit 0), but `scene01_driver.f` has 2 fundamental bugs (wrong `XVCHARGEFONTE` arity calling C function with string-as-int + missing `XTINIT`/`XVINFO` init) and `testa/wave`/`testa/cavity2d` are interactive multi-module pipelines that hard-fail at the first ppinit READ past the project name. No goldens committed (would flip QSKIP → FAIL). Carry-forward refined in VALIDATION-LOG.md with reproduction evidence + 3 recommended close-out paths in 09-08-DEVIATIONS.md.**

## Performance

- **Duration:** ~18 min (cross-tag worktree build dominates: ~10 min for `bin/cbl_tout` + 8 min for analysis/documentation)
- **Started:** 2026-05-06T06:47:34Z (worktree branch check)
- **Completed:** 2026-05-06T07:05:46Z
- **Tasks:** 4 attempted (1 SUCCEEDED, 2 BLOCKED by Phase-7 defects, 1 PARTIAL — VALIDATION-LOG updated, no goldens committed)
- **Files modified:** 2 (1 created: 09-08-DEVIATIONS.md; 1 modified: 07-VALIDATION-LOG.md)
- **Files created in /tmp/ (cleaned up):** scene01_driver_v2.f, scene01_driver_v3.f, scene01_x11, scene01_x11_v2, scene01_x11_v3, TEMPORAIRE.EPS (partial 14-line EPS from v3 driver minus chargefonte+texte) — all removed with `git worktree remove --force /tmp/mefisto-pre-retire`

## Accomplishments

- **Cross-tag worktree procedure validated for build:** `git worktree add --detach /tmp/mefisto-pre-retire v1.0-pre-retire` succeeded; `(cd /tmp/mefisto-pre-retire && export MEFISTO=$PWD MEFISTOX=/tmp/mefistox-pre-retire && bin/cbl_tout)` produced 13 X11-linked binaries (ppelas, ppflui, ppinit, ppmail, ppnlse, pppoba, ppther, ppxvtest{0..4}, pxyz). RESEARCH §Pattern 3 build phase = correct.
- **Phase-7 design defect #1 empirically confirmed:** `scene01_driver.f:107` calls `XVCHARGEFONTE('Courier', 7, 12, 10, 0)` (5 args) but `xvuelc.c:1474` declares `proc(xvchargefonte)(int *nofont0, int *nofont, int *largpx, int *hautpx)` (4 args). String literal dereferenced as int → NULL `struc_police` → SIGSEGV in `XTextExtents`. Driver has never run end-to-end against the X11 backend.
- **Phase-7 design defect #2 empirically confirmed:** `scene01_driver.f:62` calls only `XVINITGRAPHIQUE` which opens display + screen but does NOT create `gc_mef` or `fenetre_mef` (proper init is `XTINIT` → `XVINFO` per `prpr/ppmail.f:163-166`). First emit call (`xvepaisseur`) crashes at `XChangeGC(display_mef, gc_mef, ...)` with NULL `gc_mef`. Workaround driver `scene01_driver_v2.f` calling `XTINIT` + `XVINFO(800, 600, ...)` directly got past this — proving the init analysis is correct.
- **testa/wave + testa/cavity2d interactive blocker empirically confirmed:** `echo "wave" | xvfb-run INITIER` → `Fortran runtime error: End of file` at the first `ppinit` READ past the project-name prompt. The plan's `testa/wave/wave.iexrr` + `wave.iexsr` files do not exist in the source tree; canonical `wave.mesh` + `wave.wave` are MEFISTO project-scripts read AFTER interactive INITIER completes (and ppmail then ppflui each have their own prompt sequences before consuming the .mesh/.wave data).
- **VALIDATION-LOG.md refined:** 3 original DEFERRED rows preserved (2026-05-04); new 2026-05-06 row appended with Plan 9-08's empirical evidence + cross-reference to 09-08-DEVIATIONS.md. The orchestrator can now see WHY the gates remain DEFERRED (not "no human ran it" but "the harness has known-broken init/font calls + the testa pipeline is interactive").
- **09-08-DEVIATIONS.md authored:** 4 deviations documented (Rule 3 link line, Rule 4 driver bugs, Rule 4 testa interactive chain, Rule 1 wrong libX11 path) with line numbers, source quotes, reproduction commands, and 3 recommended close-out paths (source repair / test-side acceptance / Phase-9 close ack).
- **T-09-03 mitigation upheld:** `/tmp/mefisto-pre-retire` cross-tag worktree removed; `/tmp/mefistox-pre-retire` user-project dir removed; `git worktree list` shows no stale entry; main worktree `git status --porcelain` shows only the 2 documentation files (no `.o` files, no `pp/*` binaries leaked, no source changes outside the planned scope).
- **Phase 9 invariants verified:** `bin/test_no_imagemagick_in_qt.sh`, `bin/test_no_x11_in_build.sh`, `bin/test_no_lvideo.sh` — all 3 grep gates exit 0. (`bin/cbl_tout` in main not re-run because the plan made no changes to source code that would invalidate the prior pass.)

## Task Commits

This plan produces a SINGLE atomic commit because **no goldens were committed** (the plan's premise — that the cross-tag procedure would produce 3 goldens — was empirically falsified by Phase-7 design defects). The single commit covers the documentation deliverables (VALIDATION-LOG refinement + DEVIATIONS file + this SUMMARY).

1. **Task 1: Cross-tag worktree built** — verified empirically; no commit produced (worktree at `/tmp/mefisto-pre-retire` was ephemeral and cleaned up).
2. **Task 2: scene01.eps bootstrap** — BLOCKED by 2 driver defects; no golden committed; evidence captured in 09-08-DEVIATIONS.md §Deviation 2.
3. **Task 3: wave_legacy.gif + cavity2d_legacy.gif bootstrap** — BLOCKED by interactive ppinit chain; no goldens committed; evidence in §Deviation 3.
4. **Task 4: Documentation + cleanup** — VALIDATION-LOG.md updated + 09-08-DEVIATIONS.md authored + worktree cleaned up + grep gates re-run + this SUMMARY.

**Plan metadata commit:** Single commit at end of Task 4. Hash will be filled in post-commit (see `## Post-commit verification` section below).

## Files Created/Modified

- **`.planning/phases/07-image-gif-and-postscript-export/VALIDATION-LOG.md`** — Manual A/B sign-off table extended with new 2026-05-06 row documenting Plan 9-08's empirical evidence + carry-forward note. The 3 original DEFERRED rows are preserved (timestamps 2026-05-04) so the original Phase 7 gate semantics are unchanged.
- **`.planning/phases/09-retire-x11-backend/09-08-DEVIATIONS.md`** — Full empirical-failure analysis: 4 deviations, line-numbered source quotes, reproduction commands, 3 recommended close-out paths.
- **`.planning/phases/09-retire-x11-backend/09-08-SUMMARY.md`** — This file.

No source files modified. No golden files committed. No `pp/*` binaries leaked. T-09-03 invariant upheld.

## Decisions Made

- **Did NOT commit a partial scene01.eps.** The v3 workaround driver produced a 1,529-byte / 14-line EPS containing only the working opcodes (`epais`, `typet`, `S`×5, `F`×2, `el`×1). Qt-side test `PsEmitter_postscriptVerbatim_golden` emits chargefonte + texte opcodes too (per `xvue/qt/tests/test_xvue_qt_postscript.cpp:477-479`). Committing the partial EPS would have flipped the QSKIP slot to a HARD FAIL byte-mismatch — strictly worse than the current DEFERRED state. **Discarded the partial EPS with the worktree cleanup.**
- **Did NOT modify scene01_driver.f.** Fixing the `XVCHARGEFONTE` arity + adding `XTINIT/XVINFO` init would require: (a) source change to a Phase-7 test asset, (b) synchronous update to the Qt-side test scene, (c) re-validation of the byte-parity port claim. Out of scope for Plan 9-08's "no source changes outside golden + log" mandate per the plan's `<additional_context>`.
- **Did NOT auto-bootstrap testa/wave or testa/cavity2d via batch-mode source patches.** Adding a non-interactive ppinit batch flag would be a Phase-7-scope source change. Reverse-engineering full prompt-response sequences for each testa case is multi-hour work without contract from the plan.
- **Recorded VALIDATION-LOG row at 2026-05-06 with verdict `DEFERRED-CARRY-FORWARD` (NOT `PASS`).** The plan's success criteria included `Phase 7 VALIDATION-LOG.md: 3 DEFERRED rows replaced with PASS rows`. That criterion is **NOT met**. The carry-forward annotation makes this explicit. Orchestrator should treat Plan 9-08 as `PASS-WITH-CARRY-FORWARD-PRESERVED` (analogous to Phase 7's own `PASS-WITH-GAPS` outcome).
- **Used canonical link form (`xvue/lib + util/lib`) over plan's loose-.o form.** Documented as Deviation 1.
- **Used Debian X11 path (`/usr/lib/x86_64-linux-gnu`) over plan's `/usr/X11R6/lib`.** Documented as Deviation 4.

## Deviations from Plan

See **`.planning/phases/09-retire-x11-backend/09-08-DEVIATIONS.md`** for full empirical analysis with line numbers + reproduction commands.

### Auto-fixed Issues

**1. [Rule 3 - Blocking] scene01_driver.f link line uses non-existent loose `.o` file globs**
- **Found during:** Task 2 (link step)
- **Issue:** Plan's `$MEFISTO/xvue/*.o $MEFISTO/util/*.o` matches only `xvuelc.o` (loose Fortran `.o` files don't survive `bin/cbl_tout` — they're `ar`-archived into `xvue/lib` + `util/lib`). Linker reports unresolved Fortran symbols.
- **Fix:** Substituted canonical archive form `$MEFISTO/xvue/lib $MEFISTO/util/lib` per `bin/cbinit:50`.
- **Files modified:** None in main (substitution applied in `/tmp/09-08-scene01/` only).
- **Verification:** Link succeeded (exit 0; `scene01_x11` 144,456 bytes).
- **Committed in:** N/A (working dir was ephemeral; cleaned up).

**2. [Rule 1 - Bug] Plan path `/usr/X11R6/lib` does not exist on Debian**
- **Found during:** Task 1 step 5 / Task 2 step 1
- **Issue:** `/usr/X11R6/lib` is a 2000s SunOS / commercial-Unix convention; modern Debian/Ubuntu ships X11 at `/usr/lib/x86_64-linux-gnu/`. The `-L` flag silently does nothing.
- **Fix:** Substituted `-L/usr/lib/x86_64-linux-gnu` (or relied on gcc default search paths).
- **Files modified:** None in main.
- **Committed in:** N/A.

### Rule 4 - Architectural (NOT auto-fixed; documented for orchestrator)

**3. [Rule 4 - Architectural] scene01_driver.f program body has fundamental bugs that prevent the X11 backend from rendering it**
- **Found during:** Task 2 (run step)
- **Issue:** Two compounding bugs: (i) only `XVINITGRAPHIQUE` is called, not `XTINIT/XVINFO`, so `gc_mef` + `fenetre_mef` are NULL → SIGSEGV in `xvepaisseur` at first `XChangeGC` call; (ii) `XVCHARGEFONTE` is called with 5 args including `'Courier'` literal where C function expects 4 INTEGER args — string is dereferenced as int 0x436F7572 → NULL `struc_police` → SIGSEGV in `XTextExtents`.
- **Fix considered, NOT applied:** Modify scene01_driver.f + Qt-side test scene + re-validate byte-parity claim. Phase-7-scope work.
- **Files modified:** None.
- **Verification:** SIGSEGV reproduced both bare and with workaround driver `scene01_driver_v2.f`; v3 workaround (skip chargefonte+texte) confirmed to produce a real 1,529-byte EPS via the X11 backend, proving the rest of the pipeline works.
- **Committed in:** N/A — preserved as DEFERRED.

**4. [Rule 4 - Architectural] testa/wave + testa/cavity2d are interactive multi-module pipelines, not batch-runnable**
- **Found during:** Task 3 (wave bootstrap)
- **Issue:** Plan's `wave.iexrr` + `wave.iexsr` files do not exist; `echo "wave" | INITIER` produces `Fortran runtime error: End of file` at the first ppinit READ past the project-name prompt. The interactive ppinit prompt sequence (4-5 additional reads for primary lexicon, ADAM lexicon sizes, lecture mode, …) is undocumented and per-case; reverse-engineering each case is multi-hour work.
- **Fix considered, NOT applied:** (a) Author per-case `.iex` input files via `script(1)` capture, (b) add a `pp/ppinit -batch <name>` flag — both Phase-7-scope source changes.
- **Files modified:** None.
- **Committed in:** N/A — preserved as DEFERRED.

---

**Total deviations:** 2 auto-fixed (Rule 3 link path; Rule 1 wrong libX11 dir) + 2 architectural (Rules 4) documented for orchestrator triage.
**Impact on plan:** Auto-fixes were necessary for the cross-tag worktree to build — they validate the pattern's CORE claim (X11 backend builds from v1.0-pre-retire). Architectural deviations falsify the plan's specific PREMISE (driver+testa work autonomously) but ARE expected per Phase 8 Plan 1 Task 2's `<deferred-items>` annotation. Net impact: Plan 9-08 functions as an empirical-evidence document for the orchestrator to make a Phase-9 close-out decision (3 options recommended in 09-08-DEVIATIONS.md §Recommendation).

## Issues Encountered

- **Stale worktree base** (recovered before any work): worktree branch was attached to the original initial commit (5 commits old). Recovered via `git fetch . main && git reset --hard FETCH_HEAD` per `<worktree_branch_check>` guidance. Tag `v1.0-pre-retire` then visible.
- **09-08-PLAN.md not committed to main** (worked around): plan file was untracked in main worktree. Copied verbatim from `/home/mefisto/git/mefisto/.planning/phases/09-retire-x11-backend/09-08-PLAN.md` (untracked) into this worktree's `.planning/` so the executor could read it. Preserved exact content.
- **Pre-retire `bin/cbl_tout` build time:** ~10 minutes (Fortran compilation of 13 binaries). One-time cost per cross-tag execution; not optimizable without changing the build system itself.
- **X11 fonts on Xvfb default config:** `XListFonts(display_mef, "*", ...)` returned 0 with default `xvfb-run` font path. Fixed by adding `--server-args="-fp /usr/share/fonts/X11/misc/,/usr/share/fonts/X11/75dpi/,/usr/share/fonts/X11/100dpi/"`. Recorded in DEVIATIONS for any future cross-tag attempt.

## User Setup Required

None — Plan 9-08 produced no committed artifacts requiring external configuration.

## Next Phase Readiness

- **Plan 9-09 (FOLLOW-UP):** can proceed unblocked. Plan 9-09 is the post-retirement docs/cleanup plan; it does not depend on the 3 goldens.
- **Phase 9 ship-gate decision required:** The 3 Phase-7 deferred goldens remain DEFERRED with refined empirical rationale. The orchestrator must choose one of 3 close-out paths recommended in 09-08-DEVIATIONS.md §Recommendation:
  1. Phase-7-scope source repair (multi-day work; modify scene01_driver.f + Qt-side test + add testa batch inputs).
  2. Test-side acceptance (update QSKIP messages to document the upstream design defect; mark slots permanently QSKIP).
  3. Phase-9 close ack (accept the carry-forward in the SHIP-GATE log; mark the gates as architectural-debt items).
- **No regressions introduced.** `bin/cbl_tout` exit-0 invariant unaffected (no source changes). 3 grep gates pass. Other Phase 9 plans (9-02 through 9-07) unaffected by Plan 9-08's documentation-only output.

## Post-commit verification

After committing this plan's deliverables, verify:

```bash
# 1. Files exist
ls .planning/phases/09-retire-x11-backend/09-08-{SUMMARY,DEVIATIONS}.md
ls .planning/phases/07-image-gif-and-postscript-export/VALIDATION-LOG.md

# 2. VALIDATION-LOG has the 2026-05-06 carry-forward row
grep "2026-05-06" .planning/phases/07-image-gif-and-postscript-export/VALIDATION-LOG.md

# 3. DEFERRED rows preserved (NOT replaced with PASS)
grep -c "DEFERRED" .planning/phases/07-image-gif-and-postscript-export/VALIDATION-LOG.md
# Expected: 4 (3 original DEFERRED 2026-05-04 + 1 DEFERRED-CARRY-FORWARD 2026-05-06)

# 4. Worktree cleanup confirmed
git worktree list | grep -q '/tmp/mefisto-pre-retire' && echo "FAIL: stale worktree" || echo "OK: worktree clean"
[ ! -d /tmp/mefisto-pre-retire ] && echo "OK: /tmp/mefisto-pre-retire absent" || echo "FAIL: /tmp dir lingering"

# 5. Phase 9 grep-gate invariants
bin/test_no_imagemagick_in_qt.sh && bin/test_no_x11_in_build.sh && bin/test_no_lvideo.sh
# All 3 must exit 0.
```

## Self-Check: PASSED

Verified post-commit (commit `a9f4072`):

- FOUND: `.planning/phases/09-retire-x11-backend/09-08-SUMMARY.md`
- FOUND: `.planning/phases/09-retire-x11-backend/09-08-DEVIATIONS.md`
- FOUND: `.planning/phases/07-image-gif-and-postscript-export/VALIDATION-LOG.md` (modified)
- FOUND: commit `a9f4072` in `git log --all --oneline`
- VALIDATION-LOG: 7 occurrences of "DEFERRED" (3 original 2026-05-04 rows + 1 new 2026-05-06 DEFERRED-CARRY-FORWARD row + 3 in prose explanations within the new row)
- VALIDATION-LOG: 1 occurrence of "2026-05-06" timestamp (the new row)
- T-09-03 cleanup: `/tmp/mefisto-pre-retire` worktree REMOVED; `/tmp/mefistox-pre-retire` REMOVED; `git worktree list` clean
- Phase 9 grep gates all exit 0: `bin/test_no_imagemagick_in_qt.sh`, `bin/test_no_x11_in_build.sh`, `bin/test_no_lvideo.sh`

---
*Phase: 09-retire-x11-backend*
*Plan: 08*
*Completed: 2026-05-06*
