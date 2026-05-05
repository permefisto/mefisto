# Phase 9: Retire X11 backend - Context

**Gathered:** 2026-05-05
**Status:** Ready for planning (process-gated — see Phase Boundary)

<domain>
## Phase Boundary

Delete every X11 and ImageMagick code path, linker line, script, and
documentation reference, leaving Qt 6 as the single graphics backend.
Plus fold in 4 Phase-8 carry-forward items (defects + Phase-7
goldens + harness bug) so the retirement leaves a single clean post-X11
codebase.

In scope:
- RETIRE-01: delete `xvue/xvuelc.c` + `bin/ccxvue`
- RETIRE-02: remove `libX11` linker lines + `/usr/X11R6/lib64` paths from every `bin/cb*` script
- RETIRE-03: delete `bin/convertepsgif` + `xvue/video1.f` + `xvue/videofin.f` + `xvue/videonm.f` + 18+ Fortran tracer subroutines (LVIDEO pipeline) + audit `grep -rn 'convert' bin/ td/ testa/ testf/`
- RETIRE-04: update README, LISEZMOI, install scripts to list only Qt 6 runtime deps (`qt6-base`, `libqt6imageformats6-plugins`); remove `libX11-dev` + ImageMagick
- Phase-8 carry-forward #1 (matched-dim Qt recapture): wire `MEFISTO_QT_WINDOW_SIZE` (or equivalent env) so Qt backing pixmap can be sized to match X11 baseline dims; eliminates the resample-confound that dominated Phase 8 AE diffs
- Phase-8 carry-forward #2 (ppnlse_qt deadlock fix): ppnlse_qt offscreen + MEFISTO_BATCH_X11=1 deadlock at startup (no NLSER banner reached)
- Phase-8 carry-forward #3 (3 Phase-7 deferred goldens): scene01.eps + wave_legacy.gif + cavity2d_legacy.gif — bootstrap procedure under Qt-only conditions (X11 backend is gone)
- Phase-8 carry-forward #4 (harness `--out-dir` relative-path bug): canonicalize via `realpath` early in `bin/ab_sweep_phase8.sh`
- Phase-8 D-09 (CMake guard): add `verify_pp_qt_freshness` ALL target

Out of scope:
- New Qt features (Qt 6 surface stays as-is from Phase 7)
- Qt 7 migration (future milestone)
- testa/ test case retirement (cases stay; only X11 paths to invoke them go away)
- Pre-Phase-9 release tagging (handled by Rollback safety contract before phase start)

</domain>

<decisions>
## Implementation Decisions

### Deletion approach

- **D-01:** Incremental per-RETIRE-NN. 4 mechanical retirement plans
  (one per RETIRE-NN requirement) plus separate carry-forward plans.
  Easier bisect, smaller blast radius, recoverable per-plan.
- **D-02:** Plan order matches RETIRE-NN order: RETIRE-01 (xvuelc.c +
  ccxvue) → RETIRE-02 (libX11 linker lines) → RETIRE-03 (ImageMagick +
  LVIDEO) → RETIRE-04 (docs). Each plan independently buildable +
  testable; later plans depend on earlier ones (no LVIDEO removal
  before xvuelc.c removal).

### Carry-forward scope

- **D-03:** All 4 Phase-8 carry-forwards fold into Phase 9. Plans
  9-NN_A through 9-NN_D ship after the 4 RETIRE-NN plans:
    9-05: matched-dim Qt recapture (overrides #1 + #2 closure for
          the 14 CHECK cells that Phase 8 had to override)
    9-06: ppnlse_qt offscreen+BATCH_X11 deadlock fix (override #5)
    9-07: 3 Phase-7 deferred goldens bootstrap under Qt-only conditions
          — scene01.eps + wave_legacy.gif + cavity2d_legacy.gif. Phase 7
          ctest QSKIPs flip to PASS. (X11 backend is gone, so legacy-side
          baseline production must use a pre-Phase-9 git checkout if any.)
    9-08: harness `--out-dir` relative-path bug fix in
          `bin/ab_sweep_phase8.sh` + CMake `verify_pp_qt_freshness` ALL
          target (Phase-8 D-09)
- **D-04:** Carry-forward plans are NOT prerequisites for the 4
  RETIRE-NN plans. RETIRE-NN can run first; carry-forwards run last.
  Phase 9 closes only when ALL 8 plans (4 retire + 4 carry-forward)
  are signed off.

### LVIDEO retirement scope

- **D-05:** RETIRE-03 retires the entire LVIDEO pipeline alongside
  `bin/convertepsgif`. Files deleted: `xvue/video1.f`, `xvue/videofin.f`,
  `xvue/videonm.f` + the 18+ Fortran tracer subroutines that call into
  it (`flui/trvi2d.f`, `ther/trplse.f`, etc.). Single coherent
  retirement. Phase 7 README §9 explicitly identified LVIDEO as Phase
  9 RETIRE-03 deferral target.
- **D-06:** Tracer subroutines that ONLY exist for LVIDEO get deleted
  outright. Tracer subroutines that ALSO serve non-LVIDEO purposes
  (e.g., interactive draw routines) lose only their LVIDEO entry
  points; surrounding callers + Fortran COMMON blocks adjust.
- **D-07:** Animation export under Qt is provided by the Phase 7
  XvueExport `saveGifTo` ffmpeg path (already shipped). LVIDEO removal
  does NOT remove animation capability — it removes the legacy parallel
  pipeline. Users with Phase-7 Qt builds get GIFs via File→Export→GIF.

### Rollback safety

- **D-08:** Pre-Phase-9 checkpoint: maintainer runs `git tag v1.0-pre-retire`
  AND creates branch `retire-restore-point` from current main. Both
  pushed to origin. Recoverable via single `git reset --hard
  v1.0-pre-retire` command.
- **D-09:** v1.0-pre-retire tag MUST be created BEFORE Phase 9 Plan 1
  starts. The first Phase 9 plan blocks until `git tag --list
  v1.0-pre-retire` returns non-empty.
- **D-10:** During Phase 9 execution, each plan's commit chain stays
  bisectable: `bin/cbl_tout` MUST exit 0 after every commit (failed
  commits are amended/squashed before the next plan starts).

### Process gate (CRITICAL)

- **D-11:** Phase 9 starts ONLY when maintainer explicitly confirms the
  one-release-cycle A/B window has closed. The window opened
  2026-05-05 with Phase 8 ship-gate sign-off. Closure is a
  maintainer-judgment call — NOT a fixed-date deadline. The user
  signals closure via:
    - committing `## Phase 9 entry: A/B window closed YYYY-MM-DD —
      maintainer initials` to STATE.md decisions block, OR
    - explicit reply to `/gsd-progress` indicating Phase 9 is unblocked
- **D-12:** Pre-closure work allowed: `/gsd-discuss-phase 9` (this
  artifact), `/gsd-plan-phase 9 --research-phase 9` (research only),
  `/gsd-plan-phase 9` (plan creation). Disallowed pre-closure:
  `/gsd-execute-phase 9` (the actual retirement). The plan-checker
  must confirm the gate before any Wave 1 dispatch.

### Claude's Discretion

- Per-plan task breakdown (planner picks the exact split per RETIRE-NN +
  per carry-forward).
- Research scope: empirical audit of every reference to `xvuelc`, `libX11`,
  `convert`, `convertepsgif`, `video1.f` / `videofin.f` / `videonm.f`,
  ImageMagick across the tree, before planning.
- Test post-deletion: re-run all 5 BUILD-10 testa cases through the Qt-only
  pipeline; commit fresh post-Phase-9 baselines.
- Carry-forward 9-07 (Phase-7 goldens): X11 baseline production — use
  v1.0-pre-retire git checkout in a separate worktree, run
  bin/convertepsgif there, copy artifacts back, commit on main. This is
  the single cross-tag operation in Phase 9.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Roadmap and requirements

- `.planning/ROADMAP.md` §"Phase 9: Retire X11 backend" — phase goal + success criteria + RETIRE-01..04 requirement IDs
- `.planning/REQUIREMENTS.md` §RETIRE-01..04 — the 4 retirement requirements
- `.planning/PROJECT.md` — Key Decisions table; D-7 mandates parallel X11 build kept alive for one release cycle for A/B (now closing)

### Phase 7 hand-off (mandatory read)

- `.planning/phases/07-image-gif-and-postscript-export/VERIFICATION.md` §9 — 3 deferred goldens (scene01.eps, wave_legacy.gif, cavity2d_legacy.gif) bootstrap procedure; D-03 of this phase covers Phase 9 plan 9-07
- `.planning/phases/07-image-gif-and-postscript-export/07-06-SUMMARY.md` §"Plan 06 Task 3" — golden bootstrap procedure detail
- `xvue/qt/README.md` §Phase-7 — explicit LVIDEO Phase 9 RETIRE-03 deferral note (D-05 of this phase)

### Phase 8 hand-off (mandatory read)

- `.planning/phases/08-ab-validation-on-testa-subset/08-CHECKLIST.md` — 5 override decisions; overrides #1 + #2 + #5 explicitly carry into Phase 9 plans 9-05 + 9-06
- `.planning/phases/08-ab-validation-on-testa-subset/08-VALIDATION.md` — VALID-01..07 coverage; VALID-06 redefined to pan2d-only (override #3 waiver applies)
- `.planning/phases/08-ab-validation-on-testa-subset/08-CONTEXT.md` D-09 — CMake `verify_pp_qt_freshness` ALL target deferred to Phase 9 (carry-forward #4)
- `.planning/phases/08-ab-validation-on-testa-subset/08-05-SUMMARY.md` "Two harness deviations" — `--out-dir` relative-path bug (carry-forward #4)

### File targets (deletion candidates)

- `xvue/xvuelc.c` — RETIRE-01 (after one-release-cycle A/B window closes)
- `bin/ccxvue` — RETIRE-01
- `bin/cb*` (every script linking against `libX11` or `/usr/X11R6/lib64`) — RETIRE-02
- `bin/convertepsgif` — RETIRE-03
- `xvue/video1.f`, `xvue/videofin.f`, `xvue/videonm.f` — RETIRE-03 (LVIDEO)
- 18+ Fortran tracer subroutines (e.g., `flui/trvi2d.f`, `ther/trplse.f`) — RETIRE-03
- `README`, `LISEZMOI`, install scripts — RETIRE-04

### Project guardrails (CLAUDE.md)

- `CLAUDE.md` §"Working rules" — `bin/cbl_tout` MUST exit 0 after every commit; testa cases MUST keep passing; build never breaks
- `CLAUDE.md` §"Active project goals" — Qt migration end-state target

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets

- **Phase 7 XvueExport `saveGifTo` ffmpeg path** — replaces LVIDEO + bin/convertepsgif. Already shipped, A/B-validated in Phase 8.
- **`bin/cbl_tout_qt`** — the Qt-only build entry. Becomes the canonical `bin/cbl_tout` after RETIRE-02 collapses the dual-build distinction.
- **Phase 8 testa baselines** — evidence/{case}-x11.png + evidence/{case}-qt-1x.png + sweep-logs. Phase 9 doesn't touch these (validation artifacts stay; the _binaries_ producing them go away).
- **Phase 8 harness scripts** — `bin/ab_compare_pair.sh`, `bin/ab_capture_x11.sh`, `bin/ab_sweep_phase8.sh`, `bin/phase8_case_batch_map.sh`. RETIRE-NN must not break these (they may need post-retirement updates: ab_capture_x11.sh becomes obsolete; ab_sweep_phase8.sh `--mode x11` branch goes away).

### Established Patterns

- **Atomic-commit retirement** — Phase 7 EXPORT-06 grep gate established the pattern of CMake ALL targets that fail-fast on regressed scope. RETIRE-NN plans add similar grep gates for the post-retirement scope (e.g., `verify_no_xvuelc_references_in_qt_build`).
- **Per-cell sign-off matrix** — Phase 8 D-10. RETIRE-NN closure can use a similar matrix to track per-RETIRE-ID acceptance.

### Integration Points

- **Phase 8 evidence/** stays read-only. Phase 9's "post-retirement testa A/B" produces NEW captures stored under `evidence/post-retire/` (or similar) — does NOT mutate Phase 8 artifacts.
- **`bin/cbl_tout` entry-point semantics:** after RETIRE-02 the script either becomes a Qt-only build OR is replaced by `bin/cbl_tout_qt` (deletion + rename). Decision deferred to planner per RETIRE-02 task split.

</code_context>

<specifics>
## Specific Ideas

- **Pre-Phase-9 git tag = v1.0-pre-retire** + branch retire-restore-point. Single-command rollback. Created BEFORE Phase 9 Plan 1 starts (D-08, D-09).
- **8-plan structure**: 4 RETIRE-NN plans (D-02 ordering) + 4 carry-forward plans (matched-dim recapture, ppnlse_qt deadlock fix, 3 P7 goldens, harness `--out-dir` bug + CMake freshness target).
- **LVIDEO retirement is wholesale** (D-05): video1.f + videofin.f + videonm.f + 18+ tracers + bin/convertepsgif all retire under RETIRE-03. Animation capability preserved via Phase 7 XvueExport `saveGifTo`.
- **Cross-tag operation for 9-07**: produces 3 Phase-7 deferred goldens by checking out `v1.0-pre-retire` in a separate worktree, running `bin/convertepsgif` there, copying artifacts to main. Single deviation from "Phase 9 only touches Qt-only main".

</specifics>

<deferred>
## Deferred Ideas

- **Phase-9-finalize summary** that catalogs deleted file count + Fortran subroutine count. Could be auto-generated by the gsd-verifier post-execute. Not in Phase 9 scope; future cleanup metric.
- **Qt 7 migration** — Phase 9 leaves Qt 6 as the single backend; Qt 7 is a future milestone (Qt 7 GA expected 2027+).
- **Test-suite reorganization** post-retirement (testa/ pruning, removal of any X11-only tests). Phase 9 leaves testa/ structure untouched; cleanup deferred.
- **Pre-Phase-9 release tag artifact bundling** (tarball, source archive) — outside scope; maintainer's release process. v1.0-pre-retire git tag is the only Phase-9-required artifact.
- **post-retirement README badge / version bump** (v1.0 → v1.1 or v2.0 to mark Qt-only era) — release-engineering decision; not Phase 9 scope.

### Reviewed Todos (not folded)

None — no pending todos matched Phase 9 scope.

</deferred>

---

*Phase: 09-retire-x11-backend*
*Context gathered: 2026-05-05*
