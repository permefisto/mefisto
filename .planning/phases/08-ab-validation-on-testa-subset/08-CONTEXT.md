# Phase 8: A/B validation on testa subset - Context

**Gathered:** 2026-05-05
**Status:** Ready for planning

<domain>
## Phase Boundary

Verify the 5 canonical `testa/` cases (BUILD-10 baseline: pan2d, nafems_le1,
cavity2d, heat1d, nlsecu) render and behave identically on the X11 backend
versus the Qt 6 backend, including HiDPI (`QT_SCALE_FACTOR=2`) and `_OMP`
sweeps. Phase 8 produces `08-CHECKLIST.md` — the explicit ship gate that
declares v1 release-ready and opens the one-release-cycle A/B window that
gates Phase 9 (X11 retirement).

In scope:
- A/B run of all 5 cases under both backends (pp/pp* X11 and pp/pp*_qt Qt)
- Tolerance-band image diffing of canvas screenshots per (case, backend, mode) cell
- HiDPI sweep via `QT_SCALE_FACTOR=2` headless under Xvfb
- All 5 _OMP variants where they exist (mesher / elas / flui / ther / nlse)
- Bootstrap of the 3 Phase-7-deferred goldens (scene01.eps, wave_legacy.gif, cavity2d_legacy.gif) as Plan 1 prerequisite
- Color-bar spot checks on nafems_le1, heat1d, cavity2d (VALID-05)
- Font-metric spot checks on pan2d, hexahedron (VALID-06)
- 08-CHECKLIST.md sign-off matrix as v1 ship gate

Out of scope:
- Code changes inside xvue/qt/ to fix any drift — those become gap-closure plans, but only for drift the maintainer accepts as a real defect
- Retirement of xvuelc.c, bin/convertepsgif, libX11 linker lines (Phase 9 RETIRE-01..04)
- LVIDEO pipeline retirement (xvue/video1.f, videofin.f, videonm.f — Phase 9 RETIRE-03)
- New testa cases beyond the 5 canonical baseline

</domain>

<decisions>
## Implementation Decisions

### Evidence capture

- **D-01:** Side-by-side screenshots committed to repo under
  `.planning/phases/08-ab-validation-on-testa-subset/evidence/{case}-{backend}-{mode}.png`.
  Reproducible, reviewable in PR diffs, expected ~5 MB total. PNG only —
  no JPEG (lossy artifacts would defeat tolerance band).

### Pass threshold

- **D-02:** Tolerance band, NOT eyeball-only or pixel-perfect. Geometry
  drift cap: sub-2 pixels. Color drift cap: delta-E < 5. Pixel-perfect
  rejected up front (Qt anti-aliasing + font hinting differ from X11).
- **D-03:** Tolerance enforcement tool: ImageMagick `compare -metric AE
  -fuzz 5%`. The dev host already ships ImageMagick. Phase 9 RETIRE-03
  removes ImageMagick from `xvue/qt/` ONLY — not from validation tooling
  under `bin/` or `.planning/`. Reusing legacy `compare` is in-scope.

### HiDPI methodology

- **D-04:** `QT_SCALE_FACTOR=2` under `xvfb-run --auto-servernum`. Cheap,
  reproducible, runnable headless on the build host. Captures the 2x
  scaling logic but not real 4K display pixels. Real-4K eyeball check is
  NOT a Phase 8 gate — recorded as a deferred idea for the maintainer's
  ad-hoc validation.

### OMP scope

- **D-05:** Sweep all 5 `_OMP` variants where they exist on disk:
  mesher / elas / flui / ther / nlse. VALID-03's literal text says
  "ELASTICER_OMP" — that is the Phase-7-end smoke check; Phase 8 broadens
  to all 5 for full main-thread-invariant coverage. Cells where no `_OMP`
  variant exists for a module are marked `N-A` in the checklist.

### Phase 7 deferred-golden integration

- **D-06:** Plan 1 of Phase 8 bootstraps the 3 Phase-7-deferred goldens
  before any A/B sweep runs:
    1. Compile `xvue/qt/tests/golden/scene01_driver.f` against X11 +
       xvuelc.o under Xvfb, capture `TEMPORAIRE.EPS`, commit as
       `xvue/qt/tests/golden/scene01.eps`.
    2. Run `testa/wave` through X11 + `bin/convertepsgif`, commit
       `xvue/qt/tests/golden/wave_legacy.gif`.
    3. Run `testa/cavity2d` through X11 + `bin/convertepsgif`, commit
       `xvue/qt/tests/golden/cavity2d_legacy.gif`.
- **D-07:** After Plan 1 commits the goldens, re-run
  `ctest -R 'xvue_qt_(postscript|export)_tests'` and confirm all 3
  Phase-7 QSKIPs flip to PASS. This is the Phase-7-close gate folded
  into Phase 8 entry.

### Build hygiene (Phase 7 Gap-A)

- **D-08:** Plan 1 first task runs `bin/cbl_tout_qt` from a clean tree
  to refresh `pp/*_qt`. Phase 7 Gap-A documented that the original
  verifier ran `cmake --build xvue/qt/build` + `bin/cbl_tout` but never
  `bin/cbl_tout_qt`, so `pp/*_qt` could be stale relative to
  `libxvueqt.a`. A/B sweep MUST run against fresh `pp/*_qt`.
- **D-09:** Permanent CMake guard target proposed but NOT in Phase 8
  scope: a CMake `verify_pp_qt_freshness` target that fails the build
  when `libxvueqt.a` mtime > any `pp/*_qt` mtime. Tracked as deferred
  idea for Phase 9 cleanup.

### CHECKLIST.md sign-off shape

- **D-10:** Per-cell verdict matrix. Rows = 5 canonical cases. Columns =
  {X11 baseline, Qt 1x, Qt HiDPI 2x, Qt _OMP}. Each cell records:
  PASS / FAIL / N-A + maintainer initials + path-to-evidence-png +
  ImageMagick AE pixel count. v1 ship gate = every cell in {PASS, N-A}
  with maintainer signature on the bottom-row sign-off line. Any FAIL
  blocks ship until either resolved (gap-closure plan) or explicitly
  overridden with documented rationale.

### Claude's Discretion

- Per-plan task decomposition (Plan 1 = bootstrap + freshness, Plan 2..N
  = A/B sweep batches per backend mode, Plan N+1 = CHECKLIST.md
  finalize). Planner picks the exact split.
- Screenshot capture script: shell wrapper using kwin-mcp /
  xdotool / Xvfb's xwd, or a small Qt-based recorder. Planner chooses
  based on what's reproducible across backends.
- Xvfb screen resolution defaults: 1280x800 (matches Phase 7 UAT
  session). Planner can adjust per-case if a test has hard-coded
  larger window expectations.
- Tolerance-band command-line tuning (`-fuzz` percent, `-metric` choice
  AE vs PAE vs RMSE). D-02 specifies the gate; planner picks the exact
  invocation.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Roadmap and requirements

- `.planning/ROADMAP.md` §"Phase 8: A/B validation on testa subset" — phase goal + success criteria + VALID-01..07 mapping
- `.planning/REQUIREMENTS.md` §VALID-01..07 — the seven validation requirements that bind this phase
- `.planning/REQUIREMENTS.md` §BUILD-10 — locks the 5-case baseline
- `.planning/validation/BASELINE.md` — immutable list of the 5 canonical testa cases (pan2d / nafems_le1 / cavity2d / heat1d / nlsecu)
- `.planning/PROJECT.md` — Key Decisions table + project-level invariants

### Phase 7 hand-off (mandatory read)

- `.planning/phases/07-image-gif-and-postscript-export/VERIFICATION.md` §9 — three human-bootstrap items merged into Phase 8 Plan 1 (D-06)
- `.planning/phases/07-image-gif-and-postscript-export/07-UAT.md` §"Gap-A" — pp/*_qt build-hygiene gap (D-08)
- `.planning/phases/07-image-gif-and-postscript-export/VALIDATION-LOG.md` — manual A/B sign-off ledger; Phase 8 Plan 1 fills the DEFERRED rows
- `.planning/phases/07-image-gif-and-postscript-export/07-06-SUMMARY.md` §"Plan 06 Task 3" — bootstrap procedure for the 3 goldens

### Tooling

- `xvue/qt/tests/golden/scene01_driver.f` — deterministic Fortran scene driver (Phase 7 Plan 06 deliverable)
- `bin/convertepsgif` — legacy ImageMagick wrapper used to produce wave/cavity2d_legacy.gif baselines
- `bin/cbl_tout_qt` — full Qt build script that relinks pp/*_qt (Phase 7 Gap-A discovery)
- `bin/test_no_imagemagick_in_qt.sh` — scope guard: ImageMagick legitimately stays under bin/ for validation, must NOT appear under xvue/qt/

### Codebase maps

- `.planning/codebase/TESTING.md` — existing test conventions
- `.planning/codebase/ARCHITECTURE.md` — backend split between X11 (xvuelc.c) and Qt (xvue/qt/)
- `.planning/codebase/INTEGRATIONS.md` — external tool deps (ImageMagick, ffmpeg, Xvfb)

### Project guardrails (CLAUDE.md)

- `CLAUDE.md` §"Working rules" — large/visual tests are user-run; build must never break; xvuelc.c byte-identical until Phase 9
- `CLAUDE.md` §"Active project goals" — A/B window is the gate, not a date

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets

- **kwin-mcp session_start + accessibility_tree + screenshot** —
  proven in Phase 7 UAT for driving live Qt windows headlessly. Phase 8
  Plan 2+ can reuse the same session pattern to drive each (case, mode)
  cell.
- **Xvfb + xvfb-run --auto-servernum** — already used by every Qt
  ctest target. Drop-in for headless A/B runs.
- **bin/convertepsgif** — legacy X11+ImageMagick GIF pipeline. Reused
  for D-06 Plan 1 baseline production; not modified.
- **xvue/qt/tests/golden/scene01_driver.f** — Phase 7's Fortran scene
  driver, already deterministic. Plan 1 step 1 inputs.
- **ImageMagick `compare`** — already on dev host. Tolerance-band
  enforcer per D-02, D-03.

### Established Patterns

- **Bootstrap-as-Plan-1, sweep-as-Plan-N pattern** — D-06 mirrors
  Phase 7 Plan 06's "harness in autonomous portion, golden bootstrap
  in human checkpoint" split. Plan 1 of Phase 8 IS that human
  checkpoint, run head-on.
- **Per-cell verdict + maintainer sign-off** — VERIFICATION.md §9 in
  Phase 7 already used the "Required / Why human / Expected" triple.
  Phase 8 CHECKLIST.md keeps the same three columns per cell.
- **Test scope respects subsystem boundary** — EXPORT-06 grep gate
  scopes to xvue/qt/ only; Phase 8 inherits that contract (ImageMagick
  under bin/ is allowed; under xvue/qt/ is not).

### Integration Points

- **Phase 9 entry gate** — CHECKLIST.md sign-off is the literal start
  signal for the one-release-cycle A/B window. Phase 9 RETIRE-01..04
  cannot start until this file is signed.
- **CI loop** — current .planning/ workflow has no continuous CI;
  Phase 8 evidence captures (PNGs + checklist) ARE the CI artifact
  for v1.

</code_context>

<specifics>
## Specific Ideas

- **5 canonical cases are immutable** (per BUILD-10): pan2d (mesher),
  nafems_le1 (elas), cavity2d (flui), heat1d (ther), nlsecu (nlse).
  Phase 8 cannot extend or replace them.
- **HiDPI = QT_SCALE_FACTOR=2 only** under Xvfb. No real 4K hardware
  test required for the ship gate.
- **Tolerance band**: sub-2px geometry + delta-E < 5 color, enforced by
  `compare -metric AE -fuzz 5%`. AE pixel count goes into the cell.
- **OMP sweep across all 5 modules** where the `_OMP` binary exists.
  N-A cells permitted where it doesn't.
- **Phase 7 Gap-A** is fixed by D-08 (cbl_tout_qt as Plan 1 first task)
  AND noted as a permanent CMake-guard candidate (D-09, deferred).

</specifics>

<deferred>
## Deferred Ideas

- **Real-4K HiDPI eyeball check** — beyond Phase 8's QT_SCALE_FACTOR=2
  ship gate. Maintainer ad-hoc validation when 4K hardware is
  available. Not blocking.
- **CMake `verify_pp_qt_freshness` ALL target** — permanent fix for
  Phase 7 Gap-A. Belongs in Phase 9 cleanup (or earlier hot-fix
  outside Phase 8 scope).
- **CI integration** — automated Xvfb runs in GitHub Actions / GitLab
  CI. Out of scope for v1 ship; Phase 8 evidence is committed
  artifacts.
- **Tolerance-band per-case overrides** — some cases (notably
  cavity2d's velocity arrows) may legitimately need a looser fuzz than
  5%. If discovered during Plan 2+, capture as a CHECKLIST.md
  per-cell override note rather than relaxing the global gate.
- **SSIM/skimage-based diff** — alternative tolerance metric.
  Considered and rejected (D-03) for Phase 8 to avoid Python+pip
  dependency.

### Reviewed Todos (not folded)

None — no pending todos matched Phase 8 scope.

</deferred>

---

*Phase: 08-ab-validation-on-testa-subset*
*Context gathered: 2026-05-05*
