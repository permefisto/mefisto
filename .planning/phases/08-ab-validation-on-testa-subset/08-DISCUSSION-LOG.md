# Phase 8: A/B validation on testa subset - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-05-05
**Phase:** 08-ab-validation-on-testa-subset
**Areas discussed:** Evidence capture, Pass threshold, HiDPI method, OMP scope, Phase 7 deferred goldens, Build hygiene (Gap-A), CHECKLIST.md ship-gate format, Tolerance-band tool

---

## Evidence capture

| Option | Description | Selected |
|--------|-------------|----------|
| Screenshots in repo (Recommended) | PNG pairs committed under `.planning/phases/08-*/evidence/{case}-{backend}-{mode}.png` — reproducible, reviewable in PR diffs, ~5MB total. | ✓ |
| Diff images only | Run pairs through ImageMagick compare → only diff PNGs commit if non-trivial. Cheaper storage; harder to audit pass cases. | |
| Text notes only | CHECKLIST.md captures verdict + 1-line description. Maintainer eyeballs locally. No artifacts in repo. | |

**User's choice:** Screenshots in repo (Recommended).
**Notes:** Encoded as D-01 in CONTEXT.md. Path under `.planning/phases/08-ab-validation-on-testa-subset/evidence/`.

---

## Pass threshold

| Option | Description | Selected |
|--------|-------------|----------|
| Eyeball-only verdict (Recommended) | Maintainer signs off per cell. Subjective but matches CLAUDE.md "large/visual tests are user-run". Catches what users actually see. | |
| Tolerance band | Sub-2px geometry + delta-E < 5 color. Automated check via ImageMagick/skimage. Brittle on font rendering and AA. | ✓ |
| Pixel-perfect | Bit-exact match. Will fail (Qt anti-aliasing, font hinting). Rejecting up front. | |

**User's choice:** Tolerance band.
**Notes:** D-02. Eyeball-only was the recommended default but user picked tolerance band — automation preferred over per-cell judgement. Drives downstream tool choice (Tolerance-band tool area).

---

## HiDPI method

| Option | Description | Selected |
|--------|-------------|----------|
| QT_SCALE_FACTOR=2 under Xvfb (Recommended) | Cheap, reproducible, runnable headless. Captures the 2x scaling logic but not real 4K display pixels. | ✓ |
| Real 4K monitor | True HiDPI, manual maintainer step. Hardware-gated; CI cannot reproduce. | |
| Both | QT_SCALE_FACTOR=2 in CI plus one-off 4K eyeball check by maintainer. Best coverage; double the work. | |

**User's choice:** QT_SCALE_FACTOR=2 under Xvfb (Recommended).
**Notes:** D-04. Real-4K eyeball check moved to deferred ideas — maintainer ad-hoc, not a ship gate.

---

## OMP scope

| Option | Description | Selected |
|--------|-------------|----------|
| All 5 _OMP variants (Recommended) | Sweep mesher/elas/flui/ther/nlse _OMP where they exist. Proves main-thread invariant across every solver. More runs but more confidence. | ✓ |
| ELASTICER_OMP only | Strict VALID-03 read. One representative case (nafems_le1). Faster but doesn't catch per-module threading issues. | |
| Plus nlsecu_OMP | Two cases (elas + nlse) since they are the heaviest solvers. Compromise. | |

**User's choice:** All 5 _OMP variants (Recommended).
**Notes:** D-05. VALID-03's literal text mentions ELASTICER_OMP; CONTEXT.md broadens to all 5 with N-A cells where the _OMP binary doesn't exist.

---

## Tolerance-band tool

| Option | Description | Selected |
|--------|-------------|----------|
| ImageMagick compare (Recommended) | `compare -metric AE -fuzz 5%` already shipped on the dev host. Phase 9 retires ImageMagick from xvue/qt only, NOT from validation tooling. Reuse legacy. | ✓ |
| Python skimage SSIM | Adds a Python+pip dep to the build host. SSIM is structural, more forgiving than AA. Needs requirements.txt. | |
| Custom Qt diff | Write small QImage-based diff tool committed under xvue/qt/tools/. No new deps but more code to maintain. | |

**User's choice:** ImageMagick compare (Recommended).
**Notes:** D-03. Reusing legacy tool keeps the ImageMagick scope-guard contract (xvue/qt/ has no ImageMagick; bin/ + .planning/ legitimately do).

---

## Phase 7 deferred goldens

| Option | Description | Selected |
|--------|-------------|----------|
| Plan 1 prerequisite (Recommended) | First Phase 8 plan bootstraps the 3 goldens (compile scene01_driver.f, run X11+convertepsgif on testa/wave + cavity2d). Flips Phase 7 ctest QSKIPs → PASS before main A/B sweep starts. | ✓ |
| Folded into per-case A/B | wave/cavity2d goldens emerge naturally from running those testa cases through both backends. scene01 remains a separate small plan. | |
| Defer to Phase 9 | Skip golden bootstrap until X11 retirement closes the A/B window. Phase 7 ctest stays QSKIP indefinitely. | |

**User's choice:** Plan 1 prerequisite (Recommended).
**Notes:** D-06, D-07. Plan 1 procedure documented in CONTEXT.md with the three sub-steps from VERIFICATION.md §9.

---

## Build hygiene (Phase 7 Gap-A)

| Option | Description | Selected |
|--------|-------------|----------|
| First task: bin/cbl_tout_qt (Recommended) | Phase 8 Plan 1 runs cbl_tout_qt as explicit prerequisite. Guarantees pp/*_qt reflects current source. Adds CMake guard target later. | ✓ |
| CMake guard only | Add CMake target that checks pp/*_qt mtime vs libxvueqt.a; fails build if stale. Permanent fix; no Phase 8 manual step. | |
| Skip | Trust maintainer to rebuild before A/B runs. Document in CHECKLIST.md prerequisites section. | |

**User's choice:** First task: bin/cbl_tout_qt (Recommended).
**Notes:** D-08. CMake guard captured as deferred (D-09) — proper fix belongs in Phase 9 cleanup, not Phase 8.

---

## CHECKLIST.md ship-gate format

| Option | Description | Selected |
|--------|-------------|----------|
| Per-cell verdict + maintainer initials (Recommended) | Grid: rows = 5 cases, columns = {X11, Qt, Qt-OMP, Qt-HiDPI}. Each cell: PASS/FAIL/N-A + initials. v1 ship gate = all PASS or explicit override. | ✓ |
| Aggregate per-case | One row per case with overall verdict. Less granular but easier to read. | |
| Auto-pass via tolerance script | Tolerance-band script writes verdicts. Maintainer countersigns only the overall ship decision. | |

**User's choice:** Per-cell verdict + maintainer initials (Recommended).
**Notes:** D-10. Each cell carries PASS/FAIL/N-A + initials + path-to-evidence + AE pixel count from `compare`.

---

## Claude's Discretion

- Per-plan task decomposition (Plan 1 = bootstrap + freshness, Plan 2..N = A/B sweep batches per backend mode, Plan N+1 = CHECKLIST.md finalize). Planner picks the exact split.
- Screenshot capture script: shell wrapper using kwin-mcp / xdotool / Xvfb's xwd, or a small Qt-based recorder. Planner chooses based on what's reproducible.
- Xvfb screen resolution defaults: 1280x800 (matches Phase 7 UAT session).
- Tolerance-band command-line tuning (`-fuzz` percent, `-metric` choice AE vs PAE vs RMSE) — planner picks the exact invocation.

## Deferred Ideas

- Real-4K HiDPI eyeball check (beyond QT_SCALE_FACTOR=2 ship gate; maintainer ad-hoc when 4K hardware available).
- CMake `verify_pp_qt_freshness` ALL target (permanent fix for Gap-A; Phase 9 cleanup).
- CI integration (automated Xvfb runs in GitHub Actions / GitLab CI; out of scope for v1).
- Tolerance-band per-case overrides (e.g., cavity2d velocity arrows may need looser fuzz; capture as per-cell override note rather than global relaxation).
- SSIM/skimage-based diff (rejected to avoid Python+pip dependency).
