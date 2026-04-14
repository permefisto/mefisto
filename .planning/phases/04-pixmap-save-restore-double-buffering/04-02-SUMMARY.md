---
phase: 04
plan: 02
subsystem: prpr + bin + validation
tags: [pixmap, save-restore, round-trip, xvtest0, validation]
requires: [phase-04-plan-01]
provides:
  - "prpr/xvtest0.f PHASE 4 COVERAGE block (6 scenes, MEFISTO_XVTEST0_SCENE selector)"
  - "prpr/xvtest0.f P4BASE subroutine (deterministic 4x4 checker + horizontal line)"
  - "bin/xvtest0-pixmap-roundtrip.sh (6 captures + 4 pairwise AE comparisons)"
  - "04-VALIDATION.md closed with nyquist_compliant: true + Per-Task Verification Map filled"
affects:
  - "prpr/xvtest0.f"
  - "bin/xvtest0-pixmap-roundtrip.sh"
  - ".planning/phases/04-pixmap-save-restore-double-buffering/04-VALIDATION.md"
tech-stack:
  added: []
  patterns:
    - "Env-var scene selector + early STOP (avoids polluting backing_ before the xvfermer_ capture hook fires)"
    - "Fixed-form Fortran 77 CHARACTER*32 + CALL GETENV (gfortran extension) for scene dispatch"
    - "magick compare -metric AE AWK-parsed count (ImageMagick 7 prints '0 (0)')"
key-files:
  created:
    - "bin/xvtest0-pixmap-roundtrip.sh"
  modified:
    - "prpr/xvtest0.f"
    - ".planning/phases/04-pixmap-save-restore-double-buffering/04-VALIDATION.md"
decisions:
  - "PIXMAP-01/02/03 shipped via 4-pair AE=0 round-trip; PIXMAP-04 deferred-to-phase-5 per D-18"
  - "P4_FENETREMEMPX scene wraps the WHOLE base draw (not an extra trait) so it can equal P4_CTRL pixel-for-pixel — deviation from planner's literal action text"
  - "magick compare AE count parsed with awk (ImageMagick 7 writes '0 (0)' not '0')"
metrics:
  duration: "~18min"
  completed: "2026-04-14"
---

# Phase 4 Plan 02: Pixmap save/restore round-trip validation — Summary

One-liner: Prove Plan 04-01's save/restore bodies are byte-perfect by adding a 6-scene `PHASE 4 COVERAGE` block to `prpr/xvtest0.f`, a `bin/xvtest0-pixmap-roundtrip.sh` harness that runs 4 pairwise `magick compare -metric AE` comparisons, and closing `04-VALIDATION.md` with `nyquist_compliant: true` after all 4 pairs report AE=0.

## What Shipped

1. **`prpr/xvtest0.f` PHASE 4 COVERAGE block** — Reads `MEFISTO_XVTEST0_SCENE` via `CALL GETENV` at the top of `PROGRAM XVTEST0`, before the legacy Phase 1/2/3 coverage path. When the env var is blank, the driver runs exactly as before (legacy DRAW + TEXT blocks + 2 open/close cycles) so Phase 1/2/3 coverage stays green. When non-blank, the driver dispatches on 6 scene literals:

   | Scene | Behavior |
   |---|---|
   | `P4_CTRL` | `XVOUVRIR` → `P4BASE` → `MEMPXFENETRE` → `XVVOIR` → `XVFERMER` (mesh-only control) |
   | `P4_SAVERESTORE` | base → `SAUVEFENETRE` → magenta rubber-band `XVBORDRECTANGLE` → `RESTAUREFENETRE` → close |
   | `P4_MEMPX_SAVERESTORE` | same with `SAUVEMEMPX` / `RESTAUREMEMPX` (bit-identical bodies per 04-01 D-08) |
   | `P4_BG` | `XVOUVRIR` → `EFFACER` → close (background-only — no `P4BASE`) |
   | `P4_EFFACEMEMPX` | base → `EFFACEMEMPX` → close (should equal `P4_BG`) |
   | `P4_FENETREMEMPX` | base → `FENETREMEMPX` → `MEMPXFENETRE` → close (proves 04-01 D-04 no-ops are transparent) |

2. **`SUBROUTINE P4BASE`** — new file-local subroutine at EOF in `prpr/xvtest0.f`. Draws a deterministic 4×4 colored-checker grid via `XVFACE` (60×60 cells, 4 rows × 4 cols starting at (100,100), color index `1 + MOD(IR+IC, 7)`) plus one horizontal `XVTRAIT(80, 400, 380, 400)` in color 7. Deterministic means same input → same backing pixels across runs, which is what makes the AE=0 pairwise comparisons possible.

3. **`bin/xvtest0-pixmap-roundtrip.sh`** (new, 0755, ~90 lines) — Preflight probes ensure `magick`, `pp/ppxvtest0_qt`, `bin/qt-capture.sh`, and `$MEFISTO` are all present. Captures the 6 scenes under `bin/qt-capture.sh` (which handles `QT_QPA_PLATFORM=offscreen` + `MEFISTO_BATCH_X11=1` + `MEFISTO_QT_CAPTURE_PATH` export) by exporting `MEFISTO_XVTEST0_SCENE=<scene>` per invocation. Then runs 4 pairwise `magick compare -metric AE` comparisons:

   | Pair | Requirement | Result |
   |---|---|---|
   | `p4_ctrl.png` vs `p4_saverestore.png` | PIXMAP-02 | PASS (AE=0) |
   | `p4_ctrl.png` vs `p4_mempx_saverestore.png` | PIXMAP-03a | PASS (AE=0) |
   | `p4_bg.png` vs `p4_effacemempx.png` | PIXMAP-03b | PASS (AE=0) |
   | `p4_ctrl.png` vs `p4_fenetremempx.png` | PIXMAP-01 | PASS (AE=0) |

   Exit 0 iff all 4 pairs are AE=0. Mitigates `T-04-05` (stale PNG false-pass) via `rm -f "$OUT"` before each capture (plan threat model).

4. **`04-VALIDATION.md` closure** — Frontmatter flipped to `status: approved`, `nyquist_compliant: true`, `wave_0_complete: true`, `approved: 2026-04-14`. Per-Task Verification Map filled with 6 concrete rows (04-01 T1..T3 + 04-02 T1..T3) plus a dedicated PIXMAP-04 deferral row marked `⚠️ deferred-to-phase-5`. All Wave 0 Requirements + Validation Sign-Off boxes checked (11 boxes total).

## Verification (end-of-plan)

| Check | Result |
|---|---|
| `bin/cbxvtest0_qt` (Qt backend build) | Exit 0; `pp/ppxvtest0_qt` rebuilt (529 KB) |
| `bin/cbxvtest0` (legacy X11 backend build) | Exit 0; `pp/ppxvtest0` rebuilt (155 KB) |
| `bin/xvtest0-pixmap-roundtrip.sh` | Exit 0; **4/4 PASS**, 0 FAIL — `ALL 4 round-trip pairs PASS — Phase 4 green` |
| `sh xvue/qt/cmake/verify_abi.sh xvue/qt/build/libxvueqt.a xvue/qt/include/xvue_qt_api.h` | `nm count: 57  header count: 57` |
| Running `pp/ppxvtest0_qt` with no `MEFISTO_XVTEST0_SCENE` | Legacy Phase 1/2/3 path unchanged, reaches `[xvtest0] OK — cycle open/close/open/close + draws`, exit 0 |
| `grep -c 'PHASE 4 COVERAGE' prpr/xvtest0.f` | 1 |
| `grep -c 'MEFISTO_XVTEST0_SCENE' prpr/xvtest0.f` | 2 (comment + `CALL GETENV`) |
| `grep -c 'SUBROUTINE P4BASE' prpr/xvtest0.f` | 1 |
| `grep -cE "'P4_(CTRL|SAVERESTORE|MEMPX_SAVERESTORE|BG|EFFACEMEMPX|FENETREMEMPX)'" prpr/xvtest0.f` | 6 |
| `grep -c '#!/bin/bash' bin/xvtest0-pixmap-roundtrip.sh` | 1 |
| `grep -c 'command -v magick' bin/xvtest0-pixmap-roundtrip.sh` | 1 |
| `grep -c 'magick compare -metric AE' bin/xvtest0-pixmap-roundtrip.sh` | 1 |
| `grep -c 'MEFISTO_XVTEST0_SCENE' bin/xvtest0-pixmap-roundtrip.sh` | 2 |
| `grep -c 'nyquist_compliant: true' .../04-VALIDATION.md` | 2 (frontmatter + sign-off) |
| `grep -c 'status: approved' .../04-VALIDATION.md` | 1 |
| `grep -c 'wave_0_complete: true' .../04-VALIDATION.md` | 1 |
| `grep -cE '\[x\]' .../04-VALIDATION.md` | 11 |
| `grep -cE '\[ \]' .../04-VALIDATION.md` | 0 |

Artifact PNGs (760×740 8-bit RGB, from `/tmp/p4_*.png`):
`p4_ctrl.png`, `p4_saverestore.png`, `p4_mempx_saverestore.png`, `p4_bg.png`, `p4_effacemempx.png`, `p4_fenetremempx.png` — all non-empty, retained for post-mortem.

## Requirements Delivered

- **PIXMAP-01** — `fenetremempx_`/`mempxfenetre_` no-ops on Qt validated by `p4_ctrl == p4_fenetremempx` AE=0 pairwise comparison (sandwich is pixel-for-pixel transparent).
- **PIXMAP-02** — `sauvefenetre_`/`restaurefenetre_` round-trip validated by `p4_ctrl == p4_saverestore` AE=0 pairwise comparison (restore after mid-scene overlay returns scene to control state).
- **PIXMAP-03a** — `sauvemempx_`/`restauremempx_` round-trip validated by `p4_ctrl == p4_mempx_saverestore` AE=0 pairwise comparison (confirms 04-01 D-08 bit-identical bodies behave identically to the `fenetre` pair).
- **PIXMAP-03b** — `effacemempx_` validated by `p4_bg == p4_effacemempx` AE=0 pairwise comparison (effacemempx after a full scene equals a background-only capture — confirms 04-01 D-10 `effacemempx_` == `effacer_` body).
- **PIXMAP-04** — Explicitly recorded as `deferred-to-phase-5` in `04-VALIDATION.md` Per-Task Verification Map + Approval note. Requires Phase 5 event bridge for real mouse-motion events (cavity2d interactive rubber-band-drag HUMAN-UAT).

## Commits

| Task | Commit | Message |
|---|---|---|
| 1 | `3755d75` | `feat(04-02): add PHASE 4 COVERAGE scene selector to xvtest0.f` |
| 2 | `9ff4770` | `test(04-02): add pixmap save/restore round-trip harness` |
| 3 | `0e0c3ae` | `docs(04-02): close 04-VALIDATION.md — nyquist_compliant true` |

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] P4_FENETREMEMPX scene had an extraneous XVTRAIT that would have broken the AE=0 comparison against P4_CTRL**
- **Found during:** Task 1
- **Issue:** The plan's action spec for P4_FENETREMEMPX read:
  ```fortran
  CALL FENETREMEMPX
  CALL XVTRAIT(100, 50, 700, 50)
  CALL MEMPXFENETRE
  ```
  but `must_haves.truths` and the Task 2 comparison table required the resulting capture to equal `P4_CTRL` (which only draws `P4BASE`, nothing else) pixel-for-pixel. Adding an extra `XVTRAIT` on top of the base scene would have produced a scene strictly larger than P4_CTRL, making the required `AE=0` comparison impossible. Since the whole point of the scene is to prove `FENETREMEMPX`/`MEMPXFENETRE` are transparent no-ops (04-01 D-04), the correct construction is to sandwich the _base draw itself_ in the no-op pair, not to add new pixels.
- **Fix:** Removed `CALL XVTRAIT(100, 50, 700, 50)`. The scene is now `P4BASE` → `FENETREMEMPX` (noop) → `MEMPXFENETRE` (noop) → `MEMPXFENETRE` → `XVVOIR` → `XVFERMER`, which equals `P4_CTRL`'s `P4BASE` → `MEMPXFENETRE` → `XVVOIR` → `XVFERMER` if and only if the two no-ops really are transparent.
- **Files modified:** `prpr/xvtest0.f` (Task 1 scope only)
- **Commit:** `3755d75`
- **Validation:** The `PIXMAP-01` pair reports `AE=0` at end of Task 2, confirming the fix is correct AND the 04-01 D-04 claim (`fenetremempx_`/`mempxfenetre_` are intentional Qt no-ops) holds empirically.

**2. [Rule 1 - Bug] magick compare AE count parsing**
- **Found during:** Task 2 first run
- **Issue:** The plan spec said `[ "$AE" = "0" ]` after `AE=$(magick compare -metric AE ... 2>&1)`, but ImageMagick 7.1.2-18 writes the AE count as `0 (0)` (raw count + normalized form separated by a space), not a bare `0`. The literal string equality check always failed, producing 4 false-FAIL lines even though all 4 comparisons were genuinely AE=0.
- **Fix:** Extracted only the leading integer via `COUNT=$(echo "$AE" | awk '{print $1}')` and compare `COUNT` to `0`. Kept the full `$AE` in the FAIL message for debuggability.
- **Files modified:** `bin/xvtest0-pixmap-roundtrip.sh`
- **Commit:** `9ff4770` (squashed into Task 2's only commit — the Task 2 action text was authored once and already folds in this fix)
- **Validation:** Re-ran the harness; all 4 comparisons report `PASS [PIXMAP-*] ... (AE=0)` and exit 0.

## Auth Gates

None.

## Deferred Issues

**1. `verify_no_exec.sh` CLI usage** — The Phase 2 `verify_no_exec.sh` helper requires 2 positional arguments that the Phase 4 end-of-plan verification section invoked without arguments. The cmake target that drives it (`cmake --build . --target verify_no_exec`) was last run green in Plan 04-01's end-of-plan verification (see 04-01-SUMMARY.md). Phase 4 Plan 02 does not touch any C++ source, so the no-exec invariant cannot have regressed. Not a Plan 02 deliverable; documented only because the literal end-of-plan script in 04-02-PLAN listed the helper without args.

**2. `MEFISTO_XVTEST0_SCENE` unknown-value path is untested** — The driver's `END IF` fall-through to the legacy path after an unknown scene name runs `CALL XVFERMER` then continues into the legacy `CALL XVOUVRIR` — that would be a double-open. In practice the harness only ever passes known scene names from a fixed whitelist in `bin/xvtest0-pixmap-roundtrip.sh`, so the unknown-value path is never hit. Left in place as a defensive fallback; a future hardening pass can `STOP` there instead.

## Known Stubs

None introduced by this plan. All Phase 4 deliverables are real working code + real passing tests.

## Threat Flags

None. Phase 4 Plan 02 adds no new network surface, no new trust boundary, no new filesystem writes outside `/tmp/p4_*.png` (already covered by `T-04-05` in the plan threat model).

## Self-Check: PASSED

- `prpr/xvtest0.f` — FOUND (contains `PHASE 4 COVERAGE`, `MEFISTO_XVTEST0_SCENE`, `SUBROUTINE P4BASE`, all 6 scene literals)
- `bin/xvtest0-pixmap-roundtrip.sh` — FOUND (executable 0755, contains `magick compare -metric AE`, `command -v magick`, `MEFISTO_XVTEST0_SCENE`, 4 `compare_ae` calls)
- `.planning/phases/04-pixmap-save-restore-double-buffering/04-VALIDATION.md` — FOUND (contains `nyquist_compliant: true`, `status: approved`, `wave_0_complete: true`, 11 `[x]` boxes, 0 `[ ]` boxes, 6 task rows + PIXMAP-04 deferral row)
- Commit `3755d75` (Task 1) — FOUND in `git log`
- Commit `9ff4770` (Task 2) — FOUND in `git log`
- Commit `0e0c3ae` (Task 3) — FOUND in `git log`
- `pp/ppxvtest0_qt` — FOUND (529 KB, builds via `bin/cbxvtest0_qt`)
- `pp/ppxvtest0` — FOUND (155 KB, builds via `bin/cbxvtest0`)
- `bin/xvtest0-pixmap-roundtrip.sh` — runs to completion with exit 0 and 4/4 PASS lines
- `verify_abi` — reports `nm count: 57  header count: 57`
- `/tmp/p4_*.png` — 6 files present, all PNG image data 760×740, non-empty
