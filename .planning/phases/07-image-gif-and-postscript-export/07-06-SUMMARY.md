---
phase: 07-image-gif-and-postscript-export
plan: 06
subsystem: xvue/qt
tags: [qt6, validation, golden, byte-parity, grep-gate, phase-close, lvideo-deferral, autonomous-partial]

requires:
  - phase: 07-02-postscript-emitter-state-machine
    provides: PsEmitter::handleLasops body — Plan 06 wraps it in a byte-level golden compare slot.
  - phase: 07-03-postscript-emit-helpers
    provides: 12 PsEmitter helper bodies with verbatim xvuelc.c format strings — Plan 06 ships the per-primitive golden coverage and the full-scene byte-level golden harness.
  - phase: 07-04-png-jpeg-pdf-export
    provides: XvueExport static class baseline. Plan 06 does not touch XvueExport itself.
  - phase: 07-05-gif-ffmpeg-fallback
    provides: XvueExport animation surface (begin/end/captureFrame/saveGifTo). Plan 06 ships the GIF A/B test harness.

provides:
  - "xvue/qt/tests/golden/scene01_driver.f — deterministic Fortran scene driver whose X11+xvuelc.c output materializes the byte-level golden scene01.eps in Plan 06 Task 3 (human checkpoint)"
  - "Test 15 PsEmitter_postscriptVerbatim_golden — full-scene byte-for-byte compare against scene01.eps; QSKIPs cleanly until the golden lands"
  - "Test 16 PsEmitter_perPrimitive_golden — per-helper byte-substring assertions (always runs; independent of the .eps golden file)"
  - "XvueExport_gif_AB_compare_wave / _cavity2d — frame-count + readability gates against legacy X11+convertepsgif baseline GIFs; QSKIP cleanly until baselines committed"
  - "bin/test_no_imagemagick_in_qt.sh — standalone EXPORT-06 grep-gate enforcer; exit 0 / 1 / 2 (env err); auto-detects repo root if MEFISTO unset"
  - "verify_no_imagemagick_in_qt CMake ALL target — fails the build if convert/ImageMagick/magick tokens appear under xvue/qt/"
  - "xvue/qt/README.md Phase 7 section (~125 lines) — Surfaces / PsEmitter byte-for-byte contract / Y-flip rule / HiDPI math (Pitfall 7 logical-vs-physical) / Threading invariant (XVUE_QT_ASSERT_MAIN_THREAD) / GIF probe-driven dispatch / LVIDEO Phase 9 deferral / Frame caps T-07-03 / Runtime deps / Test inventory / XVUE_ANIM=1 doc"
  - "VALIDATION-LOG.md — manual A/B sign-off ledger; current state DEFERRED for testa/wave + testa/cavity2d (autonomous executor with no human terminal); automated verifications all PASS"

affects: [phase-close]

tech-stack:
  added: []
  patterns:
    - "Defer-to-human-via-QSKIP: golden-file-dependent test slots compile and ship green via QSKIP with explicit pointer-to-procedure messages. When the human commits the golden file, slots flip from QSKIP to PASS-required automatically. No silent skipping — every SKIP carries the bootstrap procedure inline."
    - "Auto-detect repo root in shell scripts: bin/test_no_imagemagick_in_qt.sh falls back to walking up from $(dirname $0) when MEFISTO is unset; usable from CMake custom_target without environment plumbing."
    - "verify_no_*-as-CMake-ALL-target: EXPORT-06 grep gate joins verify_abi / verify_no_exec / verify_shortcut_modifiers / verify_icon_source as a build-time guard. Adding the target costs ~50ms per build and catches every drift attempt."

key-files:
  created:
    - xvue/qt/tests/golden/scene01_driver.f (117 LOC — Fortran 77 fixed-form scene driver; deterministic constant input; bootstrap procedure documented in header)
    - xvue/qt/tests/golden/ (NEW directory — will hold scene01.eps + wave_legacy.gif + cavity2d_legacy.gif post-Task 3)
    - bin/test_no_imagemagick_in_qt.sh (76 LOC — executable; auto-detects repo root; allowlists Qt API tokens)
    - .planning/phases/07-image-gif-and-postscript-export/VALIDATION-LOG.md (manual A/B sign-off ledger; current state DEFERRED for the human-action steps)
    - .planning/phases/07-image-gif-and-postscript-export/07-06-SUMMARY.md (this file)
  modified:
    - xvue/qt/tests/test_xvue_qt_postscript.cpp (+171 LOC — Test 15 + Test 16 byte-level golden slots)
    - xvue/qt/tests/test_xvue_qt_export.cpp (+93 LOC, -8 LOC — XvueExport_gif_AB_compare_wave + _cavity2d slots; comment scrub for grep-gate compliance)
    - xvue/qt/CMakeLists.txt (+18 LOC — verify_no_imagemagick_in_qt ALL target)
    - xvue/qt/README.md (+125 LOC — Phase 7 section)

key-decisions:
  - "Plan 06 Task 1 ships the harness; Task 3 ships the binaries. The byte-level golden file (scene01.eps) and the GIF baselines (wave/cavity2d_legacy.gif) require running the X11 backend on real testa/ cases, which is a manual procedure (CLAUDE.md §Tests: large/long-running tests are user-run). Plan 06 ships the test slots as QSKIPs that flip to PASS-required automatically when the golden binaries land. The harness is fully tested today; only the human-A/B verdict + golden binaries are deferred."
  - "Headless-executor handling: when the autonomous agent has no attached human terminal, the manual A/B sign-off is documented as DEFERRED in VALIDATION-LOG.md rather than erroring out. The orchestrator routes this as a verification gap during phase verification. The QSKIP test slots are the orthogonal proof that the gap is visible: any future executor running ctest sees the SKIP banner and knows the gate is incomplete."
  - "EXPORT-06 grep-gate scope is xvue/qt/ ONLY (D-16/D-17). The script and CMake target are precise about this — they grep ONLY xvue/qt/, never xvue/. The LVIDEO pipeline (xvue/video1.f, videofin.f, videonm.f) and bin/convertepsgif legitimately remain in the legacy tree until Phase 9 RETIRE-03. README Phase 7 section explicitly documents this deferral so future agents do not re-discover the LVIDEO pipeline and over-delete."
  - "Allowlist policy in the grep gate: \\b(convert|ImageMagick|magick)\\b matches via word boundary anchors, then a negative-filter pipeline strips legitimate Qt API tokens (QPageSize, convertToOther, convertTo(, QString::convertTo, convertFromUtf, Qt::ConvertibleTo). The two-step pipeline keeps the gate precise without weakening it. Documentation comments in the new test files / CMake / README must avoid the literal forbidden tokens — caught + scrubbed in the commit chain (Rule 1/2 deviation)."
  - "Independent inline per-primitive coverage (Test 16): not every primitive is in scene01.eps' driver scene; the per-primitive slot inlines its own short scenes whose expected byte substrings are taken directly from the Plan 03 Format-String Parity Table. This means Test 16 always runs (no QSKIP) and provides finer-grained drift detection if Test 15 ever flags."

threat-mitigations:
  - "T-07-06 (LVIDEO scope creep): grep gate scope is xvue/qt/ only. Verified by running the gate today: PASS, exit 0. xvue/video1.f / videofin.f / videonm.f / bin/convertepsgif are NOT modified by Plan 06 (or any prior Phase 7 plan). Phase 9 RETIRE-03 owns their retirement — explicitly documented in README.md."
  - "T-07-07 (PsEmitter byte-output regression): byte-level golden compare ships in two layers — full-scene Test 15 (against scene01.eps; QSKIPs until Task 3) and per-primitive Test 16 (always runs; uses inline format-string assertions from Plan 03's Format-String Parity Table). A failing diff is a hard test failure; no manual-eye fallback."
  - "T-07-08 (xvue/xvuelc.c regression): bin/cbl_tout exits 0 on the worktree; all 13 X11 pp/* binaries produced. git diff HEAD~5 HEAD -- xvue/xvuelc.c is empty (legacy untouched)."

requirements:
  - EXPORT-03 (animated GIF — A/B test harness shipped; manual sign-off DEFERRED)
  - EXPORT-06 (no ImageMagick in xvue/qt/ — grep gate live in CMake + standalone script)

duration: ~32m (autonomous Task 1 + Task 2 + verification + SUMMARY)

build:
  qt: "cmake --build xvue/qt/build -j4 → exit 0; verify_no_imagemagick_in_qt target green"
  abi: "verify_abi → 58 symbols (unchanged)"
  x11-full: "bin/cbl_tout → exit 0; 13 pp/pp* X11 binaries built (ppelas, ppflui, ppinit, ppmail, ppnlse, pppoba, ppther, ppxvtest0..4, pxyz)"
  qt-test-summary: "xvue_qt_postscript_tests 18 PASS + 1 SKIP; xvue_qt_export_tests 16 PASS + 2 SKIP; all other Qt test binaries pass at pre-Plan-06 rate (one pre-existing testPerModuleGroupIsolation failure documented in Plan 02/03/05 SUMMARYs)"

next-steps:
  - "Plan 06 Task 3 (human checkpoint): bootstrap scene01.eps + wave_legacy.gif + cavity2d_legacy.gif from the X11 backend; eyeball-compare Qt vs X11 GIFs on testa/wave + testa/cavity2d; record verdict in VALIDATION-LOG.md; commit goldens; re-run ctest to flip the QSKIPs to PASS."
  - "Phase 9 RETIRE-03 (future): retire xvue/video1.f / videofin.f / videonm.f, bin/convertepsgif, and xvue/xvuelc.c after the one-release A/B window closes."
---

# Phase 7 Plan 06: Validation A/B + EXPORT-06 Grep Gate Summary

**Phase 7 close-out scaffolding — byte-level golden test harness, EXPORT-06 grep gate live in CMake + standalone script, README Phase 7 section ships HiDPI math + threading invariant + LVIDEO deferral, manual A/B sign-off DEFERRED for orchestrator routing (autonomous executor without human terminal).**

## Performance

- **Duration:** ~32 minutes wall (autonomous portion — Task 1 scaffold + Task 2 grep-gate + Task 2 README + automated verifications + SUMMARY)
- **Started:** 2026-05-04 (Wave 5 entry, post Plan 05 completion)
- **Tasks:** 2 of 3 autonomous; Task 3 is `checkpoint:human-verify` and is DEFERRED with documented gap (no human terminal attached)
- **Commits:** 2 (Task 1 test scaffold; Task 2 grep gate + README)
- **Files:** 5 modified, 4 created

## Accomplishments

- **EXPORT-06 deliverable shipped:** the literal `convert` / `ImageMagick` / `magick` tokens cannot appear under `xvue/qt/` without breaking the build. Two enforcement points:
  - **CMake `verify_no_imagemagick_in_qt` ALL target** runs on every Qt build (DEPENDS xvueqt; ~50ms per build).
  - **Standalone `bin/test_no_imagemagick_in_qt.sh`** runnable outside CMake (CI hook, manual sanity, etc.). Auto-detects repo root when `$MEFISTO` is unset.
  - **Allowlist policy:** Qt API tokens (`QPageSize`, `convertToOther`, `convertTo(`, `QString::convertTo`, `convertFromUtf`, `Qt::ConvertibleTo`) are explicitly excluded so legitimate `convert` substring matches do not trip the gate. Word-boundary `\b` anchors keep the grep precise.

- **Byte-level golden test harness shipped (EXPORT-04 verification gate):**
  - **`xvue/qt/tests/golden/scene01_driver.f`** — deterministic Fortran scene driver. Emits 5 lines, 2 faces, 1 ellipse, 2 text strings, 1 epaisseur change, 1 typetrait change, 1 chargefonte call. The X11 backend's TEMPORAIRE.EPS output on this driver is the byte-level reference.
  - **Test 15 `PsEmitter_postscriptVerbatim_golden`** — drives the same primitive sequence via PsEmitter direct calls, captures TEMPORAIRE.EPS, byte-compares against `scene01.eps`. QSKIPs cleanly until the golden lands.
  - **Test 16 `PsEmitter_perPrimitive_golden`** — independent of the `.eps` golden; inlines the verbatim format-string substrings from the Plan 03 Format-String Parity Table. Always runs (no QSKIP).

- **GIF visual A/B test slots scaffolded:**
  - **`XvueExport_gif_AB_compare_wave` + `XvueExport_gif_AB_compare_cavity2d`** — locate the legacy baseline GIFs under `xvue/qt/tests/golden/`, validate frame count > 0 + readability via `QImageReader`. QSKIP cleanly until baselines committed in Task 3. The full frame-count + first/last-frame md5 diff against a Qt-emitted run is documented as a manual procedure in Plan 06 Task 3 step 9.

- **Phase 7 README section landed (~125 LOC):**
  - Surface table (PNG / JPEG / PDF / GIF / Capture Animation / TEMPORAIRE.EPS triggers + outputs + source files).
  - PsEmitter byte-for-byte contract + golden-test gate pointer.
  - Y-flip-inside-PsEmitter-only rule (README_COORDS.md mandate).
  - HiDPI export math (Pitfall 7: PDF uses LOGICAL canvas dims, raster formats use BACKING physical pixels).
  - Threading invariant (`XVUE_QT_ASSERT_MAIN_THREAD()` at every public method entry).
  - GIF probe-driven dispatch (D-10 native + D-11 ffmpeg fallback).
  - LVIDEO Phase 9 deferral (full pipeline + filenames + tracer-subroutine count).
  - Frame caps (T-07-03: 100 soft + 10000 hard).
  - Runtime dependencies (ffmpeg, qt6-imageformats-plugins).
  - Test target inventory + manual A/B verification log pointer.
  - `XVUE_ANIM=1` env var documented.

- **Non-regression gates green:**
  - **X11 full build** (`bin/cbl_tout`) — exit 0; 13 `pp/pp*` X11 binaries produced (T-07-08).
  - **Qt full build** (`cmake --build xvue/qt/build`) — exit 0; libxvueqt.a + 14 test binaries linked; all 4 verify_* gates pass.
  - **ABI symbol count** — 58 (unchanged; `verify_abi.sh` exit 0).
  - **xvue/xvuelc.c byte-identical** — `git diff HEAD~5 HEAD -- xvue/xvuelc.c` empty.
  - **LVIDEO untouched** — `git diff HEAD~5 HEAD -- xvue/video1.f videofin.f videonm.f` empty.
  - **Phase 7 unit tests green** — 18 PASS + 1 SKIP for postscript; 16 PASS + 2 SKIP for export. SKIPs flip to PASS-required automatically when Task 3 commits the goldens.
  - **Other Qt test binaries** — pass at pre-Plan-06 rate; one pre-existing failure (`testPerModuleGroupIsolation` in `xvue_qt_i18n_menu_prefs_tests`) documented as out-of-scope per SCOPE BOUNDARY rule (Plan 02 / 03 / 05 SUMMARYs).

## Task Commits

1. **Task 1 — `test(07-06): scaffold byte-level golden + GIF A/B test slots`** — `9cfd437`
   - `xvue/qt/tests/golden/scene01_driver.f` (117 LOC) — Fortran 77 fixed-form deterministic scene driver. Header BOOTSTRAP NOTE block documents the human procedure to materialize `scene01.eps`.
   - `xvue/qt/tests/test_xvue_qt_postscript.cpp` — +171 LOC. Tests 15 + 16. Test 15 walks up the directory tree from cwd to find `xvue/qt/tests/golden/scene01.eps` (or uses `$MEFISTO` if set), QSKIPs cleanly with bootstrap message when absent. Test 16 is independent of any external file.
   - `xvue/qt/tests/test_xvue_qt_export.cpp` — +93 LOC. `findGoldenPath` helper + two GIF A/B slots that QSKIP cleanly until baselines land.
   - All 18 + 16 = 34 unit tests still pass; 3 new SKIPs are explicit and pointer-to-procedure.

2. **Task 2 — `feat(07-06): ship EXPORT-06 grep gate + Phase 7 README section`** — `031357f`
   - `bin/test_no_imagemagick_in_qt.sh` (NEW, 76 LOC, executable). Allowlist regex + auto-repo-root detect.
   - `xvue/qt/CMakeLists.txt` — +18 LOC. `verify_no_imagemagick_in_qt ALL` target; `${CMAKE_CURRENT_SOURCE_DIR}/../../bin/...` resolves the repo root from `xvue/qt/`.
   - `xvue/qt/README.md` — +125 LOC. Phase 7 section.
   - `xvue/qt/tests/test_xvue_qt_export.cpp` — comment scrub (Rule 1/2: Plan 06 Task 1 inadvertently introduced two new hits; the gate added in this commit catches them; comment rephrased to remove the literal forbidden tokens without changing meaning).

_Plan metadata commit will be added by the orchestrator after the wave completes._

## Files Created/Modified

### Created

- `xvue/qt/tests/golden/scene01_driver.f` — Fortran 77 scene driver. Header documents the BOOTSTRAP procedure that materializes `xvue/qt/tests/golden/scene01.eps` when run by the human in Task 3.
- `xvue/qt/tests/golden/` directory — placeholder; Task 3 commits `scene01.eps`, `wave_legacy.gif`, `cavity2d_legacy.gif` here.
- `bin/test_no_imagemagick_in_qt.sh` (executable, mode 100755) — standalone EXPORT-06 grep gate.
- `.planning/phases/07-image-gif-and-postscript-export/VALIDATION-LOG.md` — manual A/B sign-off ledger; current state DEFERRED for the human-action steps; automated gates all PASS.
- `.planning/phases/07-image-gif-and-postscript-export/07-06-SUMMARY.md` — this file.

### Modified

- `xvue/qt/tests/test_xvue_qt_postscript.cpp` — Test 15 + Test 16 added at the end of the slot list. Test 15 walks up from `origCwd_` to find the golden; Test 16 inlines its expected substrings.
- `xvue/qt/tests/test_xvue_qt_export.cpp` — `findGoldenPath` helper + `XvueExport_gif_AB_compare_wave` + `_cavity2d` slots. Comment scrub for EXPORT-06 grep-gate compliance.
- `xvue/qt/CMakeLists.txt` — `verify_no_imagemagick_in_qt ALL` target added near `verify_icon_source`. Mirrors the existing `verify_*` guard pattern.
- `xvue/qt/README.md` — Phase 7 section appended before the existing References section. References section augmented with pointers to 07-CONTEXT.md and 07-RESEARCH.md.

## Decisions Made

### Defer-to-human-via-QSKIP (Plan ⊆ contract)

The byte-level golden compare against `scene01.eps` and the GIF A/B compare against `{wave,cavity2d}_legacy.gif` REQUIRE the human-executable bootstrap procedure (running scene01_driver.f against the X11 backend; running testa/wave + testa/cavity2d through both backends; eyeballing the GIFs). CLAUDE.md §Tests is unambiguous: large/long-running tests are user-run, and visual A/B is the canonical project verification model. Plan 06 Task 1 therefore ships the harness as QSKIPs that flip to PASS-required when the human commits the goldens — no silent skipping; every SKIP carries the bootstrap-procedure pointer inline.

### Headless-executor DEFERRED handling (autonomous: false adaptation)

The plan is `autonomous: false` because Task 3 requires a human verdict. The orchestrator-provided system reminder directs the executor to "document the running-without-A/B state in VALIDATION-LOG.md and SUMMARY.md so the orchestrator can route the sign-off as a verification gap during phase verification." Plan 06 honors this: VALIDATION-LOG.md has DEFERRED rows for each of the three human-action items (scene01.eps materialization, testa/wave A/B, testa/cavity2d A/B). The QSKIP test slots are the orthogonal in-codebase signal that the gate is incomplete. When a human-attached executor runs Plan 06 Task 3 in the future, they'll commit the goldens and the SKIPs flip to PASS automatically.

### EXPORT-06 grep-gate scope is xvue/qt/ ONLY (D-16/D-17)

Per CONTEXT.md the gate is precise about scope. The grep is `grep -rn ... xvue/qt/`, never `grep -rn ... xvue/`. This is the explicit T-07-06 mitigation: `xvue/video1.f`, `xvue/videofin.f`, `xvue/videonm.f`, `bin/convertepsgif` continue to legitimately match the forbidden tokens in the legacy tree until Phase 9 RETIRE-03. The README Phase 7 section explicitly documents this deferral so future agents do not re-discover the LVIDEO pipeline and over-delete.

### Allowlist for legitimate Qt API tokens

The grep matcher is `\b(convert|ImageMagick|magick)\b` (word-boundary-anchored). The negative-filter pipeline then strips: `QPageSize`, `convertToOther`, `convertTo(`, `QString::convertTo`, `convertFromUtf`, `Qt::ConvertibleTo`. These are Qt API names that legitimately contain "convert" as a substring. The two-step approach keeps the gate precise without weakening it against any genuine ImageMagick reference.

### Comment scrub in test_xvue_qt_export.cpp (Rule 1/2 deviation)

Task 1's first draft of the GIF A/B comment block mentioned `bin/convertepsgif` and `ImageMagick` literally — those substrings would have failed the grep gate added in Task 2. The comment was rephrased in Task 2's commit to convey equivalent meaning without the literal forbidden tokens (referring to "the legacy GIF post-processor" instead). Same logic applied to the new CMake comment in `xvue/qt/CMakeLists.txt` (rephrased "ImageMagick reference" → "legacy raster-tool reference").

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1/2 — Bug + Missing critical functionality] Plan 06 Task 1 documentation accidentally introduced 2 new EXPORT-06 grep-gate hits**

- **Found during:** Task 2 — running `bash bin/test_no_imagemagick_in_qt.sh` after creating the script.
- **Issue:** The new test_xvue_qt_export.cpp comment block in Task 1 mentioned `bin/convertepsgif` and `ImageMagick` literally. The grep gate added in Task 2 caught it (FAIL with "ImageMagick references found in xvue/qt/").
- **Fix:** Rephrased the comment block to avoid the literal forbidden tokens while keeping the explanatory intent ("the legacy GIF post-processor" instead of "bin/convertepsgif ImageMagick"). Same fix applied to the new comment in `xvue/qt/CMakeLists.txt`.
- **Files modified:** `xvue/qt/tests/test_xvue_qt_export.cpp`, `xvue/qt/CMakeLists.txt`.
- **Verification:** Re-running `bash bin/test_no_imagemagick_in_qt.sh` exits 0 with PASS.
- **Committed in:** `031357f` (Task 2 — bundled with the gate-introducing commit so the grep-gate-violating comment never lands as a separate commit).

### Task-level deferral (NOT auto-fixed)

**2. [DEFERRED — checkpoint:human-verify] Task 3 manual A/B sign-off**

- **Found during:** Task 3 entry — the autonomous executor has no attached human terminal.
- **Issue:** Task 3 requires (a) compiling scene01_driver.f against the X11 backend to produce `scene01.eps`, (b) running testa/wave through both X11 and Qt backends and eyeball-comparing the resulting GIFs, (c) same for testa/cavity2d, (d) recording PASS/FAIL/PARTIAL verdicts in VALIDATION-LOG.md, (e) committing the golden binaries, (f) re-running ctest to flip the QSKIPs to PASS.
- **Resolution:** Surfaced as a verification gap in VALIDATION-LOG.md (3 DEFERRED rows) per the orchestrator's system-reminder direction for headless-executor execution of `autonomous: false` plans. The plan's automated gates (EXPORT-06 grep gate, ABI count, X11 + Qt builds, all unit tests except the 3 SKIPs that flip when goldens land) all PASS today.
- **Verification:** `bash bin/test_no_imagemagick_in_qt.sh` PASS; `bin/cbl_tout` PASS; `cmake --build xvue/qt/build` PASS; `ctest -R '^xvue_qt_(postscript|export)_tests$'` 2/2 PASS (18+1 SKIP, 16+2 SKIP); `verify_abi.sh` 58 == 58.
- **Routing:** Phase verification step (orchestrator) will catch the DEFERRED rows in VALIDATION-LOG.md and the 3 SKIPs in ctest output, and route the sign-off as a verification gap.

---

**Total deviations:** 1 auto-fixed (Rule 1/2 grep-gate documentation hit) + 1 deferred (Task 3 manual A/B sign-off, properly surfaced).

**Impact on plan:** No architectural changes; no scope creep; the EXPORT-06 grep-gate contract is honored; the byte-parity test harness ships in a state that flips from QSKIP to PASS-required automatically when the human runs Task 3.

## Issues Encountered

- **Worktree base alignment:** the system-reminder worktree-branch-check pointed at `2a2a435` (the Plan 05 docs commit). The agent's branch was already at that base, so no rewind was needed. HEAD assertion + branch namespace check passed.
- **No `$MEFISTO` env var in the worktree session:** the `bin/cbl_tout` script and the new grep gate both need a repo root path. The grep gate auto-detects via `$(dirname $0)` walk-up; for `bin/cbl_tout` the agent set `MEFISTO=$(pwd)` inline. Build-time-regenerated artifacts (`incl/homdir.inc`, `*/lib`, `xvue/xvuelc.o`) were restored via `git checkout` after the build to keep the working tree clean for the SUMMARY commit.
- **Pre-existing `testPerModuleGroupIsolation` failure** in `xvue_qt_i18n_menu_prefs_tests` — same one Plan 02 / 03 / 05 SUMMARYs document. Verified by inspection (the test asserts QSettings group isolation behavior unrelated to Plan 06). Out of scope per SCOPE BOUNDARY rule.
- **Pre-existing `-Wdangling-reference` warnings** in `xvue/qt/src/xvue_qt_ther_actions.cpp` lines 191-193 — same warnings noted in Plan 02 / 03 / 05 SUMMARYs. Out of scope.

## TDD Gate Compliance

Plan 06 frontmatter is `type: execute` (not `type: tdd`), and Task 1 was tagged `tdd="true"`. Because the test slots designed-by-contract QSKIP when the (deferred) golden file is absent, the per-task TDD cycle collapsed into a single `test(...)` commit:

- **RED gate:** `9cfd437` — `test(07-06): scaffold byte-level golden + GIF A/B test slots`. Test 15 SKIPs cleanly because `scene01.eps` does not exist; Test 16 PASSES via the inline format-string substrings; the two GIF A/B slots SKIP cleanly because the baseline GIFs do not exist. The "RED" condition (failure-to-pass-without-golden) is encoded as the QSKIP message — once the golden lands, the slot flips to PASS-required, which is the GREEN gate the human runs in Task 3.
- **GREEN gate (deferred to Task 3 human checkpoint):** the human commits `scene01.eps` + `wave_legacy.gif` + `cavity2d_legacy.gif`; the SKIPs flip to PASS; the verdict lines in VALIDATION-LOG.md replace DEFERRED with PASS / FAIL / PARTIAL.

This split is documented in the test-slot QSKIP messages themselves so a future executor reading the SKIP banner sees the bootstrap-procedure pointer in-line.

## Self-Check: PASSED

**Files verified to exist:**
- `xvue/qt/tests/golden/scene01_driver.f` — FOUND (117 LOC, contains `! Phase 7 Plan 06`, `! BOOTSTRAP NOTE`)
- `xvue/qt/tests/golden/` directory — FOUND (created)
- `bin/test_no_imagemagick_in_qt.sh` — FOUND (76 LOC, executable mode 100755)
- `xvue/qt/CMakeLists.txt` — MODIFIED, contains `add_custom_target(verify_no_imagemagick_in_qt ALL` (1) and `DEPENDS xvueqt` (multiple — including the new target's instance)
- `xvue/qt/README.md` — MODIFIED, contains: `Phase 7` (5x), `LVIDEO` (4x), `Pitfall 7` (2x), `byte-for-byte` (3x), `XVUE_QT_ASSERT_MAIN_THREAD` (1x), `XVUE_ANIM=1` (2x)
- `xvue/qt/tests/test_xvue_qt_postscript.cpp` — MODIFIED, contains slots `PsEmitter_postscriptVerbatim_golden` (1) and `PsEmitter_perPrimitive_golden` (1)
- `xvue/qt/tests/test_xvue_qt_export.cpp` — MODIFIED, contains slots `XvueExport_gif_AB_compare_wave` (1) and `XvueExport_gif_AB_compare_cavity2d` (1)
- `.planning/phases/07-image-gif-and-postscript-export/VALIDATION-LOG.md` — FOUND (DEFERRED rows for testa/wave + testa/cavity2d + scene01.eps materialization)
- `.planning/phases/07-image-gif-and-postscript-export/07-06-SUMMARY.md` — FOUND (this file)

**Commits verified:**
- `9cfd437` (test 07-06 scaffold) — FOUND in git log
- `031357f` (feat 07-06 grep gate + README) — FOUND in git log

**Build gates verified:**
- `bash bin/test_no_imagemagick_in_qt.sh` — exit 0, "EXPORT-06 PASS: no ImageMagick references in xvue/qt/"
- `cmake --build xvue/qt/build -j4` — exit 0; verify_no_imagemagick_in_qt target green
- `xvue/qt/cmake/verify_abi.sh xvue/qt/build/libxvueqt.a xvue/qt/include/xvue_qt_api.h` — exit 0, "nm count: 58 header count: 58"
- `bin/cbl_tout` — exit 0; 13 X11 pp/pp* binaries produced
- `cd xvue/qt/build && xvfb-run --auto-servernum ctest -R '^xvue_qt_(postscript|export)_tests$' --output-on-failure` — 2/2 PASS
- Direct test run: `xvue_qt_postscript_tests` 18 PASS + 1 SKIP; `xvue_qt_export_tests` 16 PASS + 2 SKIP
- `git diff HEAD~5 HEAD -- xvue/xvuelc.c` — empty (T-07-08)
- `git diff --name-only HEAD~5 HEAD -- xvue/video1.f xvue/videofin.f xvue/videonm.f` — empty (T-07-06)

---

*Phase: 07-image-gif-and-postscript-export*
*Completed (autonomous portion): 2026-05-04*
*Task 3 (manual A/B): DEFERRED — see VALIDATION-LOG.md for elevation procedure.*
