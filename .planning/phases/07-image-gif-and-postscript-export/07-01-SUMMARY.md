---
phase: 07-image-gif-and-postscript-export
plan: 01
subsystem: infra
tags: [qt6, qimagewriter, cmake, probe, ffmpeg, gif, export]

requires:
  - phase: 06.5-shared-shell-menu-bridge-dialogs (and predecessors 00..06.5)
    provides: libxvueqt.a static library, xvue/qt/CMakeLists.txt scaffolding, ABI symbol count = 58, bin/cb*_qt shell-script convention
provides:
  - Empirically-validated PROBE.md recording QImageWriter::supportedImageFormats() for this host (gif_write_supported=0)
  - Reproducible probe binary build infrastructure (xvue/qt/probes/{CMakeLists.txt, qimagewriter_probe.cpp})
  - bin/cb_probe_qt — shell wrapper that builds + runs the probe and atomically rewrites PROBE.md
  - Locked GIF strategy decision for Plan 05 — D-11 ffmpeg fallback (D-10 native path unreachable)
affects: [07-02-postscript-emitter-state-machine, 07-03-postscript-per-primitive, 07-04-png-jpeg-pdf-export, 07-05-gif-ffmpeg-fallback, 07-06-validation-ab]

tech-stack:
  added:
    - QImageWriter (Qt 6 image-format probe API)
    - QImageIOHandler::Animation option flag
    - CMake option(XVUE_QT_BUILD_PROBES OFF) gating idiom
  patterns:
    - "Standalone probe binary outside libxvueqt.a (D-09): probes/ subdirectory with own CMakeLists.txt, gated by an OFF-by-default option so canonical builds (bin/cbl_tout, bin/cbl_tout_qt) are unaffected"
    - "T-07-01 single-redirect-target shell pattern: hard-coded > \"$PROBE_MD\" path with no argv handling — probe binary writes only to stdout"
    - "Empirical-first decision-locking: research predicts outcome, plan executes the probe, PROBE.md commits the result to git so downstream plans branch on a permanent record not maintainer memory"

key-files:
  created:
    - xvue/qt/probes/qimagewriter_probe.cpp
    - xvue/qt/probes/CMakeLists.txt
    - bin/cb_probe_qt
    - .planning/phases/07-image-gif-and-postscript-export/PROBE.md
  modified:
    - xvue/qt/CMakeLists.txt (added option(XVUE_QT_BUILD_PROBES OFF) + guarded add_subdirectory(probes))
    - .gitignore (ignore xvue/qt/build/ and xvue/qt/build-probe/ CMake build directories)

key-decisions:
  - "GIF strategy locked to D-11 (ffmpeg fallback). Empirical: Qt 6.10.2 gif_write_supported=0 — the bundled GIF plugin is read-only on Debian trixie."
  - "ffmpeg 8.1-3+b1 confirmed available on the build host — no missing-dependency blocker for Plan 05."
  - "Probe binary kept OFF in canonical builds via option(XVUE_QT_BUILD_PROBES OFF). T-07-08 mitigation: bin/cbl_tout exits 0 with all 12 pp/* executables, ABI symbol count unchanged at 58."

patterns-established:
  - "Standalone-probe pattern: any future Qt-runtime-capability check (font availability, plugin presence, GL feature) follows the same probes/ subdir + option-gated CMake target + bin/cb_*_qt shell wrapper recipe"
  - "Empirical-evidence-as-frozen-decision pattern: research-flagged uncertainties get resolved by a committed empirical artifact (PROBE.md here) before downstream plans execute"

requirements-completed: [EXPORT-01]

duration: 11m20s
completed: 2026-05-04
---

# Phase 7 Plan 01: QImageWriter probe + GIF strategy lock Summary

**Empirical PROBE.md committed: gif_write_supported=0 on Qt 6.10.2 / Debian trixie — Plan 05 GIF path locked to D-11 ffmpeg fallback (8.1-3+b1 available)**

## Performance

- **Duration:** 11m20s (680 sec wall — including a 367-sec `bin/cbl_tout` non-regression build)
- **Started:** 2026-05-04T16:32:20Z
- **Completed:** 2026-05-04T16:43:41Z
- **Tasks:** 2 / 2
- **Files modified:** 6 (4 created, 2 modified)

## Accomplishments

- **EXPORT-01 deliverable shipped:** `.planning/phases/07-image-gif-and-postscript-export/PROBE.md` is committed at `ed1349c` with the empirical Qt 6.10.2 / Debian trixie probe output. Plans 02-06 can reference it via `@.planning/phases/07-image-gif-and-postscript-export/PROBE.md`.
- **Probe build infrastructure reproducible:** running `bin/cb_probe_qt` (with `$MEFISTO` exported) regenerates PROBE.md from scratch in ~10 seconds. Future maintainers can re-validate after Qt or imageformats-plugin upgrades.
- **Canonical builds untouched (T-07-08 verified):** `bin/cbl_tout` exits 0; all 12 X11 executables (`ppelas`, `ppflui`, `ppinit`, `ppmail`, `ppnlse`, `pppoba`, `ppther`, `ppxvtest0..4`, `pxyz`) build cleanly. Default Qt build (`cmake -S xvue/qt -B …`) compiles `libxvueqt.a` with no probe dependency. `verify_abi.sh` confirms ABI symbol count remains at 58.
- **GIF strategy decision locked:** `gif_write_supported=0` confirms research prediction. Plan 05 will dispatch via `QProcess::execute("ffmpeg", …)` (D-11) — no QImageWriter multi-frame path needed.

## Task Commits

Each task was committed atomically:

1. **Task 1: Probe binary + probes CMakeLists + xvue/qt/CMakeLists.txt subdirectory hook** — `2c30d9c` (feat)
2. **Task 2: bin/cb_probe_qt + run + capture PROBE.md + commit** — `ed1349c` (feat)

_Plan metadata commit will be added by the orchestrator after the wave completes._

## Files Created/Modified

- `xvue/qt/probes/qimagewriter_probe.cpp` (NEW, 25 LOC) — standalone Qt 6 probe; prints `qt_version`, `supported_write_formats`, `gif_write_supported`, `gif_animation_supported` to stdout. Pure C++17, links Qt6::Core + Qt6::Gui only.
- `xvue/qt/probes/CMakeLists.txt` (NEW) — `add_executable(qimagewriter_probe qimagewriter_probe.cpp)` with -Wall -Wextra. NOT linked into libxvueqt.a (D-09).
- `xvue/qt/CMakeLists.txt` (MODIFIED) — added `option(XVUE_QT_BUILD_PROBES "..." OFF)` guard + `if(XVUE_QT_BUILD_PROBES) add_subdirectory(probes) endif()`. Inserted between `verify_icon_source` target and `install(TARGETS xvueqt …)` to keep the file's logical flow (verify-targets, then probe-gate, then install).
- `bin/cb_probe_qt` (NEW, +x, 70 lines) — shell wrapper. Defensive `command -v ffmpeg` check (Pitfall 6). Build dir hard-coded to `$MEFISTO/xvue/qt/build-probe`. Single redirect target `> "$PROBE_MD"` (T-07-01 mitigation).
- `.planning/phases/07-image-gif-and-postscript-export/PROBE.md` (NEW, committed) — empirical record. Plans 02-06 kickoff gate.
- `.gitignore` (MODIFIED) — ignore `xvue/qt/build/` and `xvue/qt/build-probe/` CMake artifact directories.

## Empirical Probe Output

```
qt_version=6.10.2
supported_write_formats=bmp cur icns ico jfif jpeg jpg pbm pgm png ppm tif tiff wbmp webp xbm xpm
gif_write_supported=0
gif_animation_supported=0
```

ffmpeg version: `ffmpeg version 8.1-3+b1 Copyright (c) 2000-2026 the FFmpeg developers`

## Decisions Made

- **GIF strategy → D-11 ffmpeg fallback (locked).** Empirical confirmation that Qt 6.10.2's bundled GIF plugin on Debian trixie is read-only. Plan 05 will use PNG-sequence + `QProcess::execute("ffmpeg", …)` synchronous dispatch (no GUI-thread async needed for testa-scale frame counts).
- **ffmpeg version policy.** ffmpeg 8.1 is verified working — no further version pinning required at this stage. Plan 05's `QProcess::execute` invocation will pass `-y -framerate -i frame_%04d.png` arguments, which are stable across ffmpeg 4.x → 8.x.
- **Probe stays out of libxvueqt.a.** D-09 reaffirmed: standalone executable in `probes/`, built only when `-DXVUE_QT_BUILD_PROBES=ON` is explicitly passed. Default builds (`bin/cbl_tout`, `bin/cbl_tout_qt`, raw `cmake -S xvue/qt -B …`) do NOT touch the probes subdirectory — verified.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing Critical] Added .gitignore entries for CMake build directories**
- **Found during:** Task 2 (after `bin/cb_probe_qt` created `xvue/qt/build-probe/`)
- **Issue:** The new shell script creates a CMake artifact directory `xvue/qt/build-probe/` that the existing `.gitignore` did not exclude. Without an ignore rule the directory shows up as untracked on every probe re-run, polluting `git status` and risking accidental staging of build artifacts.
- **Fix:** Appended a "Qt CMake build directories" section to `.gitignore` that ignores both `xvue/qt/build/` (used by other plans / `bin/cbl_tout_qt`) and `xvue/qt/build-probe/` (this plan's probe build dir).
- **Files modified:** `.gitignore`
- **Verification:** `git status --short` after `bin/cb_probe_qt` no longer lists `xvue/qt/build-probe/` as untracked.
- **Committed in:** `ed1349c` (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 missing critical hygiene)
**Impact on plan:** Hygiene-level — no scope creep. The original plan referenced `xvue/qt/build-probe` explicitly via `bin/cb_probe_qt` but did not call out the gitignore consequence; the fix is a 4-line additive change keeping git status clean for downstream agents.

## Issues Encountered

- **Worktree base mismatch at agent start.** The worktree branch was checked out at `ac282f8` (the historical Initial-Commit) instead of `2d9e272` (the current Phase 7 planning HEAD). Worktree was sparse — no `.planning/`, no `xvue/qt/`. Resolved with `git reset --soft 2d9e272` followed by `git checkout HEAD -- .` to populate the working tree from the new HEAD. No commits were lost (no commits had been made yet). Time cost: ~2 min. Captured in the worktree-branch-check protocol — no further action needed by orchestrator.
- **Pre-existing `-Wdangling-reference` warnings in `xvue/qt/src/xvue_qt_ther_actions.cpp`** lines 191-193 surfaced when building `libxvueqt.a` to verify the default-config path. These warnings predate Plan 07-01 and are out of scope per the SCOPE BOUNDARY rule. Logged here for awareness; not blocking. (`XvueMenuFileParser::loadFor()` returns by value into a `const MenuFile&` — three sites in the ther-actions module. Phase 6.4 territory.)

## Next Phase Readiness

- **Plan 02 (PostScript state machine) — UNBLOCKED.** No dependency on PROBE.md content; can start at Wave 2 entry.
- **Plan 03 (PostScript per-primitive) — UNBLOCKED.** Same as Plan 02.
- **Plan 04 (PNG/JPEG/PDF export) — UNBLOCKED.** PROBE.md confirms `png` and `jpeg` are in `supported_write_formats`; no surprises.
- **Plan 05 (GIF ffmpeg fallback) — STRATEGY LOCKED.** PROBE.md commits `gif_write_supported=0`; Plan 05 takes the D-11 path. ffmpeg 8.1 verified available.
- **Plan 06 (validation A/B) — UNBLOCKED for the parts that don't depend on Plan 05's output.**

No blockers, no concerns.

## Self-Check: PASSED

**Files verified to exist:**
- `xvue/qt/probes/qimagewriter_probe.cpp` — FOUND
- `xvue/qt/probes/CMakeLists.txt` — FOUND
- `bin/cb_probe_qt` — FOUND, executable
- `.planning/phases/07-image-gif-and-postscript-export/PROBE.md` — FOUND, non-empty (with `supported_write_formats=`, `gif_write_supported=`, `gif_animation_supported=`, `qt_version=`, D-10/D-11 references, ffmpeg version)
- `xvue/qt/CMakeLists.txt` — FOUND, contains `option(XVUE_QT_BUILD_PROBES` AND `add_subdirectory(probes)` AND `add_library(xvueqt STATIC` (untouched library target)

**Commits verified:**
- `2c30d9c` (Task 1) — FOUND in git log
- `ed1349c` (Task 2) — FOUND in git log

**Build gates verified:**
- `cmake -S xvue/qt -B /tmp/xvueqt-probe-build -DXVUE_QT_BUILD_PROBES=ON` — exit 0
- `cmake --build /tmp/xvueqt-probe-build --target qimagewriter_probe` — exit 0
- `cmake -S xvue/qt -B /tmp/xvueqt-default && cmake --build /tmp/xvueqt-default --target xvueqt` — exit 0 (probe NOT a dependency)
- `bin/cbl_tout` (X11 + all Fortran libraries) — exit 0, 12 pp/* executables produced (T-07-08 non-regression PASSED)
- `bash xvue/qt/cmake/verify_abi.sh /tmp/xvueqt-default/libxvueqt.a xvue/qt/include/xvue_qt_api.h` — exit 0, "nm count: 58 header count: 58" (ABI unchanged)

---
*Phase: 07-image-gif-and-postscript-export*
*Completed: 2026-05-04*
