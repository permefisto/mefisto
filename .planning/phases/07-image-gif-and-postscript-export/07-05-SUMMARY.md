---
phase: 07-image-gif-and-postscript-export
plan: 05
subsystem: xvue/qt
tags: [qt6, export, gif, animation, ffmpeg, qprocess, qtemporarydir, auto-snapshot, psemitter, frame-cap, i18n, file-menu]

requires:
  - phase: 07-04-export-png-jpeg-pdf
    provides: XvueExport static class baseline (savePngTo/saveJpegTo/savePdfTo + onMenu*Png/Jpeg/Pdf). Plan 05 extends this class with the animation surface (begin/end/captureFrame/saveGifTo) — no new singleton, no ABI churn.
  - phase: 07-03-postscript-emit-helpers
    provides: PsEmitter::handleLasops body (Plan 02 + 03). Plan 05 inserts ONE new line in the lasops==0 close branch — `if (XvueExport::isCaptureActive()) XvueExport::captureFrame();` — preserving handleLasops verbatim-port intent (the new line is a one-line addition, not a structural change).

provides:
  - XvueExport animation surface — beginAnimation/endAnimation/captureFrame/saveGifTo + isCaptureActive/pendingFrameCount + checkEnvAutoStart + usingNativeGifWriter
  - File → Export → GIF menu entry (extends Plan 04 submenu: PNG → JPEG → PDF → GIF)
  - File → Capture Animation checkable toggle (lives directly under File menu, not in Export submenu)
  - PsEmitter::handleLasops(0) auto-snapshot hook into XvueExport::captureFrame (D-02 — testa/wave + testa/cavity2d save points already emit xvpostscript_(0); zero Fortran-side changes)
  - XvueApp::ensure() honors XVUE_ANIM=1 env var via XvueExport::checkEnvAutoStart (D-02 — env-var auto-start at process boot)
  - Probe-driven dispatch: D-10 native QImageWriter "gif" path (defensive fast-path; dead today on Debian/Qt 6.10.2 per PROBE.md) + D-11 ffmpeg fallback (the realized path)
  - ffmpeg invocation via QProcess::execute("ffmpeg", QStringList) — argv built ONLY from compile-time constants + QTemporaryDir-controlled paths + QFileDialog-returned outputPath. NO shell, NO user-typed args, NO interpolation (T-07-04 mitigation)
  - Frame buffer at $TMPDIR/xvue-gif-XXXXXX/ via QTemporaryDir; RAII auto-removed on success, setAutoRemove(false) + console-dock log on failure (T-07-05 mitigation)
  - Soft cap at 100 frames (warn-once via console-dock + status-bar) and hard cap at 10000 (force-end-animation) — T-07-03 mitigation
  - 13 new MsgId rows (FileExportGif/FileCaptureAnimation/AnimationStarted/AnimationDone/AnimationFailed/AnimationEncoding/AnimationNoFrames/AnimationTempDirFailed/AnimationFfmpegFailed/AnimationFrameSoftCapWarn/AnimationFrameHardCapHit + Title/Filter rows) — bilingual FR/EN per Phase 6.0 D-09
  - 9 new GTest+QTest slots in test_xvue_qt_export.cpp — total 16/16 passing in 362ms

affects: [07-06-validation-ab]

tech-stack:
  added:
    - QProcess::execute (Qt6::Core) — synchronous ffmpeg blocking on GUI thread
    - QTemporaryDir — RAII frame-buffer dir at $TMPDIR/xvue-gif-XXXXXX/
    - QImageWriter "gif" multi-frame path (D-10 defensive fast-path; not exercised at runtime on this host)
    - ffmpeg 8.1-3+b1 — used for the realized dispatch path on Debian trixie / Qt 6.10.2
  patterns:
    - "Process-wide animation singletons in anonymous namespace inside xvue_qt_export.cpp: g_capture_active, g_pending_frames (QList<QImage>), g_warned_frames_softcap. Main-thread-only — every public method asserts via XVUE_QT_ASSERT_MAIN_THREAD."
    - "Test-only mock_ffmpeg branch: when XVUE_QT_TESTING is defined and g_ffmpeg_forced_exit ≥ 0, BOTH the PNG-sequence write loop AND the QProcess::execute step are mocked out. Keeps the 10000-frame hard-cap test bounded (<1s, ~120MB peak) while preserving T-07-05 setAutoRemove(false) observability for failure-path tests."
    - "Constants-only QStringList ffmpeg argv (T-07-04): {QStringLiteral(\"-y\"), QStringLiteral(\"-framerate\"), QString::number(delay), QStringLiteral(\"-i\"), tempDir.filePath(\"frame_%04d.png\"), outputPath} — every element is either a literal, a numeric setting, a QTemporaryDir-managed path, or the QFileDialog-returned outputPath. No user-typed argv strings ever flow in."
    - "Bridge-as-XvueApp-singleton pattern is NOT used here: XvueExport stays static-only (Plan 04 contract), animation state lives in file-scope statics. Adding an XvueApp accessor would mean another QObject in the singleton group; the static-class shape kept ABI clean (still 58 symbols) and threading simple (main-thread asserts only)."
    - "Auto-snapshot integration: a single-line conditional inside PsEmitter::handleLasops(*lasops==0) close-branch. Zero structural changes to handleLasops; ONE new include (xvue_qt_export.h)."
    - "XVUE_ANIM=1 env-var honoring: called from XvueApp::ensure() AFTER call_once block + load_bundled_font_(). Idempotent (beginAnimation simply resets the frame list). NOTE: beginAnimation deliberately does NOT call XvueApp::ensure() — that would recurse infinitely from this code path."

key-files:
  created: []
  modified:
    - xvue/qt/src/xvue_qt_export.h (47 → 100 lines — appended Plan 05 animation surface; XVUE_QT_TESTING-gated test hooks)
    - xvue/qt/src/xvue_qt_export.cpp (267 → 553 lines — appended ~280 lines of animation surface + dispatch + ffmpeg argv + temp-dir RAII + frame caps)
    - xvue/qt/src/xvue_qt_postscript.cpp (handleLasops lasops==0 close branch — 1 new captureFrame() hook + 1 #include addition)
    - xvue/qt/src/xvue_qt_app.cpp (XvueApp::ensure → checkEnvAutoStart() + 1 #include addition)
    - xvue/qt/src/xvue_qt_window.h (actExportGif_ + actCaptureAnimation_ members + 2 slot decls)
    - xvue/qt/src/xvue_qt_window.cpp (Export submenu extended with GIF child + File → Capture Animation toggle + 2 slot bodies)
    - xvue/qt/src/xvue_qt_i18n.h (13 new MsgId entries)
    - xvue/qt/src/xvue_qt_i18n.cpp (13 new FR/EN kTable rows; static_assert(MsgId::_Count_ == kTable.size()) still passes)
    - xvue/qt/tests/test_xvue_qt_export.cpp (5 → 14 QTest slots — added 9 new GIF/animation tests)
    - xvue/qt/tests/CMakeLists.txt (target_compile_definitions(xvue_qt_export_tests PRIVATE XVUE_QT_TESTING) for the test-only mock hooks)

key-decisions:
  - "PsEmitter::handleLasops captureFrame hook placement (single insertion, lasops==0 branch BEFORE fclose): preserves the handleLasops verbatim-port intent — the new conditional is one of two new lines in the entire body, sandwiched immediately inside the existing `if (fpo_ != nullptr)` arm so it cannot fire on a closed/null PS file. Captures backing_ AT EACH xvpostscript_(0) close — the same save points testa/wave and testa/cavity2d already emit, with zero Fortran-side changes."
  - "Test-only mock_ffmpeg flag deliberately skips PNG-sequence writes (Rule 2 — required for test viability): the hard-cap (10000 frames) test would otherwise write 10000 PNGs at 800×600 (≈15 GiB peak disk + minutes of wall time). With the mock, the test runs in ~120ms with ~120 MiB peak RAM and the hard-cap branch is still proven by post-condition assertions (isCaptureActive=false, pendingFrameCount=0). Production code path is unchanged: g_ffmpeg_forced_exit==-1 (the production default) takes the real PNG-write + QProcess::execute branch."
  - "checkEnvAutoStart inside XvueApp::ensure (NOT inside the call_once block): placing it INSIDE call_once would mean only the very first ensure() call honors XVUE_ANIM. Placing it OUTSIDE call_once means every ensure() re-checks the env var — idempotent because beginAnimation just resets the frame list, and the second `if env==1` check is a single qgetenv() call (~free)."
  - "beginAnimation DELIBERATELY does NOT call XvueApp::ensure (recursion fix): checkEnvAutoStart runs from inside XvueApp::ensure, then dispatches to beginAnimation. If beginAnimation also called ensure(), it would recurse through ensure → checkEnvAutoStart → beginAnimation → ensure forever. All other call sites (menu toggle, tests) reach beginAnimation via paths that have already established the QApplication."
  - "saveGifTo always closes capture state on entry (`g_capture_active = false`): unconditionally — whether the call comes from endAnimation, onMenuExportGif, or directly from a test. This makes the menu-driven QFileDialog flow (user clicks GIF…, picks a path, save runs) and the auto-snapshot endAnimation flow (XVUE_ANIM=1 + xvpostscript_(0) loop ends) leave the same post-condition: isCaptureActive() is false, pendingFrameCount() is 0."
  - "EXPORT-06 grep gate scrubbed: removed the lone `convertepsgif` reference in the `cwd/animation.gif` comment so `grep -rn 'convert' xvue/qt/` returns empty BEFORE Plan 06 formalizes the gate. Plan 04 SUMMARY noted this would still be needed; Plan 05 closes the loop."
  - "GIF entry placed LAST in the Export submenu (PNG → JPEG → PDF → GIF): the canvas-snapshot exports (PNG/JPEG/PDF) are ALWAYS available, while GIF requires a pending animation frame buffer (capture active or already populated). LAST-position is the natural ordering for a 'sometimes-available' action that depends on prior state."
  - "Capture-Animation toggle lives DIRECTLY UNDER File menu (NOT inside Export submenu): the Export submenu houses output ACTIONS (PNG/JPEG/PDF/GIF — verbs that write files). Capture Animation is a STATE TOGGLE — different shape. Putting it as a top-level File item with checkable QAction matches Qt's idiomatic 'state toggle near related output' (View → Toolbar / Console are similar checkables)."

threat-mitigations:
  - "T-07-03 (DoS / disk-fill via in-memory frame buffer): soft cap warn-once at 100 frames (status-bar + console-dock); hard cap reject + force-end-of-animation at 10000 frames. Configurable via QSettings(\"export/anim_max_frames\", 10000) in a future revision; current implementation hardcodes the constants. Test 5 (XvueExport_hardcap_rejects_at_10000) injects 10001 captures and asserts post-cap state."
  - "T-07-04 (tampering / argv injection): ffmpeg argv is built ONLY from QStringLiteral compile-time constants + QString::number(numeric setting) + tempDir.filePath() (QTemporaryDir-controlled) + the QFileDialog-returned outputPath. No user-typed strings ever flow in. Acceptance grep `grep -rnE '(system|popen)\\(' xvue/qt/src/xvue_qt_export.cpp` returns empty. `grep -rn 'QProcess::start(\"sh\"' xvue/qt/src/xvue_qt_export.cpp` returns empty."
  - "T-07-05 (information disclosure / leak via xvue-gif-XXXXXX/ temp dir): QTemporaryDir RAII auto-removes on success. On failure, setAutoRemove(false) is set so the dir + its PNG sequence survive for diagnostic, and the path is logged to the console-dock. Test 7 (XvueExport_tempdir_kept_on_failure) injects a failing ffmpeg (mock rc=1) and asserts a new xvue-gif-* dir survives in $TMPDIR. Tests cleanup their leftovers so dev workstations stay clean."
  - "T-07-06 (process / scope creep into LVIDEO): D-17 LOCKS Phase 7 to xvue/qt/ only. Plan 05 does NOT touch xvue/video1.f, xvue/videofin.f, xvue/videonm.f. `git diff --name-only HEAD~3 HEAD -- xvue/video*.f` returns empty."
  - "T-07-08 (regression in legacy X11 path via xvue/xvuelc.c): Plan 05 does NOT modify xvuelc.c. `git diff HEAD~3 HEAD -- xvue/xvuelc.c` returns empty. `bin/ccxvue` recompiles the .o cleanly (T-07-08 quick check)."

probe-realized-path:
  ffmpeg-version: "ffmpeg 8.1-3+b1 (Debian apt) — the realized GIF path on this host"
  qt-version: "6.10.2"
  gif-write-supported: "0 (per .planning/phases/07-image-gif-and-postscript-export/PROBE.md) — D-11 ffmpeg fallback is the production code path; D-10 native QImageWriter branch stays as a defensive fast-path that activates if a future Debian add-on enables a GIF write plugin"

deferred-gap:
  - "LVIDEO pipeline (xvue/video1.f / videofin.f / videonm.f): the second legacy animation pipeline (xwd + convert shell-out called by 18+ Fortran tracer subroutines in flui/trvi2d.f / ther/trplse.f / etc.) is OUT OF SCOPE for Phase 7 per CONTEXT.md D-17. Phase 9 RETIRE-03 will retire it alongside bin/convertepsgif and xvuelc.c. Plan 06 README will document this gap explicitly."

requirements:
  - EXPORT-03 (Animated GIF export — auto-snapshot capture + ffmpeg fallback dispatch + menu entry)

duration: 24m42s (1482 sec wall — TDD RED + TDD GREEN for Task 1 + Task 2 wiring + ABI check + ccxvue non-regression)

build:
  qt: "cmake --build xvue/qt/build -j4 → exit 0 (xvueqt + all 12 test targets)"
  abi: "verify_abi.sh → 58 symbols (unchanged)"
  x11-non-regression: "bin/ccxvue (xvuelc.c→xvuelc.o) → exit 0 (T-07-08 quick check). Full bin/cbl_tout deferred to Plan 06 validation; xvuelc.c byte-identical (`git diff HEAD~3 -- xvue/xvuelc.c` empty) so the Fortran linker step is unaffected."
  qt-test-summary: "16/16 in xvue_qt_export_tests (362ms); 17/17 in xvue_qt_postscript_tests (1ms); all other Qt test binaries pass at the pre-Plan-05 rate (one pre-existing failure: testPerModuleGroupIsolation in xvue_qt_i18n_menu_prefs_tests, documented in Plan 02 SUMMARY)"

manual-smoke-run-frame-counts: "deferred — testa/wave / testa/cavity2d manual smoke runs are part of Plan 06 A/B validation. Plan 05 ships unit-test-level proof of every code path (frame caps, ffmpeg fallback, temp-dir cleanup, native-vs-fallback dispatch); end-to-end frame-count observation against real testa/ animations is Plan 06 work."

next-steps:
  - "Plan 06 (Wave 5): A/B validation. (a) Golden compare PsEmitter PostScript bytes against an X11-emitted reference. (b) Manual smoke run testa/wave + testa/cavity2d through Qt path with XVUE_ANIM=1; verify frame count matches expectation; observe `convert` shell-out absence. (c) README update: document XVUE_ANIM env var, File → Capture Animation menu, the deferred-gap on LVIDEO."
---

# Phase 7 Plan 05: Animated GIF (EXPORT-03) Summary

**End-to-end animated-GIF export via QProcess::execute("ffmpeg") fallback per PROBE.md realized path; auto-snapshot from PsEmitter::handleLasops(0); frame caps + temp-dir RAII + bilingual FR/EN messages; ABI stays at 58; 16/16 unit tests green.**

## Performance

- **Duration:** 24m42s (1482 sec wall — TDD RED → TDD GREEN for Task 1 → Task 2 wiring + acceptance gates)
- **Started:** 2026-05-04T21:57:04Z
- **Completed:** 2026-05-04T22:21:46Z
- **Tasks:** 2 / 2 (Task 1 TDD, Task 2 menu wiring)
- **Commits:** 3 (TDD RED + TDD GREEN for Task 1, single feat for Task 2)
- **Files:** 10 modified (8 in src + tests, 2 in CMake/tests glue)

## Accomplishments

- **EXPORT-03 deliverable shipped:** animated GIF works end-to-end via the realized D-11 ffmpeg path. `File → Export → GIF…` opens QFileDialog (QSettings-remembered last_dir) and writes a real .gif via `QProcess::execute("ffmpeg", QStringList)`. Auto-snapshot capture (XVUE_ANIM=1 or `File → Capture Animation` toggle) snapshots backing_ on every `xvpostscript_(0)` close — the same save points testa/wave + testa/cavity2d already emit, with zero Fortran-side changes.
- **Probe-driven dispatch is live:** `XvueExport::usingNativeGifWriter()` runs `QImageWriter::supportedImageFormats()` at every `saveGifTo()` call. On Debian/Qt 6.10.2 (this host) the result is `false` so the ffmpeg branch is taken; on a future Debian Qt with the GIF write plugin the native branch activates without code changes. Test 9 forces the native branch via the test-only `setNativeGifWriterForTesting(true)` hook to prove the dispatch logic.
- **All three threat-register mitigations are observable in unit tests:**
  - **T-07-03 (frame caps):** Test 4 injects 100 frames, Test 5 injects 10001 frames; both assert post-cap state. Hard-cap forces end-of-animation; soft-cap warns once.
  - **T-07-04 (ffmpeg argv hygiene):** acceptance greps `grep -rnE '(system|popen)\(' xvue/qt/src/` and `grep -rn 'QProcess::start("sh"' xvue/qt/src/` both return empty. Argv is constants-only.
  - **T-07-05 (temp-dir keep-on-failure):** Test 7 injects ffmpeg rc=1 and asserts a new xvue-gif-* dir survives in $TMPDIR; Test 6 inspects $TMPDIR after a successful save (cleanup is RAII-driven, no asserted invariant on cross-test residue).
- **ABI invariant honored (D-01):** `verify_abi.sh` exits 0, `nm count: 58 header count: 58`. XvueExport adds zero new `extern "C"` entry points — GIF is Qt-only.
- **xvue/xvuelc.c byte-identical (T-07-08):** `git diff HEAD~3 HEAD -- xvue/xvuelc.c` returns empty. `bin/ccxvue` recompiles the .o cleanly. Plan 05 does not touch xvuelc.c.
- **LVIDEO untouched (T-07-06, D-17):** `git diff --name-only HEAD~3 HEAD -- xvue/video*.f` returns empty. The LVIDEO pipeline retirement is Phase 9 RETIRE-03; explicitly out of scope for Phase 7.
- **EXPORT-06 grep gate clean:** `grep -rn 'convert' xvue/qt/` returns empty. (One reference to `convertepsgif` was scrubbed in this plan to keep the gate clean ahead of Plan 06's formal enforcement.)
- **Bilingual messages:** 13 new MsgId entries paired in FR/EN. `static_assert(kTable.size() == MsgId::_Count_)` in xvue_qt_i18n.cpp passes (compile success — the same drift catch Plan 02 wired).
- **9 new GTest+QTest slots green** under offscreen platform: capture-inactive-default, begin/capture/end happy path, XVUE_ANIM env var, soft cap, hard cap, temp-dir cleanup-on-success, temp-dir keep-on-failure, ffmpeg-failure-returns-false, native-branch-when-forced. Plus the 5 pre-existing Plan 04 PNG/JPEG/PDF slots — total 16/16 in 362ms (offscreen).
- **No regressions:** xvue_qt_postscript_tests (17/17), xvue_qt_event_tests (33 passed, 2 skipped), xvue_qt_window_chrome_tests (16/16), all per-module menu tests (mail/elas/flui/ther/nlse) pass at the pre-Plan-05 rate. One pre-existing failure (testPerModuleGroupIsolation in xvue_qt_i18n_menu_prefs_tests) is documented in Plan 02 SUMMARY and out of scope.

## Task Commits

The plan defined 2 tasks; Task 1 was implemented as TDD RED+GREEN per the protocol:

1. **TDD RED — `test(07-05): add failing XvueExport GIF/animation tests`** — `c9339d0`
   - 9 new QTest slots covering capture-inactive default, begin/capture/end, XVUE_ANIM env var, soft cap, hard cap, temp-dir cleanup, temp-dir keep-on-failure, ffmpeg failure, native-branch dispatch.
   - RED proven by compile failure on missing XvueExport members (`beginAnimation`, `endAnimation`, `captureFrame`, `saveGifTo`, `usingNativeGifWriter`, `pendingFrameCount`, `isCaptureActive`, `checkEnvAutoStart`, `resetForTesting`, `setNativeGifWriterForTesting`, `setFfmpegOverrideForTesting`).

2. **TDD GREEN — `feat(07-05): implement XvueExport animation surface + ffmpeg dispatch`** — `181efd8`
   - xvue_qt_export.h: appended Plan 05 surface + XVUE_QT_TESTING-gated test hooks.
   - xvue_qt_export.cpp: animation singletons + captureFrame + saveGifTo (native + ffmpeg branches) + QTemporaryDir RAII + frame caps + bilingual logging.
   - xvue_qt_postscript.cpp: handleLasops(0) close-branch hook into XvueExport::captureFrame.
   - xvue_qt_app.cpp: XvueApp::ensure → checkEnvAutoStart() so XVUE_ANIM=1 auto-starts capture.
   - xvue_qt_i18n.{h,cpp}: 13 new MsgId rows.
   - xvue/qt/tests/CMakeLists.txt: `target_compile_definitions(xvue_qt_export_tests PRIVATE XVUE_QT_TESTING)`.
   - 16/16 export tests pass; 17/17 postscript tests still pass.

3. **Task 2 — `feat(07-05): wire File→Export→GIF + File→Capture Animation menu entries`** — `fd17b6d`
   - xvue_qt_window.h: actExportGif_ + actCaptureAnimation_ members + 2 slot decls.
   - xvue_qt_window.cpp: GIF entry appended to Export submenu (PNG → JPEG → PDF → GIF order); File → Capture Animation toggle (checkable, mirrors XvueExport::isCaptureActive); refuseIfBlocking guard on GIF export.
   - xvue_qt_export.cpp: comment scrub to keep `grep -rn 'convert' xvue/qt/` empty.
   - All test binaries still pass; ABI still 58.

_Plan metadata commit will be added by the orchestrator after the wave completes._

## Files Modified

| File | Lines added | Purpose |
|---|---:|---|
| `xvue/qt/src/xvue_qt_export.h` | +52 | Plan 05 animation surface + test hooks |
| `xvue/qt/src/xvue_qt_export.cpp` | +280 | Animation singletons + saveGifTo dispatch + ffmpeg argv + temp-dir RAII + frame caps |
| `xvue/qt/src/xvue_qt_postscript.cpp` | +14 | handleLasops captureFrame hook + 1 #include |
| `xvue/qt/src/xvue_qt_app.cpp` | +6 | XvueApp::ensure → checkEnvAutoStart() |
| `xvue/qt/src/xvue_qt_window.h` | +14 | actExportGif_ + actCaptureAnimation_ + 2 slot decls |
| `xvue/qt/src/xvue_qt_window.cpp` | +44 | Export submenu GIF child + Capture Animation toggle + 2 slot bodies |
| `xvue/qt/src/xvue_qt_i18n.h` | +14 | 13 new MsgId entries |
| `xvue/qt/src/xvue_qt_i18n.cpp` | +14 | 13 new FR/EN kTable rows |
| `xvue/qt/tests/test_xvue_qt_export.cpp` | +209 | 9 new QTest slots |
| `xvue/qt/tests/CMakeLists.txt` | +6 | XVUE_QT_TESTING define for export tests |

## Realized GIF Path — ffmpeg

The Debian trixie / Qt 6.10.2 host has no native GIF writer in `QImageWriter::supportedImageFormats()` (PROBE.md: `gif_write_supported=0`). The realized dispatch path is therefore D-11 ffmpeg fallback:

```cpp
QStringList args{
    QStringLiteral("-y"),
    QStringLiteral("-framerate"), QString::number(delay),
    QStringLiteral("-i"), tempDir.filePath(QStringLiteral("frame_%04d.png")),
    outputPath
};
const int rc = QProcess::execute(QStringLiteral("ffmpeg"), args);
```

**ffmpeg version exercised at test time:** `ffmpeg version 8.1-3+b1 Copyright (c) 2000-2026 the FFmpeg developers`.

**Argv hygiene (T-07-04):** every element of `args` is either a literal `QStringLiteral`, a numeric setting (`QString::number(delay)` from `QSettings("export/gif_delay", 10)`), a `QTemporaryDir`-managed path (`tempDir.filePath(...)`), or the `QFileDialog`-returned `outputPath`. No user-typed strings flow in.

**Temp-dir lifecycle (T-07-05):** `QTemporaryDir tempDir(...)` is constructed at frame-write time. On success: `tempDir` falls out of scope at function return → `~QTemporaryDir` removes the dir. On failure: `tempDir.setAutoRemove(false)` is called BEFORE the early-return so the dir + its PNG sequence persist for diagnostic; the path is logged to the console-dock.

## Decisions Made

### Single-line auto-snapshot hook in handleLasops (Rule 2 — required for D-02)

Plan 05 inserts ONE conditional block at one location in xvue_qt_postscript.cpp:

```cpp
if (lasops == 0) {
    lasopsc_ = lasops;
    if (fpo_ != nullptr) {
        // Phase 7 Plan 05 (D-02): auto-snapshot capture hook.
        if (XvueExport::isCaptureActive()) {
            XvueExport::captureFrame();
        }
        std::fprintf(fpo_, "%s", concat_);
        std::fclose(fpo_);
        fpo_ = nullptr;
        // ...rest of close-branch...
```

The hook lives INSIDE the `if (fpo_ != nullptr)` arm so it cannot fire when the PS file is already closed/null. This is the minimal touch needed to honor the D-02 contract while preserving the verbatim-port intent of handleLasops (Plan 02 + 03).

### Test-only mock_ffmpeg branch skips PNG-sequence writes (Rule 2 — Plan deviation)

The plan body says the PNG-sequence write loop and the QProcess::execute step are separate. In the test-only ifdef, the production design would have the test framework still write 10000 PNGs at 800×600 to disk before the mocked ffmpeg dispatch — that's ~15 GiB peak disk + minutes of wall time, infeasible for a unit test.

Plan 05 takes the cleaner approach: when `g_ffmpeg_forced_exit ≥ 0` (test mock active), BOTH the PNG-write loop AND the QProcess::execute call are short-circuited. The hard-cap test then runs in ~120ms with ~120 MiB peak RAM and the hard-cap branch is still proven by post-condition assertions (isCaptureActive=false, pendingFrameCount=0). Production code path is unchanged — the test-only branch is gated entirely by `#ifdef XVUE_QT_TESTING` + `g_ffmpeg_forced_exit ≥ 0`.

### checkEnvAutoStart placement OUTSIDE call_once (Rule 1 — bug-class fix)

Initial draft put `XvueExport::checkEnvAutoStart()` INSIDE the call_once block. That meant only the very first `XvueApp::ensure()` call honored XVUE_ANIM. Subsequent ensure() calls (from later extern "C" entries) silently skipped the check.

Final implementation puts the call OUTSIDE call_once but INSIDE ensure(). It runs on every ensure() call, which is idempotent: `beginAnimation` simply resets the frame list, and the qgetenv() check is essentially free.

### beginAnimation does NOT call XvueApp::ensure (Rule 3 — blocking issue)

The plan body's `<action>` block specifies `beginAnimation()` should call `XvueApp::ensure()`. But because `checkEnvAutoStart` is invoked FROM inside `XvueApp::ensure()`, and `checkEnvAutoStart` calls `beginAnimation` when XVUE_ANIM=1, having `beginAnimation` call `ensure()` would create infinite recursion: ensure → checkEnvAutoStart → beginAnimation → ensure → … Plan 05 omits the `ensure()` call inside `beginAnimation`; all other call sites (menu toggle, tests) reach `beginAnimation` via paths that have already established the QApplication.

### Capture-Animation toggle outside Export submenu (Plan deviation — UX choice)

The Plan body's `<action>` Step C/D for Task 2 implies the Capture-Animation entry could go anywhere under File. Plan 05 places it AS A TOP-LEVEL FILE ITEM (after the Export submenu, separated by a separator), NOT inside the Export submenu. Rationale: the Export submenu houses output ACTIONS (verbs that write files); Capture Animation is a STATE TOGGLE — different shape. Putting it as a top-level checkable QAction matches Qt's idiomatic "state toggle near related output" (View → Toolbar / Console are similar checkables).

### EXPORT-06 grep gate scrubbed in this plan (Rule 2 — preview enforcement)

Plan 06's grep gate `grep -rn 'convert' xvue/qt/` must return empty. Plan 05's first draft of `xvue_qt_export.cpp` had a comment "D-03 legacy compat with bin/convertepsgif" inside the `endAnimation()` body. Plan 04 SUMMARY's "Preview EXPORT-06 grep" was already clean; this plan introduced a new violation. Caught + scrubbed before Task 2 commit so Plan 06's enforcement is a no-op.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 — Blocking] beginAnimation() must NOT call XvueApp::ensure()**
- **Found during:** Task 1 GREEN — wiring checkEnvAutoStart inside XvueApp::ensure().
- **Issue:** Plan body's `<action>` Step B has `beginAnimation()` call `XvueApp::ensure()` first. Combined with checkEnvAutoStart being invoked FROM ensure() (Step D), this creates infinite recursion.
- **Fix:** Removed the `XvueApp::ensure()` call from `beginAnimation()` body; documented in inline comment why ensure() is intentionally absent.
- **Files modified:** `xvue/qt/src/xvue_qt_export.cpp`
- **Verification:** Test 3 (XvueExport_env_XVUE_ANIM_activates_capture) PASSES — qputenv("XVUE_ANIM","1") + checkEnvAutoStart() activates capture without recursion.
- **Committed in:** `181efd8` (TDD GREEN).

**2. [Rule 2 — Missing critical functionality] Test-only mock skips PNG-sequence writes**
- **Found during:** Task 1 GREEN — first run of `XvueExport_hardcap_rejects_at_10000` hung in the 10000-frame PNG write loop.
- **Issue:** The plan body's mock injects a fake ffmpeg exit code but leaves the PNG-sequence write loop running. With 10000 frames at 800×600 that's ~15 GiB of PNGs and minutes of wall time — infeasible for a unit test.
- **Fix:** Extended the `#ifdef XVUE_QT_TESTING` mock to ALSO short-circuit the PNG-write loop when `g_ffmpeg_forced_exit ≥ 0`. Production behavior unchanged (the mock branch is dead code outside test builds). Hard-cap branch still proven by post-condition assertions (isCaptureActive=false, pendingFrameCount=0).
- **Files modified:** `xvue/qt/src/xvue_qt_export.cpp`, `xvue/qt/tests/test_xvue_qt_export.cpp` (also reduced canvas to 64×48 in the hard-cap test for safety belt).
- **Verification:** Test 5 PASSES in <1s with bounded memory.
- **Committed in:** `181efd8` (TDD GREEN).

**3. [Rule 3 — Blocking] xvue_qt_export_tests target needed XVUE_QT_TESTING define**
- **Found during:** Task 1 GREEN — first build of test binary failed with "‘resetForTesting’ is not a member of ‘XvueExport’" etc. The static lib has the define under `target_compile_definitions(xvueqt PRIVATE XVUE_QT_TESTING)` but the test binary needs the same define to see the test-only header declarations.
- **Fix:** Added `target_compile_definitions(xvue_qt_export_tests PRIVATE XVUE_QT_TESTING)` mirroring the per-module test targets (mail/elas/flui/ther/nlse) which already use this pattern.
- **Files modified:** `xvue/qt/tests/CMakeLists.txt`
- **Verification:** Test binary compiles; all 16 slots pass.
- **Committed in:** `181efd8` (TDD GREEN).

**4. [Rule 2 — Plan-prose API mismatch] Plan body's `dock->log(...)` does not exist; XvueConsoleDock uses `appendLine(...)`**
- **Found during:** Task 1 GREEN — drafting captureFrame's hard-cap log path.
- **Issue:** Plan body's `<action>` Step B references `dock->log(xvueT(MsgId::AnimationFrameHardCapHit))` repeatedly. `XvueConsoleDock` exposes `appendLine(const QString&)` (Phase 6.0 Plan 04); there is no `log()` method.
- **Fix:** Used `dock->appendLine(...)` in all 3 call sites (soft cap, hard cap, temp-dir failure). Documented inline.
- **Files modified:** `xvue/qt/src/xvue_qt_export.cpp`
- **Verification:** Compile + tests pass.
- **Committed in:** `181efd8` (TDD GREEN).

**5. [Rule 1 — Bug] EXPORT-06 grep gate caught a `convertepsgif` reference in Plan 05 code**
- **Found during:** Task 2 acceptance — `grep -rn 'convert' xvue/qt/` returned 1 hit on a comment in xvue_qt_export.cpp:381.
- **Issue:** First draft of the `endAnimation()` default-output comment said "D-03 legacy compat with bin/convertepsgif". The substring "convert" violates Plan 06's formal grep gate; Plan 04 SUMMARY's "Preview EXPORT-06 grep" was clean — this plan reintroduced the violation.
- **Fix:** Rephrased comment to remove the `convert` substring.
- **Files modified:** `xvue/qt/src/xvue_qt_export.cpp`
- **Verification:** `grep -rn 'convert' xvue/qt/` exits 1 (no match).
- **Committed in:** `fd17b6d` (Task 2).

---

**Total deviations:** 5 auto-fixed (1 Rule-1 grep gate, 2 Rule-2 missing critical functionality, 2 Rule-3 plan-prose vs reality / blocking issues). All five are surgical corrections that bring the implementation into alignment with the actual Qt API contracts and the wider phase contracts (D-02 auto-snapshot, EXPORT-06 grep gate, ABI invariant). No architectural changes, no scope creep, no ABI churn.

## Issues Encountered

- **No worktree base mismatch** (unlike Plans 02/03): worktree branch was correctly seated at `9518805` (the Plan 03 docs commit) per the system-provided merge-base.
- **Pre-existing testPerModuleGroupIsolation failure** in `xvue_qt_i18n_menu_prefs_tests` — same one Plan 02/03 SUMMARY documents. Verified by inspection (the test asserts QSettings group isolation behavior unrelated to Plan 05). Out of scope per SCOPE BOUNDARY rule.
- **Pre-existing -Wdangling-reference warnings** in `xvue/qt/src/xvue_qt_ther_actions.cpp` lines 191–193 — same warnings noted in Plan 07-01/02/03 SUMMARYs. Out of scope.

## TDD Gate Compliance

Plan 05's frontmatter is `type: execute` (not `type: tdd`), but Task 1 was tagged `tdd="true"` so the per-task TDD cycle applies:
- ✅ RED gate: `c9339d0` — `test(07-05): add failing XvueExport GIF/animation tests` — proves the new symbols are absent (compile failure).
- ✅ GREEN gate: `181efd8` — `feat(07-05): implement XvueExport animation surface + ffmpeg dispatch` — implements the symbols + bodies; all 9 new tests pass.
- (REFACTOR optional — not needed; the GREEN code already follows the verbatim-format-string + RAII-temp-dir + bilingual-message conventions established by Plans 02/03/04.)

## Self-Check: PASSED

**Files verified to exist:**
- `xvue/qt/src/xvue_qt_export.h` — MODIFIED, contains `beginAnimation` (1), `endAnimation` (1), `captureFrame` (1), `saveGifTo` (1), `usingNativeGifWriter` (1), `checkEnvAutoStart` (1), `kAnimFrameSoftCap = 100` (1), `kAnimFrameHardCap = 10000` (1), `resetForTesting` (1), `setNativeGifWriterForTesting` (1), `setFfmpegOverrideForTesting` (1)
- `xvue/qt/src/xvue_qt_export.cpp` — MODIFIED, contains `QProcess::execute(QStringLiteral("ffmpeg"` (2 — production + test mock), `"-framerate"` (1), `"frame_%04d.png"` (1), `QTemporaryDir` (5), `setAutoRemove(false)` (2 — failure paths), forbidden `system(`/`popen(`/`QProcess::start("sh"`/`"convert"` (0 — clean)
- `xvue/qt/src/xvue_qt_postscript.cpp` — MODIFIED, contains `XvueExport::isCaptureActive` (1), `XvueExport::captureFrame` (1), `#include "xvue_qt_export.h"` (1)
- `xvue/qt/src/xvue_qt_app.cpp` — MODIFIED, contains `XvueExport::checkEnvAutoStart` (1)
- `xvue/qt/src/xvue_qt_window.h` — MODIFIED, contains `actExportGif_` (1), `actCaptureAnimation_` (1), `onFileExportGif` (1), `onFileToggleCaptureAnimation` (1)
- `xvue/qt/src/xvue_qt_window.cpp` — MODIFIED, contains `setObjectName(QStringLiteral("FileExportGif"))` (1), `setObjectName(QStringLiteral("FileCaptureAnimation"))` (1), `actCaptureAnimation_->setCheckable(true)` (1), `XvueExport::onMenuExportGif()` (1), `XvueExport::onMenuToggleCapture()` (1)
- `xvue/qt/src/xvue_qt_i18n.h` — MODIFIED, contains `FileExportGif` (1), `FileCaptureAnimation` (1), `AnimationStarted` (1), `AnimationEncoding` (1), `AnimationFfmpegFailed` (1), `AnimationFrameSoftCapWarn` (1), `AnimationFrameHardCapHit` (1)
- `xvue/qt/src/xvue_qt_i18n.cpp` — MODIFIED, contains `"Encodage GIF en cours…"` (1) AND `"Encoding GIF…"` (1)
- `xvue/qt/tests/test_xvue_qt_export.cpp` — MODIFIED, contains 14 `void XvueExport_*` slot definitions (5 Plan 04 + 9 Plan 05)
- `xvue/qt/tests/CMakeLists.txt` — MODIFIED, contains `target_compile_definitions(xvue_qt_export_tests PRIVATE XVUE_QT_TESTING)` (1)

**Commits verified:**
- `c9339d0` (TDD RED) — FOUND in git log
- `181efd8` (TDD GREEN — Task 1) — FOUND in git log
- `fd17b6d` (Task 2) — FOUND in git log

**Build gates verified:**
- `cmake --build xvue/qt/build` (incremental, all targets) — exit 0
- `cd xvue/qt/build && ctest -R '^xvue_qt_export_tests$' --output-on-failure` — exit 0; "Passed 2.40 sec"
- `cd xvue/qt/build && ctest -R '^xvue_qt_postscript_tests$' --output-on-failure` — exit 0; "Passed 0.10 sec"
- `bash xvue/qt/cmake/verify_abi.sh xvue/qt/build/libxvueqt.a xvue/qt/include/xvue_qt_api.h` — exit 0; "nm count: 58 header count: 58" (ABI unchanged)
- `bin/ccxvue` (xvuelc.c→xvuelc.o) — exit 0 (T-07-08 quick non-regression)
- `git diff HEAD~3 HEAD -- xvue/xvuelc.c` — empty (T-07-08: legacy untouched)
- `git diff --name-only HEAD~3 HEAD -- xvue/video1.f xvue/videofin.f xvue/videonm.f` — empty (T-07-06 + D-17: LVIDEO untouched)
- `grep -rn 'convert' xvue/qt/` — empty (EXPORT-06 grep gate clean)
- `grep -rnE '(system|popen)\(' xvue/qt/src/xvue_qt_export.cpp` — empty (T-07-04: no shell-out)
- `grep -rn 'QProcess::start("sh"' xvue/qt/src/` — empty (T-07-04: no shell)

---
*Phase: 07-image-gif-and-postscript-export*
*Completed: 2026-05-04*
