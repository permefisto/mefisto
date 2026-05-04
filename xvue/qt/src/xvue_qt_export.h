// xvue/qt/src/xvue_qt_export.h
// Phase 7 Plan 04 (EXPORT-02, EXPORT-05): XvueExport — PNG/JPEG/PDF static
// helpers driven by File → Export → … menu actions. Plan 05 grows the file
// with GIF assembly + auto-snapshot beginAnimation/endAnimation.
//
// Pure-static class (no state — XvueApp owns the canvas; QSettings persists
// last_dir / jpeg_quality / gif_delay). Bridge-as-XvueApp-singleton pattern
// is NOT used here because there is no per-process state; a static class
// with `XvueApp::ensure()` calls inside each method is the right shape.
//
// Trust boundaries (07-04-PLAN.md threat model):
//   - User -> QFileDialog: Qt validates path; overwrite prompt is built-in.
//   - QImageWriter / QPdfWriter -> filesystem: writes only to dialog-returned
//     path. Disk-full / write-failure: status-bar message + console-dock log,
//     QMessageBox::critical only when menu-triggered (interactive=true).
//
// HiDPI math (Pitfall 7 from CONTEXT.md D-15):
//   - PNG/JPEG snapshot at backing pixel dims (devicePixelRatio-scaled).
//   - PDF page geometry uses LOGICAL canvas pixels (XvueCanvas::width() /
//     height(), Qt-logical) for setPageSize at 72 dpi. drawPixmap source
//     is the backing_ pixmap (raster-faithful).
#pragma once
#include <QString>

class XvueExport {
public:
    // EXPORT-02 — PNG single-frame snapshot of XvueCanvas backing.
    // Returns true on success. On failure: status-bar + console-dock log,
    // QMessageBox::critical only if interactive. The "interactive" boolean
    // is true for menu-triggered, false for unit-test/silent paths.
    static bool savePngTo(const QString& path, bool interactive = true);

    // EXPORT-02 — JPEG snapshot at QSettings("export/jpeg_quality", 90).
    static bool saveJpegTo(const QString& path, bool interactive = true);

    // EXPORT-05 — Raster PDF at canvas-native 72 dpi (Pitfall 7: page
    // geometry uses XvueCanvas::width()/height() — Qt-logical, NOT
    // backing_->width() which is HiDPI-scaled physical pixels).
    static bool savePdfTo(const QString& path, bool interactive = true);

    // Menu-action entry points: prompt via QFileDialog (QSettings-remembered
    // last_dir), call the corresponding saveXxxTo. These are the slots the
    // File → Export → … QActions connect to.
    static void onMenuExportPng();
    static void onMenuExportJpeg();
    static void onMenuExportPdf();

    // ---- Phase 7 Plan 05 (EXPORT-03) — animated GIF + auto-snapshot ----
    //
    // beginAnimation() activates the in-memory frame buffer. While active,
    // PsEmitter::handleLasops(0) (the close-PS-file branch) calls
    // captureFrame() to snapshot the current backing pixmap. Frame caps
    // (T-07-03 mitigation): warn at kAnimFrameSoftCap=100 (status-bar
    // message), force-end at kAnimFrameHardCap=10000.
    //
    // saveGifTo dispatches via probe-driven branch (D-09/D-10/D-11):
    //   - Native QImageWriter "gif" if usingNativeGifWriter() returns true
    //     (defensive runtime re-check beyond Plan 01 PROBE.md).
    //   - Otherwise ffmpeg fallback via QProcess::execute with a constants-
    //     only QStringList (T-07-04: no shell, no user-typed argv).
    // Frame buffer = QTemporaryDir at $TMPDIR/xvue-gif-XXXXXX/. On success:
    // RAII auto-removed. On failure: setAutoRemove(false) + path logged
    // to console-dock (T-07-05).
    static void beginAnimation();
    static void endAnimation();           // flushes to default cwd/animation.gif
    static bool isCaptureActive();
    static void captureFrame();           // snapshot current backing into the list
    static int  pendingFrameCount();      // for tests + telemetry

    // Menu / env-var entry points.
    static void onMenuExportGif();        // QFileDialog-prompt path
    static void onMenuToggleCapture();    // File → Capture Animation toggle (D-02)
    static void checkEnvAutoStart();      // called from XvueApp::ensure — reads XVUE_ANIM

    // Probe-driven branch selector (D-09/D-10/D-11): true iff Qt's
    // QImageWriter supports "gif" write at runtime. Defensive runtime
    // re-check; PROBE.md is the kickoff snapshot.
    static bool usingNativeGifWriter();

    // GIF assembly — sync, blocking on GUI thread (D-11). outputPath is
    // an absolute file path; interactive=true gates QMessageBox surfaces.
    static bool saveGifTo(const QString& outputPath, bool interactive = true);

    // Frame caps (T-07-03 mitigation).
    static constexpr int kAnimFrameSoftCap = 100;
    static constexpr int kAnimFrameHardCap = 10000;

#ifdef XVUE_QT_TESTING
    // Test-only hooks. Compiled in only when the library is built with
    // XVUE_QT_TESTING (set by xvue/qt/CMakeLists.txt when XVUE_QT_BUILD_TESTS
    // is ON — same gate that exposes the menu_file_parser cache reset).
    static void resetForTesting();
    static void setNativeGifWriterForTesting(bool v);
    static void setFfmpegOverrideForTesting(int forced_exit_code);
#endif
};
