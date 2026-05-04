// xvue/qt/src/xvue_qt_export.cpp
// Phase 7 Plan 04 (EXPORT-02, EXPORT-05): bodies — PNG/JPEG snapshots via
// QImage::save / QImageWriter; raster PDF via QPdfWriter at canvas-native
// 72 dpi. Pitfall 7: PDF page geometry uses XvueCanvas::width()/height()
// (Qt-LOGICAL pixels) for setPageSize; drawPixmap source is the
// HiDPI-scaled backing pixmap. Result: a 800×600 logical canvas always
// emits /MediaBox [0 0 800 600] regardless of devicePixelRatio.
//
// Failure surface (T-07-03 mitigation):
//   - Status-bar message + console-dock log on every failure path.
//   - QMessageBox::critical only when interactive=true (menu-triggered) per
//     CONTEXT.md `<discretion>` "QMessageBox::critical only when menu-triggered".
//   - Disk-full / nonexistent-dir surfaces as QImageWriter::write()==false
//     and is logged via the same path. Returns false; caller decides what next.
//
// Threading: every public method is main-thread-only via
// XVUE_QT_ASSERT_MAIN_THREAD. Mirrors the rest of xvue/qt/src/.
#include "xvue_qt_export.h"

#include "xvue_qt_api.h"           // XVUE_QT_ASSERT_MAIN_THREAD
#include "xvue_qt_app.h"
#include "xvue_qt_canvas.h"
#include "xvue_qt_console_dock.h"
#include "xvue_qt_i18n.h"
#include "xvue_qt_state.h"
#include "xvue_qt_window.h"

#include <QByteArray>
#include <QChar>
#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QImage>
#include <QImageWriter>
#include <QList>
#include <QMessageBox>
#include <QPageSize>
#include <QPainter>
#include <QPdfWriter>
#include <QPixmap>
#include <QProcess>
#include <QRect>
#include <QSettings>
#include <QSizeF>
#include <QStandardPaths>
#include <QStatusBar>
#include <QString>
#include <QStringList>
#include <QTemporaryDir>

namespace {

// Pull the current backing QPixmap pointer from XvueApp's window slot. Returns
// nullptr if the window is not yet alive (between xvfermer_ and xvinitgraphique_)
// or if the canvas has not yet had its first resizeEvent (backing_ is allocated
// lazily — Phase 2 D-04). Callers MUST nullptr-check before reading.
const QPixmap* currentBacking() {
    auto& slot = XvueApp::window_slot();
    if (!slot) return nullptr;
    XvueState* st = slot->state();
    if (!st) return nullptr;
    return st->backing_;   // raw pointer member; lifetime owned by XvueState
}

// Pull the current canvas widget. Used for LOGICAL pixel dims (Pitfall 7):
// XvueCanvas::width()/height() are Qt-logical pixels, not the HiDPI-scaled
// backing dims. PDF setPageSize uses these; drawPixmap source uses backing_.
const XvueCanvas* currentCanvas() {
    auto& slot = XvueApp::window_slot();
    if (!slot) return nullptr;
    return slot->canvas();
}

// Centralized failure surface (T-07-03):
//   - Console-dock log (Phase 6.0 D-07): always on if dock exists.
//   - Status-bar message (5s): always on if statusBar() exists.
//   - QMessageBox::critical: ONLY when interactive=true (menu-triggered).
//
// auto_snapshot path (Plan 05) will pass interactive=false so failed frame
// writes log to the console dock without freezing the auto-capture loop on
// a modal dialog.
void logFailure(const QString& msg, bool interactive) {
    if (auto& slot = XvueApp::window_slot(); slot) {
        if (auto* dock = slot->consoleDock()) {
            dock->appendLine(msg);
        }
        if (auto* sb = slot->statusBar()) {
            sb->showMessage(msg, /*ms=*/5000);
        }
    }
    if (interactive) {
        QMessageBox::critical(nullptr, xvueT(MsgId::ErrorMsgBoxTitle), msg);
    }
}

// QSettings("MEFISTO","xvue") namespace — matches Phase 6.0's QSettings usage
// and CONTEXT.md `<discretion>` documented keys (export/last_dir,
// export/jpeg_quality, export/gif_delay).
QSettings xvueSettings() {
    return QSettings(QStringLiteral("MEFISTO"), QStringLiteral("xvue"));
}

// After a successful save, remember the parent dir for next time so
// QFileDialog::getSaveFileName starts in the same directory the user
// last picked. CONTEXT.md D-03 hybrid-filename rule.
void rememberLastDir(const QString& chosen) {
    if (chosen.isEmpty()) return;
    QSettings s = xvueSettings();
    s.setValue(QStringLiteral("export/last_dir"),
               QFileInfo(chosen).absolutePath());
    s.sync();
}

QString lastDirOrHome() {
    QSettings s = xvueSettings();
    QString d = s.value(QStringLiteral("export/last_dir")).toString();
    if (!d.isEmpty()) return d;
    return QString::fromLocal8Bit(qgetenv("HOME"));
}

}  // namespace

// ============================================================================
// savePngTo — EXPORT-02 PNG snapshot.
// ============================================================================
bool XvueExport::savePngTo(const QString& path, bool interactive) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    XvueApp::ensure();

    const QPixmap* bg = currentBacking();
    if (!bg) {
        // No live canvas (window torn down or pre-first-resize). Surface as a
        // failure rather than producing an empty file — there is literally
        // nothing to save.
        logFailure(xvueT(MsgId::ExportNoCanvas), interactive);
        return false;
    }

    // QPixmap::toImage() is a pure snapshot (Qt copies the rendered surface).
    // Phase 2 D-04's long-lived state_->painter_ continues drawing on backing_
    // unaffected.
    QImage img = bg->toImage();
    if (!img.save(path, "PNG")) {
        logFailure(xvueT(MsgId::ExportPngFailed) +
                   QStringLiteral(" (%1)").arg(path),
                   interactive);
        return false;
    }
    return true;
}

// ============================================================================
// saveJpegTo — EXPORT-02 JPEG snapshot at QSettings("export/jpeg_quality", 90).
// ============================================================================
bool XvueExport::saveJpegTo(const QString& path, bool interactive) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    XvueApp::ensure();

    const QPixmap* bg = currentBacking();
    if (!bg) {
        logFailure(xvueT(MsgId::ExportNoCanvas), interactive);
        return false;
    }

    QImageWriter w(path, "jpeg");
    // Default 90 — matches CONTEXT.md `<discretion>` ("JPEG quality default.
    // Qt default (90)."). Configurable via QSettings("export/jpeg_quality").
    const int quality = xvueSettings()
        .value(QStringLiteral("export/jpeg_quality"), 90).toInt();
    w.setQuality(quality);

    if (!w.write(bg->toImage())) {
        logFailure(xvueT(MsgId::ExportJpegFailed) +
                   QStringLiteral(" (%1: %2)")
                       .arg(path, w.errorString()),
                   interactive);
        return false;
    }
    return true;
}

// ============================================================================
// savePdfTo — EXPORT-05 raster PDF at canvas-native 72 dpi.
//
// Pitfall 7 enforcement: setPageSize uses XvueCanvas::width()/height()
// (Qt-LOGICAL pixels). drawPixmap source is the HiDPI-scaled backing,
// painted into a logical-pixel rect; QPdfWriter at 72 dpi maps 1 logical
// pixel → 1 PDF point so /MediaBox emits the canvas dims directly.
// ============================================================================
bool XvueExport::savePdfTo(const QString& path, bool interactive) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    XvueApp::ensure();

    const QPixmap*    bg     = currentBacking();
    const XvueCanvas* canvas = currentCanvas();
    if (!bg || !canvas) {
        logFailure(xvueT(MsgId::ExportNoCanvas), interactive);
        return false;
    }

    // CRITICAL: logical canvas pixels — NOT backing dims. Phase 7 Pitfall 7.
    // On a HiDPI display backing_->width() == 2 * canvas->width(); using the
    // backing dims would emit /MediaBox [0 0 1600 1200] for a 800×600 canvas.
    const int logicalW = canvas->width();
    const int logicalH = canvas->height();

    QPdfWriter writer(path);
    writer.setResolution(72);   // 1 logical pixel == 1 PDF point at 72 dpi
    writer.setPageSize(QPageSize(QSizeF(logicalW, logicalH),
                                 QPageSize::Point));

    // Some QPdfWriter failure modes (e.g. unwritable path) only surface inside
    // QPainter::begin. Detect via isActive() after construction.
    QPainter p;
    if (!p.begin(&writer)) {
        logFailure(xvueT(MsgId::ExportPdfFailed) +
                   QStringLiteral(" (%1)").arg(path),
                   interactive);
        return false;
    }
    // Draw backing pixmap into the (logicalW × logicalH) page rect. Qt scales
    // the physical-pixel pixmap onto the logical-point destination — raster
    // PDF per CONTEXT.md D-13. True-vector PDF deferred.
    p.drawPixmap(QRect(0, 0, logicalW, logicalH), *bg);
    p.end();
    return true;
}

// ============================================================================
// onMenuExportPng / onMenuExportJpeg / onMenuExportPdf — File menu slots.
//   - Default dir from QSettings("export/last_dir") or $HOME.
//   - Default filename "canvas.{png,jpg,pdf}".
//   - On accept: rememberLastDir(path), then saveXxxTo(path, interactive=true).
// ============================================================================
void XvueExport::onMenuExportPng() {
    XVUE_QT_ASSERT_MAIN_THREAD();
    const QString dir  = lastDirOrHome();
    const QString seed = dir + QStringLiteral("/canvas.png");
    const QString path = QFileDialog::getSaveFileName(
        nullptr,
        xvueT(MsgId::FileExportPngDialogTitle),
        seed,
        xvueT(MsgId::FileExportPngFilter));
    if (path.isEmpty()) return;
    rememberLastDir(path);
    savePngTo(path, /*interactive=*/true);
}

void XvueExport::onMenuExportJpeg() {
    XVUE_QT_ASSERT_MAIN_THREAD();
    const QString dir  = lastDirOrHome();
    const QString seed = dir + QStringLiteral("/canvas.jpg");
    const QString path = QFileDialog::getSaveFileName(
        nullptr,
        xvueT(MsgId::FileExportJpegDialogTitle),
        seed,
        xvueT(MsgId::FileExportJpegFilter));
    if (path.isEmpty()) return;
    rememberLastDir(path);
    saveJpegTo(path, /*interactive=*/true);
}

void XvueExport::onMenuExportPdf() {
    XVUE_QT_ASSERT_MAIN_THREAD();
    const QString dir  = lastDirOrHome();
    const QString seed = dir + QStringLiteral("/canvas.pdf");
    const QString path = QFileDialog::getSaveFileName(
        nullptr,
        xvueT(MsgId::FileExportPdfDialogTitle),
        seed,
        xvueT(MsgId::FileExportPdfFilter));
    if (path.isEmpty()) return;
    rememberLastDir(path);
    savePdfTo(path, /*interactive=*/true);
}

// ============================================================================
// Phase 7 Plan 05 (EXPORT-03) — animated GIF + auto-snapshot.
//
// The animation surface adds NO new extern "C" symbols (D-01: ABI stays at
// 58). It is all C++-internal — Fortran continues to call xvpostscript_(0)
// at frame boundaries, and PsEmitter::handleLasops(0) (Plan 02) hooks into
// XvueExport::captureFrame() when isCaptureActive() is true (Plan 05 step C
// in the PLAN body).
//
// Process-wide singletons live in the anonymous namespace below. Main-thread
// only — every public method asserts via XVUE_QT_ASSERT_MAIN_THREAD.
// ============================================================================
namespace {

// ---- Animation state — process-wide singletons. ----
bool          g_capture_active        = false;
QList<QImage> g_pending_frames;
int           g_warned_frames_softcap = 0;   // sticky warn-once flag, reset on begin
#ifdef XVUE_QT_TESTING
int           g_ffmpeg_forced_exit = -1;     // -1 = real ffmpeg; else mock rc
bool          g_native_gif_forced  = false;  // force usingNativeGifWriter()=true
#endif

// QStandardPaths::TempLocation — typically /tmp on Linux. xvue-gif-XXXXXX is
// created inside it via QTemporaryDir.
QString animTempDirRoot() {
    return QStandardPaths::writableLocation(QStandardPaths::TempLocation);
}

// Emit a "GIF status" message via console-dock + status-bar with no modal.
// Mirrors logFailure() but without the QMessageBox surface — animation
// path is mostly silent on success and warn-only on cap-hits.
void logAnimStatus(const QString& msg) {
    if (auto& slot = XvueApp::window_slot(); slot) {
        if (auto* dock = slot->consoleDock()) {
            dock->appendLine(msg);
        }
        if (auto* sb = slot->statusBar()) {
            sb->showMessage(msg, /*ms=*/3000);
        }
    }
}

}  // anonymous namespace

bool XvueExport::isCaptureActive() {
    return g_capture_active;
}

int XvueExport::pendingFrameCount() {
    return static_cast<int>(g_pending_frames.size());
}

void XvueExport::beginAnimation() {
    XVUE_QT_ASSERT_MAIN_THREAD();
    // Note: deliberately do NOT call XvueApp::ensure() here — checkEnvAutoStart()
    // is invoked FROM inside XvueApp::ensure() so recursing back through
    // ensure() would trigger the call_once pre-condition guard. All other
    // callers (menu toggle, tests) reach beginAnimation() via paths that
    // have already established the QApplication.
    g_capture_active = true;
    g_pending_frames.clear();
    g_warned_frames_softcap = 0;
    logAnimStatus(xvueT(MsgId::AnimationStarted));
}

void XvueExport::captureFrame() {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!g_capture_active) return;

    const QPixmap* bg = currentBacking();
    if (!bg) return;

    // T-07-03 hard cap: reject frame, force end-of-animation.
    if (static_cast<int>(g_pending_frames.size()) >= kAnimFrameHardCap) {
        if (auto& slot = XvueApp::window_slot(); slot) {
            if (auto* dock = slot->consoleDock()) {
                dock->appendLine(xvueT(MsgId::AnimationFrameHardCapHit));
            }
        }
        endAnimation();
        return;
    }

    g_pending_frames.append(bg->toImage());

    // T-07-03 soft warn at 100 (warn-once).
    if (!g_warned_frames_softcap &&
        static_cast<int>(g_pending_frames.size()) >= kAnimFrameSoftCap) {
        g_warned_frames_softcap = 1;
        if (auto& slot = XvueApp::window_slot(); slot) {
            if (auto* dock = slot->consoleDock()) {
                dock->appendLine(xvueT(MsgId::AnimationFrameSoftCapWarn));
            }
        }
    }
}

void XvueExport::endAnimation() {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!g_capture_active && g_pending_frames.isEmpty()) return;
    g_capture_active = false;
    if (g_pending_frames.isEmpty()) return;

    // Default output: cwd/animation.gif (D-03 legacy filename compat with
    // the legacy post-processor under bin/ — see CONTEXT.md D-03).
    const QString defaultPath = QDir::current().filePath(QStringLiteral("animation.gif"));
    saveGifTo(defaultPath, /*interactive=*/false);
}

bool XvueExport::usingNativeGifWriter() {
#ifdef XVUE_QT_TESTING
    if (g_native_gif_forced) return true;
#endif
    const auto formats = QImageWriter::supportedImageFormats();
    for (const auto& f : formats) {
        if (QString::fromLatin1(f).toLower() == QStringLiteral("gif")) {
            return true;
        }
    }
    return false;
}

bool XvueExport::saveGifTo(const QString& outputPath, bool interactive) {
    XVUE_QT_ASSERT_MAIN_THREAD();

    // saveGifTo always closes capture state (whether it succeeds or fails),
    // so menu-driven QFileDialog flow + auto-snapshot endAnimation flow
    // both leave isCaptureActive() == false on return.
    g_capture_active = false;

    if (g_pending_frames.isEmpty()) {
        logFailure(xvueT(MsgId::AnimationNoFrames), interactive);
        return false;
    }

    // Status-bar feedback (D-11): "Encoding GIF…" / "Encodage GIF en cours…".
    if (auto& slot = XvueApp::window_slot(); slot) {
        if (auto* sb = slot->statusBar()) {
            sb->showMessage(xvueT(MsgId::AnimationEncoding), /*ms=*/0);
        }
    }

    bool ok = false;

    if (usingNativeGifWriter()) {
        // ----- D-10 native QImageWriter multi-frame path -----
        // Per PROBE.md this branch is dead today on Debian/Qt 6.10.2 (no GIF
        // write plugin) but kept for portability if a future Debian add-on
        // exposes one. Plan 06 unit-tests the dispatch logic.
        QImageWriter w(outputPath, "gif");
        const int delay = xvueSettings()
            .value(QStringLiteral("export/gif_delay"), 10).toInt();
        w.setText(QStringLiteral("Delay"), QString::number(delay));
        bool any_failed = false;
        for (const auto& img : g_pending_frames) {
            if (!w.write(img)) { any_failed = true; break; }
        }
        ok = !any_failed;
    } else {
        // ----- D-11 ffmpeg fallback path (the realized path on this host) -----
        // Frame buffer dir under $TMPDIR/xvue-gif-XXXXXX/ via QTemporaryDir.
        const QString tmpRoot = animTempDirRoot();
        QTemporaryDir tempDir(QDir(tmpRoot).filePath(QStringLiteral("xvue-gif-XXXXXX")));
        if (!tempDir.isValid()) {
            logFailure(xvueT(MsgId::AnimationTempDirFailed), interactive);
            g_pending_frames.clear();
            g_warned_frames_softcap = 0;
            return false;
        }

#ifdef XVUE_QT_TESTING
        // When the test-only ffmpeg override is set, both the PNG-sequence
        // write AND the QProcess::execute step are mocked out. This keeps
        // the hard-cap (10000 frame) test bounded in time and disk usage.
        // Behavior contract preserved: ok branch acts as if PNG write +
        // ffmpeg encode succeeded; failure branch sets autoRemove(false)
        // so T-07-05 is still observable in tests.
        const bool mock_ffmpeg = (g_ffmpeg_forced_exit >= 0);
#else
        const bool mock_ffmpeg = false;
#endif

        if (!mock_ffmpeg) {
            // Write PNG sequence frame_0001.png … frame_NNNN.png.
            bool any_failed = false;
            int idx = 0;
            for (const auto& img : g_pending_frames) {
                ++idx;
                const QString fname = tempDir.filePath(
                    QString("frame_%1.png").arg(idx, 4, 10, QChar('0')));
                if (!img.save(fname, "PNG")) { any_failed = true; break; }
            }
            if (any_failed) {
                tempDir.setAutoRemove(false);   // T-07-05: keep on failure
                if (auto& slot = XvueApp::window_slot(); slot) {
                    if (auto* dock = slot->consoleDock()) {
                        dock->appendLine(
                            QStringLiteral("Frame write failed; tempdir kept at ")
                            + tempDir.path());
                    }
                }
                g_pending_frames.clear();
                g_warned_frames_softcap = 0;
                return false;
            }
        }

        // ffmpeg dispatch (T-07-04): argv built from compile-time constants
        // + the QTemporaryDir-controlled path + the QFileDialog-returned
        // outputPath. NO shell, NO user-typed args, NO interpolation.
        const int delay = xvueSettings()
            .value(QStringLiteral("export/gif_delay"), 10).toInt();
        const QStringList args{
            QStringLiteral("-y"),
            QStringLiteral("-framerate"), QString::number(delay),
            QStringLiteral("-i"), tempDir.filePath(QStringLiteral("frame_%04d.png")),
            outputPath
        };
#ifdef XVUE_QT_TESTING
        const int rc = mock_ffmpeg
                       ? g_ffmpeg_forced_exit
                       : QProcess::execute(QStringLiteral("ffmpeg"), args);
#else
        (void)args;  // silence -Wunused-variable when mock_ffmpeg branch is dead
        const int rc = QProcess::execute(QStringLiteral("ffmpeg"), args);
#endif
        if (rc != 0) {
            tempDir.setAutoRemove(false);   // T-07-05: keep on failure
            const QString msg = xvueT(MsgId::AnimationFfmpegFailed)
                + QStringLiteral(" (rc=%1, tempdir=%2)")
                      .arg(rc).arg(tempDir.path());
            logFailure(msg, interactive);
            g_pending_frames.clear();
            g_warned_frames_softcap = 0;
            return false;
        }
        ok = true;
        // tempDir auto-removed by RAII destructor on success (T-07-05).
    }

    g_pending_frames.clear();
    g_warned_frames_softcap = 0;

    if (auto& slot = XvueApp::window_slot(); slot) {
        if (auto* sb = slot->statusBar()) {
            sb->showMessage(ok ? xvueT(MsgId::AnimationDone)
                               : xvueT(MsgId::AnimationFailed),
                            /*ms=*/5000);
        }
    }
    if (!ok) {
        logFailure(xvueT(MsgId::AnimationFailed), interactive);
    }
    return ok;
}

void XvueExport::onMenuExportGif() {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!isCaptureActive() && g_pending_frames.isEmpty()) {
        logFailure(xvueT(MsgId::AnimationNoFrames), /*interactive=*/true);
        return;
    }
    const QString dir  = lastDirOrHome();
    const QString seed = dir + QStringLiteral("/animation.gif");
    const QString path = QFileDialog::getSaveFileName(
        nullptr,
        xvueT(MsgId::FileExportGifDialogTitle),
        seed,
        xvueT(MsgId::FileExportGifFilter));
    if (path.isEmpty()) return;
    rememberLastDir(path);
    saveGifTo(path, /*interactive=*/true);
}

void XvueExport::onMenuToggleCapture() {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (isCaptureActive()) endAnimation();
    else                   beginAnimation();
}

void XvueExport::checkEnvAutoStart() {
    // D-02: XVUE_ANIM=1 toggles the same surface as File → Capture Animation.
    // Called once from XvueApp::ensure() so users can drive testa/wave or
    // testa/cavity2d animations purely via env var with no menu interaction.
    if (qgetenv("XVUE_ANIM").trimmed() == QByteArray("1")) {
        beginAnimation();
    }
}

#ifdef XVUE_QT_TESTING
void XvueExport::resetForTesting() {
    g_capture_active = false;
    g_pending_frames.clear();
    g_warned_frames_softcap = 0;
    g_ffmpeg_forced_exit = -1;
    g_native_gif_forced  = false;
}
void XvueExport::setNativeGifWriterForTesting(bool v) {
    g_native_gif_forced = v;
}
void XvueExport::setFfmpegOverrideForTesting(int forced_exit_code) {
    g_ffmpeg_forced_exit = forced_exit_code;
}
#endif
