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

#include <QFileDialog>
#include <QFileInfo>
#include <QImage>
#include <QImageWriter>
#include <QMessageBox>
#include <QPageSize>
#include <QPainter>
#include <QPdfWriter>
#include <QPixmap>
#include <QRect>
#include <QSettings>
#include <QSizeF>
#include <QStatusBar>
#include <QString>

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
