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
};
