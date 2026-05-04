// xvue/qt/src/xvue_qt_export.cpp
// Phase 7 Plan 04 (EXPORT-02, EXPORT-05): TDD RED stub — bodies return false
// (PNG/JPEG/PDF not yet wired). The Plan 04 GREEN commit fills the bodies
// using QImageWriter / QPdfWriter / QFileDialog. Tests in
// xvue/qt/tests/test_xvue_qt_export.cpp drive the contract.
#include "xvue_qt_export.h"
#include "xvue_qt_api.h"  // XVUE_QT_ASSERT_MAIN_THREAD

bool XvueExport::savePngTo(const QString& /*path*/, bool /*interactive*/) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    return false;
}

bool XvueExport::saveJpegTo(const QString& /*path*/, bool /*interactive*/) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    return false;
}

bool XvueExport::savePdfTo(const QString& /*path*/, bool /*interactive*/) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    return false;
}

void XvueExport::onMenuExportPng()  {
    XVUE_QT_ASSERT_MAIN_THREAD();
}
void XvueExport::onMenuExportJpeg() {
    XVUE_QT_ASSERT_MAIN_THREAD();
}
void XvueExport::onMenuExportPdf()  {
    XVUE_QT_ASSERT_MAIN_THREAD();
}
