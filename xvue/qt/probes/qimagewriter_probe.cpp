// xvue/qt/probes/qimagewriter_probe.cpp
// Phase 7 Plan 01 (EXPORT-01): probe binary that prints
// QImageWriter::supportedImageFormats() and the GIF write-support flag to
// stdout. Output is captured by bin/cb_probe_qt into PROBE.md.
//
// Per CONTEXT.md D-09: standalone executable, NOT linked into libxvueqt.a.
// Per CLAUDE.md: no Fortran here; pure Qt C++17.
#include <QGuiApplication>
#include <QImageWriter>
#include <QImageIOHandler>
#include <iostream>

int main(int argc, char** argv) {
    QGuiApplication app(argc, argv);
    std::cout << "qt_version=" << qVersion() << "\n";
    std::cout << "supported_write_formats=";
    for (const auto& f : QImageWriter::supportedImageFormats()) {
        std::cout << f.toStdString() << " ";
    }
    std::cout << "\n";
    QImageWriter probe("/tmp/qimagewriter_probe.gif", "gif");
    std::cout << "gif_write_supported=" << probe.canWrite() << "\n";
    std::cout << "gif_animation_supported="
              << probe.supportsOption(QImageIOHandler::Animation) << "\n";
    return 0;
}
