// xvue/qt/src/xvue_qt_window.cpp
// Phase 1 (D-02): 800x600 logical pixels, title "MEFISTO", central widget
// is XvueCanvas. Nothing else. Pitfall 6 applies — no chrome.
#include "xvue_qt_window.h"
#include "xvue_qt_canvas.h"

XvueWindow::XvueWindow(QWidget* parent)
    : QMainWindow(parent)
{
    setWindowTitle(QStringLiteral("MEFISTO"));
    resize(800, 600);                         // D-02: logical pixels
    canvas_ = new XvueCanvas(&state_, this);  // D-15: raw pointer into state_
    setCentralWidget(canvas_);
}

XvueWindow::~XvueWindow() = default;
