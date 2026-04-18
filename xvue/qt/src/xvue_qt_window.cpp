// xvue/qt/src/xvue_qt_window.cpp
// Phase 1 (D-02): 800x600 logical pixels, title "MEFISTO", central widget
// is XvueCanvas. Nothing else. Pitfall 6 applies — no chrome.
// Phase 5 (D-02, EVENT-01): also constructs the XvueEventBridge and installs
// it as event filter on the canvas. Bridge is a QObject child of the window
// so Qt parent-child destruction cleans it up automatically on ~XvueWindow.
#include "xvue_qt_window.h"
#include "xvue_qt_canvas.h"
#include "xvue_qt_event.h"

XvueWindow::XvueWindow(QWidget* parent)
    : QMainWindow(parent)
{
    setWindowTitle(QStringLiteral("MEFISTO"));
    resize(800, 600);                         // D-02: logical pixels
    canvas_ = new XvueCanvas(&state_, this);  // D-15: raw pointer into state_
    setCentralWidget(canvas_);

    // Phase 5 (D-02, EVENT-01). Bridge owns the event filter; lifetime
    // matches the window. Installed on the canvas so QEvent::MouseMove /
    // ButtonPress / KeyPress reach the filter before the canvas's default
    // handlers (Pitfall 3). Parent=this -> Qt parent-child destruction
    // deletes bridge_ in ~XvueWindow, no manual delete needed.
    bridge_ = new XvueEventBridge(canvas_, this);
    canvas_->installEventFilter(bridge_);
}

XvueWindow::~XvueWindow() = default;
