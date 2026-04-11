// xvue/qt/src/xvue_qt_window.h
// Phase 1 (D-02, D-04, D-15, Pitfall 6, Pitfall 8): bare QMainWindow.
// No menu bar, no toolbar, no status bar, no dock widgets — those arrive in
// Phase 6 (ROADMAP.md). XvueWindow owns the XvueState and the XvueCanvas.
#pragma once
#include <QMainWindow>
#include "xvue_qt_state.h"

class XvueCanvas;

class XvueWindow : public QMainWindow {
    Q_OBJECT  // Pitfall 8: required even with no signals/slots in Phase 1.
public:
    explicit XvueWindow(QWidget* parent = nullptr);
    ~XvueWindow() override;

    XvueState*  state()  { return &state_; }
    XvueCanvas* canvas() { return canvas_; }

private:
    XvueState     state_{};    // D-04: single-field state struct
    XvueCanvas*   canvas_ = nullptr;  // Qt-owned via setCentralWidget
};
