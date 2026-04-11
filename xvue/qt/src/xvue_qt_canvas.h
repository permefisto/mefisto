// xvue/qt/src/xvue_qt_canvas.h
// Phase 1 (D-05) + Phase 2 (D-04..D-08): central widget. paintEvent blits
// the persistent backing pixmap owned by XvueState; resizeEvent reallocates
// that backing, preserves top-left content (DRAW-09), ends+restarts the
// long-lived QPainter, and re-applies Antialiasing + pen + brush.
#pragma once
#include <QWidget>

struct XvueState;
class QResizeEvent;

class XvueCanvas : public QWidget {
    Q_OBJECT  // Pitfall 8
public:
    explicit XvueCanvas(XvueState* state, QWidget* parent = nullptr);
    ~XvueCanvas() override;

protected:
    void paintEvent(QPaintEvent* event) override;
    void resizeEvent(QResizeEvent* event) override;   // Phase 2 (D-04, D-07)

private:
    XvueState* state_ = nullptr;  // raw pointer — owned by XvueWindow (D-15)
};
