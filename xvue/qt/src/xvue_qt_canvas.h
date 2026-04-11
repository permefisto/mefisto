// xvue/qt/src/xvue_qt_canvas.h
// Phase 1 (D-05, D-15, Pitfall 8): central widget. paintEvent fills
// background with state_->background_. Phase 2 will swap this body to
// drawPixmap(0, 0, *backing_) as a one-line change.
#pragma once
#include <QWidget>

struct XvueState;

class XvueCanvas : public QWidget {
    Q_OBJECT  // Pitfall 8
public:
    explicit XvueCanvas(XvueState* state, QWidget* parent = nullptr);
    ~XvueCanvas() override;

protected:
    void paintEvent(QPaintEvent* event) override;

private:
    XvueState* state_ = nullptr;  // raw pointer — owned by XvueWindow (D-15)
};
