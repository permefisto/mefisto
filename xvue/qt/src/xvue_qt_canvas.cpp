// xvue/qt/src/xvue_qt_canvas.cpp
// Phase 1 (D-05): paintEvent body is exactly one QPainter fillRect call.
// Phase 2 will replace the body wholesale with a drawPixmap blit.
#include "xvue_qt_canvas.h"
#include "xvue_qt_state.h"
#include <QPainter>
#include <QPaintEvent>

XvueCanvas::XvueCanvas(XvueState* state, QWidget* parent)
    : QWidget(parent), state_(state)
{
}

XvueCanvas::~XvueCanvas() = default;

void XvueCanvas::paintEvent(QPaintEvent* event) {
    Q_UNUSED(event);
    if (state_) {
        QPainter(this).fillRect(rect(), state_->background_);
    }
}
