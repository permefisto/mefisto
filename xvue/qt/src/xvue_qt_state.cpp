// xvue/qt/src/xvue_qt_state.cpp
// Phase 2 (D-16, D-17, D-20): applyPen() single source of truth for pen+brush;
// destructor enforces the single-long-lived-painter invariant from DRAW-01.
#include "xvue_qt_state.h"
#include <QPainter>
#include <QPixmap>
#include <algorithm>   // std::max

void XvueState::applyPen() {
    Qt::PenStyle qt_style = Qt::SolidLine;
    int          width    = pen_width_base_;
    switch (pen_style_) {
        case 0:
            qt_style = Qt::SolidLine;
            break;
        case 1:
            qt_style = Qt::DashLine;
            break;
        default:  // type 2 or anything else -> dash, double width (D-17)
            qt_style = Qt::DashLine;
            width    = std::max(1, pen_width_base_ * 2);
            break;
    }
    pen_.setColor(foreground_);
    pen_.setWidth(width);
    pen_.setStyle(qt_style);
    brush_ = QBrush(foreground_, Qt::SolidPattern);
    if (painter_ && painter_->isActive()) {
        painter_->setPen(pen_);
        painter_->setBrush(brush_);
    }
}

XvueState::~XvueState() {
    if (painter_) {
        if (painter_->isActive()) {
            painter_->end();
        }
        delete painter_;
        painter_ = nullptr;
    }
    delete backing_;
    backing_ = nullptr;
}
