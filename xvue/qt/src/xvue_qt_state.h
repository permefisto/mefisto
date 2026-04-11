// xvue/qt/src/xvue_qt_state.h
// Phase 1 (D-04): runtime state holder; grown additively in Phase 2.
// Phase 2 (D-04, D-05, D-16, D-20): adds backing_, painter_, foreground_,
// pen_, brush_, pen_style_, pen_width_base_, applyPen() + destructor.
// NEVER reorder existing fields — background_ stays first.
#pragma once
#include <QColor>
#include <QPen>
#include <QBrush>
#include <Qt>

class QPixmap;
class QPainter;

struct XvueState {
    // Phase 1 — untouched. Matches legacy BlackPixel default
    // (xvue/xvuelc.c:935). Keeping the same default prevents visible drift in
    // the Phase 8 A/B validation.
    QColor background_ = Qt::black;

    // Phase 2 (D-21). Hardcoded white for Phase 2; Phase 3 unlocks via xvcouleurs_.
    // TODO(phase 3): make mutable via xvcouleurs_ / xvCouleursImposees_.
    QColor foreground_ = Qt::white;

    // Phase 2 (D-04, D-05). Allocated on first XvueCanvas::resizeEvent.
    // painter_ lives exactly as long as backing_.
    QPixmap*  backing_ = nullptr;
    QPainter* painter_ = nullptr;

    // Phase 2 (D-16, D-17, D-18, D-20). Rebuilt by applyPen().
    QPen     pen_;
    QBrush   brush_          = QBrush(Qt::white, Qt::SolidPattern);
    int      pen_style_      = 0;   // 0=solid, 1=dash, 2=dash-double
    int      pen_width_base_ = 0;   // 0 = cosmetic pen (1 device px)

    // Rebuild pen_ and brush_ from pen_style_/pen_width_base_/foreground_,
    // push to painter_ if the painter is active.
    void applyPen();

    // Ends painter_ (if active), deletes painter_, deletes backing_.
    ~XvueState();
};
