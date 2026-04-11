// xvue/qt/src/xvue_qt_canvas.cpp
// Phase 2 (D-06): paintEvent = drawPixmap blit ONLY. Never mutate state here.
// Phase 2 (D-04, D-07, D-08): resizeEvent reallocates backing, preserves
// content top-left, and recreates the long-lived QPainter per DRAW-01.
#include "xvue_qt_canvas.h"
#include "xvue_qt_state.h"
#include <QPainter>
#include <QPaintEvent>
#include <QResizeEvent>
#include <QPixmap>

XvueCanvas::XvueCanvas(XvueState* state, QWidget* parent)
    : QWidget(parent), state_(state)
{
    // D-04: DO NOT allocate backing_ here. Qt 6 guarantees resizeEvent fires
    // before the first paintEvent on X11/Wayland (Pitfall 6). The first
    // resizeEvent after construction performs the initial allocation.
}

XvueCanvas::~XvueCanvas() = default;

void XvueCanvas::paintEvent(QPaintEvent* event) {
    Q_UNUSED(event);
    // D-06 + Pitfall 1: paintEvent body is ONE drawPixmap call; never more.
    // Defensive null-check guards against a pathological pre-first-resize
    // paint (should not happen on X11/Wayland but costs nothing).
    if (state_ && state_->backing_) {
        QPainter(this).drawPixmap(0, 0, *state_->backing_);
    }
}

void XvueCanvas::resizeEvent(QResizeEvent* event) {
    Q_UNUSED(event);
    if (!state_) return;

    const qreal dpr = devicePixelRatioF();
    const QSize logical = size();
    const QSize device  = logical * dpr;

    // (a) end old painter BEFORE deleting its device (Pitfall 7)
    if (state_->painter_ && state_->painter_->isActive()) {
        state_->painter_->end();
    }

    // (b)(c) allocate new backing in device pixels, tag with DPR so QPainter
    //       on it takes LOGICAL coordinates (Pitfall 8).
    QPixmap* old_backing = state_->backing_;
    QPixmap* new_backing = new QPixmap(device);
    new_backing->setDevicePixelRatio(dpr);

    // (d)(e) fill background, then blit old content top-left (DRAW-09,
    //       README_RESIZE.md invariants 2 and 3). One scoped temp painter.
    {
        QPainter tmp(new_backing);
        tmp.fillRect(new_backing->rect(), state_->background_);
        if (old_backing) {
            // drawPixmap at (0, 0) in LOGICAL coordinates — Qt handles DPR.
            tmp.drawPixmap(0, 0, *old_backing);
        }
    }  // ~QPainter tmp — end() called

    // (f) delete old, (g) swap
    delete old_backing;
    state_->backing_ = new_backing;

    // (h)(i) recreate/begin the long-lived painter on the new backing
    if (!state_->painter_) {
        state_->painter_ = new QPainter();
    }
    state_->painter_->begin(new_backing);

    // (j) DRAW-08 + Pitfall 5: hints do not carry across begin()/end()
    state_->painter_->setRenderHint(QPainter::Antialiasing, true);

    // (k) re-push pen+brush through applyPen() (D-22, Pitfall 5)
    state_->applyPen();
}
