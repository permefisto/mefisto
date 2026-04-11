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
#include <cstdio>

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

    // WR-05: reorder so that every allocation that may throw happens BEFORE
    // we tear down the old painter/backing. If allocation fails mid-resize,
    // the old painter+backing stay intact and remain usable.
    QPixmap* old_backing = state_->backing_;

    QPixmap* new_backing = nullptr;
    try {
        new_backing = new QPixmap(device);
    } catch (...) {
        std::fprintf(stderr,
            "xvue-qt: resizeEvent: new QPixmap(%dx%d) failed; "
            "keeping previous backing\n",
            device.width(), device.height());
        return;  // old painter+backing untouched
    }
    new_backing->setDevicePixelRatio(dpr);

    // (d)(e) fill background, then blit old content top-left (DRAW-09,
    //       README_RESIZE.md invariants 2 and 3). One scoped temp painter.
    //       Do this BEFORE touching state_->painter_ so we can still bail.
    {
        QPainter tmp(new_backing);
        tmp.fillRect(new_backing->rect(), state_->background_);
        if (old_backing) {
            // drawPixmap at (0, 0) in LOGICAL coordinates — Qt handles DPR.
            tmp.drawPixmap(0, 0, *old_backing);
        }
    }  // ~QPainter tmp — end() called

    // Only now is it safe to tear down the old painter/backing:
    // (a) end old painter BEFORE deleting its device (Pitfall 7)
    if (state_->painter_ && state_->painter_->isActive()) {
        state_->painter_->end();
    }

    // (f) delete old, (g) swap
    delete old_backing;
    state_->backing_ = new_backing;

    // (h)(i) recreate/begin the long-lived painter on the new backing
    if (!state_->painter_) {
        try {
            state_->painter_ = new QPainter();
        } catch (...) {
            std::fprintf(stderr,
                "xvue-qt: resizeEvent: new QPainter() failed; "
                "display will freeze until next successful resize\n");
            return;
        }
    }
    if (!state_->painter_->begin(new_backing)) {
        // WR-05: loud diagnostic so a frozen display is not silent.
        std::fprintf(stderr,
            "xvue-qt: resizeEvent: QPainter::begin failed on %dx%d backing; "
            "display will freeze until next successful resize\n",
            device.width(), device.height());
        return;
    }

    // (j) DRAW-08 + Pitfall 5: hints do not carry across begin()/end()
    state_->painter_->setRenderHint(QPainter::Antialiasing, true);

    // (k) re-push pen+brush through applyPen() (D-22, Pitfall 5)
    state_->applyPen();
}
