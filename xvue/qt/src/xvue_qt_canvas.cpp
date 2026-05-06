// xvue/qt/src/xvue_qt_canvas.cpp
// Phase 2 (D-06): paintEvent = drawPixmap blit ONLY. Never mutate state here.
// Phase 2 (D-04, D-07, D-08): resizeEvent reallocates backing, preserves
// content top-left, and recreates the long-lived QPainter per DRAW-01.
// Phase 6.0 Plan 05 (UX-12, UX-13, UI-SPEC §Canvas, §Copywriting Flag #3):
// extends paintEvent with a setTransform(view_transform_) before the single
// drawPixmap (DRAW-01 Phase 6 documented extension — RESEARCH §7), plus an
// optional centered empty-state hint when state_->has_user_content_ is false.
// Adds 5 new event overrides (wheel zoom, middle-drag pan, context menu) and
// the resetView() slot. The view transform is composited at paint time only;
// the backing pixmap and Fortran draw calls remain in canonical pixel
// coordinates ("Fortran must not notice" invariant at the UX layer).
//
// MMB pan vs Phase 5 bridge-abort interaction (RESEARCH §7 Pitfall 7.2):
// Qt event order is eventFilter(canvas, event) -> canvas::eventHandler(event).
// Phase 5's XvueEventBridge::eventFilter intercepts MiddleButton in Souris
// mode as abort (notypeevent=0, nbc=2) and returns true to eat the event.
// Therefore middle-drag pan only runs when no blocking read is active
// (loop_ == nullptr). v1 behavior:
//   - User idle (no picking loop) -> middle-drag pans (normal desktop UX).
//   - User mid-picking (xvsouris_ active) -> middle-click aborts (Phase 5).
// Plan 03's A6 audit will confirm X11 parity. If X11 also pans during picking,
// Plan 03 may file a revision to require Ctrl+MMB instead.
#include "xvue_qt_canvas.h"
#include "xvue_qt_state.h"
#include "xvue_qt_app.h"
#include "xvue_qt_window.h"
// xvue_qt_menu_bridge.h: pulled in by future Plan 02-driven update of this
// file when XvueWindow::menuBridge() exists; currently the contextMenuEvent
// branch only references the window pointer.
#include "xvue_qt_i18n.h"
#include "xvue_qt_postscript.h"   // Phase 7 Plan 03: setCanvasDims hook

#include <QPainter>
#include <QPaintEvent>
#include <QResizeEvent>
#include <QWheelEvent>
#include <QMouseEvent>
#include <QContextMenuEvent>
#include <QPixmap>
#include <QMenu>
#include <QPalette>
#include <QFontMetrics>
#include <QString>

#include <cmath>
#include <cstdio>

XvueCanvas::XvueCanvas(XvueState* state, QWidget* parent)
    : QWidget(parent), state_(state)
{
    // Phase 3 D-20: dark-mode-freeze defensive guards. Remove the canvas
    // from Qt's style-system auto-fill path so even a future stylesheet in
    // Phase 6 cannot recolor the backing pixmap. WA_OpaquePaintEvent tells
    // Qt the widget will paint every pixel itself on every paintEvent.
    setAutoFillBackground(false);
    setAttribute(Qt::WA_OpaquePaintEvent, true);

    // Phase 5 (Pitfall 4). Required for QKeyEvent delivery to the canvas;
    // without StrongFocus, QTest::keyClick and real keyboard events bypass
    // the event filter the Plan 02 bridge will install here.
    setFocusPolicy(Qt::StrongFocus);

    // Phase 5 Plan 03: mouse tracking MUST be enabled so that MouseMove
    // events are delivered to the canvas when no mouse button is held.
    // The X11 xvsouris_ body reports motion with nbc=0 (no button) — Qt's
    // default widget behaviour is to only deliver MouseMove while a button
    // is pressed, which would break the xvsouris_ motion contract (and the
    // rubber-band zoom UI). This is the Qt equivalent of the X11 mesher
    // window having PointerMotionMask in its event mask (xvuelc.c opens the
    // window with pointer motion events always on).
    setMouseTracking(true);

    // D-04: DO NOT allocate backing_ here. Qt 6 guarantees resizeEvent fires
    // before the first paintEvent on X11/Wayland (Pitfall 6). The first
    // resizeEvent after construction performs the initial allocation.
}

XvueCanvas::~XvueCanvas() = default;

void XvueCanvas::paintEvent(QPaintEvent* event) {
    Q_UNUSED(event);
    // Phase 6.0 Plan 05 (UX-12) — DRAW-01 Phase 6 documented extension:
    // paintEvent body now does setTransform(view_transform_) + drawPixmap +
    // optional empty-state drawText. Still a single drawPixmap of the backing.
    // backing_ is never mutated by the transform — Fortran-side coordinates
    // remain canonical pixel coordinates.
    QPainter p(this);
    if (!state_) return;

    // Phase 9.1 fix 2026-05-06: maintainer reports blurry / doubled text on
    // canvas. Default QPainter pixmap transform uses smooth (bilinear)
    // interpolation when view_transform_ is anything other than identity
    // (HiDPI device-pixel-ratio, user zoom, fractional scale). Smooth
    // interpolation on text pixels in the backing pixmap creates ghost
    // glyphs that read as doubled / blurred. Disable so nearest-neighbour
    // scaling preserves crisp pixel-aligned rendering.
    p.setRenderHint(QPainter::SmoothPixmapTransform, false);
    p.setRenderHint(QPainter::Antialiasing, false);

    if (state_->backing_) {
        p.setTransform(state_->view_transform_);
        p.drawPixmap(0, 0, *state_->backing_);
    }

    // Reset transform for the empty-state text overlay so the hint is
    // centered in widget coords regardless of view_transform_.
    p.resetTransform();

    // Phase 6.0 Plan 05 + UI-SPEC Flag #3 (Empty-state hint).
    // Default: false at startup -> render hint until first meaningful draw.
    // Plan 06 will flip has_user_content_ in xvue_module_init_ once a module
    // owns the canvas. The Plan 01 i18n scaffold returns "" for these MsgIds;
    // Plan 02 fills the FR/EN copy.
    lastPaintDrewEmptyState_ = false;
    if (!state_->has_user_content_) {
        p.setPen(palette().color(QPalette::WindowText));
        const QRect r = rect();
        QFontMetrics fm(p.font());
        const QString heading = xvueT(MsgId::EmptyStateHeading);
        const QString body    = xvueT(MsgId::EmptyStateBody);
        const int lineH = fm.height();
        const int headW = fm.horizontalAdvance(heading);
        const int bodyW = fm.horizontalAdvance(body);
        const int headX = (r.width() - headW) / 2;
        const int bodyX = (r.width() - bodyW) / 2;
        const int centerY = r.height() / 2;
        // Always call drawText so test observables fire even with empty
        // strings (Plan 01 scaffold). Plan 02 fills real text.
        p.drawText(headX, centerY - lineH / 2, heading);
        p.drawText(bodyX, centerY + lineH,    body);
        lastPaintDrewEmptyState_ = true;
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
    state_->painter_->setRenderHint(QPainter::TextAntialiasing, true);  // D-07

    // (k) re-push pen+brush through applyPen() (D-22, Pitfall 5)
    state_->applyPen();

    // Phase 5 (Pitfall 10). The saved accrochage undo tile may point at a
    // location outside the new backing extents; invalidate it so the next
    // xvsouris2_ motion allocates a fresh 13x13 tile from the current
    // position. mempxaccro_ (the sprite template) is resolution-independent
    // and is NOT touched here.
    if (state_->accroche_undo_tile_) {
        delete state_->accroche_undo_tile_;
        state_->accroche_undo_tile_ = nullptr;
    }

    // Phase 7 Plan 03 (EXPORT-04): teach PsEmitter the current canvas
    // dimensions so its pyFlip() helper can compute ypixels-y at PS-emit
    // time. Logical pixels (NOT device pixels) match what every Qt-side
    // primitive passes through the painter (HiDPI scaling is owned by
    // the QPainter, not the Fortran caller).
    XvueApp::psEmitter().setCanvasDims(logical.width(), logical.height());
}

// ---------------------------------------------------------------------------
// Phase 6.0 Plan 05 — Qt-native gesture handlers below.
// ---------------------------------------------------------------------------

void XvueCanvas::wheelEvent(QWheelEvent* event) {
    // UX-12 + RESEARCH §7. Each wheel notch on a standard mouse delivers
    // angleDelta().y() = +/- 120. Apply 1.15x scale per notch (matches the
    // gesture density users expect from QGraphicsView wheel-zoom samples).
    if (!state_) {
        event->ignore();
        return;
    }
    const int notches = event->angleDelta().y() / 120;
    if (notches == 0) {
        event->ignore();
        return;
    }
    const qreal factor = std::pow(1.15, notches);
    QTransform t = state_->view_transform_ * QTransform().scale(factor, factor);
    // Clamp via m11 (and m22 by symmetry — equal scale factors). T-06.0-05-01
    // mitigation: refuse a zoom step that would push beyond [0.1, 10.0].
    const qreal s = t.m11();
    if (s < 0.1 || s > 10.0) {
        event->accept();
        return;   // refuse zoom step that would exceed bounds
    }
    state_->view_transform_ = t;
    update();
    event->accept();
}

void XvueCanvas::mousePressEvent(QMouseEvent* event) {
    if (event->button() == Qt::MiddleButton && state_) {
        // UX-12: middle-drag pan. We DO NOT forward to the bridge — the
        // bridge's eventFilter has already seen this press (Qt event order
        // is eventFilter -> canvas handler). When loop_ is non-null (Souris
        // mode active), the filter ate it and this handler is never reached.
        // When loop_ is null (idle), the filter passed it through and we own
        // the pan gesture.
        pan_start_            = event->pos();
        pan_anchor_transform_ = state_->view_transform_;
        pan_active_           = true;
        event->accept();
        return;
    }
    QWidget::mousePressEvent(event);
}

void XvueCanvas::mouseMoveEvent(QMouseEvent* event) {
    if (pan_active_ && (event->buttons() & Qt::MiddleButton) && state_) {
        const QPoint delta = event->pos() - pan_start_;
        // Apply translate AFTER the anchor transform so the drag delta is
        // in widget coordinates (a 5-pixel right-drag moves the view 5 px
        // right regardless of the current zoom). Math:
        //   final = anchor * translate(delta)
        // QTransform composes in column-vector order, so this lays the
        // translate "on top" of the anchored transform.
        state_->view_transform_ = pan_anchor_transform_ *
                                  QTransform().translate(delta.x(), delta.y());
        update();
        event->accept();
        return;
    }
    QWidget::mouseMoveEvent(event);
}

void XvueCanvas::mouseReleaseEvent(QMouseEvent* event) {
    if (event->button() == Qt::MiddleButton && pan_active_) {
        pan_active_ = false;
        event->accept();
        return;
    }
    QWidget::mouseReleaseEvent(event);
}

void XvueCanvas::contextMenuEvent(QContextMenuEvent* event) {
    // UI-SPEC §Canvas interaction: right-click opens a module-populated menu.
    // D-08 (Phase 5/6 modal re-entrancy guard): suppressed when blockingDepth
    // > 0 — the Fortran caller is mid-blocking-read and a menu would be a
    // re-entrancy hazard. The QAction itself is NOT disabled (per D-08 design
    // — keep it discoverable); the menu just refuses silently here.
    if (XvueApp::blockingDepth() > 0) {
        event->accept();
        return;
    }
    // Increment counter even when the menu is empty so test observability
    // captures the "menu would have been shown" event independently of the
    // populator state. Plan 07 manual sweep verifies the actual menu UI.
    ++contextMenuShownCount_;

    // Try to reach the menu bridge through the window's accessor. May be
    // nullptr in v1 (no module has populated the context menu populator yet)
    // — in that case the menu stays empty and we exit without exec().
    QMenu menu(this);
    if (auto& win = XvueApp::window_slot()) {
        // XvueWindow does not currently expose menuBridge() — Plan 02 will
        // add the accessor + setter. Until then, skip the populate step
        // (the menu stays empty and the early-return below handles it).
        // The contextMenuShownCount_ counter has already been bumped so
        // tests can verify the not-suppressed-by-blocking branch.
        (void)win;
    }
    if (!menu.isEmpty()) {
        menu.exec(event->globalPos());
    }
    event->accept();
}

void XvueCanvas::resetView() {
    // UX-12: View -> Fit to Window (Ctrl+0 in Plan 06). Restore identity and
    // schedule a repaint. Defensive against an extremely unlikely race where
    // resetView is called after canvas destruction (T-06.0-LAMBDA-01: Plan 06
    // should use QPointer<XvueCanvas> in any lambda that might outlive the
    // canvas); the state_ null-guard catches the obvious case.
    if (!state_) return;
    state_->view_transform_ = QTransform();
    update();
}
