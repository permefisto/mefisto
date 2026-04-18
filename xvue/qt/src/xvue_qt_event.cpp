// xvue/qt/src/xvue_qt_event.cpp
// Phase 5 Plan 02 Task 2. Synchronous event-filter dispatch for button/key
// events. Motion handling lands in Plan 03. xvsouris2_ accrochage lands in
// Plan 05.
#include "xvue_qt_event.h"
#include "xvue_qt_canvas.h"
#include "xvue_qt_app.h"

#include <QEvent>
#include <QEventLoop>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QCoreApplication>
#include <QCursor>
#include <QThread>

XvueEventBridge::XvueEventBridge(XvueCanvas* canvas, QObject* parent)
    : QObject(parent), canvas_(canvas) {}

XvueEventBridge::~XvueEventBridge() = default;

int XvueEventBridge::translateKey(QKeyEvent* ev) {
    if (!ev) return 0;
    // Try text() first (handles AltGr + composed chars on AZERTY).
    const QByteArray bytes = ev->text().toLatin1();
    if (!bytes.isEmpty()) {
        const unsigned char c = static_cast<unsigned char>(bytes.at(0));
        if (c != 0) return static_cast<int>(c);
    }
    // Fallback switch for control keys whose text() is empty in some layouts.
    switch (ev->key()) {
    case Qt::Key_Escape:    return 27;
    case Qt::Key_Return:
    case Qt::Key_Enter:     return 13;
    case Qt::Key_Tab:       return 9;
    case Qt::Key_At:        return 64;   // D-06 belt-and-braces
    case Qt::Key_Backspace: return 8;    // saclav.f:286
    default:                return 0;    // arrows, F-keys, modifiers — dropped
    }
}

XvueEventBridge::Result
XvueEventBridge::waitForEvent(WaitMode mode, int* items, int* pmin0) {
    Q_ASSERT(QThread::currentThread() == QCoreApplication::instance()->thread());

    BlockingDepthGuard depth_guard;  // Pitfall 6: first local, RAII

    // Save-restore every member the filter reads so nested waitForEvent()
    // calls do not clobber the outer call's state. This is the single-bridge
    // equivalent of the per-call QEventLoop stack — the filter only ever
    // sees the innermost call's loop/mode/pending via these scalars.
    QEventLoop* saved_loop         = loop_;
    WaitMode    saved_mode         = mode_;
    Result      saved_pending      = pending_;
    bool        saved_quit_pending = quit_pending_;
    int*        saved_items        = items_;
    int*        saved_pmin0        = pmin0_;

    // Reset per-call state (Pitfall 9: never trust stale quit_pending_).
    mode_         = mode;
    pending_      = Result{};
    quit_pending_ = false;
    items_        = items;
    pmin0_        = pmin0;

    QEventLoop loop;
    loop_ = &loop;
    // Plan 03 will layer motion coalescing on top of this. For Plan 02 we
    // simply run the loop until a non-motion event sets pending_ and calls
    // loop.quit() directly.
    loop.exec();

    Result result = pending_;

    // Restore the outer call's state so its filter sees its own loop/pending
    // when the stack unwinds.
    loop_         = saved_loop;
    mode_         = saved_mode;
    pending_      = saved_pending;
    quit_pending_ = saved_quit_pending;
    items_        = saved_items;
    pmin0_        = saved_pmin0;

    return result;
}

bool XvueEventBridge::eventFilter(QObject* watched, QEvent* event) {
    (void)watched;
    if (!loop_) return false;  // pass-through when not blocking

    switch (event->type()) {

    case QEvent::MouseButtonPress: {
        auto* me = static_cast<QMouseEvent*>(event);
        // Plan 02: press-only path. Plan 03 refines to full-click detection.
        pending_.notypeevent = -1;
        switch (me->button()) {
        case Qt::LeftButton:   pending_.nbc = 1; break;
        case Qt::MiddleButton: pending_.nbc = 2; break;  // X11 btn2 = abort
        case Qt::RightButton:  pending_.nbc = 3; break;
        default:               pending_.nbc = 0; break;
        }
        pending_.x = me->pos().x();
        pending_.y = me->pos().y();
        // X11 parity: middle-button historically aborts (see xvuelc.c:2272);
        // kept for full parity even though D-06 focuses on Esc/@.
        if (me->button() == Qt::MiddleButton) {
            pending_.notypeevent = 0;
            pending_.nbc         = 2;
        }
        loop_->quit();
        return true;   // eat: Fortran caller owns this event
    }

    case QEvent::MouseButtonRelease: {
        auto* me = static_cast<QMouseEvent*>(event);
        // Full-click path: after a ButtonPress we set notypeevent=1 on
        // release. In Plan 02 we simply emit release alone (Plan 03 will
        // sequence press->release into a single click when appropriate).
        pending_.notypeevent = 1;
        switch (me->button()) {
        case Qt::LeftButton:   pending_.nbc = 1; break;
        case Qt::MiddleButton: pending_.nbc = 2; break;
        case Qt::RightButton:  pending_.nbc = 3; break;
        default:               pending_.nbc = 0; break;
        }
        pending_.x = me->pos().x();
        pending_.y = me->pos().y();
        loop_->quit();
        return true;
    }

    case QEvent::KeyPress: {
        auto* ke = static_cast<QKeyEvent*>(event);
        const int code = translateKey(ke);
        // D-06: Esc (27) and @ (64) are abort.
        if (code == 27 || code == 64) {
            pending_.notypeevent = 0;
            pending_.nbc         = code;
        } else if (code != 0) {
            pending_.notypeevent = 2;
            pending_.nbc         = code;
        } else {
            // Unhandled key (F-key, arrow, modifier alone) — don't quit.
            return true;  // eat but keep waiting
        }
        // X/Y: last known canvas-local cursor position, or 0 for
        // xvpause_ (Plan 04 will pass pending_x/y through deplsouris_).
        pending_.x = canvas_ ? canvas_->mapFromGlobal(QCursor::pos()).x() : 0;
        pending_.y = canvas_ ? canvas_->mapFromGlobal(QCursor::pos()).y() : 0;
        loop_->quit();
        return true;
    }

    // Motion: Plan 03.
    case QEvent::MouseMove:
        return false;  // pass-through until Plan 03 implements coalescing

    default:
        return false;
    }
}
