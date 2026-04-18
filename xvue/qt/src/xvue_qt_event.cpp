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
#include <QTimer>

#include <cstdio>
#include <cstdlib>

XvueEventBridge::XvueEventBridge(XvueCanvas* canvas, QObject* parent)
    : QObject(parent), canvas_(canvas) {}

XvueEventBridge::~XvueEventBridge() = default;

// Plan 03: MEFISTO_XVSOURIS_DEBUG env-var cache.
//
// The env var is read once via getenv on first call and cached in a static
// local; subsequent calls are a simple bool read. Any non-empty, non-"0"
// value enables logging — "1", "true", "yes", etc. all count as on.
// Gated-off by default so production stderr is clean.
bool XvueEventBridge::debug_logging_enabled() {
    static const bool enabled = [] {
        const char* env = std::getenv("MEFISTO_XVSOURIS_DEBUG");
        if (env == nullptr)      return false;
        if (env[0] == '\0')      return false;
        if (env[0] == '0' && env[1] == '\0') return false;
        return true;
    }();
    return enabled;
}

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
    int         saved_motion_count = motion_count_;  // Plan 03
    int*        saved_items        = items_;
    int*        saved_pmin0        = pmin0_;

    // Reset per-call state (Pitfall 9: never trust stale quit_pending_).
    // Plan 03 also resets motion_count_ so the diagnostic counter reflects
    // only this invocation.
    mode_         = mode;
    pending_      = Result{};
    quit_pending_ = false;
    motion_count_ = 0;
    items_        = items;
    pmin0_        = pmin0;

    QEventLoop loop;
    loop_ = &loop;
    // Plan 03 layers motion coalescing via QTimer::singleShot(0, loop_,
    // &QEventLoop::quit) inside the MouseMove filter branch so waitForEvent
    // returns the *tail* of a motion burst — X11 XEventsQueued(QueuedAfterFlush)
    // parity. Button/key events still quit the loop synchronously.
    loop.exec();

    Result result = pending_;

    // Plan 03 diagnostic: when MEFISTO_XVSOURIS_DEBUG=1, log one line per
    // waitForEvent return with fields that let Plan 06 answer Assumption A2.
    if (debug_logging_enabled()) {
        std::fprintf(stderr,
                     "[xvsouris] mode=%d notypeevent=%d nbc=%d "
                     "x=%d y=%d motion_count=%d depth=%d\n",
                     static_cast<int>(mode),
                     result.notypeevent,
                     result.nbc,
                     result.x,
                     result.y,
                     motion_count_,
                     XvueApp::blockingDepth());
        std::fflush(stderr);
    }

    // Restore the outer call's state so its filter sees its own loop/pending
    // when the stack unwinds.
    loop_         = saved_loop;
    mode_         = saved_mode;
    pending_      = saved_pending;
    quit_pending_ = saved_quit_pending;
    motion_count_ = saved_motion_count;
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

    // Plan 03 (D-04): mouse-motion coalescing via deferred-quit timer.
    //
    // Reference semantics (xvuelc.c:2248-2263): the X11 body calls
    // XEventsQueued(QueuedAfterFlush) after each MotionNotify and only
    // returns the event if nothing else is queued behind it. The Qt
    // equivalent is: stash the (x, y) into pending_ on every MouseMove
    // (Pitfall 8: fresh every branch — never stale), arm a zero-delay
    // single-shot timer the first time (Pitfall 9: quit_pending_ was
    // already reset at the top of waitForEvent), and eat the event. The
    // timer enqueues a loop.quit() at the tail of the event queue, so any
    // motion events already queued ahead of it are dispatched (updating
    // pending_ as they go) before the timer fires. The result: loop.exec()
    // returns with the *last* coordinate pair in the burst, zero added
    // latency — identical to the X11 XEventsQueued(QueuedAfterFlush) path.
    case QEvent::MouseMove: {
        auto* me = static_cast<QMouseEvent*>(event);
        pending_.notypeevent = -2;
        pending_.nbc         = 0;      // X11 motion contract: no button carried
        pending_.x           = me->pos().x();  // Pitfall 8: fresh, not stale
        pending_.y           = me->pos().y();
        ++motion_count_;
        if (!quit_pending_) {
            quit_pending_ = true;
            QTimer::singleShot(0, loop_, &QEventLoop::quit);
        }
        return true;  // eat — we own the burst
    }

    default:
        return false;
    }
}
