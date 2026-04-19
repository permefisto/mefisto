// xvue/qt/src/xvue_qt_event.cpp
// Phase 5 Plan 02 Task 2. Synchronous event-filter dispatch for button/key
// events. Motion coalescing lands in Plan 03. xvsouris2_ accrochage lands
// in Plan 05 Task 2 (WaitMode::Souris2 branch below).
#include "xvue_qt_event.h"
#include "xvue_qt_canvas.h"
#include "xvue_qt_app.h"
#include "xvue_qt_menu_bridge.h"   // Phase 6.0 Plan 03: popChar pre-drain
#include "xvue_qt_window.h"
#include "xvue_qt_state.h"

#include <QEvent>
#include <QEventLoop>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QCoreApplication>
#include <QCursor>
#include <QPainter>
#include <QPixmap>
#include <QRect>
#include <QThread>
#include <QTimer>

#include <cstdio>
#include <cstdlib>
#include <limits>

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

// Phase 5 Plan 05 Task 2 (EVENT-03, Strategy B). File-local helpers for the
// accrochage save-restore-blit dance. All three reference the fixed 13x13
// sprite geometry from xvuelc.c:142-143 (lmempxaccro = hmempxaccro = 13)
// and treat (cx, cy) as the ITEM's canvas-local center. The sprite is
// blitted with its top-left corner at (cx - 6, cy - 6) so the black square
// border is centered on the item — byte-parity with the X11 body at
// xvuelc.c:2421-2434 (XCopyArea destination `items[pmin]-lmempxaccro/2`).
//
// Plan 05 invariant: every save_tile_under call is preceded either by a
// restore_tile on the PREVIOUS (*pmin0) location (motion-to-new-item
// case) or by a guard ensuring no tile was previously saved (*pmin0 < 0).
// This keeps accroche_undo_tile_ in lockstep with what is actually drawn
// on the backing, so the Rule-2/T-05-05-02 leak-on-repeated-motion threat
// is closed by construction.

// X11-parity nearest-item search. Layout verified against xvue/saclav.f:61
// (`CALL XVSOURIS2(MCN(MNIT), ...)`) and xvuelc.c:2397-2413:
//   items[0] = mots    (words per item, set by itemau.f:17)
//   items[1] = maxcap  (max items stockable, set by itemau.f:20)
//   items[2] = nbitem  (actual item count, set elsewhere in the mesher)
//   first item pair at items[mots], items[mots+1]
//   k-th item pair   at items[(k+1)*mots], items[(k+1)*mots+1]
//
// Returns the OFFSET `p` into items[] for the nearest (x, y) pair (mirrors
// the X11 semantic where *pmin0 holds the offset, NOT the index). -1 means
// "no item found" — empty array or nbitem==0.
//
// T-05-05-01 mitigation: clamp nbitem to a sane ceiling so a garbage
// items[2] cannot walk off the end of memory.
static int nearest_item_offset(const int* items, int cx, int cy) {
    if (!items) return -1;
    const int mots = items[0];
    if (mots < 2) return -1;   // need at least x,y per item
    int nbitem = items[2];
    if (nbitem <= 0) return -1;
    // T-05-05-01: clamp to 65536; saclav.f-driven counts are a few thousand.
    if (nbitem > 65536) nbitem = 65536;

    int dmin = std::numeric_limits<int>::max();
    int pmin = -1;
    int p    = 0;
    for (int k = 0; k < nbitem; ++k) {
        p += mots;
        const int dx = items[p]     - cx;
        const int dy = items[p + 1] - cy;
        const int d  = dx * dx + dy * dy;
        if (d < dmin) {
            dmin = d;
            pmin = p;
        }
    }
    return pmin;
}

// Save a 13x13 tile of the canvas backing centered on (cx, cy) into
// state->accroche_undo_tile_. Allocates the tile pixmap on first use and
// reuses it on subsequent saves (its size is invariant). Safe against null
// state/backing.
static void save_tile_under(XvueState* state, int cx, int cy) {
    if (!state || !state->backing_) return;
    if (!state->accroche_undo_tile_) {
        state->accroche_undo_tile_ = new QPixmap(13, 13);
    }
    // DPR must match the backing so the copyFromCanvas source aligns 1:1
    // with the sprite destination. Backing DPR can change when the window
    // moves between monitors (Phase 4 WR-02).
    state->accroche_undo_tile_->setDevicePixelRatio(
        state->backing_->devicePixelRatio());
    // The sprite's top-left on the backing is (cx - 6, cy - 6); copy that
    // 13x13 region out. QPainter::drawPixmap(target=(0,0), source=backing,
    // sourceRect=...) uses LOGICAL coordinates — Qt handles DPR scaling.
    QPainter p(state->accroche_undo_tile_);
    p.drawPixmap(QPoint(0, 0), *state->backing_,
                 QRect(cx - 6, cy - 6, 13, 13));
}

// Blit the saved undo tile back onto the canvas backing at (cx - 6, cy - 6),
// erasing the sprite. Uses the long-lived painter_ so the drawing invariants
// from Phase 2 D-05 (single active painter on backing) are preserved.
static void restore_tile(XvueState* state, int cx, int cy) {
    if (!state || !state->backing_ || !state->accroche_undo_tile_) return;
    if (!state->painter_ || !state->painter_->isActive())         return;
    state->painter_->drawPixmap(cx - 6, cy - 6, *state->accroche_undo_tile_);
}

// Blit the accrochage sprite onto the canvas backing at (cx - 6, cy - 6).
// Must be called AFTER save_tile_under so the prior content is preserved
// and can be restored on the next motion. No-op if mempxaccro_ is null
// (initaccrochage_ was never called — defensive; xvsouris2_'s guard catches
// this but the helper stays safe).
static void draw_sprite(XvueState* state, int cx, int cy) {
    if (!state || !state->backing_ || !state->mempxaccro_)        return;
    if (!state->painter_ || !state->painter_->isActive())         return;
    state->painter_->drawPixmap(cx - 6, cy - 6, *state->mempxaccro_);
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

    // Phase 6.0 Plan 03 (UX-03): menu-queue pre-drain. If a QAction handler
    // queued a synthetic lexicon character (XvueMenuBridge::queueLexicon),
    // return it as notypeevent=2 (KeyPress) BEFORE entering the nested
    // QEventLoop. Drains one byte per waitForEvent call — matches Fortran
    // SACLAV's one-char-per-read contract (xvue/saclav.f §270-333,
    // 06.0-RESEARCH §6).
    //
    // GUARD: Souris mode ONLY.
    //   - Souris2 (xvsouris2_) is the accrochage/picking path; the Fortran
    //     caller expects mouse/key events with item-tracking state — a
    //     synthetic menu char would corrupt the accrochage state machine.
    //   - Pause (xvpause_) blocks on a specific key-press gate per Phase 5
    //     D-06; menu chars must not bypass that semantic.
    //   In both alternate modes the queued char STAYS in the queue and is
    //   drained by the next plain xvsouris_ call (standard Souris mode).
    //
    // ORDERING (per Plan 03 must_haves):
    //   XVUE_QT_ASSERT_MAIN_THREAD (in xvsouris_ wrapper)
    //   -> BlockingDepthGuard ctor (above)
    //   -> save-restore + per-call-state reset (above)
    //   -> menu-queue pre-drain (THIS BLOCK)
    //   -> QEventLoop ctor (below)
    //   -> loop.exec() (below)
    //
    // NOTE: AUTOEXIT lives in the xvsouris_ wrapper (xvue_qt_api.cpp), NOT
    // here, so the wrapper's pre-drain runs BEFORE AUTOEXIT — see the
    // matching block at xvue_qt_api.cpp::xvsouris_. Keeping the pre-drain
    // ALSO here is defense in depth for direct waitForEvent() callers
    // (e.g., the Phase 5 test suite that uses XvueEventBridge directly
    // without going through xvsouris_).
    //
    // WARNING: The early-return restore of saved_* below DUPLICATES the
    // end-of-function restore block. Any future addition to per-call mutable
    // state (e.g., a new tracking field for Souris2) MUST be mirrored in
    // BOTH restore sites or this early-return path will leak stale state
    // into the outer (caller's) waitForEvent frame. The duplication is
    // localised; long-term refactor would extract a SavedWaitState helper.
    if (mode == WaitMode::Souris) {
        if (auto& win = XvueApp::window_slot()) {
            if (auto* mb = win->menuBridge()) {
                if (auto c = mb->popChar()) {
                    pending_.notypeevent = 2;
                    pending_.nbc         = static_cast<unsigned char>(*c);
                    pending_.x           = canvas_
                        ? canvas_->mapFromGlobal(QCursor::pos()).x() : 0;
                    pending_.y           = canvas_
                        ? canvas_->mapFromGlobal(QCursor::pos()).y() : 0;
                    if (debug_logging_enabled()) {
                        std::fprintf(stderr,
                            "[xvsouris] mode=%d notypeevent=2 nbc=%d "
                            "x=%d y=%d motion_count=0 depth=%d "
                            "[menu-pre-drain]\n",
                            static_cast<int>(mode),
                            pending_.nbc, pending_.x, pending_.y,
                            XvueApp::blockingDepth());
                        std::fflush(stderr);
                    }
                    Result r = pending_;
                    // ---- early-return restore (MIRROR of end-of-function) ----
                    loop_         = saved_loop;
                    mode_         = saved_mode;
                    pending_      = saved_pending;
                    quit_pending_ = saved_quit_pending;
                    motion_count_ = saved_motion_count;
                    items_        = saved_items;
                    pmin0_        = saved_pmin0;
                    return r;
                }
            }
        }
    }

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

// Plan 05 Task 2 (Strategy B). Accrochage cleanup helper: if a sprite was
// previously drawn (*pmin0_ >= 0), restore the tile under it and invalidate
// *pmin0_. Used on ButtonRelease and abort paths so the canvas is clean of
// residue when the Fortran caller resumes.
void XvueEventBridge::cleanupAccrochage() {
    if (!items_ || !pmin0_) return;
    if (*pmin0_ < 0)        return;
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* state = win->state();
    if (!state || !state->accroche_undo_tile_) return;
    const int old_cx = items_[*pmin0_];
    const int old_cy = items_[*pmin0_ + 1];
    restore_tile(state, old_cx, old_cy);
    *pmin0_ = -2;
    // Free the tile — next call's first motion will allocate a fresh 13x13.
    delete state->accroche_undo_tile_;
    state->accroche_undo_tile_ = nullptr;
    if (win->canvas()) win->canvas()->update();
}

bool XvueEventBridge::eventFilter(QObject* watched, QEvent* event) {
    (void)watched;
    if (!loop_) return false;  // pass-through when not blocking

    switch (event->type()) {

    case QEvent::MouseButtonPress: {
        auto* me = static_cast<QMouseEvent*>(event);
        // Plan 05 Task 2 (EVENT-03): Souris2 mode treats a press as the final
        // "item picked" signal. Mirrors xvuelc.c:2383-2439 where MotionNotify
        // and ButtonPress share the accrochage path and return notypeevent=5
        // (or 0 on middle-button abort).
        if (mode_ == WaitMode::Souris2) {
            int btn = 0;
            switch (me->button()) {
            case Qt::LeftButton:   btn = 1; break;
            case Qt::MiddleButton: btn = 2; break;
            case Qt::RightButton:  btn = 3; break;
            default:               btn = 0; break;
            }
            // Middle-button parity with xvsouris_: abort. saclav.f:83-86
            // remaps btn=2 to notypeevent=-1 / nbc=64 (@ abort path).
            if (me->button() == Qt::MiddleButton) {
                cleanupAccrochage();
                pending_.notypeevent = 0;
                pending_.nbc         = 2;
                pending_.x           = me->pos().x();
                pending_.y           = me->pos().y();
                loop_->quit();
                return true;
            }
            // Normal path: accrochage redraw + return notypeevent=5.
            const int cx = me->pos().x();
            const int cy = me->pos().y();
            auto& win = XvueApp::window_slot();
            auto* state = win ? win->state() : nullptr;
            const int new_offset = nearest_item_offset(items_, cx, cy);

            // Erase the previously-drawn sprite if we had one and the new
            // nearest item differs. Mirrors xvuelc.c:2415-2425.
            if (pmin0_ && state && state->accroche_undo_tile_ &&
                *pmin0_ >= 0 && *pmin0_ != new_offset) {
                const int old_cx = items_[*pmin0_];
                const int old_cy = items_[*pmin0_ + 1];
                restore_tile(state, old_cx, old_cy);
                *pmin0_ = -2;
            }

            // Draw the new sprite if an item is nearest.
            if (new_offset >= 0 && state) {
                const int icx = items_[new_offset];
                const int icy = items_[new_offset + 1];
                save_tile_under(state, icx, icy);
                draw_sprite(state, icx, icy);
                if (pmin0_) *pmin0_ = new_offset;
                if (win->canvas()) win->canvas()->update();
            }
            pending_.notypeevent = 5;
            pending_.nbc         = btn;
            pending_.x           = cx;
            pending_.y           = cy;
            loop_->quit();
            return true;
        }
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
        // Plan 05 Task 2: Souris2 mode erases the accrochage sprite on
        // release and returns notypeevent=1. Mirrors xvuelc.c:2442-2463.
        if (mode_ == WaitMode::Souris2) {
            cleanupAccrochage();
            int btn = 0;
            switch (me->button()) {
            case Qt::LeftButton:   btn = 1; break;
            case Qt::MiddleButton: btn = 2; break;
            case Qt::RightButton:  btn = 3; break;
            default:               btn = 0; break;
            }
            pending_.notypeevent = 1;
            pending_.nbc         = btn;
            pending_.x           = me->pos().x();
            pending_.y           = me->pos().y();
            loop_->quit();
            return true;
        }
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
            // Plan 05 Task 2: on abort in Souris2 mode we erase the
            // sprite first so the canvas is clean when the caller resumes.
            if (mode_ == WaitMode::Souris2) cleanupAccrochage();
            pending_.notypeevent = 0;
            pending_.nbc         = code;
        } else if (code != 0) {
            // Plan 05 Task 2: also erase in Souris2 mode on ordinary key
            // press. xvuelc.c:2465-2477 does not explicitly erase on
            // keypress but saclav.f calls xvsouris2_ in a loop, so any
            // intermediate press handler benefits from a clean canvas.
            if (mode_ == WaitMode::Souris2) cleanupAccrochage();
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
    //
    // Plan 05 Task 2: WaitMode::Souris2 overrides the motion branch. On
    // each MouseMove we run the accrochage logic (nearest-item, erase old
    // sprite, draw new sprite on the backing). Unlike the Souris branch
    // we STILL arm the deferred-quit timer so fast drags across the canvas
    // coalesce to a single return — the sprite-draw keeps the user's
    // visual feedback smooth without causing one waitForEvent return per
    // raw pixel of mouse movement.
    case QEvent::MouseMove: {
        auto* me = static_cast<QMouseEvent*>(event);
        const int cx = me->pos().x();
        const int cy = me->pos().y();

        // Phase 6.0 Plan 03 (UX-12): emit live coord signal BEFORE the bridge
        // consumes the event. Phase 5's filter eats MouseMove during blocking
        // reads — without this emit, Plan 06's status-bar coord readout would
        // freeze whenever a Fortran-side xvsouris_/xvsouris2_ is active.
        // Friendship between XvueCanvas and XvueEventBridge (declared in
        // xvue_qt_canvas.h) permits the cross-class emit without Qt 6.4+
        // "calling signal from outside class" warnings. The emit fires for
        // BOTH Souris and Souris2 modes (the canvas's coord readout is
        // mode-agnostic). Cost: one signal dispatch per move (negligible vs
        // motion event-dispatch overhead). Null-guarded against the
        // canvas-less defensive path.
        if (canvas_) {
            emit canvas_->mouseCoords(QPoint(cx, cy));
        }

        if (mode_ == WaitMode::Souris2) {
            auto& win = XvueApp::window_slot();
            auto* state = win ? win->state() : nullptr;
            const int new_offset = nearest_item_offset(items_, cx, cy);

            // Erase previous sprite if the nearest item changed.
            if (pmin0_ && state && state->accroche_undo_tile_ &&
                *pmin0_ >= 0 && *pmin0_ != new_offset) {
                const int old_cx = items_[*pmin0_];
                const int old_cy = items_[*pmin0_ + 1];
                restore_tile(state, old_cx, old_cy);
                *pmin0_ = -2;
            }

            // Draw the new sprite at the new item (if any and not already drawn).
            if (new_offset >= 0 && state &&
                (!pmin0_ || *pmin0_ != new_offset)) {
                const int icx = items_[new_offset];
                const int icy = items_[new_offset + 1];
                save_tile_under(state, icx, icy);
                draw_sprite(state, icx, icy);
                if (pmin0_) *pmin0_ = new_offset;
                if (win->canvas()) win->canvas()->update();
            }

            pending_.notypeevent = 5;
            pending_.nbc         = 0;  // motion carries no button
            pending_.x           = cx;
            pending_.y           = cy;
            ++motion_count_;
            if (!quit_pending_) {
                quit_pending_ = true;
                QTimer::singleShot(0, loop_, &QEventLoop::quit);
            }
            return true;
        }

        pending_.notypeevent = -2;
        pending_.nbc         = 0;      // X11 motion contract: no button carried
        pending_.x           = cx;      // Pitfall 8: fresh, not stale
        pending_.y           = cy;
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
