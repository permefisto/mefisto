// xvue/qt/src/xvue_qt_event.h
// Phase 5 (D-01, D-02, D-03, D-07, D-08, Pitfall 2, Pitfall 6, Pitfall 9).
// XvueEventBridge: the QObject event filter that owns the stack-local
// QEventLoop used by xvsouris_ / xvsouris2_ / xvpause_ to block the Fortran
// caller until a filtered mouse/keyboard event arrives.
//
// Thread-affinity: main thread only (SHELL-07). No atomics.
// ABI: NOT exposed via extern "C" — the four Fortran entry points call the
//      bridge through XvueApp::window_slot()->bridge() internally.
#pragma once
#include <QObject>
#include "xvue_qt_app.h"   // friend struct BlockingDepthGuard

class QEventLoop;
class QKeyEvent;
class XvueCanvas;

// Phase 5 (D-03, Pitfall 6). RAII increments XvueApp::blockingDepth_ on
// construction, decrements on destruction even via exception. Placed as the
// first local in waitForEvent().
struct BlockingDepthGuard {
    BlockingDepthGuard() noexcept {
        ++XvueApp::blockingDepth_;
        // Phase 5 (Pitfall 2). Hard ceiling against pathological Fortran
        // re-entry. Phase 6 modal dialogs legitimately push depth to 2.
        // >= 5 is a bug, not a workload.
        Q_ASSERT(XvueApp::blockingDepth_ < 4);
    }
    ~BlockingDepthGuard() noexcept {
        --XvueApp::blockingDepth_;
    }
    BlockingDepthGuard(const BlockingDepthGuard&)            = delete;
    BlockingDepthGuard& operator=(const BlockingDepthGuard&) = delete;
};

class XvueEventBridge : public QObject {
    Q_OBJECT
public:
    enum class WaitMode { Souris, Souris2, Pause };

    // Mirrors the Fortran out-parameter contract:
    //   notypeevent: -2 motion / -1 press / 0 abort / 1 full click / 2 key
    //   nbc:         button number (1/2/3) or ASCII code
    //   x, y:        canvas-local coordinates (logical px)
    struct Result {
        int notypeevent = 0;
        int nbc         = 0;
        int x           = 0;
        int y           = 0;
    };

    explicit XvueEventBridge(XvueCanvas* canvas, QObject* parent = nullptr);
    ~XvueEventBridge() override;

    // Blocks on a stack-local QEventLoop until a matching event arrives.
    // Increments XvueApp::blockingDepth_ via BlockingDepthGuard.
    // For WaitMode::Souris2, items[] is the candidate point array and
    // pmin0 is the in/out nearest-item index (-2 = unset).
    Result waitForEvent(WaitMode mode,
                        int* items = nullptr,
                        int* pmin0 = nullptr);

protected:
    bool eventFilter(QObject* watched, QEvent* event) override;

private:
    // Key → nbc translation: hybrid QKeyEvent::text()[0] + fallback switch.
    // D-07 defaults: Esc→27, Ret→13, Tab→9, @→64, Backspace→8.
    static int translateKey(QKeyEvent* ev);

    // Plan 03: MEFISTO_XVSOURIS_DEBUG env-var cache. The lookup runs once on
    // first call (static-local initialization is thread-safe in C++17 but the
    // bridge is main-thread-only anyway per SHELL-07). When enabled, the
    // diagnostic counter is written to stderr at each waitForEvent return so
    // Plan 06 can resolve Assumption A2 (whether AA_CompressHighFrequencyEvents
    // compresses before or after the filter).
    static bool debug_logging_enabled();

    XvueCanvas* canvas_        = nullptr;
    QEventLoop* loop_          = nullptr;   // non-null only inside waitForEvent
    WaitMode    mode_          = WaitMode::Souris;
    Result      pending_;                   // stash before loop.quit()
    bool        quit_pending_  = false;     // deferred-quit timer armed flag (Plan 03)
    int         motion_count_  = 0;         // Plan 03: motion events seen this waitForEvent

    // Souris2 only (Plan 05).
    int*        items_         = nullptr;
    int*        pmin0_         = nullptr;
};
