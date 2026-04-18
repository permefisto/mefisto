// xvue/qt/src/xvue_qt_event.cpp
// Phase 5 — real bodies land in Plan 02 Task 2.
#include "xvue_qt_event.h"

XvueEventBridge::XvueEventBridge(XvueCanvas* canvas, QObject* parent)
    : QObject(parent), canvas_(canvas) {}

XvueEventBridge::~XvueEventBridge() = default;

XvueEventBridge::Result XvueEventBridge::waitForEvent(WaitMode, int*, int*) {
    // Plan 02 Task 2 replaces this with the real nested QEventLoop body.
    BlockingDepthGuard depth_guard;  // still exercises the guard
    Result r;
    return r;
}

bool XvueEventBridge::eventFilter(QObject*, QEvent*) {
    return false;  // pass-through until Task 2
}

int XvueEventBridge::translateKey(QKeyEvent*) { return 0; }
