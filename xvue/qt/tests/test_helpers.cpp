// xvue/qt/tests/test_helpers.cpp
// Phase 5 Wave 0 (05-01 Task 1): minimal implementations; real canvas helpers
// land in Plan 02.
#include "test_helpers.h"
#include "xvue_qt_app.h"
#include "xvue_qt_window.h"
#include "xvue_qt_canvas.h"
#include <QCoreApplication>
#include <QEventLoop>
#include <QTimer>

namespace xvue_test {

void pumpEvents(int ms) {
    QEventLoop loop;
    QTimer::singleShot(ms, &loop, &QEventLoop::quit);
    loop.exec();
}

XvueCanvas* createTestCanvas() {
    XvueApp::ensure();
    // Placeholder: opening a real window requires the extern "C" Fortran ABI
    // path (Plan 02+). Wave 0 tests QSKIP anything that would need this.
    return nullptr;
}

void destroyTestCanvas() {
    // Placeholder — real teardown will happen in Plan 02+.
}

} // namespace xvue_test
