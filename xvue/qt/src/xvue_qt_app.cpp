// xvue/qt/src/xvue_qt_app.cpp
// Phase 1 (D-01, D-07..D-09, D-18): XvueApp implementation.
// Source: 01-RESEARCH.md §Pattern 1; 01-CONTEXT.md D-01..D-09.
#include "xvue_qt_app.h"
#include "xvue_qt_window.h"
#include <QApplication>
#include <cstdlib>

std::once_flag                     XvueApp::once_flag_;
std::unique_ptr<QApplication>      XvueApp::qapp_;
std::unique_ptr<XvueWindow>        XvueApp::window_;

void XvueApp::ensure() {
    std::call_once(once_flag_, []{
        // D-01 / Pitfall 5: static storage so Qt can safely cache argv for the
        // full process lifetime. Qt 6 keeps pointers into argv alive across
        // reopens; heap-allocating these is a foot-gun.
        static int   fake_argc = 1;
        static char  arg0[] = "mefisto";
        static char* fake_argv[] = { arg0, nullptr };
        qapp_ = std::make_unique<QApplication>(fake_argc, fake_argv);
        std::atexit(&XvueApp::teardown_atexit);
    });
}

QApplication* XvueApp::qapp() {
    return qapp_.get();
}

std::unique_ptr<XvueWindow>& XvueApp::window_slot() {
    return window_;
}

void XvueApp::teardown_atexit() {
    // D-08: the single, process-exit teardown site. Destructor-time teardown
    // is unsafe (Pitfall 5 — Qt internal static-destruction order).
    if (qapp_) {
        qapp_->quit();
    }
    window_.reset();
    qapp_.reset();
}
