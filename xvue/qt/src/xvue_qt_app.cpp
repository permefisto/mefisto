// xvue/qt/src/xvue_qt_app.cpp
// Phase 1 (D-01, D-07..D-09, D-18): XvueApp implementation.
// Source: 01-RESEARCH.md §Pattern 1; 01-CONTEXT.md D-01..D-09.
//
// 2026-04-11 debug session phase-01-xvtest0-teardown-segfault:
// D-08 as originally written ("destroy QApplication via unique_ptr reset at
// atexit") is revised. Empirically on Linux/gfortran/Qt 6, destroying the
// QApplication inside an atexit handler that runs alongside libgfortran's
// own atexit chain races Qt's internal static teardown and crashes in
// __run_exit_handlers. The canonical embedded-Qt idiom (PyQt, Qt plugin
// hosts, etc.) is to LEAK QApplication: construct once, never destroy. The
// OS reclaims memory at process exit. See Pitfall 5 and Research A4:
// "construct once, never destroy until atexit" — the "construct once" half
// is load-bearing; the "destroy at atexit" half is not.
#include "xvue_qt_app.h"
#include "xvue_qt_window.h"
#include <QApplication>
#include <QCoreApplication>
#include <QEventLoop>
#include <QFontDatabase>
#include <QStringLiteral>
#include <cstdio>
#include <cstdlib>

std::once_flag                     XvueApp::once_flag_;
std::unique_ptr<QApplication>      XvueApp::qapp_;
std::unique_ptr<XvueWindow>        XvueApp::window_;
int                                XvueApp::font_id_ = -1;

void XvueApp::load_bundled_font_()
{
    if (font_id_ >= 0) return;
    Q_INIT_RESOURCE(xvue_fonts);
    font_id_ = QFontDatabase::addApplicationFont(
        QStringLiteral(":/xvue/qt/fonts/DejaVuSansMono.ttf"));
    if (font_id_ < 0) {
        std::fprintf(stderr,
            "xvue-qt: WARNING — bundled DejaVuSansMono.ttf failed to "
            "load from Qt resource path; falling back to platform "
            "default monospace font.\n");
    }
}

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
    load_bundled_font_();
}

QApplication* XvueApp::qapp() {
    return qapp_.get();
}

std::unique_ptr<XvueWindow>& XvueApp::window_slot() {
    return window_;
}

void XvueApp::teardown_atexit() {
    // D-08 (revised 2026-04-11): the single, process-exit hook. We must NOT
    // destroy QApplication here — doing so races libgfortran's atexit chain
    // and Qt's internal static teardown, crashing in __run_exit_handlers.
    //
    // What we DO: drain the deferred-delete queue so no orphaned events
    // remain, then destroy our window_slot() while Qt is still in a
    // well-defined state, then drain once more. Finally, deliberately
    // release() — never reset() — the unique_ptr so QApplication outlives
    // this function and leaks until the OS reclaims process memory. This is
    // the canonical "embed Qt in a non-Qt main" idiom.
    //
    // Intentionally: no qapp_->quit() — exec() was never called, so quit()
    // is a no-op and only serves to confuse readers.
    if (qapp_) {
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    window_.reset();
    if (qapp_) {
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    // Documented leak. Do not change to qapp_.reset() — see file header.
    (void)qapp_.release();
}
