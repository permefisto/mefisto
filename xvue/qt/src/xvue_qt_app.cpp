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
#include "xvue_qt_menu_bridge.h"   // Phase 6.0 Plan 06: menuBridge() forwarding
#include "xvue_qt_prefs.h"          // Phase 6.0 Plan 06: colorScheme() probe
#include <QApplication>
#include <QCoreApplication>
#include <QColor>
#include <QEventLoop>
#include <QFontDatabase>
#include <QGuiApplication>
#include <QPalette>
#include <QStringLiteral>
#include <QStyle>
#include <QStyleHints>
#include <Qt>
#include <cstdio>
#include <cstdlib>

std::once_flag                     XvueApp::once_flag_;
std::unique_ptr<QApplication>      XvueApp::qapp_;
std::unique_ptr<XvueWindow>        XvueApp::window_;
int                                XvueApp::font_id_ = -1;

// Phase 5 (D-03, EVENT-08). Static zero-init is load-bearing: any call to
// blockingDepth() before the first RAII guard must observe 0.
int                                XvueApp::blockingDepth_ = 0;

int XvueApp::blockingDepth() { return blockingDepth_; }

// Phase 6.1 Plan 03 Rule 3 auto-fix: force-link the per-module
// strong-override TUs. See xvue_qt_mail_actions.cpp for the full
// rationale. Without these keepalive references, GNU ld stops at the
// weak stub in xvue_qt_api.cpp.o and never pulls the mail/elas/flui/
// ther/nlse actions TUs from libxvueqt.a. 6.2 adds elas; 6.3..6.5
// will add flui/ther/nlse alongside in alphabetical order.
extern "C" int xvue_qt_elas_actions_keepalive();
static const int g_xvue_qt_elas_actions_keepalive_ref =
    xvue_qt_elas_actions_keepalive();
extern "C" int xvue_qt_mail_actions_keepalive();
static const int g_xvue_qt_mail_actions_keepalive_ref =
    xvue_qt_mail_actions_keepalive();

void XvueApp::load_bundled_font_()
{
    if (font_id_ >= 0) return;
    Q_INIT_RESOURCE(xvue_fonts);
    // Phase 6.1 Plan 03 Rule 3 auto-fix: initialize the shared icons
    // resource (see xvue_qt_icons.qrc). Plan 02 registered 10 mail SVG
    // icons inside the existing /xvue/qt/icons prefix but never called
    // Q_INIT_RESOURCE(xvue_icons), so the static-archive resource init
    // never ran for pp*_qt or the test binaries — qt.svg warnings
    // ("Cannot open file :/xvue/qt/icons/icons/mail/*.svg") followed.
    Q_INIT_RESOURCE(xvue_icons);
    // Silence unused-var warning for the force-link refs (g++ -Wall).
    (void)g_xvue_qt_elas_actions_keepalive_ref;
    (void)g_xvue_qt_mail_actions_keepalive_ref;
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
        // Phase 5 (D-05, Pitfall 7). Defensive — default-true on X11 per Qt
        // docs; setting explicitly documents intent and protects against any
        // future platform where the default flips. Must be set BEFORE the
        // QApplication constructor so the dispatcher picks it up during init.
        QCoreApplication::setAttribute(Qt::AA_CompressHighFrequencyEvents);
        qapp_ = std::make_unique<QApplication>(fake_argc, fake_argv);

        // Phase 6.0 Plan 06 (UX-13). Connect system theme-change signal so
        // when the user flips desktop dark-mode while pp*_qt is running, we
        // re-apply the palette (only if pref == "system"). Qt 6.5+ feature;
        // on DEs without a theme plugin (XFce bare), the signal may never
        // fire — fallback: user restarts app, applyColorSchemePreference on
        // startup syncs.
        if (auto* sh = QGuiApplication::styleHints()) {
            QObject::connect(sh, &QStyleHints::colorSchemeChanged,
                qapp_.get(), [](Qt::ColorScheme){
                    // Only re-apply if pref is "system" — otherwise the user's
                    // forced light/dark override stays sticky.
                    if (XvuePrefs::colorScheme() == QStringLiteral("system")) {
                        XvueApp::applyColorSchemePreference();
                    }
                });
        }

        std::atexit(&XvueApp::teardown_atexit);
    });
    load_bundled_font_();
}

// Phase 6.0 Plan 06 (UX-03/D-02). Forwarding accessor — null-safe so call
// sites inside extern "C" entries can use it before xvinitgraphique_ has
// constructed the window slot, or after xvfermer_ has reset it.
XvueMenuBridge* XvueApp::menuBridge() {
    auto& win = window_slot();
    if (!win) return nullptr;
    return win->menuBridge();
}

// Phase 6.0 Plan 06 (UX-13, D-05). Reads XvuePrefs::colorScheme() and applies
// the matching QPalette to QApplication. Three branches:
//   - "dark"   -> hand-crafted dark palette per RESEARCH §8 (standard Qt dark
//                 palette: Window=(53,53,53), white text, etc.).
//   - "light"  -> Qt default QPalette (the constructor builds the standard
//                 light palette).
//   - "system" -> QApplication::style()->standardPalette() — the value the
//                 platform style would apply on a fresh launch. On Linux,
//                 this is what the desktop theme plugin reports.
// Idempotent: calling this multiple times with the same pref leaves the
// palette unchanged. Safe to call before any window exists (only QApplication
// global palette is touched).
void XvueApp::applyColorSchemePreference() {
    const QString pref = XvuePrefs::colorScheme();
    if (pref == QStringLiteral("dark")) {
        // RESEARCH §8 standard Qt dark palette. Constants chosen to match
        // the most common third-party "Qt dark style" implementations so
        // existing third-party styles do not collide visually.
        QPalette p;
        p.setColor(QPalette::Window,          QColor(53, 53, 53));
        p.setColor(QPalette::WindowText,      Qt::white);
        p.setColor(QPalette::Base,            QColor(42, 42, 42));
        p.setColor(QPalette::AlternateBase,   QColor(66, 66, 66));
        p.setColor(QPalette::ToolTipBase,     Qt::white);
        p.setColor(QPalette::ToolTipText,     Qt::white);
        p.setColor(QPalette::Text,            Qt::white);
        p.setColor(QPalette::Button,          QColor(53, 53, 53));
        p.setColor(QPalette::ButtonText,      Qt::white);
        p.setColor(QPalette::BrightText,      Qt::red);
        p.setColor(QPalette::Highlight,       QColor(38, 79, 120));
        p.setColor(QPalette::HighlightedText, Qt::black);
        QApplication::setPalette(p);
    } else if (pref == QStringLiteral("light")) {
        QApplication::setPalette(QPalette());   // Qt default light palette
    } else {
        // "system" or anything else (Plan 02 prefs.cpp clamps unknown values
        // to "system" on read, so the else-branch is reached only for
        // "system").
        if (auto* style = QApplication::style()) {
            QApplication::setPalette(style->standardPalette());
        }
    }
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
