// xvue/qt/tests/test_xvue_qt_window_chrome.cpp
// Phase 6.0 Plan 06: integration QTest suite covering Tasks 1+2+3 of the
// Wave 4 wiring plan.
//
// Coverage matrix:
//   Task 1 (XvueApp): testMenuBridgeAccessorNonNullAfterWindow,
//                     testDarkPaletteApplication, testLightPaletteApplication,
//                     testSystemPaletteFallback.
//   Task 2 (XvueWindow chrome): testWindowHasMenuBar, testWindowHasToolBar,
//                               testWindowHasStatusBar, testWindowHasConsoleDock,
//                               testMenuBridgeInstalledOnWindow,
//                               testCoordLabelUpdatesOnMouseCoords,
//                               testFileOpenGuardedByBlockingDepth,
//                               testFitToWindowResetsCanvas.
//   Task 3 (xvue_module_init_): testModuleInitMarksBridgeRegistered,
//                               testModuleInitFlipsHasUserContent.
//
// Custom main pattern (mirrors test_xvue_qt_dialogs_console.cpp): set
// QT_QPA_PLATFORM=offscreen, call XvueApp::ensure() to construct the
// QApplication BEFORE QTest::qExec so AA_CompressHighFrequencyEvents is set
// in time, then exec the test class.
//
// QSettings isolation (mirrors test_xvue_qt_i18n_menu_prefs.cpp): redirect
// QSettings::IniFormat path to QDir::tempPath() so the developer's
// ~/.config/LJLL/mefisto* are not scribbled on by the test run.
#include "test_helpers.h"
#include "xvue_qt_app.h"
#include "xvue_qt_canvas.h"
#include "xvue_qt_console_dock.h"
#include "xvue_qt_event.h"          // BlockingDepthGuard
#include "xvue_qt_menu_bridge.h"
#include "xvue_qt_prefs.h"
#include "xvue_qt_state.h"
#include "xvue_qt_window.h"

#include <QtTest/QtTest>
#include <QAction>
#include <QApplication>
#include <QCoreApplication>
#include <QDir>
#include <QDockWidget>
#include <QEventLoop>
#include <QLabel>
#include <QMenu>
#include <QMenuBar>
#include <QPalette>
#include <QSettings>
#include <QStatusBar>
#include <QToolBar>
#include <QTransform>

extern "C" void xvinitgraphique_(void);
extern "C" void xvfermer_(void);
extern "C" void xvue_module_init_(char*, int*);

class TestXvueQtWindowChrome : public QObject {
    Q_OBJECT

private slots:
    void initTestCase() {
        // QSettings isolation — write to QDir::tempPath() not ~/.config.
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope,
                           QDir::tempPath());
        QCoreApplication::setOrganizationName(QStringLiteral("LJLL"));
        QCoreApplication::setApplicationName(QStringLiteral("mefisto-test-p6"));
        QSettings().clear();
        // Initialize XvuePrefs so colorScheme()/recentProjects()/etc. read a
        // real backend. The "mail" group is a deterministic placeholder; the
        // window-geometry/state tests use the per-module key namespace.
        XvuePrefs::initialize("mail");
        // Stand up the production window so all subsequent tests share a
        // single XvueApp::window_slot(). Mirrors test_xvue_qt_event.cpp.
        xvinitgraphique_();
        auto& win = XvueApp::window_slot();
        QVERIFY(win != nullptr);
    }

    void cleanupTestCase() {
        xvfermer_();
    }

    // ---- Drain Qt deleteLater between tests so heap-allocated Qt children
    //      do not race the leaked-QApplication atexit (D-08).
    void cleanup() {
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }

    // ============================ Task 1: XvueApp ============================

    // testMenuBridgeAccessorNonNullAfterWindow
    void testMenuBridgeAccessorNonNullAfterWindow() {
        // After xvinitgraphique_, XvueApp::menuBridge() forwards to the
        // window's installed bridge.
        QVERIFY(XvueApp::menuBridge() != nullptr);
    }

    // testDarkPaletteApplication
    void testDarkPaletteApplication() {
        XvuePrefs::saveColorScheme(QStringLiteral("dark"));
        XvueApp::applyColorSchemePreference();
        const QColor w = QApplication::palette().color(QPalette::Window);
        QCOMPARE(w, QColor(53, 53, 53));
        // Reset for later tests.
        XvuePrefs::saveColorScheme(QStringLiteral("system"));
        XvueApp::applyColorSchemePreference();
    }

    // testLightPaletteApplication
    void testLightPaletteApplication() {
        XvuePrefs::saveColorScheme(QStringLiteral("light"));
        XvueApp::applyColorSchemePreference();
        const QColor w = QApplication::palette().color(QPalette::Window);
        const QColor defW = QPalette().color(QPalette::Window);
        QCOMPARE(w, defW);
    }

    // testSystemPaletteFallback
    void testSystemPaletteFallback() {
        XvuePrefs::saveColorScheme(QStringLiteral("system"));
        XvueApp::applyColorSchemePreference();
        // Just assert no crash and palette is valid.
        QVERIFY(QApplication::palette().color(QPalette::Window).isValid());
    }

    // ====================== Task 2: XvueWindow chrome ========================

    // testWindowHasMenuBar
    void testWindowHasMenuBar() {
        auto* win = XvueApp::window_slot().get();
        QVERIFY(win);
        QVERIFY(win->menuBar() != nullptr);
        QVERIFY(win->menuBar()->findChild<QMenu*>("File") != nullptr);
        QVERIFY(win->menuBar()->findChild<QMenu*>("View") != nullptr);
        QVERIFY(win->menuBar()->findChild<QMenu*>("Help") != nullptr);
    }

    // testWindowHasToolBar
    void testWindowHasToolBar() {
        auto* win = XvueApp::window_slot().get();
        QList<QToolBar*> tbs = win->findChildren<QToolBar*>();
        QVERIFY(tbs.size() >= 1);
        // 7 explicit actions + 2 separators = 9 actions on the toolbar.
        QVERIFY(tbs[0]->actions().size() >= 7);
    }

    // testWindowHasStatusBar
    void testWindowHasStatusBar() {
        auto* win = XvueApp::window_slot().get();
        QVERIFY(win->statusBar() != nullptr);
        // The permanent coordLabel_ widget exists as a child of the statusBar.
        QList<QLabel*> labels = win->statusBar()->findChildren<QLabel*>();
        QVERIFY(labels.size() >= 1);
    }

    // testWindowHasConsoleDock
    void testWindowHasConsoleDock() {
        auto* win = XvueApp::window_slot().get();
        QVERIFY(win->consoleDock() != nullptr);
    }

    // testMenuBridgeInstalledOnWindow
    void testMenuBridgeInstalledOnWindow() {
        auto* win = XvueApp::window_slot().get();
        QVERIFY(win->menuBridge() != nullptr);
    }

    // testCoordLabelUpdatesOnMouseCoords
    void testCoordLabelUpdatesOnMouseCoords() {
        auto* win = XvueApp::window_slot().get();
        // Invoke the slot directly via the meta-object so we exercise the
        // signal-slot path that XvueCanvas::mouseCoords would normally drive.
        // (Driving it via emit canvas->mouseCoords() would require a friend
        // declaration in XvueCanvas; the slot is the equivalent test surface.)
        QMetaObject::invokeMethod(win, "updateStatusCoords",
                                  Q_ARG(QPoint, QPoint(42, 84)));
        QList<QLabel*> labels = win->statusBar()->findChildren<QLabel*>();
        QVERIFY(!labels.isEmpty());
        const QString text = labels[0]->text();
        QVERIFY2(text.contains("42"),
                 qPrintable(QStringLiteral("status text: %1").arg(text)));
        QVERIFY2(text.contains("84"),
                 qPrintable(QStringLiteral("status text: %1").arg(text)));
    }

    // testFileOpenGuardedByBlockingDepth
    void testFileOpenGuardedByBlockingDepth() {
        auto* win = XvueApp::window_slot().get();
        // Inside this RAII scope, blockingDepth() > 0 -> File→Open's
        // refuseIfBlocking() fires and returns immediately. If the guard
        // worked, no QFileDialog ever opens (which would block the test).
        {
            BlockingDepthGuard g;
            QCOMPARE(XvueApp::blockingDepth(), 1);
            QMenu* fileMenu = win->menuBar()->findChild<QMenu*>("File");
            QVERIFY(fileMenu != nullptr);
            // First action is File→Open (per buildMenuBar order).
            QAction* openAct = fileMenu->actions().value(0);
            QVERIFY(openAct != nullptr);
            openAct->trigger();
            // If we reached this point without QFileDialog hijacking the
            // event loop, the guard is doing its job.
            QVERIFY(true);
        }
        QCOMPARE(XvueApp::blockingDepth(), 0);
    }

    // testFitToWindowResetsCanvas
    void testFitToWindowResetsCanvas() {
        auto* win = XvueApp::window_slot().get();
        // Set a non-identity transform, then trigger View → Fit (Ctrl+0).
        win->state()->view_transform_ = QTransform().scale(2.0, 2.0);
        QVERIFY(!win->state()->view_transform_.isIdentity());
        QMenu* viewMenu = win->menuBar()->findChild<QMenu*>("View");
        QVERIFY(viewMenu != nullptr);
        bool triggered = false;
        for (QAction* a : viewMenu->actions()) {
            // Match by shortcut so the test is locale-agnostic.
            if (a->shortcut() == QKeySequence("Ctrl+0")) {
                a->trigger();
                triggered = true;
                break;
            }
        }
        QVERIFY(triggered);
        QVERIFY(win->state()->view_transform_.isIdentity());
    }

    // ====================== Task 3: xvue_module_init_ ========================

    // testModuleInitMarksBridgeRegistered
    void testModuleInitMarksBridgeRegistered() {
        auto* win = XvueApp::window_slot().get();
        QVERIFY(win->menuBridge() != nullptr);
        char name[] = "mail";
        int  len    = 4;
        xvue_module_init_(name, &len);
        QVERIFY(win->menuBridge()->hasRegisteredModule());
    }

    // testModuleInitFlipsHasUserContent
    void testModuleInitFlipsHasUserContent() {
        auto* win = XvueApp::window_slot().get();
        win->state()->has_user_content_ = false;
        char name[] = "mail";
        int  len    = 4;
        xvue_module_init_(name, &len);
        QVERIFY(win->state()->has_user_content_);
    }
};

int main(int argc, char* argv[]) {
    qputenv("QT_QPA_PLATFORM", "offscreen");
    // CRITICAL: disable XvueConsoleDock::installStdoutRedirect for the QTest
    // run, otherwise xvue_module_init_'s dup2-over-STDOUT_FILENO swallows the
    // QtTest reporter's PASS/FAIL/Totals output and the suite appears to
    // truncate. Production pp*_qt invocations do not set this var.
    qputenv("XVUE_QT_DISABLE_STDOUT_REDIRECT", "1");
    XvueApp::ensure();
    TestXvueQtWindowChrome tc;
    int   qt_argc   = 1;
    char  qt_arg0[] = "xvue_qt_window_chrome_tests";
    char* qt_argv[] = { qt_arg0, nullptr };
    (void)argc;
    (void)argv;
    return QTest::qExec(&tc, qt_argc, qt_argv);
}

#include "test_xvue_qt_window_chrome.moc"
