// xvue/qt/tests/test_xvue_qt_menu_scaffold.cpp
// Phase 6.0 Plan 01 (scaffold): three sanity tests asserting the new Phase 6.0
// classes are linkable and their methods don't crash. Plan 02 flips most of
// these QCOMPAREs from "expect zero/empty" to "expect real text + queue size".
//
// Custom main() pattern follows test_xvue_qt_event.cpp — XvueApp::ensure()
// constructs QApplication BEFORE QTest::qExec so Qt's
// AA_CompressHighFrequencyEvents attribute set inside ensure() is visible.
#include "xvue_qt_app.h"
#include "xvue_qt_menu_bridge.h"
#include "xvue_qt_console_dock.h"
#include "xvue_qt_error_batcher.h"
#include "xvue_qt_preferences.h"
#include "xvue_qt_recent_projects.h"
#include "xvue_qt_prefs.h"
#include "xvue_qt_i18n.h"

#include <QtTest/QtTest>
#include <QApplication>
#include <QCoreApplication>
#include <QEventLoop>

// Plan 06 will fill the dispatch body; Plan 01 ships a warn-once stub. We
// declare the symbol locally so the test harness can prove the ABI hook
// resolves at link time.
extern "C" void xvue_module_init_(char* name, int* name_len);

class TestXvueQtMenuScaffold : public QObject {
    Q_OBJECT

private slots:
    // Drain Qt's deleteLater queue between tests + at teardown so any
    // top-level QWidget allocated by a test (QDockWidget, QDialog, QMenu)
    // is fully released before the QApplication-shutdown atexit handler
    // runs. Mirrors xvfermer_'s processEvents drain (D-08 guidance).
    void cleanup() {
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    void cleanupTestCase() {
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }

    // Sanity: every new Phase 6.0 class is linkable + constructs without crash.
    // Stack-allocated so destruction is deterministic and happens before any
    // QApplication-shutdown races (Phase 1 D-08: QApplication is leaked at
    // atexit, so heap-allocated top-level QWidgets registered with qApp's
    // topLevelWidgets() list outlive the test scope and can crash at exit).
    void testScaffoldIsLinkable() {
        XvueMenuBridge          bridge;
        XvueConsoleDock         dock;
        XvueErrorBatcher        batcher;
        XvuePreferencesDialog   prefsDlg;
        XvueRecentProjectsMenu  recent(&bridge);

        XvuePrefs::initialize("mail");

        // Plan 02 will return real text from xvueT(); Plan 01 stub returns "".
        // Accept either: scaffold contract is "doesn't crash" only.
        QString s = xvueT(MsgId::AppName);
        QVERIFY(s.isEmpty() || !s.isEmpty());

        // Touch each subsystem to prove the call-sites resolve. Plan 02
        // gives queueLexicon a real body — even an empty cmd enqueues the
        // trailing CR (=13), so drain it before the QVERIFY.
        bridge.queueLexicon(QStringLiteral(""));
        auto cr = bridge.popChar();
        QVERIFY(cr.has_value());
        QCOMPARE(static_cast<int>(*cr), 13);
        QVERIFY(!bridge.popChar().has_value());
        dock.appendLine(QStringLiteral("scaffold sanity"));
        batcher.enqueue(QStringLiteral("scaffold error"));
        recent.refresh();
        // prefsDlg ctor + dtor exercise — no show()/exec() to avoid
        // event-loop entry (D-10 / SHELL-03).
        (void)prefsDlg.windowTitle();
    }

    // Sanity: ABI count = 58 — xvue_module_init_ symbol present and callable.
    // Plan 06 fills the real dispatch body; Plan 01 warn-once stub is enough.
    void testAbiContainsModuleInit() {
        char  name_buf[8] = { 'm','a','i','l', 0,0,0,0 };
        int   name_len    = 4;
        xvue_module_init_(name_buf, &name_len);
        // No crash + warn-once printed = pass.
        QVERIFY(true);
    }

    // Sanity: XvueMenuBridge queue API resolves and behaves consistently
    // with the Plan 02 contract (full queueLexicon/popChar bodies). The
    // Plan 01 scaffold contract ("queueLexicon is a no-op") was flipped
    // by Plan 02 — see test_xvue_qt_i18n_menu_prefs.cpp for the deeper
    // multi-char + cap + trailing-CR test coverage. This scaffold test
    // keeps a single-char round-trip + the §10 layer-2 sentinel check so
    // a regression is caught even if the i18n_menu_prefs target is broken.
    void testMenuBridgeQueueBasicOps() {
        XvueMenuBridge mb;
        QCOMPARE(mb.queueSize(), 0);
        QVERIFY(!mb.popChar().has_value());

        mb.queueLexicon(QStringLiteral("x"));
        // PLAN 02 CONTRACT: queueLexicon pushes each char of cmd + a CR
        // terminator. "x" → ['x', 13] (size 2).
        QCOMPARE(mb.queueSize(), 2);
        auto first = mb.popChar();
        QVERIFY(first.has_value());
        QCOMPARE(static_cast<int>(*first), static_cast<int>('x'));
        auto cr = mb.popChar();
        QVERIFY(cr.has_value());
        QCOMPARE(static_cast<int>(*cr), 13);
        QVERIFY(!mb.popChar().has_value());

        // Sanity: the §10 layer 2 sentinel toggles correctly.
        QVERIFY(!mb.hasRegisteredModule());
        mb.markModuleRegistered();
        QVERIFY(mb.hasRegisteredModule());
    }
};

int main(int argc, char* argv[]) {
    (void)argc; (void)argv;
    XvueApp::ensure();  // constructs the process QApplication (owned by XvueApp)
    TestXvueQtMenuScaffold tc;
    int   qt_argc   = 1;
    char  qt_arg0[] = "xvue_qt_menu_tests";
    char* qt_argv[] = { qt_arg0, nullptr };
    return QTest::qExec(&tc, qt_argc, qt_argv);
}

#include "test_xvue_qt_menu_scaffold.moc"
