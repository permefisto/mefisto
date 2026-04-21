// xvue/qt/tests/test_xvue_qt_mail_menu.cpp
// Phase 6.1 Plan 03 Task 2 — D-13 mesher menu round-trip QTest cases.
//
// What this covers:
//   testMeshPointsCreate            — Mesh → Points → Create triggers "1;1;"
//   testFileImport                  — File → Import PLSVO triggers "80;"
//   testViewDrawMesh                — View → Draw PLSVO meshes triggers "10;"
//   testHelpAboutNoQueue            — shared 6.0 Help actions never queue
//   testCloseEventPushes99WhenBlocking
//                                  — D-09 closeEvent + blockingDepth() > 0
//                                    gate: push "99;" and ignore the close
//
// Pattern source:
//   - xvue/qt/tests/test_xvue_qt_event_menu_predrain.cpp (queueLexicon
//     round-trip via xvsouris_ with offscreen QPA + extern "C" decls)
//   - xvue/qt/tests/test_xvue_qt_window_chrome.cpp:185-226 (QAction::trigger)
//   - RESEARCH §Pitfall 4 (QAction::trigger, NOT postEvent+processEvents)
//   - RESEARCH §Example 1 (full 80-line template)
//
// Discipline:
//   - The production XvueMenuBridge (constructed in XvueWindow::ctor) is the
//     one registerMailActions_stub_ receives via xvue_module_init_. Our
//     xvsouris_ pre-drain reads from that same bridge, so the round-trip is
//     end-to-end without any test-only bridge injection.
//   - QAction::trigger() is synchronous on the main thread; the connected
//     lambda (see xvue_qt_mail_actions.cpp:72) calls queueLexicon directly,
//     so the chars are queued before trigger() returns.
#include "test_helpers.h"
#include "xvue_qt_app.h"
#include "xvue_qt_event.h"        // BlockingDepthGuard
#include "xvue_qt_menu_bridge.h"
#include "xvue_qt_window.h"

#include <QtTest/QtTest>
#include <QAction>
#include <QCloseEvent>
#include <QDir>
#include <QMenu>
#include <QMenuBar>
#include <QSettings>

extern "C" void xvinitgraphique_(void);
extern "C" void xvfermer_(void);
extern "C" void xvsouris_(int* notypeevent, int* nbc, int* x1, int* y1);
extern "C" void xvue_module_init_(char* name, int* name_len);

class TestXvueQtMailMenu : public QObject {
    Q_OBJECT

private:
    // Recursive helper — locate a QAction whose "lexicon" dynamic property
    // matches `lex`. Property naming matches makeLeafAction in Plan 02.
    // Walks submenus.
    static QAction* actionByLexicon(QMenu* menu, const QString& lex) {
        if (!menu) return nullptr;
        for (QAction* a : menu->actions()) {
            if (a->property("lexicon").toString() == lex) return a;
            if (a->menu()) {
                if (auto* child = actionByLexicon(a->menu(), lex)) return child;
            }
        }
        return nullptr;
    }

    // Drain N chars via xvsouris_ and assert the sequence matches `expected`.
    // Every read must return notypeevent=2 and nbc=<ASCII>.
    static void expectChars(const int* expected, int n) {
        for (int i = 0; i < n; ++i) {
            int nt = 0, nbc = 0, x = 0, y = 0;
            xvsouris_(&nt, &nbc, &x, &y);
            QCOMPARE(nt, 2);
            QCOMPARE(nbc, expected[i]);
        }
    }

private slots:
    void initTestCase() {
        // QSettings isolation: mirror test_xvue_qt_window_chrome.cpp so we do
        // not scribble ~/.config/LJLL/mefisto*.
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope,
                           QDir::tempPath());
        QCoreApplication::setOrganizationName(QStringLiteral("LJLL"));
        QCoreApplication::setApplicationName(
            QStringLiteral("mefisto-test-p6.1-mail-menu"));
        QSettings().clear();

        // Defensive: clear AUTOEXIT in case a parent shell exported it.
        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT");
        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS");

        // Construct the production window. menuBridge_ is created in the
        // ctor; registerMailActions_stub_ runs below and uses that exact
        // bridge, which is the same one xvsouris_'s pre-drain reads from.
        xvinitgraphique_();
        auto* win = XvueApp::window_slot().get();
        QVERIFY(win != nullptr);
        QVERIFY(win->menuBridge() != nullptr);

        // Wire the mail module. extern "C" signature is char*, not const
        // char*, so `name` must be a mutable buffer.
        char name[] = "mail";
        int  nlen   = 4;
        xvue_module_init_(name, &nlen);

        // Sanity: Plan 02 registered a Mesh menu at top level.
        QVERIFY(win->menuBar()->findChild<QMenu*>(QStringLiteral("Mesh"))
                != nullptr);
    }

    void cleanupTestCase() {
        // Drain whatever the tests left behind so the teardown is clean.
        auto* win = XvueApp::window_slot().get();
        if (win && win->menuBridge()) {
            while (win->menuBridge()->popChar().has_value()) {}
        }
        xvfermer_();
    }

    void init() {
        // Drain queue between tests so a leak in one test does not poison
        // the next.
        auto* win = XvueApp::window_slot().get();
        QVERIFY(win != nullptr);
        auto* mb = win->menuBridge();
        QVERIFY(mb != nullptr);
        while (mb->popChar().has_value()) {}
    }

    // ---- D-13 test 1: Mesh > Points > Create (lexicon "1;1;") ----
    // Expected chars after trigger(): '1', ';', '1', ';', 13 (CR).
    void testMeshPointsCreate() {
        auto* win = XvueApp::window_slot().get();
        auto* meshMenu = win->menuBar()->findChild<QMenu*>(
            QStringLiteral("Mesh"));
        QVERIFY(meshMenu != nullptr);
        auto* act = actionByLexicon(meshMenu, QStringLiteral("1;1;"));
        QVERIFY2(act != nullptr,
                 "Mesh > Points > Create QAction not found (lexicon=1;1;)");

        act->trigger();

        const int expected[5] = { '1', ';', '1', ';', 13 };
        expectChars(expected, 5);
    }

    // ---- D-13 test 2: File > Import PLSVO (lexicon "80;") ----
    // Expected chars: '8', '0', ';', 13.
    void testFileImport() {
        auto* win = XvueApp::window_slot().get();
        auto* fileMenu = win->menuBar()->findChild<QMenu*>(
            QStringLiteral("File"));
        QVERIFY(fileMenu != nullptr);
        auto* act = actionByLexicon(fileMenu, QStringLiteral("80;"));
        QVERIFY2(act != nullptr,
                 "File > Import PLSVO QAction not found (lexicon=80;)");

        act->trigger();

        const int expected[4] = { '8', '0', ';', 13 };
        expectChars(expected, 4);
    }

    // ---- D-13 test 3: View > Draw PLSVO meshes (lexicon "10;") ----
    // Expected chars: '1', '0', ';', 13.
    void testViewDrawMesh() {
        auto* win = XvueApp::window_slot().get();
        auto* viewMenu = win->menuBar()->findChild<QMenu*>(
            QStringLiteral("View"));
        QVERIFY(viewMenu != nullptr);
        auto* act = actionByLexicon(viewMenu, QStringLiteral("10;"));
        QVERIFY2(act != nullptr,
                 "View > Draw mesh QAction not found (lexicon=10;)");

        act->trigger();

        const int expected[4] = { '1', '0', ';', 13 };
        expectChars(expected, 4);
    }

    // ---- D-13 test 4: Help actions carry NO lexicon payload ----
    // The shared 6.0 Help > Documentation + About actions are generic; they
    // must not queue any char into the bridge.
    void testHelpAboutNoQueue() {
        auto* win = XvueApp::window_slot().get();
        auto* help = win->menuBar()->findChild<QMenu*>(QStringLiteral("Help"));
        QVERIFY(help != nullptr);

        auto* mb = win->menuBridge();
        const int before = mb->queueSize();

        for (QAction* a : help->actions()) {
            QVERIFY2(a->property("lexicon").toString().isEmpty(),
                     qPrintable(QStringLiteral(
                         "Help action \"%1\" carries lexicon payload \"%2\" "
                         "but should not (D-07 + Plan 02)")
                         .arg(a->text(),
                              a->property("lexicon").toString())));
        }

        QCOMPARE(mb->queueSize(), before);
    }

    // ---- Bonus: D-09 closeEvent pushes 99; when Fortran is blocking ----
    // Uses the BlockingDepthGuard RAII from xvue_qt_event.h (Phase 5) to
    // emulate the nested xvsouris_ state. When depth > 0, closeEvent must
    // (a) push the 4-char "99;"+CR sequence into the bridge and (b) ignore
    // the close event. When depth == 0, closeEvent must fall through to
    // QMainWindow::closeEvent and accept the close.
    //
    // closeEvent is protected, so we drive it through QWidget::close()
    // which invokes closeEvent() internally and returns true iff the event
    // was accepted (i.e., false if our D-09 guard called event->ignore()).
    void testCloseEventPushes99WhenBlocking() {
        auto* win = XvueApp::window_slot().get();
        QVERIFY(win != nullptr);
        auto* mb = win->menuBridge();
        QVERIFY(mb != nullptr);

        // Drain anything prior tests might have left.
        while (mb->popChar().has_value()) {}
        const int before = mb->queueSize();
        QCOMPARE(before, 0);

        bool closed = true;
        {
            BlockingDepthGuard g;   // depth becomes 1
            QCOMPARE(XvueApp::blockingDepth(), 1);

            // QWidget::close() calls closeEvent() synchronously and returns
            // true iff the event was accepted (not ignored).
            closed = win->close();
        }
        QCOMPARE(XvueApp::blockingDepth(), 0);

        QVERIFY2(!closed,
                 "closeEvent should call event->ignore() when "
                 "blockingDepth() > 0 (D-09)");

        // Bridge now contains ['9', '9', ';', 13].
        const int after = mb->queueSize();
        QCOMPARE(after - before, 4);

        const int expected[4] = { '9', '9', ';', 13 };
        for (int i = 0; i < 4; ++i) {
            auto ch = mb->popChar();
            QVERIFY(ch.has_value());
            QCOMPARE(static_cast<int>(static_cast<unsigned char>(*ch)),
                     expected[i]);
        }
    }
};

int main(int argc, char* argv[]) {
    (void)argc; (void)argv;
    qputenv("QT_QPA_PLATFORM", "offscreen");
    // CRITICAL: disable XvueConsoleDock::installStdoutRedirect for the
    // QTest run, otherwise xvue_module_init_'s dup2-over-STDOUT_FILENO
    // swallows QtTest's PASS/FAIL/Totals output and the suite appears to
    // hang (mirrors test_xvue_qt_window_chrome.cpp:257). Production
    // pp*_qt invocations do NOT set this var.
    qputenv("XVUE_QT_DISABLE_STDOUT_REDIRECT", "1");
    // AA_CompressHighFrequencyEvents + QApplication ownership (Phase 5 D-05).
    XvueApp::ensure();
    TestXvueQtMailMenu tc;
    int   qt_argc    = 1;
    char  qt_arg0[]  = "xvue_qt_mail_menu_tests";
    char* qt_argv[]  = { qt_arg0, nullptr };
    return QTest::qExec(&tc, qt_argc, qt_argv);
}

#include "test_xvue_qt_mail_menu.moc"
