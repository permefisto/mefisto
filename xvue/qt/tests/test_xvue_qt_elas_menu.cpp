// xvue/qt/tests/test_xvue_qt_elas_menu.cpp
// Phase 6.2 Plan 03 Task 2 — D-13 elasticity menu round-trip QTest cases.
//
// What this covers (mirrors 6.1 test_xvue_qt_mail_menu.cpp, substitute
// "mail"->"elas" and "Mesh"->"Solve"; leaf lexicons drawn from the frozen
// LEXICON-AUDIT-elas.md and matched against Plan 02's wiring in
// xvue/qt/src/xvue_qt_elas_actions.cpp):
//
//   testSolveStatic                 — Solve > Static > "Run static solve..."
//                                     fires lexicon "3;"
//   testViewDrawStress              — View > "Draw displacements & stresses..."
//                                     fires lexicon "8;"
//   testFileMSTools                 — File > SUIVI FICHIERS MS fires "72;"
//   testHelpNoQueue                 — shared 6.0 Help actions never queue
//   testCloseEventPushes99WhenBlocking
//                                   — D-09 closeEvent + blockingDepth() > 0
//                                     gate: push "99;" and ignore the close
//                                     (shared XvueWindow behaviour inherited
//                                     by ppelas_qt from 6.1 Plan 03).
//
// Pattern source:
//   - xvue/qt/tests/test_xvue_qt_mail_menu.cpp (6.1 D-13 template)
//   - RESEARCH §Pitfall 4 (QAction::trigger, NOT postEvent+processEvents)
//
// Discipline:
//   - xvue_module_init_('elas', 4) is called in initTestCase so the Plan 02
//     registerElasActions_stub_ runs against the production XvueWindow +
//     XvueMenuBridge (same pair that xvsouris_ drains from).
//   - QAction::trigger() is synchronous on the main thread; the connected
//     lambda in makeLeafAction (xvue_qt_elas_actions.cpp) calls queueLexicon
//     so the chars are queued before trigger() returns.
//   - actionByLexicon walks submenus so "3;" (inside Solve.Static) and "8;"
//     (inside View.DrawStress) are both reachable by starting from the
//     top-level Solve/View menus respectively.
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

class TestXvueQtElasMenu : public QObject {
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
        // QSettings isolation: mirror 6.1 mail-menu test so we do not
        // scribble ~/.config/LJLL/mefisto*.
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope,
                           QDir::tempPath());
        QCoreApplication::setOrganizationName(QStringLiteral("LJLL"));
        QCoreApplication::setApplicationName(
            QStringLiteral("mefisto-test-p6.2-elas-menu"));
        QSettings().clear();

        // Defensive: clear AUTOEXIT in case a parent shell exported it.
        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT");
        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS");

        // Construct the production window. menuBridge_ is created in the
        // ctor; registerElasActions_stub_ runs below and uses that exact
        // bridge, which is the same one xvsouris_'s pre-drain reads from.
        xvinitgraphique_();
        auto* win = XvueApp::window_slot().get();
        QVERIFY(win != nullptr);
        QVERIFY(win->menuBridge() != nullptr);

        // Wire the elas module. extern "C" signature is char*, not const
        // char*, so `name` must be a mutable buffer.
        char name[] = "elas";
        int  nlen   = 4;
        xvue_module_init_(name, &nlen);

        // Sanity: Plan 02 registered a Solve menu at top level (plus File,
        // View, Help inherited from shared 6.0 chrome).
        QVERIFY(win->menuBar()->findChild<QMenu*>(QStringLiteral("Solve"))
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

    // ---- D-13 test 1: Solve > Static > "Run static solve..." (lexicon "3;") ----
    // Expected chars after trigger(): '3', ';', 13 (CR appended by queueLexicon).
    void testSolveStatic() {
        auto* win = XvueApp::window_slot().get();
        auto* solveMenu = win->menuBar()->findChild<QMenu*>(
            QStringLiteral("Solve"));
        QVERIFY(solveMenu != nullptr);
        auto* act = actionByLexicon(solveMenu, QStringLiteral("3;"));
        QVERIFY2(act != nullptr,
                 "Solve > Static > static-solve QAction not found (lexicon=3;)");

        act->trigger();

        const int expected[3] = { '3', ';', 13 };
        expectChars(expected, 3);
    }

    // ---- D-13 test 2: View > "Draw displacements & stresses..." (lexicon "8;") ----
    // Expected chars: '8', ';', 13.
    void testViewDrawStress() {
        auto* win = XvueApp::window_slot().get();
        auto* viewMenu = win->menuBar()->findChild<QMenu*>(
            QStringLiteral("View"));
        QVERIFY(viewMenu != nullptr);
        auto* act = actionByLexicon(viewMenu, QStringLiteral("8;"));
        QVERIFY2(act != nullptr,
                 "View > Draw stress QAction not found (lexicon=8;)");

        act->trigger();

        const int expected[3] = { '8', ';', 13 };
        expectChars(expected, 3);
    }

    // ---- D-13 test 3: File > SUIVI FICHIERS de la MS (lexicon "72;") ----
    // Expected chars: '7', '2', ';', 13.
    void testFileMSTools() {
        auto* win = XvueApp::window_slot().get();
        auto* fileMenu = win->menuBar()->findChild<QMenu*>(
            QStringLiteral("File"));
        QVERIFY(fileMenu != nullptr);
        auto* act = actionByLexicon(fileMenu, QStringLiteral("72;"));
        QVERIFY2(act != nullptr,
                 "File > MS File Tools QAction not found (lexicon=72;)");

        act->trigger();

        const int expected[4] = { '7', '2', ';', 13 };
        expectChars(expected, 4);
    }

    // ---- D-13 test 4: Help actions carry NO lexicon payload ----
    // The shared 6.0 Help > Documentation + About actions are generic; they
    // must not queue any char into the bridge.
    void testHelpNoQueue() {
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

    // ---- Plan 04 Gap-1: menu bar order is exactly {File, Solve, View, Help} ----
    // 06.2-HUMAN-UAT.md gap 1 root cause: ensureTopLevelMenu used
    // mbar->addMenu(...) which appends. Plan 04 Task 1 routes the
    // Solve menu through QMenuBar::insertMenu(viewMenu->menuAction(),
    // solveMenu) so it lands BETWEEN File and View. This test codifies
    // the expected sequence so 6.3/6.4/6.5 cannot silently regress it.
    void testMenuOrder() {
        auto* win = XvueApp::window_slot().get();
        QVERIFY(win != nullptr);
        auto* mb = win->menuBar();
        QVERIFY(mb != nullptr);

        QStringList topLevelObjectNames;
        for (QAction* a : mb->actions()) {
            if (a->menu()) {
                topLevelObjectNames << a->menu()->objectName();
            }
        }

        QVERIFY2(topLevelObjectNames.size() >= 4,
                 qPrintable(QStringLiteral(
                     "expected at least 4 top-level menus, got %1: [%2]")
                     .arg(topLevelObjectNames.size())
                     .arg(topLevelObjectNames.join(QStringLiteral(", ")))));

        const QStringList expected = {
            QStringLiteral("File"),
            QStringLiteral("Solve"),
            QStringLiteral("View"),
            QStringLiteral("Help"),
        };
        QStringList actual = topLevelObjectNames.mid(0, 4);
        QCOMPARE(actual, expected);
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
    //
    // This covers the pp/ppelas_qt inheritance of 6.1 Plan 03 Task 2's
    // shared XvueWindow rewrite — Plan 03 6.2 MUST NOT re-edit that file.
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
                 "blockingDepth() > 0 (D-09, shared with 6.1 ppmail_qt)");

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
    // 06.2-REVIEW IN-04: Real argc/argv are intentionally NOT forwarded to
    // QTest::qExec. XvueApp::ensure() (called below) constructs
    // QApplication with a process-static fake_argv array that must outlive
    // every QApplication API call. Forwarding the real argv would create
    // a dangling pointer once main() returns. QTest CLI flags (-v2, -o,
    // -testcase) are therefore unreachable from this binary by design —
    // run with QT_QPA_PLATFORM=offscreen and accept the default verbosity.
    (void)argc; (void)argv;
    qputenv("QT_QPA_PLATFORM", "offscreen");
    // 6.1 auto-fix #3: disable XvueConsoleDock::installStdoutRedirect for the
    // QTest run, otherwise xvue_module_init_'s dup2-over-STDOUT_FILENO
    // swallows QtTest's PASS/FAIL/Totals output and the suite appears to
    // hang. Production pp*_qt invocations do NOT set this var.
    qputenv("XVUE_QT_DISABLE_STDOUT_REDIRECT", "1");
    // AA_CompressHighFrequencyEvents + QApplication ownership (Phase 5 D-05).
    XvueApp::ensure();
    TestXvueQtElasMenu tc;
    int   qt_argc    = 1;
    char  qt_arg0[]  = "xvue_qt_elas_menu_tests";
    char* qt_argv[]  = { qt_arg0, nullptr };
    return QTest::qExec(&tc, qt_argc, qt_argv);
}

#include "test_xvue_qt_elas_menu.moc"
