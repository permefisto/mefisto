// xvue/qt/tests/test_xvue_qt_nlse_menu.cpp
// Phase 6.5 Plan 03 Task 2 -- D-13 nlse menu round-trip QTest cases.
//
// What this covers (mirrors 6.4 test_xvue_qt_ther_menu.cpp, substitute
// "ther"->"nlse" and "Thermal"->"Nonlinear"; leaf lexicons drawn from the
// frozen LEXICON-AUDIT-nlse.md and matched against Plan 02's wiring in
// xvue/qt/src/xvue_qt_nlse_actions.cpp):
//
//   testNonlinearSolveImplicit       -- Nonlinear > Solve Implicit ("5;")
//   testViewDrawNlseModule           -- View      > Draw NLSE Module ("11;")
//   testFileTMSTools                 -- File      > SUIVI TMS ("71;")
//   testHelpNoQueue                  -- Help carries only audit-allowlisted
//                                       {97;} (MEFISTO Date Version) -- NOT
//                                       {98;} like ther; matches flui's
//                                       value. Codifies 6.3 Plan 03
//                                       Auto-fix Rule 1 lesson by reading the
//                                       LEXICON-AUDIT-nlse.md row 97.
//   testCloseEventPushes99WhenBlocking
//                                    -- D-09 closeEvent + blockingDepth() > 0
//                                       gate: push "99;" and ignore the close
//                                       (shared XvueWindow behaviour inherited
//                                       by ppnlse_qt from 6.1 Plan 03).
//   testMenuOrder                    -- {File, Nonlinear, View, Help} order
//                                       locked (codifies 6.2 Plan 04 fix as
//                                       automated gate, replacing the 6.2
//                                       manual A/B checkpoint).
//   testBilingualLabelsEnglish       -- EN labels routed via td/ma/debunlse
//                                       when anglais flag set; QSKIPs cleanly
//                                       on missing preconditions (codifies the
//                                       6.2 Plan 05 fix as an automated gate).
//
// Pattern source:
//   - xvue/qt/tests/test_xvue_qt_ther_menu.cpp (6.4 D-13 template)
//   - xvue/qt/tests/test_xvue_qt_flui_menu.cpp (6.3 D-13 template -- {97;})
//   - RESEARCH §Pitfall 4 (QAction::trigger, NOT postEvent+processEvents)
//
// Discipline:
//   - xvue_module_init_('nlse', 4) is called in initTestCase so the Plan 02
//     registerNlseActions_stub_ runs against the production XvueWindow +
//     XvueMenuBridge (same pair that xvsouris_ drains from).
//   - QAction::trigger() is synchronous on the main thread; the connected
//     lambda in makeLeafAction (xvue_qt_nlse_actions.cpp) calls queueLexicon
//     so the chars are queued before trigger() returns.
//   - actionByLexicon walks submenus so "5;" (inside Nonlinear) and "11;"
//     (inside View) are reachable by starting from the top-level
//     Nonlinear/View menus respectively.
#include "test_helpers.h"
#include "xvue_qt_app.h"
#include "xvue_qt_event.h"        // BlockingDepthGuard
#include "xvue_qt_i18n.h"         // xvueClearLanguageCacheForTesting
#include "xvue_qt_menu_bridge.h"
#include "xvue_qt_menu_file_parser.h"
#include "xvue_qt_window.h"

#include <QtTest/QtTest>
#include <QAction>
#include <QCloseEvent>
#include <QDir>
#include <QFile>
#include <QMenu>
#include <QMenuBar>
#include <QSet>
#include <QSettings>

extern "C" void xvinitgraphique_(void);
extern "C" void xvfermer_(void);
extern "C" void xvsouris_(int* notypeevent, int* nbc, int* x1, int* y1);
extern "C" void xvue_module_init_(char* name, int* name_len);

class TestXvueQtNlseMenu : public QObject {
    Q_OBJECT

private:
    // Recursive helper -- locate a QAction whose "lexicon" dynamic property
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
        // QSettings isolation: mirror 6.4 ther-menu test so we do not
        // scribble ~/.config/LJLL/mefisto*.
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope,
                           QDir::tempPath());
        QCoreApplication::setOrganizationName(QStringLiteral("LJLL"));
        QCoreApplication::setApplicationName(
            QStringLiteral("mefisto-test-p6.5-nlse-menu"));
        QSettings().clear();

        // Defensive: clear AUTOEXIT in case a parent shell exported it.
        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT");
        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS");

        // Construct the production window. menuBridge_ is created in the
        // ctor; registerNlseActions_stub_ runs below and uses that exact
        // bridge, which is the same one xvsouris_'s pre-drain reads from.
        xvinitgraphique_();
        auto* win = XvueApp::window_slot().get();
        QVERIFY(win != nullptr);
        QVERIFY(win->menuBridge() != nullptr);

        // Wire the nlse module. extern "C" signature is char*, not const
        // char*, so `name` must be a mutable buffer.
        char name[] = "nlse";
        int  nlen   = 4;
        xvue_module_init_(name, &nlen);

        // Sanity: Plan 02 registered a Nonlinear menu at top level (plus File,
        // View, Help inherited from shared 6.0 chrome).
        QVERIFY(win->menuBar()->findChild<QMenu*>(QStringLiteral("Nonlinear"))
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

    // ---- D-13 test 1: Nonlinear > Solve Implicit (lexicon "5;") ----
    // Expected chars after trigger(): '5', ';', 13 (CR appended by
    // queueLexicon). ppnlse.f:340-343 sets TESTNL=7 + CALL NLSEINS
    // (IMPLICIT scheme with global 2nx2n non-symmetric matrix).
    void testNonlinearSolveImplicit() {
        auto* win = XvueApp::window_slot().get();
        auto* nonlinearMenu = win->menuBar()->findChild<QMenu*>(
            QStringLiteral("Nonlinear"));
        QVERIFY(nonlinearMenu != nullptr);
        auto* act = actionByLexicon(nonlinearMenu, QStringLiteral("5;"));
        QVERIFY2(act != nullptr,
                 "Nonlinear > Solve Implicit QAction not found (lexicon=5;)");

        act->trigger();

        const int expected[3] = { '5', ';', 13 };
        expectChars(expected, 3);
    }

    // ---- D-13 test 2: View > Draw NLSE Module (lexicon "11;") ----
    // Expected chars: '1', '1', ';', 13. ppnlse.f:391-393 CALL TRTHER(KNOMOB,
    // 4, IERR) -- NTYPDESS=4 selects MODULE (the abs(U) primary visual).
    void testViewDrawNlseModule() {
        auto* win = XvueApp::window_slot().get();
        auto* viewMenu = win->menuBar()->findChild<QMenu*>(
            QStringLiteral("View"));
        QVERIFY(viewMenu != nullptr);
        auto* act = actionByLexicon(viewMenu, QStringLiteral("11;"));
        QVERIFY2(act != nullptr,
                 "View > Draw NLSE Module QAction not found (lexicon=11;)");

        act->trigger();

        const int expected[4] = { '1', '1', ';', 13 };
        expectChars(expected, 4);
    }

    // ---- D-13 test 3: File > SUIVI TMS (lexicon "71;") ----
    // Expected chars: '7', '1', ';', 13. NOTE: nlse exposes 71; as the FIRST
    // of FOUR separate top-level TMS codes (71;72;73;74;) instead of the
    // unified 70; that ther/elas/flui use. This test pins the 71; binding.
    // ppnlse.f:217 IF NMTCL.GT.70 GOTO 7001 then NMTCL=NMTCL-70=1 ->
    // label 7100 (CALL SUITMS).
    void testFileTMSTools() {
        auto* win = XvueApp::window_slot().get();
        auto* fileMenu = win->menuBar()->findChild<QMenu*>(
            QStringLiteral("File"));
        QVERIFY(fileMenu != nullptr);
        auto* act = actionByLexicon(fileMenu, QStringLiteral("71;"));
        QVERIFY2(act != nullptr,
                 "File > SUIVI TMS QAction not found (lexicon=71;)");

        act->trigger();

        const int expected[4] = { '7', '1', ';', 13 };
        expectChars(expected, 4);
    }

    // ---- D-13 test 4: Help menu carries only the audited 97; lexicon ----
    // The shared 6.0 Help > Documentation + About actions are generic; they
    // must not queue any char into the bridge. Plan 02 added a single
    // module-specific Help leaf for "97;" (Mefisto Date Version, per
    // LEXICON-AUDIT-nlse.md row "97; ... DATE de version de Mefisto;
    // Help-allowlist for Plan 03 testHelpNoQueue: {97;}"). Any OTHER
    // lexicon payload in Help is an unintended bug.
    //
    // CRITICAL: nlse's Help-allowlist is `{97;}`, MATCHING flui (6.3),
    // but DIFFERING from ther's `{98;}` (6.4). This is the 6.3 Plan 03
    // Auto-fix Rule 1 lesson: read the per-module LEXICON-AUDIT BEFORE
    // writing the test. Per LEXICON-AUDIT-nlse.md §"Help-allowlist (for
    // Plan 03 testHelpNoQueue -- explicit hand-off)": the Help menu
    // carries `{97;}`.
    void testHelpNoQueue() {
        auto* win = XvueApp::window_slot().get();
        auto* help = win->menuBar()->findChild<QMenu*>(QStringLiteral("Help"));
        QVERIFY(help != nullptr);

        auto* mb = win->menuBridge();
        const int before = mb->queueSize();

        // Allowed Help lexicon set per LEXICON-AUDIT-nlse.md row 97.
        // nlse: {97;}; flui: {97;}; ther: {98;}; mail/elas: {} (empty).
        const QSet<QString> allowedHelpLexicons = {
            QStringLiteral("97;"),
        };

        for (QAction* a : help->actions()) {
            const QString lex = a->property("lexicon").toString();
            if (lex.isEmpty()) continue;
            QVERIFY2(allowedHelpLexicons.contains(lex),
                     qPrintable(QStringLiteral(
                         "Help action \"%1\" carries unexpected lexicon "
                         "payload \"%2\" (only {97;} is audited for nlse)")
                         .arg(a->text(), lex)));
        }

        // testHelpNoQueue must NOT trigger any action -- just inspect.
        // queueSize stays unchanged.
        QCOMPARE(mb->queueSize(), before);
    }

    // ---- 6.2 Plan 04 carry-over: menu bar order is exactly
    //                              {File, Nonlinear, View, Help} ----
    // Codifies the gap-1 fix as an automated regression guard. Plan 02
    // landed ensureTopLevelMenu(insertBefore=viewMenuForAnchor) for the
    // Nonlinear menu so it lands BETWEEN File and View. This slot replaces
    // the 6.2 manual A/B menu-order checkpoint with a machine-verifiable
    // gate (Phase 6.5 is fully autonomous because of this slot).
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
            QStringLiteral("Nonlinear"),
            QStringLiteral("View"),
            QStringLiteral("Help"),
        };
        QStringList actual = topLevelObjectNames.mid(0, 4);
        QCOMPARE(actual, expected);
    }

    // ---- 6.2 Plan 05 carry-over: bilingual labels render in EN under
    //                              anglais flag ----
    // Codifies the gap-2 fix as an automated regression guard.
    // XvueMenuFileParser::loadFor("debunlse") consults td/ma/debunlse first
    // when xvueIsEnglish() returns true and the file exists. This slot
    // replaces the 6.2 manual A/B bilingual-labels checkpoint with a
    // machine-verifiable gate (Phase 6.5 is fully autonomous because of
    // this slot).
    //
    // QSKIP if the test env doesn't have the bilingual td/ma/ tree
    // -- this slot deliberately does NOT mock the filesystem because
    // mocking would not catch real-install regressions.
    void testBilingualLabelsEnglish() {
        const QByteArray mefisto = qgetenv("MEFISTO");
        if (mefisto.isEmpty()) {
            QSKIP("MEFISTO env not set -- cannot probe td/ma/ tree");
        }
        const QString home = QString::fromLocal8Bit(mefisto);
        if (!QFile::exists(home + QStringLiteral("/td/m/anglais"))) {
            QSKIP("td/m/anglais flag not present -- EN path not exercised");
        }
        if (!QFile::exists(home + QStringLiteral("/td/ma/debunlse"))) {
            QSKIP("td/ma/debunlse missing -- cannot test EN parser path");
        }

        // Drop any cached MenuFile AND the language probe from initTestCase
        // so loadFor() re-resolves both the file path and the EN/FR flag
        // now that we've verified the on-disk preconditions hold.
        // (WR-02 gap-closure: clearCacheForTesting alone did not reset
        // the xvueIsEnglish() static, which could lock to false if a
        // prior slot called it before td/m/anglais was confirmed present.)
        xvueClearLanguageCacheForTesting();
        XvueMenuFileParser::clearCacheForTesting();

        const MenuFile& mf = XvueMenuFileParser::loadFor(
            QStringLiteral("debunlse"));
        QVERIFY2(mf.ok(),
                 "td/ma/debunlse parsed empty -- file may be malformed");

        // Verify a couple of canonical EN strings from td/ma/debunlse.
        // Code 5;:   IMPLICIT i dU(t)/dt -Alfa Laplacian U(t) +N(...)U(t)=F(t,X)
        // Code 11;:  Drawing of abs(U(t,X)) NLSE WAVE MODULE
        // Code 99;:  SAVE DATA and QUIT
        // Case-insensitive substring matching so a future cosmetic edit
        // (e.g., trimmed double-space) does not break the test.
        const QString label5 = mf.label(5);
        QVERIFY2(label5.contains(QStringLiteral("IMPLICIT"),
                                 Qt::CaseInsensitive)
                 || label5.contains(QStringLiteral("LAPLACIAN"),
                                    Qt::CaseInsensitive),
                 qPrintable(QStringLiteral(
                     "Expected EN implicit-NLSE label for 5;, got: ") + label5));
        const QString label11 = mf.label(11);
        QVERIFY2(label11.contains(QStringLiteral("MODULE"),
                                  Qt::CaseInsensitive)
                 || label11.contains(QStringLiteral("WAVE"),
                                     Qt::CaseInsensitive),
                 qPrintable(QStringLiteral(
                     "Expected EN modulus label for 11;, got: ") + label11));
        const QString label99 = mf.label(99);
        QVERIFY2(label99.contains(QStringLiteral("SAVE"),
                                  Qt::CaseInsensitive)
                 || label99.contains(QStringLiteral("QUIT"),
                                     Qt::CaseInsensitive),
                 qPrintable(QStringLiteral(
                     "Expected EN save/quit label for 99;, got: ") + label99));

        // Defensive: clear both the MenuFile cache and the language probe
        // so subsequent slots in this class probe afresh.
        xvueClearLanguageCacheForTesting();
        XvueMenuFileParser::clearCacheForTesting();
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
    // This covers the pp/ppnlse_qt inheritance of 6.1 Plan 03 Task 2's
    // shared XvueWindow rewrite -- Plan 03 6.5 MUST NOT re-edit that file.
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
                 "blockingDepth() > 0 (D-09, shared with 6.1 ppmail_qt + "
                 "6.2 ppelas_qt + 6.3 ppflui_qt + 6.4 ppther_qt)");

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
    // -testcase) are therefore unreachable from this binary by design --
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
    TestXvueQtNlseMenu tc;
    int   qt_argc    = 1;
    char  qt_arg0[]  = "xvue_qt_nlse_menu_tests";
    char* qt_argv[]  = { qt_arg0, nullptr };
    return QTest::qExec(&tc, qt_argc, qt_argv);
}

#include "test_xvue_qt_nlse_menu.moc"
