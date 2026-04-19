// xvue/qt/tests/test_xvue_qt_event_menu_predrain.cpp
// Phase 6.0 Plan 03 (UX-03 / UX-12). QTest cases for the menu-queue
// pre-drain wired into XvueEventBridge::waitForEvent and the
// proc(xvsouris) wrapper. Also covers the eventFilter MouseMove
// mouseCoords emit (UX-12 — Plan 05 declared the signal + friendship,
// Plan 03 emits it).
//
// Test inventory:
//   testPreDrainSingleChar         — one char + trailing CR drain in order.
//   testPreDrainMultiChar          — "5;90;" round-trip via 6 reads.
//   testPreDrainTrailingCR         — CR is the LAST byte returned, value 13.
//   testPreDrainPreservesDepth     — early-return path still RAII-balances.
//   testPreDrainBeforeAutoexit     — menu queue beats AUTOEXIT short-circuit.
//   testPreDrainIgnoredInSouris2Mode — guard: xvsouris2_ does NOT pre-drain.
//   testPreDrainIgnoredInPauseMode — guard: xvpause_ does NOT pre-drain.
//   testCoordSignalDuringSourisMode — eventFilter emits mouseCoords on move.
//
// Design notes:
//   - Custom main() invokes XvueApp::ensure() before QTest::qExec so the
//     process QApplication is owned by XvueApp and AA_CompressHighFrequencyEvents
//     is set BEFORE the test classes construct (Phase 5 D-05/Pitfall 7).
//   - The test menu bridge is heap-allocated with parent = win (Qt
//     parent-child ownership: dies with window in cleanupTestCase).
//   - QSignalSpy uses the SIGNAL() macro form because the friend signal is
//     declared in xvue_qt_canvas.h and accessing &XvueCanvas::mouseCoords
//     from a non-friend test class requires an explicit qualified call.
//   - postMove uses QCoreApplication::postEvent — same helper pattern as
//     test_xvue_qt_event.cpp testMotionCoalescing — because QTest::mouseMove
//     uses QCursor::setPos which is a no-op on offscreen QPA.
#include "test_helpers.h"
#include "xvue_qt_app.h"
#include "xvue_qt_canvas.h"
#include "xvue_qt_event.h"
#include "xvue_qt_menu_bridge.h"
#include "xvue_qt_state.h"
#include "xvue_qt_window.h"

#include <QtTest/QtTest>
#include <QApplication>
#include <QCoreApplication>
#include <QEventLoop>
#include <QMouseEvent>
#include <QSignalSpy>
#include <QTimer>

extern "C" void xvinitgraphique_(void);
extern "C" void xvfermer_(void);
extern "C" void xvsouris_(int* notypeevent, int* nbc, int* x1, int* y1);
extern "C" void xvsouris2_(int* items, int* pmin0, int* notypeevent,
                           int* ibutton, int* x1, int* y1);
extern "C" void xvpause_(void);

class TestXvueQtEventMenuPredrain : public QObject {
    Q_OBJECT

    static void postMove(QWidget* target, QPoint p) {
        auto* ev = new QMouseEvent(QEvent::MouseMove, QPointF(p), QPointF(p),
                                   Qt::NoButton, Qt::NoButton, Qt::NoModifier);
        QCoreApplication::postEvent(target, ev);
    }

    XvueMenuBridge* mb_ = nullptr;

private slots:
    void initTestCase() {
        // Defensive: clear AUTOEXIT in case a parent shell exported it.
        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT");
        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS");

        xvinitgraphique_();
        auto* win = XvueApp::window_slot().get();
        QVERIFY(win != nullptr);
        QVERIFY(win->canvas() != nullptr);

        // Heap-allocate with parent = win so Qt parent-child ownership
        // destroys the bridge in xvfermer_'s ~XvueWindow chain.
        mb_ = new XvueMenuBridge(win);
        win->setMenuBridgeForTesting(mb_);
        QVERIFY(win->menuBridge() == mb_);
    }

    void cleanupTestCase() {
        // Drain anything tests left behind so xvfermer_ teardown is clean.
        if (mb_) {
            while (mb_->popChar().has_value()) {}
        }
        xvfermer_();
        // win destroyed -> mb_ destroyed (Qt parent-child).
        mb_ = nullptr;
    }

    void init() {
        // Drain queue between tests so a leak in one does not poison others.
        QVERIFY(mb_ != nullptr);
        while (mb_->popChar().has_value()) {}
        // Clear AUTOEXIT in case a previous test set it but failed to unset.
        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT");
        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS");
    }

    // ----- testPreDrainSingleChar -----
    // queueLexicon("a") enqueues ['a', 13]. Two consecutive xvsouris_ calls
    // must return notypeevent=2 with nbc='a' then nbc=13 (CR).
    void testPreDrainSingleChar() {
        mb_->queueLexicon(QStringLiteral("a"));
        int notype = -99, nbc = -99, x = -99, y = -99;

        xvsouris_(&notype, &nbc, &x, &y);
        QCOMPARE(notype, 2);
        QCOMPARE(nbc, (int)'a');
        QCOMPARE(XvueApp::blockingDepth(), 0);

        xvsouris_(&notype, &nbc, &x, &y);
        QCOMPARE(notype, 2);
        QCOMPARE(nbc, 13);   // trailing CR per saclav.f line-terminator contract
        QCOMPARE(XvueApp::blockingDepth(), 0);
    }

    // ----- testPreDrainMultiChar -----
    // "5;90;" -> ['5', ';', '9', '0', ';', 13]. Six reads in order.
    void testPreDrainMultiChar() {
        mb_->queueLexicon(QStringLiteral("5;90;"));
        const int expected[6] = {'5', ';', '9', '0', ';', 13};
        for (int i = 0; i < 6; ++i) {
            int notype = -99, nbc = -99, x = -99, y = -99;
            xvsouris_(&notype, &nbc, &x, &y);
            QCOMPARE(notype, 2);
            QCOMPARE(nbc, expected[i]);
        }
        QCOMPARE(XvueApp::blockingDepth(), 0);
    }

    // ----- testPreDrainTrailingCR -----
    // "99;" -> ['9', '9', ';', 13]. Fourth call must return CR=13.
    void testPreDrainTrailingCR() {
        mb_->queueLexicon(QStringLiteral("99;"));
        int notype = 0, nbc = 0, x = 0, y = 0;
        xvsouris_(&notype, &nbc, &x, &y);  // '9'
        QCOMPARE(nbc, (int)'9');
        xvsouris_(&notype, &nbc, &x, &y);  // '9'
        QCOMPARE(nbc, (int)'9');
        xvsouris_(&notype, &nbc, &x, &y);  // ';'
        QCOMPARE(nbc, (int)';');
        xvsouris_(&notype, &nbc, &x, &y);  // 13
        QCOMPARE(notype, 2);
        QCOMPARE(nbc, 13);
        QCOMPARE(XvueApp::blockingDepth(), 0);
    }

    // ----- testPreDrainPreservesDepth -----
    // The pre-drain early-return path lives in waitForEvent AFTER
    // BlockingDepthGuard ctor, so the guard's RAII still fires. Verify by
    // calling the bridge directly so the early-return branch is the one
    // exercised (the xvsouris_ wrapper short-circuits before entering
    // waitForEvent at all, which would NOT exercise the depth guard).
    void testPreDrainPreservesDepth() {
        mb_->queueLexicon(QStringLiteral("x"));
        const int before = XvueApp::blockingDepth();
        QCOMPARE(before, 0);

        auto* win = XvueApp::window_slot().get();
        QVERIFY(win != nullptr);
        auto* bridge = win->bridge();
        QVERIFY(bridge != nullptr);

        // Direct waitForEvent call exercises the in-bridge pre-drain that
        // sits between BlockingDepthGuard ctor and QEventLoop ctor.
        auto r = bridge->waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r.notypeevent, 2);
        QCOMPARE(r.nbc, (int)'x');

        const int after = XvueApp::blockingDepth();
        QCOMPARE(after, 0);   // RAII fires cleanly even on early return
    }

    // ----- testPreDrainBeforeAutoexit -----
    // queueLexicon("x") + AUTOEXIT armed. Pre-drain must win:
    //   call 1: pre-drain returns 'x'.
    //   call 2: pre-drain returns 13 (the trailing CR queueLexicon added).
    //   call 3: queue empty -> AUTOEXIT fires -> notypeevent=2, nbc=' '.
    void testPreDrainBeforeAutoexit() {
        qputenv("MEFISTO_XVSOURIS_AUTOEXIT", "1");
        // Tiny delay so AUTOEXIT does not slow the test suite.
        qputenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS", "10");

        mb_->queueLexicon(QStringLiteral("x"));

        int notype = 0, nbc = 0, x = 0, y = 0;
        // Call 1: pre-drain -> 'x'.
        xvsouris_(&notype, &nbc, &x, &y);
        QCOMPARE(notype, 2);
        QCOMPARE(nbc, (int)'x');

        // Call 2: pre-drain -> CR.
        xvsouris_(&notype, &nbc, &x, &y);
        QCOMPARE(notype, 2);
        QCOMPARE(nbc, 13);

        // Call 3: queue empty -> AUTOEXIT synthesizes SPACE.
        xvsouris_(&notype, &nbc, &x, &y);
        QCOMPARE(notype, 2);
        QCOMPARE(nbc, (int)' ');

        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT");
        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS");
    }

    // ----- testPreDrainIgnoredInSouris2Mode -----
    // xvsouris2_ does NOT pre-drain (the bridge's mode == WaitMode::Souris
    // guard rejects Souris2). Strategy:
    //   - queueLexicon("x") so the menu queue holds ['x', 13].
    //   - Arm AUTOEXIT so xvsouris2_ returns deterministically without
    //     blocking on a real event.
    //   - Call xvsouris2_; it should follow the AUTOEXIT path (returns
    //     notypeevent=2, ibutton=' ') and the queue must STILL hold ['x', 13]
    //     afterwards (verified by a follow-up plain xvsouris_ that drains 'x').
    void testPreDrainIgnoredInSouris2Mode() {
        qputenv("MEFISTO_XVSOURIS_AUTOEXIT", "1");
        qputenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS", "10");

        mb_->queueLexicon(QStringLiteral("x"));
        QCOMPARE(mb_->queueSize(), 2);   // ['x', 13]

        int items[6] = { 3, 1, 1,   50, 50, 1 };
        int pmin0 = -2;
        int notype = 0, ibutton = 0, x = 0, y = 0;
        xvsouris2_(items, &pmin0, &notype, &ibutton, &x, &y);

        // AUTOEXIT path returns notypeevent=2 with ibutton=' ' (space).
        QCOMPARE(notype, 2);
        QCOMPARE(ibutton, (int)' ');
        // Queue still intact — Souris2 did NOT consume the menu char.
        QCOMPARE(mb_->queueSize(), 2);

        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT");
        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS");

        // Follow-up Souris read drains 'x' to prove the queue was preserved.
        notype = 0; int nbc = 0; x = 0; y = 0;
        xvsouris_(&notype, &nbc, &x, &y);
        QCOMPARE(notype, 2);
        QCOMPARE(nbc, (int)'x');
        // Drain CR too so the next test starts clean.
        xvsouris_(&notype, &nbc, &x, &y);
        QCOMPARE(nbc, 13);
    }

    // ----- testPreDrainIgnoredInPauseMode -----
    // xvpause_ does NOT pre-drain. xvpause_ has no out-params so we verify
    // the queue STAYS intact across a xvpause_ call. AUTOEXIT armed so the
    // call returns deterministically without blocking.
    void testPreDrainIgnoredInPauseMode() {
        qputenv("MEFISTO_XVSOURIS_AUTOEXIT", "1");
        qputenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS", "10");

        mb_->queueLexicon(QStringLiteral("x"));
        QCOMPARE(mb_->queueSize(), 2);

        xvpause_();   // AUTOEXIT short-circuits; no menu drain.
        QCOMPARE(mb_->queueSize(), 2);   // queue untouched

        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT");
        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS");

        // Drain so the next test starts clean.
        int notype = 0, nbc = 0, x = 0, y = 0;
        xvsouris_(&notype, &nbc, &x, &y);  // 'x'
        QCOMPARE(nbc, (int)'x');
        xvsouris_(&notype, &nbc, &x, &y);  // CR
        QCOMPARE(nbc, 13);
    }

    // ----- testCoordSignalDuringSourisMode -----
    // The eventFilter MouseMove branch emits canvas_->mouseCoords(QPoint).
    // Verify by:
    //   - Connect QSignalSpy to mouseCoords.
    //   - Schedule a MouseMove at (42, 84) via QTimer (delivered after the
    //     bridge enters waitForEvent).
    //   - Call xvsouris_ (which enters the bridge's nested loop). The move
    //     fires the deferred-quit timer; waitForEvent returns notypeevent=-2.
    //   - Spy must have at least one entry, and the last point recorded must
    //     match (42, 84).
    void testCoordSignalDuringSourisMode() {
        auto* win = XvueApp::window_slot().get();
        QVERIFY(win != nullptr);
        auto* canvas = win->canvas();
        QVERIFY(canvas != nullptr);

        QSignalSpy spy(canvas, SIGNAL(mouseCoords(QPoint)));
        QVERIFY(spy.isValid());

        QTimer::singleShot(20, [canvas]{ postMove(canvas, QPoint(42, 84)); });

        int notype = 0, nbc = 0, x = 0, y = 0;
        xvsouris_(&notype, &nbc, &x, &y);

        QCOMPARE(notype, -2);                // motion return
        QCOMPARE(x, 42);
        QCOMPARE(y, 84);
        QVERIFY(spy.count() >= 1);

        // Read the last spy entry — coalescing may produce a single emit
        // for a single posted move; if multiple, the last one is the tail.
        const QList<QVariant> last = spy.last();
        QCOMPARE(last.size(), 1);
        QCOMPARE(last.first().toPoint(), QPoint(42, 84));
        QCOMPARE(XvueApp::blockingDepth(), 0);
    }
};

// Custom main: XvueApp::ensure() owns the QApplication so the
// AA_CompressHighFrequencyEvents attribute is set BEFORE the test class
// constructs (Phase 5 D-05). Mirrors test_xvue_qt_event.cpp main().
int main(int argc, char* argv[]) {
    (void)argc; (void)argv;
    XvueApp::ensure();
    TestXvueQtEventMenuPredrain tc;
    int    qt_argc   = 1;
    char   qt_arg0[] = "xvue_qt_event_menu_predrain_tests";
    char*  qt_argv[] = { qt_arg0, nullptr };
    return QTest::qExec(&tc, qt_argc, qt_argv);
}

#include "test_xvue_qt_event_menu_predrain.moc"
