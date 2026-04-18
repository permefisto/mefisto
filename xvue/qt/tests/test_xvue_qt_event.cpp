// xvue/qt/tests/test_xvue_qt_event.cpp
// Phase 5 Wave 0 (05-01 Task 1): QTest skeleton covering EVENT-01..EVENT-08.
// Plan 02 Task 2: replaces EVENT-08 QSKIPs with real RAII/nested-depth assertions
// and the EVENT-01 non-motion key/button rows. testBridgeInstallation is
// intentionally still QSKIPped here — Task 3 wires XvueWindow::bridge() and
// implements that test + testReopenCreatesFreshBridge.
#include "test_helpers.h"
#include "xvue_qt_app.h"
#include "xvue_qt_canvas.h"
#include "xvue_qt_event.h"
#include "xvue_qt_state.h"
#include "xvue_qt_window.h"

#include <QtTest/QtTest>
#include <QApplication>
#include <QCoreApplication>
#include <QEventLoop>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QPixmap>
#include <QTimer>

// Forward declarations for Fortran-ABI entry points used to stand up a real
// XvueWindow under the test. These are the byte-identical entries the real
// Fortran callers use; Task 3 adds XvueWindow::bridge() and testBridgeInstallation
// will then verify the bridge pointer. For Task 2 we instantiate the bridge
// manually on the existing canvas so the filter+loop plumbing is exercised
// without depending on the Task 3 wiring.
extern "C" void xvinitgraphique_(void);
extern "C" void xvfermer_(void);

class TestXvueQtEvent : public QObject {
    Q_OBJECT

    // Task 2 helper: synthesize a KeyPress event and post it to the canvas
    // via QCoreApplication::postEvent. QTest::keyClick is not used because
    // it requires the canvas to have active keyboard focus and the offscreen
    // QPA may not raise a focus-in event in time; postEvent goes straight
    // to the canvas's event queue where the bridge filter intercepts it.
    static void postKey(QWidget* target, int qtKey, const QString& text) {
        auto* ev = new QKeyEvent(QEvent::KeyPress, qtKey, Qt::NoModifier, text);
        QCoreApplication::postEvent(target, ev);
    }

    static void postButtonPress(QWidget* target, Qt::MouseButton b, QPoint p) {
        auto* ev = new QMouseEvent(QEvent::MouseButtonPress, QPointF(p), QPointF(p),
                                   b, b, Qt::NoModifier);
        QCoreApplication::postEvent(target, ev);
    }

    static void postButtonRelease(QWidget* target, Qt::MouseButton b, QPoint p) {
        auto* ev = new QMouseEvent(QEvent::MouseButtonRelease, QPointF(p), QPointF(p),
                                   b, Qt::NoButton, Qt::NoModifier);
        QCoreApplication::postEvent(target, ev);
    }

private slots:
    void initTestCase() {
        // XvueApp::ensure() was called by main() before qExec. Now stand up a
        // real window+canvas so the bridge has a valid target to filter.
        xvinitgraphique_();
        auto& win = XvueApp::window_slot();
        QVERIFY(win != nullptr);
        QVERIFY(win->canvas() != nullptr);
    }

    void cleanupTestCase() {
        // Tear down the window between test invocations so we do not leak the
        // process-wide window_slot into downstream test runners.
        xvfermer_();
    }

    // ---- EVENT-01: bridge installation ----
    void testBridgeInstallation() {
        // After initTestCase's xvinitgraphique_, XvueWindow should expose a
        // non-null bridge pointer and the bridge must be installed on the
        // canvas so keyboard/mouse events reach its filter. Drive a Space
        // keypress through the window's bridge and assert it dispatches.
        auto& win = XvueApp::window_slot();
        QVERIFY(win != nullptr);
        QVERIFY(win->canvas() != nullptr);
        QVERIFY(win->bridge() != nullptr);
        // Drive a space key and verify waitForEvent returns the expected
        // translation. No local bridge this time — we use the window's.
        auto* canvas = win->canvas();
        auto* bridge = win->bridge();
        QTimer::singleShot(10, [&]{ postKey(canvas, Qt::Key_Space, QStringLiteral(" ")); });
        auto r = bridge->waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r.notypeevent, 2);
        QCOMPARE(r.nbc, 32);
        QCOMPARE(XvueApp::blockingDepth(), 0);
    }

    // Phase 5 Plan 02 Task 3: xvfermer_ -> xvinitgraphique_ must yield a
    // brand-new bridge instance on the brand-new canvas; the previous
    // bridge was Qt-parent-owned by the old window and destroyed with it.
    // The pointer-identity checks are "best effort" — the heap allocator
    // CAN reuse the same address after free+new, in which case the test
    // falls back to the functional check (the new bridge must work).
    void testReopenCreatesFreshBridge() {
        auto& win = XvueApp::window_slot();
        QVERIFY(win != nullptr);
        auto* first_bridge = win->bridge();
        auto* first_canvas = win->canvas();
        QVERIFY(first_bridge != nullptr);
        QVERIFY(first_canvas != nullptr);

        xvfermer_();
        QVERIFY(XvueApp::window_slot() == nullptr);

        xvinitgraphique_();
        auto& win2 = XvueApp::window_slot();
        QVERIFY(win2 != nullptr);
        auto* second_bridge = win2->bridge();
        auto* second_canvas = win2->canvas();
        QVERIFY(second_bridge != nullptr);
        QVERIFY(second_canvas != nullptr);

        // Functional check: the fresh bridge dispatches events correctly.
        // This is the load-bearing invariant (the addresses may or may not
        // coincide depending on the heap allocator).
        QTimer::singleShot(10, [&]{ postKey(second_canvas, Qt::Key_Escape, QString()); });
        auto r = second_bridge->waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r.notypeevent, 0);
        QCOMPARE(r.nbc, 27);
        QCOMPARE(XvueApp::blockingDepth(), 0);
    }

    // ---- EVENT-02: xvsouris_ ----
    void testXvsourisMotion()             { QSKIP("Plan 03: motion coalescing not yet wired"); }

    void testXvsourisButtonPress() {
        auto* canvas = XvueApp::window_slot()->canvas();
        XvueEventBridge bridge(canvas);
        canvas->installEventFilter(&bridge);
        // ButtonPress alone: bridge returns notypeevent=-1, nbc=1 (left).
        QTimer::singleShot(0, [&]{ postButtonPress(canvas, Qt::LeftButton, QPoint(10, 20)); });
        auto r = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r.notypeevent, -1);
        QCOMPARE(r.nbc, 1);
        QCOMPARE(r.x, 10);
        QCOMPARE(r.y, 20);
        QCOMPARE(XvueApp::blockingDepth(), 0);
        canvas->removeEventFilter(&bridge);
    }

    void testXvsourisButtonRelease() {
        auto* canvas = XvueApp::window_slot()->canvas();
        XvueEventBridge bridge(canvas);
        canvas->installEventFilter(&bridge);
        // ButtonRelease alone: bridge returns notypeevent=1 (full click), nbc=1.
        QTimer::singleShot(0, [&]{ postButtonRelease(canvas, Qt::LeftButton, QPoint(30, 40)); });
        auto r = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r.notypeevent, 1);
        QCOMPARE(r.nbc, 1);
        QCOMPARE(r.x, 30);
        QCOMPARE(r.y, 40);
        QCOMPARE(XvueApp::blockingDepth(), 0);
        canvas->removeEventFilter(&bridge);
    }

    void testXvsourisKeyPress() {
        // Alias for "a space keypress returns notypeevent=2, nbc=32".
        auto* canvas = XvueApp::window_slot()->canvas();
        XvueEventBridge bridge(canvas);
        canvas->installEventFilter(&bridge);
        QTimer::singleShot(0, [&]{ postKey(canvas, Qt::Key_Space, QStringLiteral(" ")); });
        auto r = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r.notypeevent, 2);
        QCOMPARE(r.nbc, 32);
        QCOMPARE(XvueApp::blockingDepth(), 0);
        canvas->removeEventFilter(&bridge);
    }

    void testXvsourisKeyPressSpace() {
        // Plan 02 acceptance name: QTest::keyClick(Key_Space) returns nbc=32.
        auto* canvas = XvueApp::window_slot()->canvas();
        XvueEventBridge bridge(canvas);
        canvas->installEventFilter(&bridge);
        QTimer::singleShot(0, [&]{ postKey(canvas, Qt::Key_Space, QStringLiteral(" ")); });
        auto r = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r.notypeevent, 2);
        QCOMPARE(r.nbc, 32);
        canvas->removeEventFilter(&bridge);
    }

    void testXvsourisEscapeAbort() {
        // D-06 parity: Esc (27) -> notypeevent=0 (ABANDON).
        auto* canvas = XvueApp::window_slot()->canvas();
        XvueEventBridge bridge(canvas);
        canvas->installEventFilter(&bridge);
        QTimer::singleShot(0, [&]{ postKey(canvas, Qt::Key_Escape, QString()); });
        auto r = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r.notypeevent, 0);
        QCOMPARE(r.nbc, 27);
        canvas->removeEventFilter(&bridge);
    }

    void testXvsourisAtSignAbort() {
        // D-06 parity: '@' (64) -> notypeevent=0 (ABANDON) for French tutorials.
        auto* canvas = XvueApp::window_slot()->canvas();
        XvueEventBridge bridge(canvas);
        canvas->installEventFilter(&bridge);
        QTimer::singleShot(0, [&]{ postKey(canvas, Qt::Key_At, QStringLiteral("@")); });
        auto r = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r.notypeevent, 0);
        QCOMPARE(r.nbc, 64);
        canvas->removeEventFilter(&bridge);
    }

    void testXvsourisReturnKey() {
        // D-07: Return -> nbc=13 (CR), notypeevent=2.
        auto* canvas = XvueApp::window_slot()->canvas();
        XvueEventBridge bridge(canvas);
        canvas->installEventFilter(&bridge);
        // Empty text() forces the fallback switch to kick in for Key_Return;
        // QTest::keyClick normally fills text() with \r but we want to verify
        // the D-07 fallback also works.
        QTimer::singleShot(0, [&]{ postKey(canvas, Qt::Key_Return, QString()); });
        auto r = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r.notypeevent, 2);
        QCOMPARE(r.nbc, 13);
        canvas->removeEventFilter(&bridge);
    }

    // ---- EVENT-03: xvsouris2_ ----
    void testXvsouris2Accrochage()        { QSKIP("Plan 05: xvsouris2_ accrochage not yet wired"); }

    // ---- EVENT-04: xvpause_ ----
    void testXvpauseReturnsOnKey()        { QSKIP("Plan 04: xvpause_ ABI wiring not yet"); }
    void testXvpauseReturnsOnMouseClick() { QSKIP("Plan 04: xvpause_ ABI wiring not yet"); }
    void testXvpauseAutoexit()            { QSKIP("Plan 04: AUTOEXIT extension pending"); }

    // ---- EVENT-05: deplsouris_ ----
    void testDeplsourisNonBlocking()      { QSKIP("Plan 04: deplsouris_ not yet wired"); }

    // ---- EVENT-06: initaccrochage_ ----
    void testInitaccrochage()             { QSKIP("Plan 05: initaccrochage_ not yet wired"); }

    // ---- EVENT-07: motion coalescing ----
    void testMotionCoalescing()           { QSKIP("Plan 03: motion coalescing pending"); }

    // ---- EVENT-08: blocking depth ----
    void testBlockingDepthRAII() {
        // RAII guard: depth increments on construction, decrements on
        // destruction even when an exception unwinds. Exercise via a scoped
        // block.
        QCOMPARE(XvueApp::blockingDepth(), 0);
        {
            BlockingDepthGuard g1;
            QCOMPARE(XvueApp::blockingDepth(), 1);
        }
        QCOMPARE(XvueApp::blockingDepth(), 0);

        // Real waitForEvent call: ends with depth back at 0.
        auto* canvas = XvueApp::window_slot()->canvas();
        XvueEventBridge bridge(canvas);
        canvas->installEventFilter(&bridge);
        QTimer::singleShot(0, [&]{ postKey(canvas, Qt::Key_Escape, QString()); });
        auto r = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r.notypeevent, 0);
        QCOMPARE(r.nbc, 27);
        QCOMPARE(XvueApp::blockingDepth(), 0);
        canvas->removeEventFilter(&bridge);
    }

    void testBlockingDepthNested() {
        // Nested: inner waitForEvent runs while outer waitForEvent is still
        // blocking. Depth must reach 2 mid-stream and return to 0 after both
        // resolve.
        auto* canvas = XvueApp::window_slot()->canvas();
        XvueEventBridge bridge(canvas);
        canvas->installEventFilter(&bridge);

        int observed_outer_depth  = -1;  // inside outer's event callback
        int observed_nested_depth = -1;  // inside inner's event callback
        int observed_between      = -1;  // between inner return and outer quit

        // When the outer loop starts pumping events, this fires first.
        // Note: small delays (not 0) for both timers — empirically xvfb-run
        // can starve 0-ms timers when the outer waitForEvent is nested inside
        // the QTest runner's own event pump. 5 ms is enough to let both loops
        // enter their exec() cleanly and the timer callbacks fire reliably.
        QTimer::singleShot(10, [&]{
            observed_outer_depth = XvueApp::blockingDepth();  // expect 1

            // Schedule an event that fires INSIDE the inner waitForEvent
            // loop so we can observe depth from there.
            QTimer::singleShot(10, [&]{
                observed_nested_depth = XvueApp::blockingDepth();  // expect 2
                // Post Esc to quit the inner loop.
                postKey(canvas, Qt::Key_Escape, QString());
            });

            // Enter the nested waitForEvent.
            auto inner = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
            QCOMPARE(inner.notypeevent, 0);
            QCOMPARE(inner.nbc, 27);

            observed_between = XvueApp::blockingDepth();  // expect 1

            // Post Esc to quit the outer loop.
            postKey(canvas, Qt::Key_Escape, QString());
        });

        auto outer = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(outer.notypeevent, 0);
        QCOMPARE(outer.nbc, 27);
        QCOMPARE(observed_outer_depth,  1);
        QCOMPARE(observed_nested_depth, 2);
        QCOMPARE(observed_between,      1);
        QCOMPARE(XvueApp::blockingDepth(), 0);
        canvas->removeEventFilter(&bridge);
    }

    // Phase 5 Wave 0 Task 2: XvueApp::blockingDepth() must return 0 in a
    // fresh process with no nested waitForEvent() calls active.
    void testBlockingDepthAccessorZero() {
        QCOMPARE(XvueApp::blockingDepth(), 0);
    }

    // Phase 5 Wave 0 Task 3 (D-08, Pitfall 10): fresh XvueState must have
    // both accrochage pixmap pointers nullptr (they are lazily allocated).
    void testXvueStateAccrochageFieldsNull() {
        XvueState s;
        QCOMPARE(s.mempxaccro_, static_cast<QPixmap*>(nullptr));
        QCOMPARE(s.accroche_undo_tile_, static_cast<QPixmap*>(nullptr));
    }

    // Phase 5 Wave 0 Task 3 (Pitfall 4): XvueCanvas must have StrongFocus.
    // Full canvas construction requires window/state wiring Plan 02 will
    // handle; this stub records the intent in the matrix.
    void testCanvasFocusPolicyStrong() {
        auto* canvas = XvueApp::window_slot()->canvas();
        QVERIFY(canvas != nullptr);
        QCOMPARE(canvas->focusPolicy(), Qt::StrongFocus);
    }

    // Phase 5 Wave 0 Task 2 (D-05, Pitfall 7): AA_CompressHighFrequencyEvents
    // must be observable once XvueApp::ensure() has constructed the
    // QApplication. main() below calls XvueApp::ensure() before QTest::qExec,
    // so the attribute set inside the ensure() call_once lambda is visible
    // here regardless of X11's default-true behavior.
    void testCompressHighFrequencyEventsSet() {
        QVERIFY(QCoreApplication::testAttribute(Qt::AA_CompressHighFrequencyEvents));
    }
};

// Phase 5 Wave 0 Task 2: custom main instead of QTEST_MAIN so that
// XvueApp::ensure() owns QApplication construction (and therefore gets to
// set Qt::AA_CompressHighFrequencyEvents BEFORE the ctor runs — D-05).
// Using QTEST_MAIN would construct a second QApplication and trip Qt's
// "only one QApplication instance" assertion.
int main(int argc, char* argv[]) {
    (void)argc; (void)argv;
    XvueApp::ensure();  // constructs the process QApplication (owned by XvueApp)
    TestXvueQtEvent tc;
    // Hand qExec a minimal argv — QApplication is already alive.
    int    qt_argc     = 1;
    char   qt_arg0[]   = "xvue_qt_event_tests";
    char*  qt_argv[]   = { qt_arg0, nullptr };
    return QTest::qExec(&tc, qt_argc, qt_argv);
}

#include "test_xvue_qt_event.moc"
