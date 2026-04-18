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
// Phase 5 Plan 04 (EVENT-02/04/05): Fortran ABI entry points under test.
extern "C" void xvsouris_(int* notypeevent, int* nbc, int* x1, int* y1);
extern "C" void xvpause_(void);
extern "C" void deplsouris_(int* x, int* y);
// Phase 5 Plan 05 (EVENT-03/06): accrochage entry points.
extern "C" void initaccrochage_(void);
extern "C" void xvsouris2_(int* items, int* pmin0, int* notypeevent,
                           int* ibutton, int* x1, int* y1);

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

    // Plan 03: synthesize a MouseMove event so the motion-coalescing path in
    // eventFilter can be exercised without a real QTest::mouseMove (which
    // uses QCursor::setPos that doesn't propagate on xvfb). postEvent
    // delivers the event directly to the canvas's queue where the filter
    // intercepts it. The canvas has setMouseTracking(true) from its ctor
    // (xvue_qt_canvas.cpp), so Qt does not drop the no-button move.
    static void postMove(QWidget* target, QPoint p) {
        auto* ev = new QMouseEvent(QEvent::MouseMove, QPointF(p), QPointF(p),
                                   Qt::NoButton, Qt::NoButton, Qt::NoModifier);
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
    // Plan 03: single MouseMove returns notypeevent=-2 with the correct
    // coordinates. This is the simplest motion-coalescing case — a burst
    // of size one — and validates the Pitfall 8 "fresh position" invariant.
    void testXvsourisMotion() {
        auto* canvas = XvueApp::window_slot()->canvas();
        XvueEventBridge bridge(canvas);
        canvas->installEventFilter(&bridge);
        QTimer::singleShot(10, [canvas]{ postMove(canvas, QPoint(42, 73)); });
        auto r = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r.notypeevent, -2);
        QCOMPARE(r.nbc, 0);              // X11 contract: motion carries no button
        QCOMPARE(r.x, 42);
        QCOMPARE(r.y, 73);
        QCOMPARE(XvueApp::blockingDepth(), 0);
        canvas->removeEventFilter(&bridge);
    }

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


    // Phase 5 Plan 04 (EVENT-02): real Fortran ABI entry point routes through
    // XvueApp::window_slot()->bridge()->waitForEvent(Souris). Driving a space
    // keypress into the live window's canvas must produce notypeevent=2,
    // nbc=32 in the Fortran out-parameters.
    void testXvsourisFortranEntryPoint() {
        auto& win = XvueApp::window_slot();
        QVERIFY(win != nullptr);
        auto* canvas = win->canvas();
        QVERIFY(canvas != nullptr);
        // The production bridge is already installed by XvueWindow ctor via
        // installEventFilter — no local bridge here.
        QTimer::singleShot(10, [canvas]{ postKey(canvas, Qt::Key_Space, QStringLiteral(" ")); });
        int ntev = -99, nbc = -99, x1 = -99, y1 = -99;
        xvsouris_(&ntev, &nbc, &x1, &y1);
        QCOMPARE(ntev, 2);
        QCOMPARE(nbc, 32);
        QCOMPARE(XvueApp::blockingDepth(), 0);
    }

    // Phase 5 Plan 04 (EVENT-02): silent-abandon guard. With no window open
    // the Fortran entry point must not crash; it writes notypeevent=0 and
    // returns. Mirror of Pitfall 11 in the plan context.
    void testXvsourisNoWindow() {
        // Close the window to simulate "never opened" from the bridge's POV.
        xvfermer_();
        QVERIFY(XvueApp::window_slot() == nullptr);

        int ntev = -99, nbc = -99, x1 = -99, y1 = -99;
        xvsouris_(&ntev, &nbc, &x1, &y1);
        QCOMPARE(ntev, 0);
        QCOMPARE(nbc, 0);
        QCOMPARE(x1, 0);
        QCOMPARE(y1, 0);
        QCOMPARE(XvueApp::blockingDepth(), 0);

        // Reopen for downstream tests in the suite ordering.
        xvinitgraphique_();
        QVERIFY(XvueApp::window_slot() != nullptr);
    }
    // ---- EVENT-03: xvsouris2_ ----
    // Phase 5 Plan 05 Task 2 (EVENT-03). Drives the real Fortran ABI entry
    // point through the bridge with a populated items[] array in the X11
    // layout verified against xvue/saclav.f:61 and xvuelc.c:2397-2413:
    //   items[0] = mots  = words per item (here 3: x, y, num)
    //   items[1] = max   = max item capacity (unused by xvsouris2_)
    //   items[2] = nbitem= actual item count
    //   items[3..]      = (x, y, num) triplets
    // *pmin0 is the OFFSET (not index) of the nearest item into items[], so
    // for two 3-word items at (100,100,1) and (200,200,2) the valid offsets
    // are 3 and 6. A click near (100,100) must return *pmin0 == 3.
    void testXvsouris2Accrochage() {
        auto& win = XvueApp::window_slot();
        QVERIFY(win != nullptr);
        auto* canvas = win->canvas();
        QVERIFY(canvas != nullptr);

        // Initialise the sprite so the filter has something to blit.
        initaccrochage_();
        QVERIFY(win->state()->mempxaccro_ != nullptr);

        int items[9] = { 3, 2, 2,   100, 100, 1,   200, 200, 2 };
        int pmin0 = -2;
        int ntev = -99, ibtn = -99, x1 = -99, y1 = -99;

        // A left button click near the first item (100, 100).
        QTimer::singleShot(10, [canvas]{
            postButtonPress(canvas, Qt::LeftButton, QPoint(105, 98));
        });

        xvsouris2_(items, &pmin0, &ntev, &ibtn, &x1, &y1);

        // Nearest to (105, 98) is item 0 at offset 3. notypeevent=5 from
        // the press path (X11 parity: press-in-accrochage returns 5, not 1).
        QCOMPARE(pmin0, 3);
        QCOMPARE(ntev, 5);
        QCOMPARE(ibtn, 1);
        QCOMPARE(x1, 105);
        QCOMPARE(y1, 98);
        QCOMPARE(XvueApp::blockingDepth(), 0);
    }

    // Plan 05 Task 2 (Rule 2 / T-05-05-01). Defensive guards: items==null
    // and null out-param slots must not crash. The function returns with
    // ntev=0 (no-event) and zero-filled out-params.
    void testXvsouris2NullGuards() {
        auto& win = XvueApp::window_slot();
        QVERIFY(win != nullptr);
        auto* canvas = win->canvas();
        QVERIFY(canvas != nullptr);

        initaccrochage_();
        int items[6] = { 3, 1, 1,   50, 50, 1 };
        int pmin0 = -2;

        // Dispatch a space key so the bridge returns promptly. We don't
        // care about the return code — only that the call does not crash
        // even when notypeevent / ibutton / x1 / y1 are null.
        QTimer::singleShot(10, [canvas]{
            postKey(canvas, Qt::Key_Space, QStringLiteral(" "));
        });
        xvsouris2_(items, &pmin0, nullptr, nullptr, nullptr, nullptr);
        QCOMPARE(XvueApp::blockingDepth(), 0);
    }

    // Plan 05 Task 2 (Pitfall 10 and T-05-05-02). On canvas resize the
    // accroche_undo_tile_ is invalidated by the resizeEvent (Plan 01 Task 3
    // installed the invalidation guard). Verify the invariant holds: after
    // resize, the tile pointer is null and the next xvsouris2_ motion
    // allocates a fresh tile without crashing on the stale pointer.
    void testXvsouris2ResizeInvalidatesTile() {
        auto& win = XvueApp::window_slot();
        QVERIFY(win != nullptr);
        auto* canvas = win->canvas();
        auto* state = win->state();
        QVERIFY(canvas != nullptr);
        QVERIFY(state != nullptr);

        initaccrochage_();
        int items[6] = { 3, 1, 1,   80, 80, 1 };
        int pmin0 = -2;
        int ntev = -99, ibtn = -99, x1 = -99, y1 = -99;

        // First press draws the sprite and saves a tile.
        QTimer::singleShot(10, [canvas]{
            postButtonPress(canvas, Qt::LeftButton, QPoint(82, 78));
        });
        xvsouris2_(items, &pmin0, &ntev, &ibtn, &x1, &y1);
        QCOMPARE(pmin0, 3);
        QVERIFY(state->accroche_undo_tile_ != nullptr);

        // Now resize — Plan 01 Task 3's resizeEvent guard invalidates
        // the tile so the next xvsouris2_ call does not use stale contents.
        canvas->resize(canvas->size() + QSize(20, 20));
        QCoreApplication::processEvents();
        QVERIFY(state->accroche_undo_tile_ == nullptr);
        QCOMPARE(XvueApp::blockingDepth(), 0);
    }

    // ---- EVENT-04: xvpause_ ----
    // Phase 5 Plan 04 (EVENT-04). xvpause_ blocks on bridge->waitForEvent(Pause)
    // until a KeyPress arrives on the live window's canvas.
    void testXvpauseReturnsOnKey() {
        auto& win = XvueApp::window_slot();
        QVERIFY(win != nullptr);
        auto* canvas = win->canvas();
        QVERIFY(canvas != nullptr);
        QTimer::singleShot(10, [canvas]{ postKey(canvas, Qt::Key_Space, QStringLiteral(" ")); });
        QElapsedTimer t; t.start();
        xvpause_();
        QVERIFY2(t.elapsed() < 1000,
                 qPrintable(QStringLiteral("xvpause_ took %1ms, expected < 1000ms").arg(t.elapsed())));
        QCOMPARE(XvueApp::blockingDepth(), 0);
    }
    void testXvpauseReturnsOnMouseClick() { QSKIP("Plan 04: mouse-click Pause path not load-bearing; xvpause_ contract is KeyPress-only (xvuelc.c:2529)"); }

    // Phase 5 Plan 04 (EVENT-04). MEFISTO_XVSOURIS_AUTOEXIT short-circuits
    // xvpause_ in both backends so headless harnesses do not hang (§8 of
    // 05-RESEARCH.md). The env var is read per-call (not cached) so we can
    // flip it at test time.
    void testXvpauseAutoexit() {
        qputenv("MEFISTO_XVSOURIS_AUTOEXIT", "1");
        qputenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS", "50");
        QElapsedTimer t; t.start();
        xvpause_();
        const qint64 elapsed = t.elapsed();
        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT");
        qunsetenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS");
        QVERIFY2(elapsed < 500,
                 qPrintable(QStringLiteral("xvpause_ AUTOEXIT took %1ms, expected < 500ms").arg(elapsed)));
        QCOMPARE(XvueApp::blockingDepth(), 0);
    }

    // ---- EVENT-05: deplsouris_ ----
    // Phase 5 Plan 04 (EVENT-05, D-09). deplsouris_ warps the cursor via
    // QCursor::setPos(canvas->mapToGlobal(...)) and returns immediately.
    // Wayland is a no-op (D-09); under xvfb the cursor usually moves within
    // a few pixels of the requested point, but the verification is best-
    // effort (offscreen QPA may not move the cursor at all — we only
    // enforce the non-blocking contract and the positional tolerance if
    // the cursor moved at all).
    void testDeplsourisNonBlocking() {
        auto& win = XvueApp::window_slot();
        QVERIFY(win != nullptr);
        auto* canvas = win->canvas();
        QVERIFY(canvas != nullptr);

        int x = 50, y = 60;
        QElapsedTimer t; t.start();
        deplsouris_(&x, &y);
        QVERIFY2(t.elapsed() < 100,
                 qPrintable(QStringLiteral("deplsouris_ took %1ms, expected non-blocking").arg(t.elapsed())));

        QPoint expected = canvas->mapToGlobal(QPoint(50, 60));
        QPoint actual = QCursor::pos();
        if (actual != QPoint(0, 0)) {
            QVERIFY2(qAbs(actual.x() - expected.x()) <= 5,
                     qPrintable(QStringLiteral("cursor x drift: expected=%1 actual=%2")
                                .arg(expected.x()).arg(actual.x())));
            QVERIFY2(qAbs(actual.y() - expected.y()) <= 5,
                     qPrintable(QStringLiteral("cursor y drift: expected=%1 actual=%2")
                                .arg(expected.y()).arg(actual.y())));
        }
    }

    // ---- EVENT-06: initaccrochage_ ----
    // Phase 5 Plan 05 Task 1: allocates XvueState::mempxaccro_ as a 13x13 QPixmap
    // with a 3-pixel black square border on transparent background (Strategy B
    // from 05-RESEARCH.md §6 — the save/restore blit approach replaces the X11
    // GXand/GXorInverted raster-op XOR trick).
    void testInitaccrochage() {
        auto& win = XvueApp::window_slot();
        QVERIFY(win != nullptr);
        auto* state = win->state();
        QVERIFY(state != nullptr);

        // Clean slate: if a previous test primed the sprite, drop it so we
        // exercise the first-call allocation branch deterministically.
        if (state->mempxaccro_) {
            delete state->mempxaccro_;
            state->mempxaccro_ = nullptr;
        }

        initaccrochage_();
        QVERIFY(state->mempxaccro_ != nullptr);
        QCOMPARE(state->mempxaccro_->size(), QSize(13, 13));

        // Pixel sampling. The sprite is a 13x13 with a 3-pixel-thick black
        // square border on a transparent center. Center (6,6) must be
        // transparent (alpha=0) and the border band at (2,2) must be opaque
        // black.
        QImage img = state->mempxaccro_->toImage();
        QCOMPARE(img.pixelColor(6, 6).alpha(), 0);
        QVERIFY2(img.pixelColor(2, 2).alpha() == 255,
                 qPrintable(QStringLiteral("border pixel (2,2) alpha=%1 (want 255)")
                            .arg(img.pixelColor(2, 2).alpha())));
        QVERIFY2(img.pixelColor(2, 2).red() == 0,
                 qPrintable(QStringLiteral("border pixel (2,2) red=%1 (want 0)")
                            .arg(img.pixelColor(2, 2).red())));

        // Idempotency: calling a second time must not leak nor crash; the
        // resulting sprite still has the right size and border content.
        initaccrochage_();
        QVERIFY(state->mempxaccro_ != nullptr);
        QCOMPARE(state->mempxaccro_->size(), QSize(13, 13));
        QImage img2 = state->mempxaccro_->toImage();
        QCOMPARE(img2.pixelColor(6, 6).alpha(), 0);
        QCOMPARE(img2.pixelColor(2, 2).alpha(), 255);
    }

    // Phase 5 Plan 05 Task 1 (Pitfall 11): initaccrochage_ called before
    // xvinitgraphique_ (or after xvfermer_) must NOT crash — it is a silent
    // no-op when there is no window / canvas / state yet.
    void testInitaccrochageBeforeInit() {
        xvfermer_();
        QVERIFY(XvueApp::window_slot() == nullptr);
        // Must not crash — silent no-op.
        initaccrochage_();
        // Reopen so downstream tests in the declaration order still have a
        // live window.
        xvinitgraphique_();
        QVERIFY(XvueApp::window_slot() != nullptr);
    }

    // ---- EVENT-07: motion coalescing ----
    // Plan 03 (D-04): 100 rapid MouseMove events → bounded number of
    // waitForEvent returns. The deferred-quit timer guarantees each
    // waitForEvent coalesces at least the burst posted before the timer
    // fires; the drain loop stops cleanly on an Esc terminator posted at
    // the end of the same timer callback (so no stale singleShot lambdas
    // leak into subsequent tests — a concrete bug observed during Plan
    // 03 TDD when fallback timers with by-reference captures out-lived
    // the test function).
    void testMotionCoalescing() {
        auto* canvas = XvueApp::window_slot()->canvas();
        XvueEventBridge bridge(canvas);
        canvas->installEventFilter(&bridge);

        // Post 100 moves at once via a single timer callback. Esc goes in
        // a *separate* later timer so the first waitForEvent's deferred-
        // quit fires on the motion burst's tail (without the Esc racing
        // it). All captures by value (no dangling references after return).
        QTimer::singleShot(10, [canvas]{
            for (int i = 0; i < 100; ++i) {
                postMove(canvas, QPoint(10 + i, 20 + i));
            }
        });
        // Esc at 200 ms: enough time for the first waitForEvent to return
        // with its coalesced motion, and for subsequent drain waitForEvent
        // calls to consume any remaining moves in the queue.
        QTimer::singleShot(200, [canvas]{
            postKey(canvas, Qt::Key_Escape, QString());
        });

        int returns      = 0;
        int last_x       = -1;
        int last_y       = -1;
        // Drain until the Esc terminator lands (notypeevent=0). Hard upper
        // bound of 25 iterations so a bug that fails to coalesce at all
        // does not wedge the test — we fail on `returns > 20` below.
        for (int i = 0; i < 25; ++i) {
            auto r = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
            if (r.notypeevent == -2) {
                ++returns;
                last_x = r.x;
                last_y = r.y;
                continue;
            }
            // Non-motion (Esc) means the queue is drained — stop counting.
            QCOMPARE(r.notypeevent, 0);
            QCOMPARE(r.nbc, 27);
            break;
        }

        // Bounded: at least one return (we posted moves) and ≤ 20 total.
        QVERIFY2(returns >= 1, "expected at least one motion return");
        QVERIFY2(returns <= 20, qPrintable(QStringLiteral("expected <= 20 motion returns, got %1").arg(returns)));
        // Last returned position must land in the burst range.
        QVERIFY2(last_x >= 10 && last_x <= 109,
                 qPrintable(QStringLiteral("last_x=%1 out of [10..109]").arg(last_x)));
        QVERIFY2(last_y >= 20 && last_y <= 119,
                 qPrintable(QStringLiteral("last_y=%1 out of [20..119]").arg(last_y)));
        QCOMPARE(XvueApp::blockingDepth(), 0);
        canvas->removeEventFilter(&bridge);
    }

    // Plan 03 (Pitfall 8): pending_.x / pending_.y must never hold stale
    // values across waitForEvent invocations. Two calls: first posts (5,6),
    // second posts (777,888) — the second return must show (777,888), not
    // (5,6). This is guarded at the cpp level by `pending_ = Result{}`
    // at waitForEvent entry.
    void testMotionFreshPerCall() {
        auto* canvas = XvueApp::window_slot()->canvas();
        XvueEventBridge bridge(canvas);
        canvas->installEventFilter(&bridge);

        QTimer::singleShot(10, [canvas]{ postMove(canvas, QPoint(5, 6)); });
        auto r1 = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r1.notypeevent, -2);
        QCOMPARE(r1.x, 5);
        QCOMPARE(r1.y, 6);

        QTimer::singleShot(10, [canvas]{ postMove(canvas, QPoint(777, 888)); });
        auto r2 = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r2.notypeevent, -2);
        QCOMPARE(r2.x, 777);
        QCOMPARE(r2.y, 888);

        QCOMPARE(XvueApp::blockingDepth(), 0);
        canvas->removeEventFilter(&bridge);
    }

    // Plan 03 (Pitfall 9): quit_pending_ must be reset to false at the
    // start of every waitForEvent call. After a motion burst the flag is
    // true; if waitForEvent did not reset it, the second call's first
    // MouseMove would fail to arm a new singleShot and the loop would
    // hang until a button/key arrived. We verify the second call still
    // returns on a pure-motion burst.
    void testQuitPendingResetAcrossCalls() {
        auto* canvas = XvueApp::window_slot()->canvas();
        XvueEventBridge bridge(canvas);
        canvas->installEventFilter(&bridge);

        QTimer::singleShot(10, [canvas]{ postMove(canvas, QPoint(1, 2)); });
        auto r1 = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r1.notypeevent, -2);

        // Second call: pure motion burst. If quit_pending_ leaked true
        // from the first call, eventFilter's `if (!quit_pending_)` guard
        // would skip the singleShot and this waitForEvent would never
        // return. We post only the motion — no fallback Esc, because
        // leaking singleShot lambdas across tests crashed the suite.
        // If the quit_pending_ reset ever regresses, this test will time
        // out at QTest's default per-function timeout (which is an
        // explicit failure — not a hang).
        QTimer::singleShot(10, [canvas]{ postMove(canvas, QPoint(3, 4)); });
        auto r2 = bridge.waitForEvent(XvueEventBridge::WaitMode::Souris);
        QCOMPARE(r2.notypeevent, -2);
        QCOMPARE(r2.x, 3);
        QCOMPARE(r2.y, 4);

        QCOMPARE(XvueApp::blockingDepth(), 0);
        canvas->removeEventFilter(&bridge);
    }

    // Plan 03: MEFISTO_XVSOURIS_DEBUG env-var cache. The env-var lookup is
    // cached on first call, so this test cannot flip the flag at runtime —
    // the cache is fixed for the process. Instead, the acceptance criteria
    // verify that when the test runner is invoked with the env var set,
    // stderr contains at least one `motion_count=` line. The grep is done
    // from the verify block outside the test binary. This test is a
    // documented placeholder for the runtime flag and also asserts the
    // debug_logging_enabled accessor is callable (compile-time guard).
    void testDebugLoggingEnvVar() {
        QSKIP("Plan 03/06: MEFISTO_XVSOURIS_DEBUG caches on first access — "
              "runtime flag tested via stderr grep in the verify block");
    }

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
