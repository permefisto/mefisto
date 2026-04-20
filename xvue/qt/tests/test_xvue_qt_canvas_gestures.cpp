// xvue/qt/tests/test_xvue_qt_canvas_gestures.cpp
// Phase 6.0 Plan 05 Task 2 (UX-12, UX-13, UI-SPEC §Canvas, Flag #3):
// 10 QTest cases covering the new XvueCanvas event handlers + view transform
// + empty-state paint branch.
//
// Test strategy: instantiate XvueCanvas directly with a stack-allocated
// XvueState (no XvueWindow parent needed for wheel/pan/paintEvent/resetView).
// This keeps the suite isolated from Plan 03's XvueEventBridge and Plan 02's
// XvueWindow::menuBridge() accessor — Plan 05 only verifies the canvas-side
// behavior. The "menu populated by XvueMenuBridge" branch is covered by
// Plan 04's tests + the Plan 07 manual sweep.
//
// Custom main() pattern follows test_xvue_qt_event.cpp / test_xvue_qt_menu_*:
// XvueApp::ensure() constructs the QApplication BEFORE QTest::qExec so the
// AA_CompressHighFrequencyEvents attribute set inside ensure() is visible
// (D-05 / Pitfall 7 from Phase 5).
#include "test_helpers.h"
#include "xvue_qt_canvas.h"
#include "xvue_qt_state.h"
#include "xvue_qt_app.h"
#include "xvue_qt_event.h"   // BlockingDepthGuard

#include <QtTest/QtTest>
#include <QApplication>
#include <QCoreApplication>
#include <QEventLoop>
#include <QWheelEvent>
#include <QMouseEvent>
#include <QContextMenuEvent>
#include <QPixmap>
#include <QPainter>
#include <QPoint>
#include <QPointF>
#include <QTransform>

class TestXvueQtCanvasGestures : public QObject {
    Q_OBJECT

    // ---- helpers (mirror test_xvue_qt_event.cpp's postKey / postMove pattern)
    //
    // Use sendEvent (synchronous) rather than postEvent + processEvents, because
    // processEvents(ExcludeUserInputEvents) drops the wheel/mouse events the
    // tests are trying to deliver — that was the root cause of the flake on
    // testWheelZoomIn documented in deferred-items.md. sendEvent guarantees
    // synchronous delivery via QCoreApplication::notify, bypassing the queue
    // and the user-input filter, which makes assertions deterministic.
    static void postWheel(QWidget* target, int angleDeltaY) {
        QWheelEvent ev(
            QPointF(10, 10), QPointF(10, 10),
            QPoint(),          // pixelDelta — unused for wheel-step model
            QPoint(0, angleDeltaY),
            Qt::NoButton, Qt::NoModifier,
            Qt::NoScrollPhase, /*inverted*/ false);
        QCoreApplication::sendEvent(target, &ev);
    }
    static void postMousePress(QWidget* target, Qt::MouseButton b, QPoint p) {
        QMouseEvent ev(QEvent::MouseButtonPress, QPointF(p), QPointF(p),
                       b, b, Qt::NoModifier);
        QCoreApplication::sendEvent(target, &ev);
    }
    static void postMouseMove(QWidget* target, Qt::MouseButtons held, QPoint p) {
        QMouseEvent ev(QEvent::MouseMove, QPointF(p), QPointF(p),
                       Qt::NoButton, held, Qt::NoModifier);
        QCoreApplication::sendEvent(target, &ev);
    }
    static void postMouseRelease(QWidget* target, Qt::MouseButton b, QPoint p) {
        QMouseEvent ev(QEvent::MouseButtonRelease, QPointF(p), QPointF(p),
                       b, Qt::NoButton, Qt::NoModifier);
        QCoreApplication::sendEvent(target, &ev);
    }
    static void postContextMenu(QWidget* target, QPoint p) {
        QContextMenuEvent ev(QContextMenuEvent::Mouse, p, p);
        QCoreApplication::sendEvent(target, &ev);
    }
    static void drainEvents() {
        // Kept for backward compatibility — sendEvent is synchronous, but a
        // few tests still rely on a deleteLater drain between iterations.
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }

    // Stack-allocated state + heap-allocated canvas (owned by us). The canvas
    // is created in initTestCase and destroyed in cleanupTestCase. Plan 01
    // (D-08): heap-allocated top-level QWidgets registered with qApp's
    // topLevelWidgets() outlive their delete and can crash at exit. The
    // canvas here has parent=nullptr so it IS a top-level widget — but we
    // explicitly destroy it BEFORE QApplication's atexit (cleanupTestCase
    // runs before main returns) and we drain deleteLater between tests.
    XvueState state_;
    XvueCanvas* canvas_ = nullptr;

private slots:
    void initTestCase() {
        canvas_ = new XvueCanvas(&state_);
        canvas_->resize(400, 300);
        // Allocate a backing pixmap so paintEvent exercises the drawPixmap
        // branch (the empty-state branch runs after, regardless).
        state_.backing_ = new QPixmap(400, 300);
        state_.backing_->fill(Qt::black);
    }

    void cleanupTestCase() {
        delete canvas_;
        canvas_ = nullptr;
        // state_'s dtor cleans backing_.
        // Drain deleteLater so any queued widget destruction completes before
        // QApplication shutdown (mirrors test_xvue_qt_menu_scaffold.cpp).
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }

    void cleanup() {
        // Drain deleteLater between tests too.
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }

    void init() {
        // Reset per-test mutable state on the shared canvas + XvueState.
        state_.view_transform_         = QTransform();
        state_.has_user_content_       = false;
        canvas_->contextMenuShownCount_   = 0;
        canvas_->lastPaintDrewEmptyState_ = false;
    }

    // ---- Wheel zoom tests ----
    void testWheelZoomIn() {
        postWheel(canvas_, 120);
        drainEvents();
        QVERIFY(state_.view_transform_.m11() > 1.0);
        // Single notch -> exactly 1.15 (within fuzzy tolerance).
        QVERIFY(qFuzzyCompare(state_.view_transform_.m11(), 1.15));
    }

    void testWheelZoomOut() {
        postWheel(canvas_, -120);
        drainEvents();
        QVERIFY(state_.view_transform_.m11() < 1.0);
        QVERIFY(state_.view_transform_.m11() > 0.0);
    }

    void testWheelClamping() {
        // 30 zoom-in steps would be 1.15^30 ≈ 66.2x without the clamp;
        // m11 must stay <= 10.0.
        for (int i = 0; i < 30; ++i) postWheel(canvas_, 120);
        drainEvents();
        QVERIFY(state_.view_transform_.m11() <= 10.0);
        // Then zoom out 60 steps to push past 0.1; m11 must stay >= 0.1.
        for (int i = 0; i < 60; ++i) postWheel(canvas_, -120);
        drainEvents();
        QVERIFY(state_.view_transform_.m11() >= 0.1);
    }

    // ---- Middle-drag pan ----
    void testMiddleDragPan() {
        postMousePress(canvas_, Qt::MiddleButton, QPoint(10, 10));
        postMouseMove(canvas_, Qt::MiddleButton, QPoint(15, 20));
        postMouseRelease(canvas_, Qt::MiddleButton, QPoint(15, 20));
        drainEvents();
        // Drag delta = (5, 10) in widget coordinates; anchor was identity
        // so final transform should be a pure translate.
        QCOMPARE(state_.view_transform_.dx(), 5.0);
        QCOMPARE(state_.view_transform_.dy(), 10.0);
    }

    // Verify pan-press is accepted (i.e., NOT propagated to default handler).
    // We don't have a way to inspect QMouseEvent::isAccepted() from the
    // postEvent path, so we use an indirect proof: a press-only (no move,
    // no release) leaves pan_active_ = true, but state_.view_transform_
    // remains identity. A subsequent move with NO MMB held does NOT update
    // the transform either (no buttons & MiddleButton).
    void testMiddleDragNoFortranDispatch() {
        postMousePress(canvas_, Qt::MiddleButton, QPoint(10, 10));
        drainEvents();
        // After press alone, the view transform is still identity (no drag).
        QVERIFY(state_.view_transform_.isIdentity());
        // Post a move with NO buttons held — must be ignored by the pan path.
        postMouseMove(canvas_, Qt::NoButton, QPoint(50, 50));
        drainEvents();
        QVERIFY(state_.view_transform_.isIdentity());
        // Cleanup: release the middle button so subsequent tests start clean.
        postMouseRelease(canvas_, Qt::MiddleButton, QPoint(10, 10));
        drainEvents();
    }

    // ---- Context menu tests ----
    void testContextMenuSuppressedBlocking() {
        {
            BlockingDepthGuard guard;   // bumps blockingDepth_ to 1
            postContextMenu(canvas_, QPoint(50, 50));
            drainEvents();
        }   // guard goes out of scope -> blockingDepth_ back to 0
        QCOMPARE(canvas_->contextMenuShownCount_, 0);
    }

    void testContextMenuShownWhenIdle() {
        QCOMPARE(XvueApp::blockingDepth(), 0);
        postContextMenu(canvas_, QPoint(50, 50));
        drainEvents();
        // Counter bumps even though the menu is empty (no populator wired
        // in Plan 05 — Plan 02 wires XvueMenuBridge; Plan 04 covers the
        // populated-menu path).
        QCOMPARE(canvas_->contextMenuShownCount_, 1);
    }

    // ---- resetView ----
    void testResetViewIdentity() {
        postWheel(canvas_, 120);
        drainEvents();
        QVERIFY(!state_.view_transform_.isIdentity());
        canvas_->resetView();
        QVERIFY(state_.view_transform_.isIdentity());
    }

    // ---- Empty-state ----
    //
    // QWidget::repaint() is a no-op on hidden widgets (and the canvas here is
    // a top-level widget that's never shown), which was the root cause of the
    // testEmptyStateRendersText flake documented in deferred-items.md. Use
    // QWidget::grab() instead — it forces a full paintEvent into an offscreen
    // pixmap regardless of visibility, making the test deterministic without
    // requiring the widget to be shown (cleaner under offscreen QPA).
    void testEmptyStateRendersText() {
        state_.has_user_content_ = false;
        (void)canvas_->grab();   // forces synchronous paintEvent on hidden widget
        QVERIFY(canvas_->lastPaintDrewEmptyState_);
    }

    void testEmptyStateClearAfterContent() {
        state_.has_user_content_ = true;
        (void)canvas_->grab();
        QVERIFY(!canvas_->lastPaintDrewEmptyState_);
    }
};

int main(int argc, char* argv[]) {
    (void)argc; (void)argv;
    XvueApp::ensure();   // constructs the process QApplication (D-05/Pitfall 7)
    TestXvueQtCanvasGestures tc;
    int   qt_argc   = 1;
    char  qt_arg0[] = "xvue_qt_canvas_gestures_tests";
    char* qt_argv[] = { qt_arg0, nullptr };
    return QTest::qExec(&tc, qt_argc, qt_argv);
}

#include "test_xvue_qt_canvas_gestures.moc"
