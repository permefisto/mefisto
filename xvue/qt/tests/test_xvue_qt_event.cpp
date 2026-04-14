// xvue/qt/tests/test_xvue_qt_event.cpp
// Phase 5 Wave 0 (05-01 Task 1): QTest skeleton covering EVENT-01..EVENT-08.
// All event-bridge test bodies currently QSKIP; each downstream plan replaces
// a subset with real assertions as the corresponding bridge behavior lands.
#include "test_helpers.h"

#include <QtTest/QtTest>
#include <QApplication>

class TestXvueQtEvent : public QObject {
    Q_OBJECT
private slots:
    void initTestCase() {
        // XvueApp::ensure() will allocate the QApplication if needed — the
        // first helper call that touches XvueApp does so.
    }

    // ---- EVENT-01: bridge installation ----
    void testBridgeInstallation()         { QSKIP("Plan 02: XvueEventBridge not yet implemented"); }

    // ---- EVENT-02: xvsouris_ ----
    void testXvsourisMotion()             { QSKIP("Plan 02/03: motion not yet wired"); }
    void testXvsourisButtonPress()        { QSKIP("Plan 02: xvsouris_ not yet wired"); }
    void testXvsourisButtonRelease()      { QSKIP("Plan 02: xvsouris_ not yet wired"); }
    void testXvsourisKeyPress()           { QSKIP("Plan 02: xvsouris_ not yet wired"); }
    void testXvsourisEscapeAbort()        { QSKIP("Plan 04: xvsouris_ body not yet"); }
    void testXvsourisAtSignAbort()        { QSKIP("Plan 04: xvsouris_ body not yet"); }

    // ---- EVENT-03: xvsouris2_ ----
    void testXvsouris2Accrochage()        { QSKIP("Plan 05: xvsouris2_ not yet wired"); }

    // ---- EVENT-04: xvpause_ ----
    void testXvpauseReturnsOnKey()        { QSKIP("Plan 04: xvpause_ not yet wired"); }
    void testXvpauseReturnsOnMouseClick() { QSKIP("Plan 04: xvpause_ not yet wired"); }
    void testXvpauseAutoexit()            { QSKIP("Plan 04: AUTOEXIT extension pending"); }

    // ---- EVENT-05: deplsouris_ ----
    void testDeplsourisNonBlocking()      { QSKIP("Plan 04: deplsouris_ not yet wired"); }

    // ---- EVENT-06: initaccrochage_ ----
    void testInitaccrochage()             { QSKIP("Plan 05: initaccrochage_ not yet wired"); }

    // ---- EVENT-07: motion coalescing ----
    void testMotionCoalescing()           { QSKIP("Plan 03: motion coalescing pending"); }

    // ---- EVENT-08: blocking depth ----
    void testBlockingDepthRAII()          { QSKIP("Plan 02: RAII guard pending"); }
    void testBlockingDepthNested()        { QSKIP("Plan 02: nested guard pending"); }

    // Smoke test — Wave 0 Task 2 lands the probe symbol. Before Task 2 this
    // test fails at link time (unresolved symbol), so Task 1 acceptance
    // tolerates a single known link/runtime failure here.
    void testBlockingDepthAccessorZero() {
        extern int xvue_qt_test_blocking_depth_probe();  // forward decl
        QCOMPARE(xvue_qt_test_blocking_depth_probe(), 0);
    }
};

// Provide a weak probe symbol so Task 1 can link the test binary even before
// XvueApp::blockingDepth() lands in Task 2. Task 2 replaces this file's body
// with a real #include "xvue_qt_app.h" call.
int xvue_qt_test_blocking_depth_probe() { return 0; }

QTEST_MAIN(TestXvueQtEvent)
#include "test_xvue_qt_event.moc"
