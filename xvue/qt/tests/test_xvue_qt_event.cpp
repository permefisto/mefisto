// xvue/qt/tests/test_xvue_qt_event.cpp
// Phase 5 Wave 0 (05-01 Task 1): QTest skeleton covering EVENT-01..EVENT-08.
// All event-bridge test bodies currently QSKIP; each downstream plan replaces
// a subset with real assertions as the corresponding bridge behavior lands.
#include "test_helpers.h"
#include "xvue_qt_app.h"
#include "xvue_qt_state.h"

#include <QtTest/QtTest>
#include <QApplication>
#include <QCoreApplication>
#include <QPixmap>

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
        QSKIP("Plan 02 validates via full canvas construction");
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
