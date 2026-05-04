// xvue/qt/tests/test_xvue_qt_postscript.cpp
// Phase 7 Plan 02 (EXPORT-04): unit tests for PsEmitter::handleLasops state
// machine + the pyFlip helper that Plan 03's per-primitive emit methods will
// rely on. Six slots cover the full xvuelc.c:1187-1304 dispatch matrix
// (open/close, mode-100 erase, chaine dispatch, negative-erase,
// modepsc==0 chaine guard, Y-flip helper).
//
// All slots run under QT_QPA_PLATFORM=offscreen via xvfb-run --auto-servernum
// (per ctest invocation in xvue/qt/tests/CMakeLists.txt) so the QApplication
// constructed by XvueApp::ensure() does not require a real X display.
//
// Each test slot chdirs into a per-test QTemporaryDir BEFORE calling
// handleLasops so the legacy `TEMPORAIRE.EPS` filename writes to a scratch
// directory, not the project tree. The dir is restored in cleanup() so the
// next slot starts clean.
#include "xvue_qt_app.h"
#include "xvue_qt_postscript.h"

#include <QtTest/QtTest>
#include <QCoreApplication>
#include <QDir>
#include <QEventLoop>
#include <QTemporaryDir>

class TestXvueQtPostscript : public QObject {
    Q_OBJECT

private:
    QString origCwd_;
    QTemporaryDir tmpDir_;

private slots:
    void initTestCase() {
        // Capture original cwd so we can restore in cleanup() — every slot
        // chdirs into its own QTemporaryDir to keep TEMPORAIRE.EPS writes
        // out of the project tree.
        origCwd_ = QDir::currentPath();
        QVERIFY(tmpDir_.isValid());
    }

    void cleanupTestCase() {
        QDir::setCurrent(origCwd_);
    }

    void cleanup() {
        // Drain Qt deleteLater between tests so Qt children do not race the
        // leaked-QApplication atexit (D-08 in xvue_qt_app.cpp).
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
        QDir::setCurrent(origCwd_);
    }

    void init() {
        // Each test gets its own scratch subdirectory so TEMPORAIRE.EPS
        // writes never collide between slots. QTest::currentTestFunction()
        // gives us a unique-per-slot name.
        QString sub = tmpDir_.path() + "/" +
                      QString::fromUtf8(QTest::currentTestFunction());
        QDir().mkpath(sub);
        QVERIFY(QDir::setCurrent(sub));
    }

    // Test 1 — state machine, lasopsc==0 close path.
    // Construct PsEmitter, call handleLasops(1) — fpo opens, lasopsc_==1.
    // Then call handleLasops(0) — fpo flushed and closed, concat_[0]=='\0'.
    void PsEmitter_handleLasops_open_close() {
        PsEmitter pe;
        QCOMPARE(pe.lasopsc(), 0);
        QVERIFY(pe.fpoForTesting() == nullptr);

        pe.handleLasops(1);
        QCOMPARE(pe.lasopsc(), 1);
        QVERIFY(pe.fpoForTesting() != nullptr);
        // TEMPORAIRE.EPS now exists in cwd.
        QVERIFY(QFile::exists("TEMPORAIRE.EPS"));

        pe.handleLasops(0);
        QCOMPARE(pe.lasopsc(), 0);
        QVERIFY(pe.fpoForTesting() == nullptr);
        // concat is reset to empty string after close.
        QCOMPARE(pe.concatForTesting()[0], '\0');
    }

    // Test 2 — state machine, mode 100 erase reopens TEMPORAIRE.EPS and
    // resets the iTe/iFa/ity/iep/iPo/ire/iRe/iel/iEl/iFP macro counters.
    // handleLasops(1), then handleLasops(101) — fpo_ closed/reopened;
    // lasopsc_ ends at 1 (101-100).
    void PsEmitter_handleLasops_mode100_reset() {
        PsEmitter pe;
        pe.handleLasops(1);
        QCOMPARE(pe.lasopsc(), 1);
        FILE* before = pe.fpoForTesting();
        QVERIFY(before != nullptr);

        pe.handleLasops(101);
        // After mode-100 fall-through, lasopsc_ == 101 - 100 == 1.
        QCOMPARE(pe.lasopsc(), 1);
        // The fpo_ pointer was closed and reopened — should be non-null
        // again (file rewritten under the same name).
        QVERIFY(pe.fpoForTesting() != nullptr);
        QVERIFY(QFile::exists("TEMPORAIRE.EPS"));
    }

    // Test 3 — chaine dispatch, lasopsc==4. After handleLasops(4) under
    // modepsc_==1, chaine_[0] is allocated by the calloc loop (since
    // modepsc_!=0). lasopsc_==4 selects the chaine path for downstream
    // emit helpers (Plan 03 territory; Plan 02 only verifies the
    // dispatch field).
    void PsEmitter_handleLasops_chaine_dispatch() {
        PsEmitter pe;
        pe.setModepscForTesting(1);   // enable the chaine calloc path
        pe.handleLasops(4);
        QCOMPARE(pe.lasopsc(), 4);
        // chaine_[0] is allocated when modepsc_!=0 inside the open path.
        QVERIFY(pe.chaineForTesting()[0] != nullptr);
    }

    // Test 4 — negative-erase, lasopsc<-1. handleLasops(1) opens the file,
    // then handleLasops(-5) clears chaine_[1] (since -lasopsc-4 ==
    // -(-5)-4 == 1). modepsc_ must be 1 so the chaine[] buffers exist.
    void PsEmitter_handleLasops_negative_erase() {
        PsEmitter pe;
        pe.setModepscForTesting(1);
        pe.handleLasops(1);
        QVERIFY(pe.chaineForTesting()[1] != nullptr);
        // Pre-poison chaine_[1] so we can detect the clobber.
        pe.chaineForTesting()[1][0] = 'X';
        QCOMPARE(pe.chaineForTesting()[1][0], 'X');

        pe.handleLasops(-5);
        QCOMPARE(pe.lasopsc(), -5);
        // The first byte of chaine_[1] must be cleared by the negative-erase
        // branch (idx = -lasopsc-4 = -(-5)-4 = 1).
        QCOMPARE(pe.chaineForTesting()[1][0], '\0');
    }

    // Test 5 — modepsc==0 chaine guard. With modepsc_==0 the calloc loop
    // is skipped, so chaine_[0] stays nullptr after handleLasops(1).
    void PsEmitter_handleLasops_modepsc_zero_skips_chaine() {
        PsEmitter pe;
        // modepsc_ default is 0 — do NOT call setModepscForTesting.
        pe.handleLasops(1);
        QCOMPARE(pe.lasopsc(), 1);
        QVERIFY(pe.chaineForTesting()[0] == nullptr);
    }

    // Test 6 — pyFlip helper. Asserts that with ypixels_==600,
    // pyFlip(20) == 580. Nails down the flip-INSIDE-PsEmitter contract
    // before Plan 03 wires it. README_COORDS.md mandate.
    void PsEmitter_pyFlip_yields_ypixels_minus_y() {
        PsEmitter pe;
        pe.setCanvasDims(800, 600);
        QCOMPARE(pe.pyFlip(20), 580);
        QCOMPARE(pe.pyFlip(0), 600);
        QCOMPARE(pe.pyFlip(600), 0);
    }
};

int main(int argc, char* argv[]) {
    qputenv("QT_QPA_PLATFORM", "offscreen");
    XvueApp::ensure();
    TestXvueQtPostscript tc;
    int   qt_argc   = 1;
    char  qt_arg0[] = "xvue_qt_postscript_tests";
    char* qt_argv[] = { qt_arg0, nullptr };
    (void)argc;
    (void)argv;
    return QTest::qExec(&tc, qt_argc, qt_argv);
}

#include "test_xvue_qt_postscript.moc"
