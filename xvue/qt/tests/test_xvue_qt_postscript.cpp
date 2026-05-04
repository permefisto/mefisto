// xvue/qt/tests/test_xvue_qt_postscript.cpp
// Phase 7 Plan 02 (EXPORT-04): unit tests for PsEmitter::handleLasops state
// machine + the pyFlip helper.
// Phase 7 Plan 03: extends with 8 per-primitive byte-output tests covering
// line / face / epaisseur / typetrait / chargefonte (Courier + fallback) /
// Y-flip-not-double-applied / inactive-no-emit.
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
#include "xvue_qt_api.h"          // MefistoPoint (Plan 03 face() test)
#include "xvue_qt_postscript.h"

#include <QtTest/QtTest>
#include <QCoreApplication>
#include <QDir>
#include <QEventLoop>
#include <QFile>
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
    //
    // VERBATIM semantics from xvuelc.c:1286 — `lasopsc = lasopsc - 100`
    // subtracts 100 from the EXISTING lasopsc_ (the file-static), NOT
    // from the incoming lasops parameter. The legacy callers (effacer at
    // xvuelc.c:1414, 1435) update the file-static lasopsc to (100 + old)
    // BEFORE the recursive call, so 100+1 - 100 == 1 yields the expected
    // "post-erase" state. From a test that calls handleLasops(101)
    // directly without that caller-side mutation, lasopsc_ goes
    // 1 - 100 == -99, which is the verbatim-correct outcome.
    //
    // What the test verifies:
    //  - mode-100 branch was entered (lasopsc_ ended at -99, not 101)
    //  - fpo_ was closed and reopened (TEMPORAIRE.EPS still exists)
    //  - 100-subtract executed (sentinel: -99 == 1 - 100)
    void PsEmitter_handleLasops_mode100_reset() {
        PsEmitter pe;
        pe.handleLasops(1);
        QCOMPARE(pe.lasopsc(), 1);
        FILE* before = pe.fpoForTesting();
        QVERIFY(before != nullptr);

        pe.handleLasops(101);
        // Verbatim: lasopsc_ (was 1) - 100 == -99. NOT 101-100=1.
        QCOMPARE(pe.lasopsc(), -99);
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

    // ===== Plan 03 per-primitive byte-output tests =====
    //
    // Helper: read TEMPORAIRE.EPS contents AFTER the file pointer is closed
    // (handleLasops(0) flushes + closes). Plan 03 emit helpers use buffered
    // fprintf so callers MUST close before reading.
    //
    // Note: handleLasops(0) flushes the pending `concat_` buffer to fpo_
    // before closing — so for line() where the segment-merging logic stages
    // bytes in concat_ rather than writing them immediately, the close is
    // what makes them visible on disk. This mirrors xvuelc.c's :1248 flush.

    // Test 7 — line emit (S opcode). Asserts:
    //   * 6-wide right-justified ints with leading spaces ("    10    580")
    //   * Y-flipped coords (ypixels=600, y=20 -> 580; y=40 -> 560)
    //   * trailing " S\n" segment opcode
    // Source-of-truth format string: xvuelc.c:1954
    //   "%6i %6i %6i %6i %3i %4.2f %4.2f %4.2f %4.2f S\n"
    void PsEmitter_line_emits_S_with_yflip() {
        PsEmitter pe;
        pe.setCanvasDims(800, 600);
        pe.setCourgbForTesting(1.00f, 0.50f, 0.25f);
        pe.setCounbForTesting(0.75f);
        pe.handleLasops(1);
        pe.line(10, 20, 30, 40);
        pe.handleLasops(0);

        QFile f("TEMPORAIRE.EPS");
        QVERIFY(f.open(QIODevice::ReadOnly));
        const QByteArray contents = f.readAll();
        QVERIFY2(contents.contains("    10    580    30    560"),
                 QByteArray("expected 6-wide ints with Y-flipped 20->580 / 40->560: ")
                 + contents.left(200));
        QVERIFY2(contents.contains(" S\n"),
                 QByteArray("expected trailing ' S\\n' segment opcode: ")
                 + contents.left(200));
        QVERIFY2(contents.contains("1.00 0.50 0.25 0.75 S\n"),
                 QByteArray("expected courgb+counb verbatim format: ")
                 + contents.left(200));
    }

    // Test 7b — line emit, counb_==-1 branch (xvuelc.c:1958 "0.00 S\n").
    void PsEmitter_line_emits_S_counb_default() {
        PsEmitter pe;
        pe.setCanvasDims(800, 600);
        pe.setCourgbForTesting(1.00f, 0.50f, 0.25f);
        pe.setCounbForTesting(-1.0f);
        pe.handleLasops(1);
        pe.line(10, 20, 30, 40);
        pe.handleLasops(0);

        QFile f("TEMPORAIRE.EPS");
        QVERIFY(f.open(QIODevice::ReadOnly));
        const QByteArray contents = f.readAll();
        QVERIFY2(contents.contains("1.00 0.50 0.25 0.00 S\n"),
                 QByteArray("expected counb_==-1 branch '0.00 S\\n': ")
                 + contents.left(200));
    }

    // Test 8 — epaisseur emits "%2i epais\n" verbatim from xvuelc.c:1895.
    void PsEmitter_epaisseur_emits_2i_epais() {
        PsEmitter pe;
        pe.setCanvasDims(800, 600);
        pe.handleLasops(1);
        pe.epaisseur(3);
        pe.handleLasops(0);

        QFile f("TEMPORAIRE.EPS");
        QVERIFY(f.open(QIODevice::ReadOnly));
        const QByteArray contents = f.readAll();
        QVERIFY2(contents.contains(" 3 epais\n"),
                 QByteArray("expected ' 3 epais\\n' (2-wide int): ")
                 + contents.left(200));
    }

    // Test 9 — typetrait emits "%2i typet\n" verbatim from xvuelc.c:1856.
    // Plan body claimed "typtr"; xvuelc.c source emits "typet".
    void PsEmitter_typetrait_emits_2i_typet() {
        PsEmitter pe;
        pe.setCanvasDims(800, 600);
        pe.handleLasops(1);
        pe.typetrait(2);
        pe.handleLasops(0);

        QFile f("TEMPORAIRE.EPS");
        QVERIFY(f.open(QIODevice::ReadOnly));
        const QByteArray contents = f.readAll();
        QVERIFY2(contents.contains(" 2 typet\n"),
                 QByteArray("expected ' 2 typet\\n' (xvuelc.c:1856 verbatim): ")
                 + contents.left(200));
    }

    // Test 10 — face emit (F close-and-fill). xvuelc.c:1799 emits the
    // polygon as "x y " coordinate pairs followed by
    // "<n> <r> <g> <b> <a> F\n" close-and-fill. With ypixels=100,
    // y=20->80, y=40->60, y=60->40 (Y-flipped).
    void PsEmitter_face_emit_F_close() {
        PsEmitter pe;
        pe.setCanvasDims(200, 100);
        pe.setCourgbForTesting(1.00f, 0.50f, 0.25f);
        pe.setCounbForTesting(-1.0f);
        pe.handleLasops(1);
        MefistoPoint pts[3] = {{10, 20}, {30, 40}, {50, 60}};
        pe.face(pts, 3);
        pe.handleLasops(0);

        QFile f("TEMPORAIRE.EPS");
        QVERIFY(f.open(QIODevice::ReadOnly));
        const QByteArray contents = f.readAll();
        QVERIFY2(contents.contains("    10     80"),
                 QByteArray("expected first vertex (10, 80) Y-flipped: ")
                 + contents.left(200));
        QVERIFY2(contents.contains(" F\n"),
                 QByteArray("expected trailing ' F\\n' close-and-fill opcode: ")
                 + contents.left(200));
    }

    // Test 11 — chargefonte D-08 mapping: "Courier" -> "/Courier".
    // Output uses xvuelc.c:1553 "%d %d %s charge\n" format (NOT setfont).
    void PsEmitter_chargefonte_courier_maps_to_PSCourier() {
        PsEmitter pe;
        pe.setCanvasDims(800, 600);
        pe.handleLasops(1);
        pe.chargefonte("Courier", 12, 10, 3, false, false);
        pe.handleLasops(0);

        QFile f("TEMPORAIRE.EPS");
        QVERIFY(f.open(QIODevice::ReadOnly));
        const QByteArray contents = f.readAll();
        QVERIFY2(contents.contains("/Courier"),
                 QByteArray("expected '/Courier' PS family name: ")
                 + contents.left(200));
        QVERIFY2(contents.contains("charge\n"),
                 QByteArray("expected 'charge\\n' verb (xvuelc.c:1553 verbatim): ")
                 + contents.left(200));
    }

    // Test 12 — chargefonte D-08 fallback: unknown family -> "/Helvetica".
    // Per CONTEXT.md D-08, DejaVuSansMono (the bundled Phase 3 font) and any
    // other Qt family not in the 4-entry mapping table fall through to
    // /Helvetica with a warn-once stderr line.
    void PsEmitter_chargefonte_unknown_falls_back_to_Helvetica() {
        PsEmitter pe;
        pe.setCanvasDims(800, 600);
        pe.handleLasops(1);
        pe.chargefonte("DejaVu Sans Mono", 12, 10, 3, false, false);
        pe.handleLasops(0);

        QFile f("TEMPORAIRE.EPS");
        QVERIFY(f.open(QIODevice::ReadOnly));
        const QByteArray contents = f.readAll();
        QVERIFY2(contents.contains("/Helvetica"),
                 QByteArray("expected '/Helvetica' fallback per D-08: ")
                 + contents.left(200));
        QVERIFY2(contents.contains("charge\n"),
                 QByteArray("expected 'charge\\n' verb on fallback path too: ")
                 + contents.left(200));
    }

    // Test 13 — Y-flip applied exactly once (Pitfall 2 / README_COORDS.md).
    // pe.line(10, 20, 30, 40) with ypixels=600 must emit Y values 580/560,
    // not 420/440 (would happen if the caller AND helper both flipped).
    void PsEmitter_yflip_not_double_applied() {
        PsEmitter pe;
        pe.setCanvasDims(800, 600);
        pe.setCourgbForTesting(0.0f, 0.0f, 0.0f);
        pe.setCounbForTesting(-1.0f);
        pe.handleLasops(1);
        pe.line(10, 20, 30, 40);
        pe.handleLasops(0);

        QFile f("TEMPORAIRE.EPS");
        QVERIFY(f.open(QIODevice::ReadOnly));
        const QByteArray contents = f.readAll();
        QVERIFY2(contents.contains("   580"),
                 QByteArray("expected single Y-flip 20->580: ")
                 + contents.left(200));
        QVERIFY2(contents.contains("   560"),
                 QByteArray("expected single Y-flip 40->560: ")
                 + contents.left(200));
        QVERIFY2(!contents.contains("   420"),
                 QByteArray("FAIL: Y-flip applied twice (saw 420 = 600-20-160?): ")
                 + contents.left(200));
    }

    // Test 14 — inactive-no-emit guard. With handleLasops never called
    // (lasopsc_==0), every emit helper must early-return without touching
    // fpo_ or concat_. Mirrors xvuelc.c's `if (lasopsc > 0)` gate.
    void PsEmitter_inactive_no_emit() {
        PsEmitter pe;
        pe.setCanvasDims(800, 600);
        QVERIFY(!pe.active());

        // No handleLasops(1) call — fpo_ stays nullptr, lasopsc_ stays 0.
        pe.line(0, 0, 1, 1);
        pe.epaisseur(3);
        pe.typetrait(1);
        MefistoPoint pts[2] = {{0, 0}, {1, 1}};
        pe.face(pts, 2);
        pe.chargefonte("Courier", 12, 10, 3, false, false);

        // None of those should have opened a file.
        QVERIFY(pe.fpoForTesting() == nullptr);
        // concat_ should still be zero-initialized.
        QCOMPARE(pe.concatForTesting()[0], '\0');
        // No TEMPORAIRE.EPS on disk either.
        QVERIFY(!QFile::exists("TEMPORAIRE.EPS"));
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
