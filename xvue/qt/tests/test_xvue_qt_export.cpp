// xvue/qt/tests/test_xvue_qt_export.cpp
// Phase 7 Plan 04 (EXPORT-02, EXPORT-05): unit tests for XvueExport.
//
// Coverage matrix (from 07-04-PLAN.md `<behavior>` block):
//   1. PNG round-trip — known-marker pixel survives toImage()->save()->load().
//   2. JPEG round-trip + default quality 90 — file written, decodable.
//   3. JPEG quality from QSettings — quality 50 yields smaller file than 95.
//   4. PDF page geometry — 800x600 logical canvas → /MediaBox [0 0 800 600]
//      at 72 dpi (Pitfall 7 from CONTEXT.md D-15: setPageSize uses LOGICAL
//      pixels — XvueCanvas::width()/height() — NOT backing_->width()).
//   5. Write-failure handling — savePngTo("/nonexistent/dir/foo.png") returns
//      false without crashing; interactive=false silences the QMessageBox.
//
// All slots run under QT_QPA_PLATFORM=offscreen via xvfb-run --auto-servernum
// so QApplication does not need a real X display. Each slot uses its own
// QTemporaryDir so output paths never collide across tests.
//
// Window lifecycle: xvinitgraphique_() in initTestCase stands up the
// production XvueWindow + XvueCanvas; xvfermer_() in cleanupTestCase tears
// it down. Tests resize the canvas to the well-known 800x600 logical size
// before exporting so the PDF-geometry assertion has a deterministic target.
#include "xvue_qt_app.h"
#include "xvue_qt_canvas.h"
#include "xvue_qt_export.h"
#include "xvue_qt_state.h"
#include "xvue_qt_window.h"

#include <QtTest/QtTest>
#include <QByteArray>
#include <QCoreApplication>
#include <QDir>
#include <QEventLoop>
#include <QFile>
#include <QFileInfo>
#include <QImage>
#include <QImageReader>
#include <QPainter>
#include <QPixmap>
#include <QSettings>
#include <QTemporaryDir>

extern "C" void xvinitgraphique_(void);
extern "C" void xvfermer_(void);

class TestXvueQtExport : public QObject {
    Q_OBJECT

private:
    QTemporaryDir tmpDir_;

    // Resize canvas to a known logical size and paint a deterministic pattern
    // into the backing pixmap so the round-trip tests have something to read
    // back. 800x600 is the canonical PDF-geometry target (07-04-PLAN.md
    // acceptance criterion: /MediaBox [0 0 800 600]).
    void prepareCanvas800x600() {
        auto& slot = XvueApp::window_slot();
        QVERIFY(slot);
        auto* canvas = slot->canvas();
        QVERIFY(canvas);
        // Forcing logical 800x600 — Qt fires resizeEvent which reallocates
        // the backing pixmap to (800 * dpr) x (600 * dpr).
        canvas->resize(800, 600);
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
        // Paint a known marker into backing_ so PNG/JPEG round-trip can
        // verify the bytes survived. XvueState::painter_ is a LONG-LIVED
        // QPainter on backing_ (Phase 2 D-04, D-05). To paint into the
        // backing for test setup we must reuse it — opening a second
        // QPainter on the same QPixmap is illegal in Qt and surfaces as
        // `QPainter::begin: A paint device can only be painted by one
        // painter at a time` at runtime.
        XvueState* st = slot->state();
        QVERIFY(st->backing_ != nullptr);
        QVERIFY(st->painter_ != nullptr && st->painter_->isActive());
        st->painter_->fillRect(st->backing_->rect(), Qt::black);
        // Marker rectangle at backing-coords (10,10)..(60,60) so PNG/JPEG
        // round-trip detects it in the saved bytes.
        st->painter_->fillRect(QRect(10, 10, 50, 50), QColor(255, 0, 0));
    }

private slots:
    void initTestCase() {
        QVERIFY(tmpDir_.isValid());
        // Isolate QSettings so test-set values (e.g. export/jpeg_quality=50)
        // do not bleed into the developer's ~/.config and so that a fresh
        // run starts clean.
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope,
                           tmpDir_.path());
        QCoreApplication::setOrganizationName(QStringLiteral("MEFISTO"));
        QCoreApplication::setApplicationName(QStringLiteral("xvue-test-7-04"));
        QSettings(QStringLiteral("MEFISTO"), QStringLiteral("xvue")).clear();
        // Stand up the production window so canvas + state + backing exist.
        xvinitgraphique_();
        auto& win = XvueApp::window_slot();
        QVERIFY(win != nullptr);
        QVERIFY(win->canvas() != nullptr);
    }

    void cleanupTestCase() {
        xvfermer_();
    }

    void cleanup() {
        // Drain Qt deleteLater between tests so Qt children do not race the
        // leaked-QApplication atexit (D-08 in xvue_qt_app.cpp).
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }

    // ---- Test 1: PNG round-trip ----------------------------------------
    // Paint a known marker pixel into backing_, savePngTo tmpfile, reload
    // via QImage::load, assert the marker pixel matches within 1 bit/channel
    // drift. PNG is lossless so exact match is the right contract.
    void XvueExport_png_roundTrip() {
        prepareCanvas800x600();
        const QString path = tmpDir_.path() + "/png_roundtrip.png";

        QVERIFY(XvueExport::savePngTo(path, /*interactive=*/false));
        QVERIFY(QFile::exists(path));
        QVERIFY(QFileInfo(path).size() > 0);

        QImage reloaded;
        QVERIFY(reloaded.load(path, "PNG"));
        QVERIFY(!reloaded.isNull());
        // Marker pixel — devicePixelRatio scales the backing, so the marker
        // we painted at backing-coords (10,10) lands at pixel (10,10) in the
        // saved PNG (saved bytes are physical pixels).
        const QColor px = reloaded.pixelColor(20, 20);
        QCOMPARE(px.red(),   255);
        QCOMPARE(px.green(),   0);
        QCOMPARE(px.blue(),    0);
    }

    // ---- Test 2: JPEG default quality 90 -------------------------------
    // Default quality is 90 (CONTEXT.md `<discretion>`). Verify the file is
    // written, has nonzero size, and QImageReader can read it back.
    void XvueExport_jpeg_default_quality_90() {
        prepareCanvas800x600();
        // Ensure no quality override is set so the default branch fires.
        QSettings(QStringLiteral("MEFISTO"), QStringLiteral("xvue"))
            .remove(QStringLiteral("export/jpeg_quality"));
        const QString path = tmpDir_.path() + "/jpeg_default.jpg";

        QVERIFY(XvueExport::saveJpegTo(path, /*interactive=*/false));
        QVERIFY(QFile::exists(path));
        QVERIFY(QFileInfo(path).size() > 0);

        QImageReader r(path);
        QVERIFY(r.canRead());
        QImage img = r.read();
        QVERIFY(!img.isNull());
    }

    // ---- Test 3: JPEG quality from QSettings ---------------------------
    // Save the same backing once at quality 50 and once at quality 95.
    // The quality-50 file must be smaller than the quality-95 file.
    // Sanity-checks that QSettings("export/jpeg_quality") drives
    // QImageWriter::setQuality.
    void XvueExport_jpeg_quality_from_QSettings() {
        prepareCanvas800x600();
        const QString lo = tmpDir_.path() + "/jpeg_q50.jpg";
        const QString hi = tmpDir_.path() + "/jpeg_q95.jpg";

        QSettings s(QStringLiteral("MEFISTO"), QStringLiteral("xvue"));

        s.setValue(QStringLiteral("export/jpeg_quality"), 50);
        s.sync();
        QVERIFY(XvueExport::saveJpegTo(lo, /*interactive=*/false));

        s.setValue(QStringLiteral("export/jpeg_quality"), 95);
        s.sync();
        QVERIFY(XvueExport::saveJpegTo(hi, /*interactive=*/false));

        const qint64 loSize = QFileInfo(lo).size();
        const qint64 hiSize = QFileInfo(hi).size();
        QVERIFY2(loSize > 0 && hiSize > 0,
                 qPrintable(QStringLiteral("loSize=%1 hiSize=%2")
                                .arg(loSize).arg(hiSize)));
        QVERIFY2(loSize < hiSize,
                 qPrintable(QStringLiteral("expected q50 < q95: loSize=%1 hiSize=%2")
                                .arg(loSize).arg(hiSize)));

        // Reset for next test.
        s.remove(QStringLiteral("export/jpeg_quality"));
        s.sync();
    }

    // ---- Test 4: PDF page geometry (Pitfall 7) -------------------------
    // 800x600 logical canvas at 72 dpi → MediaBox encodes 800 × 600 in PDF
    // points. Reads the produced .pdf bytes and matches one of the spacing/
    // numeric-format variants Qt's QPdfWriter is observed to emit:
    //   /MediaBox [0 0 800 600]              (integer literals)
    //   /MediaBox [ 0 0 800 600 ]            (integer literals, padded)
    //   /MediaBox [0 0 800.000000 600.000000] (Qt 6.10 float-format actual)
    // The last form is what the running Qt 6.10.2 build emits at the
    // moment; we accept any of the three so the test survives Qt format
    // tweaks across versions but still locks the page-geometry contract.
    void XvueExport_pdf_geometry_72dpi_logical() {
        prepareCanvas800x600();
        const QString path = tmpDir_.path() + "/geometry_800x600.pdf";

        QVERIFY(XvueExport::savePdfTo(path, /*interactive=*/false));
        QVERIFY(QFile::exists(path));
        QVERIFY(QFileInfo(path).size() > 0);

        QFile f(path);
        QVERIFY(f.open(QIODevice::ReadOnly));
        const QByteArray bytes = f.readAll();
        f.close();

        // Accept either integer or float numeric format (Qt 6.10 emits floats).
        const bool tightInt    = bytes.contains("/MediaBox [0 0 800 600]");
        const bool spacedInt   = bytes.contains("/MediaBox [ 0 0 800 600 ]");
        const bool tightFloat  = bytes.contains("/MediaBox [0 0 800.000000 600.000000]");
        const bool spacedFloat = bytes.contains("/MediaBox [ 0 0 800.000000 600.000000 ]");
        const bool found = tightInt || spacedInt || tightFloat || spacedFloat;
        if (!found) {
            // Surface a leading byte-encoded dump so the diagnostic is useful.
            const int idx = bytes.indexOf("/MediaBox");
            QFAIL(qPrintable(QStringLiteral(
                "expected /MediaBox encoding 800 × 600 (any int/float form); "
                "idx=%1; first 64 bytes after MediaBox: %2")
                .arg(idx)
                .arg(QString::fromLatin1(
                     bytes.mid(idx >= 0 ? idx : 0, 64).toPercentEncoding()))));
        }
        QVERIFY(found);
    }

    // ---- Test 5: write-failure returns false ---------------------------
    // savePngTo to a directory that does not exist must return false without
    // crashing. interactive=false silences the QMessageBox so the test
    // doesn't deadlock waiting for user input.
    void XvueExport_writeFailure_returns_false() {
        prepareCanvas800x600();
        const QString badPath = "/nonexistent_dir_for_xvue_test/foo.png";

        QVERIFY(!XvueExport::savePngTo(badPath, /*interactive=*/false));
        QVERIFY(!QFile::exists(badPath));
    }
};

int main(int argc, char* argv[]) {
    qputenv("QT_QPA_PLATFORM", "offscreen");
    // CRITICAL (mirrors test_xvue_qt_window_chrome.cpp): disable
    // XvueConsoleDock stdout-redirect so QtTest's PASS/FAIL output reaches
    // the actual stdout fd.
    qputenv("XVUE_QT_DISABLE_STDOUT_REDIRECT", "1");
    XvueApp::ensure();
    TestXvueQtExport tc;
    int   qt_argc   = 1;
    char  qt_arg0[] = "xvue_qt_export_tests";
    char* qt_argv[] = { qt_arg0, nullptr };
    (void)argc;
    (void)argv;
    return QTest::qExec(&tc, qt_argc, qt_argv);
}

#include "test_xvue_qt_export.moc"
