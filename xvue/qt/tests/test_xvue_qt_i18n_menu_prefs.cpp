// xvue/qt/tests/test_xvue_qt_i18n_menu_prefs.cpp
// Phase 6.0 Plan 02 (06.0-02): QTest suite covering the three filled-in
// foundation classes (i18n / menu-bridge / prefs).
//
// Coverage:
//   i18n         testBilingualNonEmptyEveryEntry, testAppNameIdentical
//   menu-bridge  testQueueMultiChar, testQueueTrailingCR, testQueueEmpty,
//                testQueueCapEnforced, testQueueContextPopulator,
//                testMarkModuleRegistered
//   prefs        testQSettingsRoundTripGeometry, testQSettingsRoundTripState,
//                testQSettingsRoundTripConsoleVisible,
//                testQSettingsRoundTripColorScheme,
//                testRecentProjectsCap10, testRecentProjectsDedupe,
//                testRecentProjectsClear, testCorruptSettingsDefault,
//                testPerModuleGroupIsolation
//
// Note on language flip: xvueIsEnglish() caches its result via a C++17
// static-local on first call, so the tests in this binary cannot toggle
// FR/EN mid-process. Per the plan, language-flip is verified manually on
// the UI-SPEC sign-off; here we just verify the table is fully populated
// (every MsgId returns a non-empty string) and that AppName is identical
// across both languages.
//
// QSettings isolation: per RESEARCH §Pitfall 12.3 we re-route the QSettings
// UserScope path to QDir::tempPath() in initTestCase so the test does not
// touch the real $HOME/.config/LJLL/mefisto.conf. Each test calls clear()
// in init() to start from a known empty state.
#include "xvue_qt_i18n.h"
#include "xvue_qt_menu_bridge.h"
#include "xvue_qt_prefs.h"

#include <QtTest/QtTest>
#include <QApplication>
#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QMenu>
#include <QSettings>
#include <QStringList>
#include <QTemporaryDir>

class TestXvueQtI18nMenuPrefs : public QObject {
    Q_OBJECT

private slots:
    void initTestCase() {
        // Isolate from the user's real config (RESEARCH Pitfall 12.3).
        // QSettings::setPath only affects the format/scope pair; we use
        // IniFormat/UserScope everywhere in xvue_qt_prefs.cpp.
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope,
                           QDir::tempPath());
        QCoreApplication::setOrganizationName(QStringLiteral("LJLL"));
        QCoreApplication::setApplicationName(QStringLiteral("mefisto"));
    }

    void init() {
        // Per-test reset of the test QSettings file so tests do not leak
        // state into each other.
        QSettings s;
        s.clear();
        s.sync();
    }

    void cleanup() {
        // Drain Qt's deleteLater queue between tests to keep
        // QApplication-shutdown atexit happy (mirrors menu_scaffold).
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }

    void cleanupTestCase() {
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }

    // ---------------------------------------------------------------------
    // i18n tests
    // ---------------------------------------------------------------------
    void testBilingualNonEmptyEveryEntry() {
        for (int i = 0; i < static_cast<int>(MsgId::_Count_); ++i) {
            // Disambiguate from QObject::tr (both now in scope via Q_OBJECT).
            const char* s = ::tr(static_cast<MsgId>(i));
            QVERIFY2(s != nullptr && s[0] != '\0',
                     qPrintable(QStringLiteral("empty string at MsgId=%1").arg(i)));
        }
    }

    void testAppNameIdentical() {
        // AppName is "MEFISTO" in both languages.
        QCOMPARE(QString::fromUtf8(::tr(MsgId::AppName)),
                 QStringLiteral("MEFISTO"));
    }

    // ---------------------------------------------------------------------
    // menu-bridge tests
    // ---------------------------------------------------------------------
    void testQueueMultiChar() {
        XvueMenuBridge mb;
        QCOMPARE(mb.queueSize(), 0);
        mb.queueLexicon(QStringLiteral("5;90;"));
        // 5 chars + trailing CR = 6 elements.
        QCOMPARE(mb.queueSize(), 6);
        const char expected[] = {'5', ';', '9', '0', ';', 13};
        for (char c : expected) {
            auto got = mb.popChar();
            QVERIFY(got.has_value());
            QCOMPARE(static_cast<int>(*got), static_cast<int>(c));
        }
        // 7th pop returns nullopt.
        QVERIFY(!mb.popChar().has_value());
    }

    void testQueueTrailingCR() {
        XvueMenuBridge mb;
        mb.queueLexicon(QStringLiteral("99;"));
        // 3 chars + CR = 4 elements.
        QCOMPARE(mb.queueSize(), 4);
        // Skip '9','9',';'.
        mb.popChar();
        mb.popChar();
        mb.popChar();
        // Last element MUST be 13 (Fortran SACLAV line terminator).
        auto cr = mb.popChar();
        QVERIFY(cr.has_value());
        QCOMPARE(static_cast<int>(*cr), 13);
    }

    void testQueueEmpty() {
        XvueMenuBridge mb;
        QVERIFY(!mb.popChar().has_value());
    }

    void testQueueCapEnforced() {
        // T-06.0-BRIDGE-01: queue cap at kMaxQueueSize=10000. Each "x" push
        // contributes 2 elements ('x' + CR), so 6000 pushes target 12000
        // elements — must cap at 10000.
        XvueMenuBridge mb;
        for (int i = 0; i < 6000; ++i) {
            mb.queueLexicon(QStringLiteral("x"));
        }
        QVERIFY(mb.queueSize() <= XvueMenuBridge::kMaxQueueSize);
        // We expect to have stopped accepting pushes once the next push
        // would breach the cap. With cap=10000, accepting 5000 pushes
        // (10000 elements) is the limit; the 5001st push (would-be 10002)
        // is rejected wholesale. Final size should therefore be exactly
        // 10000.
        QCOMPARE(mb.queueSize(), XvueMenuBridge::kMaxQueueSize);
    }

    void testQueueContextPopulator() {
        XvueMenuBridge mb;
        mb.setContextMenuPopulator([](QMenu* m){ m->addAction(QStringLiteral("fake")); });
        QMenu menu;
        mb.populateContextMenu(&menu);
        QCOMPARE(menu.actions().size(), 1);
    }

    void testMarkModuleRegistered() {
        XvueMenuBridge mb;
        QVERIFY(!mb.hasRegisteredModule());
        mb.markModuleRegistered();
        QVERIFY(mb.hasRegisteredModule());
    }

    // ---------------------------------------------------------------------
    // prefs tests
    // ---------------------------------------------------------------------
    void testQSettingsRoundTripGeometry() {
        XvuePrefs::initialize("mail");
        const QByteArray blob("geom-blob-1");
        XvuePrefs::saveWindowGeometry(blob);
        QCOMPARE(XvuePrefs::windowGeometry(), blob);
    }

    void testQSettingsRoundTripState() {
        XvuePrefs::initialize("mail");
        const QByteArray blob("state-blob-1");
        XvuePrefs::saveWindowState(blob);
        QCOMPARE(XvuePrefs::windowState(), blob);
    }

    void testQSettingsRoundTripConsoleVisible() {
        XvuePrefs::initialize("mail");
        XvuePrefs::saveConsoleVisible(false);
        QVERIFY(!XvuePrefs::consoleVisible());
        XvuePrefs::saveConsoleVisible(true);
        QVERIFY(XvuePrefs::consoleVisible());
    }

    void testQSettingsRoundTripColorScheme() {
        XvuePrefs::initialize("mail");
        XvuePrefs::saveColorScheme(QStringLiteral("dark"));
        QCOMPARE(XvuePrefs::colorScheme(), QStringLiteral("dark"));
        // Invalid value rejected silently — current value unchanged.
        XvuePrefs::saveColorScheme(QStringLiteral("invalid-value"));
        QCOMPARE(XvuePrefs::colorScheme(), QStringLiteral("dark"));
        XvuePrefs::saveColorScheme(QStringLiteral("system"));
        QCOMPARE(XvuePrefs::colorScheme(), QStringLiteral("system"));
    }

    void testRecentProjectsCap10() {
        XvuePrefs::initialize("mail");
        // Create 11 real temp dirs so D-06 prune-on-read does not strip
        // them as missing.
        QTemporaryDir tmp;
        QVERIFY(tmp.isValid());
        QStringList made;
        for (int i = 0; i < 11; ++i) {
            const QString p = tmp.filePath(QStringLiteral("proj%1").arg(i));
            QDir().mkpath(p);
            made << p;
            XvuePrefs::pushRecentProject(p);
        }
        const QStringList list = XvuePrefs::recentProjects();
        QCOMPARE(list.size(), 10);
        // First-pushed (oldest) is evicted.
        QVERIFY(!list.contains(made.first()));
        // Last-pushed (newest) is at the head.
        QCOMPARE(list.first(), made.last());
    }

    void testRecentProjectsDedupe() {
        XvuePrefs::initialize("mail");
        QTemporaryDir tmp;
        QVERIFY(tmp.isValid());
        const QString a = tmp.filePath(QStringLiteral("a"));
        const QString b = tmp.filePath(QStringLiteral("b"));
        QDir().mkpath(a);
        QDir().mkpath(b);
        XvuePrefs::pushRecentProject(a);
        XvuePrefs::pushRecentProject(b);
        XvuePrefs::pushRecentProject(a);   // dedupe + re-promote to head
        const QStringList list = XvuePrefs::recentProjects();
        QCOMPARE(list.size(), 2);
        QCOMPARE(list[0], a);
        QCOMPARE(list[1], b);
    }

    void testRecentProjectsClear() {
        XvuePrefs::initialize("mail");
        QTemporaryDir tmp;
        QVERIFY(tmp.isValid());
        const QString a = tmp.filePath(QStringLiteral("a"));
        QDir().mkpath(a);
        XvuePrefs::pushRecentProject(a);
        XvuePrefs::clearRecentProjects();
        QVERIFY(XvuePrefs::recentProjects().isEmpty());
    }

    void testCorruptSettingsDefault() {
        XvuePrefs::initialize("mail");
        // T-06.0-SETTINGS-01 contract is: "do not crash on type-mismatched
        // or corrupt value". Two paths to verify:
        //
        //  1. Inject a non-ByteArray value under the geometry key.
        //     QSettings + IniFormat persists EVERY value as a string and
        //     lets QVariant::toByteArray coerce on read; the round-trip
        //     is lossy but never throws. The call returning normally is
        //     the actual mitigation.
        //
        //  2. (Stronger) saveColorScheme silently rejects an invalid
        //     value; colorScheme() then clamps any unknown value back to
        //     "system". This exercises the explicit clamp branch.
        //
        // Path 1: inject corrupt value, prove getter returns a value.
        QSettings s;
        s.setValue(QStringLiteral("mail/window/geometry"), 42);
        s.sync();
        const QByteArray got = XvuePrefs::windowGeometry();
        // Do not assert the exact bytes — Qt INI coercion of int → bytes
        // varies (4-byte BE integer vs ASCII "42" vs Qt @Variant blob
        // depending on Qt version + INI format flags). Accept any size,
        // including the "@Variant(...)" Qt-encoded blob; the threat is
        // mitigated as long as the call returns without throwing.
        Q_UNUSED(got);

        // Path 2: corrupt the color-scheme key by writing a garbage
        // string, prove the getter clamps to "system".
        s.setValue(QStringLiteral("ui/color-scheme"), QStringLiteral("?garbage?"));
        s.sync();
        QCOMPARE(XvuePrefs::colorScheme(), QStringLiteral("system"));
    }

    void testPerModuleGroupIsolation() {
        // Per UI-SPEC §Window Chrome: window/* keys are namespaced by
        // module ("mail/window/console-visible" vs "elas/window/console-
        // visible"). saveConsoleVisible(true) under "mail" should NOT
        // affect "elas".
        XvuePrefs::initialize("mail");
        XvuePrefs::saveConsoleVisible(true);
        XvuePrefs::initialize("elas");
        // For elas, no value yet; consoleVisible(false) should return
        // the explicit fallback (false) — proving "mail"'s true did not
        // leak across the per-module key namespace.
        QVERIFY(!XvuePrefs::consoleVisible(false));
        // Re-initialize back to "mail" — the saved true is still there.
        XvuePrefs::initialize("mail");
        QVERIFY(XvuePrefs::consoleVisible(false));
    }
};

int main(int argc, char** argv) {
    // QApplication (not QCoreApplication) so future tests that touch
    // QWidget-derived dialogs can share this binary. Offscreen platform
    // keeps xvfb-run optional but consistent with Phase 5.
    qputenv("QT_QPA_PLATFORM", "offscreen");
    QApplication app(argc, argv);
    TestXvueQtI18nMenuPrefs suite;
    return QTest::qExec(&suite, argc, argv);
}

#include "test_xvue_qt_i18n_menu_prefs.moc"
