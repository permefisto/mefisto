// xvue/qt/src/xvue_qt_prefs.h
// Phase 6.0 Plan 01 (scaffold): XvuePrefs static API per RESEARCH §4. Plan 02
// fills the QSettings-backed bodies (org=Mefisto, app=$module).
#pragma once
#include <QByteArray>
#include <QString>
#include <QStringList>

class XvuePrefs {
public:
    // Called once per pp*_qt process by Plan 06's xvue_module_init_ (or by
    // tests directly). Plan 02 sets the QCoreApplication org/app names so
    // QSettings finds the right INI path.
    static void initialize(const char* moduleName);

    static QByteArray windowGeometry();
    static void       saveWindowGeometry(const QByteArray&);

    static QByteArray windowState();
    static void       saveWindowState(const QByteArray&);

    static bool consoleVisible(bool fallback = true);
    static void saveConsoleVisible(bool);

    static QStringList recentProjects();
    static void        pushRecentProject(const QString& absPath);
    static void        clearRecentProjects();

    // "system" | "light" | "dark"  per UI-SPEC §Color scheme.
    static QString colorScheme();
    static void    saveColorScheme(const QString&);

private:
    static QString moduleGroup_;
};
