// xvue/qt/src/xvue_qt_prefs.cpp
// Phase 6.0 Plan 01 (scaffold): inert stub bodies. Plan 02 wires QSettings.
#include "xvue_qt_prefs.h"

QString XvuePrefs::moduleGroup_;

void XvuePrefs::initialize(const char* moduleName) {
    // Plan 02 will also call QCoreApplication::setOrganizationName("Mefisto")
    // and ::setApplicationName(moduleName) so QSettings picks the right INI.
    moduleGroup_ = QString::fromUtf8(moduleName ? moduleName : "");
}

QByteArray XvuePrefs::windowGeometry() { return {}; }
void       XvuePrefs::saveWindowGeometry(const QByteArray&) {}

QByteArray XvuePrefs::windowState() { return {}; }
void       XvuePrefs::saveWindowState(const QByteArray&) {}

bool XvuePrefs::consoleVisible(bool fallback) { return fallback; }
void XvuePrefs::saveConsoleVisible(bool) {}

QStringList XvuePrefs::recentProjects() { return {}; }
void        XvuePrefs::pushRecentProject(const QString&) {}
void        XvuePrefs::clearRecentProjects() {}

QString XvuePrefs::colorScheme() { return QStringLiteral("system"); }
void    XvuePrefs::saveColorScheme(const QString&) {}
