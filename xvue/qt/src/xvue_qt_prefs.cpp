// xvue/qt/src/xvue_qt_prefs.cpp
// Phase 6.0 Plan 02 (06.0-02): QSettings-backed implementation of the 12
// static methods declared in xvue_qt_prefs.h. Replaces the Plan 01 stub.
//
// Storage layout (06.0-UI-SPEC §Window Chrome & Layout — locked):
//   Organization  : "LJLL"        (QCoreApplication::organizationName)
//   Application   : "mefisto"     (QCoreApplication::applicationName)
//   Format        : INI / UserScope ($HOME/.config/LJLL/mefisto.conf)
//
//   Per-module keys (one INI section per pp*_qt executable):
//     <module>/window/geometry       QByteArray (saveGeometry blob)
//     <module>/window/state          QByteArray (saveState blob)
//     <module>/window/console-visible bool
//     <module>/recent-projects       QStringList (LRU, capped at 10)
//
//   Shared keys (single value across all modules — UI-SPEC §Color scheme):
//     ui/color-scheme                "system" | "light" | "dark"
//
// T-06.0-SETTINGS-01 mitigation: every getter is wrapped in try/catch with
// a safe default. QSettings::value().toByteArray()/toBool() coerce silently
// on type mismatch already, but the try/catch is documented defensive depth
// (some Qt builds have been observed to throw on read of corrupt INI).
//
// D-06 (06.0-CONTEXT.md): recentProjects() prunes paths that no longer
// exist on disk. Pruning happens at READ time so the QSettings file does
// not silently grow forever; the on-disk list is rewritten on the next
// pushRecentProject().
#include "xvue_qt_prefs.h"

#include <QCoreApplication>
#include <QFileInfo>
#include <QSettings>

QString XvuePrefs::moduleGroup_;

namespace {

// File-local helpers. Free-function `scoped` takes the group as a parameter
// so it does not need access to XvuePrefs::moduleGroup_ (private static).
QString scoped(const QString& group, const QString& key) {
    if (group.isEmpty()) {
        return key;
    }
    return group + QStringLiteral("/") + key;
}

// All getters/setters share a single QSettings construction pattern. Per
// Qt docs, constructing QSettings is cheap (no heap I/O on UNIX with
// IniFormat — the file is mmap'd lazily); re-instantiating per call keeps
// us free of cross-call cache invalidation issues.
QSettings settings() {
    return QSettings(QSettings::IniFormat, QSettings::UserScope,
                     QStringLiteral("LJLL"), QStringLiteral("mefisto"));
}

}  // namespace

void XvuePrefs::initialize(const char* moduleName) {
    // Idempotent — safe to call from every pp*_qt main(). Sets
    // QCoreApplication::organizationName + applicationName so the rest of
    // Qt (QSettings default ctor, QStandardPaths, etc.) finds the right
    // location. We also seed our own moduleGroup_ for per-module keys.
    QCoreApplication::setOrganizationName(QStringLiteral("LJLL"));
    QCoreApplication::setApplicationName(QStringLiteral("mefisto"));
    moduleGroup_ = QString::fromLatin1(moduleName ? moduleName : "");
}

QByteArray XvuePrefs::windowGeometry() {
    try {
        return settings()
            .value(scoped(moduleGroup_, QStringLiteral("window/geometry")), QByteArray())
            .toByteArray();
    } catch (...) {
        return QByteArray();
    }
}

void XvuePrefs::saveWindowGeometry(const QByteArray& b) {
    settings().setValue(scoped(moduleGroup_, QStringLiteral("window/geometry")), b);
}

QByteArray XvuePrefs::windowState() {
    try {
        return settings()
            .value(scoped(moduleGroup_, QStringLiteral("window/state")), QByteArray())
            .toByteArray();
    } catch (...) {
        return QByteArray();
    }
}

void XvuePrefs::saveWindowState(const QByteArray& b) {
    settings().setValue(scoped(moduleGroup_, QStringLiteral("window/state")), b);
}

bool XvuePrefs::consoleVisible(bool fallback) {
    try {
        return settings()
            .value(scoped(moduleGroup_, QStringLiteral("window/console-visible")), fallback)
            .toBool();
    } catch (...) {
        return fallback;
    }
}

void XvuePrefs::saveConsoleVisible(bool v) {
    settings().setValue(scoped(moduleGroup_, QStringLiteral("window/console-visible")), v);
}

QStringList XvuePrefs::recentProjects() {
    try {
        const QStringList raw = settings()
            .value(scoped(moduleGroup_, QStringLiteral("recent-projects")), QStringList())
            .toStringList();
        // D-06: prune missing paths silently. We do not rewrite the INI
        // here — the next pushRecentProject() call rewrites with the
        // pruned list. Read is hot-path; write is cold.
        QStringList pruned;
        pruned.reserve(raw.size());
        for (const QString& p : raw) {
            if (!p.isEmpty() && QFileInfo(p).exists()) {
                pruned << p;
            }
        }
        return pruned;
    } catch (...) {
        return QStringList();
    }
}

void XvuePrefs::pushRecentProject(const QString& absPath) {
    if (absPath.isEmpty()) {
        return;
    }
    QStringList list = recentProjects();   // already pruned
    list.removeAll(absPath);                // dedupe (case-sensitive — paths are)
    list.prepend(absPath);                  // most-recent at head
    while (list.size() > 10) {              // D-06 cap
        list.removeLast();
    }
    settings().setValue(scoped(moduleGroup_, QStringLiteral("recent-projects")), list);
}

void XvuePrefs::clearRecentProjects() {
    settings().setValue(scoped(moduleGroup_, QStringLiteral("recent-projects")),
                        QStringList());
}

QString XvuePrefs::colorScheme() {
    // ui/color-scheme is SHARED across all modules per UI-SPEC §Color
    // scheme — one global appearance preference, not a per-pp*_qt setting.
    try {
        const QString v = settings()
            .value(QStringLiteral("ui/color-scheme"), QStringLiteral("system"))
            .toString();
        // T-06.0-SETTINGS-01 — clamp unknown/corrupt to "system".
        if (v == QLatin1String("system") ||
            v == QLatin1String("light")  ||
            v == QLatin1String("dark")) {
            return v;
        }
        return QStringLiteral("system");
    } catch (...) {
        return QStringLiteral("system");
    }
}

void XvuePrefs::saveColorScheme(const QString& s) {
    // T-06.0-SETTINGS-01 — only valid values are persisted; unknown values
    // are rejected silently (caller bug, not user-facing error).
    if (s == QLatin1String("system") ||
        s == QLatin1String("light")  ||
        s == QLatin1String("dark")) {
        settings().setValue(QStringLiteral("ui/color-scheme"), s);
    }
}
