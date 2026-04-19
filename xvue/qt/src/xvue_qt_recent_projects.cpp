// xvue/qt/src/xvue_qt_recent_projects.cpp
// Phase 6.0 Plan 04: File → Recent Projects submenu population (D-06).
//
// refresh() rebuilds the menu from XvuePrefs::recentProjects() (which Plan 02
// wires to QSettings + path-existence prune):
//   - One QAction per recent path; label = basename, tooltip = full path
//   - Triggering an action emits openProjectRequested(absPath) so the host
//     window (Plan 06) can fork a new pp*_qt with that project
//   - Empty list shows a disabled "(none)" placeholder so the submenu is
//     never empty (Qt would otherwise hide it on some platforms)
//   - Always-visible "Clear Recent" action at the bottom guards a
//     destructive-confirm QMessageBox per UI-SPEC §Destructive actions
//
// T-06.0-LAMBDA-01: each QAction is a child of `this` (QMenu::addAction) so
// the menu's destruction destroys the action before any queued lambda could
// fire on a dead `this`. Safe.
#include "xvue_qt_recent_projects.h"

#include "xvue_qt_i18n.h"
#include "xvue_qt_menu_bridge.h"
#include "xvue_qt_prefs.h"

#include <QAction>
#include <QFileInfo>
#include <QMessageBox>

XvueRecentProjectsMenu::XvueRecentProjectsMenu(XvueMenuBridge* bridge,
                                               QWidget* parent)
    : QMenu(parent), bridge_(bridge)
{
    setTitle(xvueT(MsgId::FileRecentSubmenu));
    refresh();
}

XvueRecentProjectsMenu::~XvueRecentProjectsMenu() = default;

void XvueRecentProjectsMenu::refresh() {
    clear();   // remove all existing QActions (ownership = this menu)

    const QStringList recent = XvuePrefs::recentProjects();
    if (recent.isEmpty()) {
        // Disabled "(none)" placeholder. Translated label not yet in MsgId
        // table — Qt's tr() returns the literal until Plan 02 adds an entry.
        QAction* empty = addAction(tr("(none)"));
        empty->setEnabled(false);
    } else {
        for (const QString& path : recent) {
            const QString label = QFileInfo(path).fileName();
            QAction* a = addAction(label);
            a->setToolTip(path);
            // Capture path by value so each action carries its own copy.
            connect(a, &QAction::triggered, this, [this, path]{
                emit openProjectRequested(path);
            });
        }
    }

    addSeparator();
    QAction* clearAct = addAction(xvueT(MsgId::FileRecentClear));
    connect(clearAct, &QAction::triggered,
            this,     &XvueRecentProjectsMenu::onClearRequested);
}

void XvueRecentProjectsMenu::onClearRequested() {
    // UI-SPEC §Destructive actions: confirm before clearing.
    // QMessageBox::question with explicit Ok+Cancel matches the destructive-
    // confirm template (default = Cancel = safe).
    const auto ret = QMessageBox::question(
        this->parentWidget(),
        xvueT(MsgId::FileRecentClearConfirmTitle),
        xvueT(MsgId::DestructiveConfirmBodyGeneric),
        QMessageBox::Ok | QMessageBox::Cancel,
        QMessageBox::Cancel);
    if (ret == QMessageBox::Ok) {
        XvuePrefs::clearRecentProjects();
        refresh();
    }
}
