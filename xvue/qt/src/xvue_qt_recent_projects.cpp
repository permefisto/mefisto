// xvue/qt/src/xvue_qt_recent_projects.cpp
// Phase 6.0 Plan 01 (scaffold): no-op refresh + onClearRequested.
#include "xvue_qt_recent_projects.h"

XvueRecentProjectsMenu::XvueRecentProjectsMenu(XvueMenuBridge* bridge,
                                               QWidget* parent)
    : QMenu(parent), bridge_(bridge) {}

XvueRecentProjectsMenu::~XvueRecentProjectsMenu() = default;

void XvueRecentProjectsMenu::refresh() {
    // Plan 04 fills: clear actions, iterate XvuePrefs::recentProjects(),
    // build QAction per path, append "Clear Recent" with confirm dialog.
}

void XvueRecentProjectsMenu::onClearRequested() {
    // Plan 04 fills: confirm + XvuePrefs::clearRecentProjects() + refresh().
}
