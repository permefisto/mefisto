// xvue/qt/src/xvue_qt_recent_projects.h
// Phase 6.0 Plan 01 (scaffold): "Open Recent" submenu populator. Plan 04
// fills refresh() with QAction-per-path + Clear Recent + D-06 startup prune.
#pragma once
#include <QMenu>

class XvueMenuBridge;

class XvueRecentProjectsMenu : public QMenu {
    Q_OBJECT
public:
    explicit XvueRecentProjectsMenu(XvueMenuBridge* bridge, QWidget* parent = nullptr);
    ~XvueRecentProjectsMenu() override;

    // Rebuild actions from XvuePrefs::recentProjects(). Plan 04 fills.
    void refresh();

signals:
    void openProjectRequested(const QString& absPath);

private slots:
    void onClearRequested();

private:
    XvueMenuBridge* bridge_ = nullptr;
};
