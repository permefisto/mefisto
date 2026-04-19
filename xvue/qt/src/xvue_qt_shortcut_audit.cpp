// xvue/qt/src/xvue_qt_shortcut_audit.cpp
// Phase 6.0 Plan 01 (scaffold): logs the QAction child count on demand.
// Plan 02/04 fills the actual collision-detection loop.
#include "xvue_qt_shortcut_audit.h"

#include <QAction>
#include <QDebug>
#include <QWidget>

void XvueShortcutAudit::validateNoCollisions(QWidget* win) {
    if (!win) return;
    auto actions = win->findChildren<QAction*>();
    qInfo() << "XvueShortcutAudit: scanned" << actions.size() << "QAction(s)";
}
