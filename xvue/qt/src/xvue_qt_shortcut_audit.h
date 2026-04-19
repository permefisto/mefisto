// xvue/qt/src/xvue_qt_shortcut_audit.h
// Phase 6.0 Plan 01 (scaffold): runtime shortcut-collision check. Each
// per-module registerXxxActions(win) calls validateNoCollisions(win) at the
// end so QKeySequence dupes Q_ASSERT loudly during testing (RESEARCH §10
// layer 3).
#pragma once

class QWidget;

class XvueShortcutAudit {
public:
    // Walks win->findChildren<QAction*>() and Q_ASSERTs no collisions.
    // Plan 02/04 may extend with allowlist for module-scoped reuse.
    static void validateNoCollisions(QWidget* win);
};
