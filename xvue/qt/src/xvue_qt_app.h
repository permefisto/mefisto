// xvue/qt/src/xvue_qt_app.h
// Phase 1 (D-01, D-07, D-08, D-09): XvueApp owns the process-lifetime
// QApplication and the (lazily allocated, reopen-able) XvueWindow slot.
// Source: 01-RESEARCH.md §Pattern 1; 01-CONTEXT.md D-01..D-09.
#pragma once
#include <memory>
#include <mutex>

class QApplication;
class XvueWindow;
class XvueMenuBridge;   // Phase 6.0 Plan 06: menuBridge() forwarding accessor

class XvueApp {
public:
    // Idempotent; safe to call from any extern "C" entry point as the first
    // statement (D-18 makes this a precondition for XVUE_QT_ASSERT_MAIN_THREAD).
    static void ensure();

    // Non-null after ensure() has run. Use for qApp-level introspection.
    static QApplication* qapp();

    // Reference to the unique_ptr slot that holds the live XvueWindow (may be
    // null between xvfermer_ and the next xvinitgraphique_). Callers must not
    // cache the returned reference across xvfermer_ calls.
    static std::unique_ptr<XvueWindow>& window_slot();

    static int font_id_;
    static void load_bundled_font_();

    // Phase 5 (D-03, EVENT-08). Re-entrancy counter for nested waitForEvent().
    // Main-thread-only (SHELL-07) — no atomics. Phase 6 modal dialogs query
    // this via blockingDepth() > 0 to refuse QDialog::exec() re-entry.
    static int blockingDepth();

    // Phase 6.0 Plan 06 additions.
    // Convenience accessor — forwards to window_slot()->menuBridge() with
    // nullptr guards (returns nullptr if window_slot is empty or bridge not
    // installed). Mirrors the blockingDepth() static style so call sites
    // inside extern "C" entries stay consistent.
    static XvueMenuBridge* menuBridge();

    // UX-13 (D-05). Reads XvuePrefs::colorScheme() and applies the palette.
    // Call on startup after ensure() and whenever XvuePreferencesDialog
    // commits a new value. Idempotent. The QStyleHints::colorSchemeChanged
    // signal connected inside ensure() invokes this when the user flips
    // desktop dark-mode while pp*_qt is running (Qt 6.5+ feature; on DEs
    // without a theme plugin the signal may never fire — fallback: user
    // restarts app, applyColorSchemePreference on startup syncs).
    static void applyColorSchemePreference();

    // RAII guard (declared in Plan 02's xvue_qt_event.h). Friending keeps the
    // counter strictly internal — only the guard may touch blockingDepth_.
    friend struct BlockingDepthGuard;

private:
    static std::once_flag                     once_flag_;
    static std::unique_ptr<QApplication>      qapp_;
    static std::unique_ptr<XvueWindow>        window_;
    static void teardown_atexit();

    // Phase 5 (D-03). Process-wide; main thread only (SHELL-07). Initialized
    // to 0 in xvue_qt_app.cpp — static zero-init is also a safety net.
    static int blockingDepth_;
};
