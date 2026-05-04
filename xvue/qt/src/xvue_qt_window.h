// xvue/qt/src/xvue_qt_window.h
// Phase 1 (D-02, D-04, D-15, Pitfall 6, Pitfall 8): bare QMainWindow.
// Phase 5 (D-02, EVENT-01): also owns the XvueEventBridge installed as
// event filter on the canvas; lifetime matches the window.
// Phase 6.0 Plan 03 (UX-03): adds menuBridge() accessor + setMenuBridgeForTesting()
// helper. Plan 06 wires the production bridge in the ctor.
// Phase 6.0 Plan 06 (UX-01, UX-04, UX-06, UX-07, UX-09, UX-13): grows from
// the bare QMainWindow into the full Phase 6.0 shell — menu bar (File/View/
// Help with 13 shared QActions), toolbar, status bar with a permanent coord
// label, console dock at the bottom, menu bridge installed in the ctor,
// and QSettings-backed geometry/state/console-visible round-trip via
// XvuePrefs in closeEvent + ctor restoration. Modal-guard helper
// refuseIfBlocking() shows a 3000ms status-bar message when blockingDepth>0.
#pragma once
#include <QMainWindow>
#include "xvue_qt_state.h"

class XvueCanvas;
class XvueEventBridge;
class XvueMenuBridge;          // Phase 6.0 Plan 03 fwd decl
class XvueConsoleDock;         // Phase 6.0 Plan 06 fwd decl
class XvueRecentProjectsMenu;  // Phase 6.0 Plan 06 fwd decl
class QLabel;
class QToolBar;
class QAction;
class QMenu;
class QPoint;

class XvueWindow : public QMainWindow {
    Q_OBJECT  // Pitfall 8: required even with no signals/slots in Phase 1.
public:
    explicit XvueWindow(QWidget* parent = nullptr);
    ~XvueWindow() override;

    XvueState*       state()  { return &state_; }
    XvueCanvas*      canvas() { return canvas_; }
    XvueEventBridge* bridge() { return bridge_; }  // non-null after construction

    // Phase 6.0 Plan 03 (UX-03): menu bridge accessor. Plan 06 wires the real
    // XvueMenuBridge in the ctor so this returns non-null after construction.
    // The pre-drain block in xvue_qt_event.cpp::waitForEvent and
    // xvue_qt_api.cpp::xvsouris_ both null-check, so a future regression that
    // leaves menuBridge_ unset will degrade gracefully rather than crash.
    XvueMenuBridge*  menuBridge() { return menuBridge_; }

    // Phase 6.0 Plan 03 (test-only): inject a bridge BEFORE Plan 06 wires
    // the production one. Ownership: caller retains; we do NOT reparent.
    // Plan 06 ctor replaces this wiring by parenting a fresh bridge to the
    // window. Production code never invokes this setter — it exists solely
    // so Plan 03's pre-drain QTest cases (testPreDrainSingleChar et al.)
    // can exercise the new code path.
    void setMenuBridgeForTesting(XvueMenuBridge* mb) { menuBridge_ = mb; }

    // Phase 6.0 Plan 06: discoverable accessors for tests + 6.1..6.5 modules.
    XvueConsoleDock* consoleDock() { return consoleDock_; }
    QToolBar*        toolBar()     { return toolBar_; }

public slots:
    // Phase 6.0 Plan 06: connected to XvueCanvas::mouseCoords so the
    // permanent status-bar label tracks the cursor in canvas coordinates.
    void updateStatusCoords(QPoint p);

protected:
    // Phase 6.0 Plan 06: persists geometry/state/console-visible to XvuePrefs.
    void closeEvent(QCloseEvent* event) override;

private:
    // Phase 6.0 Plan 06: chrome-construction helpers, called from ctor in
    // the order menu -> toolbar -> status bar -> console dock.
    void buildMenuBar();
    void buildToolBar();
    void buildStatusBar();
    void buildConsoleDock();

    // Phase 6.0 Plan 06: shared QAction handlers (UI-SPEC §6.0 menu contents).
    void onFileOpen();
    void onFileSave();
    void onFileSaveAs();
    // Phase 7 Plan 04 (EXPORT-02, EXPORT-05): File → Export → submenu slots.
    // Each delegates to XvueExport::onMenuExport{Png,Jpeg,Pdf}, which prompts
    // via QFileDialog (QSettings-remembered last_dir) and then writes the
    // canvas backing pixmap to the user-selected path. Plan 05 will add
    // onFileExportGif and onFileCaptureAnimation alongside these.
    void onFileExportPng();
    void onFileExportJpeg();
    void onFileExportPdf();
    void onFileQuit();
    void onViewPreferences();
    void onViewFit();
    void onHelpDocumentation();
    void onHelpAbout();
    void onRecentProjectRequested(const QString& path);

    // Phase 6.0 Plan 06 (D-08): modal-guard helper. Returns true when
    // XvueApp::blockingDepth() > 0 and shows the 3000ms refuse message.
    // Caller uses the return value to early-exit the QAction handler.
    bool refuseIfBlocking();

    XvueState        state_{};    // D-04: single-field state struct
    XvueCanvas*      canvas_     = nullptr;  // Qt-owned via setCentralWidget
    XvueEventBridge* bridge_     = nullptr;  // Phase 5 (D-02): Qt-owned
    XvueMenuBridge*  menuBridge_ = nullptr;  // Phase 6.0 Plan 06: Qt-owned
    XvueConsoleDock* consoleDock_ = nullptr; // Phase 6.0 Plan 06: Qt-owned
    XvueRecentProjectsMenu* recentMenu_ = nullptr;  // child of File menu

    // Status-bar permanent widget (UI-SPEC §Status bar coord-label slot).
    QLabel*          coordLabel_  = nullptr;
    QToolBar*        toolBar_     = nullptr;

    // Phase 6.0 Plan 06: 13 shared QActions (UI-SPEC §6.0 menu contents).
    // Kept as members so tests can discover them and 6.1..6.5 module-action
    // registrations can re-target their visibility/enabled state.
    QAction* actOpen_          = nullptr;
    QAction* actSave_          = nullptr;
    QAction* actSaveAs_        = nullptr;
    // Phase 7 Plan 04 (EXPORT-02, EXPORT-05): replaces the 6.0 placeholder
    // QAction* actExport_ with a real submenu + 3 child QActions. Plan 05
    // will append GIF + Capture-Animation entries to exportMenu_.
    QMenu*   exportMenu_       = nullptr;
    QAction* actExportPng_     = nullptr;
    QAction* actExportJpeg_    = nullptr;
    QAction* actExportPdf_     = nullptr;
    QAction* actQuit_          = nullptr;
    QAction* actToolbarToggle_ = nullptr;
    QAction* actConsoleToggle_ = nullptr;
    QAction* actZoomIn_        = nullptr;
    QAction* actZoomOut_       = nullptr;
    QAction* actFit_           = nullptr;
    QAction* actPreferences_   = nullptr;
    QAction* actDocumentation_ = nullptr;
    QAction* actAbout_         = nullptr;

    // UI-SPEC Flag #5: warn once if $MEFISTOX is unset on first File→Open.
    bool             mefistoxWarned_ = false;
};
