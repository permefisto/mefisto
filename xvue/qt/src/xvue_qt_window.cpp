// xvue/qt/src/xvue_qt_window.cpp
// Phase 1 (D-02): central widget = XvueCanvas; minimal QMainWindow.
// Phase 5 (D-02, EVENT-01): also constructs the XvueEventBridge.
// Phase 6.0 Plan 06 (UX-01, UX-04, UX-06, UX-07, UX-09, UX-13):
//   - menu bar with {File, View, Help} containing 13 shared QActions
//   - toolbar with Open/Save/separator/ZoomIn/ZoomOut/Fit/separator/Console-toggle
//   - status bar with permanent coord label
//   - console dock at the bottom (XvueConsoleDock)
//   - menu bridge installed and visible via XvueWindow::menuBridge()
//   - geometry/state/console-visible round-trip via XvuePrefs
//   - color-scheme palette applied at construction (UX-13)
//   - modal-guard helper refuseIfBlocking() for D-08 modal re-entry
//   - QFileDialog defaults to $MEFISTOX (UI-SPEC Flag #5 warn-once)
//   - Recent Projects submenu (UX-06, D-06)
//   - File→Quit destructive-confirm (UI-SPEC §Destructive actions)
#include "xvue_qt_window.h"
#include "xvue_qt_canvas.h"
#include "xvue_qt_event.h"
#include "xvue_qt_menu_bridge.h"
#include "xvue_qt_console_dock.h"
#include "xvue_qt_recent_projects.h"
#include "xvue_qt_about_dialog.h"
#include "xvue_qt_preferences.h"
#include "xvue_qt_prefs.h"
#include "xvue_qt_i18n.h"
#include "xvue_qt_app.h"
#include "xvue_qt_export.h"   // Phase 7 Plan 04: File → Export submenu wiring

#include <QAction>
#include <QApplication>
#include <QCloseEvent>
#include <QCoreApplication>
#include <QDesktopServices>
#include <QDialog>
#include <QFileDialog>
#include <QKeySequence>
#include <QLabel>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QStatusBar>
#include <QStyle>
#include <QToolBar>
#include <QUrl>

XvueWindow::XvueWindow(QWidget* parent)
    : QMainWindow(parent)
{
    setWindowTitle(QStringLiteral("MEFISTO"));
    resize(1024, 768);                            // UI-SPEC default first-launch

    canvas_ = new XvueCanvas(&state_, this);     // D-15: raw pointer into state_
    setCentralWidget(canvas_);

    // Phase 5 (D-02, EVENT-01). Bridge owns the event filter; lifetime
    // matches the window. Parent=this -> Qt parent-child destruction deletes
    // bridge_ in ~XvueWindow.
    bridge_ = new XvueEventBridge(canvas_, this);
    canvas_->installEventFilter(bridge_);

    // Phase 6.0 Plan 06 (UX-03): install the production menu bridge. Parent=this
    // so Qt destruction deletes it in ~XvueWindow. The Plan 03 test-only setter
    // (setMenuBridgeForTesting) remains usable for unit tests that need a
    // fresh bridge — they should construct one with parent=nullptr to avoid
    // double-destruction, then setMenuBridgeForTesting(theirBridge).
    menuBridge_ = new XvueMenuBridge(this);

    // Build chrome in order: menu -> toolbar -> status bar -> console dock.
    // Order matters because buildToolBar() reuses QActions created in
    // buildMenuBar() (Qt aliasing) and buildConsoleDock() connects the
    // View → Console toggle QAction created in buildMenuBar().
    buildMenuBar();
    buildToolBar();
    buildStatusBar();
    buildConsoleDock();

    // Connect canvas mouseCoords signal (Plan 05) to the status-bar coord
    // label. Plan 03's eventFilter MouseMove emits mouseCoords during the
    // event-bridge dispatch, so the coord label updates whether the user
    // is idle or mid-picking.
    connect(canvas_, &XvueCanvas::mouseCoords,
            this,    &XvueWindow::updateStatusCoords);

    // Restore persisted geometry + state (UX-07). On first launch (no saved
    // bytes) the resize(1024, 768) and default state above stay in effect.
    const QByteArray geom = XvuePrefs::windowGeometry();
    if (!geom.isEmpty()) restoreGeometry(geom);
    const QByteArray st = XvuePrefs::windowState();
    if (!st.isEmpty()) restoreState(st);

    // Apply color-scheme preference (UX-13). XvueApp::ensure() must have run
    // already (Phase 1 D-09 -- XvueWindow is only created via the extern "C"
    // entry points, all of which call XvueApp::ensure() first), so the
    // QApplication instance is guaranteed alive.
    XvueApp::applyColorSchemePreference();
}

XvueWindow::~XvueWindow() = default;

// ---- closeEvent: persist geometry/state/console-visible (UX-07) ----
// Phase 6.1 D-09: when Fortran is blocked on xvsouris_ (blockingDepth()
// > 0), intercept the window close by pushing "99;" through the menu
// bridge so Fortran's LIMTCL save path (prpr/ppmail.f label 9900 ->
// ARRET(0) -> STOP) drives the real exit. If blockingDepth() == 0 we
// fall through to QMainWindow::closeEvent (batch mode, atexit teardown,
// or Fortran already exited — RESEARCH Pitfall 3 "stuck window"
// mitigation).
void XvueWindow::closeEvent(QCloseEvent* event) {
    XvuePrefs::saveWindowGeometry(saveGeometry());
    XvuePrefs::saveWindowState(saveState());
    if (consoleDock_) {
        XvuePrefs::saveConsoleVisible(consoleDock_->isVisible());
    }

    if (menuBridge_ && XvueApp::blockingDepth() > 0) {
        if (consoleDock_) {
            consoleDock_->appendLine(QStringLiteral("[menu] 99;"));
        }
        menuBridge_->queueLexicon(QStringLiteral("99;"));
        event->ignore();
        return;
    }

    QMainWindow::closeEvent(event);
}

// ---- Modal-guard helper (D-08) ----
// Returns true when a blocking read is active (XvueApp::blockingDepth() > 0)
// and shows the 3000ms refuse message; the caller uses the return to early-
// exit the QAction handler. The QAction itself is NOT disabled — the guard
// lives in the handler so the menu item remains clickable; the message
// explains why nothing opened.
bool XvueWindow::refuseIfBlocking() {
    if (XvueApp::blockingDepth() > 0) {
        statusBar()->showMessage(xvueT(MsgId::ModalRefuse), 3000);
        return true;
    }
    return false;
}

// ---- buildMenuBar: 13 shared QActions in {File, View, Help} ----
void XvueWindow::buildMenuBar() {
    auto* mb = menuBar();

    // File menu --------------------------------------------------------------
    auto* fileMenu = mb->addMenu(xvueT(MsgId::FileMenuTitle));
    fileMenu->setObjectName("File");

    actOpen_ = new QAction(xvueT(MsgId::FileOpen), this);
    actOpen_->setShortcut(QKeySequence::Open);
    actOpen_->setIcon(style()->standardIcon(QStyle::SP_DialogOpenButton));
    actOpen_->setStatusTip(xvueT(MsgId::FileOpenTip));
    connect(actOpen_, &QAction::triggered, this, &XvueWindow::onFileOpen);
    fileMenu->addAction(actOpen_);

    actSave_ = new QAction(xvueT(MsgId::FileSave), this);
    actSave_->setShortcut(QKeySequence::Save);
    actSave_->setIcon(style()->standardIcon(QStyle::SP_DialogSaveButton));
    actSave_->setStatusTip(xvueT(MsgId::FileSaveTip));
    connect(actSave_, &QAction::triggered, this, &XvueWindow::onFileSave);
    fileMenu->addAction(actSave_);

    actSaveAs_ = new QAction(xvueT(MsgId::FileSaveAs), this);
    actSaveAs_->setShortcut(QKeySequence::SaveAs);
    actSaveAs_->setIcon(style()->standardIcon(QStyle::SP_DialogSaveButton));
    connect(actSaveAs_, &QAction::triggered, this, &XvueWindow::onFileSaveAs);
    fileMenu->addAction(actSaveAs_);

    fileMenu->addSeparator();

    // Recent Projects submenu (UX-06, D-06). Constructed and immediately
    // refreshed from XvuePrefs::recentProjects(). 6.1+ overrides the
    // openProjectRequested handler if it needs module-specific lexicon.
    recentMenu_ = new XvueRecentProjectsMenu(menuBridge_, this);
    fileMenu->addMenu(recentMenu_);
    connect(recentMenu_, &XvueRecentProjectsMenu::openProjectRequested,
            this,        &XvueWindow::onRecentProjectRequested);

    fileMenu->addSeparator();

    // Phase 7 Plan 04 (EXPORT-02, EXPORT-05, D-04, D-13, D-14): File → Export
    // → submenu. Per CONTEXT.md D-04 the submenu is added once in the shared
    // shell builder; all five pp*_qt executables inherit it automatically. Per
    // D-13 the PDF child is raster (QPdfWriter::drawPixmap on backing_); per
    // D-14 PDF triggering is Qt-File-menu-only (no Fortran ABI extension —
    // ABI stays at 58). Plan 05 will append GIF + Capture-Animation entries
    // after the PDF child.
    exportMenu_ = fileMenu->addMenu(xvueT(MsgId::FileExport));
    exportMenu_->setObjectName(QStringLiteral("FileExport"));
    exportMenu_->setIcon(style()->standardIcon(QStyle::SP_ArrowDown));

    actExportPng_ = exportMenu_->addAction(xvueT(MsgId::FileExportPng));
    actExportPng_->setObjectName(QStringLiteral("FileExportPng"));
    connect(actExportPng_, &QAction::triggered,
            this,          &XvueWindow::onFileExportPng);

    actExportJpeg_ = exportMenu_->addAction(xvueT(MsgId::FileExportJpeg));
    actExportJpeg_->setObjectName(QStringLiteral("FileExportJpeg"));
    connect(actExportJpeg_, &QAction::triggered,
            this,           &XvueWindow::onFileExportJpeg);

    actExportPdf_ = exportMenu_->addAction(xvueT(MsgId::FileExportPdf));
    actExportPdf_->setObjectName(QStringLiteral("FileExportPdf"));
    connect(actExportPdf_, &QAction::triggered,
            this,          &XvueWindow::onFileExportPdf);

    fileMenu->addSeparator();

    actQuit_ = new QAction(xvueT(MsgId::FileQuit), this);
    actQuit_->setShortcut(QKeySequence::Quit);
    actQuit_->setIcon(style()->standardIcon(QStyle::SP_DialogCloseButton));
    actQuit_->setStatusTip(xvueT(MsgId::FileQuitTip));
    connect(actQuit_, &QAction::triggered, this, &XvueWindow::onFileQuit);
    fileMenu->addAction(actQuit_);

    // View menu --------------------------------------------------------------
    auto* viewMenu = mb->addMenu(xvueT(MsgId::ViewMenuTitle));
    viewMenu->setObjectName("View");

    actToolbarToggle_ = new QAction(xvueT(MsgId::ViewToolbar), this);
    actToolbarToggle_->setCheckable(true);
    actToolbarToggle_->setChecked(true);
    viewMenu->addAction(actToolbarToggle_);
    // Toolbar visibility connect lives in buildToolBar() (toolBar_ exists then).

    actConsoleToggle_ = new QAction(xvueT(MsgId::ViewConsole), this);
    actConsoleToggle_->setShortcut(QKeySequence("F9"));
    actConsoleToggle_->setCheckable(true);
    actConsoleToggle_->setChecked(true);
    viewMenu->addAction(actConsoleToggle_);
    // Console dock connect lives in buildConsoleDock() (consoleDock_ exists then).

    viewMenu->addSeparator();

    actZoomIn_ = new QAction(xvueT(MsgId::ViewZoomIn), this);
    actZoomIn_->setShortcut(QKeySequence::ZoomIn);
    actZoomIn_->setIcon(style()->standardIcon(QStyle::SP_ArrowUp));
    // 6.0 placeholder; 6.1+ may queue module-specific zoom lexicon.
    viewMenu->addAction(actZoomIn_);

    actZoomOut_ = new QAction(xvueT(MsgId::ViewZoomOut), this);
    actZoomOut_->setShortcut(QKeySequence::ZoomOut);
    actZoomOut_->setIcon(style()->standardIcon(QStyle::SP_ArrowDown));
    viewMenu->addAction(actZoomOut_);

    actFit_ = new QAction(xvueT(MsgId::ViewFit), this);
    actFit_->setShortcut(QKeySequence("Ctrl+0"));
    actFit_->setIcon(style()->standardIcon(QStyle::SP_BrowserReload));
    connect(actFit_, &QAction::triggered, this, &XvueWindow::onViewFit);
    viewMenu->addAction(actFit_);

    viewMenu->addSeparator();

    actPreferences_ = new QAction(xvueT(MsgId::ViewPreferences), this);
    actPreferences_->setShortcut(QKeySequence("Ctrl+,"));
    connect(actPreferences_, &QAction::triggered,
            this,            &XvueWindow::onViewPreferences);
    viewMenu->addAction(actPreferences_);

    // Help menu --------------------------------------------------------------
    auto* helpMenu = mb->addMenu(xvueT(MsgId::HelpMenuTitle));
    helpMenu->setObjectName("Help");

    actDocumentation_ = new QAction(xvueT(MsgId::HelpDocumentation), this);
    actDocumentation_->setShortcut(QKeySequence("F1"));
    actDocumentation_->setIcon(style()->standardIcon(QStyle::SP_DialogHelpButton));
    connect(actDocumentation_, &QAction::triggered,
            this,              &XvueWindow::onHelpDocumentation);
    helpMenu->addAction(actDocumentation_);

    actAbout_ = new QAction(xvueT(MsgId::HelpAbout), this);
    actAbout_->setIcon(style()->standardIcon(QStyle::SP_MessageBoxInformation));
    connect(actAbout_, &QAction::triggered, this, &XvueWindow::onHelpAbout);
    helpMenu->addAction(actAbout_);
}

// ---- buildToolBar: shares QActions with the menu bar via Qt aliasing ----
void XvueWindow::buildToolBar() {
    toolBar_ = addToolBar(QStringLiteral("main"));
    toolBar_->setObjectName(QStringLiteral("MainToolBar"));  // for saveState()
    toolBar_->setIconSize(QSize(22, 22));

    // 7 visible items + 2 separators = 9 actions (separator counts as a
    // QAction in Qt). UI-SPEC §Toolbar lists Open / Save / Sep / ZoomIn /
    // ZoomOut / Fit / Sep / Console-toggle. Adding to a toolbar references
    // the same QAction object created in buildMenuBar() — clicks on either
    // surface trigger the same handler.
    toolBar_->addAction(actOpen_);
    toolBar_->addAction(actSave_);
    toolBar_->addSeparator();
    toolBar_->addAction(actZoomIn_);
    toolBar_->addAction(actZoomOut_);
    toolBar_->addAction(actFit_);
    toolBar_->addSeparator();
    toolBar_->addAction(actConsoleToggle_);

    // Wire toolbar-toggle QAction <-> toolbar visibility (View menu's
    // toolbar checkbox now controls the toolbar).
    connect(actToolbarToggle_, &QAction::toggled,
            toolBar_,          &QToolBar::setVisible);
}

// ---- buildStatusBar: permanent coord label ----
void XvueWindow::buildStatusBar() {
    auto* sb = statusBar();
    coordLabel_ = new QLabel(sb);
    sb->addPermanentWidget(coordLabel_);
    updateStatusCoords(QPoint(0, 0));   // initial text
}

// ---- buildConsoleDock: docked at bottom, default-visible per D-07 ----
void XvueWindow::buildConsoleDock() {
    consoleDock_ = new XvueConsoleDock(this);
    consoleDock_->setObjectName("ConsoleDock");
    addDockWidget(Qt::BottomDockWidgetArea, consoleDock_);
    consoleDock_->resize(600, 140);
    // D-07 default visible on first launch; QSettings overrides on subsequent
    // launches via the consoleVisible(true)-fallback parameter.
    consoleDock_->setVisible(XvuePrefs::consoleVisible(true));

    // Wire console-toggle QAction <-> dock visibility.
    connect(actConsoleToggle_, &QAction::toggled,
            consoleDock_,      &QDockWidget::setVisible);
    // Keep the action check-state in sync when user closes via the dock's
    // X button (or any other path that toggles visibility).
    connect(consoleDock_,      &QDockWidget::visibilityChanged,
            actConsoleToggle_, &QAction::setChecked);
}

// ---- updateStatusCoords slot ----
void XvueWindow::updateStatusCoords(QPoint p) {
    if (coordLabel_) {
        coordLabel_->setText(
            xvueT(MsgId::StatusCoordFormat).arg(p.x()).arg(p.y()));
    }
}

// ============================================================================
// QAction handlers
// ============================================================================

void XvueWindow::onFileOpen() {
    if (refuseIfBlocking()) return;

    // UX-06: start in $MEFISTOX if set; warn once if unset (UI-SPEC Flag #5).
    const QByteArray mfx = qgetenv("MEFISTOX");
    if (mfx.isEmpty() && !mefistoxWarned_) {
        mefistoxWarned_ = true;
        statusBar()->showMessage(
            QStringLiteral("MEFISTOX env var unset \xE2\x80\x94 defaulting to $HOME"),
            3000);
    }
    const QString startDir = mfx.isEmpty()
        ? QString::fromLocal8Bit(qgetenv("HOME"))
        : QString::fromLocal8Bit(mfx);

    const QString path = QFileDialog::getExistingDirectory(
        this, xvueT(MsgId::OpenProjectDialogTitle), startDir);
    if (path.isEmpty()) return;
    XvuePrefs::pushRecentProject(path);
    if (recentMenu_) recentMenu_->refresh();
    // 6.1..6.5 may override this handler to queue the module-specific open
    // lexicon. Plan 06 leaves a breadcrumb here for that wiring.
}

void XvueWindow::onFileSave() {
    if (refuseIfBlocking()) return;
    // 6.0 placeholder: 6.1..6.5 wire module-specific save lexicon. The
    // shared action infrastructure is ready; the lexicon strings are
    // per-module so 6.0 emits nothing.
}

void XvueWindow::onFileSaveAs() {
    if (refuseIfBlocking()) return;
    // 6.0 placeholder — same rationale as onFileSave.
}

// ---- Phase 7 Plan 04: File → Export slot bodies ----
// Each delegates to the matching XvueExport static helper, which prompts via
// QFileDialog (QSettings-remembered last_dir) and writes the canvas backing
// pixmap to the user-selected path. refuseIfBlocking() gates so a Fortran-
// blocked picking loop cannot accidentally fire a QFileDialog re-entry.
void XvueWindow::onFileExportPng()  {
    if (refuseIfBlocking()) return;
    XvueExport::onMenuExportPng();
}
void XvueWindow::onFileExportJpeg() {
    if (refuseIfBlocking()) return;
    XvueExport::onMenuExportJpeg();
}
void XvueWindow::onFileExportPdf()  {
    if (refuseIfBlocking()) return;
    XvueExport::onMenuExportPdf();
}

void XvueWindow::onFileQuit() {
    // Phase 6.1 D-09 silent dispatch: push "99;" through the menu bridge
    // so Fortran's LIMTCL save path drains the typed-equivalent. Matches
    // the typed-lexicon `99;` which is also unconfirmed — consistent
    // click/type semantics per CONTEXT D-09. Fortran's prpr/ppmail.f
    // label 9900 calls ARRET(0) → STOP → clean process exit; no
    // QMessageBox::question confirm, no QCoreApplication::quit() bypass.
    if (menuBridge_) {
        if (consoleDock_) {
            consoleDock_->appendLine(QStringLiteral("[menu] 99;"));
        }
        menuBridge_->queueLexicon(QStringLiteral("99;"));
    }
    // Deliberately no close() / QCoreApplication::quit() — Fortran owns
    // the exit. The saclav.f LIMTCL dispatch will pick up "99;" on its
    // next xvsouris_ iteration and run the SAUVEGARDE path.
}

void XvueWindow::onViewPreferences() {
    if (refuseIfBlocking()) return;
    XvuePreferencesDialog dlg(this);
    if (dlg.exec() == QDialog::Accepted) {
        // Re-apply palette in case color-scheme changed.
        XvueApp::applyColorSchemePreference();
        // Console dock visibility persistence is handled in closeEvent;
        // no immediate refresh needed.
    }
}

void XvueWindow::onViewFit() {
    if (canvas_) canvas_->resetView();
}

void XvueWindow::onHelpDocumentation() {
    // UX-09: launch documentation via QDesktopServices. $MEFISTO/doc/ is the
    // canonical location (active doc symlinked to doca/ or docf/ at install
    // time per CLAUDE.md). If $MEFISTO is unset, surface a status message
    // rather than opening a relative URL.
    const QByteArray home = qgetenv("MEFISTO");
    if (home.isEmpty()) {
        statusBar()->showMessage(QStringLiteral("MEFISTO env var unset"), 3000);
        return;
    }
    const QString url = QStringLiteral("file://") +
                        QString::fromLocal8Bit(home) + QStringLiteral("/doc/");
    QDesktopServices::openUrl(QUrl(url));
}

void XvueWindow::onHelpAbout() {
    if (refuseIfBlocking()) return;
    XvueAboutDialog::show(this);
}

void XvueWindow::onRecentProjectRequested(const QString& path) {
    if (refuseIfBlocking()) return;
    // Re-push the path so it bubbles to the top of the LRU list.
    XvuePrefs::pushRecentProject(path);
    if (recentMenu_) recentMenu_->refresh();
    // 6.1+ queues the module-specific open lexicon here.
}
