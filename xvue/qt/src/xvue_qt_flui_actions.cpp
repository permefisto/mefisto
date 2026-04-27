// xvue/qt/src/xvue_qt_flui_actions.cpp -- Phase 6.3 Plan 02
//
// STRONG-symbol body of `registerFluiActions_stub_(XvueWindow*, XvueMenuBridge*)`.
// The weak counterpart at xvue/qt/src/xvue_qt_api.cpp:256 is displaced at
// link time when this TU is present in libxvueqt.a -- the dispatch at
// xvue_module_init_ then calls this body whenever `mod == "flui"`.
//
// Composition contract (locked by LEXICON-AUDIT-flui.md Plan 01 freeze):
//   top-level menus = {File, Fluid, View, Help} (Fluid replaces Solve from 6.2)
//   File menu : append 70; (99; Save&Quit stays on shared 6.0 actQuit_
//               which Plan 03 rewires via queueLexicon)
//   Fluid menu : flat leaves 1; plus Input/Heat/Stokes/NavierStokes submenus
//                expanding cl_flui / methrecr leaves; Parameters submenu (20;, 60;)
//   View menu : extend with 10; (DrawResults nested vitepr2d/3d/trflevit/etc) and 19;
//   Help menu : append 97;
//   toolbar   : top-5 = {2;, 5;, 7;, 10;, 99;} (99; via shared 6.0 actQuit_)
//   flat QActions for leaves; menu headers are QMenu (not QAction)
//   every QAction triggered lambda echoes `[menu] <lexicon>` to the console
//   dock BEFORE queuing the chars (mail D-07 pattern, mirrored by elas)
//
// Security posture (see 06.3-02-PLAN.md <threat_model>):
//   T-6.3-ABI-DRIFT      -- mitigated by verify_abi.sh (count stays 58)
//   T-6.3-FORCE-LINK     -- mitigated (xvue_qt_flui_actions_keepalive export
//                           referenced from xvue_qt_app.cpp static init)
//   T-6.3-MENU-ORDER     -- mitigated (Q_ASSERT_X on viewMenuForAnchor +
//                           ensureTopLevelMenu(insertBefore=...))
//   T-6.3-PARALLEL-BUILD -- mitigated by serial cbl_tout_qt + cbl_tout sequence
#include "xvue_qt_window.h"
#include "xvue_qt_menu_bridge.h"
#include "xvue_qt_console_dock.h"
#include "xvue_qt_menu_file_parser.h"

#include <QAction>
#include <QIcon>
#include <QMenu>
#include <QMenuBar>
#include <QObject>
#include <QStyle>
#include <QToolBar>

namespace {

// Builds a leaf QAction, wires [menu] echo + queueLexicon, returns it.
// Mirrors the 6.1 mail / 6.2 elas helper verbatim -- same semantics,
// same parent ownership, same dynamic-property test hook.
QAction* makeLeafAction(XvueWindow* win,
                        XvueMenuBridge* mb,
                        const QString& label,
                        const QString& lexicon,
                        const QString& iconResPath = QString(),
                        QStyle::StandardPixmap sp = QStyle::SP_CustomBase)
{
    auto* a = new QAction(label, win);
    if (!iconResPath.isEmpty()) {
        a->setIcon(QIcon(iconResPath));
    } else if (sp != QStyle::SP_CustomBase) {
        a->setIcon(win->style()->standardIcon(sp));
    }
    // Plan 03 test hook -- `actionByLexicon(menu, lex)` walks the QAction
    // tree filtering on this dynamic property.
    a->setProperty("lexicon", lexicon);

    QObject::connect(a, &QAction::triggered, [win, mb, lexicon]() {
        // Echo BEFORE queue (6.1 D-07): appendLine writes directly to
        // the QPlainTextEdit, not through the stdout redirect.
        if (auto* dock = win->consoleDock()) {
            dock->appendLine(QStringLiteral("[menu] ") + lexicon);
        }
        if (mb) mb->queueLexicon(lexicon);
    });
    return a;
}

// Helper for top-level menu creation with a stable objectName so QTest
// cases can locate the Fluid menu via `mbar->findChild<QMenu*>("Fluid")`.
// Idempotent -- repeated registration is safe.
//
// 6.2 Plan 04 generalization: when `insertBefore` is non-null, route the
// new menu through QMenuBar::insertMenu(insertBefore->menuAction(), m) so
// it lands BEFORE that anchor menu (used to position Fluid between File
// and View, matching {File, Fluid, View, Help}).
QMenu* ensureTopLevelMenu(QMenuBar* mbar,
                          const QString& title,
                          const QString& objectName,
                          QMenu* insertBefore = nullptr)
{
    if (auto* existing = mbar->findChild<QMenu*>(objectName)) {
        return existing;
    }
    auto* m = new QMenu(title, mbar);
    m->setObjectName(objectName);
    if (insertBefore) {
        mbar->insertMenu(insertBefore->menuAction(), m);
    } else {
        mbar->addMenu(m);
    }
    return m;
}

// Scan `menu` (top level + one submenu level) for a QAction whose
// `lexicon` dynamic property matches `lex`, and append it to `tb`.
// Appends the EXISTING QAction (never clones) so menu and toolbar
// share a single instance.
// NOTE: depth limit is two levels (top + one submenu). Actions nested
// more than one submenu level deep will silently not be found. All
// toolbar-targeted lexicons here (2;, 5;, 7;, 10;) are within that depth.
void addToolbarByLexicon(QToolBar* tb, QMenu* menu, const QString& lex)
{
    if (!tb || !menu) return;
    for (QAction* a : menu->actions()) {
        if (a->property("lexicon").toString() == lex) {
            tb->addAction(a);
            return;
        }
        if (a->menu()) {
            for (QAction* sub : a->menu()->actions()) {
                if (sub->property("lexicon").toString() == lex) {
                    tb->addAction(sub);
                    return;
                }
            }
        }
    }
}

} // namespace

// Force-link keepalive (Phase 6.1 Plan 03 Rule 3 auto-fix, reused per
// HANDOFF.json reusable_patterns_for_6_2_through_6_5[1]). GNU ld's
// archive-search semantics only pull an archive member when it resolves
// an otherwise-undefined reference. The weak `registerFluiActions_stub_`
// definition in xvue_qt_api.cpp.o already resolves the dispatch-site
// reference, so without this keepalive the linker NEVER pulls
// xvue_qt_flui_actions.cpp.o from libxvueqt.a -- the strong override
// below would silently be excluded and the warn-once stub ships in
// pp/ppflui_qt. XvueApp's load_bundled_font_() references this symbol
// via a static ref so the TU is always pulled.
extern "C" int xvue_qt_flui_actions_keepalive() {
    return 1;
}

// STRONG C-linkage symbol -- displaces the weak warn-once stub in
// xvue_qt_api.cpp. The dispatch at xvue_module_init_ is module-name-gated,
// so this body ONLY fires when xvue_module_init_ is called with
// name == "flui".
extern "C" void registerFluiActions_stub_(XvueWindow* win,
                                          XvueMenuBridge* mb)
{
    if (!win || !mb) return;
    // One-shot guard: registerFluiActions_stub_ must only be called once per
    // XvueWindow lifetime. xvue_module_init_ is the sole intended call site.
    // Without this guard a second call would add duplicate separators, QActions,
    // and toolbar buttons because only ensureTopLevelMenu (the Fluid-menu block)
    // is itself idempotent; the File/View extensions and toolbar block are not.
    //
    // IN-02 note: this guard is process-scoped (static bool), which is safe
    // under the current single-window-per-process invariant enforced by XvueApp.
    // If a future phase ever introduces window recycling, replace with a per-window
    // flag (e.g. win->setProperty("fluiRegistered", true)) to avoid the second
    // window silently receiving no flui actions.
    // registerMailActions_stub_ / registerElasActions_stub_ carry the same
    // guard with the same assumption.
    static bool registered = false;
    if (registered) return;
    registered = true;

    // Parse bilingual labels once per sub-menu (cached internally by
    // XvueMenuFileParser -- shared across modules). Parser falls back
    // to "N;" if a file is missing or malformed.
    const MenuFile& debuflui = XvueMenuFileParser::loadFor(QStringLiteral("debuflui"));
    const MenuFile& cl_flui  = XvueMenuFileParser::loadFor(QStringLiteral("cl_flui"));
    const MenuFile& methrecr = XvueMenuFileParser::loadFor(QStringLiteral("methrecr"));
    const MenuFile& nbvitpre = XvueMenuFileParser::loadFor(QStringLiteral("nbvitpre"));
    const MenuFile& vitepr2d = XvueMenuFileParser::loadFor(QStringLiteral("vitepr2d"));
    const MenuFile& vitepr3d = XvueMenuFileParser::loadFor(QStringLiteral("vitepr3d"));
    const MenuFile& trflevit = XvueMenuFileParser::loadFor(QStringLiteral("trflevit"));

    auto* mbar = win->menuBar();
    auto* tb   = win->toolBar();

    // ------------------------------------------------------------------
    // File menu -- extended via findChild to avoid creating a duplicate
    // top-level File menu. The 6.0 actQuit_ stays put; Plan 03 rewires
    // its handler to queueLexicon("99;").
    // ------------------------------------------------------------------
    if (auto* fileMenu = mbar->findChild<QMenu*>(QStringLiteral("File"))) {
        fileMenu->addSeparator();
        // 70; -- MANAGE TMS Files Units (debuflui leaf)
        fileMenu->addAction(makeLeafAction(win, mb,
            debuflui.label(70),
            QStringLiteral("70;"),
            QString(), QStyle::SP_DirIcon));
    }

    // ------------------------------------------------------------------
    // Fluid menu -- NEW top-level menu per UI-SPEC "Fluid / Fluide"
    // Per-Module Conformance Contract. Uses ensureTopLevelMenu, which is
    // itself idempotent (returns an existing QMenu for the same objectName),
    // but the outer registerFluiActions_stub_ is guarded above to run once.
    //
    // 6.2 Plan 04 Gap-1 (carried forward): look up the View menu so we
    // can insert Fluid BEFORE it, yielding the {File, Fluid, View, Help}
    // sequence locked by ROADMAP Phase 6.3 + 06.0 Per-Module Conformance
    // Contract. 06.0 chrome added File/View/Help in that index order;
    // this line repositions Fluid between File and View.
    // ------------------------------------------------------------------
    auto* viewMenuForAnchor = mbar->findChild<QMenu*>(QStringLiteral("View"));
    // 6.2 WR-03 gap-closure (mirrored): assert that the View anchor is
    // present so a future chrome-init reordering fails loudly in debug
    // builds rather than silently producing {File, Help, Fluid} instead
    // of {File, Fluid, View, Help}.
    Q_ASSERT_X(viewMenuForAnchor != nullptr,
               "registerFluiActions_stub_",
               "View menu not found — Fluid will be appended after Help "
               "(expected {File, Fluid, View, Help} order); "
               "chrome init order may have changed.");
    auto* fluidMenu = ensureTopLevelMenu(mbar,
        QObject::tr("&Fluid"),
        QStringLiteral("Fluid"),
        viewMenuForAnchor);

    // Fluid > Object name (leaf 1;)
    fluidMenu->addAction(makeLeafAction(win, mb,
        debuflui.label(1), QStringLiteral("1;"),
        QStringLiteral(":/xvue/qt/icons/icons/flui/flui-object.svg")));

    // Fluid > Input physical data (2; parent) -- cl_flui submenu holds
    // 2;1; .. 2;5; leaves; only 2;1; 2;2; 2;3; are qaction=yes per audit.
    auto* inputMenu = fluidMenu->addMenu(debuflui.label(2));
    inputMenu->setObjectName(QStringLiteral("Fluid.Input"));
    inputMenu->setIcon(QIcon(
        QStringLiteral(":/xvue/qt/icons/icons/flui/flui-input-physics.svg")));
    // Parent itself is also a runnable lexicon for users who want "just 2;"
    inputMenu->addAction(makeLeafAction(win, mb,
        QObject::tr("Open input physical data dialog..."),
        QStringLiteral("2;"),
        QStringLiteral(":/xvue/qt/icons/icons/flui/flui-input-physics.svg")));
    inputMenu->addSeparator();
    inputMenu->addAction(makeLeafAction(win, mb,
        cl_flui.label(1), QStringLiteral("2;1;")));
    inputMenu->addAction(makeLeafAction(win, mb,
        cl_flui.label(2), QStringLiteral("2;2;")));
    inputMenu->addAction(makeLeafAction(win, mb,
        cl_flui.label(3), QStringLiteral("2;3;")));

    // Fluid > Input heat data (3; flat leaf)
    fluidMenu->addAction(makeLeafAction(win, mb,
        debuflui.label(3), QStringLiteral("3;"),
        QStringLiteral(":/xvue/qt/icons/icons/flui/flui-input-heat.svg")));

    // Fluid > Steady Stokes (5; parent) -- methrecr leaves under it
    auto* stokesSteadyMenu = fluidMenu->addMenu(debuflui.label(5));
    stokesSteadyMenu->setObjectName(QStringLiteral("Fluid.StokesSteady"));
    stokesSteadyMenu->setIcon(QIcon(
        QStringLiteral(":/xvue/qt/icons/icons/flui/solve-stokes.svg")));
    stokesSteadyMenu->addAction(makeLeafAction(win, mb,
        QObject::tr("Run steady Stokes..."),
        QStringLiteral("5;"),
        QStringLiteral(":/xvue/qt/icons/icons/flui/solve-stokes.svg")));
    stokesSteadyMenu->addSeparator();
    // 5;81;1; -- methrecr leaf 1 (Crout factorization), synthetic prefix
    // per audit D-04. Real Fortran dispatch reads methrecr.label(N) within
    // STOKESTA's CHMERESO call.
    stokesSteadyMenu->addAction(makeLeafAction(win, mb,
        methrecr.label(1), QStringLiteral("5;81;1;")));
    // 5;81;2; -- methrecr leaf 2 (Conjugate Gradient)
    stokesSteadyMenu->addAction(makeLeafAction(win, mb,
        methrecr.label(2), QStringLiteral("5;81;2;")));

    // Fluid > Unsteady Stokes (6; flat leaf -- icon REUSED from 5;)
    fluidMenu->addAction(makeLeafAction(win, mb,
        debuflui.label(6), QStringLiteral("6;"),
        QStringLiteral(":/xvue/qt/icons/icons/flui/solve-stokes.svg")));

    // Fluid > Implicit Navier-Stokes (7; parent) -- canonical NS solver
    auto* nsImplicitMenu = fluidMenu->addMenu(debuflui.label(7));
    nsImplicitMenu->setObjectName(QStringLiteral("Fluid.NSImplicit"));
    nsImplicitMenu->setIcon(QIcon(
        QStringLiteral(":/xvue/qt/icons/icons/flui/solve-navier-stokes.svg")));
    nsImplicitMenu->addAction(makeLeafAction(win, mb,
        QObject::tr("Run implicit Navier-Stokes..."),
        QStringLiteral("7;"),
        QStringLiteral(":/xvue/qt/icons/icons/flui/solve-navier-stokes.svg")));

    // Fluid > Fractional Navier-Stokes (8; flat leaf -- icon REUSED)
    fluidMenu->addAction(makeLeafAction(win, mb,
        debuflui.label(8), QStringLiteral("8;"),
        QStringLiteral(":/xvue/qt/icons/icons/flui/solve-navier-stokes.svg")));

    // Fluid > Boussinesq Heat Navier-Stokes (9; flat leaf -- icon REUSED)
    fluidMenu->addAction(makeLeafAction(win, mb,
        debuflui.label(9), QStringLiteral("9;"),
        QStringLiteral(":/xvue/qt/icons/icons/flui/solve-navier-stokes.svg")));

    // Fluid > Parameters submenu (audit's codes 20; 60;)
    auto* paramsMenu = fluidMenu->addMenu(QObject::tr("&Parameters"));
    paramsMenu->setObjectName(QStringLiteral("Fluid.Parameters"));
    paramsMenu->addAction(makeLeafAction(win, mb,
        debuflui.label(20), QStringLiteral("20;")));
    paramsMenu->addAction(makeLeafAction(win, mb,
        debuflui.label(60), QStringLiteral("60;"),
        QString(), QStyle::SP_ComputerIcon));

    // ------------------------------------------------------------------
    // View menu -- extended via findChild (Pitfall 7 inherited from 6.1).
    // Adds 10; (Draw results, deepest reach via TFLUIDE) and 19; (Draw mesh).
    // ------------------------------------------------------------------
    if (auto* viewMenu = mbar->findChild<QMenu*>(QStringLiteral("View"))) {
        viewMenu->addSeparator();

        // 10; Draw velocities/pressures/temperatures -- TFLUIDE driver.
        // QMenu parent + parent-shortcut leaf + selected qaction=yes leaves
        // from nbvitpre / vitepr2d / vitepr3d / trflevit.
        auto* drawResultsMenu = viewMenu->addMenu(debuflui.label(10));
        drawResultsMenu->setObjectName(QStringLiteral("View.DrawResults"));
        drawResultsMenu->setIcon(QIcon(
            QStringLiteral(":/xvue/qt/icons/icons/flui/draw-velocity.svg")));
        drawResultsMenu->addAction(makeLeafAction(win, mb,
            QObject::tr("Draw results..."),
            QStringLiteral("10;"),
            QStringLiteral(":/xvue/qt/icons/icons/flui/draw-velocity.svg")));
        drawResultsMenu->addSeparator();

        // 10;1;{1,2}; -- nbvitpre leaves (numbers vs. times for selection)
        drawResultsMenu->addAction(makeLeafAction(win, mb,
            nbvitpre.label(1), QStringLiteral("10;1;1;")));
        drawResultsMenu->addAction(makeLeafAction(win, mb,
            nbvitpre.label(2), QStringLiteral("10;1;2;")));

        // View.DrawResults.2D submenu -- vitepr2d leaves (qaction=yes subset)
        auto* draw2dMenu = drawResultsMenu->addMenu(QObject::tr("2D drawings"));
        draw2dMenu->setObjectName(QStringLiteral("View.DrawResults.2D"));
        draw2dMenu->addAction(makeLeafAction(win, mb,
            vitepr2d.label(1), QStringLiteral("10;20;1;")));
        draw2dMenu->addAction(makeLeafAction(win, mb,
            vitepr2d.label(2), QStringLiteral("10;20;2;")));
        draw2dMenu->addAction(makeLeafAction(win, mb,
            vitepr2d.label(3), QStringLiteral("10;20;3;")));
        draw2dMenu->addAction(makeLeafAction(win, mb,
            vitepr2d.label(5), QStringLiteral("10;20;5;")));
        draw2dMenu->addAction(makeLeafAction(win, mb,
            vitepr2d.label(9), QStringLiteral("10;20;9;")));
        draw2dMenu->addAction(makeLeafAction(win, mb,
            vitepr2d.label(10), QStringLiteral("10;20;10;")));
        draw2dMenu->addAction(makeLeafAction(win, mb,
            vitepr2d.label(90), QStringLiteral("10;20;90;")));

        // View.DrawResults.3D submenu -- vitepr3d leaves (qaction=yes subset)
        auto* draw3dMenu = drawResultsMenu->addMenu(QObject::tr("3D drawings"));
        draw3dMenu->setObjectName(QStringLiteral("View.DrawResults.3D"));
        draw3dMenu->addAction(makeLeafAction(win, mb,
            vitepr3d.label(1), QStringLiteral("10;30;1;")));
        draw3dMenu->addAction(makeLeafAction(win, mb,
            vitepr3d.label(2), QStringLiteral("10;30;2;")));
        draw3dMenu->addAction(makeLeafAction(win, mb,
            vitepr3d.label(3), QStringLiteral("10;30;3;")));
        draw3dMenu->addAction(makeLeafAction(win, mb,
            vitepr3d.label(5), QStringLiteral("10;30;5;")));
        draw3dMenu->addAction(makeLeafAction(win, mb,
            vitepr3d.label(10), QStringLiteral("10;30;10;")));
        draw3dMenu->addAction(makeLeafAction(win, mb,
            vitepr3d.label(90), QStringLiteral("10;30;90;")));

        // View.DrawResults.Arrows submenu -- trflevit leaves
        // (velocity-arrow scale + execute trigger)
        auto* arrowsMenu = drawResultsMenu->addMenu(
            QObject::tr("Velocity arrows"));
        arrowsMenu->setObjectName(QStringLiteral("View.DrawResults.Arrows"));
        arrowsMenu->addAction(makeLeafAction(win, mb,
            trflevit.label(2), QStringLiteral("10;80;2;")));
        arrowsMenu->addAction(makeLeafAction(win, mb,
            trflevit.label(3), QStringLiteral("10;80;3;")));
        arrowsMenu->addAction(makeLeafAction(win, mb,
            trflevit.label(90), QStringLiteral("10;80;90;")));

        // 19; Draw PLSVO mesh (REUSES 6.1 mail mesh-draw.svg -- shared
        // qrc prefix resolves icons/mail/ from anywhere).
        viewMenu->addAction(makeLeafAction(win, mb,
            debuflui.label(19), QStringLiteral("19;"),
            QStringLiteral(":/xvue/qt/icons/icons/mail/mesh-draw.svg")));
    }

    // ------------------------------------------------------------------
    // Help menu -- extended via findChild for the version banner.
    // ------------------------------------------------------------------
    if (auto* helpMenu = mbar->findChild<QMenu*>(QStringLiteral("Help"))) {
        helpMenu->addSeparator();
        helpMenu->addAction(makeLeafAction(win, mb,
            debuflui.label(97), QStringLiteral("97;"),
            QString(), QStyle::SP_MessageBoxInformation));
    }

    // ------------------------------------------------------------------
    // Toolbar top-5 per LEXICON-AUDIT-flui.md toolbar=yes rows:
    // {2;, 5;, 7;, 10;, 99;}. Four are wired here by locating the
    // existing QActions on fluidMenu / viewMenu via addToolbarByLexicon.
    // 99; Save&Quit is owned by the shared 6.0 actQuit_ action which
    // Plan 03 rewires via queueLexicon("99;"); no toolbar entry here.
    // ------------------------------------------------------------------
    tb->addSeparator();
    addToolbarByLexicon(tb, fluidMenu, QStringLiteral("2;"));   // Input data
    addToolbarByLexicon(tb, fluidMenu, QStringLiteral("5;"));   // Steady Stokes
    addToolbarByLexicon(tb, fluidMenu, QStringLiteral("7;"));   // Implicit NS
    if (auto* viewMenu = mbar->findChild<QMenu*>(QStringLiteral("View"))) {
        addToolbarByLexicon(tb, viewMenu, QStringLiteral("10;"));  // Draw results
    }
}
