// xvue/qt/src/xvue_qt_ther_actions.cpp -- Phase 6.4 Plan 02
//
// STRONG-symbol body of `registerTherActions_stub_(XvueWindow*, XvueMenuBridge*)`.
// The weak counterpart at xvue/qt/src/xvue_qt_api.cpp:263 is displaced at
// link time when this TU is present in libxvueqt.a -- the dispatch at
// xvue_module_init_ then calls this body whenever `mod == "ther"`.
//
// Composition contract (locked by LEXICON-AUDIT-ther.md Plan 01 freeze):
//   top-level menus = {File, Thermal, View, Help} (Thermal replaces Fluid
//                      from 6.3 / Solve from 6.2 / Mesh from 6.1)
//   File menu  : append 0; (Number of INPUT DATA GAME) + 70; (MS File Tools)
//                 -- 99; Save&Quit stays on shared 6.0 actQuit_ (Plan 03 wires
//                    queueLexicon)
//   Thermal menu : flat 1; (Object), 2; HEAT INPUT submenu (cl_ther leaves),
//                  3; STEADY HEAT submenu (resothst + methreso), 4; UNSTEADY
//                  HEAT submenu (resothin + methreso), 5; WAVE solver flat,
//                  6; LINEAR EIGEN flat, 12; nonlinear elliptic flat,
//                  19; POLY EIGEN flat (DEAD CODE per ppther.f:451 -- still
//                  wired for menu completeness), Parameters submenu
//                  (20; precision, 38; window pixels, 39; background)
//   View menu : extend with 7; (Draw eigenvectors), 8; (Draw temperatures
//               + flux) submenu expanding tempgrad/tractemp/tracflux/tracgrad
//               leaves, 9; (Draw wave), 10; (Draw mesh -- REUSES 6.1 mail
//               mesh-draw.svg), 11; (Draw iso-surfaces) submenu expanding
//               valisoth leaves, 13; Energies, 16; L1 norm
//   Help menu : append 98; (Mefisto Date Version -- NOTE: 98; not 97; like
//               flui; Plan 03's testHelpNoQueue allowlist is {98;})
//   toolbar   : top-5 = {2;, 3;, 4;, 8;, 99;} (99; via shared 6.0 actQuit_)
//   flat QActions for leaves; menu headers are QMenu (not QAction)
//   every QAction triggered lambda echoes `[menu] <lexicon>` to the console
//   dock BEFORE queuing the chars (mail D-07 pattern, mirrored by elas/flui)
//
// Security posture (see 06.4-02-PLAN.md <threat_model>):
//   T-6.4-ABI-DRIFT      -- mitigated by verify_abi.sh (count stays 58)
//   T-6.4-FORCE-LINK     -- mitigated (xvue_qt_ther_actions_keepalive export
//                           referenced from xvue_qt_app.cpp static init)
//   T-6.4-MENU-ORDER     -- mitigated (Q_ASSERT_X on viewMenuForAnchor +
//                           ensureTopLevelMenu(insertBefore=...))
//   T-6.4-PARALLEL-BUILD -- mitigated by serial cbl_tout_qt + cbl_tout sequence
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
// Mirrors the 6.1 mail / 6.2 elas / 6.3 flui helper verbatim -- same
// semantics, same parent ownership, same dynamic-property test hook.
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
// cases can locate the Thermal menu via `mbar->findChild<QMenu*>("Thermal")`.
// Idempotent -- repeated registration is safe.
//
// 6.2 Plan 04 generalization: when `insertBefore` is non-null, route the
// new menu through QMenuBar::insertMenu(insertBefore->menuAction(), m) so
// it lands BEFORE that anchor menu (used to position Thermal between File
// and View, matching {File, Thermal, View, Help}).
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
// toolbar-targeted lexicons here (2;, 3;, 4;, 8;) are within that depth.
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
// an otherwise-undefined reference. The weak `registerTherActions_stub_`
// definition in xvue_qt_api.cpp.o already resolves the dispatch-site
// reference, so without this keepalive the linker NEVER pulls
// xvue_qt_ther_actions.cpp.o from libxvueqt.a -- the strong override
// below would silently be excluded and the warn-once stub ships in
// pp/ppther_qt. XvueApp's load_bundled_font_() references this symbol
// via a static ref so the TU is always pulled.
extern "C" int xvue_qt_ther_actions_keepalive() {
    return 1;
}

// STRONG C-linkage symbol -- displaces the weak warn-once stub in
// xvue_qt_api.cpp. The dispatch at xvue_module_init_ is module-name-gated,
// so this body ONLY fires when xvue_module_init_ is called with
// name == "ther".
extern "C" void registerTherActions_stub_(XvueWindow* win,
                                          XvueMenuBridge* mb)
{
    if (!win || !mb) return;
    // One-shot guard: registerTherActions_stub_ must only be called once per
    // XvueWindow lifetime. xvue_module_init_ is the sole intended call site.
    // Without this guard a second call would add duplicate separators, QActions,
    // and toolbar buttons because only ensureTopLevelMenu (the Thermal-menu block)
    // is itself idempotent; the File/View extensions and toolbar block are not.
    //
    // IN-02 note: this guard is process-scoped (static bool), which is safe
    // under the current single-window-per-process invariant enforced by XvueApp.
    // If a future phase ever introduces window recycling, replace with a per-window
    // flag (e.g. win->setProperty("therRegistered", true)) to avoid the second
    // window silently receiving no ther actions.
    // registerMailActions_stub_ / registerElasActions_stub_ /
    // registerFluiActions_stub_ carry the same guard with the same assumption.
    static bool registered = false;
    if (registered) return;
    registered = true;

    // Parse bilingual labels once per sub-menu (cached internally by
    // XvueMenuFileParser -- shared across modules). Parser falls back
    // to "N;" if a file is missing or malformed.
    const MenuFile& debuther = XvueMenuFileParser::loadFor(QStringLiteral("debuther"));
    const MenuFile& cl_ther  = XvueMenuFileParser::loadFor(QStringLiteral("cl_ther"));
    const MenuFile& methreso = XvueMenuFileParser::loadFor(QStringLiteral("methreso"));
    const MenuFile& resothst = XvueMenuFileParser::loadFor(QStringLiteral("resothst"));
    const MenuFile& resothin = XvueMenuFileParser::loadFor(QStringLiteral("resothin"));
    const MenuFile& scheonde = XvueMenuFileParser::loadFor(QStringLiteral("scheonde"));
    const MenuFile& tempgrad = XvueMenuFileParser::loadFor(QStringLiteral("tempgrad"));
    const MenuFile& tractemp = XvueMenuFileParser::loadFor(QStringLiteral("tractemp"));
    const MenuFile& tractem1 = XvueMenuFileParser::loadFor(QStringLiteral("tractem1"));
    const MenuFile& tracflux = XvueMenuFileParser::loadFor(QStringLiteral("tracflux"));
    const MenuFile& tracgrad = XvueMenuFileParser::loadFor(QStringLiteral("tracgrad"));
    const MenuFile& valisoth = XvueMenuFileParser::loadFor(QStringLiteral("valisoth"));

    auto* mbar = win->menuBar();
    auto* tb   = win->toolBar();

    // ------------------------------------------------------------------
    // File menu -- extended via findChild to avoid creating a duplicate
    // top-level File menu. The 6.0 actQuit_ stays put; Plan 03 rewires
    // its handler to queueLexicon("99;").
    // ------------------------------------------------------------------
    if (auto* fileMenu = mbar->findChild<QMenu*>(QStringLiteral("File"))) {
        fileMenu->addSeparator();
        // 0; -- Number of INPUT DATA GAME (debuther leaf; niche parameter)
        fileMenu->addAction(makeLeafAction(win, mb,
            debuther.label(0),
            QStringLiteral("0;"),
            QString(), QStyle::SP_FileIcon));
        // 70; -- MANAGEMENT TMS Files Units (debuther leaf; bilingual parser
        // falls through to FR per 6.2 Plan 05 fix because td/ma/debuther
        // line 21 EN translation is missing)
        fileMenu->addAction(makeLeafAction(win, mb,
            debuther.label(70),
            QStringLiteral("70;"),
            QString(), QStyle::SP_DirIcon));
    }

    // ------------------------------------------------------------------
    // Thermal menu -- NEW top-level menu per UI-SPEC "Thermal / Thermique"
    // Per-Module Conformance Contract. Uses ensureTopLevelMenu, which is
    // itself idempotent (returns an existing QMenu for the same objectName),
    // but the outer registerTherActions_stub_ is guarded above to run once.
    //
    // 6.2 Plan 04 Gap-1 (carried forward): look up the View menu so we
    // can insert Thermal BEFORE it, yielding the {File, Thermal, View, Help}
    // sequence locked by ROADMAP Phase 6.4 line 213 + 06.0 Per-Module
    // Conformance Contract. 06.0 chrome added File/View/Help in that index
    // order; this line repositions Thermal between File and View.
    // ------------------------------------------------------------------
    auto* viewMenuForAnchor = mbar->findChild<QMenu*>(QStringLiteral("View"));
    // 6.2 WR-03 gap-closure (mirrored): assert that the View anchor is
    // present so a future chrome-init reordering fails loudly in debug
    // builds rather than silently producing {File, Help, Thermal} instead
    // of {File, Thermal, View, Help}.
    Q_ASSERT_X(viewMenuForAnchor != nullptr,
               "registerTherActions_stub_",
               "View menu not found — Thermal will be appended after Help "
               "(expected {File, Thermal, View, Help} order); "
               "chrome init order may have changed.");
    auto* thermalMenu = ensureTopLevelMenu(mbar,
        QObject::tr("&Thermal"),
        QStringLiteral("Thermal"),
        viewMenuForAnchor);

    // Thermal > Object name (leaf 1;)
    thermalMenu->addAction(makeLeafAction(win, mb,
        debuther.label(1), QStringLiteral("1;"),
        QStringLiteral(":/xvue/qt/icons/icons/ther/ther-object.svg")));

    // Thermal > Heat transfer input data (2; parent) -- cl_ther submenu holds
    // 2;1; .. 2;9; leaves; 2;1; 2;2; 2;3; 2;8; are qaction=yes per audit.
    auto* inputMenu = thermalMenu->addMenu(debuther.label(2));
    inputMenu->setObjectName(QStringLiteral("Thermal.Input"));
    inputMenu->setIcon(QIcon(
        QStringLiteral(":/xvue/qt/icons/icons/ther/ther-input-heat.svg")));
    // Parent itself is also a runnable lexicon for users who want "just 2;"
    inputMenu->addAction(makeLeafAction(win, mb,
        QObject::tr("Open heat transfer input dialog..."),
        QStringLiteral("2;"),
        QStringLiteral(":/xvue/qt/icons/icons/ther/ther-input-heat.svg")));
    inputMenu->addSeparator();
    inputMenu->addAction(makeLeafAction(win, mb,
        cl_ther.label(1), QStringLiteral("2;1;")));
    inputMenu->addAction(makeLeafAction(win, mb,
        cl_ther.label(2), QStringLiteral("2;2;")));
    inputMenu->addAction(makeLeafAction(win, mb,
        cl_ther.label(3), QStringLiteral("2;3;")));
    inputMenu->addAction(makeLeafAction(win, mb,
        cl_ther.label(8), QStringLiteral("2;8;")));

    // Thermal > Steady heat (3; parent) -- resothst leaves + methreso
    // factorisation choices under it.
    auto* steadyMenu = thermalMenu->addMenu(debuther.label(3));
    steadyMenu->setObjectName(QStringLiteral("Thermal.Steady"));
    steadyMenu->setIcon(QIcon(
        QStringLiteral(":/xvue/qt/icons/icons/ther/solve-heat-steady.svg")));
    steadyMenu->addAction(makeLeafAction(win, mb,
        QObject::tr("Run steady heat solver..."),
        QStringLiteral("3;"),
        QStringLiteral(":/xvue/qt/icons/icons/ther/solve-heat-steady.svg")));
    steadyMenu->addSeparator();
    // 3;1; -- coefficient choice (linear)
    steadyMenu->addAction(makeLeafAction(win, mb,
        resothst.label(1), QStringLiteral("3;1;")));
    // 3;2; -- coefficient choice (T-dependent)
    steadyMenu->addAction(makeLeafAction(win, mb,
        resothst.label(2), QStringLiteral("3;2;")));
    // 3;81;1; -- methreso leaf 1 (Cholesky/Crout direct), synthetic prefix
    steadyMenu->addAction(makeLeafAction(win, mb,
        methreso.label(1), QStringLiteral("3;81;1;")));
    // 3;81;2; -- methreso leaf 2 (Conjugate Gradient)
    steadyMenu->addAction(makeLeafAction(win, mb,
        methreso.label(2), QStringLiteral("3;81;2;")));

    // Thermal > Unsteady heat (4; parent) -- resothin leaves + methreso
    auto* unsteadyMenu = thermalMenu->addMenu(debuther.label(4));
    unsteadyMenu->setObjectName(QStringLiteral("Thermal.Unsteady"));
    unsteadyMenu->setIcon(QIcon(
        QStringLiteral(":/xvue/qt/icons/icons/ther/solve-heat-unsteady.svg")));
    unsteadyMenu->addAction(makeLeafAction(win, mb,
        QObject::tr("Run unsteady heat solver..."),
        QStringLiteral("4;"),
        QStringLiteral(":/xvue/qt/icons/icons/ther/solve-heat-unsteady.svg")));
    unsteadyMenu->addSeparator();
    // 4;1; -- linear-time-dependent unsteady (canonical)
    unsteadyMenu->addAction(makeLeafAction(win, mb,
        resothin.label(1), QStringLiteral("4;1;")));
    // 4;2; -- nonlinear-source unsteady
    unsteadyMenu->addAction(makeLeafAction(win, mb,
        resothin.label(2), QStringLiteral("4;2;")));

    // Thermal > 2D wave solver (5; parent) -- scheonde leaves under it
    auto* waveMenu = thermalMenu->addMenu(debuther.label(5));
    waveMenu->setObjectName(QStringLiteral("Thermal.Wave"));
    waveMenu->setIcon(QIcon(
        QStringLiteral(":/xvue/qt/icons/icons/ther/solve-wave.svg")));
    waveMenu->addAction(makeLeafAction(win, mb,
        QObject::tr("Run 2D wave solver..."),
        QStringLiteral("5;"),
        QStringLiteral(":/xvue/qt/icons/icons/ther/solve-wave.svg")));
    waveMenu->addSeparator();
    // 5;81;1; -- scheonde leaf 1 (Newmark implicit)
    waveMenu->addAction(makeLeafAction(win, mb,
        scheonde.label(1), QStringLiteral("5;81;1;")));

    // Thermal > Eigenvalue linear solver (6; flat leaf -- REUSES 6.2 elas
    // solve-eigen.svg)
    thermalMenu->addAction(makeLeafAction(win, mb,
        debuther.label(6), QStringLiteral("6;"),
        QStringLiteral(":/xvue/qt/icons/icons/elas/solve-eigen.svg")));

    // Thermal > Nonlinear elliptic (12; flat leaf, no icon -- research-grade)
    thermalMenu->addAction(makeLeafAction(win, mb,
        debuther.label(12), QStringLiteral("12;"),
        QString()));

    // Thermal > Polynomial eigenvalue (19; flat leaf -- REUSES solve-eigen.svg
    // even though ppther.f:451 dispatch is COMMENTED OUT (DEAD CODE at runtime).
    // Wired for menu completeness; user clicks fire the lexicon but the Fortran
    // side is a no-op.)
    thermalMenu->addAction(makeLeafAction(win, mb,
        debuther.label(19), QStringLiteral("19;"),
        QStringLiteral(":/xvue/qt/icons/icons/elas/solve-eigen.svg")));

    // Thermal > Parameters submenu (audit's codes 20; 38; 39;)
    auto* paramsMenu = thermalMenu->addMenu(QObject::tr("&Parameters"));
    paramsMenu->setObjectName(QStringLiteral("Thermal.Parameters"));
    paramsMenu->addAction(makeLeafAction(win, mb,
        debuther.label(20), QStringLiteral("20;"),
        QString(), QStyle::SP_DialogApplyButton));
    paramsMenu->addAction(makeLeafAction(win, mb,
        debuther.label(38), QStringLiteral("38;"),
        QString(), QStyle::SP_ComputerIcon));
    paramsMenu->addAction(makeLeafAction(win, mb,
        debuther.label(39), QStringLiteral("39;")));

    // ------------------------------------------------------------------
    // View menu -- extended via findChild (Pitfall 7 inherited from 6.1).
    // Adds 7; (Draw eigenvectors), 8; (Draw temperatures + flux) submenu,
    // 9; (Draw wave), 10; (Draw mesh), 11; (Draw iso-surfaces) submenu,
    // 13; Energies, 16; L1 norm.
    // ------------------------------------------------------------------
    if (auto* viewMenu = mbar->findChild<QMenu*>(QStringLiteral("View"))) {
        viewMenu->addSeparator();

        // 7; Draw eigenvectors -- REUSES 6.2 elas solve-eigen.svg (same family)
        viewMenu->addAction(makeLeafAction(win, mb,
            debuther.label(7),
            QStringLiteral("7;"),
            QStringLiteral(":/xvue/qt/icons/icons/elas/solve-eigen.svg")));

        // 8; Draw temperatures + flux -- TRTHER driver. QMenu parent +
        // parent-shortcut leaf + qaction=yes leaves from tempgrad / tractemp
        // / tractem1 / tracflux / tracgrad.
        auto* drawTempMenu = viewMenu->addMenu(debuther.label(8));
        drawTempMenu->setObjectName(QStringLiteral("View.DrawTemperature"));
        drawTempMenu->setIcon(QIcon(
            QStringLiteral(":/xvue/qt/icons/icons/ther/draw-temperature.svg")));
        drawTempMenu->addAction(makeLeafAction(win, mb,
            QObject::tr("Draw temperatures + flux..."),
            QStringLiteral("8;"),
            QStringLiteral(":/xvue/qt/icons/icons/ther/draw-temperature.svg")));
        drawTempMenu->addSeparator();

        // 7;1; .. 7;4; -- tempgrad leaves (canonical drawing dispatcher)
        drawTempMenu->addAction(makeLeafAction(win, mb,
            tempgrad.label(1), QStringLiteral("7;1;")));
        drawTempMenu->addAction(makeLeafAction(win, mb,
            tempgrad.label(2), QStringLiteral("7;2;")));
        drawTempMenu->addAction(makeLeafAction(win, mb,
            tempgrad.label(3), QStringLiteral("7;3;")));
        drawTempMenu->addAction(makeLeafAction(win, mb,
            tempgrad.label(4), QStringLiteral("7;4;")));

        // View.DrawTemperature.Tractemp submenu -- 2D-temperature canonical
        // (tractemp leaves; qaction=yes subset 1; 2;)
        auto* tractempMenu = drawTempMenu->addMenu(QObject::tr("2D temperatures"));
        tractempMenu->setObjectName(QStringLiteral("View.DrawTemperature.2D"));
        tractempMenu->addAction(makeLeafAction(win, mb,
            tractemp.label(1), QStringLiteral("8;20;1;")));
        tractempMenu->addAction(makeLeafAction(win, mb,
            tractemp.label(2), QStringLiteral("8;20;2;")));

        // View.DrawTemperature.Tractem1 submenu -- 1D-time temperature
        // (tractem1 leaves; qaction=yes subset 1; 2;)
        auto* tractem1Menu = drawTempMenu->addMenu(QObject::tr("1D-time temperatures"));
        tractem1Menu->setObjectName(QStringLiteral("View.DrawTemperature.1DTime"));
        tractem1Menu->addAction(makeLeafAction(win, mb,
            tractem1.label(1), QStringLiteral("8;30;1;")));
        tractem1Menu->addAction(makeLeafAction(win, mb,
            tractem1.label(2), QStringLiteral("8;30;2;")));

        // View.DrawTemperature.Flux submenu -- normal flux drawing
        // (tracflux leaves; qaction=yes 1; 15;)
        auto* fluxMenu = drawTempMenu->addMenu(QObject::tr("Normal flux"));
        fluxMenu->setObjectName(QStringLiteral("View.DrawTemperature.Flux"));
        fluxMenu->addAction(makeLeafAction(win, mb,
            tracflux.label(1), QStringLiteral("8;60;1;")));
        fluxMenu->addAction(makeLeafAction(win, mb,
            tracflux.label(15), QStringLiteral("8;60;15;")));

        // View.DrawTemperature.Gradient submenu -- gradient drawing
        // (tracgrad leaves; qaction=yes 1; 15;)
        auto* gradMenu = drawTempMenu->addMenu(QObject::tr("Gradient"));
        gradMenu->setObjectName(QStringLiteral("View.DrawTemperature.Gradient"));
        gradMenu->addAction(makeLeafAction(win, mb,
            tracgrad.label(1), QStringLiteral("8;70;1;")));
        gradMenu->addAction(makeLeafAction(win, mb,
            tracgrad.label(15), QStringLiteral("8;70;15;")));

        // 9; Draw wave -- REUSES the new ther solve-wave.svg
        viewMenu->addAction(makeLeafAction(win, mb,
            debuther.label(9), QStringLiteral("9;"),
            QStringLiteral(":/xvue/qt/icons/icons/ther/solve-wave.svg")));

        // 10; Draw PLSVO mesh (REUSES 6.1 mail mesh-draw.svg -- shared
        // qrc prefix resolves icons/mail/ from anywhere).
        viewMenu->addAction(makeLeafAction(win, mb,
            debuther.label(10), QStringLiteral("10;"),
            QStringLiteral(":/xvue/qt/icons/icons/mail/mesh-draw.svg")));

        // 11; Draw iso-surfaces -- valisoth leaves under it (qaction=yes 2;)
        auto* drawIsoMenu = viewMenu->addMenu(debuther.label(11));
        drawIsoMenu->setObjectName(QStringLiteral("View.DrawIsoSurfaces"));
        drawIsoMenu->setIcon(QIcon(
            QStringLiteral(":/xvue/qt/icons/icons/ther/draw-iso-surface.svg")));
        drawIsoMenu->addAction(makeLeafAction(win, mb,
            QObject::tr("Draw iso-surfaces..."),
            QStringLiteral("11;"),
            QStringLiteral(":/xvue/qt/icons/icons/ther/draw-iso-surface.svg")));
        drawIsoMenu->addSeparator();
        // 11;2; -- regular iso between min and max (canonical)
        drawIsoMenu->addAction(makeLeafAction(win, mb,
            valisoth.label(2), QStringLiteral("11;2;")));

        // 13; Energies of solutions (no icon -- niche post-processing)
        viewMenu->addAction(makeLeafAction(win, mb,
            debuther.label(13), QStringLiteral("13;"),
            QString()));

        // 16; L1 norm (no icon -- niche post-processing)
        viewMenu->addAction(makeLeafAction(win, mb,
            debuther.label(16), QStringLiteral("16;"),
            QString()));
    }

    // ------------------------------------------------------------------
    // Help menu -- extended via findChild for the version banner.
    // CRITICAL: ther uses code 98; (NOT 97; like flui). Plan 03's
    // testHelpNoQueue allowlist is {98;} (per Auto-fix Rule 1 lesson
    // from 6.3 -- per-module Help-allowlist drawn from LEXICON-AUDIT).
    // ------------------------------------------------------------------
    if (auto* helpMenu = mbar->findChild<QMenu*>(QStringLiteral("Help"))) {
        helpMenu->addSeparator();
        helpMenu->addAction(makeLeafAction(win, mb,
            debuther.label(98), QStringLiteral("98;"),
            QString(), QStyle::SP_MessageBoxInformation));
    }

    // ------------------------------------------------------------------
    // Toolbar top-5 per LEXICON-AUDIT-ther.md toolbar=yes rows:
    // {2;, 3;, 4;, 8;, 99;}. Four are wired here by locating the
    // existing QActions on thermalMenu / viewMenu via addToolbarByLexicon.
    // 99; Save&Quit is owned by the shared 6.0 actQuit_ action which
    // Plan 03 rewires via queueLexicon("99;"); no toolbar entry here.
    // ------------------------------------------------------------------
    tb->addSeparator();
    addToolbarByLexicon(tb, thermalMenu, QStringLiteral("2;"));   // Heat input
    addToolbarByLexicon(tb, thermalMenu, QStringLiteral("3;"));   // Steady heat
    addToolbarByLexicon(tb, thermalMenu, QStringLiteral("4;"));   // Unsteady heat
    if (auto* viewMenu = mbar->findChild<QMenu*>(QStringLiteral("View"))) {
        addToolbarByLexicon(tb, viewMenu, QStringLiteral("8;"));  // Draw temperatures
    }
}
