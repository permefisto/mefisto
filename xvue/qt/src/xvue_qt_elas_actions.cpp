// xvue/qt/src/xvue_qt_elas_actions.cpp -- Phase 6.2 Plan 02
//
// STRONG-symbol body of `registerElasActions_stub_(XvueWindow*, XvueMenuBridge*)`.
// The weak counterpart at xvue/qt/src/xvue_qt_api.cpp:249 is displaced at
// link time when this TU is present in libxvueqt.a -- the dispatch at
// xvue_module_init_ then calls this body whenever `mod == "elas"`.
//
// Composition contract (locked by LEXICON-AUDIT-elas.md Task 2 checkpoint):
//   top-level menus = {File, Solve, View, Help} (Solve replaces Mesh from 6.1)
//   File menu : append 72;, 73;, 74; (99; Save&Quit stays on shared 6.0 actQuit_
//               which 6.2 Plan 03 rewires via queueLexicon)
//   Solve menu : flat leaves 1;..8; at top level plus Parameters submenu
//                (20;, 38;, 39;) and Input/Static/Unsteady/Eigen submenus
//                expanding cl_elas / resoelst / resoelin /
//                methreso / therelas / methvvpr leaves
//   View menu : extend with 7; 8; 10; 71;
//   toolbar   : top-5 = {2;, 3;, 8;, 10;, 99;} (99; via shared 6.0 actQuit_)
//   flat QActions for leaves; menu headers are QMenu (not QAction)
//   every QAction triggered lambda echoes `[menu] <lexicon>` to the console
//   dock BEFORE queuing the chars (mail D-07 pattern)
//
// Security posture (see 06.2-02-PLAN.md <threat_model>):
//   T-6.2-BRIDGE-FLOOD   -- mitigated by XvueMenuBridge::kMaxQueueSize=10000
//   T-6.2-SVG-XXE        -- accepted (Qt 6.5+ disables external-entity SVG;
//                           icons compiled into qrc, no runtime file read)
//   T-6.2-PATH-INJECT    -- mitigated (parser `name` is hard-coded literals)
//   T-6.2-WEAK-OVERRIDE  -- mitigated (keepalive force-link; strong wins)
//   T-6.2-KEEPALIVE-MISS -- mitigated (xvue_qt_elas_actions_keepalive export
//                           referenced from xvue_qt_app.cpp static init)
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
// Mirrors the 6.1 mail helper verbatim -- same semantics, same parent
// ownership, same dynamic-property test hook.
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
// cases can locate the Solve menu via `mbar->findChild<QMenu*>("Solve")`.
// Idempotent -- repeated registration is safe.
//
// Plan 04 Gap-1: when `insertBefore` is non-null, route the new menu
// through QMenuBar::insertMenu(insertBefore->menuAction(), m) so it lands
// BEFORE that anchor menu (used to position Solve between File and View).
// When `insertBefore` is null, fall back to addMenu (append). We construct
// the QMenu(title, mbar) directly because addMenu(QString) creates AND
// appends in one call — there is no addMenu overload that creates without
// appending. The QMenu parent is `mbar` so Qt's parent-child tree handles
// deletion just like the addMenu(QString) variant.
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
// more than one submenu level deep will silently not be found. All five
// lexicons currently requested by this module are within that depth.
// If a future caller targets a three-level-deep action, use a recursive
// helper (see actionByLexicon in the test files for a depth-unlimited model).
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
// an otherwise-undefined reference. The weak `registerElasActions_stub_`
// definition in xvue_qt_api.cpp.o already resolves the dispatch-site
// reference, so without this keepalive the linker NEVER pulls
// xvue_qt_elas_actions.cpp.o from libxvueqt.a -- the strong override
// below would silently be excluded and the warn-once stub ships in
// pp/ppelas_qt. XvueApp's load_bundled_font_() references this symbol
// via a static ref so the TU is always pulled.
extern "C" int xvue_qt_elas_actions_keepalive() {
    return 1;
}

// STRONG C-linkage symbol -- displaces the weak warn-once stub in
// xvue_qt_api.cpp. The dispatch at xvue_module_init_ is module-name-gated,
// so this body ONLY fires when xvue_module_init_ is called with
// name == "elas".
extern "C" void registerElasActions_stub_(XvueWindow* win,
                                          XvueMenuBridge* mb)
{
    if (!win || !mb) return;
    // One-shot guard: registerElasActions_stub_ must only be called once per
    // XvueWindow lifetime. xvue_module_init_ is the sole intended call site.
    // Without this guard a second call would add duplicate separators, QActions,
    // and toolbar buttons because only ensureTopLevelMenu (the Solve-menu block)
    // is itself idempotent; the File/View extensions and toolbar block are not.
    //
    // IN-02 note: this guard is process-scoped (static bool), which is safe
    // under the current single-window-per-process invariant enforced by XvueApp.
    // If Phase 6.3+ ever introduces window recycling, replace with a per-window
    // flag (e.g. win->setProperty("elasRegistered", true)) to avoid the second
    // window silently receiving no elas actions.
    // registerMailActions_stub_ carries the same guard with the same assumption.
    static bool registered = false;
    if (registered) return;
    registered = true;

    // Parse bilingual labels once per sub-menu (cached internally by
    // XvueMenuFileParser -- shared across modules). Parser falls back
    // to "N;" if a file is missing or malformed.
    const MenuFile& debuelas = XvueMenuFileParser::loadFor(QStringLiteral("debuelas"));
    const MenuFile& cl_elas  = XvueMenuFileParser::loadFor(QStringLiteral("cl_elas"));
    const MenuFile& resoelst = XvueMenuFileParser::loadFor(QStringLiteral("resoelst"));
    const MenuFile& resoelin = XvueMenuFileParser::loadFor(QStringLiteral("resoelin"));
    const MenuFile& methreso = XvueMenuFileParser::loadFor(QStringLiteral("methreso"));
    const MenuFile& methvvpr = XvueMenuFileParser::loadFor(QStringLiteral("methvvpr"));
    const MenuFile& therelas = XvueMenuFileParser::loadFor(QStringLiteral("therelas"));
    const MenuFile& contdepl = XvueMenuFileParser::loadFor(QStringLiteral("contdepl"));
    const MenuFile& traccont = XvueMenuFileParser::loadFor(QStringLiteral("traccont"));
    const MenuFile& tracdepl = XvueMenuFileParser::loadFor(QStringLiteral("tracdepl"));
    const MenuFile& valzone  = XvueMenuFileParser::loadFor(QStringLiteral("valzone"));

    auto* mbar = win->menuBar();
    auto* tb   = win->toolBar();

    // ------------------------------------------------------------------
    // File menu -- extended via findChild to avoid creating a duplicate
    // top-level File menu. The 6.0 actQuit_ stays put; Plan 03 rewires
    // its handler to queueLexicon("99;").
    // ------------------------------------------------------------------
    if (auto* fileMenu = mbar->findChild<QMenu*>(QStringLiteral("File"))) {
        fileMenu->addSeparator();
        fileMenu->addAction(makeLeafAction(win, mb,
            debuelas.label(72),
            QStringLiteral("72;"),
            QString(), QStyle::SP_DirIcon));
        fileMenu->addAction(makeLeafAction(win, mb,
            debuelas.label(73),
            QStringLiteral("73;"),
            QString(), QStyle::SP_DirIcon));
        fileMenu->addAction(makeLeafAction(win, mb,
            debuelas.label(74),
            QStringLiteral("74;"),
            QString(), QStyle::SP_TrashIcon));
    }

    // ------------------------------------------------------------------
    // Solve menu -- NEW top-level menu per UI-SPEC "Solve / Calcul"
    // Per-Module Conformance Contract. Uses ensureTopLevelMenu, which is
    // itself idempotent (returns an existing QMenu for the same objectName),
    // but the outer registerElasActions_stub_ is guarded above to run once.
    //
    // Plan 04 Gap-1: look up the View menu so we can insert Solve BEFORE
    // it, yielding the {File, Solve, View, Help} sequence locked by
    // ROADMAP Phase 6.2 goal + 06.0 Per-Module Conformance Contract.
    // 06.0 chrome added File/View/Help in that index order; this line
    // repositions Solve between File and View.
    // ------------------------------------------------------------------
    auto* viewMenuForAnchor = mbar->findChild<QMenu*>(QStringLiteral("View"));
    // WR-03 gap-closure: assert that the View anchor is present so a future
    // chrome-init reordering fails loudly in debug builds rather than
    // silently producing {File, Help, Solve} instead of {File, Solve, View, Help}.
    Q_ASSERT_X(viewMenuForAnchor != nullptr,
               "registerElasActions_stub_",
               "View menu not found — Solve will be appended after Help "
               "(expected {File, Solve, View, Help} order); "
               "chrome init order may have changed.");
    auto* solveMenu = ensureTopLevelMenu(mbar,
        QObject::tr("&Solve"),
        QStringLiteral("Solve"),
        viewMenuForAnchor);

    // Solve > Object name (leaf 1;)
    solveMenu->addAction(makeLeafAction(win, mb,
        debuelas.label(1), QStringLiteral("1;"),
        QStringLiteral(":/xvue/qt/icons/icons/elas/elas-object.svg")));

    // Solve > Input (2; parent) -- cl_elas submenu holds 2;1; 2;2; leaves
    auto* inputMenu = solveMenu->addMenu(debuelas.label(2));
    inputMenu->setObjectName(QStringLiteral("Solve.Input"));
    inputMenu->setIcon(QIcon(QStringLiteral(":/xvue/qt/icons/icons/elas/elas-input.svg")));
    // Parent itself is also a runnable lexicon for users who want "just 2;"
    inputMenu->addAction(makeLeafAction(win, mb,
        QObject::tr("Open input dialog..."), QStringLiteral("2;"),
        QStringLiteral(":/xvue/qt/icons/icons/elas/elas-input.svg")));
    inputMenu->addSeparator();
    inputMenu->addAction(makeLeafAction(win, mb,
        cl_elas.label(1), QStringLiteral("2;1;")));
    inputMenu->addAction(makeLeafAction(win, mb,
        cl_elas.label(2), QStringLiteral("2;2;")));

    // Solve > Static (3; parent) -- resoelst/methreso/therelas leaves
    auto* staticMenu = solveMenu->addMenu(debuelas.label(3));
    staticMenu->setObjectName(QStringLiteral("Solve.Static"));
    staticMenu->setIcon(QIcon(QStringLiteral(":/xvue/qt/icons/icons/elas/solve-static.svg")));
    staticMenu->addAction(makeLeafAction(win, mb,
        QObject::tr("Run static solve..."), QStringLiteral("3;"),
        QStringLiteral(":/xvue/qt/icons/icons/elas/solve-static.svg")));
    staticMenu->addSeparator();
    // 3;1; -- resoelst leaf (Coefficients INDEPENDANTS des DEPLACEMENTS)
    staticMenu->addAction(makeLeafAction(win, mb,
        resoelst.label(1), QStringLiteral("3;1;")));
    // 3;81; -- methreso leaf (synthetic prefix per audit D-04);
    // cited lexicon as 3;81; but the real Fortran dispatch reads
    // methreso.label(1). Audit's synthetic-prefix discipline keeps
    // every QAction's lexicon unique so the action-tree walker can
    // find them by property.
    staticMenu->addAction(makeLeafAction(win, mb,
        methreso.label(1), QStringLiteral("3;81;")));
    // 3;91; -- therelas leaf (synthetic prefix; ELASTICITY ONLY)
    staticMenu->addAction(makeLeafAction(win, mb,
        therelas.label(1), QStringLiteral("3;91;")));

    // Solve > Unsteady (4; parent) -- resoelin leaf 4;1;
    auto* unsteadyMenu = solveMenu->addMenu(debuelas.label(4));
    unsteadyMenu->setObjectName(QStringLiteral("Solve.Unsteady"));
    unsteadyMenu->setIcon(QIcon(QStringLiteral(":/xvue/qt/icons/icons/elas/solve-dynamic.svg")));
    unsteadyMenu->addAction(makeLeafAction(win, mb,
        QObject::tr("Run unsteady solve..."), QStringLiteral("4;"),
        QStringLiteral(":/xvue/qt/icons/icons/elas/solve-dynamic.svg")));
    unsteadyMenu->addSeparator();
    unsteadyMenu->addAction(makeLeafAction(win, mb,
        resoelin.label(1), QStringLiteral("4;1;")));

    // Solve > Eigenmodes (6; parent) -- methvvpr leaf 6;1;
    auto* eigenMenu = solveMenu->addMenu(debuelas.label(6));
    eigenMenu->setObjectName(QStringLiteral("Solve.Eigen"));
    eigenMenu->setIcon(QIcon(QStringLiteral(":/xvue/qt/icons/icons/elas/solve-eigen.svg")));
    eigenMenu->addAction(makeLeafAction(win, mb,
        QObject::tr("Run eigenmodes solve..."), QStringLiteral("6;"),
        QStringLiteral(":/xvue/qt/icons/icons/elas/solve-eigen.svg")));
    eigenMenu->addSeparator();
    eigenMenu->addAction(makeLeafAction(win, mb,
        methvvpr.label(1), QStringLiteral("6;1;")));

    // Solve > Parameters submenu (audit's codes 20; 38; 39;)
    auto* paramsMenu = solveMenu->addMenu(QObject::tr("&Parameters"));
    paramsMenu->setObjectName(QStringLiteral("Solve.Parameters"));
    paramsMenu->addAction(makeLeafAction(win, mb,
        debuelas.label(20), QStringLiteral("20;")));
    paramsMenu->addAction(makeLeafAction(win, mb,
        debuelas.label(38), QStringLiteral("38;")));
    paramsMenu->addAction(makeLeafAction(win, mb,
        debuelas.label(39), QStringLiteral("39;")));

    // ------------------------------------------------------------------
    // View menu -- extended via findChild (Pitfall 7).
    // Adds 7;, 8;+contdepl+traccont+tracdepl+valzone, 10;, 71;
    // ------------------------------------------------------------------
    if (auto* viewMenu = mbar->findChild<QMenu*>(QStringLiteral("View"))) {
        viewMenu->addSeparator();
        viewMenu->addAction(makeLeafAction(win, mb,
            debuelas.label(7), QStringLiteral("7;"),
            QStringLiteral(":/xvue/qt/icons/icons/elas/draw-modes.svg")));

        // 8; Draw Disp+Stress parent -- contdepl submenu with leaves
        // 8;0; 8;1; 8;2; 8;3; 8;6; plus traccont/tracdepl/valzone
        // second-level leaves for the qaction=yes rows.
        auto* drawStressMenu = viewMenu->addMenu(debuelas.label(8));
        drawStressMenu->setObjectName(QStringLiteral("View.DrawStress"));
        drawStressMenu->setIcon(QIcon(
            QStringLiteral(":/xvue/qt/icons/icons/elas/draw-stress.svg")));
        drawStressMenu->addAction(makeLeafAction(win, mb,
            QObject::tr("Draw displacements && stresses..."),
            QStringLiteral("8;"),
            QStringLiteral(":/xvue/qt/icons/icons/elas/draw-stress.svg")));
        drawStressMenu->addSeparator();
        // 8;0; case number prompt (contdepl entry 0)
        drawStressMenu->addAction(makeLeafAction(win, mb,
            contdepl.label(0), QStringLiteral("8;0;")));

        // Stresses submenu (8;1; + traccont leaves 8;1;1; 8;1;15;)
        auto* stressesMenu = drawStressMenu->addMenu(contdepl.label(1));
        stressesMenu->setObjectName(QStringLiteral("View.DrawStress.Stresses"));
        stressesMenu->addAction(makeLeafAction(win, mb,
            QObject::tr("Open stresses dialog..."),
            QStringLiteral("8;1;")));
        stressesMenu->addSeparator();
        stressesMenu->addAction(makeLeafAction(win, mb,
            traccont.label(1), QStringLiteral("8;1;1;")));
        stressesMenu->addAction(makeLeafAction(win, mb,
            traccont.label(15), QStringLiteral("8;1;15;")));

        // Displacements submenu (8;2; + tracdepl leaves 8;2;1; 8;2;2; 8;2;90;)
        auto* dispsMenu = drawStressMenu->addMenu(contdepl.label(2));
        dispsMenu->setObjectName(QStringLiteral("View.DrawStress.Displacements"));
        dispsMenu->addAction(makeLeafAction(win, mb,
            QObject::tr("Open displacements dialog..."),
            QStringLiteral("8;2;")));
        dispsMenu->addSeparator();
        dispsMenu->addAction(makeLeafAction(win, mb,
            tracdepl.label(1), QStringLiteral("8;2;1;")));
        dispsMenu->addAction(makeLeafAction(win, mb,
            tracdepl.label(2), QStringLiteral("8;2;2;")));
        dispsMenu->addAction(makeLeafAction(win, mb,
            tracdepl.label(90), QStringLiteral("8;2;90;")));

        // Von Mises submenu (8;3; + valzone leaf 8;3;1;)
        auto* vmMenu = drawStressMenu->addMenu(contdepl.label(3));
        vmMenu->setObjectName(QStringLiteral("View.DrawStress.VonMises"));
        vmMenu->addAction(makeLeafAction(win, mb,
            QObject::tr("Open von Mises dialog..."),
            QStringLiteral("8;3;")));
        vmMenu->addSeparator();
        vmMenu->addAction(makeLeafAction(win, mb,
            valzone.label(1), QStringLiteral("8;3;1;")));

        // 8;6; eigenmode movement (contdepl leaf 6)
        drawStressMenu->addAction(makeLeafAction(win, mb,
            contdepl.label(6), QStringLiteral("8;6;")));

        // 10; Draw mesh (reuses 6.1 mail mesh-draw.svg -- shared qrc prefix)
        viewMenu->addAction(makeLeafAction(win, mb,
            debuelas.label(10), QStringLiteral("10;"),
            QStringLiteral(":/xvue/qt/icons/icons/mail/mesh-draw.svg")));

        // 71; TMS tools
        viewMenu->addAction(makeLeafAction(win, mb,
            debuelas.label(71), QStringLiteral("71;"),
            QString(), QStyle::SP_BrowserReload));
    }

    // ------------------------------------------------------------------
    // Toolbar top-5 per LEXICON-AUDIT-elas.md toolbar=yes rows:
    // {2;, 3;, 8;, 10;, 99;}. Four are wired here by locating the
    // existing QActions on solveMenu / viewMenu via addToolbarByLexicon.
    // 99; Save&Quit is owned by the shared 6.0 actQuit_ action which
    // Plan 03 rewires via queueLexicon("99;"); no toolbar entry here.
    // ------------------------------------------------------------------
    tb->addSeparator();
    addToolbarByLexicon(tb, solveMenu, QStringLiteral("2;"));   // Input data
    addToolbarByLexicon(tb, solveMenu, QStringLiteral("3;"));   // Static solve
    if (auto* viewMenu = mbar->findChild<QMenu*>(QStringLiteral("View"))) {
        addToolbarByLexicon(tb, viewMenu, QStringLiteral("8;"));   // Draw stress
        addToolbarByLexicon(tb, viewMenu, QStringLiteral("10;"));  // Draw mesh
    }
}
