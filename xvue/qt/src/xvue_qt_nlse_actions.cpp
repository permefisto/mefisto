// xvue/qt/src/xvue_qt_nlse_actions.cpp -- Phase 6.5 Plan 02
//
// STRONG-symbol body of `registerNlseActions_stub_(XvueWindow*, XvueMenuBridge*)`.
// The weak counterpart at xvue/qt/src/xvue_qt_api.cpp:270 is displaced at
// link time when this TU is present in libxvueqt.a -- the dispatch at
// xvue_module_init_ then calls this body whenever `mod == "nlse"`.
//
// Composition contract (locked by LEXICON-AUDIT-nlse.md Plan 01 freeze):
//   top-level menus = {File, Nonlinear, View, Help} (Nonlinear replaces
//                      Thermal from 6.4 / Fluid from 6.3 / Solve from 6.2 /
//                      Mesh from 6.1)
//   File menu  : append 71; (TMS Tools), 72; (MS File Tools), 73; (I/O
//                  Units), 74; (Delete TMS) -- 99; Save&Quit stays on shared
//                  6.0 actQuit_ (Plan 03 wires queueLexicon)
//   Nonlinear menu : flat 1; (Object), 2; INPUT submenu (cl_nlse leaves +
//                    coefnlse leaves under synthetic 2;91; offset),
//                    5; Solve Implicit flat, 6; Solve i-time flat,
//                    Parameters submenu (20; precision, 38; window pixels,
//                    39; background)
//   View menu : extend with 11; (Draw NLSE Module -- primary visual),
//               12; (Draw Real Part), 13; (Draw Imag Part), 14; (Stop-test
//               diagnostic), 15; (Max|U|-vs-time), 16; (Error-vs-time),
//               19; (Draw mesh -- REUSES 6.1 mail mesh-draw.svg)
//   Help menu : append 97; (MEFISTO VERSION DATE -- NOTE: 97; matches flui;
//               NOT 98; like ther; Plan 03's testHelpNoQueue allowlist
//               is {97;} per Auto-fix Rule 1 lesson from 6.3)
//   toolbar   : top-5 = {2;, 5;, 11;, 19;, 99;} (99; via shared 6.0 actQuit_)
//   flat QActions for leaves; menu headers are QMenu (not QAction)
//   every QAction triggered lambda echoes `[menu] <lexicon>` to the console
//   dock BEFORE queuing the chars (mail D-07 pattern, mirrored by
//   elas/flui/ther)
//
// Security posture (see 06.5-02-PLAN.md <threat_model>):
//   T-6.5-ABI-DRIFT      -- mitigated by verify_abi.sh (count stays 58)
//   T-6.5-FORCE-LINK     -- mitigated (xvue_qt_nlse_actions_keepalive export
//                           referenced from xvue_qt_app.cpp static init,
//                           positioned alphabetically BETWEEN mail and ther
//                           since 'mail' < 'nlse' < 'ther')
//   T-6.5-MENU-ORDER     -- mitigated (Q_ASSERT_X on viewMenuForAnchor +
//                           ensureTopLevelMenu(insertBefore=...))
//   T-6.5-PARALLEL-BUILD -- mitigated by serial cbl_tout_qt + cbl_tout sequence
//   T-6.5-XML-PARSE      -- mitigated by ASCII hyphens in qrc comment headers
//
// Lessons applied from 6.4 REVIEW (WR-01 / IN-04 / IN-05):
//   WR-01: codes 11; 12; 13; all dispatch through TRTHER but with DIFFERENT
//          NTYPDESS values (4=module, 5=real, 6=imaginary per ppnlse.f:
//          392, 399, 406). For nlse the View > Draw NLSE * entries are FLAT
//          leaves only -- the downstream tempgrad/tractem* cascade is owned
//          by 6.4 ther audit (compressed proxy rows in LEXICON-AUDIT-nlse.md)
//          and not re-expanded here. No prefix copy-paste hazard since we do
//          not nest sub-leaves under 11/12/13 in this TU.
//   IN-04: every qaction=yes leaf in LEXICON-AUDIT-nlse.md is wired below.
//          Verified one-to-one against the 20 qaction=yes rows:
//          top-level: 1; 2; 5; 6; 11; 12; 13; 14; 15; 16; 19; (11)
//          cl_nlse:   2;2; 2;3; 2;8;                          (3)
//          coefnlse:  2;91;1; 2;91;2; 2;91;4; 2;91;6; 2;91;8; (5)
//          ABI total: 11+3+5 = 19 qaction=yes nlse leaves (audit reports
//          20 -- 99; is the 20th, owned by shared 6.0 actQuit_).
//   IN-05: ppnlse.f:215-222 dispatch table -- codes 8; 9; 17; 18; fall
//          through to label 30 (loop-back). NONE of these are in the
//          qaction=yes set, so there is no dead-code-without-affordance
//          hazard here. (Checked: audit lists no qaction=yes leaf with
//          a dead dispatch.)
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
// Mirrors the 6.1 mail / 6.2 elas / 6.3 flui / 6.4 ther helper verbatim --
// same semantics, same parent ownership, same dynamic-property test hook.
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
// cases can locate the Nonlinear menu via
// `mbar->findChild<QMenu*>("Nonlinear")`.  Idempotent -- repeated
// registration is safe.
//
// 6.2 Plan 04 generalization: when `insertBefore` is non-null, route the
// new menu through QMenuBar::insertMenu(insertBefore->menuAction(), m) so
// it lands BEFORE that anchor menu (used to position Nonlinear between
// File and View, matching {File, Nonlinear, View, Help}).
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
// NOTE: depth limit is two levels (top + one submenu). All toolbar-targeted
// lexicons here (2;, 5;, 11;, 19;) are within that depth.
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
// an otherwise-undefined reference. The weak `registerNlseActions_stub_`
// definition in xvue_qt_api.cpp.o already resolves the dispatch-site
// reference, so without this keepalive the linker NEVER pulls
// xvue_qt_nlse_actions.cpp.o from libxvueqt.a -- the strong override
// below would silently be excluded and the warn-once stub ships in
// pp/ppnlse_qt. XvueApp's load_bundled_font_() references this symbol
// via a static ref so the TU is always pulled.
extern "C" int xvue_qt_nlse_actions_keepalive() {
    return 1;
}

// STRONG C-linkage symbol -- displaces the weak warn-once stub in
// xvue_qt_api.cpp. The dispatch at xvue_module_init_ is module-name-gated,
// so this body ONLY fires when xvue_module_init_ is called with
// name == "nlse".
extern "C" void registerNlseActions_stub_(XvueWindow* win,
                                          XvueMenuBridge* mb)
{
    if (!win || !mb) return;
    // One-shot guard: registerNlseActions_stub_ must only be called once per
    // XvueWindow lifetime. xvue_module_init_ is the sole intended call site.
    // Without this guard a second call would add duplicate separators, QActions,
    // and toolbar buttons because only ensureTopLevelMenu (the Nonlinear-menu
    // block) is itself idempotent; the File/View extensions and toolbar block
    // are not.
    //
    // IN-02 note: this guard is process-scoped (static bool), which is safe
    // under the current single-window-per-process invariant enforced by XvueApp.
    // If a future phase ever introduces window recycling, replace with a
    // per-window flag (e.g. win->setProperty("nlseRegistered", true)) to avoid
    // the second window silently receiving no nlse actions.
    // registerMailActions_stub_ / registerElasActions_stub_ /
    // registerFluiActions_stub_ / registerTherActions_stub_ carry the same
    // guard with the same assumption.
    static bool registered = false;
    if (registered) return;
    registered = true;

    // Parse bilingual labels once per sub-menu (cached internally by
    // XvueMenuFileParser -- shared across modules). Parser falls back
    // to "N;" if a file is missing or malformed.
    const MenuFile& debunlse = XvueMenuFileParser::loadFor(QStringLiteral("debunlse"));
    const MenuFile& cl_nlse  = XvueMenuFileParser::loadFor(QStringLiteral("cl_nlse"));
    const MenuFile& coefnlse = XvueMenuFileParser::loadFor(QStringLiteral("coefnlse"));

    auto* mbar = win->menuBar();
    auto* tb   = win->toolBar();

    // ------------------------------------------------------------------
    // File menu -- extended via findChild to avoid creating a duplicate
    // top-level File menu. The 6.0 actQuit_ stays put; Plan 03 rewires
    // its handler to queueLexicon("99;").
    //
    // nlse exposes 71;/72;/73;/74; as FOUR separate top-level codes for
    // SUITMS / SUIFMS / SUIVES / TUER (matches 6.1 mail's pattern -- UNLIKE
    // 6.2/6.3/6.4 which use combined 70; entry into managmef).
    // ------------------------------------------------------------------
    if (auto* fileMenu = mbar->findChild<QMenu*>(QStringLiteral("File"))) {
        fileMenu->addSeparator();
        // 71; -- SUIVI TMS (TMS administration)
        fileMenu->addAction(makeLeafAction(win, mb,
            debunlse.label(71),
            QStringLiteral("71;"),
            QString(), QStyle::SP_FileDialogDetailedView));
        // 72; -- SUIVI FILES MS (MS file administration)
        fileMenu->addAction(makeLeafAction(win, mb,
            debunlse.label(72),
            QStringLiteral("72;"),
            QString(), QStyle::SP_DirIcon));
        // 73; -- UNITES I/O (I/O unit administration)
        fileMenu->addAction(makeLeafAction(win, mb,
            debunlse.label(73),
            QStringLiteral("73;"),
            QString(), QStyle::SP_DriveHDIcon));
        // 74; -- TUER PLSVO (destructive op -- TMS deletion)
        fileMenu->addAction(makeLeafAction(win, mb,
            debunlse.label(74),
            QStringLiteral("74;"),
            QString(), QStyle::SP_TrashIcon));
    }

    // ------------------------------------------------------------------
    // Nonlinear menu -- NEW top-level menu per UI-SPEC "Nonlinear / Non Linéaire"
    // Per-Module Conformance Contract (per ROADMAP Phase 6.5 line 220 which
    // overrides the older 06.0-UI-SPEC table at line 245 listing
    // "Solve / Calcul" -- see audit's Menu-taxonomy section).
    //
    // ensureTopLevelMenu is itself idempotent (returns an existing QMenu for
    // the same objectName), but the outer registerNlseActions_stub_ is guarded
    // above to run once.
    //
    // 6.2 Plan 04 Gap-1 (carried forward): look up the View menu so we
    // can insert Nonlinear BEFORE it, yielding the {File, Nonlinear, View,
    // Help} sequence locked by ROADMAP Phase 6.5 + 06.0 Per-Module
    // Conformance Contract. 06.0 chrome added File/View/Help in that index
    // order; this line repositions Nonlinear between File and View.
    // ------------------------------------------------------------------
    auto* viewMenuForAnchor = mbar->findChild<QMenu*>(QStringLiteral("View"));
    // 6.2 WR-03 gap-closure (mirrored): assert that the View anchor is
    // present so a future chrome-init reordering fails loudly in debug
    // builds rather than silently producing {File, Help, Nonlinear} instead
    // of {File, Nonlinear, View, Help}.
    Q_ASSERT_X(viewMenuForAnchor != nullptr,
               "registerNlseActions_stub_",
               "View menu not found - Nonlinear will be appended after Help "
               "(expected {File, Nonlinear, View, Help} order); "
               "chrome init order may have changed.");
    auto* nonlinearMenu = ensureTopLevelMenu(mbar,
        QObject::tr("&Nonlinear"),
        QStringLiteral("Nonlinear"),
        viewMenuForAnchor);

    // Nonlinear > Object name (leaf 1;) -- mandatory first command per
    // ppnlse.f:227-302 (every NLSE simulation begins by selecting the object).
    nonlinearMenu->addAction(makeLeafAction(win, mb,
        debunlse.label(1), QStringLiteral("1;"),
        QStringLiteral(":/xvue/qt/icons/icons/nlse/nlse-object.svg")));

    // Nonlinear > NLSE input data (2; parent) -- DFNLSE prompts coefnlse
    // (via SDDEF2) + cl_nlse (LIMTCL). cl_nlse leaves 2;2; 2;3; 2;8; are
    // qaction=yes per audit. coefnlse leaves use synthetic prefix 2;91;<N>;
    // to avoid collision with cl_nlse leaves (cl_nlse has indices 2,3,8;
    // coefnlse has indices 1,2,4,6,8,9 -- collision on 2 and 8 forced
    // the synthetic 2;91; offset; see D-6.5-01-03).
    auto* inputMenu = nonlinearMenu->addMenu(debunlse.label(2));
    inputMenu->setObjectName(QStringLiteral("Nonlinear.Input"));
    inputMenu->setIcon(QIcon(
        QStringLiteral(":/xvue/qt/icons/icons/nlse/nlse-input.svg")));
    // Parent itself is also a runnable lexicon for users who want "just 2;"
    inputMenu->addAction(makeLeafAction(win, mb,
        QObject::tr("Open NLSE input dialog..."),
        QStringLiteral("2;"),
        QStringLiteral(":/xvue/qt/icons/icons/nlse/nlse-input.svg")));
    inputMenu->addSeparator();
    // cl_nlse leaves -- Boundary Condition complex wave inputs
    // 2;2; -- Fgamma FLUX Complexe IMPOSE (Neumann BC)
    inputMenu->addAction(makeLeafAction(win, mb,
        cl_nlse.label(2), QStringLiteral("2;2;")));
    // 2;3; -- ONDE_D VALEUR Complexe IMPOSEE (Dirichlet BC -- canonical)
    inputMenu->addAction(makeLeafAction(win, mb,
        cl_nlse.label(3), QStringLiteral("2;3;")));
    // 2;8; -- ONDE_0 INITIAL Complexe IMPOSEE (initial BC value)
    inputMenu->addAction(makeLeafAction(win, mb,
        cl_nlse.label(8), QStringLiteral("2;8;")));
    inputMenu->addSeparator();
    // coefnlse leaves -- material parameters & forcing & initial condition
    // Synthetic prefix 2;91;<N>; per D-6.5-01-03; runtime dispatch at
    // ther/dfnlse.f:213 SDDEF2 with 'coefnlse' arg under debunlse code 2;.
    // 2;91;1; -- Rho Mass Density (mandatory)
    inputMenu->addAction(makeLeafAction(win, mb,
        coefnlse.label(1), QStringLiteral("2;91;1;")));
    // 2;91;2; -- Alfa Laplacian Coefficient (mandatory)
    inputMenu->addAction(makeLeafAction(win, mb,
        coefnlse.label(2), QStringLiteral("2;91;2;")));
    // 2;91;4; -- FOmega Complex Force Second Member
    inputMenu->addAction(makeLeafAction(win, mb,
        coefnlse.label(4), QStringLiteral("2;91;4;")));
    // 2;91;6; -- Beta NLSE non-linearity coefficient (mandatory)
    inputMenu->addAction(makeLeafAction(win, mb,
        coefnlse.label(6), QStringLiteral("2;91;6;")));
    // 2;91;8; -- U(t0,X) Initial Complex Wave Value
    inputMenu->addAction(makeLeafAction(win, mb,
        coefnlse.label(8), QStringLiteral("2;91;8;")));

    // Nonlinear > Solve Implicit (5; flat leaf) -- IMPLICIT scheme with
    // global 2nx2n non-symmetric matrix; ppnlse.f:340-343 sets TESTNL=7 +
    // CALL NLSEINS. Canonical default solver run.
    nonlinearMenu->addAction(makeLeafAction(win, mb,
        debunlse.label(5), QStringLiteral("5;"),
        QStringLiteral(":/xvue/qt/icons/icons/nlse/solve-nlse-implicit.svg")));

    // Nonlinear > Solve i-time / Gross-Pitaevskii (6; flat leaf, no icon --
    // specialty workflow); ppnlse.f:365-368 sets TESTNL=9 + CALL NLSEINS
    // (semi-implicit nxn scheme). Currently no icon -- user may add
    // solve-nlse-itime.svg in a polish-pass.
    nonlinearMenu->addAction(makeLeafAction(win, mb,
        debunlse.label(6), QStringLiteral("6;"),
        QString()));

    // Nonlinear > Parameters submenu (audit's codes 20; 38; 39;)
    // freq=low forces qaction=no per Rule 8 -- but audit's design wires
    // these via menu-only path (bypassing LIMTCL) so they are accessible
    // via the Parameters submenu without contributing to the
    // registerNlseActions_stub_ "20 qaction=yes" set. See audit notes:
    // 38; -> XVPXFE direct call; 39; -> couleurs prompt.
    auto* paramsMenu = nonlinearMenu->addMenu(QObject::tr("&Parameters"));
    paramsMenu->setObjectName(QStringLiteral("Nonlinear.Parameters"));
    paramsMenu->addAction(makeLeafAction(win, mb,
        debunlse.label(20), QStringLiteral("20;"),
        QString(), QStyle::SP_DialogApplyButton));
    paramsMenu->addAction(makeLeafAction(win, mb,
        debunlse.label(38), QStringLiteral("38;"),
        QString(), QStyle::SP_ComputerIcon));
    paramsMenu->addAction(makeLeafAction(win, mb,
        debunlse.label(39), QStringLiteral("39;")));

    // ------------------------------------------------------------------
    // View menu -- extended via findChild (Pitfall 7 inherited from 6.1).
    // Adds 11; (Draw NLSE Module -- primary visual), 12; (Draw Real Part),
    // 13; (Draw Imag Part), 14; (Stop-test diagnostic), 15; (Max|U|-vs-time
    // diagnostic), 16; (Error-vs-time diagnostic), 19; (Draw mesh -- REUSE
    // 6.1 mail mesh-draw.svg).
    //
    // WR-01 lesson: codes 11/12/13 share TRTHER dispatch but with different
    // NTYPDESS values (4/5/6 per ppnlse.f:392/399/406). The downstream
    // tempgrad/tractem*/tracflux/tracgrad cascade is owned by the 6.4 ther
    // audit (compressed proxy rows in LEXICON-AUDIT-nlse.md), so we wire
    // these as FLAT leaves only -- no nested sub-menus. This avoids the
    // "8;1; vs 7;1; prefix copy-paste" hazard entirely.
    // ------------------------------------------------------------------
    if (auto* viewMenu = mbar->findChild<QMenu*>(QStringLiteral("View"))) {
        viewMenu->addSeparator();

        // 11; Draw NLSE Module -- primary visual; ppnlse.f:391-393
        // CALL TRTHER(KNOMOB, 4, IERR) -- NTYPDESS=4 selects MODULE.
        viewMenu->addAction(makeLeafAction(win, mb,
            debunlse.label(11), QStringLiteral("11;"),
            QStringLiteral(":/xvue/qt/icons/icons/nlse/draw-nlse-modulus.svg")));

        // 12; Draw NLSE Real Part -- ppnlse.f:398-400 CALL TRTHER(KNOMOB,
        // 5, IERR) -- NTYPDESS=5 selects REAL part. Icon shared with 13;
        // since Re/Im U components are visually paired.
        viewMenu->addAction(makeLeafAction(win, mb,
            debunlse.label(12), QStringLiteral("12;"),
            QStringLiteral(":/xvue/qt/icons/icons/nlse/draw-nlse-component.svg")));

        // 13; Draw NLSE Imag Part -- ppnlse.f:405-407 CALL TRTHER(KNOMOB,
        // 6, IERR) -- NTYPDESS=6 selects IMAGINARY part. Icon shared with 12;.
        viewMenu->addAction(makeLeafAction(win, mb,
            debunlse.label(13), QStringLiteral("13;"),
            QStringLiteral(":/xvue/qt/icons/icons/nlse/draw-nlse-component.svg")));

        // 14; Draw STOP TEST values -- ppnlse.f:411-414 CALL TRNLSETST after
        // LXLXOU. NLSE-primary routine; reaches no further LIMTCL sub-menus
        // (verified via grep "LIMTCL" ther/trnlsetst.f -> 0 hits). No icon
        // (post-iteration convergence diagnostic).
        viewMenu->addAction(makeLeafAction(win, mb,
            debunlse.label(14), QStringLiteral("14;"),
            QString()));

        // 15; Draw Max|U(Node)|(Time) -- ppnlse.f:418-421 CALL TRNLSEMXU
        // after LXLXOU. NLSE-primary routine; post-iteration time-series
        // diagnostic. No icon.
        viewMenu->addAction(makeLeafAction(win, mb,
            debunlse.label(15), QStringLiteral("15;"),
            QString()));

        // 16; Draw ERROR(Time) -- ppnlse.f:425-428 CALL TRNLSERR after
        // LXLXOU. NLSE-primary routine; convergence diagnostic. No icon.
        viewMenu->addAction(makeLeafAction(win, mb,
            debunlse.label(16), QStringLiteral("16;"),
            QString()));

        // 19; Draw PLSVO mesh -- ppnlse.f:433-434 CALL TRMAIL.
        // REUSES 6.1 mail mesh-draw.svg via the shared qrc prefix
        // (resolves icons/mail/ from anywhere). Same reuse pattern as
        // 6.2 elas code 10;, 6.3 flui code 19;, and 6.4 ther code 10;.
        viewMenu->addAction(makeLeafAction(win, mb,
            debunlse.label(19), QStringLiteral("19;"),
            QStringLiteral(":/xvue/qt/icons/icons/mail/mesh-draw.svg")));
    }

    // ------------------------------------------------------------------
    // Help menu -- extended via findChild for the version banner.
    // CRITICAL: nlse uses code 97; (matches flui; NOT 98; like ther).
    // Plan 03's testHelpNoQueue allowlist is {97;} per Auto-fix Rule 1
    // lesson from 6.3 -- per-module Help-allowlist drawn from
    // LEXICON-AUDIT, NOT inherited from the previous-module template.
    // Verified against td/m/debunlse line 20 and td/ma/debunlse line 20.
    // ------------------------------------------------------------------
    if (auto* helpMenu = mbar->findChild<QMenu*>(QStringLiteral("Help"))) {
        helpMenu->addSeparator();
        helpMenu->addAction(makeLeafAction(win, mb,
            debunlse.label(97), QStringLiteral("97;"),
            QString(), QStyle::SP_MessageBoxInformation));
    }

    // ------------------------------------------------------------------
    // Toolbar top-5 per LEXICON-AUDIT-nlse.md toolbar=yes rows:
    // {2;, 5;, 11;, 19;, 99;}. Four are wired here by locating the
    // existing QActions on nonlinearMenu / viewMenu via addToolbarByLexicon.
    // 99; Save&Quit is owned by the shared 6.0 actQuit_ action which
    // Plan 03 rewires via queueLexicon("99;"); no toolbar entry here.
    // ------------------------------------------------------------------
    tb->addSeparator();
    addToolbarByLexicon(tb, nonlinearMenu, QStringLiteral("2;"));   // NLSE Input
    addToolbarByLexicon(tb, nonlinearMenu, QStringLiteral("5;"));   // Solve Implicit
    if (auto* viewMenu = mbar->findChild<QMenu*>(QStringLiteral("View"))) {
        addToolbarByLexicon(tb, viewMenu, QStringLiteral("11;"));   // Draw NLSE Module
        addToolbarByLexicon(tb, viewMenu, QStringLiteral("19;"));   // Draw mesh
    }
}
