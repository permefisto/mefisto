// xvue/qt/src/xvue_qt_mail_actions.cpp — Phase 6.1 Plan 02
//
// STRONG-symbol body of `registerMailActions_stub_(XvueWindow*, XvueMenuBridge*)`.
// The weak counterpart at xvue_qt_api.cpp is displaced at link time when this
// TU is present in libxvueqt.a — the dispatch at xvue_module_init_ then calls
// this body whenever `mod == "mail"`.
//
// Composition contract (locked by CONTEXT decisions):
//   D-01  top-level menus = {File, Mesh, View, Help}
//   D-02  File menu: append 80; Import, 90; Export, 70; Manage TMS
//         (99; Save&Quit stays on the existing 6.0 actQuit_, Plan 03 rewires
//         its handler through queueLexicon("99;"))
//   D-03  Mesh menu: flat leaves 1;..7; at the top level; Parameters submenu
//         for 0;, 11;, 20;, 21;; Points submenu holds the 3 saisi_pt leaves
//   D-04  View menu: append 10; Draw mesh, 19; Min/MAX XYZ, 60; Manage window
//   D-05  flat QActions for leaves only; menu headers are QMenu (not QAction)
//   D-06  toolbar: top-5 QActions — {1;1;, 2;, 3;, 10;, 99;} per LEXICON-AUDIT
//   D-07  every QAction triggered lambda echoes `[menu] <lexicon>` to the
//         console dock BEFORE queuing the chars (Pitfall 6: appendLine goes
//         directly to the QPlainTextEdit, not through stdout redirect)
//   D-12  QAction labels read at runtime from $MEFISTO/td/m/<menu> files via
//         XvueMenuFileParser (the LIMTCL source of truth); fallback label
//         is the lexicon-path string itself
//
// Security posture (see 06.1-02-PLAN.md <threat_model>):
//   T-6.1-BRIDGE-FLOOD — mitigated by XvueMenuBridge::kMaxQueueSize=10000
//   T-6.1-SVG-XXE       — accepted (Qt 6.5+ disables external-entity SVG)
//   T-6.1-PATH-INJECT   — mitigated (parser input is hard-coded literals)
//   T-6.1-WEAK-OVERRIDE — accepted (module-gated dispatch — Pitfall 8)
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
// Parent = `win` so Qt's parent-child tree handles deletion. The lexicon
// string is stashed on the QAction via a dynamic property so Plan 03's
// QTest cases can locate actions by payload (RESEARCH §Example 1).
//
// `iconResPath` — Qt resource path (":/xvue/qt/icons/icons/mail/...svg")
//                 OR empty to skip custom icons.
// `sp`          — QStyle::StandardPixmap fallback when `iconResPath` is
//                 empty; pass QStyle::SP_CustomBase to skip icons entirely.
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

    // Plan 03 test hook — `actionByLexicon(menu, lex)` walks the QAction
    // tree filtering on this dynamic property.
    a->setProperty("lexicon", lexicon);

    QObject::connect(a, &QAction::triggered, [win, mb, lexicon]() {
        // D-07: echo BEFORE queue. appendLine writes directly to the
        // QPlainTextEdit (RESEARCH Pitfall 6 / Assumption A5), so we
        // don't recurse through the stdout redirect.
        if (auto* dock = win->consoleDock()) {
            dock->appendLine(QStringLiteral("[menu] ") + lexicon);
        }
        if (mb) mb->queueLexicon(lexicon);
    });
    return a;
}

// Helper for top-level menu creation with a stable objectName so Plan 03
// QTest cases can locate the Mesh menu deterministically via
// `mbar->findChild<QMenu*>("Mesh")`. Idempotent — if a menu with the
// objectName already exists, return it (Pitfall 7 avoid duplicate menus).
//
// 6.2 Plan 04 Gap-1 (mail co-fix): when `insertBefore` is non-null, route
// the new menu through QMenuBar::insertMenu(insertBefore->menuAction(), m)
// so it lands BEFORE that anchor menu (used to position Mesh between File
// and View). When `insertBefore` is null, fall back to addMenu (append).
// We construct the QMenu(title, mbar) directly because addMenu(QString)
// creates AND appends in one call — there is no addMenu overload that
// creates without appending. The QMenu parent is `mbar` so Qt's parent-
// child tree handles deletion just like the addMenu(QString) variant.
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

// Scan `menu` (and one level of its submenus) for a QAction whose
// `lexicon` dynamic property matches `lex`, and append it to `tb`.
// We append the EXISTING QAction — never clone — so the menu and the
// toolbar share one QAction instance (RESEARCH §Don't Hand-Roll).
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

// Force-link keepalive (Phase 6.1 Plan 03 Rule 3 auto-fix). GNU ld's
// archive-search semantics only pull an archive member when it resolves
// an otherwise-undefined reference. The weak `registerMailActions_stub_`
// definition in xvue_qt_api.cpp.o already resolves the dispatch-site
// reference, so without this keepalive the linker NEVER pulls
// xvue_qt_mail_actions.cpp.o from libxvueqt.a — the strong override
// below would silently be excluded and the mesher would run with the
// empty warn-once stub. XvueApp::ensure() references this symbol via
// `xvue_qt_mail_actions_keepalive()` so the TU is always pulled.
extern "C" int xvue_qt_mail_actions_keepalive() {
    return 1;
}

// STRONG C-linkage symbol — displaces the weak warn-once stub in
// xvue_qt_api.cpp per RESEARCH §Pattern 1. The dispatch at
// xvue_module_init_ is module-name-gated (Pitfall 8 verdict), so this
// body ONLY fires when xvue_module_init_ is called with name == "mail".
extern "C" void registerMailActions_stub_(XvueWindow* win,
                                          XvueMenuBridge* mb)
{
    if (!win || !mb) return;

    // D-12: parse bilingual labels once per sub-menu (cached internally).
    // Explicitly load the sub-menus referenced by the QActions below.
    // The parser falls back to "N;" if a file is missing or malformed
    // (handles the td/m/modifm2d `80;` typo — RESEARCH Pitfall 5).
    const MenuFile& debut    = XvueMenuFileParser::loadFor(QStringLiteral("debut"));
    const MenuFile& saisipt  = XvueMenuFileParser::loadFor(QStringLiteral("saisi_pt"));
    const MenuFile& importer = XvueMenuFileParser::loadFor(QStringLiteral("importer"));
    const MenuFile& exporter_= XvueMenuFileParser::loadFor(QStringLiteral("exporter"));
    const MenuFile& managmef = XvueMenuFileParser::loadFor(QStringLiteral("managmef"));
    const MenuFile& managfen = XvueMenuFileParser::loadFor(QStringLiteral("managfen"));

    auto* mbar = win->menuBar();
    auto* tb   = win->toolBar();

    // ------------------------------------------------------------------
    // File menu (D-02) — extended via findChild to AVOID creating a
    // duplicate top-level File menu (RESEARCH Pitfall 7). The existing
    // 6.0 actQuit_ stays put — Plan 03 rewires its handler to
    // queueLexicon("99;").
    // ------------------------------------------------------------------
    if (auto* fileMenu = mbar->findChild<QMenu*>(QStringLiteral("File"))) {
        fileMenu->addSeparator();
        fileMenu->addAction(makeLeafAction(win, mb,
            importer.title().isEmpty()
                ? QStringLiteral("Import PLSVO mesh")
                : importer.title(),
            QStringLiteral("80;"),
            QString(), QStyle::SP_DialogOpenButton));
        fileMenu->addAction(makeLeafAction(win, mb,
            exporter_.title().isEmpty()
                ? QStringLiteral("Export PLSVO mesh")
                : exporter_.title(),
            QStringLiteral("90;"),
            QString(), QStyle::SP_DialogSaveButton));
        fileMenu->addAction(makeLeafAction(win, mb,
            managmef.title().isEmpty()
                ? QStringLiteral("Manage TMS file units")
                : managmef.title(),
            QStringLiteral("70;"),
            QString(), QStyle::SP_DirIcon));
    }

    // ------------------------------------------------------------------
    // Mesh menu (D-03) — NEW top-level menu, added via ensureTopLevelMenu
    // so the same register function is safe to call more than once per
    // window (defensive). The top-level geometry entries are flat
    // QActions except for Points, which has a 3-leaf submenu matching
    // the saisi_pt sub-menu in td/m/.
    //
    // 6.2 Plan 04 Gap-1 (mail co-fix): look up the View menu so we can
    // insert Mesh BEFORE it, yielding the {File, Mesh, View, Help}
    // sequence locked by ROADMAP Phase 6.1 goal + 06.0 Per-Module
    // Conformance Contract. 06.0 chrome added File/View/Help in that
    // index order; this line repositions Mesh between File and View.
    // ------------------------------------------------------------------
    auto* viewMenuForAnchor = mbar->findChild<QMenu*>(QStringLiteral("View"));
    auto* meshMenu = ensureTopLevelMenu(mbar,
        debut.title().isEmpty()
            ? QStringLiteral("&Mesh")
            : debut.title(),
        QStringLiteral("Mesh"),
        viewMenuForAnchor);

    // Mesh > Points submenu (saisi_pt — 3 leaves for point-input methods).
    auto* pointsMenu = meshMenu->addMenu(debut.label(1));
    pointsMenu->setObjectName(QStringLiteral("Mesh.Points"));
    pointsMenu->addAction(makeLeafAction(win, mb,
        saisipt.label(1), QStringLiteral("1;1;"),
        QStringLiteral(":/xvue/qt/icons/icons/mail/mesh-point-add.svg")));
    pointsMenu->addAction(makeLeafAction(win, mb,
        saisipt.label(2), QStringLiteral("1;2;"),
        QStringLiteral(":/xvue/qt/icons/icons/mail/mesh-point-add.svg")));
    pointsMenu->addAction(makeLeafAction(win, mb,
        saisipt.label(3), QStringLiteral("1;3;"),
        QStringLiteral(":/xvue/qt/icons/icons/mail/mesh-point-add.svg")));

    // Mesh > top-level leaves 2;..7; (D-03 "each top-level geometry
    // entry is flat"). Users type the leaf option after clicking these.
    meshMenu->addAction(makeLeafAction(win, mb,
        debut.label(2), QStringLiteral("2;"),
        QStringLiteral(":/xvue/qt/icons/icons/mail/mesh-line-add.svg")));
    meshMenu->addAction(makeLeafAction(win, mb,
        debut.label(3), QStringLiteral("3;"),
        QStringLiteral(":/xvue/qt/icons/icons/mail/mesh-surface-add.svg")));
    meshMenu->addAction(makeLeafAction(win, mb,
        debut.label(4), QStringLiteral("4;"),
        QStringLiteral(":/xvue/qt/icons/icons/mail/mesh-volume-add.svg")));
    meshMenu->addAction(makeLeafAction(win, mb,
        debut.label(5), QStringLiteral("5;"),
        QStringLiteral(":/xvue/qt/icons/icons/mail/mesh-object.svg")));
    meshMenu->addAction(makeLeafAction(win, mb,
        debut.label(6), QStringLiteral("6;"),
        QStringLiteral(":/xvue/qt/icons/icons/mail/mesh-transform.svg")));
    meshMenu->addAction(makeLeafAction(win, mb,
        debut.label(7), QStringLiteral("7;"),
        QStringLiteral(":/xvue/qt/icons/icons/mail/mesh-interpolate.svg")));

    // Mesh > Parameters submenu (D-03) — 0;, 11;, 20;, 21; live here
    // rather than at the top level to avoid cluttering the geometry row.
    auto* paramsMenu = meshMenu->addMenu(QStringLiteral("Parameters"));
    paramsMenu->setObjectName(QStringLiteral("Mesh.Parameters"));
    paramsMenu->addAction(makeLeafAction(win, mb,
        debut.label(0),  QStringLiteral("0;")));
    paramsMenu->addAction(makeLeafAction(win, mb,
        debut.label(11), QStringLiteral("11;")));
    paramsMenu->addAction(makeLeafAction(win, mb,
        debut.label(20), QStringLiteral("20;")));
    paramsMenu->addAction(makeLeafAction(win, mb,
        debut.label(21), QStringLiteral("21;")));

    // ------------------------------------------------------------------
    // View menu (D-04) — extended via findChild (Pitfall 7).
    // ------------------------------------------------------------------
    if (auto* viewMenu = mbar->findChild<QMenu*>(QStringLiteral("View"))) {
        viewMenu->addSeparator();
        viewMenu->addAction(makeLeafAction(win, mb,
            debut.label(10), QStringLiteral("10;"),
            QStringLiteral(":/xvue/qt/icons/icons/mail/mesh-draw.svg")));
        viewMenu->addAction(makeLeafAction(win, mb,
            debut.label(19), QStringLiteral("19;"),
            QString(), QStyle::SP_ArrowRight));
        viewMenu->addAction(makeLeafAction(win, mb,
            managfen.title().isEmpty()
                ? debut.label(60)
                : managfen.title(),
            QStringLiteral("60;"),
            QString(), QStyle::SP_BrowserReload));
    }

    // ------------------------------------------------------------------
    // Toolbar top-5 (D-06) — locate the QActions we already created on
    // the mesh / view menus and append them after a separator. This
    // shares ONE QAction between menu and toolbar per RESEARCH §Don't
    // Hand-Roll. The 99; Save&Quit is the shared 6.0 actQuit_ action —
    // Plan 03 rewires it through queueLexicon("99;") AND adds it to the
    // toolbar; here we do NOT add it (so this plan's toolbar count = 4;
    // Plan 03 brings it to 5 as locked by LEXICON-AUDIT Top-5).
    //
    // NOTE: the LEXICON-AUDIT "Top-5 Toolbar" lists {1;1;, 2;, 3;, 10;,
    // 99;}. This plan adds only the first four; Plan 03's
    // onFileQuit rewrite is what finally binds `99;` to the existing
    // actQuit_ (toolbar placement is up to Plan 03's handler change).
    // ------------------------------------------------------------------
    tb->addSeparator();
    addToolbarByLexicon(tb, meshMenu, QStringLiteral("1;1;"));
    addToolbarByLexicon(tb, meshMenu, QStringLiteral("2;"));
    addToolbarByLexicon(tb, meshMenu, QStringLiteral("3;"));
    if (auto* viewMenu = mbar->findChild<QMenu*>(QStringLiteral("View"))) {
        addToolbarByLexicon(tb, viewMenu, QStringLiteral("10;"));
    }
}
