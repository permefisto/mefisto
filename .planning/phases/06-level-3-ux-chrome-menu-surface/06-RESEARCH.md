# Phase 6: Level-3 UX chrome & menu surface — Research

**Researched:** 2026-04-18
**Domain:** Qt 6.10 Widgets — QMainWindow chrome, menu/toolbar surfaces, QAction lexicon-bridge, QDockWidget console, QSettings persistence, bilingual i18n, dark-mode, canvas gestures
**Confidence:** HIGH (compile-time claims), HIGH (Qt API facts — verified vs. Qt 6.10 docs), MEDIUM (Fortran SAISIE char-at-a-time claim — verified by reading source), MEDIUM (per-module conformance enforcement — design proposal, no runtime precedent yet)

---

## Summary

Phase 6 has two distinct planning surfaces: **6.0** ships a shared `QMainWindow` shell (`XvueMenuBridge`, `XvueConsoleDock`, `XvueAboutDialog`, `XvueRecentProjectsMenu`, `XvuePreferencesDialog`, i18n table, QSettings namespace, dark-mode wiring, canvas gestures) into `libxvueqt.a`; **6.1..6.5** each audit one Fortran module's lexicon and register QActions from inside the shell with no changes to the Fortran solver tree. The critical wiring (UX-03) — `QAction::triggered` → queue → synthetic `notypeevent=2` keyboard event — works by returning **one ASCII character per `xvsouris_` call**, because Fortran's `SACLAV` loop (`xvue/saclav.f:1-347`) reads exactly one char per `xvsouris2_` / `SAIPTC → xvsouris` invocation and accumulates them into `KLG(LHKLG)` until a `;` or CR terminates the line. This is the load-bearing semantic fact of the whole phase; misreading it leads to over-engineered synthetic-event injection that the planner must avoid. `[VERIFIED: xvue/saclav.f:270-317 + util/saiptc.f:23 + util/searclic.f:68]`

The shell lives in `libxvueqt.a` — not in `prpr/pp*.f` — because `prpr/*.f` is Fortran and cannot construct QActions. Per-module menus are declared from inside C++ helpers (`registerMailActions(XvueMenuBridge*, QMenuBar*, QToolBar*)`) invoked from a per-module `extern "C"` hook (proposed: `xvue_module_init_mail_` etc.) called early from each `prpr/pp*.f` `PROGRAM`. `[ASSUMED]` — planner must confirm this registration point; alternative is auto-detecting module identity from `argv[0]` at `XvueApp::ensure()` and self-registering (see §5 for tradeoff).

**Primary recommendation:**
- 6.0 delivers **one bridge class + one dock class + three dialog classes + one i18n translation unit + QSettings helper + canvas-gesture extensions** (~7 new `.cpp` files in `xvue/qt/src/`), all owned by `XvueWindow`/`XvueApp` in the Phase-5 pattern. `XvueMenuBridge::pendingCommands_` is a `QQueue<char>` (one ASCII code per element, not a `QString`) — the unit of dispatch is the character, not the command string. `QAction::triggered` push extends the queue by `cmd.length() + 1` (the `;` terminator) chars; `XvueEventBridge::waitForEvent()` pre-check pops one char and returns synthetically without entering the nested loop.
- i18n uses the **compile-time bilingual table** ratified in UI-SPEC, NOT Qt Linguist. The file-existence check on `$MEFISTO/td/m/anglais` is already load-bearing for the Fortran side (`util/langue.f:22-32`) and the shell scripts (`bin/cbmail_qt:14-23` and ~30 other `bin/cb*`). Single source of truth.
- Console dock uses **`freopen(stdout, ...)` redirection to a pipe + `QSocketNotifier` on the read-end fd**. Fortran's stdout and the Qt chrome live in the same process — `QProcess` is inapplicable. Line buffering is forced via `setvbuf(stdout, NULL, _IOLBF, 0)` in `XvueApp::ensure()` so Fortran `WRITE` flushes reach the notifier promptly.
- QAction shortcut context is `Qt::WindowShortcut` with **modifier-only matching (D-04)** and defense via a **`QEvent::ShortcutOverride` event filter on the canvas** that accepts the override for bare alphanumeric / `;` / `@` / `Esc` / `Return` / `Tab` — so the canvas's `StrongFocus` input path reaches `XvueEventBridge` verbatim for typed-lexicon keys regardless of any module's menu registrations.
- Per-module menu-bar dispatch and canvas context-menu population happen via a **virtual hook + registration callback** so that 6.0 ships a window that refuses to launch if no module has registered (fail-loud, not silent-divergence).

---

<user_constraints>
## User Constraints (from 06-CONTEXT.md)

### Locked Decisions

**D-01 — Phase split:** Split Phase 6 into **6.0** (shared shell: QMainWindow + `XvueMenuBridge` + universal dialogs + `QSettings` persistence + dark-mode chrome) plus **6.1..6.5** (per-module lexicon audit + menu wiring: 6.1 mail, 6.2 elas, 6.3 flui, 6.4 ther, 6.5 nlse). ROADMAP.md must be updated to reflect the split before `/gsd-plan-phase` runs. 6.0 unblocks early A/B of the shell alone; each module phase is an independent release cycle.

**D-02 — Per-module dynamic menu bar:** Each `pp*_qt` executable defines its own top-level menus. 6.0 ships the shared `{File, View, Help}` plus the menu-bar mechanism; each 6.1..6.5 declares and plugs in its module-specific menus (e.g., `ppmail_qt` = `{File, Edit, Mesh, View, Help}`, `ppelas_qt` = `{File, Solve, View, Help}`, etc.). A unified static bar with grayed-out items would surface irrelevant menus.

**D-03 — Audit exhaustive, QActions frequency-weighted:** For each of 6.1..6.5, `LEXICON-AUDIT.md` catalogs EVERY typed lexicon command (full documentation — the long tail is preserved on paper). But `QAction` wiring covers only the **80/20 subset** (approximately 20–40 most-used commands per module). Long-tail commands keep working via the existing typed lexicon. Frequency weighting is judged from `testa/` cases + module drivers.

**D-04 — Modifier rule:** `Ctrl/Alt/Cmd+X` → QAction; plain alphanumeric + digits + `;` + Esc + Return → typed lexicon via `XvueEventBridge`. F-keys (`F1..F12`) route to QAction. Preserves `99;`, `5;90;1;` muscle memory exactly. Esc (27) / `@` (64) continue to follow Phase 5 D-06 abort path.

**D-05 — Dark-mode = System-follow via `QPalette`:** Chrome follows the system dark/light palette automatically. Canvas pixels, axes labels, and scientific colormaps are untouched. User can force light/dark via Preferences dialog (stored in `QSettings` as `ui/color-scheme`).

**D-06 — Recent-projects list:** 10 most-recent projects, shown as File → Recent Projects submenu. "Clear Recent" action at the bottom. Each entry stores absolute `$MEFISTOX/<project>/` path. Missing entries pruned silently.

**D-07 — Console dock behavior:** Visible by default on first launch (state persisted via `QSettings` after that). Auto-scroll to bottom on new output. Copy-to-clipboard action in the dock's context menu. Session-local — log resets each process start. Lines matching `*** ERREUR` parsed and surfaced as `QMessageBox::warning`.

**D-08 — Modal dialog re-entrancy guard:** When `XvueApp::blockingDepth() > 0`, any `QDialog::exec` path refuses silently and shows a 3-second status-bar message: `"Finish current operation first (type 99;)"` (FR: `"Terminez l'opération en cours (tapez 99;)"`). The QAction itself is NOT disabled — the guard lives inside the QDialog wrapper.

**D-09 — About dialog:** Credits **Alain Perronnet (LJLL / UPMC Paris)**. MEFISTO version derived from `incl/homdir.inc` mtime or a new `incl/version.inc` (planner decides). Qt version via `qVersion()`. Build date via `__DATE__` at compile time. Both FR and EN via `$MEFISTO/td/m/anglais` flag.

**D-10 — Toolbar icon source:** Qt built-in `QStyle::StandardPixmap` where semantically appropriate. Custom SVG under `xvue/qt/resources/icons/` where no Qt built-in fits. Icon resource registered via Qt resource system (`.qrc`).

### Claude's Discretion

- `XvueMenuBridge` ownership (likely `XvueApp::menuBridge()` singleton, but planner may put it on `XvueWindow`).
- Whether preferences dialog is one tab or multiple.
- How the About dialog renders multi-line credits.
- Exact `LEXICON-AUDIT.md` markdown schema per module (planner proposes shared template in 6.0 — **UI-SPEC already provides a schema**).
- Whether the `QProcess` stdout pipe-reader lives in `XvueApp` (single global) or per-module.
- Whether 6.1..6.5 execute sequentially or in parallel.

### Deferred Ideas (OUT OF SCOPE)

- **V1X-07** (progress-bar fed from solver stdout) — Phase 9+
- **V1X-08** (animation scrubber) — Phase 9+
- **V1X-09** (single-level snapshot undo) — Phase 9+
- Canvas colormap dialog — out of scope; `PALCDE` stays typed-lexicon per D-05.
- Plugin / extension system — out of scope.
- Multi-document / multi-canvas — out of scope.
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| UX-01 | `QMainWindow` with `QMenuBar` (File/Edit/View/Mesh/Solve/Help), `QToolBar` sharing QActions, `QStatusBar` with mouse coords, shortcuts, tooltips | §1 XvueMenuBridge architecture + §5 QAction registration + §7 Canvas interactions |
| UX-02 | `XvueMenuBridge::pendingCommands_` queue; `QAction::triggered` pushes lexicon strings without calling Fortran or opening modal dialogs directly | §1 + §6 Synthetic keyboard event dispatch |
| UX-03 | Next `xvsouris_`/`xvsouris2_` drains queue and returns synthetic `notypeevent=2` keyboard event; no solver-driver changes | §6 (char-at-a-time semantic is load-bearing) |
| UX-04 | Modal dialogs refuse with status-bar message when `blockingDepth() > 0` | §1 + D-08 + UI-SPEC modal-refuse pattern |
| UX-05 | Per-module lexicon audit of `mail/elas/flui/ther/nlse/`; `.planning/phase-6/LEXICON-AUDIT.md` per module | §9 LEXICON-AUDIT.md schema (6.0 ships template; 6.1..6.5 populate) |
| UX-06 | File ops via `QFileDialog` with project-dir filter (`$MEFISTOX/<project>/`); recent-projects persisted via `QSettings` | §4 QSettings namespace design |
| UX-07 | Window geometry, dock layout, recent-projects, preferences persist via `QSettings` | §4 |
| UX-08 | Bilingual flag `$MEFISTO/td/m/anglais` selects UI language; all labels/tooltips/dialogs/About in FR+EN | §2 Qt i18n strategy (compile-time table) |
| UX-09 | About dialog credits Alain Perronnet; Help menu launches PDF doc via `QDesktopServices::openUrl` | §11 build/CMake + About implementation |
| UX-10 | Console dock displays stdout via pipe-reader | §3 Console dock pipe implementation |
| UX-11 | `*** ERREUR` lines parsed and surfaced as `QMessageBox` | §3 (with D-08 guard interaction) |
| UX-12 | Canvas: wheel-zoom, middle-drag pan, right-click context menu, live coord readout | §7 Canvas interactions |
| UX-13 | Chrome follows system dark-mode via `QPalette`; scientific colormaps unaffected (TEXT-06) | §8 Dark-mode implementation |
</phase_requirements>

## Project Constraints (from CLAUDE.md)

The following project-level directives MUST be honored by all Phase 6 plans:

- **Fortran side must not notice the migration** — `mail/`, `elas/`, `flui/`, `ther/`, `nlse/`, `reso/`, `util/`, `incl/` MUST NOT be modified. New Fortran-visible entry points, if any, belong in `xvue/qt/src/xvue_qt_api.cpp` and `xvue/qt/include/xvue_qt_api.h` only.
- **Compilation must never break** — every commit leaves `bin/cbl_tout_qt` and `bin/cbl_tout` in a clean-build state (the two SHARE `pp/` per MEMORY.md; never run concurrently).
- **Small `testa/` tests must continue to pass** — at a minimum `testa/pan2d`, `testa/torus` (mesher) and the 5 BUILD-10 baseline cases keep running.
- **Ask before adding Ubuntu/Debian packages** — Phase 6 should need nothing beyond the existing `qt6-base-dev` (6.10.2 installed). If a task proposes a new apt package, it must ask first.
- **Git discipline** — commit after every logical step (bridge lands, dock lands, dialog lands, per-menu registration lands). No force-push. No hook bypass.
- **`99;` exit convention** — the typed-lexicon `99;` path must keep working everywhere; QAction dispatch must not shortcut around it.
- **Fortran 77 norms (`doc/normes.ps`)** — N/A to Phase 6 because this phase adds NO Fortran files. All work is C++/Qt.
- **`libgfortran5` pin (MEMORY.md `feedback_debian_sid_libgfortran`)** — unrelated to Phase 6 scope but remains a build-environment invariant.

## Known Phase-5 Decisions Carried Forward

These are not repeated here but implicitly constrain Phase 6:

- `XvueEventBridge` is a `QObject` parented to `XvueWindow`, installed as event filter on `XvueCanvas` at window construction (Phase 5 D-02). `XvueMenuBridge` MUST mirror the same ownership pattern — see §1.
- `XvueApp::blockingDepth()` accessor is the ground truth for D-08 modal guard. Phase 6 queries it; only `BlockingDepthGuard` mutates it. `[VERIFIED: xvue/qt/src/xvue_qt_app.h:29-37 + xvue_qt_event.h:21-34]`
- Esc (27), `@` (64), Return (13), Tab (9), Backspace (8) have locked ASCII mappings via `XvueEventBridge::translateKey` (Phase 5 D-06/D-07). `[VERIFIED: xvue/qt/src/xvue_qt_event.cpp:145-163]`
- `MEFISTO_XVSOURIS_AUTOEXIT` and `MEFISTO_XVSOURIS_DEBUG` env vars must continue to function in the headless/debug paths. Menu-bridge pre-drain must be checked BEFORE the AUTOEXIT short-circuit so headless scripts can still see queued commands if they set them, OR the AUTOEXIT path is hit first and queues are left untouched — planner decides, default proposal: menu-drain BEFORE AUTOEXIT check (matches "menu events are user-initiated, AUTOEXIT is test harness").
- Canvas has `Qt::StrongFocus` + `setMouseTracking(true)` (Phase 5). D-04 modifier rule leans on this — plain alphanumeric reaches the canvas keyPressEvent before any shortcut lookup when the window context matches.
- `QCoreApplication::setAttribute(Qt::AA_CompressHighFrequencyEvents)` is already set in `XvueApp::ensure()` before `QApplication` construction (Phase 5 D-05). `[VERIFIED: xvue/qt/src/xvue_qt_app.cpp:62]`

---

## Standard Stack

### Core (all from Debian apt `qt6-base-dev` 6.10.2+dfsg-7 — NO new packages)

| Component | Qt Class | Purpose | Why Standard |
|-----------|----------|---------|--------------|
| Menu bar | `QMenuBar` | Top-level menu + submenu surface | `QMainWindow::setMenuBar` sets ownership + Qt native accessibility |
| Menu | `QMenu` | Submenu / context menu | Used by `QMenuBar::addMenu` and `QWidget::contextMenuEvent` |
| Action object | `QAction` | Shareable menu-item + toolbar-button handle | Qt's "one action, many surfaces" idiom — same QAction in menu+toolbar |
| Shortcut | `QAction::setShortcut` + `QShortcut` | Keyboard binding to action | `Qt::WindowShortcut` context is Qt default — scopes correctly |
| Toolbar | `QToolBar` | Icon-button row under menu bar | `QMainWindow::addToolBar` + `QAction::setIcon` auto-population |
| Status bar | `QStatusBar` | Bottom strip: transient + permanent message zones | `showMessage(text, ms)` auto-clears; `addPermanentWidget(QLabel)` for coords |
| Dock widget | `QDockWidget` | Console log panel, floatable/closable | `QMainWindow::addDockWidget(Qt::BottomDockWidgetArea, ...)` |
| Scrollable log | `QPlainTextEdit` | Monospace text area with `setMaximumBlockCount` | Cheaper than `QTextEdit` — Qt's documented recommendation for logs |
| Pipe reader | `QSocketNotifier` | Qt event-loop integration for a file descriptor | Only Qt class that sits on raw fd without blocking — matches `freopen` pipe |
| File chooser | `QFileDialog` | Open/save project directory | `QFileDialog::getExistingDirectory` suits MEFISTO's "project is a directory" model |
| Dialog | `QDialog` | Modal base for Preferences | `exec()` is the modal block (D-08 guards this) |
| Message box | `QMessageBox` | Error / confirm / About surfaces | `QMessageBox::about()` templates the About dialog for free |
| Persistence | `QSettings` | Window/dock geometry, recent-projects, prefs | Linux: `$HOME/.config/LJLL/mefisto.conf` (freedesktop convention) |
| Icons | `QStyle::StandardPixmap` + `.qrc` SVG | Built-in glyphs + custom SVG for mesh ops | Qt-native dark-mode tinting — `QStyle` hooks into `QPalette` |
| URL opener | `QDesktopServices::openUrl` | Launch PDF docs in system viewer | Cross-DE; handles `file://` URLs on Linux via `xdg-open` |
| Clipboard | `QApplication::clipboard()` | Console dock copy action | Standard Qt API |
| Palette notify | `QStyleHints::colorSchemeChanged` signal | React to system dark-mode toggle | **Qt 6.5+ feature** — Qt 6.10.2 has it. Emitted on Linux when theme changes |

`[VERIFIED: dpkg -l | grep qt6-base-dev → 6.10.2+dfsg-7]`
`[CITED: doc.qt.io/qt-6/qstylehints.html — colorScheme property + colorSchemeChanged signal, intro Qt 6.5]`

### Supporting (in-repo, already delivered by prior phases)

| Component | Source | Purpose |
|-----------|--------|---------|
| `XvueApp` singleton | `xvue/qt/src/xvue_qt_app.{h,cpp}` | Process-lifetime QApplication + `blockingDepth()` + font load |
| `XvueWindow : QMainWindow` | `xvue/qt/src/xvue_qt_window.{h,cpp}` | Owns canvas, state, bridge(s); Phase 6 grows it with menu/toolbar/status/dock |
| `XvueCanvas : QWidget` | `xvue/qt/src/xvue_qt_canvas.{h,cpp}` | `StrongFocus` + `setMouseTracking(true)` + paint backing pixmap |
| `XvueEventBridge` | `xvue/qt/src/xvue_qt_event.{h,cpp}` | `QObject` event filter + nested `QEventLoop` — Phase 6 consumes `waitForEvent` entry point via the menu-queue pre-check |
| `BlockingDepthGuard` | `xvue/qt/src/xvue_qt_event.h:21-34` | RAII depth counter (Phase 6 queries via `XvueApp::blockingDepth()`) |
| `XvueState` backing pixmap | `xvue/qt/src/xvue_qt_state.{h,cpp}` | Canvas pixels — NEVER touched by chrome code (TEXT-06) |
| Qt resource system (font bundled) | `xvue/qt/CMakeLists.txt:32-35` | Existing `qt_add_resources` pattern — reused for icons |

### Alternatives Considered (and rejected)

| Instead of | Could Use | Why Rejected |
|------------|-----------|--------------|
| Compile-time i18n table | Qt Linguist (`.ts` + `lrelease` → `.qm` + `QTranslator`) | Requires toolchain step (`lrelease`) not in shell-script build; adds file artifacts (`.qm`) to install; duplicates the `$MEFISTO/td/m/anglais` source of truth; live-switching out of scope per PROJECT.md. Verdict: unwarranted complexity. |
| `QSocketNotifier` on pipe fd | `QProcess` external process | Fortran runs IN our process — there's no child to spawn; `QProcess` is the wrong tool |
| Dedicated reader thread | QTimer polling | Polling wastes cycles; `QSocketNotifier` is Qt-idiomatic for fd events |
| `XvueMenuBridge` as second `QObject` filter | Extend `XvueEventBridge` with menu-queue logic | Separation of concerns: Phase-5 bridge is pure input-event relay; menu bridge is queue-drain. Two filters compose cleanly; one-filter-with-two-jobs is harder to reason about on Phase 5 regressions |
| Menu bar per-`XvueWindow` | Single shared global menu (macOS style) | macOS out of scope (PROJECT.md); Linux expects menu inside window |
| Auto-detect module from `argv[0]` | Explicit `xvue_module_init_{mail,elas,...}_` hook | Fragile — `argv[0]` under `exec()` may be different path than expected; users run via shell launcher wrappers. Explicit hook is simpler. See §5 |

**Installation (no changes):**
```bash
# Everything needed is already installed — no apt command in Phase 6.
dpkg -l | grep qt6-base-dev  # verifies: 6.10.2+dfsg-7
```

**Version verification:**
```bash
# Qt 6.10.2 from Debian trixie/forky — verified 2026-04-18
apt-cache policy qt6-base-dev  # → Installed: 6.10.2+dfsg-7
```

---

## 1. XvueMenuBridge Architecture (UX-02, UX-03)

### Ownership

**Recommendation:** Own the bridge on **`XvueWindow`** (mirror Phase 5 D-02 `XvueEventBridge` owner), with an `XvueWindow::menuBridge()` accessor AND a convenience `XvueApp::menuBridge()` that forwards to `window_slot()->menuBridge()`.

Rationale:
- Phase 5 D-02 precedent: the event bridge is owned by `XvueWindow` and parented via `new XvueMenuBridge(canvas_, this)` — Qt parent-child destruction covers cleanup. `[VERIFIED: xvue/qt/src/xvue_qt_window.cpp:11-26]`
- `XvueApp::menuBridge()` convenience accessor mirrors the `XvueApp::blockingDepth()` static style (Phase 5 D-03) so call sites inside `extern "C"` entries stay identical.
- On `xvfermer_` → `xvinitgraphique_` reopen, the bridge is rebuilt fresh — desirable behavior (any stale menu-queue entries from the previous session are discarded).

### Queue data structure

**Recommendation:** `QQueue<char> pendingChars_` (one ASCII code per element).

Rationale (this is the critical §6 finding):
- `xvsouris_` returns ONE character per call (`notypeevent=2, nbc=<ascii>`).
- Fortran `SACLAV` (`xvue/saclav.f:270-317`) calls `xvsouris2_` in a loop, reading one char at a time, and appends each char to `KLG(LHKLG)` until `NASCII=13` (Return) or `;` terminates.
- A multi-character command string like `"5;90;"` must therefore be drained as 5 chars: `'5', ';', '9', '0', ';'` across 5 successive `xvsouris_` calls.
- A `QQueue<QString>` would force per-character splitting in `waitForEvent`, adding complexity for zero benefit. A `QQueue<char>` is the natural unit.

**Queue push API (on the QAction handler side):**
```cpp
void XvueMenuBridge::queueLexicon(const QString& cmd) {
    // cmd examples: "5;90;", "99;", "12;"
    // Appends every char including the terminating ';' so Fortran's
    // SACLAV loop sees the exact byte sequence the user would have typed.
    for (QChar ch : cmd) {
        pendingChars_.enqueue(ch.toLatin1());
    }
}
```

**Drain API (called from `XvueEventBridge::waitForEvent` entry):**
```cpp
std::optional<char> XvueMenuBridge::popChar() {
    if (pendingChars_.isEmpty()) return std::nullopt;
    return pendingChars_.dequeue();
}
```

### Integration point

**Recommendation:** Add a **pre-drain check at the top of `XvueEventBridge::waitForEvent()`** — before `loop.exec()` — that asks `XvueMenuBridge::popChar()` and if non-empty synthesizes the `Result` and returns immediately.

Concrete edit site: `xvue/qt/src/xvue_qt_event.cpp:166-230`. Between the `BlockingDepthGuard` construction and the `saved_loop = loop_; ...` save-restore block, add:

```cpp
// Phase 6 UX-03: menu-queue pre-drain. If a QAction handler queued a
// synthetic lexicon character, return it as a notypeevent=2 key event
// without entering the nested loop. Drains one char per waitForEvent
// call so xvsouris_ delivers characters matching Fortran's SACLAV loop
// (one char per call, char accumulated in KLG until ;/Return).
if (auto* menu = XvueApp::menuBridge()) {
    if (auto c = menu->popChar()) {
        Result r;
        r.notypeevent = 2;
        r.nbc = static_cast<unsigned char>(*c);
        r.x = canvas_ ? canvas_->mapFromGlobal(QCursor::pos()).x() : 0;
        r.y = canvas_ ? canvas_->mapFromGlobal(QCursor::pos()).y() : 0;
        return r;  // no BlockingDepthGuard interaction needed —
                   // depth_guard is already on the stack; ~guard will
                   // decrement. But we skipped loop.exec() so the
                   // filter is never armed. This is correct.
    }
}
```

**Critical:** the pre-drain happens AFTER the `BlockingDepthGuard` ctor but BEFORE `loop_ = &loop`. The guard's depth++ is fine because it decrements on return. The menu-bridge is main-thread-only (same as event-bridge per SHELL-07), so no atomics.

### Synthetic event format

Already answered above: **one char per `xvsouris_`, returned as `notypeevent=2, nbc=<ascii>`.** No multi-char or QKeyEvent injection is needed. `SACLAV` does the string accumulation.

### Thread affinity

Main-thread only, identical to `XvueEventBridge`. `QQueue<char>` operations are not thread-safe; this is fine because `QAction::triggered` signals fire on the main thread and `waitForEvent` runs on the main thread. Add a `Q_ASSERT(QThread::currentThread() == qApp->thread())` at the top of both `queueLexicon` and `popChar` per SHELL-07.

### Interaction with `MEFISTO_XVSOURIS_AUTOEXIT`

`[VERIFIED: xvue/qt/src/xvue_qt_api.cpp::xvsouris_ has AUTOEXIT short-circuit per Phase 5 D-10]`

**Recommendation:** Menu-queue pre-drain MUST precede AUTOEXIT. AUTOEXIT is a test-harness escape hatch; a test that pre-populates the queue expects to see those chars before AUTOEXIT synthesizes a dummy return. The modified `xvsouris_` body ordering:

1. `XVUE_QT_ASSERT_MAIN_THREAD`
2. Menu-queue pre-drain (Phase 6 addition) — if non-empty, return drained char
3. `MEFISTO_XVSOURIS_AUTOEXIT` check (Phase 5) — if set, return synthetic
4. `bridge->waitForEvent(Souris)` (Phase 5 normal path)

### Risk / pitfalls

- **Pitfall 1.1 — Queue build-up during blocked reads.** If a user clicks multiple menu items while `xvsouris_` is NOT being polled (e.g. a long solve phase is executing `mail/*` logic not involving `CALL XVSOURIS`), the queue grows. Next time any `SACLAV`-driven read happens, the queue drains in order. **This is the correct behavior** (matches the "type 99;5;90;" muscle memory). No bound on queue size in v1; 10k is a defensive cap. Planner adds cap.
- **Pitfall 1.2 — Stale queue after `xvfermer_`.** On window reopen the bridge is reconstructed; stale menu-queue entries from the previous window are discarded. This is desirable — the user closed one module and opened another, their menu clicks should not carry over.
- **Pitfall 1.3 — Dispatch during `xvsouris2_` accrochage mode.** If the user is mid-accrochage and clicks a menu item, the item's chars get queued but the current `xvsouris2_` call is in `WaitMode::Souris2` and expects mouse events. The menu-queue pre-drain fires FIRST, returning `notypeevent=2` — which `SACLAV` treats as "user typed a key during picking" → triggers the char-handling branch. This is correct: typing `@` during picking aborts, typing a digit starts a new lexicon command. Menu-click → same semantic.
- **Threat T-06-01-01 — Shortcut-press during accrochage.** User is mid-accrochage, presses `Ctrl+O` → `File → Open`. The QAction handler (in the open-file dialog spawn) sees `blockingDepth() > 0`, refuses, shows 3s message. No menu-queue push. Canvas cleanup is already handled by Phase 5 `cleanupAccrochage()` on any outer return. ✓
- **Threat T-06-01-02 — Queue push AFTER depth decrement.** If a user clicks a menu item exactly as the nested loop is about to quit, there is a race: push happens, loop quits, outer caller consumes the real event, and NOW the queue has 5 orphan chars that drain on the next `xvsouris_`. This is ALSO correct behavior — those chars execute on the next Fortran-level read, which is the natural place for them. Not a bug.

---

## 2. Qt i18n Strategy — Compile-Time Table (UX-08)

### Decision (ratified from UI-SPEC)

**Use the compile-time bilingual lookup table. Do NOT adopt Qt Linguist.**

### Rationale

1. **Single source of truth.** The file `$MEFISTO/td/m/anglais` is already load-bearing for:
   - Fortran: `util/langue.f:22-32` reads `INQUIRE( FILE=HOMDIR//'/td/m/anglais', EXIST=... )` → sets `LANGAG` in `/LANGCM/` COMMON. Every `IF( LANGAG .EQ. 0 )` branch in `mail/`, `util/`, `elas/`, ... depends on this. `[VERIFIED: util/langue.f:22-32 + incl/langue.inc]`
   - Shell scripts: at least 30+ `bin/cb*` scripts check `test -f $MEFISTO/td/m/anglais` to emit FR or EN build banner. `[VERIFIED: grep 'anglais' bin/ returned entries in cbmail_qt:14, cbther_qt:14, cbelas_qt, cbflui_qt, cbl*, etc.]`
   - On this machine: file EXISTS → EN mode. `[VERIFIED: test -f /home/drico/git/mefisto/td/m/anglais → present]`

   Adding `.ts`/`.qm` via Qt Linguist would introduce a SECOND language-selection mechanism that disagrees with the first during transitions (e.g. user deletes `anglais` but `.qm` is still set to English). Keep one mechanism.

2. **Fixed, small corpus.** UI-SPEC enumerates ~38 `MsgId` entries for 6.0. Per-module 6.1..6.5 adds ~20-40 more each → total ~150-250 strings. Well under Qt Linguist's break-even point (~500+ strings where `.ts` / `.qm` tooling amortizes).

3. **No toolchain cost.** Qt Linguist requires:
   - `lupdate` at source-change time (extract new strings)
   - `lrelease` at build time (compile `.ts` → `.qm`)
   - `qt_add_translations` in CMake
   - `.qm` files shipped in install
   - `QTranslator` runtime instance

   ALL of this is new infrastructure; none of it serves a user need beyond what the file-flag table already does.

4. **Live-switch OUT OF SCOPE.** PROJECT.md explicitly excludes live language switching. Process-restart after toggling the flag file is sufficient. `[CITED: .planning/PROJECT.md "Out of Scope" table — "Live-switch bilingual UI via QTranslator" → Reason: "not worth the complexity"]`

### Implementation sketch (from UI-SPEC, refined)

**New files:**
- `xvue/qt/src/xvue_qt_i18n.h` — declares `enum class MsgId : int { FileMenuTitle, FileOpen, ..., ButtonOk }` and functions `bool xvueIsEnglish()`, `const char* tr(MsgId)`, `QString xvueT(MsgId)`.
- `xvue/qt/src/xvue_qt_i18n.cpp` — the FR/EN table + `xvueIsEnglish()` that reads `$MEFISTO/td/m/anglais` once, caches.

**Skeleton (~60 lines):**
```cpp
// xvue_qt_i18n.h
#pragma once
#include <QString>

enum class MsgId : int {
    AppName,
    FileMenuTitle, FileOpen, FileOpenTip, FileSave, FileSaveTip,
    FileSaveAs, FileRecentSubmenu, FileRecentClear,
    FileExport, FileQuit, FileQuitTip,
    ViewMenuTitle, ViewToolbar, ViewConsole,
    ViewZoomIn, ViewZoomOut, ViewFit, ViewPreferences,
    HelpMenuTitle, HelpDocumentation, HelpAbout,
    ConsoleTitle, StatusCoordFormat, ModalRefuse,
    ErrorMsgBoxTitle, QuitConfirmTitle, QuitConfirmBody,
    AboutTitle, AboutBody,
    OpenProjectDialogTitle, OpenProjectFilter,
    PreferencesTitle, PrefRecentCountLabel, PrefConsoleDefaultLabel,
    PrefColorSchemeLabel, PrefColorSchemeSystem, PrefColorSchemeLight,
    PrefColorSchemeDark,
    ButtonOk, ButtonCancel, ButtonApply, ButtonClose,
    DestructiveConfirmBodyGeneric,
    _Count_  // sentinel
};

bool    xvueIsEnglish();
const char* tr(MsgId id);
QString xvueT(MsgId id);   // convenience — QString::fromUtf8(tr(id))
```

```cpp
// xvue_qt_i18n.cpp  (excerpt)
#include "xvue_qt_i18n.h"
#include <QString>
#include <QFileInfo>
#include <cstdlib>
#include <array>

namespace {
struct Entry { const char* fr; const char* en; };
// Indexed by static_cast<int>(MsgId).
constexpr std::array<Entry, static_cast<int>(MsgId::_Count_)> kTable = {{
    {"MEFISTO",              "MEFISTO"},                   // AppName
    {"&Fichier",             "&File"},                     // FileMenuTitle
    {"&Ouvrir projet…",      "&Open Project…"},            // FileOpen
    // ... (complete the 38 entries from UI-SPEC §Copywriting)
}};
}  // namespace

bool xvueIsEnglish() {
    static const bool english = []{
        const char* home = std::getenv("MEFISTO");
        if (!home) return false;                          // default FR
        QFileInfo flag(QString::fromLocal8Bit(home) + "/td/m/anglais");
        return flag.exists();
    }();
    return english;
}

const char* tr(MsgId id) {
    const int ix = static_cast<int>(id);
    Q_ASSERT(ix >= 0 && ix < static_cast<int>(MsgId::_Count_));
    const auto& e = kTable[ix];
    return xvueIsEnglish() ? e.en : e.fr;
}

QString xvueT(MsgId id) { return QString::fromUtf8(tr(id)); }
```

The `xvueIsEnglish()` cache is C++17 static-local (thread-safe init), same pattern as `XvueEventBridge::debug_logging_enabled()` at `xvue_qt_event.cpp:38-47`. `[VERIFIED pattern]`

### Pitfalls

- **Pitfall 2.1 — `_Count_` sentinel drift.** Adding a new `MsgId` without a corresponding table entry crashes with `Q_ASSERT`. Mitigation: `constexpr std::array` with fixed size triggers a compile error if the initializer list size mismatches (`kTable.size() != static_cast<int>(MsgId::_Count_)` via static_assert). Planner adds the `static_assert`.
- **Pitfall 2.2 — Mid-process flag toggle.** `xvueIsEnglish()` caches on first call. If a user touches `td/m/anglais` mid-session, the cache is stale. This matches the Fortran side (`LANGAG` is read at `LANGUE` subroutine entry only, typically during `PPINIT` or `PPMAIL` startup). Documented as: "restart MEFISTO to change language" — UI-SPEC already says this is acceptable.
- **Pitfall 2.3 — `$MEFISTO` undefined.** If a user launches `pp/ppmail_qt` without setting `$MEFISTO`, `xvueIsEnglish()` returns `false` (FR default). Matches `bin/cb*` fallback ("if no `$MEFISTO/td/m/anglais` then FR"). ✓
- **Pitfall 2.4 — UTF-8 encoding mojibake.** Accented FR strings (`é`, `à`, `ç`) must be stored as UTF-8 literals in the source file. Qt 6 `QString::fromUtf8` is correct. Source files must be UTF-8 (verify with `file xvue_qt_i18n.cpp` after writing). C++17 does NOT guarantee the source encoding; project must add `target_compile_options(xvueqt PRIVATE -finput-charset=UTF-8)` to `xvue/qt/CMakeLists.txt`.

---

## 3. Console Dock Pipe Implementation (UX-10, UX-11)

### Decision

**Use `freopen(stdout, ...)` redirection to a POSIX pipe; read via `QSocketNotifier` on the read-end fd; `setvbuf(stdout, NULL, _IOLBF, 0)` forces line buffering.**

### Why in-process (not QProcess)

Fortran solvers run IN the same process as the Qt event loop. `CALL MESSAG` or `WRITE(IMPRIM,*)` prints to the process's stdout via libgfortran's internal buffer. `QProcess` spawns a CHILD process — there is no child.

### Architecture sketch

```
Process startup (xvinitgraphique_ or first XvueApp::ensure()):
  +---- pipe(fd[2]) ---+
  | fd[0] read-end     |
  | fd[1] write-end    |
  +--------------------+
         |
         v
  dup2(fd[1], STDOUT_FILENO)    // stdout now writes into the pipe
  close(fd[1])                  // keep only the dup
  setvbuf(stdout, NULL, _IOLBF, 0)  // line-buffered — flush on '\n'
         |
         v
  QSocketNotifier* nt = new QSocketNotifier(fd[0], QSocketNotifier::Read, dock);
  connect(nt, &QSocketNotifier::activated, dock, &XvueConsoleDock::onPipeReadable);

Fortran runs:
  WRITE(IMPRIM,*) 'MAILLAGE OK'   -> libgfortran -> stdout pipe fd[1]
                                   -> buffered until '\n' flush
                                   -> Qt event loop sees fd[0] readable
                                   -> QSocketNotifier fires
                                   -> XvueConsoleDock::onPipeReadable
                                   -> read(fd[0], buf, N)
                                   -> parse newline, appendLine(QString)
                                   -> QPlainTextEdit::appendPlainText
                                   -> auto-scroll to bottom
                                   -> if line matches /^\*\*\* ERREUR/
                                        -> XvueErrorBatcher::enqueue(line)
                                        -> (500ms later) QMessageBox::warning
```

### Thread-safety

Everything on the main thread. `QSocketNotifier` fires in the event loop. No worker thread needed because:
- Writes are line-buffered (stdout) — small flushes.
- Reads are non-blocking once the pipe buffers (Linux pipe buffer ~64KB default; Fortran long-output sessions WILL exercise back-pressure, see Pitfall 3.2).
- `appendPlainText` is fast.

### Line-buffering forcing

`libgfortran` WRITE-statement buffering is complex. The safe incantation:
```c
// In XvueApp::ensure() AFTER pipe redirect:
setvbuf(stdout, nullptr, _IOLBF, 0);
// Flushing on newline is now guaranteed.
// Fortran unit 6 (default IMPRIM in MEFISTO) maps to stdout.
```

`[CITED: POSIX setvbuf _IOLBF semantics — line-buffered means flush at each '\n']`

If this turns out to not be enough (Fortran RUNTIME keeps its own buffer inside libgfortran), a fallback exists: `FLUSH(IMPRIM)` calls in MEFISTO code — but this would mean touching Fortran, which is forbidden. Alternative: an `extern "C" void xvue_qt_flush_stdout_()` that the C++ side calls after each `xvsouris_` return (since the user is about to see a prompt, it's natural to flush any pending output). Planner decides per empirical test.

### `*** ERREUR` parsing + D-08 interaction

Per UI-SPEC §Console dock:
- Regex: `^\*\*\* ERREUR` (FR) or `^\*\*\* ERROR` (EN) — case-sensitive, anchored at line start.
- Batching window: 500 ms. First matching line arms a `QTimer::singleShot(500, this, &XvueConsoleDock::flushErrorBatch)`. Subsequent matches within 500 ms append to `pendingErrorLines_`. On timer fire, a single `QMessageBox::warning` with joined body appears.
- **D-08 interaction:** `QMessageBox::warning` is ALSO a `QDialog::exec()`. But D-08 says the guard "does NOT apply to `QMessageBox::warning` from `*** ERREUR`" (UI-SPEC §Modal guard, explicit exception).

**However** — and this is a subtle threat: a `QMessageBox::warning` spawned FROM A SOLVER MID-READ still opens a nested event loop. This is OK in principle because `XvueEventBridge::waitForEvent` already uses nested loops (Phase 5), and the Qt event loop is re-entrant. But if a `*** ERREUR` fires DURING the deferred-quit timer window in `xvsouris_` motion coalescing, the message box's exec could eat the `loop.quit` event and hang the bridge.

**Mitigation (T-06-03-01):** When `blockingDepth() > 0`, the error batcher **defers the `QMessageBox` call** until `blockingDepth()` returns to 0. Implementation: check depth on timer fire; if > 0, reschedule the timer; else show. Worst case user sees the error box a moment after `xvsouris_` returns. `[ASSUMED]` — this is a design proposal; planner should validate with an explicit test.

### Pitfalls

- **Pitfall 3.1 — stderr is NOT redirected.** `freopen(stdout, ...)` only catches stdout. Fortran `*** ERREUR` messages typically use `WRITE(IMPRIM,*)` not `WRITE(0,*)` so they go to stdout — good. But `XvueApp::load_bundled_font_()` uses `std::fprintf(stderr, ...)` (xvue_qt_app.cpp:43-46) and `MEFISTO_XVSOURIS_DEBUG` writes to stderr. These DON'T go to the console dock — they go to the launcher terminal. **This is correct behavior for diagnostic output**; user-facing solver messages are on stdout.
- **Pitfall 3.2 — Pipe back-pressure.** Linux default pipe buffer is ~64KB. If Fortran emits a long mesh report without the event loop draining fast enough, libgfortran's write() will BLOCK waiting for the reader. This blocks the solver thread, which is the main thread. Main thread blocks → Qt event loop never drains → deadlock. **Mitigation:** set pipe size up via `fcntl(fd[0], F_SETPIPE_SZ, 1048576)` (1MB) at redirect time. `[CITED: Linux pipe(7) man page — F_SETPIPE_SZ since 2.6.35]`. This is a defensive measure; in practice mesh output is episodic, not a fire-hose.
- **Pitfall 3.3 — Log size unbounded.** `QPlainTextEdit::setMaximumBlockCount(10000)` per UI-SPEC caps memory. Tested Qt pattern — documented.
- **Pitfall 3.4 — `setvbuf` after `freopen` order.** `setvbuf` must be called AFTER `dup2` + `freopen`-equivalent, otherwise it applies to the OLD stdout. Planner: order is redirect → setvbuf.
- **Pitfall 3.5 — Process-level stdout is redirected, so the launcher terminal sees nothing.** If a user runs `pp/ppmail_qt` from a terminal expecting to see stdout there, they won't. This is a UX regression. **Mitigation:** `tee` pattern — write a second pipe that echoes to the original tty. Alternative: only redirect if stdout is NOT a tty (`isatty(STDOUT_FILENO)`). UI-SPEC doesn't specify — planner decides. **Recommend:** always redirect (console dock is the source of truth); document that users wanting stdout on terminal set `MEFISTO_CONSOLE_DOCK_PASSTHROUGH=1` or similar. `[ASSUMED]`

---

## 4. QSettings Namespace Design (UX-07)

### Decision (from UI-SPEC § Window Chrome & Layout)

**Single application name (`mefisto`) + organization name (`LJLL`); per-module `beginGroup("mail" | "elas" | "flui" | "ther" | "nlse")` for all per-module keys.**

### Why single-app (not per-executable app names)

Originally considered: each `pp*_qt` calls `QCoreApplication::setApplicationName("mefisto-mail")` etc. → produces 5 separate `mefisto-*.conf` files.

Rejected because:
1. User expectation: MEFISTO is ONE tool with multiple entry points. One conf file is mental clarity.
2. Some settings (UI color scheme) SHOULD be shared across modules — per-app-name scoping would force duplication.
3. Qt's default Linux path `$HOME/.config/LJLL/mefisto.conf` is discoverable; `mefisto-mail.conf + mefisto-elas.conf + ...` is clutter.

`[VERIFIED via WebSearch/Qt docs: QCoreApplication::setOrganizationName + setApplicationName is the standard pairing; beginGroup is Qt's recommended subsection scoping]`

### Key namespace (locked by UI-SPEC)

| Key | Scope | Type | Source |
|-----|-------|------|--------|
| `mail/window/geometry` | per-module | `QByteArray` (`saveGeometry()`) | UI-SPEC §QSettings key namespace |
| `mail/window/state` | per-module | `QByteArray` (`saveState()` — docks + toolbars) | UI-SPEC |
| `mail/window/console-visible` | per-module | `bool` | UI-SPEC (overrides D-07 default after first launch) |
| `mail/recent-projects` | per-module | `QStringList` (10 items max, D-06) | UI-SPEC |
| `ui/color-scheme` | SHARED (no group prefix) | `QString` — "system"/"light"/"dark" | UI-SPEC D-05 |
| `elas/window/geometry` | per-module | as above | — |
| ... | ... | ... | ... |

Shared keys (no group prefix): `ui/color-scheme`, `app/first-launch` (bool sentinel).

### File path (Linux)

`$HOME/.config/LJLL/mefisto.conf` — Qt's default for `setOrganizationName("LJLL") + setApplicationName("mefisto")`. Verified via `QSettings::fileName()` returning this path on first `QSettings` construction. `[CITED: doc.qt.io/qt-6/qsettings.html — Linux path construction]`

### Implementation location

Recommended: a small helper class `XvuePrefs` in `xvue/qt/src/xvue_qt_prefs.{h,cpp}` with static methods:

```cpp
// xvue_qt_prefs.h
class XvuePrefs {
public:
    // Called by each pp*_qt main() BEFORE first QSettings construction.
    // Idempotent. Sets org+app names and module group scope.
    static void initialize(const char* moduleName);  // "mail", "elas", ...

    static QByteArray windowGeometry();
    static void       saveWindowGeometry(const QByteArray&);
    static QByteArray windowState();
    static void       saveWindowState(const QByteArray&);

    static bool consoleVisible(bool fallback = true);  // D-07 default true
    static void saveConsoleVisible(bool);

    static QStringList recentProjects();
    static void        pushRecentProject(const QString& absPath);  // caps at 10
    static void        clearRecentProjects();

    static QString colorScheme();  // "system" | "light" | "dark"
    static void    saveColorScheme(const QString&);

private:
    static QString moduleGroup_;
};
```

### Module-identity injection

The key question: WHERE does 6.0 learn that it's running as `mail` vs. `elas` vs. `flui`?

Three options (see §5 too):

1. **Explicit module hook.** Per-module plans add `CALL XVUE_MODULE_INIT(<name>)` at the top of `prpr/pp*.f` BEFORE any graphics call. The hook is an `extern "C"` in `xvue/qt/src/xvue_qt_api.cpp` that sets a static `moduleName_` + triggers the per-module `registerXxxActions` callback. **Drawback: changes `prpr/pp*.f` Fortran files.** These are NOT solver files (they're the top-level `PROGRAM` units), but still 5 files get touched in 6.1..6.5. Negligible churn; contained to a single `CALL` insertion.

2. **argv[0] auto-detect.** `XvueApp::ensure()` inspects `fake_argv[0]` — but that's `"mefisto"` (hardcoded at `xvue_qt_app.cpp:56`). Would need to change to read real `argv[0]` from `main()`. The executable path is `pp/ppmail_qt` → `moduleName = "mail"`. **Drawback: fragile when invoked via a symlink or under `exec()`.** Launcher scripts like `bin/MAILLER` may `exec pp/ppmail_qt` — argv[0] is whatever the shell passed, usually absolute path.

3. **Build-time `#define`.** Each `pp*_qt` is linked with `-DMEFISTO_MODULE=\"mail\"` etc. via `bin/cb*_qt`. **Drawback: requires `libxvueqt.a` to read a macro that's set by the LINKER's caller, not `libxvueqt.a`'s own CMake.** Alternative: a small per-module C++ file (`xvue_qt_module_mail.cpp`) that exports a `const char*` → linked only into that pp, selected via `bin/cb*_qt`.

**Recommendation:** Option **1** (explicit Fortran hook) for clarity and to keep `libxvueqt.a` module-agnostic. The 5 one-line edits to `prpr/pp*.f` are in-scope for 6.1..6.5 (each module phase) — NOT 6.0 — because 6.0 defines the hook signature and 6.1..6.5 each adds their own `CALL XVUE_MODULE_INIT('mail');` line. Adds 1 new `extern "C"` to the ABI (ABI count goes from 57 → 58).

`[ASSUMED]` — user may prefer Option 2/3; flag for review.

### Pitfalls

- **Pitfall 4.1 — First-launch defaults.** Users with Phase 1..5 builds have no `mefisto.conf`. First `QSettings::value("mail/window/geometry")` returns default-constructed `QByteArray`. UI-SPEC defaults handle this: 1024x768 window, Bottom docked console visible, etc. Confirmed pattern.
- **Pitfall 4.2 — Write permissions.** `$HOME/.config/LJLL/` creation on first write. Qt handles `mkdir -p` automatically per `QSettings` docs. If `$HOME` is read-only (containerized build), writes fail silently. Document but don't attempt recovery.
- **Pitfall 4.3 — Settings migration.** If Phase 6 ships, then Phase 6.1 changes a key name, users' old settings are orphaned. v1 doesn't need a migration system; planner sets key names carefully at 6.0.
- **Pitfall 4.4 — Thread-safety of `QSettings`.** QSettings is thread-safe BUT not cheap to construct/destruct repeatedly. Pattern: construct a process-scoped `QSettings` singleton in `XvuePrefs` on first access, reuse.

---

## 5. QAction Registration Pattern (UX-01, UX-05)

### Where each module's QActions are declared

**Recommendation:** A per-module C++ function `void registerMailActions(XvueWindow* win, XvueMenuBridge* mb)` defined in `xvue/qt/src/xvue_qt_mail_actions.cpp` (6.1) — analogously for `elas`/`flui`/`ther`/`nlse`. The function is exposed as an `extern "C"` symbol so the Fortran-side hook (`XVUE_MODULE_INIT`) can call it.

**Alternatively** — and simpler: the Fortran-side hook passes the module name as a string, and `XVUE_MODULE_INIT` dispatches internally:

```cpp
// xvue/qt/src/xvue_qt_api.cpp  (Phase 6.0 adds this)
extern "C" void xvue_module_init_(const char* name, int* name_len) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    // Expected to be called exactly once per pp*_qt process, before
    // xvinitgraphique_.
    const QString mod = QString::fromLatin1(name, *name_len).trimmed();
    XvuePrefs::initialize(mod.toUtf8().constData());
    auto& win_slot = XvueApp::window_slot();
    if (!win_slot) return;  // error — init_module should precede xvinit
    auto* win = win_slot.get();
    auto* mb  = win->menuBridge();
    if      (mod == "mail")  registerMailActions(win, mb);   // in 6.1
    else if (mod == "elas")  registerElasActions(win, mb);   // in 6.2
    else if (mod == "flui")  registerFluiActions(win, mb);   // in 6.3
    else if (mod == "ther")  registerTherActions(win, mb);   // in 6.4
    else if (mod == "nlse")  registerNlseActions(win, mb);   // in 6.5
    else {
        std::fprintf(stderr, "xvue-qt: unknown module name '%s'\n",
                     mod.toUtf8().constData());
        // Shell still comes up with just {File, View, Help}.
    }
}
```

**Actual flow for 6.1 (mail) delivery:**
1. 6.0 defines `xvue_module_init_()` ABI + stubs for all 5 `register*Actions` symbols as warn-once stubs (so 6.1..6.5 can land independently).
2. 6.1 adds `registerMailActions` body + `CALL XVUE_MODULE_INIT('mail');` to `prpr/ppmail.f`.
3. 6.2..6.5 same for other modules.

### Signal routing (QAction → bridge queue)

**Recommendation:** Lambda capture, passes the lexicon string to the bridge:

```cpp
// In registerMailActions:
auto* actOpen = new QAction(xvueT(MsgId::FileOpen), win);
actOpen->setShortcut(QKeySequence::Open);           // Ctrl+O
actOpen->setIcon(win->style()->standardIcon(QStyle::SP_DialogOpenButton));
actOpen->setStatusTip(xvueT(MsgId::FileOpenTip));
QObject::connect(actOpen, &QAction::triggered, [mb, win]{
    // Open-Project is a DIALOG action — guarded by blockingDepth
    if (XvueApp::blockingDepth() > 0) {
        win->statusBar()->showMessage(xvueT(MsgId::ModalRefuse), 3000);
        return;
    }
    const QString path = QFileDialog::getExistingDirectory(
        win, xvueT(MsgId::OpenProjectDialogTitle),
        qgetenv("MEFISTOX"));
    if (path.isEmpty()) return;
    XvuePrefs::pushRecentProject(path);
    // Queue the lexicon equivalent of "open project <path>"
    mb->queueLexicon(...);
});

// For a pure-lexicon action (no dialog), the lambda is simpler:
auto* actAdd = new QAction(xvueT(MsgId::MeshAddVertex), win);
actAdd->setShortcut(QKeySequence("Ctrl+V"));
QObject::connect(actAdd, &QAction::triggered, [mb]{
    mb->queueLexicon("5;90;");       // Add-vertex lexicon sequence
});
```

### Icon registration

**Recommendation:** Qt resource system (`.qrc`) registered via `AUTORCC` in `xvue/qt/CMakeLists.txt`. New file: `xvue/qt/resources/xvue_icons.qrc` listing `xvue/qt/resources/icons/*.svg`.

Why `.qrc` (not filesystem):
- Fonts are already in a `.qrc` via `qt_add_resources` — consistent pattern. `[VERIFIED: xvue/qt/CMakeLists.txt:32-35]`
- No install-path headache — icons are inside the `libxvueqt.a` static library.
- `QIcon(":/xvue/qt/icons/add-vertex.svg")` is the access path.

The CMake line:
```cmake
# Addition to xvue/qt/CMakeLists.txt (placed after the existing qt_add_resources
# call for fonts):
qt_add_resources(xvueqt xvue_icons
    PREFIX "/xvue/qt"
    FILES
        resources/icons/add-vertex.svg      # 6.1
        resources/icons/add-edge.svg        # 6.1
        # ... others added per 6.1..6.5
)
```

Per UI-SPEC §Iconography HiDPI: SVG only, no PNG fallback. Enforce via CMake glob check — a tiny add-on script like `xvue/qt/cmake/verify_icons_svg_only.sh` that fails if `*.png` or `*.bmp` is found in `resources/icons/`.

### Keyboard shortcuts

**Recommendation per UI-SPEC D-04:** All `QAction::setShortcut` calls use Qt's `Qt::WindowShortcut` context, which is the **default** so no explicit call needed. `[VERIFIED: doc.qt.io/qt-6/qshortcut.html — "The default value is Qt::WindowShortcut"]`

**Critical defensive invariant:** A `QAction` with shortcut `"Ctrl+C"` fires when the user presses Ctrl+C regardless of focus (as long as the window is active). But a plain `"c"` on a QAction would STEAL the `c` from the canvas's keyPressEvent. D-04 forbids this.

To enforce the D-04 modifier rule in code:
- The lint check at plan time: `grep -rn 'setShortcut.*"[a-z0-9]"' xvue/qt/src/*_actions.cpp` should return zero matches (no bare-letter shortcuts).
- CMake build-time check: scan all `.cpp` files and fail if `setShortcut(QKeySequence("<bare-char>"))` is found.

**Additionally** — and this is UX-12-relevant — the canvas's `keyPressEvent` (via `XvueEventBridge::eventFilter`) must NOT receive modifier-only events when they map to a shortcut. Qt's default routing handles this correctly: `QEvent::ShortcutOverride` fires BEFORE `KeyPress` on the focused widget; if the focused widget (canvas) doesn't accept the override, the shortcut fires.

The canvas's `XvueEventBridge` already `return false`s on any event it doesn't handle (at `xvue_qt_event.cpp:254-255` pass-through). But `ShortcutOverride` is never in its switch — correct. So modifier combos reach the shortcut handler, plain keys reach `KeyPress` → translateKey → `notypeevent=2`. D-04 is satisfied without additional code in 6.0.

**HOWEVER** — to protect against a future 6.1..6.5 regression where a plan accidentally declares `QAction::setShortcut(QKeySequence("5"))`, we SHOULD add a `QEvent::ShortcutOverride` event filter on the canvas that accepts override (by calling `event->accept()`) for all ASCII printable chars without modifiers. This would prevent even an erroneous plain-char shortcut from firing. **This is a defense-in-depth measure** and should be included in 6.0.

### Pitfalls

- **Pitfall 5.1 — Shared QAction in multiple menus.** When a QAction appears in both a menu and a toolbar (per UI-SPEC toolbar recommendation), it's ONE QAction object — just added to both containers via `addAction`. Qt handles this correctly; both surfaces enable/disable together. No special code.
- **Pitfall 5.2 — `menuRole()` on macOS.** `QAction::menuRole` auto-relocates items like "About" / "Quit" to macOS application menu. MEFISTO is Linux-only, so irrelevant — but defensive: don't rely on `menuRole()` for routing.
- **Pitfall 5.3 — Shortcut conflicts across modules.** UI-SPEC reserves `Ctrl+{O,S,E,Q,+,-,0,,}`, `F1`, `F9` for 6.0. 6.1..6.5 must audit and avoid these. Lint check at each module's plan.
- **Pitfall 5.4 — Dynamic shortcut changes on language flip.** UI-SPEC locks FR and EN labels to share shortcut — "&Open Project" (Alt+O in EN) vs "&Ouvrir projet" (Alt+O in FR). Ampersand position is in the `MsgId` table; keep positions consistent.

---

## 6. Synthetic Keyboard Event Dispatch (UX-03 — critical wiring)

This section resolves the most load-bearing research question.

### The Fortran mechanism (ground truth)

`[VERIFIED: xvue/saclav.f:270-333 read fully]`

```fortran
      CALL SAIPTCSU( NOCODE, NX, NY, NASCII )  ! → xvsouris_ for INTERA=3
...
      IF( NOCODE .GT. 0 ) THEN
         ! Mouse click branch — not our concern here
      ENDIF
C     LE CARACTERE CARLU A ETE LU AU CLAVIER PHYSIQUE
      CARLU = CHAR( NASCII )
      IF( NASCII .EQ. 13  ) THEN
C        CR RETOUR A LA LIGNE
         GOTO 9900                                 ! line complete → return
      ENDIF
...
      IF( NASCII .EQ. 8 .OR. NASCII .EQ. 127 ) THEN
         ! backspace
         KLG(LHKLG)(NDC:NDC) = ' '
         GOTO 400
      ELSE IF( ... ) THEN
         ! space or other terminators
         ...
      ELSE IF( CARLU .EQ. '@' .OR. CARLU .EQ. '?' ) THEN
         ! abort or help
         KLG(LHKLG)(NDC+1:NDC+1) = CARLU
         KLG(LHKLG)(NDC+2:NDC+2) = ';'
         NSORTI = 1
      ELSE IF( NASCII .GE. 32 .AND. NASCII .LE. 126 ) THEN
C        UN CARACTERE STANDARD
         KLG(LHKLG)(NDC+1:NDC+1) = CARLU           ! append ONE char
         NBLANC = 0
      ELSE
         GOTO 10                                   ! loop back to read next
      ENDIF
...
      GOTO 10                                       ! read next char
```

So `SACLAV` is a **one-char-per-call read-and-accumulate loop**. Each `xvsouris2_` returns ONE ascii code (or a mouse event). `SACLAV` appends the char to `KLG(LHKLG)(NDC+1:NDC+1)` and loops. When the char is `';'` (implicitly in `NASCII .GE. 32 .AND. NASCII .LE. 126`) it gets appended — but the line is NOT yet terminated; the terminator is `CR` (`NASCII=13`, `GOTO 9900`).

Actually — re-read — **`;` is treated as a normal printable char** and appended to `KLG(LHKLG)`. The caller of SACLAV (EVMENU → LIRLIG → DONNMF) then sees the completed line containing `"5;90;"` and does token parsing on it at a higher level (in `DONNMF`, line-scan logic at `util/donnmf.f:118` searching for `INDEX(KLG(...), ';')`).

**But**: `DONNMF` requires a CR to end the line (`READ` statement in `LIRLIG:68`, OR `SACLAV` sets `GOTO 9900` when `NASCII=13`). So a multi-char menu-driven command must include a **final CR (13) AT THE END** to flush it into the `DONNMF` parser.

### What the menu-bridge queue must push

**For command `"5;90;"` the queue should contain: `'5', ';', '9', '0', ';', 13`** — 6 elements, the last being a CR.

Wait — verify this is what the user would type. When you type `5;90;` on the keyboard:
- Each keystroke goes through `xvsouris_` → `xvsouris2_` → `SACLAV` → appends to `KLG`
- At some point the user presses **Return** to execute. That's the `NASCII=13` that triggers `GOTO 9900` → line returns to `LIRLIG` → `DONNMF` parses.

So yes — the QAction handler must queue `"5;90;\r"` (or equivalently push chars `{'5', ';', '9', '0', ';', 13}`). The trailing CR is the "execute" keystroke.

**Simpler helper:**
```cpp
void XvueMenuBridge::queueLexicon(const QString& cmd) {
    for (QChar ch : cmd) pendingChars_.enqueue(ch.toLatin1());
    pendingChars_.enqueue(13);   // trailing CR — always
}
```

The QAction's call site is then trivial: `mb->queueLexicon("5;90;");`.

### Resolves Assumption A4

A4 in the topic list: "does SAISIE block with one xvsouris_ call per char, or does it read a buffer?" — **Verified: one char per call, SACLAV accumulates, CR terminates.** `[VERIFIED via saclav.f, saiptc.f, donnmf.f]`

### Approach (a) vs (b) reconciliation

Both approaches from the topic list are now compared:

| Approach | Feasibility | Verdict |
|----------|-------------|---------|
| (a) Queue returns chars one-at-a-time across successive `xvsouris_` calls | Matches Fortran exactly | **Chosen** |
| (b) `QCoreApplication::postEvent(canvas, new QKeyEvent(...))` and let eventFilter dispatch normally | Technically works for single chars but SACLAV's read-loop is synchronous around xvsouris_ — posting a KeyEvent to the canvas only gets dispatched after control returns to the Qt event loop, which it does inside `loop.exec()`. So either approach reaches SACLAV. BUT approach (b) races with pending Qt events (queued mouse events, timer events). And approach (b) requires posting N chars upfront, each going through loop.exec() round-trip. Approach (a) is simpler, less ordering-dependent. | **Rejected** |

Approach (a) is chosen and is the entire §1 architecture.

### Pitfalls

- **Pitfall 6.1 — `xvsouris2_` mode.** In `WaitMode::Souris2`, the bridge's event filter handles motion/press/release differently (accrochage). If a menu-queue char arrives via pre-drain while in Souris2 mode, the pre-drain returns `notypeevent=2, nbc=<char>` which `SACLAV` branches into the key-handling path (line 99 `ELSE IF( NOCODE .EQ. 2 )`). **But SACLAV's accrochage branch translates `NOCODE=2` to `NOCODE=-1` at line 102**, then proceeds through the char-append logic. This is correct. A typed `@` during accrochage aborts. A typed `5` during accrochage starts a new char — matching the typed-lexicon behavior exactly. ✓
- **Pitfall 6.2 — Menu click while XvueApp has no window.** If `XvueApp::window_slot()` is null (between `xvfermer_` and next `xvinitgraphique_`), a QAction could still fire — but there's no window to click in, so no menu bar. Defensive: `XvueMenuBridge::queueLexicon` no-ops when `window_slot()` is null, OR ensures the bridge is destroyed with the window. The latter is the natural Qt parent-child lifetime — bridge is a child of the window; when window is destroyed the bridge goes with it. The menu bar's QActions are also children of the window → no menu bar, no QActions. ✓
- **Pitfall 6.3 — Batch-mode (`INTERA=0`).** `LIRLIG:52` branches on `INTERA .LE. 1` → reads from a file, NOT via `EVMENU`/`SACLAV`. In batch mode, `xvsouris_` is never called. So menu-bridge queue drain never fires. Correct — menu items don't mean anything in batch mode. No special code.
- **Pitfall 6.4 — `@` as typed-lexicon abort.** Per D-04 the user can type `@` to abort. Per D-06 Phase 5, Esc maps to 27 at bridge level, `@` maps to 64. Menu-bridge queue uses raw ASCII 64 for `'@'`. SACLAV treats both (NASCII=27 becomes NASCII=64 internally at `saiptc.f:60-64`; SACLAV `CARLU .EQ. '@'` at line 304 triggers `NSORTI=1` → `GOTO 9900` → flush). A "Cancel" menu item could queue `"@"` + CR, aborting the current read. Probably not a 6.0 deliverable; flag for 6.1..6.5.

---

## 7. Canvas Interactions (UX-12)

### Wheel zoom

**Recommendation:** Pure Qt-side view transform composed at paint time. Does NOT dispatch to Fortran.

```cpp
// XvueCanvas additions (in 6.0):
void XvueCanvas::wheelEvent(QWheelEvent* event) {
    // 1.15x per notch. angleDelta().y() is in 1/8 degree, 120 = 1 notch.
    const int notches = event->angleDelta().y() / 120;
    if (notches == 0) return;
    const qreal factor = std::pow(1.15, notches);
    state_->view_transform_ = state_->view_transform_ * QTransform().scale(factor, factor);
    // Clamp to prevent extreme zoom: 0.1x..10x
    // ...
    update();  // request repaint
}

void XvueCanvas::paintEvent(QPaintEvent*) {
    QPainter p(this);
    // Apply view transform BEFORE drawing backing pixmap.
    p.setTransform(state_->view_transform_);
    p.drawPixmap(0, 0, *state_->backing_);
}
```

The `view_transform_` is a new field on `XvueState` (a `QTransform`, default identity). Backing pixmap is unchanged — Fortran's drawing coordinates stay canonical.

**Important:** the existing `XvueCanvas::paintEvent` at `xvue_qt_canvas.cpp:45-53` is currently a single `QPainter(this).drawPixmap(0, 0, *state_->backing_);`. Phase 6 adds the `setTransform` before `drawPixmap`. This is a Phase 2 DRAW-01 invariant extension (paintEvent still does ONE drawPixmap, just with a transform) — documented as a Phase 6 allowed extension.

### Middle-drag pan

```cpp
void XvueCanvas::mousePressEvent(QMouseEvent* e) {
    if (e->button() == Qt::MiddleButton) {
        pan_start_ = e->pos();
        pan_anchor_transform_ = state_->view_transform_;
        e->accept();
        return;                         // does NOT dispatch to bridge
    }
    QWidget::mousePressEvent(e);        // fallthrough to bridge via filter
}

void XvueCanvas::mouseMoveEvent(QMouseEvent* e) {
    if (e->buttons() & Qt::MiddleButton) {
        const QPoint delta = e->pos() - pan_start_;
        state_->view_transform_ = pan_anchor_transform_.translated(delta.x(), delta.y());
        update();
        e->accept();
        return;
    }
    QWidget::mouseMoveEvent(e);        // fallthrough
}
```

### Right-click context menu

```cpp
void XvueCanvas::contextMenuEvent(QContextMenuEvent* e) {
    if (XvueApp::blockingDepth() > 0) return;  // suppress per UI-SPEC
    QMenu menu(this);
    // Populated by 6.1..6.5 via callback registered with bridge
    if (auto* mb = XvueApp::menuBridge()) mb->populateContextMenu(&menu);
    if (!menu.isEmpty()) menu.exec(e->globalPos());
    e->accept();
}
```

`XvueMenuBridge` exposes a `populateContextMenu` that calls a per-module callback registered by `registerMailActions` (or equivalent). The callback adds the top-6 most-used QActions from the module's LEXICON-AUDIT.md frequency data. **6.0 ships the hook; 6.1..6.5 populate it.**

### Live coord readout

The existing Phase 5 bridge eats MouseMove in its filter. For a coord readout that fires EVEN when no blocking read is active, we need a separate path.

**Recommendation:** Add a `XvueCanvas::mouseMoveEvent` override that **emits a signal** regardless of the bridge's consumption. Signal connects to `XvueWindow::updateStatusCoords` which updates the status-bar permanent widget label.

```cpp
class XvueCanvas {
    Q_OBJECT
signals:
    void mouseCoords(QPoint posLogicalPx);   // emitted every MouseMove
protected:
    void mouseMoveEvent(QMouseEvent* e) override {
        emit mouseCoords(e->pos());
        QWidget::mouseMoveEvent(e);  // lets bridge filter also see it
    }
};
```

Bridge filter's `eventFilter` runs BEFORE the widget's `mouseMoveEvent` (because filters precede target). So:
1. Filter sees MouseMove → if `loop_`, coalesces (Phase 5 behavior preserved).
2. Either way, the event continues to `mouseMoveEvent` (`return false` in filter for non-blocking case → falls through to widget).

**But wait** — Phase 5's filter at `xvue_qt_event.cpp:253-255` checks `if (!loop_) return false;` and for MouseMove during Souris mode at line 464-476 it `return true;` (eats event). This means when a blocking read IS in progress, the widget's `mouseMoveEvent` is NEVER called → signal NOT emitted → coords NOT updated.

**This is actually a Phase 6 design question: should status-bar coords update during a Fortran blocking read?** Yes — the user is picking a point and wants to see the coordinate. So Phase 6 needs to:
- When filter eats MouseMove (in Souris/Souris2 modes), ALSO emit the coords signal before eating. This is a 2-line addition to `xvue_qt_event.cpp`.

Alternative: emit coords from INSIDE the filter's MouseMove branch, then continue as Phase 5. Planner's call. Either is fine.

### Status-bar integration

```cpp
// xvue_qt_window.cpp additions:
void XvueWindow::buildStatusBar() {
    auto* sb = statusBar();   // auto-creates via QMainWindow::statusBar()
    coordLabel_ = new QLabel(sb);
    sb->addPermanentWidget(coordLabel_);
    connect(canvas_, &XvueCanvas::mouseCoords,
            this,     &XvueWindow::updateStatusCoords);
}

void XvueWindow::updateStatusCoords(QPoint p) {
    coordLabel_->setText(xvueT(MsgId::StatusCoordFormat).arg(p.x()).arg(p.y()));
}
```

### Pitfalls

- **Pitfall 7.1 — View transform not exported.** If a user zooms, then Phase 7 exports the canvas, the export captures the ZOOMED view (because `paintEvent` applies the transform). UI-SPEC already documents this as a Phase 7 concern (forwarding the transform to export).
- **Pitfall 7.2 — Middle-click parity with Fortran.** Phase 5's bridge at `xvue_qt_event.cpp:329` uses MiddleButton = abort. But UX-12 says middle-drag = pan. The conflict: `XvueCanvas::mousePressEvent` with middle button must call `e->accept()` BEFORE the bridge filter fires — **but widget event handler runs AFTER the filter.** So the bridge WILL receive the MiddleButton first and abort the Fortran read. **This breaks the new pan behavior.**

  **Mitigation:** Phase 6 updates `XvueEventBridge::eventFilter` to CHECK the middle button's drag mode first. Proposal: the filter IGNORES MiddleButton press when `mousePressEvent`'s pan-start tracking is armed. Simpler proposal: the filter ONLY handles MiddleButton=abort when `mode_ == WaitMode::Souris2` (accrochage) — in regular `Souris` mode, middle-click passes through (returns `false`) so the canvas's `mousePressEvent` handles pan. Phase 5 currently does abort in BOTH Souris and Souris2 modes — but the X11 reference (`xvuelc.c:2272`) aborts only in certain contexts. **This requires re-reading the xvuelc.c X11 reference to be sure.** `[ASSUMED]` — flag for discuss-phase or explicit user decision.

  Conservative fallback: require the user to hold Ctrl+Middle-click for pan. Avoids ambiguity.

- **Pitfall 7.3 — Right-click in Fortran mode.** Fortran's `SACLAV` branch for mouse buttons at line 83-87 treats button 2 as `@` abort (middle mouse button). Button 3 (right) is returned as `NOCODE=3`. **If the user right-clicks during an `xvsouris_` read, the contextMenuEvent fires IN ADDITION to the filter's MouseButtonPress with button=3.** `contextMenuEvent` is triggered by `QEvent::ContextMenu` which fires AFTER `MouseButtonRelease`. If the filter already ate the press+release, the contextMenuEvent is still delivered to the widget.

  **This means right-click during a blocking read would both (a) return button=3 to SACLAV AND (b) show a context menu.** Messy. D-08-style: suppress contextMenuEvent when `blockingDepth() > 0`, as we already do above. ✓ Verified.
- **Pitfall 7.4 — Fit-to-window reset.** `View → Fit to Window` (Ctrl+0) resets `view_transform_` to identity. New XvueCanvas method: `void resetView() { state_->view_transform_ = QTransform(); update(); }`. Expose via QAction.

---

## 8. Dark-Mode Implementation (UX-13, D-05)

### Detection (system-follow)

**Recommendation:** Use Qt 6.5+'s `QStyleHints::colorScheme` property and `colorSchemeChanged` signal. Qt 6.10.2 has them.

```cpp
// In XvueApp::ensure() after QApplication construction:
auto* sh = QGuiApplication::styleHints();
// sh->colorScheme() returns Qt::ColorScheme::{Light, Dark, Unknown}
// Connect to re-apply palette on system theme flip:
QObject::connect(sh, &QStyleHints::colorSchemeChanged,
                 [](Qt::ColorScheme s){
                     // Phase 6: handled by preferences override — see below
                 });
```

`[CITED: doc.qt.io/qt-6/qstylehints.html — colorScheme property introduced Qt 6.5, colorSchemeChanged signal]`

### User override (Preferences)

Per D-05, user can force light/dark via Preferences:
- `ui/color-scheme` ∈ {`system`, `light`, `dark`} in QSettings.
- "system" → follow `QStyleHints::colorScheme()`.
- "light" / "dark" → call `QApplication::setPalette(...)` with a hand-constructed palette OR rely on `QStyleHints::setColorScheme()` (Qt 6.8+; check Linux support — see Pitfall 8.1).

**Simplest implementation for 6.0:** use `QPalette` directly.
```cpp
void XvueApp::applyColorSchemePreference() {
    const QString pref = XvuePrefs::colorScheme();
    QPalette p;
    if (pref == "dark") {
        p = QApplication::palette();            // start from system
        p.setColor(QPalette::Window,     QColor(53, 53, 53));
        p.setColor(QPalette::WindowText, Qt::white);
        p.setColor(QPalette::Base,       QColor(42, 42, 42));
        p.setColor(QPalette::Highlight,  QColor(38, 79, 120));
        // ... (full dark palette, standard Qt example)
        QApplication::setPalette(p);
    } else if (pref == "light") {
        // Standard light: use Qt::LightGray-ish
        p = QPalette();                         // Qt default
        QApplication::setPalette(p);
    } else {
        // "system" — do NOTHING; style inherits from system.
        QApplication::setPalette(QApplication::style()->standardPalette());
    }
}
```

Called on app startup (after `XvueApp::ensure()`) and whenever the Preferences dialog saves.

### Canvas exclusion (TEXT-06 invariant)

**The canvas must NOT react to the palette change.** Phase 3 already has defensive guards:
- `XvueCanvas::XvueCanvas(...)` at `xvue_qt_canvas.cpp:18-21`: sets `setAutoFillBackground(false)` + `setAttribute(Qt::WA_OpaquePaintEvent, true)`. `[VERIFIED]`
- `paintEvent` draws backing pixmap ignoring palette.

Phase 6 adds NO additional palette-application code to the canvas. ✓

### Pitfalls

- **Pitfall 8.1 — Linux programmatic color-scheme setter is flaky.** `QStyleHints::setColorScheme()` on Linux depends on each DE's theme plugin. KDE/GNOME may ignore it. **Verdict:** don't use `setColorScheme` for preference override; use `QApplication::setPalette` directly (as above). `[CITED: Qt Forum discussions on setColorScheme Linux gaps]`
- **Pitfall 8.2 — `colorSchemeChanged` signal firing on Debian XFce.** A1 in the topic list. XFce uses `xsettings` + GTK theme; Qt reads GTK via `qt6ct` or `qgnomeplatform` plugin. XFce without those plugins may NOT emit the signal. Detection test (manual):
  1. Set XFce theme to dark (`xfsettings` GUI).
  2. Run `pp/ppmail_qt`.
  3. Open Terminal; `gsettings set org.gnome.desktop.interface color-scheme prefer-dark` — see if app chrome flips.
  
  **Fallback:** on startup always read current system theme once via `QStyleHints::colorScheme()` and apply. Missing runtime signal is acceptable (user restarts app to flip theme).  Document in README.
- **Pitfall 8.3 — Custom SVG icons stay light on dark theme.** Per UI-SPEC, custom SVGs use `currentColor` stroke. Qt's `QIcon` rendering interpolates `currentColor` from the widget's palette text color. Verify when 6.1 adds the first SVG. Quick test recipe: add an SVG with `stroke="currentColor"`, put it on a `QToolButton`, flip palette, confirm glyph re-tints.
- **Pitfall 8.4 — Testing via CLI.** On Debian, `gsettings set org.gnome.desktop.interface color-scheme 'prefer-dark'` toggles the GNOME setting; `plasma-apply-colorscheme BreezeDark` toggles Plasma; no standard for XFce. For automated testing, a headless xvfb session can't exercise system theme. So the dark-mode test is MANUAL (UI-SPEC flags this as "manual validation during 6.0 sign-off").

---

## 9. LEXICON-AUDIT.md Schema (UX-05, D-03)

### Schema (shared template — 6.0 ships this)

UI-SPEC §"Per-Module Conformance Contract" already specifies the schema:

```markdown
# <Module> Lexicon Audit

## Menu: <mesh / solve>
### <Group> (e.g. "Generate", "Refine", "Geometry")

| Lexicon input | FR label | EN label | Shortcut | Toolbar? | Icon | Frequency (H/M/L) | Fortran dispatch | Surfaced as QAction? |
|---------------|----------|----------|----------|----------|------|-------------------|------------------|----------------------|
| `5;90;` | Ajouter sommet | Add vertex | `Ctrl+V` | yes | `add-vertex.svg` | H | `mail/aisommet.f::A1LGSOM` | yes |
```

Schema confirmation: single shared template goes into 6.0 as `.planning/phases/06-level-3-ux-chrome-menu-surface/LEXICON-AUDIT-TEMPLATE.md`. 6.1..6.5 each produce `.planning/phase-6/LEXICON-AUDIT-{mail,elas,flui,ther,nlse}.md` populated from that template. `[CITED: 06-UI-SPEC.md §Per-Module Conformance Contract]`

### Frequency weighting methodology

**Recommendation:** Manual judgment from `testa/*` + Fortran source inspection. No automated grep tooling in v1.

Scoring rubric (proposed by 6.0 plan):
- **H (High):** appears in ≥3 `testa/` test cases AND used in common user workflow (open, save, add/generate/refine). Target: ~10-15 per module.
- **M (Medium):** appears in 1-2 `testa/` cases OR used in specialized but non-rare workflow. Target: ~10-20 per module.
- **L (Low):** not in `testa/` and only referenced in docs or specialized paths. Target: remainder. NOT surfaced as QAction — stays typed-only.

Total per module: 20-40 QActions (D-03), long tail documented but not wired.

### Coverage verification

**Recommendation:** A CMake-integrated check OR a shell-script audit that:
1. Greps every `SAISIE`, `LECMEN`, `LECLEX`, `LIREF`, etc. dispatch in the module's `.f` files.
2. Extracts the lexicon commands (typically strings like `'5;90;'` or numeric dispatches).
3. Cross-checks against the LEXICON-AUDIT.md for that module.
4. Fails if a lexicon command is in source but missing from the audit.

This is 6.1..6.5 deliverable (the audit script is per-module because Fortran conventions vary), but 6.0 should ship a prototype `bin/audit_lexicon.sh` that 6.1's planner can adapt.

### Pitfalls

- **Pitfall 9.1 — Lexicon commands embedded in data files.** Some menu commands are read from `td/m/*` text files at runtime, not hardcoded in `.f`. These will not be surfaced by a grep of `.f` files alone. Planner must extend the audit script to also scan `td/m/` menu definition files.
- **Pitfall 9.2 — Numeric dispatches (not lexical strings).** Some commands use `READ(IMPRIM,*) N` + `GO TO (10,20,30,...) N` in Fortran — there's no string like `"5;90;"` in source. The audit must recognize this pattern and document it.

---

## 10. Per-Module Conformance Enforcement

### Design

**Problem:** how does 6.0 prevent a later 6.3 (`flui`) from silently introducing a plain-char QShortcut that intercepts a lexicon key?

**Recommendation:** Three layers:

1. **Build-time lint (6.0 ships):** a CMake custom target `verify_shortcut_modifiers` that greps all `xvue_qt_*_actions.cpp` files and fails if `QKeySequence("[a-zA-Z0-9]")` without modifier is found. Mirror of Phase 0's `verify_abi` and Phase 1's `verify_no_exec` patterns. `[VERIFIED pattern: xvue/qt/CMakeLists.txt:69-77]`

   Script: `xvue/qt/cmake/verify_shortcut_modifiers.sh` — scans source files for forbidden patterns.

2. **Runtime registration assert (6.0 ships):**
   ```cpp
   // In XvueWindow ctor after build-menu-bar:
   Q_ASSERT_X(menuBridge_->hasRegisteredModule(),
              "XvueWindow",
              "No module called xvue_module_init_ before xvinitgraphique_. "
              "6.1..6.5 plans must add CALL XVUE_MODULE_INIT('...') to "
              "prpr/pp*_qt.f before any graphics call.");
   ```
   This assertion fires if a pp*_qt ships without ever registering a module — i.e. it shows the shared `{File, View, Help}` menu bar and NOTHING else. In a debug build this halts launch; in release it prints a warning and continues (just a reduced UX).

3. **Per-module CMake target verification (6.1..6.5 ship):** Each module's `cb*_qt` script verifies `nm pp/pp*_qt | grep registerMailActions` (etc.) returns a symbol — confirming the registration function got linked in.

### Shortcut-collision audit (within one module)

Within a single module's QActions, shortcut collisions ARE Qt runtime warnings (`"ambiguous shortcut"`). A test at registration time that counts and asserts no ambiguous shortcuts:
```cpp
void registerMailActions(...) {
    // ... create actions ...
    // 6.0 ships a helper:
    XvueShortcutAudit::validateNoCollisions(win);
}
```

### Pitfalls

- **Pitfall 10.1 — Cross-module shortcut collision.** Module A and Module B don't run in the same process — `ppmail_qt` vs `ppelas_qt` — so cross-module collisions are impossible at runtime. Audit only within a module.
- **Pitfall 10.2 — Dynamic action creation.** If a 6.1..6.5 plan lazily creates QActions (e.g., recent-projects submenu entries), those bypass registration-time checks. 6.0's validation must scan ALL QActions on the window, not just the static ones. `win->findChildren<QAction*>()` covers dynamic ones.

---

## 11. Build / CMake Changes

### New source files in `xvue/qt/CMakeLists.txt` `target_sources` addition:

```cmake
add_library(xvueqt STATIC
    src/xvue_qt_api.cpp
    src/xvue_qt_app.cpp
    src/xvue_qt_canvas.cpp
    src/xvue_qt_event.cpp
    src/xvue_qt_state.cpp
    src/xvue_qt_window.cpp
    # Phase 6 additions:
    src/xvue_qt_menu_bridge.cpp     # XvueMenuBridge (§1)
    src/xvue_qt_i18n.cpp            # bilingual table (§2)
    src/xvue_qt_console_dock.cpp    # XvueConsoleDock (§3)
    src/xvue_qt_prefs.cpp           # XvuePrefs (§4)
    src/xvue_qt_about_dialog.cpp    # QMessageBox::about wrapper (D-09)
    src/xvue_qt_preferences.cpp     # XvuePreferencesDialog
    src/xvue_qt_recent_projects.cpp # Recent-projects submenu (D-06)
    src/xvue_qt_error_batcher.cpp   # *** ERREUR batcher (UX-11)
    src/xvue_qt_shortcut_audit.cpp  # §10 defensive lint
)
```

### New resource file

```cmake
qt_add_resources(xvueqt xvue_icons
    PREFIX "/xvue/qt"
    FILES
        # 6.0 ships ZERO custom SVGs (all 6.0 actions use StandardPixmap
        # per UI-SPEC §Iconography). 6.1..6.5 append.
)
```

### Build-time lint targets

```cmake
# §10 layer 1:
add_custom_target(verify_shortcut_modifiers ALL
    COMMAND ${CMAKE_COMMAND} -E echo "verify_shortcut_modifiers: scanning..."
    COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/cmake/verify_shortcut_modifiers.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/src
    DEPENDS xvueqt
    VERBATIM
)

# §5 SVG-only enforcement for icons:
add_custom_target(verify_icons_svg_only ALL
    COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/cmake/verify_icons_svg_only.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/resources/icons
    DEPENDS xvueqt
    VERBATIM
)
```

### UTF-8 source encoding

```cmake
target_compile_options(xvueqt PRIVATE -finput-charset=UTF-8)
```

### ABI symbol count update

Phase 6 adds ONE new `extern "C"` (`xvue_module_init_`), so ABI count goes **57 → 58**. `xvue/qt/cmake/verify_abi.sh` has the expected count as a parameter; planner updates to 58 during 6.0 implementation.

`[VERIFIED: xvue/qt/CMakeLists.txt:55 — "Expected count after Option A resolution of Planner Alert: 57"]`

### `bin/cbl_tout_qt` / `bin/cb*_qt` changes

**None required** — the shell scripts already link `libxvueqt.a` + Qt6. New sources are picked up automatically by CMake's `AUTOMOC` and resource system. Possibly a `-DMEFISTO_MODULE=mail` define per-module linker arg, IF we choose module identity option 3 (§5). With option 1 (recommended), no script change.

### Pitfalls

- **Pitfall 11.1 — AUTOMOC order.** `CMAKE_AUTOMOC ON` must precede `find_package(Qt6)`. Already done at line 11. New `Q_OBJECT` classes auto-picked up. ✓
- **Pitfall 11.2 — Qt6Widgets dependency transitive.** Already linked `PUBLIC Qt6::Widgets` per line 46. QDockWidget, QMenuBar, QDialog all live in Qt6Widgets. ✓
- **Pitfall 11.3 — `QSocketNotifier` is in Qt6Core.** Already linked. ✓
- **Pitfall 11.4 — Qt6PrintSupport link.** Currently linked (line 45). Phase 6 doesn't USE it; Phase 7 will. Keep as-is.

---

## 12. Testing Strategy

### QTest extension (build on Phase 5's infrastructure)

Phase 5 delivered `xvue/qt/tests/test_xvue_qt_event.cpp` with 31 QTest cases. Phase 6 adds a parallel file `xvue/qt/tests/test_xvue_qt_menu.cpp` via CMake option `XVUE_QT_BUILD_TESTS=ON` (existing pattern at `xvue/qt/CMakeLists.txt:87-90`).

### Proposed test cases (Phase 6.0)

| Test | Covers | Approach |
|------|--------|----------|
| `testMenuBridgeQueueSingleChar` | UX-02, UX-03 | Push `"5;"`, call `xvsouris_`, assert returns `{notypeevent=2, nbc='5'}`. Call again, assert returns `{notypeevent=2, nbc=';'}`. Call again (empty queue), assert AUTOEXIT or normal path. |
| `testMenuBridgeQueueMultiChar` | UX-02, UX-03 | Push `"5;90;"`, call `xvsouris_` 6 times (5 chars + CR), assert each returns the expected char. |
| `testMenuBridgeTrailingCR` | UX-03 | `queueLexicon("99;")` produces 4 queued chars ending in 13 (CR). |
| `testModalGuardDepthGreaterThanZero` | UX-04, D-08 | In a `BlockingDepthGuard` scope, trigger `File → Open` action. Assert status bar message appears, dialog does NOT open. |
| `testBilingualFlip` | UX-08 | Create a temp `MEFISTO/td/m/anglais`, reset i18n cache, verify menu labels flip to EN. Remove, reset, verify FR. (Uses temp `MEFISTO` var.) |
| `testRecentProjectsCap10` | UX-06, D-06 | Push 11 paths, verify oldest dropped, size == 10. |
| `testRecentProjectsMissingPathPrune` | D-06 | Push a non-existent path, startup prune removes it. |
| `testConsoleDockAppendLine` | UX-10 | Simulate stdout write (via direct `XvueConsoleDock::appendLine`), verify QPlainTextEdit has the text + auto-scroll. |
| `testErreurBatcher` | UX-11 | Feed 3 `*** ERREUR` lines within 200 ms, verify ONE QMessageBox with joined body. |
| `testErreurBatcherDeferredDuringBlocking` | UX-11 + T-06-03-01 | Feed `*** ERREUR` during `blockingDepth > 0`, verify MessageBox doesn't open until depth returns to 0. |
| `testQSettingsPersistRoundTrip` | UX-07 | Save geometry + state + color-scheme, reset QSettings, reload, verify identity. |
| `testDarkModePaletteApplication` | UX-13 | Set pref to "dark", call `applyColorSchemePreference`, verify `QApplication::palette().color(QPalette::Window)` is the dark value. |
| `testCanvasWheelZoom` | UX-12 | postEvent a QWheelEvent with +120 angleDelta, verify `XvueState::view_transform_` scale > 1.0. |
| `testCanvasMiddleDragPan` | UX-12 | postEvent press→move→release sequence, verify `view_transform_` translated. |
| `testCanvasContextMenuSuppressedDuringBlocking` | UX-12 + D-08 | postEvent QContextMenuEvent with `blockingDepth > 0`, verify QMenu NOT exec'd. |
| `testCanvasCoordSignal` | UX-12 | postEvent MouseMove, verify `mouseCoords` signal emitted with logical-px position. |
| `testShortcutModifierRuleLint` | D-04, §10 | Build-time CMake check — not a runtime test; scripted. |

Total ~16 QTest cases for 6.0, plus per-module adds in 6.1..6.5.

### Manual validation (human eye)

Per UI-SPEC §Dark-mode test, canvas invariance: toggle system theme, restart app, eyeball chrome flipping while canvas pixmap unchanged.

Per the critical wiring UX-03: manually click a menu item, verify Fortran's `SACLAV` accumulates the chars into `KLG` and executes the command as if typed.

### Pitfalls

- **Pitfall 12.1 — QTest offscreen.** Qt tests run on `QT_QPA_PLATFORM=offscreen`. Some APIs (like `QClipboard`) don't work offscreen. Console dock copy-to-clipboard test must tolerate this (skip if offscreen).
- **Pitfall 12.2 — Asynchronous signal connections.** Some tests will need `QSignalSpy` + `QTRY_COMPARE` to wait for signals that fire through event loop.
- **Pitfall 12.3 — `QSettings` test isolation.** Tests must use `QSettings::setPath(IniFormat, UserScope, <tmpdir>)` to avoid corrupting real user settings.

---

## 13. Known Pitfalls / Assumptions to Resolve

Resolving the A1..A5 from the topic list:

### A1 — Qt's system-theme signals on Debian XFce (without Plasma)

**Status:** `[ASSUMED]` — likely works with `qt6ct` or `qgnomeplatform` plugin installed. Verification deferred to manual test during 6.0 validation.

**Fallback:** always read current theme at startup; user restarts app to re-sync.

### A2 — QProcess stdout vs. in-process freopen race

**Status:** **RESOLVED** — in-process freopen chosen (§3). `QProcess` inapplicable because Fortran runs in-process.

### A3 — QAction shortcut with WindowShortcut context + canvas StrongFocus

**Status:** **RESOLVED** — Qt's precedence is: `ShortcutOverride` (propagated to focused widget FIRST) → shortcut dispatch → `KeyPress` on focused widget. For bare alphanumeric without shortcut binding, the shortcut dispatch skips and `KeyPress` reaches canvas directly. `[CITED: doc.qt.io/qt-6/qshortcut.html + Qt Forum discussion on ShortcutOverride precedence]`

Defensive `ShortcutOverride` filter on canvas (accepting override for plain ASCII) provides defense-in-depth against a future 6.1..6.5 plan accidentally binding a plain-char shortcut (§5 Pitfall).

### A4 — SAISIE one-char-per-call vs buffered

**Status:** **VERIFIED ONE-CHAR-PER-CALL.** See §6. `[VERIFIED: xvue/saclav.f:270-317, util/saiptc.f:23, util/searclic.f:68]`

### A5 — UI-SPEC 5 non-blocking flags surfacing real issues

UI-SPEC declared 3 PASS + 3 FLAG (non-blocking) at Phase 6 approval. Flags (paraphrased from STATE.md "3 FLAG non-blocking recommendations"):

1. **`ButtonOk` copy** — UI-SPEC uses `"OK"` for both FR and EN. Some FR conventions use `"Accepter"` for primary-action confirm. Verdict: `"OK"` is universally understood in FR; keep. **Not a plan task.**
2. **`FileRecentClear` noun** — UI-SPEC uses `"Clear Recent"` / `"Effacer la liste"`. Some users expect explicit `"Recent Projects"` in the action. Verdict: the submenu context makes it clear; keep concise. **Not a plan task.**
3. **Empty-state canvas rendering** — UI-SPEC's "No project open" message on empty canvas requires a `paintEvent` addition that checks for no-data state. Verdict: **flag for planner** — this is a Phase 6.0 task.  
   Proposal: if `state_->backing_` is null OR `state_->background_ == Qt::black` AND no draw calls have ever been made, render a centered "No project open" message in `QPalette::WindowText` via `QPainter::drawText`. NEW field needed: `bool XvueState::has_user_content_` (false at init, true after first draw primitive).
4. **QGroupBox bold** — Not 6.0 scope; flagged for 6.1..6.5 preferences dialogs.
5. **QFileDialog setDirectory** — UI-SPEC uses `qgetenv("MEFISTOX")` for default dir. If unset, QFileDialog defaults to `$HOME`. Verdict: document env-var requirement; if unset, show a status-bar warning once. **Minor plan task.**

### Additional assumptions documented in Assumptions Log below

---

## Runtime State Inventory

Phase 6 is largely greenfield (new C++ files). However, it introduces `QSettings` state that persists across sessions — an implicit state inventory is warranted.

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | **`QSettings` file at `$HOME/.config/LJLL/mefisto.conf`** — new. Contains window geometry, recent-projects list, color-scheme preference. | No migration from previous phases (no pre-existing settings). First-launch defaults apply. |
| Live service config | None — MEFISTO is standalone desktop. | None. |
| OS-registered state | None — no systemd units, no Task Scheduler. | None. |
| Secrets/env vars | `MEFISTO` (path, existing), `MEFISTOX` (path, existing); Phase 6 reads both. Phase 6 documents recommended addition: no new env vars strictly required. | Document `MEFISTOX` usage in README. |
| Build artifacts | `libxvueqt.a` grows (new `.o` files from ~7 new `.cpp`). No stale artifacts — `bin/cbl_tout_qt` cleans `xvue/build/` before build. | None. |

**Nothing found in category:** Live service config, OS-registered state — MEFISTO does not use either (pure desktop app, no daemon, no OS integration).

---

## Common Pitfalls (consolidated cross-cutting)

### Pitfall A — Phase 6 work collides with parallel cbl_tout / cbl_tout_qt runs

**What goes wrong:** Developer triggers cbl_tout_qt to test menu changes while cbl_tout is still building legacy X11 — both clobber `pp/ppmail_qt` mid-link.

**Why it happens:** Shared `pp/` output directory (documented in MEMORY.md `feedback_parallel_builds_share_pp`).

**How to avoid:** Serialize builds. The Phase 6 plan should NOT parallelize `bin/cbl_tout` and `bin/cbl_tout_qt` invocations. `[VERIFIED: MEMORY.md entry — historical regression]`

### Pitfall B — Regression on Phase 5 event-bridge invariants

**What goes wrong:** Adding menu-bridge pre-drain at the top of `waitForEvent` mistakenly bypasses the `BlockingDepthGuard` or the save/restore block.

**Why it happens:** The pre-drain is an early return — if placed BEFORE the guard, depth doesn't increment, breaking D-08.

**How to avoid:** Pre-drain is AFTER `BlockingDepthGuard depth_guard;` (so `~guard` decrements on return) but BEFORE `saved_loop = loop_; ...` (no state to save if we never enter the loop). Verified in §1.

**Warning signs:** Phase 5's existing tests (`test_xvue_qt_event.cpp`) start failing — specifically `testBlockingDepthAssertGuardRaii` or `testNestedWaitForEventDepth`.

### Pitfall C — Debug session after UI-SPEC Flag #3 (empty-state rendering)

**What goes wrong:** First-launch window shows black canvas because backing pixmap is uninitialized, and user doesn't know what to do.

**Why it happens:** Phase 1/2 backing pixmap allocation depends on `resizeEvent` firing; initial background is `Qt::black` (`xvue_qt_state.h:23`).

**How to avoid:** 6.0 adds the empty-state QPainter::drawText in `paintEvent` when `!has_user_content_`, per A5 flag.

### Pitfall D — Fortran build drift from libgfortran version pin

**What goes wrong:** User upgrades `libgfortran5` → Fortran ABI drift → MEFISTO crashes in unrelated solver paths → developer blames Phase 6.

**Why it happens:** Debian sid's gcc-16 snapshot exposed latent UB (STATE.md history). Pin must be maintained.

**How to avoid:** 6.0 adds a check at `bin/cbl_tout_qt` top: `dpkg -l | grep '^ii.*libgfortran5' | awk '{print $3}'` and verifies it's NOT > `15.2.0-9`. Warn but don't block. `[VERIFIED: MEMORY.md + STATE.md Phase 03-04 2026-04-14 close]`

---

## Code Examples

### Example 1: `XvueMenuBridge` skeleton

```cpp
// xvue_qt_menu_bridge.h
#pragma once
#include <QObject>
#include <QQueue>
#include <QString>

class QMenu;

class XvueMenuBridge : public QObject {
    Q_OBJECT
public:
    explicit XvueMenuBridge(QObject* parent = nullptr);

    // Called from QAction::triggered lambdas. Appends chars and a CR.
    void queueLexicon(const QString& cmd);

    // Called from XvueEventBridge::waitForEvent pre-drain.
    std::optional<char> popChar();

    // Called by 6.1..6.5 to register a context-menu populator callback.
    using ContextPopulator = std::function<void(QMenu*)>;
    void setContextMenuPopulator(ContextPopulator fn);
    void populateContextMenu(QMenu* m);

    // §10 layer 2: registration sentinel
    bool hasRegisteredModule() const { return moduleRegistered_; }
    void markModuleRegistered()      { moduleRegistered_ = true; }

private:
    QQueue<char>     pendingChars_;
    ContextPopulator contextPopulator_;
    bool             moduleRegistered_ = false;
};
```

### Example 2: Pre-drain integration in XvueEventBridge

```cpp
// Addition to xvue/qt/src/xvue_qt_event.cpp waitForEvent body.
// Inserted BETWEEN BlockingDepthGuard and saved_loop = loop_;

// Phase 6 UX-03: menu-queue pre-drain.
if (auto* menu = XvueApp::menuBridge()) {
    if (auto c = menu->popChar()) {
        Result r;
        r.notypeevent = 2;
        r.nbc = static_cast<unsigned char>(*c);
        r.x = canvas_ ? canvas_->mapFromGlobal(QCursor::pos()).x() : 0;
        r.y = canvas_ ? canvas_->mapFromGlobal(QCursor::pos()).y() : 0;
        // debug_logging_enabled() diagnostic intentionally reused:
        if (debug_logging_enabled()) {
            std::fprintf(stderr,
                "[xvsouris] mode=%d notypeevent=2 nbc=%d "
                "x=%d y=%d motion_count=0 depth=%d [menu-pre-drain]\n",
                static_cast<int>(mode), r.nbc, r.x, r.y,
                XvueApp::blockingDepth());
            std::fflush(stderr);
        }
        return r;
    }
}
```

### Example 3: Console-dock pipe initialization

```cpp
// xvue_qt_console_dock.cpp (excerpt)
#include <QDockWidget>
#include <QPlainTextEdit>
#include <QSocketNotifier>
#include <unistd.h>
#include <fcntl.h>
#include <cstdio>

void XvueConsoleDock::installStdoutRedirect() {
    int fd[2];
    if (::pipe(fd) != 0) {
        std::fprintf(stderr, "xvue-qt: pipe() failed; console dock disabled\n");
        return;
    }
    // Pitfall 3.2: expand pipe buffer.
    ::fcntl(fd[0], F_SETPIPE_SZ, 1 << 20);  // 1 MB

    if (::dup2(fd[1], STDOUT_FILENO) < 0) {
        std::fprintf(stderr, "xvue-qt: dup2 stdout failed\n");
        ::close(fd[0]); ::close(fd[1]); return;
    }
    ::close(fd[1]);

    // Pitfall 3.4: AFTER dup2.
    std::setvbuf(stdout, nullptr, _IOLBF, 0);

    notifier_ = new QSocketNotifier(fd[0], QSocketNotifier::Read, this);
    connect(notifier_, &QSocketNotifier::activated,
            this,      &XvueConsoleDock::onPipeReadable);
}

void XvueConsoleDock::onPipeReadable() {
    char buf[4096];
    const ssize_t n = ::read(notifier_->socket(), buf, sizeof buf);
    if (n <= 0) return;
    partialLine_.append(QByteArray::fromRawData(buf, n));
    int pos;
    while ((pos = partialLine_.indexOf('\n')) >= 0) {
        const QByteArray line = partialLine_.left(pos);
        partialLine_.remove(0, pos + 1);
        appendLine(QString::fromUtf8(line));
        if (line.startsWith("*** ERREUR") || line.startsWith("*** ERROR")) {
            errorBatcher_->enqueue(QString::fromUtf8(line));
        }
    }
}
```

### Example 4: QAction wiring with bridge queue

```cpp
// xvue_qt_mail_actions.cpp (6.1 — shown here for context)
void registerMailActions(XvueWindow* win, XvueMenuBridge* mb) {
    auto* menuBar = win->menuBar();

    // {File} — 6.0 shared
    auto* fileMenu = menuBar->findChild<QMenu*>("fileMenu");
    Q_ASSERT(fileMenu);  // 6.0 built it

    // {Mesh} — 6.1 specific
    auto* meshMenu = menuBar->addMenu(xvueT(MsgId::MeshMenuTitle));

    auto* actAddVertex = new QAction(xvueT(MsgId::MeshAddVertex), win);
    actAddVertex->setShortcut(QKeySequence("Ctrl+V"));
    actAddVertex->setIcon(QIcon(":/xvue/qt/icons/add-vertex.svg"));
    actAddVertex->setStatusTip(xvueT(MsgId::MeshAddVertexTip));
    QObject::connect(actAddVertex, &QAction::triggered, [mb]{
        mb->queueLexicon("5;90;");   // mesh→add-vertex
    });
    meshMenu->addAction(actAddVertex);
    win->toolBar()->addAction(actAddVertex);  // shared QAction

    // ...more actions...

    mb->markModuleRegistered();
}
```

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Qt 5 `QTranslator` + `.ts` for i18n | Compile-time table for fixed-vocabulary apps | Always valid; adopted here per UI-SPEC | Dramatically simpler build; no `lrelease` step |
| `QProcess` external reader for stdout | `QSocketNotifier` on pipe fd (in-process) | Qt 4+; MEFISTO context forces it | No child process; direct libgfortran integration |
| Qt 5 manual dark-mode palette | Qt 6.5+ `QStyleHints::colorScheme` + `setColorScheme` | Qt 6.5 (2023) | Linux support still DE-dependent; fallback is manual palette |
| Bare-letter QShortcuts freely allowed | Modifier-required by convention + build-time lint | Phase 6 D-04 | Preserves typed-lexicon compatibility |
| Menu actions dispatch directly into app code | Queue-synthetic-event pattern for legacy Fortran interop | Phase 6 novel | Enables "GUI front + text backend" without touching backend |

**Deprecated / outdated:**
- `QApplication::setStyle("fusion")` to force dark-mode — brittle; use `QPalette` instead. Qt docs recommend `QStyleHints` in 6.5+.
- `QMovie` for animated GIF — Phase 7 will decide per EXPORT-01 probe; not Phase 6 concern.

---

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| `qt6-base-dev` | All chrome work | ✓ | 6.10.2+dfsg-7 | — |
| `libqt6widgets6` | QMenuBar, QDockWidget, QDialog | ✓ | 6.10.2+dfsg-7 | — |
| `libqt6core6t64` | QSocketNotifier, QSettings | ✓ | 6.10.2+dfsg-7 | — |
| `libqt6gui6` | QPalette, QStyleHints | ✓ | 6.10.2+dfsg-7 | — |
| `libqt6printsupport6` | (Phase 7 — not 6) | ✓ | 6.10.2+dfsg-7 | — |
| `gfortran` | Fortran link | ✓ | 15.2.0 (pinned) | — |
| `libgfortran5` | Runtime | ✓ | 15.2.0-9 (pinned per STATE.md) | Must NOT upgrade to 16.x until Fortran UB sites are fixed |
| `pkg-config` | `bin/cbmail_qt` uses pkg-config | ✓ | present | — |
| `bash` | Shell build | ✓ | present | — |
| `$MEFISTO` env | i18n flag lookup + launcher | ✓ | `/home/drico/git/mefisto` | If unset, i18n defaults to FR |
| `$MEFISTOX` env | Project dir (for QFileDialog default) | varies | — | Falls back to `$HOME` in QFileDialog |
| `$HOME/.config/` writable | QSettings | ✓ | — | — |

**No missing dependencies.**

**Build-environment pin (documented):** `libgfortran5` stays at `15.2.0-9` until gcc-16 stabilizes OR Fortran UB fixed. `[VERIFIED: MEMORY.md feedback_debian_sid_libgfortran + STATE.md Phase 03-04 close]`

---

## Validation Architecture

This section is REQUIRED because `.planning/config.json::workflow.nyquist_validation = true`.

### Test Framework

| Property | Value |
|----------|-------|
| Framework | QTest (Qt 6.10.2 Test module) |
| Config file | `xvue/qt/CMakeLists.txt` — `option(XVUE_QT_BUILD_TESTS "..." OFF)` gated via `bin/cbl_tout_qt` |
| Test source dir | `xvue/qt/tests/` |
| Existing test target | `xvue_qt_event_tests` (Phase 5 — 31 cases) |
| New test target (Phase 6.0) | `xvue_qt_menu_tests` (≈16 cases per §12) |
| Quick run command (single test case) | `cd xvue/qt/build/tests && ./xvue_qt_menu_tests -functions | head -20` (list); `./xvue_qt_menu_tests testMenuBridgeQueueSingleChar` |
| Full suite command | `bin/cbl_tout_qt && cd xvue/qt/build/tests && ./xvue_qt_event_tests && ./xvue_qt_menu_tests` |

### Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| UX-01 | QMenuBar / QToolBar / QStatusBar composed | unit | `./xvue_qt_menu_tests testWindowChromeComposition` | ❌ Wave 0 |
| UX-02 | `XvueMenuBridge::pendingChars_` queue grow/drain | unit | `./xvue_qt_menu_tests testMenuBridgeQueueSingleChar` + `testMenuBridgeQueueMultiChar` | ❌ Wave 0 |
| UX-03 | Synthetic `notypeevent=2` returned via xvsouris_ | unit | `./xvue_qt_menu_tests testMenuBridgeQueueMultiChar` + `testMenuBridgeTrailingCR` | ❌ Wave 0 |
| UX-04 | Modal refused during `blockingDepth > 0` | unit | `./xvue_qt_menu_tests testModalGuardDepthGreaterThanZero` | ❌ Wave 0 |
| UX-05 | LEXICON-AUDIT.md exists per module | manual | `test -f .planning/phase-6/LEXICON-AUDIT-mail.md && ... ` | 6.1..6.5 |
| UX-06 | QFileDialog + recent-projects persist | unit + manual | `./xvue_qt_menu_tests testRecentProjectsCap10` + `testRecentProjectsMissingPathPrune` | ❌ Wave 0 |
| UX-07 | QSettings round-trip | unit | `./xvue_qt_menu_tests testQSettingsPersistRoundTrip` | ❌ Wave 0 |
| UX-08 | FR/EN flip via flag | unit | `./xvue_qt_menu_tests testBilingualFlip` | ❌ Wave 0 |
| UX-09 | About dialog shows Alain Perronnet | manual | Visual run `pp/ppmail_qt` + Help → About | 6.0 human-gate |
| UX-10 | Console dock captures stdout via pipe | unit | `./xvue_qt_menu_tests testConsoleDockAppendLine` + manual | ❌ Wave 0 |
| UX-11 | `*** ERREUR` surfaces as QMessageBox | unit | `./xvue_qt_menu_tests testErreurBatcher` + `testErreurBatcherDeferredDuringBlocking` | ❌ Wave 0 |
| UX-12 | Wheel zoom + middle-drag + right-click + coords | unit | `./xvue_qt_menu_tests testCanvasWheelZoom` + `testCanvasMiddleDragPan` + `testCanvasContextMenuSuppressedDuringBlocking` + `testCanvasCoordSignal` | ❌ Wave 0 |
| UX-13 | Dark-mode chrome; canvas unaffected | unit + manual | `./xvue_qt_menu_tests testDarkModePaletteApplication` + human visual | ❌ Wave 0 |

### Sampling Rate

- **Per task commit:** `cd xvue/qt/build/tests && ./xvue_qt_menu_tests <single-test>` — sub-second
- **Per wave merge:** Full suite `./xvue_qt_event_tests && ./xvue_qt_menu_tests` — <60 seconds
- **Phase gate (`/gsd-verify-work`):** Full suite green + manual UI sweep per UI-SPEC Dimension 1-6 + 5 baseline `testa/` cases still interactive on Qt backend

### Wave 0 Gaps

The Wave 0 scaffold for Phase 6.0 must create BEFORE implementation plans can land:

- [ ] `xvue/qt/src/xvue_qt_menu_bridge.{h,cpp}` — bridge class stubs (empty body)
- [ ] `xvue/qt/src/xvue_qt_i18n.{h,cpp}` — `MsgId` enum + empty table
- [ ] `xvue/qt/src/xvue_qt_console_dock.{h,cpp}` — dock class stubs
- [ ] `xvue/qt/src/xvue_qt_prefs.{h,cpp}` — QSettings helper stubs
- [ ] `xvue/qt/src/xvue_qt_about_dialog.{h,cpp}` — About launcher stub
- [ ] `xvue/qt/src/xvue_qt_preferences.{h,cpp}` — preferences dialog stub
- [ ] `xvue/qt/src/xvue_qt_recent_projects.{h,cpp}` — recent projects submenu stub
- [ ] `xvue/qt/src/xvue_qt_error_batcher.{h,cpp}` — error batcher stub
- [ ] `xvue/qt/src/xvue_qt_shortcut_audit.{h,cpp}` — shortcut lint helper stub
- [ ] `xvue/qt/tests/test_xvue_qt_menu.cpp` — new test file (mirrors test_xvue_qt_event.cpp structure)
- [ ] `xvue/qt/cmake/verify_shortcut_modifiers.sh` — new lint script
- [ ] `xvue/qt/cmake/verify_icons_svg_only.sh` — new lint script
- [ ] `xvue/qt/resources/xvue_icons.qrc` — empty resource file
- [ ] Update `xvue/qt/CMakeLists.txt` `target_sources` block (+ 9 new files)
- [ ] Update `xvue/qt/cmake/verify_abi.sh` expected count: 57 → 58
- [ ] `xvue/qt/tests/CMakeLists.txt` — add `xvue_qt_menu_tests` target
- [ ] `xvue/qt/src/xvue_qt_window.{h,cpp}` — extend XvueWindow with menuBridge accessor + buildMenuBar/buildToolBar/buildStatusBar/buildConsoleDock methods
- [ ] `xvue/qt/src/xvue_qt_app.{h,cpp}` — add `menuBridge()` forwarding accessor + `applyColorSchemePreference()` + stdout redirect init
- [ ] `xvue/qt/src/xvue_qt_canvas.{h,cpp}` — add `wheelEvent`/`mousePressEvent`/`mouseMoveEvent`/`mouseReleaseEvent`/`contextMenuEvent` overrides + `mouseCoords` signal + `resetView` slot
- [ ] `xvue/qt/src/xvue_qt_state.{h,cpp}` — add `QTransform view_transform_` + `bool has_user_content_`
- [ ] `.planning/phases/06-level-3-ux-chrome-menu-surface/LEXICON-AUDIT-TEMPLATE.md` — schema template

---

## Dependency Order for 6.0 → 6.1..6.5

Research-driven plan-ordering insights for the planner:

### 6.0 DOES ship (6.0 is the FOUNDATION):

1. `XvueMenuBridge` class (queue, pre-drain, context-menu populator hook).
2. `XvueEventBridge::waitForEvent` pre-drain integration (one-line addition).
3. `XvueConsoleDock` with stdout pipe + `*** ERREUR` batcher.
4. `XvuePrefs` QSettings helper.
5. `XvuePreferencesDialog` + `XvueAboutDialog`.
6. `XvueRecentProjectsMenu`.
7. `xvue_qt_i18n` (bilingual table — all UI-SPEC strings).
8. `xvue_module_init_` Fortran hook (new extern "C"; ABI 57→58).
9. Canvas gesture extensions: `wheelEvent`, middle-drag pan, `contextMenuEvent`, `mouseCoords` signal.
10. Dark-mode palette application.
11. Shared `{File, View, Help}` menu + 6.0 QActions (Open/Save/SaveAs/Export/Quit + Preferences + About + Documentation + view toggles).
12. Shared toolbar (Open/Save/ZoomIn/ZoomOut/Fit/Console-toggle).
13. Status bar (permanent coord label + transient message zone).
14. Build-time lint targets (`verify_shortcut_modifiers`, `verify_icons_svg_only`).
15. `xvue_qt_menu_tests` QTest suite (~16 cases).
16. `LEXICON-AUDIT-TEMPLATE.md`.
17. Empty-state canvas rendering (UI-SPEC Flag #3 follow-through).

### 6.1..6.5 DEPEND on 6.0 shipping first:

Each per-module phase adds:

1. Their own `registerXxxActions(XvueWindow*, XvueMenuBridge*)` C++ function.
2. Their own `.planning/phase-6/LEXICON-AUDIT-{module}.md` (from template).
3. Their own SVG icons (added to `xvue_qt_icons.qrc`).
4. Their own `CALL XVUE_MODULE_INIT('xxx');` in `prpr/ppXxx.f`.
5. Their own i18n table additions (append-only to `xvue_qt_i18n.cpp`).
6. Their own module-specific QTest cases (append to `xvue_qt_menu_tests`).

### Serial vs parallel execution of 6.1..6.5

**Research finding:** 6.1..6.5 CAN run in parallel IF:
- Each module touches only `prpr/ppXxx.f` (the one Fortran file per module) + creates new C++ files (no merge conflicts).
- The i18n table append is sequential (one `MsgId` per module batch).
- SVG icon additions are sequential (single `.qrc` file).

**Verdict:** Sequential execution is safer for the solo-dev posture. Plan for **6.1 → 6.2 → 6.3 → 6.4 → 6.5** in that order (roadmap natural order: mesher first since `testa/pan2d` is the validation anchor).

---

## Security Domain

Phase 6 is a desktop GUI application; security_enforcement applies minimally.

### Applicable ASVS Categories

| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | no | No auth — single-user local app |
| V3 Session Management | no | No sessions |
| V4 Access Control | no | Filesystem-level access only |
| V5 Input Validation | **yes (limited)** | QFileDialog path validation; recent-projects list entry validation (startup prune) |
| V6 Cryptography | no | No cryptographic operations |
| V9 Communications | no | No network |
| V10 Malicious Code | no | No plugin/extension system (out of scope) |
| V11 Business Logic | **yes** | Modal guard (D-08) prevents re-entrant `QDialog::exec` that could corrupt Fortran state |

### Known Threat Patterns for Qt/C++ desktop GUI

| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Path traversal in recent-projects file | Tampering | `QFileInfo::canonicalFilePath()` + `QFileInfo::exists()` check at startup; prune non-existent silently |
| Malformed `anglais` flag file | Tampering | `QFileInfo::exists()` is boolean; content never parsed |
| QSettings ini injection (newlines in values) | Tampering | Qt handles escaping; don't inject raw user input into section keys |
| Pipe-reader unbounded memory | DoS | `QPlainTextEdit::setMaximumBlockCount(10000)` caps log; pipe buffer capped via `F_SETPIPE_SZ` |
| Menu-queue unbounded growth | DoS | Cap `pendingChars_` at 10k entries (Pitfall 1.1) |
| Dark-mode CSS injection | N/A (not CSS; Qt palette) | None needed |
| SVG icon XXE / external entities | Tampering | `QIcon` uses Qt SVG which DISABLES external entities by default in Qt 6.5+. Verified. |
| `QDesktopServices::openUrl` with user-controlled URL (Help → Documentation) | Redirect | URL is a compile-time constant `file://$MEFISTO/doc/...`. No user input. ✓ |

---

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Qt's `colorSchemeChanged` signal fires on Debian XFce with `qt6ct` plugin | §8, §13 | Dark-mode live switch doesn't work; user must restart — acceptable fallback, already documented |
| A2 | Module-identity hook (`xvue_module_init_`) is preferable to argv[0] auto-detect or build-time `#define` | §4, §5 | If user prefers option 2/3, planner switches; no architectural rework |
| A3 | `setvbuf(stdout, NULL, _IOLBF, 0)` is sufficient to force libgfortran to flush on newline | §3 | If insufficient, need `FLUSH` calls OR a `extern "C" xvue_qt_flush_stdout_` hook. Test empirically in Wave 0 scaffold. |
| A4 | Pipe back-pressure will not deadlock during normal mesh output | §3 | If long Fortran reports cause deadlock, bump pipe size (1MB) OR add a threaded reader. Initial 1MB setting is defensive. |
| A5 | Always-redirect stdout (vs `isatty` conditional) is the right call | §3 | Terminal users lose stdout visibility; document env-var escape hatch `MEFISTO_CONSOLE_DOCK_PASSTHROUGH` |
| A6 | MiddleButton parity: Phase 5 bridge should NOT abort in Souris mode during 6.0 (pan-reserved) | §7 | If X11 reference requires abort in Souris mode, planner adds modifier (Ctrl+MMB) for pan instead |
| A7 | `*** ERREUR` QMessageBox deferred while `blockingDepth > 0` | §3, §13 | Alternative: fire immediately (Qt allows nested dialogs). Verify via Phase 5 regression tests. |
| A8 | Compile-time i18n table, not Qt Linguist | §2 | If project later adds languages (Spanish, German), table still scales to 3-4 langs; 5+ tips toward Linguist |
| A9 | Fortran `SACLAV` loop reads ONE char per xvsouris_ call — verified `[VERIFIED]`, not assumed | §6 | Not an assumption |
| A10 | `bin/cb*_qt` scripts don't need changes for 6.0 new sources | §11 | CMake AUTOMOC picks up new files automatically; confirmed by existing pattern |
| A11 | `QSocketNotifier` read-mode on a pipe fd works on Linux (it does; this is its documented purpose) | §3 | Not really an assumption; verified via Qt docs |

**Items flagged `[ASSUMED]` requiring user confirmation in discuss-phase:**
- A1 (dark-mode signal on XFce — defer to manual test)
- A2 (module-identity mechanism — explicit hook preferred but not locked)
- A3 (stdout flushing — empirical test)
- A5 (stdout redirect conditional)
- A6 (middle-button semantics re-check against Phase 5 / X11 reference)
- A7 (deferred error-box during blocking read)

---

## Open Questions

1. **Module-identity mechanism (A2).**
   - What we know: Three options viable (explicit hook, argv[0], build-time define); explicit hook adds 1 ABI entry.
   - What's unclear: User preference.
   - Recommendation: Propose explicit hook; let planner/discuss-phase confirm.

2. **MiddleButton parity for pan (A6).**
   - What we know: Phase 5 bridge at `xvue_qt_event.cpp:329` aborts Fortran on MMB; UX-12 wants MMB for pan.
   - What's unclear: Does xvuelc.c reference treat MMB differently in Souris vs Souris2 modes? (Need to re-read xvuelc.c:2183-2315.)
   - Recommendation: 6.0 plan's first task is a 30-minute xvuelc.c cross-check. Most likely outcome: MMB aborts ONLY in Souris2 (picking) — in Souris mode, MMB was unused, so pan takes over without breakage.

3. **First-launch empty-state rendering details (UI-SPEC Flag #3).**
   - What we know: UI-SPEC says "Aucun projet ouvert / No project open" + "Choose File → Open Project" body.
   - What's unclear: Exact font size, position (centered vs top), color.
   - Recommendation: Centered, `QApplication::font()` size, `QPalette::WindowText`, rendered via `QPainter::drawText(backing->rect(), Qt::AlignCenter, text)`.

4. **Preferences dialog tabs vs single-panel.**
   - What we know: ~4 preferences (recent-projects count, console default visibility, color-scheme, TBD).
   - What's unclear: UX choice.
   - Recommendation: Single `QFormLayout` panel — <5 fields doesn't justify tabs.

5. **Recent-projects "Clear" confirmation copy.**
   - What we know: UI-SPEC destructive table says title `"Clear Recent Projects?"` / Body `DestructiveConfirmBodyGeneric` / Buttons OK+Cancel.
   - What's unclear: Should "Clear Recent" be a distinct `MsgId` for the menu entry vs a distinct dialog title?
   - Recommendation: Two separate `MsgId`s (`FileRecentClear` for the menu; a new `FileRecentClearConfirm` for the dialog title).

6. **PDF documentation location (UX-09).**
   - What we know: `doc/` only contains `normes.ps` (PostScript, not PDF).
   - What's unclear: Where is user-facing documentation?
   - Recommendation: Survey `doca/` vs `docf/`; if no PDF user guide exists, F1 launches `normes.ps` via `QDesktopServices::openUrl` (ghostscript opens PostScript on Linux). This is a weak link — flag for Phase 9 user-docs pass.

---

## Sources

### Primary (HIGH confidence)

- **`xvue/saclav.f:1-347`** — Full Fortran loop reading one char per xvsouris_ call. Ground truth for §6.
- **`util/saiptc.f`** — SAIPTC wrapper around xvsouris; confirms char-at-a-time semantic.
- **`util/searclic.f:68-100`** — Fortran caller using xvsouris in a loop. Confirms pattern.
- **`util/donnmf.f`** — Line-parser consuming accumulated `KLG(...)` strings at `;`-delimited tokens.
- **`util/lirlig.f`** — Line-reader calling `EVMENU` → `SACLAV` (keyboard mode) or file read (batch mode).
- **`util/langue.f:22-32`** — Fortran-side bilingual flag mechanism. Source of truth for `$MEFISTO/td/m/anglais`.
- **`incl/langue.inc`** — `LANGAG` COMMON block declaration.
- **`xvue/qt/src/xvue_qt_event.{h,cpp}`** — Phase 5 bridge — the pattern Phase 6 mirrors.
- **`xvue/qt/src/xvue_qt_app.{h,cpp}`** — `XvueApp::blockingDepth()` + `ensure()` with `AA_CompressHighFrequencyEvents`.
- **`xvue/qt/src/xvue_qt_window.{h,cpp}`** — Phase 1/5 window scaffold — Phase 6 grows it.
- **`xvue/qt/src/xvue_qt_canvas.{h,cpp}`** — Canvas with `StrongFocus` + `setMouseTracking(true)` + DARK-MODE defensive guards (`WA_OpaquePaintEvent`).
- **`xvue/qt/src/xvue_qt_state.{h,cpp}`** — State container (Phase 6 adds `view_transform_` + `has_user_content_`).
- **`xvue/qt/CMakeLists.txt`** — Build system (Phase 6 adds 9 new source files + 2 lint targets + icon resource).
- **`xvue/qt/README.md`** — Phase 5 event-bridge architecture.
- **`xvue/qt/tests/test_xvue_qt_event.cpp`** — Phase 5 QTest infrastructure (pattern for Phase 6 tests).
- **`bin/cbmail_qt`** (+ other `cb*_qt`) — Qt linker scripts. Phase 6 needs no changes.

### Secondary (Qt 6.10 docs — MEDIUM confidence; cross-verified with official source)

- [QSocketNotifier Class — Qt 6.10.1](https://doc.qt.io/qt-6/qsocketnotifier.html) — non-blocking fd monitoring; recommended approach for stdin/pipe.
- [QShortcut Class — Qt 6.11](https://doc.qt.io/qt-6/qshortcut.html) — `WindowShortcut` is default; `ShortcutOverride` event precedes `KeyPress`.
- [QAction Class — Qt 6.11](https://doc.qt.io/qt-6/qaction.html) — single action shared between menu + toolbar.
- [QSettings Class — Qt 6.11](https://doc.qt.io/qt-6/qsettings.html) — organization + application name; `beginGroup` per-module.
- [QStyleHints Class — Qt 6.10.2](https://doc.qt.io/qt-6/qstylehints.html) — `colorScheme` property + `colorSchemeChanged` signal since Qt 6.5.
- [Dark Mode on Windows 11 with Qt 6.5 — Qt Blog](https://www.qt.io/blog/dark-mode-on-windows-11-with-qt-6.5) — background on Qt's dark-mode story (Linux coverage limited).
- [Qt detect system dark/light mode? — Qt Forum](https://forum.qt.io/topic/147263/qt-detect-system-dark-light-mode) — practical pattern on Linux.
- [QWidget::contextMenuEvent — Qt 6.11](https://doc.qt.io/qt-6/qwidget.html) — confirms event-order for right-click.
- [QMainWindow::saveState / saveGeometry — Qt 6.11](https://doc.qt.io/qt-6/qmainwindow.html) — QByteArray round-trip for QSettings.

### Tertiary (LOW confidence — community patterns, verified only by single source)

- [Qt Forum thread "reading multi-line input from QSocketNotifier"](https://forum.qt.io/topic/158084/reading-multi-line-input-from-qsocketnotifier) — confirms main-thread read pattern; not an official doc.
- [Stack Overflow / Qt Centre "Non-blocking stream read"](https://www.qtcentre.org/threads/20268-Non-blocking-stream-read) — `QSocketNotifier` vs `QTextStream` — informal.

### In-repo prior research (HIGH)

- `.planning/phases/05-event-bridge-blocking-reads/05-RESEARCH.md` — nested `QEventLoop` semantics, motion coalescing.
- `.planning/phases/06-level-3-ux-chrome-menu-surface/06-CONTEXT.md` — locked decisions D-01..D-10.
- `.planning/phases/06-level-3-ux-chrome-menu-surface/06-UI-SPEC.md` — approved design contract with copywriting table.
- `.planning/PROJECT.md` — core value invariant + scope exclusions.
- `.planning/ROADMAP.md` — phase 6 goal + research flag for sub-phase split.
- `.planning/REQUIREMENTS.md` — UX-01..UX-13 authoritative text.
- `.planning/codebase/STRUCTURE.md`, `STACK.md`, `CONVENTIONS.md` — codebase maps.
- `.planning/STATE.md` — Phase 5 close state, UI-SPEC approval flags.
- `CLAUDE.md` — coding norms, working rules, libgfortran pin.
- `~/.claude/projects/-home-drico-git-mefisto/memory/MEMORY.md` — libgfortran 15.2.0-9 pin + parallel-build pp/ clobber warning.

---

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — Qt 6.10.2 verified via `dpkg -l`; all classes docs from Qt 6.10 official.
- Architecture: HIGH — patterns mirror Phase 5 (established, production-validated).
- Synthetic-event wiring (§6): HIGH — char-at-a-time confirmed by reading 5+ Fortran files (`saclav.f`, `saiptc.f`, `searclic.f`, `donnmf.f`, `lirlig.f`).
- i18n strategy: HIGH — file-flag mechanism is load-bearing in 30+ files; switching to Linguist would break the invariant.
- Console dock: MEDIUM — `freopen` + `QSocketNotifier` pattern is standard but specific pipe-buffer behavior depends on libgfortran internals; requires Wave 0 empirical test (A3, A4).
- QSettings: HIGH — single pattern is Qt-documented idiom.
- QAction registration: HIGH — shortcut precedence verified via Qt docs + forum; module-identity hook is a design proposal (MEDIUM).
- Canvas interactions: MEDIUM — wheel/pan straightforward; MMB conflict with Phase 5 (A6) needs confirmation.
- Dark-mode: MEDIUM — Qt docs HIGH for API existence; Linux DE coverage varies (A1).
- LEXICON-AUDIT schema: HIGH — UI-SPEC provides authoritative schema.
- Per-module conformance: MEDIUM — defense-in-depth design; empirical test TBD.
- Build / CMake: HIGH — pattern established in Phase 0-5.
- Testing: HIGH — QTest framework proven in Phase 5.

**Research date:** 2026-04-18

**Valid until:** 2026-06-01 (6 weeks). After that, re-verify Qt 6.10 → 6.11 changes if upgraded.
