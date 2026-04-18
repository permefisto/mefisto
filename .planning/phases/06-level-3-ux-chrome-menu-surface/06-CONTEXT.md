# Phase 6: Level-3 UX chrome & menu surface - Context

**Gathered:** 2026-04-19
**Status:** Ready for planning (requires ROADMAP split: 6 → 6.0 + 6.1..6.5 before `/gsd-plan-phase`)

<domain>
## Phase Boundary

Full Qt `QMainWindow` shell — menu bar, toolbar, status bar, dock widgets, dialogs, HiDPI and dark-mode chrome — dispatching user actions into the existing Fortran text-lexicon via the `QAction::triggered` → `XvueMenuBridge::pendingCommands_` queue → next `xvsouris_` drains and returns synthetic `notypeevent=2` keyboard event pattern.

**Invariant:** zero changes to solver drivers in `mail/`, `elas/`, `flui/`, `ther/`, `nlse/`, `reso/`, `util/`. The typed lexicon (`99;`, `5;90;`, …) remains the always-available fallback — this phase is additive UX, not a replacement.

</domain>

<decisions>
## Implementation Decisions

### Phase structure (meta)
- **D-01:** Split Phase 6 into **6.0** (shared shell: QMainWindow + `XvueMenuBridge` + universal dialogs + `QSettings` persistence + dark-mode chrome) plus **6.1..6.5** (per-module lexicon audit + menu wiring: 6.1 mail, 6.2 elas, 6.3 flui, 6.4 ther, 6.5 nlse). ROADMAP.md must be updated to reflect the split before `/gsd-plan-phase` runs. 6.0 unblocks early A/B of the shell alone; each module phase is an independent release cycle.

### Menu taxonomy
- **D-02:** **Per-module dynamic menu bar.** Each `pp*_qt` executable defines its own top-level menus. 6.0 ships the shared `{File, View, Help}` plus the menu-bar mechanism; each 6.1..6.5 declares and plugs in its module-specific menus (e.g., `ppmail_qt` = `{File, Edit, Mesh, View, Help}`, `ppelas_qt` = `{File, Solve, View, Help}`, etc.). Rationale: a mesher is genuinely a different tool than a solver; a unified static bar with grayed-out items would surface irrelevant menus. `QMenuBar::addMenu` happens in the module's own `main` / construction path after the shared shell is built.

### LEXICON-AUDIT scope
- **D-03:** **Audit exhaustive, QActions frequency-weighted.** For each of 6.1..6.5, `LEXICON-AUDIT.md` catalogs EVERY typed lexicon command (full documentation — the long tail is preserved on paper). But `QAction` wiring covers only the **80/20 subset** (approximately 20–40 most-used commands per module). Long-tail commands keep working via the existing typed lexicon — nothing regresses, just not GUI-surfaced. Frequency weighting is judged from `testa/` cases + module drivers; planner decides the exact cutoff per module.

### Shortcuts + typed-lexicon coexistence
- **D-04:** **Modifier rule.** `Ctrl/Alt/Cmd+X` → QAction; plain alphanumeric + digits + `;` + Esc + Return → typed lexicon via `XvueEventBridge`. F-keys (`F1..F12`) route to QAction. Preserves `99;`, `5;90;1;` muscle memory exactly. Esc (27) / `@` (64) continue to follow Phase 5 D-06 abort path. Implemented by scoping QShortcuts to the `QMainWindow` (not widget-local) while the canvas keeps `StrongFocus` — Qt routes non-modifier keys directly to the focused canvas before any shortcut lookup, and modifier combos match the shortcut table first.

### Dark-mode scope (Claude's Discretion default — confirmed)
- **D-05:** **System-follow via `QPalette`.** Chrome (menus / toolbars / status bar / dock widgets / dialogs) follows the system dark/light palette automatically. Canvas pixels, axes labels, and scientific colormaps are untouched — the canvas ignores `QPalette` and continues using the palette loaded via `PALCDE` / `XvCouleursImposees`. User can force light/dark via Preferences dialog (stored in `QSettings` as `ui/color-scheme`).

### Recent-projects list (Claude's Discretion default — confirmed)
- **D-06:** 10 most-recent projects, shown as File → Recent Projects submenu. "Clear Recent" action at the bottom. Each entry stores absolute `$MEFISTOX/<project>/` path. Startup check validates the path still exists; missing entries are pruned silently.

### Console dock behavior (Claude's Discretion default — confirmed)
- **D-07:** Visible by default on first launch (state persisted via `QSettings` after that). Auto-scroll to bottom on new output. Copy-to-clipboard action in the dock's context menu. Session-local — log resets each process start (no cross-session persistence). Lines matching `*** ERREUR` parsed and surfaced as `QMessageBox::warning` with the matched line as message body.

### Modal dialog re-entrancy guard (Claude's Discretion default — confirmed)
- **D-08:** When `XvueApp::blockingDepth() > 0` (inside a nested `xvsouris_`), any `QDialog::exec` path refuses silently and shows a 3-second status-bar message: `"Finish current operation first (type 99;)"` (FR: `"Terminez l'opération en cours (tapez 99;)"`). The QAction itself is NOT disabled — the guard lives inside the QDialog wrapper so the menu item remains clickable; the message explains why nothing opened. Inherits Phase 5 D-03 `blockingDepth()` accessor.

### About dialog (Claude's Discretion default — confirmed)
- **D-09:** Credits **Alain Perronnet (LJLL / UPMC Paris)**. MEFISTO version derived from `incl/homdir.inc` mtime or a new `incl/version.inc` (planner decides). Qt version via `qVersion()`. Build date via `__DATE__` at compile time. Both FR and EN text selectable via the existing `$MEFISTO/td/m/anglais` flag (Phase UX-08).

### Toolbar icon source (Claude's Discretion default — confirmed)
- **D-10:** Qt built-in `QStyle::StandardPixmap` where semantically appropriate (`SP_DialogSaveButton`, `SP_DialogOpenButton`, `SP_MediaPlay`, `SP_BrowserReload`, etc.). Custom SVG under `xvue/qt/resources/icons/` where no Qt built-in fits (mesh-specific actions: add-vertex, add-edge, delete-face, etc.). Icon resource registered via Qt resource system (`.qrc`).

### Claude's Discretion
- `XvueMenuBridge` ownership (likely `XvueApp::menuBridge()` static singleton, but planner may put it on `XvueWindow` if that's cleaner — mirrors Phase 5 D-02 for `XvueEventBridge`).
- Whether preferences dialog is one tab or multiple (recent-projects count, console visibility default, dark-mode override can all fit on one tab).
- How the About dialog renders multi-line credits (plain label vs `QLabel` with HTML).
- Exact `LEXICON-AUDIT.md` markdown schema per module (planner will propose a shared template in 6.0).
- Whether the `QProcess` stdout pipe-reader lives in `XvueApp` (single global) or is instantiated per-module in `pp*_qt::main`.
- Whether 6.1..6.5 execute sequentially or can run in parallel (depends on whether the module audits touch shared files — likely need some serialization).

### Folded Todos

None — no pending todos matched Phase 6 scope.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase 6 scope
- `.planning/ROADMAP.md` §"Phase 6: Level-3 UX chrome & menu surface" — Goal, depends-on, 13 requirements (UX-01..UX-13), 5 success criteria, research flag for sub-phase split
- `.planning/REQUIREMENTS.md` — UX-01..UX-13 full text (QMenuBar, XvueMenuBridge, QDockWidget console, QSettings persistence, etc.)
- `.planning/PROJECT.md` — Core value invariant ("Fortran must not notice the migration"), bilingual FR/EN constraint, single-developer no-CI posture

### Phase 5 decisions carried forward
- `.planning/phases/05-event-bridge-blocking-reads/05-CONTEXT.md` — D-02 (bridge-as-QObject pattern for menu bridge), D-03 (`blockingDepth()` accessor for modal guard), D-06 (Esc=27, @=64 abort path), D-07 (key translation contract), D-09 (Wayland `QCursor::setPos` caveat), D-10 (`MEFISTO_XVSOURIS_AUTOEXIT` preserved for headless test)
- `.planning/phases/05-event-bridge-blocking-reads/05-RESEARCH.md` — Nested `QEventLoop` semantics, event-filter dispatch patterns
- `xvue/qt/README.md` — Phase 5 event-bridge architecture, D-09 Wayland caveat, `MEFISTO_XVSOURIS_DEBUG` diagnostic counter

### Codebase maps
- `.planning/codebase/STRUCTURE.md` — Module layout (`mail/`, `elas/`, `flui/`, `ther/`, `nlse/`, `util/`, `xvue/`, `pp/`)
- `.planning/codebase/CONVENTIONS.md` — Fortran 77 fixed-form norms, `include` file discipline, `CALL` → `extern "C"` wrapper pattern
- `.planning/codebase/STACK.md` — Qt 6 from apt, CMake owns `xvue/` only, shell-script build for everything else
- `.planning/codebase/INTEGRATIONS.md` — X11/ImageMagick/gfortran toolchain
- `doc/normes.ps` — MEFISTO project coding norms (PostScript — `evince` / `gs` to view) — Fortran norm compliance required for any Fortran edits (N/A for this phase; UI work is C++/Qt only)

### Lexicon / menu system (to be audited per-module in 6.1..6.5)
- `td/m/anglais` — Bilingual FR/EN flag file (exists → English menus; absent → French)
- `td/m/` (mesh menu files), `td/ma/` (English mesh menus), `td/mf/` (French mesh menus)
- `mail/*.f`, `elas/*.f`, `flui/*.f`, `ther/*.f`, `nlse/*.f` — Each module's interactive drivers contain the `CALL SAISIE` / `CALL LECMEN` lexicon-dispatch code that the audit enumerates

### Existing Qt shell
- `xvue/qt/src/xvue_qt_window.h` / `.cpp` — `XvueWindow : QMainWindow` (shell) from Phase 1; `bridge()` accessor from Phase 5
- `xvue/qt/src/xvue_qt_app.h` / `.cpp` — `XvueApp` singleton with `window_slot()`, `blockingDepth()`, `AA_CompressHighFrequencyEvents` setup
- `xvue/qt/src/xvue_qt_canvas.h` / `.cpp` — `XvueCanvas` with `StrongFocus`, `setMouseTracking(true)`, resize invalidation

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **`XvueApp` singleton** (`xvue/qt/src/xvue_qt_app.{h,cpp}`): `ensure()` / `window_slot()` / `blockingDepth()` / `qapp()`. Adding `menuBridge()` accessor is an obvious symmetric extension for 6.0's `XvueMenuBridge`.
- **`XvueWindow : QMainWindow`** (`xvue/qt/src/xvue_qt_window.{h,cpp}`): Already a `QMainWindow`. Adding `QMenuBar`, `QToolBar`, `QStatusBar`, `QDockWidget` is standard Qt API on top — no structural change.
- **Bridge-as-QObject pattern** (`xvue/qt/src/xvue_qt_event.{h,cpp}`): Phase 5's `XvueEventBridge` established the pattern — `XvueMenuBridge` mirrors it (QObject, event queue, accessor from window or app).
- **`BlockingDepthGuard`** (`xvue/qt/src/xvue_qt_event.h`): RAII guard available for any future nested call inside menu handlers.
- **`XvueState` backing pixmap** (`xvue/qt/src/xvue_qt_state.{h,cpp}`): `saved_canvas_` / `mempxaccro_` / `accroche_undo_tile_` — menu operations that trigger a redraw must respect the state's owned pixmaps.

### Established Patterns
- **`extern "C"` ABI wall** at `xvue/qt/src/xvue_qt_api.cpp` — Fortran sees byte-identical signatures from `xvue/qt/include/xvue_qt_api.h`. Menu bridge is **not** exposed through `extern "C"` (it's internal to Qt side); it interacts with Fortran only by feeding synthetic events into `xvsouris_`'s return value.
- **Warn-once stub pattern** for unimplemented entry points (kept for any xvue entry points Phase 6 doesn't wire — shouldn't be any if roadmap is accurate).
- **Per-module executable build** via `bin/cbmail`, `bin/cbelas`, etc. Each 6.1..6.5 modifies its own module's build path to link the shared shell + the module's menu declarations.
- **Bilingual flag via file existence**: `$MEFISTO/td/m/anglais` → English. No new mechanism needed; `XvueMenuBridge` reads this flag at construction and populates menu labels accordingly.
- **X11/Qt dual-build** (`bin/cbl_tout` for X11, `bin/cbl_tout_qt` for Qt — share `pp/` so never run in parallel). Phase 6 changes land on the Qt side; X11 backend stays unchanged.

### Integration Points
- **`XvueWindow` ctor** (`xvue/qt/src/xvue_qt_window.cpp`): After canvas is created and `installEventFilter(bridge)` is called, add `setMenuBar(new QMenuBar(this))`, `addToolBar(new QToolBar(this))`, `setStatusBar(new QStatusBar(this))`, `addDockWidget(Qt::BottomDockWidgetArea, new XvueConsoleDock(this))`. Module-specific menu additions happen in each `pp*_qt::main` between `XvueApp::ensure()` and event loop entry.
- **`XvueEventBridge::eventFilter` / `waitForEvent`**: Before returning normally, check if `XvueMenuBridge::pendingCommands_` is non-empty. If so, pop the head command and return it as a synthetic `notypeevent=2` keyboard event sequence (one character at a time, or wire-format encoded — planner decides). This is THE critical wiring point — UX-03.
- **`QProcess` pipe-reader** (UX-10): `XvueApp::ensure()` is a reasonable home for the stdout pipe setup IF the console is one-per-process. Per-module option: each `pp*_qt::main` creates its own `QProcess` (but since these are the Fortran solvers themselves running in-process, stdout is captured via `freopen(stdout, ...)` redirection to a pipe, not an external `QProcess`).
- **`QSettings` key namespace**: Use `QCoreApplication::setOrganizationName("LJLL")` + `setApplicationName("mefisto")` + per-module section (`QSettings::beginGroup("mail")`). Avoids key collisions across `pp*_qt` executables sharing the same ini file.

</code_context>

<specifics>
## Specific Ideas

- **"Feels like a proper desktop app"** — the user approved the move away from the typed-lexicon-only UX in the initial project scoping. Level-3 means real menus, real toolbars, real dialogs — not a cosmetic veneer over the same text prompt.
- **Bilingual must be frictionless** — existing MEFISTO users split between French and English. The `td/m/anglais` flag mechanism is load-bearing and must not break; FR is the default when the flag is absent (French origin of the project).
- **Keyboard muscle memory is sacred** — Phase 5 D-06 established that `Esc`=27 and `@`=64 must continue to reach Fortran as abort codes. `99;` / `5;90;1;` sequences must continue to work verbatim. Any shortcut scheme that intercepts plain characters is rejected.

</specifics>

<deferred>
## Deferred Ideas

- **V1X-07** (progress-bar fed by parsing iteration/residual lines from solver stdout) — Phase 9 or future; UX-10 ships the pipe-reader that V1X-07 would plug into later.
- **V1X-08** (animation scrubber for time-stepping solvers `flui`, `ther`, `nlse`) — Phase 9 or future.
- **V1X-09** (single-level snapshot undo for mesh edits) — Phase 9 or future.
- **Dialog for canvas colormap customization** — out of scope; scientific colormaps (`PALCDE`) stay managed via the typed lexicon per D-05.
- **Plugin / extension system for custom menu items** — out of scope; single-developer project doesn't need it.
- **Multi-document / multi-canvas support** — out of scope; each `pp*_qt` is single-window by design.

### Reviewed Todos (not folded)

None — no pending todos were surfaced for Phase 6 by `todo match-phase`.

</deferred>

---

*Phase: 06-level-3-ux-chrome-menu-surface*
*Context gathered: 2026-04-19*
