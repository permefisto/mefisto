# Feature Research

**Domain:** Qt 6 desktop GUI for a Fortran scientific computing package (interactive 2D/3D mesher + FE solvers + on-screen visualization)
**Researched:** 2026-04-10
**Confidence:** HIGH for Qt table stakes and anti-features (well-established Qt/desktop conventions); MEDIUM for differentiators (competitive reference points Gmsh / ParaView / Salome are authoritative but feature parity trade-offs depend on MEFISTO's specific workflows)
**Reference points:** Gmsh 4.x (mesher GUI), ParaView 5.x (post-processing), Salome 9.x (integrated FE workbench), FreeCAD FEM workbench

---

## Context recap (so the table below makes sense)

- **Core Value (non-negotiable):** every existing MEFISTO workflow keeps working; Fortran solver code in `mail/`/`elas/`/`flui/`/`ther/`/`nlse/`/`reso/`/`util/` stays bit-identical. The Fortran side must not notice the migration.
- **Graphics surface:** ~60 `extern "C"` entry points in `xvue/xvuelc.c` (drawing primitives, color/font, pixmap double-buffering, mouse events, PostScript export, interactive input loop). Those entry points MUST be preserved by name and signature.
- **Level 3 UX modernization:** the text-lexicon menu (typed commands like `99;` dispatched from a graphics-window text area) is replaced by a real `QMainWindow` + `QMenuBar` + `QToolBar` + `QDialog` surface — but the dispatch target is the SAME Fortran subroutines that the typed commands reach today.
- **"Maps to existing entry point"** below means: the feature is already expressed in `xvuelc.c`'s C surface and just gets a Qt backend. **"Net-new"** means: the feature exists only above the C surface (pure Qt chrome) or requires additional plumbing between Qt and the Fortran wrappers in `xvue/*.f`.

---

## Feature Landscape

### Table Stakes (Users Expect These from Any Modern Qt Desktop App in 2025/2026)

These are the features a modern desktop user simply assumes. They are cheap individually, together they set the "this is not a 1990s X11 window" bar.

| Feature | Why Expected | Complexity | xvuelc mapping | Notes |
|---------|--------------|------------|----------------|-------|
| `QMainWindow` shell with central widget | Any Qt app uses this pattern; gives dockable surroundings for free | S | Net-new | Replaces the raw `XCreateWindow` + text-input-area layout of the current `xvuelc.c`. Central widget is the canvas that hosts the Qt reimplementation of the Xlib drawing primitives. |
| Real `QMenuBar` with File / Edit / View / Mesh / Solve / Help | Level 3 modernization requirement from PROJECT.md; replaces text lexicon | M | Net-new (dispatches to existing Fortran subs) | The hard part is auditing the text lexicon per module (`mail/`, `elas/`, `flui/`, `ther/`, `nlse/`) and deciding which commands get a menu item vs. an action vs. a dialog. Dispatch target is unchanged. |
| `QToolBar` with icons for the most-used actions | Matches user expectation for a discoverable, clickable interface | S | Net-new | Same action objects as the menu bar (use `QAction` once, attach to both). Icons can start with Qt's built-in theme icons; bespoke icons are v1.x. |
| Keyboard shortcuts (`QKeySequence`) on every `QAction` | Power users expect Ctrl+S, Ctrl+Z, F1, Esc-to-cancel-picking, etc. | S | Net-new | Standard bindings (`QKeySequence::Save`, `Open`, `Quit`, `HelpContents`). Module-specific shortcuts (mesh refine, run solver) get custom bindings documented in Help. |
| Standard file dialogs (`QFileDialog`) for project open / save / export | Replaces typed paths in text lexicon; respects XDG recent-files | S | Partial (xvuelc has no file dialog — today it reads from `$MEFISTOX/<project>`) | Dialogs must still land the user in `$MEFISTOX/<project>/` because the Fortran side reads/writes files relative to that directory via its `util/` helpers — don't break that invariant. |
| About dialog (`QMessageBox::aboutQt` + custom about) | Table stakes; also the place to show MEFISTO's version, the Qt version, and credit Alain Perronnet / LJLL | S | Net-new | Include a link to `doc/normes.ps` or its rendering, and the MEFISTO version string. |
| Status bar (`QStatusBar`) with permanent + transient slots | Users expect progress, current mode, hover coordinates in the bottom bar | S | Net-new | Permanent slots: current project, current module, current language. Transient: mouse coordinates in mesh space, last typed command, last error. |
| HiDPI / fractional scaling support | Required on any modern laptop / 4K display | S | Free with Qt 6 | Qt 6 enables HiDPI by default; just make sure the custom-painted canvas uses `devicePixelRatioF()` for line widths and font metrics so meshes don't look like hairlines on 4K. |
| Dark-mode awareness (follow system palette) | Expected on GNOME / KDE with dark themes; Qt 6 picks it up automatically | S | Partial — needs care in `xvuelc` color mapping | The tricky part: MEFISTO's colormap (`xvStockeRGBtoColormap`, `xvCouleursImposees`) encodes physical meaning (stress magnitude, temperature, etc.). Scientific colors must NOT invert with the system theme; only chrome (menus, toolbars, dialog backgrounds) follows the palette. Document this as a design rule. |
| Clipboard integration for text fields and rendered views | `QApplication::clipboard()`; users expect Ctrl+C on any widget | S | Net-new | Copy-as-image of the current canvas view into the clipboard is a cheap win. |
| Window state persistence (`QSettings` — geometry, dock positions, last project) | Reopening the app to its previous state is an expected baseline | S | Net-new | Store under `~/.config/MEFISTO/` following the XDG base-dir spec. Keyed by module (MAILLER layout is different from ELASTICER layout). |
| Language selection at runtime (French / English) | Already exists via `$MEFISTO/td/m/anglais`; must be preserved | S | Maps to existing | Expose as a `Settings → Language` menu item that flips the same flag the launcher scripts check today. Must not break bilingual identifiers in error messages coming up from the Fortran side. |
| Print / export to PDF via `QPrinter` | Standard desktop expectation; replaces MEFISTO's EPS-only export | S | Partial — replaces `xvpostscript` entry point | `QPrinter::setOutputFormat(PdfFormat)` gives PDF for free. Keep the `xvpostscript` entry-point name/signature, swap the implementation to use `QPrinter`. |
| PNG / JPEG export via `QImage::save()` | Needed because `testa/` cases compare plots visually | S | Partial — `xvuelc.c` currently emits EPS only | Already in PROJECT.md Active requirements ("Qt takes over image export"). Output paths stay under `$MEFISTOX/<project>/` to match what the Fortran side expects. |
| Animated GIF export via `QMovie` / `QImageWriter` | Already exists (via `convert` shell-out); must be preserved after dropping ImageMagick | S | Maps to existing (replaces `bin/convertepsgif`) | Listed in PROJECT.md Active. `QImageWriter` supports GIF with loop count in Qt 6.10. |
| Undo framework (`QUndoStack`) for interactive mesh edits | A modern GUI without Ctrl+Z feels broken | **L** | **Net-new — heavy caveat** | The Fortran mesher is not designed around a command pattern — undo would require either (a) snapshotting mesh state between operations (cheap Qt-side, expensive in memory for large meshes) or (b) reverse operations per mesh command (requires Fortran-side changes → **violates Core Value**). **Recommendation: defer to v1.x**, start with a single-level "undo last mesh modification" using a snapshot-of-the-scratch-files approach that needs no Fortran changes. Flag this in PITFALLS. |
| Progress reporting for long solver runs (`QProgressDialog` or status-bar `QProgressBar`) | Users expect feedback on long operations | M | Net-new | The Fortran solvers today print progress to stdout. Option A: parse stdout from the solver-driving goroutine in C++. Option B: add a tiny `extern "C" void xvprogress(int pct)` entry point callable from Fortran — BUT that requires edits to solver files → **out of scope**. **Recommendation: Option A** (pipe-reader in the Qt side, parses the existing `WRITE(*,*)` diagnostics). Fragile but zero-Fortran-change. |
| Cancel-long-operation button | Paired with progress dialog; expected when a progress bar is visible | M | Net-new with caveats | True cancellation requires cooperative checks in the Fortran solver, which the solver doesn't have. Realistic v1: "Interrupt" sends SIGINT to the solver process — but MEFISTO treats Ctrl+C as dangerous (`CLAUDE.md`: "never Ctrl-C, leaves project files inconsistent"). **Recommendation: omit in v1, document as a known limitation**, add a scary confirmation dialog for "Force stop" in v1.x. |
| Context menus on canvas (`QContextMenuEvent`) | Right-click is expected to show a menu | S | Net-new | Context-sensitive: different menu when right-clicking a node vs. an edge vs. empty space. Uses the same `QAction` objects as the main menu bar. |
| Tooltips on all actions and toolbar buttons | `QAction::setToolTip` — one-line per action | S | Net-new | Bilingual: tooltip strings come from the same resource bundle used for the `anglais` language flag. |
| "What's this?" help (`Qt::WhatsThisRole`) on dialog fields | Modern desktop convention for field-level help | S | Net-new | Optional — can be v1.x if time-boxed. |
| Error dialogs with actionable messages (`QMessageBox::critical`) | Replaces `WRITE(*,*)` diagnostics dumped to the terminal | M | Partial | Requires capturing stdout/stderr from the Fortran side (same pipe-reader used for progress). Parse `ERREUR`/`ERROR` lines, surface as Qt dialogs. Keep terminal output as a fallback log. |
| Log / console dock widget | Replaces the terminal scrollback; expected in any modern dev-tool-like app | M | Net-new | Dockable `QDockWidget` containing a `QPlainTextEdit` fed by the stdout-pipe reader. Searchable. User can detach or hide. Doubles as the "console" for any residual text-lexicon commands the user wants to type by hand during the transition. |
| Recent projects list (`File → Open Recent`) | Standard desktop convention | S | Net-new | Stored in `QSettings`. Entries point to directories under `$MEFISTOX/`. |
| Drag-and-drop project directories onto the window | Bonus expectation on modern desktops; cheap with `QMimeData` | S | Net-new | v1.x candidate, not v1. |
| Help → User Manual opens `doc/` HTML or PDF | Users expect F1 to open documentation | S | Net-new | MEFISTO has `doc/`, `doca/`, `docf/` (PostScript/HTML bilingual). F1 launches the system PDF viewer on the right file based on current language. |
| Crash / unexpected-exit detection on startup | "Recover last project?" prompt on startup is standard in engineering tools | M | Net-new | Implemented by writing a lockfile in `$MEFISTOX/<project>/` on load and clearing it on clean exit. On next launch, if lockfile exists, offer to recover from the scratch files. Zero Fortran-side changes. Bonus: partially mitigates the "STOP loses project state" concern from CONCERNS.md. |

**Table-stakes total rough cost estimate:** ~60% of the Qt migration effort is just table stakes chrome. None of it is individually hard; the cost is in the breadth and in the per-module lexicon audit.

---

### Differentiators (Scientific-FE-Specific Features That Set MEFISTO Qt Apart From the X11 Version)

These are features that the current `xvue/` cannot do well or at all, and where Qt 6 (via `QPainter`, `QOpenGLWidget`, `QGraphicsView`) enables a real UX improvement. They map to what Gmsh / ParaView / Salome users already expect from a modern FE GUI.

| Feature | Value Proposition | Complexity | xvuelc mapping | Notes |
|---------|-------------------|------------|----------------|-------|
| Smooth, anti-aliased 2D rendering via `QPainter` with `Antialiasing` hint | Eliminates the pixelated look of the X11 build immediately | S | Maps to existing drawing primitives | Free win: set `QPainter::Antialiasing` when reimplementing the primitives. Noticeable to users the first time they see a mesh. |
| Pan / zoom / rotate with mouse wheel + middle-drag (canvas navigation) | Gmsh/ParaView-standard interaction; X11 build has clunky keyboard-driven navigation | M | Net-new on Qt side; passes the resulting view matrix to existing drawing primitives | Keep the Fortran view-state unchanged — transform happens in `QPainter` via `QTransform`. No solver-side change. |
| Coordinate readout in status bar (mouse-over mesh-space coordinates) | Modern expectation in any CAD/mesh tool; replaces having to type `?` at the text lexicon | S | Net-new | Status bar slot; updates on `mouseMoveEvent`. Coordinates are computed by inverting the current `QTransform`. |
| Interactive picking feedback (hover highlights, selection rubber-band) | Makes mesh editing discoverable; the X11 build has hit-testing but no hover feedback | M | Partial — hit-testing exists via `xvsouris`/`xvsouris2`/`deplsouris` entry points; hover visualization is net-new | The existing entry points return an index; the Qt side can decorate that index (outline the triangle, show its ID in the status bar) without touching the Fortran. |
| Mesh-quality color overlay (aspect ratio, skew, Jacobian) | Directly addresses mesher heuristics flagged in CONCERNS.md (`mail/a1teqm.f`) | M | Net-new color mapping on top of existing mesh drawing | Quality metric values must come from the Fortran side (existing routines in `mail/` compute them). Qt maps to colormap and renders. If the existing Fortran subs don't expose the metric per-element, this is v1.x. |
| Solution result color bars with log-scale toggle | Scientific visualization baseline; X11 build has linear-only | M | Net-new (`QWidget` colorbar next to the canvas) | Log-scale toggle is a `QAction` that remaps the color LUT Qt-side. No Fortran change. |
| Multiple viewports / split view (`QSplitter` with two canvas widgets) | Compare mesh vs. solution, or two solver outputs side-by-side | M | Net-new | Each canvas is a separate instance of the Qt reimplementation of the X11 drawing surface. Requires xvuelc state to be instance-scoped, not global — **check `xvuelc.c` for globals**, this is a known risk. |
| Solution-over-time animation scrubber (`QSlider` + Play button) | Standard ParaView feature; MEFISTO has time-stepping in `flui/`, `ther/`, `nlse/` already writing frames to scratch files | M | Net-new Qt chrome on top of existing frame-reading routines | v1.x candidate. Fortran already writes per-timestep output; the scrubber is pure Qt that re-invokes the existing display routines for the selected frame index. |
| Clip planes / section cuts for 3D views | Table stakes for 3D FE visualization in ParaView/Salome; differentiator for MEFISTO | **L** | Net-new | Only meaningful if the rendering path goes via `QOpenGLWidget`. With pure `QPainter` 2D projection (what `xvuelc.c` effectively does today), clip planes require rewriting the display path. **Defer to v2.** |
| Command log / history panel showing the lexicon commands equivalent to each GUI action | Killer feature for transition: users learning the Qt UI can see what the old text command would have been | M | Net-new | Doubles as documentation and builds trust with long-time users. Each `QAction::triggered` logs its equivalent lexicon command to the console dock. |
| Replay-command / scripted-session capture | Record the command log, save as a `.mef` script, replay on another project | M | Net-new | Replay just feeds the saved commands into the same dispatch layer that the GUI actions use. No Fortran change because the text lexicon already parses these commands. Note: this is NOT a full scripting language — see anti-features. |
| Dark-mode-safe scientific colormaps (viridis, plasma, turbo) | ParaView-standard; replaces the MEFISTO rainbow colormap which is known to distort perception | S | Maps to existing `xvStockeRGBtoColormap`/`xvCouleursImposees` — just add new palette presets | Adds a `View → Colormap → {Rainbow / Viridis / Plasma / Grayscale}` menu. Existing rainbow stays as the default for `testa/` A/B comparison. |
| Mesh statistics panel (element count, quality histogram, min/max edge length) | Every modern mesher has this; MEFISTO computes the values internally but doesn't surface them nicely | M | Partial — stats exist in `mail/` already, wrapped in text output | Dockable `QDockWidget` with a `QTableView` + a tiny `QChartView` histogram. Reads values by parsing the existing Fortran text output OR by adding a read-only query entry point (if the Fortran already has one in `util/`). |
| Per-module workspace layouts (`QMainWindow::saveState`) | A user in FLUIDER wants different dock visibility than in MAILLER | S | Net-new | Free once the docks exist; Qt persists per-module layouts in `QSettings`. |
| Bilingual UI with live-switch (no restart) | The current launcher-script language flag requires restart | S | Net-new | Store strings in `QTranslator` files generated from the existing `td/m/anglais` vs. French resource. Low user priority but removes a long-standing papercut. |

---

### Anti-Features (Commonly Requested for "Modern" GUIs, Explicitly Out of Scope Here)

These are features that sound modern but either (a) violate the Core Value by requiring Fortran-side changes, (b) balloon the scope beyond what a single developer can finish, or (c) are explicitly deferred in PROJECT.md's Out of Scope section.

| Feature | Why Requested | Why Problematic for xvue-qt | Alternative |
|---------|---------------|-----------------------------|-------------|
| **Custom theming engine / skin system** | "Make it look cool" | Adds a whole styling subsystem; Qt already follows the system palette for free; blocks dark-mode-follow-system; maintenance burden forever | Use system `QStyle` + `QPalette`. Dark mode is inherited from the OS. Zero code. |
| **Plugin system (dynamic-loaded modules)** | "Let users extend it" | No existing use case; MEFISTO has zero extension points in Fortran; would require a C++ plugin ABI that drifts with Qt versions; single developer cannot maintain an ABI | If extension is ever needed, users modify the Fortran sources and rebuild — consistent with the project's tradition. |
| **Qt Quick / QML UI** | "QML is the modern way" | QML is the modern way for *mobile/embedded/touch*. For a dense desktop engineering app with dozens of menu commands, **Qt Widgets is the correct choice** — better accessibility, richer widgets, proven on Gmsh / FreeCAD / QGIS / KDE. QML would force custom implementations of trees, tables, docks that Widgets provides free. | Use Qt Widgets (`QMainWindow`, `QDockWidget`, `QTableView`, `QTreeView`, `QMenuBar`). This is already in PROJECT.md constraints (Qt 6 Widgets from Debian apt). |
| **Embedded scripting language (Python / Lua / QtScript)** | "Let users automate workflows" | Requires a stable C++ binding layer *and* Fortran-callable bridges; QtScript is deprecated; Python binding (PyQt/PySide) would triple the dependency surface; single-developer maintenance is unrealistic | Command log / replay (listed as a differentiator) provides "scripted" workflows using the existing text lexicon — **no new language, no new runtime**. |
| **Cloud sync / remote project storage** | "Modern apps sync to the cloud" | MEFISTO projects are directory trees with scratch files holding paged mesh/matrix data — the format is implicit in the Fortran code (CONCERNS.md). Any cloud sync has to re-upload multi-GB scratch files on every save. Zero user demand. Out of scope in PROJECT.md. | Users can rsync `$MEFISTOX/` themselves. Not our problem. |
| **Real-time collaborative editing** | "Google-Docs for meshes" | Requires a CRDT layer over a Fortran program that holds state in `COMMON` blocks; effectively impossible without rewriting the solver core — violates Core Value. | N/A. Not a use case for a single-user scientific tool. |
| **Auto-updater / in-app update** | "Users should get fixes automatically" | MEFISTO is distributed as source; there are no binary releases; the build embeds the install path via `incl/homdir.inc` generated at build time (CONCERNS.md) — auto-replacing binaries would break that. | Rely on distro packaging (Debian apt) or `git pull && bin/cbl_tout`. |
| **Crash-reporter uploading dumps to a server** | "Breakpad is modern" | Requires a server, privacy policy, data-handling policy, maintenance. Single developer. Zero budget. Users are scientists, not telemetry donors. | Local crash dump: on unexpected exit, write the console-log dock contents to `$MEFISTOX/<project>/crash.log` and tell the user where to find it. |
| **Telemetry / analytics** | "We need to know how users use it" | Same privacy and maintenance concerns as crash-reporter uploading. Scientific users will refuse. | None. Do not collect data. |
| **Integrated code editor for .f files** | "Edit Fortran inside the GUI" | Orthogonal to the graphics layer. Users have `vim`/`emacs`/`vscode`. Would balloon scope by an order of magnitude. | External editor. |
| **Built-in CAD / geometry kernel (OpenCASCADE)** | "Salome has it" | Adds a ~100 MB dependency, a whole second API surface, and a second build system. MEFISTO's existing geometry input (typed points / lines / surfaces in the MAILLER lexicon) is the scope. | Keep MAILLER's existing geometry definition as-is. If users want OpenCASCADE, they can use Salome. |
| **Python-based post-processing (like ParaView-Python)** | "ParaView has it" | Python bridge is a giant dependency and a new FFI surface. No demand in MEFISTO's user base. | Export to VTK/CSV via a new `QAction` (see v1.x candidates) if users want to post-process in ParaView. |
| **Full retheming of `bin/INITIER`, `bin/MAILLER`, etc., launcher scripts as Qt dialogs** | "Everything should be Qt" | Explicitly Out of Scope in PROJECT.md — "Only the graphics windows they spawn become Qt". | Launcher scripts stay bash. |
| **Rewriting the existing text lexicon parser in C++** | "Clean up legacy" | The lexicon lives in `util/` and is called from every solver module. Removing it touches thousands of call sites — violates Core Value. | Keep the lexicon alive as the dispatch backbone; `QAction::triggered` calls the same subroutines the lexicon would have called. The console dock also lets users type lexicon commands directly during transition. |
| **Plugin-based colormap system** | "Users want custom colormaps" | Over-engineered for the actual requirement (4 presets: rainbow, viridis, plasma, grayscale). | Hard-code the 4 presets. Done. |
| **Multi-document interface (MDI) with multiple projects open simultaneously** | "Modern apps support tabs" | Each MEFISTO executable is one-project-per-process by design (launchers `cd` into the project dir first, Fortran reads/writes relative paths). MDI would require making the Fortran side project-aware *per-window* — deep change, violates Core Value. | One project per window, launch a second window for a second project. Same as today. |
| **Integrated test runner for `testa/` / `testf/`** | "Run tests from the GUI" | `testa/` cases are interactive; they require human eyeballing per PROJECT.md and TESTING.md. A "run all tests" button is misleading. | Out of scope. Tests stay manual. |

---

## Feature Dependencies

```
QMainWindow shell
    |
    +--requires--> Qt reimplementation of xvuelc.c entry points (from STACK.md)
    |                   |
    |                   +--requires--> xvuelc state made instance-scoped (not global)
    |
    +--enables----> QMenuBar / QToolBar / QAction framework
    |                   |
    |                   +--requires--> per-module lexicon audit (mail, elas, flui, ther, nlse)
    |                   |
    |                   +--enables----> Keyboard shortcuts
    |                   +--enables----> Context menus
    |                   +--enables----> Command log / replay
    |
    +--enables----> QStatusBar
    |                   +--enables----> Coordinate readout
    |                   +--enables----> Progress reporting
    |
    +--enables----> QDockWidget framework
                        +--enables----> Console / log dock
                        |                   +--requires--> stdout pipe-reader from Fortran
                        |                   +--enables----> Error dialogs (parsed from stdout)
                        +--enables----> Mesh statistics panel
                        +--enables----> Solution animation scrubber
                        +--enables----> Workspace layout persistence

QPainter with Antialiasing
    +--enables----> Pan/zoom/rotate via QTransform
    |                   +--enables----> Coordinate readout in status bar
    |                   +--enables----> Interactive picking feedback
    |                   +--requires--> Existing xvsouris/xvsouris2 entry points preserved
    +--enables----> Mesh-quality color overlay
    +--enables----> Scientific colormaps (viridis/plasma)
    +--enables----> Color bars with log-scale toggle

QFileDialog + QSettings
    +--enables----> Recent projects list
    +--enables----> Window state persistence
    +--enables----> Per-module layout persistence

QPrinter (PDF)   ┐
QImage           ├--replaces--> xvpostscript entry point
QImageWriter/GIF ┘                    (implementation only — name preserved)

Undo framework (QUndoStack)
    +--conflicts-with--> Core Value (if true undo is pursued)
    +--alternative----> Scratch-file snapshot before mesh-modifying actions (no Fortran change)
```

### Dependency Notes

- **`xvuelc.c` state must become instance-scoped before split-view / multi-viewport work.** Today's C file likely has module-level globals (drawing context, current color, current font). The Qt reimplementation should encapsulate these in a C++ class, with the `extern "C"` entry points dispatching to a current-instance pointer. **This is a v1 prerequisite, not a feature** — surfaces in STACK.md, called out here because it blocks the multi-viewport differentiator.

- **Stdout pipe-reader is a keystone.** It enables progress reporting, error dialogs, and the console/log dock — three features that together cover most of the "this feels like a modern app" impression. One implementation unlocks three features; build it first.

- **Per-module lexicon audit is the long pole.** The menu-bar / toolbar / dialog work is only as good as the audit that enumerates which lexicon commands deserve which UI surface. This must be scheduled as its own phase (per PROJECT.md's Level 3 modernization), and should probably be broken into one sub-phase per solver module (`mail`, `elas`, `flui`, `ther`, `nlse`).

- **Undo conflicts with Core Value if pursued naively.** The only no-Fortran-change undo path is scratch-file snapshotting, which is space-expensive. Mark as v1.x and document the approach in PITFALLS.

- **Animation scrubber depends on existing time-stepping output.** `flui`, `ther`, `nlse` already emit per-timestep results (confirm during phase 1). The scrubber is pure Qt chrome around whatever frame-selection mechanism the text lexicon exposes today. Cheap if the pre-req holds; expensive if the Fortran side currently only re-runs from scratch.

- **Command log enables replay enables light "scripting" — without adding a scripting language.** This is the answer to anyone who asks for embedded Python: they already have it, it's called the text lexicon, and the console dock exposes it.

---

## MVP Definition

The MVP answer must honor the Core Value: **every workflow that works today keeps working**, plus the minimum Qt chrome that makes the result feel like a 2026 desktop app rather than an X11 window with a Qt skin.

### Launch With (v1 — the xvue-qt migration proper)

Must-haves, grouped by dependency order.

**Backend (no user-visible chrome, but gates everything):**
- [ ] Qt reimplementation of all ~60 `extern "C"` entry points in `xvuelc.c` (PROJECT.md Active) — **P1, existing entry points, complexity L in aggregate, S per entry point**
- [ ] xvuelc state encapsulated in a C++ class, instance-scoped — **P1, enables multi-viewport later, S**
- [ ] Anti-aliased 2D rendering via `QPainter` — free with the above — **P1, S**
- [ ] Stdout pipe-reader for solver process — **P1, M, keystone**

**Chrome (the Level 3 modernization):**
- [ ] `QMainWindow` shell with central canvas widget — **P1, S**
- [ ] `QMenuBar` with File / Edit / View / Mesh / Solve / Help — **P1, M (per-module lexicon audit is the cost)**
- [ ] `QToolBar` sharing `QAction` objects with the menu bar — **P1, S**
- [ ] Keyboard shortcuts on every `QAction` — **P1, S**
- [ ] `QFileDialog` for project open / save / export — **P1, S**
- [ ] `QStatusBar` with permanent + transient slots, including mouse coordinates — **P1, S**
- [ ] About dialog with MEFISTO version, Qt version, Alain Perronnet / LJLL credit — **P1, S**
- [ ] Context menus on canvas — **P1, S**
- [ ] Tooltips on all actions and toolbar buttons — **P1, S**
- [ ] Error dialogs parsed from Fortran stdout via the pipe-reader — **P1, M**
- [ ] Console / log dock widget (hideable) — **P1, M**
- [ ] Recent projects list in File menu — **P1, S**
- [ ] Window state persistence via `QSettings` — **P1, S**
- [ ] Language selection at runtime (preserves existing French/English split) — **P1, S**
- [ ] Help → User Manual (opens `doc/` PDF in system viewer) — **P1, S**
- [ ] HiDPI handling on the canvas (`devicePixelRatioF`) — **P1, S**
- [ ] Dark-mode aware chrome; scientific colors frozen regardless of theme — **P1, S**

**Export replacements (already in PROJECT.md Active):**
- [ ] PNG / JPEG export via `QImage::save()` — **P1, S, replaces xvuelc EPS-only**
- [ ] PDF export via `QPrinter` — **P1, S, reuses xvpostscript entry point name**
- [ ] Animated GIF via `QMovie` / `QImageWriter`, dropping ImageMagick dep — **P1, S**

**Navigation (cheap differentiator that users will notice immediately):**
- [ ] Pan / zoom / rotate via mouse wheel and middle-drag — **P1, M**
- [ ] Coordinate readout in status bar — **P1, S (depends on navigation)**
- [ ] Interactive picking feedback (hover highlight) — **P1, M (reuses existing xvsouris entry points)**

**Validation:**
- [ ] Visual A/B on the `testa/` subset (PROJECT.md Active) — **P1, process not feature**

### Add After Validation (v1.x — deferred but planned)

- [ ] Crash-recovery lockfile + startup "Recover last project?" prompt — **P2, M, mitigates CONCERNS.md STOP issue**
- [ ] Command log dock showing the lexicon-equivalent of each GUI action — **P2, M, huge for user trust**
- [ ] Command replay / save-and-replay session — **P2, M, depends on command log**
- [ ] Scientific colormap presets (viridis / plasma / grayscale), rainbow kept as default for A/B — **P2, S**
- [ ] Color bar widget with log-scale toggle — **P2, M**
- [ ] Mesh statistics dock (element count, min/max edge length, quality histogram) — **P2, M**
- [ ] Per-module saved workspace layouts — **P2, S**
- [ ] Drag-and-drop project directories onto the main window — **P2, S**
- [ ] Context-sensitive right-click menus differentiated by clicked entity — **P2, S**
- [ ] Progress reporting via `QProgressDialog` parsed from solver stdout — **P2, M**
- [ ] Single-level snapshot-based undo for mesh edits — **P2, L, scratch-file snapshot approach**
- [ ] Solution animation scrubber for time-stepping solvers (`flui`, `ther`, `nlse`) — **P2, M, depends on frame-selection mechanism in current lexicon**

### Future Consideration (v2+)

- [ ] Multiple viewports / split view — depends on instance-scoped `xvuelc` state being rock-solid — **P3, M**
- [ ] Mesh-quality color overlay — depends on per-element quality metrics being exposed from `mail/` without Fortran changes — **P3, M**
- [ ] "Force stop" interrupt for runaway solver runs (with scary confirmation) — **P3, M**
- [ ] Live-switch bilingual UI without restart via `QTranslator` — **P3, S**
- [ ] VTK / CSV export for ParaView post-processing — **P3, S**
- [ ] 3D clip planes / section cuts — requires `QOpenGLWidget` rewrite of the display path — **P3, L, only if users ask**
- [ ] Multi-level undo — only if an alternative to scratch-file snapshots emerges — **P3, L**
- [ ] True cooperative solver cancellation — requires Fortran-side changes, explicitly out of Core Value — **P3, never unless Core Value is renegotiated**

---

## Feature Prioritization Matrix

| Feature | User Value | Implementation Cost | Priority |
|---------|------------|---------------------|----------|
| Qt reimplementation of `xvuelc.c` entry points | HIGH (mandatory) | HIGH | **P1** |
| `QMainWindow` + `QMenuBar` + `QToolBar` + `QAction` | HIGH | MEDIUM | **P1** |
| Stdout pipe-reader for Fortran process | HIGH (enables 3 features) | MEDIUM | **P1** |
| `QFileDialog` / recent projects / `QSettings` | MEDIUM | LOW | **P1** |
| Console / log dock | HIGH | MEDIUM | **P1** |
| PNG/PDF/GIF export via Qt | HIGH (replaces ImageMagick) | LOW | **P1** |
| Pan / zoom / rotate navigation | HIGH | MEDIUM | **P1** |
| Interactive picking feedback | HIGH | MEDIUM | **P1** |
| HiDPI + dark-mode chrome | MEDIUM | LOW | **P1** |
| Context menus | MEDIUM | LOW | **P1** |
| Error dialogs from parsed stdout | HIGH | MEDIUM | **P1** |
| Crash-recovery lockfile | MEDIUM | MEDIUM | **P2** |
| Command log / replay | HIGH | MEDIUM | **P2** |
| Scientific colormaps (viridis/plasma) | MEDIUM | LOW | **P2** |
| Color bar + log-scale toggle | MEDIUM | MEDIUM | **P2** |
| Mesh statistics dock | MEDIUM | MEDIUM | **P2** |
| Progress reporting | MEDIUM | MEDIUM | **P2** |
| Snapshot-based undo (single-level) | MEDIUM | HIGH | **P2** |
| Animation scrubber for time-stepping | MEDIUM | MEDIUM | **P2** |
| Multiple viewports / split view | LOW | MEDIUM | **P3** |
| Mesh-quality color overlay | MEDIUM | MEDIUM | **P3** |
| 3D clip planes | LOW (unless 3D use grows) | HIGH | **P3** |
| Live-switch bilingual UI | LOW | LOW | **P3** |
| Custom theming / plugins / scripting language | LOW | HIGH | **anti-feature** |
| Cloud sync / collaboration / telemetry | LOW | HIGH | **anti-feature** |

**Priority key:**
- P1: Must have for xvue-qt v1 landing (the Qt migration with Level 3 modernization)
- P2: Should have in v1.x (post-A/B validation)
- P3: Future consideration (v2+ or never)

---

## Competitor Feature Analysis

| Feature | Gmsh 4.x | ParaView 5.x | Salome 9.x | MEFISTO Qt (our approach) |
|---------|----------|--------------|------------|---------------------------|
| GUI framework | FLTK | Qt Widgets + VTK | Qt Widgets + PyQt + OpenCASCADE | **Qt 6 Widgets** (matches ParaView / Salome) |
| Menu bar + toolbars | yes | yes | yes | **yes (v1, P1)** |
| Embedded scripting | Gmsh scripting language (.geo) | Python + Trace | Python + YACS | **No — command log / replay of existing text lexicon instead** (anti-feature) |
| Pan / zoom / rotate | yes | yes | yes | **yes (v1, P1)** |
| Picking feedback | yes | yes | yes | **yes (v1, P1)** — reuses existing xvsouris |
| Mesh quality visualization | yes (element quality view) | via filter | yes | **v2 (P3)** — depends on Fortran-side metric exposure |
| Solution animation | limited | yes (scrubber + playback) | yes | **v1.x (P2)** — scrubber only, no playback speed control in v1 |
| Color bars | yes | yes | yes | **v1.x (P2)** |
| Scientific colormaps | limited | full (viridis, plasma, turbo, etc.) | inherits ParaView | **v1.x (P2)** — 4 presets hard-coded |
| Clip planes (3D) | yes | yes | yes | **v2+ (P3)** — requires `QOpenGLWidget` rewrite |
| Multiple viewports | yes | yes | yes | **v2 (P3)** — requires instance-scoped xvuelc state |
| Undo/redo | limited | yes (for pipeline) | yes | **v1.x (P2)** — snapshot-based, single-level |
| Dark mode | partial | yes | yes | **v1 (P1)** — chrome only, scientific colors frozen |
| HiDPI | yes | yes | yes | **v1 (P1)** — free with Qt 6 |
| Command log | no | yes (trace) | yes (Python trace) | **v1.x (P2)** — lexicon-command log |
| File export PDF | EPS | yes | yes | **v1 (P1)** — `QPrinter` |
| Image export PNG | yes | yes | yes | **v1 (P1)** — `QImage::save()` |
| Animated GIF | no | yes | yes | **v1 (P1)** — `QImageWriter` loop |
| Plugin system | limited | yes | yes (Python) | **No** (anti-feature) |

**Takeaway:** MEFISTO Qt should aim squarely at **ParaView-class chrome** (menu bar, toolbar, docks, status bar, HiDPI, dark mode, `QImage`-based export, scientific colormaps) while explicitly rejecting **ParaView/Salome's scripting language ambitions** because those require Fortran-side bridges that violate the Core Value. The closest "spiritual" peer for scope is actually **Gmsh** — a focused mesh+solver GUI with strong interactive editing — which is why Gmsh's navigation and picking UX are the best reference points for MEFISTO's v1.

---

## Sources

- PROJECT.md (in-repo) — Core Value, constraints, Out of Scope section, Level 3 modernization target
- .planning/codebase/ARCHITECTURE.md (in-repo) — xvue/ layer, 60 extern "C" entry points, lexicon dispatch pattern
- .planning/codebase/CONCERNS.md (in-repo) — STOP on fatal errors losing project state, X11 fragility, no structured logging, GIF export via ImageMagick shell-out
- CLAUDE.md (in-repo) — working rules, `99;` save-exit convention, Qt migration as active goal
- Gmsh project (https://gmsh.info/) — reference for mesher GUI interaction patterns
- ParaView project (https://www.paraview.org/) — reference for scientific visualization chrome (color bars, animation scrubber, clip planes)
- Salome platform (https://www.salome-platform.org/) — reference for integrated Qt-based FE workbench
- FreeCAD FEM workbench — reference for Qt Widgets + Fortran/C++ backend integration in a desktop engineering app
- Qt 6.10 documentation (https://doc.qt.io/qt-6/) — `QMainWindow`, `QAction`, `QPainter`, `QImageWriter`, `QPrinter`, `QSettings`, `QUndoStack`, `QTranslator`, HiDPI handling, dark-mode / palette following

**Confidence calibration:**
- **HIGH** for Qt Widgets table-stakes features and anti-features (these are well-established desktop conventions and Qt 6 provides them directly).
- **HIGH** for the anti-features list (each is grounded in an explicit PROJECT.md constraint or a Core Value violation).
- **MEDIUM** for the differentiators and the competitor-feature comparison — the reference points (Gmsh, ParaView, Salome) are authoritative, but feature parity/trade-off judgments depend on MEFISTO workflows that will only be fully understood during the per-module lexicon audit.
- **MEDIUM** for undo / progress reporting / cancellation feasibility — the "no Fortran changes" constraint forces workarounds (scratch-file snapshots, stdout pipe-reader) whose practical cost will only be known when prototyped.

---

*Feature research for: Qt 6 GUI replacing MEFISTO's X11/Motif graphics layer*
*Researched: 2026-04-10*
