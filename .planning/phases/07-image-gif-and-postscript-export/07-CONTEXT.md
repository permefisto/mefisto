# Phase 7: Image, GIF, and PostScript export - Context

**Gathered:** 2026-05-01
**Status:** Ready for planning

<domain>
## Phase Boundary

Replace the X11/ImageMagick export stack with a Qt-native one:

- **PNG / JPEG single-frame snapshots** of the current canvas via `QImageWriter` / `QImage::save`, triggered from the Qt File menu only.
- **Animated GIF** of a frame sequence via probe-driven dispatch (`QImageWriter` GIF plugin if `EXPORT-01` PROBE.md confirms support; otherwise PNG-sequence + `ffmpeg` via `QProcess`). Captured automatically by intercepting `xvpostscript_` save points so testa/wave and testa/cavity2d need zero solver changes.
- **PostScript** via the hand-rolled `~120-line xvpostscript_` body ported verbatim into a new `PsEmitter` class in `xvue/qt/src/xvue_qt_postscript.{h,cpp}`, with the ~150 per-primitive `if(lasopsc>0) fprintf(fpo,...)` branches scattered through `xvuelc.c` reproduced as byte-identical helper calls.
- **PDF** as an additive bonus — raster PDF via `QPdfWriter::drawPixmap(backing_)` from a File menu action.

**Invariants:**

- Fortran ABI stays at 58 symbols. `xvpostscript_(int *lasops)` remains the ONLY Fortran-driven export entry point ("Fortran must not notice"). PNG/JPEG/GIF/PDF are user-triggered Qt menu actions (or env-var-toggled auto-capture for animation).
- Phase 7 removes ImageMagick **from the Qt code path only**. `bin/convertepsgif` stays alive for the X11 backend during the one-release-cycle A/B window. Full removal across `bin/`, `td/`, `testa/`, `testf/` is RETIRE-03 in Phase 9.
- PostScript Y-flip (`ypixels - y`) preserved verbatim at emit time per `xvue/README_COORDS.md` Phase 7 row.
- Single-backing pixmap (Phase 2 D-04, Phase 4 D-04): `backing_` is the snapshot source for every raster export.
- Bilingual FR/EN status / dialog text (Phase 6.0 D-09 pattern).

</domain>

<decisions>
## Implementation Decisions

### Export ABI & trigger surface

- **D-01:** **No new Fortran ABI for image/GIF/PDF.** PNG, JPEG, GIF, and PDF are triggered exclusively via the Qt File menu (or, for animation, via env-var/menu toggle that Qt-side intercepts). `xvpostscript_(int *lasops)` stays the ONLY Fortran-callable export entry point. ABI symbol count unchanged at 58. Honors the project-wide "Fortran must not notice the migration" invariant.
- **D-02:** **Auto-snapshot animation capture.** When `XVUE_ANIM=1` env var is set OR the File → "Capture animation" menu toggle is active, `PsEmitter` snapshots `backing_` to an in-memory `QImage` list on every `xvpostscript_(*lasops==0)` close (the existing "PS save complete" trigger that testa/wave and testa/cavity2d already emit between frames). The frame list is flushed to `animation.gif` on the next `xvpostscript_(0)` close after the toggle is disabled, or via explicit File → "End animation" menu. Zero changes to solver/driver code; reuses existing PS-save points in `testa/`.
- **D-03:** **Hybrid filename.** Menu-triggered exports prompt via `QFileDialog::getSaveFileName` with `QSettings("xvue/export/last_dir")` remembered across runs. Auto-snapshot frames during animation capture use fixed legacy-compatible names (`zfxy0001.png`, `zfxy0002.png`, …) in `cwd` — mirrors `bin/convertepsgif`'s `zfxy0*.eps` glob convention so external scripts that watch this prefix still work. Final GIF name: `animation.gif` in `cwd` (overwrites without warning, matching legacy convertepsgif behavior).
- **D-04:** **File menu wired in shared shell.** `File → Export → {PNG…, JPEG…, GIF…, PDF…}` entries added once in the shared `xvue_qt_app` File-menu builder (Phase 6.0 module). All five `pp*_qt` executables inherit the entries automatically. No per-module audit/wiring needed in 6.1..6.5 follow-ups; no changes to per-module `register<Mod>Actions_stub_` keepalives.

### PostScript port granularity

- **D-05:** **EXPORT-04 'verbatim ~120 lines' = `xvpostscript_` function body only.** Lines 1187–1304 of `xvuelc.c` (open/close `TEMPORAIRE.EPS`, `lasops` dispatch, `concat` flush, mode 100/101/102 erase) are ported byte-for-byte into `PsEmitter::handleLasops(int lasops)`. The ~150 in-primitive `if(lasopsc>0) fprintf(fpo,...)` branches scattered through `xvchargefonte_`, `xvtrait_`, `xvtraitcouleur_`, `xvface_`, `xvfaceisocouleur_`, `xvflpt_`, `xvellipse_`, `xvfond_`, `effacer`, etc., are PS-emit helpers — ported as separate methods on `PsEmitter` with the SAME format strings so output bytes match.
- **D-06:** **`PsEmitter` helper class.** New `xvue/qt/src/xvue_qt_postscript.{h,cpp}`. Public methods mirror xvuelc.c emit sites: `PsEmitter::line(x1,y1,x2,y2,rgb,epais)`, `PsEmitter::face(coords[],n,col_pol,col_fac)`, `PsEmitter::text(s,x,y,rgb,font)`, `PsEmitter::loadFont(qfont,ascent,width,mp)`, `PsEmitter::ellipse(...)`, `PsEmitter::flpt(...)`, `PsEmitter::fond(rgb)`, `PsEmitter::clear()`. Each method's body uses the EXACT format strings from the corresponding `xvuelc.c` `fprintf(fpo,...)` site (preserves byte-identical PS output). README_COORDS.md `ypixels - y` flip applied inside each method, never in the caller. Qt primitives in `xvue/qt/src/xvue_qt_api.cpp` call the matching `PsEmitter::*` after their `QPainter::*` calls (when `psEmitter().active()`).
- **D-07:** **PS state encapsulated in `PsEmitter` instance.** All ~15 file-statics from `xvuelc.c` (`lasopsc`, `modepsc`, `fpo`, `concat[1024]` + `nbrcon`, `chaine[8]` + `longchaine[8]`, `fontcour`, `courgb[3]`, `iTe`, `iFa`, `ity`, `iep`, `iPo`, `ire`, `iRe`, `iel`, `iEl`, `iFP`, `menu`) become non-static members of `PsEmitter`. Single global instance owned by `XvueApp`: `XvueApp::psEmitter()` accessor. Mirrors Phase 5 D-02 (`XvueEventBridge`) and Phase 6.0 D-XX (`XvueMenuBridge`) bridge-as-XvueApp-singleton pattern.
- **D-08:** **PS font-name translation = hardcoded mapping table.** The X11 font-name string-bash in `xvchargefonte_` (lines ~1500–1568, parsing `-foundry-family-weight-slant-…` BDF font names) is dead code on the Qt path (no `XLoadQueryFont`). `PsEmitter::loadFont(QFont qf, …)` uses a 4-entry lookup table:
  ```
  QFont::family() → PS family
  ───────────────────────────
  "Courier"                → "/Courier"
  "Times"                  → "/Times-Roman"
  "Helvetica"              → "/Helvetica"
  "New Century Schoolbook" → "/NewCenturySchlbk"
  ```
  plus Italic/Oblique/Bold modifier suffix logic ported verbatim from xvuelc.c lines ~1530–1565. These four families are what MEFISTO actually loads (see `xvchargefonte` family table); other Qt fonts fall back to `/Helvetica` with a warn-once log entry.

### GIF capture & fallback

- **D-09:** **Standalone `QImageWriter` probe executable.** New `xvue/qt/probes/qimagewriter_probe.cpp` (10–20 lines: prints `QImageWriter::supportedImageFormats()` joined to stdout) built via new `bin/cb_probe_qt` script and CMake target. Run once at Phase 7 kickoff; stdout redirected to `.planning/phases/07-image-gif-and-postscript-export/PROBE.md` (EXPORT-01). PROBE.md commit is the kickoff gate. Defensive runtime re-check at first GIF write inside `PsEmitter::beginAnimation()` — if probe-confirmed support evaporates at runtime (plugin missing post-install), fallback path activates with a console-dock warning.
- **D-10:** **Native GIF path (probe-positive).** When `QImageWriter::supportedImageFormats()` lists `"gif"`, single `QImageWriter` opened with `setFormat("gif")`, each captured `QImage` written as a frame via the multi-frame extension (`setText("Delay","10")` for 10/100 sec = 100 ms ≈ 10 fps). No extra deps beyond `qt6-imageformats-plugins` (already in PROJECT.md install list).
- **D-11:** **ffmpeg fallback (probe-negative).** Write PNG sequence (`frame_0001.png`, `frame_0002.png`, …) into the temp dir from D-12, then `QProcess::execute("ffmpeg", {"-y", "-framerate","10", "-i","frame_%04d.png", "out.gif"})` blocking on the GUI thread (status-bar shows `"Encoding GIF…"` / FR `"Encodage GIF en cours…"`). testa/wave and testa/cavity2d animations are short (seconds-scale), so sync is acceptable for v1. ffmpeg is the **only** allowed external tool — ImageMagick `convert` is forbidden anywhere in `xvue/qt/`. Verified at phase end by `grep -rn 'convert' xvue/qt/` returning empty (EXPORT-06).
- **D-12:** **Frame buffer = mkstemp-pattern temp dir.** `QStandardPaths::writableLocation(QStandardPaths::TempLocation) + "/xvue-gif-XXXXXX/"` (Qt has no native mkdtemp — emulate with `QTemporaryDir`). Frames deleted on successful GIF write; kept on failure with the path logged to the console dock for diagnosis (Phase 6.0 console-dock surface). Project `cwd` stays clean except for the final `animation.gif`.

### PDF entry point

- **D-13:** **Raster PDF v1.** `QPdfWriter` + `QPainter::drawPixmap(0, 0, backing_)`. ~10 lines total. PDF page contains the canvas pixmap as embedded raster image. Fast, simple, fits EXPORT-05 "additive bonus" wording. True-vector PDF (which would require capturing every `QPainter` call into a display list and replaying onto `QPdfWriter`) is deferred to v2.
- **D-14:** **PDF triggered via Qt File menu only.** `File → Export → PDF…`. No Fortran ABI extension. Consistent with D-01.
- **D-15:** **Canvas-native page geometry at 72 dpi.** `QPdfWriter::setResolution(72)` + `setPageSize(QPageSize(QSizeF(xpixels, ypixels), QPageSize::Point))`. PDF page aspect matches canvas exactly — no fit-to-A4 distortion. xmin/ymin/xmax/ymax bbox written into PDF metadata mirrors the legacy PostScript header `%%BoundingBox` / DSC contract.

### ImageMagick scope (locked)

- **D-16:** **Phase 7 removes ImageMagick from the Qt code path only.** `bin/convertepsgif` and `convert` shell-outs from `td/`, `testa/`, `testf/` stay alive for the legacy X11 backend during the one-release-cycle A/B window. Full audit + removal is RETIRE-03 in Phase 9. Phase 7 verification: `grep -rn 'convert' xvue/qt/` returns empty (or only the PDF format keyword inside `QPageSize` enum text). Prevents downstream agents from over-deleting in Phase 7 and from blocking Phase 7 close on non-Qt convertepsgif callers that legitimately stay until Phase 9.

- **D-17:** **LVIDEO pipeline (`xvue/video1.f`, `xvue/videofin.f`, `xvue/videonm.f`) is OUT OF SCOPE for Phase 7.** The 07-RESEARCH.md inventory uncovered a second legacy animation pipeline independent of `xvpostscript_` and `bin/convertepsgif`: the LVIDEO `xwd` + `convert` shell-out called by 18+ Fortran tracer subroutines (`flui/trvi2d.f`, `ther/trplse.f`, `ther/trisot.f`, others). It lives in `xvue/*.f` (the legacy backend), NOT in `xvue/qt/`, so D-16's "ImageMagick from Qt code path only" wording covers it: Phase 7 leaves it alone. Full removal / Qt re-port is Phase 9 RETIRE-03 alongside `bin/convertepsgif` and `xvuelc.c`. EXPORT-06 grep gate for Phase 7 is `grep -rn 'convert' xvue/qt/` (empty), NOT `grep -rn 'convert' xvue/` (which would still legitimately match LVIDEO until Phase 9). This decision is recorded explicitly so future agents don't re-discover the LVIDEO pipeline and try to bundle it into Phase 7.

### Claude's Discretion

- **PsEmitter threading.** GUI thread only — `XVUE_QT_ASSERT_MAIN_THREAD()` at the top of every public method, matching the rest of `xvue/qt/src/xvue_qt_api.cpp`. `QImageWriter` and `QProcess::execute` called synchronously from main thread. Async / progress-dialog path deferred until evidence that testa/wave produces hundreds of frames (current observation: tens, not hundreds).
- **Default GIF frame delay.** 10/100 sec (≈ 10 fps). Matches typical scientific-viz cadence and the "convertepsgif → ImageMagick default 0.04 sec/frame" legacy timing roughly. Adjustable via `QSettings("xvue/export/gif_delay")` if testa/wave ends up looking too fast or slow.
- **JPEG quality default.** Qt default (90). Configurable via `QSettings("xvue/export/jpeg_quality")` if a user complains about file size or artifacts.
- **HiDPI export resolution.** `backing_` is already `devicePixelRatioF`-scaled. PNG/JPEG snapshot exports at backing pixel dims (i.e., `2*` logical for HiDPI displays). PDF stays at 72 dpi canvas-native logical pixels (vector-equivalent). Documented in `xvue/qt/README.md` Phase 7 section.
- **Disk-full / write-failure handling.** Status-bar message + console-dock log (Phase 6.0 surfaces). `QMessageBox::critical` only when menu-triggered (user expects a response); env-var-driven auto-snapshot path stays silent on failure but logs to console dock.
- **Filename-collision policy.** Menu-triggered: `QFileDialog` returns full path; if user picks an existing name, `QFileDialog` itself prompts for overwrite (Qt default). Auto-snapshot fixed names (`zfxy0NNN.png`) overwrite without warning — matches legacy convertepsgif behavior, and the only callers are testa/ scripts that expect deterministic filenames.
- **`xvpostscript_` body in `xvue_qt_api.cpp`.** Stays as the ABI-frozen entry; immediately delegates to `XvueApp::psEmitter().handleLasops(*lasops)`. Single dispatch line; PsEmitter owns all real logic.
- **Probe build wiring.** New `bin/cb_probe_qt` script (5–10 lines) + `xvue/qt/probes/CMakeLists.txt` target. PROBE.md committed at Phase 7 kickoff before the first plan executes (per EXPORT-01).

### Folded Todos

None — no pending todos matched Phase 7 scope.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase 7 scope

- `.planning/ROADMAP.md` §"Phase 7: Image, GIF, and PostScript export" (lines 230–240) — Goal, dependencies (Phase 6), 6 requirements (EXPORT-01..06), 5 success criteria, research-flag for GIF probe.
- `.planning/REQUIREMENTS.md` §"Export — image, GIF, and PostScript" (lines 88–95) — EXPORT-01..06 full text. Note EXPORT-04 wording: "moving the existing hand-rolled PostScript emitter (~120 lines of fprintf) **verbatim** from xvuelc.c; no switch to QPrinter PostScript output". EXPORT-05: PDF as bonus via QPrinter::PdfFormat NOT by modifying xvpostscript_.
- `.planning/PROJECT.md` §"Constraints" (lines 93–104) — Fortran ABI invariant, Linux-x86_64-only, no CI, manual A/B validation, single developer.
- `.planning/PROJECT.md` §"Key Decisions" (lines 106–121) — "Qt takes over image/GIF export", "Replace xvue/xvuelc.c entirely (not coexist forever)", "Big-bang Qt landing + one-release A/B window".

### PostScript invariant (Phase 0 audit, write-once until Phase 9)

- `xvue/README_COORDS.md` §"Where Y IS flipped (PostScript export only — Phase 7 concern)" (lines 38–52) — `ps_y = ypixels - screen_y` formula. Direct evidence at xvuelc.c:1895, 1932, 1953, 1966.
- `xvue/README_COORDS.md` §"Implications for the Qt migration" Phase 7 row (line 64) — **PRESERVE the `ypixels - y` flip verbatim** when porting; do NOT clean up by storing flipped coordinates.
- `xvue/README_COORDS.md` §"Anti-patterns" (lines 67–76) — three named anti-patterns to avoid in PsEmitter implementation.

### xvuelc.c source — exact line ranges to port

- `xvue/xvuelc.c:170-179` — file-static PS state declarations (`lasopsc`, `modepsc`, `chaine[]`, `longchaine[]`, `fpo`, `concat`, `nbrcon`, `fontcour`, `courgb[]`, `iTe`/`iFa`/`ity`/`iep`/`iPo`/`ire`/`iRe`/`iel`/`iEl`/`iFP`, `menu`). All become non-static `PsEmitter` members.
- `xvue/xvuelc.c:1187-1304` — `xvpostscript_(int *lasops)` function body (~120 lines, EXPORT-04 'verbatim' scope). Ports verbatim into `PsEmitter::handleLasops(int)`.
- `xvue/xvuelc.c:1500-1568` — `xvchargefonte_` PS-emit branch (X11 font-name string-bash). Replaced by `PsEmitter::loadFont(QFont, …)` with hardcoded mapping table per D-08.
- `xvue/xvuelc.c:1740-1820` — `xvtrait_` line-drawing PS-emit branch (Y-flip applied here at lines 1751, 1809). Ports as `PsEmitter::line(...)`.
- `xvue/xvuelc.c:1845-1975` — line-drawing block with explicit `ypixels - *y1` flips at 1895, 1932, 1953, 1966. Verbatim format strings → `PsEmitter::line()` body.
- `xvue/xvuelc.c:2050-2150` — `xvface_` / `xvfaceisocouleur_` PS branches. Port to `PsEmitter::face()`.
- `xvue/xvuelc.c:2580-2780` — `xvflpt_` / `xvellipse_` / additional primitive PS branches. Port to `PsEmitter::flpt()`, `PsEmitter::ellipse()`.
- `xvue/xvuelc.c:1414, 1435, 3030, 3033` — `proc(xvpostscript)(&lasopsc)` call sites within `effacer` and other primitives. These become `XvueApp::psEmitter().handleLasops(modeChange)` on the Qt side.
- `xvue/xvuelc.c:3140-3680` — PostScript header / dictionary / footer fprintf block (DSC comments, font ops, color palette). Verbatim format strings preserved.

### Existing Qt scaffolding (must integrate, not replace)

- `xvue/qt/src/xvue_qt_api.cpp:607-614` — current `xvpostscript_` warn-once stub. Replaced in Phase 7 with single-line dispatch to `XvueApp::psEmitter().handleLasops(*lasops)`.
- `xvue/qt/src/xvue_qt_app.{h,cpp}` — `XvueApp` singleton (Phase 1). New accessor: `XvueApp::psEmitter()` mirroring existing bridge accessors.
- `xvue/qt/src/xvue_qt_canvas.{h,cpp}` — `XvueCanvas::backing_` is the `QPixmap` snapshot source for image/GIF/PDF (Phase 2 D-04 single-backing).
- `xvue/qt/src/xvue_qt_window.{h,cpp}` — `XvueWindow` File menu construction site (Phase 6.0). New File → Export submenu added here.
- `xvue/qt/CMakeLists.txt` — adds `xvue_qt_postscript.cpp` to library sources; new `xvue/qt/probes/` subdirectory + CMake target for the QImageWriter probe.

### Prior CONTEXT decisions carried into Phase 7

- `.planning/phases/01-window-shell-xvueapp-xvuewindow-xvuecanvas/01-CONTEXT.md` — `XvueApp` singleton ownership pattern (singleton via `XvueApp::ensure()`).
- `.planning/phases/02-drawing-primitives-backing-pixmap/02-CONTEXT.md` D-04 — single-backing collapse: `backing_` is the visible surface, no separate fenetre/mempx distinction. Snapshot directly off `backing_`.
- `.planning/phases/04-pixmap-save-restore-double-buffering/04-CONTEXT.md` D-04 — same single-backing reaffirmed under solver tracer load.
- `.planning/phases/05-event-bridge-blocking-reads/05-CONTEXT.md` D-02 — bridge-as-QObject-singleton-on-XvueApp pattern. PsEmitter follows this shape.
- `.planning/phases/06.0-shared-shell-menu-bridge-dialogs/06.0-CONTEXT.md` D-04 — modifier-keyed Ctrl/Alt/Cmd routing via QShortcuts; D-07 console dock for diagnostics; D-09 bilingual FR/EN dialog text.

### Codebase maps

- `.planning/codebase/STRUCTURE.md` — `xvue/qt/` layout (src, tests, probes), pp/ executables, bin/ shell-build scripts.
- `.planning/codebase/CONVENTIONS.md` — `extern "C"` `proc(name)` macro pattern for Fortran-callable symbols, `XVUE_QT_ASSERT_MAIN_THREAD()` invariant on every public Qt API.
- `.planning/codebase/STACK.md` — Qt 6.10 from Debian apt, `qt6-imageformats-plugins` for GIF (probe-confirmed), CMake owns only `xvue/`, ffmpeg available from Debian apt.
- `.planning/codebase/CONCERNS.md` §"Graphics layer — X11/Motif" — `libX11-dev` and `convert` are required at runtime today; Phase 7 starts dropping the latter from the Qt path. §"Numerical / algorithmic fragility" — preserve numerical output (PS bytes are an output users depend on for paper figures).

### Legacy ImageMagick reference (informational only)

- `bin/convertepsgif` (single line: `convert -rotate 90 zfxy0*.eps -extent 980x550 cyl53zfxy.gif`) — the legacy post-processor Phase 7's GIF path replaces. The 90° rotation + 980×550 extent are case-specific to the cyl53 example and NOT universal; the new Qt path emits canvas-native dims and the testa/wave + testa/cavity2d cases must be re-checked to confirm whether per-case post-rotation is needed (researcher should grep `td/`, `testa/`, `testf/` for additional `convert` invocations beyond bin/convertepsgif).
- `bin.lnx64/convertepsgif` — duplicate of bin/convertepsgif; both stay alive until Phase 9 RETIRE-03.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets

- **`XvueApp` singleton** (`xvue/qt/src/xvue_qt_app.{h,cpp}`) — already hosts `XvueEventBridge` and `XvueMenuBridge` accessors; add `XvueApp::psEmitter()` following the same pattern.
- **`XvueCanvas::backing_`** (`xvue/qt/src/xvue_qt_canvas.{h,cpp}`) — `QPixmap` member already covers single-backing snapshot semantics. `QImageWriter` accepts `backing_.toImage()` directly.
- **`XvueWindow` File menu** (`xvue/qt/src/xvue_qt_window.{h,cpp}`) — Phase 6.0 shared shell already builds the File menu. Add Export submenu in the same builder.
- **Console dock surface** (Phase 6.0) — already wired for diagnostic logging; PsEmitter and ffmpeg fallback errors land here.
- **`XVUE_QT_ASSERT_MAIN_THREAD()` macro** — inherit from existing API contract; assert at every PsEmitter public method entry.
- **`bin/cb_*_qt` shell-build scripts** — pattern to follow for new `bin/cb_probe_qt`.

### Established Patterns

- **Bridge-as-singleton-on-XvueApp** (Phase 5 D-02, Phase 6 D-XX) — PsEmitter follows the same shape: QObject derivative owned by XvueApp, accessed via `XvueApp::psEmitter()`.
- **`proc(name)` `extern "C"` macro** (Phase 0 ABI freeze) — `xvpostscript_` stays as `proc(xvpostscript)`; only the body changes (now: dispatch one-liner).
- **Bilingual messages** (Phase 6.0 D-09) — every user-facing status-bar / console / message-box string has FR + EN variants selected via `langage_` flag. Apply to "Encoding GIF…" / "Encodage GIF…", "Animation capture started" / "Capture d'animation activée", etc.
- **Y-flip at emit time only** (Phase 0 audit, README_COORDS.md) — never store flipped Y; flip inside PsEmitter::line/face/text only.
- **CMake target = static library** (Phase 0 D-04) — PsEmitter sources land inside `libxvueqt.a`; no separate library.

### Integration Points

- **`xvue/qt/src/xvue_qt_api.cpp:607` `xvpostscript_` body** — replaced by single dispatch line.
- **Every Qt drawing primitive in `xvue/qt/src/xvue_qt_api.cpp`** that has a counterpart with `if(lasopsc>0) fprintf(fpo,...)` in `xvuelc.c` — adds an `XvueApp::psEmitter().<verb>(...)` call after its `QPainter::*` call. Touched primitives (researcher should produce the exact list from the xvuelc.c line ranges in canonical_refs): `xvtrait_`, `xvtraitcouleur_`, `xvface_`, `xvfaceisocouleur_`, `xvflpt_`, `xvellipse_`, `xvfond_`, `xvchargefonte_`, plus the `effacer` PS-erase trigger.
- **Phase 6.0 File menu builder in `XvueWindow::buildFileMenu`** — append Export submenu.
- **`xvue/qt/CMakeLists.txt`** — add `src/xvue_qt_postscript.cpp` to source list; add `add_subdirectory(probes)` target.
- **New `bin/cb_probe_qt` shell-build script** — wraps the CMake probes target build + run; redirects stdout into `.planning/phases/07-image-gif-and-postscript-export/PROBE.md`.

</code_context>

<specifics>
## Specific Ideas

- **`bin/convertepsgif` 90° rotation + 980×550 extent are case-specific, not universal.** The single line `convert -rotate 90 zfxy0*.eps -extent 980x550 cyl53zfxy.gif` is for the cyl53 example. testa/wave and testa/cavity2d may have their own per-case post-processing (or none). Researcher must `grep -rn 'convert\|convertepsgif' td/ testa/ testf/` to enumerate all callers and document per-case transform needs in RESEARCH.md before the GIF emit code is written. The Qt path emits canvas-native by default; per-case transforms (if any) become explicit options in the auto-snapshot configuration.
- **`xvpostscript_` is called recursively from `effacer` (xvuelc.c:1414, 1435) and possibly other primitives** with `lasopsc + 100`. The recursive call pattern triggers the "mode 100/101/102 erase" branch. PsEmitter::handleLasops must preserve this re-entrancy.
- **`chaine[8]` buffer is a per-menu accumulator** for "draw into this stored menu region" (lasopsc 3..10 routes into `chaine[lasopsc-4]`). Used so that menu redraws can replay PS output without re-emitting from scratch. Preserve buffer semantics; don't replace with `QString` unless byte-identical concatenation is verified.
- **`fpo` (FILE*) opens TEMPORAIRE.EPS unconditionally on `lasops>0 && lasopsc==0`** (the "start capturing" transition). Failure to open exits the process (`exit(1)`). On the Qt side this becomes a `QMessageBox::critical` + `XvueApp::quit()` — preserve the abort semantics, not the `exit(1)` strategy.

</specifics>

<deferred>
## Deferred Ideas

- **True-vector PDF.** Phase 7 ships raster-only PDF (D-13). True-vector PDF requires capturing every `QPainter` call into a display list and replaying onto `QPdfWriter`. Defer to v2; revisit if scientific users request vector PDFs for paper insertion.
- **Async GIF encoding for very large animations.** D-11 specifies sync `QProcess::execute`. Async + `QProgressDialog` becomes worth it only if testa/wave / testa/cavity2d are observed to produce hundreds of frames. Defer until evidence.
- **`QPageSetupDialog` for PDF page sizing.** D-15 picks canvas-native at 72 dpi as the only option. User-driven page-size dialog deferred until requested.
- **Multi-window `PsEmitter` instances.** D-07 puts state on a single XvueApp-owned instance. MEFISTO uses one window at a time; multi-window with independent PS captures is deferred (would move state to XvueState per-window).
- **`SVG` export.** Not in REQUIREMENTS.md; not requested. v2 candidate if vector PDF lands.
- **EPS variant separate from .ps.** Legacy emits `TEMPORAIRE.EPS`. Phase 7 emits the same byte stream under that name; users can rename to `.eps` or `.ps` freely. No separate emitter needed.
- **GIF compression / palette tuning.** Default Qt GIF plugin output is acceptable for testa/ visual A/B. Custom palette / dithering deferred.

### Reviewed Todos (not folded)

None — no pending todos surfaced for Phase 7 review.

</deferred>

---

*Phase: 07-image-gif-and-postscript-export*
*Context gathered: 2026-05-01*
