# Phase 3: Text, fonts, colormap - Context

**Gathered:** 2026-04-11
**Status:** Ready for planning

<domain>
## Phase Boundary

Text rendering and the indexed-palette colormap faithfully mirror the X11 backend, with scientific colormaps frozen against system dark-mode so color-encoded physical data stays accurate. Phase 3 delivers TEXT-01 through TEXT-06: bundled reproducible font loading, `QFontMetrics`-driven text extents, `QPainter::drawText` rendering into the Phase 2 backing pixmap, the `std::array<QColor, 256>` indexed palette mirroring X11 colormap semantics, the full `xvcouleur_`/`xvCouleursImposees_`/`xvStockeRGBtoColormap_`/`xvColormapToRGB_`/`xvrecuprgbdec_`/`xvactivervb_` entry-point family, and a documented dark-mode-freeze invariant.

**In scope:**
- Ship one bundled TTF (`xvue/qt/fonts/DejaVuSansMono.ttf`) and load it via `QFontDatabase::addApplicationFont` at `XvueApp` construction (one-shot, idempotent). Legacy integer font numbers map to point sizes in a hardcoded table.
- Real implementations of: `xvchargefonte_` (TEXT-01), `xvnbpixeltexte_` (TEXT-02), `xvtexte_` + `xvftexte_` (TEXT-03), `xvcouleur_` (TEXT-04), `xvCouleursImposees_` + `xvStockeRGBtoColormap_` + `xvColormapToRGB_` + `xvrecuprgbdec_` + `xvactivervb_` (TEXT-05).
- Grow `XvueState` additively: add the static palette arrays (`float red[256]`, `float green[256]`, `float blue[256]`, and a lazily-computed `QColor palette_cache_[256]`) plus a `QFont current_font_` + `int current_font_size_idx_` pair and a `QFontMetrics current_metrics_` cache.
- Retire the Phase 2 hardcoded `foreground_ = Qt::white`: `xvcouleur_` writes into `XvueState::foreground_` via `palette_cache_[icolor]` and calls `applyPen()` (which Phase 2 already routes through `setBrush` per 02/D-20).
- Retire the Phase 2 minimal 2-entry `xvfond_` palette hack (currently at `xvue/qt/src/xvue_qt_api.cpp:337-358`) and route it through the real `palette_cache_`.
- Freeze the palette against `QEvent::PaletteChange`: the canvas's backing pixmap and painter never source colors from `qApp->palette()`. QPalette remains free to affect window chrome that later phases (6/UX-13) add.
- Populate the currently-zeroed `xvinfo_` palette/font outputs at `xvue_qt_api.cpp:189-223`: `maxfonts`, `nbfonts`, `namefonts[][256]`, `visuclass`, and the palette range ints that reference the real colormap.
- Extend `prpr/xvtest0.f` with a text + color coverage section exercising every Phase 3 entry point (pattern: same "add after existing holds" approach as the Phase 2 extension).
- Re-enable the deferred A/B gate against `prpr/xvtest1.f`–`xvtest4.f` and the 5-case `testa/` driver set (Phase 2 explicitly deferred these because they call `xvtexte_`/`xvchargefonte_`/`xvcouleurs_` which were warn-once stubs — Phase 3 is the catch-up gate).

**Explicitly out of scope:**
- Pixmap save/restore slots — Phase 4.
- Event bridge / mouse / keyboard — Phase 5.
- Menu bar, toolbar, dock widgets, `QMainWindow` chrome, dark-mode QPalette wiring — Phase 6. Phase 3's dark-mode-freeze work is a defensive posture only: it guarantees the backing pixmap *ignores* QPalette changes. The actual palette-change plumbing that Phase 6 needs is not wired up in Phase 3.
- Image / GIF / PostScript export — Phase 7. `lasopsc`/`fpo`/`concat[]` blocks from the legacy text/color bodies are dropped at translation time, same rule as Phase 2 (02/D-26).
- Any dependency on real input events — `xvsouris_`/`xvpause_`/`deplsouris_` stay warn-once stubs.
- `CR-01`/`WR-01`/`WR-02` latent color-divergence follow-ups from Phase 2 are **implicitly closed** by Phase 3: once Phase 3 introduces independent pen/brush colors, the existing rectangle/polygon stroke-vs-fill dispatcher (added in commit `6741c68`) starts producing visible correct output. No additional Phase 3 task needed for these; the Phase 3 A/B gate against `xvtest1..4` is the verification.

</domain>

<decisions>
## Implementation Decisions

### Font strategy

- **D-01:** A single monospace TTF is bundled in `xvue/qt/fonts/DejaVuSansMono.ttf`. License is Bitstream Vera / DejaVu (Bitstream Vera Fonts Copyright with free-software addendum — trivially redistributable, already in Debian `fonts-dejavu-core`). The file is committed to the repository so `bin/cbl_tout_qt` does not depend on any system font package.
- **D-02:** `XvueApp` constructor calls `QFontDatabase::addApplicationFont(":/xvue/qt/fonts/DejaVuSansMono.ttf")` (loaded via Qt resource system — CMake `qt_add_resources` in Phase 3 plan) exactly once. The returned font ID is stored in a `static int XvueApp::font_id_ = -1` guard; subsequent constructions skip the load. This is safe across `xvfermer_`/`xvinitgraphique_` cycles because `XvueApp` is deliberately leaked at atexit (01/D-08) — there is exactly one `XvueApp` per process.
- **D-03:** Legacy X11 integer font numbers map to a **hardcoded point-size table** rather than distinct families:
  | legacy `nofont` | Qt point size |
  |-----------------|---------------|
  | 0               | 8             |
  | 1               | 10            |
  | 2               | 12            |
  | 3               | 14            |
  | 4               | 16            |
  | 5               | 18            |
  | 6               | 20            |
  | 7               | 24            |
  | 8               | 28            |
  | ≥9              | 32            |
  The table lives as `constexpr int XvueState::kFontSizes[]` in `xvue/qt/src/xvue_qt_state.h`. `nbpolices_qt = sizeof(kFontSizes)/sizeof(int)`. `listfonts_qt[i]` is `"DejaVu Sans Mono 8pt"`-style display names used by `xvinfo_`.
- **D-04:** `XvueState::current_font_` is a `QFont("DejaVu Sans Mono", pt)` rebuilt every time `xvchargefonte_` runs. `XvueState::current_metrics_` is a `QFontMetrics` constructed from `current_font_` and cached. `xvchargefonte_`'s body: validate `nofont0` and `nofont` bounds, set `current_font_size_idx_ = *nofont`, rebuild `current_font_` with `kFontSizes[*nofont]`, call `painter_->setFont(current_font_)`, rebuild `current_metrics_`, write the metric outputs: `*largpx = current_metrics_.horizontalAdvance('0')` (match legacy "width of digit '0'" heuristic) and `*hautpx = current_metrics_.height()`.
- **D-05:** `xvnbpixeltexte_` body: convert the Fortran character buffer to `QString::fromLatin1(texte, *length)` (legacy uses `XDrawString` which is 8-bit Latin-1), then `*nbpxla = current_metrics_.horizontalAdvance(qstr)` and `*nbpxha = current_metrics_.height()`. No locale handling — matches legacy semantics.

### Text rendering

- **D-06:** `xvtexte_` and `xvftexte_` are **semantically identical** in Phase 3 under the single-backing model (same pattern as Phase 2's `xvtrait_`/`xvftrait_` collapse per 02/D-09). Both bodies: build a `QString::fromLatin1(string, *length)`, call `painter_->setFont(current_font_)` (cheap Qt no-op if already set), call `painter_->drawText(*x1, *y1, qstr)` — the **baseline form**, byte-for-byte matching legacy `XDrawString(display, ?, gc, *x1, *y1, string, *length)`. Followed by the standard D-01 epilogue from Phase 2 (`canvas_->update(); processEvents(ExcludeUserInputEvents)`). The two symbols stay exported separately for ABI preservation; a header comment on both references this D-06.
- **D-07:** Text rendering uses `QPainter::Antialiasing` (inherited from 02/D-22 — set once after every `painter_->begin()`) **plus** `QPainter::TextAntialiasing` which must be set alongside it. Phase 3 updates the `applyPen()` / painter-recreation path (currently at `xvue/qt/src/xvue_qt_canvas.cpp` resizeEvent) to add `painter_->setRenderHint(QPainter::TextAntialiasing, true);` right after the existing antialiasing line. This is a one-line additive change to the Phase 2 code.
- **D-08:** No attempt to preserve legacy baseline-offset quirks (X11's `XDrawString` sometimes mishandles descenders). `QPainter::drawText(x, y, str)` is the authoritative baseline; if `testa/nafems_le1` or `testa/pan2d` shows a 1–2 pixel Y drift vs X11, we document it as an acceptable drift in the Phase 3 A/B report and move on. TEXT-02 requires "no clipping or overlap", not "byte-identical Y".

### Palette storage and lifetime

- **D-09:** The palette lives as **`static` members of `XvueState`**, not instance members:
  ```cpp
  class XvueState {
    ...
    // Process-lifetime palette (matches legacy xvuelc.c static arrays at :100-104).
    static float red[256];
    static float green[256];
    static float blue[256];
    static unsigned long norgb[256];           // retained for xvrecuprgbdec_ compatibility
    static QColor palette_cache_[256];          // lazily rebuilt when r/g/b change
    static bool palette_cache_dirty_[256];      // per-index dirty flag
    ...
  };
  ```
  Class-scope `static` means the palette survives the `XvueWindow → XvueState` destruction that `xvfermer_` triggers, exactly like legacy `xvuelc.c` (whose palette arrays are file-scope static, also process-lifetime). On `xvinitgraphique_` (second call and later), the new `XvueState` instance sees the palette the previous session left behind. No explicit save/restore, no `XvueApp` ownership move — the static storage class handles it.
- **D-10:** `MAX_PALETTE = 256` — matches legacy `CMAPSIZE 256` (xvue/xvuelc.c:85). Not a gray area; pinned here so downstream grep finds it.
- **D-11:** Palette is initialized to an X11-default-like state exactly once, at the **first** `XvueState` construction in the process (guarded by a `static bool palette_initialized_ = false` inside the ctor). Default values: index 0 = black (0,0,0), index 1 = white (1,1,1), indices 2-15 = the standard X11 named colors from legacy `xvuelc.c:colorcell_defs` init path, indices 16-255 = black (caller must populate via `xvStockeRGBtoColormap_`). The initial default mirrors what `XAllocColor` produced on a freshly-opened X11 display.
- **D-12:** Storage format is **float `[0..1]`** in the `red/green/blue` arrays, matching the Fortran ABI at `xvrecuprgbdec_` / `xvStockeRGBtoColormap_` / `xvactivervb_` exactly. `QColor` values are computed lazily via `QColor::fromRgbF(r[i], g[i], b[i])` and cached in `palette_cache_[i]` on first access after a `palette_cache_dirty_[i] = true` flag flip. `xvrecuprgbdec_` literally returns `r[*i]`/`g[*i]`/`b[*i]` by value — no conversion, no allocation. `xvStockeRGBtoColormap_` writes into `r[i]`/`g[i]`/`b[i]` and sets `palette_cache_dirty_[i] = true`.
- **D-13 [AMENDED 2026-04-11 post-research]:** `norgb[256]` is **dropped entirely**. Research confirmed via exhaustive `.f/.F/.inc` grep that nothing reads or writes `norgb` — zero callers, no ABI exposure, no X11 Pixel indirection needed in the Qt model. Do NOT add the field to `XvueState`. Original D-13 "identity stub fallback" is superseded.

### Color entry-point bodies

- **D-14:** `xvcouleur_` body: bounds-check `*icolor` against `[0, 255]`, ensure `palette_cache_[*icolor]` is rebuilt if dirty (`QColor::fromRgbF(red[i], green[i], blue[i])`), set `state_->foreground_ = palette_cache_[*icolor]`, call `applyPen()` (which syncs both pen and brush per 02/D-20). No painter begin/end, no flush — this is a state-change operation, subsequent drawing primitives pick up the new color through their own epilogue.
- **D-15 [AMENDED 2026-04-11 post-research]:** `xvCouleursImposees` has **zero Fortran callers** (exhaustive grep) and was deliberately excluded from the Phase 0 public header (00-02-SUMMARY.md:43). Implement as a **file-static C++ helper** `static void xvCouleursImposees_helper(...)` inside `xvue_qt_api.cpp`, **not** as an `extern "C"` ABI symbol. Body still writes an RGB-triple array into `red[]`/`green[]`/`blue[]` starting at index 0 and marks indices dirty, but it is called only from other internal helpers (or not at all if no internal caller emerges). Keeps the 57-symbol ABI count intact per 00/D-33 literal parity.
- **D-16 [AMENDED 2026-04-11 post-research]:** `xvStockeRGBtoColormap` and `xvColormapToRGB` likewise have **zero Fortran callers** and are Phase-0 excluded. Implement both as **file-static C++ helpers** (`static void xvStockeRGBtoColormap_helper(int idx, float r, float g, float b)` and `static void xvColormapToRGB_helper(int idx, float* r, float* g, float* b)`) inside `xvue_qt_api.cpp`. Bodies are unchanged from the original D-16 (bounds-check + float store/load + dirty flag flip) but they are NOT exported. 57-symbol ABI preserved.
- **D-17 [AMENDED 2026-04-11 post-research]:** `xvactivervb_(nbcells, r, g, b)` is a **bulk palette-load**, NOT a transient foreground setter. The legacy body at `xvuelc.c:1072-1116` copies `r/g/b[0..nbcells-1]` into the static `red/green/blue[]` arrays and marks those indices dirty; the only live Fortran caller `xvue/palcde.f:619` passes `NDCOUL+1` cells — a full palette refresh. Phase 3 body: bounds-clamp `*nbcells` against `MAX_PALETTE`, copy `red[i] = r[i]; green[i] = g[i]; blue[i] = b[i]; palette_cache_dirty_[i] = true;` for `i ∈ [0, nbcells)`. Does NOT touch `state_->foreground_`, does NOT call `applyPen()`, does NOT begin/end the painter. This supersedes the original D-17 which described it as a transient one-shot write — that semantics would silently drop every color past index 0.
- **D-18:** `xvrecuprgbdec_(nbcolor, r, g, b)` body: copy `red[0..*nbcolor-1]`, `green[0..*nbcolor-1]`, `blue[0..*nbcolor-1]` into the three Fortran output arrays. Bounds-clamp `*nbcolor` against `MAX_PALETTE`. No palette cache touch — this is a read-only snapshot.

### Dark-mode freeze

- **D-19:** The dark-mode-freeze invariant is enforced by **construction**, not by event interception: `XvueCanvas` paints via `QPainter(this).drawPixmap(0, 0, *state_->backing_)` (01/D-05 + 02/D-05) — the backing pixmap's pixels are absolute RGB values constructed via `QColor::fromRgbF(r, g, b)` from the float arrays. Since no canvas-level `QColor` is ever sourced from `qApp->palette()` or `this->palette()`, a `QEvent::PaletteChange` has no path to modify backing-pixmap pixels. The invariant is: **`grep -rE 'qApp->palette|->palette\(\)' xvue/qt/src/xvue_qt_canvas.* xvue/qt/src/xvue_qt_api.*` returns zero matches.** Plan 03-04 (or wherever the verification task lands) adds this grep to `verify_no_exec` alongside the existing `lasopsc` guard.
- **D-20:** Defensive secondary guard: `XvueCanvas::setAutoFillBackground(false)` and `setAttribute(Qt::WA_OpaquePaintEvent, true)` are set in the `XvueCanvas` constructor. This removes the canvas from Qt's default style-system auto-fill path, so even if Phase 6 later adds a stylesheet that would recolor a widget's background, the canvas stays pixel-exact. Both attributes are additive to the Phase 1 canvas construction — touch the constructor once in a Phase 3 plan task.
- **D-21:** No `QEvent::PaletteChange` interception is added. The two construction-level guards (D-19 + D-20) are sufficient. Option B from discussion (override `changeEvent`) is explicitly rejected as belt-and-suspenders that adds code surface without closing any failure mode the construction-level guards miss. If Phase 6 discovers a leak, Phase 6 owns the fix.

### `xvinfo_` real implementation (TEXT-04/05 supporting)

- **D-22:** `xvinfo_` currently zeroes out `maxfonts`, `nbfonts`, `namefonts`, `visuclass`, and the palette range ints (`xvue_qt_api.cpp:189-223`). Phase 3 replaces the zero-fill with real values:
  - `*maxfonts = nbpolices_qt` (from D-03 table)
  - `*nbfonts  = nbpolices_qt`
  - `namefonts[k]` filled with `listfonts_qt[k]` strings, clamped to `*maxfonts`
  - `nbchar[k] = strlen(namefonts[k])`
  - `*visuclass = 4` (TrueColor — matches a 24-bit display, the only thing Qt 6 supports)
  - Palette range ints (legacy `ix`/`iy` params — check xvuelc.c:xvinfo body for their meaning at plan time)
- **D-23:** Font-name strings in `namefonts` are stable `constexpr const char*` — no dynamic allocation, no `malloc` / `XListFonts` equivalent. Matches legacy X11 semantic where `XListFonts` returns a server-allocated array the caller never frees (MEFISTO leaks it).

### Test driver extension (D-36 pattern)

- **D-24:** `prpr/xvtest0.f` gets a Phase 3 coverage section added **after** the Phase 2 two-hold structure (currently at :33-77 with two `SLEEP` holds). New section: (a) load font 2 (12pt) via `xvchargefonte_`, (b) populate palette indices 0..7 with the legacy-default 8-color block via `xvStockeRGBtoColormap_`, (c) set each of indices 0..7 as foreground via `xvcouleur_` and draw one short line + one text label per color (so the visual checkpoint can see "8 distinct colored lines with labels"), (d) exercise `xvactivervb_` by drawing one line in a custom RGB triple not in the palette, (e) exercise `xvnbpixeltexte_` by drawing a text label whose width was computed via `xvnbpixeltexte_` (visual check: no clipping beyond the computed box). Uses the existing third hold — no new `SLEEP` unless the existing holds are too short.
- **D-25:** The Phase 3 visual checkpoint (analog to 02-04 Task 2) must verify: (a) all 8 colored lines render in distinct colors matching the palette, (b) text labels are readable, (c) `xvactivervb_` line renders in the exact RGB triple passed, (d) `xvnbpixeltexte_`-computed bounding box matches the rendered text, (e) dragging the window resizes without color drift (DRAW-09 + TEXT-06 regression guard). Plan 03-04 owns this checkpoint.

### Legacy A/B catch-up gate

- **D-26:** Phase 3 re-enables the `prpr/xvtest1.f`–`xvtest4.f` and 5-case `testa/` A/B re-run that Phase 2 explicitly deferred (see 02/domain "Explicitly out of scope" para 7). This is the **hard gate** for Phase 3 completion: those drivers call every DRAW-01..09 + TEXT-01..06 entry point and producing byte-comparable output against the legacy X11 backend (`bin/cbl_tout` path) is the empirical proof that Phase 3 is done. One plan task in Phase 3 is dedicated to this — it's not a side-check of another plan.
- **D-27:** A/B output comparison is **visual by human eye**, not pixel-exact. The AA text renderer in Qt will produce subpixel differences from X11 core fonts even at identical point sizes. The success criterion is: "no clipping, no overlap, no missing glyphs, no miscolored regions, no missing geometry, no coordinate drift > 2 pixels." A pixel-diff tool is not Phase 3 scope; if the human reviewer wants one, it goes on the Phase 8 (A/B validation) checklist.
- **D-28:** A/B runs use `QT_QPA_PLATFORM=xcb` (native X11 path, not `offscreen`) on the Linux dev box. This keeps the Qt rendering path identical to what a real user would see and exercises the same compositor path as the legacy X11 backend. Pure headless rendering (`offscreen`) is used only for CI-style exit-code checks.

### Claude's Discretion

- Exact layout of the `xvue/qt/fonts/` directory and CMake integration — planner picks between `qt_add_resources` (binary-embedded TTF, smaller runtime dependency surface) and an explicit `xvue/qt/fonts` install path (simpler but adds a runtime file). `qt_add_resources` is recommended.
- Exact name of the new `03-VALIDATION.md` Per-Task Verification Map rows — follows the Phase 2 template.
- Whether the Phase 3 palette-init default (D-11, "X11 named-color mimicry") reads an embedded constant table in `xvue_qt_state.cpp` or computes from `QColor::fromString("black"/"white"/"red"/…)`. The former is faster and less magical; pick one at plan time.
- Whether `xvcouleurs_` (legacy plural, called during some `xvinitgraphique_` paths) gets a real body in Phase 3 or stays a warn-once shim — depends on planner's call-site audit. If any live caller exists, real body. Otherwise, minimal shim that routes to `xvCouleursImposees_` with a fixed index range.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Requirements and roadmap
- `.planning/ROADMAP.md` §"Phase 3: Text, fonts, colormap" — goal, dependencies, success criteria (goal statement + 4 numbered SCs)
- `.planning/REQUIREMENTS.md` §TEXT-01..TEXT-06 — acceptance criteria per entry-point family
- `.planning/REQUIREMENTS.md` §UX-13 — dark-mode chrome delegation (forward reference, Phase 6)

### Prior phase locked decisions
- `.planning/phases/00-build-skeleton-abi-stubs/00-CONTEXT.md` §decisions (D-33 literal ABI parity)
- `.planning/phases/01-window-shell-xvueapp-xvuewindow-xvuecanvas/01-CONTEXT.md` §decisions (D-01 exposure pump, D-05 paintEvent swap, D-08 QApplication lifetime leak)
- `.planning/phases/02-drawing-primitives-backing-pixmap/02-CONTEXT.md` §decisions (D-01 epilogue, D-05 long-lived painter, D-09 symbol-collapse pattern, D-14 xvfond hack to retire, D-20 brush sync, D-21 foreground_ placeholder, D-22 antialiasing reapply, D-26 PS-state drop rule)
- `.planning/phases/02-drawing-primitives-backing-pixmap/02-SUMMARY.md` files (01..04) — what actually shipped in Phase 2 and any deviations from 02-CONTEXT
- `.planning/phases/02-drawing-primitives-backing-pixmap/02-REVIEW-FIX.md` — CR-01/WR-01/WR-02 fixes that Phase 3 inherits (rectangle/polygon stroke+fill dispatch, null-guard symmetry); Phase 3 must not regress these

### Legacy reference (byte-for-byte parity target)
- `xvue/xvuelc.c:85` — `#define CMAPSIZE 256` (palette size source of truth)
- `xvue/xvuelc.c:100-104` — `static float red/green/blue[CMAPSIZE]` + `static unsigned long norgb[CMAPSIZE]` (storage layout mirror)
- `xvue/xvuelc.c:126` — `static XFontStruct *struc_police` (current-font state mirror)
- `xvue/xvuelc.c:131-132` — `int nbpolices` + `char **listfonts` (font metadata mirror)
- `xvue/xvuelc.c:481-500` — `XQueryColor` → float `r/g/b` conversion (ABI layout confirmation)
- `xvue/xvuelc.c:525-560` — `xvStockeRGBtoColormap` body (canonical float→storage path)
- `xvue/xvuelc.c:974-1025` — `xvinfo_` body (font-list enumeration reference)
- `xvue/xvuelc.c:1119-1175` — `xvcouleur_` body (indexed color activation reference)
- `xvue/xvuelc.c:1463-1580` — `xvchargefonte_` body (X11 font load path; Phase 3 replaces with bundled-TTF path but must match return-value semantics)
- `xvue/xvuelc.c:1650-1700` — `xvtexte_` / `xvftexte_` bodies (`XDrawString` baseline semantics target)
- `xvue/xvuelc.c:2050-2070` — `xvfacetraits_` ncf/nca reference (Phase 2 TODO being closed in Phase 3 once palette lands)

### Current Qt stubs to replace
- `xvue/qt/src/xvue_qt_api.cpp:189-223` — `xvinfo_` palette/font zero-out (Phase 3 populates per D-22/D-23)
- `xvue/qt/src/xvue_qt_api.cpp:228-245` — `xvrecuprgbdec_` warn-once stub (Phase 3 D-18)
- `xvue/qt/src/xvue_qt_api.cpp:246-280` — `xvcouleur_` warn-once stub (Phase 3 D-14)
- `xvue/qt/src/xvue_qt_api.cpp:337-358` — `xvfond_` Phase 2 hack (Phase 3 retires, routes through real `palette_cache_[]`)
- `xvue/qt/src/xvue_qt_api.cpp:380-420` — `xvchargefonte_` / `xvnbpixeltexte_` / `xvtexte_` / `xvftexte_` stubs (Phase 3 D-04/D-05/D-06)
- `xvue/qt/src/xvue_qt_state.h` + `.cpp` — grow additively per D-03 / D-09 / D-10 / D-12

### Project rules
- `CLAUDE.md` §"Compilation must never break" — Qt and legacy X11 builds both green after every Phase 3 change
- `CLAUDE.md` §"Programming norms" — `doc/normes.ps` (fixed-form Fortran 77 + comment style) for the `xvtest0.f` extension
- `doc/normes.ps` — authoritative norms (PostScript; view with `evince`)

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `XvueState` (xvue/qt/src/xvue_qt_state.{h,cpp}) — already owns `QPainter* painter_`, `QPixmap* backing_`, `QColor foreground_`, `QPen pen_`, `QBrush brush_`, `applyPen()` helper. Phase 3 grows it additively.
- `XvueCanvas::resizeEvent` — already reasserts `QPainter::Antialiasing` after every painter begin (02/D-22). Phase 3 adds `QPainter::TextAntialiasing` as a one-line sibling.
- Rectangle/polygon stroke+fill dispatcher (shipped in commit `6741c68`) — already handles outline vs fill vs both via enum dispatch. Phase 3 doesn't touch this code path; it benefits automatically once `foreground_` stops being hardcoded white.
- `warn_once` + `register_stub_name` macro (xvue/qt/src/xvue_qt_api.cpp Phase 0 infra) — Phase 3 *removes* warn-once lines for every TEXT-01..06 entry point the way Phase 2 removed them for DRAW-01..09. Do NOT leave any stubs in place once an entry point is implemented.
- `verify_no_exec` CMake target — already greps for PS-state tokens (`lasopsc`/`courgb`/`ypixels`/`iFa`/`iRe`/`iEl`/`iel`). Phase 3 adds the D-19 grep alongside: `qApp->palette|->palette\(\)` must return zero matches in `xvue_qt_canvas.*` and `xvue_qt_api.*`.
- `bin/cbxvtest0_qt` — the Phase 2-extended build script rebuilds `prpr/xvtest0.f` into `pp/ppxvtest0_qt` with the Qt backend. Phase 3 driver extension reuses this script unchanged.

### Established Patterns
- **Literal ABI parity** — every public `void proc(xvfoo)_(…)` keeps the byte-identical signature from `xvue/xvuelc.c` (00/D-33). Phase 3's five new color entry points and three new text entry points must be diffed byte-for-byte against the legacy header before the body is written.
- **Symbol-collapse** — `xvtrait_` / `xvftrait_` share a body (02/D-09); `xvrectangle_` / `xvbordrectangle_` / `xvfrectangle_` / `xvfbordrectangle_` share an enum-dispatched body (02/D-13 + CR-01 fix). Phase 3 applies the same rule to `xvtexte_` / `xvftexte_` (D-06) and to any color entry points that turn out to be synonyms under the Qt model.
- **Null-guard symmetry** — WR-03 from Phase 2 review established that every pointer parameter gets a null check at entry, no exceptions. Phase 3 entry points must do the same or the verifier will flag them.
- **No PS-recording side effects** — `lasopsc`, `fpo`, `concat[]`, `ire`/`iel`/`iFa`/`iRe`/`iEl` blocks from legacy bodies are dropped at translation time (02/D-26). `xvchargefonte_` and `xvtexte_` both have these blocks in the legacy — strip them.
- **State-change entry points don't flush** — `xvcouleur_`, `xvtypetrait_`, `xvepaisseur_` set state and return; they do not call `update()` or `processEvents()` (02/D-01). Only drawing primitives flush. Applies to all Phase 3 color/font setters.
- **Drawing entry points flush via D-01 epilogue** — `xvtexte_` / `xvftexte_` are drawing primitives and must end with `canvas_->update(); QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);`.

### Integration Points
- `XvueApp` constructor — one-shot `addApplicationFont` call (D-02)
- `XvueState` header — add palette statics + font state (D-09 / D-04)
- `XvueState` constructor — palette default-init guard (D-11)
- `XvueState::applyPen()` — already reads `foreground_`, no Phase 3 change needed; `xvcouleur_` writes `foreground_` and calls `applyPen()` (D-14)
- `XvueCanvas::resizeEvent` — add `TextAntialiasing` render hint (D-07)
- `XvueCanvas` constructor — add `setAutoFillBackground(false)` + `Qt::WA_OpaquePaintEvent` (D-20)
- `xvue_qt_api.cpp` — eight real bodies replace eight warn-once stubs; `xvfond_` Phase 2 hack retired; `xvinfo_` zero-fill replaced
- CMake `xvue/qt/CMakeLists.txt` — `qt_add_resources` declaration for the bundled TTF (D-02)
- `bin/verify_no_exec` (or equivalent CMake target) — add the `qApp->palette` grep (D-19)
- `prpr/xvtest0.f` — extended draw-coverage section with colors + text (D-24)

</code_context>

<specifics>
## Specific Ideas

- "DejaVu Sans Mono is the safest" — bundled mono TTF per D-01/D-02. No Liberation Mono, no Noto, no system fallback.
- "palette survive open/close" — satisfied by class-scope `static` storage (D-09), matching legacy `xvuelc.c` file-scope static arrays byte-for-byte.
- "Raw RGB, never QPalette" — construction-level guard only (D-19 + D-20). No `changeEvent` interception.
- "Use the baseline form `drawText(x, y, str)` for literal parity" — D-06, exactly matches `XDrawString` semantics.
- "Preserve Fortran float arrays" — storage format is float `[0..1]` matching the Fortran ABI (D-12), with a lazy `QColor` cache rebuilt only on write.
- Phase 3 is where the deferred `xvtest1..4` + `testa/` A/B runs finally close — this is the hard gate for phase completion (D-26/D-27), not a side-verification.
- Dark-mode-freeze is a construction-level invariant, not a runtime event handler. One grep in `verify_no_exec` is sufficient proof (D-19).

</specifics>

<deferred>
## Deferred Ideas

- **Pixel-diff A/B tooling** — a scripted pixel-by-pixel comparison between X11 and Qt backing pixmaps is tempting but belongs in Phase 8 (A/B validation on testa subset), not Phase 3. Phase 3 uses human visual review with a ≤2-pixel coordinate tolerance and a "no clipping/overlap/missing geometry" rubric (D-27).
- **Multiple font families** — rejected in D-03. Legacy X11 font numbers become point sizes on one mono family. If the mesher turns out to depend on a proportional font for any label, it's a Phase 6 chrome concern, not Phase 3.
- **Font metric Y-drift precision** — the ≤2-pixel tolerance in D-08 / D-27 defers pixel-exact baseline reproduction to Phase 8 if anyone cares.
- **Stylesheet / theme support** — Phase 6 concern. Phase 3 only guarantees the canvas ignores `qApp->palette()`.
- **Palette serialization** — if a user wants to save a custom colormap to disk and reload it next session, that's a new capability for a future phase (probably Phase 6 alongside user preferences). Phase 3 keeps the palette in process-lifetime static storage only.
- **Internationalized text (UTF-8 / complex scripts)** — legacy uses `XDrawString` which is 8-bit Latin-1. Phase 3 matches this via `QString::fromLatin1`. Full UTF-8 support is a rewrite that belongs in a future chrome phase, not Phase 3.
- **`norgb` usage decision** — D-13 defers the final call (drop vs identity-stub) to planner RESEARCH when the call-site audit lands. Not a gray area for the user.

</deferred>

---

*Phase: 03-text-fonts-colormap*
*Context gathered: 2026-04-11*
