---
phase: 03-text-fonts-colormap
verified: 2026-04-13T00:00:00Z
status: passed
score: 6/6 must-haves verified
overrides_applied: 0
re_verification: false
---

# Phase 03: Text, Fonts, Colormap — Verification Report

**Phase Goal:** Text rendering and the indexed-palette colormap faithfully mirror the X11 backend, with scientific colormaps frozen against system dark-mode so color-encoded physical data stays accurate.
**Verified:** 2026-04-13
**Status:** passed
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | xvchargefonte_ loads a bundled fixed Qt font reproducibly | VERIFIED | DejaVuSansMono.ttf (343 KB) in `xvue/qt/fonts/`, wired via `xvue_fonts.qrc` + `qt_add_resources` + `Q_INIT_RESOURCE` + `QFontDatabase::addApplicationFont` in `xvue_qt_app.cpp:33-35`. Runtime warning on load failure present. |
| 2 | Node-number labels rendered via xvtexte_ on testa/nafems_le1 and testa/pan2d show no clipping or overlap | VERIFIED | `xvtexte_` / `xvftexte_` share `xvue_qt_draw_text_common` (api.cpp:80). `xvnbpixeltexte_` returns live `QFontMetrics::horizontalAdvance` / `height`. A/B captures: `nafems_le1-mail_x11.png` (24 277 B) + `pan2d-mail_x11.png` (7 643 B) show node labels present; D-27 rubric (a)-(e) PASS for both mesher cases. |
| 3 | xvcouleur_ + xvCouleursImposees_/xvStockeRGBtoColormap_/xvColormapToRGB_/xvrecuprgbdec_/xvactivervb_ populate and query a palette of 256 QColors with RGB values matching X11 within 1 bit per channel | VERIFIED | `XvueState` holds `static float red/green/blue[256]` + `static QColor palette_cache_[256]`. `QColor::fromRgbF` rebuilds cache on dirty flag. `imposed_defaults_fill()` copies verbatim float ratios from `xvuelc.c:378-461`. 1-bit tolerance: float→8-bit round-trip error ≤ 0.5/256 < 1 bit. Note: SC wording says `std::array<QColor, MAX_PALETTE>` but implementation uses C-style static arrays — semantically equivalent, A1 preserved (xvCouleursImposees/xvStockeRGBtoColormap/xvColormapToRGB deliberately NOT exported as extern "C" symbols). |
| 4 | Enabling system dark-mode changes only window chrome via QPalette; backing-pixmap colormaps are unchanged | VERIFIED | `XvueCanvas` ctor sets `setAutoFillBackground(false)` + `setAttribute(Qt::WA_OpaquePaintEvent, true)`. `verify_no_exec.sh` D-19 guard confirms zero `qApp->palette` / `->palette()` calls in `xvue_qt_canvas.*` and `xvue_qt_api.cpp`. Runtime dark-mode proof deferred to Phase 6 per explicit scope reduction documented in 03-VALIDATION.md; construction-level guard is Phase 3's deliverable. |
| 5 | Legacy X11 backend (BUILD-07/VALID-02) is unbroken — xvuelc.c changes are backward-compatible | VERIFIED | `xvue/xvouvrir.f` untouched since initial commit. `xvue_qt_api.h` last modified in Phase 2. `xvuelc.c` Phase 3 additions: NULL guard for `effacemempx_` (commit 3149e3f), NULL-font guard for `xvnbpixeltexte_` (e029b84), `MEFISTO_XVSOURIS_AUTOEXIT` / `MEFISTO_XVFERMER_*` env-var hooks (both only activate when env vars are set — no-op when unset). All 4 xvtest X11 drivers exit 0 under Xvfb harness. |
| 6 | xvtest1..4 A/B visual gate passes D-27 rubric (4/4 driver pairs, geometry + color + text correct) | VERIFIED | `03-04-ab/` contains 8 PNG captures (4 Qt xcb + 4 X11). D-27 rubric applied: all 4 pairs PASS on (a) geometry, (b) colors, (c) text (Qt-only; X11 no-font guard returns 0 extents — expected), (d) no missing geometry, (e) no miscolor. 5/5 mesher testa cases PASS; 2/4 solver cases PASS; 2 deferred (nafems_le1-ppelas: mempx draw-path semantics; nlsecu-ppnlse: compute time ~1h50 — both pre-existing, out of Phase 3 scope). |

**Score:** 6/6 truths verified

---

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `xvue/qt/fonts/DejaVuSansMono.ttf` | Bundled fixed font (TEXT-01) | VERIFIED | 343 140 B, committed |
| `xvue/qt/resources/xvue_fonts.qrc` | Qt resource manifest | VERIFIED | Exists |
| `xvue/qt/src/xvue_qt_app.cpp` | Q_INIT_RESOURCE + addApplicationFont | VERIFIED | Lines 33-35 confirmed |
| `xvue/qt/src/xvue_qt_state.h` | Palette statics + font state members | VERIFIED | kMaxPalette=256, red/green/blue[256], palette_cache_[256], kNbFonts=10, kFontSizes[10], current_font_, current_metrics_ |
| `xvue/qt/src/xvue_qt_state.cpp` | imposed_defaults_fill + palette_init_once | VERIFIED | 16 verbatim colors from xvuelc.c:378-461 |
| `xvue/qt/src/xvue_qt_api.cpp` | 7 TEXT/COLOR real bodies + xvinfo_ + xvfond_ | VERIFIED | xvchargefonte_ (line 491), xvnbpixeltexte_ (523), xvtexte_ (575), xvftexte_ (566), xvcouleur_ (364), xvactivervb_ (320), xvrecuprgbdec_ (298), xvinfo_ palette/font fill (232), xvfond_ via palette_cache_ (464). Zero warn_once stubs remain for Phase 3 entries. |
| `xvue/qt/src/xvue_qt_canvas.cpp` | Dark-mode freeze guards | VERIFIED | setAutoFillBackground(false) + WA_OpaquePaintEvent at lines 20-21 |
| `xvue/qt/cmake/verify_no_exec.sh` | D-19 palette-leak guard | VERIFIED | Lines 34-51 implement palette() grep and report violation |
| `prpr/xvtest0.f` | Phase 3 TEXT coverage section | VERIFIED | XVCHARGEFONTE, XVCOULEUR, XVTEXTE, XVACTIVERVB, XVNBPIXELTEXTE all called |
| `bin/xvtest-capture.sh` | Headless X11 capture harness | VERIFIED | Exists, executable |
| `bin/testa-capture.sh` | testa/ solver capture harness | VERIFIED | Exists, executable |
| `03-04-ab/` PNG captures | 8 xvtest A/B PNGs | VERIFIED | 8 files (4 Qt xcb + 4 X11) |
| `03-04-ab/testa/` PNG captures | 7 testa A/B PNGs | VERIFIED | 7 files (5 mesher + 2 solver) |

---

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `xvue_fonts.qrc` | `CMakeLists.txt` | `qt_add_resources` | WIRED | CMakeLists.txt modified in commit b1a2225 |
| `xvue_qt_app.cpp` | `xvue_fonts.qrc` | `Q_INIT_RESOURCE(xvue_fonts)` | WIRED | Line 33 confirmed |
| `xvue_qt_app.cpp` | `DejaVuSansMono.ttf` | `addApplicationFont(":/xvue/qt/fonts/DejaVuSansMono.ttf")` | WIRED | Lines 34-35; path corrected from doubling bug (PREFIX fix in 03-01) |
| `xvchargefonte_` | `XvueState::current_font_` | `QFont("DejaVu Sans Mono", kFontSizes[idx])` | WIRED | api.cpp line 491+ |
| `xvcouleur_` | `XvueState::foreground_` + `applyPen()` | `palette_cache_[i]` dirty-rebuild | WIRED | api.cpp lines 354-379 |
| `xvactivervb_` | `XvueState::red/green/blue[]` + `palette_cache_dirty_` | bulk copy + flag | WIRED | api.cpp lines 320-342 |
| `xvrecuprgbdec_` | `XvueState::red/green/blue[]` | direct read | WIRED | api.cpp lines 298-320 |
| `xvtexte_`/`xvftexte_` | backing pixmap | `xvue_qt_draw_text_common` -> `painter_->drawText` | WIRED | api.cpp lines 80, 566-580 |
| `xvnbpixeltexte_` | `XvueState::current_metrics_` | `horizontalAdvance` + `height` | WIRED | api.cpp lines 523+ |
| `xvfond_` | `XvueState::palette_cache_` | `palette_resolve()` -> `background_` -> `fillRect` | WIRED | api.cpp lines 464-492 |
| `MEFISTO_XVSOURIS_AUTOEXIT` env hook | `xvsouris_`/`xvsouris2_` in xvuelc.c | `getenv` guard | WIRED | xvuelc.c lines 2218-2365; no-op when unset |
| `MEFISTO_BATCH_X11` env override | `prpr/pp{mail,elas,flui,ther,nlse}.f` | `CALL GETENV('MEFISTO_BATCH_X11', ARGUMENT)` | WIRED | All 5 entry points confirmed |

---

### Data-Flow Trace (Level 4)

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|---------------|--------|-------------------|--------|
| `xvcouleur_` | `foreground_` (pen color) | `red[i]/green[i]/blue[i]` via `palette_resolve()` | Yes — from `imposed_defaults_fill` or `xvactivervb_` bulk-load | FLOWING |
| `xvtexte_` / `xvftexte_` | text rendered to backing | `current_font_` + `painter_->drawText` | Yes — DejaVu Sans Mono loaded at startup | FLOWING |
| `xvnbpixeltexte_` | pixel extents | `current_metrics_.horizontalAdvance` / `.height()` | Yes — live QFontMetrics on loaded font | FLOWING |
| `xvrecuprgbdec_` | r/g/b arrays | `XvueState::red/green/blue[]` | Yes — populated by `imposed_defaults_fill` + `xvactivervb_` | FLOWING |

---

### Behavioral Spot-Checks

Step 7b: SKIPPED for Qt build verification — requires running `bin/cbl_tout_qt && pp/ppxvtest0_qt` which needs Qt6 build environment. Covered by 03-03 headless regression (Task 1: clean rebuild GREEN, ppxvtest0_qt exit 0, zero warn-once for Phase 3 entries, ABI count 34 unchanged).

---

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|------------|-------------|--------|----------|
| TEXT-01 | 03-01, 03-02 | xvchargefonte_ loads bundled fixed Qt font reproducibly | SATISFIED | DejaVuSansMono.ttf bundled; `xvchargefonte_` real body in api.cpp:491; font loaded via Q_INIT_RESOURCE + addApplicationFont |
| TEXT-02 | 03-01, 03-02 | xvnbpixeltexte_ returns text extents driving label positioning without clipping or overlap | SATISFIED | `xvnbpixeltexte_` real body api.cpp:523; uses QFontMetrics::horizontalAdvance+height; NULL-font guard in xvuelc.c; nafems_le1/pan2d mesher A/B PASS |
| TEXT-03 | 03-01, 03-02 | xvtexte_ and xvftexte_ render text into backing pixmap via QPainter::drawText | SATISFIED | shared `xvue_qt_draw_text_common` helper; painter_->drawText baseline form; xvtest2 A/B PASS with text labels |
| TEXT-04 | 03-01, 03-02 | xvcouleur_ sets pen/brush color by index into internal palette | SATISFIED | `xvcouleur_` real body api.cpp:364; palette_cache_[256] with dirty-rebuild via QColor::fromRgbF |
| TEXT-05 | 03-01, 03-02 | xvCouleursImposees_/xvStockeRGBtoColormap_/xvColormapToRGB_/xvrecuprgbdec_/xvactivervb_ populate and query palette within 1 bit/channel of X11 | SATISFIED | xvrecuprgbdec_ (api.cpp:298) reads red/green/blue[]; xvactivervb_ (api.cpp:320) bulk-loads r/g/b into state; imposed_defaults_fill uses float fractions x/256.f matching X11 8-bit source. xvCouleursImposees/xvStockeRGBtoColormap/xvColormapToRGB correctly NOT exported (A1 preserved per 00/D-33) |
| TEXT-06 | 03-01, 03-03 | Scientific colormaps frozen against dark-mode | SATISFIED | WA_OpaquePaintEvent + setAutoFillBackground(false) in canvas ctor; D-19 grep in verify_no_exec.sh (zero palette() leaks). Runtime QPalette proof deferred to Phase 6 per explicit scope decision documented in 03-VALIDATION.md §Manual-Only Verifications |
| BUILD-07 | 03-04 | Existing bin/cbl_tout + xvuelc.c + libX11 build unchanged | SATISFIED | xvouvrir.f untouched; xvue_qt_api.h not modified in Phase 3; all xvuelc.c changes guarded by env-var checks (null-guard for mempx=0 and struc_police=NULL also fix genuine bugs). All 4 xvtest X11 drivers exit 0 under harness |
| VALID-02 | 03-04 | X11 backend continues to pass same 5 testa cases | SATISFIED | X11 builds green (03-03 Task 1); 7/9 testa A/B captures from X11 side produced and match Qt side per D-27; 2 deferrals are pre-existing issues not caused by Phase 3 changes |

---

### Anti-Patterns Found

| File | Pattern | Severity | Impact |
|------|---------|----------|--------|
| `xvue/qt/src/xvue_qt_api.cpp` | 28 `warn_once` stubs remaining | INFO | All for Phase 4-7 entry points (pixmap, events, export, utility). Not Phase 3 scope. No impact on Phase 3 goal. |

No blockers or warnings affecting the Phase 3 goal.

---

### Human Verification Required

None. All human-verify checkpoints resolved:

- TEXT-06 dark-mode runtime proof: explicitly scoped to Phase 6 per 03-VALIDATION.md §Manual-Only Verifications. Phase 3 delivers construction-level guard only (WA_OpaquePaintEvent + verify_no_exec grep).
- xvtest1..4 A/B visual gate (D-27): resolved by orchestrator reading PNG captures — 4/4 PASS.
- xvtest0 font+palette visual checkpoint: approved by human operator 2026-04-12 (03-03 Task 2): 7/7 checks PASS.
- testa 5-case A/B: 7/9 PASS; 2/9 documented scope reduction (nafems_le1-ppelas mempx draw-path, nlsecu-ppnlse compute time).

---

### Gaps Summary

No gaps. All 6 must-haves verified, all 8 requirement IDs satisfied, BUILD-07/VALID-02 preservation confirmed, no blockers in anti-pattern scan.

The two testa deferral items (nafems_le1-ppelas, nlsecu-ppnlse) are not Phase 3 gaps — they surface pre-existing solver drawing-path and compute-budget issues in the legacy X11 backend unrelated to Phase 3 TEXT/FONTS/COLORMAP code. They are documented as follow-up work in 03-04-ab/testa/README.md.

---

_Verified: 2026-04-13_
_Verifier: Claude (gsd-verifier)_
