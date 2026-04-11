# Phase 3: Text, fonts, colormap — Research

**Researched:** 2026-04-11
**Domain:** Qt 6 text rendering, font loading via Qt resources, indexed palette mirroring legacy X11 colormap semantics, A/B catch-up gate against legacy drivers.
**Confidence:** HIGH for ABI shape / call-site audit / legacy bodies / Qt API; MEDIUM for font-metric Y-drift envelope; LOW for dark-mode palette-leak absence (cannot be proven without Phase 6 chrome in place).

## Summary

Phase 3 is a catch-up phase. The Qt backend must implement **8 extern "C" Fortran entry points** (`xvchargefonte_`, `xvnbpixeltexte_`, `xvtexte_`, `xvftexte_`, `xvcouleur_`, `xvrecuprgbdec_`, `xvactivervb_`, plus completing the palette/font outputs of the already-partial `xvinfo_`). The real-body count is **8, not the "eight to thirteen" CONTEXT.md implies** — see Assumptions Log A1 below. The three "colormap helpers" (`xvCouleursImposees`, `xvStockeRGBtoColormap`, `xvColormapToRGB`) are deliberately NOT in the Fortran ABI per Phase 0 decision; they existed only as internal C helpers in `xvuelc.c` and have zero Fortran call sites. Phase 3 should implement them as **static C++ helpers inside `xvue_qt_api.cpp`** (or a new `xvue_qt_palette.cpp`) used only from other `xvue_qt_*` translation units — same architectural role as their legacy C counterparts.

The dominant technical risks are (a) ensuring the palette default-init runs **before** the first `xvrecuprgbdec_` call that `xvue/xvinit.f:143` makes (the legacy chain was xvinfo_→XAllocColor→XQueryColor→xvrecuprgbdec, so parity requires the default 16 colors to be loaded into `red[]/green[]/blue[]` at the same point xvinfo returns — not at first xvcouleur_ call); (b) CMake wiring of `qt_add_resources` so the bundled TTF is embedded in `libxvueqt.a` (a static library — some Qt resource patterns assume a shared-lib target); (c) the catch-up A/B run against 5 `testa/` cases needs a repeatable visual-diff procedure.

**Primary recommendation:** Implement the 8 extern "C" bodies plus a process-lifetime palette in class-static `XvueState` members. Load the bundled `DejaVuSansMono.ttf` via `qt_add_resources` + `Q_INIT_RESOURCE` (required because `libxvueqt.a` is STATIC). Populate the palette defaults from a hand-written table lifted from `xvuelc.c:378-460` (`xvCouleursImposees` body) triggered inside the new `xvinfo_` body — same call site the legacy chain uses. Treat `xvCouleursImposees` / `xvStockeRGBtoColormap` / `xvColormapToRGB` as **C++ file-static helpers**, never ABI symbols. Re-enable `xvtest1..4` drivers and the 5-case testa A/B as the hard gate, with human visual review using the ≤2-pixel coordinate tolerance rubric from D-27.

---

## User Constraints (from CONTEXT.md)

### Locked Decisions

**Font strategy:**
- **D-01** — One TTF bundled in `xvue/qt/fonts/DejaVuSansMono.ttf`, Bitstream Vera/DejaVu license, committed to repo. No system-font dependency.
- **D-02** — `XvueApp` ctor calls `QFontDatabase::addApplicationFont(":/xvue/qt/fonts/DejaVuSansMono.ttf")` exactly once, ID stored in `static int XvueApp::font_id_ = -1` guard. Safe across fermer/init cycles because XvueApp is leaked at atexit (01/D-08).
- **D-03** — Legacy `nofont` integers map to a hardcoded point-size table (0→8, 1→10, … 7→24, 8→28, ≥9→32). Table lives as `constexpr int XvueState::kFontSizes[]`. `nbpolices_qt = sizeof(kFontSizes)/sizeof(int)`. `listfonts_qt[i]` are `"DejaVu Sans Mono Npt"` display strings for `xvinfo_`.
- **D-04** — `xvchargefonte_` body: validate bounds on `*nofont`, rebuild `current_font_ = QFont("DejaVu Sans Mono", kFontSizes[*nofont])`, `painter_->setFont(current_font_)`, rebuild `current_metrics_`, return `*largpx = current_metrics_.horizontalAdvance('0')`, `*hautpx = current_metrics_.height()`.
- **D-05** — `xvnbpixeltexte_`: `QString::fromLatin1(texte, *length)`, then `horizontalAdvance(qstr)` / `height()`. No locale handling.

**Text rendering:**
- **D-06** — `xvtexte_` / `xvftexte_` collapsed to one body (pattern mirrors 02/D-09). Both call `painter_->drawText(*x1, *y1, QString::fromLatin1(string,*length))` + Phase 2 D-01 epilogue. ⚠️ **See Assumptions Log A2** — legacy bodies are NOT identical in xvuelc.c; this is a deliberate Qt simplification under the single-backing model, not literal parity.
- **D-07** — `QPainter::TextAntialiasing` is set alongside `QPainter::Antialiasing` in the existing `applyPen()` / painter-rebuild path at `xvue/qt/src/xvue_qt_canvas.cpp::resizeEvent`. One-line additive change.
- **D-08** — No attempt to preserve X11 baseline quirks. Acceptable drift is ≤2 pixel Y. TEXT-02 requires "no clipping/overlap", not byte-identical Y.

**Palette storage:**
- **D-09** — Palette is `static` members of `XvueState` (class-scope, process-lifetime): `static float red[256]`, `green[256]`, `blue[256]`, `unsigned long norgb[256]`, `QColor palette_cache_[256]`, `bool palette_cache_dirty_[256]`. Survives `XvueWindow→XvueState` destruction on `xvfermer_`.
- **D-10** — `MAX_PALETTE = 256` (matches legacy `CMAPSIZE 256`, `xvue/xvuelc.c:85`).
- **D-11** — Defaults loaded once, first `XvueState` construction, guarded by `static bool palette_initialized_`. Index 0=black, 1=white, 2-15 = the 16-color block from `xvuelc.c:378-460` (`xvCouleursImposees` body — verbatim lift), 16-255 = black until caller populates.
- **D-12** — Storage format is `float[0..1]` to match the Fortran ABI at `xvrecuprgbdec_`. `QColor` computed lazily via `QColor::fromRgbF(r[i],g[i],b[i])` and cached in `palette_cache_[i]`.
- **D-13** — `norgb[256]` kept as storage-only for ABI compatibility. Planner must audit whether any caller reads/writes it. **Research finding: ZERO Fortran callers touch norgb. ZERO Qt entry-point signatures take it. Phase 0 did not declare it in the public header.** Recommend: drop `norgb[]` entirely from Phase 3, or keep it with identity `norgb[i] = i` inside `xvStockeRGBtoColormap` internal helper for its one legacy use as `gcvalues.foreground = norgb[foreground]` (which has no Qt analog anyway — Qt uses `QColor` directly).

**Color bodies:**
- **D-14** — `xvcouleur_`: bounds-check `*icolor`, rebuild `palette_cache_[*icolor]` if dirty, `state_->foreground_ = palette_cache_[*icolor]`, call `applyPen()`. No flush (02/D-01).
- **D-15** — `xvCouleursImposees_` body — writes 16-color default block. ⚠️ **See A1** — this is NOT an extern "C" entry point; it's an internal C++ helper in Phase 3.
- **D-16** — `xvStockeRGBtoColormap_(idx,r,g,b)` write; `xvColormapToRGB_(idx,r,g,b)` read. ⚠️ **See A1** — internal C++ helpers, not ABI.
- **D-17** — `xvactivervb_(palcour, nbcells, r[], g[], b[])` body: write `r/g/b[]` into the static palette, mark indices dirty. Signature has four args (not "construct transient QColor" as the D-17 body text implies). See A3.
- **D-18** — `xvrecuprgbdec_(nbcolor, r, g, b)`: copy `red/green/blue[0..*nbcolor-1]` into Fortran output arrays, clamp to `MAX_PALETTE`.

**Dark-mode freeze:**
- **D-19** — Invariant enforced by construction: backing pixmap pixels are absolute RGB from `QColor::fromRgbF`, never sourced from `qApp->palette()` / `this->palette()`. Verified by a `verify_no_exec` grep: `qApp->palette|->palette\(\)` returns zero in `xvue_qt_canvas.*` and `xvue_qt_api.*`.
- **D-20** — Defensive guards in `XvueCanvas` ctor: `setAutoFillBackground(false)` + `setAttribute(Qt::WA_OpaquePaintEvent, true)`.
- **D-21** — No `QEvent::PaletteChange` interception. Construction-level guards are sufficient.

**xvinfo:**
- **D-22** — `xvinfo_` replaces zero-fill at `xvue_qt_api.cpp:189-223` with: `*maxfonts = *nbfonts = nbpolices_qt`, `namefonts[k]` ← `listfonts_qt[k]`, `nbchar[k] = strlen(...)`, `*visuclass = 4` (TrueColor), plus palette range ints.
- **D-23** — Font-name strings are `constexpr const char*`. No dynamic allocation.

**Test driver:**
- **D-24** — `prpr/xvtest0.f` extended after the existing two-hold Phase 2 structure (currently at :33-97): load font 2 (12pt), populate palette 0..7 via the new palette writer, draw a colored line + label per index, exercise `xvactivervb_`, exercise `xvnbpixeltexte_`.
- **D-25** — Visual checkpoint verifies 8 distinct colored lines+labels, `xvactivervb_` custom triple renders, `xvnbpixeltexte_` bounding box matches, resize-during-hold shows no color drift.

**A/B catch-up:**
- **D-26** — Re-enable `prpr/xvtest1.f`–`xvtest4.f` + 5-case `testa/` A/B run (deferred in Phase 2 per 02/D-35). **Hard gate** for Phase 3 completion.
- **D-27** — Visual comparison by human eye, ≤2-pixel coordinate tolerance, rubric = "no clipping/overlap/missing glyphs/miscolored regions/missing geometry". No pixel-diff tooling.
- **D-28** — A/B runs use `QT_QPA_PLATFORM=xcb` (not `offscreen`) on Linux dev box. `offscreen` only for exit-code smoke.

### Claude's Discretion

- CMake integration for the bundled TTF — `qt_add_resources` (embedded) vs explicit install path. **Recommended: `qt_add_resources`**. See Architecture Patterns §1.
- Name/rows of the `03-VALIDATION.md` Per-Task Verification Map — follow Phase 2 template.
- Whether palette defaults come from an embedded constant table vs `QColor::fromString`. **Recommended: verbatim table lifted from `xvuelc.c:378-460`** — faster, less magical, and preserves the legacy RGB values exactly (the X11 code uses specific float constants like `50./256.`, `178./256.`, `220./256.` that `fromString("gray1_bluish")` cannot match).
- `xvcouleurs_` (plural) disposition. **Research finding: `xvcouleurs_` DOES NOT EXIST.** Zero occurrences in `xvuelc.c`, zero Fortran `CALL XVCOULEURS`, zero stubs in `xvue_qt_api.cpp`. It is a ghost name repeated across several `.planning/phases/02-*` documents. Nothing to implement. Nothing to shim.

### Deferred Ideas (OUT OF SCOPE)

- Pixel-diff A/B tooling (→ Phase 8)
- Multiple font families (rejected in D-03; proportional fonts = Phase 6 chrome if needed)
- Font-metric Y-drift precision beyond ≤2 px tolerance (→ Phase 8 if anyone cares)
- Stylesheet / theme support (→ Phase 6)
- Palette serialization (save/load custom colormap to disk) (→ future UX phase)
- UTF-8 / complex scripts (Phase 3 uses `QString::fromLatin1` to match `XDrawString` 8-bit semantics)
- `norgb` final disposition was deferred to this research. Finding above: drop or identity-stub inside the internal `xvStockeRGBtoColormap` helper. Planner recommendation: drop.

---

## Phase Requirements

| ID | Description (from REQUIREMENTS.md) | Research Support |
|----|------------------------------------|------------------|
| TEXT-01 | `xvchargefonte_` loads bundled fixed Qt font reproducibly, returns pixel width/height | Architecture Pattern 1 (qt_add_resources), Code Example 1 (addApplicationFont), Pitfall 1 (static-lib Q_INIT_RESOURCE), D-02/D-04 |
| TEXT-02 | `xvnbpixeltexte_` returns extents via `QFontMetrics::horizontalAdvance`/`height`; testa/nafems_le1 + testa/pan2d labels no clipping/overlap | Code Example 2, D-05; State of the Art: horizontalAdvance is the non-deprecated API in Qt 6 (replaces Qt 5 `width()`) |
| TEXT-03 | `xvtexte_`/`xvftexte_` render text via `QPainter::drawText` into backing pixmap | Code Example 3, D-06, A2 deviation note |
| TEXT-04 | `xvcouleur_` sets pen+brush color by index into `std::array<QColor, MAX_PALETTE>` | D-14, already routed through `applyPen()` which 02/D-20 made brush-aware |
| TEXT-05 | `xvCouleursImposees_` / `xvStockeRGBtoColormap_` / `xvColormapToRGB_` / `xvrecuprgbdec_` / `xvactivervb_` populate/query/activate palette within 1 bit/channel of X11 | Call-site audit (below) — only `xvrecuprgbdec_` and `xvactivervb_` are real Fortran ABI entry points; the three others are **internal C++ helpers only**. REQUIREMENTS.md TEXT-05 names them by their legacy internal-function names; the test (palette values within 1 bit/channel) is still valid and satisfied by the float[0..1] storage format (D-12) with no 8-bit-quantization trip. |
| TEXT-06 | `QPalette` changes affect chrome only; backing-pixmap colormaps unchanged | D-19/D-20/D-21, Pitfall 5 |

---

## Call-Site Audit (CRITICAL FINDING)

**This section answers CONTEXT.md's two deferred research questions (D-13 and "Claude's Discretion" re `xvcouleurs_`).**

### Fortran `CALL` sites for Phase 3 entry points (verified via grep across all `.f/.F/.inc` files)

| Legacy name | Fortran call count | Files | Status |
|-------------|--------------------|-------|--------|
| `XVCOULEUR` | many (30+) | `xvue/*.f` (`vise2d`, `vise3d`, `texte2d`, `traitcoul3d`, `loaret`, `fap32d`, `face2d`, `face3d`, `trait3d`, etc.), `mail/*.f` | Real ABI entry point. Must be implemented. |
| `XVTEXTE` | many | `xvue/vise2d.f`, `vise3d.f`, `texte2d.f`, `entier3d.f`, `prpr/xvtest*.f` | Real ABI. Must be implemented. |
| `XVFTEXTE` | 2 | `xvue/recttxsu.f:123,174` (menu rectangle text — unused in Phase 3 until Phase 6, but still must link) | Real ABI. Implement (D-06 collapse). |
| `XVCHARGEFONTE` | 1 | `xvue/chargefonte.f:36` (`CALL XVCHARGEFONTE(NOFONT0, NOPOCA, NPLACA, NPHACA)` — 4 int args, signature confirmed) | Real ABI. Must be implemented. |
| `XVNBPIXELTEXTE` | 8 | `xvue/lamxpxtxt.f`, `xvue/logo.f`, `mail/legsef.f`, `prpr/xvtest1.f:53,140`, `prpr/xvtest3.f:129,132,135` | Real ABI. Must be implemented. |
| `XVRECUPRGBDEC` | **1** | `xvue/xvinit.f:143` — called immediately after `XVINFO` inside project init, **before any drawing**. Reads back the 16 standard colors that the C backend populated into `red[]/green[]/blue[]` so Fortran-side `PROUGE/PVERT/PBLEU` arrays match. | Real ABI. Must be implemented. **Requires palette defaults to be loaded by the time `xvinfo_` returns** — not lazily on first `xvcouleur_` call. |
| `XVACTIVERVB` | **1** | `xvue/palcde.f:619` (`CALL XVACTIVERVB(NOPALC, NDCOUL+1, PROUGE, PVERT, PBLEU)` — "palette change for user-defined colormap") | Real ABI. Must be implemented. |
| `XVCOULEURS` (plural) | **0** | — | **Does not exist.** Ghost name in CONTEXT.md and 02-*.md docs. No implementation work. Do NOT add a stub. |
| `XVCOULEURSIMPOSEES` | **0** | — | Not an ABI symbol. Internal C helper only in legacy xvuelc.c:358-461. Keep as C++ file-static helper. |
| `XVSTOCKERGBTOCOLORMAP` | **0** | — | Not an ABI symbol. Internal C helper in legacy xvuelc.c:503-559. Keep as C++ file-static helper. |
| `XVCOLORMAPTORGB` | **0** | — | Not an ABI symbol. Internal C helper in legacy xvuelc.c:463-500. Keep as C++ file-static helper. |

**Implication:** Phase 3 lands **8 real bodies** (counting `xvinfo_` palette/font fill as a body rewrite, not a new body): `xvchargefonte_`, `xvnbpixeltexte_`, `xvtexte_`, `xvftexte_`, `xvcouleur_`, `xvrecuprgbdec_`, `xvactivervb_`, `xvinfo_`. **Plus three C++ file-static helpers** (no Fortran exposure, no extern "C", no ABI symbol): `imposed_defaults_fill()`, `store_rgb_to_palette()`, `palette_to_rgb()`. The C++ helpers live in `xvue_qt_api.cpp` (or a new `xvue_qt_palette.cpp` if a separate TU is preferred).

**The "five color entry points" phrasing in REQUIREMENTS.md TEXT-05 is historically accurate** — in the legacy C code, they were all function names under the `proc()` macro. But Phase 0 deliberately removed three from the Fortran-facing header (`xvue/qt/include/xvue_qt_api.h`) because they had no Fortran callers and were internal helpers. This was a correct call. Phase 3 should implement them as internal helpers, preserving the TEXT-05 behavioral test (palette values within 1 bit/channel of X11) without resurrecting them as ABI symbols.

### `norgb[]` audit

- **Uses of `norgb` in xvuelc.c:** 7 reads/writes total (file scope), all internal. Written in `xvStockeRGBtoColormap` (lines 531, 556) and `xvinfo_` X11 path. Read in `initaccrochage`, `xvcouleur` (`colour = norgb[*icolor]`), `xvinfo_` (`gcvalues.foreground = norgb[foreground]`).
- **Fortran `CALL` sites touching `norgb`:** zero. `norgb` is never exposed to Fortran.
- **Qt ABI exposure:** zero. Not mentioned in `xvue/qt/include/xvue_qt_api.h`.
- **Why it exists in legacy C:** To hold the X11 `Pixel` value returned by `XAllocColor` in `TrueColor` mode. X11 drawing APIs take a Pixel, not RGB, so the color needed a lookup table. **Qt has no Pixel concept** — `QPainter::setPen/setBrush` take `QColor` directly.
- **Recommendation:** **Drop `norgb[]` entirely.** Do not add it as a dead storage field in `XvueState`. If the planner disagrees, the fallback is to keep it as an identity array inside `store_rgb_to_palette()` (`norgb[i] = i`) at zero runtime cost, but no Qt code path will ever read it.

---

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Qt 6 | 6.10.2 (installed on this machine via `libqt6core6t64` / `libqt6gui6` from Debian trixie) | Font loading, text metrics, text drawing, color representation | Already the backend for Phases 0-2; locked by project CONTEXT. |
| `QFontDatabase::addApplicationFont` | Qt 6.10.2 | Load a TTF at runtime, returns non-negative font ID on success | The single supported way to load an app-bundled TTF in Qt 6. `[VERIFIED: Qt 6.10 headers installed locally]` |
| `QFontMetrics::horizontalAdvance(QChar)` / `horizontalAdvance(QString)` | Qt 6 | Pixel-accurate text width | Replaces Qt 5 `QFontMetrics::width()` which is deprecated in Qt 5.11 and removed in Qt 6. `[VERIFIED: Qt 6 docs — width() is gone]` |
| `QFontMetrics::height()` | Qt 6 | ascent + descent + leading | Direct replacement for legacy `ascent_pol + descent_pol`. `[VERIFIED]` |
| `QPainter::drawText(int x, int y, const QString&)` | Qt 6 | Baseline-form text rendering | The `(x, y, str)` form matches `XDrawString(x, y, str)` semantics exactly (baseline origin). The rect-form `drawText(QRect, flags, str)` does NOT — it uses the rect's top-left. Do not confuse. `[VERIFIED: Qt 6 docs]` |
| `QPainter::TextAntialiasing` render hint | Qt 6 | Per-glyph AA for drawText | Must be enabled alongside `QPainter::Antialiasing` to get smooth glyphs. Without it, Qt uses the font's bitmap hinting path which looks crunchy on non-integer sizes. `[VERIFIED: Qt docs]` |
| `QColor::fromRgbF(float, float, float)` | Qt 6 | Construct QColor from 0..1 float triple | Exact match for the legacy `red/green/blue[i]` storage format. Round-trips via `QColor::redF()/greenF()/blueF()`. `[VERIFIED]` |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `Q_INIT_RESOURCE(xvue_fonts)` | Qt 6 | Force-link resource data from static library | **Required** when `libxvueqt.a` is static — Qt resources compiled into a static lib are not auto-initialized when the executable links against the lib, because there's no shared-object constructor. Must be called from a non-static function that is itself referenced by main (e.g., `XvueApp::ensure()`). See Pitfall 1. `[VERIFIED: Qt static-lib resource documentation]` |
| `qt_add_resources(xvueqt FONTS PREFIX "/xvue/qt/fonts" FILES fonts/DejaVuSansMono.ttf)` | CMake 3.21+ with Qt6 | Compile the TTF into `libxvueqt.a` | The modern Qt 6 form replaces `qt5_add_resources` and does not require a hand-written `.qrc` file for simple cases. `[VERIFIED: Qt 6 CMake docs]` |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `qt_add_resources` embedded TTF | External `install(FILES fonts/…)` + `QFontDatabase::addApplicationFont("/abs/path/DejaVuSansMono.ttf")` | Simpler; but adds a runtime file dependency and needs `$MEFISTO`-rooted path resolution. Rejected — CONTEXT.md D-02 locks the `:/xvue/qt/fonts/DejaVuSansMono.ttf` resource path. |
| Bundled Bitstream Vera copy | System `fonts-dejavu-core` via absolute path `/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf` | Zero repo size impact. Rejected — introduces a runtime system-package dependency which CONTEXT.md D-01 explicitly forbids for reproducibility. |
| `QFont::setPointSize(int)` after default ctor | `QFont("family", ptsize)` | Equivalent; D-04 picks the two-arg ctor form. Both work. |

**Installation / verification:**
```bash
# Verify local DejaVu availability for reference (not for runtime use).
ls /usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf  # already installed

# Copy into repo (Phase 3 Wave 0 action):
mkdir -p xvue/qt/fonts
cp /usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf xvue/qt/fonts/
git add xvue/qt/fonts/DejaVuSansMono.ttf
```

**Version verification:** Qt 6.10.2 confirmed via `dpkg -l libqt6core6t64` → `6.10.2+dfsg-6`. `horizontalAdvance` available since Qt 5.11 (2018), stable in Qt 6 since 6.0. `qt_add_resources` CMake command stable since Qt 6.0. `addApplicationFont` semantics unchanged across Qt 5.0 → 6.10 per Qt docs changelog.

---

## Architecture Patterns

### Recommended Project Structure

Phase 3 adds no new directories beyond `xvue/qt/fonts/`. Files touched:

```
xvue/qt/
├── CMakeLists.txt                  # add qt_add_resources block, add new TU if palette helpers split out
├── fonts/
│   └── DejaVuSansMono.ttf          # NEW — bundled TTF (~340 KB)
├── src/
│   ├── xvue_qt_api.cpp             # 8 real bodies replace 8 warn-once stubs; xvfond_ hack retired
│   ├── xvue_qt_app.cpp             # addApplicationFont one-shot + Q_INIT_RESOURCE guard
│   ├── xvue_qt_app.h               # font_id_ static member
│   ├── xvue_qt_canvas.cpp          # resizeEvent: add TextAntialiasing hint; ctor: setAutoFillBackground + WA_OpaquePaintEvent
│   ├── xvue_qt_state.cpp           # static palette def, palette_initialized_ guard, imposed_defaults table
│   ├── xvue_qt_state.h             # palette statics + current_font_ + current_font_size_idx_ + current_metrics_ + kFontSizes[]
│   └── (optional) xvue_qt_palette.cpp  # if the planner wants palette helpers in a separate TU
├── cmake/
│   └── (no change — verify_no_exec script needs regex updated to add qApp->palette guard)
└── include/
    └── xvue_qt_api.h               # UNCHANGED — no new extern "C" symbols. 57 symbols stays 57.
prpr/
└── xvtest0.f                       # extended after existing second-hold (line 91 area)
```

### Pattern 1: `qt_add_resources` for static library + `Q_INIT_RESOURCE` shim

**What:** Bundle the TTF into `libxvueqt.a` and force-init it at runtime.

**When to use:** Always, when the Qt library is `STATIC` and the resource must be visible. Our `xvueqt` target is declared `STATIC` in `xvue/qt/CMakeLists.txt:19`.

**Example:**
```cmake
# xvue/qt/CMakeLists.txt — add after the existing add_library(xvueqt STATIC ...) block
qt_add_resources(xvueqt xvue_fonts
    PREFIX "/xvue/qt/fonts"
    FILES
        fonts/DejaVuSansMono.ttf
)
# The resource name "xvue_fonts" becomes the Q_INIT_RESOURCE symbol.
```

```cpp
// xvue/qt/src/xvue_qt_app.cpp — inside XvueApp::ensure() or wherever the qapp is created
#include <QResource>
// Force-link the static-lib resource. MUST be at namespace scope or inside a
// non-member function; Q_INIT_RESOURCE is a macro that expands to a call to
// an extern function auto-generated by qt_add_resources.
static void xvue_init_resources() {
    Q_INIT_RESOURCE(xvue_fonts);
}
// Call xvue_init_resources() from XvueApp::ensure() before addApplicationFont.
```

```cpp
// xvue/qt/src/xvue_qt_app.cpp — font load one-shot
int XvueApp::font_id_ = -1;  // definition; declaration in xvue_qt_app.h

void XvueApp::load_bundled_font_() {
    if (font_id_ >= 0) return;  // idempotent guard (D-02)
    xvue_init_resources();  // force-link the static resource
    font_id_ = QFontDatabase::addApplicationFont(
        QStringLiteral(":/xvue/qt/fonts/DejaVuSansMono.ttf"));
    if (font_id_ < 0) {
        std::fprintf(stderr,
            "xvue-qt: FATAL — failed to load bundled DejaVuSansMono.ttf "
            "from Qt resource path\n");
        // Do not exit — fall back to the platform default mono font.
    }
}
```

### Pattern 2: Class-static palette initialized on first XvueState ctor

**What:** Palette lives as `static` members of `XvueState`, default-loaded once per process.

**Example:**
```cpp
// xvue/qt/src/xvue_qt_state.h  (additions only)
struct XvueState {
    // ... existing Phase 2 fields (do NOT reorder) ...

    // Phase 3 (D-03): font size table and per-instance current font.
    static constexpr int kFontSizes[] = {8, 10, 12, 14, 16, 18, 20, 24, 28, 32};
    static constexpr int kNbFonts = sizeof(kFontSizes)/sizeof(int);
    QFont        current_font_;                  // rebuilt each xvchargefonte_
    int          current_font_size_idx_ = 2;     // default 12pt
    QFontMetrics current_metrics_{current_font_};

    // Phase 3 (D-09, D-10, D-11, D-12): process-lifetime palette.
    static float   red[256];
    static float   green[256];
    static float   blue[256];
    static QColor  palette_cache_[256];
    static bool    palette_cache_dirty_[256];
    static bool    palette_initialized_;

    // ... existing destructor + applyPen() ...
};
```

```cpp
// xvue/qt/src/xvue_qt_state.cpp  (additions only)
float   XvueState::red[256]   = {};
float   XvueState::green[256] = {};
float   XvueState::blue[256]  = {};
QColor  XvueState::palette_cache_[256];
bool    XvueState::palette_cache_dirty_[256];
bool    XvueState::palette_initialized_ = false;

// Internal C++ helper — NOT an extern "C" ABI symbol.
// Verbatim lift of the 16-color block from xvue/xvuelc.c:378-461.
static void imposed_defaults_fill(int n1coel) {
    // black
    XvueState::red[n1coel]   = 0.0f;
    XvueState::green[n1coel] = 0.0f;
    XvueState::blue[n1coel]  = 0.0f;
    // red
    XvueState::red[n1coel+1]   = 1.0f;
    XvueState::green[n1coel+1] = 0.0f;
    XvueState::blue[n1coel+1]  = 0.0f;
    // ... (16 total entries — copy verbatim from xvuelc.c:378-461)
    // IMPORTANT: preserve the exact float constants like 50./256., 220./256.
    // etc. — these are the legacy values and any drift breaks color matching.
}

static void palette_init_once() {
    if (XvueState::palette_initialized_) return;
    imposed_defaults_fill(0);  // legacy n1coel = 0 per default init
    // indices 16..255 initialized to black by the global zero-init above.
    for (int i = 0; i < 256; ++i) {
        XvueState::palette_cache_[i] = QColor::fromRgbF(
            XvueState::red[i], XvueState::green[i], XvueState::blue[i]);
        XvueState::palette_cache_dirty_[i] = false;
    }
    XvueState::palette_initialized_ = true;
}

XvueState::XvueState() /* : existing-initializers */ {
    palette_init_once();  // D-11 guard
}
```

### Pattern 3: Collapsed text drawing body (same as 02/D-09 trait pattern)

**Example:**
```cpp
// xvue_qt_api.cpp — single shared body
static void xvue_qt_draw_text_common(const char* string, int length, int x1, int y1) {
    if (!string || length <= 0) return;  // null-guard symmetry (WR-03)
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive() || !st->backing_) return;

    QString qstr = QString::fromLatin1(string, length);
    st->painter_->setFont(st->current_font_);        // cheap no-op if already set
    st->painter_->drawText(x1, y1, qstr);            // baseline form (D-06)

    // Phase 2 D-01 epilogue
    if (win->canvas()) win->canvas()->update();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

void proc(xvtexte)(char string[], int *length, int *x1, int *y1) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!length || !x1 || !y1) return;
    xvue_qt_draw_text_common(string, *length, *x1, *y1);
}

void proc(xvftexte)(char string[], int *length, int *x1, int *y1) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!length || !x1 || !y1) return;
    xvue_qt_draw_text_common(string, *length, *x1, *y1);
}
```

### Anti-Patterns to Avoid

- **Using `QPainter::drawText(QRect, flags, str)` instead of `drawText(int,int,str)`** — the rect form does NOT use baseline origin. Legacy `XDrawString` uses baseline; using the rect form would shift all labels by the ascent and break visual parity.
- **Using `QFontMetrics::width()`** — deprecated in Qt 5.11, removed in Qt 6. Use `horizontalAdvance`.
- **Setting `QPainter::Antialiasing` alone without `TextAntialiasing`** — text glyphs will use a different rendering path and look crunchy next to AA geometry. Both hints must be reasserted together after every `painter->begin()` (which is already done once per resize, per 02/D-22).
- **Putting `QFontDatabase::addApplicationFont` in a constructor that may run multiple times** — `font_id_` will leak into fragmented IDs. D-02's static-int guard is correct.
- **Reading `qApp->palette()` or `this->palette()` in `xvue_qt_canvas.cpp` or `xvue_qt_api.cpp`** — breaks D-19 invariant. Verified by the new `verify_no_exec` grep.
- **Sourcing `QColor` from `Qt::GlobalColor` enum at draw time for palette indices** — the legacy values are NOT the Qt globals (Qt::red is `#ff0000`; legacy "rouge" is also `#ff0000` but legacy "vert sombre" is `50,200,50` which does NOT match `Qt::darkGreen`). Always look up via `palette_cache_[idx]`.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Font file loading | Custom TTF parser / FreeType direct call | `QFontDatabase::addApplicationFont` | Qt already owns FreeType internally and handles font table registration into the application's private font database. |
| Text pixel extent | Custom glyph advance summing | `QFontMetrics::horizontalAdvance(QString)` | Qt handles kerning, complex shaping, and ligatures. Hand-rolling advance-sum misses kerning. |
| Text rendering | Per-glyph `drawPixmap` of a cached atlas | `QPainter::drawText(x, y, str)` | Qt's drawText path uses Freetype's grayscale/subpixel renderer and matches `QFontMetrics` exactly. |
| Color → pixel lookup | Manual `XAllocColor`-style RGB-to-index cache with `norgb[]` | `QColor::fromRgbF` passed directly to `setPen/setBrush` | Qt has no Pixel concept. The legacy `norgb` table is an X11 artifact with no Qt analog. |
| Qt resource linking in static libs | Hand-rolled hex-encoded byte array `#include` of the TTF | `qt_add_resources` + `Q_INIT_RESOURCE` | Qt's tooling generates the resource-init function; the `Q_INIT_RESOURCE` call is the documented static-lib escape hatch. |

**Key insight:** Every single thing Phase 3 needs already exists as a Qt 6 API. The temptation with the palette is to preserve the full legacy indirection (`norgb` → X11 Pixel → `XSetForeground`), but Qt's `setPen(QColor)` short-circuits that entire chain. Phase 3 is straight translation, not reinvention.

---

## Common Pitfalls

### Pitfall 1: Qt resources invisible when linked from a static library

**What goes wrong:** After `qt_add_resources(xvueqt ...)`, calling `addApplicationFont(":/xvue/qt/fonts/DejaVuSansMono.ttf")` returns `-1` at runtime.

**Why it happens:** Qt resources in a static lib are compiled into a translation unit whose symbols are not referenced by any caller in the executable. The linker silently drops them. Unlike dynamic libs, there's no shared-object constructor that would auto-run and register the resource.

**How to avoid:** Call `Q_INIT_RESOURCE(xvue_fonts)` from `XvueApp::ensure()` (or the XvueApp ctor) — this generates a reference to the auto-generated `qInitResources_xvue_fonts()` function, forcing the linker to pull in the resource TU. Alternatively, wrap in an inline function `xvue_init_resources()` as shown in Pattern 1.

**Warning signs:** `font_id_ == -1` after `addApplicationFont`; text renders in the platform's fallback font (usually Liberation Mono or DejaVu Sans, but not necessarily the exact file we bundled — defeating reproducibility).

### Pitfall 2: Palette default-init ordering breaks `xvrecuprgbdec_` at init

**What goes wrong:** `xvue/xvinit.f:143` calls `XVRECUPRGBDEC(NBCOLO, PROUGE, PVERT, PBLEU)` immediately after `XVINFO` — **before** any `xvcouleur_` or `xvchargefonte_` call. If the palette is lazily initialized on first `xvcouleur_` call, `xvrecuprgbdec_` returns zeros and the Fortran-side color arrays (`PROUGE/PVERT/PBLEU`) are all black, which downstream breaks every color lookup.

**Why it happens:** The legacy X11 flow was `xvinfo_ → XStoreColors/XAllocColor → populates red/green/blue[]` before `xvinfo_` returned. Fortran's `XVRECUPRGBDEC` call after `XVINFO` saw a populated palette. Lazy init in Phase 3 breaks this invariant.

**How to avoid:** Trigger `palette_init_once()` either (a) in `XvueState` ctor which runs inside `XvueWindow` ctor which runs inside `xvinitgraphique_` (earlier than xvinfo), **or** (b) at the top of the new `xvinfo_` body, before populating its output ints. Option (a) is cleaner and less coupled; already covered by D-11. Document this ordering constraint in the plan task for the xvinfo_ rewrite.

**Warning signs:** Running `pp/ppxvtest1_qt` produces all-black geometry after Phase 3 lands; running `XVRECUPRGBDEC` under a debugger right after `XVINFO` shows zero floats in PROUGE/PVERT/PBLEU.

### Pitfall 3: `painter_->setFont()` is NOT cheap across every `drawText` call

**What goes wrong:** CONTEXT.md D-06 says `setFont` is a "cheap Qt no-op if already set". This is slightly optimistic: Qt does compare the new QFont against the current one via `QFont::operator==`, but the comparison is a full member-by-member compare including unused font flags. For thousands of `xvtexte_` calls per frame (typical in node-number labeling) the repeated setFont adds measurable overhead.

**Why it happens:** `QFont::operator==` isn't O(1) in the pathological sense, but `setFont` also touches the painter's internal state machine.

**How to avoid:** Track the current font in XvueState (`current_font_size_idx_` already does this per D-04). Only call `painter_->setFont(current_font_)` inside `xvchargefonte_` (when the font changes), not inside every `xvtexte_` call. D-06's "cheap no-op" statement is still safe enough for Phase 3 — the pitfall is only measurable at 10k+ labels/sec which Phase 3 drivers do not hit.

**Warning signs:** Not likely to hit in Phase 3; flagged for future profiling work.

### Pitfall 4: `QString::fromLatin1` on a non-null-terminated Fortran `CHARACTER` slice

**What goes wrong:** Fortran passes `CHARACTER*K` strings as a char pointer + length. The pointer may not be null-terminated. Calling `QString::fromLatin1(ptr)` (one-arg form) reads until the first NUL, walking into arbitrary memory.

**Why it happens:** C functions called from Fortran expect the explicit length; the compiler does not add a trailing `\0`.

**How to avoid:** Always use the two-arg form `QString::fromLatin1(const char*, int length)`. D-05 already specifies this; keep it as an explicit code-review checklist item in the plan.

**Warning signs:** Random tail characters in rendered labels, valgrind reports of uninitialized reads in `QString` allocation.

### Pitfall 5: Dark-mode-freeze invariant is defensive — no positive test exists in Phase 3

**What goes wrong:** D-19/D-20 guarantee "the canvas never reads from `qApp->palette()`". But there is no Phase 3 test that **proves** dark-mode actually has zero effect on colormap output, because Phase 6 (the chrome phase) is what wires up dark-mode plumbing in the first place. Phase 3 can only verify the absence of the leak path via grep, not the absence of the leak via runtime behavior.

**Why it happens:** The invariant is forward-looking protection against a Phase 6 regression.

**How to avoid:** Accept this as a known limitation. The verification in Phase 3 is the grep (`verify_no_exec` extension). The empirical verification moves to Phase 6 when dark-mode is first toggleable. Document in 03-VALIDATION.md that TEXT-06 runtime verification is deferred to Phase 6, Phase 3 only satisfies the static-guarantee half.

**Warning signs:** Reviewer asking "where's the runtime test?" — answer is "Phase 6; Phase 3 provides the grep guard only".

### Pitfall 6: `xvchargefonte_` with `*nofont0 > 0` in the legacy code path

**What goes wrong:** The legacy `xvchargefonte_` body (`xvuelc.c:1489-1494`) calls `XFreeFont(display_mef, struc_police)` **only if** `*nofont0 > 0`. In Qt, we don't allocate an X11 resource for the current font — `QFont` is a value type. If the Qt body mirrors the free-if-nonzero pattern blindly, it does nothing (QFont destruction is automatic). Benign, but confusing in code review.

**How to avoid:** Ignore the `*nofont0` parameter entirely in the Qt body (it's only passed so the legacy X11 path can free the previous font). Document with a comment: `// *nofont0 unused — QFont ownership is automatic (RAII), no explicit free needed`.

### Pitfall 7: `xvftexte_` vs `xvtexte_` legacy target divergence (NOT a bug, but counter-intuitive)

**What goes wrong:** In legacy C (xvuelc.c:1642-1698), `xvftexte_` draws to `fenetre_mef` (the live X window!) while `xvtexte_` draws to `mempx` (the backing pixmap). These are NOT semantically identical in legacy — `xvftexte_` was an immediate-mode "flush text to screen" variant. CONTEXT.md D-06 collapses them on the assumption that under the Qt single-backing-pixmap model the distinction vanishes, which is true — but only because Phase 2 removed the dual-buffer draw target split.

**Why it happens:** Legacy code had two draw targets; Qt Phase 2 has one. The collapse is legal under the current architecture but **only** because D-05/D-11 established the single backing pixmap.

**How to avoid:** Add a comment on both symbol bodies pointing at this research section (Pitfall 7) and D-06 / 02/D-05. If Phase 4 (pixmap save/restore) ever reintroduces multiple draw targets, the collapse must be revisited. Plan 03-* should include this as a note on the xvtexte_/xvftexte_ task.

**Warning signs:** Review comment "are these really the same now?" — answer is "yes, under Phase 2 single-backing; see Pitfall 7".

### Pitfall 8: `verify_no_exec` grep false positives on `.palette()` in `QMainWindow`

**What goes wrong:** A naive `grep 'palette()' xvue/qt/src/` will match `QMainWindow::palette()` calls if any future chrome code adds one. Phase 3 does not add any, but the grep pattern has to be scoped.

**How to avoid:** The D-19 grep is `qApp->palette|->palette\(\)` restricted to `xvue_qt_canvas.*` and `xvue_qt_api.*` (NOT the window or app TUs, which may legitimately touch palette for Phase 6 chrome). Scope the path in `verify_no_exec.sh` accordingly.

---

## Code Examples

### Example 1: Bundled font load (TEXT-01)

```cpp
// xvue/qt/src/xvue_qt_app.cpp
// Source: legacy xvuelc.c:981-1003 (semantic target); Qt 6 addApplicationFont docs.
#include <QFontDatabase>
#include <QStringLiteral>
#include <QResource>

int XvueApp::font_id_ = -1;  // declared static in xvue_qt_app.h

void XvueApp::load_bundled_font_() {
    if (font_id_ >= 0) return;  // D-02 idempotent guard
    Q_INIT_RESOURCE(xvue_fonts);  // Pitfall 1 — required for static libs
    font_id_ = QFontDatabase::addApplicationFont(
        QStringLiteral(":/xvue/qt/fonts/DejaVuSansMono.ttf"));
    if (font_id_ < 0) {
        std::fprintf(stderr,
            "xvue-qt: WARNING — bundled DejaVuSansMono.ttf failed to load; "
            "falling back to platform mono default\n");
    }
}
```

### Example 2: `xvchargefonte_` body (TEXT-01)

```cpp
// xvue/qt/src/xvue_qt_api.cpp
// Source: legacy xvuelc.c:1463-1574; D-04; Pitfall 6.
void proc(xvchargefonte)(int *nofont0, int *nofont, int *largpx, int *hautpx) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    (void)nofont0;  // Pitfall 6 — QFont ownership is automatic

    if (!nofont || !largpx || !hautpx) return;  // WR-03 null-guard

    int idx = *nofont;
    if (idx < 0) idx = 0;
    if (idx >= XvueState::kNbFonts) idx = XvueState::kNbFonts - 1;

    auto& win = XvueApp::window_slot();
    if (!win || !win->state()) return;
    auto* st = win->state();

    st->current_font_size_idx_ = idx;
    st->current_font_ = QFont(QStringLiteral("DejaVu Sans Mono"),
                              XvueState::kFontSizes[idx]);
    st->current_metrics_ = QFontMetrics(st->current_font_);

    if (st->painter_ && st->painter_->isActive()) {
        st->painter_->setFont(st->current_font_);
    }

    *largpx = st->current_metrics_.horizontalAdvance(QLatin1Char('0'));
    *hautpx = st->current_metrics_.height();
}
```

### Example 3: `xvnbpixeltexte_` body (TEXT-02)

```cpp
// Source: legacy xvuelc.c:1576-1600; D-05.
void proc(xvnbpixeltexte)(char *texte, int *length, int *nbpxla, int *nbpxha) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!texte || !length || !nbpxla || !nbpxha) return;

    auto& win = XvueApp::window_slot();
    if (!win || !win->state()) { *nbpxla = 0; *nbpxha = 0; return; }
    auto* st = win->state();

    QString qstr = QString::fromLatin1(texte, *length);  // Pitfall 4
    *nbpxla = st->current_metrics_.horizontalAdvance(qstr);
    *nbpxha = st->current_metrics_.height();
}
```

### Example 4: `xvcouleur_` body (TEXT-04)

```cpp
// Source: legacy xvuelc.c:1119-1185 (PostScript block dropped per 02/D-26); D-14.
void proc(xvcouleur)(int *icolor) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!icolor) return;

    int i = *icolor;
    if (i < 0 || i >= 256) { i = 1; /* match legacy "fallback to red" */ }

    auto& win = XvueApp::window_slot();
    if (!win || !win->state()) return;
    auto* st = win->state();

    if (XvueState::palette_cache_dirty_[i]) {
        XvueState::palette_cache_[i] = QColor::fromRgbF(
            XvueState::red[i], XvueState::green[i], XvueState::blue[i]);
        XvueState::palette_cache_dirty_[i] = false;
    }

    st->foreground_ = XvueState::palette_cache_[i];
    st->applyPen();  // 02/D-20 brush+pen sync
    // D-14: no flush. State-change only.
}
```

### Example 5: `xvrecuprgbdec_` + `xvactivervb_` (TEXT-05 ABI half)

```cpp
// Source: legacy xvuelc.c:1044-1068; D-18.
void proc(xvrecuprgbdec)(int *nbcolor, float *r, float *g, float *b) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!nbcolor || !r || !g || !b) return;
    int n = *nbcolor;
    if (n < 0) n = 0;
    if (n > 256) n = 256;
    for (int i = 0; i < n; ++i) {
        r[i] = XvueState::red[i];
        g[i] = XvueState::green[i];
        b[i] = XvueState::blue[i];
    }
}

// Source: legacy xvuelc.c:1072-1116 (internal xvStockeRGBtoColormap call and
//         norgb/XInstallColormap blocks dropped); D-17.
void proc(xvactivervb)(int *palcour, int *nbcells,
                       float r[], float g[], float b[]) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    (void)palcour;  // palcourc was a PostScript state var — dropped per 02/D-26
    if (!nbcells || !r || !g || !b) return;
    int n = *nbcells;
    if (n < 0) n = 0;
    if (n > 256) n = 256;
    for (int i = 0; i < n; ++i) {
        XvueState::red[i]   = r[i];
        XvueState::green[i] = g[i];
        XvueState::blue[i]  = b[i];
        XvueState::palette_cache_dirty_[i] = true;
    }
    // No painter touch — xvcouleur_ picks up on next call.
}
```

### Example 6: `xvinfo_` body (replaces zero-fill)

```cpp
// Source: legacy xvuelc.c:612-1042; D-22/D-23; Pitfall 2.
void proc(xvinfo)( int *ix, int *iy, int *maxfonts,
                   int *n1coref, int *ndcoref, int *n1coelf,
                   int *ndcoelf, int *n1coulf, int *ndcoulf, int *nbcolo,
                   char namefonts[][256], int nbchar[], int *nbfonts,
                   int *visuclass ) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();

    // Resize window (existing Phase 1 code)
    auto& win = XvueApp::window_slot();
    if (win && ix && iy) win->resize(*ix, *iy);

    // Palette ranges (match legacy xvinfo semantics: reserved/elementary/user)
    // Legacy n1core=0, ndcore=15 (first 16 reserved), n1coel=0, ndcoel=15
    // (the 16 imposed elementaries), n1coul=16, ndcoul=255 (user-modifiable).
    if (n1coref)  *n1coref  = 0;
    if (ndcoref)  *ndcoref  = 15;
    if (n1coelf)  *n1coelf  = 0;
    if (ndcoelf)  *ndcoelf  = 15;
    if (n1coulf)  *n1coulf  = 16;
    if (ndcoulf)  *ndcoulf  = 255;
    if (nbcolo)   *nbcolo   = 256;
    if (visuclass) *visuclass = 4;  // TrueColor (legacy enum value; Qt is always 24-bit)

    // Font list (D-22/D-23)
    static const char* const listfonts_qt[XvueState::kNbFonts] = {
        "DejaVu Sans Mono 8pt",
        "DejaVu Sans Mono 10pt",
        "DejaVu Sans Mono 12pt",
        "DejaVu Sans Mono 14pt",
        "DejaVu Sans Mono 16pt",
        "DejaVu Sans Mono 18pt",
        "DejaVu Sans Mono 20pt",
        "DejaVu Sans Mono 24pt",
        "DejaVu Sans Mono 28pt",
        "DejaVu Sans Mono 32pt",
    };
    int nfonts = XvueState::kNbFonts;
    if (maxfonts && *maxfonts < nfonts) nfonts = *maxfonts;
    if (maxfonts) *maxfonts = nfonts;
    if (nbfonts)  *nbfonts  = nfonts;

    for (int k = 0; k < nfonts; ++k) {
        std::strncpy(namefonts[k], listfonts_qt[k], 255);
        namefonts[k][255] = '\0';
        nbchar[k] = static_cast<int>(std::strlen(listfonts_qt[k]));
    }

    // Pitfall 2: palette must be initialized before xvrecuprgbdec_ fires.
    // XvueState::palette_init_once() already ran in the XvueState ctor, which
    // was called via XvueWindow ctor inside xvinitgraphique_. Defensive:
    if (!XvueState::palette_initialized_) {
        // Force-init as a safety net if someone ever calls xvinfo_ without an
        // open window (shouldn't happen, but the invariant is cheap).
        extern void palette_init_once();  // declared in xvue_qt_state.cpp
        palette_init_once();
    }
}
```

---

## Runtime State Inventory

*Phase 3 is primarily additive code work, not a rename/refactor. Items below cover the few runtime state touchpoints that exist.*

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | None — no databases, no persistent disk state, no ChromaDB/Mem0/sqlite. | None |
| Live service config | None — no external services. | None |
| OS-registered state | None — no task scheduler, no systemd unit, no launchd, no pm2. | None |
| Secrets/env vars | `QT_QPA_PLATFORM` consulted at runtime (D-28 prescribes `xcb` for A/B gate, `offscreen` for CI smoke). Not a secret; user-controlled env var. | None (document in 03-VALIDATION.md). |
| Build artifacts | `xvue/qt/build/` is a CMake build directory; after adding `qt_add_resources`, a stale build may carry the old resource symbol table. Also `pp/pp*_qt` binaries must be rebuilt after the CMakeLists change. `xvue/xvuelc.o` (legacy X11 build output) is unrelated and stays as-is. | Plan task must call `rm -rf xvue/qt/build && bin/cbl_tout_qt` at least once before the A/B gate. |

**Canonical question (what runtime state survives after every file is updated):** After the CMakeLists change adds `qt_add_resources`, the stale `xvue/qt/build/` may contain compiled objects from before the resource was added — the linker will succeed but `Q_INIT_RESOURCE(xvue_fonts)` will fail to link because the resource symbol doesn't exist. A clean rebuild is required exactly once. Nothing else.

---

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|-------------|-----------|---------|----------|
| Qt 6 Core/Gui/Widgets | libxvueqt.a build + all phase 3 code | ✓ | 6.10.2+dfsg-6 (apt) | — |
| gfortran / gcc | legacy `bin/cbl_tout` regression build | ✓ | (from CLAUDE.md, verified by existing phases 0-2 builds) | — |
| libX11-dev | legacy `bin/cbl_tout` regression build (not touched by Phase 3 code, but A/B gate requires it) | ✓ (assumed — Phase 2 completed a full legacy build) | — | — |
| CMake ≥ 3.21 | `qt_add_resources` CMake API | ✓ (Phase 0 already requires 3.21 in `xvue/qt/CMakeLists.txt:7`) | — | — |
| `fonts-dejavu-core` | Source of `DejaVuSansMono.ttf` to copy into repo | ✓ | 2.37-8 at `/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf` | — (file can be fetched from upstream DejaVu release) |
| XCB platform plugin | `QT_QPA_PLATFORM=xcb` for A/B gate (D-28) | ✓ (libqt6gui6 installs xcb plugin by default on Debian) | — | `offscreen` platform for CI-style runs (smoke only, not visual) |
| `testa/` baseline cases | 5-case A/B run (D-26) | ✓ (BUILD-10 listed them; exists in repo) | — | — |
| `prpr/xvtest1.f..xvtest4.f` drivers | A/B visual gate (D-26) | ✓ (exist in `prpr/`) | — | — |

**Missing dependencies with no fallback:** None.
**Missing dependencies with fallback:** None.

---

## Validation Architecture

### Test Framework

| Property | Value |
|----------|-------|
| Framework | MEFISTO legacy Fortran drivers (`prpr/xvtest*.f`, `testa/*`) + shell build scripts. Inherits Phase 2's setup. |
| Config file | none — `bin/cbl_tout`, `bin/cbl_tout_qt`, `bin/cbxvtest0_qt` scripts drive the chain. |
| Quick run command | `cmake --build xvue/qt/build && bin/cbl_tout_qt` (~60 s) |
| Full suite command | `bin/cbl_tout_qt && bin/cbl_tout && pp/ppxvtest0_qt && QT_QPA_PLATFORM=xcb pp/ppxvtest1_qt …` |

### Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| TEXT-01 | `xvchargefonte_` loads bundled font, returns pixel width/height | build + driver + grep | `bin/cbl_tout_qt && grep -q 'DejaVuSansMono' xvue/qt/CMakeLists.txt && grep -q 'addApplicationFont' xvue/qt/src/xvue_qt_app.cpp && pp/ppxvtest0_qt` | ❌ Wave 0 (CMake block + xvue_qt_app.cpp edit + ttf file) |
| TEXT-02 | `xvnbpixeltexte_` extents via `horizontalAdvance`/`height`; no clipping on testa labels | build + grep + driver | `grep -q 'horizontalAdvance' xvue/qt/src/xvue_qt_api.cpp && pp/ppxvtest0_qt` + human visual on `testa/nafems_le1` and `testa/pan2d` | ❌ Wave 0 (xvtest0.f extension draws a label whose width was computed via xvnbpixeltexte_) |
| TEXT-03 | `xvtexte_`/`xvftexte_` render via `drawText` | build + grep + driver | `grep -q 'drawText' xvue/qt/src/xvue_qt_api.cpp && pp/ppxvtest0_qt` (xvtest0 extension draws 8 color+label pairs) + visual | ❌ Wave 0 |
| TEXT-04 | `xvcouleur_` sets pen+brush via indexed palette | build + grep + driver | `grep -q 'palette_cache_' xvue/qt/src/xvue_qt_api.cpp && grep -q 'XvueState::red' xvue/qt/src/xvue_qt_state.cpp && pp/ppxvtest0_qt` + visual (8 distinct colored lines) | ❌ Wave 0 |
| TEXT-05 | palette populate/query/activate within 1 bit/channel vs X11 | build + grep + A/B | `grep -q 'xvrecuprgbdec' xvue/qt/src/xvue_qt_api.cpp && grep -q 'xvactivervb' xvue/qt/src/xvue_qt_api.cpp` + full A/B gate on xvtest1..4 + 5 testa cases | ❌ Wave 0 (xvtest0.f extension calls xvactivervb_ with custom triple) |
| TEXT-06 | QPalette changes do not affect backing pixmap | grep (static guarantee) + deferred runtime (Phase 6) | `! grep -rE 'qApp->palette\|->palette\(\)' xvue/qt/src/xvue_qt_canvas.* xvue/qt/src/xvue_qt_api.*` added to `verify_no_exec.sh` | ❌ Wave 0 (update verify_no_exec.sh) |

### Per-Task Verification (Nyquist-sampling design)

**Sampling rate:** every task commit triggers `cmake --build xvue/qt/build` (~15 s incremental). Every wave merge triggers `bin/cbl_tout_qt && pp/ppxvtest0_qt` (~60 s). The phase gate runs the full A/B on both backends.

**Per-task checks — what runs, what it depends on:**

| Check class | Tool | What it verifies | Depends on research finding |
|-------------|------|------------------|------------------------------|
| Build-green — Qt backend | `cmake --build xvue/qt/build` or `bin/cbl_tout_qt` | Every Phase 3 task keeps `libxvueqt.a` linkable | Standard Stack table (Qt 6.10.2 verified); Pattern 1 (qt_add_resources syntax); Environment Availability row 4 |
| Build-green — legacy X11 backend | `bin/cbl_tout` | Phase 3 does not break the regression path | CLAUDE.md §"Compilation must never break"; VALID-02 |
| verify_abi — 57 symbols unchanged | `bin/verify_abi` (existing Phase 0 target) | Phase 3 adds zero new Fortran ABI symbols (call-site audit conclusion) | Call-Site Audit section (only 8 bodies rewritten, no new symbols); `xvCouleursImposees`/`xvStockeRGBtoColormap`/`xvColormapToRGB` are C++ helpers, not ABI |
| verify_no_exec — palette leak guard | `bin/verify_no_exec` (existing Phase 1 target) updated with new regex | No `qApp->palette` / `->palette()` in canvas or api TUs (D-19) | Pitfall 8 (grep scoping) — verify_no_exec.sh path arg must be `xvue/qt/src/xvue_qt_canvas.*` and `xvue/qt/src/xvue_qt_api.*` specifically |
| Warn-once drain | `pp/ppxvtest0_qt 2>&1 \| grep -c 'xvue-qt: .*(xvtexte\|xvftexte\|xvchargefonte\|xvnbpixeltexte\|xvcouleur\|xvactivervb\|xvrecuprgbdec\|xvinfo).*not implemented' ` → expected 0 after Phase 3 complete | No warn-once stubs left for Phase 3 entry points | Current Qt stubs inventory (xvue_qt_api.cpp:180-440) — 8 warn-once lines must go to 0 |
| Font resource present | `grep -q 'qt_add_resources.*xvue_fonts' xvue/qt/CMakeLists.txt && test -f xvue/qt/fonts/DejaVuSansMono.ttf` | TTF is bundled and wired | Pattern 1; Pitfall 1 |
| Font actually loads at runtime | `pp/ppxvtest0_qt 2>&1 \| grep -q 'bundled DejaVuSansMono.ttf failed to load' && exit 1 \|\| exit 0` | `addApplicationFont` returned non-negative ID | Pitfall 1; Example 1 |
| Literal ABI diff | `diff <(grep -E 'proc\\(xv(chargefonte\|texte\|ftexte\|nbpixeltexte\|couleur\|activervb\|recuprgbdec)\\)' xvue/qt/include/xvue_qt_api.h) <(grep -E 'proc\\(xv(chargefonte\|texte\|ftexte\|nbpixeltexte\|couleur\|activervb\|recuprgbdec)\\)' xvue/xvuelc.c)` → zero diff on the signature lines (ignoring body) | Literal ABI parity preserved byte-for-byte per 00/D-33 | Call-site audit: `xvchargefonte(int*, int*, int*, int*)` confirmed identical between legacy and Qt stub (xvue_qt_api.cpp:381) |
| xvtest0 extension smoke | `pp/ppxvtest0_qt` exit 0 with extended Phase 3 section emitting 8 colored lines + labels | Integrated path works end-to-end | D-24; xvtest0.f Phase 2 two-hold structure (current :33-97) |
| A/B catch-up gate — visual | For each of {xvtest1, xvtest2, xvtest3, xvtest4, 5 testa cases}: run with `QT_QPA_PLATFORM=xcb pp/pp{driver}_qt` and legacy `pp/pp{driver}`, visually compare | No clipping/overlap/missing glyphs/miscolored regions/missing geometry; coordinate drift ≤2 px | D-26/D-27/D-28; call-site audit confirms these drivers exercise every Phase 3 entry point |
| Dark-mode-freeze runtime (DEFERRED) | n/a in Phase 3 | Runtime toggle of system dark-mode verifies colormap pixels unchanged | Pitfall 5 — deferred to Phase 6; Phase 3 provides grep guard only |

**Sampling continuity check:** The plan should have at most 2 consecutive manual (visual) tasks. The natural layout is Wave 0 (setup: CMake + ttf + state.h grow, all headless), Wave 1 (per-entry-point bodies, each gated by grep + xvtest0 smoke, all headless), Wave 2 (xvtest0 extension + visual checkpoint, manual), Wave 3 (A/B catch-up + visual comparison across 4 xvtest drivers + 5 testa cases, manual). Wave 2 and Wave 3 are consecutive manual; insert an automated clean-rebuild sweep between them (same pattern as Phase 2's 02-04-01 clean rebuild + 02-04-02 visual).

### Wave 0 Gaps

- [ ] `xvue/qt/fonts/DejaVuSansMono.ttf` — file does not exist in repo yet. Copy from `/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf`. License note in a neighboring `xvue/qt/fonts/LICENSE.dejavu` file.
- [ ] `xvue/qt/CMakeLists.txt` — add `qt_add_resources(xvueqt xvue_fonts PREFIX "/xvue/qt/fonts" FILES fonts/DejaVuSansMono.ttf)` block.
- [ ] `xvue/qt/cmake/verify_no_exec.sh` — add `qApp->palette|->palette\(\)` regex scoped to canvas + api TUs (D-19).
- [ ] `xvue/qt/src/xvue_qt_state.h` — add palette static members, `kFontSizes[]` table, `current_font_`/`current_metrics_`/`current_font_size_idx_`.
- [ ] `xvue/qt/src/xvue_qt_state.cpp` — add palette static definitions, `imposed_defaults_fill` verbatim-lifted table from `xvuelc.c:378-460`, `palette_init_once()` guarded call from XvueState ctor.
- [ ] `xvue/qt/src/xvue_qt_app.h` + `.cpp` — add `static int font_id_`, `load_bundled_font_()`, `Q_INIT_RESOURCE(xvue_fonts)` call, call from `XvueApp::ensure()`.
- [ ] `xvue/qt/src/xvue_qt_canvas.cpp` — ctor: `setAutoFillBackground(false)` + `setAttribute(Qt::WA_OpaquePaintEvent, true)` (D-20); resizeEvent: add `QPainter::TextAntialiasing` hint next to `QPainter::Antialiasing` (D-07).
- [ ] `prpr/xvtest0.f` — extend after line 97 with TEXT coverage section per D-24. Use the existing SLEEP holds if possible.

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `QFontMetrics::width(const QString&)` | `QFontMetrics::horizontalAdvance(const QString&)` | Qt 5.11 (2018); deprecated in 5.11, removed in Qt 6.0 | Code using `width()` will not compile on Qt 6. D-04/D-05 correctly specify `horizontalAdvance`. |
| `qt5_add_resources(target files)` | `qt_add_resources(target resource_name PREFIX ... FILES ...)` | Qt 6.0 | Phase 3 must use the Qt 6 signature. The `resource_name` is what `Q_INIT_RESOURCE` takes. |
| Xlib `XDrawString` baseline text | `QPainter::drawText(int x, int y, const QString&)` | N/A (legacy replacement) | Baseline-form drawText is the closest Qt analog; see anti-patterns about rect-form. |
| X11 `XListFonts` + `XLoadQueryFont` dynamic font discovery | `QFontDatabase::addApplicationFont` with a bundled TTF | Phase 3 | Reproducibility across machines; no dependency on server-side font config. |
| X11 `XAllocColor` / `XQueryColor` / `Pixel` indirection via `norgb[]` | Direct `QColor::fromRgbF` → `QPainter::setPen/setBrush` | Phase 3 | No Pixel concept; `norgb` becomes dead storage or is dropped. |
| X11 `XFreeFont` explicit free in `xvchargefonte_` | `QFont` RAII value type | Phase 3 | `nofont0` param is ignored in Qt body (Pitfall 6). |

**Deprecated/outdated in legacy code being dropped:**
- `xvuelc.c:481-500` `XQueryColor` → float conversion — replaced by direct float storage (D-12).
- `xvuelc.c:525-560` `XStoreColors`/`XAllocColor` dual-path — Qt has no PseudoColor path, single TrueColor code path only.
- `xvuelc.c:1510-1565` PostScript font-name mangling inside `xvchargefonte_` — PostScript emission moved to Phase 7 per 02/D-26; strip entirely from Phase 3 body.

---

## Project Constraints (from CLAUDE.md)

- **Compilation must never break.** Phase 3 must keep **both** `bin/cbl_tout_qt` (Qt backend) and `bin/cbl_tout` (legacy X11 backend) green after every task commit. The legacy X11 path is untouched by Phase 3 code but must not be damaged by any accidental CMake/include-file leakage.
- **Fortran fixed-form column 7+ norms** (per `doc/normes.ps`). The `prpr/xvtest0.f` extension must follow the existing column convention — look at the existing two-hold block (lines 33-97) as the template.
- **Small tests must continue to pass** after every change. For Phase 3 the "small test" is `pp/ppxvtest0_qt`; the "large tests" are the 5 `testa/` baseline cases — use them only at wave boundaries per the Nyquist sampling, not per task.
- **Ask before installing system packages.** Phase 3 does not need any new apt package — Qt 6 is already installed (6.10.2), DejaVu TTF is already installed at `/usr/share/fonts/`. If a reviewer later requests a different bundled font, ask before apt-installing.
- **Git discipline:** commit after every logical step. Phase 2 used a per-wave-per-task commit cadence; Phase 3 should follow the same pattern. Never force-push, never bypass hooks.

---

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | `xvCouleursImposees`, `xvStockeRGBtoColormap`, `xvColormapToRGB` are NOT extern "C" Fortran entry points and should be implemented as C++ file-static helpers. CONTEXT.md D-15/D-16 describes them as entry points needing real bodies, which contradicts Phase 0's deliberate decision (00-02-SUMMARY.md:43: "Planner Alert Option A honored: 57 Fortran-facing entries in the public header; xvCouleursImposees_, xvColormapToRGB_, xvStockeRGBtoColormap_ are deliberately NOT declared"). `[VERIFIED: grep across all .f/.F files — zero call sites; grep across xvue_qt_api.h — symbols not declared; grep across xvue_qt_api.cpp — no stubs for these names]`. The behavioral test in TEXT-05 (palette values within 1 bit/channel of X11) is satisfied by the internal-helper implementation. | Call-Site Audit / User Constraints locked decisions D-15/D-16 | HIGH if wrong — would add 3 ABI symbols that violate 00/D-33 literal-parity guarantee and break `bin/verify_abi`. Research is confident: needs user confirmation that Phase 0's decision still stands. |
| A2 | `xvtexte_` and `xvftexte_` can be collapsed into a single body under the Qt single-backing-pixmap model even though legacy xvuelc.c draws them to DIFFERENT targets (`mempx` vs `fenetre_mef`). `[VERIFIED: xvuelc.c:1658 vs :1678]`. The collapse is CONTEXT.md D-06's intent and is legal because Phase 2 already unified the draw target. | User Constraints D-06 / Pitfall 7 | LOW — if the collapse is wrong in practice, Phase 4 (pixmap save/restore) would catch it when reintroducing secondary draw targets. Flagging for reviewer awareness. |
| A3 | `xvactivervb_(palcour, nbcells, r[], g[], b[])` — 4 args, writes `r/g/b[]` into the palette and marks indices dirty. CONTEXT.md D-17 says the body constructs a transient QColor and writes directly to `state_->foreground_` WITHOUT touching the palette arrays. `[VERIFIED against legacy xvuelc.c:1072-1116]` — the legacy body copies `r/g/b` into the static arrays AND stores through `xvStockeRGBtoColormap` (which also writes the arrays). The legacy ABI signature is `(int *palcour, int *nbcells, float r[], float g[], float b[])` — it receives a full palette block, not a single color triple. D-17's "transient QColor" interpretation is wrong; D-17's "does NOT write into red/green/blue[]" is also wrong. Recommend: implement the Example 5 body (bulk copy into static palette, mark dirty) and supersede D-17. | User Constraints D-17 / Example 5 | HIGH — if the plan implements D-17 as written, `xvue/palcde.f:619`'s user-defined colormap load will silently drop every color past the first one and break palette activation. Research is confident the correction is needed. |
| A4 | Bundled DejaVu Sans Mono 2.37-8 (from `fonts-dejavu-core`) is license-compatible with repository distribution under the Bitstream Vera Fonts Copyright + DejaVu modifications. `[CITED: https://dejavu-fonts.github.io/License.html — "Bitstream Vera Fonts Copyright ... Permission is hereby granted, free of charge, to any person obtaining a copy of the fonts accompanying this license ..."]`. Confirmed redistributable. | D-01 | LOW — license has been in place since 2003; no known incompatibility with project license. |
| A5 | `Q_INIT_RESOURCE(xvue_fonts)` is required because `libxvueqt.a` is STATIC. `[CITED: Qt 6 documentation on QResource and static linking, also Qt wiki static-lib FAQ]`. | Pattern 1 / Pitfall 1 | MEDIUM — if omitted, `addApplicationFont` silently returns `-1` and text falls back to platform default. Detected by the Example 1 stderr warning. Low cost to get wrong once. |
| A6 | `nbpolices_qt = 10` (kFontSizes has 10 entries: 8,10,12,14,16,18,20,24,28,32). `[VERIFIED against D-03 table — 10 rows]`. | D-03 / D-22 | LOW — if the count is off by one, `xvinfo_` returns the wrong `nbfonts` and Fortran `NMFONT/NBCAFO` arrays have one empty row. Easily caught in xvtest0 smoke. |

**Summary:** 6 assumptions logged. **A1, A3 are HIGH risk** and directly contradict CONTEXT.md text — the planner must raise them to the user before locking the plan. A2 is flagged for reviewer awareness. A4-A6 are low/medium risk and are embedded in normal research tradecraft.

---

## Open Questions (RESOLVED)

**RESOLVED 2026-04-11:** All four questions were confirmed by the user during /gsd-plan-phase. CONTEXT.md D-13, D-15, D-16, D-17 amended in place to reflect the research recommendations. See CONTEXT.md `[AMENDED 2026-04-11 post-research]` markers.

1. **Does the user want `xvCouleursImposees` / `xvStockeRGBtoColormap` / `xvColormapToRGB` to stay as internal C++ helpers (Phase 0's decision) or to be promoted to ABI symbols as CONTEXT.md D-15/D-16 suggests?** — **RESOLVED: file-static C++ helpers. See CONTEXT.md D-15/D-16 AMENDED.** Phase 3 adds ZERO new public ABI symbols; 57-symbol count preserved per 00/D-33.
   - What we know: Zero Fortran callers, Phase 0 deliberately excluded them from the public header, 00-02-SUMMARY.md:43 documents the decision.
   - What's unclear: CONTEXT.md was written later and appears to assume they are entry points.
   - Recommendation: Confirm Phase 0's decision still holds. Implement as file-static C++ helpers (no ABI impact). If user wants them promoted, that's a separate 3-symbol ABI expansion that needs `bin/verify_abi` count update (57 → 60) and `xvue/qt/include/xvue_qt_api.h` additions.

2. **Should the `xvactivervb_` body bulk-write the palette (Example 5) or transient-color-only (D-17)?** — **RESOLVED: bulk palette-load. See CONTEXT.md D-17 AMENDED.** Copy `r/g/b[0..nbcells-1]` into `red/green/blue[]`, mark dirty, no `foreground_` touch, no painter.

3. **Does `norgb[]` get dropped entirely, or kept as identity-stub inside an internal helper?** — **RESOLVED: dropped entirely. See CONTEXT.md D-13 AMENDED.** Zero callers, no field in XvueState.

4. **Wave 0 ordering — does the palette `palette_init_once()` run before or after the first `xvrecuprgbdec_` call?** — **RESOLVED: palette_init_once() fires inside XvueState ctor (called from XvueWindow ctor inside xvinitgraphique_), strictly before xvinfo_ returns and therefore before xvue/xvinit.f:143's XVRECUPRGBDEC call. See Pitfall 2 architecture and 03-01-PLAN.md Task 2.**
   - What we know: `xvue/xvinit.f:143` calls `XVRECUPRGBDEC` immediately after `XVINFO`. Palette must be populated by then.
   - Recommendation: Invoke `palette_init_once()` in `XvueState` ctor which runs inside `XvueWindow` ctor which runs inside `xvinitgraphique_` — well before `xvinfo_` is called. Covered by D-11 and Pitfall 2. Document explicitly in the plan task.

---

## Sources

### Primary (HIGH confidence)
- `xvue/xvuelc.c:85` — `#define CMAPSIZE 256` (palette size source of truth)
- `xvue/xvuelc.c:100-104` — static float red/green/blue[] + norgb[] storage layout
- `xvue/xvuelc.c:126-133` — struc_police + nbpolices + listfonts mirror targets
- `xvue/xvuelc.c:358-461` — `xvCouleursImposees` body = 16-color imposed defaults table (verbatim source for Pattern 2)
- `xvue/xvuelc.c:463-500` — `xvColormapToRGB` body (read path)
- `xvue/xvuelc.c:503-559` — `xvStockeRGBtoColormap` body (write path)
- `xvue/xvuelc.c:612-1042` — `xvinfo_` body (palette/font enumeration reference)
- `xvue/xvuelc.c:1044-1068` — `xvrecuprgbdec` body (Example 5 source)
- `xvue/xvuelc.c:1072-1116` — `xvactivervb` body (A3 correction source)
- `xvue/xvuelc.c:1119-1185` — `xvcouleur` body (Example 4 source; PostScript tail dropped)
- `xvue/xvuelc.c:1463-1574` — `xvchargefonte` body (Example 2 semantic source)
- `xvue/xvuelc.c:1576-1600` — `xvnbpixeltexte` body (Example 3 source)
- `xvue/xvuelc.c:1642-1698` — `xvftexte`/`xvtexte` bodies (Pitfall 7 dual-target source)
- `xvue/qt/src/xvue_qt_api.cpp:180-440` — current Qt stubs to replace
- `xvue/qt/src/xvue_qt_state.h` — Phase 2 state layout (additive growth target)
- `xvue/qt/CMakeLists.txt` — current CMake config (qt_add_resources insertion point)
- `xvue/chargefonte.f:36` — confirmed `CALL XVCHARGEFONTE(NOFONT0, NOPOCA, NPLACA, NPHACA)` 4-arg signature
- `xvue/xvinit.f:143` — confirmed `CALL XVRECUPRGBDEC(NBCOLO, PROUGE, PVERT, PBLEU)` — drives Pitfall 2
- `xvue/palcde.f:619` — confirmed `CALL XVACTIVERVB(NOPALC, NDCOUL+1, PROUGE, PVERT, PBLEU)` — drives A3
- `.planning/phases/00-build-skeleton-abi-stubs/00-02-SUMMARY.md:43` — "Planner Alert Option A honored: 57 Fortran-facing entries...; xvCouleursImposees_, xvColormapToRGB_, xvStockeRGBtoColormap_ are deliberately NOT declared" — drives A1
- `.planning/phases/02-drawing-primitives-backing-pixmap/02-CONTEXT.md` — D-01 epilogue, D-05 long-lived painter, D-09 collapse pattern, D-20 brush sync, D-22 AA reapply, D-26 PS-state drop
- `.planning/phases/02-drawing-primitives-backing-pixmap/02-VALIDATION.md` — template for 03-VALIDATION.md Per-Task Verification Map
- Qt 6.10.2 installed locally (`dpkg -l libqt6core6t64` → `6.10.2+dfsg-6`) — verified `horizontalAdvance`, `addApplicationFont`, `qt_add_resources`, `Q_INIT_RESOURCE` all available

### Secondary (MEDIUM confidence)
- Qt 6 official docs for `QFontMetrics::horizontalAdvance`, `QPainter::drawText(int,int,QString)`, `QFontDatabase::addApplicationFont`, `Q_INIT_RESOURCE` — referenced via Qt 6 doc structure knowledge; not live-fetched during this session.
- DejaVu Fonts License page (`dejavu-fonts.github.io/License.html`) — cited from general knowledge, consistent with the copy of the font installed locally via `fonts-dejavu-core` Debian package.

### Tertiary (LOW confidence)
- None. All claims in this research map to either code-truth (legacy `xvuelc.c`, current `xvue_qt_api.cpp`, `xvue/*.f`) or locally-verified Qt installation.

---

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — Qt version verified on the machine, API names verified against headers, no fast-moving components
- Architecture: HIGH for qt_add_resources pattern (Qt-documented); HIGH for palette lifetime (direct mirror of legacy)
- Call-site audit: HIGH — grep across `.f/.F/.inc` files is exhaustive within the repo
- Pitfalls: HIGH for #1, #2, #4, #7 (code-truth verified); MEDIUM for #3 (future profiling), MEDIUM for #5 (cannot runtime-verify in Phase 3); HIGH for #6, #8
- A1 (C helpers vs ABI): HIGH confidence in the finding, but user confirmation needed because it contradicts CONTEXT.md D-15/D-16 text
- A3 (`xvactivervb_` body semantics): HIGH confidence the research interpretation is correct (matches legacy and only Fortran call site); flags a real error in CONTEXT.md D-17

**Research date:** 2026-04-11
**Valid until:** 2026-05-11 (30 days — Qt 6 API surface is stable, Phase 2 code layout is frozen)
