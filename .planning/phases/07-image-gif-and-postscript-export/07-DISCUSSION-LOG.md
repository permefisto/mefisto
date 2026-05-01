# Phase 7: Image, GIF, and PostScript export - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in 07-CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-05-01
**Phase:** 07-image-gif-and-postscript-export
**Areas discussed:** ABI + trigger + filename, PS port granularity, GIF capture & fallback, PDF entry point

---

## Gray Areas Selected

| Option | Description | Selected |
|--------|-------------|----------|
| ABI + trigger + filename | New xvexport_ vs per-format vs Qt-only menu; frame feed; filename source | ✓ |
| PS port granularity | Verbatim core vs full state machine; PsEmitter helper class; globals scope; font name translation | ✓ |
| GIF capture & fallback | Probe location; native vs ffmpeg path; intermediate frame buffer | ✓ |
| PDF entry point | Raster vs vector; trigger; page geometry | ✓ |

---

## ABI + trigger + filename

### Q1: How should PNG/JPEG/GIF/PDF export be invoked from Fortran solver code (testa/wave drives multi-frame GIF today)?

| Option | Description | Selected |
|--------|-------------|----------|
| No new Fortran ABI | Image/GIF/PDF only via Qt File menu; xvpostscript_ stays the ONLY Fortran-driven export. Honors 'Fortran must not notice' invariant. testa/wave keeps emitting PS as today; Qt-side captures parallel snapshots. | ✓ |
| Add xvexport_(fmt,...) | Single new Fortran symbol with format selector. Solver code can opt in to direct image emission. Requires ABI extension. | |
| Per-format ABI symbols | Add xvpng_, xvjpeg_, xvgif_, xvpdf_ as separate Fortran-callable entry points. | |

**User's choice:** No new Fortran ABI

### Q2: How does animated GIF capture frames from testa/wave?

| Option | Description | Selected |
|--------|-------------|----------|
| Auto-snapshot on PS save | Qt intercepts every xvpostscript_(0) close as implicit 'frame done' and snapshots backing_ to in-memory QImage list. Zero solver changes. Flushed via env var XVUE_ANIM=1 or File menu toggle. | ✓ |
| New xvanimframe_ ABI | Add explicit xvanimstart_/xvanimframe_/xvanimend_ Fortran symbols. testa drivers must be patched. | |
| Defer to Phase 8 | Phase 7 ships PNG/JPEG/PS/PDF only. Animated GIF deferred. EXPORT-03 slips. | |

**User's choice:** Auto-snapshot on PS save

### Q3: Where do exported files land?

| Option | Description | Selected |
|--------|-------------|----------|
| Hybrid menu/auto | Menu → QFileDialog. Fortran-triggered → fixed legacy-compatible names (zfxy0NNN.png) in cwd. Preserves testa automation. | ✓ |
| Always QFileDialog | Every export prompts. Breaks testa/ automation. | |
| QSettings last-dir | Remember last-used-directory across runs. Auto-numbered filenames inside it. | |

**User's choice:** Hybrid menu/auto

### Q4: Where does the File → Export… menu wiring live?

| Option | Description | Selected |
|--------|-------------|----------|
| Shared shell (xvue_qt_app) | Export submenu added once in shared File menu (Phase 6.0). All five pp*_qt get entries automatically. No per-module audit. | ✓ |
| Per-module register | Each pp*_qt's register<Mod>Actions_stub_ adds module-specific Export entries. | |
| Skip menu in Phase 7 | Fortran/auto-snapshot only. File menu wiring deferred. | |

**User's choice:** Shared shell (xvue_qt_app)

---

## PS port granularity

### Q1: What scope does EXPORT-04 'preserve verbatim ~120 lines' actually cover?

| Option | Description | Selected |
|--------|-------------|----------|
| Core state machine only | Just xvpostscript_ function body (lines 1187–1304, ~120 lines). The ~150 'if(lasopsc>0) fprintf(fpo,...)' branches scattered in primitives are PS-emit helpers, ported separately. | ✓ |
| Full PS subsystem | All 95 lasopsc references + 150 fprintf sites + globals copied verbatim across all of xvuelc.c. | |
| Function body + primitive tails | xvpostscript_ verbatim, plus each Qt primitive gets verbatim fprintf-tail copied from xvuelc.c counterpart. Middle ground. | |

**User's choice:** Core state machine only

### Q2: How are the per-primitive PS-emit branches reproduced on the Qt side?

| Option | Description | Selected |
|--------|-------------|----------|
| Shared PsEmitter helper | New xvue_qt_postscript.{h,cpp} with PsEmitter::line/face/text/font/etc. wraps EXACT format strings from xvuelc.c. Single source of PS truth. | ✓ |
| Verbatim fprintf in each primitive | Copy fprintf-tail from xvuelc.c into each Qt primitive function. | |
| Defer primitive PS to Phase 8 | Phase 7 ports xvpostscript_ open/close machinery only. Primitives emit nothing. | |

**User's choice:** Shared PsEmitter helper

### Q3: Where does PS state live (lasopsc, modepsc, fpo, concat, chaine[8], fontcour, courgb, iTe/iFa/ity/…)?

| Option | Description | Selected |
|--------|-------------|----------|
| PsEmitter class instance | State as PsEmitter members; single global instance owned by XvueApp (XvueApp::psEmitter()). Mirrors Phase 5/6 bridge pattern. | ✓ |
| File-static in xvue_qt_postscript.cpp | Literal `static int lasopsc;` etc. mirroring xvuelc.c. Externs called from primitives. | |
| Per-window in XvueState | PS state attached to XvueState (per-XvueWindow). | |

**User's choice:** PsEmitter class instance

### Q4: How are X11 font names translated to PostScript font names?

| Option | Description | Selected |
|--------|-------------|----------|
| Hardcoded mapping table | X11 font-name parsing dead on Qt path. Replace with mapping QFont::family() → PS family + Italic/Oblique modifier. | ✓ |
| Port verbatim | Copy X11 font-name string-bash into PsEmitter::loadFont. Won't fire on Qt path but stays literally faithful. | |
| Always /Courier | PS output uses /Courier exclusively. Crude. | |

**User's choice:** Hardcoded mapping table

---

## GIF capture & fallback

### Q1: Where does the QImageWriter::supportedImageFormats() probe live (EXPORT-01)?

| Option | Description | Selected |
|--------|-------------|----------|
| Standalone probe exec | New xvue/qt/probes/qimagewriter_probe.cpp built by CMake. Run once at kickoff; output → PROBE.md. Matches EXPORT-01 wording verbatim. | ✓ |
| First-call inline probe | No standalone probe. PROBE.md hand-written; runtime check at first GIF write. | |

**User's choice:** Standalone probe exec

### Q2: If QImageWriter supports GIF, what's the emit path?

| Option | Description | Selected |
|--------|-------------|----------|
| QImageWriter multi-frame | Single QImageWriter opened in 'gif' format, append each captured QImage as a frame. No extra deps beyond qt6-imageformats-plugins. | ✓ |
| Always ffmpeg | Single code path even with GIF writer available. Adds ffmpeg as soft runtime dep. | |
| QMovie-based | Heavier API designed for GIF reading. Wrong fit. | |

**User's choice:** QImageWriter multi-frame

### Q3: If QImageWriter does NOT support GIF, what's the fallback?

| Option | Description | Selected |
|--------|-------------|----------|
| ffmpeg sync via QProcess | Write PNG sequence to temp dir, then QProcess::execute('ffmpeg', …) blocking on UI. testa/ animations are short. | ✓ |
| ffmpeg async + progress | QProcess started async with QProgressDialog. More wiring. | |
| Manual: PNG sequence only | Emit PNG sequence + README. No auto-conversion. EXPORT-03 risk. | |

**User's choice:** ffmpeg sync via QProcess

### Q4: Where do intermediate PNG frames live during GIF capture?

| Option | Description | Selected |
|--------|-------------|----------|
| Temp dir, auto-cleanup | QStandardPaths::TempLocation/xvue-gif-XXXXXX/ (mkstemp pattern). Cleaned on success; kept on failure. | ✓ |
| cwd zfxy0NNN.png | Match legacy convertepsgif zfxy0*.eps glob naming in cwd. | |
| In-memory + ffmpeg pipe | QImage list in RAM → QProcess::write() into ffmpeg image2pipe. Fastest, more code. | |

**User's choice:** Temp dir, auto-cleanup

---

## PDF entry point

### Q1: Raster PDF (embed canvas pixmap) or true-vector PDF?

| Option | Description | Selected |
|--------|-------------|----------|
| Raster PDF v1 | QPdfWriter + QPainter::drawPixmap(backing_). ~10 lines. Fits 'additive bonus' wording. Vector deferred to v2. | ✓ |
| True-vector PDF v1 | Capture every QPainter call as paint command, replay onto QPdfWriter. Large refactor. | |
| PS→PDF via ghostscript | Shell-out gs ps2pdf. Adds ghostscript runtime dep. | |

**User's choice:** Raster PDF v1

### Q2: How is PDF triggered?

| Option | Description | Selected |
|--------|-------------|----------|
| Qt File menu only | File → Export → PDF…. No Fortran ABI. Consistent with D-01. | ✓ |
| New xvpdf_ Fortran ABI | Extend ABI to 59 symbols. Breaks 'Fortran must not notice' for additive bonus. | |

**User's choice:** Qt File menu only

### Q3: PDF page geometry?

| Option | Description | Selected |
|--------|-------------|----------|
| Canvas-native at 72 dpi | QPdfWriter::setResolution(72) + setPageSize(QPageSize(QSizeF(xpixels, ypixels), Points)). Page matches canvas aspect exactly. | ✓ |
| A4 fit-to-page | Standard A4, auto-rotate landscape. Distorts aspect. | |
| User dialog picks | QPageSetupDialog. Overkill v1. | |

**User's choice:** Canvas-native at 72 dpi

---

## Summary

All 4 areas resolved with the recommended pick on every question (15 questions across 4 areas, 15/15 recommended). The Phase 7 design is conservative, additive, and preserves both the Fortran ABI (no new symbols) and the PostScript byte-output contract (verbatim format strings inside PsEmitter helpers).

## Claude's Discretion

Captured inline in 07-CONTEXT.md `<decisions>` section under "Claude's Discretion":
- PsEmitter threading (GUI-thread only)
- Default GIF frame delay (~10 fps)
- JPEG quality default (Qt 90)
- HiDPI export resolution behavior
- Disk-full / write-failure handling surface
- Filename-collision policy
- xvpostscript_ body becomes a single dispatch line
- Probe build wiring (cb_probe_qt + CMake target)

## Deferred Ideas

Captured in 07-CONTEXT.md `<deferred>` section:
- True-vector PDF (v2)
- Async GIF encoding (until evidence of hundreds of frames)
- QPageSetupDialog for PDF
- Multi-window PsEmitter instances
- SVG export
- GIF compression / palette tuning
