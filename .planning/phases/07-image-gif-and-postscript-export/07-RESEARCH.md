# Phase 7: Image, GIF, and PostScript export — Research

**Researched:** 2026-05-01
**Domain:** Qt 6 image / PDF / PostScript / GIF export wired into a frozen Fortran ABI; replacement of `bin/convertepsgif` ImageMagick shell-out by ffmpeg + native Qt writers.
**Confidence:** HIGH on the LOCKED decisions; HIGH on the GIF probe outcome (verified empirically on this machine); MEDIUM on the per-primitive PS-emit catalog (xvuelc.c lines verified, but a re-read during PsEmitter implementation is mandatory because the format strings vary by primitive); HIGH on the runtime-state inventory (existing VIDEO1/VIDEOFIN xwd+convert pipeline discovered — must be addressed in plan).

## Summary

Phase 7 replaces three loosely-coupled legacy export paths (the `xvpostscript_` hand-rolled emitter, the `bin/convertepsgif` ImageMagick post-processor, and the `xvue/video1.f`+`xvue/videofin.f` runtime xwd+convert shell-out) with three Qt-native ones (a verbatim `PsEmitter` C++ class, an ffmpeg-fallback animated-GIF path, and a `QFileDialog`-driven PNG/JPEG/PDF surface in the shared File menu). The Fortran ABI stays at 58 — no new entry points; `xvpostscript_(int*)` remains the sole Fortran-driven export trigger.

The critical research finding is empirically verified: **Qt 6.10.2 on Debian trixie does NOT support GIF writing.** `QImageWriter::supportedImageFormats()` returns `{bmp,cur,icns,ico,jfif,jpeg,jpg,pbm,pgm,png,ppm,tif,tiff,wbmp,webp,xbm,xpm}` — `gif` is conspicuously absent. `libqgif.so` is installed (it ships with `libqt6gui6`), but it is a read-only QImageReader plugin. Direct probe: `QImageWriter::write` for `setFormat("gif")` returns `false` with `errorString()=="Unsupported image format"`. **The PROBE.md (EXPORT-01) outcome is pre-determined: ffmpeg fallback (D-11) is the only viable path on the target system.** This does not invalidate the probe — it must still run at kickoff for documentation and to detect any future Debian QtImageFormats add-on.

A second critical finding: there is already a Fortran-side animation-capture system (`xvue/video1.f` / `xvue/videofin.f` driving `LVIDEO` flag in `incl/trvari.inc`) that uses `xwd` + ImageMagick `convert` shell-outs. It is **independent** of `xvpostscript_` and is called from 18+ tracer subroutines (`flui/trvi2d.f`, `ther/trplse.f`, `ther/trisot.f`, …). The `xvue_qt_postscript`-snapshot interception in CONTEXT D-02 does NOT capture these — it captures the EPS-emit save points. Both paths must be addressed: D-02's auto-snapshot covers the EPS-driven workflow; the existing VIDEO1/VIDEOFIN Fortran routines must be redirected to a Qt-native equivalent (or kept routing through `convert` until Phase 9 RETIRE-03). This is a runtime-state-inventory item.

**Primary recommendation:** Execute Phase 7 in 6 plans. Plan 01 (probe + PROBE.md) gates everything. Plans 02 (PsEmitter skeleton + xvpostscript_ verbatim port) and 04 (PNG/JPEG/PDF File menu actions) can run in parallel. Plan 03 (per-primitive PsEmitter helpers + emit-site wiring across xvue_qt_api.cpp) follows Plan 02. Plan 05 (animated GIF + ffmpeg dispatch + auto-snapshot + VIDEO1/VIDEOFIN reroute) follows Plan 04. Plan 06 (test/validation suite + golden EPS samples + EXPORT-06 grep gate) closes. Total ABI delta: 0. Total LOC delta: ~700 added (PsEmitter + dispatchers + tests), ~120 deleted (Qt stub bodies replaced).

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**Export ABI & trigger surface**
- **D-01:** No new Fortran ABI for image/GIF/PDF. PNG, JPEG, GIF, and PDF are triggered exclusively via the Qt File menu (or, for animation, via env-var/menu toggle that Qt-side intercepts). `xvpostscript_(int *lasops)` stays the ONLY Fortran-callable export entry point. ABI symbol count unchanged at 58.
- **D-02:** Auto-snapshot animation capture. When `XVUE_ANIM=1` env var is set OR the File → "Capture animation" menu toggle is active, `PsEmitter` snapshots `backing_` to an in-memory `QImage` list on every `xvpostscript_(*lasops==0)` close. Frame list flushed to `animation.gif` on next `xvpostscript_(0)` close after toggle disabled, or via explicit File → "End animation" menu.
- **D-03:** Hybrid filename. Menu-triggered exports prompt via `QFileDialog::getSaveFileName` with `QSettings("xvue/export/last_dir")` remembered across runs. Auto-snapshot frames during animation capture use fixed legacy-compatible names (`zfxy0001.png`, `zfxy0002.png`, …) in cwd. Final GIF name: `animation.gif` in cwd (overwrites without warning).
- **D-04:** File menu wired in shared shell. `File → Export → {PNG…, JPEG…, GIF…, PDF…}` entries added once in the shared `xvue_qt_app` File-menu builder (Phase 6.0 module). All five `pp*_qt` executables inherit automatically. No per-module audit/wiring needed.

**PostScript port granularity**
- **D-05:** EXPORT-04 'verbatim ~120 lines' = `xvpostscript_` function body only (xvuelc.c:1187–1304). Ports byte-for-byte into `PsEmitter::handleLasops(int lasops)`. The ~150 in-primitive `if(lasopsc>0) fprintf(fpo,...)` branches are PS-emit helpers ported as separate methods on `PsEmitter` with the SAME format strings.
- **D-06:** `PsEmitter` helper class. New `xvue/qt/src/xvue_qt_postscript.{h,cpp}`. Public methods mirror xvuelc.c emit sites: `line()`, `face()`, `text()`, `loadFont()`, `ellipse()`, `flpt()`, `fond()`, `clear()`. Each method's body uses the EXACT format strings from the corresponding xvuelc.c `fprintf(fpo,...)` site. README_COORDS.md `ypixels - y` flip applied inside each method, never in the caller. Qt primitives in `xvue/qt/src/xvue_qt_api.cpp` call the matching `PsEmitter::*` after their `QPainter::*` calls (when `psEmitter().active()`).
- **D-07:** PS state encapsulated in `PsEmitter` instance. All ~15 file-statics from xvuelc.c (`lasopsc`, `modepsc`, `fpo`, `concat[1024]`+`nbrcon`, `chaine[8]`+`longchaine[8]`, `fontcour`, `courgb[3]`, `iTe`, `iFa`, `ity`, `iep`, `iPo`, `ire`, `iRe`, `iel`, `iEl`, `iFP`, `menu`) become non-static members. Single global instance owned by `XvueApp`: `XvueApp::psEmitter()` accessor.
- **D-08:** PS font-name translation = hardcoded mapping table. The X11 font-name string-bash in `xvchargefonte_` is dead code on the Qt path. `PsEmitter::loadFont(QFont qf, …)` uses a 4-entry lookup table:
  ```
  "Courier"                → "/Courier"
  "Times"                  → "/Times-Roman"
  "Helvetica"              → "/Helvetica"
  "New Century Schoolbook" → "/NewCenturySchlbk"
  ```
  Plus Italic/Oblique/Bold modifier suffix logic ported verbatim from xvuelc.c lines ~1530–1565. Other Qt fonts fall back to `/Helvetica` with warn-once log entry.

**GIF capture & fallback**
- **D-09:** Standalone `QImageWriter` probe executable. New `xvue/qt/probes/qimagewriter_probe.cpp` built via new `bin/cb_probe_qt` script and CMake target. Run once at Phase 7 kickoff; stdout redirected to PROBE.md. Defensive runtime re-check at first GIF write inside `PsEmitter::beginAnimation()`.
- **D-10:** Native GIF path (probe-positive). When `QImageWriter::supportedImageFormats()` lists `"gif"`, single `QImageWriter` opened with `setFormat("gif")`, each captured `QImage` written as a frame via the multi-frame extension (`setText("Delay","10")`).
- **D-11:** ffmpeg fallback (probe-negative). Write PNG sequence into the temp dir from D-12, then `QProcess::execute("ffmpeg", {"-y", "-framerate","10", "-i","frame_%04d.png", "out.gif"})` blocking on the GUI thread. `convert` (ImageMagick) is forbidden anywhere in `xvue/qt/`. EXPORT-06 verified by `grep -rn 'convert' xvue/qt/` returning empty.
- **D-12:** Frame buffer = mkstemp-pattern temp dir. `QStandardPaths::writableLocation(QStandardPaths::TempLocation) + "/xvue-gif-XXXXXX/"` via `QTemporaryDir`. Frames deleted on successful GIF write; kept on failure with path logged to console dock.

**PDF entry point**
- **D-13:** Raster PDF v1. `QPdfWriter` + `QPainter::drawPixmap(0, 0, backing_)`. ~10 lines total. True-vector PDF deferred to v2.
- **D-14:** PDF triggered via Qt File menu only. No Fortran ABI extension.
- **D-15:** Canvas-native page geometry at 72 dpi. `QPdfWriter::setResolution(72)` + `setPageSize(QPageSize(QSizeF(xpixels, ypixels), QPageSize::Point))`. PDF page aspect matches canvas exactly — no fit-to-A4 distortion.

**ImageMagick scope**
- **D-16:** Phase 7 removes ImageMagick from the Qt code path only. `bin/convertepsgif` and `convert` shell-outs from `td/`, `testa/`, `testf/` stay alive for the legacy X11 backend during the one-release-cycle A/B window. Full audit + removal is RETIRE-03 in Phase 9. Verification: `grep -rn 'convert' xvue/qt/` returns empty.

### Claude's Discretion

- PsEmitter threading: GUI thread only — `XVUE_QT_ASSERT_MAIN_THREAD()` at the top of every public method. `QImageWriter` and `QProcess::execute` called synchronously from main thread. Async / progress-dialog path deferred.
- Default GIF frame delay: 10/100 sec (≈ 10 fps). Adjustable via `QSettings("xvue/export/gif_delay")`.
- JPEG quality default: Qt default (90). Configurable via `QSettings("xvue/export/jpeg_quality")`.
- HiDPI export resolution: PNG/JPEG snapshot exports at backing pixel dims (i.e., `2*` logical for HiDPI). PDF stays at 72 dpi canvas-native logical pixels.
- Disk-full / write-failure handling: status-bar message + console-dock log. `QMessageBox::critical` only when menu-triggered.
- Filename-collision policy: menu-triggered uses `QFileDialog`'s overwrite prompt; auto-snapshot fixed names overwrite without warning.
- `xvpostscript_` body in `xvue_qt_api.cpp`: stays as the ABI-frozen entry; immediately delegates to `XvueApp::psEmitter().handleLasops(*lasops)`.
- Probe build wiring: new `bin/cb_probe_qt` script + `xvue/qt/probes/CMakeLists.txt` target. PROBE.md committed at Phase 7 kickoff before the first plan executes.

### Deferred Ideas (OUT OF SCOPE)

- **True-vector PDF.** Phase 7 ships raster-only PDF (D-13). True-vector PDF requires capturing every `QPainter` call into a display list and replaying onto `QPdfWriter`. Defer to v2.
- **Async GIF encoding for very large animations.** Defer until evidence testa/wave / testa/cavity2d produce hundreds of frames.
- **`QPageSetupDialog` for PDF page sizing.** Defer until requested.
- **Multi-window `PsEmitter` instances.** Single XvueApp-owned instance suffices for v1.
- **SVG export.** Not in REQUIREMENTS.md.
- **EPS variant separate from .ps.** Same byte stream under `TEMPORAIRE.EPS`; users can rename freely.
- **GIF compression / palette tuning.** Default Qt GIF plugin output (or ffmpeg defaults) acceptable for testa/ visual A/B.
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| EXPORT-01 | Probe `QImageWriter::supportedImageFormats()` at kickoff; result recorded in `PROBE.md`; chosen GIF strategy documented (`QImageWriter` loop OR PNG-sequence + ffmpeg — never ImageMagick). | §"Empirical probe outcome" — verified Qt 6.10.2 on Debian trixie does NOT write GIF; ffmpeg fallback path (D-11) is the only viable route. Probe still required by spec; PROBE.md content is `gif: ABSENT` confirmed. |
| EXPORT-02 | `XvueExport` writes PNG and JPEG single-frame snapshots of the current canvas via `QImageWriter` / `QImage::save`. | §"Standard Stack" `QImageWriter`; §"Code Examples" `XvueExport::savePng` / `saveJpeg`. `backing_.toImage().save(...)` is the canonical 3-line pattern. |
| EXPORT-03 | `XvueExport` writes animated GIFs from a sequence of frames matching the legacy `bin/convertepsgif` visual output on `testa/wave` and `testa/cavity2d`. | §"GIF strategy" — ffmpeg via `QProcess::execute` with `-framerate 10`. Visual parity verified by frame-hash + frame-count comparison vs legacy convertepsgif (golden file in `testa/wave/expected/`). §"Runtime State Inventory" item 2 — there are TWO legacy GIF pipelines: `bin/convertepsgif` (single line, cyl53-specific) AND `xvue/video1.f`+`xvue/videofin.f` (runtime xwd+convert called from 18+ tracer subroutines). Both must be addressed. |
| EXPORT-04 | `xvpostscript_` body moved verbatim from `xvuelc.c` into `xvue/qt/src/xvue_qt_postscript.cpp` — no switch to `QPrinter` PostScript output. | §"PostScript per-primitive emit catalog" — exact line ranges (1187–1304 verified ~117 lines body; 14 per-primitive `if(lasopsc>0)` gates; 151 total `fprintf(fpo, ...)` sites). PsEmitter port preserves the format strings byte-for-byte; `ypixels-y` flip preserved at emit time. |
| EXPORT-05 | PDF export added as a bonus via `QPdfWriter` (raster path), NOT by modifying `xvpostscript_`. | §"Standard Stack" `QPdfWriter`; §"Code Examples" 72-dpi canvas-native page geometry pattern. New File menu action; ~10 LOC. |
| EXPORT-06 | After EXPORT-01..05 ship, ImageMagick is no longer invoked by any Qt-backend code path. Verify: `grep -rn 'convert' xvue/qt/` returns empty (or only the PDF format keyword). | §"Common Pitfalls" Pitfall 5 — phase-end grep gate; Plan 06 owns it. NOTE: `bin/convertepsgif`, `xvue/video1.f`, `xvue/videofin.f` are NOT in `xvue/qt/` — they remain alive until Phase 9 RETIRE-03. The grep gate is scoped to `xvue/qt/` ONLY (locked by D-16). |
</phase_requirements>

## Project Constraints (from CLAUDE.md)

These directives are extracted from `CLAUDE.md` and have the same authority as locked decisions. The planner must verify Phase 7 plans comply:

- **`bin/cbl_tout` must succeed before any commit** — full build canonical gate. Phase 7 adds new sources to `xvue/qt/CMakeLists.txt` (`xvue_qt_postscript.cpp`) and a new probes subdirectory; both must compile under the existing flags (Qt6::Widgets, Qt6::Gui, Qt6::PrintSupport, `-Wall -Wextra -fno-openmp`).
- **The X11 build (`bin/cbl_tout`) must keep producing identical `pp/pp*` executables** — `xvue/xvuelc.c` is NOT modified. PsEmitter is Qt-only; the legacy emitter stays in xvuelc.c untouched.
- **`testa/` smallest test cases must continue to pass after every change** — `testa/wave` and `testa/cavity2d` are the canonical Phase 7 fixtures (REQUIREMENTS.md EXPORT-03).
- **Ask before installing system packages** — Phase 7 wants to use `ffmpeg` (already installed at `/usr/bin/ffmpeg` 8.1-3+b1, verified) and the existing Qt6 GIF plugin (already installed via `libqt6gui6`). No new apt installs needed. Plan 01 task should `command -v ffmpeg || ASK-USER` defensively.
- **Fortran 77 fixed-form discipline** — Phase 7 does NOT add Fortran code (D-01: no new ABI). The existing `xvue/video1.f` / `xvue/videofin.f` routines may be touched if Plan 05 reroutes them to call `xvue_module_init`-registered Qt sinks — if so, fixed-form columns 1-5 / 6-continuation / 7-72 must be honored. Strong recommendation: leave `video1.f` / `videofin.f` alone in Phase 7 (they are X11-side legacy, slated for Phase 9). Document the gap explicitly.
- **Never force-push; never bypass hooks.** N/A for Phase 7 mechanics.
- **Commit after every logical step where rollback would be useful.** Phase 7's plan boundaries (PROBE.md commit, PsEmitter skeleton landed, per-primitive emit site wired, GIF dispatcher landed, golden tests passing) are natural rollback points.
- **Programming norms in `doc/normes.ps`** — applies to Fortran only; the Qt/C++ side follows the Qt 6 idiomatic style established in Phases 0-6.5 (`XVUE_QT_ASSERT_MAIN_THREAD`, `warn_once`, snake_case proc names, member trailing underscore).

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Qt6::Gui | 6.10.2 [VERIFIED: `pkg-config --modversion Qt6Gui` on host] | `QImageWriter`, `QImage`, `QPdfWriter`, `QPainter`, `QPixmap` | Already linked into `libxvueqt.a`. PNG/JPEG/PDF writers are part of Qt 6 base — no extra apt package needed. `QPdfWriter` ships in `Qt6::Gui` (NOT `Qt6::PrintSupport`) [CITED: doc.qt.io/qt-6/qpdfwriter.html]. |
| Qt6::Core | 6.10.2 [VERIFIED] | `QProcess`, `QStandardPaths`, `QTemporaryDir`, `QSettings`, `QFile` | Already linked. `QProcess::execute` is the canonical sync `wait+exit-code` invoker for ffmpeg. `QTemporaryDir` (mkdtemp wrapper) handles the D-12 `xvue-gif-XXXXXX/` scratch directory. |
| Qt6::Widgets | 6.10.2 [VERIFIED] | `QFileDialog::getSaveFileName`, `QMessageBox::critical`, `QStatusBar` | Already linked. File menu actions and progress feedback. |
| ffmpeg | 8.1-3+b1 [VERIFIED: `ffmpeg -version` on host] | PNG-sequence → animated GIF encoder | The D-11 fallback. Available at `/usr/bin/ffmpeg`. Standard Debian package. Replaces ImageMagick `convert` for the Qt path — produces identical (or better-compressed) GIFs from a `frame_%04d.png` glob. |

**Verified versions** (registry probe 2026-05-01):
- Qt 6.10.2: `pkg-config --modversion Qt6Core` returned `6.10.2`. Debian package `qt6-base-dev:amd64` version `6.10.2+dfsg-7`.
- ffmpeg 8.1-3+b1: shipped as `/usr/bin/ffmpeg` on Debian trixie.

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Qt6::Test | 6.10.2 | QTest slots for PsEmitter golden-output and PNG/PDF writer regression | Plan 06 — new `xvue_qt_export_tests` and `xvue_qt_postscript_tests` binaries follow the Phase 6.x pattern (one `add_executable` per concern in `xvue/qt/tests/CMakeLists.txt`). |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `QPdfWriter` raster (D-13) | `QPdfWriter` vector via display-list capture | True-vector PDF means hooking every QPainter call, replaying on the PDF QPainter. Out of scope for v1 (Deferred). |
| `QPdfWriter` (Qt6::Gui) | `QPrinter` + `QPrinter::PdfFormat` (Qt6::PrintSupport) | REQUIREMENTS.md EXPORT-05 mentions `QPrinter::PdfFormat` but `QPdfWriter` is the modern Qt-6 idiom; same output, simpler API, no need for `Qt6::PrintSupport` (already linked anyway). Either is fine; pick `QPdfWriter` for clarity. |
| `ffmpeg` for GIF | `gifski` / `magick` / `lottie` / `QtGifImage` (3rd-party) | ffmpeg is already in the dependency chain, ImageMagick is being dropped, and gifski/QtGifImage need new apt packages or vendoring. ffmpeg wins on zero-deps and known-good dithering. |
| `QImageWriter` multi-frame GIF (D-10) | ffmpeg fallback (D-11) | **Empirically forced:** Qt 6.10.2 on Debian trixie does NOT export GIF. `QImageWriter::supportedImageFormats()` does not include `gif`. The native path remains in code as a defensive fast-path for future Debian QtImageFormats add-on, but D-11 is the realized path. |
| Hand-rolled mkdtemp | `QTemporaryDir` | `QTemporaryDir` provides RAII cleanup + unique name generation with a single line. No reason to hand-roll. |
| `system("ffmpeg ...")` | `QProcess::execute` | `QProcess` properly escapes argv, captures exit code, lets us inspect stderr on failure; `system()` shells out and is shell-injection-prone. CONVENTIONS already established (`QProcess` is the chosen primitive). |

**Installation** (no apt installs needed — all already installed):
```bash
# Already present (verified):
dpkg -l qt6-base-dev      # 6.10.2+dfsg-7
which ffmpeg              # /usr/bin/ffmpeg (8.1-3+b1)
ls /usr/lib/x86_64-linux-gnu/qt6/plugins/imageformats/libqgif.so  # exists, libqt6gui6
```

**Version verification command pattern for the planner:**
```bash
pkg-config --modversion Qt6Core    # confirm Qt 6.10+
ffmpeg -version | head -1          # confirm ffmpeg 4+ (any modern version OK)
```

## Architecture Patterns

### Recommended Project Structure

```
xvue/qt/
├── src/
│   ├── xvue_qt_api.cpp                  ← Plan 02: xvpostscript_ delegates to PsEmitter (1 line)
│   ├── xvue_qt_api.cpp (other 13-15 emit-site primitives)  ← Plan 03: emit-site calls into PsEmitter
│   ├── xvue_qt_app.{h,cpp}              ← Plan 02: psEmitter() accessor (mirrors menuBridge, eventBridge)
│   ├── xvue_qt_postscript.{h,cpp}       ← NEW (Plan 02): PsEmitter class
│   ├── xvue_qt_export.{h,cpp}           ← NEW (Plan 04): XvueExport (PNG/JPEG/PDF/GIF dispatcher)
│   ├── xvue_qt_window.cpp               ← Plan 04: File → Export submenu builder
│   └── ... (existing untouched)
├── probes/                              ← NEW (Plan 01)
│   ├── CMakeLists.txt
│   └── qimagewriter_probe.cpp           ← 10-20 LOC; prints supportedImageFormats()
├── tests/
│   ├── test_xvue_qt_postscript.cpp      ← NEW (Plan 06): golden-EPS byte-level compare
│   ├── test_xvue_qt_export.cpp          ← NEW (Plan 06): PNG/JPEG/PDF/GIF round-trip
│   └── ... (existing untouched)
└── CMakeLists.txt                       ← Plans 01-05 each add one or two source lines
```

```
bin/
├── cb_probe_qt                          ← NEW (Plan 01): wraps CMake probes target build + run
└── ... (existing cbl_tout_qt etc. untouched)
```

```
.planning/phases/07-image-gif-and-postscript-export/
├── 07-CONTEXT.md                         ← already gathered
├── 07-RESEARCH.md                        ← THIS FILE
└── PROBE.md                              ← NEW (Plan 01): output of qimagewriter_probe; commit gates everything
```

### Pattern 1: Bridge-as-XvueApp-singleton (carried from Phase 5 D-02 / Phase 6.0)

**What:** Stateful Qt-side helper class is a `QObject` member on `XvueApp`, accessed via accessor.
**When to use:** Any helper that owns mutable state and lives for the QApplication lifetime.
**Example pattern (PsEmitter, mirrors `XvueEventBridge` / `XvueMenuBridge`):**
```cpp
// Source: existing pattern in xvue/qt/src/xvue_qt_app.{h,cpp}
class XvueApp : public QObject {
public:
    static XvueApp& ensure();
    XvueEventBridge& eventBridge();   // Phase 5 [VERIFIED: README.md]
    XvueMenuBridge&  menuBridge();    // Phase 6.0 [VERIFIED: README.md]
    PsEmitter&       psEmitter();     // Phase 7 NEW
private:
    std::unique_ptr<PsEmitter> ps_emitter_;
};
```

### Pattern 2: ABI entry as one-line dispatcher

**What:** Frozen `extern "C"` Fortran entry points stay tiny; their body is a single `XvueApp::ensure(); XVUE_QT_ASSERT_MAIN_THREAD(); SOMETHING.method(args);` invocation.
**When to use:** Any `proc(name)` whose semantics live in a stateful helper.
**Example:**
```cpp
// Source: per D-05/D-06/D-07 + existing Phase 6.0 dispatch pattern in xvue_qt_api.cpp
void proc(xvpostscript)(int *lasops) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    XvueApp::ensure().psEmitter().handleLasops(*lasops);
}
```

### Pattern 3: Per-primitive PS-emit branch

**What:** After every `QPainter::*` call in `xvue_qt_api.cpp` that has a counterpart `if(lasopsc>0) fprintf(fpo,...)` branch in xvuelc.c, append the equivalent `psEmitter().<verb>(...)` call.
**When to use:** Plan 03 — for each primitive enumerated in §"PostScript per-primitive emit catalog".
**Example (xvtrait_):**
```cpp
// Source: PsEmitter::line corresponds to xvuelc.c:1942-1985 verbatim
void proc(xvtrait)(int *x1, int *y1, int *x2, int *y2) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto* st = XvueApp::ensure().window_slot()->canvas()->state();
    st->painter()->drawLine(*x1, *y1, *x2, *y2);
    XvueApp::ensure().psEmitter().line(*x1, *y1, *x2, *y2);  // NEW
}
```
The `PsEmitter::line` method internally checks `if (!active()) return;` and applies `ypixels - y` flip.

### Pattern 4: ffmpeg dispatch via QProcess::execute

**What:** Synchronous shell-out with no shell interpolation, exit code checked, stderr captured.
**When to use:** Plan 05 GIF assembly.
**Example:**
```cpp
// Source: idiomatic Qt6 — QProcess::execute returns int (0 = success)
QStringList args{"-y", "-framerate", "10",
                 "-i", tempDir.filePath("frame_%04d.png"),
                 outputPath};
int exit_code = QProcess::execute("ffmpeg", args);
if (exit_code != 0) {
    XvueApp::ensure().consoleDock().log("ffmpeg failed: exit " + QString::number(exit_code));
    QMessageBox::critical(...);
}
```

### Anti-Patterns to Avoid

- **Storing flipped Y in PsEmitter state.** The `ypixels - y` flip MUST happen at emit time inside each `PsEmitter::line/face/text` body. Never pre-flip and store. (xvue/README_COORDS.md is the write-once normative source.)
- **Calling `system("convert ...")` from xvue_qt_api.cpp.** EXPORT-06 verification is `grep -rn 'convert' xvue/qt/`; any C-string `"convert"` literal in this tree fails the gate.
- **Calling `exit(1)` on `fopen` failure.** xvuelc.c does this in two places (lines 1234, 1282). Qt port replaces with `QMessageBox::critical` + `XvueApp::ensure().quitGracefully()` so unsaved project state isn't lost. (Documented in CONTEXT specifics.)
- **Releasing `chaine[i]` with `*chaine[i] = '\0'` then `free(chaine[i])`.** xvuelc.c:1252 has this pattern (clobbers freed memory). Port faithfully but document — subtle bug; preserving for byte-parity is acceptable since the freed memory isn't read after.
- **Calling `XApplication::exec()` or `QApplication::exec()` from any export path.** Already enforced by the SHELL-03 verify-no-exec gate (Phase 1 Plan 01). Phase 7 just inherits.
- **Letting auto-snapshot frames pile up unbounded.** `XVUE_ANIM=1` mode opens a frame buffer; without an upper-bound check, a long-running solver could fill `/tmp` and crash. Document in Plan 05 — soft cap at, e.g., 10000 frames with a console-dock warning at 1000.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Animated-GIF encoding | A C++ LZW encoder | ffmpeg via `QProcess::execute` | Decades of optimization, color-quantization, and dithering are baked into ffmpeg. Hand-rolled LZW is a multi-week project and would be slower + worse-looking. |
| Temporary directory creation | `mkdtemp` + manual `rm -rf` | `QTemporaryDir` | RAII cleanup. Cross-platform handling. Already in Qt6::Core. |
| PNG / JPEG encoding | libpng / libjpeg-turbo direct calls | `QImage::save("path", "PNG")` | One line vs hundreds. Qt's `QImageWriter` is what `QImage::save` calls under the hood — well-tested. |
| PDF generation | hand-rolling PDF byte stream | `QPdfWriter` + `QPainter` | PDF is a structured format with object streams, xref tables, byte offsets. Hand-rolling means recreating Qt's 5000+ LOC PDF engine. Use `QPdfWriter`. |
| Multi-frame GIF probe at runtime | a hand-rolled magic-byte sniff | `QImageWriter::supportedImageFormats()` | Qt's plugin loader is the source of truth on this system. The list is exactly what `QImageWriter::write` will accept. |
| Hardcoded ffmpeg path | `/usr/bin/ffmpeg` literal | `"ffmpeg"` (relies on PATH) | More portable across distros (and across the maintainer's local installs). Defensive `command -v ffmpeg` check at first GIF write surfaces a clean error to console dock. |
| Hand-rolling argv string concatenation for ffmpeg | `system("ffmpeg -y -i ...")` | `QProcess::execute("ffmpeg", QStringList{...})` | No shell interpolation, no escaping bugs, exit code captured. |

**Key insight:** Phase 7's value is in the *exact-byte preservation* of the PostScript output (D-05, D-06: the format strings ARE the contract because end-users have papers and pipelines built on the .eps byte format). Everything else (PNG, JPEG, PDF, GIF) should be the most idiomatic Qt one-liner possible — the implementation surface area beyond `PsEmitter` is intentionally minimal.

## Runtime State Inventory

> Phase 7 is a refactor / migration phase that crosses backend lines. The grep audit on `xvue_qt_api.cpp` finds Qt-side code, but several runtime systems on the X11 + Fortran side embed `convert`, `xwd`, `TEMPORAIRE.EPS`, and `LVIDEO` references that the planner needs to know about explicitly.

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| **Stored data** | `TEMPORAIRE.EPS` written into project `cwd` by `xvpostscript_(*lasops>0 && lasopsc==0)` (xvuelc.c:1228, 1276). Same filename used by Qt path (D-05 verbatim port). `TEMPORAIRE.QUA` quality histogram side-file at xvuelc.c:1222, 1269. | Code edit only: PsEmitter writes the same filenames. No data migration. **Document:** if a user has both X11 and Qt MEFISTO running in same dir, they will collide on `TEMPORAIRE.EPS`. Out of scope for v1 (single backend per session is the model). |
| **Live service config** | None — MEFISTO has no external service config (n8n / Datadog / Tailscale / etc.). | None. |
| **OS-registered state** | None — no Windows Task Scheduler, no systemd unit, no pm2, no launchd plists. MEFISTO is run interactively via `INITIER` / `MAILLER` / `ELASTICER` shell launchers. | None. |
| **Secrets/env vars** | New env var introduced by Phase 7: `XVUE_ANIM=1` (D-02 auto-snapshot enable). Existing related env vars from Phase 5: `MEFISTO_XVSOURIS_AUTOEXIT`, `MEFISTO_XVSOURIS_DEBUG`, `MEFISTO_XVFERMER_READY_FILE`, `MEFISTO_XVFERMER_HOLD_MS`, `MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS`. No SOPS / .env / vault rename. | Document XVUE_ANIM in `xvue/qt/README.md` Phase 7 section as an additive feature. No code-rename hazard since this is a new key. |
| **Build artifacts / installed packages** | `pp/pp*_qt` rebuild required (Phase 7 adds `xvue_qt_postscript.cpp` + `xvue_qt_export.cpp` to `libxvueqt.a`). `bin/cbl_tout_qt` is the canonical full rebuild — already invokes a clean `xvue/build/` step. The new `bin/cb_probe_qt` script must be marked executable (`chmod 755`) to match repo convention (every script in `bin/` has `+x`). | Plan 01 includes `chmod 755 bin/cb_probe_qt`. The probes/ CMake target is built only on demand (NOT by `cbl_tout_qt`); `bin/cb_probe_qt` invokes it standalone. |
| **Critical newly-discovered: existing GIF capture in solvers** | `xvue/video1.f`, `xvue/videofin.f`, `xvue/videonm.f` shell out to `xwd` + `convert` from 18 Fortran tracer subroutines (`flui/trvi2d.f`, `flui/trvi3d.f`, `flui/parpartr.f`, `flui/tttsupa2d.f`, `ther/trplse.f`, `ther/trisot.f` x2, `ther/trlldr.f`, `ther/trnlserr.f`, `ther/trnlsetst.f`, `ther/trnlsemxu.f`, `ther/trsosf.f`, `ther/trso1so.f`, `ther/trzont.f`, `flui/trco2d.f`, … — full set found by `grep -rln 'VIDEO1' flui/ ther/`). Triggered by `LVIDEO` flag in COMMON `/TRVAPS/` (`incl/trvari.inc`). Produces `<NMPROJ>_<name>NNNN.jpg` per frame, then `convert *.jpg -delay 100 ...gif` as final assembly. **This is a SECOND legacy GIF pipeline distinct from `bin/convertepsgif`.** | **Phase 7 scope decision (must be in PLAN 05):** *Do NOT touch* `video1.f` / `videofin.f` in Phase 7. They are X11-side legacy that stays alive during the A/B window and is retired in Phase 9 (RETIRE-03). Phase 7's auto-snapshot path (D-02) is independent — it captures EPS save points via `xvpostscript_`, not the LVIDEO/`xwd` capture points. Document the gap in PLAN 05 + `xvue/qt/README.md`: Qt-backend solvers will produce no GIF from `LVIDEO=1` paths until Phase 9 wires them through Qt; meanwhile, the X11 backend continues to handle these. **OR** if the planner judges this too big a hole, Plan 05 can additionally provide a runtime `XVUE_QT_VIDEO_HOOK` env var that intercepts `xwd` → `QPixmap::save(jpg)` and `convert` → `ffmpeg` — but this expands scope. Default recommendation: defer. |
| **`bin/convertepsgif` callers in td/testa/testf** | Empirically searched `grep -rn 'convert\|convertepsgif' td/ testa/ testf/` returned **zero matches**. The single `convert -rotate 90 zfxy0*.eps -extent 980x550 cyl53zfxy.gif` line in `bin/convertepsgif` and `bin.lnx64/convertepsgif` is the only caller. Researcher could not find any `cyl53` directory under testa/ or testf/ either. | The planner can safely state: bin/convertepsgif is "orphaned" / cyl53-specific, NOT used by any tracked test. Phase 7 needs to faithfully replicate the rotate+extent post-processing ONLY if/when a user calls the script manually. The Qt-side GIF path emits canvas-native dims; if cyl53-style 90°-rotation is wanted, document as a per-case CLI option (`xvue_export gif --rotate 90`) — out of scope for v1 unless requested. |
| **`testa/wave` and `testa/cavity2d` GIF expectations** | `testa/wave/` has only `wave.mesh` (1.9K) and `wave.wave` (2.6K). `testa/cavity2d/` has only `.meshbf`/`.meshth`/`.stoke56cr` files. **Neither directory contains an existing reference .gif file.** No hand-rolled "expected" GIF byte-stream in-tree. EXPORT-03 visual A/B vs `bin/convertepsgif` must therefore be done by *running* the legacy pipeline on these tests (X11 backend) and comparing to the Qt pipeline output by eye + frame-count + frame-hash — there is no committed golden file. | Plan 06 must commit a baseline `testa/wave/expected_legacy.gif` and `testa/cavity2d/expected_legacy.gif` produced from the X11 backend BEFORE Phase 7's emitter lands, as the goldens for visual A/B. This is a one-time bootstrap step in Plan 06 (or Plan 01 — earlier is safer). |

**Nothing found in:** Live service config, OS-registered state. Verified by `find /etc/systemd /etc/cron.d -name 'mefisto*'` (none) and absence of any Windows Task Scheduler / pm2 references in the repo.

## Common Pitfalls

### Pitfall 1: PostScript byte-format drift via `printf` width specifiers

**What goes wrong:** A "harmless" tidy-up of `%6i %6i %4.2f %4.2f %4.2f` to `%i %i %f %f %f` changes the byte stream. Users have downstream PostScript-postprocessing scripts (TeX inclusion, lpr filters, `gs` renders) keyed to the column widths.
**Why it happens:** Modern Qt/C++ developers reflexively prefer `QString::number` or `<<`, which produce different formatting than `printf`.
**How to avoid:** PsEmitter MUST use `std::fprintf(fpo_, ...)` literally (xvuelc.c port), NOT `QTextStream` or `<<`. Format strings ARE the contract. Plan 03 task description must explicitly forbid format-string changes. Plan 06 byte-level golden compare catches any drift.
**Warning signs:** Any code review comment suggesting "we could simplify this with QString::number" — push back firmly.

### Pitfall 2: `ypixels - y` Y-flip applied twice or in the wrong layer

**What goes wrong:** PsEmitter::line caller already flips Y (e.g., the `xvtrait_` Qt body pre-computes `int ypiof = state->ypixels() - *y1;` and passes that), so PsEmitter flips a second time → wrong PostScript Y. Or Qt-side QPainter accidentally consumes a pre-flipped Y → screen mesh upside-down.
**Why it happens:** Coordinate transforms across module boundaries are notoriously easy to misalign.
**How to avoid:** README_COORDS.md is the write-once normative source: **flip ONLY in `PsEmitter::line/face/text` bodies, NEVER in callers**. The mantra: "store and pass canvas Y-down (matching Qt and Fortran); flip at PostScript emit time." Plan 02 PsEmitter unit test asserts: given `(x=10, y=20)` with `ypixels=600`, the emitted line contains `580` (= 600-20).
**Warning signs:** A `psEmitter().line(x1, ypixels - y1, ...)` call site → bug.

### Pitfall 3: `chaine[lasopsc-4]` buffer underflow when lasopsc < 4

**What goes wrong:** xvuelc.c:1568, 1754, 1812, 1861, 1898, 2088, 2614, 2740, … all dispatch on `lasopsc < 3 ? fprintf(fpo,...) : sprintf(&chaine[lasopsc-4][...], ...)`. The `<3` branch writes to file; the `>=3` branch writes to a stored menu buffer. **If lasopsc takes any value `0..2` after the gate `lasopsc > 0` is true (i.e., `lasopsc==1` or `lasopsc==2`), the negative-index `chaine[1-4]` would corrupt memory** — but the dispatch correctly routes those to the file branch. lasopsc 3..10 routes to chaine[lasopsc-4] = chaine[-1..6]. There IS a chaine[-1] buffer underflow when `lasopsc==3`!
**Why it happens:** Subtle off-by-one in original C code. Looking at the code again: `chaine[lasopsc-4]` for `lasopsc==4` writes to `chaine[0]`. So lasopsc>=4 is safe. But `lasopsc==3` triggers `chaine[-1]`. **Audit:** xvuelc.c:1298 caps `lasopsc = min(*lasops, 11)` and the `>= 3` arm is gated on `lasopsc < 3` being false, so `lasopsc==3` does take the chaine path. This may be a latent bug in xvuelc.c — but for verbatim port (D-05), preserve faithfully and document.
**How to avoid:** Plan 02 reads xvuelc.c:1187-1304 with surgical care, ports the dispatch table as-is, and adds a `[[maybe_unused]]` comment noting the chaine[-1] question. Phase 7 is NOT a bug-fix phase for legacy emitter quirks — fidelity wins (project value: scientific reproducibility).
**Warning signs:** Address-sanitizer red on Plan 02 unit tests at `lasopsc==3`. If it's already broken in xvuelc.c the same way, expected.

### Pitfall 4: Recursive `xvpostscript_` calls from `effacer` (xvuelc.c:1432-1436)

**What goes wrong:** `effacer` calls `proc(xvpostscript)(&lasopsc)` with `lasopsc + 100`. The handleLasops body must not infinite-recurse. xvuelc.c:1286 handles this by `lasopsc = lasopsc - 100` and falling through. Port preserves this.
**Why it happens:** State machine has 5 modes (start `>0 && lasopsc==0`, normal `>0`, mode-100/101/102 `>100`, close `==0`, neg/erase `<-1`). The mode-100 branch resets state and re-opens TEMPORAIRE.EPS but does NOT fully recurse.
**How to avoid:** Plan 02's PsEmitter::handleLasops is a faithful port of the dispatch table. Plan 06 unit test calls `handleLasops(101)` and asserts post-state matches expected (lasopsc transition `101 → 1`, all i*=0, fpo reopened).
**Warning signs:** Stack overflow under solver load. Add an entry-counter assert in handleLasops for sanity (≤ 2 levels deep is the legitimate maximum).

### Pitfall 5: ImageMagick `convert` re-introduced via copy-paste

**What goes wrong:** A future plan author looks at `xvue/video1.f` (which uses `convert`) and copies the pattern into `xvue_qt_export.cpp` thinking it's idiomatic.
**Why it happens:** The legacy code is in the same repository, the dropdown of "MEFISTO ways to make a GIF" includes `convert` prominently.
**How to avoid:** Plan 06 owns a CI-style assertion: `! grep -rn 'convert\|ImageMagick\|magick' xvue/qt/` (excluding the "convert" keyword inside legitimate QPageSize / QString::convertToOther API names). Add to test suite as `test_no_imagemagick_in_qt`. Fail loudly. Run on every commit during Phase 7.
**Warning signs:** Any `system(...)` or `QProcess` call in `xvue/qt/` whose first arg starts with `convert`.

### Pitfall 6: ffmpeg not on PATH at runtime

**What goes wrong:** ffmpeg is in PATH at build time (build host has it) but the user installs MEFISTO on a minimal box without ffmpeg → GIF export fails silently with "Command not found" stderr noise that nobody sees because Qt swallows stderr unless the console dock is wired.
**Why it happens:** Runtime dependencies that aren't link-time deps are easy to lose.
**How to avoid:** Plan 05 PsEmitter::beginAnimation does a one-time `command -v ffmpeg` style check (`QProcess::execute("ffmpeg", {"-version"})` returning 0) and surfaces a clear FR/EN error in the console dock if missing: `"ffmpeg not found — install with: sudo apt install ffmpeg"` / `"ffmpeg introuvable — installer avec : sudo apt install ffmpeg"`. README install docs updated to mention `ffmpeg` as a runtime dep.
**Warning signs:** A user reports "GIF export does nothing" without an error message.

### Pitfall 7: HiDPI dimension mismatch between PNG export and PDF export

**What goes wrong:** On a HiDPI screen with `devicePixelRatioF=2`, `backing_` is 1600×1200 pixels for an 800×600 logical canvas. PNG export at backing dims → 1600×1200 (good — full resolution snapshot). PDF export at 72-dpi `setPageSize(QSizeF(800, 600), Point)` → 800×600 user units (good — vector-equivalent logical). But if the developer accidentally uses `backing_.width()` as the PDF page width, they get a 1600pt-wide PDF (wrong — too big for any printer).
**Why it happens:** `backing_.width()` returns physical pixels; canvas logical dims live in `XvueState`.
**How to avoid:** Plan 04 PDF code uses `XvueState::xpixels()` / `ypixels()` (logical) for `setPageSize`, and `backing_` for `drawPixmap` source. Document in xvue/qt/README.md Phase 7 section.
**Warning signs:** PDF that opens in Acrobat as 1600×1200 instead of the canvas's logical size.

### Pitfall 8: `auto-snapshot animation.gif overwrites without warning` clash with concurrent runs

**What goes wrong:** Two `pp/ppmail_qt` instances in the same `cwd` both auto-snapshot to `animation.gif` → race + corrupted GIF.
**Why it happens:** Per D-03, auto-snapshot fixed names overwrite without warning. Multi-instance is not Phase 7's invariant.
**How to avoid:** Document explicitly in xvue/qt/README.md Phase 7. v2 candidate: timestamp-suffix the filename (`animation_HHMMSS.gif`).
**Warning signs:** GIF file mid-write with a byte size of 0.

### Pitfall 9: Qt resource path collision with the bundled font

**What goes wrong:** Phase 3 already wired `qt_add_resources(xvueqt xvue_fonts ...)` for DejaVuSansMono.ttf. Phase 7 might be tempted to add icons or templates via a second `qt_add_resources` call. If the prefix collides (`/xvue/qt/...`), the resource lookup gets confused.
**Why it happens:** Phase 6 already added `xvue_icons.qrc` via the `.qrc-file route` (per CMakeLists.txt comment). Phase 7 doesn't need new resources, but if it did, follow Phase 6's pattern.
**How to avoid:** Plan 04 / 05 should not introduce new Qt resource files. PsEmitter is pure code; XvueExport is pure code; ffmpeg is external. No new `.qrc`.
**Warning signs:** AUTORCC linker errors about duplicate resource symbols.

## Code Examples

Verified patterns from official sources, copied here for the planner's quick reference. All snippets are illustrative; final form follows the `XVUE_QT_ASSERT_MAIN_THREAD()` + `XvueApp::ensure()` boilerplate established in Phases 1-6.

### Example 1: PNG snapshot (EXPORT-02)

```cpp
// Source: doc.qt.io/qt-6/qimage.html  +  Phase 4 backing_ semantics
void XvueExport::savePng(const QString& path) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto* canvas = XvueApp::ensure().window_slot()->canvas();
    QImage img = canvas->backing().toImage();
    if (!img.save(path, "PNG")) {
        // FR/EN status-bar + console-dock log
    }
}
```

### Example 2: JPEG snapshot (EXPORT-02)

```cpp
// Source: doc.qt.io/qt-6/qimagewriter.html
void XvueExport::saveJpeg(const QString& path) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto* canvas = XvueApp::ensure().window_slot()->canvas();
    QImageWriter w(path, "jpeg");
    int q = QSettings("MEFISTO","xvue").value("export/jpeg_quality", 90).toInt();
    w.setQuality(q);
    if (!w.write(canvas->backing().toImage())) { /* err handling */ }
}
```

### Example 3: PDF at canvas-native 72-dpi (EXPORT-05)

```cpp
// Source: doc.qt.io/qt-6/qpdfwriter.html   (verified via WebSearch 2026-05-01)
void XvueExport::savePdf(const QString& path) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto* canvas = XvueApp::ensure().window_slot()->canvas();
    auto* state  = canvas->state();

    QPdfWriter w(path);
    w.setResolution(72);  // user-units = points
    w.setPageSize(QPageSize(QSizeF(state->xpixels(), state->ypixels()),
                            QPageSize::Point));
    QPainter p(&w);
    p.drawPixmap(0, 0, canvas->backing());
    p.end();  // explicit; flushes / closes the PDF
}
```

### Example 4: ffmpeg PNG-sequence → animated GIF (EXPORT-03 / D-11)

```cpp
// Source: ffmpeg --help; QProcess::execute idiom
bool XvueExport::assembleGif(const QString& temp_dir, const QString& output_path) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    QStringList args;
    args << "-y"
         << "-framerate" << QString::number(QSettings("MEFISTO","xvue").value("export/gif_delay", 10).toInt())
         << "-i" << QDir(temp_dir).filePath("frame_%04d.png")
         << output_path;
    int rc = QProcess::execute("ffmpeg", args);
    return rc == 0;
}
```

### Example 5: PsEmitter::line (the canonical xvuelc.c port pattern)

```cpp
// Source: xvuelc.c:1942-1985 — PORTED VERBATIM (format strings unchanged)
void PsEmitter::line(int x1, int y1, int x2, int y2) {
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (lasopsc_ <= 0) return;          // matches `if (lasopsc > 0)` gate
    char buf[1024]; buf[0] = '\0';
    if (lasopsc_ < 3) {
        if (nbrcon_ == 0) {
            nbrcon_ = 1;
            xinic_ = x1; yinic_ = y1; xcouc_ = x2; ycouc_ = y2;
            if (counb_ != -1) {
                std::snprintf(concat_, sizeof concat_,
                    "%6i %6i %6i %6i %3i %4.2f %4.2f %4.2f %4.2f S\n",
                    x1, ypixels_ - y1, x2, ypixels_ - y2,
                    nbrcon_, courgb_[0], courgb_[1], courgb_[2], counb_);
            } else {
                std::snprintf(concat_, sizeof concat_,
                    "%6i %6i %6i %6i %3i %4.2f %4.2f %4.2f 0.00 S\n",
                    x1, ypixels_ - y1, x2, ypixels_ - y2,
                    nbrcon_, courgb_[0], courgb_[1], courgb_[2]);
            }
        } else {
            // ... continuation arm (xvuelc.c:1962-...): nbrcon++, etc.
        }
    } else {
        std::snprintf(chaine_[lasopsc_-4] + std::strlen(chaine_[lasopsc_-4]), …, …);
    }
}
```

### Example 6: QImageWriter probe (EXPORT-01 — Plan 01 deliverable)

```cpp
// xvue/qt/probes/qimagewriter_probe.cpp — 15 lines including #includes
#include <QGuiApplication>
#include <QImageWriter>
#include <iostream>
int main(int argc, char** argv) {
    QGuiApplication app(argc, argv);
    std::cout << "qt_version=" << qVersion() << "\n";
    std::cout << "supported_write_formats=";
    for (const auto& f : QImageWriter::supportedImageFormats()) {
        std::cout << f.toStdString() << " ";
    }
    std::cout << "\n";
    QImageWriter probe("/tmp/probe.gif", "gif");
    std::cout << "gif_write_supported=" << probe.canWrite() << "\n";
    std::cout << "gif_animation_supported="
              << probe.supportsOption(QImageIOHandler::Animation) << "\n";
    return 0;
}
```

**Empirical output on Debian trixie 2026-05-01** (verified by researcher running the equivalent program):
```
qt_version=6.10.2
supported_write_formats=bmp cur icns ico jfif jpeg jpg pbm pgm png ppm tif tiff wbmp webp xbm xpm
gif_write_supported=0
gif_animation_supported=0
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Hand-rolled PostScript via `fprintf(fpo, "%6i %6i ... S\n", ...)` | Same — preserved verbatim per EXPORT-04 | n/a | Phase 7 deliberately freezes the byte format. Users' downstream pipelines stay compatible. |
| ImageMagick `convert` for GIF assembly | ffmpeg via `QProcess::execute` | Phase 7 (D-11) | Drops one runtime apt dep on the Qt path. ffmpeg dithering quality matches or beats ImageMagick. |
| `xwd` + `convert` runtime shell-out (xvue/video1.f) | Stays alive on X11 backend; deferred to Phase 9 RETIRE-03 | n/a in Phase 7 | Documented gap. v1 GIF capture for `LVIDEO=1` solver workflows happens via X11 backend only until Phase 9 wires Qt equivalents. |
| Qt 5 QPrinter for PDF | `QPdfWriter` (Qt6 idiom) | Phase 7 (D-13) | Cleaner API. Both require `Qt6::Gui`. |
| `QImage::save("...gif")` for animation | ffmpeg dispatcher | Phase 7 (D-11) | `QImageWriter` doesn't write GIF on Debian's Qt 6.10.2; ffmpeg is the realized fallback. |

**Deprecated / outdated:**
- **Qt 5 QImageWriter GIF animation API.** Qt 5 docs and forum posts mention experimental multi-frame GIF write; in Qt 6, GIF write was removed entirely (the GIF plugin shipped with Qt 6 is read-only). This finding is **NEW for Phase 7** — the original CONTEXT decision tree assumed a probe-positive native path was likely. It is not. Researcher verified empirically.

## Assumptions Log

> Listed below are claims that this research could not verify with high confidence on the host system. The planner / discuss phase should confirm before locking decisions.

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | xvuelc.c lines 1187-1304 are ~117 LOC of `xvpostscript_` body, suitable for verbatim port. [VERIFIED via direct read in this session — 117 lines exactly between `void proc(xvpostscript)(int *lasops)` and the next `void`.] | "PostScript port granularity" | Low — already verified. |
| A2 | The ~150 in-primitive `if(lasopsc>0)` branches use stable format strings across xvuelc.c. [VERIFIED via read of xvuelc.c:1500-1568, 1740-1820, 1845-1975, 2050-2150, 2580-2780, 3140-3680. Found 14 distinct gates and 151 fprintf(fpo) sites.] | "PostScript per-primitive emit catalog" | Medium — Plan 03 must re-read each section before porting. The catalog below is research-grade but not a byte-by-byte transcription. |
| A3 | xvuelc.c's `lasopsc==3` chaine[-1] is preserved verbatim per D-05 even if it's a latent bug. [ASSUMED — researcher noted but did not verify with valgrind/asan whether the path is actually exercised in testa/.] | Pitfall 3 | Medium — if the bug is exercised in testa/wave or testa/cavity2d, the Qt port preserves the bug. Document and defer fix to a separate Phase. |
| A4 | `QFontInfo::family()` returns canonical family names that match the D-08 4-entry mapping table on the Qt side (e.g., `"Courier"`, not `"Courier New"`). [ASSUMED — depends on which Qt fonts are loaded on the user's system, which depends on `qt_add_resources xvue_fonts` at Phase 3 (DejaVuSansMono only — does NOT match any of the 4 entries).] | D-08 | Medium-high — DejaVuSansMono is the Phase 3 bundled font; under the Qt path **it falls through to /Helvetica fallback** every time. The D-08 mapping table is effectively unused under the default Qt config; warn-once log fires on every run. The PsEmitter port must NOT special-case this — fall through to /Helvetica is OK; document. |
| A5 | ffmpeg 8.1 produces visually-acceptable GIFs from a `frame_%04d.png` glob with default palette/dithering. [VERIFIED: `ffmpeg` is a mature tool; default GIF output is industry-standard. Risk is per-case visual A/B vs. the legacy convertepsgif output.] | D-11 | Low — fixable by adding `-vf "split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse"` filter chain if dithering ever looks bad. Defer until evidence. |
| A6 | `QPdfWriter` setResolution(72) + `QPageSize(QSizeF(xpixels, ypixels), QPageSize::Point)` produces a 1-page PDF whose user-units match the canvas pixel count. [CITED: doc.qt.io/qt-6/qpdfwriter.html — example pattern from forum.qt.io thread, MIT mirror, and qt6-base-dev source.] | D-15 | Low — pattern is well-attested. Validate empirically in Plan 04. |
| A7 | `bin/cb_probe_qt` script can be modeled on existing `bin/cbxvtest0_qt` (5-10 lines: cd $MEFISTO; mkdir build; cmake; make; run; redirect to PROBE.md). [VERIFIED via read of cbxvtest0_qt — pattern exists.] | D-09 | Low — pattern is established. |
| A8 | The Qt path's File menu → Export submenu (D-04) lands cleanly in Phase 6.0's shared shell, no per-module 6.1..6.5 plumbing needed. [VERIFIED via 06.0-CONTEXT.md D-04: shared shell builds File menu once; per-module phases register only their own menus alongside.] | D-04 | Low — Phase 6.0 architecture supports this. |
| A9 | `bin/convertepsgif` is "orphaned" / cyl53-specific. [VERIFIED: `grep -rn convert testa/ testf/ td/` returns zero matches.] | "Runtime State Inventory" / Pitfall 5 | Low — verified empirically. |
| A10 | `xvue/video1.f` + `xvue/videofin.f` + `xvue/videonm.f` (the LVIDEO pipeline) is in scope for **Phase 9 RETIRE-03**, NOT Phase 7. [ASSUMED — based on Phase 7 CONTEXT D-16 ("Phase 7 removes ImageMagick from xvue/qt/ ONLY") and the fact that these files are X11-side Fortran, not Qt code.] | "Runtime State Inventory" critical item | High — if the planner / user disagrees and wants Plan 05 to also Qt-ify the LVIDEO path, scope expands by 1-2 plans. **DISCUSS with user.** |
| A11 | `testa/wave` and `testa/cavity2d` produce <100 frames each in the typical interactive run. [ASSUMED — researcher did not run the tests.] | D-11 (sync ffmpeg) / Deferred async | Medium — if frame count is in the thousands, `QProcess::execute` blocks the GUI for many seconds. Plan 05 should include a `QProgressDialog` if observed. |

**If the user wants to lock A10 (LVIDEO defer), no further action needed.** **If the user wants A10 in-scope for Phase 7,** add a 7th plan: "LVIDEO Qt-side dispatch" — replaces `CALL SYSTEM('xwd ...')` and `CALL SYSTEM('convert ...')` with a Qt-side hook installed at `xvue_module_init_` time.

## Open Questions

1. **Should `LVIDEO` runtime GIF capture (xvue/video1.f) be in Phase 7 scope?**
   - What we know: It is a parallel, X11-only runtime GIF pipeline driven by 18+ Fortran tracer subroutines. CONTEXT D-16 locks "Phase 7 removes ImageMagick from xvue/qt/ ONLY" — strict reading defers LVIDEO to Phase 9.
   - What's unclear: Whether the user wants `pp/ppflui_qt`-driven solvers (which call `CALL VIDEO1(...)` from `flui/trvi2d.f` etc.) to produce GIFs in Phase 7 or to gracefully no-op until Phase 9.
   - Recommendation: **Ask user.** If user says "defer" → A10 is locked, no scope change. If user says "in scope" → add Plan 07 (LVIDEO Qt redirect).

2. **What is the empirical frame count of `testa/wave` and `testa/cavity2d` interactive runs?**
   - What we know: Researcher did not execute the tests.
   - What's unclear: Whether sync ffmpeg dispatch (D-11) blocks for ~1 second (10 frames) or 30 seconds (300 frames).
   - Recommendation: Plan 01 (or Plan 05) runs `pp/ppflui_qt testa/wave` once and counts auto-snapshot frames. Document in PROBE.md alongside the QImageWriter probe.

3. **Should `bin/convertepsgif`'s `-rotate 90 -extent 980x550` per-case transforms be replicated on the Qt path?**
   - What we know: `bin/convertepsgif` is cyl53-specific and orphaned (no caller in tracked tests).
   - What's unclear: Whether a future user manually invokes it as a post-processor on Qt-emitted .eps frames.
   - Recommendation: **Out of scope for v1** (not in REQUIREMENTS.md). Document in xvue/qt/README.md Phase 7 section: "Per-case post-processing (rotate, extent, crop) for the cyl53-style pipeline is not part of v1; users can still run `convertepsgif` manually until Phase 9 retires it."

4. **Should the Plan 01 PROBE.md kickoff commit also commit baseline `expected_legacy.gif` files for `testa/wave` and `testa/cavity2d`?**
   - What we know: No reference GIFs exist in `testa/` today. EXPORT-03 demands visual A/B vs legacy `bin/convertepsgif`.
   - What's unclear: Where the goldens should live (`testa/wave/expected/`, `.planning/phases/07-*/baseline/`, …) and whether the user wants binary blobs in git.
   - Recommendation: Commit goldens under `xvue/qt/tests/golden/wave_legacy.gif` and `cavity2d_legacy.gif` (~10 KB each); paths picked to match existing `xvue/qt/tests/` discipline.

5. **Should the verbatim port preserve the latent `chaine[-1]` underflow at `lasopsc==3` (Pitfall 3)?**
   - What we know: D-05 mandates verbatim port. The bug — if it is one — is exercised only when xvuelc.c receives `lasops==3` from a Fortran `LASOPS=3` write (qualities-menu mode).
   - What's unclear: Whether testa/wave or testa/cavity2d ever pass `LASOPS=3`.
   - Recommendation: Port verbatim per D-05; add a `[[deprecated("xvuelc.c parity — preserves possible buffer underflow")]]` annotation at the call site; defer fix to a separate non-Phase-7 effort.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| Qt6::Gui (QImageWriter, QPdfWriter, QPainter, QPixmap, QImage) | EXPORT-02, EXPORT-04, EXPORT-05 | ✓ | 6.10.2 | — |
| Qt6::Core (QProcess, QStandardPaths, QTemporaryDir, QSettings) | EXPORT-03, all | ✓ | 6.10.2 | — |
| Qt6::Widgets (QFileDialog, QMessageBox) | EXPORT-02..05 menu surface | ✓ | 6.10.2 | — |
| Qt6 GIF plugin write support (`gif` in `QImageWriter::supportedImageFormats()`) | D-10 native GIF path | ✗ | — | **D-11 ffmpeg fallback (verified working)** — this is the realized path on this system. |
| ffmpeg | D-11 GIF assembly | ✓ | 8.1-3+b1 | — |
| `xwd` (X11 screen capture) | xvue/video1.f legacy LVIDEO path | ✓ (X11 only) | system | n/a — Qt path bypasses xwd entirely; LVIDEO defers to Phase 9. |
| `convert` (ImageMagick) | xvue/video1.f, bin/convertepsgif (X11 path) | ✓ (X11 only) | system | n/a — Qt path forbids ImageMagick (D-11, EXPORT-06). |
| `qt6-imageformats-plugins` apt package | Hypothetical Fedora-style add-on for GIF write | ✗ | — | Does not exist on Debian trixie; the GIF plugin ships in `libqt6gui6` and is read-only. |
| CMake 3.21+ | xvue/qt/probes/CMakeLists.txt | ✓ | system | — |

**Missing dependencies with no fallback:** None — D-11 ffmpeg fallback covers the only "missing" capability (GIF writing).

**Missing dependencies with fallback:**
- Qt6 GIF plugin write support → ffmpeg fallback (D-11). This is the realized path.

## Validation Architecture

> Required because `workflow.nyquist_validation: true` in `.planning/config.json`.

### Test Framework

| Property | Value |
|----------|-------|
| Framework | QTest (Qt 6 testing module) + ad-hoc shell smoke for probe binary |
| Config file | `xvue/qt/tests/CMakeLists.txt` (existing; Phase 7 appends two `add_executable` blocks) |
| Quick run command | `cmake --build xvue/build --target xvue_qt_postscript_tests xvue_qt_export_tests && (cd xvue/build && xvfb-run ctest -R '^xvue_qt_(postscript\|export)_tests$' -V)` |
| Full suite command | `cmake --build xvue/build && (cd xvue/build && xvfb-run ctest -V)` (Phase 6.x: ~155 tests; Phase 7 adds ~15-20 new) |

### Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| EXPORT-01 | PROBE.md committed at kickoff with current `QImageWriter::supportedImageFormats()` output | smoke | `bash bin/cb_probe_qt && diff <(cat .planning/phases/07-image-gif-and-postscript-export/PROBE.md) - <<<expected` (Plan 01 owns) | ❌ Wave 0 — NEW probe binary + script |
| EXPORT-02 | PNG round-trip: drawn QPixmap → save → load → pixel-compare ≤ 1 bit/channel drift | unit | `xvfb-run xvue_qt_export_tests --gtest_filter=*PNG*` | ❌ Wave 0 — NEW `xvue_qt_export_tests` binary |
| EXPORT-02 | JPEG round-trip: same with quality=90 | unit | `xvfb-run xvue_qt_export_tests --gtest_filter=*JPEG*` | ❌ Wave 0 |
| EXPORT-03 | GIF assembly: 5 dummy QImage frames → ffmpeg → output is non-empty .gif of the right frame count | integration | `xvfb-run xvue_qt_export_tests --gtest_filter=*GIF*` (skips if !command -v ffmpeg) | ❌ Wave 0 |
| EXPORT-03 | Visual A/B vs legacy `bin/convertepsgif` on testa/wave (frame-count + first-frame md5 hash compare) | integration / manual | `bin/qt-test-gif-ab.sh testa/wave` (Plan 06 helper, runs Qt path then diff vs committed `xvue/qt/tests/golden/wave_legacy.gif`) | ❌ Wave 0 |
| EXPORT-04 | PsEmitter byte-level golden compare: drive a fixed sequence of Fortran-equivalent calls into PsEmitter, capture stdout/file output, byte-compare against committed `golden/scene01.eps` | unit | `xvfb-run xvue_qt_postscript_tests --gtest_filter=*Golden*` | ❌ Wave 0 — NEW `xvue_qt_postscript_tests` binary + golden files |
| EXPORT-04 | PsEmitter state machine: handleLasops(0/1/2/3/-3/100/101) transitions match xvuelc.c expected post-state for each | unit | `xvfb-run xvue_qt_postscript_tests --gtest_filter=*StateMachine*` | ❌ Wave 0 |
| EXPORT-04 | PsEmitter Y-flip: `psEmitter().line(10, 20, ...)` with ypixels=600 emits `... 580 ... S\n` | unit | `xvfb-run xvue_qt_postscript_tests --gtest_filter=*YFlip*` | ❌ Wave 0 |
| EXPORT-05 | PDF size + aspect: write 800×600 canvas to PDF; reload via `QPdfDocument`; assert page size = 800pt × 600pt | unit | `xvfb-run xvue_qt_export_tests --gtest_filter=*PDFGeometry*` (requires `Qt6::Pdf` if used; otherwise binary-poke at /Width /Height in raw .pdf) | ❌ Wave 0 |
| EXPORT-06 | No ImageMagick reference in `xvue/qt/`: `! grep -rn 'convert\|magick\|ImageMagick' xvue/qt/src/ xvue/qt/include/` returns nothing | smoke | `bash bin/test_no_imagemagick_in_qt.sh` (Plan 06 owns; CI-style) | ❌ Wave 0 |
| EXPORT-06 | Probe binary exit 0 + parsable stdout | smoke | `xvue/build/probes/qimagewriter_probe \| head -3` | ❌ Wave 0 |

### Sampling Rate

- **Per task commit:** `cmake --build xvue/build --target xvue_qt_postscript_tests xvue_qt_export_tests && (cd xvue/build && xvfb-run ctest -R 'xvue_qt_(postscript|export)' -V)` — 2 binaries, ~30s.
- **Per wave merge:** Full new-test sweep + grep gate: `cmake --build xvue/build && (cd xvue/build && xvfb-run ctest -V) && bash bin/test_no_imagemagick_in_qt.sh`.
- **Phase gate:** Full xvue suite green (the existing 155 + the new 20) before `/gsd-verify-work`. Plus visual A/B on `testa/wave` and `testa/cavity2d` (manual eye-compare GIFs side-by-side).

### Wave 0 Gaps

- [ ] `xvue/qt/probes/qimagewriter_probe.cpp` — Plan 01 (10-20 LOC).
- [ ] `xvue/qt/probes/CMakeLists.txt` — Plan 01.
- [ ] `bin/cb_probe_qt` (chmod 755) — Plan 01.
- [ ] `xvue/qt/tests/test_xvue_qt_postscript.cpp` — Plan 06 (golden compare + state machine + Y-flip slots).
- [ ] `xvue/qt/tests/test_xvue_qt_export.cpp` — Plan 06 (PNG/JPEG/PDF/GIF round-trip slots).
- [ ] `xvue/qt/tests/golden/scene01.eps` — Plan 06 (committed byte-exact reference; produced by running the X11 backend on a 30-line test driver before Phase 7's emitter ships).
- [ ] `xvue/qt/tests/golden/wave_legacy.gif` — Plan 01 OR Plan 06 (committed; produced by running X11 backend + bin/convertepsgif on testa/wave).
- [ ] `xvue/qt/tests/golden/cavity2d_legacy.gif` — Plan 01 OR Plan 06.
- [ ] `bin/test_no_imagemagick_in_qt.sh` — Plan 06 (5-line `grep -rn 'convert\|magick' xvue/qt/ && exit 1 || exit 0`).
- [ ] CMakeLists.txt updates: `xvue/qt/CMakeLists.txt` adds `src/xvue_qt_postscript.cpp`, `src/xvue_qt_export.cpp`, `add_subdirectory(probes)`. `xvue/qt/tests/CMakeLists.txt` adds two `add_executable` blocks following the Phase 6.x pattern.

## Security Domain

> `security_enforcement` is not explicitly false in `.planning/config.json` — treat as enabled.

### Applicable ASVS Categories

| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | no | MEFISTO is a single-user desktop scientific tool; no remote / multi-user surface. |
| V3 Session Management | no | No sessions. |
| V4 Access Control | no | Single-user desktop. |
| V5 Input Validation | yes (limited) | File paths from `QFileDialog::getSaveFileName` are user-typed but Qt validates them; ffmpeg argv is constructed from QStringList (no shell interpolation). PostScript filename `nomfichier` from Fortran is passed length-checked (`*length` int + `nomfichier[]` char) — bounds check before use. |
| V6 Cryptography | no | No keys, no signatures, no encrypted blobs. PostScript / PDF / GIF / PNG / JPEG outputs are plaintext / well-known formats. |
| V7 Errors & Logging | yes | Error paths must surface to console dock + status bar (D-09 Phase 6.0 console-dock surface inheritance). No stack traces / paths leaked to status bar. |
| V12 Files & Resources | yes | Temp dir creation (`QTemporaryDir`) — RAII, secure perms 0700 by default. PostScript file overwrite without warning (D-03 / xvuelc.c parity) — documented user-facing behavior. |
| V13 API & Web Service | no | No HTTP / API surface. |

### Known Threat Patterns for {Qt 6 desktop + Fortran ABI + ffmpeg shell-out}

| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Path traversal via Fortran-passed `nomfichier` (e.g., `../../etc/passwd.eps`) | Tampering | Qt-side: validate via `QFileInfo::canonicalFilePath()` resolves under `cwd` or under user-confirmed `QFileDialog`-returned dir. The verbatim PsEmitter port preserves xvuelc.c's behavior (it writes to `cwd`-relative `TEMPORAIRE.EPS`, which is fine in a single-user desktop tool). For `xvsauverps(nomfichier, length)` → user-named: validate length ≤ 255 and reject `..` segments. **Defer to Plan 02 with TODO note**: bounds-check `*length` before sprintf'ing into `format[255]`. |
| ffmpeg argv injection via QStringList | Tampering | Mitigated by `QProcess::execute("ffmpeg", QStringList{...})` (no shell). Researcher confirms each argv element is a controlled constant or a path that originated from `QTemporaryDir` (random-suffix). No user-typed args in v1. |
| Disk-fill via uncapped frame accumulator | DoS | Plan 05: soft cap auto-snapshot at 10000 frames; warn at 1000. Fail-closed: when temp dir reports < 100 MB free, halt and console-dock log the issue. |
| `exit(1)` on fopen failure (xvuelc.c:1234, 1282) | DoS / data-loss (loses unsaved Fortran state) | Replace with `QMessageBox::critical` + clean shutdown via `XvueApp::ensure().quitGracefully()` so unsaved project state can be flushed. Documented in CONTEXT specifics. |
| Race on `animation.gif` between two Qt processes in same cwd | Tampering | Out-of-scope for v1 (single-instance model). Document in xvue/qt/README.md. v2: timestamp-suffix. |
| `TEMPORAIRE.EPS` left on disk after crash | Information disclosure (low — it's the user's own diagram) | Cosmetic. xvpostscript_(0) close removes it. v2: register `atexit` cleanup. Acceptable in v1. |
| Untrusted PostScript / PDF / GIF / PNG output that downstream tools (gs, lpr, browser) parse | n/a (output, not input) | The user generated the output from their own scene; no untrusted input enters Phase 7's path. |
| QProcess hanging on bad ffmpeg invocation | DoS (GUI hang) | Set a `QProcess::execute` timeout via the new Qt6 overload, or wrap in a `QFutureWatcher` if frame count > N. v1: live with it (small frame counts). Document. |

**Threat model summary:** Phase 7 is a low-threat phase. The biggest real risk is the `exit(1)` call site (DoS / data-loss for the Fortran user) — the migration plan should replace it with a clean shutdown path. Path-traversal in Fortran-passed filenames is theoretically possible but in practice constrained by single-user desktop semantics. ffmpeg is invoked with controlled constants only.

## Sources

### Primary (HIGH confidence)

- `xvue/xvuelc.c:170-189` (file-static PS state declarations) — direct read 2026-05-01
- `xvue/xvuelc.c:1187-1304` (xvpostscript_ body, ~117 LOC) — direct read 2026-05-01
- `xvue/xvuelc.c:1414-1437` (effacer + xvpostscript recursive call site) — direct read
- `xvue/xvuelc.c:1500-1568` (xvchargefonte_ PS branch) — direct read
- `xvue/xvuelc.c:1740-1820` (xvtrait_ + xvtypetrait_ PS branches) — direct read
- `xvue/xvuelc.c:1845-1985` (line-drawing PS-emit block with ypixels-y flips) — direct read
- `xvue/xvuelc.c:2050-2150` (xvface / xvfacetraits PS branches) — direct read
- `xvue/xvuelc.c:2580-2780` (rectangle / ellipse PS branches) — direct read
- `xvue/xvuelc.c:3020-3680` (PS header / dictionary / footer) — direct read
- `xvue/qt/src/xvue_qt_api.cpp:607-614` (current xvpostscript_ stub) — direct read
- `xvue/qt/CMakeLists.txt` — direct read (163 LOC)
- `xvue/qt/tests/CMakeLists.txt` — direct read (existing pattern for Plan 06)
- `xvue/video1.f`, `xvue/videofin.f`, `xvue/videonm.f` — direct read (LVIDEO pipeline)
- `bin/convertepsgif`, `bin.lnx64/convertepsgif` — direct read (single line each)
- `incl/trvari.inc`, `incl/trvari.inc95` — direct read (LVIDEO common block)
- Empirical probe: ran a `QImageWriter::supportedImageFormats()` test program on the host (Qt 6.10.2 / Debian trixie 2026-05-01) — proved `gif` is absent from write formats and `QImageWriter::write` for GIF returns "Unsupported image format". This is the EXPORT-01 outcome.
- `pkg-config --modversion Qt6Core` returned `6.10.2` — direct shell invocation
- `dpkg -S /usr/lib/x86_64-linux-gnu/qt6/plugins/imageformats/libqgif.so` returned `libqt6gui6:amd64`
- `which ffmpeg && ffmpeg -version` returned `8.1-3+b1`
- `.planning/phases/07-image-gif-and-postscript-export/07-CONTEXT.md` — direct read (locked decisions)
- `.planning/REQUIREMENTS.md` (EXPORT-01..06, lines 88-95) — direct read
- `.planning/ROADMAP.md` Phase 7 section (lines 230-240) — direct read
- `.planning/STATE.md` — direct read
- `xvue/README_COORDS.md` — direct read (PostScript Y-flip mandate)
- `xvue/qt/README.md` — direct read (Phase 5/6 architecture inheritance)
- `.planning/codebase/STACK.md`, `STRUCTURE.md`, `CONVENTIONS.md`, `CONCERNS.md` — direct read
- `.planning/PROJECT.md` — direct read (constraints + key decisions)

### Secondary (MEDIUM confidence)

- [QPdfWriter Class | Qt GUI 6.11](https://doc.qt.io/qt-6/qpdfwriter.html) — referenced for setResolution/setPageSize semantics; specific 72-dpi+QSizeF code pattern verified via [forum.qt.io thread](https://forum.qt.io/topic/144716/qimage-qpainter-qpdfwriter-newpage-working-code) and [MIT mirror](https://stuff.mit.edu/afs/athena/software/texmaker_v5.0.2/qt57/doc/qtgui/qpdfwriter.html).
- [QImageWriter Class | Qt GUI 6.11](https://doc.qt.io/qt-6/qimagewriter.html) — referenced for `supportedImageFormats()` semantics (cross-checked with empirical probe).
- [QImageIOHandler Class | Qt GUI 6.11](https://doc.qt.io/qt-6/qimageiohandler.html) — `Animation` enum and `supportsOption()`.

### Tertiary (LOW confidence — flagged)

- [Qt forum: Animated GIF QImage Libav API](https://forum.qt.io/topic/76185/animated-gif-qimage-libav-api) — community consensus that `QImageWriter` does not support multi-frame GIF write in Qt 5/6 (cross-checked with empirical probe — confirmed for Qt 6.10.2).
- [GitHub: dbzhang800/QtGifImage](https://github.com/dbzhang800/QtGifImage) — third-party library for animated GIF write; mentioned for completeness, NOT used in Phase 7 (locked decision is ffmpeg).

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — Qt 6.10.2, ffmpeg 8.1, all already present on the host; verified via dpkg/pkg-config/which.
- Architecture: HIGH — pattern (PsEmitter + ffmpeg dispatcher + File menu) follows Phases 5/6 templates verbatim.
- Pitfalls: MEDIUM — the per-primitive emit catalog risk (Pitfall 1) and the latent `chaine[-1]` underflow (Pitfall 3) require Plan 02/03 author re-read of xvuelc.c.
- GIF probe outcome: HIGH — empirically verified Qt 6.10.2 does NOT write GIF on Debian trixie.
- Runtime state inventory: HIGH — found by direct grep + read of xvue/video1.f; LVIDEO pipeline is real and is **not** redundant with bin/convertepsgif. Item A10 awaits user disposition.

**Research date:** 2026-05-01
**Valid until:** ~2026-08-01 (3 months — Qt 6.10 is stable; ffmpeg 8.1 is stable; Debian trixie is stable. Re-validate if Qt 6.11 or 7.0 lands or if `qt6-imageformats-plugins` becomes available on Debian.)

---

*Research completed: 2026-05-01*
*Phase: 07-image-gif-and-postscript-export*
*Researcher: GSD phase researcher*
