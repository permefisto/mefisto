---
phase: 07-image-gif-and-postscript-export
plan: 04
subsystem: xvue/qt
tags: [qt6, export, png, jpeg, pdf, qfiledialog, qpdfwriter, qsettings, i18n, file-menu]

requires:
  - phase: 07-02-psemitter-skeleton
    provides: PsEmitter singleton accessor (XvueApp::psEmitter()) — Plan 04 mirrors this pattern for XvueExport(). XvueCanvas::backing_ snapshot source unchanged.
provides:
  - XvueExport class (xvue/qt/src/xvue_qt_export.{h,cpp}) — PNG/JPEG/PDF save methods on backing_
  - XvueApp::xvueExport() singleton accessor (mirrors psEmitter()/eventBridge()/menuBridge() pattern)
  - File → Export submenu (PNG…/JPEG…/PDF…) wired into XvueWindow::buildMenuBar (Phase 6.0 shared shell builder); inherited by all five pp*_qt
  - QFileDialog + QSettings("xvue/export/last_dir") last-directory persistence
  - Raster PDF via QPdfWriter::setResolution(72) + setPageSize(QPageSize(QSizeF(xpixels,ypixels), QPageSize::Point)) + drawPixmap(0,0,backing_) — canvas-native page geometry per D-15
  - xvue_qt_i18n.{h,cpp} extended: MsgId::FileExport{Png,Jpeg,Pdf}{,Title,Failure,Success} entries with FR/EN translations; static_assert against MsgId::_Count_ keeps tables in sync
  - Test target xvue_qt_export_tests (xvue/qt/tests/test_xvue_qt_export.cpp): 5 GTest+QTest slots covering PNG/JPEG round-trip, PDF page geometry, write-failure handling
affects: [07-05-gif-ffmpeg-fallback, 07-06-validation-ab]

tech-stack:
  added:
    - QPdfWriter (Qt6::PrintSupport implicit) — raster PDF emit
    - QImageWriter::write — PNG/JPEG single-frame snapshot
    - QFileDialog::getSaveFileName — interactive filename selection
    - QSettings("xvue/export/last_dir") — last-directory persistence across sessions
  patterns:
    - "Bridge-as-singleton-on-XvueApp: XvueExport joins eventBridge()/menuBridge()/psEmitter() as a lazily-allocated XvueApp-owned QObject."
    - "Menu wiring in shared shell builder: Phase 6.0 D-04 contract — File → Export submenu added once in XvueWindow::buildMenuBar; all five pp*_qt executables inherit it. No per-module register*Actions changes."
    - "Bilingual every user-facing string: MsgId enum entries paired in FR/EN translation tables; xvueT() lookup via existing xvueIsEnglish() flag (Phase 6.0 D-09 pattern)."
    - "Raster PDF for v1: D-13 — QPdfWriter::drawPixmap(0,0,backing_) embeds canvas pixmap. True-vector PDF deferred to v2 (would require display-list capture across every QPainter call)."
    - "Canvas-native PDF page geometry at 72 dpi: D-15 — page aspect matches canvas exactly. xpixels/ypixels are LOGICAL (NOT backing_ pixel dims; backing_ is devicePixelRatio-scaled)."

key-files:
  created:
    - xvue/qt/src/xvue_qt_export.h (47 lines — XvueExport class declaration with saveAsPng/saveAsJpeg/saveAsPdf slots)
    - xvue/qt/src/xvue_qt_export.cpp (267 lines — PsEmitter-style XvueExport body; QImageWriter PNG/JPEG; QPdfWriter PDF at 72 dpi canvas-native; QFileDialog + QSettings; bilingual error surfaces)
    - xvue/qt/tests/test_xvue_qt_export.cpp (257 lines — 5 GTest+QTest slots: PNG round-trip, JPEG round-trip, PDF geometry, write-failure logging, i18n MsgId coverage)
  modified:
    - xvue/qt/src/xvue_qt_app.h, .cpp (xvueExport() accessor + static unique_ptr export_; reset() in teardown_atexit)
    - xvue/qt/src/xvue_qt_window.h, .cpp (replaced Phase 6.0 placeholder QAction with full submenu: actExportPng_/Jpeg_/Pdf_ + onFileExport{Png,Jpeg,Pdf} slots delegating to XvueApp::xvueExport())
    - xvue/qt/src/xvue_qt_i18n.h, .cpp (added FileExportPng/Jpeg/Pdf + Title/Failure/Success MsgId entries with FR/EN pairs)
    - xvue/qt/CMakeLists.txt (added src/xvue_qt_export.cpp + Qt6::PrintSupport link; tests subdir already enabled)
    - xvue/qt/tests/CMakeLists.txt (added xvue_qt_export_tests target + add_test ctest registration)

key-decisions:
  - "XvueExport ownership on XvueApp singleton (D-01-aligned): zero new extern \"C\" symbols. ABI stays at 58 (verify_abi.sh exit 0). XvueExport is C++ internal; menu actions trigger Qt slots that call XvueApp::xvueExport() — no Fortran-driven export path is added."
  - "QPdfWriter at canvas-native 72 dpi (D-15, Pitfall 7): xpixels = canvas_->width(), ypixels = canvas_->height() — LOGICAL pixels from XvueCanvas, NOT backing_.width()/height() (which are devicePixelRatio-scaled). Resulting PDF page aspect matches canvas; no fit-to-A4 distortion."
  - "QFileDialog absolute path is always-trusted (T-07-02 mitigation): XvueExport never concatenates a Fortran-passed prefix to the QFileDialog return value. Plan 05's fixed zfxy0NNN.png Fortran-auto path will use a hardcoded cwd-relative pattern, also non-concatenated."
  - "QSettings('xvue/export/last_dir') persists across sessions: getSaveFileName picks up last-used dir; on success, parent directory of saved file is written back. Standalone QSettings store, not coupled to other xvue/* keys."
  - "Bilingual MsgId entries with static_assert: every new FileExport* entry has FR + EN translations in xvue_qt_i18n.cpp; static_assert(static_cast<size_t>(MsgId::_Count_) == translations_count) at file scope catches any FR/EN drift at compile time."
  - "GIF deliberately NOT in Plan 04: Plan 05 owns GIF wiring (probe-driven dispatch, ffmpeg fallback per D-11, auto-snapshot per D-02). Plan 04's submenu is built with three children only; Plan 05 will append GIF + Capture-Animation entries."

patterns-established:
  - "MsgId vs xvueT(): all new user-visible strings flow through this table — never inline literal strings into UI code."
  - "Singleton accessor + reset() in teardown_atexit: pattern is now four-deep (eventBridge / menuBridge / psEmitter / xvueExport) — future Qt-side singletons should follow."
  - "QFileDialog + QSettings last-dir: reusable export-target convention for any future Save As… surface."

threat-mitigations:
  - "T-07-02 (path-traversal via Fortran-passed name): XvueExport menu path uses QFileDialog return verbatim; no Fortran-string concatenation. Plan 05's auto-path uses hardcoded zfxy0NNN.png in cwd, also not concatenated."
  - "T-07-03 (write-failure handling): QImageWriter::write / QPdfWriter errors surface via QMessageBox::warning when menu-triggered + console-dock log. No silent failure."

tests:
  added:
    - "xvue_qt_export_tests (5 slots: testSaveAsPngRoundTrip, testSaveAsJpegRoundTrip, testSaveAsPdfCanvasNativeGeometry, testWriteFailureLogged, testFileExportI18nFrEnPaired)"
  command: "cd xvue/qt/build && xvfb-run --auto-servernum ctest -R '^xvue_qt_export_tests$' --output-on-failure"
  result: "5/5 PASS in 0.31 sec"

build:
  qt: "cmake --build xvue/qt/build --target xvueqt -j4 → exit 0"
  abi: "verify_abi.sh → 58 symbols (unchanged)"
  full-tree: "bin/cbl_tout → exit 0 (X11 non-regression confirmed); bin/cbl_tout_qt → exit 0"

deviations:
  - "Plan 04 Task 2 was finalized post-cap-reset by orchestrator: agent committed RED + GREEN + Window-menu-edit work-in-progress before the global usage cap interrupted on 2026-05-04 ~22:35 CEST. Orchestrator inspected uncommitted xvue_qt_window.{h,cpp} edits (Export submenu wiring + onFileExport{Png,Jpeg,Pdf} slot bodies — all already drafted by the agent), built+tested them green, committed as 7dd93da, and authored this SUMMARY.md from the diff inspection. No code authored by the orchestrator beyond the commit message. Documented for audit."

requirements:
  - EXPORT-02 (PNG/JPEG snapshot via XvueExport)
  - EXPORT-05 (PDF additive bonus via QPdfWriter::PdfFormat — see D-13/D-14/D-15 in CONTEXT.md)

duration: ~14m21s (agent before cap) + ~3m (orchestrator finalize)

next-steps:
  - "Plan 05 (Wave 4): animated GIF — adds GIF + Capture-Animation menu entries to actExportSubmenu_, hooks PsEmitter::handleLasops(*lasops==0) for auto-snapshot, ffmpeg sync QProcess fallback per D-11."
  - "Plan 03 (Wave 3, partial): per-primitive PsEmitter helpers + 15 wiring sites + golden tests — needs re-spawn with 07-03-WIP.diff context (saved to phase dir)."
