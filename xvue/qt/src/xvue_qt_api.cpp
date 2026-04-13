// xvue/qt/src/xvue_qt_api.cpp
// Phase 0: warn-once no-op stubs for every Fortran-facing extern "C" entry
// point in xvue_qt_api.h. Each stub prints a single stderr line on first
// call, then is silent for the rest of the process.
//
// Per D-18: do NOT abort, do NOT set error flags, do NOT touch any Qt object.
// Non-void returns use the simplest safe default (nullptr, 0.0).
//
// Per D-17: warn-once diagnostic pattern is a per-function `static bool warned`.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include "xvue_qt_api.h"
#include "xvue_qt_app.h"
#include "xvue_qt_window.h"
#include "xvue_qt_canvas.h"
#include "xvue_qt_state.h"
#include <QApplication>
#include <QCoreApplication>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QFont>
#include <QFontMetrics>
#include <QGuiApplication>
#include <QLatin1Char>
#include <QPainter>
#include <QPixmap>
#include <QPoint>
#include <QPolygon>
#include <QRect>
#include <QScreen>
#include <QString>
#include <QWindow>
#include <vector>

namespace {

inline void warn_once(bool &flag, const char *name) {
    if (!flag) {
        std::fprintf(stderr, "xvue-qt: stub %s not implemented yet\n", name);
        flag = true;
    }
}

// File-local helper (D-13). Not part of the Fortran ABI.
// All four xv*rectangle_ entry points route through this one implementation.
//
// CR-01 (2026-04-11): legacy X11 distinguishes outline-only
// (XDrawRectangle: xvfbordrectangle_/xvbordrectangle_) from fill-only
// (XFillRectangle: xvfrectangle_/xvrectangle_). QPainter::drawRect strokes
// AND fills in one call, so the four symbols must dispatch on an explicit
// mode and use Qt::NoBrush / fillRect to keep the legacy semantics when
// pen and brush colors diverge (Phase 3).
enum class RectMode { Outline, Fill };

inline void xvue_qt_draw_rect_common(int x, int y, int w, int h, RectMode mode) {
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;
    const QRect r(x, y, w, h);
    if (mode == RectMode::Outline) {
        // XDrawRectangle: stroke only, no fill.
        QBrush saved = st->painter_->brush();
        st->painter_->setBrush(Qt::NoBrush);
        st->painter_->drawRect(r);
        st->painter_->setBrush(saved);
    } else {
        // XFillRectangle: fill only, no outline.
        st->painter_->fillRect(r, st->painter_->brush());
    }
    if (win->canvas()) win->canvas()->update();
}

// Phase 3 D-06 + Pitfall 7: xvtexte_ and xvftexte_ share this body under
// the Phase 2 single-backing model (02/D-05). In legacy xvuelc.c the two
// drew to DIFFERENT targets (mempx at :1658 vs fenetre_mef at :1678); the
// collapse here is legal ONLY because Phase 2 unified the draw target.
// If a future phase reintroduces multiple draw targets, revisit.
inline void xvue_qt_draw_text_common(const char* string, int length,
                                     int x1, int y1) {
    if (!string || length <= 0) return;              // WR-03 null-guard
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive() || !st->backing_) return;

    // Pitfall 4: TWO-ARG fromLatin1 — the legacy Fortran strings are NOT
    // NUL-terminated; explicit length is mandatory.
    QString qstr = QString::fromLatin1(string, length);
    st->painter_->setFont(st->current_font_);        // cheap no-op if unchanged
    st->painter_->drawText(x1, y1, qstr);            // baseline form (D-06)

    // Phase 2 D-01 epilogue: drawing primitives flush.
    if (win->canvas()) win->canvas()->update();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

} // anonymous namespace

extern "C" {

// ---- 1. languemefisto_ ----
void proc(languemefisto)(int *langue) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "languemefisto_");
    (void)langue;
}

// ---- 2. dctnmc_ (returns void*) ----
void *proc(dctnmc)(int *nboctets) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "dctnmc_");
    (void)nboctets;
    return nullptr;
}

// ---- 3. dstnmc_ ----
void proc(dstnmc)(void *mcoctets) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "dstnmc_");
    (void)mcoctets;
}

// ---- 4. nomrepmefisto_ ----
void proc(nomrepmefisto)(char *chaine, int *size) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "nomrepmefisto_");
    (void)chaine; (void)size;
}

// ---- 5. xvinitgraphique_ ----
void proc(xvinitgraphique)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();

    auto& win = XvueApp::window_slot();
    if (!win) {
        // D-07: lazy (re)allocation. QApplication is NOT recreated — call_once
        // guarantees exactly-one. D-02: title "MEFISTO", 800x600 set in ctor.
        win = std::make_unique<XvueWindow>();
    }
    win->show();
    win->raise();
    win->activateWindow();

    // D-01 revised 2026-04-11 (debug phase-01-xvtest0-teardown-segfault):
    // the original "exactly one processEvents" rule is insufficient on X11
    // to actually realize a window (MapRequest, ConfigureNotify, Expose
    // each need their own trip through the event loop). We pump in a
    // bounded loop until QWindow::isExposed() is true or the 2 s timeout
    // fires (so a broken display cannot hang the Fortran main). The
    // ExcludeUserInputEvents flag remains mandatory (Pitfall 4).
    QElapsedTimer t;
    t.start();
    const int timeout_ms = 2000;
    QWindow* wh = win->windowHandle();
    while (t.elapsed() < timeout_ms) {
        QCoreApplication::processEvents(
            QEventLoop::ExcludeUserInputEvents, 20);
        wh = win->windowHandle();  // may become non-null after first pump
        if (wh && wh->isExposed()) break;
    }
}

// ---- 6. xtinit_ ----
// Phase 03.1 Pitfall 1 fix: promoted from warn-once stub to a real
// body so that XVOUVRIR -> XTINIT -> (XVINIT chain) opens a Qt window
// on the Qt backend the same way XVINITGRAPHIQUE does on the xvtest0
// code path. Idempotency is provided by window_slot()'s "if (!win)"
// lazy alloc — calling XVINITGRAPHIQUE after XTINIT (or in any order)
// does NOT double-create. XvueApp::ensure() is guarded by call_once
// per D-07. See .planning/phases/03.1-.../03.1-RESEARCH.md §Pitfall 1.
//
// The legacy X11 backend (xvue/xvuelc.c) retains its own independent
// xtinit_ body — this fix is Qt-side only, preserving BUILD-07
// bit-identity on the legacy path.
void proc(xtinit)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    // Forward to xvinitgraphique_'s body: lazy-alloc window_slot,
    // show/raise/activate, pump the bounded exposure loop. This is
    // idempotent — a subsequent xvinitgraphique_ call from any other
    // Fortran path will reuse the existing XvueWindow via the
    // "if (!win)" guard inside xvinitgraphique_.
    proc(xvinitgraphique)();
}

// ---- 7. xvpxecran_ ----
void proc(xvpxecran)(int *xp, int *yp) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    // D-16: logical pixels from the primary screen, callable before
    // xvinitgraphique_ (D-18). Multi-monitor awareness deferred.
    QScreen* s = QGuiApplication::primaryScreen();
    if (s && xp && yp) {
        *xp = s->size().width();
        *yp = s->size().height();
    }
}

// ---- 8. xvmmecran_ ----
void proc(xvmmecran)(int *xmm, int *ymm) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    // D-17: physical millimetres from QScreen::physicalSize(). Qt 6 guarantees
    // mm regardless of DPR.
    QScreen* s = QGuiApplication::primaryScreen();
    if (s && xmm && ymm) {
        QSizeF mm = s->physicalSize();
        *xmm = static_cast<int>(mm.width()  + 0.5);
        *ymm = static_cast<int>(mm.height() + 0.5);
    }
}

// ---- 9. initaccrochage_ ----
void proc(initaccrochage)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "initaccrochage_");
}

// ---- 10. xvinfo_ ----
void proc(xvinfo)( int *ix, int *iy, int *maxfonts,
                   int *n1coref, int *ndcoref, int *n1coelf,
                   int *ndcoelf, int *n1coulf, int *ndcoulf, int *nbcolo,
                   char namefonts[][256], int nbchar[], int *nbfonts, int *visuclass ) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();

    // D-03: window resize portion from Phase 1 stays.
    auto& win = XvueApp::window_slot();
    if (win && ix && iy) {
        win->resize(*ix, *iy);
    }

    // Phase 3 D-22/D-23: real font enumeration + palette ranges.
    // Pitfall 2 ordering: palette_init_once() already ran inside XvueState
    // ctor (via XvueWindow ctor via xvinitgraphique_), so any caller that
    // reads red/green/blue immediately after xvinfo_ returns (e.g.
    // xvue/xvinit.f:143 calling XVRECUPRGBDEC) sees populated data.
    static const char* const kListFonts[XvueState::kNbFonts] = {
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
    if (maxfonts && *maxfonts > 0 && *maxfonts < nfonts) nfonts = *maxfonts;
    // Phase 03.1 Plan 03 — Rule 1 bug fix:
    //   Do NOT write back to *maxfonts. The sole live callers
    //   (xvue/xvinit.f:84, xvue/xvpxfe.f:65) pass the Fortran PARAMETER
    //   constant MAXFONTS (=512 from incl/xvfontes.inc). gfortran -O may
    //   place PARAMETER-sourced actual arguments in read-only memory, so
    //   writing through the pointer crashes with SIGSEGV at the first
    //   xvtest1..4 run. The legacy C backend (xvue/xvuelc.c:612-1042)
    //   never writes *maxfonts — it is a pure input/capacity parameter.
    //   Matching that contract is both correct and safe.
    if (nbfonts)  *nbfonts  = nfonts;
    if (namefonts && nbchar) {
        for (int k = 0; k < nfonts; ++k) {
            std::strncpy(namefonts[k], kListFonts[k], 255);
            namefonts[k][255] = '\0';
            nbchar[k] = static_cast<int>(std::strlen(kListFonts[k]));
        }
    }

    if (visuclass) *visuclass = 4;   // TrueColor — Qt is always 24-bit

    // Palette ranges (legacy xvinfo semantics — ref xvuelc.c:612-1042):
    //   [0..15]   = imposed defaults
    //   [16..255] = user-modifiable
    if (n1coref)  *n1coref  = 0;
    if (ndcoref)  *ndcoref  = 15;
    if (n1coelf)  *n1coelf  = 0;
    if (ndcoelf)  *ndcoelf  = 15;
    if (n1coulf)  *n1coulf  = 16;
    if (ndcoulf)  *ndcoulf  = 255;
    if (nbcolo)   *nbcolo   = 256;
}

// ---- 11. xvrecuprgbdec_ (Phase 3 D-18, TEXT-05) ----
// Read-only snapshot of the process-lifetime palette into Fortran arrays.
// Pitfall 2: palette_init_once() runs inside XvueState ctor, which is
// inside XvueWindow ctor, which is called from xvinitgraphique_ — strictly
// before xvue/xvinit.f:143's XVRECUPRGBDEC call (the only live caller per
// the call-site audit in 03-RESEARCH.md).
void proc(xvrecuprgbdec)(int *nbcolor, float *r, float *g, float *b) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!nbcolor || !r || !g || !b) return;            // WR-03

    int n = *nbcolor;
    if (n < 0) n = 0;
    if (n > XvueState::kMaxPalette) n = XvueState::kMaxPalette;

    for (int i = 0; i < n; ++i) {
        r[i] = XvueState::red[i];
        g[i] = XvueState::green[i];
        b[i] = XvueState::blue[i];
    }
}

// ---- 12. xvactivervb_ (Phase 3 D-17 AMENDED per A3, TEXT-05) ----
// Bulk palette load. Research correction over the original D-17 text:
// the single live caller xvue/palcde.f:619 passes NDCOUL+1 cells, so this
// is a full user-palette refresh, NOT a transient one-shot color write.
// Matches legacy xvuelc.c:1072-1116 bulk-copy semantics.
void proc(xvactivervb)( int *palcour, int *nbcells,
                        float r[], float g[], float b[] ) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    (void)palcour;  // palcourc was a PostScript state var, dropped (02/D-26)
    if (!nbcells || !r || !g || !b) return;            // WR-03

    int n = *nbcells;
    if (n < 0) n = 0;
    if (n > XvueState::kMaxPalette) n = XvueState::kMaxPalette;

    for (int i = 0; i < n; ++i) {
        XvueState::red[i]   = r[i];
        XvueState::green[i] = g[i];
        XvueState::blue[i]  = b[i];
        XvueState::palette_cache_dirty_[i] = true;
    }
    // No painter touch, no flush — xvcouleur_ picks up new values via the
    // dirty-flag rebuild on its next call.
}

// IN-02: shared palette-index clamp + dirty-cache rebuild used by
// xvcouleur_, xvfond_, and apply_palette_foreground. Returns the resolved
// QColor; callers decide whether to write foreground_/background_ and
// whether to call applyPen(). File-local (static) by design — all call
// sites are in this translation unit.
static const QColor& palette_resolve(int icolor, int fallback) {
    int i = icolor;
    if (i < 0 || i >= XvueState::kMaxPalette) i = fallback;
    if (XvueState::palette_cache_dirty_[i]) {
        XvueState::palette_cache_[i] = QColor::fromRgbF(
            XvueState::red[i],
            XvueState::green[i],
            XvueState::blue[i]);
        XvueState::palette_cache_dirty_[i] = false;
    }
    return XvueState::palette_cache_[i];
}

// ---- 13. xvcouleur_ (Phase 3 D-14, TEXT-04) ----
// State-change entry: install palette_cache_[i] as current foreground.
// NO flush — drawing primitives run their own D-01 epilogue.
void proc(xvcouleur)(int *icolor) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!icolor) return;                                // WR-03

    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st) return;

    // IN-02: legacy "fallback to red" (index 1) via shared helper.
    st->foreground_ = palette_resolve(*icolor, 1);
    st->applyPen();  // 02/D-20 syncs brush + pen from foreground_
    // D-14: no update(), no processEvents — state-change only.
}

// ---- 14. xvpostscript_ ----
void proc(xvpostscript)(int *lasops) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvpostscript_");
    (void)lasops;
}

// ---- 15. fenetremempx_ ----
void proc(fenetremempx)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "fenetremempx_");
}

// ---- 16. mempxfenetre_ ----
void proc(mempxfenetre)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "mempxfenetre_");
}

// ---- 17. sauvefenetre_ ----
void proc(sauvefenetre)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "sauvefenetre_");
}

// ---- 18. restaurefenetre_ ----
void proc(restaurefenetre)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "restaurefenetre_");
}

// ---- 19. sauvemempx_ ----
void proc(sauvemempx)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "sauvemempx_");
}

// ---- 20. restauremempx_ ----
void proc(restauremempx)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "restauremempx_");
}

// ---- 21. effacemempx_ ----
void proc(effacemempx)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "effacemempx_");
}

// ---- 22. effacer_ (D-15) ----
void proc(effacer)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (st && st->painter_ && st->painter_->isActive() && st->backing_) {
        st->painter_->fillRect(st->backing_->rect(), st->background_);
    }
    if (win->canvas()) {
        win->canvas()->update();
    }
    // WR-02: deferred flush -- xvvoir_/xvpause_ pump the event loop.
}

// ---- 23. xvfond_ (Phase 3: retired Phase 2 two-entry hack) ----
// Phase 2 used a hardcoded {black, white} 2-entry palette (D-14). Phase 3
// routes the color source through the real palette_cache_ the same way
// xvcouleur_ does, while preserving the Phase 2 D-24 fillRect mechanism.
void proc(xvfond)(int *icolor) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!icolor) return;                               // WR-03

    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st) return;

    // IN-02: xvfond_ defaults to palette index 0 (black) via shared helper.
    // D-15: update XvueState::background_ through the live window, repaint.
    // Phase 2 D-24: re-fill the backing with the new background and flush.
    st->background_ = palette_resolve(*icolor, 0);
    if (st->painter_ && st->painter_->isActive() && st->backing_) {
        st->painter_->fillRect(st->backing_->rect(), st->background_);
    }
    if (win->canvas()) {
        win->canvas()->update();
    }
    // WR-02: deferred flush -- xvvoir_/xvpause_ pump the event loop.
}

// ---- 24. xvchargefonte_ (Phase 3 D-04, TEXT-01) ----
// Bundled DejaVu Sans Mono at XvueState::kFontSizes[*nofont]. QFont is RAII
// (Pitfall 6) so *nofont0 (the legacy "previous font to free") is ignored.
void proc(xvchargefonte)(int *nofont0, int *nofont, int *largpx, int *hautpx) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    (void)nofont0;  // Pitfall 6: QFont is RAII, no explicit free needed

    if (!nofont || !largpx || !hautpx) return;       // WR-03

    int idx = *nofont;
    if (idx < 0) idx = 0;
    if (idx >= XvueState::kNbFonts) idx = XvueState::kNbFonts - 1;

    auto& win = XvueApp::window_slot();
    if (!win) { *largpx = 0; *hautpx = 0; return; }
    auto* st = win->state();
    if (!st) { *largpx = 0; *hautpx = 0; return; }

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

// ---- 25. xvnbpixeltexte_ (Phase 3 D-05, TEXT-02) ----
void proc(xvnbpixeltexte)(char *texte, int *length, int *nbpxla, int *nbpxha) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!texte || !length || !nbpxla || !nbpxha) return;   // WR-03

    auto& win = XvueApp::window_slot();
    if (!win) { *nbpxla = 0; *nbpxha = 0; return; }
    auto* st = win->state();
    if (!st) { *nbpxla = 0; *nbpxha = 0; return; }

    // Pitfall 4: TWO-ARG fromLatin1 — explicit length from Fortran.
    QString qstr = QString::fromLatin1(texte, *length);
    *nbpxla = st->current_metrics_.horizontalAdvance(qstr);
    *nbpxha = st->current_metrics_.height();
}

// ---- 26. xvfermer_ ----
void proc(xvfermer)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    // Headless-test capture hook — matches the legacy xvuelc.c contract
    // (cf. xvue/xvuelc.c xvfermer_). Three optional env vars:
    //   MEFISTO_XVFERMER_READY_FILE  — path to touch right before
    //       destroying the window so an external screenshot tool can
    //       grab the final rendered state (root-window capture path).
    //   MEFISTO_XVFERMER_HOLD_MS     — milliseconds to sleep after
    //       flushing + touching the sentinel, before the window
    //       is destroyed. Default 0 (no hold).
    //   MEFISTO_QT_CAPTURE_PATH      — if set, call QWidget::grab()
    //       on the canvas and save directly to this PNG path before
    //       destroying the window. Works under QT_QPA_PLATFORM=offscreen
    //       with NO X11 dependency, so it does not need xcb-cursor0 or
    //       an Xvfb server — ideal for headless CI and for cases where
    //       the xcb platform plugin fails to initialise.
    // All three default to "no-op" when unset; interactive behavior is
    // unchanged. Qt-side counterpart of the legacy hook (commit e029b84)
    // with the extra in-process grab() path for full headless support.
    {
        auto& win = XvueApp::window_slot();
        if (win && win->canvas()) {
            win->canvas()->update();
        }
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);

        const char* qt_capture_path = std::getenv("MEFISTO_QT_CAPTURE_PATH");
        if (qt_capture_path && qt_capture_path[0] != '\0') {
            // In-process screenshot of the XvueState backing pixmap —
            // this is the authoritative surface every Phase 2 drawing
            // primitive paints on (see xvue_qt_canvas.cpp:35-36:
            // "paintEvent = drawPixmap blit ONLY" from state_->backing_).
            // Capturing it directly bypasses the widget paint pipeline,
            // which under QT_QPA_PLATFORM=offscreen does not always
            // trigger a full repaint between update() and grab(), and
            // guarantees the captured PNG matches exactly what the
            // drawing primitives produced up to this point. Also ends
            // any active painter session first so the pixmap is in a
            // consistent readable state.
            XvueState* st = (win ? win->state() : nullptr);
            if (st && st->backing_) {
                if (st->painter_ && st->painter_->isActive()) {
                    st->painter_->end();
                }
                if (!st->backing_->save(
                        QString::fromLatin1(qt_capture_path), "PNG")) {
                    std::fprintf(stderr,
                        "xvue-qt: failed to save backing to "
                        "MEFISTO_QT_CAPTURE_PATH=%s\n",
                        qt_capture_path);
                }
            } else if (win && win->canvas()) {
                // Fallback to widget grab if no backing is available
                // (e.g., the driver never called xvinfo_ to allocate it).
                QPixmap shot = win->canvas()->grab();
                if (!shot.save(
                        QString::fromLatin1(qt_capture_path), "PNG")) {
                    std::fprintf(stderr,
                        "xvue-qt: failed to save widget grab to "
                        "MEFISTO_QT_CAPTURE_PATH=%s\n",
                        qt_capture_path);
                }
            }
        }

        const char* ready_path = std::getenv("MEFISTO_XVFERMER_READY_FILE");
        const char* hold_env   = std::getenv("MEFISTO_XVFERMER_HOLD_MS");
        int hold_ms = 0;
        if (hold_env && hold_env[0] != '\0') {
            int hv = std::atoi(hold_env);
            if (hv >= 0 && hv <= 60000) hold_ms = hv;
        }
        if (ready_path && ready_path[0] != '\0') {
            FILE* f = std::fopen(ready_path, "w");
            if (f) { std::fputs("ready\n", f); std::fclose(f); }
        }
        if (hold_ms > 0) {
            // Pump the event loop during the hold so paint events
            // delivered to the canvas actually reach the X server.
            QElapsedTimer t; t.start();
            while (t.elapsed() < hold_ms) {
                QCoreApplication::processEvents(
                    QEventLoop::ExcludeUserInputEvents, 50);
                usleep(10 * 1000);
            }
        }
    }
    // D-06: destroy the window only. Do NOT touch qApp, do NOT call
    // qApp->quit(). The QApplication lives until the atexit handler.
    XvueApp::window_slot().reset();
    // D-06 addendum 2026-04-11 (debug phase-01-xvtest0-teardown-segfault):
    // drain DeferredDelete events queued by the widget-hierarchy teardown
    // immediately, while Qt is in a well-defined state. Otherwise these
    // events survive until the atexit handler and risk running against a
    // torn-down Qt. ExcludeUserInputEvents per Pitfall 4.
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

// ---- 27. xvpxfenetre_ (D-23) ----
void proc(xvpxfenetre)(int *x, int *y) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y) return;
    auto& win = XvueApp::window_slot();
    if (!win || !win->canvas()) { *x = 0; *y = 0; return; }
    *x = win->canvas()->width();   // logical pixels (SHELL-06)
    *y = win->canvas()->height();
}

// ---- 28. xvftexte_ (Phase 3 D-06, TEXT-03) ----
// Shares body with xvtexte_ under Phase 2 single-backing (02/D-05, Pitfall 7).
void proc(xvftexte)(char string[], int *length, int *x1, int *y1) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!string || !length || !x1 || !y1) return;      // WR-03
    xvue_qt_draw_text_common(string, *length, *x1, *y1);
}

// ---- 29. xvtexte_ (Phase 3 D-06, TEXT-03) ----
void proc(xvtexte)(char string[], int *length, int *x1, int *y1) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!string || !length || !x1 || !y1) return;      // WR-03
    xvue_qt_draw_text_common(string, *length, *x1, *y1);
}

// ---- 30. xvface_ (D-11, DRAW-03) ----
// Legacy: XFillPolygon -- fill only, no edge stroke (WR-01).
void proc(xvface)(int *n, MefistoPoint *pts) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;
    if (!n || *n < 3 || !pts) return;

    QPolygon poly;
    poly.reserve(*n);
    for (int i = 0; i < *n; ++i) {
        poly << QPoint(pts[i].x, pts[i].y);
    }
    // WR-01: drawPolygon strokes the outline with the current pen. Legacy
    // XFillPolygon only fills. Suppress the pen for the fill step so the
    // two colors stay independent once Phase 3 unlocks pen/brush divergence.
    QPen saved_pen = st->painter_->pen();
    st->painter_->setPen(Qt::NoPen);
    st->painter_->drawPolygon(poly, Qt::OddEvenFill);  // auto 1<->n close
    st->painter_->setPen(saved_pen);

    if (win->canvas()) win->canvas()->update();
}

// ---- 31. xvtypetrait_ (D-17, D-19, DRAW-06) ----
void proc(xvtypetrait)(int *ptype) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!ptype) return;   // WR-03
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st) return;
    st->pen_style_ = *ptype;
    st->applyPen();   // applyPen() self-gates on painter_->isActive()
}

// ---- 32. xvepaisseur_ (D-18, DRAW-06) ----
void proc(xvepaisseur)(int *pepais) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!pepais) return;   // WR-03
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st) return;
    st->pen_width_base_ = *pepais;
    st->applyPen();
}

// ---- 34. xvtrait_ (D-09, DRAW-02) ----
void proc(xvtrait)(int *x1, int *y1, int *x2, int *y2) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x1 || !y1 || !x2 || !y2) return;   // WR-03
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;
    st->painter_->drawLine(*x1, *y1, *x2, *y2);
    if (win->canvas()) win->canvas()->update();
    // WR-02: deferred flush -- xvvoir_/xvpause_ pump the event loop.
}

// ---- 33. xvftrait_ (D-09 -- semantically identical to xvtrait_ under
// the single-backing model; legacy window-vs-mempx split obsolete) ----
void proc(xvftrait)(int *x1, int *y1, int *x2, int *y2) {
    proc(xvtrait)(x1, y1, x2, y2);
}

// ---- 35. xvtraits_ (D-10, DRAW-02) ----
void proc(xvtraits)(int *nbpoints, MefistoPoint *points) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;
    if (!nbpoints || *nbpoints < 2 || !points) return;

    constexpr int STACK_LIMIT = 128;
    if (*nbpoints <= STACK_LIMIT) {
        QPoint qpts[STACK_LIMIT];
        for (int i = 0; i < *nbpoints; ++i) {
            qpts[i] = QPoint(points[i].x, points[i].y);
        }
        st->painter_->drawPolyline(qpts, *nbpoints);
    } else {
        std::vector<QPoint> qpts;
        qpts.reserve(*nbpoints);
        for (int i = 0; i < *nbpoints; ++i) {
            qpts.emplace_back(points[i].x, points[i].y);
        }
        st->painter_->drawPolyline(qpts.data(),
                                   static_cast<int>(qpts.size()));
    }
    if (win->canvas()) win->canvas()->update();
    // WR-02: deferred flush -- xvvoir_/xvpause_ pump the event loop.
}

// Internal helper — inline palette install for callers that already hold
// a resolved XvueState*. Mirrors proc(xvcouleur) body but skips
// XvueApp::ensure()/thread-assert/window-slot re-lookup. IN-02: delegates
// clamp + dirty-cache rebuild to the shared palette_resolve() helper.
static void apply_palette_foreground(XvueState* st, int icolor) {
    st->foreground_ = palette_resolve(icolor, 1);  // legacy fallback to red
    st->applyPen();  // rebuilds pen from pen_style_ + foreground_, brush from foreground_
}

// ---- 36. xvfacetraits_ (D-12, DRAW-03; 02.1/D-1: honor ncf/nca) ----
void proc(xvfacetraits)(int *ncf, int *nca, int *n, MefistoPoint *pts) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!ncf || !nca || !n || *n < 3 || !pts) return;   // WR-03
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;

    QPolygon poly;
    poly.reserve(*n);
    for (int i = 0; i < *n; ++i) {
        poly << QPoint(pts[i].x, pts[i].y);
    }

    // 02.1/D-1: mirror xvue/xvuelc.c:2055,2064 — install ncf, fill; install
    // nca, outline. Dashed style set by xvtypetrait_ survives the outline
    // step because pen_style_ is threaded through XvueState::applyPen()
    // (see xvue_qt_state.cpp:136-142), which rebuilds pen_ from
    // pen_style_ + foreground_ on every apply_palette_foreground() call.
    // WR-01: fill must not stroke (setPen(Qt::NoPen)); outline must not
    // re-fill (setBrush(Qt::NoBrush)).
    apply_palette_foreground(st, *ncf);
    st->painter_->setPen(Qt::NoPen);
    st->painter_->drawPolygon(poly, Qt::OddEvenFill);   // fill (D-12 order)
    // IN-01: no pen restore between fill and outline — the next
    // apply_palette_foreground(st, *nca) call invokes applyPen() which
    // unconditionally rebuilds pen_ from pen_style_ + foreground_, so any
    // intermediate setPen(Qt::NoPen) is immediately superseded.

    apply_palette_foreground(st, *nca);                 // rebuilds pen with nca color, preserves dash style
    st->painter_->setBrush(Qt::NoBrush);
    st->painter_->drawPolygon(poly);                    // outline (D-12 order) — uses current pen (incl dash)
    // IN-03 / WR-04: legacy leaves foreground_ at nca (we match that by
    // not touching foreground_). The outline step above left the painter
    // brush at Qt::NoBrush; while any subsequent xvcouleur_/applyPen()
    // would refresh it, WR-04 tightens the post-condition by restoring
    // the brush here so primitives that read the live painter brush
    // directly (e.g. xvface_) observe the canonical post-ncf->nca final
    // state without depending on an intervening state-change call.
    st->painter_->setBrush(QBrush(st->foreground_, Qt::SolidPattern));

    if (win->canvas()) win->canvas()->update();
}

// ---- 37. xvsouris_ ----
// Headless-test short-circuit: when MEFISTO_XVSOURIS_AUTOEXIT is set,
// return a synthetic SPACE keypress so driver loops of the form
//   IF (NOTYEV .NE. 2) GOTO 100
// exit on first iteration. Pairs with the MEFISTO_XVFERMER_* capture
// hook above and matches the legacy xvuelc.c contract.
//
// When the env var is not set, the Qt-side behavior stays the stub it
// was before (returns no event — the driver's input loop is expected
// to be handled by the Qt event pump). Interactive paths are never
// regressed.
void proc(xvsouris)(int *notypeevent, int *nbc, int *x1, int *y1) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    const char* autoexit = std::getenv("MEFISTO_XVSOURIS_AUTOEXIT");
    if (autoexit && autoexit[0] != '\0') {
        int delay_ms = 500;
        const char* delay_env = std::getenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS");
        if (delay_env && delay_env[0] != '\0') {
            int d = std::atoi(delay_env);
            if (d >= 0 && d <= 60000) delay_ms = d;
        }
        // Paint the canvas and pump events so the window is up to date.
        auto& win = XvueApp::window_slot();
        if (win && win->canvas()) win->canvas()->update();
        QElapsedTimer t; t.start();
        while (t.elapsed() < delay_ms) {
            QCoreApplication::processEvents(
                QEventLoop::ExcludeUserInputEvents, 25);
            usleep(5 * 1000);
        }
        if (notypeevent) *notypeevent = 2;
        if (nbc)         *nbc         = ' ';
        if (x1)          *x1          = 0;
        if (y1)          *y1          = 0;
        return;
    }
    static bool warned = false;
    warn_once(warned, "xvsouris_");
    (void)notypeevent; (void)nbc; (void)x1; (void)y1;
}

// ---- 38. xvsouris2_ ----
void proc(xvsouris2)(int *items, int *pmin0, int *notypeevent, int *ibutton, int *x1, int *y1) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    const char* autoexit = std::getenv("MEFISTO_XVSOURIS_AUTOEXIT");
    if (autoexit && autoexit[0] != '\0') {
        int delay_ms = 500;
        const char* delay_env = std::getenv("MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS");
        if (delay_env && delay_env[0] != '\0') {
            int d = std::atoi(delay_env);
            if (d >= 0 && d <= 60000) delay_ms = d;
        }
        auto& win = XvueApp::window_slot();
        if (win && win->canvas()) win->canvas()->update();
        QElapsedTimer t; t.start();
        while (t.elapsed() < delay_ms) {
            QCoreApplication::processEvents(
                QEventLoop::ExcludeUserInputEvents, 25);
            usleep(5 * 1000);
        }
        if (notypeevent) *notypeevent = 2;
        if (ibutton)     *ibutton     = ' ';
        if (x1)          *x1          = 0;
        if (y1)          *y1          = 0;
        (void)items; (void)pmin0;
        return;
    }
    static bool warned = false;
    warn_once(warned, "xvsouris2_");
    (void)items; (void)pmin0; (void)notypeevent; (void)ibutton; (void)x1; (void)y1;
}

// ---- 39. deplsouris_ ----
void proc(deplsouris)(int *x, int *y) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "deplsouris_");
    (void)x; (void)y;
}

// ---- 40. xvvoir_ (D-02) ----
void proc(xvvoir)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    auto& win = XvueApp::window_slot();
    if (win && win->canvas()) {
        win->canvas()->update();
    }
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

// ---- 41. xvpause_ ----
void proc(xvpause)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvpause_");
}

// ---- 42. xvfbordrectangle_ (D-13, DRAW-04) ----
// Legacy: XDrawRectangle -- outline only.
void proc(xvfbordrectangle)(int *x, int *y, int *width, int *height) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y || !width || !height) return;   // WR-03
    xvue_qt_draw_rect_common(*x, *y, *width, *height, RectMode::Outline);
}

// ---- 43. xvbordrectangle_ (D-13, DRAW-04) ----
// Legacy: XDrawRectangle -- outline only.
void proc(xvbordrectangle)(int *x, int *y, int *width, int *height) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y || !width || !height) return;   // WR-03
    xvue_qt_draw_rect_common(*x, *y, *width, *height, RectMode::Outline);
}

// ---- 44. xvfrectangle_ (D-13, DRAW-04) ----
// Legacy: XFillRectangle -- fill only.
void proc(xvfrectangle)(int *x, int *y, int *width, int *height) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y || !width || !height) return;   // WR-03
    xvue_qt_draw_rect_common(*x, *y, *width, *height, RectMode::Fill);
}

// ---- 45. xvrectangle_ (D-13, DRAW-04) ----
// Legacy: XFillRectangle -- fill only.
void proc(xvrectangle)(int *x, int *y, int *width, int *height) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y || !width || !height) return;   // WR-03
    xvue_qt_draw_rect_common(*x, *y, *width, *height, RectMode::Fill);
}

// ---- 46. xvbordarcellipse_ (D-14, RESEARCH Q1: drawArc, DRAW-05) ----
// Legacy xvuelc.c:2554 -- uses XDrawArc, outline only, float* angles in degrees.
// X11 x64 (1/64 deg) -> Qt x16 (1/16 deg). NOT x64.
void proc(xvbordarcellipse)(int *x, int *y, int *width, int *height,
                            float *angle1, float *angle2) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y || !width || !height || !angle1 || !angle2) return;  // WR-03
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;

    const QRect bbox(*x - *width, *y - *height,
                     *width * 2,  *height * 2);
    const int start_16 = static_cast<int>(*angle1 * 16.0f);
    const int span_16  = static_cast<int>(*angle2 * 16.0f);

    st->painter_->drawArc(bbox, start_16, span_16);  // outline -- matches XDrawArc

    if (win->canvas()) win->canvas()->update();
    // WR-02: deferred flush -- xvvoir_/xvpause_ pump the event loop.
}

// ---- 47. xvarcellipse_ (D-14, RESEARCH Q1 CORRECTION: drawPie, DRAW-05) ----
// Legacy xvuelc.c:2616 -- uses XFillArc (filled pie slice to center).
// Qt's equivalent of XFillArc is drawPie, NOT drawArc. drawArc would only
// stroke the curve; XFillArc fills the pie wedge. See 02-RESEARCH.md Q1.
void proc(xvarcellipse)(int *x, int *y, int *width, int *height,
                        float *angle1, float *angle2) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!x || !y || !width || !height || !angle1 || !angle2) return;  // WR-03
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;

    const QRect bbox(*x - *width, *y - *height,
                     *width * 2,  *height * 2);
    const int start_16 = static_cast<int>(*angle1 * 16.0f);
    const int span_16  = static_cast<int>(*angle2 * 16.0f);

    st->painter_->drawPie(bbox, start_16, span_16);  // filled wedge -- matches XFillArc

    if (win->canvas()) win->canvas()->update();
    // WR-02: deferred flush -- xvvoir_/xvpause_ pump the event loop.
}

// ---- 48. tempscpu_ ----
void proc(tempscpu)(double *tclock) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "tempscpu_");
    (void)tclock;
}

// ---- 49. secondes1970_ ----
void proc(secondes1970)(double *secondes) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "secondes1970_");
    (void)secondes;
}

// ---- 50. secondes1969_ ----
void proc(secondes1969)(double *secondes) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "secondes1969_");
    (void)secondes;
}

// ---- 51. nomordinateurhote_ ----
void proc(nomordinateurhote)(char *host, int *nbcar) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "nomordinateurhote_");
    (void)host; (void)nbcar;
}

// ---- 52. ladate_ ----
void proc(ladate)(int *a, int *m, int *j) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "ladate_");
    (void)a; (void)m; (void)j;
}

// ---- 53. heureminuteseconde_ ----
void proc(heureminuteseconde)(int *h, int *m, int *s, int *millis) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "heureminuteseconde_");
    (void)h; (void)m; (void)s; (void)millis;
}

// ---- 54. valvarenv_ ----
void proc(valvarenv)( char *nom, int *lval_admis,
                      char *val, int *lval_trouve ) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "valvarenv_");
    (void)nom; (void)lval_admis; (void)val; (void)lval_trouve;
}

// ---- 55. xvinitierps_ ----
void proc(xvinitierps)(int *modeps) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvinitierps_");
    (void)modeps;
}

// ---- 56. xvimprimerps_ ----
void proc(xvimprimerps)(char nomfichier[], int *length) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvimprimerps_");
    (void)nomfichier; (void)length;
}

// ---- 57. xvsauverps_ ----
void proc(xvsauverps)(char nomfichier[], int *length) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvsauverps_");
    (void)nomfichier; (void)length;
}

} // extern "C"
