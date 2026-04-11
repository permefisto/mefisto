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
#include "xvue_qt_api.h"
#include "xvue_qt_app.h"
#include "xvue_qt_window.h"
#include "xvue_qt_canvas.h"
#include "xvue_qt_state.h"
#include <QApplication>
#include <QCoreApplication>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QGuiApplication>
#include <QPainter>
#include <QPixmap>
#include <QPoint>
#include <QPolygon>
#include <QRect>
#include <QScreen>
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
void proc(xtinit)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xtinit_");
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

    // D-03: Phase 1 partial. Resize the window if it exists; everything else
    // (palette ranges, font tables, visual class) is zeroed out and the
    // warn-once line is kept verbatim from the Phase 0 pattern. Phase 3 will
    // replace the palette/font branch with real colormap plumbing.
    auto& win = XvueApp::window_slot();
    if (win && ix && iy) {
        win->resize(*ix, *iy);
    }

    // Zero palette/font outputs so Fortran callers see deterministic values.
    if (maxfonts)  *maxfonts  = 0;
    if (n1coref)   *n1coref   = 0;
    if (ndcoref)   *ndcoref   = 0;
    if (n1coelf)   *n1coelf   = 0;
    if (ndcoelf)   *ndcoelf   = 0;
    if (n1coulf)   *n1coulf   = 0;
    if (ndcoulf)   *ndcoulf   = 0;
    if (nbcolo)    *nbcolo    = 0;
    if (nbfonts)   *nbfonts   = 0;
    if (visuclass) *visuclass = 0;
    (void)namefonts;
    (void)nbchar;

    static bool warned_xvinfo_partial = false;
    if (!warned_xvinfo_partial) {
        std::fprintf(stderr,
            "xvue-qt: stub xvinfo_ palette outputs not implemented yet\n");
        warned_xvinfo_partial = true;
    }
}

// ---- 11. xvrecuprgbdec_ ----
void proc(xvrecuprgbdec)(int *nbcolor, float *r, float *g, float *b) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvrecuprgbdec_");
    (void)nbcolor; (void)r; (void)g; (void)b;
}

// ---- 12. xvactivervb_ ----
void proc(xvactivervb)( int *palcour, int *nbcells,
                        float r[], float g[], float b[] ) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvactivervb_");
    (void)palcour; (void)nbcells; (void)r; (void)g; (void)b;
}

// ---- 13. xvcouleur_ ----
void proc(xvcouleur)(int *icolor) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvcouleur_");
    (void)icolor;
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
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

// ---- 23. xvfond_ ----
void proc(xvfond)(int *icolor) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    if (!icolor) return;

    // D-14: Phase 1 has no palette yet. Minimal 2-entry mapping matches the
    // legacy X11 BlackPixel/WhitePixel convention (xvuelc.c:935 et seq.).
    QColor chosen = Qt::black;
    if (*icolor == 0) {
        chosen = Qt::black;
    } else if (*icolor == 1) {
        chosen = Qt::white;
    } else {
        static bool warned_xvfond_range = false;
        if (!warned_xvfond_range) {
            std::fprintf(stderr,
                "xvue-qt: xvfond_ palette index %d out of Phase 1 range "
                "(Phase 3 will add full colormap)\n", *icolor);
            warned_xvfond_range = true;
        }
        chosen = Qt::black;
    }

    // D-15: update XvueState::background_ through the live window, schedule
    // repaint. With no open window, xvfond_ is a no-op past the warn-once line
    // (Phase 1 XvueState is owned by XvueWindow).
    // Phase 2 (D-24): re-fill the backing with the new background and flush.
    auto& win = XvueApp::window_slot();
    if (win) {
        auto* st = win->state();
        if (st) {
            st->background_ = chosen;
            if (st->painter_ && st->painter_->isActive() && st->backing_) {
                st->painter_->fillRect(st->backing_->rect(), st->background_);
            }
        }
        if (win->canvas()) {
            win->canvas()->update();
        }
    }
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

// ---- 24. xvchargefonte_ ----
void proc(xvchargefonte)(int *nofont0, int *nofont, int *largpx, int *hautpx) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvchargefonte_");
    (void)nofont0; (void)nofont; (void)largpx; (void)hautpx;
}

// ---- 25. xvnbpixeltexte_ ----
void proc(xvnbpixeltexte)(char *texte, int *length, int *nbpxla, int *nbpxha) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvnbpixeltexte_");
    (void)texte; (void)length; (void)nbpxla; (void)nbpxha;
}

// ---- 26. xvfermer_ ----
void proc(xvfermer)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
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

// ---- 28. xvftexte_ ----
void proc(xvftexte)(char string[], int *length, int *x1, int *y1) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvftexte_");
    (void)string; (void)length; (void)x1; (void)y1;
}

// ---- 29. xvtexte_ ----
void proc(xvtexte)(char string[], int *length, int *x1, int *y1) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvtexte_");
    (void)string; (void)length; (void)x1; (void)y1;
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
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;
    st->painter_->drawLine(*x1, *y1, *x2, *y2);
    if (win->canvas()) win->canvas()->update();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
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
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

// ---- 36. xvfacetraits_ (D-12, DRAW-03) ----
void proc(xvfacetraits)(int *ncf, int *nca, int *n, MefistoPoint *pts) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    (void)ncf; (void)nca;  // TODO(phase 3): honor fill/edge color indices
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
    // WR-01: split fill and outline cleanly. Fill must not stroke the edge
    // (setPen(Qt::NoPen)); outline must not re-fill (setBrush(Qt::NoBrush)).
    // Without this, once Phase 3 differentiates ncf/nca colors, the fill
    // step would paint the edge in the fill color and the outline step
    // would repaint with the edge color, wasting work and producing a
    // thicker edge line than legacy.
    QPen   saved_pen   = st->painter_->pen();
    QBrush saved_brush = st->painter_->brush();
    st->painter_->setPen(Qt::NoPen);
    st->painter_->drawPolygon(poly, Qt::OddEvenFill);  // fill (D-12 order)
    st->painter_->setPen(saved_pen);
    st->painter_->setBrush(Qt::NoBrush);
    st->painter_->drawPolygon(poly);                   // outline (D-12 order)
    st->painter_->setBrush(saved_brush);

    if (win->canvas()) win->canvas()->update();
}

// ---- 37. xvsouris_ ----
void proc(xvsouris)(int *notypeevent, int *nbc, int *x1, int *y1) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    static bool warned = false;
    warn_once(warned, "xvsouris_");
    (void)notypeevent; (void)nbc; (void)x1; (void)y1;
}

// ---- 38. xvsouris2_ ----
void proc(xvsouris2)(int *items, int *pmin0, int *notypeevent, int *ibutton, int *x1, int *y1) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
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
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

// ---- 47. xvarcellipse_ (D-14, RESEARCH Q1 CORRECTION: drawPie, DRAW-05) ----
// Legacy xvuelc.c:2616 -- uses XFillArc (filled pie slice to center).
// Qt's equivalent of XFillArc is drawPie, NOT drawArc. drawArc would only
// stroke the curve; XFillArc fills the pie wedge. See 02-RESEARCH.md Q1.
void proc(xvarcellipse)(int *x, int *y, int *width, int *height,
                        float *angle1, float *angle2) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
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
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
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
