// xvue/qt/src/xvue_qt_state.cpp
// Phase 2 (D-16, D-17, D-20): applyPen() single source of truth for pen+brush;
// destructor enforces the single-long-lived-painter invariant from DRAW-01.
#include "xvue_qt_state.h"
#include <QPainter>
#include <QPixmap>
#include <algorithm>   // std::max

float   XvueState::red[XvueState::kMaxPalette]             = {};
float   XvueState::green[XvueState::kMaxPalette]           = {};
float   XvueState::blue[XvueState::kMaxPalette]            = {};
QColor  XvueState::palette_cache_[XvueState::kMaxPalette];
bool    XvueState::palette_cache_dirty_[XvueState::kMaxPalette] = {};
bool    XvueState::palette_initialized_                    = false;

// Verbatim lift of xvue/xvuelc.c:378-461 (xvCouleursImposees) with n1coel=0.
static void imposed_defaults_fill()
{
    /* noir    */
    XvueState::red  [0] = 0.0f;
    XvueState::green[0] = 0.0f;
    XvueState::blue [0] = 0.0f;

    /* rouge   */
    XvueState::red  [1] = 1.0f;
    XvueState::green[1] = 0.0f;
    XvueState::blue [1] = 0.0f;

    /* vert sombre */
    XvueState::red  [2] =  50.f/256.f;
    XvueState::green[2] = 200.f/256.f;
    XvueState::blue [2] =  50.f/256.f;

    /* bleu    */
    XvueState::red  [3] = 0.0f;
    XvueState::green[3] = 0.0f;
    XvueState::blue [3] = 1.0f;

    /* cyan    */
    XvueState::red  [4] = 0.0f;
    XvueState::green[4] = 0.8f;
    XvueState::blue [4] = 1.0f;

    /* jaune   */
    XvueState::red  [5] = 1.0f;
    XvueState::green[5] = 1.0f;
    XvueState::blue [5] = 0.0f;

    /* magenta */
    XvueState::red  [6] = 1.0f;
    XvueState::green[6] = 0.0f;
    XvueState::blue [6] = 1.0f;

    /* blanc   */
    XvueState::red  [7] = 1.0f;
    XvueState::green[7] = 1.0f;
    XvueState::blue [7] = 1.0f;

    /* gris1 sombre bleute */
    XvueState::red  [8] =  80.f/256.f;
    XvueState::green[8] =  80.f/256.f;
    XvueState::blue [8] = 100.f/256.f;

    /* gris2 moyen bleute */
    XvueState::red  [9] = 150.f/256.f;
    XvueState::green[9] = 150.f/256.f;
    XvueState::blue [9] = 178.f/256.f;

    /* gris3 clair bleute */
    XvueState::red  [10] = 220.f/256.f;
    XvueState::green[10] = 220.f/256.f;
    XvueState::blue [10] = 256.f/256.f;

    /* peachpuff  remplace  beige */
    XvueState::red  [11] = 256.f/256.f;
    XvueState::green[11] = 218.f/256.f;
    XvueState::blue [11] = 185.f/256.f;

    /* orange   */
    XvueState::red  [12] = 1.0f;
    XvueState::green[12] = 0.5f;
    XvueState::blue [12] = 0.0f;

    /* saumon */
    XvueState::red  [13] = 250.f/256.f;
    XvueState::green[13] = 128.f/256.f;
    XvueState::blue [13] = 114.f/256.f;

    /* rose */
    XvueState::red  [14] = 1.0f;
    XvueState::green[14] = 190.f/256.f;
    XvueState::blue [14] = 206.f/256.f;

    /* turquoise */
    XvueState::red  [15] =  74.f/256.f;
    XvueState::green[15] = 250.f/256.f;
    XvueState::blue [15] = 160.f/256.f;
}

static void palette_init_once()
{
    if (XvueState::palette_initialized_) return;

    imposed_defaults_fill();

    for (int i = 0; i < XvueState::kMaxPalette; ++i) {
        XvueState::palette_cache_[i] = QColor::fromRgbF(
            XvueState::red[i],
            XvueState::green[i],
            XvueState::blue[i]);
        XvueState::palette_cache_dirty_[i] = false;
    }

    XvueState::palette_initialized_ = true;
}

XvueState::XvueState() {
    palette_init_once();
}

void XvueState::applyPen() {
    Qt::PenStyle qt_style = Qt::SolidLine;
    int          width    = pen_width_base_;
    switch (pen_style_) {
        case 0:
            qt_style = Qt::SolidLine;
            break;
        case 1:
            qt_style = Qt::DashLine;
            break;
        default:  // type 2 or anything else -> dash, double width (D-17)
            qt_style = Qt::DashLine;
            width    = std::max(1, pen_width_base_ * 2);
            break;
    }
    pen_.setColor(foreground_);
    pen_.setWidth(width);
    pen_.setStyle(qt_style);
    brush_ = QBrush(foreground_, Qt::SolidPattern);
    if (painter_ && painter_->isActive()) {
        painter_->setPen(pen_);
        painter_->setBrush(brush_);
    }
}

XvueState::~XvueState() {
    // Phase 4 D-03: delete saved_canvas_ BEFORE painter_ and backing_.
    // It is never the active QPainter target so ordering is safe, but
    // matching the "least-entangled first" pattern keeps future audits honest.
    delete saved_canvas_;
    saved_canvas_ = nullptr;

    // Phase 5 (D-08): mirror the saved_canvas_ pattern for the two
    // accrochage pixmaps. Both are orthogonal to painter_/backing_ so the
    // order relative to those two does not matter.
    delete mempxaccro_;
    mempxaccro_ = nullptr;
    delete accroche_undo_tile_;
    accroche_undo_tile_ = nullptr;

    if (painter_) {
        if (painter_->isActive()) {
            painter_->end();
        }
        delete painter_;
        painter_ = nullptr;
    }
    delete backing_;
    backing_ = nullptr;
}
