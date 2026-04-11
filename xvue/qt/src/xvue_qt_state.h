// xvue/qt/src/xvue_qt_state.h
// Phase 1 (D-04): runtime state holder. In Phase 1 this struct has exactly
// one field: the background colour used by XvueCanvas::paintEvent. Phase 2
// will add pen/brush/line width and a long-lived QPainter*; Phase 3 will add
// the indexed palette; Phase 4 will add the XvuePixmapStack. Fields are added
// additively — never rename or reorder existing fields.
#pragma once
#include <QColor>
#include <Qt>

struct XvueState {
    // D-04: default matches legacy X11 attributes.background_pixel = BlackPixel
    // (xvue/xvuelc.c:935). Keeping the same default prevents visible drift in
    // the Phase 8 A/B validation.
    QColor background_ = Qt::black;
};
