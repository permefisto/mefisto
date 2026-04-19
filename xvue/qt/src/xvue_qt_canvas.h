// xvue/qt/src/xvue_qt_canvas.h
// Phase 1 (D-05) + Phase 2 (D-04..D-08): central widget. paintEvent blits
// the persistent backing pixmap owned by XvueState; resizeEvent reallocates
// that backing, preserves top-left content (DRAW-09), ends+restarts the
// long-lived QPainter, and re-applies Antialiasing + pen + brush.
// Phase 6.0 Plan 05 (UX-12, UX-13, UI-SPEC §Canvas, Flag #3): adds Qt-native
// gestures (wheel zoom, middle-drag pan, right-click context menu),
// mouseCoords signal (Plan 03 emits, Plan 06 connects), resetView() slot
// (View → Fit), and an empty-state hint rendered when state_->has_user_content_
// is false. View transform is composited at paintEvent time only — backing
// pixmap and Fortran draw calls stay in canonical pixel coordinates ("Fortran
// must not notice" invariant at the UX layer).
#pragma once
#include <QWidget>
#include <QPoint>
#include <QTransform>

struct XvueState;
class QResizeEvent;
class QWheelEvent;
class QMouseEvent;
class QContextMenuEvent;
class XvueEventBridge;   // Phase 6.0 Plan 03 friendship (forward decl)

class XvueCanvas : public QWidget {
    Q_OBJECT  // Pitfall 8
    // Phase 6.0 Plan 03 friendship: allow XvueEventBridge to emit mouseCoords
    // from inside its eventFilter MouseMove branch without Qt 6.4+ "calling
    // signal from outside class" warnings.
    friend class XvueEventBridge;
public:
    explicit XvueCanvas(XvueState* state, QWidget* parent = nullptr);
    ~XvueCanvas() override;

    // Phase 6.0 Plan 05 (UX-12): reset viewport to identity transform + repaint.
    // Wired to View → Fit to Window (Ctrl+0) in Plan 06.
    void resetView();

    // Phase 6.0 Plan 05 — test-only observability hooks. Not part of the
    // production contract; documented as such so downstream code does not
    // become accidental consumers. Used by xvue_qt_canvas_gestures_tests
    // to assert behavior that is otherwise hard to verify (menu-shown
    // counter, empty-state branch entered).
    int  contextMenuShownCount_ = 0;
    bool lastPaintDrewEmptyState_ = false;

signals:
    // Phase 6.0 Plan 05 (UX-12): emitted on every MouseMove (no-button tracking
    // enabled in ctor since Phase 5). Plan 06 connects to status bar coord label.
    // Plan 03 will emit this from XvueEventBridge::eventFilter (the friend
    // declaration above permits the cross-class emit).
    void mouseCoords(QPoint posLogicalPx);

protected:
    void paintEvent(QPaintEvent* event) override;
    void resizeEvent(QResizeEvent* event) override;   // Phase 2 (D-04, D-07)

    // Phase 6.0 Plan 05 additions (UX-12, UX-13, UI-SPEC §Canvas interaction):
    void wheelEvent(QWheelEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;
    void contextMenuEvent(QContextMenuEvent* event) override;

private:
    XvueState* state_ = nullptr;  // raw pointer — owned by XvueWindow (D-15)

    // Phase 6.0 Plan 05 — middle-button pan tracking. pan_anchor_transform_
    // captures the view_transform_ at press time so a continuous drag
    // produces a single linear translation rather than accumulating drift.
    QPoint     pan_start_            { };
    QTransform pan_anchor_transform_ { };
    bool       pan_active_           = false;
};
