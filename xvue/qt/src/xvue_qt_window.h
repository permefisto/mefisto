// xvue/qt/src/xvue_qt_window.h
// Phase 1 (D-02, D-04, D-15, Pitfall 6, Pitfall 8): bare QMainWindow.
// No menu bar, no toolbar, no status bar, no dock widgets — those arrive in
// Phase 6 (ROADMAP.md). XvueWindow owns the XvueState and the XvueCanvas.
// Phase 5 (D-02, EVENT-01): XvueWindow also owns the XvueEventBridge that
// is installed as event filter on the canvas; lifetime matches the window.
// Phase 6.0 Plan 03 (UX-03): adds menuBridge() accessor + setMenuBridgeForTesting()
// helper. The accessor returns nullptr in Plan 03-era builds (no production
// menu bridge wired yet). Plan 06 ctors a real XvueMenuBridge in the window
// constructor and the accessor returns it. The setter is a TEST-ONLY hook
// that lets Plan 03's pre-drain unit tests inject a bridge before Plan 06
// lands. Production code never calls the setter.
#pragma once
#include <QMainWindow>
#include "xvue_qt_state.h"

class XvueCanvas;
class XvueEventBridge;
class XvueMenuBridge;   // Phase 6.0 Plan 03 fwd decl (avoid header cycle)

class XvueWindow : public QMainWindow {
    Q_OBJECT  // Pitfall 8: required even with no signals/slots in Phase 1.
public:
    explicit XvueWindow(QWidget* parent = nullptr);
    ~XvueWindow() override;

    XvueState*       state()  { return &state_; }
    XvueCanvas*      canvas() { return canvas_; }
    XvueEventBridge* bridge() { return bridge_; }  // non-null after construction

    // Phase 6.0 Plan 03 (UX-03): menu bridge accessor. Returns nullptr in
    // Plan 03-era builds; Plan 06 wires the real XvueMenuBridge in the ctor
    // and this accessor will then return a non-null pointer. Callers MUST
    // null-check (the pre-drain block in xvue_qt_event.cpp::waitForEvent
    // and xvue_qt_api.cpp::xvsouris_ both do).
    XvueMenuBridge*  menuBridge() { return menuBridge_; }

    // Phase 6.0 Plan 03 (test-only): inject a bridge BEFORE Plan 06 wires
    // the production one. Ownership: caller retains; we do NOT reparent.
    // Plan 06 ctor replaces this wiring by parenting a fresh bridge to the
    // window. Production code never invokes this setter — it exists solely
    // so Plan 03's pre-drain QTest cases (testPreDrainSingleChar et al.) can
    // exercise the new code path before the production wiring lands.
    void setMenuBridgeForTesting(XvueMenuBridge* mb) { menuBridge_ = mb; }

private:
    XvueState        state_{};    // D-04: single-field state struct
    XvueCanvas*      canvas_     = nullptr;  // Qt-owned via setCentralWidget
    XvueEventBridge* bridge_     = nullptr;  // Phase 5 (D-02): Qt-owned via parent=this
    XvueMenuBridge*  menuBridge_ = nullptr;  // Phase 6.0 Plan 03; Plan 06 assigns
};
