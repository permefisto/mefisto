// xvue/qt/tests/test_helpers.h
// Phase 5 Wave 0 (05-01 Task 1): test helpers for driving nested QEventLoop
// from QTest bodies. Wave 0 only provides the pump helper and placeholder
// canvas stubs; Plan 02+ will wire through the full Fortran ABI.
#pragma once
#include <QPoint>

class QWidget;
class XvueCanvas;

namespace xvue_test {

// Spins the Qt event loop for a short duration so that posted events reach
// the bridge before the test inspects its state.
void pumpEvents(int ms = 10);

// Placeholder for Plan 02+: will open an offscreen XvueCanvas with the
// application state installed and return its pointer. Wave 0 returns nullptr.
XvueCanvas* createTestCanvas();
void        destroyTestCanvas();

} // namespace xvue_test
