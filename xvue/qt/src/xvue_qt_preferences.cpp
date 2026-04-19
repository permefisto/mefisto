// xvue/qt/src/xvue_qt_preferences.cpp
// Phase 6.0 Plan 01 (scaffold): empty dialog body. Plan 04 fills the form.
#include "xvue_qt_preferences.h"

XvuePreferencesDialog::XvuePreferencesDialog(QWidget* parent)
    : QDialog(parent)
{
    // Plan 04 will set windowTitle to xvueT(MsgId::PreferencesTitle) and
    // build the form (color-scheme combo, console-visible check, recent-count
    // spin). Plan 01 leaves it empty so the symbol exists for tests.
    setWindowTitle(QString());
}

XvuePreferencesDialog::~XvuePreferencesDialog() = default;

void XvuePreferencesDialog::onAccept() {
    // Plan 04 fills: write XvuePrefs::saveColorScheme/saveConsoleVisible/etc.
    accept();
}
