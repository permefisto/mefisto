// xvue/qt/src/xvue_qt_preferences.cpp
// Phase 6.0 Plan 04: single-panel Preferences dialog (UX D-05/D-06/D-07).
//
// Layout (UI-SPEC §Preferences):
//   QFormLayout with three rows:
//     1. Color scheme    : QComboBox (System / Light / Dark)
//     2. Console default : QCheckBox (visible at startup)
//     3. Recent count    : QSpinBox (1..20) — v1 placeholder, disabled
//   Bottom: QDialogButtonBox(Ok | Cancel)
//
// On Ok we round-trip into XvuePrefs (Plan 02 wires the QSettings backend).
// On Cancel we discard via QDialog::reject — no XvuePrefs writes.
//
// Plan-04 known limitation: the recent-count spin box is disabled because the
// 10-entry cap currently lives as a compile-time constant inside
// XvuePrefs::pushRecentProject. Wiring it through is deferred to a future
// plan; we ship the control today for UI-SPEC layout conformance.
#include "xvue_qt_preferences.h"

#include "xvue_qt_i18n.h"
#include "xvue_qt_prefs.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QPushButton>
#include <QSpinBox>
#include <QVBoxLayout>

XvuePreferencesDialog::XvuePreferencesDialog(QWidget* parent)
    : QDialog(parent)
{
    setWindowTitle(xvueT(MsgId::PreferencesTitle));

    // UI-SPEC §Spacing scale: inner margin = 8 px (sm).
    auto* outer = new QVBoxLayout(this);
    outer->setContentsMargins(8, 8, 8, 8);

    auto* form = new QFormLayout;

    // ---- Row 1: color scheme combo ------------------------------------------
    colorSchemeCombo_ = new QComboBox(this);
    colorSchemeCombo_->addItem(xvueT(MsgId::PrefColorSchemeSystem),
                               QStringLiteral("system"));
    colorSchemeCombo_->addItem(xvueT(MsgId::PrefColorSchemeLight),
                               QStringLiteral("light"));
    colorSchemeCombo_->addItem(xvueT(MsgId::PrefColorSchemeDark),
                               QStringLiteral("dark"));
    {
        const QString cur = XvuePrefs::colorScheme();
        const int     ix  = colorSchemeCombo_->findData(cur);
        colorSchemeCombo_->setCurrentIndex(ix < 0 ? 0 : ix);
    }
    form->addRow(xvueT(MsgId::PrefColorSchemeLabel), colorSchemeCombo_);

    // ---- Row 2: console default visibility ----------------------------------
    consoleVisibleCheck_ = new QCheckBox(this);
    consoleVisibleCheck_->setChecked(XvuePrefs::consoleVisible(true));
    form->addRow(xvueT(MsgId::PrefConsoleDefaultLabel), consoleVisibleCheck_);

    // ---- Row 3: recent-projects count (v1 placeholder, disabled) ------------
    recentCountSpin_ = new QSpinBox(this);
    recentCountSpin_->setRange(1, 20);
    recentCountSpin_->setValue(10);
    recentCountSpin_->setEnabled(false);   // v1: not yet wired through
    form->addRow(xvueT(MsgId::PrefRecentCountLabel), recentCountSpin_);

    outer->addLayout(form);

    // ---- Buttons ------------------------------------------------------------
    auto* buttons = new QDialogButtonBox(
        QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
    if (auto* ok = buttons->button(QDialogButtonBox::Ok)) {
        ok->setText(xvueT(MsgId::ButtonOk));
    }
    if (auto* cancel = buttons->button(QDialogButtonBox::Cancel)) {
        cancel->setText(xvueT(MsgId::ButtonCancel));
    }
    connect(buttons, &QDialogButtonBox::accepted,
            this,    &XvuePreferencesDialog::onAccept);
    connect(buttons, &QDialogButtonBox::rejected,
            this,    &QDialog::reject);
    outer->addWidget(buttons);
}

XvuePreferencesDialog::~XvuePreferencesDialog() = default;

void XvuePreferencesDialog::onAccept() {
    XvuePrefs::saveColorScheme(colorSchemeCombo_->currentData().toString());
    XvuePrefs::saveConsoleVisible(consoleVisibleCheck_->isChecked());
    // recentCountSpin_->value() is intentionally NOT persisted (v1 placeholder).
    accept();
}
