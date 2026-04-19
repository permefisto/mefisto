// xvue/qt/src/xvue_qt_preferences.h
// Phase 6.0 Plan 01 (scaffold): Single-panel preferences dialog declaration.
// Plan 04 fills the form layout (recent count, console default, color scheme)
// and wires accept() to XvuePrefs::save*.
#pragma once
#include <QDialog>

class QComboBox;
class QSpinBox;
class QCheckBox;

class XvuePreferencesDialog : public QDialog {
    Q_OBJECT
public:
    explicit XvuePreferencesDialog(QWidget* parent = nullptr);
    ~XvuePreferencesDialog() override;

private slots:
    void onAccept();

private:
    QComboBox* colorSchemeCombo_    = nullptr;
    QCheckBox* consoleVisibleCheck_ = nullptr;
    QSpinBox*  recentCountSpin_     = nullptr;
};
