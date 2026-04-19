// xvue/qt/src/xvue_qt_about_dialog.h
// Phase 6.0 Plan 01 (scaffold): static launcher wrapping QMessageBox::about.
// Per D-09 + UX-09, no custom QDialog — uses native About style.
#pragma once

class QWidget;

class XvueAboutDialog {
public:
    // Plan 04 fills: QMessageBox::about(parent, AboutTitle, AboutBody) with
    // the i18n-resolved strings from xvueT(MsgId::*).
    static void show(QWidget* parent);
};
