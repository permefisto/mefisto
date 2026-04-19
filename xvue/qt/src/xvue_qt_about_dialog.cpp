// xvue/qt/src/xvue_qt_about_dialog.cpp
// Phase 6.0 Plan 04: native QMessageBox::about wrapper for UX-09 (D-09).
//
// Body composition (UI-SPEC §Copywriting AboutBody):
//   "MEFISTO %1\n(MEsh and Finite element Software TO…)\n\n"
//   "Auteur : Alain Perronnet\nLaboratoire Jacques-Louis Lions (LJLL), UPMC Paris\n\n"
//   "Qt %2\nCompilé le %3"
//
// Argument substitutions (positional %1/%2/%3):
//   %1 = MEFISTO_VERSION (compile-time macro, see below)
//   %2 = qVersion()      (Qt runtime version)
//   %3 = __DATE__        (compile-time build date)
//
// Plan 06 may centralize MEFISTO_VERSION into a versions header — for now we
// declare a local fallback so this file can land independently of Plan 06.
#include "xvue_qt_about_dialog.h"

#include "xvue_qt_i18n.h"

#include <QMessageBox>
#include <QString>
#include <QtGlobal>     // qVersion()

#ifndef MEFISTO_VERSION
#define MEFISTO_VERSION "v1.0-qt6"   // Plan 06 may override via -D on the cmd line
#endif

void XvueAboutDialog::show(QWidget* parent) {
    const QString title = xvueT(MsgId::AboutTitle);
    const QString body  = xvueT(MsgId::AboutBody)
                              .arg(QStringLiteral(MEFISTO_VERSION))
                              .arg(QString::fromLatin1(qVersion()))
                              .arg(QString::fromLatin1(__DATE__));
    // QMessageBox::about uses the platform's native About style; on most
    // desktops the body is auto-formatted with the application icon. Until
    // Plan 02 fills the i18n table, title and body are empty strings — the
    // dialog still constructs and dismisses cleanly.
    QMessageBox::about(parent, title, body);
}
