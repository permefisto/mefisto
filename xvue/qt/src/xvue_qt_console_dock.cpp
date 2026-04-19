// xvue/qt/src/xvue_qt_console_dock.cpp
// Phase 6.0 Plan 01 (scaffold): minimal QPlainTextEdit-in-QDockWidget body.
// Plan 04 fills installStdoutRedirect/onPipeReadable.
#include "xvue_qt_console_dock.h"

#include <QFontDatabase>
#include <QPlainTextEdit>
#include <QTextCursor>

XvueConsoleDock::XvueConsoleDock(QWidget* parent)
    : QDockWidget(parent)
{
    // Plan 04 will set windowTitle to xvueT(MsgId::ConsoleTitle); empty for now.
    setWindowTitle(QString());

    textEdit_ = new QPlainTextEdit(this);
    textEdit_->setReadOnly(true);
    textEdit_->setMaximumBlockCount(10000);
    textEdit_->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
    setWidget(textEdit_);
}

XvueConsoleDock::~XvueConsoleDock() = default;

void XvueConsoleDock::installStdoutRedirect() {
    // Plan 04 fills: pipe()/dup2()/F_SETPIPE_SZ(1<<20)/QSocketNotifier.
}

void XvueConsoleDock::appendLine(const QString& line) {
    if (!textEdit_) return;
    textEdit_->appendPlainText(line);
    textEdit_->moveCursor(QTextCursor::End);
}

void XvueConsoleDock::onPipeReadable() {
    // Plan 04 fills: read(readFd_, ...) + dispatch lines vs. error batcher.
}
