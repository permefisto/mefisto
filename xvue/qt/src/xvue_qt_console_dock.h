// xvue/qt/src/xvue_qt_console_dock.h
// Phase 6.0 Plan 01 (scaffold): QDockWidget hosting the redirected stdout/
// stderr console output. Plan 04 fills installStdoutRedirect with the real
// pipe()/dup2()/F_SETPIPE_SZ(1MB) + QSocketNotifier wiring per RESEARCH §3.
#pragma once
#include <QByteArray>
#include <QDockWidget>

class QPlainTextEdit;
class QSocketNotifier;
class XvueErrorBatcher;

class XvueConsoleDock : public QDockWidget {
    Q_OBJECT
public:
    explicit XvueConsoleDock(QWidget* parent = nullptr);
    ~XvueConsoleDock() override;

    // Plan 04 fills: pipe() + dup2 + setvbuf + F_SETPIPE_SZ(1MB) + QSocketNotifier.
    void installStdoutRedirect();

    // Public for test access; appends a line + triggers auto-scroll.
    void appendLine(const QString& line);

private slots:
    void onPipeReadable();

private:
    QPlainTextEdit*   textEdit_     = nullptr;
    QSocketNotifier*  notifier_     = nullptr;
    XvueErrorBatcher* errorBatcher_ = nullptr;
    QByteArray        partialLine_;
    int               readFd_       = -1;
};
