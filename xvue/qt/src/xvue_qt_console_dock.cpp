// xvue/qt/src/xvue_qt_console_dock.cpp
// Phase 6.0 Plan 04: real body for the stdout-redirected console dock.
//
// Architecture (per UI-SPEC §Console dock + RESEARCH §3):
//   POSIX pipe()  →  dup2(fd[1], STDOUT_FILENO)  →  Fortran WRITE(IMPRIM,*) flows
//                                                   into fd[1].
//   QSocketNotifier(fd[0], Read)                  →  on activate, drain bytes,
//                                                   split on '\n', forward each
//                                                   complete line to textEdit_
//                                                   (and *** ERREUR lines to
//                                                   the XvueErrorBatcher).
//
// Pitfall mitigations from 06.0-RESEARCH.md §3:
//   3.2  F_SETPIPE_SZ (1<<20) expands the pipe buffer from 64 KB to 1 MB so
//        long Fortran mesh reports do not back-pressure the main thread.
//   3.3  setMaximumBlockCount(10000) caps log memory.
//   3.4  setvbuf(_IOLBF) is called AFTER dup2 so it applies to the new stdout.
//   3.5  We always redirect — users wanting tty stdout will need a future
//        passthrough flag (out of Plan 04 scope).
#include "xvue_qt_console_dock.h"

#include "xvue_qt_error_batcher.h"
#include "xvue_qt_i18n.h"

#include <QFontDatabase>
#include <QPlainTextEdit>
#include <QSocketNotifier>
#include <QTextCursor>

#include <cstdio>
#include <cstdlib>
#include <fcntl.h>
#include <sys/types.h>
#include <unistd.h>

XvueConsoleDock::XvueConsoleDock(QWidget* parent)
    : QDockWidget(parent)
{
    // Plan 02 fills the FR/EN table; until then xvueT() returns "" and the
    // window keeps an empty title (Plan 06 will sync title on locale change).
    setWindowTitle(xvueT(MsgId::ConsoleTitle));

    textEdit_ = new QPlainTextEdit(this);
    textEdit_->setReadOnly(true);
    textEdit_->setMaximumBlockCount(10000);   // UI-SPEC cap (Pitfall 3.3)
    textEdit_->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
    setWidget(textEdit_);

    // T-06.0-04-01 — error batcher coalesces *** ERREUR cascades into ONE
    // QMessageBox per 500 ms window.
    errorBatcher_ = new XvueErrorBatcher(this);

    // UI-SPEC §Dock area policy: top area disallowed (would clash with toolbar).
    setAllowedAreas(Qt::BottomDockWidgetArea
                  | Qt::LeftDockWidgetArea
                  | Qt::RightDockWidgetArea);
}

XvueConsoleDock::~XvueConsoleDock() {
    if (notifier_) {
        notifier_->setEnabled(false);
    }
    // readFd_ is intentionally NOT closed: closing it would invalidate the
    // dup'd STDOUT_FILENO that Fortran is still writing to. The pipe lives
    // for the process lifetime — same discipline as Phase 1 D-08 (leak the
    // QApplication instead of racing libgfortran's atexit).
}

void XvueConsoleDock::installStdoutRedirect() {
    // Idempotent guard — only redirect once per process. A second call on the
    // same dock would leak the previous pipe and confuse the notifier.
    if (readFd_ >= 0) return;

    // Phase 6.0 Plan 06 test-isolation knob: when XVUE_QT_DISABLE_STDOUT_REDIRECT
    // is set in the environment, skip the dup2 over STDOUT_FILENO. This is
    // mandatory for QTest-driven test binaries that exercise xvue_module_init_
    // — otherwise QtTest's reporter writes to a pipe the test never reads
    // and the suite output disappears mid-run. Production pp*_qt invocations
    // never set this var; the redirect is the default behavior there.
    const char* disable = std::getenv("XVUE_QT_DISABLE_STDOUT_REDIRECT");
    if (disable && *disable && *disable != '0') {
        return;
    }

    int fd[2];
    if (::pipe(fd) != 0) {
        std::fprintf(stderr, "xvue-qt: pipe() failed; console dock disabled\n");
        return;
    }
    // T-06.0-CONSOLE-01 mitigation (RESEARCH Pitfall 3.2): expand pipe buffer
    // from default ~64 KB to 1 MB so long Fortran reports do not backpressure
    // the main thread. F_SETPIPE_SZ may silently clamp to /proc/sys/fs/pipe-
    // max-size; we ignore the return value and let the kernel apply what it can.
    (void)::fcntl(fd[0], F_SETPIPE_SZ, 1 << 20);   // 1 MB

    // Dup the write end onto STDOUT_FILENO so Fortran WRITE(IMPRIM,*) and any
    // C printf/fprintf to stdout flow into the pipe.
    if (::dup2(fd[1], STDOUT_FILENO) < 0) {
        std::fprintf(stderr, "xvue-qt: dup2 stdout failed; console dock disabled\n");
        ::close(fd[0]);
        ::close(fd[1]);
        return;
    }
    // Original write-end fd no longer needed — STDOUT_FILENO now refers to
    // the pipe's write side.
    ::close(fd[1]);

    // RESEARCH Pitfall 3.4: setvbuf MUST be called AFTER dup2 so it applies
    // to the new stdout file descriptor. _IOLBF flushes on each '\n' which
    // matches Fortran's WRITE(*,*) semantics — gives the dock real-time output.
    std::setvbuf(stdout, nullptr, _IOLBF, 0);

    readFd_   = fd[0];
    notifier_ = new QSocketNotifier(readFd_, QSocketNotifier::Read, this);
    connect(notifier_, &QSocketNotifier::activated,
            this,      &XvueConsoleDock::onPipeReadable);
}

void XvueConsoleDock::appendLine(const QString& line) {
    if (!textEdit_) return;
    textEdit_->appendPlainText(line);
    // Auto-scroll to bottom on every append (UI-SPEC §Console dock behavior).
    QTextCursor c = textEdit_->textCursor();
    c.movePosition(QTextCursor::End);
    textEdit_->setTextCursor(c);
}

void XvueConsoleDock::onPipeReadable() {
    if (readFd_ < 0) return;
    char    buf[4096];
    ssize_t n;
    // Drain until read() would block (n < requested) or returns 0/EOF.
    // QSocketNotifier is level-triggered; we still drain in a loop to amortize
    // the per-event cost when libgfortran flushes a long block.
    while ((n = ::read(readFd_, buf, sizeof buf)) > 0) {
        partialLine_.append(buf, static_cast<int>(n));
        drainPartialLineBuffer();
        if (n < static_cast<ssize_t>(sizeof buf)) break;   // pipe drained
    }
}

void XvueConsoleDock::feedRawBytesForTest(const QByteArray& bytes) {
    // Test-only path: bypass the OS pipe entirely so unit tests can drive the
    // line-splitting + *** ERREUR detection deterministically.
    partialLine_.append(bytes);
    drainPartialLineBuffer();
}

void XvueConsoleDock::drainPartialLineBuffer() {
    int pos;
    while ((pos = partialLine_.indexOf('\n')) >= 0) {
        const QByteArray rawLine = partialLine_.left(pos);
        partialLine_.remove(0, pos + 1);
        const QString line = QString::fromUtf8(rawLine);
        appendLine(line);
        // UI-SPEC §Console dock: lines matching ^\*\*\* ERREUR or
        // ^\*\*\* ERROR forward to the batcher for QMessageBox::warning.
        if (line.startsWith(QStringLiteral("*** ERREUR")) ||
            line.startsWith(QStringLiteral("*** ERROR"))) {
            if (errorBatcher_) {
                errorBatcher_->enqueue(line);
            }
        }
    }
}
