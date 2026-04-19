// xvue/qt/src/xvue_qt_error_batcher.cpp
// Phase 6.0 Plan 04: 500 ms cascade-batching for *** ERREUR alerts.
//
// Policy (UI-SPEC §Console dock + threat T-06.0-04-01):
//   N error lines arriving within kBatchWindowMs collapse into ONE
//   QMessageBox::warning whose body is the lines joined by '\n'. This caps
//   user-facing dialog spam at ~2 boxes/second worst case.
//
// D-08 / T-06-03-01 mitigation:
//   When XvueApp::blockingDepth() > 0 (Fortran is currently inside a nested
//   xvsouris_/xvsouris2_/xvpause_ wait), opening a modal QMessageBox would
//   re-enter the QEventLoop and risk eating the deferred-quit timer Phase 5
//   uses for motion coalescing. We defer the flush by re-arming the timer
//   for another window and retrying.
#include "xvue_qt_error_batcher.h"

#include "xvue_qt_app.h"
#include "xvue_qt_i18n.h"

#include <QCoreApplication>
#include <QMessageBox>
#include <QTimer>

XvueErrorBatcher::XvueErrorBatcher(QObject* parent)
    : QObject(parent), timer_(new QTimer(this))
{
    timer_->setSingleShot(true);
    timer_->setInterval(kBatchWindowMs);
    connect(timer_, &QTimer::timeout, this, &XvueErrorBatcher::flushErrorBatch);
}

XvueErrorBatcher::~XvueErrorBatcher() = default;

void XvueErrorBatcher::enqueue(const QString& errorLine) {
    pendingErrorLines_ << errorLine;
    // Arm the batching window on the FIRST line of a cascade. Subsequent
    // enqueue() calls inside the window simply append; the timer keeps its
    // original deadline. T-06.0-04-01 — bounds dialogs at 1 per kBatchWindowMs.
    if (!timer_->isActive()) {
        timer_->start();
    }
}

void XvueErrorBatcher::flushErrorBatch() {
    if (pendingErrorLines_.isEmpty()) return;

    // T-06-03-01 mitigation: defer the QMessageBox while a Fortran blocking
    // read is in progress. Re-arming the timer keeps the batch alive without
    // dropping any error lines.
    if (XvueApp::blockingDepth() > 0) {
        timer_->start();
        return;
    }

    const QString title = xvueT(MsgId::ErrorMsgBoxTitle);
    const QString body  = pendingErrorLines_.join(QChar::fromLatin1('\n'));
    pendingErrorLines_.clear();

    // Fire the QMessageBox via QueuedConnection so we exit the timer slot
    // before the modal exec begins. This avoids subtle reentry where the
    // modal would be opened from inside QTimer::timeout dispatch.
    QMetaObject::invokeMethod(qApp, [title, body]{
        QMessageBox::warning(nullptr, title, body);
    }, Qt::QueuedConnection);
}
