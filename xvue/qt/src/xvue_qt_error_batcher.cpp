// xvue/qt/src/xvue_qt_error_batcher.cpp
// Phase 6.0 Plan 01 (scaffold): single-shot timer member; bodies are no-ops.
// Plan 04 wires enqueue/flush to a QMessageBox + D-08 deferral.
#include "xvue_qt_error_batcher.h"

#include <QTimer>

XvueErrorBatcher::XvueErrorBatcher(QObject* parent)
    : QObject(parent), timer_(new QTimer(this))
{
    timer_->setSingleShot(true);
    connect(timer_, &QTimer::timeout, this, &XvueErrorBatcher::flushErrorBatch);
}

XvueErrorBatcher::~XvueErrorBatcher() = default;

void XvueErrorBatcher::enqueue(const QString& errorLine) {
    // Plan 04 fills: pendingErrorLines_.append(errorLine) + timer_->start(kBatchWindowMs).
    (void)errorLine;
}

void XvueErrorBatcher::flushErrorBatch() {
    // Plan 04 fills: QMessageBox::warning with the joined batch.
}
