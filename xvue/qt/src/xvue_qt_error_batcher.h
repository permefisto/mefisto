// xvue/qt/src/xvue_qt_error_batcher.h
// Phase 6.0 Plan 01 (scaffold): batches Fortran-side error lines into a
// single QMessageBox per UI-SPEC §Console dock (500 ms window). Plan 04 fills
// the QTimer + D-08 deferral logic.
#pragma once
#include <QObject>
#include <QStringList>

class QTimer;

class XvueErrorBatcher : public QObject {
    Q_OBJECT
public:
    explicit XvueErrorBatcher(QObject* parent = nullptr);
    ~XvueErrorBatcher() override;

    void enqueue(const QString& errorLine);

    // UI-SPEC §Console dock — coalesce window in milliseconds.
    static constexpr int kBatchWindowMs = 500;

private slots:
    void flushErrorBatch();

private:
    QStringList pendingErrorLines_;
    QTimer*     timer_ = nullptr;
};
