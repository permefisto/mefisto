// xvue/qt/src/xvue_qt_menu_bridge.cpp
// Phase 6.0 Plan 02 (06.0-02): real bodies for queueLexicon/popChar.
// Replaces the Plan 01 no-op scaffold.
//
// Why each char + trailing CR (=13)? Verified against xvue/saclav.f §270-333:
// SACLAV's lexer expects one byte at a time and uses CR as the line
// terminator that flushes its KLG accumulator. Anything else breaks the
// lexicon's parse loop. UX-02 / UX-03 wiring (Plan 03 plugs popChar into
// XvueEventBridge::waitForEvent's pre-drain).
//
// Threading: every public mutator/reader lives on the main thread. The
// bridge is parented under XvueWindow which itself is main-thread; QAction
// triggers are delivered on the main thread; XvueEventBridge::waitForEvent
// runs on the main thread (Phase 5 D-02). We assert this with
// QThread::currentThread() == qApp->thread() — a violation indicates a
// downstream caller has gone off-thread and must be fixed at the source,
// not papered over with a mutex.
//
// Threat T-06.0-BRIDGE-01: kMaxQueueSize=10000 cap. A user mash-clicking
// a menu item while Fortran is stuck in a long compute would otherwise
// grow the queue without bound. Past the cap we silently drop new pushes;
// a filling queue is a symptom of a stuck Fortran loop, not a user error,
// and the dropped events would be unrecoverable garbage anyway.
#include "xvue_qt_menu_bridge.h"

#include <QCoreApplication>
#include <QMenu>
#include <QThread>

#include <utility>

XvueMenuBridge::XvueMenuBridge(QObject* parent) : QObject(parent) {}

XvueMenuBridge::~XvueMenuBridge() = default;

void XvueMenuBridge::queueLexicon(const QString& cmd) {
    // SHELL-07 — main-thread-only. The bridge is a QObject parented under
    // XvueWindow (main-thread). Crossing threads here would race popChar()
    // and the Fortran-side reader.
    Q_ASSERT(QThread::currentThread() == qApp->thread());

    // T-06.0-BRIDGE-01 cap. Cmd contributes cmd.length() bytes + 1 CR; if
    // accepting the push would breach kMaxQueueSize, drop the entire push
    // (atomic — never half-write a command, otherwise the trailing CR may
    // separate from its prefix and confuse SACLAV's parser).
    const int wouldBe = pendingChars_.size() + cmd.length() + 1;
    if (wouldBe > kMaxQueueSize) {
        return;
    }

    // 06.0-RESEARCH §6 + saclav.f §270-333: each char goes as one
    // QQueue<char> element; trailing CR (=13) is the Fortran line
    // terminator that flushes SACLAV's KLG(LHKLG) accumulator.
    for (QChar ch : cmd) {
        pendingChars_.enqueue(ch.toLatin1());
    }
    pendingChars_.enqueue(static_cast<char>(13));
}

std::optional<char> XvueMenuBridge::popChar() {
    Q_ASSERT(QThread::currentThread() == qApp->thread());
    if (pendingChars_.isEmpty()) {
        return std::nullopt;
    }
    return pendingChars_.dequeue();
}

void XvueMenuBridge::setContextMenuPopulator(ContextPopulator fn) {
    contextPopulator_ = std::move(fn);
}

void XvueMenuBridge::populateContextMenu(QMenu* m) {
    if (m && contextPopulator_) {
        contextPopulator_(m);
    }
}
