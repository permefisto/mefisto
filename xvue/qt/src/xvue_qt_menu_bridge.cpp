// xvue/qt/src/xvue_qt_menu_bridge.cpp
// Phase 6.0 Plan 01 (scaffold): stub bodies for XvueMenuBridge. Plan 02
// replaces queueLexicon/popChar with the real text-lexicon plumbing.
#include "xvue_qt_menu_bridge.h"
#include <QMenu>

XvueMenuBridge::XvueMenuBridge(QObject* parent) : QObject(parent) {}

XvueMenuBridge::~XvueMenuBridge() = default;

void XvueMenuBridge::queueLexicon(const QString& cmd) {
    // Plan 02 will fill: push each char of cmd + trailing CR (13).
    (void)cmd;
}

std::optional<char> XvueMenuBridge::popChar() {
    // Plan 02 will fill: pop the front of pendingChars_.
    if (pendingChars_.isEmpty()) return std::nullopt;
    char c = pendingChars_.dequeue();
    return c;
}

void XvueMenuBridge::setContextMenuPopulator(ContextPopulator fn) {
    contextPopulator_ = std::move(fn);
}

void XvueMenuBridge::populateContextMenu(QMenu* m) {
    if (contextPopulator_ && m) contextPopulator_(m);
}
