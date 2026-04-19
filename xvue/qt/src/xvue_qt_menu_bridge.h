// xvue/qt/src/xvue_qt_menu_bridge.h
// Phase 6.0 Plan 01 (scaffold): XvueMenuBridge — translates QAction triggers
// into the Fortran-side text lexicon by enqueueing characters that
// xvsouris_/xvsouris2_ can dequeue. Plan 02 fills queueLexicon/popChar with
// real implementations; this file declares the API only.
//
// §10 layer 2: also tracks "module registered" sentinel so Plan 06's
// xvue_module_init_ can verify per-module action wiring before the menu bar
// is shown.
#pragma once
#include <QObject>
#include <QQueue>
#include <QString>
#include <functional>
#include <optional>

class QMenu;

class XvueMenuBridge : public QObject {
    Q_OBJECT
public:
    explicit XvueMenuBridge(QObject* parent = nullptr);
    ~XvueMenuBridge() override;

    // Plan 02 fills body: push each char of cmd + trailing CR (13).
    void queueLexicon(const QString& cmd);
    // Plan 02 fills body: dequeue one ASCII code; nullopt if empty.
    std::optional<char> popChar();

    // Hooks populated by 6.1..6.5 per UI-SPEC §Canvas right-click.
    using ContextPopulator = std::function<void(QMenu*)>;
    void setContextMenuPopulator(ContextPopulator fn);
    void populateContextMenu(QMenu* m);

    // §10 layer 2 — registration sentinel.
    bool hasRegisteredModule() const { return moduleRegistered_; }
    void markModuleRegistered()      { moduleRegistered_ = true; }

    // T-06.0-BRIDGE-01 bound (Plan 02 will enforce against this).
    static constexpr int kMaxQueueSize = 10000;
    int queueSize() const { return pendingChars_.size(); }

private:
    QQueue<char>     pendingChars_;
    ContextPopulator contextPopulator_;
    bool             moduleRegistered_ = false;
};
