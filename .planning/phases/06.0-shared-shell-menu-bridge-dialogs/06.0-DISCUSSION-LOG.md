# Phase 6: Level-3 UX chrome & menu surface - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in 06-CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-04-19
**Phase:** 06-level-3-ux-chrome-menu-surface
**Mode:** discuss (interactive)
**Areas discussed:** Phase split, Menu taxonomy, LEXICON-AUDIT scope, Shortcuts + typed-lexicon coexistence

---

## Gray-area selection

The orchestrator surfaced 4 gray areas based on Phase 6's ROADMAP goal, UX-01..UX-13, and Phase 5 carry-forward decisions. Locked architecture (XvueMenuBridge queue, blockingDepth modal guard, QFileDialog, QSettings, bilingual via `td/m/anglais`, `QDockWidget` console, `*** ERREUR` parsing, wheel-zoom / middle-pan / right-click / live coords, `QPalette` dark-mode) was not discussed — already decided by requirement text.

**User selected:** all 4 (Phase split, Menu taxonomy, LEXICON-AUDIT scope, Shortcuts + typed-lexicon).

---

## Phase split

| Option | Description | Selected |
|--------|-------------|----------|
| Monolithic Phase 6 | One discuss/plan/execute cycle ships all 13 UX requirements. 8-12 plans across 4-5 waves. Highest coherence, longest elapsed time. | |
| Split: 6.0 + 6.1..6.5 (Recommended) | 6.0 ships shell + XvueMenuBridge + universal dialogs + persistence + dark-mode (~5 plans). 6.1..6.5 each ship that module's lexicon audit + QAction wiring (~2-3 plans each). | ✓ |
| Hybrid: one phase, waves per module | Single Phase 6 but plans organized as Wave 0 (shell+bridge) then Waves 1-5 (one per module). No phase-boundary overhead but shell-rework contamination risk. | |

**User's choice:** Split: 6.0 + 6.1..6.5 (Recommended)
**Notes:** Triggers ROADMAP.md update before `/gsd-plan-phase` — Phase 6 entry must be replaced with 6 new entries.

---

## Menu taxonomy

| Option | Description | Selected |
|--------|-------------|----------|
| Unified static (File, Edit, View, Mesh, Solve, Help) | One fixed top-level taxonomy for all modules; items enabled/disabled per active module. Consistent UX, simpler routing. | |
| Per-module dynamic bar (Recommended) | Each pp*_qt defines its own menu bar. Matches mental model (mesher ≠ solver). 6.0 ships File/View/Help mechanism; each 6.1..6.5 declares module menus. | ✓ |
| Hybrid shared + module menu | 6.0 ships fixed {File, Edit, View, Help}; each 6.1..6.5 contributes ONE module-specific menu. Consistent shell + module identity. | |

**User's choice:** Per-module dynamic bar (Recommended)
**Notes:** Different `pp*_qt` genuinely have different concerns; grayed-out menus would be visual noise.

---

## LEXICON-AUDIT scope

| Option | Description | Selected |
|--------|-------------|----------|
| Exhaustive QActions | Every typed command gets a QAction + menu entry + toolbar icon + keyboard shortcut. Highest ceremony, scales with lexicon size. | |
| Pragmatic: audit exhaustive, QActions frequency-weighted (Recommended) | LEXICON-AUDIT.md catalogs every typed command (full documentation). QActions wired only for 80/20 subset (~20-40 per module). Long-tail keeps working via typed lexicon. | ✓ |
| Minimal (only top-10 per module) | QActions for ~10 most-critical commands only. Fastest, least complete. | |

**User's choice:** Pragmatic: audit exhaustive, QActions frequency-weighted (Recommended)
**Notes:** Preserves MEFISTO's zero-regression invariant — typed lexicon is the safety net, QActions are the modern UX layer on top.

---

## Shortcuts + typed-lexicon coexistence

| Option | Description | Selected |
|--------|-------------|----------|
| Modifier rule: Ctrl/Alt/Cmd → QAction, plain chars → lexicon (Recommended) | Modifier combos fire QActions; plain alphanumeric + `;` + digits + Esc/Return go to typed lexicon via XvueEventBridge. F-keys route to QAction. Standard CUA conventions. | ✓ |
| Context-sensitive (focus-aware) | `QShortcut::setContext` disables global shortcuts when canvas has focus. More flexible but complex edge cases. | |
| Canvas-first: plain keys always → lexicon, modifiers → QAction, no F-keys as QActions | Strictest preservation of typed UX. Simplest to reason about. | |

**User's choice:** Modifier rule (Recommended)
**Notes:** Phase 5 D-06 (Esc=27, @=64 abort) continues unchanged. `99;` / `5;90;1;` muscle memory preserved exactly.

---

## Claude's Discretion defaults (presented, user confirmed)

| Default | Description |
|---------|-------------|
| Dark-mode scope | System-follow via `QPalette`, chrome only; canvas untouched. User-toggle via Preferences (QSettings-backed). |
| Recent-projects list | 10 most-recent, File → Recent Projects submenu, "Clear Recent" action. |
| Console dock | Visible by default, auto-scroll, copy-to-clipboard, session-local. |
| Modal refuse UX | 3-second status-bar message `"Finish current operation first (type 99;)"` when `blockingDepth() > 0`. |
| About dialog | Credits Alain Perronnet / LJLL / UPMC Paris; MEFISTO version + Qt version + build date. |
| Toolbar icons | Qt built-in `QStyle::StandardPixmap` where applicable; custom SVG in `xvue/qt/resources/icons/` where none fits. |

**User's choice:** "I'm ready for context" — no overrides requested.

---

## Deferred Ideas

Surfaced from scope-creep detection during analysis (not discussion):
- V1X-07 progress-bar fed by solver stdout parsing → Phase 9 or future
- V1X-08 animation scrubber for time-stepping solvers → Phase 9 or future
- V1X-09 snapshot undo for mesh edits → Phase 9 or future
- Canvas colormap customization dialog → out of scope
- Plugin / extension system → out of scope
- Multi-document / multi-canvas → out of scope

---

*Phase: 06-level-3-ux-chrome-menu-surface*
*Discussion log written: 2026-04-19*
