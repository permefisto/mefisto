# Phase 1: Window shell - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions captured in `01-CONTEXT.md` — this log preserves the back-and-forth.

**Date:** 2026-04-11
**Phase:** 01-window-shell-xvueapp-xvuewindow-xvuecanvas
**Mode:** discuss (delegated via "think and take the best approach")
**Areas discussed:** Window open timing & sizing; xvfermer_ and reopen semantics; Assertions & exec() ban enforcement; xvfond_ + SHELL-04 source + test driver (grouped)

## Prior Context Loaded

- `.planning/PROJECT.md` — core value "Fortran must not notice the change"; Qt 6; CMake scoped to xvue/; parallel X11 build
- `.planning/REQUIREMENTS.md` §SHELL-01..07 — 7 requirements mapped to Phase 1
- `.planning/ROADMAP.md` §Phase 1 — goal, success criteria, depends-on Phase 0
- `.planning/phases/00-build-skeleton-abi-stubs/00-CONTEXT.md` — Phase 0 decisions carried forward (D-01/D-02 xvue/qt/ layout, D-04 single ABI header, D-05 trailing-underscore macro, D-08/D-09 shell-script pattern, D-12 verify_abi CMake target precedent, D-17 warn-once stub pattern, Claude discretion for `XVUE_QT_ASSERT_MAIN_THREAD()` skeleton)
- `.planning/research/ARCHITECTURE.md` — XvueApp/Window/Canvas/State component split, `std::call_once` singleton recipe, atexit teardown
- `.planning/research/PITFALLS.md` — Pitfall 5 (QApplication double-init), Pitfall 6 (exec() inverts control flow), Pitfall 7 (processEvents starvation), Pitfall 8 (nested QDialog::exec re-entrancy)
- `.planning/codebase/STRUCTURE.md`, `.planning/codebase/CONVENTIONS.md` — bin/ shell-script and prpr/ driver conventions
- Direct reads of `xvue/xvuelc.c` lines 286–303 (xvinitgraphique_), 319–334 (xvpxecran_), 337–356 (xvmmecran_), 612–980 (xvinfo_ incl. XCreateWindow at ~943), 935 (BlackPixel default), 1434 (xvfond_), 1602 (xvfermer_)

## Gray Areas Offered

1. **Window open timing & sizing** — what xvinitgraphique_ does, what xvinfo_ does in Phase 1, initial window size
2. **xvfermer_ + reopen semantics** — hide vs close vs delete
3. **Assertions & exec() ban enforcement** — SHELL-07 scope, SHELL-03 mechanism
4. **xvfond_ + SHELL-04 source + test driver** — grouped loose ends

## User Response

User selected "Other" with note: **"think and take the best approach"** — delegated all four areas to Claude with authority to pick recommended defaults and lock them as decisions.

## Resolutions (Claude-chosen per delegation)

### Area 1 — Window open timing & sizing

| Decision | Chosen | Alternatives considered | Rationale |
|----------|--------|------------------------|-----------|
| D-01 | xvinitgraphique_ does app+window+show+processEvents(ExcludeUserInputEvents) | (a) app only, defer window to xvinfo_; (b) app+window, no show until first paint | ROADMAP success criterion 1 requires a visible window after xvinitgraphique_; ExcludeUserInputEvents flag avoids re-entrancy (Pitfalls 6, 8) |
| D-02 | Initial size 800×600, title "MEFISTO" | (a) QScreen * 0.8; (b) 1024×768; (c) no size until xvinfo_ | Sentinel size fits any monitor, distinct from common defaults, easy to spot as Phase-1 state |
| D-03 | xvinfo_ is partial real (resize only), palette outputs stay warn-once | (a) pure stub; (b) full real implementation | Fortran callers drive window sizing through xvinfo_; palette plumbing is Phase 3 work; partial impl keeps the hook alive without importing Phase 3 scope |
| D-04, D-05 | Minimal XvueState { QColor background_ = Qt::black; }; XvueCanvas::paintEvent fills with it | (a) QPalette on widget; (b) no state struct, hardcoded black | Clean additive path to Phase 2 (adds painter/pen/brush); one-line paintEvent swap in Phase 2 |

### Area 2 — xvfermer_ and reopen semantics

| Decision | Chosen | Alternatives considered | Rationale |
|----------|--------|------------------------|-----------|
| D-06 | xvfermer_ = window_.reset() (destroy XvueWindow, keep QApplication) | (a) window_->hide(); (b) window_->close() with WA_DeleteOnClose; (c) keep window, reset state only | Research (ARCHITECTURE.md §Singleton discipline) recommends "re-creates a new window inside the still-alive QApplication"; matches mesher real-world flow; synchronous teardown is observable |
| D-07 | Second xvinitgraphique_ detects null window and allocates fresh one | — | Natural consequence of D-06 |
| D-08 | std::atexit handler installed on first XvueApp::ensure() | (a) destructor at static-dtor time; (b) manual call from main | Qt embedding idiom; avoids static-destruction-order hostility |
| D-09 | std::call_once guards QApplication only; plain `if (!window_)` guards window | (a) call_once for both | call_once is one-shot by design; window must be reconstructable 0-to-N times |

### Area 3 — Assertions & exec() ban enforcement

| Decision | Chosen | Alternatives considered | Rationale |
|----------|--------|------------------------|-----------|
| D-10 | CMake custom target `verify_no_exec` post-build grep-and-fail | (a) git pre-commit hook; (b) both | Follows Phase 0 D-12 verify_abi precedent exactly; portable, unbypassable, single-developer friendly |
| D-11 | No git pre-commit hook | — | Would require install coordination; can be skipped with --no-verify; CMake guard is strictly stronger |
| D-12 | Retrofit XVUE_QT_ASSERT_MAIN_THREAD() into all ~60 existing stubs in xvue_qt_api.cpp | (a) Phase 1 entries only, add retroactively per phase; (b) only on real impls | Literal compliance with SHELL-07 "every extern C entry point"; single-file bulk edit; removes Phase 2+ onboarding friction |
| D-13 | Macro includes null-guard `if (qApp)` to handle the very-first-call window | — | First xvinitgraphique_ call enters with qApp == nullptr; null-guard lets the macro be the first statement uniformly |

### Area 4 — xvfond_, SHELL-04, test driver

| Decision | Chosen | Alternatives considered | Rationale |
|----------|--------|------------------------|-----------|
| D-14 | xvfond_ stores in XvueState.background_ via minimal 2-entry hardcoded table (0→black, 1→white) | (a) QPalette on widget; (b) warn-once no-op until Phase 3 | Matches legacy X11 BlackPixel/WhitePixel convention at xvuelc.c:935; disappears cleanly in Phase 3 |
| D-15 | XvueCanvas holds raw ptr to XvueState owned by XvueWindow; xvfond_ calls canvas_->update() | — | Automatic observation of next paint; consistent with Phase 2's state-lives-in-window model |
| D-16 | xvpxecran_ → QGuiApplication::primaryScreen()->size() | (a) window->screen(); (b) conditional on window existence | Deterministic, works before window creation, Phase 1 explicitly defers multi-monitor |
| D-17 | xvmmecran_ → primaryScreen()->physicalSize() | (a) size()/logicalDotsPerInch() math | Qt 6 guarantees physicalSize() in mm; no manual unit conversion |
| D-18 | Every real Phase 1 entry point starts with XvueApp::ensure() before XVUE_QT_ASSERT_MAIN_THREAD() | — | Guarantees qApp non-null for the assertion and for primaryScreen() access |
| D-19 | New prpr/xvtest0.f: CALL XVINITGRAPHIQUE; CALL XVFERMER; CALL XVINITGRAPHIQUE; CALL XVFERMER; STOP | (a) C++ test main; (b) reuse pp/ppmail_qt with lexicon input | Smallest thing exercising the full stack (gfortran link + extern C + reopen cycle); follows existing prpr/xvtest{1..4}.f convention reserved for Phase 2 drawing drivers |
| D-20 | New bin/cbxvtest0_qt cloned from a Phase 0 variant; bin/cbl_tout_qt NOT modified | (a) add ppxvtest0_qt to cbl_tout_qt; (b) make it a test-only CMake target | Preserves Phase 0 D-08's "big script builds five canonical modules" scope anchor |
| D-21 | Phase 1 validation gate: pp/ppxvtest0_qt runs, window visibly appears twice; 5 testa/ baseline cases NOT re-run until Phase 2 | (a) Full 5-case A/B sweep at end of Phase 1 | Nothing in Phase 1 changes legacy X11 behavior and the Qt backend has no drawing yet; the A/B sweep would be noise |

## Discretion Items Left Open

Recorded in the `<decisions>` §"Claude's Discretion" block of `01-CONTEXT.md`:

- Exact implementation of `XvueApp::ensure()` (internal organization)
- `XvueWindow` parent class bareness (empty menu bar? status bar?)
- `xvtest0.f` exact line count (optional diagnostic `PRINT *` calls)
- `verify_no_exec` exact CMake invocation form
- Header/source file fanout inside `xvue/qt/` (follow research's per-component split literally, or collapse)
- `QT_SCALE_FACTOR` shell var vs. explicit `Qt::AA_EnableHighDpiScaling` code call

## Deferred Ideas

All listed in `<deferred>` of `01-CONTEXT.md`: multi-monitor, QPalette background, git pre-commit hook, ppxvtest0_qt in cbl_tout_qt, defensive HiDPI attribute, extended XvueState fields, WA_DeleteOnClose+close() idiom, Phase 1 5-case A/B sweep.

## Scope Creep

None. All four areas stayed within the Phase 1 domain boundary (window lifecycle + screen metrics + background color + build-time guards). Drawing primitives, pixmap state, events, menus, and export were repeatedly identified as belonging to Phases 2–7 and kept out.
