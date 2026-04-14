# Phase 4: Pixmap save/restore (double-buffering) - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions captured in 04-CONTEXT.md — this log preserves the analysis.

**Date:** 2026-04-14
**Phase:** 04-pixmap-save-restore-double-buffering
**Mode:** discuss (interactive)
**Areas analyzed:** Slot model, Flip ops semantics, Resize handling, mempxaccro placement, Validation strategy

## Domain Boundary Established

Off-screen pixmap slot supports the mesher's rubber-band drag workflow without flicker, reproducing the X11 double-buffering semantics on the Qt backend. Delivers 7 ABI entry points (`fenetremempx_`, `mempxfenetre_`, `sauvefenetre_`, `restaurefenetre_`, `sauvemempx_`, `restauremempx_`, `effacemempx_`) plus `XvueState::saved_canvas_` allocation/lifecycle, validated by a synthetic save/draw/restore round-trip test (no Phase 5 event bridge dependency).

## Carried Forward From Earlier Phases

- **Phase 2 D-04** — `XvueState::backing_` is the single visible surface; drawing primitives target it; paintEvent is one drawPixmap blit. This is THE load-bearing invariant for Phase 4 — it forces the no-op decision on `fenetremempx_`/`mempxfenetre_` (D-04) and the bit-identical bodies for `effacemempx_`/`effacer_` (D-10).
- **Phase 2 D-05** — long-lived `painter_` on `backing_`. Phase 4 helpers must NOT swap painter targets; the save direction uses a temporary scoped painter on `saved_canvas_`.
- **Phase 2 D-06** — DPR-aware backing. Phase 4 inherits the same setDevicePixelRatio() pattern for `saved_canvas_`.
- **Phase 2 D-09** — `xvtrait_` / `xvftrait_` collapsed. Phase 2 explicitly noted that Phase 4 might need to revisit this; the answer is **no** — Phase 4 confirms the collapse holds.
- **Phase 2 D-32 / D-33 / D-34** — `XVUE_QT_ASSERT_MAIN_THREAD`, signature-literal-copy, `verify_abi` + `verify_no_exec` invariants. All preserved.
- **Phase 03-04 close** — empirical validation that calling `MEMPXFENETRE` from Fortran tracers produces correct displays on Qt even when the symbol is a warn_once stub. This empirically backs Phase 4 D-04's no-op decision.

## Gray Areas Identified

1. **Slot model** — How many off-screen slots?
2. **Flip ops** — What should `fenetremempx_` / `mempxfenetre_` do?
3. **Resize handling** — What if canvas resizes between save and restore?
4. **mempxaccro placement** — Phase 4 or Phase 5?
5. **Validation strategy** — Interactive cavity2d test or synthetic round-trip?

## User Decisions

### Slot model
- **Question:** How many off-screen slots should `XvueState` own, and how should they be named?
- **Choice:** "1 slot, fixed name (Recommended)"
- **Rationale captured in D-01:** Legacy `xvuelc.c` only allocates one `mempxsauvfen`, no Fortran caller does nested save/restore, lowest complexity & blast radius.

### Flip ops semantics
- **Question:** What should `fenetremempx_` and `mempxfenetre_` actually do on Qt?
- **Choice:** "True no-op + suppress warn-once (Recommended)"
- **Rationale captured in D-04, D-05, D-06:** Qt's backing IS the visible surface (Phase 2 D-04). Phase 03-04 empirically validated that no behavior change is needed; only warn-once diagnostic noise is removed.

### Resize handling
- **Question:** What happens if the canvas resizes between `sauvefenetre_` and `restaurefenetre_`?
- **Choice:** "Reallocate slot to current size on restore (Recommended)"
- **Rationale captured in D-12, D-13, D-14:** Pragmatic — legacy X11 `XCopyArea` also clips silently when sizes differ. Lazy invalidation at restore time avoids coupling Phase 4 to Phase 2's resize logic.

### mempxaccro placement
- **Question:** The legacy code also has `mempxaccro` — where does it belong?
- **Choice:** "Defer to Phase 5 — picking/event bridge (Recommended)"
- **Rationale captured in domain boundary "out of scope" + deferred ideas:** `mempxaccro` is consumed by `xvsouris2_` which is Phase 5's symbol. Keeps Phase 4 narrow.

### Validation strategy
- **Question:** How should Phase 4 validate the rubber-band scenario given the SC-4 dependency on Phase 5?
- **Choice:** "Synthetic xvtest harness driving the API directly (Recommended)"
- **Rationale captured in D-15, D-16, D-17, D-18:** Round-trip pixel equality test — draw mesh, save, draw overlay, restore, compare to mesh-only capture. Reuses Phase 03-04's existing `bin/xvtest-capture.sh` and `bin/qt-capture.sh` harnesses. SC-4 (interactive cavity2d) explicitly deferred to Phase 5 with pointer.

## Corrections Made

None — all 5 questions answered with the recommended option on the first pass. The user chose the conservative, smallest-blast-radius path across the entire phase.

## Scope Creep Avoided

- The conversation surfaced `mempxaccro` as a related-but-different symbol. Rather than expand Phase 4 to cover it, it was deferred to Phase 5 where its consumer (`xvsouris2_`) lives.
- The conversation surfaced N-slot stacks and a real second draw target as alternatives. Both were considered and rejected at decision time, then captured in `<deferred>` so the rationale is not lost.
- ROADMAP SC-4 (interactive cavity2d) was explicitly NOT pulled into Phase 4 scope despite being listed as a Phase 4 success criterion — Phase 5 is the natural home for anything that needs real mouse motion.

## External Research

None — all questions resolved from codebase analysis (xvue/xvuelc.c legacy bodies, xvue/qt/src/xvue_qt_api.cpp current stubs, Phase 2 / 03 / 03-04 prior CONTEXT.md decisions). No library-specific knowledge gaps required Context7 or web search.

## Notes for Planner

- Phase 4 has a **tiny blast radius**: 7 entry-point bodies, 1 new `XvueState` field, 1 destructor cleanup line, 1 extension of `prpr/xvtest0.f`. No new files, no CMakeLists changes, no Fortran solver changes.
- The validation harness (D-15/D-16/D-17) reuses Phase 03-04 infrastructure; the planner should NOT introduce a new test framework.
- SC-4 is deferred to Phase 5 with explicit documentation. The Per-Task Verification Map should record SC-4 as `deferred-to-phase-5`, not `green` and not `red`.
