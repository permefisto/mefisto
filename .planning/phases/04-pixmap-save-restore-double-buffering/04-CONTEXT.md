# Phase 4: Pixmap save/restore (double-buffering) - Context

**Gathered:** 2026-04-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Named off-screen pixmap slot supports the mesher's rubber-band drag workflow without flicker, reproducing the X11 double-buffering semantics on the Qt backend. Delivers the 7 ABI entry points for save/restore/clear of off-screen surfaces, validated by a synthetic save/draw/restore round-trip test. The **interactive** cavity2d mesher rubber-band test belongs in Phase 5 (event bridge); Phase 4 ships the API and a headless A/B-comparable validation.

In scope: `fenetremempx_`, `mempxfenetre_`, `sauvefenetre_`, `restaurefenetre_`, `sauvemempx_`, `restauremempx_`, `effacemempx_` real bodies; `XvueState::saved_canvas_` allocation; resize-handling on the slot; round-trip validation harness.

Out of scope: `mempxaccro` (snap-handle pixmap consumed by `xvsouris2_`) — belongs to Phase 5 picking. Interactive `testa/cavity2d` rubber-band drag — belongs to Phase 5 event bridge. Any change to `xvtrait_` / `xvftrait_` semantics — Phase 2 D-09 already collapsed them.

</domain>

<decisions>
## Implementation Decisions

### Slot model

- **D-01:** `XvueState` gains a single `QPixmap* saved_canvas_ = nullptr` field (raw pointer, manual lifecycle, matches Phase 2 D-04 ownership style). One slot is sufficient — the legacy `xvue/xvuelc.c` only allocates one `mempxsauvfen` (xvuelc.c:138 + xvuelc.c:1036) and no Fortran caller in `mail/`, `elas/`, `flui/`, `ther/`, `nlse/`, `util/`, `reso/`, or `prpr/` performs nested save/restore. Both `sauvefenetre_`/`restaurefenetre_` and `sauvemempx_`/`restauremempx_` operate on this same slot.
- **D-02:** The slot is allocated lazily on first use (first call to `sauvefenetre_` or `sauvemempx_`), not at `XvueState` construction. Avoids paying the QPixmap allocation cost for solver runs that never enter the rubber-band code path. Allocation is `new QPixmap(backing_->size())` with `setDevicePixelRatio(backing_->devicePixelRatio())` to match HiDPI semantics inherited from D-06 of Phase 2.
- **D-03:** The slot is destroyed in `XvueState::~XvueState()`, before `painter_` and `backing_` deletion. No risk of dangling references because the slot is never the active QPainter target.

### `fenetremempx_` / `mempxfenetre_` — flip ops

- **D-04:** Both are **true no-ops**: empty function bodies, no `update()`, no `processEvents()`, no warn_once. Qt's `XvueState::backing_` IS the visible surface (Phase 2 D-04, paintEvent = single drawPixmap blit). There is no fenetre/mempx distinction to flip on the Qt backend, so any caller that invokes the pair is satisfied by the backing already showing what it wanted.
- **D-05:** Each empty body carries a one-line comment pointing to Phase 2 D-04 and Phase 4 D-04 so a future reader doesn't mistake the no-op for a missing implementation. Format: `// D-04 (Phase 4): no-op — backing_ is the visible surface (Phase 2 D-04 single-backing collapse)`.
- **D-06:** The Phase 03-04 close-out empirically validated that calling `MEMPXFENETRE` from Fortran tracers (e.g., `elas/trelas.f:281` added during the Phase 3 reopen) produces correct displays on Qt even when the symbol is a warn-once stub. Phase 4 just upgrades the warn-once stub to a documented intentional no-op — no behavior change, only diagnostic noise removed.

### `sauvefenetre_` / `restaurefenetre_` / `sauvemempx_` / `restauremempx_` — save/restore

- **D-07:** All four entry points share two helper functions in `xvue_qt_api.cpp` (file-local, `static inline`):
  - `xvue_qt_save_to_slot()` — ensures `saved_canvas_` exists and matches `backing_` size, then `QPainter::drawPixmap(0, 0, *backing_)` from a temporary painter on `saved_canvas_`. The active long-lived `painter_` on `backing_` is NOT touched (it stays alive on `backing_` per Phase 2 D-05).
  - `xvue_qt_restore_from_slot()` — if `saved_canvas_` is null or the wrong size, log a stderr warning (`xvue-qt: restore_from_slot: no slot or size mismatch`) and return. Otherwise blit `saved_canvas_` → `backing_` via the active `painter_->drawPixmap(0, 0, *saved_canvas_)`, then schedule a `canvas_->update()` (no `processEvents` — match Phase 2 D-01 epilogue).
- **D-08:** `sauvefenetre_` and `sauvemempx_` are bit-identical bodies wrapping `xvue_qt_save_to_slot()`. `restaurefenetre_` and `restauremempx_` are bit-identical bodies wrapping `xvue_qt_restore_from_slot()`. The four ABI symbols are kept distinct (no symbol consolidation) to preserve `verify_abi` count and the legacy ABI surface.
- **D-09:** No `XVUE_QT_ASSERT_MAIN_THREAD` exemption — every entry point body still starts with the assertion (Phase 2 D-32 invariant).

### `effacemempx_` — clear

- **D-10:** `effacemempx_` body fills `backing_` with `state_->background_` via the active `painter_->fillRect(backing_->rect(), background_)`, followed by `canvas_->update()` (Phase 2 D-01 epilogue without `processEvents` — same lighter-weight epilogue as save/restore in D-07). Functionally identical to the existing Qt `effacer_` (Phase 2 D-15) because Qt has no separate mempx surface; the two symbols stay distinct for ABI preservation but share the same body.
- **D-11:** A one-line comment notes the equivalence: `// D-10 (Phase 4): same body as effacer_ — Qt has no separate mempx surface (Phase 2 D-04, D-15)`.

### Resize handling

- **D-12:** If `backing_->size()` differs from `saved_canvas_->size()` at restore time, the restore is a **no-op with stderr warning**. Mirrors the legacy X11 undefined behavior (XCopyArea silently clips when source/destination sizes differ) but makes it auditable. Pragmatic — the typical lifetime of a save/restore pair is shorter than a user resize, so the size-mismatch path is rare.
- **D-13:** No proactive invalidation on `XvueCanvas::resizeEvent`. The slot stays around until the next save (which reallocates if needed) or destruction. Adding resize-driven invalidation would couple Phase 4 to Phase 2's resize logic; the size check at restore time is sufficient.
- **D-14:** When a save reallocates because the size changed, the old `saved_canvas_` is `delete`-ed and a new one allocated. No attempt to preserve content across the size change; that's not what the legacy did either.

### Validation strategy

- **D-15:** Phase 4 ships a **synthetic save/restore round-trip test** that exercises the API headlessly without the Phase 5 event bridge. Approach:
  - Extend `prpr/xvtest0.f` with a Phase 4 coverage section (or add a sibling `xvtest0.f` block, planner's call):
    1. `xvinitgraphique`, draw a known mesh using existing Phase 2 primitives (e.g. a 4×4 checkered grid via `xvface_`)
    2. `sauvefenetre`
    3. Draw a "rubber band" overlay (e.g., a magenta rectangle via `xvbordrectangle_`) on top
    4. `restaurefenetre`
    5. Capture backing
  - The captured pixmap should be **bit-identical** to a control capture taken right after step 1 (mesh-only, no rubber-band overlay). The harness compares the two PNGs via `cmp` or a small pixel-diff helper.
- **D-16:** The same round-trip is repeated for `sauvemempx_`/`restauremempx_` and for `effacemempx_` (the latter clears to background, so the comparison is "after effacemempx == backing right after xvinitgraphique"). Three sub-tests, all headless, all integrated into a new `bin/xvtest0-pixmap-roundtrip.sh` harness or extended `bin/cbxvtest0_qt`.
- **D-17:** A/B comparison against the legacy X11 backend uses the existing `bin/xvtest-capture.sh` infrastructure from Phase 03-04. Each round-trip produces a Qt PNG and an X11 PNG; the orchestrator reads both via the `Read` tool to confirm visual equality.
- **D-18:** ROADMAP Success Criterion 4 (interactive cavity2d rubber-band test) is **explicitly deferred to Phase 5**. Phase 4's validation map records SC-4 as `deferred-to-phase-5` rather than `green`, with a comment pointing to Phase 5's planned cavity2d HUMAN-UAT entry. Phase 4 still closes because SC-1, SC-2, SC-3 are all covered by D-15/D-16/D-17.

### Claude's Discretion

- **Sanity-driver vehicle.** Whether Phase 4 extends `prpr/xvtest0.f` with a `PHASE 4 COVERAGE` block, adds a sibling Fortran driver, or introduces a small C++ test harness under `xvue/qt/tests/` is the planner's call. Phase 2 D-36 left the same choice open and chose the Fortran extension; Phase 4 should match for consistency unless a strong reason emerges.
- **Pixel-diff helper.** D-15 says "compare via `cmp` or a small pixel-diff helper". `cmp` works on the raw PNG bytes (which is brittle if PNG metadata varies) and `convert -metric AE` from ImageMagick is a more robust option. Either is fine; if the planner picks `convert`, ensure the harness probes for ImageMagick at the start and bails cleanly if absent.
- **Save-side helper layout.** D-07 specifies file-local `static inline` helpers in `xvue_qt_api.cpp`. If the helpers grow beyond ~15 lines each, splitting them out into `xvue_qt_state.cpp` is acceptable; the public API contract (which 7 ABI entry points exist) is fixed by the requirements.
- **Lazy allocation safety.** D-02 says lazy allocation. If the planner finds that a fully-eager allocation in `XvueState::XvueState()` simplifies error handling significantly, that's an acceptable tradeoff — the cost is one ~6KB QPixmap per process, which is negligible. The decision is laziness for memory hygiene, not correctness.
- **Stderr message wording.** D-12's "no slot or size mismatch" warning text is illustrative. The planner may pick wording that matches the existing diagnostic style in `xvue_qt_api.cpp`.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing Phase 4.**

### Phase scope anchors
- `.planning/ROADMAP.md` §"Phase 4: Pixmap save/restore (double-buffering)" — phase boundary, goal, depends-on (Phase 2), success criteria 1–4
- `.planning/REQUIREMENTS.md` §"Pixmap save/restore" PIXMAP-01 through PIXMAP-04 — the 4 requirements this phase delivers
- `.planning/REQUIREMENTS.md` Requirements traceability matrix lines 197–200 — PIXMAP-01..04 → Phase 4 mapping

### Legacy semantics being replicated
- `xvue/xvuelc.c:1307–1428` — the 7 legacy entry-point bodies (`fenetremempx_`, `mempxfenetre_`, `sauvefenetre_`, `restaurefenetre_`, `sauvemempx_`, `restauremempx_`, `effacemempx_`). Read these literally before implementing — they document the X11 `XCopyArea` semantics that the Qt bodies must reproduce or intentionally collapse.
- `xvue/xvuelc.c:138` (declaration) and `xvuelc.c:1036` (allocation) — the `mempxsauvfen` slot is a single `Pixmap`, not an array or stack. Confirms D-01's "1 slot, fixed name" choice.
- `xvue/xvuelc.c:1395–1428` — `effacemempx_` body, including the NULL-display batch-mode guard added in commit 3149e3f (Phase 03 prework). Phase 4's Qt body has no equivalent guard need because `painter_->isActive()` already covers the not-yet-allocated case.

### Phase 2 invariants (load-bearing for Phase 4)
- `.planning/phases/02-drawing-primitives-backing-pixmap/02-CONTEXT.md` D-04 — XvueState owns one QPixmap `backing_`; Phase 4 D-01 explicitly mirrors the same ownership style for `saved_canvas_`
- `.planning/phases/02-drawing-primitives-backing-pixmap/02-CONTEXT.md` D-05 — `painter_` lives continuously on `backing_`; Phase 4 must NOT swap targets, only blit FROM backing_ TO slot or from slot back to backing_
- `.planning/phases/02-drawing-primitives-backing-pixmap/02-CONTEXT.md` D-06 — DPR-aware backing; Phase 4 D-02 inherits the same DPR for the slot
- `.planning/phases/02-drawing-primitives-backing-pixmap/02-CONTEXT.md` D-09 — `xvtrait_` / `xvftrait_` semantically collapsed; confirms Phase 4 has nothing to do for those two (they don't appear in Phase 4's symbol set anyway)
- `.planning/phases/02-drawing-primitives-backing-pixmap/02-CONTEXT.md` D-15 — `effacer_` body; Phase 4 D-10 explicitly notes that `effacemempx_` is bit-identical
- `.planning/phases/02-drawing-primitives-backing-pixmap/02-CONTEXT.md` D-32, D-33, D-34 — `XVUE_QT_ASSERT_MAIN_THREAD` invariant + ABI signature literal-copy + `verify_abi` and `verify_no_exec` cmake targets must stay green

### Phase 03-04 capture infrastructure (reused by Phase 4 validation)
- `bin/xvtest-capture.sh` — Xvfb + sentinel + import harness for legacy X11 captures
- `bin/qt-capture.sh` — `QT_QPA_PLATFORM=offscreen` + backing-pixmap grab harness for Qt
- `xvue/qt/src/xvue_qt_api.cpp:543–640` — `xvfermer_` `MEFISTO_QT_CAPTURE_PATH` hook used by qt-capture.sh; Phase 4 round-trip tests reuse this path

### Phase 5 boundary handoff (informs Phase 4 deferrals)
- `.planning/ROADMAP.md` §"Phase 5: Event bridge & blocking reads" — explicit Phase 5 scope for `xvsouris_` / `xvsouris2_` and the cavity2d interactive HUMAN-UAT. Phase 4 D-18 references this for the SC-4 deferral.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `XvueState` (xvue/qt/src/xvue_qt_state.h, .cpp) — already owns `backing_`, `painter_`, `background_`, `foreground_`, `pen_`, `brush_`. Phase 4 adds exactly one new field (`saved_canvas_`) plus its destructor cleanup. The applyPen/foreground/background machinery is untouched.
- `XvueCanvas::resizeEvent` (xvue/qt/src/xvue_qt_canvas.cpp:40–113) — already handles backing_ reallocation + content preservation via top-left blit. Phase 4 D-13 deliberately does NOT hook into this path; the slot is invalidated lazily at restore time instead.
- `xvue_qt_draw_rect_common` (xvue/qt/src/xvue_qt_api.cpp:58–75) — the existing rect-draw helper pattern is the model for the new `xvue_qt_save_to_slot` / `xvue_qt_restore_from_slot` file-local helpers (D-07).
- `bin/xvtest-capture.sh` + `bin/qt-capture.sh` (Phase 03-04) — the round-trip validation harness in D-15/D-16 reuses these unchanged.
- `verify_abi` cmake target — already counts 57 ABI symbols; Phase 4 changes 7 stub bodies to real bodies, leaving the symbol count unchanged.

### Established Patterns
- **Single-backing collapse** (Phase 2 D-04) — the most load-bearing Phase 2 invariant for Phase 4. Drawing primitives target backing_; paintEvent is one drawPixmap blit. Phase 4 D-04 (no-op flip ops) and D-10 (effacemempx == effacer) both inherit directly from this.
- **Long-lived QPainter** (Phase 2 D-05) — `painter_` stays bound to `backing_`. Phase 4 helpers use a TEMPORARY QPainter on `saved_canvas_` for the save direction (because `painter_` cannot be active on two targets at once); the restore direction uses the existing active `painter_` on `backing_`.
- **HiDPI via setDevicePixelRatio** (Phase 2 D-06) — Phase 4 D-02 inherits the same DPR pattern when allocating `saved_canvas_`.
- **`XVUE_QT_ASSERT_MAIN_THREAD` first line in every body** (Phase 2 D-32) — Phase 4 D-09 keeps this invariant.
- **One-line `// D-NN comment` pointing to the deciding ADR** — used throughout Phase 1, 2, 3 for non-obvious choices. Phase 4 D-05 and D-11 follow the same pattern.

### Integration Points
- `xvue/qt/src/xvue_qt_api.cpp` — the 7 entry-point bodies (currently warn_once stubs at lines 394, 403, 411, 419, 427, 435, 443). Phase 4 replaces all 7 in place. No new file is created.
- `xvue/qt/src/xvue_qt_state.h` — add `QPixmap* saved_canvas_ = nullptr;` field. One-line addition near the existing `backing_` declaration.
- `xvue/qt/src/xvue_qt_state.cpp` — extend the destructor to delete `saved_canvas_`. One-line addition.
- `prpr/xvtest0.f` — extend with a `PHASE 4 COVERAGE` section per D-15. Coordinated with the Phase 4 validation harness.
- `bin/cbxvtest0_qt` — already builds xvtest0; no changes needed unless the planner picks a sibling driver vehicle (then a new cb script is needed).
- No changes to `bin/cbl_tout`, `bin/cbl_tout_qt`, `xvue/qt/CMakeLists.txt`, or any Fortran solver source. Phase 4's blast radius is intentionally tiny.

</code_context>

<specifics>
## Specific Ideas

- "Lowest complexity, smallest blast radius" — preferred posture, restated from D-01 picking single-slot over N-slots.
- "Match Phase 2 D-04's single-backing collapse" — the architectural anchor that drives D-04 (no-op flip) and D-10 (effacemempx == effacer).
- "Explicit deferral with documented pointer to Phase 5" — D-18 over silently-deferred SC-4.
- "Round-trip pixel equality" as the validation primitive (D-15) — reuses Phase 03-04's existing capture infrastructure rather than introducing a new test framework.

</specifics>

<deferred>
## Deferred Ideas

- **`mempxaccro` (snap-handle pixmap)** — the 13×13 pixmap consumed by `xvsouris2_` at `xvuelc.c:2415–2450` for interactive picking snap indicators. Belongs to Phase 5 (event bridge & blocking reads) which owns `xvsouris_` / `xvsouris2_`. See ROADMAP.md Phase 5 scope.
- **Interactive cavity2d rubber-band-drag HUMAN-UAT** — ROADMAP SC-4. Requires the Phase 5 event bridge to dispatch real mouse-motion events into the Fortran mesher. Phase 4 D-18 explicitly defers this with a pointer to Phase 5.
- **N-slot named pixmap stack** — rejected at D-01 in favor of single slot. Reconsider if a future caller (likely Phase 6 or 7) needs nested save/restore that a single slot can't cover. The decision is reversible — adding more slots is mostly additive.
- **Real second draw target (`mempx_`)** — rejected at D-04. Would require reverting Phase 2 D-04's single-backing collapse and updating every drawing primitive to target `mempx_` instead of `backing_`. The Phase 03-04 close empirically validated that no caller actually needs this on Qt, so it stays deferred indefinitely unless a real failure case surfaces.
- **Stale Fortran latent UB audit** — separate hardening item carried over from Phase 03-04 close (uninitialized `TPSINI` in `ther/thed1t.f` and similar). Out of Phase 4 scope; tracked under the libgfortran5 hold deferred item.

</deferred>

---

*Phase: 04-pixmap-save-restore-double-buffering*
*Context gathered: 2026-04-14*
