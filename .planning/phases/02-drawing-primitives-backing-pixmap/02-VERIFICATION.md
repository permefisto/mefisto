---
phase: 02-drawing-primitives-backing-pixmap
verified: 2026-04-11T20:48:00Z
status: passed
score: 9/9 DRAW requirements verified; 5/5 ROADMAP success criteria verified (SC-1 scope-limited to headless smoke + xvtest0.f, full xvtest1..4 parity explicitly deferred to Phase 3 per D-35)
overrides_applied: 0
re_verification:
  previous_status: none
  previous_score: n/a
  gaps_closed: []
  gaps_remaining: []
  regressions: []
deferred:
  - truth: "Full visual parity vs X11 on prpr/xvtest1.f..xvtest4.f (ROADMAP SC-1 literal reading)"
    addressed_in: "Phase 3 (text, fonts, colormap)"
    evidence: "02-CONTEXT.md D-35, 02-DISCUSSION-LOG.md:136, 02-VALIDATION.md:31 — xvtest1..4 call xvtexte_/xvchargefonte_/xvcouleur_ which remain warn-once stubs in Phase 2; output would be uninterpretable. Deferral documented in phase context BEFORE execution; ROADMAP SC-1 scope narrowed to the xvtest0.f driver + backing/painter invariants."
  - truth: "Full testa/ 5-case A/B re-run (related to SC-1)"
    addressed_in: "Phase 8 (A/B validation on testa subset)"
    evidence: "ROADMAP Phase 8 explicitly owns this; Phase 2 never claimed it."
---

# Phase 02: Drawing Primitives & Backing Pixmap — Verification Report

**Phase Goal (ROADMAP):** "All pure synchronous drawing primitives render into a persistent off-screen `QPixmap` via one long-lived `QPainter`, matching the X11 backend visually on `prpr/xvtest1.f`–`xvtest4.f`."

**Verified:** 2026-04-11
**Status:** passed
**Re-verification:** No — initial verification.

## Scope Note

The ROADMAP goal literally mentions `prpr/xvtest1.f`–`xvtest4.f`. Per `02-CONTEXT.md` D-35 and `02-DISCUSSION-LOG.md:136`, the phase scope was explicitly narrowed at planning time: those drivers call `xvtexte_` / `xvchargefonte_` / `xvcouleur_` which remain warn-once stubs until Phase 3, so their output would be uninterpretable as a parity gate. Phase 2 substitutes `prpr/xvtest0.f` (extended in plan 02-01) as the driver-level smoke gate and defers full `xvtest1..4` parity to Phase 3. This verification respects that documented scope narrowing.

## Goal Achievement

### Observable Truths (merged ROADMAP SC + must_haves)

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Running Qt backend produces lines, polylines, filled polygons, rectangles, and ellipse arcs visually indistinguishable from X11 on extended `xvtest0.f` driver (ROADMAP SC-1, scoped) | VERIFIED | `02-04-SUMMARY.md` Task 2 table: all 7 human checks PASS on live Qt 6 display 2026-04-11; headless `pp/ppxvtest0_qt` still exits 0 with zero DRAW warn-once lines (re-run now: `EXIT=0`). Full `xvtest1..4` parity deferred to Phase 3 per D-35. |
| 2 | `XvueState` holds exactly one `QPainter*` bound to the backing `QPixmap` for widget lifetime; `paintEvent` body is a single `drawPixmap(0, 0, backing)` call (ROADMAP SC-2 / DRAW-01) | VERIFIED | `xvue_qt_state.h:27-28` declares `QPixmap* backing_`, `QPainter* painter_`; `xvue_qt_canvas.cpp:22-30` paintEvent body is exactly `QPainter(this).drawPixmap(0, 0, *state_->backing_)` guarded by null-check; grep `painter_->begin` returns exactly ONE hit at `xvue_qt_canvas.cpp:70` (inside resizeEvent) — single long-lived painter invariant grep-verified. |
| 3 | `xvface_`/`xvfacetraits_` work unchanged via `MefistoPoint { short x; short y; }` byte-identical to `XPoint` (ROADMAP SC-3 / DRAW-03) | VERIFIED | `xvue_qt_api.cpp:425-443` (xvface_) and `:519-539` (xvfacetraits_) dereference `MefistoPoint*` into `QPolygon`; `xvue/qt/MEFISTO_POINT_AUDIT.md` certifies 29/29 live Fortran callers use `INTEGER*2 (2,N)` — static_assert at header lines 30-35 enforces 4-byte layout. |
| 4 | Pen style (`xvtypetrait_`), pen width (`xvepaisseur_`), `effacer_`/`xvvoir_`/`xvpxfenetre_`, and `Antialiasing` enabled by default (ROADMAP SC-4 / DRAW-06, DRAW-07, DRAW-08) | VERIFIED | `xvue_qt_api.cpp:446-467` pen writers store into `st->pen_style_`/`pen_width_base_` + `applyPen()`; `:303-316` effacer_ fillRect + flush epilogue; `:569-577` xvvoir_ = update+processEvents; `:396-404` xvpxfenetre_ returns canvas width/height. `xvue_qt_canvas.cpp:73` sets `QPainter::Antialiasing, true` after every painter begin (Pitfall 5 handled). |
| 5 | Resizing the canvas reallocates the backing pixmap preserving prior content per documented convention (ROADMAP SC-5 / DRAW-09) | VERIFIED | `xvue_qt_canvas.cpp:32-77` resizeEvent implements full D-07 sequence: end old painter → allocate DPR-aware new backing → scoped `QPainter tmp` fillRect(background) + drawPixmap(0,0,*old) → delete old → swap → begin painter on new backing → reapply Antialiasing + applyPen(). `xvue/qt/README_RESIZE.md` documents the top-left / no-scale / no-center convention. Human visual checkpoint Check 5 PASS 2026-04-11: drag-resize preserved geometry at top-left, new area filled black. |

**Score:** 5/5 ROADMAP success criteria verified; 9/9 DRAW requirements closed.

### Deferred Items

| # | Item | Addressed In | Evidence |
|---|------|--------------|----------|
| 1 | Full `xvtest1..4` visual parity vs X11 | Phase 3 | `02-CONTEXT.md` D-35, `02-RESEARCH.md:87`, `02-DISCUSSION-LOG.md:136` — blocked on text/colormap entry points |
| 2 | `testa/` 5-case A/B re-run | Phase 8 | `ROADMAP.md:23` owns A/B validation as Phase 8 |
| 3 | `xvcouleur_` / independent pen+brush colors | Phase 3 | `xvue_qt_state.h:23` TODO(phase 3) on `foreground_`; latent CR-01 / WR-01 are blocked by this |

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `xvue/qt/src/xvue_qt_state.h` | `backing_`, `painter_`, `foreground_`, `pen_`, `brush_`, `pen_style_`, `pen_width_base_`, `applyPen()`, `~XvueState()` | VERIFIED | All 9 declarations present (lines 19-41); `background_` remains first field (D-04 invariant); forward decls for `QPixmap`/`QPainter` used correctly. |
| `xvue/qt/src/xvue_qt_state.cpp` | `applyPen()` body + `~XvueState()` body | VERIFIED | `applyPen()` lines 9-32 implements 3-case switch (solid / dash / dash-double with max(1, 2×) width) and self-gates on `painter_->isActive()`; `~XvueState()` lines 34-44 ends painter if active, deletes painter and backing, nulls both. |
| `xvue/qt/src/xvue_qt_canvas.cpp` | `paintEvent` = drawPixmap blit, `resizeEvent` = D-07 sequence | VERIFIED | `paintEvent` lines 22-30 is exactly one drawPixmap call guarded by `state_ && state_->backing_`; `resizeEvent` lines 32-77 implements D-07 steps a-k in documented order. |
| `xvue/qt/src/xvue_qt_api.cpp` | 4 Wave-1 bodies (effacer_, xvvoir_, xvpxfenetre_, xvfond_ ext) + 13 Wave-2 bodies + 4 rect shims through helper | VERIFIED | All 17 real bodies present; 1 helper `xvue_qt_draw_rect_common` in anonymous namespace (lines 42-50); `drawPie` ≠ `drawArc` split correctly applied (lines 632, 656) with `× 16.0f` angle conversion (NOT ×64). |
| `xvue/qt/MEFISTO_POINT_AUDIT.md` | 29/29 callers INTEGER*2 | VERIFIED | Exists; contains "INTEGER*2", "ABI-safe"; zero INTEGER*4 offenders. |
| `xvue/qt/README_RESIZE.md` | Top-left anchor convention | VERIFIED | Exists; contains "top-left". |
| `prpr/xvtest0.f` | Extended DRAW-01..09 coverage section, SLEEP(15/10/3) holds | VERIFIED | Plan 02-01 Task 4 (`b556037`) added draw-coverage section; plan 02-04 Task 0 (`a82247b`) bumped SLEEPs to 15/10/3 for human checkpoint usability. |
| `.planning/phases/02-.../02-VALIDATION.md` | `status: complete`, `nyquist_compliant: true`, `wave_0_complete: true`, dated approval | VERIFIED | Frontmatter lines 4-8 all match (`complete` / `true` / `true`); closed commit `c84741b`. |
| `pp/ppxvtest0_qt` | Built, exit 0, zero DRAW warn-once lines | VERIFIED | Binary 136512 bytes rebuilt 20:29; fresh run during this verification: EXIT=0, zero `xvue-qt: stub xv…_` lines; reopen cycle clean. |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|----|--------|---------|
| `xvue_qt_canvas.cpp::resizeEvent` | `xvue_qt_state.h::backing_`/`painter_`/`applyPen` | `state_->backing_` reassignment + `painter_->begin/end` + `applyPen()` | WIRED | Lines 41-76: end → allocate → blit → swap → begin → setRenderHint → applyPen(). All state fields touched; single `painter_->begin` call site. |
| `xvue_qt_api.cpp::effacer_` | `xvue_qt_state.h::backing_`/`painter_`/`background_` | `painter_->fillRect(backing_->rect(), background_)` | WIRED | Line 310 — fillRect call present and guarded. |
| `xvue_qt_api.cpp::xvpxfenetre_` | `xvue_qt_canvas.h` | `canvas()->width()` / `height()` | WIRED | Lines 402-403 with logical-pixel comment referencing SHELL-06. |
| `xvue_qt_api.cpp::{xvtypetrait_, xvepaisseur_}` | `xvue_qt_state.h::applyPen` | `st->applyPen()` | WIRED | Lines 454, 466 — both call applyPen(); applyPen self-gates on `painter_->isActive()` so state persists across reopen (Wave 2 decision). |
| `xvue_qt_api.cpp::xvtrait_` / `xvftrait_` | `xvue_qt_state.h::painter_` | `st->painter_->drawLine` | WIRED | `xvtrait_` line 477; `xvftrait_` delegates via `proc(xvtrait)(x1,y1,x2,y2)` line 485. |
| `xvue_qt_api.cpp::xvarcellipse_` | `xvue_qt_state.h::painter_` | `st->painter_->drawPie` | WIRED | Line 656, `× 16.0f` angle conversion matches Qt 1/16° convention. |
| `xvue_qt_api.cpp::xvbordarcellipse_` | `xvue_qt_state.h::painter_` | `st->painter_->drawArc` | WIRED | Line 632, drawArc (outline) ≠ drawPie (xvarcellipse_). RESEARCH Q1 correction correctly applied. |
| `xvue_qt_api.cpp::{4 × xv*rectangle_}` | helper | `xvue_qt_draw_rect_common(*x,*y,*width,*height)` | WIRED (with latent CR-01 flaw, see Anti-Patterns) | All 4 entry points call the helper (lines 591, 598, 605, 612); helper calls `drawRect` line 47. See CR-01 — not a DRAW-04 goal failure under Phase 2's single-color regime but a documented latent bug. |

### Data-Flow Trace (Level 4)

Drawing primitives do not render "fetched data"; their inputs are Fortran-passed integers/floats, so the classic data-source check does not apply. The equivalent check here is: does an invocation of each primitive reach the backing pixmap?

| Primitive | Drives | Flows to backing | Status |
|-----------|--------|------------------|--------|
| xvtrait_ / xvftrait_ / xvtraits_ | `painter_->drawLine` / `drawPolyline` on `backing_` via long-lived painter | Yes | FLOWING |
| xvface_ / xvfacetraits_ | `painter_->drawPolygon` on `backing_` | Yes | FLOWING |
| 4 × rectangles | `painter_->drawRect` on `backing_` | Yes | FLOWING (goal-level); CR-01 flags latent color-divergence semantics |
| xvarcellipse_ / xvbordarcellipse_ | `painter_->drawPie` / `drawArc` on `backing_` | Yes | FLOWING |
| effacer_ | `painter_->fillRect(backing_->rect(), background_)` | Yes | FLOWING |
| xvfond_ | updates `background_` + fillRect on `backing_` | Yes | FLOWING |
| paintEvent | `QPainter(this).drawPixmap(0, 0, *backing_)` surfaces backing to widget | Yes | FLOWING |
| resizeEvent | allocate → fillRect(bg) → drawPixmap(0,0,*old) → swap | Yes | FLOWING |

Human visual checkpoint (`02-04-SUMMARY.md` table lines 113-120) confirms the full pipeline renders to the real display: Check 1 (geometry present) PASS, Check 3 (antialiasing smooth) PASS, Check 5 (resize preserves top-left) PASS. This is the strongest possible Level-4 evidence for an interactive graphics phase.

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| `pp/ppxvtest0_qt` runs, exits 0, no DRAW warn-once lines | `pp/ppxvtest0_qt` | EXIT=0, banner printed, two open/close cycles, zero `xvue-qt: stub xv…_` lines | PASS |
| DRAW-XX symbol distinctness (4 rect symbols survive, not collapsed in ABI) | `verify_abi` reported 57 symbols (per 02-04 SUMMARY, re-confirmed in recent builds) | 57/57 | PASS |
| Single `painter_->begin` call site (DRAW-01 invariant) | grep `painter_->begin xvue/qt/src` | 1 hit (canvas.cpp:70) | PASS |
| Zero forbidden PS-recording tokens survive | grep `lasopsc|courgb|counb|ypixels xvue_qt_api.cpp` | 0 hits | PASS |
| Angle convention enforced (×16, not ×64) | Plan 02-03 SUMMARY self-check: `grep '\* 64\|\*64'` | 0 hits | PASS |
| Build artifact present | `ls pp/ppxvtest0_qt` | 136512 bytes, rebuilt 20:29 local | PASS |
| Audit artifact present + ABI-safe | `ls xvue/qt/MEFISTO_POINT_AUDIT.md` + grep | OK | PASS |
| Resize convention doc present | `ls xvue/qt/README_RESIZE.md` + grep "top-left" | OK | PASS |

### Requirements Coverage

| Req | Plan(s) | Description | Status | Evidence |
|-----|---------|-------------|--------|----------|
| DRAW-01 | 02-02 | Single long-lived QPainter + paintEvent = drawPixmap blit | SATISFIED | `xvue_qt_state.h:27-28` + `xvue_qt_canvas.cpp:22-30`, `:70` (single begin site) |
| DRAW-02 | 02-03 | Lines / polylines via drawLine / drawPolyline | SATISFIED | `xvue_qt_api.cpp:470-516` (xvtrait/xvftrait/xvtraits) |
| DRAW-03 | 02-01, 02-03 | Polygons via drawPolygon + MefistoPoint ABI | SATISFIED | `xvue_qt_api.cpp:425-443` (xvface_) + `:519-539` (xvfacetraits_) + `MEFISTO_POINT_AUDIT.md` |
| DRAW-04 | 02-03 | 4 rectangle variants | SATISFIED (with latent CR-01) | `xvue_qt_api.cpp:587-613` — 4 distinct symbols route through `xvue_qt_draw_rect_common`. Goal-level behavior verified visually; CR-01 latent-bug noted below. |
| DRAW-05 | 02-01, 02-03 | Ellipse arcs via drawArc + drawPie, ×16 conversion, float* signatures | SATISFIED | `xvue_qt_api.cpp:618-636` (drawArc outline), `:642-660` (drawPie filled) — the RESEARCH Q1 correction is literally in the code with documented comments |
| DRAW-06 | 02-03 | Pen style + width via xvtypetrait_/xvepaisseur_ + applyPen | SATISFIED | `xvue_qt_api.cpp:446-467` + `xvue_qt_state.cpp:9-32` — 3-case switch with max(1, width×2) for type-2 |
| DRAW-07 | 02-02 | effacer_ / xvvoir_ / xvpxfenetre_ | SATISFIED | `xvue_qt_api.cpp:303-316` + `:569-577` + `:396-404` |
| DRAW-08 | 02-02 | Antialiasing default on the painter | SATISFIED | `xvue_qt_canvas.cpp:73` setRenderHint(Antialiasing, true) re-applied after every begin() (Pitfall 5) |
| DRAW-09 | 02-01, 02-02 | Resize preserves content top-left via documented convention | SATISFIED | `xvue_qt_canvas.cpp:32-77` D-07 sequence + `README_RESIZE.md` convention doc + human checkpoint Check 5 PASS |

No orphaned requirements: all 9 DRAW-XX IDs declared in PLAN frontmatter, traced to implementation, matched in REQUIREMENTS.md §Drawing, and cross-checked against ROADMAP Phase 2 scope.

### Anti-Patterns Found

| File | Line(s) | Pattern | Severity | Impact |
|------|---------|---------|----------|--------|
| `xvue/qt/src/xvue_qt_api.cpp` | 42-50, 587-613 | **CR-01**: 4 rect symbols collapse via `drawRect` which strokes+fills; legacy splits `XDrawRectangle` (outline-only) vs `XFillRectangle` (fill-only) | Warning (latent) | Invisible today because `applyPen()` sets both pen and brush from `foreground_ = Qt::white`; will produce wrong output once Phase 3 unlocks independent colors. Not a Phase 2 goal failure — ROADMAP SC-4 says "pen style / pen width / effacer / xvvoir / xvpxfenetre behave identically" which the helper does in the current single-color regime. Documented in 02-REVIEW.md with full fix sketch. Should be tracked for Phase 3 closure. |
| `xvue/qt/src/xvue_qt_api.cpp` | 439, 534-535 | **WR-01**: `drawPolygon` with default pen+brush strokes an implicit outline on `xvface_`; mirrors CR-01 | Warning (latent) | Same color-regime latency. Phase 3 fix. |
| `xvue/qt/src/xvue_qt_api.cpp` | 49, 315, 359, 392, 442, 479, 515, 538, 576, 635, 659 | **WR-02**: per-primitive `processEvents` enables re-entrant resize mid-batch | Warning | The single-long-lived-painter invariant remains intact because resize ends+rebegins the painter before the next primitive. No observed misbehavior. Phase 3 should consider deferring the flush to `xvvoir_`. |
| `xvue/qt/src/xvue_qt_api.cpp` | 453, 465, 477, 591-612, 627-651 | **WR-03**: primitive entry points dereference `int*` / `float*` without null guards (inconsistent with `xvpxfenetre_`) | Warning | Fortran callers never pass null in practice. Defensive-consistency fix only. |
| `xvue/qt/src/xvue_qt_canvas.cpp` | 41-73 | **WR-05**: resizeEvent is not strongly exception-safe — if `new QPixmap` or `new QPainter` throws after the old painter has been ended, state is inconsistent | Warning | `std::bad_alloc` on pixmap allocation is essentially fatal in this codebase; not a Phase 2 blocker. |
| `prpr/xvtest0.f` | 72, 91, 97 | **IN-02**: ~28 s cumulative sleep per run; no CI-friendly mode | Info | Intentional (plan 02-04 Task 0) for the human visual checkpoint. |

**None of these are Phase 2 goal failures.** CR-01 and WR-01 are latent color-regime bugs that are invisible until Phase 3 differentiates pen and brush colors. Phase 2 ships a **single-color** painter by design (`foreground_ = Qt::white` with TODO(phase 3) marker), so under its own scope every primitive renders correctly and visual parity with xvtest0.f is confirmed by the human operator.

### Human Verification Required

**None outstanding.** The single non-automated gate (Wave 3 visual checkpoint) was executed and approved 2026-04-11 by the human operator on a live Qt 6 display session: all 7 checks PASS (`02-04-SUMMARY.md` table lines 113-120 + Validation Sign-Off approval line). Checks 4/6/7 also covered headlessly by the automated clean-rebuild sweep. This report does NOT require a new human pass.

### Gaps Summary

No gaps. Every ROADMAP success criterion for Phase 2 is satisfied under its documented scope. Every DRAW-01..09 requirement has real code, real wiring, and real visible behavior confirmed both by a headless smoke run (`pp/ppxvtest0_qt` exit 0, zero warn-once lines) and by a dated human visual checkpoint. The single `painter_->begin` call site grep-verifies the DRAW-01 single-long-lived-painter invariant. The resize preserve convention is both documented (`README_RESIZE.md`) and enforced by the resizeEvent D-07 sequence. The arc `drawPie`/`drawArc` split and the ×16 (not ×64) angle conversion — the one non-mechanical port in the phase — are correctly applied with explanatory comments.

The latent CR-01 / WR-01 color-regime bugs identified in `02-REVIEW.md` are documented, don't affect Phase 2 goal achievement (single-color regime by design), and are naturally scoped to Phase 3 when `xvcouleur_` unlocks independent pen/brush colors. They should be tracked on the Phase 3 entry checklist.

Full `xvtest1..4` visual parity (the literal ROADMAP SC-1 text) is explicitly deferred to Phase 3 per `02-CONTEXT.md` D-35, ratified at planning time before execution, because those drivers depend on text/colormap entry points not yet implemented. Phase 2 substitutes `xvtest0.f` + human checkpoint as the equivalent gate under its own scope.

## Recommendation

**Close Phase 2.** Advance to Phase 3 (text, fonts, colormap) with the following entries on its follow-up checklist:
1. Fix CR-01 (rectangle fill vs outline semantics) once `xvcouleur_` differentiates pen/brush colors.
2. Fix WR-01 (`xvface_` implicit-outline stroke) as part of the same color-divergence work.
3. Consider deferring per-primitive `processEvents` to `xvvoir_` (WR-02) once the first batching-sensitive Fortran caller surfaces.
4. Optional hardening: WR-03 null guards, WR-05 exception-safety in resizeEvent.
5. Complete the literal ROADMAP SC-1 by running `prpr/xvtest1..4` against the Qt backend with Phase 3's text+palette in place (this is already listed in Phase 3's scope per D-35).

---

_Verified: 2026-04-11T20:48:00Z_
_Verifier: Claude (gsd-verifier)_
