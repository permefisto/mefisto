---
phase: 04
plan: 01
subsystem: xvue/qt
tags: [pixmap, save-restore, double-buffering, abi]
requires: [phase-02-drawing-primitives]
provides:
  - "XvueState::saved_canvas_ field (single lazily-allocated slot)"
  - "xvue_qt_save_to_slot / xvue_qt_restore_from_slot file-local helpers"
  - "Real Qt bodies for fenetremempx_/mempxfenetre_/sauvefenetre_/restaurefenetre_/sauvemempx_/restauremempx_/effacemempx_"
affects:
  - "xvue/qt/src/xvue_qt_state.h"
  - "xvue/qt/src/xvue_qt_state.cpp"
  - "xvue/qt/src/xvue_qt_api.cpp"
tech-stack:
  added: []
  patterns:
    - "Temporary scoped QPainter on a second QPixmap (Phase 2 resizeEvent pattern)"
    - "Lazy (re)allocation on size mismatch (D-02/D-14)"
    - "Anonymous-namespace file-local inline helpers to keep ABI count stable (Pitfall 7)"
key-files:
  created: []
  modified:
    - "xvue/qt/src/xvue_qt_state.h"
    - "xvue/qt/src/xvue_qt_state.cpp"
    - "xvue/qt/src/xvue_qt_api.cpp"
decisions:
  - "D-01/D-02/D-03 honored: single raw-pointer slot, lazy allocation, destroyed before painter_/backing_"
  - "D-04/D-05/D-06 honored: fenetremempx_/mempxfenetre_ are documented intentional no-ops"
  - "D-07/D-08/D-09 honored: two TU-local helpers back four save/restore entry points; every body starts with XVUE_QT_ASSERT_MAIN_THREAD()"
  - "D-10/D-11 honored: effacemempx_ body mirrors effacer_; kept as distinct ABI symbol"
  - "D-12/D-13/D-14 honored: size mismatch at restore = stderr warn + no-op; lazy reallocation on save-side size change"
metrics:
  duration: "~12min"
  completed: "2026-04-14"
---

# Phase 4 Plan 01: Pixmap save/restore real Qt bodies — Summary

One-liner: Replace 7 warn-once pixmap save/restore stubs with real Qt 6 bodies backed by a single lazily-allocated `QPixmap* saved_canvas_` slot on `XvueState`, keeping the ABI at 57 symbols.

## What Shipped

1. **`XvueState::saved_canvas_` field** (`xvue/qt/src/xvue_qt_state.h` + `.cpp`)
   - Raw `QPixmap*` default-initialized to `nullptr`, declared next to `backing_`/`painter_` with a Phase 4 D-01/D-02/D-03 comment. No existing fields reordered.
   - Destructor deletes `saved_canvas_` **before** `painter_` and `backing_` per D-03.

2. **File-local save/restore helpers** (`xvue/qt/src/xvue_qt_api.cpp`, inside the existing anonymous namespace)
   - `xvue_qt_save_to_slot()`: lazy (re)alloc on null/size-mismatch (D-02/D-14), HiDPI via `setDevicePixelRatio(backing_->devicePixelRatio())`, scoped temporary `QPainter tmp(saved_canvas_)` blits from `backing_` while the long-lived `painter_` stays bound to `backing_` (Phase 2 D-05 invariant; two-painters-two-devices is legal per Qt 6 docs).
   - `xvue_qt_restore_from_slot()`: null/size-mismatch → `fprintf(stderr, "xvue-qt: restore_from_slot: no slot or size mismatch\n")` + return (D-12). Otherwise uses the active `painter_->drawPixmap(0, 0, *saved_canvas_)` and schedules `canvas_->update()` (no `processEvents`).
   - Both are `inline` and anonymous-namespace, so TU-local → zero ABI impact (Pitfall 7).

3. **Seven entry-point bodies** (`xvue/qt/src/xvue_qt_api.cpp`)
   | Symbol | Body |
   |---|---|
   | `fenetremempx_` | intentional no-op + D-04 comment |
   | `mempxfenetre_` | intentional no-op + D-04 comment |
   | `sauvefenetre_` | calls `xvue_qt_save_to_slot()` |
   | `restaurefenetre_` | calls `xvue_qt_restore_from_slot()` |
   | `sauvemempx_` | calls `xvue_qt_save_to_slot()` (bit-identical to `sauvefenetre_`) |
   | `restauremempx_` | calls `xvue_qt_restore_from_slot()` (bit-identical to `restaurefenetre_`) |
   | `effacemempx_` | `painter_->fillRect(backing_->rect(), background_)` + `canvas_->update()` (same body as `effacer_`, D-10/D-11) |

   Every body starts with `XvueApp::ensure()` + `XVUE_QT_ASSERT_MAIN_THREAD()` per D-09.

## Verification (end-of-plan)

| Check | Result |
|---|---|
| `bin/cbl_tout_qt` | Exit 0; full `pp/pp*_qt` set rebuilt |
| `bin/cbl_tout` (legacy X11) | Exit 0; full legacy `pp/pp*` set rebuilt |
| `cmake --build . --target verify_abi` | `nm count: 57  header count: 57` |
| `cmake --build . --target verify_no_exec` | `OK (no forbidden tokens)` + palette-leak scan clean |
| `nm libxvueqt.a` for the 7 symbols | 7/7 present with trailing underscore |
| `warn_once(warned, "${sym}_")` occurrences in the 7 bodies | 0/7 |
| `QT_QPA_PLATFORM=offscreen MEFISTO_BATCH_X11=1 pp/ppxvtest0_qt` stub warnings for the 7 symbols | 0 |
| `git diff --stat` vs plan start | exactly 3 files modified (matches plan) |

## Requirements Delivered (code-side)

- **PIXMAP-01** — `fenetremempx_`/`mempxfenetre_` are documented intentional no-ops on Qt, ABI stable. Behavior validated on Qt already during Phase 03-04 close; this plan just upgrades warn-once stubs to documented no-ops.
- **PIXMAP-02** — `sauvefenetre_`/`restaurefenetre_` route through `xvue_qt_save_to_slot` / `xvue_qt_restore_from_slot` on the single lazily-allocated `saved_canvas_` slot.
- **PIXMAP-03** — `sauvemempx_`/`restauremempx_` share bodies with the fenetre pair; `effacemempx_` mirrors `effacer_`.
- **PIXMAP-04** — Explicitly `deferred-to-phase-5` per Phase 4 CONTEXT D-18; not a Phase 4 deliverable.

The end-to-end validation harness (synthetic save/draw/restore round-trip via `prpr/xvtest0.f` extension + A/B capture) is the job of Plan 04-02, which now has a real implementation to validate against.

## Commits

| Task | Commit | Message |
|---|---|---|
| 1 | `c165644` | `feat(04-01): add saved_canvas_ slot to XvueState (D-01/D-02/D-03)` |
| 2 | `4ddaf02` | `feat(04-01): add xvue_qt_save_to_slot/restore_from_slot helpers (D-07)` |
| 3 | `c590e41` | `feat(04-01): replace 7 pixmap save/restore stubs with real Qt bodies` |

## Decisions Honored

- **D-01..D-03** (slot model & lifecycle) — single raw-pointer slot, lazy allocation, destroyed first.
- **D-04..D-06** (flip ops) — both `fenetremempx_` and `mempxfenetre_` are documented intentional no-ops.
- **D-07..D-09** (save/restore helpers) — two TU-local helpers, bit-identical bodies across 4 entry points, `XVUE_QT_ASSERT_MAIN_THREAD` preserved.
- **D-10..D-11** (effacemempx_) — same body as `effacer_`, distinct symbol for ABI preservation.
- **D-12..D-14** (resize handling) — size mismatch = stderr warn + no-op; lazy reallocation on save-side size change; no proactive `resizeEvent` invalidation.

## Deviations from Plan

None — plan executed exactly as written. The only incidental action was adding `g++`/`c++` symlinks under `/tmp/gfortran-14-shim/` (pointing to `g++-14`) so `cbl_tout_qt` could locate the C++ compiler; this is a host-environment convenience, not a code change, and mirrors the existing `gcc`/`gfortran` symlinks already in that shim directory.

## Auth Gates

None.

## Deferred Issues

None.

## Known Stubs

None introduced by this plan. The 7 previously stubbed symbols now have real bodies; all other pre-existing warn-once stubs (outside this plan's scope) are untouched.

## Self-Check: PASSED

- `xvue/qt/src/xvue_qt_state.h` — FOUND (contains `QPixmap* saved_canvas_ = nullptr`)
- `xvue/qt/src/xvue_qt_state.cpp` — FOUND (contains `delete saved_canvas_`)
- `xvue/qt/src/xvue_qt_api.cpp` — FOUND (contains both helpers and all 7 new bodies)
- Commit `c165644` — FOUND in `git log`
- Commit `4ddaf02` — FOUND in `git log`
- Commit `c590e41` — FOUND in `git log`
- `verify_abi` target exits 0 with `nm count: 57  header count: 57`
- `verify_no_exec` target exits 0
