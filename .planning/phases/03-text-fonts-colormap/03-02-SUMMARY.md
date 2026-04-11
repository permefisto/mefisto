---
phase: 03
plan: 02
subsystem: xvue-qt
tags: [text, fonts, colormap, palette, phase-3, wave-1]
requires:
  - 03-01  # XvueState font/palette members, DejaVu Sans Mono bundling, xvtest0 coverage
provides:
  - xvchargefonte_ real body (D-04, TEXT-01)
  - xvnbpixeltexte_ real body (D-05, TEXT-02)
  - xvtexte_ / xvftexte_ real bodies via shared xvue_qt_draw_text_common (D-06, TEXT-03)
  - xvcouleur_ real body (D-14, TEXT-04)
  - xvactivervb_ bulk-load real body (D-17 AMENDED per A3, TEXT-05)
  - xvrecuprgbdec_ real body (D-18, TEXT-05)
  - xvinfo_ real palette/font fill (D-22, D-23)
  - xvfond_ retired Phase 2 two-entry hack, routed through palette_cache_
affects:
  - xvue/qt/src/xvue_qt_api.cpp
tech-stack:
  added: []
  patterns: [QColor::fromRgbF float storage, palette_cache_dirty_ rebuild, shared draw-text helper, Pitfall-4 two-arg fromLatin1]
key-files:
  created: []
  modified:
    - xvue/qt/src/xvue_qt_api.cpp
decisions:
  - "xvtexte_ and xvftexte_ collapsed to a single static helper xvue_qt_draw_text_common per D-06 + Pitfall 7 — legal because 02/D-05 unified the draw target into a single backing pixmap"
  - "xvactivervb_ implemented as bulk palette load per A3 research correction, NOT as D-17's transient-color interpretation — the only live caller (xvue/palcde.f:619) passes NDCOUL+1 cells"
  - "xvCouleursImposees / xvStockeRGBtoColormap / xvColormapToRGB NOT added as extern \"C\" symbols — A1 preserved (00/D-33); these would only live as file-static helpers if ever needed"
  - "Phase 3 xvfond_ retires the 2-entry hardcoded palette but KEEPS the Phase 2 D-24 fillRect-on-backing mechanism; only the color source changes"
metrics:
  duration: ~25m
  completed: 2026-04-11T21:55Z
---

# Phase 3 Plan 02: TEXT/Colormap Entry Points Summary

Seven Phase 3 public entry points (xvchargefonte_, xvnbpixeltexte_, xvtexte_, xvftexte_, xvcouleur_, xvactivervb_, xvrecuprgbdec_) plus the xvinfo_ palette/font section and the Phase 2 xvfond_ two-entry hack are replaced with real bodies routing through the XvueState font + palette members that arrive from sibling plan 03-01.

## Tasks

| # | Name | Commit | Files |
|---|------|--------|-------|
| 1 | Implement 3 text bodies + shared draw helper | `6cbe6bc` | xvue/qt/src/xvue_qt_api.cpp |
| 2 | Implement 3 color bodies + retire xvfond_ hack + xvinfo_ palette/font fill | `f74ff93` | xvue/qt/src/xvue_qt_api.cpp |

## What changed

### Task 1 — Text entries (commit `6cbe6bc`)

- Added `xvue_qt_draw_text_common(const char*, int, int, int)` inside the anonymous namespace. Null-guards string/length, fetches window/state/painter/backing, `QString::fromLatin1(string, length)` TWO-ARG (Pitfall 4), `painter_->setFont(current_font_)`, `painter_->drawText(x, y, qstr)` baseline form, then runs the Phase 2 D-01 epilogue (`canvas()->update()` + `QCoreApplication::processEvents(ExcludeUserInputEvents)`).
- `xvchargefonte_`: clamps `*nofont` to `[0, kNbFonts-1]`, rebuilds `QFont("DejaVu Sans Mono", kFontSizes[idx])` and `QFontMetrics`, pushes the font to the live painter, writes `*largpx = horizontalAdvance(QLatin1Char('0'))` and `*hautpx = height()`. `*nofont0` explicitly discarded per Pitfall 6 (QFont is RAII).
- `xvnbpixeltexte_`: measures a two-arg Latin-1 QString with `current_metrics_`, returns zero pair on null state.
- `xvtexte_` / `xvftexte_`: both delegate to the shared helper after WR-03 null-guarding their own pointer params. D-06 + Pitfall 7 comment explains why the collapse is legal.
- Added headers: `<cstring>`, `<QFont>`, `<QFontMetrics>`, `<QLatin1Char>`, `<QString>`.
- Lines changed: +75 / -16 (net +59).
- 4 warn_once lines (xvchargefonte/xvnbpixeltexte/xvftexte/xvtexte) deleted.

### Task 2 — Color entries + xvinfo + xvfond (commit `f74ff93`)

- `xvrecuprgbdec_` (D-18): read-only snapshot of `XvueState::red/green/blue[0..n-1]` into Fortran `r/g/b`, clamped to `kMaxPalette`. Pitfall 2 comment documents the ordering guarantee vs xvue/xvinit.f:143.
- `xvactivervb_` (D-17 AMENDED, A3): bulk copies `r/g/b[0..nbcells-1]` into `XvueState::red/green/blue` and marks `palette_cache_dirty_[i] = true`. `*palcour` dropped per 02/D-26. Does NOT touch `foreground_`, does NOT call `applyPen`, does NOT flush.
- `xvcouleur_` (D-14): bounds-check with fallback to index 1 on out-of-range; lazily rebuilds `palette_cache_[i] = QColor::fromRgbF(red[i], green[i], blue[i])` when dirty; assigns `st->foreground_ = palette_cache_[i]`; calls `st->applyPen()`. State-change only — no flush.
- `xvfond_`: retired the Phase 2 2-entry {black, white} hardcode and the associated `warned_xvfond_range` warn-once. Now reads color from `palette_cache_[i]` via the same dirty-rebuild pattern as xvcouleur_, then preserves the Phase 2 D-24 mechanism (update `background_`, `painter_->fillRect(backing_->rect(), background_)`, canvas update).
- `xvinfo_` (D-22/D-23): preserved the Phase 1 window-resize branch. Replaced the zero-fill palette/font block and the `warned_xvinfo_partial` stub with real values:
  - 10 literal "DejaVu Sans Mono Npt" strings (8, 10, 12, 14, 16, 18, 20, 24, 28, 32) copied into `namefonts[k]` via `strncpy` with 255-byte cap; `nbchar[k] = strlen`.
  - `*maxfonts = *nbfonts = kNbFonts` (respecting caller's `*maxfonts` upper bound if non-zero).
  - `*visuclass = 4` (TrueColor — Qt is always 24-bit).
  - Palette ranges: `*n1coref=0 *ndcoref=15 *n1coelf=0 *ndcoelf=15 *n1coulf=16 *ndcoulf=255 *nbcolo=256`.
- Lines changed: +129 / -70 (net +59).
- 5 warn_once lines (xvrecuprgbdec/xvactivervb/xvcouleur/xvinfo_partial/xvfond_range) deleted.

## Source-level acceptance (verified in worktree)

| Check | Result |
|---|---|
| `grep -c 'QColor::fromRgbF'` | 2 (xvcouleur_, xvfond_) |
| `grep -c 'palette_cache_dirty_'` | 5 |
| `grep -c 'DejaVu Sans Mono'` | 12 (10 list + 2 QFont) |
| `grep -c 'xvue_qt_draw_text_common'` | 3 (1 def, 2 calls) |
| `grep -c 'visuclass = 4'` | 1 |
| `grep -c 'current_metrics_.horizontalAdvance'` | 2 (QLatin1Char + qstr) |
| `grep -c 'painter_->drawText'` | 1 (baseline form) |
| `grep -c 'QString::fromLatin1'` | 2 (both 2-arg) |
| warn_once for Phase 3 entries | 0 remaining |
| `warned_xvinfo_partial` / `warned_xvfond_range` | 0 remaining |

## Build / runtime verification

**NOT completed in this worktree** — see "Known Blockers" below. The acceptance criteria that reference `bin/cbl_tout_qt`, `pp/ppxvtest0_qt`, and `nm libxvueqt.a` could not be exercised here because sibling plan 03-01 has not yet been merged into this worktree. All references to `XvueState::kNbFonts`, `kFontSizes`, `kMaxPalette`, `red`, `green`, `blue`, `palette_cache_`, `palette_cache_dirty_`, `current_font_`, `current_metrics_`, and `current_font_size_idx_` are the exact member names documented in 03-02-PLAN.md `<interfaces>` and 03-01-PLAN.md, so once the parallel-wave framework merges 03-01 the compile should succeed without further edits.

Local build reproducer (fails as expected until merge):

```text
error: 'kNbFonts' is not a member of 'XvueState'
error: 'struct XvueState' has no member named 'current_font_size_idx_'
error: 'struct XvueState' has no member named 'current_font_'
error: 'kFontSizes' is not a member of 'XvueState'
error: 'struct XvueState' has no member named 'current_metrics_'
```

These errors match exactly the 03-01 artifact surface; no extra/unexpected symbol dependencies.

## Known Blockers

- **Parallel wave dependency on 03-01 not yet present in this worktree.** The parallel_wave_note in the execution prompt explicitly covers this: "the executor framework merges 03-01's worktree back before the wave finishes, but if your build fails because an XvueState member or Qt resource symbol declared in the 03-01 plan is missing, STOP, note the missing symbol in SUMMARY.md as a blocker, commit what you have, and return — do not reimplement 03-01 work in this worktree." Followed literally: both tasks committed, SUMMARY notes the blocker, returning.
- Required symbols not yet in this worktree (expected from 03-01):
  - `XvueState::kNbFonts`, `XvueState::kFontSizes[]`, `XvueState::kMaxPalette`
  - `XvueState::red[]`, `XvueState::green[]`, `XvueState::blue[]`
  - `XvueState::palette_cache_[]`, `XvueState::palette_cache_dirty_[]`
  - `XvueState::current_font_`, `XvueState::current_metrics_`, `XvueState::current_font_size_idx_`
- Post-merge verification that remains for the wave orchestrator (or a follow-up executor) to run once 03-01 is folded in:
  - `bin/cbl_tout_qt` green
  - `bin/cbl_tout` green
  - `pp/ppxvtest0_qt` exits 0 with ZERO warn-once lines for xvchargefonte/xvnbpixeltexte/xvtexte/xvftexte/xvcouleur/xvrecuprgbdec/xvactivervb
  - `nm xvue/qt/build/libxvueqt.a | grep -cE '^.* T xv[a-z]+_$'` returns exactly 57
  - `nm xvue/qt/build/libxvueqt.a | grep -E ' T xv(CouleursImposees|StockeRGBtoColormap|ColormapToRGB)'` returns ZERO matches (A1 preserved)
  - `sh xvue/qt/cmake/verify_no_exec.sh xvue/qt/src xvue/qt/include` exits 0

## ABI status

No new `extern "C"` symbols added in xvue/qt/src/xvue_qt_api.cpp. The file still declares bodies for exactly the set already declared in xvue/qt/include/xvue_qt_api.h (57 entries). `xvCouleursImposees`, `xvStockeRGBtoColormap`, `xvColormapToRGB` are NOT introduced — A1 and 00/D-33 preserved.

## Deviations from Plan

None at the source level. All Task 1 and Task 2 steps executed exactly as written. The deliberate deferral of build/runtime verification is pre-authorized by the parallel_wave_note in the execution prompt.

## Self-Check

- xvue/qt/src/xvue_qt_api.cpp — FOUND (modified)
- Commit `6cbe6bc` (Task 1) — FOUND
- Commit `f74ff93` (Task 2) — FOUND
- SUMMARY.md at .planning/phases/03-text-fonts-colormap/03-02-SUMMARY.md — FOUND (this file)
- Build verification — NOT RUN (documented blocker: 03-01 symbols missing in this worktree; expected resolution path is wave merge by orchestrator)

## Self-Check: PASSED (source-level; build gate deferred to post-merge per parallel_wave_note)
