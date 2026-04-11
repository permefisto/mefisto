# Phase 3: Text, fonts, colormap - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions captured in 03-CONTEXT.md — this log preserves the analysis.

**Date:** 2026-04-11
**Phase:** 03-text-fonts-colormap
**Mode:** discuss (interactive, standard)

## Gray Areas Presented

| # | Area | Key question |
|---|------|--------------|
| 1 | Bundled font strategy (TEXT-01) | Single TTF + size table vs small font bundle vs system+fallback |
| 2 | Palette size + storage location | `std::array<QColor, 256>` in XvueState vs singleton in XvueApp; survive open/close? |
| 3 | Dark-mode freeze mechanism (TEXT-06) | Raw-RGB-never-QPalette vs `changeEvent` intercept vs `WA_OpaquePaintEvent` |
| 4 | Text positioning semantics (TEXT-03) | Baseline `drawText(x, y, str)` vs `QRect` form |
| 5 | RGB storage precision | Preserve Fortran `float[256]` vs precomputed `QColor[256]` |

## User Responses

1. **Bundled font:** Option A (single TTF multiple sizes) with **DejaVu Sans Mono**.
2. **Palette:** Option A (`std::array<QColor, 256>` in `XvueState`) — **and palette must survive open/close**.
3. **Dark-mode freeze:** Option A — raw RGB, never QPalette.
4. **Text positioning:** Baseline form `drawText(x, y, str)` for literal X11 `XDrawString` parity.
5. **RGB precision:** Preserve Fortran float arrays.

## Claude Resolution Note

User picked (A) "palette lives in XvueState" AND "palette survives open/close". `XvueState` is owned by `XvueWindow` (Phase 1 D-05) and is destroyed on `xvfermer_`, so instance members would not persist.

**Resolution:** Palette arrays become **class-scope `static` members** of `XvueState` (D-09). This:
- Satisfies "in XvueState" (code locality preserved, grep finds it in the expected place)
- Satisfies "survive open/close" (process-lifetime storage)
- Matches legacy `xvue/xvuelc.c:100-104` byte-for-byte (file-scope static arrays)
- Requires no ownership move to `XvueApp` or singleton indirection

Flagged explicitly in D-09 so the planner cannot miss the distinction between instance and class storage.

## Scope-Creep Detected

None. All user input was within the TEXT-01..TEXT-06 phase boundary.

## Prior Context Applied

- Literal ABI parity rule from 00/D-33 — extended to the 8 new Phase 3 entry points
- Symbol-collapse pattern from 02/D-09 — extended to `xvtexte_`/`xvftexte_`
- Long-lived painter from 02/D-05 — text rendering uses the same painter, no new begin/end
- `foreground_ = Qt::white` placeholder from 02/D-21 — Phase 3 is where it gets retired
- `xvfond_` minimal 2-entry hack from 02/D-14 — Phase 3 is where it gets retired
- Antialiasing reassertion from 02/D-22 — Phase 3 adds `TextAntialiasing` as a one-line sibling
- PS-state drop rule from 02/D-26 — applied to `xvchargefonte_` / `xvtexte_` / color entry points
- Null-guard symmetry rule (from Phase 2 WR-03 fix, commit `b9b37d6`) — extended to all Phase 3 entries
- Phase 2 deferred A/B gate against `xvtest1..4` + `testa/` — Phase 3 is where it closes (D-26)
- Phase 2 latent CR-01/WR-01/WR-02 fixes (commit `6741c68` / `fe52db4`) — implicitly verified once Phase 3 introduces independent pen/brush colors

## Codebase Scout Findings (applied in context)

- `xvue/xvuelc.c:85` — `CMAPSIZE 256` → D-10 pinning
- `xvue/xvuelc.c:100-104` — file-scope `static float red/green/blue[CMAPSIZE]` + `norgb[CMAPSIZE]` → D-09 / D-12 / D-13 storage mirror
- `xvue/xvuelc.c:131-132` — `nbpolices` + `listfonts` → D-03 / D-22 font metadata mirror
- `xvue/xvuelc.c:1463-1580` — legacy `xvchargefonte_` uses `XLoadQueryFont` + dynamic `XListFonts` → **not reproducible**, Phase 3 replaces with bundled-TTF path (D-01/D-02/D-04)
- `xvue/xvuelc.c:1650-1700` — `xvtexte_` / `xvftexte_` use `XDrawString(display, ?, gc, x, y, …)` → baseline form required (D-06 / D-08)
- `xvue/qt/src/xvue_qt_api.cpp:189-223` — current zero-fill `xvinfo_` → D-22/D-23 real values
- `xvue/qt/src/xvue_qt_api.cpp:337-358` — Phase 2 `xvfond_` 2-entry hack → retired in Phase 3
- `xvue/qt/src/xvue_qt_api.cpp:380-420` — current text/font stubs to replace

## No Follow-Up Rounds Needed

All 5 gray areas resolved in a single round. No correction loop.
