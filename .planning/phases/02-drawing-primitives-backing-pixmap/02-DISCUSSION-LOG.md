# Phase 2: Drawing primitives & backing pixmap - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the reasoning.

**Date:** 2026-04-11
**Phase:** 02-drawing-primitives-backing-pixmap
**Mode:** discuss (claude-driven, user deferred to Claude's judgment)
**Areas discussed:** Paint flush strategy, `xvftrait_`/`xvtrait_` routing, Pen style #2 mapping, Resize preserve convention, `effacer_` semantics, PS side-effect hooks, Phase 2 test driver, Baseline `testa/` A/B resume, `MefistoPoint` header placement, plus Qt-6-idiomatic bonus decisions (HiDPI backing, foreground color bridge, `QPainter` resize lifecycle).

## User directive

> "I don't really know, make it the more Qt compatible way and modern. think thoroughly test will be deferred"

Claude interpreted this as: (a) resolve every gray area toward the modern Qt-6-idiomatic approach when there's no strong legacy reason otherwise, (b) defer the test-driver and baseline `testa/` A/B decisions, (c) think thoroughly before writing decisions (no shallow defaults).

## Areas resolved

### Paint flush strategy

**Options considered:**
| Option | Description | Selected |
|--------|-------------|----------|
| Per-primitive pump | Every draw entry point ends with `canvas_->update() + processEvents(ExcludeUserInputEvents)` | ✓ |
| `xvvoir_`-only pump | Primitives just mutate the backing; flush deferred until the next `xvvoir_` / `effacer_` call | |
| Hybrid | Pump for "visible" primitives but not state-setters | |

**Rationale:** Matches Pitfall 6's HIGH-confidence recommendation. Legacy `xvtest*.f` callers have most `xvvoir_` calls commented out; a pump-per-primitive strategy keeps them visually correct without requiring a Fortran-side change. `ExcludeUserInputEvents` preserves the Phase 1 D-01 invariant (no input can reach Fortran through the pump). The overhead is measurable but small (one event-loop drain per primitive) and Phase 2 is not in a hot-path performance regime.

### `xvftrait_` vs `xvtrait_` routing

**Options considered:**
| Option | Description | Selected |
|--------|-------------|----------|
| Both to backing | Semantically identical bodies; legacy distinction obsolete | ✓ |
| `xvftrait_` to a separate overlay QPixmap | Preserve "immediate, not backed" semantic via a second surface | |
| `xvftrait_` via `QPainter(canvas_)` direct | Draw directly on the widget, bypassing the backing | |

**Rationale:** The legacy distinction existed because X11 had two surfaces (window + mempx) and some call sites wanted immediate-only drawing that wouldn't survive `fenetremempx`. Under the single-backing Qt model there is no second surface; `paintEvent` always blits the backing; any draw that doesn't touch the backing would vanish on the next repaint. Making the two symbols semantically identical is the only internally consistent option. Keeping both symbols (not merging) preserves the Fortran ABI.

### Pen style #2 mapping

**Options considered:**
| Option | Description | Selected |
|--------|-------------|----------|
| Preserve X11 semantic | `Qt::DashLine` with effective width = `base * 2` | ✓ |
| Use `Qt::DashDotLine` | Map 0/1/2 → Solid/Dash/DashDot, no width doubling | |
| Use `Qt::CustomDashPattern` | Hand-craft a dash pattern that matches X11 LineDoubleDash visually | |

**Rationale:** The legacy comment at `xvuelc.c:1769-1770` says explicitly "LIGNE TIRETEE D'EPAISSEUR DOUBLE" — callers relying on style 2 expect a visually heavier line. Mapping to `Qt::DashDotLine` would change the visual weight and the dash pattern, producing a two-way drift from X11. Width-doubling under the hood preserves the visual contract and keeps the Qt pen-style enum usage clean.

### Resize preserve convention

**Options considered:**
| Option | Description | Selected |
|--------|-------------|----------|
| Top-left sub-blit | Old content anchored at (0, 0); new area filled with background | ✓ |
| Center anchor | Old content centered in new backing; margins filled with background | |
| Scale to fit | `QPainter::drawPixmap` with source/target rect stretching | |
| Clear to background | Discard old content entirely on resize | |

**Rationale:** DRAW-09 says "optional preserve via sub-blit (documented convention)". Top-left is the convention every GUI toolkit defaults to; it maps directly to `QPainter::drawPixmap(0, 0, old)` with no arithmetic; it matches Fortran callers' coordinate intuition (origin at top-left, growing right and down). Centering introduces integer division edge cases (odd delta); scaling distorts the content and breaks the pixel-exact semantics scientific plotting relies on; clearing loses information.

### `effacer_` semantics

**Options considered:**
| Option | Description | Selected |
|--------|-------------|----------|
| Backing fillRect + update + pump | `painter_->fillRect(backing->rect(), background_)` + D-01 epilogue | ✓ |
| Preserve legacy dual-surface copy | Also call the (Phase 4) pixmap-slot copy mechanism | |
| Clear without pump | Backing fillRect only; rely on caller to flush | |

**Rationale:** The legacy dual-surface coordination exists because X11 `XClearWindow` clears the window surface but not `mempx`, so the legacy body also does `fenetremempx` to sync them. Under single-backing Qt the two surfaces are the same surface; the copy is definitionally unnecessary. Following the D-01 epilogue policy for symmetry makes the entry point behave like every other "visible mutation" entry point.

### PostScript side-effect hooks

**Options considered:**
| Option | Description | Selected |
|--------|-------------|----------|
| Strip entirely | No `lasopsc`-style hooks in Qt primitives; PS export is a Phase 7 clean-slate job | ✓ |
| Preserve the hook structure | Port the `if (lasopsc > 0)` blocks 1:1, writing to a Qt-side `fpo`/`concat` | |
| Hybrid | Keep a per-primitive callback slot for Phase 7 to attach a recorder | |

**Rationale:** The legacy pattern couples drawing primitives to PS recording via shared global state (`lasopsc`, `courgb`, `counb`, `fpo`, `concat[]`, `ypixels`, etc.). Porting that to the Qt backend would recreate 500+ lines of state-machine code whose only purpose is a feature that Phase 7 is going to rewrite anyway. Phase 7 can use `QPdfWriter` or a direct-emit PS writer with its own `QPainter` — no need to share state with the interactive drawing path. Stripping the hooks is a load-bearing decision for Phase 7's design freedom (D-26/D-27/D-28).

### Phase 2 test driver

**Options considered:**
| Option | Description | Selected |
|--------|-------------|----------|
| Defer exact form to planner, lock coverage minimum | D-36: planner picks vehicle, D-36 lists required coverage | ✓ |
| Extend `xvtest0.f` directly | Single-file change to Phase 1's shell driver | |
| New `xvtest0_draw.f` | Sibling driver in `prpr/` | |
| C++ unit test | Harness under `xvue/qt/tests/` | |

**Rationale:** User directive: "test will be deferred". Claude honored this by not locking the vehicle. D-36 constrains the *coverage* (primitives that must be exercised) as the Phase 2 exit gate, without dictating where the test lives. The three vehicles are all viable and the decision is a planner call under Claude's Discretion.

### Baseline `testa/` A/B resume

**Options considered:**
| Option | Description | Selected |
|--------|-------------|----------|
| Defer to Phase 3 | All 5 baseline cases need text + colormap; partial run is uninterpretable | ✓ |
| Run partial now | Accept that text/color paths produce warn-once output in Phase 2 | |
| Run a subset without text | Pick cases that don't exercise `xvtexte_` / `xvcouleurs_` (likely none) | |

**Rationale:** Phase 1 D-21 noted the A/B resume would land "at Phase 2", but closer inspection reveals every `testa/` baseline case calls `xvtexte_` and `xvcouleurs_` at least once. Phase 2 leaves those on the warn-once path; running the A/B would produce output with missing labels and wrong colors and the diff would be dominated by the expected warnings, not by real drawing regressions. Phase 3 (text + colormap) is the natural earliest point where the A/B is interpretable.

### `MefistoPoint` header placement & Fortran call-site audit

**Options considered:**
| Option | Description | Selected |
|--------|-------------|----------|
| In `xvue_qt_api.h` + mandatory audit deliverable | ABI-visible struct + `MEFISTO_POINT_AUDIT.md` grep artifact | ✓ |
| Internal header | Struct hidden in `xvue/qt/src/`, audit optional | |
| Reuse Xlib `XPoint` via `#include <X11/Xlib.h>` | Avoids declaring a new type | |

**Rationale:** The struct is part of the `extern "C"` ABI (the Fortran-side callers pass arrays that match this layout) so it belongs in the single Fortran-facing header alongside the function prototypes (Phase 0 D-04). Reusing Xlib `XPoint` would drag `libX11` back into the Qt backend, which the xvue-qt project exists to eliminate. Making the audit a first-class deliverable (not just a planner checklist item) turns Pitfall 4 from a one-off concern into a persistent, re-runnable artifact (D-31).

## Claude-bonus decisions (Qt-6-idiomatic)

Decisions Claude made beyond the presented gray areas, applying the "modern Qt-compatible" directive:

- **HiDPI backing** — Backing pixmap sized at `widget.size() * devicePixelRatioF()` with `setDevicePixelRatio(dpr)`; `paintEvent` uses `drawPixmap(0, 0, *backing)` (Qt 6 DPR-aware idiom). D-06.
- **Foreground color bridge** — Phase 2 hardcodes `XvueState::foreground_ = Qt::white` to avoid leaking palette state across the Phase 2/3 boundary; pen and brush both read from it; Phase 3 takes ownership. D-21.
- **`applyPen()` helper** — Single-source-of-truth for pen state rebuilding. Mutations to `pen_style_` or `pen_width_base_` both flow through it. D-16, D-17, D-19.
- **`QPainter` resize lifecycle** — `resizeEvent` is the only place `painter_->end()` and `begin()` happen; D-07 specifies the exact 8-step sequence to prevent partial-state bugs.
- **Polygon fill rule** — `Qt::OddEvenFill` matches Xlib `Complex+CoordModeOrigin` for the geometries MEFISTO produces; `WindingFill` alternative noted in deferred ideas if an A/B case surfaces. D-11.
- **Arc angle convention** — D-14 flags the conversion as a "literal diff against `xvuelc.c`" planner task so that whatever angular domain the legacy code expects at the Fortran ABI is preserved.
- **Cosmetic pen** — Qt width=0 "cosmetic pen (1 device pixel)" matches X11 `line_width=0` "thin line"; no special case. D-18.

## Deferred Ideas

Captured in `02-CONTEXT.md` `<deferred>` section. Summary:

- Full `xvtest1–4` parity + 5-case `testa/` A/B → Phase 3
- `ncf` / `nca` color indices in `xvfacetraits_` → Phase 3
- Border-color parameter in `xvfrectangle_` family → Phase 3
- Color-coded output generally → Phase 3
- PostScript export → Phase 7
- Pixmap save/restore slots → Phase 4
- Mouse/keyboard/event delivery → Phase 5
- `Qt::WindingFill` alternative → Phase 3 A/B gate
- Alternate arc-angle conversion → planner diff
- Multi-pixmap compositing → Phase 4+
- C++ unit-test framework introduction → planner discretion
- `xvue/qt/README_RESIZE.md` location → planner discretion

## Reviewed Todos (not folded)

None — no pending todos matched Phase 2 scope at init time.
