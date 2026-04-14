# Phase 4: Pixmap save/restore (double-buffering) — Research

**Researched:** 2026-04-13
**Domain:** Qt 6 QPainter/QPixmap lifecycle, off-screen slot semantics, round-trip pixel validation
**Confidence:** HIGH (Qt 6 semantics verified against doc.qt.io; codebase facts verified by direct read)

## Summary

Phase 4 ports 7 ABI entry points (`fenetremempx_`, `mempxfenetre_`, `sauvefenetre_`, `restaurefenetre_`, `sauvemempx_`, `restauremempx_`, `effacemempx_`) from warn-once stubs to real Qt 6 bodies, backed by a single lazily-allocated `QPixmap* saved_canvas_` field on `XvueState`. CONTEXT.md locked every design decision; this research answers five specific implementation questions the user asked and produces a Validation Architecture section for the planner.

**Key verified facts:**
1. Qt 6 explicitly permits **two concurrent `QPainter`s on two different `QPaintDevice`s** — Phase 4 D-07's "temporary scoped painter on `saved_canvas_` while the long-lived `painter_` stays bound to `backing_`" is correct and has no documented gotcha. [CITED: doc.qt.io/qt-6/qpainter.html]
2. `QPainter::drawPixmap(x, y, src)` coordinates are **logical**, not device pixels. When source and destination DPRs match, "it is drawn directly onto the device with no additional transformation applied." [CITED: doc.qt.io/qt-6/qpainter.html#drawPixmap]
3. `QPixmap` is **implicitly shared / copy-on-write**; `QPixmap::copy()` is an explicit deep copy; `swap()` is O(1). [CITED: doc.qt.io/qt-6/qpixmap.html]
4. `QPixmap::save(path, "PNG")` — used by the existing `xvfermer_` capture hook — is **not documented as deterministic**; Qt's docs are silent on whether it injects a `tIME` chunk. Naive `cmp` on two PNG outputs is therefore fragile. Use `compare -metric AE` from ImageMagick for pixel-level equality. [CITED: doc.qt.io/qt-6/qpixmap.html#save; imagemagick.org/script/compare.php]
5. `QPixmap::fill()` has an **explicit warning**: "The effect of this function is undefined when the pixmap is being painted on." This directly applies to `effacemempx_` (D-10) — it MUST go through `painter_->fillRect(...)`, not `backing_->fill(...)`. [CITED: doc.qt.io/qt-6/qpixmap.html#fill]

**Primary recommendation:** Keep D-07 as written (temporary scoped painter on `saved_canvas_` for save, active `painter_->drawPixmap` for restore). The three Qt invariants above all align; no amendment needed. Use `compare -metric AE` (ImageMagick) as the round-trip validation primitive — `cmp` is a failure mode waiting to happen.

## User Constraints (from CONTEXT.md)

### Locked Decisions

- **D-01:** Single `QPixmap* saved_canvas_ = nullptr` field on `XvueState`; raw pointer, manual lifecycle. One slot is sufficient — legacy `xvuelc.c` only allocates one `mempxsauvfen` (xvuelc.c:138 + :1036).
- **D-02:** Lazy allocation on first `sauvefenetre_` / `sauvemempx_`. `new QPixmap(backing_->size())` with `setDevicePixelRatio(backing_->devicePixelRatio())`.
- **D-03:** Destroyed in `~XvueState()` before `painter_` and `backing_`.
- **D-04:** `fenetremempx_` and `mempxfenetre_` are **true no-ops** (empty bodies, no `update()`, no `processEvents()`, no warn_once).
- **D-05:** Each empty body carries a one-line `// D-04 (Phase 4): no-op — backing_ is the visible surface (Phase 2 D-04)` comment.
- **D-06:** Phase 03-04 close empirically validated MEMPXFENETRE-from-Fortran-tracer on Qt even when stubbed — behavior is unchanged, only diagnostic noise removed.
- **D-07:** Two file-local `static inline` helpers in `xvue_qt_api.cpp`:
  - `xvue_qt_save_to_slot()` — ensures slot exists and size-matches; `QPainter::drawPixmap(0, 0, *backing_)` via a **temporary scoped painter on `saved_canvas_`**; active long-lived `painter_` is NOT touched.
  - `xvue_qt_restore_from_slot()` — null/size-mismatch → stderr warn + return; else `painter_->drawPixmap(0, 0, *saved_canvas_)` + `canvas_->update()` (no `processEvents`).
- **D-08:** `sauvefenetre_`/`sauvemempx_` share body (save helper); `restaurefenetre_`/`restauremempx_` share body (restore helper). All 4 ABI symbols stay distinct — no consolidation (preserves `verify_abi = 57`).
- **D-09:** Every body still starts with `XVUE_QT_ASSERT_MAIN_THREAD()` (Phase 2 D-32 invariant).
- **D-10:** `effacemempx_` body = `painter_->fillRect(backing_->rect(), background_)` + `canvas_->update()`. Functionally identical to `effacer_` (Phase 2 D-15) since Qt has no separate mempx surface.
- **D-11:** `// D-10 (Phase 4): same body as effacer_ — Qt has no separate mempx surface (Phase 2 D-04, D-15)`.
- **D-12:** Size mismatch at restore → **no-op with stderr warning** (mirrors X11 `XCopyArea` undefined behavior but auditable).
- **D-13:** No proactive invalidation on `resizeEvent`; lazy size check at restore time only.
- **D-14:** Save that needs a resize `delete`s old `saved_canvas_` and allocates fresh. No content preservation across size change.
- **D-15..D-17:** Validation = synthetic save/draw/restore round-trip via `prpr/xvtest0.f` extension + existing `bin/qt-capture.sh` + `bin/xvtest-capture.sh`. Three sub-tests (sauvefenetre/restaurefenetre, sauvemempx/restauremempx, effacemempx). A/B comparison against legacy X11 backend.
- **D-18:** ROADMAP SC-4 (interactive cavity2d rubber-band) is **explicitly deferred to Phase 5** — Phase 4 validation map records SC-4 as `deferred-to-phase-5`.

### Claude's Discretion

- **Sanity-driver vehicle:** extend `prpr/xvtest0.f` with PHASE 4 COVERAGE block (matches Phase 2 D-36 and Phase 3 — consistency wins).
- **Pixel-diff helper:** `cmp` (byte compare) or `convert -metric AE` (pixel compare). Research recommends **`magick compare -metric AE`** (see Pitfall 3).
- **Save-side helper layout:** file-local `static inline` in `xvue_qt_api.cpp`; split to `xvue_qt_state.cpp` only if helpers exceed ~15 lines.
- **Lazy vs eager allocation of `saved_canvas_`:** lazy preferred; eager acceptable if it simplifies error handling.
- **Stderr message wording:** planner's call to match existing `xvue_qt_api.cpp` diagnostic style.

### Deferred Ideas (OUT OF SCOPE)

- `mempxaccro` (snap-handle pixmap) → Phase 5 (owns `xvsouris_`/`xvsouris2_`).
- Interactive cavity2d rubber-band-drag HUMAN-UAT → Phase 5 event bridge.
- N-slot named pixmap stack → rejected at D-01; reversible if a future phase needs it.
- Real second draw target `mempx_` → rejected at D-04 (would revert Phase 2 D-04 single-backing collapse).
- Stale Fortran latent UB audit → tracked separately under libgfortran5 hold.

## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| PIXMAP-01 | `fenetremempx_` / `mempxfenetre_` copy between backing and off-screen slot | D-04 says both are no-ops on Qt (single-backing collapse). Research Concern 1 confirms Qt 6 permits the temporary scoped painter used when a real copy is needed (save/restore path). |
| PIXMAP-02 | `sauvefenetre_` / `restaurefenetre_` save/restore full canvas to/from named slot | D-07 save/restore helpers. Research Concern 1 + 4 confirm temporary scoped `QPainter` on `saved_canvas_` is the right approach; implicit sharing (Concern 3) doesn't hurt this path. |
| PIXMAP-03 | `sauvemempx_` / `restauremempx_` / `effacemempx_` secondary slots | D-08 same-body consolidation. D-10 `effacemempx_` = `painter_->fillRect(..., background_)`. Research Concern 4b confirms `painter_->fillRect` is required (QPixmap::fill is undefined while painted on). |
| PIXMAP-04 | Rubber-band-drag flicker-free (cavity2d interactive) | **Deferred to Phase 5** per D-18 (requires event bridge). Phase 4 still ships the API + headless round-trip validation (SC-1, SC-2, SC-3 green). |

## Project Constraints (from CLAUDE.md)

- **Compilation must never break:** after every change `bin/cbl_tout_qt && bin/cbl_tout` must stay green.
- **Fortran 77 fixed-form:** `prpr/xvtest0.f` extension follows columns-7+ fixed-form rules (verify trailing blanks, `C` in column 1 for comments, `CALL` starting at column 7).
- **ABI stability:** `verify_abi` target counts 57 symbols; Phase 4 changes 7 stub bodies in place, adds zero symbols.
- **libgfortran5 pinned** to 15.2.0-9 (`apt-mark hold`); builds must use the `/tmp/gfortran-14-shim` PATH to resolve `gcc`/`cc`/`gfortran` to `-14` variants.
- **Never bypass hooks, never force-push, never try to install Ubuntu packages.**
- **Asking before acting:** if ImageMagick or `cmp` are not present, ask — do not work around.
- **MEFISTO normes** (`doc/normes.ps`): naming, comment style, fixed-form column layout.
- **Interactive exit:** never Ctrl-C; return to main menu and type `99;`.
- **`XVUE_QT_ASSERT_MAIN_THREAD()` first line of every ABI entry point body** (Phase 2 D-32 — also restated as Phase 4 D-09).

## Standard Stack

No new libraries. Phase 4 runs on the exact stack already used by Phases 1–3:

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Qt 6 Widgets/Gui/Core | system (currently 6.8.x+ on Debian sid) | `QPixmap`, `QPainter`, `QWidget` | Already in use since Phase 1; chosen at init over Qt 5 |
| gfortran | 14 via `/tmp/gfortran-14-shim` | Fortran 77 compilation of `xvtest0.f` and solver sources | Project-wide pin (see STATE.md libgfortran5 deferred item) |
| ImageMagick | system (`/usr/bin/convert`, `/usr/bin/magick`) | Pixel-diff helper for round-trip validation | Already used by `bin/xvtest-capture.sh`; presence verified at research time |
| Xvfb | system | X11 headless for legacy A/B capture | Already used by `bin/xvtest-capture.sh` |

**Installation:** none — everything required is already installed or pinned. No `apt install`, no `cmake`/pkg-config additions.

## Runtime State Inventory

Phase 4 is an in-place replacement of 7 stub function bodies plus one field on `XvueState`. It is NOT a rename or migration; still, for auditability:

| Category | Items | Action |
|----------|-------|--------|
| Stored data | None — no datastores, no persisted project files touched | None |
| Live service config | None — no services | None |
| OS-registered state | None | None |
| Secrets/env vars | Phase 4 adds zero env vars. Reuses `MEFISTO_QT_CAPTURE_PATH`, `MEFISTO_XVSOURIS_AUTOEXIT`, `MEFISTO_XVFERMER_READY_FILE`, `MEFISTO_XVFERMER_HOLD_MS` from Phase 03-04 unchanged | None |
| Build artifacts | `pp/ppxvtest0_qt` rebuilt via existing `bin/cbxvtest0_qt`; `xvue/qt/build/libxvueqt.*` rebuilt via `bin/cbl_tout_qt`. No new artifact names, no stale egg-info. | Standard `bin/cbl_tout_qt` rebuild |

## Architecture Patterns

### Pattern 1: File-local `static inline` helper in `xvue_qt_api.cpp`
**What:** Two helpers (`xvue_qt_save_to_slot`, `xvue_qt_restore_from_slot`) share bodies across the four save/restore ABI entry points.
**When:** Matches the existing `xvue_qt_draw_rect_common` pattern (xvue_qt_api.cpp:58–75).
**Example (reference, NOT the Phase 4 code):**
```cpp
inline void xvue_qt_draw_rect_common(int x, int y, int w, int h, RectMode mode) {
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;
    // ... body ...
    if (win->canvas()) win->canvas()->update();
}
```

### Pattern 2: Temporary scoped `QPainter` on a second `QPixmap`
**What:** Paint to a pixmap that is NOT the active long-lived painter's target by constructing a short-lived `QPainter` in a block scope.
**When:** Save direction (`saved_canvas_ = backing_`). Qt 6 allows this because the restriction is **one painter per device**, not one painter per process. Two devices can each have their own active painter.
**Reference (from `xvue_qt_canvas.cpp:69–75` — resize path, proven pattern):**
```cpp
{
    QPainter tmp(new_backing);
    tmp.fillRect(new_backing->rect(), state_->background_);
    if (old_backing) {
        tmp.drawPixmap(0, 0, *old_backing);
    }
}  // ~QPainter tmp — end() called automatically
```
Phase 4's `xvue_qt_save_to_slot()` follows this exact shape, swapping the role of source and destination.

### Pattern 3: Restore via active long-lived painter
**What:** For the restore direction, `painter_` is already active on `backing_`; just blit.
**When:** `restaurefenetre_` / `restauremempx_` — the destination is `backing_`, which already has `painter_` bound to it.
**Shape:**
```cpp
// Preconditions asserted: st, st->painter_ active, st->backing_, st->saved_canvas_ present and size-matched
st->painter_->drawPixmap(0, 0, *st->saved_canvas_);
win->canvas()->update();
// No processEvents — light epilogue, matches Phase 2 D-01 lighter variant
```

### Pattern 4: Extend `prpr/xvtest0.f` with a new coverage block
**What:** Add a PHASE 4 COVERAGE section after the TEXT coverage block, before the second `XVFERMER`.
**When:** Matches Phase 2 D-36 draw-coverage section and Phase 3 D-24 text-coverage section. Consistency-driven; planner should not deviate without strong reason.

### Anti-Patterns to Avoid

- **Calling `backing_->fill(background_)` in `effacemempx_`:** Qt docs explicitly say `QPixmap::fill()` is undefined while the pixmap is being painted on, and `painter_` is permanently bound to `backing_`. Must use `painter_->fillRect(backing_->rect(), background_)` instead.
- **`QPixmap::swap()` between `backing_` and `saved_canvas_`:** would break the invariant that `painter_` lives on `backing_`. Phase 2 D-05 says `painter_` stays bound to `backing_` for the whole session — swapping the underlying pixel data silently works (implicit sharing) but swapping the objects themselves breaks the painter's device pointer.
- **Two `QPainter`s on the same `QPixmap`:** Qt 6 explicitly forbids this and `painter->begin(sameDevice)` fails. Not applicable here (Phase 4 only paints on two *different* devices concurrently) but worth naming.
- **`cmp` byte-compare on two `QPixmap::save(*.png)` outputs:** Qt's docs are silent on whether PNG metadata (tIME chunk, software tag) is deterministic. Use `magick compare -metric AE` instead.
- **Hooking `resizeEvent` to invalidate `saved_canvas_`:** D-13 explicitly forbids this. The size check at restore time is sufficient and keeps Phase 4 decoupled from Phase 2's resize logic.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Pixel-diff between two PNGs | Byte-level `cmp` + pray PNG metadata is stable | `magick compare -metric AE a.png b.png null:` (exit 0 if identical pixels, AE count on stderr) | Qt's `QPixmap::save("PNG")` determinism is not documented; `compare -metric AE` ignores PNG chunks and works on decoded pixel data |
| Deep copy of a QPixmap | `for (y) for (x) setPixel()` | `QPainter::drawPixmap(0,0,src)` OR `QPixmap::copy()` | Both are O(pixels) but hardware-accelerated; the first is in-place into a preallocated slot (Phase 4 D-07); the second allocates |
| Headless screenshot of Qt canvas | External screenshot tool | `MEFISTO_QT_CAPTURE_PATH` + existing `xvfermer_` in-process `backing_->save()` hook | Already exists (xvue_qt_api.cpp:570–606), already works under `QT_QPA_PLATFORM=offscreen`, no Xvfb needed |
| Xvfb orchestration for X11 A/B | Custom Xvfb bring-up | `bin/xvtest-capture.sh` + `MEFISTO_XVFERMER_READY_FILE` sentinel | Battle-tested in Phase 03-04; reuse unchanged |

**Key insight:** Phase 4's tooling surface is **zero new scripts**. Everything the validation harness needs already exists.

## Research Concern Deep-Dives

### Concern 1: Two `QPainter`s on two different `QPixmap`s (D-07 correctness)

**Question:** Can two QPainters coexist if they target different QPixmaps?

**Answer:** YES — Qt 6 documentation explicitly permits this. [CITED: doc.qt.io/qt-6/qpainter.html]

The restriction is phrased "A paint device can only be painted by one painter at a time" — the constraint is per-device, not per-process. Qt's documented example of failure is literally `painter->begin(myWidget); painter2->begin(myWidget)` — the same device twice.

For Phase 4, the active long-lived `painter_` is bound to `backing_`. Inside `xvue_qt_save_to_slot()`, a temporary `QPainter tmp(saved_canvas_)` begins on a different device. This is legal. The `resizeEvent` code at `xvue_qt_canvas.cpp:69–75` already does exactly this (temporary painter on `new_backing` while `painter_` is still bound to `old_backing`); that path has been green since Phase 2 landed. D-07 is correct as written.

**Subtlety:** `drawPixmap(x, y, src)` reads pixels from `src`; it does not need a painter on `src`. Qt docs: "A source QPixmap can be safely drawn by a QPainter even if no other painter is actively modifying it" — and by symmetry, even if no painter is active on it at all. The save direction (temporary on slot, blit from backing_) is safe because nothing tries to paint on `backing_` during the temporary's lifetime; the restore direction (active `painter_` on `backing_` blits from `saved_canvas_`) is safe because `saved_canvas_` has no painter.

### Concern 2: devicePixelRatio round-trip integrity

**Question:** Does a save/restore of a QPixmap with DPR=2 preserve all pixels with no DPI drift? Is the `(0, 0)` in `drawPixmap` logical or device?

**Answer:**
- `QPainter::drawPixmap` coordinates are **logical**, not device pixels. Qt: "we specify points using logical coordinates which then are converted into the physical coordinates of the paint device." [CITED: doc.qt.io/qt-6/qpainter.html]
- When source and destination DPRs match: "Should [the source DPR] match the value of the underlying QPaintDevice, it is drawn directly onto the device with no additional transformation applied." [CITED: doc.qt.io/qt-6/qpainter.html#drawPixmap]
- D-02 allocates the slot with `setDevicePixelRatio(backing_->devicePixelRatio())` — identical DPR. Round-trip is lossless.

**However**, Qt's official docs are **silent** on whether `QPixmap::copy()` preserves `devicePixelRatio`. Phase 4 D-07 does NOT use `copy()` — it uses `drawPixmap` on a preallocated slot whose DPR we set explicitly in D-02. So this silence does not affect us. Worth noting as a general pitfall: anywhere else in Qt code that uses `copy()` should explicitly re-apply `setDevicePixelRatio` to be safe.

**Action for D-02:** keep the explicit `setDevicePixelRatio` call. It's correct AND insulates us from the (unspecified) `copy()` behavior if a future reader ever reaches for `copy()` as an alternative.

### Concern 3: `QPixmap::copy()` vs `QPainter::drawPixmap()` vs `QPixmap::swap()`

**Verified facts:**
- `QPixmap` is implicitly shared / copy-on-write. "QPixmap objects can be passed around by value." [CITED: doc.qt.io/qt-6/qpixmap.html]
- `copy()` returns a **deep copy** ("Returns a deep copy of the subset of the pixmap").
- `swap()` "is very fast and never fails" — O(1) pointer swap.
- `operator=` has both copy and move variants; both respect implicit sharing.

**Alternatives evaluated for the save direction:**

| Alt | Mechanism | Allocation | Verdict |
|-----|-----------|-----------|---------|
| A | `saved_canvas_ = new QPixmap(backing_->copy())` each save | Yes, deep copy | Simpler (no temp painter), but allocates per save call. Per-call overhead small but adds a `new`/`delete` churn on every rubber-band-motion frame in Phase 5. |
| B (D-07) | Preallocated slot + temp `QPainter` + `drawPixmap(0,0,*backing_)` | No (reuses slot; deep-copies pixel data in place) | Chosen. Zero churn. Proven shape (resizeEvent uses same pattern). |
| C | `*saved_canvas_ = *backing_` (copy assignment) | No (COW — first write to either triggers detach) | **Dangerous**: because the next draw to `backing_` through `painter_` would trigger QPixmap COW detachment on `backing_`, but `painter_` is already bound to the pre-detach device — undefined. Rejected. |
| D | `saved_canvas_->swap(*backing_)` | No, O(1) | **Breaks D-05 invariant**: `painter_` would end up bound to what is now the saved slot. Rejected. |

**Winning approach: B (= D-07 current).** No amendment to D-07 needed. Rationale cites the Qt invariants, not just intuition.

### Concern 4a: `QPixmap::save` PNG determinism

**Finding:** Qt's docs **do not promise** that `QPixmap::save(path, "PNG")` produces byte-identical output for byte-identical input. [CITED: doc.qt.io/qt-6/qpixmap.html#save — "silent on determinism"]

**Impact:** The existing `xvfermer_` capture hook (xvue_qt_api.cpp:588–594) uses exactly this call. If Qt's PNG writer injects a `tIME` chunk or varies the zlib compression level across runs, two otherwise-identical runs produce PNGs that differ by a few bytes even though the pixels are identical.

**Validated mitigation:** ImageMagick's `compare -metric AE` operates on decoded pixel data. From imagemagick.org/script/compare.php: exit "0 if the images are similar" — AE counts absolute-error pixels and is zero for pixel-identical inputs regardless of PNG chunk differences. [CITED: imagemagick.org/script/compare.php]

**Recommended invocation (for the planner's harness):**
```bash
# Returns 0 on pixel-identical; prints AE count on stderr
magick compare -metric AE "$A" "$B" null: 2> /tmp/ae_count
# $A and $B identical-pixel -> exit 0, /tmp/ae_count == "0"
```

**Fallback:** `cmp -s "$A" "$B"` as an opportunistic fast path with AE as the authoritative tiebreaker. But recommendation is: use AE as primary; `cmp` is a false-fail waiting to happen.

### Concern 4b: `QPixmap::fill()` undefined-while-painted-on

**Finding:** Qt explicitly documents: "The effect of this function is undefined when the pixmap is being painted on." [CITED: doc.qt.io/qt-6/qpixmap.html#fill]

**Impact on D-10:** `effacemempx_` clears the canvas to background. Since `painter_` is permanently bound to `backing_` (Phase 2 D-05), calling `backing_->fill(background_)` is in the explicitly-undefined zone. D-10 already says "fill via the active `painter_->fillRect(backing_->rect(), background_)`" — this is correct. Phase 2's `effacer_` body does the same thing (confirming the shared-body claim in D-11).

**Cross-check against Phase 2 `effacer_`:** look at the existing `effacer_` body in `xvue_qt_api.cpp` ~ line 451 (the stub just above says `// ---- 22. effacer_ (D-15) ----`) — D-10's body MUST be bit-identical except for the D-11 comment. The planner should verify this with a `diff` between the proposed `effacemempx_` body and the Phase 2 `effacer_` body.

### Concern 5: Validation harness mechanics

**Q1: Is extending `prpr/xvtest0.f` the lowest-friction path?**

YES. Verified by direct read of `prpr/xvtest0.f` (168 lines, already has SHELL, DRAW, and TEXT coverage sections at lines 46–100, 102–142). A PHASE 4 COVERAGE block inserts naturally after the TEXT block (before line 144 `CALL MEMPXFENETRE`). No new driver binary → no new `bin/cb*` script → no `cbl_tout_qt` wiring changes → blast radius stays tiny, matching the Phase 03-04 convention.

**Q2: Is round-trip pixel equality the right validation primitive?**

YES with a caveat: **pixel equality, not byte equality**. See Concern 4a above.

**Q3: Do the existing harnesses need modification?**

NO. Verified by read:
- `bin/qt-capture.sh` (lines 1–80) already handles the Qt offscreen case; env contract is `MEFISTO_QT_CAPTURE_PATH` → backing pixmap saved via `xvfermer_` hook. Phase 4's round-trip test can take **multiple** captures from a single driver run by re-invoking `xvfermer_` — but that's heavy. Simpler: run the driver 2× with different env-var branches inside the Fortran source selecting which scene to draw. Or: run the driver once and have it write TWO captures via an auxiliary sibling hook. **Recommendation for planner:** one tiny wrapper script `bin/xvtest0-pixmap-roundtrip.sh` that invokes `pp/ppxvtest0_qt` multiple times with different scene-selector env vars, captures N PNGs, then runs `magick compare -metric AE` pairwise. No modification of `qt-capture.sh` or `xvtest-capture.sh`.
- `bin/xvtest-capture.sh` (lines 1–140) handles the X11 legacy capture unchanged.
- `bin/cbxvtest0_qt` (lines 1–68) rebuilds the Qt driver; already PATH-aware and reuses `MEFISTO`/shell env. No changes needed.

**Capture reality check:** The Qt capture hook at `xvue_qt_api.cpp:570–606` calls `st->painter_->end()` **before** `backing_->save()`. This means whatever state `backing_` is in when `xvfermer_` is called is what gets captured. For the round-trip test, the sequence `sauvefenetre → draw overlay → restaurefenetre → xvfermer` puts `backing_` back into its pre-overlay state — exactly what we want to compare against a control capture made from `draw-base-scene → xvfermer`.

## Code Examples

### Example 1: Save helper (temporary scoped painter)
```cpp
// File-local, xvue_qt_api.cpp — NOT the Fortran ABI
static inline void xvue_qt_save_to_slot() {
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->backing_) return;

    // Allocate or resize the slot (D-02, D-14)
    if (!st->saved_canvas_ || st->saved_canvas_->size() != st->backing_->size()) {
        delete st->saved_canvas_;
        st->saved_canvas_ = new QPixmap(st->backing_->size());
        st->saved_canvas_->setDevicePixelRatio(st->backing_->devicePixelRatio());
    }

    // D-07: temporary scoped painter on the slot; active painter_ on backing_ is NOT touched.
    // Qt 6 permits two concurrent QPainters on two different QPaintDevices
    // (doc.qt.io/qt-6/qpainter.html — "one painter per device").
    {
        QPainter tmp(st->saved_canvas_);
        tmp.drawPixmap(0, 0, *st->backing_);
    }  // ~QPainter tmp — end() called automatically
}
```

### Example 2: Restore helper (reuse active painter)
```cpp
static inline void xvue_qt_restore_from_slot() {
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (!st || !st->painter_ || !st->painter_->isActive()) return;

    // D-12: size mismatch or missing slot → stderr warn + no-op
    if (!st->saved_canvas_ || st->saved_canvas_->size() != st->backing_->size()) {
        std::fprintf(stderr,
            "xvue-qt: restore_from_slot: no slot or size mismatch\n");
        return;
    }

    // D-07: blit slot → backing_ via the active long-lived painter_.
    st->painter_->drawPixmap(0, 0, *st->saved_canvas_);
    if (win->canvas()) win->canvas()->update();
    // No processEvents — Phase 2 D-01 lighter epilogue.
}
```

### Example 3: Fortran 77 PHASE 4 COVERAGE block (sketch)
```fortran
C     -- begin PIXMAP coverage section (Phase 4, D-15) --------------
C     Round-trip: save full scene, draw overlay, restore, capture
C     should match a control capture made with no overlay.
      CALL XVCOULEUR(7)
      CALL XVTRAIT(50, 500, 750, 500)
      CALL SAUVEFENETRE
C     magenta rubber-band overlay
      CALL XVCOULEUR(6)
      CALL XVBORDRECTANGLE(200, 520, 400, 40)
      CALL RESTAUREFENETRE
      CALL MEMPXFENETRE
      CALL XVVOIR
C     -- end PIXMAP coverage section ------------------------------
```
(The planner will split this into the three sub-tests of D-16.)

## Common Pitfalls

### Pitfall 1: `QPixmap::fill()` while `painter_` is active
**What goes wrong:** `backing_->fill(bg)` in `effacemempx_` is explicitly undefined in Qt 6 docs when a painter is bound to the pixmap.
**Why it happens:** Natural instinct reaches for `fill()` first; the doc warning is in a different page section than the `QPixmap::fill` overview.
**How to avoid:** Always clear via `painter_->fillRect(rect, color)` when `painter_` is the long-lived one. D-10 already says so; this pitfall is for defense in depth during code review.
**Warning sign:** Any `state_->backing_->fill(` in `xvue_qt_api.cpp` or `xvue_qt_state.cpp` outside of an allocation path.

### Pitfall 2: Calling `painter_->begin()` on a device with another active painter
**What goes wrong:** Qt's `QPainter::begin` fails with the printed warning "QPainter::begin: A paint device can only be painted by one painter at a time." The target draw is silently dropped.
**Why it happens:** If a future maintainer tries to "reuse `painter_` on the slot" for the save direction instead of using a temporary, they'd either have to `end()` and re-`begin()` (breaking the Phase 2 D-05 invariant) or hit this error.
**How to avoid:** Always use a **temporary scoped** `QPainter tmp(slot)` for the save direction, exactly as `resizeEvent` at `xvue_qt_canvas.cpp:69–75` does.
**Warning sign:** Any call to `state_->painter_->begin(...)` added during Phase 4.

### Pitfall 3: PNG byte-compare as validation primitive
**What goes wrong:** `cmp -s qt.png qt_roundtrip.png` passes on one build and fails on another because Qt's `QImageWriter` backend varies the PNG `tIME` chunk or zlib compression level between runs.
**Why it happens:** Qt docs don't promise `QPixmap::save("PNG")` is deterministic. We verified the docs are silent.
**How to avoid:** Use `magick compare -metric AE A B null:` — exit 0 iff pixel-identical, independent of PNG metadata.
**Warning sign:** Any use of `cmp` or `diff` on PNG files in the Phase 4 harness without a pixel-level fallback.

### Pitfall 4: Size mismatch at restore after a window resize
**What goes wrong:** User resizes the window between `sauvefenetre_` and `restaurefenetre_`. `backing_` is the new size; `saved_canvas_` is the old size. Legacy X11 `XCopyArea` silently clips; the Qt version should stderr-warn and bail.
**Why it happens:** D-13 deliberately does not invalidate the slot on resize.
**How to avoid:** The D-12 size check at the top of `xvue_qt_restore_from_slot()` is the only line of defense. It must be present and it must return early without calling `drawPixmap`.
**Warning sign:** Any Phase 4 code path that calls `painter_->drawPixmap(0, 0, *saved_canvas_)` without first verifying `saved_canvas_->size() == backing_->size()`.

### Pitfall 5: Accidentally invoking implicit `QPixmap` sharing + COW detach mid-paint
**What goes wrong:** If someone replaces D-07's `drawPixmap` with `*saved_canvas_ = *backing_` (copy-assign), the next draw through `painter_` to `backing_` triggers QPixmap COW detachment of `backing_`, but `painter_->device()` still points to the pre-detach pointer. Undefined behavior.
**Why it happens:** Implicit sharing looks like "free deep copy" to someone who hasn't read Phase 2 D-05.
**How to avoid:** Strict `drawPixmap` only. Never use `operator=` or `swap` between `saved_canvas_` and `backing_`.
**Warning sign:** Any `*state_->saved_canvas_ = ...` or `saved_canvas_->swap(...)` in Phase 4 code.

### Pitfall 6: Forgetting `XVUE_QT_ASSERT_MAIN_THREAD()` on the two no-op bodies
**What goes wrong:** D-04 (`fenetremempx_`, `mempxfenetre_`) says the bodies are empty. "Empty" here means "no business logic" — but Phase 2 D-32 says every ABI entry point must still assert main-thread first. D-09 restates this explicitly.
**Why it happens:** Natural reading of "true no-op" → literally empty function body.
**How to avoid:** Read D-09 together with D-04; the body is `XvueApp::ensure(); XVUE_QT_ASSERT_MAIN_THREAD(); /* D-04 comment */` — two lines of scaffolding, zero lines of logic.
**Warning sign:** A `fenetremempx_` body with fewer than 3 lines (assert + comment + braces).

### Pitfall 7: `verify_abi` count drift
**What goes wrong:** Changing 7 stub bodies to real bodies should be symbol-count-neutral, but if a helper (e.g., `xvue_qt_save_to_slot`) is accidentally declared `extern "C"` or missing `static inline`, it shows up in the symbol list and `verify_abi` fails with count 58.
**Why it happens:** Copy-paste from an existing ABI entry point template.
**How to avoid:** Both helpers MUST be `static inline` (file-local, no external linkage). D-07 already specifies this; the planner should include an explicit verify step in the task DoD.
**Warning sign:** CMake `verify_abi` target reports `expected 57, got 58`.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| Qt 6 Widgets/Gui/Core | `xvue/qt/build/libxvueqt` | Assumed (Phase 1–3 green) | system | — |
| gfortran 14 | `xvtest0.f` compile | Via `/tmp/gfortran-14-shim` | 14 | — |
| libgfortran5 (pinned 15.2.0-9) | Runtime | Verified at research time (`dpkg -l`) | `15.2.0-9` (apt-mark hold) | — |
| ImageMagick `convert` + `magick` | Validation harness pixel diff | Verified: `/usr/bin/convert`, `/usr/bin/magick` | system | `cmp` opportunistic (weaker) |
| `cmp` | Validation harness (opportunistic) | Verified: `/usr/bin/cmp` | system | — |
| Xvfb | Legacy X11 A/B capture via `bin/xvtest-capture.sh` | Assumed (Phase 03-04 uses it) | system | — |
| xcb-cursor0 | NOT required — Qt offscreen path bypasses it | — | — | — |

**Missing dependencies with no fallback:** None.
**Missing dependencies with fallback:** None — everything Phase 4 needs is installed.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Shell drivers + visual A/B captures + pixel-diff (matches Phase 03 convention) |
| Config file | None (shell-based; extends `prpr/xvtest0.f`) |
| Quick run command | `bin/cbl_tout_qt && pp/ppxvtest0_qt` (headless smoke via `QT_QPA_PLATFORM=offscreen`) |
| Full suite command | `bin/xvtest0-pixmap-roundtrip.sh` (Phase 4 creates this tiny wrapper) which invokes `bin/qt-capture.sh` N times + `magick compare -metric AE` pairwise + `bin/xvtest-capture.sh` for legacy X11 A/B |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| PIXMAP-01 | `fenetremempx_` / `mempxfenetre_` no-ops do not break the scene (single-backing collapse) | smoke | `bin/qt-capture.sh pp/ppxvtest0_qt /tmp/p4_fenetremempx.png` — assert exit 0 + PNG written | ❌ Wave 0 extends `prpr/xvtest0.f` with `CALL FENETREMEMPX` + `CALL MEMPXFENETRE` |
| PIXMAP-02a | `sauvefenetre_` + `restaurefenetre_` round-trip pixel-identical to control | integration | `magick compare -metric AE /tmp/p4_ctrl.png /tmp/p4_saverestore.png null:` — exit 0 | ❌ Wave 0 extends `xvtest0.f` + adds `bin/xvtest0-pixmap-roundtrip.sh` |
| PIXMAP-02b | Same test against legacy X11 backend | integration | `bin/xvtest-capture.sh pp/ppxvtest0 /tmp/p4_x11.png` + visual A/B vs Qt capture | ❌ Wave 0 — same harness |
| PIXMAP-03a | `sauvemempx_` / `restauremempx_` round-trip pixel-identical to control | integration | `magick compare -metric AE /tmp/p4_ctrl.png /tmp/p4_mempx_saverestore.png null:` — exit 0 | ❌ Wave 0 |
| PIXMAP-03b | `effacemempx_` clears to background | integration | `magick compare -metric AE /tmp/p4_bg.png /tmp/p4_effacemempx.png null:` — exit 0; `/tmp/p4_bg.png` = capture right after `xvinitgraphique` + `effacer` | ❌ Wave 0 |
| PIXMAP-03c | `effacemempx_` body byte-identical (or functionally identical) to `effacer_` body | unit (code review) | Manual diff of the two function bodies in `xvue_qt_api.cpp` + comment cross-reference | Built into Task DoD |
| PIXMAP-04 | Interactive cavity2d rubber-band-drag no-flicker | **DEFERRED to Phase 5** | — | Validation map records `deferred-to-phase-5` per D-18 |

### Sampling Rate
- **Per task commit:** `bin/cbl_tout_qt && pp/ppxvtest0_qt` (rebuild Qt lib + xvtest0_qt, run once headless under `QT_QPA_PLATFORM=offscreen`) — ~15–25 s
- **Per wave merge:** `bin/xvtest0-pixmap-roundtrip.sh` (full round-trip with pixel-diff) + `bin/cbl_tout` (full legacy X11 build must stay green per CLAUDE.md rule "Compilation must never break") — ~45–60 s
- **Phase gate before `/gsd-verify-work`:** Full suite green + `verify_abi` target reports 57 + `verify_no_exec` target green + A/B capture pair visually reviewed by orchestrator using the `Read` tool on both PNGs (D-17)

### Wave 0 Gaps
- [ ] `prpr/xvtest0.f` — extend with PHASE 4 COVERAGE block (D-15 scene + 3 sub-tests of D-16, bracketed by `XVBEGIN`/`XVEND`-style comment banners matching existing DRAW and TEXT coverage sections). Planner owns the exact sequencing of save/draw-overlay/restore/capture between the sub-tests.
- [ ] `bin/xvtest0-pixmap-roundtrip.sh` — new ~50-line wrapper: invokes `pp/ppxvtest0_qt` under `bin/qt-capture.sh` for each sub-test scene (via a `MEFISTO_PHASE4_SCENE={ctrl,saverestore,mempx_roundtrip,effacemempx}` env-var selector the Fortran driver reads via `GETENV`), then pairwise `magick compare -metric AE`. Exit code 0 iff all AE counts are 0.
- [ ] ImageMagick probe at script start: `command -v magick >/dev/null || { echo "magick not found"; exit 2; }` — already verified installed but keep the guard per CLAUDE.md "ask before acting."
- [ ] `bin/cbxvtest0_qt` — **no changes needed** unless the planner picks a sibling driver vehicle over the `xvtest0.f` extension (planner's call; research recommends extension for consistency with Phases 2, 3).

*(If no gaps surface during planning: "None — existing test infrastructure covers all phase requirements" but two items above are firm from the research.)*

## Plan Implications

Each research finding maps back to a CONTEXT.md decision it **confirms** (keeps as-is) or **amends** (surfaces a change the planner should apply):

| Finding | CONTEXT.md decision | Verdict |
|---------|---------------------|---------|
| Qt 6 permits two QPainters on two different devices | D-07 (temp scoped painter on slot) | **Confirms** — no change. Add a citation to `doc.qt.io/qt-6/qpainter.html` in the D-07 comment when the planner writes the helper. |
| `drawPixmap` coordinates are logical; matching DPR draws directly with no transformation | D-02 (`setDevicePixelRatio` on slot) | **Confirms** — no change. |
| `QPixmap::copy()` preservation of DPR is undocumented | — (D-07 uses drawPixmap, not copy) | **Confirms** D-07 is the right choice for a different reason than originally stated. |
| `QPixmap` is implicitly shared; `operator=`/`swap` are not safe substitutes for D-07 | D-07 (drawPixmap via temp painter) | **Confirms** — no change. Pitfall 5 documents why alternatives are rejected. |
| `QPixmap::fill()` is undefined while painted on | D-10 (`effacemempx_` = `painter_->fillRect`) | **Confirms** — no change. Pitfall 1 explains why. |
| `QPixmap::save` PNG output is not documented as deterministic | Claude's Discretion ("cmp or a small pixel-diff helper") | **Amends** — research recommends `magick compare -metric AE` as primary, not `cmp`. CONTEXT.md's wording allowed either; research picks AE. Planner should codify. |
| `bin/qt-capture.sh` + `bin/xvtest-capture.sh` + `xvfermer_` hook are unmodified-ready | D-15, D-16, D-17 | **Confirms** — no harness changes. |
| `prpr/xvtest0.f` extension is the lowest-friction driver vehicle | Claude's Discretion ("xvtest0.f extension, sibling driver, or tests/ C++ harness") | **Confirms** — research picks xvtest0.f extension to match Phases 2, 3. |
| `XVUE_QT_ASSERT_MAIN_THREAD` is required on the no-op bodies too | D-09 | **Confirms** — Pitfall 6 flags this as a common misreading. |
| `verify_abi = 57` stability requires helpers to be `static inline` | D-07 (file-local helpers) | **Confirms** — Pitfall 7 makes this a DoD verification. |
| Capture happens at `xvfermer_` time and ends the painter before saving | D-15/D-16 round-trip design | **Confirms** — restore direction must leave `backing_` in the expected final state before `xvfermer_`; the planner's xvtest0.f sequence must order `RESTAUREFENETRE` before the terminal `XVFERMER`. |

**No decision amendments required.** One Claude's Discretion slot (pixel-diff helper) is now resolved: use `magick compare -metric AE`.

## State of the Art

| Old approach | Current approach | When changed | Impact |
|--------------|------------------|--------------|--------|
| `XCopyArea(src, dst, 0, 0, w, h, 0, 0)` | `QPainter::drawPixmap(0, 0, *src)` via a scoped painter on dst | Phase 2 (single-backing collapse) | Destination is always `backing_`; source is `backing_` (save) or `saved_canvas_` (restore). |
| `XFillRectangle(disp, pixmap, gc, 0, 0, w, h)` | `painter_->fillRect(backing_->rect(), background_)` | Phase 2 D-15 | Background color comes from `XvueState::background_` (applyPen-free path). |
| `XSetFunction(disp, gc, GXcopy)` before copy | nothing | Phase 2 | Qt's `drawPixmap` has no raster-op mode selector; source replaces dest by default (no XOR/AND/OR mode). Phase 4 doesn't need it — the legacy bodies only used `GXcopy` anyway. |
| `Pixmap mempxsauvfen = XCreatePixmap(...)` allocated eagerly in `xvinfo_` | `QPixmap* saved_canvas_ = nullptr` allocated lazily on first save | Phase 4 D-02 | Saves one ~6KB allocation for solver runs that never save. Reversible (D-02 discretion says eager is acceptable). |

**Deprecated/outdated:** none relevant to Phase 4.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | ImageMagick `magick compare -metric AE` ignores PNG metadata (tIME chunk) when comparing pixel data | Concern 4a, Validation Architecture | Low — AE operates on decoded pixel planes by definition, but imagemagick.org/script/compare.php doesn't literally state "metadata is ignored." Empirically will be verified on first Wave 0 run. If wrong, fall back to converting both PNGs to PPM first via `magick A.png A.ppm` before diff. |
| A2 | The existing Qt offscreen capture hook (`MEFISTO_QT_CAPTURE_PATH`) produces a deterministic pixmap for deterministic drawing input — i.e., two runs of the same Fortran driver produce pixel-identical `backing_` state | Validation Architecture, Concern 4a | Medium — the only nondeterminism sources are the PNG encoder (handled by using AE comparison) and font rendering (mitigated by Phase 3's frozen DejaVu Sans Mono choice). If wrong, the round-trip comparison becomes unreliable and we need a control-vs-experiment comparison rather than two separate runs. |
| A3 | A Fortran driver can read an env var via `CALL GETENV('MEFISTO_PHASE4_SCENE', var)` (gfortran intrinsic) to branch scene selection | Validation Architecture "Wave 0 Gaps" bullet 2 | Low — `GETENV` is a well-established gfortran extension; used elsewhere in MEFISTO (search `grep -rn GETENV` confirms). If wrong, alternative is to build N sibling drivers (more cb scripts). |
| A4 | Qt's `QPixmap::save` PNG output may vary across runs | Concern 4a | Low — even if it turns out to be deterministic in practice on this Qt version, using AE comparison is strictly safer than `cmp`. No downside to the conservative choice. |
| A5 | The legacy X11 A/B capture will reach the exact same pixel state as Qt for the round-trip test | Validation Architecture PIXMAP-02b | Medium — legacy and Qt font rendering still differ at subpixel level (Phase 3 close notes). D-17 uses visual A/B (human/orchestrator read of PNGs), not strict pixel equality, precisely to tolerate this. Research confirms D-17's choice. |

## Open Questions

1. **Does the current `xvfermer_` in-process capture hook allow MULTIPLE captures from one driver run, or is it single-shot?**
   - What we know: the hook fires once at `xvfermer_` time; `xvtest0.f` calls `XVFERMER` twice (first + reopen cycle). The hook currently saves to `MEFISTO_QT_CAPTURE_PATH` on EACH call, overwriting. Two captures possible but only the last is kept.
   - What's unclear: whether the Phase 4 harness should (a) run the driver N times with different env selectors, or (b) extend the hook to write sequenced files.
   - Recommendation: (a) — zero changes to `xvfermer_`, one tiny shell wrapper. The planner can pick (b) if extensibility matters, but (a) is the shortest path.

2. **Should the Phase 4 PHASE 4 COVERAGE block in `xvtest0.f` run BEFORE or AFTER the existing DRAW+TEXT coverage blocks?**
   - What we know: drawing state from DRAW+TEXT accumulates in `backing_` before the PHASE 4 block would execute. The round-trip control capture must therefore include those primitives in its expected pixels, which complicates comparison.
   - Recommendation: run PHASE 4 coverage FIRST (right after the first `XVOUVRIR` and before `XVTYPETRAIT(0)` at line 48), so the base scene is a clean black canvas. Then the save/overlay/restore primitives operate on an empty background. Followed by `CALL EFFACER` to reset before the existing DRAW+TEXT blocks run. The planner decides the final ordering.

3. **Is `MEFISTO_PHASE4_SCENE` the right env-var name, or does an existing convention apply?**
   - What we know: existing env vars are `MEFISTO_*` prefixed (MEFISTO_QT_CAPTURE_PATH, MEFISTO_XVSOURIS_AUTOEXIT, MEFISTO_XVFERMER_READY_FILE).
   - Recommendation: stick with `MEFISTO_` prefix. Exact suffix is planner's call — `MEFISTO_XVTEST0_SCENE` is more searchable since it ties the var to the specific driver.

## Pitfalls

*(Consolidated view — all seven pitfalls are documented individually above in "Common Pitfalls". Listed here for scan-ability.)*

1. **`QPixmap::fill()` while `painter_` is active** — use `painter_->fillRect(...)` instead.
2. **`painter_->begin()` on a device with another active painter** — use scoped temporary `QPainter tmp(otherDevice)`.
3. **`cmp` byte-compare on two PNG outputs** — use `magick compare -metric AE A B null:` instead.
4. **Size mismatch at restore after a window resize** — D-12 stderr-warn + early return.
5. **`operator=` or `swap` between `saved_canvas_` and `backing_`** — breaks Phase 2 D-05 invariant; use `drawPixmap` only.
6. **Empty no-op bodies missing `XVUE_QT_ASSERT_MAIN_THREAD`** — D-09 is firm; body has 3 lines (assert + comment + empty braces), not zero.
7. **`verify_abi` count drift** — helpers MUST be `static inline`; any extern linkage breaks the 57-count invariant.

## Sources

### Primary (HIGH confidence — cited, not assumed)
- `https://doc.qt.io/qt-6/qpainter.html` — "A paint device can only be painted by one painter at a time"; per-device constraint, not per-process
- `https://doc.qt.io/qt-6/qpainter.html#drawPixmap` — logical coordinates; matching DPR draws directly without transformation
- `https://doc.qt.io/qt-6/qpixmap.html` — implicitly shared / copy-on-write; `copy()` is deep; `swap()` is fast
- `https://doc.qt.io/qt-6/qpixmap.html#fill` — "The effect of this function is undefined when the pixmap is being painted on"
- `https://doc.qt.io/qt-6/qpixmap.html#save` — documentation silent on PNG determinism (a negative finding, valuable)
- `https://imagemagick.org/script/compare.php` — `compare` exit status conventions (0 = similar)

### Secondary (verified against local codebase)
- `/home/drico/git/mefisto/xvue/xvuelc.c:1307–1428` — legacy 7 entry-point bodies (literal read)
- `/home/drico/git/mefisto/xvue/xvuelc.c:138, 1036` — single `mempxsauvfen` allocation (D-01 support)
- `/home/drico/git/mefisto/xvue/qt/src/xvue_qt_api.cpp:394–447` — current 7 warn-once stubs (literal read)
- `/home/drico/git/mefisto/xvue/qt/src/xvue_qt_api.cpp:56–75` — `xvue_qt_draw_rect_common` pattern (Phase 4 helper template)
- `/home/drico/git/mefisto/xvue/qt/src/xvue_qt_api.cpp:540–640` — `xvfermer_` capture hook (validation harness anchor)
- `/home/drico/git/mefisto/xvue/qt/src/xvue_qt_canvas.cpp:30–113` — `paintEvent` + `resizeEvent` (temp-painter pattern reference at 69–75)
- `/home/drico/git/mefisto/xvue/qt/src/xvue_qt_state.{h,cpp}` — `XvueState` structure + destructor ordering (D-03 constraint)
- `/home/drico/git/mefisto/prpr/xvtest0.f` — 168 lines, existing DRAW + TEXT coverage sections (Phase 4 extension target)
- `/home/drico/git/mefisto/bin/qt-capture.sh` — Qt offscreen capture harness (80 lines, unchanged)
- `/home/drico/git/mefisto/bin/xvtest-capture.sh` — legacy X11 Xvfb capture harness (140 lines, unchanged)
- `/home/drico/git/mefisto/bin/cbxvtest0_qt` — Qt driver build script (68 lines, unchanged)
- Local probe: `which convert magick cmp` → all three installed; `dpkg -l libgfortran5` → `15.2.0-9` pinned

### Tertiary (assumed — flagged in Assumptions Log)
- ImageMagick AE metric ignores PNG metadata in practice (A1)
- Qt offscreen capture is run-to-run deterministic for the same input (A2)
- gfortran `GETENV` intrinsic works here (A3)
- Legacy X11 and Qt produce visually-equivalent-but-not-pixel-identical captures (A5)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — no new dependencies, everything verified at research time
- Architecture (D-07 temp painter pattern): HIGH — verified against Qt 6 official docs + proven by the existing `resizeEvent` pattern at `xvue_qt_canvas.cpp:69–75`
- Pitfalls: HIGH — 5 of 7 are direct citations from doc.qt.io; 2 are direct readings of CONTEXT.md Phase 2 invariants
- Validation Architecture: HIGH — 100% reuse of existing Phase 03-04 harnesses, verified by direct read
- Assumptions (A1–A5): explicitly flagged; none are load-bearing for correctness, only for ergonomics

**Research date:** 2026-04-13
**Valid until:** Phase 4 execution completes OR Qt 6 moves to a major version bump (Qt 7) — whichever is sooner. All citations are against `doc.qt.io/qt-6/`, which tracks the current 6.x minor release; the documented invariants (one painter per device, drawPixmap logical coordinates, fill-while-painted-on undefined) are Qt 6 architectural, not release-specific.
