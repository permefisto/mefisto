---
phase: 04-pixmap-save-restore-double-buffering
reviewed: 2026-04-14T00:00:00Z
depth: standard
files_reviewed: 5
files_reviewed_list:
  - bin/xvtest0-pixmap-roundtrip.sh
  - prpr/xvtest0.f
  - xvue/qt/src/xvue_qt_api.cpp
  - xvue/qt/src/xvue_qt_state.cpp
  - xvue/qt/src/xvue_qt_state.h
findings:
  critical: 0
  warning: 2
  info: 5
  total: 7
status: issues_found
---

# Phase 4: Code Review Report

**Reviewed:** 2026-04-14
**Depth:** standard
**Files Reviewed:** 5
**Status:** issues_found

## Summary

Phase 4 adds QPixmap-backed save/restore to the Qt backend via a single slot
(`XvueState::saved_canvas_`), routes `sauvefenetre_`/`sauvemempx_` and
`restaurefenetre_`/`restauremempx_` through two file-local helpers, and keeps
`fenetremempx_`/`mempxfenetre_` as documented no-ops (single-backing collapse
from Phase 2 D-04). The Fortran headless driver `xvtest0.f` gains a
`MEFISTO_XVTEST0_SCENE`-gated coverage block, and a bash harness drives six
captures through `bin/qt-capture.sh` and four pairwise `magick compare -metric AE`
checks.

Overall the design is sound: the scoped temporary `QPainter` on `saved_canvas_`
respects the Phase 2 D-05 single-long-lived-painter invariant, the destructor
deletes `saved_canvas_` before the active painter (D-03 "least-entangled
first"), size-mismatch on restore is an auditable stderr no-op (D-12), and the
Fortran scene dispatch correctly proves each of PIXMAP-01/02/03a/03b against
the same base scene. No critical bugs, no security issues. Two warnings are
concrete correctness risks in the test harness; five info items are code-quality
observations that don't block Phase 4 green.

## Warnings

### WR-01: Shell AE parser collapses multi-line stderr into a broken scalar

**File:** `bin/xvtest0-pixmap-roundtrip.sh:65-67`
**Issue:** `magick compare -metric AE` is invoked with `2>&1 || true`, so `$AE`
captures any progress/warning lines that IM7 may prepend (e.g. delegate
warnings, `@warning/...` lines from policy.xml, or the IM7 "0 (0)" normalized
form on multi-line output). `COUNT=$(echo "$AE" | awk '{print $1}')` then
prints `$1` for **every** input line — if stderr has any warning line at all,
`$COUNT` becomes a multi-line string and `[ "$COUNT" = "0" ]` silently fails
even when AE is zero, producing a false FAIL. The `|| true` also hides IM7
exit code 2 (real error) as if it were exit 1 (images differ), so a missing
PNG will be reported as a pixel diff rather than a harness failure.
**Fix:**
```bash
# Separate stdout/stderr; take only the first line of the metric output;
# distinguish "differ" (rc=1) from "error" (rc>=2).
AE=$(magick compare -metric AE "$A" "$B" null: 2>&1 >/dev/null)
rc=$?
if [ "$rc" -ge 2 ]; then
    echo "FAIL  [${REQ}]  ${LABEL}  (magick error rc=${rc}: ${AE})"
    FAILED=$((FAILED + 1))
    return
fi
COUNT=$(printf '%s\n' "$AE" | awk 'NR==1 {print $1; exit}')
if [ "$COUNT" = "0" ]; then
    ...
```

### WR-02: `saved_canvas_` devicePixelRatio not refreshed on slot reuse

**File:** `xvue/qt/src/xvue_qt_api.cpp:94-98`
**Issue:** The lazy-realloc branch (`!st->saved_canvas_ || size mismatch`) sets
`setDevicePixelRatio(backing_->devicePixelRatio())`, but the reuse path (size
matches) does not. If the window is dragged between monitors with different
DPR while keeping logical size constant (Qt can change `backing_->devicePixelRatio()`
without changing `backing_->size()` — size is in device pixels on Qt 6), a
subsequent `sauvefenetre_` copies device-pixel content into a slot whose DPR
no longer matches, and `restaurefenetre_`'s `drawPixmap(0,0,...)` will then
scale incorrectly. The round-trip test never exercises DPR changes, so the
current PIXMAP-02 pass does not cover this. Low-likelihood under xvtest0 but
a real correctness bug under interactive use (Phase 5 cavity2d rubber-band).
**Fix:**
```cpp
if (!st->saved_canvas_ || st->saved_canvas_->size() != st->backing_->size()) {
    delete st->saved_canvas_;
    st->saved_canvas_ = new QPixmap(st->backing_->size());
}
// Always refresh DPR — cheap and idempotent.
st->saved_canvas_->setDevicePixelRatio(st->backing_->devicePixelRatio());
```

## Info

### IN-01: `effacemempx_` and `effacer_` are byte-identical bodies

**File:** `xvue/qt/src/xvue_qt_api.cpp:498-524`
**Issue:** Comment at 494-497 acknowledges the duplication is intentional for
ABI preservation (D-08), but the bodies will drift if either is edited in
isolation (e.g. a future D-15 update to `effacer_` that does not replicate to
`effacemempx_`). PIXMAP-03b would then silently regress.
**Fix:** Factor into a file-local helper in the anonymous namespace, same
pattern as `xvue_qt_save_to_slot` / `xvue_qt_draw_rect_common`:
```cpp
inline void xvue_qt_clear_backing() {
    auto& win = XvueApp::window_slot();
    if (!win) return;
    auto* st = win->state();
    if (st && st->painter_ && st->painter_->isActive() && st->backing_) {
        st->painter_->fillRect(st->backing_->rect(), st->background_);
    }
    if (win->canvas()) win->canvas()->update();
}
```
Both entry points become one line. ABI count stays at 57.

### IN-02: Scene-name comparison is case-sensitive but script always uppercases

**File:** `prpr/xvtest0.f:53-112`
**Issue:** `IF (SCENE .EQ. 'P4_CTRL')` etc. require exact-case input. A user
invoking `MEFISTO_XVTEST0_SCENE=p4_ctrl ./pp/ppxvtest0_qt` manually will fall
through to the unknown-scene path and then into the legacy Phase 1/2 code path.
The harness always uppercases, so this is only a hand-run footgun. Not a bug,
but consider normalizing once after `GETENV` for defensive UX.
**Fix:** Document the requirement in the comment block at line 44-50, or
uppercase `SCENE` locally before dispatch.

### IN-03: Bash script uses `tr` subshell for ASCII lowercasing

**File:** `bin/xvtest0-pixmap-roundtrip.sh:42`
**Issue:** `LC=$(echo "$SCENE" | tr 'A-Z' 'a-z')` forks two processes per scene.
The shebang is `#!/bin/bash`, so the parameter-expansion form is available.
**Fix:** `LC="${SCENE,,}"`

### IN-04: Hand-rolled `compare_ae` uses `local` without `typeset`/POSIX fallback

**File:** `bin/xvtest0-pixmap-roundtrip.sh:59-73`
**Issue:** `local` is bash-only, which matches the shebang, but the variable
`FAILED` is modified from inside the function via `$((FAILED + 1))` — that
works only because `FAILED` is not declared `local`. Readers may assume it is.
Declare the mutation explicitly or use a dedicated accumulator function.
**Fix:** Add a comment at line 58 noting `FAILED` is intentionally global, or
return a non-zero rc from `compare_ae` and aggregate at the call site.

### IN-05: `xvue_qt_save_to_slot` silently no-ops when backing_ is null

**File:** `xvue/qt/src/xvue_qt_api.cpp:86-107`
**Issue:** Unlike `restore_from_slot` (which warns on size mismatch, D-12),
`save_to_slot` returns silently if any precondition fails (null painter,
inactive painter, null backing_). If xvtest0 ever drifts so that `sauvefenetre_`
runs before `xvinitgraphique_` has allocated `backing_` (e.g. a future scene
reorders calls), the save will no-op, the subsequent restore will hit the
size-mismatch branch, and the test will print a confusing "no slot or size
mismatch" stderr line with no indication that the save never happened.
**Fix:** Add a symmetric stderr warning on the save precondition failure:
```cpp
if (!st || !st->painter_ || !st->painter_->isActive() || !st->backing_) {
    std::fprintf(stderr, "xvue-qt: save_to_slot: painter/backing not ready\n");
    return;
}
```
Cheap diagnostic, zero cost on the happy path.

---

_Reviewed: 2026-04-14_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
