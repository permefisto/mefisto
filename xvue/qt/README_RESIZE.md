# xvue-qt — Canvas resize-preserve convention (DRAW-09)

**Authority:** 02-CONTEXT.md §D-07, §D-08; DRAW-09 (REQUIREMENTS.md).

When `XvueCanvas::resizeEvent` fires, the backing `QPixmap` is
reallocated and the previous drawing is preserved with a **top-left
anchor, no scaling, no centering, no clip compensation** rule:

1. The new backing is allocated at
   `widget.size() * devicePixelRatioF()` device pixels and its DPR is
   set via `setDevicePixelRatio(dpr)` so `QPainter` on it accepts
   logical coordinates unchanged (SHELL-06 invariant).
2. The new backing is first filled with `state_->background_` (the
   Phase 1 default is `Qt::black`).
3. The old backing's pixels are copied into the new backing via a
   single `QPainter(new_backing).drawPixmap(0, 0, *old_backing)` call.
   This anchors the old content at the top-left corner.
4. On **grow**, freshly uncovered area to the right and bottom shows
   through as `state_->background_` (step 2).
5. On **shrink**, content outside the new rectangle is clipped by
   Qt's `drawPixmap` automatically — no manual clip math is done.

No alternative anchor (center, scale-to-fit, etc.) is ever applied in
Phase 2. If a future phase wants a different behavior, it must edit
both the implementation in `xvue/qt/src/xvue_qt_canvas.cpp`
(`XvueCanvas::resizeEvent`) and this document in the same commit.

The long-lived `QPainter*` (D-05) is ended on the old backing and
re-`begin()`-ed on the new backing in the same sequence; render hints
and pen/brush are re-applied via `XvueState::applyPen()` because
Qt does not inherit them across `begin()`/`end()` cycles.

## References

- `.planning/phases/02-drawing-primitives-backing-pixmap/02-CONTEXT.md` §D-07, §D-08
- `xvue/qt/src/xvue_qt_canvas.cpp` (`XvueCanvas::resizeEvent`)
- `xvue/README_COORDS.md` (Phase 0 sibling artifact — Y-axis convention)
- REQUIREMENTS.md DRAW-09
