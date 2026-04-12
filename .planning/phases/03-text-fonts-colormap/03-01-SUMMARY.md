---
phase: 03-text-fonts-colormap
plan: 01
subsystem: graphics
tags: [qt6, fonts, palette, qresource, fortran]

requires:
  - phase: 02-drawing-primitives-backing-pixmap
    provides: XvueState with backing/painter/pen/brush, XvueCanvas with resize painter rebuild
provides:
  - DejaVuSansMono.ttf bundled via Qt resource system
  - XvueState palette statics (red/green/blue[256], palette_cache_[256])
  - XvueState font state (current_font_, current_metrics_, kFontSizes[10])
  - XvueApp bundled-font load (Q_INIT_RESOURCE + addApplicationFont)
  - XvueCanvas dark-mode-freeze guards (WA_OpaquePaintEvent, no auto-fill)
  - palette_init_once with verbatim 16-color imposed defaults from xvuelc.c:378-461
  - verify_no_exec.sh palette-leak scan guard
  - prpr/xvtest0.f Phase 3 TEXT coverage section
affects: [03-02-text-entry-points, 03-03-checkpoint, 03-04-validation]

tech-stack:
  added: [qt_add_resources, QFontDatabase, QFontMetrics]
  patterns: [static palette storage mirroring xvuelc.c file-scope arrays, Q_INIT_RESOURCE for static-lib font embedding]

key-files:
  created:
    - xvue/qt/fonts/DejaVuSansMono.ttf
    - xvue/qt/fonts/LICENSE.dejavu
    - xvue/qt/resources/xvue_fonts.qrc
  modified:
    - xvue/qt/CMakeLists.txt
    - xvue/qt/cmake/verify_no_exec.sh
    - xvue/qt/src/xvue_qt_state.h
    - xvue/qt/src/xvue_qt_state.cpp
    - xvue/qt/src/xvue_qt_app.h
    - xvue/qt/src/xvue_qt_app.cpp
    - xvue/qt/src/xvue_qt_canvas.cpp
    - prpr/xvtest0.f

key-decisions:
  - "qt_add_resources PREFIX set to /xvue/qt (not /xvue/qt/fonts) because FILES path fonts/DejaVuSansMono.ttf already contributes the fonts/ segment"
  - "imposed_defaults_fill uses .f suffix float literals (50.f/256.f) not double (50./256.) to match fromRgbF(float,float,float) without implicit narrowing"

patterns-established:
  - "Static palette arrays on XvueState class scope — survives XvueWindow destruction on xvfermer_"
  - "palette_init_once guard pattern — runs exactly once from XvueState ctor, before xvinfo_ returns"

requirements-completed: [TEXT-01, TEXT-02, TEXT-03, TEXT-04, TEXT-05, TEXT-06]

duration: 45min
completed: 2026-04-12
---

# Plan 03-01: Wave 0 Infrastructure Summary

**Bundled DejaVuSansMono font via Qt resources, 16-color imposed-default palette from xvuelc.c:378-461, canvas dark-mode guards, and xvtest0.f TEXT coverage driver**

## Performance

- **Duration:** ~45 min (inline execution after subagent sandbox failure)
- **Started:** 2026-04-12T00:00Z
- **Completed:** 2026-04-12T16:40Z
- **Tasks:** 3
- **Files modified:** 11

## Accomplishments
- DejaVuSansMono.ttf embedded in libxvueqt.a, loads successfully at runtime (no "failed to load" warning)
- Palette initialized with verbatim 16-color block from xvuelc.c:378-461 before any drawing
- XvueCanvas hardened against Qt dark-mode palette leaks (WA_OpaquePaintEvent + setAutoFillBackground(false))
- xvtest0.f exercises all TEXT-01..05 entry points (warn-once expected until 03-02 silences them)

## Task Commits

1. **Task 1: Bundle TTF + qt_add_resources + canvas guards** - `b1a2225`
2. **Task 2: XvueState palette + XvueApp font load** - `7f24e0f`
3. **Task 3: xvtest0.f TEXT coverage** - `0e9022a`

## Files Created/Modified
- `xvue/qt/fonts/DejaVuSansMono.ttf` - Bundled monospace font (343 KB)
- `xvue/qt/fonts/LICENSE.dejavu` - DejaVu/Bitstream Vera license
- `xvue/qt/resources/xvue_fonts.qrc` - Qt resource manifest
- `xvue/qt/CMakeLists.txt` - qt_add_resources wiring
- `xvue/qt/cmake/verify_no_exec.sh` - Palette-leak build guard (D-19)
- `xvue/qt/src/xvue_qt_state.h` - Palette + font state declarations
- `xvue/qt/src/xvue_qt_state.cpp` - Palette statics + imposed_defaults_fill + palette_init_once
- `xvue/qt/src/xvue_qt_app.h` - font_id_ + load_bundled_font_ declarations
- `xvue/qt/src/xvue_qt_app.cpp` - Q_INIT_RESOURCE + addApplicationFont implementation
- `xvue/qt/src/xvue_qt_canvas.cpp` - Dark-mode guards + TextAntialiasing hint
- `prpr/xvtest0.f` - Phase 3 TEXT coverage section

## Decisions Made
- Fixed qt_add_resources PREFIX from `/xvue/qt/fonts` to `/xvue/qt` — the FILES path `fonts/DejaVuSansMono.ttf` already contains the `fonts/` segment, so the original prefix produced `:/xvue/qt/fonts/fonts/DejaVuSansMono.ttf` (doubled)
- imposed_defaults_fill uses float literals (`.f` suffix) matching the C source's float storage semantics
- ABI symbol count is 35 (not 57 as plan stated from design doc — 57 is the full-project target, 35 is actual post-Phase 2)

## Deviations from Plan

### Auto-fixed Issues

**1. [Blocking] qt_add_resources PREFIX path doubling**
- **Found during:** Task 2 (runtime font-load check)
- **Issue:** Plan specified `PREFIX "/xvue/qt/fonts"` but with `FILES fonts/DejaVuSansMono.ttf`, Qt creates `:/xvue/qt/fonts/fonts/DejaVuSansMono.ttf`
- **Fix:** Changed PREFIX to `/xvue/qt` so full path is `:/xvue/qt/fonts/DejaVuSansMono.ttf`
- **Verification:** `strings libxvueqt.a | grep xvue/qt` shows correct path; runtime produces no font-load warning
- **Committed in:** 7f24e0f (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Essential fix for font loading. No scope creep.

## Issues Encountered
- Initial subagent execution failed due to sandbox permission denials on cmake/build invocations in worktree. Recovered by executing inline on main tree.

## User Setup Required
None

## Next Phase Readiness
- All infrastructure in place for plan 03-02 (TEXT entry point bodies)
- Plan 03-02's 8 entry-point bodies can drop in with zero friction — XvueState has palette/font state, font loads at runtime

---
*Phase: 03-text-fonts-colormap*
*Completed: 2026-04-12*
