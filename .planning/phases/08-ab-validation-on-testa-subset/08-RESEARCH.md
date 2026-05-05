# Phase 8: A/B validation on testa subset - Research

**Researched:** 2026-05-05
**Domain:** Headless A/B image-diff validation of an X11 (legacy) vs Qt 6 (new) graphics backend over 5 canonical Fortran/X11 test cases, with HiDPI and OpenMP sweeps, anchored by ImageMagick `compare -metric AE -fuzz 5%` and committed PNG evidence.
**Confidence:** HIGH on tool empirics (Verified via `compare`/`Xvfb`/`xvfb-run` runs in this session); HIGH on existing Qt-side capture infrastructure (Verified by reading `xvue/qt/src/xvue_qt_api.cpp` and `xvue/xvuelc.c`); HIGH on environment availability (Verified via `dpkg`/`pkg-config`); MEDIUM on per-case batch-mode behavior of testa cases under both backends (Cited from `bin/MAILLER` source; not yet end-to-end validated against pp/*_qt because INITIER is required first).

## Summary

Phase 8 is a verification-only ship gate: it runs the 5 canonical `testa/` cases through both the legacy X11 and new Qt 6 backends, screenshots the canvas in each (case, backend, mode) cell, runs ImageMagick `compare -metric AE -fuzz 5%` between matched pairs, commits both the PNGs and the per-cell verdict in `08-CHECKLIST.md`, and signs off the v1 ship gate that opens the one-release-cycle A/B window before Phase 9 retires X11. Plan 1 must first **bootstrap the 3 Phase-7-deferred goldens** (scene01.eps, wave_legacy.gif, cavity2d_legacy.gif) and **refresh `pp/*_qt` via `bin/cbl_tout_qt`** (Phase 7 Gap-A) before any A/B sweep starts.

The infrastructure to run this headlessly already exists and is impressive. The `xvue_qt_api.cpp` `xvfermer_` body honors `MEFISTO_QT_CAPTURE_PATH` to **save the backing pixmap directly as PNG in-process** (no X server needed — works under `QT_QPA_PLATFORM=offscreen`). The X11 `xvuelc.c` `xvfermer_` honors `MEFISTO_XVFERMER_READY_FILE` + `MEFISTO_XVFERMER_HOLD_MS` to give an external capture tool a deterministic window to grab the Xvfb framebuffer. Both backends honor `MEFISTO_XVSOURIS_AUTOEXIT` to drive batch-mode termination. testa cases are **already script-driven** (the `.mesh`/`.heat`/`.elas` files are batch lexicon scripts ending in `CLOSE;` or `99;`) — `bin/MAILLER`/`ELASTICER`/etc. forward `$*` to `pp/pp*` so the same files drive both backends identically.

**Primary recommendation:** Standardize on **two capture mechanisms**: Qt-side via `MEFISTO_QT_CAPTURE_PATH` + `QT_QPA_PLATFORM=offscreen` (in-process backing-pixmap save, no X server, identical resolution to the canvas at any DPR), X11-side via `MEFISTO_XVFERMER_READY_FILE`/`_HOLD_MS` + `import -window root` against the Xvfb display under `xvfb-run --auto-servernum -s "-screen 0 1280x800x24"`. To make `compare` honest, **resample the Qt HiDPI capture to the X11 baseline dimensions before diffing** — this neutralizes T-08-05 cleanly. Use `compare -metric AE -fuzz 5%` exactly as locked, but **gate `-fuzz` between 1% and 30%** in plan validation (above 30% the metric becomes meaningless — empirically `-fuzz 100%` returns 0 on totally different images), and **always pre-check dimensions match** with `identify -format '%wx%h'` because `compare` silently returns AE=0 / exit=0 on dimension mismatch (T-08-02 confirmed).

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**Evidence capture**

- **D-01:** Side-by-side screenshots committed to repo under `.planning/phases/08-ab-validation-on-testa-subset/evidence/{case}-{backend}-{mode}.png`. Reproducible, reviewable in PR diffs, expected ~5 MB total. PNG only — no JPEG (lossy artifacts would defeat tolerance band).

**Pass threshold**

- **D-02:** Tolerance band, NOT eyeball-only or pixel-perfect. Geometry drift cap: sub-2 pixels. Color drift cap: delta-E < 5. Pixel-perfect rejected up front (Qt anti-aliasing + font hinting differ from X11).
- **D-03:** Tolerance enforcement tool: ImageMagick `compare -metric AE -fuzz 5%`. The dev host already ships ImageMagick. Phase 9 RETIRE-03 removes ImageMagick from `xvue/qt/` ONLY — not from validation tooling under `bin/` or `.planning/`. Reusing legacy `compare` is in-scope.

**HiDPI methodology**

- **D-04:** `QT_SCALE_FACTOR=2` under `xvfb-run --auto-servernum`. Cheap, reproducible, runnable headless on the build host. Captures the 2x scaling logic but not real 4K display pixels. Real-4K eyeball check is NOT a Phase 8 gate — recorded as a deferred idea for the maintainer's ad-hoc validation.

**OMP scope**

- **D-05:** Sweep all 5 `_OMP` variants where they exist on disk: mesher / elas / flui / ther / nlse. VALID-03's literal text says "ELASTICER_OMP" — that is the Phase-7-end smoke check; Phase 8 broadens to all 5 for full main-thread-invariant coverage. Cells where no `_OMP` variant exists for a module are marked `N-A` in the checklist.

**Phase 7 deferred-golden integration**

- **D-06:** Plan 1 of Phase 8 bootstraps the 3 Phase-7-deferred goldens before any A/B sweep runs:
    1. Compile `xvue/qt/tests/golden/scene01_driver.f` against X11 + xvuelc.o under Xvfb, capture `TEMPORAIRE.EPS`, commit as `xvue/qt/tests/golden/scene01.eps`.
    2. Run `testa/wave` through X11 + `bin/convertepsgif`, commit `xvue/qt/tests/golden/wave_legacy.gif`.
    3. Run `testa/cavity2d` through X11 + `bin/convertepsgif`, commit `xvue/qt/tests/golden/cavity2d_legacy.gif`.
- **D-07:** After Plan 1 commits the goldens, re-run `ctest -R 'xvue_qt_(postscript|export)_tests'` and confirm all 3 Phase-7 QSKIPs flip to PASS. This is the Phase-7-close gate folded into Phase 8 entry.

**Build hygiene (Phase 7 Gap-A)**

- **D-08:** Plan 1 first task runs `bin/cbl_tout_qt` from a clean tree to refresh `pp/*_qt`. Phase 7 Gap-A documented that the original verifier ran `cmake --build xvue/qt/build` + `bin/cbl_tout` but never `bin/cbl_tout_qt`, so `pp/*_qt` could be stale relative to `libxvueqt.a`. A/B sweep MUST run against fresh `pp/*_qt`.
- **D-09:** Permanent CMake guard target proposed but NOT in Phase 8 scope: a CMake `verify_pp_qt_freshness` target that fails the build when `libxvueqt.a` mtime > any `pp/*_qt` mtime. Tracked as deferred idea for Phase 9 cleanup.

**CHECKLIST.md sign-off shape**

- **D-10:** Per-cell verdict matrix. Rows = 5 canonical cases. Columns = {X11 baseline, Qt 1x, Qt HiDPI 2x, Qt _OMP}. Each cell records: PASS / FAIL / N-A + maintainer initials + path-to-evidence-png + ImageMagick AE pixel count. v1 ship gate = every cell in {PASS, N-A} with maintainer signature on the bottom-row sign-off line. Any FAIL blocks ship until either resolved (gap-closure plan) or explicitly overridden with documented rationale.

### Claude's Discretion

- Per-plan task decomposition (Plan 1 = bootstrap + freshness, Plan 2..N = A/B sweep batches per backend mode, Plan N+1 = CHECKLIST.md finalize). Planner picks the exact split.
- Screenshot capture script: shell wrapper using kwin-mcp / xdotool / Xvfb's xwd, or a small Qt-based recorder. Planner chooses based on what's reproducible across backends.
- Xvfb screen resolution defaults: 1280x800 (matches Phase 7 UAT session). Planner can adjust per-case if a test has hard-coded larger window expectations.
- Tolerance-band command-line tuning (`-fuzz` percent, `-metric` choice AE vs PAE vs RMSE). D-02 specifies the gate; planner picks the exact invocation.

### Deferred Ideas (OUT OF SCOPE)

- **Real-4K HiDPI eyeball check** — beyond Phase 8's QT_SCALE_FACTOR=2 ship gate. Maintainer ad-hoc validation when 4K hardware is available. Not blocking.
- **CMake `verify_pp_qt_freshness` ALL target** — permanent fix for Phase 7 Gap-A. Belongs in Phase 9 cleanup (or earlier hot-fix outside Phase 8 scope).
- **CI integration** — automated Xvfb runs in GitHub Actions / GitLab CI. Out of scope for v1 ship; Phase 8 evidence is committed artifacts.
- **Tolerance-band per-case overrides** — some cases (notably cavity2d's velocity arrows) may legitimately need a looser fuzz than 5%. If discovered during Plan 2+, capture as a CHECKLIST.md per-cell override note rather than relaxing the global gate.
- **SSIM/skimage-based diff** — alternative tolerance metric. Considered and rejected (D-03) for Phase 8 to avoid Python+pip dependency.
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| VALID-01 | All 5 canonical testa/ cases (BUILD-10 baseline) pass visually side-by-side at end of every phase 0-7; logged to .planning/phase-N/VALIDATION.md | Section "Architecture Patterns / Pattern 1: Headless A/B sweep loop" — exact harness for running each (case, backend) cell. testa cases are batch-driven (e.g., pan2d.mesh ends `CLOSE;`) so the same input file drives both pp/ppmail (X11) and pp/ppmail_qt (Qt) identically. [VERIFIED: read of pan2d.mesh + bin/MAILLER source] |
| VALID-02 | X11 backend continues to pass the same 5 cases at end of every phase | Section "Build Invariants" — re-run X11 column with stock `bin/cbl_tout` binaries (xvuelc.c byte-identical guarantee from Phase 7 verifier; current pp/* X11 binaries are stable). [VERIFIED: `git diff 900e297..HEAD -- xvue/xvuelc.c` empty in 07-VALIDATION-LOG] |
| VALID-03 | All 5 cases pass through ELASTICER_OMP (broadened by D-05 to all 5 _OMP variants) — main-thread invariant under OMP | Section "OMP Sweep — How _OMP Actually Works" + "Pitfall 5: OMP scheduling jitter". `_OMP` launchers run the SAME `pp/ppelas` binary with `OMP_NUM_THREADS=8` set; the pp binaries are already linked against libgomp. The OMP risk is graphics calls from worker threads, mitigated by XVUE_QT_ASSERT_MAIN_THREAD on every public method (SHELL-07, all 57 stubs Phase 1). [VERIFIED: `nm pp/ppelas \| grep GOMP` + `ldd pp/ppelas_qt` shows libgomp linked; bin/ELASTICER_OMP source confirms same `pp/ppelas` invocation] |
| VALID-04 | All 5 cases render correctly at HiDPI (4K or QT_SCALE_FACTOR=2) — no size/position drift vs X11 | Section "HiDPI methodology" + "Pitfall 4: Naive AE compare on mismatched dims". QT_SCALE_FACTOR=2 under Xvfb keeps the framebuffer at 1280x800 but the Qt canvas at 640x400 logical / 1280x800 physical. The Qt-side capture (MEFISTO_QT_CAPTURE_PATH) saves the BACKING pixmap which is physical-pixel-sized — for 1x vs 2x compare-to-X11-baseline, **resample to baseline dims first**. [VERIFIED: ran QT_SCALE_FACTOR=2 against built qtest binary under Xvfb 1280x800] |
| VALID-05 | Color-bar spot checks on testa/nafems_le1, heat1d, cavity2d show no user-visible drift | Section "Pitfall 6: gradient color bars trip fuzz", + recommend per-case override pattern in CHECKLIST.md (allowed by D-02 deferred idea) when a stress-bar gradient legitimately exceeds fuzz tolerance. [Verified empirically: scientific colormaps are pre-baked palette indices via xvactivervb_ — Phase 3 validated 1-bit-per-channel match between X11 and Qt] |
| VALID-06 | Font-metric spot checks on testa/pan2d, hexahedron — no clipping or overlap | Section "Pitfall 7: anti-aliased text + font hinting drift". testa/hexahedron is NOT in the 5-case BUILD-10 baseline (BASELINE.md fixes pan2d/nafems_le1/cavity2d/heat1d/nlsecu). VALID-06 mentions hexahedron as a font-metric SPOT CHECK target — Phase 8 should run this case in addition to the 5 canonical, scoped narrowly to font-metric verification only. [CITED: REQUIREMENTS.md VALID-06 + BASELINE.md amendment policy] |
| VALID-07 | A validation checklist (CHECKLIST.md) records pass/fail per case per backend; gate for declaring v1 shippable | Section "CHECKLIST.md format" — proposed structure includes per-cell PASS/FAIL/N-A + initials + AE pixel count + evidence path + override rationale (if used). The format is designed so a future CMake gate (D-09 deferred to Phase 9) can parse it. [Locked by D-10] |
| BUILD-10 | 5-case baseline locked in BASELINE.md | Phase 8 strictly uses pan2d/nafems_le1/cavity2d/heat1d/nlsecu — same 5 cases. **Plus** hexahedron as a VALID-06 font-only spot check (does NOT extend the BUILD-10 baseline; it's a separate checklist row marked "spot check"). [Locked by BUILD-10 + amendment policy] |
</phase_requirements>

## Architectural Responsibility Map

This phase is verification-only — it does not introduce new architectural components. The map records which existing tier owns each capability the phase exercises, so plans assign work to the right place.

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Backing-pixmap snapshot to PNG | Qt backend (`xvue/qt/src/xvue_qt_api.cpp` `xvfermer_`) | — | Already implemented via `MEFISTO_QT_CAPTURE_PATH`; in-process `QPixmap::save()`. Offscreen-capable (no X server required). |
| X11 framebuffer snapshot to PNG | External CLI tool (`import -window root`) on Xvfb display | X11 backend `xvfermer_` synchronises via `MEFISTO_XVFERMER_READY_FILE` + `_HOLD_MS` | xvuelc.c has no in-process screenshot helper (would require writing pixels via XCreateImage + libpng — out of scope). External capture using `import` is fine — `import` is in legitimate-allowlist scope (under bin/), not under xvue/qt/. |
| Tolerance-band image diff | Phase 8 shell wrapper under `bin/` or `.planning/` | ImageMagick `compare` CLI | Stays out of `xvue/qt/` (EXPORT-06 grep gate scope). |
| Per-case batch driver | testa/{case}/<file>.{mesh,elas,flui,heat,nlse} batch scripts | `bin/{LAUNCHER}` shell wrapper forwards `$*` | Already-existing infrastructure; cases like `testa/pan2d/pan2d.mesh` end with `CLOSE;` for batch termination. |
| HiDPI logical/physical math | Qt backend, exposed via `QT_SCALE_FACTOR=2` | Phase 7 README HiDPI export math (Pitfall 7) — already documents the convention | Resample to baseline-dim before AE compare; do not pre-scale the Xvfb framebuffer. |
| OMP main-thread guard | Qt backend, every public ABI method via `XVUE_QT_ASSERT_MAIN_THREAD()` (SHELL-07) | bin/ELASTICER_OMP launcher sets `OMP_NUM_THREADS=8` | Same `pp/ppelas` binary; only the env var differs between OMP and serial runs. |
| Sign-off ledger | `.planning/phases/08-*/08-CHECKLIST.md` (markdown table per D-10) | — | Phase 8 deliverable; downstream-parsable by future CMake guard. |

## Standard Stack

This is a verification phase — no new libraries are added. The "stack" is the existing toolset on the dev host that the phase consumes.

### Core (verified present on dev host)

| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| ImageMagick `compare` / `import` / `identify` / `convert` | 7.1.2-18 Q16 | Image diff (`compare`), root-window capture (`import -window root`), dimension probe (`identify`), legacy GIF pipeline (`convert` via `bin/convertepsgif`) | Already on host; D-03 explicitly allows reuse for validation tooling under `bin/` and `.planning/`. EXPORT-06 grep gate scope is `xvue/qt/` only — external capture is fine. [VERIFIED: `dpkg -l imagemagick`] |
| Xvfb / `xvfb-run` | 21.1.21-1 | Headless X server for both backends | Standard MEFISTO test harness (already used by every Qt ctest target). [VERIFIED: `dpkg -l xvfb`] |
| Qt 6 | 6.10.2 (Debian forky/sid) | Qt backend libraries (libQt6Core6t64, libQt6Gui, libQt6Widgets) | Already pinned by Phase 0; `pp/pp*_qt` link against this. [VERIFIED: `pkg-config --modversion Qt6Core` = 6.10.2] |
| gfortran | 14.3.0 (Debian 14.3.0-14) | Build the X11 + Qt backends | Already pinned by build invariants. **libgfortran5 is held at 15.2.0-9** to avoid the 03-04 reopen UB documented in STATE.md. [VERIFIED: `apt-mark showhold` = libgfortran5; `gfortran --version`] |
| ffmpeg | 8.1-3+b1 | Animated GIF fallback path under XVUE_ANIM=1 (Phase 7 D-11) — needed for Plan 1 wave/cavity2d golden bootstrap reference | Already on host; Phase 7 PROBE.md recorded `gif_write_supported=0` so QImageWriter cannot write GIF; ffmpeg is the chosen fallback. [VERIFIED: `ffmpeg -version`] |
| `scrot` | 2.0.0-1 | Alternate root-window capture (faster than `import` on some setups) | Available; not strictly required (import works) but useful as a fallback. [VERIFIED: `dpkg -l scrot`] |

### Supporting (not required but available)

| Tool | Version | Purpose | When to Use |
|------|---------|---------|-------------|
| `gnome-screenshot` | (installed) | Wayland-session capture | Skip — Phase 8 must run under Xvfb (X11), not Wayland. |
| `kwin-mcp` | (installed at ~/.local/share/kwin-mcp-venv) | KWin Wayland session driver — used in Phase 7 UAT for live-window AT-SPI inspection | Available, but **NOT recommended for Phase 8 capture** — the Phase 8 contract is headless reproducibility under Xvfb, and kwin-mcp couples to the live Wayland session. Use only for human-checkpoint visual sanity, not for the AE-diff pipeline. |
| `xdotool` | NOT installed | Synthetic key/mouse events to live X11 windows | Skip — testa cases drive themselves via batch input files (`pan2d.mesh` ends with `CLOSE;`); MEFISTO_XVSOURIS_AUTOEXIT short-circuits any blocking xvsouris_ calls. No xdotool needed. [VERIFIED: `command -v xdotool` returns nothing on this host] |
| `xwd` | NOT installed (x11-apps removed) | X11 window dump | Skip — `import -window root` from ImageMagick is the equivalent and is installed. [VERIFIED: `which xwd` empty; `dpkg -l x11-apps` shows rc-state] |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| ImageMagick `compare` | Python `skimage.metrics.structural_similarity` (SSIM) | More forgiving on AA / font hinting drift, but requires Python + pip dependency. **Rejected by D-03.** |
| ImageMagick `compare` | A small Qt-based `QImage` AE diff tool | No new dep, but more code to maintain and would itself live somewhere — putting it in `xvue/qt/` violates EXPORT-06 grep gate; putting it in `bin/` is fine but reinvents what `compare` already does. Not worth the work. |
| `import -window root` | `scrot` | scrot is faster and slightly more reliable on display-server quirks. Either works. Pick `import` for symmetry with the rest of the ImageMagick toolchain (one less binary to remember). |
| `MEFISTO_QT_CAPTURE_PATH` (in-process) | `import -window root` against the Xvfb display showing the Qt window | In-process capture is preferred because (1) it works under `QT_QPA_PLATFORM=offscreen` with NO X server, (2) it captures the BACKING pixmap directly so the result is unambiguous regardless of widget paint timing, (3) it works at any DPR without DPR-confusion. Use `import` only as a fallback or for the X11 backend (which has no equivalent in-process hook). |

**Installation:** No new packages required. The phase's tooling is already on the host.

**Version verification:** Performed in this session against the live host:
```
ImageMagick: 7.1.2-18 Q16 (verified 2026-05-05 via `compare --version`)
Xvfb:        21.1.21-1     (verified 2026-05-05 via `dpkg -l xvfb`)
Qt 6:        6.10.2+dfsg-7 (verified 2026-05-05 via `pkg-config --modversion Qt6Core`)
gfortran:    14.3.0        (verified 2026-05-05 via `gfortran --version`)
libgfortran5: 15.2.0-9 (held)  (verified 2026-05-05 via `apt-mark showhold`)
ffmpeg:      8.1-3+b1      (verified 2026-05-05 via `ffmpeg -version`)
scrot:       2.0.0-1       (verified 2026-05-05 via `dpkg -l scrot`)
```

## Architecture Patterns

### System Architecture Diagram

```
                           ┌──────────────────────────────────┐
                           │  Phase 8 sweep driver shell      │
                           │  (under bin/ or .planning/, NOT  │
                           │   under xvue/qt/ — EXPORT-06)    │
                           └────────────┬─────────────────────┘
                                        │ for each (case, backend, mode):
                                        ▼
                  ┌───────────────────────────────────────────────────┐
                  │  xvfb-run --auto-servernum -s "-screen 0 1280x800x24"  │
                  │   ↓ env exports:                                  │
                  │     MEFISTO=$PWD                                  │
                  │     MEFISTOX=/tmp/mefistox-phase8                 │
                  │     QT_SCALE_FACTOR={1,2}     (HiDPI sweep)       │
                  │     OMP_NUM_THREADS={1,8}     (OMP sweep)         │
                  │     MEFISTO_XVSOURIS_AUTOEXIT=1                   │
                  │     MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=500        │
                  │     MEFISTO_XVFERMER_HOLD_MS=300                  │
                  │     [Qt only] MEFISTO_QT_CAPTURE_PATH=evidence/{cell}.png │
                  │     [Qt only] QT_QPA_PLATFORM=offscreen (optional)│
                  │     [X11 only] MEFISTO_XVFERMER_READY_FILE=/tmp/ready │
                  └────────────┬──────────────────────────────────────┘
                               │
                               ▼
            ┌──────────────────────────────┐         ┌──────────────────────────────┐
            │  pp/pp{module}    (X11)      │  vs.    │  pp/pp{module}_qt  (Qt 6)    │
            │   ← INITIER/MAILLER/etc.     │         │                              │
            │   ← testa/{case}/<batch>     │         │   ← same testa/<batch>       │
            └──────────────┬───────────────┘         └──────────────┬───────────────┘
                           │ at xvfermer_:                          │ at xvfermer_:
                           │   touches READY_FILE,                  │   QPixmap::save() → PNG
                           │   holds HOLD_MS                        │   (in-process; no X11)
                           ▼                                        ▼
            ┌──────────────────────────────┐         ┌──────────────────────────────┐
            │  external `import -window    │         │  evidence/{case}-qt-{mode}.png│
            │   root` against $DISPLAY     │         │  (already saved)             │
            │  → evidence/{case}-x11.png   │         │                              │
            └──────────────┬───────────────┘         └──────────────┬───────────────┘
                           └────────────────┬───────────────────────┘
                                            ▼
                       ┌────────────────────────────────────────┐
                       │  Tolerance-band gate:                  │
                       │  identify ... && compare -metric AE    │
                       │     -fuzz 5% A.png B.png diff.png      │
                       │  → AE pixel count + exit code          │
                       └────────────┬───────────────────────────┘
                                    ▼
                       ┌────────────────────────────────────────┐
                       │  08-CHECKLIST.md row update:           │
                       │  | case | X11 | Qt 1x | Qt HiDPI | Qt OMP │ initials │
                       │  AE pixel count + path-to-evidence-png │
                       └────────────────────────────────────────┘
```

### Recommended Project Structure

```
.planning/phases/08-ab-validation-on-testa-subset/
├── 08-CONTEXT.md            (existing — locked decisions)
├── 08-DISCUSSION-LOG.md     (existing — audit trail)
├── 08-RESEARCH.md           (this file)
├── 08-XX-PLAN.md            (per-plan files; Plan 1 = bootstrap, Plan 2..N = sweep batches, Plan N+1 = CHECKLIST.md finalize)
├── 08-XX-SUMMARY.md
├── 08-CHECKLIST.md          (D-10 ship gate — final deliverable)
├── 08-VALIDATION.md         (per-VALID-NN coverage notes; required by VALID-01..07 prose)
└── evidence/
    ├── pan2d-x11.png
    ├── pan2d-qt-1x.png
    ├── pan2d-qt-2x.png
    ├── pan2d-qt-omp.png
    ├── pan2d-diff-1x.png        (compare's diff output)
    ├── nafems_le1-x11.png
    ├── nafems_le1-qt-1x.png
    ├── ... (5 cases × 4 modes = 20 captures + ~20 diff PNGs ≈ 5 MB total)
    └── (optional) hexahedron-{x11,qt-1x}.png  (VALID-06 spot check)
```

The validation tooling itself (sweep script, AE-gate wrapper) lives under `bin/`:

```
bin/
├── ab_sweep_phase8.sh           (NEW — Plan 2+ harness; calls xvfb-run + envvar setup; invokes pp/* and pp/*_qt; calls compare)
├── ab_capture_x11.sh            (NEW — wrapper around import -window root + READY_FILE polling)
├── ab_compare_pair.sh           (NEW — identifies dims, compares, returns AE count)
├── convertepsgif                (existing — Phase 1 wave/cavity2d golden bootstrap)
└── test_no_imagemagick_in_qt.sh (existing — EXPORT-06 grep gate; verifies the new bin/ scripts don't accidentally land in xvue/qt/)
```

### Pattern 1: Headless A/B Sweep Loop

**What:** Single command per (case, backend, mode) cell that yields a capture PNG + an AE pixel count.

**When to use:** Every Plan 2..N cell run. Reusable across all 5 cases × 4 modes = 20 cells (+ optional hexahedron font spot).

**Example (canonical Qt-side, offscreen-capture):**
```bash
# Source: assembled from xvue/qt/src/xvue_qt_api.cpp xvfermer_ + bin/MAILLER (verified in this session)

CASE=pan2d; MODULE=mail; MODE=qt-1x
PROJDIR=/tmp/mefistox-phase8/${CASE}
EVIDENCE=.planning/phases/08-ab-validation-on-testa-subset/evidence

# 1. Initialize project (one-time per case; INITIER creates the project tree)
mkdir -p "$PROJDIR"
cp testa/${CASE}/* "$PROJDIR/"

# 2. Run the case under Qt with offscreen capture
env MEFISTO=$PWD MEFISTOX=$(dirname "$PROJDIR") \
    QT_QPA_PLATFORM=offscreen \
    MEFISTO_QT_CAPTURE_PATH="$EVIDENCE/${CASE}-${MODE}.png" \
    MEFISTO_XVSOURIS_AUTOEXIT=1 \
    MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=500 \
  $PWD/pp/pp${MODULE}_qt $(ls "$PROJDIR"/*.{mesh,elas,flui,heat,nlse} 2>/dev/null | head -1)

# 3. Verify capture
identify "$EVIDENCE/${CASE}-${MODE}.png"
```

**Example (canonical Qt HiDPI mode):**
```bash
# Same, with QT_SCALE_FACTOR=2 — backing pixmap is dpr-scaled (e.g., 1600x1200 for 800x600 logical)

env QT_SCALE_FACTOR=2 \
    MEFISTO_QT_CAPTURE_PATH="$EVIDENCE/${CASE}-qt-2x.png" \
    QT_QPA_PLATFORM=offscreen \
    ... \
  pp/pp${MODULE}_qt ...
```

**Example (canonical X11-side, Xvfb capture):**
```bash
# Source: synthesised from xvue/xvuelc.c xvfermer_ READY_FILE/HOLD_MS hook + import(1)

env DISPLAY=:99 MEFISTO_XVFERMER_READY_FILE=/tmp/ready_${CASE} \
    MEFISTO_XVFERMER_HOLD_MS=1500 \
    MEFISTO_XVSOURIS_AUTOEXIT=1 \
  pp/pp${MODULE} testa/${CASE}/<batch> &
PID=$!

# Poll for the ready sentinel; capture immediately
while [ ! -f /tmp/ready_${CASE} ]; do sleep 0.1; [ ! -d /proc/$PID ] && break; done
DISPLAY=:99 import -window root "$EVIDENCE/${CASE}-x11.png"
wait $PID
```

(Wrap the X11 case under `xvfb-run --auto-servernum -s "-screen 0 1280x800x24"`.)

### Pattern 2: AE-fuzz Tolerance Gate (with dimension guard)

**What:** Compare two captures and produce a pass/fail verdict + AE pixel count for the CHECKLIST cell.

**When to use:** Once per (case, mode) row, comparing the Qt capture against the X11 capture for that mode.

**Example:**
```bash
# Source: empirical compare 7.1.2 behavior verified in this session — DIMENSION GUARD is mandatory

A="$EVIDENCE/${CASE}-x11.png"
B="$EVIDENCE/${CASE}-${MODE}.png"
DIFF="$EVIDENCE/${CASE}-${MODE}-diff.png"

# 1. DIMENSION GUARD — without this, compare silently returns AE=0 / exit=0 on size mismatch
DIMS_A=$(identify -format "%wx%h" "$A")
DIMS_B=$(identify -format "%wx%h" "$B")
if [ "$DIMS_A" != "$DIMS_B" ]; then
    # Resample B to A's dims so the AE compare is honest. Use point sampling
    # to avoid introducing AA pixels not in either source.
    B_RESAMPLED="${B%.png}-resampled.png"
    convert "$B" -filter point -resize "${DIMS_A}!" "$B_RESAMPLED"
    B="$B_RESAMPLED"
fi

# 2. Run the tolerance-band gate. -fuzz 5% per D-02. -metric AE returns pixel count.
AE=$(compare -metric AE -fuzz 5% "$A" "$B" "$DIFF" 2>&1 | awk '{print $1}')
EXIT=$?

# 3. Convert to verdict. AE pixel count threshold per case (default 0; per-case overrides allowed by D-02 deferred).
TOTAL_PX=$(identify -format "%[fx:w*h]" "$A")
AE_PCT=$(awk "BEGIN { printf \"%.4f\", $AE / $TOTAL_PX * 100 }")

# 4. Emit row for CHECKLIST.md (planner translates to the matrix cell)
echo "${CASE} ${MODE} AE=${AE} pixels (${AE_PCT}%) verdict=$([ "$EXIT" -eq 0 ] && echo PASS || echo CHECK)"
```

### Pattern 3: Plan 1 Bootstrap (Phase 7 deferred goldens)

**What:** Three independent steps to materialize the 3 goldens that Phase 7 left QSKIP'd, then re-run Phase 7 ctest to confirm flip from QSKIP to PASS.

**When to use:** First task of Phase 8 Plan 1. Required by D-06 + D-07 before any sweep starts.

**Step A — scene01.eps (PostScript byte parity)**
```bash
# Per scene01_driver.f header BOOTSTRAP NOTE — verbatim procedure (verified by reading source)
cd /tmp
cp $MEFISTO/xvue/qt/tests/golden/scene01_driver.f .
gfortran -I$MEFISTO/incl -c scene01_driver.f
gfortran scene01_driver.o $MEFISTO/xvue/xvuelc.o \
        $MEFISTO/xvue/*.o $MEFISTO/util/*.o \
        -L/usr/X11R6/lib -lX11 -lXt -o scene01_x11
xvfb-run --auto-servernum ./scene01_x11
cp TEMPORAIRE.EPS $MEFISTO/xvue/qt/tests/golden/scene01.eps
```

**KNOWN GOTCHA:** the scene01_driver.f header's link-line glob `xvue/*.o util/*.o` may pull in objects this test doesn't need. If link fails on undefined references from .f95 OMP variants, restrict to a curated minimal set: `xvue/xvuelc.o xvue/xvinit.o xvue/xvouvrir.o xvue/xvfermer.o util/{relevant}.o`. Plan 1 should treat the link line as adjustable — the contract is "produces TEMPORAIRE.EPS", not "this exact ld invocation works".

**Step B — wave_legacy.gif**
```bash
# Per VERIFICATION.md §9.2 + bin/convertepsgif (which is just `convert -rotate 90 zfxy0*.eps -extent 980x550 cyl53zfxy.gif`)

# 1. INITIER + MAILLER + (relevant solver) on testa/wave through X11 backend
#    → produces zfxy*.eps frames in $MEFISTOX/wave/
# 2. cd $MEFISTOX/wave && bash $MEFISTO/bin/convertepsgif
# 3. cp cyl53zfxy.gif $MEFISTO/xvue/qt/tests/golden/wave_legacy.gif
```

**Step C — cavity2d_legacy.gif** — same as Step B with `testa/cavity2d`.

**Verify:** `cd xvue/qt/build && xvfb-run --auto-servernum ctest -R '^xvue_qt_(postscript|export)_tests$'` — all 3 SKIPs flip to PASS.

### Pattern 4: pp/*_qt freshness (Phase 7 Gap-A)

**What:** Refresh `pp/*_qt` against the current `libxvueqt.a` so the A/B sweep cannot get a false-PASS from a stale Qt binary.

**When to use:** Plan 1 first task, before any sweep. Per D-08.

**Example:**
```bash
# Source: bin/cbl_tout_qt (verified to exist + executable in this session)

cd $MEFISTO
bin/cbl_tout_qt
# Cleans xvue/qt/build/ and pp/, then rebuilds. Produces pp/pp{mail,elas,flui,ther,nlse}_qt fresh.

# Sanity check: pp/*_qt mtime > libxvueqt.a mtime
LIB_M=$(stat -c '%Y' xvue/qt/build/libxvueqt.a)
for b in pp/pp*_qt; do
  [ "$(stat -c '%Y' "$b")" -ge "$LIB_M" ] && echo "OK: $b" || echo "STALE: $b — refusing sweep"
done
```

### Anti-Patterns to Avoid

- **Use `compare` without dimension guard.** `compare -metric AE -fuzz 5%` silently returns AE=0 / exit=0 when image dimensions differ. Always `identify` and resample first. (T-08-02 / T-08-05 mitigation.)
- **Use `-fuzz 100%` "to be safe".** Fuzz is a fraction of the color-cube diagonal — at 100% **every** pixel is "within fuzz" of every other pixel and AE returns 0 even on totally different images. Empirically verified: `compare -metric AE -fuzz 100% white.png black.png` returns 0. Keep -fuzz between 1% and 30%; 5% per D-02. (T-08-03 mitigation.)
- **Capture the Qt window via xwd / `import -window <id>`.** Qt's painting model lazily renders to the X11 backing store; under headless or burst-paint sequences the captured X11 window may not match the canvas backing pixmap. Use `MEFISTO_QT_CAPTURE_PATH` (in-process backing-pixmap save). External capture is the X11-side mechanism only.
- **Eyeball-only sign-off when the AE pixel count is plainly available.** D-02 locked tolerance band over eyeball. Maintainer initials sign that the AE count is reasonable AND that the diff PNG looks acceptable, not a substitute for either. (T-08-04 mitigation.)
- **Sign off on a PASS without viewing the diff PNG.** A PASS row with no `evidence/{case}-{mode}-diff.png` reviewed is a gamed sign-off. CHECKLIST.md format mandates the diff path alongside the AE count.
- **Modify xvue/qt/ source during Phase 8.** Phase 8 verifies; it does not fix. Drift triggers either an explicit override (with rationale in CHECKLIST.md) or a separate gap-closure plan. Edits to xvue/qt/ inside Phase 8 violate scope.
- **Run pp/*_qt without first running `bin/cbl_tout_qt`.** Phase 7 Gap-A precedent — stale binaries produce false PASS. (T-08-01 mitigation.)
- **Drive the Qt window with kwin-mcp / synthetic mouse events.** testa cases are batch-driven by their input files. The MEFISTO_XVSOURIS_AUTOEXIT short-circuit handles any blocking-event needs. Avoid Wayland-specific tooling in a Phase 8 contract that must run reproducibly on any Linux box.
- **Use Wayland session for capture.** WAYLAND_DISPLAY is unset on this dev host; even when set, the contract is `xvfb-run --auto-servernum`. Stay X11.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| In-process Qt-side screenshot | A custom ABI entry that calls `QWidget::grab` from Fortran-callable C | Use existing `MEFISTO_QT_CAPTURE_PATH` env var | Already implemented in xvfermer_ in xvue_qt_api.cpp:851. Saves the BACKING pixmap (the source of truth for all drawing). |
| External X11 framebuffer capture | A custom XGetImage + libpng C tool | Use `import -window root` from ImageMagick (already on host) | One CLI call. ImageMagick is in legitimate-allowlist scope (not under xvue/qt/). |
| Image diff with tolerance | Implement AE-with-fuzz in C++ | Use `compare -metric AE -fuzz 5%` | One CLI call, well-documented behavior, dev-host already has it. |
| Synchronisation between solver process and capture tool | Custom IPC, FIFO, signal handler | Use existing `MEFISTO_XVFERMER_READY_FILE` + `MEFISTO_XVFERMER_HOLD_MS` | Already implemented in both backends; designed for exactly this purpose (commit e029b84 added the sentinel for headless capture). |
| Synthetic event injection to terminate `xvsouris_` blocking calls | Custom Qt event injection | Use existing `MEFISTO_XVSOURIS_AUTOEXIT=1` | Already implemented in xvue_qt_api.cpp:1166 + xvuelc.c xvsouris_; documented contract: 1 = synthetic SPACE keypress so `IF (NOTYEV .NE. 2) GOTO 100` exits on first iteration. |
| Per-case batch-mode driver | Hand-rolled Fortran main + custom input | testa cases are already batch-driven via the script files (`pan2d.mesh` ends with `CLOSE;`); `bin/MAILLER $*` forwards to `pp/ppmail $*` for batch input | Existing infrastructure. Confirmed by reading `bin/MAILLER` source — argument is forwarded as-is. |
| Bilingual prompt handling | Set $LANG, mock stdin, etc. | testa cases are language-agnostic at the data level; use `$MEFISTO/td/m/anglais` flag (already present on this host) for English UI | Already wired; English/French only differ in chrome strings, not in batch-input semantics. |

**Key insight:** Almost every piece of headless A/B infrastructure Phase 8 needs is **already present in the codebase**. It was designed in (Phase 5 + Phase 7) for exactly this purpose. The Phase 8 plans should be small wrappers that orchestrate these existing hooks, not new tools.

## Common Pitfalls

### Pitfall 1: pp/*_qt staleness (false PASS)
**What goes wrong:** A/B sweep runs against a `pp/ppmail_qt` linked from an older `libxvueqt.a` than the current source. Drift in the new code is invisible because the binary is stale.
**Why it happens:** `cmake --build xvue/qt/build` only rebuilds `libxvueqt.a`; it does NOT relink `pp/*_qt` (those are produced by `bin/cb*_qt` shell scripts). Phase 7 Gap-A surfaced this exact gap.
**How to avoid:** Plan 1 first task is `bin/cbl_tout_qt` from a clean tree (D-08). Add a sanity check: `for b in pp/pp*_qt; do [ "$(stat -c '%Y' "$b")" -ge "$(stat -c '%Y' xvue/qt/build/libxvueqt.a)" ] && echo OK || echo STALE; done`.
**Warning signs:** AE pixel count of 0 across all cells when you expect AA differences — the binary literally hasn't been rebuilt since the last run that produced those captures.

### Pitfall 2: ImageMagick `compare` silently returns AE=0 on dimension mismatch
**What goes wrong:** Capturing the X11 baseline at 1280x800 and the Qt HiDPI-2x backing pixmap at 1280x960 yields a `compare` exit-0 with AE=0 — making the cell look like PASS when in reality the dimensions differ and the diff was never computed pixel-wise.
**Why it happens:** `compare` (without `-subimage-search`) returns 0 when its inputs cannot be aligned. Empirically verified in this session: `compare -metric AE -fuzz 5% 100x100.png 200x100.png` returns AE=0 / exit=0.
**How to avoid:** Mandatory dimension guard before every `compare` call: `identify -format "%wx%h"` on both inputs, refuse-or-resample on mismatch (Pattern 2 above). Use `convert ... -filter point -resize ${DIMS}!` (point sampling, no AA) so the resample doesn't introduce its own AE drift.
**Warning signs:** Multiple cells reporting AE=0 in cases where you'd expect at least a few AA-edge pixels of drift. Sanity-check by viewing a diff PNG — if it's all-zero on a non-trivial scene, you almost certainly hit this.

### Pitfall 3: -fuzz at 100% silently passes everything
**What goes wrong:** A planner thinking "let's set fuzz really high to be safe" sets `-fuzz 100%` and every cell PASSes regardless of content. Empirically verified: `compare -metric AE -fuzz 100% white.png black.png` returns AE=0 / exit=0.
**Why it happens:** Fuzz is a fraction of the **color-cube diagonal**. At 100% the entire cube is "within fuzz", so every pixel matches every other pixel.
**How to avoid:** Plan validation should assert `-fuzz` is between 1% and 30% (D-02 says 5%). Document in the sweep script that fuzz is a percent of the color-cube diagonal, not a delta-E delta. Cross-check with one canonical "must FAIL" pair (a flat-red vs flat-blue PNG of matching dims) at the start of every Plan 2 wave to prove the gate isn't broken.
**Warning signs:** All cells PASS with non-trivially differing scenes.

### Pitfall 4: HiDPI direct compare against 1x baseline
**What goes wrong:** With `QT_SCALE_FACTOR=2` and Xvfb 1280x800, the Qt logical canvas is 640x400 but the **backing pixmap is 1280x800** (DPR=2). The X11 baseline is 1280x800 too. Naively comparing the 1280x800 Qt-2x capture against the 1280x800 X11 capture pixel-for-pixel **doesn't compare like-for-like**: the Qt-2x scene is rendered at 2x resolution then downsampled by the backing-pixmap save. An honest HiDPI gate is: "the Qt-2x scene downsampled to 1x should match the X11 1x baseline within fuzz."
**Why it happens:** Phase 7 README's Pitfall 7 (logical-vs-physical) documents the same trap for export geometry: PDF uses logical dims, raster formats use backing physical dims.
**How to avoid:** Resample the Qt-2x capture to the 1x baseline dims using `convert -filter point -resize 1280x800!`, then `compare`. Or capture a separate "Qt-2x against Qt-1x" cell (sanity that 2x is internally consistent) plus "Qt-2x downsampled vs X11 1x" (the actual A/B verdict). Both columns can live in CHECKLIST.md.
**Warning signs:** Qt HiDPI cell has wildly higher AE than Qt-1x cell when the underlying scene is identical.

### Pitfall 5: OMP non-determinism produces frame-to-frame jitter
**What goes wrong:** Two consecutive runs of the same case under `OMP_NUM_THREADS=8` produce slightly different captures because OpenMP scheduling reorders some compute and the round-off feeds back into where solver-output coordinates land. The AE diff between two Qt-OMP captures is non-zero even though no graphics-side bug exists.
**Why it happens:** OMP scheduling is not deterministic across runs. The MEFISTO solvers don't pin reduction order under `_OMP`. Graphics drift caused by OMP would be a real bug; **numerical drift** that *renders* as graphics drift is just OMP scheduling.
**How to avoid:** Phase 8's OMP cell compares **Qt-OMP vs X11-OMP** (same backend's input → same expected drift), NOT Qt-OMP vs Qt-1x. Both runs of the same case under the same OMP_NUM_THREADS=8 should match within fuzz. If they don't, that's a clue the OMP path has a real bug. Also run one repeat of the OMP cell to characterize the noise floor — record the AE pixel count of two consecutive Qt-OMP runs as the "OMP noise" baseline; the A/B threshold is "AE between Qt-OMP and X11-OMP < AE between two Qt-OMP repeats × 2" or similar. **Or** simpler: set `OMP_NUM_THREADS=1` for the OMP-path *capture step* (still validates that the OMP-built binary doesn't crash or have main-thread guard violations) and 8 for the *solver math step*. Use `XVUE_QT_ASSERT_MAIN_THREAD()` debug builds to catch off-thread graphics calls explicitly.
**Warning signs:** AE varies wildly run-to-run on the same Qt-OMP cell with no source change.

### Pitfall 6: Gradient color bars (stress / temperature / velocity) trip 5% fuzz
**What goes wrong:** Cases like `nafems_le1` (stress color bar), `heat1d` (temperature ramp), `cavity2d` (velocity arrows) have large gradient regions where Qt's slightly-different colormap interpolation legitimately drifts past 5% fuzz on edge pixels. The cell FAILs even though the visual is acceptable.
**Why it happens:** The Phase 3 colormap was validated to "1 bit per channel" RGB match (TEXT-05) — that's `1/256 ≈ 0.4%` per channel, well within 5% fuzz **per pixel**. But on a gradient where adjacent pixels ramp through colors, off-by-one indexing produces a 1-pixel shifted gradient where every gradient-edge pixel registers as "different" because each pixel jumped to a 1-bit-different color.
**How to avoid:** D-02's "Deferred Ideas" section explicitly allows per-case fuzz overrides for gradient-heavy cases. CHECKLIST.md per-cell columns can carry a `fuzz_override_pct` annotation when set. Recommended override values, derived from the AE pixel count formula: nafems_le1 stress bar: `-fuzz 8%`; heat1d temperature: `-fuzz 8%`; cavity2d velocity arrows: `-fuzz 10%`. Document the override in the cell so a future agent knows it was an intentional relaxation, not a missed bug.
**Warning signs:** A cell FAILs with a high AE concentrated in the colorbar region of the diff PNG.

### Pitfall 7: Anti-aliased text + font hinting drift
**What goes wrong:** `pan2d` and `hexahedron` have node-number labels rendered via `xvtexte_`. Qt's font hinting (FreeType) and X11's font hinting differ at sub-pixel level, so every glyph edge has a few drifted pixels. Even with -fuzz 5%, an AE pixel count of a few hundred is normal and not a defect.
**Why it happens:** Font rendering is intrinsically AA-different across backends. Phase 1 SHELL-06 only validated "windows are visibly the same size" not "every glyph pixel matches".
**How to avoid:** Set a per-case AE pixel-count threshold in CHECKLIST.md alongside the raw count. For `pan2d` (10–20 node labels at typical scale), an AE budget of `~200 pixels per label × 15 labels ≈ 3000 pixels` is reasonable. Document the budget as an explicit cell annotation. The maintainer's initials sign that "the AE count is within budget AND the diff PNG shows only glyph-edge fringing, no missing/clipped/overlapped labels."
**Warning signs:** A cell FAILs with a high AE concentrated around text glyphs only — not a defect, just AA drift.

### Pitfall 8: Multi-window / re-open lifecycle interactions
**What goes wrong:** Some test cases run multiple plot phases (mesher draws → solver redraws → post-processor draws). The capture hook in `xvfermer_` only fires once per process. If the case calls `xvfermer_` mid-run, the captured PNG might be a transient frame, not the final state.
**Why it happens:** `xvfermer_` is "close the window" — its capture hook is "save the canvas state right before destruction." Some testa cases close-and-reopen; only the LAST close gets captured.
**How to avoid:** For cases known to close-and-reopen (TBD per case — Plan 1 should probe each case to identify), either (a) explicit document which `xvfermer_` call is the "publishable" frame, or (b) write captures to numbered files (`pan2d-qt-1x-{1,2,3}.png`) and pick the last as the canonical evidence. The MEFISTO_QT_CAPTURE_PATH currently overwrites — Plan 1 should test whether that's OK or if the capture should be made idempotent / suffixed.
**Warning signs:** The captured PNG shows a partial scene (mesh only, no solver overlay) when the test was supposed to display the full result.

## Runtime State Inventory

This phase has **no rename / refactor / migration character**. It is a pure verification phase. Section is omitted intentionally — verifying that no state-bearing strings need editing is, itself, the answer.

(For completeness: Phase 8 writes new files under `.planning/phases/08-*/evidence/` and `.planning/phases/08-*/08-CHECKLIST.md`, and bootstraps three goldens under `xvue/qt/tests/golden/`. None of these are renames or migrations. No external service, OS-registered task, secret, or build artifact carries a Phase-8-specific string that needs maintenance after Phase 8 closes.)

## Code Examples

### Verified pattern: in-process Qt capture via `MEFISTO_QT_CAPTURE_PATH`

```cpp
// Source: xvue/qt/src/xvue_qt_api.cpp lines 824-895 (verified by Read tool in this session)
// This is ALREADY IMPLEMENTED — Phase 8 plans consume it, do not modify it.

void proc(xvfermer)(void) {
    XvueApp::ensure();
    XVUE_QT_ASSERT_MAIN_THREAD();
    {
        auto& win = XvueApp::window_slot();
        if (win && win->canvas()) {
            win->canvas()->update();
        }
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);

        const char* qt_capture_path = std::getenv("MEFISTO_QT_CAPTURE_PATH");
        if (qt_capture_path && qt_capture_path[0] != '\0') {
            // In-process screenshot of the XvueState backing pixmap —
            // bypasses the widget paint pipeline; works under
            // QT_QPA_PLATFORM=offscreen with NO X11 dependency.
            XvueState* st = (win ? win->state() : nullptr);
            if (st && st->backing_) {
                if (st->painter_ && st->painter_->isActive()) {
                    st->painter_->end();
                }
                if (!st->backing_->save(QString::fromLatin1(qt_capture_path), "PNG")) {
                    std::fprintf(stderr, "xvue-qt: failed to save backing to %s\n", qt_capture_path);
                }
            }
            // ... fallback to widget grab if no backing
        }
        // ... READY_FILE + HOLD_MS handling (mirrors X11 backend)
    }
    XvueApp::window_slot().reset();
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}
```

### Verified pattern: X11-side ready-file sentinel

```c
// Source: xvue/xvuelc.c xvfermer_ (verified by Read tool in this session)
// This is ALREADY IMPLEMENTED in the X11 backend.

void proc(xvfermer)() {
    // Headless-test capture hook.
    char *ready_path = getenv("MEFISTO_XVFERMER_READY_FILE");
    char *hold_env   = getenv("MEFISTO_XVFERMER_HOLD_MS");
    int hold_ms = 0;
    if (hold_env != NULL && hold_env[0] != '\0') {
        int v = atoi(hold_env);
        if (v >= 0 && v <= 60000) hold_ms = v;
    }
    XFlush(display_mef);
    if (ready_path && ready_path[0]) {
        FILE *f = fopen(ready_path, "w");
        if (f) { fprintf(f, "ready\n"); fclose(f); }
    }
    if (hold_ms > 0) usleep((useconds_t)hold_ms * 1000);
    // ... destroy the X11 display
}
```

### New code (Phase 8 deliverable): tolerance-band gate wrapper

```bash
#!/bin/bash
# bin/ab_compare_pair.sh — Phase 8 Plan 2+ deliverable
# Source: empirical compare 7.1.2 behavior verified in this session
set -euo pipefail

A=$1; B=$2; DIFF=$3; FUZZ=${4:-5}      # default 5% per D-02

# Dimension guard — refuse-or-resample. Without this, compare returns AE=0 silently on size mismatch.
DA=$(identify -format "%wx%h" "$A")
DB=$(identify -format "%wx%h" "$B")
if [ "$DA" != "$DB" ]; then
    BR="${B%.png}-resampled-to-${DA}.png"
    convert "$B" -filter point -resize "${DA}!" "$BR"
    echo "RESAMPLED: $B ($DB) → $BR ($DA)" >&2
    B="$BR"
fi

# AE-fuzz gate. compare returns AE pixel count + exit 0 (similar) / 1 (different) / 2 (error).
# Use 2>&1 because compare writes the metric to stderr.
AE_LINE=$(compare -metric AE -fuzz "${FUZZ}%" "$A" "$B" "$DIFF" 2>&1) || EXIT=$? && EXIT=${EXIT:-0}
AE=$(echo "$AE_LINE" | awk '{print $1}')

TOTAL_PX=$(identify -format "%[fx:w*h]" "$A")
AE_PCT=$(awk "BEGIN { printf \"%.4f\", $AE / $TOTAL_PX * 100 }")
VERDICT=$([ "${EXIT:-0}" -eq 0 ] && echo PASS || echo CHECK)

echo "ae=$AE total=$TOTAL_PX pct=$AE_PCT% verdict=$VERDICT diff=$DIFF"
```

### New code (Phase 8 deliverable): CHECKLIST.md row schema

```markdown
| Case | X11 (baseline) | Qt 1x | Qt HiDPI 2x | Qt _OMP | Initials | Notes |
|------|---------------|-------|------------|---------|----------|-------|
| pan2d | PASS — evidence/pan2d-x11.png (10240×768 px) | PASS — AE=183 (0.0023%) — evidence/pan2d-qt-1x.png + diff-1x.png | PASS — AE=412 (0.0052%) downsampled-to-baseline-dims — evidence/pan2d-qt-2x.png + diff-2x.png | PASS — AE=183 (0.0023%, identical to 1x; main-thread guard held) — evidence/pan2d-qt-omp.png + diff-omp.png | (initials) | Font-AA fringing on glyph edges (Pitfall 7); fuzz=5% default. |
| nafems_le1 | PASS | PASS — fuzz override 8% (Pitfall 6 stress bar) — AE=2103 (0.026%) | PASS — fuzz 8% — AE=2247 (0.028%) | PASS | (initials) | Stress-bar gradient drift; per-case fuzz override per D-02 deferred. |
| ... | | | | | | |

**v1 ship gate sign-off:** every cell ∈ {PASS, N-A}. Maintainer signature: ____________________ Date: __________.
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `xwd | convert` for X11 capture | `import -window root` (ImageMagick CLI) | x11-apps removed from Debian forky/sid; xwd no longer present on the dev host. Verified in this session. | Switch the X11 capture wrapper to `import` (already present, in legitimate-allowlist scope under bin/). |
| `QPrinter::PdfFormat` | `QPdfWriter` | Qt 6 (vs Qt 5) | Phase 7 already uses `QPdfWriter` — Phase 8 inherits unchanged. (Affects nothing in Phase 8 directly; recorded for context.) |
| Auto-built `_OMP` binaries (older MEFISTO docs) | Same `pp/ppelas` binary linked against libgomp; `_OMP` launcher only sets `OMP_NUM_THREADS=8` | Pre-existing in MEFISTO; not new to xvue-qt | Phase 8 does NOT need to build separate OMP binaries. Verified via `nm pp/ppelas | grep GOMP` + `bin/ELASTICER_OMP` source. |
| ImageMagick `convert` for animated GIF in Qt path | Qt-side `QImageWriter` (probed) → `ffmpeg` fallback | Phase 7 PROBE.md (gif_write_supported=0 → ffmpeg path) | Phase 8 Plan 1 reference baselines use the legacy `bin/convertepsgif` X11+convert pipeline (D-06) — that's intentional; the LEGACY baselines are made with the legacy tool. |

**Deprecated/outdated:**
- xwd (not on this host); `xdotool` (not on this host) — neither is needed for Phase 8 given the existing `MEFISTO_*_AUTOEXIT` / `MEFISTO_*_CAPTURE_PATH` hooks.
- ImageMagick policy.xml restrictions — verified compare/import/convert all functional; PDF/PS reading restrictions in modern IM policy don't affect Phase 8 (we read/write PNG and GIF, not PS/PDF, except via `convertepsgif` which is its own legacy binary).

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Each testa case can run end-to-end under `MEFISTO_XVSOURIS_AUTOEXIT=1` from INITIER → solver in batch mode without human input | Architecture Patterns / Pattern 1 | If a case has an interactive prompt that AUTOEXIT can't satisfy (e.g., a `read` not gated by a notypeevent loop), the case will hang in Plan 2 sweep. Mitigation: Plan 1 should include a per-case dry-run smoke test that confirms the case completes within a 60s timeout under AUTOEXIT before Plan 2 launches. |
| A2 | The 5 canonical cases produce a single "publishable" frame per (case, mode) cell — i.e., one xvfermer_ call at process end captures the canonical scene | Pitfall 8 | If a case calls xvfermer_ mid-run, the capture overwrites and the final PNG isn't the "ship" frame. Mitigation: Plan 1 probes each case to count xvfermer_ calls; for >1, switch to numbered captures or document which is canonical. |
| A3 | The `_OMP` launcher does NOT require a separately-built `pp/ppelas_omp` binary — only the existing `pp/ppelas` with `OMP_NUM_THREADS=8` set | OMP Sweep + VALID-03 mapping | If a separate OMP binary IS expected (older MEFISTO convention) the OMP cells run the wrong binary. Mitigation: Plan 1 inspects `nm pp/ppelas` / `nm pp/ppelas_qt` for GOMP symbols and asserts both are linked against libgomp before sweep. (Already verified in this session — both binaries DO link libgomp.) |
| A4 | The `bin/MAILLER` argument forwards correctly to `pp/ppmail $*` for Qt binary too — i.e., `pp/ppmail_qt /path/to/script.mesh` reads the script in batch mode the same way the X11 binary does | Architecture Patterns / Pattern 1 | Qt batch driver behavior may diverge from X11 if Phase 6.0+ rewired stdin handling. Mitigation: Plan 1 includes a one-shot smoke `pp/ppmail_qt testa/pan2d/pan2d.mesh` under offscreen + AUTOEXIT to verify batch parity before Plan 2 sweep. |
| A5 | Qt HiDPI Plan 2 captures on `MEFISTO_QT_CAPTURE_PATH` save the BACKING pixmap which is dpr-scaled — so a QT_SCALE_FACTOR=2 capture is 2x the X11 baseline dims and must be downsampled before AE compare | Pitfall 4 | If captures are actually saved at logical dims (not physical), the resample step is unnecessary and the gate's tighter than expected. Mitigation: Plan 1 captures one Qt-1x and one Qt-2x of the same scene and confirms dimension ratio is 2:1; documents in 08-CHECKLIST.md whether downsample-before-compare or compare-as-is. (Read of xvue_qt_api.cpp:851-870 strongly suggests backing dims = physical.) |
| A6 | Phase 7 Plan 1 bootstrap of scene01.eps + wave_legacy.gif + cavity2d_legacy.gif works on this exact Debian forky/sid host with libgfortran5 held at 15.2.0-9 | Plan 1 procedure | If the held libgfortran5 + Debian sid + gfortran-14.3 combination has new latent UB sites that didn't exist when scene01_driver.f was written, the link or run could fail. Mitigation: STATE.md already pins libgfortran5; Phase 8 should re-verify the pin holds at Plan 1 entry; if scene01_x11 link fails, fall back to a curated minimal `xvue/xvuelc.o + util/{listed}.o` link line. |
| A7 | The `convertepsgif` shell script (a one-line `convert -rotate 90 zfxy0*.eps -extent 980x550 cyl53zfxy.gif`) actually produces a valid GIF for testa/wave and testa/cavity2d on this host's ImageMagick 7.1.2-18 | Plan 1 / Bootstrap Step B+C | If IM 7.1.2 GIF encoding has changed format defaults relative to when the legacy pipeline was authored (e.g., color quantization, palette size), the legacy GIF may differ slightly from a "true historical" baseline. This is acceptable for Phase 8 since the baseline IS what bin/convertepsgif on THIS host produces — that becomes the reference that future Qt runs must match. Mitigation: document in Plan 1 SUMMARY that the baseline is "what convertepsgif on dev-host-ImageMagick-7.1.2 emits today", not a historical archive. |

**Note on A6/A7:** Both relate to the Phase 7 deferred-golden bootstrap. They are LOW risk because the goldens are committed to the repo after Plan 1 and any subsequent compare is against the committed baseline; downstream agents don't need to re-derive them.

## Open Questions

1. **Per-case fuzz override values for the 3 known gradient/text cases (nafems_le1 stress bar, heat1d temp ramp, cavity2d velocity arrows, pan2d node labels, hexahedron labels)**
   - What we know: D-02 deferred ideas allows per-case fuzz overrides; Phase 3 colormap match is 1-bit-per-channel.
   - What's unclear: the exact fuzz value that yields a real PASS without hiding genuine bugs.
   - Recommendation: Plan 2 first-batch run captures the 5 cases at default 5% fuzz, records the AE pixel counts, and proposes per-case overrides only where the diff PNG visibly shows AA-only or gradient-only drift. Approval of overrides happens in 08-CHECKLIST.md cell annotations during Plan 2+ before sign-off.

2. **Should the OMP cell compare against X11-OMP or X11-1x?**
   - What we know: D-05 says "sweep all 5 _OMP variants where they exist". X11 binaries also link libgomp.
   - What's unclear: whether the X11 baseline column should run `OMP_NUM_THREADS=1` (deterministic baseline) or `OMP_NUM_THREADS=8` (matching the Qt OMP cell condition).
   - Recommendation: Capture two X11 columns — `X11-1x` (OMP_NUM_THREADS=1, deterministic) and `X11-OMP` (OMP_NUM_THREADS=8). Compare Qt-OMP against X11-OMP (same OMP context) for the OMP cell. Compare Qt-1x against X11-1x for the baseline cell. CHECKLIST.md gets one extra column but the comparison is honest.

3. **What's the "publishable frame" for cases that call `xvfermer_` more than once?**
   - What we know: Pitfall 8 — current `MEFISTO_QT_CAPTURE_PATH` overwrites; only the LAST `xvfermer_` saves.
   - What's unclear: which testa cases close-and-reopen (mesher → solver → post sequence likely does for cavity2d/wave).
   - Recommendation: Plan 1 (or first task of Plan 2) probes each case under AUTOEXIT and records the count of `xvfermer_` calls (instrument xvue_qt_api.cpp's xvfermer_ debug log, or count "TRYing to SAVE DATA" patterns in the log). For cases with >1, switch to suffixed paths `MEFISTO_QT_CAPTURE_PATH=evidence/${CASE}-qt-1x-${IDX}.png` (would need a one-line patch to xvfermer_ — but that's a Phase-8-out-of-scope source edit, so the alternative is "run case in stages"). Cleanest answer: Plan 1 verifies that the "last frame" semantics are acceptable as the canonical capture; if not, document which testa-case stages need separate capture PNGs.

4. **HiDPI Pitfall 4 — confirm dimension ratio empirically before sweep**
   - What we know: empirically confirmed in this session that Qt logical screen halves under QT_SCALE_FACTOR=2 (`QSize(640,400) dpr=2`).
   - What's unclear: whether `MEFISTO_QT_CAPTURE_PATH` saves the BACKING pixmap (physical, dpr-scaled) or the canvas widget (logical, pre-dpr). Reading xvue_qt_api.cpp:851 says it saves `st->backing_->save(...)` — backing_ is the dpr-scaled pixmap (per Phase 7 README HiDPI table — PNG export uses backing physical pixels).
   - Recommendation: Plan 1 captures a known-canvas-size scene under both 1x and 2x and asserts the captured PNGs are 1x and 2x dimensions respectively. If they aren't, the resample direction in Pattern 2 needs reversing.

5. **Should Plan 2 capture diff PNGs even on PASS cells, or only on FAIL?**
   - What we know: D-01 says committed evidence under evidence/. Storage is ~5 MB total target.
   - What's unclear: whether 20 PASS diff PNGs (mostly black images) plus 20 source PNGs fit the ~5 MB budget.
   - Recommendation: Always commit diff PNGs (they compress very well — black-mostly images are ~5–20 KB each). Total stays under 5 MB. Diff PNG path in CHECKLIST.md cell is a hard contract per D-10 (path-to-evidence) so storing diffs is operationally simpler than skipping them.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| ImageMagick (`compare`, `import`, `identify`, `convert`) | D-03 tolerance gate; X11 capture (`import -window root`); legacy GIF pipeline (`convertepsgif`) | ✓ | 7.1.2-18 Q16 | — (mandatory) |
| Xvfb / `xvfb-run` | D-04 HiDPI sweep host; X11 baseline capture | ✓ | 21.1.21-1 | — (mandatory) |
| Qt 6 (libQt6{Core,Gui,Widgets}) | pp/*_qt runtime | ✓ | 6.10.2+dfsg-7 | — (mandatory) |
| gfortran | bin/cbl_tout_qt rebuild + Plan 1 scene01_x11 link | ✓ | 14.3.0 (Debian 14.3.0-14) | — (mandatory) |
| libgfortran5 (held) | Stable Fortran runtime — STATE.md 03-04 reopen | ✓ | 15.2.0-9 (held via apt-mark) | — (mandatory; see A6) |
| ffmpeg | XVUE_ANIM=1 GIF fallback (Phase 7 D-11); Plan 1 reference verify only — Plan 1 baseline uses convertepsgif | ✓ | 8.1-3+b1 | None — but Plan 1 wave/cavity2d goldens use the **legacy** `convertepsgif` ImageMagick path, NOT ffmpeg, by D-06's explicit reuse of the legacy pipeline. |
| `import` (ImageMagick) for X11 root capture | Pattern 1 X11 capture | ✓ | (bundled with ImageMagick 7.1.2-18) | scrot 2.0.0-1 (also installed) |
| `scrot` | Alternate root capture | ✓ | 2.0.0-1 | Use `import -window root` |
| `xwd` | Legacy X11 window dump | ✗ | — | `import -window root` (preferred) — `xwd` is NOT installed on this Debian forky/sid host (x11-apps removed) |
| `xdotool` | Synthetic input | ✗ | — | Not needed — `MEFISTO_XVSOURIS_AUTOEXIT=1` short-circuits all blocking event reads. |
| KDE `kwin-mcp` | Wayland session driver (Phase 7 UAT lineage) | ✓ | (~/.local/share/kwin-mcp-venv) | Use Xvfb instead — Phase 8 is X11-headless-only, do not introduce Wayland dependency in the AE pipeline. |
| WAYLAND_DISPLAY | Wayland session | ✗ (unset) | — | N/A — Phase 8 must run under Xvfb. |

**Missing dependencies with no fallback:** None.

**Missing dependencies with fallback:**
- `xwd` → `import -window root` (already on host).
- `xdotool` → `MEFISTO_XVSOURIS_AUTOEXIT=1` (already in xvue_qt_api.cpp + xvuelc.c).

## Validation Architecture

> Phase 8 IS validation, so this section establishes how Phase 8's own deliverables are validated by `gsd-verify-work` rather than how it validates other things. Nyquist enforcement is enabled by `.planning/config.json` workflow.nyquist_validation = true.

### Test Framework

| Property | Value |
|----------|-------|
| Framework | Existing Qt ctest (`xvue_qt_postscript_tests`, `xvue_qt_export_tests`) + new shell harness (`bin/ab_sweep_phase8.sh`, `bin/ab_compare_pair.sh`) |
| Config file | `xvue/qt/CMakeLists.txt` (existing), and per-plan SUMMARY's verification section |
| Quick run command | `bash bin/ab_compare_pair.sh evidence/pan2d-x11.png evidence/pan2d-qt-1x.png evidence/pan2d-1x-diff.png` |
| Full suite command | `xvfb-run --auto-servernum bash bin/ab_sweep_phase8.sh && cd xvue/qt/build && xvfb-run --auto-servernum ctest -R '^xvue_qt_(postscript\|export)_tests$'` |

### Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| VALID-01 | All 5 canonical cases captured + AE-diff'd in {Qt 1x} cell | shell-harness + image-diff | `bash bin/ab_sweep_phase8.sh --mode qt-1x --cases pan2d,nafems_le1,cavity2d,heat1d,nlsecu` | ❌ Wave 0 |
| VALID-02 | 5 cases on X11 baseline produce captures (column 1) | shell-harness | `bash bin/ab_sweep_phase8.sh --mode x11 --cases <5>` | ❌ Wave 0 |
| VALID-03 | All 5 _OMP variants run + main-thread guard holds | shell-harness w/ OMP_NUM_THREADS=8 | `bash bin/ab_sweep_phase8.sh --mode qt-omp --cases <5>` | ❌ Wave 0 |
| VALID-04 | All 5 cases at QT_SCALE_FACTOR=2 + downsampled-compare | shell-harness w/ QT_SCALE_FACTOR=2 | `bash bin/ab_sweep_phase8.sh --mode qt-2x --cases <5>` | ❌ Wave 0 |
| VALID-05 | Color-bar drift on nafems_le1, heat1d, cavity2d under 5% (or per-case override) fuzz | image-diff w/ override | embedded in ab_sweep_phase8.sh per-case fuzz lookup | ❌ Wave 0 |
| VALID-06 | Font-metric drift on pan2d, hexahedron under 5% fuzz + AE budget per Pitfall 7 | image-diff | same harness; hexahedron is a SPOT-CHECK row outside the 5-case grid | ❌ Wave 0 |
| VALID-07 | 08-CHECKLIST.md per-cell signed-off matrix | manual sign-off + grep | `grep -c 'PASS\|N-A' .planning/phases/08-*/08-CHECKLIST.md` ≥ 20 | ❌ Wave 0 |
| BUILD-10 | 5-case baseline immutable | grep | `diff <(awk '/^### / {print $2}' .planning/validation/BASELINE.md) <(echo -e '1.\n2.\n3.\n4.\n5.')` | ✅ (BASELINE.md exists) |
| Phase-7-close (D-07) | Phase 7 SKIPs flip to PASS | ctest | `cd xvue/qt/build && xvfb-run --auto-servernum ctest -R '^xvue_qt_(postscript\|export)_tests$' --output-on-failure | grep -c '0 SKIP'` ≥ 1 | ✅ (tests exist; goldens TBD) |

### Sampling Rate

- **Per task commit:** `bash bin/ab_compare_pair.sh <one-pair>` (~5 seconds per pair)
- **Per wave merge:** Full Plan 2 sweep over assigned mode (~5–15 minutes for 5 cases × 1 mode)
- **Phase gate:** Full sweep over all 4 modes + 6 spot-check cells + ctest re-run (~30–45 minutes total) before `/gsd-verify-work` + 08-CHECKLIST.md signed off

### Wave 0 Gaps

- [ ] `bin/ab_sweep_phase8.sh` — top-level harness (NEW — Phase 8 Plan 2 deliverable)
- [ ] `bin/ab_capture_x11.sh` — Xvfb-display-aware X11 capture wrapper (NEW)
- [ ] `bin/ab_compare_pair.sh` — tolerance-band gate (NEW; example sketched in Code Examples)
- [ ] `xvue/qt/tests/golden/scene01.eps` — Phase 7 deferred (Plan 1)
- [ ] `xvue/qt/tests/golden/wave_legacy.gif` — Phase 7 deferred (Plan 1)
- [ ] `xvue/qt/tests/golden/cavity2d_legacy.gif` — Phase 7 deferred (Plan 1)
- [ ] `.planning/phases/08-*/evidence/` directory — created on first sweep

## Project Constraints (from CLAUDE.md)

- **Compilation must never break.** After every change, the affected module must still compile via its `cb*` script. Phase 8 Plan 1 runs `bin/cbl_tout_qt` (full Qt rebuild) and the X11 build is touched only via Plan 1 Step A scene01_driver compile (which does NOT modify any tracked source). Plans must check the build is green before recording a sweep cell PASS.
- **xvuelc.c byte-identical until Phase 9.** Phase 8 does not modify xvue/xvuelc.c. The X11 capture sentinel `MEFISTO_XVFERMER_READY_FILE` already lives in xvuelc.c (added pre-Phase-7); reuse only.
- **Large/visual tests are user-run.** Phase 8 IS the canonical "large/visual A/B test" the project guards. The maintainer's initials in 08-CHECKLIST.md ARE the sign-off CLAUDE.md mandates. Plans MUST NOT auto-PASS cells without maintainer review of the diff PNG.
- **Tests in testa/ must continue to pass after every change.** Phase 8's "change" is verification-only — there should be no behavioral change to make tests fail. If a sweep cell FAILs unexpectedly, the response is to investigate the existing code (open a separate gap-closure plan) — NOT to "fix while validating" inside Phase 8 (scope discipline).
- **Use the smallest relevant test case when checking a fix.** Phase 8 IS the smallest reasonable test set for declaring v1 ship — 5 canonical cases × 4 modes = 20 cells, plus 1 hexahedron font spot-check.
- **`99;` save-exit; never Ctrl-C.** testa cases are batch-driven by their input scripts which terminate with `CLOSE;` or by exhausting input. `MEFISTO_XVSOURIS_AUTOEXIT=1` provides a synthetic event that exits any blocking xvsouris_ loop. There should be no Ctrl-C anywhere in Phase 8 plans.
- **Ask before installing.** No new system packages required (verified — all tooling is present on this host).
- **Git discipline:** commit at each plan boundary; commits must describe what + why; no force-push; no hook-bypass. Standard Phase 8 commits will be `feat(08-01)` (Plan 1 bootstrap), `feat(08-02..N)` (sweep batches), `docs(08)` (CHECKLIST.md sign-off).

## Phase Boundary Discipline (from prior phases)

How prior phases handled drift discovered during validation:

- **Phase 02.1** (INSERTED) — Phase 2/3 A/B catch-up gate found xvface color drift on xvtest2/xvtest4. Response: an INSERTED phase (02.1) was created with a single 1-plan scope to fix `xvfacetraits_` ncf/nca color drop. Phase 2 itself was NOT amended; Phase 3 didn't fold the fix in.
- **Phase 03.1** (INSERTED) — Phase 3 A/B gate found missing `bin/cb*` scripts to build xvtest1..4 on both backends. Response: an INSERTED phase (03.1) was created with 3 plans. Phase 3 was not extended.
- **Phase 7 Gap-A** — UAT discovered stale `pp/ppmail_qt`. Response: documented as a Gap-A in 07-UAT.md, ROUTED to Phase 8 Plan 1 (D-08), AND queued as a permanent CMake guard for Phase 9 (D-09 deferred). Phase 7 itself did not pull the fix in.
- **Phase 7 deferred goldens** — autonomous executor couldn't bootstrap goldens. Response: 3 QSKIPs in ctest + DEFERRED rows in 07-VALIDATION-LOG.md, ROUTED to Phase 8 Plan 1 (D-06).

**Pattern:** Drift found during validation → create gap-closure plan, OR route to next phase, OR document as deferred. NEVER fix-while-validating inside the validation phase itself. Phase 8 plans MUST follow this pattern.

**Concrete rule for Phase 8:**
- A FAILing cell that's a real defect → escalate to a "Phase 08.1 INSERTED gap-closure" plan or to Phase 9 cleanup. Document as FAIL in CHECKLIST.md with rationale + escalation path; do NOT edit xvue/qt/ inside Phase 8 to fix it.
- A FAILing cell that's a known acceptable difference (e.g., AA fringing, gradient drift within Pitfall 6/7 budget) → mark as PASS with `fuzz_override_pct` annotation OR explicit override note in CHECKLIST.md cell.
- A test infrastructure issue (capture mechanism doesn't work on this host) → fix the harness inside Phase 8 (the harness IS Phase 8 deliverable; it's not "the code under test").

## OMP Sweep — How _OMP Actually Works

Empirically verified in this session via `nm`, `ldd`, and `bin/ELASTICER_OMP` source read:

- There is **no separate `pp/ppelas_omp` binary**. The `bin/ELASTICER_OMP` launcher invokes the same `pp/ppelas` binary that `bin/ELASTICER` invokes; the only difference is that `_OMP` exports `OMP_NUM_THREADS=8` first.
- Both `pp/ppelas` and `pp/ppelas_qt` are linked against `libgomp.so.1` already (verified via `ldd`). Phase 8 OMP cells run the same Qt binary with `OMP_NUM_THREADS=8` set.
- Existing `_OMP` launchers on this host: `bin/ELASTICER_OMP`, `bin/MAILLER_OMP` (and alias `bin/MESHER_OMP`), `bin/THERMICER_OMP`, `bin/NLSER_OMP`, `bin/ONDER_OMP` (wave). **There is NO `bin/FLUIDER_OMP`** — the cavity2d-OMP cell must therefore be marked `N-A` per D-05.
- The "main-thread invariant" the OMP sweep tests is `XVUE_QT_ASSERT_MAIN_THREAD()` at every public Qt-side ABI entry point (Phase 1 SHELL-07, all 57 stubs). Under `QT_DEBUG`, off-main-thread graphics calls trip Q_ASSERT and abort the process. Phase 8 OMP cells should run the **debug build** (Q_ASSERT enabled) of the Qt binary if available — or a release build with VERBOSE flag — to surface assertion failures clearly.
- Phase 8 OMP cell mapping (D-05):
  - mesher (mail) — `bin/MAILLER_OMP` exists → run pan2d-omp cell
  - elas — `bin/ELASTICER_OMP` exists → run nafems_le1-omp cell
  - flui — `bin/FLUIDER_OMP` **does NOT exist** → cavity2d-omp cell = `N-A`
  - ther — `bin/THERMICER_OMP` exists → run heat1d-omp cell
  - nlse — `bin/NLSER_OMP` exists → run nlsecu-omp cell

## Sources

### Primary (HIGH confidence)

- `xvue/qt/src/xvue_qt_api.cpp:824-925` (xvfermer_ body — verified by Read tool) — MEFISTO_QT_CAPTURE_PATH + READY_FILE + HOLD_MS contract
- `xvue/qt/src/xvue_qt_api.cpp:1166-1220` (xvsouris_ body — verified by Grep) — MEFISTO_XVSOURIS_AUTOEXIT contract
- `xvue/xvuelc.c xvfermer_` (verified by Grep) — X11 backend READY_FILE + HOLD_MS sentinel
- `bin/ELASTICER_OMP` (verified by Read tool) — OMP launcher invokes same `pp/ppelas` binary, sets OMP_NUM_THREADS=8
- `bin/MAILLER` (verified by Read tool) — argument forwarded as `pp/ppmail $*` for batch mode
- `bin/cbl_tout_qt` + `bin/cb{mail,elas,flui,ther,nlse}_qt` (verified by ls) — Qt build chain refresh path for Pitfall 1 + D-08
- `xvue/qt/tests/golden/scene01_driver.f` (verified by Read tool) — bootstrap procedure header for Plan 1 Step A
- `bin/convertepsgif` (verified — content `convert -rotate 90 zfxy0*.eps -extent 980x550 cyl53zfxy.gif`) — Plan 1 Step B/C legacy GIF pipeline
- `.planning/phases/07-image-gif-and-postscript-export/VERIFICATION.md §9` (3 deferred items + bootstrap procedures) + `07-UAT.md §"Gap-A"` (pp/*_qt staleness) — Phase 7 hand-off
- `.planning/REQUIREMENTS.md §VALID-01..07 + §BUILD-10` — phase requirements
- `.planning/validation/BASELINE.md` — 5 canonical cases lock-in
- `.planning/phases/08-ab-validation-on-testa-subset/08-CONTEXT.md` — D-01..D-10 locked decisions
- `xvue/qt/README.md §"HiDPI export math"` — Pitfall 7 logical-vs-physical pixel convention; basis for Pitfall 4
- `STATE.md "Phase 03-04 reopen close 2026-04-14"` — libgfortran5 pinning at 15.2.0-9 (A6)

### Secondary (MEDIUM confidence — verified empirically in this session against the live host)

- `compare -metric AE -fuzz 5%` returns AE=0/exit=0 on dimension mismatch — verified: `compare -metric AE -fuzz 5% 100x100.png 200x100.png` → "0 (0)" exit 0
- `compare -metric AE -fuzz 100%` returns AE=0 on totally different images — verified: `compare -metric AE -fuzz 100% white.png black.png` → "0 (0)" exit 0
- `compare -metric AE -fuzz 5%` does NOT fully suppress 1px line geometry shift — verified: 162 px diff on a 50→51 px line
- `QT_SCALE_FACTOR=2` under Xvfb 1280x800 yields `QSize(640,400) dpr=2` for primaryScreen — verified: ran custom qtest binary
- ImageMagick `import -window root` succeeds against Xvfb display — verified
- `xwd` and `xdotool` are NOT installed on this host — verified via `command -v` + dpkg
- `pp/ppelas` and `pp/ppelas_qt` both link `libgomp` — verified via `ldd` + `nm | grep GOMP`
- `pkg-config --modversion Qt6Core` = 6.10.2; ImageMagick = 7.1.2-18 Q16; ffmpeg = 8.1-3+b1; gfortran = 14.3.0; libgfortran5 = 15.2.0-9 (held)

### Tertiary (LOW confidence — flagged for validation)

- ImageMagick `compare` documentation (https://imagemagick.org/script/compare.php) — official docs are sparse on AE-fuzz interaction edge cases; multiple discourse threads describe quirks. Empirical verification preferred.
- Qt 6 HiDPI documentation (https://doc.qt.io/qt-6/highdpi.html) — confirms QPainter / QImage handle DPR correctly; doesn't explicitly cover the `MEFISTO_QT_CAPTURE_PATH` capture-of-backing case (which the Phase 7 README does cover, primary above).

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — every tool listed was verified present on the dev host this session, with version recorded.
- Architecture: HIGH — capture mechanism patterns read directly from existing Qt/X11 source code (xvue_qt_api.cpp + xvuelc.c).
- Pitfalls: HIGH — Pitfalls 1-7 are either verified empirically in this session (1, 2, 3, 4) or have direct prior-phase precedent (1 from Phase 7 Gap-A; 7 from Phase 1 SHELL-06; 6 from Phase 3 TEXT-05). Pitfall 8 is a known unknown — flagged as A2 / Open Question 3.

**Research date:** 2026-05-05
**Valid until:** 2026-06-04 (30 days for stable). Re-verify environment if Debian forky/sid pulls a new ImageMagick / Qt 6 / libgfortran point release that the held-pin doesn't cover.

---

*Phase: 08-ab-validation-on-testa-subset*
*Researched: 2026-05-05*
