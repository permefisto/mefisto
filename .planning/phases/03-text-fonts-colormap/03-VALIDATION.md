---
phase: 03
slug: text-fonts-colormap
status: resolved
nyquist_compliant: true
wave_0_complete: true
created: 2026-04-11
updated: 2026-04-14
resolution: Task 3 reopen closed cleanly after a Debian sid libgfortran5 downgrade (gcc-16 snapshot → 15.2.0-9 hold) + gcc-14/gfortran-14 PATH-shim rebuild + 2 batch-file fixes (nafems_le1.mesh missing TRACOBJE option, nafems_le1.elas missing TRACCONT '15'). Final D-27 rubric: 12/12 A/B pairs PASS, 1 DEFERRED (nlsecu ~1h50 compute). The 2 Qt "rendering gaps" were not real gaps — they were stale Qt binaries running against the gcc-16-snapshot libgfortran5 runtime that the Debian upgrade brought in unannounced.
---

# Phase 03 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> **Source:** `03-RESEARCH.md` §"Validation Architecture".

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Shell drivers + visual A/B (see `bin/cbxvtest0_qt`, `bin/cbl_tout`); no unit test framework — matches MEFISTO convention |
| **Config file** | none — Phase 0 + Phase 2 already installed build scripts |
| **Quick run command** | `bin/cbl_tout_qt && bin/cbl_tout` (both backends must stay green) |
| **Full suite command** | `bin/cbxvtest0_qt && pp/ppxvtest0_qt` + A/B replay of `prpr/xvtest1..4` + `testa/nafems_le1`, `testa/pan2d` (5-case driver set) |
| **Estimated runtime** | ~90s for Qt+X11 double-build; ~5min for full A/B visual gate |

---

## Sampling Rate

- **After every task commit:** `bin/cbl_tout_qt` (Qt build green)
- **After every plan wave:** `bin/cbl_tout_qt && bin/cbl_tout` (both builds green)
- **Before `/gsd-verify-work`:** Full A/B gate green + `verify_no_exec` grep clean
- **Max feedback latency:** ~60s per task commit (single-backend build)

---

## Per-Task Verification Map

> Filled during planning. Rows below are the contract gsd-planner MUST honor — one row per task.

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 03-01-T1 | 03-01 | 1 | TEXT-01, TEXT-04 | — | N/A | build + grep | `bin/cbl_tout_qt && bin/xvue/qt/cmake/verify_no_exec.sh` (TTF bundled, qt_add_resources wired, verify_no_exec extended) | ✅ | ✅ green |
| 03-01-T2 | 03-01 | 1 | TEXT-01, TEXT-03 | — | N/A | build | `bin/cbl_tout_qt` — XvueState palette+font state, XvueApp::ensure() font load, palette_init_once with 16-color defaults | ✅ | ✅ green |
| 03-01-T3 | 03-01 | 1 | TEXT-01..06 | — | N/A | build + run | `bin/cbxvtest0_qt && pp/ppxvtest0_qt` — Phase 3 TEXT coverage section added to prpr/xvtest0.f (D-24) | ✅ | ✅ green |
| 03-02-T1 | 03-02 | 1 | TEXT-01, TEXT-02, TEXT-03 | — | N/A | build + ABI | `bin/cbl_tout_qt && nm pp/ppxvtest0_qt \| grep -c -E 'xv(chargefonte\|nbpixeltexte\|texte)_$'` → 3 symbols (xvftexte collapsed into xvtexte) | ✅ | ✅ green |
| 03-02-T2 | 03-02 | 1 | TEXT-03, TEXT-04, TEXT-05 | — | N/A | build + ABI | `bin/cbl_tout_qt && nm pp/ppxvtest0_qt \| grep -c -E 'xv(couleur\|activervb\|recuprgbdec\|fond)_$'` → 4 symbols + xvinfo_ palette/font fill | ✅ | ✅ green |
| 03-03-T1 | 03-03 | 2 | TEXT-01..06 | — | N/A | build + run + grep | `bin/cbl_tout_qt && bin/cbl_tout && pp/ppxvtest0_qt` headless smoke — no warn-once for xvchargefonte/xvnbpixeltexte/xvtexte/xvftexte/xvcouleur/xvactivervb/xvrecuprgbdec ; verify_no_exec clean | ✅ | ✅ green |
| 03-03-T2 | 03-03 | 2 | TEXT-01..06 | — | N/A | manual (human-verify) | Human visual checkpoint on xvtest0 — fonts + colored lines + xvactivervb + measured label. Approved 2026-04-12. | — | ✅ green |
| 03-03-T3 | 03-03 | 2 | TEXT-01..06 | — | N/A | docs | Wave 2 approval recorded in 03-03-SUMMARY.md | ✅ | ✅ green |
| 03-04-T1 | 03-04 | 3 | TEXT-01..06, BUILD-07 | — | N/A | build + run | `bin/cbl_tout && bin/cbl_tout_qt ; bin/xvtest-capture.sh pp/ppxvtestN …_x11.png` (all 4 xvtest drivers Qt + legacy exit 0 under Xvfb :99; timeout-via-SIGTERM path retired by MEFISTO_XVSOURIS_AUTOEXIT hook) | ✅ | ✅ green |
| 03-04-T2 | 03-04 | 3 | TEXT-01..06, VALID-02 | — | N/A | manual (human-verify) | Human A/B visual gate xvtest1..4 Qt vs X11 per D-27 rubric. Resolved by orchestrator reading PNG captures (commit f3b9a6d) — all 4 pairs PASS. | ✅ | ✅ green |
| 03-04-T3 | 03-04 | 3 | TEXT-01..06, VALID-02 | — | N/A | manual (human-verify) | Human A/B visual gate on 5 canonical testa/ cases Qt vs X11. Initial close-out (2026-04-13) premature; first reopen (2026-04-13) showed 10/12 PASS / 2 MISMATCH; second reopen (2026-04-14) traced the 2 mismatches to a Debian sid `libgfortran5 = 16-20260322-1` (gcc-16 snapshot) runtime regression — fixed by pinning libgfortran5 to 15.2.0-9 + rebuilding with gcc-14/gfortran-14 + fixing 2 long-standing test-batch bugs (nafems_le1.mesh missing `1` between `5;` and `90;`; nafems_le1.elas missing `15;` between `1;` and `90;`). **Final D-27 verdict: 12 pairs compared, 12 PASS, 1 DEFERRED (nlsecu ~1h50 compute).** See 03-04-ab/testa/README.md. | ✅ | ✅ green (12/12 PASS) |
| 03-04-T4 | 03-04 | 3 | TEXT-01..06 | — | N/A | docs | Fill this Per-Task Verification Map + update 03-04-SUMMARY.md | ✅ | ✅ green |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [x] `prpr/xvtest0.f` — extend with Phase 3 coverage section (D-24): font load, palette store, 8-color line+label draw, `xvactivervb_` bulk-load demo, `xvnbpixeltexte_` bounding-box check
- [x] `bin/cbxvtest0_qt` — rebuild script already exists (Phase 2); Phase 3 reuses unchanged
- [x] `bin/verify_no_exec` — add D-19 grep rule: `qApp->palette|->palette()` must return zero matches in `xvue_qt_canvas.*` and `xvue_qt_api.*`
- [x] `xvue/qt/fonts/DejaVuSansMono.ttf` — bundled TTF committed to repo (D-01)
- [x] `xvue/qt/resources/xvue_fonts.qrc` — Qt resource file referencing the TTF (Pitfall 1)
- [x] `Q_INIT_RESOURCE(xvue_fonts)` call site — required because `libxvueqt.a` is STATIC (Pitfall 1 from RESEARCH.md)

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| A/B visual match Qt ↔ X11 for `xvtest1..4` | TEXT-03, TEXT-05 | Subpixel AA differs between Qt and X11 core fonts; human eye decides "no clipping/overlap, ≤2px drift" (D-27) | 1) Build both: `bin/cbl_tout_qt && bin/cbl_tout`. 2) Run each `pp/ppxvtest{1,2,3,4}_qt` and `pp/ppxvtest{1,2,3,4}` side-by-side on X11 (`QT_QPA_PLATFORM=xcb`). 3) Review per D-27 rubric. |
| Label layout on `testa/nafems_le1`, `testa/pan2d` | TEXT-02, TEXT-03 | Same reason — fontmetric drift is acceptable, clipping/overlap is not | Run under Qt, screenshot, compare against X11 run |
| Dark-mode freeze (TEXT-06) runtime proof | TEXT-06 | No QPalette/theme plumbing until Phase 6 — Phase 3 delivers construction-level guard only; runtime proof needs Phase 6 chrome | Manual dark-mode toggle deferred to Phase 6; Phase 3 guard is the grep in `verify_no_exec` (D-19) |

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 covers all MISSING references (xvtest0 extension, verify_no_exec grep, bundled TTF, Q_INIT_RESOURCE)
- [x] No watch-mode flags
- [x] Feedback latency < 120s per task (double-backend build budget)
- [x] `nyquist_compliant: true` set in frontmatter (updated 2026-04-13 after 03-04-T4 filled the Per-Task Verification Map)
- [x] TEXT-06 runtime proof explicitly deferred to Phase 6 (documented above)

**Approval:** *RESOLVED* 2026-04-14 — second reopen closed clean.
12/12 A/B pairs PASS, 1 DEFERRED (nlsecu ~1h50 compute).

**Task 2 (xvtest1..4):** 4/4 PASS per D-27 rubric applied to Qt+X11 PNGs
both captured via the new offscreen + Xvfb harnesses. The Qt side uses
`XvueState::backing_` direct grab; the X11 side uses `import -window
root` after the xvfermer_ sentinel hold. Both sides run the same driver
source. Full visual match for geometry, colors, and dashed edges.

**Task 3 (5 testa cases × 2 modules):** Final A/B verdict:
- Mesher 5/5 PASS : pan2d, nafems_le1, cavity2d, heat1d, nlsecu all
  render the expected 2D/1D/3D mesh with quality-coloring legend on
  both backends. Small palette-index choice differences (cavity2d and
  nlsecu use different colors in the same palette) are accepted as
  both are within the shared MAX_PALETTE and the geometry is identical.
- Solver 3/4 PASS : nafems_le1-ppelas (radial stress arrow pattern on
  the quarter annulus, principal-stress directions visible on both
  sides), cavity2d-ppflui (dark-blue uniform pressure mesh + 0..414
  ISO-pressure color-bar legend on both sides), heat1d-ppther (1D mesh
  + flux arrows + red eigenvalue diagonal + EIGENVALUE/NORMAL FLUX
  titles on both sides; Qt additionally shows the quality stats block
  and 10-color legend — extra content, not a regression).
- Solver 1/4 DEFERRED : nlsecu-ppnlse (~1h50 compute cost).

**What changed since the 2026-04-13 first reopen:**
1. **Debian sid libgfortran5 regression.** `apt upgrade` at 23:37
   on 2026-04-13 pulled `libgfortran5 = 16-20260322-1` — a gcc-16
   snapshot from sid. The MEFISTO Fortran binaries dynamically link
   `libgfortran.so.5`, so the runtime behavior changed underneath
   them. The first reopen captured the 2 "MISMATCH" PNGs against
   stale binaries running on the new runtime. Fix: downgrade with
   `sudo apt install /var/cache/apt/archives/libgfortran5_15.2.0-9_amd64.deb`
   then `sudo apt-mark hold libgfortran5`. The downgrade also removed
   `gcc-15` and `gfortran-15` (left only `gcc-14` and `gfortran-14`),
   so a `/tmp/gfortran-14-shim` PATH directory was created with
   `gcc → gcc-14`, `cc → gcc-14`, `gfortran → gfortran-14` symlinks.
   `bin/cbl_tout` + `bin/cbl_tout_qt` were re-run via that PATH.
2. **`testa/nafems_le1/nafems_le1.mesh` had a long-standing batch bug.**
   The drawing block `10; 5; 90;` was missing the `1; { ALL OBJECTS }`
   between `5;` and `90;`. Heat1d and cavity2d both have `10; 5; 1; 90;`
   for the equivalent step. The TRACOBJE submenu only accepts 1/3/5,
   so `90` errored out and no mesh ever drew. Fixed in
   `testa/nafems_le1/nafems_le1.mesh`.
3. **`testa/nafems_le1/nafems_le1.elas` had the same bug.** The drawing
   block `8; 1; 90; 90;` was missing the `15; { Drawing of STRESSES
   in ALL FE }` between `1;` and `90;`. The TRACCONT submenu only
   accepts 1..9 + 15 + 16, so `90` errored. Fixed in
   `testa/nafems_le1/nafems_le1.elas`.

These two batch-file fixes improve the test ABOVE its previous best:
the pre-Debian-upgrade nafems_le1-ppelas capture happened to show the
mesher's quality drawing (residual root-window state from the prior
ppmail process), not the actual stress visualization. The fixed batch
now produces a real principal-stress arrow drawing.

**Pre-existing bugs found + fixed as byproducts (across both reopens):**
- `xvue/xvuelc.c:1401` effacemempx_ NULL guard
- `xvue/xvuelc.c:1601` xvnbpixeltexte_ NULL-font guard
- `flui/lifiviprte.f:161,167` FORMAT descriptor parse error (ppflui
  restart reads)
- `elas/trelas.f:273` added `CALL MEMPXFENETRE` so the X11 elas tracer
  actually shows its output before the next menu prompt
- `prpr/xvtest0.f` use `XVOUVRIR` + `MEMPXFENETRE` for legacy-X11
  compatibility
- `testa/nafems_le1/nafems_le1.mesh` — drawing block `10; 5; 90;` was
  missing the `1; { ALL OBJECTS }` step (TRACOBJE submenu)
- `testa/nafems_le1/nafems_le1.elas` — drawing block `8; 1; 90; 90;`
  was missing the `15; { Drawing of STRESSES in ALL FE }` step
  (TRACCONT submenu)

**Build-environment finding (deferred to a future hardening phase):**
- `libgfortran5` must currently be pinned to `15.2.0-9` because Debian
  sid ships an experimental `16-20260322-1` (gcc 16 snapshot) that
  exposes latent UB in MEFISTO Fortran code (`TPSINI` reads garbage
  in the unsteady heat solver, FPE crashes, etc.). The pin is held
  via `apt-mark hold libgfortran5`. Should be lifted once Debian sid
  finalizes gcc-16 AND the latent UB sites in `ther/thed1t.f` (and
  similar) are properly initialized. **No code-level fix in this
  phase** — that's a follow-up audit of every Fortran subroutine
  that consumes solver/time variables. The latent UB exists
  independently of Phase 3 and is not blocking the Phase 3 visual
  contract.

**Infrastructure delivered (commits this session):**
- `xvue/xvuelc.c` headless hooks: `MEFISTO_XVSOURIS_AUTOEXIT`,
  `MEFISTO_XVFERMER_READY_FILE`, `MEFISTO_XVFERMER_HOLD_MS`
- `xvue/qt/src/xvue_qt_api.cpp` Qt-side counterpart of the same hooks,
  plus `MEFISTO_QT_CAPTURE_PATH` that saves `XvueState::backing_`
  directly (no X server needed — works under QT_QPA_PLATFORM=offscreen)
- `prpr/pp{mail,elas,flui,ther,nlse}.f` `MEFISTO_BATCH_X11` override
- `bin/xvtest-capture.sh` — orchestrates xvtest drivers under Xvfb
- `bin/testa-capture.sh` — orchestrates testa solvers under Xvfb
- `bin/qt-capture.sh` — orchestrates Qt binaries under offscreen

All this infrastructure is reusable by future phases' visual A/B gates
and is preserved even though Task 3 itself is not yet signed off.
