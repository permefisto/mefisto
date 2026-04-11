---
phase: 03
slug: text-fonts-colormap
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-04-11
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
| 03-XX-XX | XX | N | TEXT-0X | — | N/A (local rendering) | build | `bin/cbl_tout_qt` | ✅ | ⬜ pending |
| 03-XX-XX | XX | N | TEXT-0X | — | N/A | grep | `grep -nE 'qApp->palette\|->palette\(\)' xvue/qt/src/xvue_qt_canvas.* xvue/qt/src/xvue_qt_api.* ; test $? -eq 1` | ✅ | ⬜ pending |
| 03-XX-XX | XX | N | TEXT-01..06 | — | N/A | ABI diff | `nm pp/ppxvtest0_qt \| grep -E 'xv(couleur\|chargefonte\|texte\|ftexte\|recuprgbdec\|activervb\|nbpixeltexte)_'` returns exactly 7 symbols | ✅ | ⬜ pending |
| 03-XX-XX | XX | N | TEXT-01..06 | — | N/A | run | `pp/ppxvtest0_qt` — xvtest0 extended driver completes with no warn-once messages | ✅ W0 | ⬜ pending |
| 03-XX-XX | XX | N | TEXT-01..06 | — | N/A | A/B visual | Run `pp/ppxvtest1..4` with Qt and X11 backends, human review vs D-27 rubric | ✅ | ⬜ pending |
| 03-XX-XX | XX | N | TEXT-03, TEXT-04 | — | N/A | A/B visual | `pp/ppelas testa/nafems_le1 && pp/ppelas testa/pan2d` under Qt — label layout matches X11 within D-08 tolerance | ✅ | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `prpr/xvtest0.f` — extend with Phase 3 coverage section (D-24): font load, palette store, 8-color line+label draw, `xvactivervb_` bulk-load demo, `xvnbpixeltexte_` bounding-box check
- [ ] `bin/cbxvtest0_qt` — rebuild script already exists (Phase 2); Phase 3 reuses unchanged
- [ ] `bin/verify_no_exec` — add D-19 grep rule: `qApp->palette|->palette()` must return zero matches in `xvue_qt_canvas.*` and `xvue_qt_api.*`
- [ ] `xvue/qt/fonts/DejaVuSansMono.ttf` — bundled TTF committed to repo (D-01)
- [ ] `xvue/qt/resources/xvue_fonts.qrc` — Qt resource file referencing the TTF (Pitfall 1)
- [ ] `Q_INIT_RESOURCE(xvue_fonts)` call site — required because `libxvueqt.a` is STATIC (Pitfall 1 from RESEARCH.md)

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| A/B visual match Qt ↔ X11 for `xvtest1..4` | TEXT-03, TEXT-05 | Subpixel AA differs between Qt and X11 core fonts; human eye decides "no clipping/overlap, ≤2px drift" (D-27) | 1) Build both: `bin/cbl_tout_qt && bin/cbl_tout`. 2) Run each `pp/ppxvtest{1,2,3,4}_qt` and `pp/ppxvtest{1,2,3,4}` side-by-side on X11 (`QT_QPA_PLATFORM=xcb`). 3) Review per D-27 rubric. |
| Label layout on `testa/nafems_le1`, `testa/pan2d` | TEXT-02, TEXT-03 | Same reason — fontmetric drift is acceptable, clipping/overlap is not | Run under Qt, screenshot, compare against X11 run |
| Dark-mode freeze (TEXT-06) runtime proof | TEXT-06 | No QPalette/theme plumbing until Phase 6 — Phase 3 delivers construction-level guard only; runtime proof needs Phase 6 chrome | Manual dark-mode toggle deferred to Phase 6; Phase 3 guard is the grep in `verify_no_exec` (D-19) |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references (xvtest0 extension, verify_no_exec grep, bundled TTF, Q_INIT_RESOURCE)
- [ ] No watch-mode flags
- [ ] Feedback latency < 120s per task (double-backend build budget)
- [ ] `nyquist_compliant: true` set in frontmatter after planner fills the Per-Task Verification Map
- [ ] TEXT-06 runtime proof explicitly deferred to Phase 6 (documented above)

**Approval:** pending
