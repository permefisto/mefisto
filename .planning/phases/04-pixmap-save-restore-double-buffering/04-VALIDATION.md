---
phase: 04
slug: pixmap-save-restore-double-buffering
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-04-14
---

# Phase 04 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> **Source:** `04-RESEARCH.md` §"Validation Architecture".

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Shell drivers + visual A/B captures + pixel-diff (matches Phase 03 convention) |
| **Config file** | none — extends `prpr/xvtest0.f` and adds `bin/xvtest0-pixmap-roundtrip.sh` |
| **Quick run command** | `bin/cbl_tout_qt && pp/ppxvtest0_qt` (headless smoke via `QT_QPA_PLATFORM=offscreen`) |
| **Full suite command** | `bin/xvtest0-pixmap-roundtrip.sh` (Wave 0 creates this — invokes `bin/qt-capture.sh` N times + `magick compare -metric AE` pairwise + `bin/xvtest-capture.sh` for legacy X11 A/B) |
| **Estimated runtime** | ~30s for full suite, ~15–25s for quick run |

---

## Sampling Rate

- **After every task commit:** `bin/cbl_tout_qt && pp/ppxvtest0_qt` (rebuild Qt lib + headless xvtest0_qt run)
- **After every plan wave:** `bin/xvtest0-pixmap-roundtrip.sh` (round-trip pixel-diff) + `bin/cbl_tout` (legacy X11 build must stay green per CLAUDE.md "Compilation must never break")
- **Before `/gsd-verify-work`:** Full suite green + `verify_abi` reports 57 + `verify_no_exec` green + Qt+X11 A/B capture pair visually reviewed via `Read` tool (D-17)
- **Max feedback latency:** ~60s per wave merge

---

## Per-Task Verification Map

> Filled during planning. Rows below are the contract gsd-planner MUST honor — one row per task.

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| _planner-fills_ | _planner_ | _planner_ | PIXMAP-01 / 02 / 03 / 04 | — / N/A | N/A (read-only pixmap blits) | smoke / integration | _per Validation Architecture map below_ | ❌ Wave 0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

### Requirement → Test Map (from RESEARCH.md)

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| PIXMAP-01 | `fenetremempx_` / `mempxfenetre_` no-ops do not break the scene (single-backing collapse) | smoke | `bin/qt-capture.sh pp/ppxvtest0_qt /tmp/p4_fenetremempx.png` — assert exit 0 + PNG written | ❌ Wave 0 extends `prpr/xvtest0.f` with `CALL FENETREMEMPX` + `CALL MEMPXFENETRE` |
| PIXMAP-02a | `sauvefenetre_` + `restaurefenetre_` round-trip pixel-identical to control | integration | `magick compare -metric AE /tmp/p4_ctrl.png /tmp/p4_saverestore.png null:` — exit 0 | ❌ Wave 0 extends `xvtest0.f` + adds `bin/xvtest0-pixmap-roundtrip.sh` |
| PIXMAP-02b | Same test against legacy X11 backend | integration | `bin/xvtest-capture.sh pp/ppxvtest0 /tmp/p4_x11.png` + visual A/B vs Qt capture | ❌ Wave 0 — same harness |
| PIXMAP-03a | `sauvemempx_` / `restauremempx_` round-trip pixel-identical to control | integration | `magick compare -metric AE /tmp/p4_ctrl.png /tmp/p4_mempx_saverestore.png null:` — exit 0 | ❌ Wave 0 |
| PIXMAP-03b | `effacemempx_` clears to background | integration | `magick compare -metric AE /tmp/p4_bg.png /tmp/p4_effacemempx.png null:` — exit 0; `/tmp/p4_bg.png` = capture right after `xvinitgraphique` + `effacer` | ❌ Wave 0 |
| PIXMAP-03c | `effacemempx_` body byte-identical (or functionally identical) to `effacer_` body | unit (code review) | Manual diff of the two function bodies in `xvue_qt_api.cpp` + comment cross-reference | Built into Task DoD |
| PIXMAP-04 | Interactive cavity2d rubber-band-drag no-flicker | **DEFERRED to Phase 5** | — | Records `deferred-to-phase-5` per Phase 4 D-18 |

---

## Wave 0 Requirements

- [ ] `prpr/xvtest0.f` — extend with `PHASE 4 COVERAGE` block (D-15 scene + 3 sub-tests of D-16, bracketed by comment banners matching existing DRAW/TEXT coverage sections). Reads `MEFISTO_XVTEST0_SCENE` (or planner-chosen equivalent) env var via `CALL GETENV` for scene selection.
- [ ] `bin/xvtest0-pixmap-roundtrip.sh` — new ~50-line wrapper: invokes `pp/ppxvtest0_qt` under `bin/qt-capture.sh` for each sub-test scene with the env-var selector, then pairwise `magick compare -metric AE`. Exit code 0 iff all AE counts are 0.
- [ ] ImageMagick probe at script start: `command -v magick >/dev/null || { echo "magick not found"; exit 2; }` — already verified installed but keep the guard per CLAUDE.md "ask before acting."
- [ ] No changes to `bin/cbxvtest0_qt`, `bin/cbl_tout`, `bin/cbl_tout_qt`, or `xvue/qt/CMakeLists.txt` — Phase 4 reuses Phase 03-04 capture infrastructure unchanged.

*If no gaps surface during planning: "None — existing test infrastructure covers all phase requirements."*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Qt+X11 A/B visual match for the 3 round-trip captures | PIXMAP-02b | Subpixel font rendering and palette index choice differ between Qt's bundled DejaVu Sans Mono and X11 system fonts (Phase 03 close confirmed); strict pixel equality across backends is not the goal — visual match is | 1) Build both: `bin/cbl_tout_qt && bin/cbl_tout` 2) Run `bin/xvtest0-pixmap-roundtrip.sh` 3) Orchestrator reads each generated PNG via the `Read` tool and applies the D-27 rubric (geometry, colors, text, no missing geometry, no miscolored regions) |
| Interactive cavity2d rubber-band-drag no-flicker | PIXMAP-04 | Requires real mouse-motion events from the Phase 5 event bridge; cannot be validated headlessly in Phase 4 | **Deferred to Phase 5** — Phase 4 ships save/restore API + headless round-trip; Phase 5 owns the interactive HUMAN-UAT |

*All other phase behaviors have automated verification via the round-trip pixel-diff harness.*

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references (`prpr/xvtest0.f` PHASE 4 block, `bin/xvtest0-pixmap-roundtrip.sh`)
- [ ] No watch-mode flags
- [ ] Feedback latency < 60s per wave merge
- [ ] `nyquist_compliant: true` set in frontmatter (flipped after planner fills the Per-Task Verification Map and gsd-plan-checker confirms Dimension 8 green)
- [ ] PIXMAP-04 explicitly recorded as `deferred-to-phase-5`, not `green` and not `red`

**Approval:** *pending* — flips to *approved YYYY-MM-DD* after gsd-plan-checker passes the plan(s).
