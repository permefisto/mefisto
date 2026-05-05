# Phase 9: Retire X11 backend - Research

**Researched:** 2026-05-05
**Domain:** Build-script + Fortran source pruning, multi-architecture removal, plus 4 carry-forward defects
**Confidence:** HIGH (empirically audited; main-tree greps cross-referenced with Phase 7/8 evidence)

## Summary

Phase 9 retires the X11 / ImageMagick / LVIDEO graphics paths from MEFISTO and folds in the four
Phase-8 carry-forward items (matched-dim Qt recapture, ppnlse_qt offscreen+BATCH_X11 deadlock,
3 deferred Phase-7 goldens, harness `--out-dir` bug + CMake `verify_pp_qt_freshness` ALL target).
The work is mechanical — every deletion candidate is identifiable by an explicit grep — but the
**LVIDEO scope is wider than CONTEXT.md spelled out** and the **Debian package name in
REQUIREMENTS.md/CONTEXT.md is wrong**. Both are flagged below.

**Empirical totals (main tree, worktrees excluded):**
- 32 `bin/cb*` + `bin/xvt*` + `bin/Makefile*` scripts contain `-lX11` and/or `/usr/X11R6/lib*`
  linker references. Most are deletion candidates; ~10 are partial-edit candidates (already
  carry both X11 and Qt linker lines as siblings in `bin/cb*_qt`).
- `xvue/xvuelc.c` (3 749 lines, 156.7 KB) + `bin/ccxvue` (the wrapper) → RETIRE-01.
- `xvue/video1.f`, `xvue/videofin.f`, `xvue/videonm.f` (3 files, 211 lines total) + `bin/convertepsgif`
  (1 line: `convert -rotate 90 zfxy0*.eps -extent 980x550 cyl53zfxy.gif`) → RETIRE-03.
- **12** Fortran tracer files in `flui/` + `ther/` + `util/` contain embedded
  `IF (LVIDEO .NE. 0) THEN ... CALL VIDEO* ...` blocks. Each tracer has OTHER drawing logic
  surrounding its LVIDEO block — none is "LVIDEO-only". (CONTEXT.md D-05 said "18+"; the
  empirical figure is **12**, all selective edits per D-06.)
- 6 `xvue/*.f` files reference `LVIDEO` directly (`coudef.f` initializer, `vise2d.f` + `vise3d.f`
  menu setters, `video1.f` + `videofin.f` self-references, plus the COMMON declaration in
  `incl/trvari.inc`). Plus 4 menu-lexicon entries in `td/m/visee[23]d` + `td/ma/visee[23]d`.
- 3 README/LISEZMOI pairs (top-level `README` + `LISEZMOI`, plus duplicates under `bin/`)
  reference `libX11-dev` + `ImageMagick`. **`bin/README` and `bin/LISEZMOI` are byte-identical
  to the top-level versions** (identical greps, identical line numbers) — appear to be the
  same content shipped twice for the binary-distribution tarball.
- `bin/png2eps` and `bin/png2jpg` are TWO ADDITIONAL `convert`-shellout scripts NOT mentioned
  in CONTEXT.md or REQUIREMENTS.md §RETIRE-03; flagged below as a research finding.

**Primary recommendation:** Run RETIRE-NN in the order CONTEXT.md D-02 specifies
(01 → 02 → 03 → 04). Within RETIRE-02, refactor `bin/cbl_tout` to call `bin/cbl_tout_qt`
internals, then delete `bin/cbl_tout_qt` (rename via merge, not copy + rename — see §9).
For carry-forward 9-05, add a `MEFISTO_QT_WINDOW_SIZE=WxH` env hook in `xvinitgraphique_`
and `xvinfo_` (single-file edit in `xvue/qt/src/xvue_qt_api.cpp`). For 9-06, the deadlock
sits in `prpr/ppnlse.f`'s long-running TIME=20 loop; the actual fix is a Qt event-pump
gap during `tttsupa2d.f` iteration, NOT a Qt-startup hang — the "10 log lines no banner"
symptom is a misreading. See §6.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Delete `xvuelc.c` + `ccxvue` | Source files / build wrappers | — | Single C source + its compile script; nothing depends on the .o except legacy `bin/cb*` linker lines |
| Strip X11 linker lines | `bin/cb*` shell scripts | `bin/Makefile*` (legacy IBM/HP relics) | All linker invocations live in shell scripts; CMakeLists.txt for Qt has no X11 deps |
| Delete LVIDEO Fortran | `xvue/video[1,fin,nm].f` | `flui/`, `ther/`, `util/` tracers (selective edits) | The 3 SUBROUTINEs are self-contained; the 12 callers each gate via `IF (LVIDEO .NE. 0) THEN` and need surgery, not deletion |
| Strip ImageMagick from docs | Top-level `README`/`LISEZMOI` | `bin/README`/`bin/LISEZMOI` (duplicates), `CLAUDE.md` deps section | Doc scope is small, all references co-located with the libX11-dev sentence |
| Matched-dim Qt recapture (carry 9-05) | `xvue/qt/src/xvue_qt_api.cpp` (xvinitgraphique_, xvinfo_) | `xvue/qt/src/xvue_qt_window.cpp` (resize policy) | Window+canvas size negotiation lives in xvue/qt/; harness already exports needed env vars |
| ppnlse_qt deadlock fix (carry 9-06) | `xvue/qt/src/xvue_qt_event.cpp` OR a Fortran event-pump call inside `flui/tttsupa2d.f` | `prpr/ppnlse.f` (timeout policy) | Symptom is "no NLSER banner"; root cause is Qt event-loop starvation during tight Fortran iteration loops, NOT a startup hang (see §6) |
| 3 Phase-7 goldens (carry 9-07) | Cross-tag worktree (`v1.0-pre-retire`) | golden artifacts committed on main | Single deviation from "Phase 9 only touches Qt-only main"; confined to Plan 9-07 |
| `--out-dir` realpath fix (carry 9-08) | `bin/ab_sweep_phase8.sh` (single-line edit) | — | Pure shell script edit |
| `verify_pp_qt_freshness` (carry 9-08) | `xvue/qt/CMakeLists.txt` ALL target | `bin/cbl_tout_qt` invocation | Mirrors existing `verify_no_imagemagick_in_qt` pattern |

## Standard Stack

This is a deletion-and-cleanup phase. There is no "standard stack to install" — the relevant
libraries already exist on the dev host. The table below documents what stays vs goes.

### Tools that REMAIN after Phase 9

| Tool | Version | Purpose | Why It Stays |
|------|---------|---------|--------------|
| Qt 6 (qt6-base-dev) | 6.10.2+dfsg-7 [VERIFIED: `apt-cache show qt6-base-dev`] | Single graphics backend post-retirement | Phase 0–7 deliverable; canonical post-Phase-9 graphics |
| gfortran | 14.3.0 [VERIFIED: `gfortran --version`] | Fortran compiler | All Fortran sources stay |
| gcc | (system default) | C compiler — no longer linking xvuelc.c | Still needed for `-lstdc++` indirect link via Qt |
| ImageMagick `compare` / `convert` | 7.1.2-18 [VERIFIED: `convert --version`] | Phase 8 A/B validation tooling under `bin/ab_compare_pair.sh` | Per Phase 8 D-03 + Phase 7 EXPORT-06: ImageMagick stays for VALIDATION; only the `xvue/qt/` runtime forbids it |
| ffmpeg | 8.1-3+b1 [VERIFIED: `ffmpeg -version`] | Qt animated GIF backend (XvueExport saveGifTo) | Phase 7 EXPORT-03 deliverable; replaces LVIDEO+convertepsgif |
| Xvfb | (xvfb-run script available) [VERIFIED: `xvfb-run --help`] | Headless rendering for tests | Used by Qt ctest + 9-07 cross-tag golden bootstrap |
| cmake | 3.31.6 [VERIFIED: `cmake --version`] | xvue/qt/ build driver | Stays |

### Tools that GO after Phase 9

| Tool | Where Referenced | Action |
|------|------------------|--------|
| `libX11-dev` apt package | `README:83-87`, `LISEZMOI:84-89`, `bin/README:86-89`, `bin/LISEZMOI:88-91`, `CLAUDE.md:39` | Remove from doc dependency lists [VERIFIED: grep] |
| `-lX11` / `-lXt` linker flags | 32 `bin/cb*` and `bin/xvt*` scripts | Remove or delete script per RETIRE-02 [VERIFIED: grep] |
| `/usr/X11R6/lib64` and `/usr/X11R6/lib` paths | Same 32 scripts | Remove. Note: this directory **does not exist on Debian sid** — the linker path was already silently being ignored by gcc, which falls back to `/usr/lib/x86_64-linux-gnu/libX11.so` [VERIFIED: `ls /usr/X11R6/lib64/` returns only `modules/`] |
| `xvue/xvuelc.c` (3749 lines) + `xvue/xvuelc.o` | Source tree | Delete [VERIFIED: `ls -la xvue/xvuelc.c`] |
| `bin/ccxvue` (44-line C compile wrapper) | bin/ | Delete [VERIFIED: `cat bin/ccxvue`] |
| `xvue/video1.f`, `xvue/videofin.f`, `xvue/videonm.f` | xvue/ | Delete [VERIFIED: `wc -l`] |
| `bin/convertepsgif` (1-line wrapper) | bin/ | Delete [VERIFIED: `cat bin/convertepsgif`] |
| `bin/png2eps` + `bin/png2jpg` (additional `convert` shell-outs) | bin/ | **NEW FINDING** — flag for planner; both invoke `convert`. Decision needed: delete (consistent with RETIRE-03) or migrate to a Qt/ffmpeg equivalent (out of scope for v1)? |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Delete `bin/cbl_tout_qt`, rename `cbl_tout` to be Qt-only | Keep both, edit `cbl_tout` to call `cbl_tout_qt` body inline | Single-name convention is cleaner; downstream tooling (CLAUDE.md, MEMORY.md, `xvue/qt/README.md:24`) already uses `bin/cbl_tout` as the canonical name; Phase 7 + 8 hand-off docs reference `bin/cbl_tout_qt` for build-hygiene gates — those references must be updated in lockstep |
| Selective LVIDEO surgery in 12 tracers | Bulk-delete tracers + remove menu options | Tracers do non-LVIDEO work too (mesh/contour drawing); bulk delete would break interactive `visee2d`/`visee3d` workflows. Selective surgery preserves D-06 contract |
| In-process `MEFISTO_QT_WINDOW_SIZE` env (9-05) | Hardcode 1280x800 | Environment knob is reversible; hardcoding regresses interactive UX (T-09-05) |

**Installation:** No new packages. Dev host already has all dependencies (verified above).

**Version verification (already done):**
- Qt 6: `pkg-config --modversion Qt6Widgets` → `6.10.2` [VERIFIED 2026-05-05]
- ffmpeg: `ffmpeg -version` → `8.1-3+b1` [VERIFIED 2026-05-05]
- cmake: `cmake --version` → `3.31.6` [VERIFIED 2026-05-05]

## Architecture Patterns

### System Architecture Diagram

Pre-Phase-9 (current):

```
                     bin/cbl_tout                bin/cbl_tout_qt
                          │                            │
                          ▼                            ▼
                  bin/ccxvue → xvue/xvuelc.o   cmake build → libxvueqt.a
                          │                            │
                          ▼                            ▼
       bin/cb{init,mail,elas,flui,ther,nlse}    bin/cb{mail,elas,flui,ther,nlse}_qt
                  (-lX11 -L/usr/X11R6/lib64)       (-lxvueqt -lQt6Widgets -lQt6Gui ...)
                          │                            │
                          ▼                            ▼
                pp/pp{init,mail,elas,...}      pp/pp{mail,elas,...}_qt
                          │                            │
                          ▼                            ▼
                 5 testa/ baseline cases      5 testa/ baseline cases
                  (X11 visual A/B)              (Qt offscreen capture)
```

Post-Phase-9 (target):

```
                            bin/cbl_tout
                                 │  (Qt-only; what was cbl_tout_qt)
                                 ▼
                        cmake build → libxvueqt.a
                                 │
                                 ▼
              bin/cb{mail,elas,flui,ther,nlse}     (no _qt suffix; renamed in lockstep)
                                 │
                                 ▼
                        pp/pp{mail,elas,...}        (no _qt suffix)
                                 │
                                 ▼
                       5 testa/ baseline cases
                              (Qt only)
```

### Recommended Project Structure (post-Phase-9)

Tree shape stays unchanged. Files removed:

```
xvue/
├── (every *.f file: STAYS — these are the Fortran wrappers, not X11 code)
├── xvuelc.c         # DELETED (RETIRE-01)
├── xvuelc.o         # DELETED (build artifact, .gitignored anyway)
├── video1.f         # DELETED (RETIRE-03)
├── videofin.f       # DELETED (RETIRE-03)
├── videonm.f        # DELETED (RETIRE-03)
├── coudef.f         # EDITED  (drop LVIDEO=0 init at line 13)
├── vise2d.f         # EDITED  (drop label 8800 at line 782 + GOTO at line 39)
├── vise3d.f         # EDITED  (drop label 8800 at line 798 + GOTO)
└── qt/              # ALL STAYS — sole graphics backend

incl/
└── trvari.inc       # EDITED  (drop LVIDEO from COMMON / TRVAPS / at line 9)

flui/
├── parpartr.f       # EDITED  (delete VIDEO* call block; preserve KNOMFGIF declaration if used elsewhere)
├── trvi2d.f         # EDITED
├── trvi3d.f         # EDITED
└── tttsupa2d.f      # EDITED

ther/
├── trisot.f         # EDITED
├── trlldr.f         # EDITED
├── trplse.f         # EDITED
├── trso1so.f        # EDITED
├── trzont.f         # EDITED
└── trztxy.f         # EDITED

util/
├── trtable.f        # EDITED  (delete IF LVIDEO block at lines 160+)
└── trtables.f       # EDITED  (similar)

td/m/
├── visee2d          # EDITED  (drop ',88: '''CONSTRUIRE un FILM VIDEO.gif''')
└── visee3d          # EDITED

td/ma/
├── visee2d          # EDITED  (drop ',88: '''Make a VIDEO.gif''')
└── visee3d          # EDITED

bin/
├── cbl_tout         # REWRITTEN to be the Qt-only build (absorbs cbl_tout_qt body)
├── cbl_tout_qt      # DELETED after merge into cbl_tout (D-09 below)
├── ccxvue           # DELETED (RETIRE-01)
├── convertepsgif    # DELETED (RETIRE-03)
├── png2eps          # DELETED (RETIRE-03 NEW FINDING)
├── png2jpg          # DELETED (RETIRE-03 NEW FINDING)
├── cb{init,mail,elas,flui,ther,nlse}            # DELETED (RETIRE-02 — replaced by _qt counterparts)
├── cb{mail,elas,flui,ther,nlse}_qt              # RENAMED to drop _qt suffix
├── cbxvtest{0,1,2,3,4}                           # DELETED (RETIRE-02 — _qt counterparts stay, then renamed)
├── cbxvtest{0,1,2,3,4}_qt                        # RENAMED to drop _qt suffix
├── cbpoba                                         # EDITED (drop -lX11 / /usr/X11R6 paths)
├── cbpara, cbgpara, cbgparaddd                    # EDITED or DELETED (legacy parallel/MPI; not in active build chain — confirm with maintainer)
├── cbadap, cbgadap, cbpppl, cbonde,              # All have X11 linker lines; status unclear (legacy)
│   cbbrezfort{2d,3d}, cbg{init,init1,mail,...}    # — likely DELETE (none invoked from cbl_tout)
├── xvt{1,2,3,4}                                   # DELETED (legacy X11-test launchers)
├── ajax.f                                         # Inspect — file in bin/ tagged in xvuelc grep but content TBD
├── avecMOTIF                                      # Legacy comment-only file referencing X11; DELETE
├── instal_2disk.{hp,src}                          # Legacy HP installer; DELETE
├── Makefile, MakefileIBM, MakefileMefisto         # Legacy IBM/HP/early-Linux Makefiles. DELETE.
├── ab_*.sh, qt-capture.sh, testa-capture.sh,     # STAY — Phase 8 harness
│   xvtest-capture.sh, xvtest0-pixmap-roundtrip.sh # STAY — uses ImageMagick LEGITIMATELY
└── test_no_imagemagick_in_qt.sh                  # STAY — keeps the EXPORT-06 grep gate alive
```

### Pattern 1: Mechanical-Deletion-with-Atomic-Commit

**What:** Each RETIRE-NN plan is one or more atomic commits, each of which leaves
`bin/cbl_tout` (post-rename, the Qt-only build) green. Failure to compile = revert that commit.
**When to use:** All 4 RETIRE-NN plans.
**Example pattern (no specific code site):**
- Step 1: stage deletions for the file(s) targeted by this commit.
- Step 2: run `bin/cbl_tout` (post-rename: just the Qt build); confirm exit 0.
- Step 3: run the 5 BUILD-10 testa cases through Phase-8 harness; confirm 5/5 captures land.
- Step 4: commit. If any step fails, `git restore -SW .` and try smaller scope.

### Pattern 2: Selective-LVIDEO-Block-Excision (D-06)

**What:** For each tracer file in `flui/`, `ther/`, `util/`:
1. Locate the `IF (LVIDEO .NE. 0) THEN ... ENDIF` block (or pattern variant
   `IF (LVIDEO .NE. 0) THEN <one CALL VIDEO* line>` without explicit `ENDIF`).
2. Delete the block and any local-variable declarations that become unused.
3. Verify: per-file Fortran compile via `bin/cbl_tous_f` (`echo flui | bin/cbl_tous_f` etc.).
4. Verify: surrounding tracer logic still produces correct mesh/contour output via
   targeted testa case (e.g., `flui/trvi2d.f` is exercised by `cavity2d`).

**When to use:** RETIRE-03 LVIDEO scope across the 12 tracer files.

### Pattern 3: Cross-Tag Worktree for Legacy Asset Production (Plan 9-07)

**What:** Plan 9-07 must produce 3 Phase-7-deferred goldens
(`scene01.eps`, `wave_legacy.gif`, `cavity2d_legacy.gif`) AFTER X11 is gone from main.
Procedure:
```bash
# In a separate location, NOT inside main worktree
git worktree add /tmp/mefisto-pre-retire v1.0-pre-retire
cd /tmp/mefisto-pre-retire
export MEFISTO=$PWD MEFISTOX=/tmp/mefistox-pre-retire
bin/cbl_tout    # Builds X11 binaries from the pre-Phase-9 tree

# Bootstrap scene01.eps (per Phase 7 VERIFICATION.md §9 procedure):
gfortran -I$MEFISTO/incl -c xvue/qt/tests/golden/scene01_driver.f
gfortran scene01_driver.o $MEFISTO/xvue/xvuelc.o $MEFISTO/xvue/*.o \
    $MEFISTO/util/*.o -L/usr/X11R6/lib -lX11 -lXt -o scene01_x11
xvfb-run --auto-servernum ./scene01_x11
cp TEMPORAIRE.EPS /home/mefisto/git/mefisto/xvue/qt/tests/golden/scene01.eps

# Bootstrap wave_legacy.gif + cavity2d_legacy.gif via legacy bin/convertepsgif chain:
# (See Phase 7 VERIFICATION.md §9 for canonical procedure.)

cd /home/mefisto/git/mefisto
git add xvue/qt/tests/golden/{scene01.eps,wave_legacy.gif,cavity2d_legacy.gif}
git commit -m "test(09-07): commit 3 Phase-7-deferred goldens via cross-tag bootstrap"
git worktree remove /tmp/mefisto-pre-retire
```

**When to use:** Once and only once, in Plan 9-07. Pitfall T-09-03 mitigation: ensure no
intermediate build artifacts (e.g., `.o` files, `pp/` binaries) leak from the worktree into
the main checkout via `git status` clean check pre-commit.

### Pattern 4: Env-Knob Window Sizing for Matched-Dim Recapture (Plan 9-05)

**What:** Add a single env-var hook to `xvue/qt/src/xvue_qt_api.cpp::xvinitgraphique_`
(or `xvinfo_`, before `win->resize(*ix, *iy)` at line 486):

```cpp
// Source: empirical analysis of xvue_qt_api.cpp:476-487 + Phase 8 evidence sweep-log-x11.md
// Fortran callers pass *ix = LAPXFE = 800, *iy = LHPXFE = 600 (hard-coded in prpr/pp*.f
// when EXISTNF != 0). Window chrome (menubar+toolbar+statusbar+console dock) eats into the
// canvas, so the captured backing pixmap lands at e.g. 760x442 not 800x600. To match the
// X11 baseline's 1280x800 (Xvfb root-window screenshot), allow override via env.
const char* qt_window_size = std::getenv("MEFISTO_QT_WINDOW_SIZE");
if (qt_window_size && qt_window_size[0]) {
    int w, h;
    if (std::sscanf(qt_window_size, "%dx%d", &w, &h) == 2 && w > 0 && h > 0) {
        // Resize the CANVAS, not the window — chrome adapts around it.
        // Easiest: call win->canvas()->setFixedSize(w, h); then win->adjustSize();
        // so the window grows to accommodate. Phase 8 baseline = 1280x800.
        if (win && win->canvas()) {
            win->canvas()->setFixedSize(w, h);
            win->adjustSize();
        }
    }
}
```

**When to use:** Plan 9-05 (matched-dim recapture). Once shipped, harness sets
`MEFISTO_QT_WINDOW_SIZE=1280x800` in the qt-1x dispatch path; AE compare runs against
1280x800 X11 baseline without resample dim-mismatch confound; CHECK cells flip toward PASS.

**Pitfall:** `setFixedSize` will lock interactive resize. Use `setMinimumSize +
setMaximumSize` to lock only when env is set; clear constraints on `xvfermer_`. Or only
honor the env when `MEFISTO_BATCH_X11=1` AND `QT_QPA_PLATFORM=offscreen` (the harness
context).

### Anti-Patterns to Avoid

- **Bulk-delete LVIDEO tracers:** Per D-06, the 12 tracer files do non-LVIDEO drawing too.
  Bulk delete = breaking interactive `visee[23]d` workflow on testa cases that don't reach
  Phase 8's 5-baseline subset. Use per-block excision.

- **Hardcoding 1280x800 in xvue/qt/:** Regresses interactive UX (T-09-05). The env-knob
  pattern is reversible per process — interactive sessions stay at the existing 1024x768
  default; harness flips to 1280x800 only when env is set.

- **Renaming `cbl_tout_qt` by `mv` instead of `git mv`:** Loses git history continuity for
  the Qt build script. Use `git mv bin/cbl_tout_qt bin/cbl_tout` AFTER deleting the legacy
  `cbl_tout` in the same commit.

- **Bumping ABI count from 58 to 59 to add a window-size entry point:** T-09-06 mitigation —
  ABI must stay at 58 symbols (Phase 6.0 D-XX locked this). Use an env var, NOT a new
  `extern "C"` entry. (Phase 8 verifier ran `verify_abi.sh` → "nm count: 58 header count: 58
  exit 0" — that gate stays green only if Plan 9-05 doesn't add symbols.)

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Animated GIF from frame sequence | New ffmpeg wrapper | `XvueExport::saveGifTo` (Phase 7 EXPORT-03 deliverable) | Already shipped, A/B-validated in Phase 8, replaces LVIDEO+convertepsgif by design (D-07) |
| PostScript byte-output | New emitter | `PsEmitter::handleLasops` (Phase 7 EXPORT-04 deliverable) | Verbatim port from xvuelc.c:1187-1304, byte-parity-tested |
| Image format probing for build-time targets | Bash test script | `xvue/qt/probes/qimagewriter_probe.cpp` (Phase 7 EXPORT-01) | Records output to PROBE.md, decides ffmpeg vs QImageWriter strategy |
| Test that `xvue/qt/` doesn't gain X11 code post-retirement | New CMake target | Existing `verify_no_imagemagick_in_qt` pattern at `xvue/qt/CMakeLists.txt:154-170` (DEPENDS xvueqt) | Already proven in Phase 7; clone the pattern with the new grep regex |
| `pp/*_qt` mtime check | Inline `find` | Cmake `add_custom_target(verify_pp_qt_freshness ALL ... DEPENDS xvueqt)` | Uniform ALL-target failure mode; same shape as `verify_no_exec`, `verify_shortcut_modifiers`, `verify_icon_source`, `verify_no_imagemagick_in_qt` |
| `--out-dir` realpath sanitization | Wrap pushd | Single `realpath -m` invocation early in script | Standard idiom; bash-portable |

**Key insight:** Phase 9 inherits a mature toolchain. The "don't hand-roll" rule mostly says
"reuse the Phase 6/7 infrastructure" — animated GIF, PostScript, image-format probes,
build-time grep gates — all already shipped and battle-tested.

## Runtime State Inventory

This phase is a deletion + refactor of a SOURCE TREE — it does NOT mutate stored data,
external services, or OS-registered processes. The inventory below is the explicit answer
to each category from the rename/refactor checklist:

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | **None** — MEFISTO does not use any database. The only persistent state is `$MEFISTOX/<project>/` (user project files) which contain mesh data, NOT codebase identifiers; no string substitution needed there. Verified by absence of any DB driver / ORM / SQL / Mem0 / ChromaDB / Redis reference in `find . -name 'CLAUDE.md' -o -name '*.cpp' -o -name '*.f' \| xargs grep -l 'redis\|mongo\|sqlite' 2>/dev/null`. | None |
| Live service config | **None** — MEFISTO is a single-process desktop app, no n8n / Datadog / Cloudflare / external service. Single-machine scientific computing tool. Verified by absence of service references in `bin/`, `xvue/`, `prpr/`. | None |
| OS-registered state | **None** — no systemd unit, no launchd plist, no Task Scheduler entry, no pm2 saved config. The `INITIER`/`MAILLER`/`ELASTICER` launchers in `bin/` are bash scripts run interactively from a shell; they don't register with the OS. Verified by `find /etc/systemd /lib/systemd -name '*mefisto*' 2>/dev/null` (returns nothing) and `crontab -l 2>/dev/null` (no MEFISTO refs). | None |
| Secrets / env vars | The codebase reads `MEFISTO`, `MEFISTOX`, `CDPATH`, `MEFISTO_BATCH_X11`, `MEFISTO_QT_CAPTURE_PATH`, `MEFISTO_XVFERMER_READY_FILE`, `MEFISTO_XVFERMER_HOLD_MS`, `MEFISTO_XVSOURIS_AUTOEXIT`, `MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS`. Phase 9 carry-forward 9-05 ADDS `MEFISTO_QT_WINDOW_SIZE`. None of these change semantics under retirement; the additions are backward-compatible (default = no-op when unset). | Update `xvue/qt/README.md` to document the new env var. No existing secret/config gets renamed. |
| Build artifacts / installed packages | **3 items.** (1) `xvue/xvuelc.o` — stale C object for retired source; deletion in RETIRE-01 covers it. (2) `pp/pp{init,mail,elas,flui,ther,nlse,xvtest0,xvtest1,xvtest2,xvtest3,xvtest4}` — the X11-built executables. After RETIRE-02 these become orphans (no producer); they get clobbered when `bin/cbl_tout` (post-rename) builds. Recommend explicit `rm -f $MEFISTO/pp/pp{init,xvtest0,xvtest1,xvtest2,xvtest3,xvtest4}` step at the start of RETIRE-02 since these names will not be regenerated. (3) `xvue/qt/build/` — CMake build dir; `bin/cbl_tout_qt` already does `rm -rf` on it before each build (line 51), so no action needed. | rm explicit step in RETIRE-02 |

**Nothing found in categories 1–3:** Stated explicitly. MEFISTO's codebase is self-contained
scientific software with no external service dependencies. Confirms the retirement is purely
a source-tree refactor.

## Common Pitfalls

### Pitfall 1: Ghost references to xvuelc symbols in Fortran COMMON blocks (T-09-01)

**What goes wrong:** A Fortran source file references a symbol that gets resolved at link
time only via `xvuelc.o`. When `xvuelc.o` goes away, the link breaks.

**Why it happens:** `xvuelc.c` contains ~60 `extern "C"` symbols. Fortran callers reference
them by Fortran name (lowercase + trailing underscore). The Qt build provides identical
symbols via `libxvueqt.a`. Phase 8 verified ABI count = 58 on both sides.

**Empirical verification (already done):**
- `nm xvue/lib | wc -l` shows ~224 Fortran wrappers; none collide with xvuelc symbol names.
- `grep -rn 'xvuelc' xvue/*.f xvue/*.inc incl/*.inc` returns 0 matches in Fortran sources
  (only xvue/xvuelc.c references its own filename in a comment line).
- Plan 9-01's verify step: after `rm xvue/xvuelc.c xvue/xvuelc.o`, run
  `bin/cbl_tout` (post-rename) and confirm clean link. CMake's xvueqt target already
  provides every symbol via `verify_abi.sh`.

**How to avoid:** Plan 9-01 must include a post-deletion rebuild + 5-case testa run before
commit.

**Warning signs:** Linker error `undefined reference to xvtrait_` (or similar) → the Qt
side is not picking up the symbol; check `nm xvue/qt/build/libxvueqt.a | grep xvtrait_`.

### Pitfall 2: Tracer subroutine that ALSO serves non-LVIDEO purpose (T-09-02)

**What goes wrong:** Plan deletes a whole tracer file thinking it's LVIDEO-only; mesher /
solver loses an interactive `visee[23]d`-driven contour-plot or mesh-overlay capability.

**Why it happens:** The 12 tracer files (`flui/parpartr.f`, `flui/trvi[23]d.f`,
`flui/tttsupa2d.f`, `ther/tr{isot,lldr,plse,so1so,zont,ztxy}.f`, `util/trtable[s].f`) each
contain BOTH `IF (LVIDEO .NE. 0) THEN CALL VIDEO* ... ENDIF` blocks AND lots of other
drawing logic. Per-file CLOC distribution: video-callout block ≈ 5 lines; surrounding
tracer body ≈ 200–1000 lines.

**How to avoid:** Per-block excision (Pattern 2 above). Verify each file compiles cleanly
post-edit AND its principal tracer behavior survives via spot-check on cavity2d/heat1d/etc.
testa cases.

**Warning signs:** Test case suddenly shows blank canvas where the X11 path showed mesh
overlays → an embedded drawing call got deleted along with the LVIDEO block.

### Pitfall 3: Cross-tag golden production leaks artifacts to main worktree (T-09-03)

**What goes wrong:** Plan 9-07 builds X11 binaries in `/tmp/mefisto-pre-retire` worktree.
On accidental `cd` back to main without `git worktree remove`, the new build artifacts
(.o files, pp/* binaries) might be confused with main-tree files — or worse, the user
runs `bin/cbl_tout` in main while it still has the X11 path, producing now-defunct
artifacts.

**How to avoid:**
1. Always use a NON-`$MEFISTO` directory for the cross-tag worktree (e.g., `/tmp/...`).
2. Set `MEFISTO=/tmp/mefisto-pre-retire` for the duration of the worktree session;
   don't leak the env var.
3. After committing the 3 goldens, `git worktree remove /tmp/mefisto-pre-retire`.
4. Pre-commit: run `git status` in main; only the 3 golden files should be staged. No
   `.o`, no `pp/*` binaries.

### Pitfall 4: Renaming `cbl_tout_qt` to `cbl_tout` breaks downstream tooling (T-09-04)

**What goes wrong:** Various docs / CI / tooling reference `bin/cbl_tout_qt` as the Qt
build entry. Phase 7 + 8 plans already do (e.g., `xvue/qt/README.md:24`,
`08-CONTEXT.md` D-08). After the rename, those references rot.

**How to avoid:**
1. Run `grep -rn 'cbl_tout_qt' . --include='*.md' --include='*.txt' --include='*.sh' --include='*.cpp'`
   BEFORE the rename (this is a planning task in Plan 9-02 or 9-08).
2. Update each reference in the same commit.
3. Add a CMake-build-time grep gate (similar to `verify_no_imagemagick_in_qt`) that
   fails if `cbl_tout_qt` reappears in `xvue/qt/`. (Optional; the codebase will not
   regress unless someone manually re-introduces it.)

### Pitfall 5: Matched-dim env knob accidentally regresses interactive UX (T-09-05)

**What goes wrong:** Plan 9-05's `MEFISTO_QT_WINDOW_SIZE` env knob, if implemented as
`setFixedSize`, locks the canvas at 1280x800 even in interactive sessions if the user
has the env var set in their shell.

**How to avoid:**
1. Gate the env-var honoring on EITHER `MEFISTO_BATCH_X11=1` OR `QT_QPA_PLATFORM=offscreen`
   (the headless contexts).
2. Use `setMinimumSize + setMaximumSize` (NOT `setFixedSize`) so the constraint can be
   cleared on `xvfermer_`.
3. Document the env var as "headless-only" in `xvue/qt/README.md`.

### Pitfall 6: ppnlse_qt deadlock fix changes the ABI (T-09-06)

**What goes wrong:** Plan 9-06 adds a new `extern "C"` entry point to pump the Qt event
loop from inside Fortran's iteration loop. ABI count goes from 58 to 59. Phase 8's
`verify_abi.sh` ALL target fails on next build.

**How to avoid:** Two safe approaches:
1. **C-only fix:** Add the event-pump call inside an EXISTING ABI entry that's already
   called every iteration (e.g., `xvvoir_`, `mempxfenetre_`). No ABI change.
2. **Add the entry but update verify_abi.sh in same commit:** Treat the new entry as a
   deliberate ABI v2 bump from 58 to 59. Update header + verify script + COORDS doc in
   lockstep. Phase 6.0 already did this once when bumping 57→58 for `xvue_module_init_`.

Recommend (1) for minimum invasiveness. Most likely call site:
`xvue/qt/src/xvue_qt_api.cpp::xvvoir_` (already pumps `processEvents` per its semantics).

### Pitfall 7: Documentation update misses one of the 4 README/LISEZMOI files (T-09-07)

**What goes wrong:** Plan 9-04 updates top-level `README` + `LISEZMOI` but forgets that
identical content is duplicated under `bin/README` + `bin/LISEZMOI` (binary-distribution
tarball assets).

**How to avoid:**
- Run `grep -ln 'libX11-dev\|imagemagick\|ImageMagick\|convert.*animated' README LISEZMOI bin/README bin/LISEZMOI install.bash bin/install.bash CLAUDE.md`
  → expect FIVE matches: `README`, `LISEZMOI`, `bin/README`, `bin/LISEZMOI`, `CLAUDE.md`
  (line 39).
- Plan 9-04 task list must enumerate all 5.

### Pitfall 8: NEW FINDING — `bin/png2eps` and `bin/png2jpg` are silent ImageMagick consumers

**What goes wrong:** Phase 9 retires "ImageMagick `convert`" but leaves two helper scripts
that invoke `convert` for ad-hoc image-format conversion (`png2eps`, `png2jpg`). They are
not in `bin/cbl_tout` build chain and not in the Phase-8 harness, but they ARE in
`bin/test_no_imagemagick_in_qt.sh` exclusion scope (`xvue/qt/` only). After Plan 9-04 says
"ImageMagick removed from install docs", these scripts will silently break on machines
that no longer have ImageMagick.

**How to avoid:** Plan 9-04 must explicitly choose ONE of:
- (a) Delete `bin/png2eps` + `bin/png2jpg` (no Qt-equivalent in v1; defer to v1.x).
- (b) Update doc to say "ImageMagick OPTIONAL — required only for `bin/png2{eps,jpg}`
  helpers" (preserves the helpers but contradicts CONTEXT.md "ImageMagick is dropped from
  the install documentation").
- (c) Replace with Qt-native `XvueExport::saveJpegTo` shell-out via a new `bin/png2jpg_qt`
  (over-scope for Phase 9).

Recommend (a) for cleanest retirement.

### Pitfall 9: NEW FINDING — Debian Qt6 image-formats package name in REQUIREMENTS.md / CONTEXT.md is WRONG

**What goes wrong:** REQUIREMENTS.md §RETIRE-04 line 112 + 09-CONTEXT.md "In scope"
list `libqt6imageformats6-plugins`. Empirical Debian package name is
`qt6-image-formats-plugins` (3 hyphens, no `6` suffix).

**Verification:** `apt-cache search '^libqt6imageformats'` → 0 results.
`apt-cache search '^qt6-image-formats'` → `qt6-image-formats-plugins - Qt 6 Image Formats plugins`. [VERIFIED 2026-05-05].

**How to avoid:** Plan 9-04 must update doc to use the correct name. Recommended
dependency line: `qt6-base-dev qt6-image-formats-plugins`. Note: `qt6-base-dev` already
provides `Qt6Widgets`, `Qt6Gui`, `Qt6Core`, `Qt6PrintSupport`, `Qt6Test` — only the
extra image-format plugins (TIFF, WebP, etc.) need the second package.

## Code Examples

Verified patterns from official sources:

### Example 1: Tracer file LVIDEO block excision (`util/trtable.f`)

Before (lines 154–179, currently in main):
```fortran
C     COPIE DE MEMPX DANS FENETRE
      CALL MEMPXFENETRE

C     POUR VIDER LE BUFFER DE X11
      CALL XVVOIR

      IF( LVIDEO .NE. 0 ) THEN

C     A PARTIR DE LA PIXMAP de la FENETRE X11 ACTUELLE CREATION DU FICHIER.xwd
      L = NUDCNB( NMFVIDEO )
      print *
      print *, 'xwd -xy -name Mefisto -out '//NMFVIDEO(1:L)//'.xwd',
     %        ' EST EXECUTE'
      CALL SYSTEM('xwd -xy -name Mefisto -out '//NMFVIDEO(1:L)//'.xwd')

C     CONVERSION DU FICHIER.xwd en le FICHIER.jpg
      LL = NUDCNB( NMPROJ )
      print *,    'convert ' // NMFVIDEO(1:L) // '.xwd' // ' '
     %                       // NMPROJ(1:LL)  // '_' // NMFVIDEO(1:L) //
     %                     '.jpg  EST EXECUTE'
      CALL SYSTEM('convert ' // NMFVIDEO(1:L) // '.xwd' // ' '
     %                       // NMPROJ(1:LL)  // '_' // NMFVIDEO(1:L) //
     %                     '.jpg')

C     DESTRUCTION DU FICHIER .xwd
      CALL SYSTEM( 'rm -Rf ' // NMFVIDEO(1:L) // '.xwd' )
```

After (lines 154–158, post-Plan 9-03):
```fortran
C     COPIE DE MEMPX DANS FENETRE
      CALL MEMPXFENETRE

C     POUR VIDER LE BUFFER DE X11
      CALL XVVOIR
```

Source: `util/trtable.f:150-180` [VERIFIED via `Read util/trtable.f:150-180`].

### Example 2: CMake `verify_pp_qt_freshness ALL` target

Add to `xvue/qt/CMakeLists.txt` after the `verify_no_imagemagick_in_qt` target
(around line 170):

```cmake
# Phase 9 carry-forward 9-08 (Phase 8 D-09): freshness guard. After cmake builds
# libxvueqt.a, every pp/pp*_qt executable must have mtime ≥ libxvueqt.a mtime.
# If any pp/pp*_qt is older than libxvueqt.a, the maintainer forgot to re-link
# and the next test sweep will use a stale binary.
add_custom_target(verify_pp_qt_freshness ALL
    COMMAND ${CMAKE_COMMAND} -E echo "verify_pp_qt_freshness: scanning..."
    COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/cmake/verify_pp_qt_freshness.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/build/libxvueqt.a
        ${CMAKE_CURRENT_SOURCE_DIR}/../../pp
    DEPENDS xvueqt
    VERBATIM
)
```

With companion `xvue/qt/cmake/verify_pp_qt_freshness.sh`:

```bash
#!/usr/bin/env sh
# Phase 9 Plan 9-08 (Phase 8 D-09): Fail if libxvueqt.a is newer than ANY pp/pp*_qt.
# After Phase 9 RETIRE-02 the suffix _qt is dropped — this script must be updated
# in lockstep when the rename lands.
set -eu
LIB=$1
PPDIR=$2

if [ ! -f "$LIB" ]; then
    echo "verify_pp_qt_freshness: libxvueqt.a not found at $LIB" >&2
    exit 1
fi

LIB_MTIME=$(stat -c '%Y' "$LIB")
EXIT=0
# Glob: post-RETIRE-02 the wildcard becomes pp/pp{init,mail,elas,flui,ther,nlse}.
for binary in "$PPDIR"/pp*_qt; do
    [ ! -f "$binary" ] && continue
    BIN_MTIME=$(stat -c '%Y' "$binary")
    if [ "$BIN_MTIME" -lt "$LIB_MTIME" ]; then
        echo "FAIL: $binary mtime ($BIN_MTIME) < libxvueqt.a mtime ($LIB_MTIME)" >&2
        EXIT=1
    else
        echo "OK: $binary mtime ($BIN_MTIME) >= libxvueqt.a mtime ($LIB_MTIME)"
    fi
done
exit $EXIT
```

Source: pattern cloned from `verify_no_imagemagick_in_qt` at
`xvue/qt/CMakeLists.txt:154-170` [VERIFIED via `Read xvue/qt/CMakeLists.txt:154-170`].

### Example 3: Harness `--out-dir` realpath sanitization (`bin/ab_sweep_phase8.sh`)

Before (lines 41–80, current):
```bash
OUT_DIR=".planning/phases/08-ab-validation-on-testa-subset/evidence"
# ... arg parse ...
mkdir -p "$OUT_DIR"
[ -z "$EVIDENCE_LOG" ] && EVIDENCE_LOG="${OUT_DIR}/sweep-log-${MODE}.md"
```

(Then later line 121 `pushd "$PROJDIR"` makes `OUT_DIR` resolve under PROJDIR if it was
relative.)

After (one-line insertion after argument parsing, around line 78):
```bash
OUT_DIR=$(realpath -m "$OUT_DIR")    # canonicalize BEFORE pushd into PROJDIR
mkdir -p "$OUT_DIR"
[ -z "$EVIDENCE_LOG" ] && EVIDENCE_LOG="${OUT_DIR}/sweep-log-${MODE}.md"
```

Note `-m` flag: do NOT require existence (the directory may not exist yet, gets created
by `mkdir -p` next line).

Source: `bin/ab_sweep_phase8.sh:42, 78, 121` [VERIFIED via `Read`].

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `xvue/xvuelc.c` Xlib backend | Qt 6 via `xvue/qt/` | 2026-04 → 2026-05 (Phases 0–7) | Single graphics backend; modern HiDPI + dark-mode support |
| LVIDEO `xwd → convert → .gif` pipeline | Phase 7 `XvueExport::saveGifTo` (ffmpeg-based) | Phase 7 (2026-05) | Cleaner; no temp `.xwd`/`.jpg` intermediates; ffmpeg is more portable than ImageMagick |
| Hand-rolled PostScript via `fprintf` in `xvuelc.c` | `PsEmitter::handleLasops` (verbatim port to C++) | Phase 7 EXPORT-04 | Byte-parity preserved; same `.eps` output format |
| Two parallel build paths (X11 + Qt) | Single Qt-only path | Phase 9 (this work) | Drops 32 `bin/cb*` X11 scripts; halves build complexity |

**Deprecated/outdated:**
- `bin/Makefile`, `bin/MakefileIBM`, `bin/MakefileMefisto` — legacy build files referencing
  HP/IBM systems and `/usr/lib/X11R5/`. Not used in any active build chain. Recommend
  delete in RETIRE-02 (or earlier).
- `bin/avecMOTIF` — comment-only file; X11/Motif relic.
- `bin/instal_2disk.{hp,src}` — HP-UX dual-disk installer; predates Linux port.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | The 12 LVIDEO-caller tracer files in flui/ + ther/ + util/ all use `IF (LVIDEO .NE. 0) THEN` as a guard, not as gating logic that affects non-LVIDEO behavior | §Common Pitfall 2, Pattern 2 | If a tracer's surrounding logic relies on a side effect from inside the LVIDEO block (e.g., `LVIDEO=0` reset that gates a later branch), excision breaks the file. Mitigation: Plan 9-03 task verifies file-by-file via per-file compile + targeted testa case |
| A2 | `bin/png2eps` + `bin/png2jpg` are dead code, not invoked from any active workflow | §Pitfall 8 | Maintainer might use them ad-hoc on shell command line; deletion would remove a long-standing utility. Mitigation: ask maintainer in Plan 9-04 |
| A3 | Plan 9-06 root cause for ppnlse_qt is event-loop starvation during long-running tttsupa2d-style iteration loops, NOT a true startup deadlock | §Architectural Map, §6 below | If actual root cause is different (e.g., a Qt6-specific QApplication-init race under offscreen+OMP), the recommended fix won't apply. **HIGH RISK** — needs hands-on debugging during Plan 9-06 task 1 (a diagnose-then-fix split is recommended) |
| A4 | The Debian package `qt6-image-formats-plugins` is the maintained successor to `libqt6imageformats6-plugins` and provides identical functionality | §Pitfall 9 | Could be that the maintainer wants the PDF-format plugin specifically (`qt6-image-formats-plugin-pdf` is a separate sub-package). Mitigation: Plan 9-04 task lists BOTH as candidates and asks maintainer |
| A5 | Cross-tag golden production in Plan 9-07 will succeed (the Phase 7 procedure has been validated by the verifier) | §Pattern 3, Plan 9-07 | If `gfortran scene01_driver.o $MEFISTO/xvue/xvuelc.o $MEFISTO/xvue/*.o $MEFISTO/util/*.o -L/usr/X11R6/lib -lX11 -lXt -o scene01_x11` fails to link (missing transitive deps), the golden cannot be produced. Phase 8 Plan 1 Task 2 (per 08-CHECKLIST.md) noted "scene01_driver.f link cannot be solved without invasive Fortran build-infra modification" — **this is a known partial blocker**. Recommendation: Plan 9-07 reproduces the link command exactly as documented in `xvue/qt/tests/golden/scene01_driver.f:36-43` and accepts that the fix may be additional `*.o` files |
| A6 | Setting MEFISTO_QT_WINDOW_SIZE in xvinitgraphique_ via `canvas->setFixedSize + win->adjustSize` will make the captured backing pixmap match the env-specified size | §Pattern 4 | Qt's window-chrome accounting may still steal a few pixels, or the setMinimum/Maximum approach may not work as expected under offscreen platform. Mitigation: Plan 9-05 includes a verify task that captures with the env set and confirms the resulting PNG dim |

**If this table is empty:** Not empty — six items need maintainer confirmation or empirical
verification during plan execution.

## Open Questions

1. **Is Phase 9 the right time to delete legacy `bin/cbg*`, `bin/cbpara*`, `bin/cbpoba`,
   `bin/cbonde`, `bin/cbbrezfort{2d,3d}`, `bin/Makefile*` scripts?**
   - What we know: They contain X11 linker lines (RETIRE-02 scope by name pattern). They
     are NOT invoked from `bin/cbl_tout` or `bin/cbl_tout_qt`. They appear to be
     legacy parallel/MPI/HP-UX build scripts.
   - What's unclear: Whether the maintainer still uses any of them ad-hoc.
   - Recommendation: Plan 9-02 task lists them for maintainer review; default action is
     delete (RETIRE-02 is meant to be aggressive).

2. **Should `bin/cbl_tout` post-rename be the merged Qt build, or should we delete the
   legacy `cbl_tout` and rename `cbl_tout_qt` to `cbl_tout`?**
   - What we know: D-12 leaves the choice to the planner. Both approaches are valid.
   - Recommendation: `git rm bin/cbl_tout && git mv bin/cbl_tout_qt bin/cbl_tout` in a single
     commit. Cleaner git history; semantically identical to a rewrite. Update all docs in
     same commit. (See §9 below for full reasoning.)

3. **Plan 9-06 (ppnlse_qt deadlock): is the symptom "10 log lines, no NLSER banner" really
   a startup deadlock, or is it the long-running TIME=20 loop just appearing hung?**
   - What we know: 08-smoke-probes.md documents nlsecu canonical TIME=20 (2000 steps)
     running ~hour-scale on the dev hardware. 60s timeout = expected fail. Phase 8 Plan 1
     SUMMARY says "deadlock at startup, no NLSER banner reached even at 240s" — but ppnlse.f
     code shows the banner is printed AT line 117/120/121 BEFORE any solver code. If 10 log
     lines come out with no banner, the banner sequence is being skipped — that suggests a
     Fortran-side branching issue (INTERA mismatch?) rather than a Qt-startup deadlock.
   - What's unclear: Whether the actual hang is in Qt event-loop init or in a loop earlier
     than the banner.
   - Recommendation: Plan 9-06 starts with a DIAGNOSE task (capture stack trace via
     `gdb -p $(pgrep ppnlse_qt)` after 30s). If stack is in `tttsupa2d.f` iteration, root
     cause is event-loop starvation (fix: pump events). If stack is in QApplication init,
     root cause is true deadlock (fix: serialise call_once differently).

4. **Should Plan 9-08's `verify_pp_qt_freshness` be a CMake ALL target (build-time) or a
   ctest test target (test-time)?**
   - What we know: Phase 8 D-09 says "ALL target". `verify_no_imagemagick_in_qt` is also
     ALL. But the `xvueqt` library is built BEFORE `pp/pp*_qt` link (since
     `bin/cbl_tout_qt` runs `cmake --build` first then the bash `cb*_qt` scripts), so an
     ALL target would always fail when libxvueqt.a is newer than the (still un-linked) pp/pp*.
   - What's unclear: Whether the freshness check needs to be deferred to AFTER the `pp/pp*`
     link step.
   - Recommendation: Make it a custom target invoked at the END of `bin/cbl_tout` (post-rename),
     after the per-module bash scripts have produced `pp/pp*`. Run via
     `cmake --build xvue/qt/build --target verify_pp_qt_freshness` as the very last step
     before the final `echo "All executables created"` line. This keeps the spirit of D-09
     (auto-fail on stale binaries) without breaking the build-order chicken-and-egg.

5. **For Plan 9-05 matched-dim: does setting `setMinimumSize + setMaximumSize` propagate
   to the backing-pixmap allocation in `XvueCanvas::resizeEvent`?**
   - What we know: `XvueCanvas::resizeEvent` (xvue_qt_canvas.cpp:126) reallocates
     `state_->backing_` from `size() * dpr`, where `size()` is the canvas widget's logical
     size. Min/Max constraints DO affect canvas size after `win->adjustSize()`.
   - What's unclear: Order of operations — is `xvinitgraphique_` window-show + min/max set
     guaranteed to fire `resizeEvent` exactly once with the target dim?
   - Recommendation: Plan 9-05 verify task captures the qt-1x PNG, confirms `identify` reports
     1280x800 (the env-specified dim), iterates if not.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| gfortran | All Fortran compiles incl. RETIRE-NN verify steps | ✓ | 14.3.0 (Debian) | — |
| gcc | Indirect via cmake; Fortran wrappers don't compile C anymore post-Phase-9 | ✓ | (system) | — |
| Qt 6 (qt6-base-dev) | xvue/qt/ build | ✓ | 6.10.2+dfsg-7 | — |
| ffmpeg | Phase 7 GIF export (already shipped); RETIRE-03 retains capability | ✓ | 8.1-3+b1 | — |
| ImageMagick | Phase 8 A/B harness (`bin/ab_compare_pair.sh`); STAYS in tooling per Phase 7 EXPORT-06 / Phase 8 D-03 | ✓ | 7.1.2-18 (Q16) | — |
| Xvfb (`xvfb-run`) | Plan 9-07 cross-tag golden bootstrap | ✓ | (apt-supplied) | Real X11 display in interactive session (slower) |
| cmake | xvue/qt/ build | ✓ | 3.31.6 | — |
| `realpath` | Plan 9-08 harness fix | ✓ | (coreutils) | `readlink -f` on systems without realpath; bash-portable |
| Git tag `v1.0-pre-retire` | Plan 9-07 cross-tag worktree | ✗ (current `git tag --list` returns empty in this dev tree) | — | **Maintainer must create the tag BEFORE Plan 9-07 starts.** This is the literal D-09 condition |
| libgfortran5 (15.2.0-9 hold) | Pre-existing UB-trip mitigation per Phase 03-04 reopen 2026-04-14 | ✓ | 15.2.0-9 (`apt-mark hold` per STATE.md) | None — version pin is mandatory until Fortran UB sites are fixed (out of Phase 9 scope) |

**Missing dependencies with no fallback:**
- `git tag v1.0-pre-retire` — **gating dependency for Plan 9-07**. Per CONTEXT.md D-08
  the tag is created BEFORE Plan 9-01 starts (rollback safety). Plan 9-07 simply
  consumes it. If the maintainer forgets to tag before /gsd-execute-phase fires, Plan 9-01
  blocks per D-09. The plan-checker should verify this.

**Missing dependencies with fallback:**
- None.

## Validation Architecture

`workflow.nyquist_validation` is `true` in `.planning/config.json` [VERIFIED 2026-05-05].

### Test Framework

| Property | Value |
|----------|-------|
| Framework | Qt Test (QTest framework, Qt6 6.10.2) for unit; bash + AUTOEXIT smoke + ImageMagick AE compare for integration |
| Config file | `xvue/qt/tests/CMakeLists.txt` (Qt Test); `bin/ab_sweep_phase8.sh` (integration sweep harness) |
| Quick run command | `cd xvue/qt/build && ctest -R '^xvue_qt_(postscript\|export\|menu_tests)$'` |
| Full suite command | `bin/cbl_tout && cd xvue/qt/build && xvfb-run --auto-servernum ctest && bin/ab_sweep_phase8.sh --mode qt-1x --cases pan2d,nafems_le1,cavity2d,heat1d,nlsecu` |

### Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| RETIRE-01 | xvuelc.c + ccxvue deleted; full Qt-only build green | smoke | `bin/cbl_tout && ls -la pp/pp*_qt` | ✅ (build script + binaries exist as gates) |
| RETIRE-01 verify | No xvuelc symbol references remain in Fortran link path | unit | `! nm xvue/qt/build/libxvueqt.a \| grep ' U xvtest_x11_only_symbol'` (negative-grep on hypothetical xvuelc-only entry) | ✅ (verify_abi.sh in xvue/qt/cmake/) |
| RETIRE-02 | No `-lX11` / `/usr/X11R6` in any active build script | unit | `! grep -r 'lX11\|X11R6' bin/cb* bin/cbl_tout` | ❌ (Wave 0 — add `bin/test_no_x11_in_build.sh` modeled on `bin/test_no_imagemagick_in_qt.sh`) |
| RETIRE-03 | LVIDEO + `convert` shell-outs gone from active source paths | unit | `! grep -rn 'CALL VIDEO' flui/ ther/ util/ xvue/` AND `! grep -rn 'CALL SYSTEM.*convert' xvue/ flui/ ther/ util/` | ❌ (Wave 0 — add a new grep-gate script `bin/test_no_lvideo.sh`) |
| RETIRE-04 | Documentation lists only Qt deps | smoke / human verify | grep `'libX11-dev\|imagemagick'` returns 0 hits across {README, LISEZMOI, bin/README, bin/LISEZMOI, CLAUDE.md} | ✅ (grep is the test) |
| Carry 9-05 | Qt 1x captures land at the env-specified WxH | integration | `MEFISTO_QT_WINDOW_SIZE=1280x800 bin/ab_sweep_phase8.sh --mode qt-1x --cases pan2d --out-dir /tmp/p9-05 && identify -format '%wx%h' /tmp/p9-05/pan2d-qt-1x.png` returns `1280x800` | ✅ (harness already exists) |
| Carry 9-06 | nlsecu Qt 1x sweep produces a non-truncated capture (NLSER banner reached, frame from solver step) | integration | `bin/ab_sweep_phase8.sh --mode qt-1x --cases nlsecu --out-dir /tmp/p9-06` exits within 600s timeout AND output PNG > 1KB AND contains an NLSER-frame marker (manual inspection or future grep on EPS dump) | partial — needs custom timeout extension and manual inspection step |
| Carry 9-07 | 3 Phase-7-deferred goldens land; ctest QSKIPs flip to PASS | unit | `ls xvue/qt/tests/golden/{scene01.eps,wave_legacy.gif,cavity2d_legacy.gif}` AND `ctest -R '^xvue_qt_(postscript\|export)_tests$' --output-on-failure` exits 0 with 0 SKIP | ✅ (test slots exist; QSKIP→PASS auto-flip) |
| Carry 9-08 (`--out-dir` fix) | Relative `--out-dir` resolves correctly | unit | `cd /tmp && bin/ab_sweep_phase8.sh --mode qt-1x --cases pan2d --out-dir ./out --smoke-only && ls /tmp/out/pan2d-qt-1x.png` | ✅ |
| Carry 9-08 (CMake guard) | `verify_pp_qt_freshness` fails when libxvueqt.a is newer than any pp/pp*_qt | unit | `touch xvue/qt/build/libxvueqt.a && cmake --build xvue/qt/build --target verify_pp_qt_freshness` should exit non-zero | ❌ (Wave 0 — add the target + script) |

### Sampling Rate

- **Per task commit:** `bin/cbl_tout` (incremental rebuild — full Qt-only build with the
  CMake guards), takes ~30 seconds. CMake ALL-target gates run automatically.
- **Per wave merge:** Above + `xvfb-run ctest` (Qt unit tests), takes ~3 minutes.
- **Phase gate:** Above + full Phase 8 harness sweep (`bin/ab_sweep_phase8.sh --mode qt-1x`
  on all 5 cases), takes ~10 minutes plus matched-dim verification.

### Wave 0 Gaps

- [ ] `bin/test_no_x11_in_build.sh` — RETIRE-02 grep gate, mirrors `bin/test_no_imagemagick_in_qt.sh`
- [ ] `bin/test_no_lvideo.sh` — RETIRE-03 grep gate (covers `CALL VIDEO*` and `CALL SYSTEM('convert')`)
- [ ] `xvue/qt/cmake/verify_pp_qt_freshness.sh` + corresponding `add_custom_target` in CMakeLists.txt — Carry 9-08
- [ ] Pattern-extension target in `xvue/qt/CMakeLists.txt` to call the two new test scripts (so they fail the build, not just stand alone)
- [ ] Update `bin/test_no_imagemagick_in_qt.sh` to remove its scope comment about LVIDEO being legitimate (it no longer is post-Phase-9)

## Security Domain

`security_enforcement` is not explicitly set in `.planning/config.json`. Per the research
template's "absent = enabled" rule, applying ASVS scoping below.

### Applicable ASVS Categories

| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | no | MEFISTO is a desktop scientific tool, no auth |
| V3 Session Management | no | No sessions |
| V4 Access Control | no | Filesystem-level; user runs as themselves |
| V5 Input Validation | partial | Phase 9 carry-forward 9-05 reads `MEFISTO_QT_WINDOW_SIZE` env — must `sscanf("%dx%d")` and bound-check (`> 0 && < 8192`) to avoid 0-byte alloc or negative dim |
| V6 Cryptography | no | No crypto |

### Known Threat Patterns for {Fortran + Qt + bash retirement}

| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Env-var injection via MEFISTO_QT_WINDOW_SIZE causing integer overflow on `setFixedSize` | Tampering | bounded `sscanf` parse; clamp to [200, 4096] each dim; ignore malformed input silently |
| Env-var injection via MEFISTO_BATCH_X11 already exists and is HONORED — Phase 9 doesn't add new attack surface | — | already battle-tested through Phase 8 |
| Cross-tag worktree leaks artifacts to main on git operation mistake | Tampering (build artifacts) | Pattern 3's `git status` clean check pre-commit |
| `realpath -m` accepts non-existent paths and could canonicalize a path under `/tmp/` to inadvertently overwrite system files | Privilege | The `OUT_DIR` is always under user-owned scope (`$MEFISTOX/`-relative or `.planning/`); no risk of writing system locations. Default value is `.planning/...` which is in-repo |
| Fortran `CALL SYSTEM('convert ...')` shell-injection via NMFVIDEO/NOMFGIF buffer | Tampering | Phase 9 deletes these CALL SYSTEM sites entirely (RETIRE-03) — the threat goes away |

## Project Constraints (from CLAUDE.md)

Extracted from `/home/mefisto/git/mefisto/CLAUDE.md` (project-level instructions):

- **Build never breaks (CRITICAL):** "After every change, verify that the affected module
  still compiles with its `cb*` script. Before committing, the full build (`bin/cbl_tout`)
  must succeed." Phase 9 must NOT leave the tree in an unbuildable state at any commit
  boundary. Each RETIRE-NN plan's last task is a full-build verify.
- **testa cases must pass:** "Small tests in `testa/` or `testf/` must continue to pass
  after every change." 5 BUILD-10 cases are the canonical baseline; verify after every plan.
- **Large/long tests delegated to maintainer:** "For large or long-running tests, ask the
  user to run them before declaring a change complete." Applies to nlsecu (Plan 9-06) full
  TIME=20 run; truncated capture acceptable for autonomous runs.
- **Ask before installing system packages:** "If a C include or an Ubuntu system package
  is missing to compile, ask the user to install it." Phase 9 doesn't install anything new
  (deletion phase) — but Plan 9-04 doc update may shift dep list; if the maintainer's machine
  doesn't have `qt6-image-formats-plugins` installed yet, ask.
- **Commit after every logical step:** "Commit after every logical step where rolling back
  would be useful (working build, passing test, completed sub-fix…)." 8-plan structure with
  multiple commits per plan supports this.
- **Never force-push, never bypass hooks:** standard git discipline.
- **Programming norms:** `doc/normes.ps` — naming + comment + fixed-form Fortran column
  layout. Plan 9-03 LVIDEO surgery in flui/ + ther/ + util/ tracers must respect column 7
  start.
- **Active project goal:** "Bug fixes: identify and correct existing bugs without altering
  the overall behaviour." Carry-forward 9-06 ppnlse_qt deadlock fits this; carry-forward
  9-05 matched-dim recapture is closer to "feature parity" than "bug" — but acceptable per
  Phase 8 override #1 explicitly deferring to Phase 9.
- **Qt migration end state:** "the X11/Motif graphical layer (`xvue/`) will eventually be
  replaced by Qt." Phase 9 IS that replacement reaching its end state.

The planner MUST verify each plan's task sequencing respects "build never breaks" — i.e.,
every task either leaves the tree green or includes the make-green step in the same commit.

## Sources

### Primary (HIGH confidence)

- `/home/mefisto/git/mefisto/CLAUDE.md` — project guardrails [VERIFIED via Read]
- `/home/mefisto/git/mefisto/.planning/REQUIREMENTS.md §RETIRE-01..04` — 4 requirements [VERIFIED via Read]
- `/home/mefisto/git/mefisto/.planning/phases/09-retire-x11-backend/09-CONTEXT.md` — 12 locked decisions [VERIFIED via Read]
- `/home/mefisto/git/mefisto/.planning/phases/07-image-gif-and-postscript-export/VERIFICATION.md §9` — 3 deferred goldens bootstrap procedure [VERIFIED via Read]
- `/home/mefisto/git/mefisto/.planning/phases/08-ab-validation-on-testa-subset/08-CHECKLIST.md` — 5 override decisions; #1 + #2 + #5 carry into Phase 9 [VERIFIED via Read]
- `/home/mefisto/git/mefisto/.planning/phases/08-ab-validation-on-testa-subset/08-VALIDATION.md` — VALID-01..07 coverage [VERIFIED via Read]
- `/home/mefisto/git/mefisto/.planning/phases/08-ab-validation-on-testa-subset/08-CONTEXT.md` D-09 — verify_pp_qt_freshness deferred to Phase 9 [VERIFIED via Read]
- `/home/mefisto/git/mefisto/xvue/xvuelc.c:38, 1640-1680` — legacy backend size and BATCH_X11 path [VERIFIED via grep + sed]
- `/home/mefisto/git/mefisto/xvue/qt/src/xvue_qt_api.cpp:336-487, 822-910` — Qt xvinitgraphique_ + xvfermer_ + xvinfo_ window-size logic [VERIFIED via Read]
- `/home/mefisto/git/mefisto/xvue/qt/src/xvue_qt_canvas.cpp:126-209` — backing pixmap resize logic [VERIFIED via Read]
- `/home/mefisto/git/mefisto/bin/ab_sweep_phase8.sh:42-218` — `--out-dir` bug location [VERIFIED via Read]
- `/home/mefisto/git/mefisto/bin/cbl_tout` and `bin/cbl_tout_qt` — build entry points [VERIFIED via Read]
- `/home/mefisto/git/mefisto/incl/trvari.inc:9` — LVIDEO COMMON declaration [VERIFIED via Read]
- `/home/mefisto/git/mefisto/xvue/video1.f, videofin.f, videonm.f` — full file contents [VERIFIED via Read]
- `/home/mefisto/git/mefisto/prpr/ppnlse.f:90-145, prpr/ppmail.f:100-170` — BATCH_X11 + INTERA dispatch [VERIFIED via Read]
- `/home/mefisto/git/mefisto/xvue/qt/CMakeLists.txt:115-170` — existing CMake ALL-target patterns [VERIFIED via Read]
- `/home/mefisto/git/mefisto/xvue/qt/tests/golden/scene01_driver.f:33-46` — bootstrap procedure [VERIFIED via Read]

### Secondary (MEDIUM confidence)

- Empirical dep verification on dev host: `apt-cache show qt6-base-dev`,
  `pkg-config --modversion Qt6Widgets`, `ffmpeg -version`, `cmake --version`,
  `gfortran --version`, `convert --version` — all returned current versions on the dev
  host where this research ran [VERIFIED 2026-05-05]
- `apt-cache search '^qt6-image-formats'` returns `qt6-image-formats-plugins`, contradicting
  REQUIREMENTS.md `libqt6imageformats6-plugins` claim [VERIFIED 2026-05-05]

### Tertiary (LOW confidence)

- Assumption A3 (ppnlse_qt deadlock root cause = event-loop starvation, NOT startup
  deadlock): inferred from comparing 08-smoke-probes.md log output with the Fortran code
  flow in prpr/ppnlse.f. Not empirically debugged. **Plan 9-06 must include a diagnose
  task before implementing a fix.**

## Metadata

**Confidence breakdown:**
- RETIRE-01 (xvuelc + ccxvue): **HIGH** — single C source, single bash wrapper, zero Fortran-side ghost references
- RETIRE-02 (libX11 linker lines + 32 bin/cb* scripts): **HIGH** — every reference greppable; deletion targets clear
- RETIRE-03 (ImageMagick + LVIDEO): **MEDIUM-HIGH** — 12 tracer files need surgery, not bulk delete (D-06); Plan 9-03 needs careful per-file review. Two new findings (`png2eps`/`png2jpg`) require maintainer decision
- RETIRE-04 (docs): **HIGH** — 5 file locations enumerated; pkg name correction empirically verified
- Carry-forward 9-05 (matched-dim recapture): **MEDIUM** — env-knob pattern is straightforward; Plan must confirm `setMinimum/Maximum + adjustSize` produces the expected backing-pixmap size empirically
- Carry-forward 9-06 (ppnlse_qt deadlock): **LOW-MEDIUM** — root cause not empirically confirmed; A3 is the highest-risk assumption in this research. Plan 9-06 needs a diagnose task before fix
- Carry-forward 9-07 (3 Phase-7 goldens): **MEDIUM-HIGH** — procedure documented in Phase 7 VERIFICATION.md §9 but Phase 8 Plan 1 Task 2 noted "scene01_driver.f link cannot be solved without invasive Fortran build-infra modification". The cross-tag worktree approach should sidestep that, but verify pre-execution
- Carry-forward 9-08 (harness fix + CMake guard): **HIGH** — both are mechanical; verify_pp_qt_freshness pattern proven; realpath fix is one line

**Research date:** 2026-05-05

**Valid until:** 2026-06-05 (30 days, stable phase). Re-verify dependency versions before Plan 9-04 if anything in `apt-cache` shifts.
