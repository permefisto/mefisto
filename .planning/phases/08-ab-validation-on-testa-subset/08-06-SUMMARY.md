---
phase: 08-ab-validation-on-testa-subset
plan: 06
subsystem: validation
tags: [phase8, valid-05, valid-06, spot-check, colorbar, font-metric, pitfall-6, pitfall-7, hexahedron, gap-documented]

requires:
  - phase: 08-ab-validation-on-testa-subset (Plan 02)
    provides: 5 X11 baseline PNG captures at 1280x800 (LEFT operands for compares)
  - phase: 08-ab-validation-on-testa-subset (Plan 03)
    provides: 5 Qt 1x PNG captures at 760x442 (RIGHT operands for compares) + per-cell pre-shaped CHECKLIST scaffolds
  - phase: 08-ab-validation-on-testa-subset (Wave-2 merge)
    provides: evidence/wave-merge-ae.md (post-merge AE re-run for all 10 Wave-2 cells)
provides:
  - .planning/phases/08-ab-validation-on-testa-subset/evidence/colorbar-nafems_le1.md (Pitfall 6 spot check, full-frame, fuzz=8 recommended)
  - .planning/phases/08-ab-validation-on-testa-subset/evidence/colorbar-heat1d.md (Pitfall 6 spot check, full-frame, fuzz=8 recommended)
  - .planning/phases/08-ab-validation-on-testa-subset/evidence/colorbar-cavity2d.md (Pitfall 6 spot check, full-frame + cropped right-edge legend, fuzz=10 recommended)
  - .planning/phases/08-ab-validation-on-testa-subset/evidence/font-pan2d.md (Pitfall 7 AA-fringe analysis, AE budget=3000 confounded by dim-mismatch)
  - .planning/phases/08-ab-validation-on-testa-subset/evidence/font-hexahedron.md (GAP-branch report; all 5 literal section headings preserved per WARNING-3 iter2)
  - .planning/phases/08-ab-validation-on-testa-subset/evidence/hexahedron-probe.md (alphabetical-order probe trail across {ln, ob, pt, sf, vl} per WARNING-5 iter2; Outcome=GAP)
  - .planning/phases/08-ab-validation-on-testa-subset/evidence/spot-check-summary.md (5-row aggregated table for Plan 07 ingest)
  - .planning/phases/08-ab-validation-on-testa-subset/08-VALIDATION.md (## Open escalations entry for hexahedron — exceptional Plan 06 modification per step 12c)
affects: [08-07-checklist-finalize]

tech-stack:
  added:
    - 7 new evidence markdown files (3 colorbar reports + 2 font reports + 1 probe report + 1 summary; total ~32 KB)
    - 1 new validation tracking file (08-VALIDATION.md with open-escalations section)
  patterns:
    - "Region-of-interest crop pattern: convert -crop WxH+X+Y +repage source.png region.png; bin/ab_compare_pair.sh region-A region-B diff fuzz — applied empirically only when an actual gradient region was visible in the X11 capture (cavity2d only)."
    - "Probe-then-document pattern for non-baseline cases (BLOCKER #1 / MAJOR-9 mitigation): probe alphabetical-order candidates with workspace re-seed between each (rm -rf workspace && cp + INITIER); 60-s per-candidate timeout; 5-min total budget; meaningful-content acceptance test (not just file-size > 0); document GAP outcome with all 5 literal section headings preserved per WARNING-3 iter2."
    - "Resample-dim-mismatch as dominant-AE source pattern: every Wave-2 cell + every Plan 06 Pitfall 6/7 cell shows ~14-53 % AE driven primarily by the 760x442 → 1280x800 nearest-neighbour resample mandated by ab_compare_pair.sh dim-guard; the planner-prescribed Pitfall 6/7 overrides + budgets cannot be applied cleanly to the resample-confounded numbers. Plan 7 + Phase 9 carry the matched-dim recapture as a deferred idea."

key-files:
  created:
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/colorbar-nafems_le1.md
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/colorbar-heat1d.md
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/colorbar-cavity2d.md
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/font-pan2d.md
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/font-hexahedron.md
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/hexahedron-probe.md
    - .planning/phases/08-ab-validation-on-testa-subset/evidence/spot-check-summary.md
    - .planning/phases/08-ab-validation-on-testa-subset/08-VALIDATION.md
  modified:
    - (none — Phase 8 verify-only contract honored; no source-tree mutation; xvue/xvuelc.c + xvue/qt/src/ byte-identical)

key-decisions:
  - "Full-frame compare for nafems_le1 + heat1d (no per-region crop): empirically determined that nafems_le1 has only 7 unique colors (discrete contour bands, no gradient bar) and heat1d has only 8 unique colors (1D thermal-curves plot, no temperature ramp). The plan's per-case crop coordinates (right-edge ~80 px wide) cannot be cleanly applied because no gradient colorbar is visible in either X11 capture. Recorded as 'full-frame — no smooth-gradient bar visible' in each report's Region section. Pitfall 6 hypothesis (gradient drift) is empirically DISCONFIRMED for both cases."
  - "cavity2d split measurement: full-frame (40.14 % AE) AND cropped right-edge legend at x=1100-1260 (160x800 region, 21.94 % AE). Cropped-region AE LOWER than full-frame disconfirms Pitfall 6 as principal failure mode even on the only case where a real visible colorbar exists. The plan's Pitfall 6 fuzz=10 override remains the recommended cell annotation for forward compatibility (matched-dim recapture would benefit from it)."
  - "hexahedron probe outcome = GAP: all 5 alphabetical candidates {ln, ob, pt, sf, vl} probed; 4 crashed at unknown-name lookup with 2045-byte all-peach background captures (xvfermer_ fired at crash time before any drawing primitive painted the mesh); 1 (`pt`) hung at the 60-s timeout (no rendering trigger; expects interactive menu pick). The 5 files are interlocking building blocks (pt → ln → sf → vl → ob), NOT independent single-file batch entry points. Legacy MAILLER protocol's single-batch invocation pattern does not apply. Plan 06 invokes the GAP branch (step 12) per WARNING-3 iter2 (preserve all 5 literal section headings) and step 12c (file ## Open escalations entry in 08-VALIDATION.md — exceptional Plan 06 modification documented in the file's owner-note)."
  - "'Meaningful capture' acceptance test: probe step 9's 'capture_size > 0 AND exit ∈ {0, 124-with-capture}' criterion, taken literally, would mark all 4 crash candidates (ln/ob/sf/vl, all 2045 B) as CAPTURED — even though they show only the peach background. The Plan 02 'Color-content sanity check' rule (a capture with only 1-2 grayscale entries indicates xvfermer_ fired before any data was painted) is the corrective interpretation: all 5 candidates fail meaningful-content. Hexahedron-probe.md documents this interpretation explicitly under §Outcome to forestall a future agent thinking the 2045-B captures qualify."
  - "Resample-dim-mismatch as dominant signal across all 4 CHECK cells: nafems_le1 (40.31 %), heat1d (13.99 %), cavity2d (40.14 %), pan2d (52.81 %). Same root cause as the wave-merge-ae.md consolidated 10-cell finding: the 760x442 Qt offscreen backing pixmap → 1280x800 X11 baseline resample produces ~1.68 px nearest-neighbour edge offsets that count as drift. The Pitfall 6/7 fuzz overrides + AE budgets cannot be applied cleanly to the resample-confounded numbers; Plan 7 cell annotations pair the override values with `dim-mismatch=resample-to-1280x800` notes. A matched-dim Qt recapture (1280x800 in-process backing pixmap if MEFISTO_QT_WINDOW_SIZE or equivalent is supported) is recorded as a Phase 9 deferred-idea."

patterns-established:
  - "Probe-first-document-on-failure pattern (BLOCKER #1 / MAJOR-9 mitigation): when a non-BUILD-10 spot-check target lacks an obvious headless driver, probe alphabetical-order candidates with strict per-candidate + total budgets; document GAP outcome with all 5 literal section headings preserved (per WARNING-3 iter2) so verify gates count cleanly; file an escalation entry in 08-VALIDATION.md with maintainer-decision options (a) CONTEXT.md amendment naming canonical driver / (b) waiver-by-rationale."
  - "Pitfall-6/7-disconfirmed-via-empirical-region-survey: instead of trusting the planner-prescribed crop coordinates, run a unique-color count strip-by-strip survey to verify a gradient region actually exists. Three of three colorbar-cases (nafems_le1, heat1d, cavity2d) found Pitfall 6 disconfirmed as the principal failure mode; the resample-dim-mismatch dominates instead."
  - "Co-existing 'reserved' override + 'unapplicable' status pattern: reports record fuzz_override values AND AE budgets even when those values cannot be applied cleanly to the current resample-confounded measurement. Forward compatibility: matched-dim recaptures will be able to apply them directly without re-running the spot-check analysis."

requirements-completed: [VALID-05, VALID-06]

duration: ~10min
completed: 2026-05-05
---

# Phase 8 Plan 06: VALID-05 + VALID-06 Spot-Check Sweep Summary

**5 spot-check reports filed (3 Pitfall-6 colorbar + 2 Pitfall-7 font) plus 1 hexahedron probe report establishing GAP outcome; aggregated 5-row table for Plan 07 ingest; verdict roll-up = 0 PASS / 4 CHECK / 0 FAIL / 1 GAP-DOCUMENTED.**

## Performance

- **Duration:** ~10 min wall-clock
- **Started:** 2026-05-05T15:21:00Z (approx)
- **Completed:** 2026-05-05T15:30:00Z (approx)
- **Tasks:** 3/3 (Task 1 colorbar, Task 2 font + probe, Task 3 summary)
- **Files created:** 8 (7 evidence reports + 1 validation file)
- **Files modified:** 0 (Phase 8 verify-only contract honored; xvuelc.c + xvue/qt/src/ byte-identical)

## Accomplishments

- **VALID-05 (3/3 colorbar cases):** filed colorbar-nafems_le1.md /
  colorbar-heat1d.md / colorbar-cavity2d.md, each with Region + AE
  measurements + Hypothesis + Recommended fuzz_override_pct + Verdict
  sections. Empirical finding: Pitfall 6 hypothesis (gradient color bar
  tripping 5 % fuzz) is **disconfirmed for all 3 cases** — nafems_le1 + heat1d
  have no smooth-gradient bar (7-8 unique colors total); cavity2d has a real
  right-edge legend but its cropped-region AE (21.94 %) is LOWER than
  full-frame AE (40.14 %), disconfirming Pitfall 6 as principal failure
  mode. Recommended fuzz overrides (8 / 8 / 10) retained for forward
  compatibility with future matched-dim recaptures.

- **VALID-06 Part A (pan2d Pitfall 7):** filed font-pan2d.md with empirical
  confirmation of Pitfall 7 (3545 unique colors in Qt-1x vs 9 in X11 = 394×
  AA explosion). The pan2d cell's 540 804 px AE (52.81 %) is dominated by
  the resample-dim-mismatch artefact common to all 5 Wave-2 cells; Pitfall 7's
  3000-px budget is confounded by the resample.

- **VALID-06 Part B (hexahedron probe):** filed hexahedron-probe.md with
  Outcome=GAP after probing all 5 alphabetical candidates (per WARNING-5
  iter2): `ln`/`ob`/`sf`/`vl` crashed at unknown-name lookup (CRASH); `pt`
  hung at 60-s timeout. The 5 files are interlocking building blocks
  (pt → ln → sf → vl → ob), NOT independent batch entry points. Filed
  font-hexahedron.md with all 5 LITERAL section headings preserved per
  WARNING-3 iter2 (Captures / Diff statistics / AE budget / Hypothesis /
  Verdict — bodies = "N/A — see hexahedron-probe.md"). Filed
  08-VALIDATION.md `## Open escalations` entry per step 12c (the only
  Phase 8 case where Plan 06 modifies this file; the file's owner-note
  documents the exception).

- **Aggregated spot-check-summary.md** for Plan 07 ingest: 5-row table +
  VALID-05 / VALID-06 roll-ups + per-main-case CHECKLIST.md cell
  annotations + separate hexahedron spot-check-row guidance + Outcome
  line. Plan 07 ingests this verbatim into 08-CHECKLIST.md.

- **Phase 8 invariants verified:** `git diff -- xvue/xvuelc.c` empty;
  `git diff -- xvue/qt/src/` empty. No source-tree mutation. Per-task
  atomic commits (3) honoring the parallel-execution contract (no
  STATE.md / ROADMAP.md modification — Wave 3 solo execution).

## Task Commits

Each task was committed atomically:

1. **Task 1: VALID-05 colorbar spot-check reports for 3 cases** — `f507f5e` (feat)
   - 3 reports (colorbar-nafems_le1.md / colorbar-heat1d.md / colorbar-cavity2d.md)
   - 5 sections each (Region / AE measurements / Hypothesis / Recommended fuzz_override_pct / Verdict)
   - 358 insertions

2. **Task 2: VALID-06 font-metric reports + hexahedron probe (GAP)** — `2f942dc` (feat)
   - font-pan2d.md (5 sections, Pitfall 7 confirmed qualitatively)
   - font-hexahedron.md (5 LITERAL sections, GAP-branch per WARNING-3 iter2)
   - hexahedron-probe.md (Outcome=GAP, 5 candidates probed alphabetically per WARNING-5 iter2)
   - 08-VALIDATION.md (## Open escalations entry — exceptional Plan 06 modification per step 12c)

3. **Task 3: Aggregated spot-check verdicts for Plan 07 ingest** — `6e8db45` (feat)
   - spot-check-summary.md (5 sections, 5 aggregated rows)
   - VALID-05 + VALID-06 roll-ups + per-main-case Notes annotations + separate hexahedron spot-check-row guidance

_Plan metadata commit will be added after this SUMMARY is staged (per
docs(08-06): finalize plan summary pattern; Wave 3 solo execution does
NOT modify STATE.md or ROADMAP.md per the parallel-execution contract)._

## Files Created/Modified

### Created (8)

| Path                                                                                              | Lines | SHA-256 (first 16) | Purpose                                                                            |
| ------------------------------------------------------------------------------------------------- | ----- | ------------------ | ---------------------------------------------------------------------------------- |
| `.planning/phases/08-ab-validation-on-testa-subset/evidence/colorbar-nafems_le1.md`               | 118   | ea0fe79e96e51e57   | Pitfall 6 spot check, fuzz=8 recommended                                           |
| `.planning/phases/08-ab-validation-on-testa-subset/evidence/colorbar-heat1d.md`                   | 112   | 1a7ef0f103eddcba   | Pitfall 6 spot check, fuzz=8 recommended                                           |
| `.planning/phases/08-ab-validation-on-testa-subset/evidence/colorbar-cavity2d.md`                 | 128   | ed618d467c889d2e   | Pitfall 6 spot check, fuzz=10 recommended; full-frame + cropped right-edge legend  |
| `.planning/phases/08-ab-validation-on-testa-subset/evidence/font-pan2d.md`                        | 154   | c17019f214971dd1   | Pitfall 7 AA-fringe spot check; Pitfall 7 confirmed qualitatively                  |
| `.planning/phases/08-ab-validation-on-testa-subset/evidence/font-hexahedron.md`                   | 59    | a0fc9f6be574e8c8   | GAP-branch font report; 5 literal section headings preserved (WARNING-3 iter2)     |
| `.planning/phases/08-ab-validation-on-testa-subset/evidence/hexahedron-probe.md`                  | 241   | da83cc496130d574   | Alphabetical-order probe trail (WARNING-5 iter2); Outcome=GAP                      |
| `.planning/phases/08-ab-validation-on-testa-subset/evidence/spot-check-summary.md`                | 128   | 30aaa76517f74c45   | 5-row aggregated table + roll-ups + Plan 07 hand-off                               |
| `.planning/phases/08-ab-validation-on-testa-subset/08-VALIDATION.md`                              | 80    | 7fe61eb59d2179e7   | ## Open escalations section (exceptional Plan 06 file modification per step 12c)   |

### Modified (0)

Phase 8 verify-only contract honored — no source files in the project's
working tree (`xvue/`, `prpr/`, `util/`, `mail/`, `elas/`, `flui/`, `ther/`,
`nlse/`, `bin/`) edited:

- `xvue/xvuelc.c` byte-identical (`git diff -- xvue/xvuelc.c` → empty).
- `xvue/qt/src/` no working-tree changes (`git diff -- xvue/qt/src/` → empty).

## Decisions Made

See `key-decisions:` in frontmatter (5 decisions). The two most load-bearing:

1. **Full-frame compare for nafems_le1 + heat1d (Pitfall 6 disconfirmed)**
   — empirical strip-by-strip unique-color survey showed neither case has a
   smooth-gradient colorbar in the X11 capture. The plan's region-crop
   instructions are inapplicable; reports record this finding explicitly
   in each Region section.

2. **hexahedron probe outcome = GAP (BLOCKER #1 / MAJOR-9 mitigation)**
   — all 5 alphabetical candidates fail the meaningful-capture acceptance
   test (4 crash, 1 hangs); the 5 files are interlocking mesh-definition
   building blocks, not single-file batch drivers. GAP branch invoked per
   step 12; font-hexahedron.md preserves all 5 literal section headings
   per WARNING-3 iter2; 08-VALIDATION.md `## Open escalations` filed per
   step 12c.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 — Blocking issue] X11 + Qt mesher binaries (`pp/pp{mail,init,mail_qt}`) absent from worktree pp/**

- **Found during:** Task 2 step 0 (hexahedron probe pre-flight) — after the
  initial worktree branch check and recovery via `git fetch . main &&
  git reset --hard FETCH_HEAD`.
- **Issue:** Worktree's `pp/` tree contained only `pp/pxyz` (the only
  tracked binary; rest are gitignored per `.gitignore` line `pp/pp*`).
  Plan 06 cannot run the hexahedron probe without `pp/ppmail_qt` and
  `pp/ppinit`.
- **Fix:** Copied `pp/pp{mail,init,mail_qt}` from the main clone
  (`/home/mefisto/git/mefisto/pp/`) into the worktree's `pp/`. Justified
  because (a) the main clone's HEAD is at the same commit as the
  worktree's reset-target, so the binaries match the source tree
  byte-for-byte; (b) full `bin/cbl_tout_qt` rebuild would take ~7 min
  to produce identical binaries; (c) Plan 02 SUMMARY documents the
  same recovery pattern (Deviation #3).
- **Files modified:** `pp/pp{mail,init,mail_qt}` (all gitignored, no
  commit impact).
- **Verification:** `git status --short pp/` returned empty after copy
  (gitignore honored); the probe ran end-to-end against the copied
  binaries.
- **Committed in:** N/A (gitignored, no commit needed).

**2. [Rule 2 — Missing critical interpretation] 'Meaningful capture' acceptance test added to probe step 9**

- **Found during:** Task 2 step 0 (hexahedron probe outcome interpretation).
- **Issue:** The plan's probe step 9 acceptance criterion (`capture_size >
  0 AND exit ∈ {0, 124-with-capture}`), taken literally, would mark all
  4 crash candidates (ln/ob/sf/vl, all 2045 B) as CAPTURED — even though
  they show only the canonical Mefisto peach background `#FFD9B8` (1
  unique color, no mesh, no labels, no data). Promoting any of these to
  the CAPTURE branch would produce a false positive and a fabricated
  hexahedron capture — exactly the BLOCKER #1 / T-08-35 spoofing risk
  the alphabetical-order probe is designed to prevent.
- **Fix:** Added the Plan 02 §"Color-content sanity check" rule as the
  corrective acceptance test: a capture with only 1-2 grayscale (or 1
  single-background) entries indicates xvfermer_ fired before any data
  was painted into the backing pixmap. Documented this interpretation
  explicitly in hexahedron-probe.md §Outcome to forestall a future agent
  thinking the 2045-B captures qualify.
- **Files modified:** none (interpretation only; no source-tree change).
- **Verification:** Each candidate's outcome row in hexahedron-probe.md
  records the capture content explicitly ("1-color peach background only")
  alongside the size; a future agent reading the table can verify the
  meaningful-content test was applied.
- **Committed in:** `2f942dc` (Task 2 commit).

**3. [Rule 2 — Missing critical step] INITIER workspace-seed step added before each probe**

- **Found during:** Task 2 step 0 first probe attempt — `pp/ppmail_qt ln`
  ran but immediately CRASHed because the workspace lacked the ms10–ms15
  super-files that pp/ppinit creates.
- **Issue:** The plan's probe recipe (step 8) sets up the workspace via
  `cp -r testa/hexahedron/* $MEFISTOX/hexahedron/` but does NOT include
  the `echo hexahedron | pp/ppinit` seeding step. Without ppinit, every
  probe candidate would CRASH on missing ms-files, masking the real
  signal (whether the candidate file is a viable batch driver).
- **Fix:** Added `(cd $MEFISTOX/hexahedron && echo hexahedron |
  pp/ppinit) >/dev/null` after each workspace re-seed and before the
  pp/ppmail_qt invocation. This matches the documented Mefisto workflow
  (1: INITIER → 2: MAILLER) and the contract documented in
  bin/phase8_case_batch_map.sh §"WORKSPACE DISCIPLINE". Documented in
  hexahedron-probe.md §"Probe attempts" preamble.
- **Files modified:** none (probe-script local to executor; no source-tree
  change).
- **Verification:** After the fix, the 5 candidates produce DIFFERENT
  failure modes (4 CRASH, 1 HANG) — confirming the probe is now sensitive
  to the actual file-content viability, not blocked by the missing ms-file
  setup.
- **Committed in:** `2f942dc` (Task 2 commit; the resulting probe trail
  in hexahedron-probe.md reflects the corrected setup).

---

**Total deviations:** 3 auto-fixed (1 blocking, 2 missing-critical, 0
architectural). All resolved at the executor level (binary copy + probe
interpretation + workspace seed); no harness mutation, no source-tree
mutation, no testa/ mutation.

**Impact on plan:** All auto-fixes essential to honor the plan's verify
gates. The binary copy is the documented recovery pattern (Plan 02 + 01
SUMMARY both record it). The meaningful-capture test is a Rule-2 mitigation
of T-08-35 (spoofing risk via fabricated hexahedron capture) — without it,
a future agent could mistakenly proceed to the CAPTURE branch with any of
the 4 crash captures. The INITIER seed step is a Rule-2 missing-critical
addition that aligns the probe with the documented workspace-discipline
contract.

## Issues Encountered

- **Worktree at fresh main HEAD (no recovery needed despite initial
  diagnosis):** The orchestrator's per-agent branch was created at the
  same commit as main (`859941e`); `git fetch . main && git reset --hard
  FETCH_HEAD` was a no-op. The bash agent's CWD reset between commands
  produced an apparent "FATAL: HEAD on main" once (when a `cd
  /home/mefisto/git/mefisto` command moved the shell to the main worktree).
  Recovered by re-issuing absolute path commands in the per-agent worktree
  CWD. Three colorbar reports were initially written to the main worktree's
  filesystem path (because the Write tool used absolute path
  `/home/mefisto/git/mefisto/.planning/...`); fixed by `mv` into the per-agent
  worktree. No commit impact.

- **5-min probe budget with bash command timeouts:** running a 60-s
  per-candidate probe sequentially across 5 candidates risked exceeding
  the bash tool's 120-s default timeout. Mitigated by running the 4
  remaining candidates (after the first probe) in a single 360-s
  bash invocation with explicit per-candidate elapsed-time tracking and
  early-exit guard at 300 s.

## User Setup Required

None — no external service configuration, no credentials. The hexahedron
probe ran headless on the build host's existing Xvfb + Qt 6 (`offscreen`
QPA platform) installation.

## Next Phase Readiness

**Plan 07 (CHECKLIST.md finalize) unblocked.** Plan 07 ingests
spot-check-summary.md verbatim and:

1. Encodes the 4 CHECK verdicts into per-cell Notes annotations (per the
   table in §"Hand-off to Plan 07"):
   - pan2d: `pitfall-7-confirmed; ae-budget=3000-confounded-by-resample`
   - nafems_le1: `pitfall-6-disconfirmed; fuzz_override_pct=8;
     dim-mismatch=resample-to-1280x800`
   - cavity2d: `pitfall-6-secondary; fuzz_override_pct=10; ...;
     colorbar-region-AE-pct=21.94`
   - heat1d: `pitfall-6-disconfirmed; fuzz_override_pct=8;
     dim-mismatch=resample-to-1280x800`
2. Adds a separate "Spot-check rows" section (NOT in the canonical 5-case
   verdict matrix per T-08-23 mitigation) with the row:
   `| hexahedron | mail | GAP | escalation per evidence/hexahedron-probe.md |`
3. Carries the `08-VALIDATION.md ## Open escalations` entry into the v1
   ship-gate sign-off matrix as a maintainer-decision item: maintainer's
   sign-off either accepts the gap with rationale-(b) ("pan2d coverage is
   sufficient for AA-drift gate") or refuses and demands gap-closure-(a)
   (CONTEXT.md amendment naming canonical hexahedron driver).

**Phase 9 deferred-idea recorded:** matched-dim Qt recapture (1280x800
in-process backing pixmap if `MEFISTO_QT_WINDOW_SIZE` or equivalent is
supported) to retire the resample-dim-mismatch confound and apply the
Pitfall 6/7 fuzz overrides + AE budgets cleanly to like-for-like compares.

## Self-Check: PASSED

All claims verified after writing this summary:

- `evidence/colorbar-nafems_le1.md`: FOUND (118 lines, sha256 `ea0fe79e...`),
  5 sections present, fuzz=5 / fuzz=8 AE rows present.
- `evidence/colorbar-heat1d.md`: FOUND (112 lines, sha256 `1a7ef0f1...`),
  5 sections present.
- `evidence/colorbar-cavity2d.md`: FOUND (128 lines, sha256 `ed618d46...`),
  5 sections present, full-frame + cropped colorbar tables both present.
- `evidence/font-pan2d.md`: FOUND (154 lines, sha256 `c17019f2...`), 5
  sections present, 3545-color empirical finding recorded.
- `evidence/font-hexahedron.md`: FOUND (59 lines, sha256 `a0fc9f6b...`), 5
  LITERAL section headings present (Captures / Diff statistics / AE budget
  / Hypothesis / Verdict), GAP-DOCUMENTED verdict text present.
- `evidence/hexahedron-probe.md`: FOUND (241 lines, sha256 `da83cc49...`),
  5 candidates recorded in alphabetical order, Outcome=GAP line present.
- `evidence/spot-check-summary.md`: FOUND (128 lines, sha256 `30aaa765...`),
  5-row Aggregated Verdicts table + 5 sections + Outcome line present.
- `08-VALIDATION.md`: FOUND (80 lines, sha256 `7fe61eb5...`), `## Open
  escalations (filed by Plan 06 Task 2)` heading present, hexahedron entry
  present.
- Commits `f507f5e` / `2f942dc` / `6e8db45` FOUND in `git log --oneline`.
- Phase 8 invariants verified: `git diff -- xvue/xvuelc.c` empty;
  `git diff -- xvue/qt/src/` empty.
- Verify gates: Task 1 (5 sections × 3 colorbar reports), Task 2 (5
  literal sections × 2 font reports + Outcome line in probe + GAP
  escalation reachable), Task 3 (5 sections + 5-row Aggregated Verdicts
  table) all PASS.

---

*Phase: 08-ab-validation-on-testa-subset*
*Plan: 06 (VALID-05 + VALID-06 spot-check sweep)*
*Completed: 2026-05-05*
