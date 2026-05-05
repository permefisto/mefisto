# Phase 8 Validation Notes

**Owner:** Plan 07 (CHECKLIST.md finalize) — except for the
`## Open escalations` section below, which is filed by Plan 06 Task 2 as
documented in `08-06-PLAN.md` step 12c (the single exceptional case where
Plan 06 modifies this file). Plan 07 will preserve the escalation entry
verbatim and integrate it into the v1 ship-gate sign-off matrix as a
maintainer-decision item.

## Open escalations (filed by Plan 06 Task 2)

### VALID-06 hexahedron font-metric spot check — GAP

**Filed:** 2026-05-05 by Phase 8 Plan 06 Task 2 (probe-then-document branch
per BLOCKER #1 / MAJOR-9 iter1 review; WARNING-3 / WARNING-5 iter2 review).

**Probe report:** `evidence/hexahedron-probe.md`
**Spot-check report:** `evidence/font-hexahedron.md`

**Summary:** VALID-06 hexahedron font-metric spot check could not be
completed because no headless batch driver was found in the time budget
of Plan 06. `testa/hexahedron/` contains 5 ASCII batch script files (`ln`,
`ob`, `pt`, `sf`, `vl`) with no recognized extension. All 5 candidates
were probed in alphabetical order under
`pp/ppmail_qt <D>` + `QT_QPA_PLATFORM=offscreen` + `MEFISTO_QT_CAPTURE_PATH`
+ `MEFISTO_XVSOURIS_AUTOEXIT=1` + `timeout 60`:

- `ln`, `ob`, `sf`, `vl`: `STOP Sorry, CRASH of MEFISTO software` at
  unknown-name lookup (each file references named entities defined in
  EARLIER stages of the construction pipeline; no single file is
  self-sufficient). 2045-byte all-peach background capture at crash time
  (xvfermer_ fired before any drawing primitive painted the mesh).
- `pt`: process HUNG at 60-s timeout (file has no rendering trigger;
  expects the next interactive menu pick to draw). 0-byte capture
  (xvfermer_ never fired).

**Architectural reading:** The 5 files are NOT independent batch entry
points — they are interlocking building blocks of a single hexahedron
mesh definition (pt → ln → sf → vl → ob). The legacy MAILLER protocol's
single-batch invocation pattern (`pp/ppmail batchfile`) does not apply
here.

**Maintainer decision required.** Maintainer must EITHER:

  **(a)** Issue a CONTEXT.md amendment naming the canonical driver. This
  could be:
  - A pointer to a wrapper script (e.g., `bin/hexahedron_drive.sh`) that
    issues the menu picks via stdin pipe to `pp/ppmail_qt` after seeding
    the workspace.
  - A documented stdin-script that drives the interactive menu through
    pt → ln → sf → vl → ob → drawings (close the workflow).
  - A new `pp/ppmail_qt` multi-file batch protocol that accepts a
    sequence of batch files (e.g., `pp/ppmail_qt pt ln sf vl ob`).

  **OR**

  **(b)** Waive VALID-06 hexahedron coverage with documented rationale
  (e.g., "pan2d coverage is sufficient for the AA-drift gate; hexahedron
  is too marginal to require headless A/B coverage given the maintenance
  cost of building a multi-file driver protocol"). The waiver should be
  recorded in `REQUIREMENTS.md` against VALID-06 with the maintainer's
  initials and date.

**Phase 7 sign-off is gated on this resolution.** Plan 07 will read this
escalation entry and emit a `hexahedron | mail | GAP | escalation per
evidence/hexahedron-probe.md` row in the CHECKLIST.md "Spot-check rows"
section (NOT the canonical 5-case grid — per T-08-23 mitigation,
hexahedron remains a spot-check-only target). The maintainer's Plan 07
sign-off either accepts the gap with rationale-(b) or refuses and demands
gap-closure-(a).

**See also:**
- `evidence/hexahedron-probe.md` — full empirical probe trail (5 candidates
  tried in alphabetical order with per-candidate exit code, capture size,
  log lines, crash markers).
- `evidence/font-hexahedron.md` — GAP-branch spot-check report (preserves
  all 5 literal section headings per WARNING-3 iter2 with N/A bodies
  pointing at hexahedron-probe.md).
- `08-06-PLAN.md` step 12c — the contract authorizing this escalation
  filing.
