# Font-metric spot check — hexahedron (VALID-06; GAP-DOCUMENTED per BLOCKER #1)

**Date:** 2026-05-05
**Plan:** Phase 8 Plan 06, Task 2 (Part B) — GAP branch (step 12).
**Probe report:** `evidence/hexahedron-probe.md` (5 candidates probed in
alphabetical order per WARNING-5 iter2; all 5 NOT-VIABLE).
**Escalation:** `08-VALIDATION.md ## Open escalations` (filed by this plan).

**Important context:** hexahedron is named in VALID-06 as a font-metric
spot-check target but is NOT in the BUILD-10 baseline (the 5 canonical
testa cases: pan2d, nafems_le1, cavity2d, heat1d, nlsecu). The iter1
review (BLOCKER #1, MAJOR-9) flagged that `testa/hexahedron/` contains 5
ASCII batch script files (`ln`, `ob`, `pt`, `sf`, `vl`) with no obvious
single-file batch driver. Task 2 step 0 probed empirically (alphabetical
order, 5-min budget); the probe outcome was **GAP** — no candidate
produced a meaningful end-to-end capture. This file is the GAP-branch
report; per WARNING-3 iter2, **all 5 literal section headings are
preserved below even when the body is "N/A — see hexahedron-probe.md"**.

## Captures

N/A — see hexahedron-probe.md for probe details. No captures produced
because no headless driver was found in the 5-min probe budget across
{ln, ob, pt, sf, vl} alphabetical order.

(Detail: `ln`, `ob`, `sf`, `vl` produced 2045-byte all-peach background
captures at CRASH time — xvfermer_ fired before any drawing primitive
painted the mesh. `pt` produced size=0 because the process hung at the
60-s timeout — `pt` has no rendering trigger. Per the Plan 02 color-content
sanity check rule, none of the 5 captures is a meaningful frame.)

## Diff statistics

N/A — see hexahedron-probe.md for probe details. (No captures → no diff
to compute.)

## AE budget

N/A — see hexahedron-probe.md for probe details. (Reserved budget would
have been 1600 px per Pitfall 7 — 8 hexahedron vertex labels × 200 px
per label — applied to the eventual capture once the maintainer resolves
the gap.)

## Hypothesis

N/A — see hexahedron-probe.md for probe details. testa/hexahedron's
driver chain may require interactive-menu input that AUTOEXIT cannot
satisfy; or the case may need a CONTEXT.md amendment naming the canonical
driver explicitly.

(Architectural note: the 5 files are interlocking building blocks
[pt → ln → sf → vl → ob] of a single hexahedron mesh definition, not
independent single-file batch entry points. The legacy MAILLER protocol's
single-batch invocation `pp/ppmail batchfile` does not apply.)

## Verdict

GAP-DOCUMENTED — VALID-06 hexahedron coverage UNRESOLVED, escalation filed
in 08-VALIDATION.md `## Open escalations` section.
