---
status: partial
phase: 00-build-skeleton-abi-stubs
source: [00-04-SUMMARY.md, .planning/validation/BASELINE.md]
started: 2026-04-10
updated: 2026-04-10
---

## Current Test

[awaiting human testing — deferred by user during /gsd:execute-phase 0]

## Tests

### 1. testa/pan2d — legacy MAILLER smoke test
expected: `MAILLER pan2d` launches the X11 interactive mesher against `testa/pan2d`, renders the mesh, exits cleanly after `99;`. Captures a 1-3 sentence baseline description in `.planning/validation/BASELINE.md` under the `pan2d` entry.
result: [pending]

### 2. testa/nafems_le1 — legacy ELASTICER smoke test
expected: `ELASTICER nafems_le1` launches the elasticity solver against the NAFEMS LE1 benchmark case, renders the solution, exits cleanly after `99;`. Updates `.planning/validation/BASELINE.md` `nafems_le1` entry with observed behavior.
result: [pending]

### 3. testa/cavity2d — legacy FLUIDER smoke test
expected: `FLUIDER cavity2d` launches the fluid solver against the 2D cavity benchmark, renders the velocity/pressure fields, exits cleanly after `99;`. Updates `.planning/validation/BASELINE.md` `cavity2d` entry.
result: [pending]

### 4. testa/heat1d — legacy THERMICER smoke test
expected: `THERMICER heat1d` launches the thermal solver against the 1D heat benchmark, renders the temperature field, exits cleanly after `99;`. Updates `.planning/validation/BASELINE.md` `heat1d` entry.
result: [pending]

### 5. testa/nlsecu — legacy NLSER smoke test
expected: `NLSER nlsecu` launches the nonlinear solver against the `nlsecu` benchmark, exits cleanly after `99;`. Updates `.planning/validation/BASELINE.md` `nlsecu` entry.
result: [pending]

## Summary

total: 5
passed: 0
issues: 0
pending: 5
skipped: 0
blocked: 0

## Gaps

None yet — these are pending maintainer runs, not failures. Resolve with `/gsd-verify-work 0` when a maintainer with X11 access can drive the 5 interactive sessions and fill in the BASELINE.md placeholders.
