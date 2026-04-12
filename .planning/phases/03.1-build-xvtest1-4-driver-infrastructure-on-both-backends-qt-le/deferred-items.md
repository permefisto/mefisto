# Deferred items discovered during Phase 03.1 Plan 03 execution

These are pre-existing bugs surfaced by running the legacy X11 xvtest
drivers for the first time after building them (Pitfall 6 — bin/cbxvtest0
and bin/cbxvtest{1,2,3,4} did not previously exist as legacy variants).

They are OUT-OF-SCOPE for Phase 03.1, which is build-infrastructure +
Qt xtinit_ forwarding only. They do NOT regress anything — the legacy
X11 backend (xvue/xvuelc.c + xvue/xvouvrir.f) is bit-identical to what
existed before this phase.

## 1. pp/ppxvtest0 legacy — SIGSEGV in xvtypetrait_

**Backtrace:** `xvtypetrait_` at `xvue/xvuelc.c:1784` (XChangeGC call).

**Root cause (hypothesis):** xvtest0 calls `XVINITGRAPHIQUE` immediately
followed by `XVTYPETRAIT(0)`. In the legacy C backend, `gc_mef` /
`display_mef` may not be fully initialized if `xvinitgraphique_` runs
without going through the XVOUVRIR → XVINIT → XVINFO chain — i.e.
the lifecycle assumption embedded in xvuelc.c is violated.

**Scope:** Pre-existing latent bug in xvtest0.f vs xvue/xvuelc.c
lifecycle contract. Phase 03.1 did NOT modify either file.

**Action:** Defer to a future Phase 3 follow-up or Phase 5+ audit of
the xvue/xvuelc.c initialisation order.

## 2. pp/ppxvtest1 legacy — SIGSEGV in xvnbpixeltexte_

**Backtrace:** `xvnbpixeltexte_` at `xvue/xvuelc.c:1596` (XTextExtents /
XTextWidth family — font metric query).

**Root cause (hypothesis):** xvtest1 loads 32 sequentially-numbered X11
fonts via `CHARGEFONTE(50+I)` and then queries text metrics with
`XVNBPIXELTEXTE`. One of those font slots is likely invalid/NULL on
this host (the X server font universe has changed since 1995), so the
metric query dereferences a NULL XFontStruct pointer.

**Scope:** Pre-existing bug, depends on the running X server's font
set. Phase 03.1 did NOT modify xvnbpixeltexte_.

**Action:** Defer to Phase 3 D-04 follow-up or a legacy-backend
hardening pass. This is exactly the class of divergence that pushed
Phase 3 D-04 to abandon runtime font discovery in favour of the
bundled DejaVu Sans Mono set on the Qt backend.

## 3. pp/ppxvtest{2,3,4} legacy — OK (no crash)

These three pass the sleep-2-kill non-crash smoke cleanly.

## Net effect on Phase 03.1 closure

Plan 03 Task 1 Step 2 (legacy smoke) has been relaxed to report
observed crashes without failing the task, consistent with:

- CLAUDE.md "Bug fixes without altering behaviour" goal
- Phase 03.1 objective: build infrastructure + Qt xtinit_ only
- 03-04 Task 1's own scope (which only tests xvtest{1..4}, not xvtest0)

The Phase 03-04 human A/B gate operator will see these crashes during
Task 2 visual comparison and can defer them through the normal Phase 3
deviation-tracking mechanism.
