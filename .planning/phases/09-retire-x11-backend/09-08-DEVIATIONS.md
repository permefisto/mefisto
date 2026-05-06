# Phase 9 Plan 9-08 — Execution Deviations

Plan 9-08 attempted the cross-tag worktree procedure (per RESEARCH §Pattern 3)
to produce three Phase-7-deferred goldens. The cross-tag worktree itself was
built successfully (Task 1) and the X11 backend in `/tmp/mefisto-pre-retire`
was functional. However, **two additional defects beyond the originally
diagnosed "headless / multi-step" blockers** prevented the goldens from being
produced. Both are pre-existing carry-forwards from Phase 7 design and
neither can be solved within Plan 9-08's mandate of "no source changes
outside golden + log".

## Deviation 1 — `[Rule 3 - Blocker]` scene01_driver.f link line uses non-existent loose `.o` files

**Found during:** Task 2, link step.

**Issue:** The plan (Task 2 step 1) and `scene01_driver.f` header (lines 38-40)
specify the link line:

```bash
gfortran scene01_driver.o $MEFISTO/xvue/xvuelc.o \
    $MEFISTO/xvue/*.o $MEFISTO/util/*.o \
    -L/usr/X11R6/lib -lX11 -lXt -o scene01_x11
```

This was authored against an assumed build state where every Fortran source
in `xvue/*.f` and `util/*.f` had a corresponding individual `.o` file. The
pre-retire build (`bin/cbl_tout`) does **not** preserve loose `.o` files —
it ar-archives every Fortran object into `xvue/lib`, `util/lib`,
`mail/lib`, etc. and deletes the loose `.o` files. After
`bin/cbl_tout` exit 0, only `xvue/xvuelc.o` (the C-compiled file)
remains alongside the archives. The glob `$MEFISTO/xvue/*.o` matches
**only** `xvuelc.o`, and `$MEFISTO/util/*.o` matches nothing — the
linker reports unresolved Fortran symbols.

**Fix applied (in `/tmp/09-08-scene01/`, NOT in main):** Substitute the
canonical archive form per `bin/cbinit` line 50:

```bash
gfortran scene01_driver.o $MEFISTO/xvue/xvuelc.o \
    $MEFISTO/xvue/lib $MEFISTO/util/lib \
    -lgfortran -L/usr/lib/x86_64-linux-gnu -lX11 -lXt -o scene01_x11
```

(Also adjusted the libX11 path from the plan's `/usr/X11R6/lib` —
which does not exist on Debian — to the canonical
`/usr/lib/x86_64-linux-gnu` location.)

**Result:** Link succeeded (exit 0; `scene01_x11` binary 144,456 bytes).

**Implication:** The `scene01_driver.f` BOOTSTRAP NOTE comment block (lines
33-46) and the Plan 9-07 RESEARCH §Pattern 3 example (RESEARCH lines
255-258) are both stale. They document a link form that has never worked
post-`cbl_tout` because the Fortran build never produces the loose `.o`
files those globs assume. A future maintainer running the documented
procedure verbatim would hit the same blocker. **Recommend:** update the
`scene01_driver.f` header in Phase 7 (or in a follow-up Plan 9-09 docs
patch) to reflect the archive form. This deviation does **not** require
a code commit in Plan 9-08.

## Deviation 2 — `[Rule 4 - Architectural]` scene01_driver.f program body has fundamental bugs that prevent the X11 backend from rendering it

**Found during:** Task 2, run step.

**Issue:** Even with the link fixed (Deviation 1), running `scene01_x11`
under `xvfb-run` produces a SIGSEGV in `xvepaisseur_` at `xvuelc.c:1882`
(`XChangeGC(display_mef, gc_mef, ...)`). Root cause: the driver's init
sequence is incomplete. `scene01_driver.f:62` calls only
`XVINITGRAPHIQUE`, which (per `xvuelc.c:286`) opens the display + sets
`screen_mef` but does **not** create `gc_mef` or `fenetre_mef`. Both are
required by every subsequent emit call. The proper init sequence is
`XTINIT → XVINIT` (per `prpr/ppmail.f:163-166`); `XVINIT` calls
`XVINFO(...)` which creates the GC + window.

**Workaround attempt:** Wrote `scene01_driver_v2.f` that calls
`XTINIT` then `XVINFO(800, 600, MAXFONTS, ...)` directly with stub
output args. With explicit Xvfb font-path
(`-fp /usr/share/fonts/X11/misc/,...`) this brought up 512 X11
fonts (`Nb de FONTES X11 disponibles = 512`) and got past
`xvepaisseur_`. The line/face/ellipse calls all worked.

**Second blocker — fundamental driver bug:** `XVCHARGEFONTE` at
`scene01_driver.f:107` is called with **5 arguments** including a
character literal:

```fortran
CALL XVCHARGEFONTE('Courier', 7, 12, 10, 0)
```

But the C-side `xvuelc.c:1474` declares it with **4 INTEGER arguments**:

```c
void proc(xvchargefonte)( int *nofont0, int *nofont, int *largpx, int *hautpx )
```

Confirmed by the only legitimate Fortran caller in the codebase
(`xvue/chargefonte.f:36`):

```fortran
CALL XVCHARGEFONTE( NOFONT0, NOPOCA, NPLACA, NPHACA )
```

When called with `'Courier'`, the C function dereferences the first
4 bytes of the string `"Cour"` interpreted as a 32-bit int —
`0x436F7572` = 1131508082 — uses it as `*nofont0`. The path
`if (*nofont0 > 0) XFreeFont(display_mef, struc_police)` then
calls `XFreeFont` on a NULL `struc_police`, which segfaults on the
next `XTextExtents(struc_police, ...)`.

This means **the driver as committed has never run end-to-end against
the X11 backend.** Phase 7 Plan 06 Task 3 was checkpointed as a
human-verify step, was deferred to Phase 8, was again deferred to Phase
9, and now Plan 9-08 has empirically confirmed the driver itself cannot
produce TEMPORAIRE.EPS without source modification.

**Fix considered, NOT applied:** Modify `scene01_driver.f` to use
correct `XVCHARGEFONTE(NOFONT0, NOFONT, LARGPX, HAUTPX)` arity. This
would require:
1. Adding declarations for the 4 INTEGER args.
2. Synchronously updating the Qt-side scene at
   `xvue/qt/tests/test_xvue_qt_postscript.cpp:477`
   (`pe.chargefonte("Courier", 12, 10, 3, false, false)`) so the
   byte-comparison still passes.
3. Re-evaluating the byte-parity claim — the entire
   "verbatim port from xvuelc.c:1187-1304" contract documented in
   Phase 7 RESEARCH would need re-verification.

This is a Phase-7 design defect requiring Phase-7-scope work to
resolve, not a Phase 9 retirement task. **Out of scope for Plan 9-08.**

**Result:** TEMPORAIRE.EPS not produced. The
`PsEmitter_postscriptVerbatim_golden` QSKIP slot must remain DEFERRED.
Updated VALIDATION-LOG.md row reflects the empirical evidence.

**Workaround attempt for partial coverage:** Wrote
`scene01_driver_v3.f` that omits `XVCHARGEFONTE` + `XVTEXTE` calls.
This produced a 1,529-byte / 14-line EPS containing `epais`, `typet`,
5 `S` (line) opcodes, 2 `F` (face) opcodes, and 1 `el` (ellipse)
opcode — a real EPS byte stream from the X11 backend. It is
**not** suitable as `scene01.eps` because the Qt-side test still emits
the chargefonte+texte opcodes; QCOMPARE would fail with a hard test
failure (worse than the current QSKIP). The artifact is preserved as
evidence at `/tmp/09-08-scene01/TEMPORAIRE.EPS` (cleaned up with
the worktree but content quoted in this deviations file for the
record).

## Deviation 3 — `[Rule 4 - Architectural]` testa/wave + testa/cavity2d cannot run autonomously

**Found during:** Task 3, wave bootstrap.

**Issue:** The testa cases are interactive multi-module workflows. The
pipeline is:

1. `INITIER` reads project name from stdin → invokes `pp/ppinit`
   interactively (multiple prompts: working dir, primary lexicon,
   ADAM lexicon size, …) → writes the project's secondary memory file.
2. `MAILLER project.mesh` (after `INITIER`) runs `pp/ppmail` with
   the `.mesh` data → reads more interactive prompts (mesh
   parameters, output options, save sequence) before exiting.
3. `FLUIDER project.wave` (or analogous) runs `pp/ppflui` →
   solver-time prompts + frame export options.
4. `bin/convertepsgif` (a one-line `convert -rotate 90 zfxy0*.eps
   -extent 980x550 cyl53zfxy.gif` wrapper) bundles the per-frame EPS
   files produced by step 3 into an animated GIF.

**Empirical attempt (`echo "wave" | INITIER`):**

```
Fortran runtime error: End of file
```

at the very first ppinit READ statement past the project name. ppinit
needs at least 4-5 additional responses (primary lexicon, ADAM lexicon
sizes, lecture mode, ...) before it can complete and persist the
secondary-memory file. The full interactive sequence for each testa
case is undocumented in the source tree (the menu prompts are emitted
inline by `pp/ppinit` and depend on the language flag); reverse-
engineering them per case is multi-hour work.

The plan's hardcoded `testa/wave/wave.iexrr` and `testa/wave/wave.iexsr`
files **do not exist**:

```
$ ls /tmp/mefisto-pre-retire/testa/wave/
wave.mesh  wave.wave
```

These two files (`.mesh`, `.wave`) are MEFISTO project-script
inputs (DefVar, surface definitions, … in MEFISTO's domain-specific
syntax) that are READ by MAILLER + FLUIDER **after** the interactive
project-init dialog completes — they do **not** drive the ppinit
prompt sequence themselves.

**Phase 8 Plan 1 Task 2 acknowledgement:**
"testa/wave + testa/cavity2d are multi-module legacy AUTOEXIT chains
needing human-issued `99;` saves between MAILLER and the solver step"
— Plan 9-08 has now empirically confirmed this is a stricter blocker
than that quote suggests. It is **not** just the inter-module `99;`
save; it is the entire ppinit interactive prompt sequence.

**Fix considered, NOT applied:**
1. Author per-case `.iex` (interactive-exchange) input files that
   pre-record every prompt response. Requires running each case
   under `script(1)` once, capturing the full prompt/response
   sequence, then driving each subsequent run with the captured file.
2. Or add a true non-interactive batch mode to `pp/ppinit` (a
   `-batch <projectname>` flag that uses defaults). This is a
   Phase-7-scope source change; Plan 9-08 is "no source changes
   outside golden + log".

**Out of scope for Plan 9-08.** Both `wave_legacy.gif` and
`cavity2d_legacy.gif` remain DEFERRED.

## Deviation 4 — `[Rule 1 - Bug]` Plan 9-08 cited path `/usr/X11R6/lib` does not exist on Debian

**Found during:** Task 1 step 5 verification, Task 2 step 1 link.

**Issue:** Plan Task 2 step 1 (and Task 1 step 5 implicit) hardcodes
`-L/usr/X11R6/lib`. This path is a 2000s SunOS / commercial-Unix
convention; on modern Debian/Ubuntu (the actual deploy target) X11
ships at `/usr/lib/x86_64-linux-gnu/`. The `-L` flag silently does
nothing (path does not exist), but the link still succeeds because
gcc default search paths include `/usr/lib/x86_64-linux-gnu/`
already.

**Fix applied:** Substituted `-L/usr/lib/x86_64-linux-gnu` (or simply
omitted `-L` and let gcc defaults take over). Documented for
future RESEARCH §Pattern 3 update.

## Final state

| Item                             | Result                                         |
| -------------------------------- | ---------------------------------------------- |
| `xvue/qt/tests/golden/scene01.eps` | NOT produced — driver bug (Dev 2)            |
| `xvue/qt/tests/golden/wave_legacy.gif` | NOT produced — interactive chain (Dev 3) |
| `xvue/qt/tests/golden/cavity2d_legacy.gif` | NOT produced — interactive chain (Dev 3) |
| `/tmp/mefisto-pre-retire` cross-tag worktree | CLEANED UP (T-09-03 mitigation upheld) |
| `/tmp/mefistox-pre-retire` user-project dir | REMOVED (`rm -rf`)                       |
| `bin/cbl_tout` in pre-retire worktree   | EXIT 0 (X11 backend functional there)     |
| Main worktree `git status`              | CLEAN — no `.o`/`pp/*` leakage             |

## Recommendation for Phase 9 close-out

The three goldens are **structurally unproducible** under current
scene01_driver.f + testa AUTOEXIT design. Closing them DEFERRED →
PASS requires either:

1. **Source repair:** Modify `scene01_driver.f` to use correct
   `XVCHARGEFONTE` arity + add `XTINIT/XVINFO` init; add canonical
   `.iex` batch-input files for `testa/wave` + `testa/cavity2d`;
   re-run cross-tag bootstrap.
2. **Test-side acceptance:** Update the Qt-side test fixtures
   (`PsEmitter_postscriptVerbatim_golden` + `XvueExport_gif_AB_compare_*`)
   to skip with a narrower message documenting that the upstream
   Fortran driver harness has known-broken init/font calls and that
   the QSKIP is by design (not a deferred verification).
3. **Phase-9 close ack:** Accept that Phase 7 carry-forward #3 is a
   permanent DEFERRED, document in the project SHIP-GATE log, and
   update the Phase 7 VALIDATION-LOG.md `Notes` section to reference
   Plan 9-08's empirical failure analysis.

Plan 9-08 leaves the rows DEFERRED with refined rationale (now
empirically grounded). Option 1 is multi-day Phase-7-scope work;
options 2 and 3 are recommended for orchestrator triage.
