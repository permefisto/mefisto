# hexahedron headless-driver probe (BLOCKER #1 mitigation)

**Date:** 2026-05-05
**Plan:** Phase 8 Plan 06, Task 2 step 0 (per BLOCKER #1 / MAJOR-9 iter1 review;
WARNING-3 / WARNING-5 iter2 review).
**Probe order:** alphabetical across {`ln`, `ob`, `pt`, `sf`, `vl`} per
WARNING-5 iter2 (drops the iter1 "vl-first" preference that lacked source
grounding).
**Total budget:** 5 minutes.
**Per-candidate budget:** 60 seconds (`timeout 60`).
**Capture mechanism:** `MEFISTO_QT_CAPTURE_PATH` (in-process Qt offscreen
backing pixmap → PNG via `xvfermer_`).
**Driver invocation:** `pp/ppmail_qt <D>` (with `<D>` ∈ {ln, ob, pt, sf, vl}).
**Workspace:** `/tmp/mefistox-phase8-spot/hexahedron/` (re-seeded with
`echo hexahedron | pp/ppinit` before each probe).

## ls + file output for testa/hexahedron/

```
$ ls -la testa/hexahedron/
-rwxrwxr-x 1 mefisto mefisto 399  ln
-rwxrwxr-x 1 mefisto mefisto 728  ob
-rwxrwxr-x 1 mefisto mefisto 341  pt
-rwxrwxr-x 1 mefisto mefisto 300  sf
-rwxrwxr-x 1 mefisto mefisto 196  vl

$ file testa/hexahedron/{ln,ob,pt,sf,vl}
testa/hexahedron/ln: ASCII text
testa/hexahedron/ob: ASCII text
testa/hexahedron/pt: ASCII text
testa/hexahedron/sf: ASCII text
testa/hexahedron/vl: ASCII text
```

5 ASCII batch script files; NO recognized extension; total size 1964 bytes.
Names are abbreviations of mesh-element types:
- `pt` = points (defines vertex coords A0…D0, A1…D1, DA1)
- `ln` = lines (defines edges A0B0, B0C0, …, D1A1, DA1; references points)
- `sf` = surfaces (defines faces ABCD0…DADA; references lines)
- `vl` = volumes (defines hexa volume; references surfaces)
- `ob` = objects (composes the entire hexahedron object; references all)

**Architectural reading:** The 5 files are NOT independent batch entry points
— they are **interlocking building blocks** of a single hexahedron mesh
definition. Each file references named entities defined in earlier-stage
files. The natural construction order is `pt → ln → sf → vl → ob`, which
is what the interactive Mefisto-MAILLER menu orchestrates.

## head -10 of {ln, ob, pt, sf, vl}

### `ln`
```
@;
@;
2; {LINES}
 A0B0; i; 2; nx; 1; A0; B0;
 B0C0; i; 2; ny; 1; B0; C0;
 C0D0; i; 2; nx; 1; C0; D0;
 D0A0; i; 2; ny; 1; D0; A0;

 A1B1; i; 2; nx; 1; A1; B1;
 B1C1; i; 2; ny; 1; B1; C1;
```

### `ob`
```
@;
@;
5; { OBJECTS parts of the object }
 hexahedron; 0; { classic 2D ou 3D }
  7; { number of Points Lines Surfaces Volumes Objects of the object }
 v; hexahedron; { POINT or LINE or SURFACE or VOLUME or OBJECT? }
 s; ABCD0;        { name of type PLSVO and name of PLSVO}
 s; ABCD1;
 s; ABAB;
 s; BCBC;
```

### `pt`
```
DefVar nx, ny, nz;
nx=8;
ny=6;
nz=4;


1; {POINTS}
 1;
 A0; i; 1;  0;  0;  0;
```

### `sf`
```
@;
@;
3; {SURFACES}
ABCD0; i; 1; A0B0; B0C0; C0D0; D0A0; 1;
ABCD1; i; 1; A1B1; B1C1; C1D1; D1A1; 1;
ABAB;  i; 1; A0B0; B0B1; A1B1; A0A1; 1;
BCBC;  i; 1; B0C0; C0C1; B1C1; B0B1; 1;
CDCD;  i; 1; C0D0; D0D1; C1D1; C0C1; 1;
DADA;  i; 1; D0A0; A0A1; D1A1; D0D1; 1;
@;
```

### `vl`
```
@;
@;
4; { VOLUMES materials of the object }
  hexa;
  i;{ in the lexicon ~>TRANSFO }
  1;{ transfinite hexahedron }
ABCD0; ABCD1; ABAB; BCBC; CDCD; DADA;
@;

10; {drawings}
```

## Probe attempts

Probes run in alphabetical order with workspace re-seeded between each
(WARNING-5 iter2): `rm -rf $MEFISTOX/hexahedron && mkdir -p
$MEFISTOX/hexahedron && cp -r testa/hexahedron/* $MEFISTOX/hexahedron/`,
then `(cd $MEFISTOX/hexahedron && echo hexahedron | pp/ppinit)` to seed the
ms10–ms15 super-files, then the probe invocation below:

```
(cd $MEFISTOX/hexahedron && \
  env QT_QPA_PLATFORM=offscreen \
      MEFISTO_QT_CAPTURE_PATH=/tmp/_hex_probe_${D}.png \
      MEFISTO_BATCH_X11=1 \
      MEFISTO_XVSOURIS_AUTOEXIT=1 \
      MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=500 \
      timeout 60 $MEFISTO/pp/ppmail_qt ${D})
```

| # | Candidate | Exit | timeout 60s | Capture size | Capture content                          | Log lines | Crash/Error markers | Outcome      |
|---|-----------|------|-------------|--------------|------------------------------------------|-----------|---------------------|--------------|
| 1 | `ln`      | 0    | no (≈1s)    | 2045 B       | 1-color peach background only            | 54        | `STOP Sorry, CRASH` | NOT-VIABLE   |
| 2 | `ob`      | 0    | no (≈1s)    | 2045 B       | 1-color peach background only            | 54        | `STOP Sorry, CRASH` | NOT-VIABLE   |
| 3 | `pt`      | 124  | YES (60s)   | 0 B          | (no capture — process hung at 60s timeout) | 51        | (no error/crash)    | NOT-VIABLE   |
| 4 | `sf`      | 0    | no (≈1s)    | 2045 B       | 1-color peach background only            | 54        | `STOP Sorry, CRASH` | NOT-VIABLE   |
| 5 | `vl`      | 0    | no (≈1s)    | 2045 B       | 1-color peach background only            | 54        | `STOP Sorry, CRASH` | NOT-VIABLE   |

**Total elapsed:** 61 seconds (well under the 5-min budget; 60s consumed by
the `pt` timeout, others returned in <1s with crash).

### Per-candidate detail

**`ln` (line definitions):** `pp/ppmail_qt ln` opens the Qt window
(760x442), parses argv[1]="ln", enters BATCH mode, reads the `ln` file.
The file references named points `A0`, `B0`, etc. that are NOT yet defined
in the workspace (because the points file `pt` was never sourced). The
mesher emits `STOP Sorry, CRASH of MEFISTO software` after the unknown-name
lookup fails. The 2045-byte PNG is the empty Qt canvas at crash time
(only `#FFD9B8` peach background — confirmed by `convert -unique-colors`).

**`ob` (objects):** Same crash signature as `ln`. The `ob` file references
surface names `ABCD0`, `ABCD1`, etc. that are not defined; CRASH at
unknown-name lookup; 2045-byte all-peach capture.

**`pt` (points):** `pp/ppmail_qt pt` opens the Qt window, parses argv[1]="pt",
processes the `DefVar nx, ny, nz` and the 9 points definitions, hits the
`19; {min MAX of XYZ of PLSVO}` block (which sets bounding box only — NO
graphical render call), then **HANGS** waiting for the next interactive
menu pick. The 60-second timeout fires; `xvfermer_` is never reached;
no PNG is written (size=0).

**`sf` (surfaces):** Same crash signature as `ln`/`ob`. References lines
`A0B0`, `B0C0`, etc. that are not defined; CRASH at unknown-name lookup;
2045-byte all-peach capture.

**`vl` (volumes):** Same crash signature as `ln`/`ob`/`sf`. References
surfaces `ABCD0`, `ABCD1`, etc. that are not defined; CRASH at unknown-name
lookup; 2045-byte all-peach capture.

### Sample log tail (representative — `ln` candidate)

```
 Screen PIXELS     WIDTH= 800  HEIGHT= 800
 Window PIXELS     WIDTH= 760  HEIGHT= 442
 X-CONVERSION MM=>PIXELS =   3.94088674
 RESERVED COLORS from    0 to   15
 ALLOWED  COLORS from   16 to  255
 RECOVERY of           10  X11-CHARACTER FONTS

xvue-qt: stub nomrepmefisto_ not implemented yet
xvue-qt: stub secondes1970_ not implemented yet
STOP Sorry, CRASH of MEFISTO software
QThreadStorage: entry 1 destroyed before end of thread 0x562decdfb050
QThreadStorage: entry 0 destroyed before end of thread 0x562decdfb050
```

The `xvue-qt: stub nomrepmefisto_ not implemented yet` line is unrelated
to the crash — it is a Phase-7-known list of unimplemented Qt-side helper
stubs that the X11 backend implements but Qt does not. The actual crash
trigger is the unknown-name lookup, occurring deep in the Fortran
mesher's PLSVO directory after the stub messages.

## Outcome

**Outcome = GAP**

No alphabetical candidate {`ln`, `ob`, `pt`, `sf`, `vl`} produces a
**meaningful** end-to-end capture under the documented headless probe
contract:

- `ln`, `ob`, `sf`, `vl`: CRASH at unknown-name lookup (because each
  references entities defined in EARLIER stages of the construction
  pipeline; no single file is self-sufficient).
- `pt`: HANGS at 60-s timeout (because `pt` has no rendering trigger;
  it expects the next menu pick to draw).

**Determining factor:** The 5 files are interlocking parts of a single
hexahedron mesh definition (pt → ln → sf → vl → ob), not independent
single-file batch entry points. The legacy MAILLER protocol's
single-batch invocation pattern (`pp/ppmail batchfile`) does not
apply: there is no canonical batch driver in `testa/hexahedron/`.

**Acceptance test for "non-empty" capture (per probe step 9 logic):**
- All 4 of {ln, ob, sf, vl} produce capture_size > 0 (2045 B), but
  the content is a 1-color peach background — NO mesh, NO labels, NO
  data — failing the Plan 02 §"Color-content sanity check" rule:
  *"a capture with only [...] 1-2 grayscale entries [...] indicates
  xvfermer_ never fired the XCopyArea-to-fenetre path"*. The 2045 B
  captures are the equivalent: xvfermer_ fired AT crash time, but
  before any drawing primitive painted into the backing pixmap.
- `pt` produces capture_size = 0 (no PNG written; process hung).

By the meaningful-content criterion, no candidate qualifies. **Probe
outcome = GAP**, branching to step 12 (document-the-gap path).

## Hand-off

The hexahedron VALID-06 row escalates to maintainer review per
`08-VALIDATION.md ## Open escalations` (filed by Plan 06 Task 2). Plan 7
ingests the escalation and either accepts the gap (waiver-by-rationale)
or refuses and demands gap-closure plan (e.g., a CONTEXT.md amendment
naming the canonical aggregator file or adding a `pp/ppmail_qt`
multi-file batch protocol).
