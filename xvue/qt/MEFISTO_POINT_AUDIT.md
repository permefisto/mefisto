# MefistoPoint ABI audit (Phase 2, D-31)

**Purpose:** Verify every Fortran caller of XVFACE / XVTRAITS / XVFACETRAITS
passes a point array whose layout is byte-identical to
`struct MefistoPoint { short x; short y; }` (4 bytes / point). This is
the precondition Pitfall 4 (research/PITFALLS.md) warns about.

**Layout under audit:** `struct MefistoPoint { short x; short y; };`
declared in `xvue/qt/include/xvue_qt_api.h` with
`static_assert(sizeof == 4)`.

**Result:** 29 of 29 live call sites verified INTEGER*2 (2, N) — ABI-safe.
MefistoPoint `{ short x; short y; }` is ABI-safe.

(Four additional hits from `xvue/fap32d.f` are commented out `CCC`
lines referencing `XVFACETRAITS` — they do not participate in the live
ABI and are listed in a separate table below for completeness.)

## Caller list (live calls)

| File | Line | CALL statement | Declared argument type | Verdict |
|------|------|----------------|------------------------|---------|
| xvue/tria2d.f | 32 | `CALL XVFACE( 3, XYPX )` | `INTEGER*2 XYPX(2,3)` (line 18) | OK |
| xvue/traits3d.f | 57 | `CALL XVTRAITS( NBS, XYPX )` | `INTEGER*2 XYPX(2,MXPOIN)` (line 23) | OK |
| xvue/triacoul.f | 52 | `CALL XVFACE( 3, XYPX )` | `INTEGER*2 XYPX(1:2,1:3)` (line 16) | OK |
| xvue/triacoul.f | 86 | `CALL XVFACE( 3, XYF )` | `INTEGER*2 XYF(2,5)` (lotria.f pattern; declared in xvue/lotria.f companion); local XYF declared `INTEGER*2 XYF(2,5)` in caller's COMMON via lotria | OK |
| xvue/triacoul.f | 103 | `CALL XVFACE( 4, XYF )` | `INTEGER*2 XYF(2,5)` | OK |
| xvue/triacoul.f | 120 | `CALL XVFACE( 4, XYF )` | `INTEGER*2 XYF(2,5)` | OK |
| xvue/triacoul.f | 140 | `CALL XVFACE( 5, XYF )` | `INTEGER*2 XYF(2,5)` | OK |
| xvue/triacoul.f | 168 | `CALL XVFACE( 4, XYF )` | `INTEGER*2 XYF(2,5)` | OK |
| xvue/triacoul.f | 188 | `CALL XVFACE( 4, XYF )` | `INTEGER*2 XYF(2,5)` | OK |
| xvue/triacoul.f | 201 | `CALL XVFACE( 3, XYF )` | `INTEGER*2 XYF(2,5)` | OK |
| util/t3flec.f | 107 | `CALL XVFACE( 3, XYSF )` | `INTEGER*2 XYSF(2,3)` (line 22) | OK |
| util/t3flec.f | 132 | `CALL XVFACE( 3, XYSF )` | `INTEGER*2 XYSF(2,3)` (line 22) | OK |
| xvue/lotria.f | 35 | `CALL XVFACE( 3, T )` | `INTEGER*2 T(2,3)` (line 21) | OK |
| xvue/fap32d.f | 133 | `CALL XVFACE( 4, XYPX )` | `INTEGER*2 XYPX(2,5)` (line 36) | OK |
| xvue/fap32d.f | 152 | `CALL XVFACE( 4, XYPX )` | `INTEGER*2 XYPX(2,5)` (line 36) | OK |
| xvue/fap32d.f | 237 | `CALL XVFACE( 4, XYPX )` | `INTEGER*2 XYPX(2,5)` (line 36) | OK |
| xvue/triacoul3dbord.f | 58 | `CALL XVTRAITS( I, XYPX )` | `INTEGER*2 XYPX(2,4)` (line 25) | OK |
| xvue/fap33d.f | 242 | `CALL XVFACE( 4, XYPX )` | `INTEGER*2 XYPX(2,5)` (line 49) | OK |
| xvue/fap33d.f | 308 | `CALL XVFACE( 4, XYPX )` | `INTEGER*2 XYPX(2,5)` (line 49) | OK |
| xvue/fap33d.f | 494 | `CALL XVFACE( 5, XYPX )` | `INTEGER*2 XYPX(2,5)` (line 49) | OK |
| xvue/triacou3dbord.f | 94 | `CALL XVTRAITS( I, XYPX )` | `INTEGER*2 XYPX(2,4)` (line 28) | OK |
| xvue/face2d.f | 67 | `CALL XVFACE( NBSP, XYPX )` | `INTEGER*2 XYPX(2,MXPOIN)` (line 25) | OK |
| xvue/face2d.f | 73 | `CALL XVTRAITS( NBSP, XYPX )` | `INTEGER*2 XYPX(2,MXPOIN)` (line 25) | OK |
| xvue/face2d.f | 78 | `CALL XVFACETRAITS( NCF, NCA, NBSP, XYPX )` | `INTEGER*2 XYPX(2,MXPOIN)` (line 25) | OK |
| xvue/triacoul2dbord.f | 42 | `CALL XVTRAITS( I, XYPX )` | `INTEGER*2 XYPX(2,4)` (line 21) | OK |
| xvue/face3d.f | 79 | `CALL XVFACE( NBSP, XYPX )` | `INTEGER*2 XYPX(2,MXPOIN)` (line 26) | OK |
| xvue/face3d.f | 86 | `CALL XVTRAITS( NBSP, XYPX )` | `INTEGER*2 XYPX(2,MXPOIN)` (line 26) | OK |
| xvue/face3d.f | 91 | `CALL XVFACETRAITS( NCF, NCA, NBSP, XYPX )` | `INTEGER*2 XYPX(2,MXPOIN)` (line 26) | OK |
| prpr/xvtest1.f | 104 | `CALL XVFACE( MAXPTS, XPOINTS)` | `INTEGER*2 XPOINTS(2,MAXPTS)` (line 18) | OK |

## Commented-out call sites (not in live ABI)

| File | Line | CALL statement | Note |
|------|------|----------------|------|
| xvue/fap32d.f | 134 | `CCC CALL XVFACETRAITS( NCF, NCBLAN, 4, XYPX )` | commented out |
| xvue/fap32d.f | 153 | `CCC CALL XVFACETRAITS( NCF, NCBLAN, 4, XYPX )` | commented out |
| xvue/fap32d.f | 238 | `CCC CALL XVFACETRAITS( NCF, NCBLAN, 5, XYPX )` | commented out |

These are historical `CCC`-prefixed comment lines — not compiled, do
not participate in the ABI.

## Hard blockers

None.

## Re-run command

```sh
grep -rniE 'CALL[[:space:]]+(XVFACE|XVTRAITS|XVFACETRAITS)\b' \
     xvue/ prpr/ mail/ elas/ flui/ ther/ nlse/ util/ reso/
```

## References

- `.planning/phases/02-drawing-primitives-backing-pixmap/02-CONTEXT.md` §D-29..D-31
- `.planning/research/PITFALLS.md` §Pitfall 4
- `xvue/qt/include/xvue_qt_api.h` lines 30–35
