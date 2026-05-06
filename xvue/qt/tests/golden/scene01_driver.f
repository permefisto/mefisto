! Phase 7 Plan 06: byte-level golden producer for PsEmitter port verification.
!
! This driver is a deterministic, constant-input scene whose PostScript
! output (TEMPORAIRE.EPS) becomes the byte-level reference against which
! Plan 03's PsEmitter port is compared.
!
! BOOTSTRAP NOTE: this driver is compiled + run by the human in Plan 06
! Task 3 (the checkpoint:human-verify task). The resulting TEMPORAIRE.EPS
! is then committed as xvue/qt/tests/golden/scene01.eps. Until that
! happens, the byte-level golden test slot QSKIPs cleanly.
!
! Deterministic-input contract (per 07-RESEARCH.md "Validation Architecture"):
!   - No random numbers, no time-dependent calls, no env-var reads.
!   - All coordinates and parameters are integer / float literals.
!   - Canvas dims are explicitly set via xvinitgraphique-equivalent so
!     pyFlip(y) has a known ypixels at emit time.
!
! Scene contents (stable; if this changes the golden MUST be regenerated):
!   - 1 epaisseur change (PS opcode "epais")
!   - 1 typetrait change (PS opcode "typet")
!   - 1 fond / clear-erase setup (mode 100 reset path stays untouched here
!     because xvfond emits no PS bytes per Plan 03 SUMMARY)
!   - 5 line segments (PS opcode "S" — covers Y-flip, counb!=-1)
!   - 2 filled faces (PS opcode "F")
!   - 1 ellipse arc (PS opcode "el")
!   - 2 text strings (PS opcode "T")
!   - 1 chargefonte call (PS opcode "charge")
!
! Total emit-helper coverage: 8 of the 12 helpers (the remaining 4 —
! traits/facetraits/bordrectangle/rectangle/arcellipse — are covered by
! per-primitive tests in test_xvue_qt_postscript.cpp; the full-scene golden
! exercises the most-used ones).
!
! How to materialize the golden (must use v1.0-pre-retire git tag — main
! tree post-Phase-9 has X11 backend retired; xvuelc.c is gone):
!   $ git worktree add /tmp/mefisto-pre-retire v1.0-pre-retire
!   $ cd /tmp/mefisto-pre-retire
!   $ MEFISTO=$PWD bin/cbl_tout                # build legacy X11 binaries
!   $ cp $MEFISTO/xvue/qt/tests/golden/scene01_driver.f /tmp/
!   $ cd /tmp
!   $ gfortran -I/tmp/mefisto-pre-retire/incl -c scene01_driver.f
!   $ gfortran scene01_driver.o /tmp/mefisto-pre-retire/xvue/xvuelc.o \
!         /tmp/mefisto-pre-retire/xvue/lib /tmp/mefisto-pre-retire/util/lib \
!         -L/usr/X11R6/lib -lX11 -lXt -o scene01_x11
!   $ MEFISTO_BATCH_X11=1 MEFISTO_XVSOURIS_AUTOEXIT=1 \
!         xvfb-run --auto-servernum ./scene01_x11
!   $ cp TEMPORAIRE.EPS $MEFISTO_MAIN/xvue/qt/tests/golden/scene01.eps
!   $ cd $MEFISTO_MAIN && git add xvue/qt/tests/golden/scene01.eps
!   $ git worktree remove /tmp/mefisto-pre-retire
!
! Once committed, ctest -R '^xvue_qt_postscript_tests$' flips the
! PsEmitter_postscriptVerbatim_golden slot from QSKIP to PASS-required.

      PROGRAM SCENE01
      IMPLICIT NONE

!     Local arrays for face coordinates (3-vertex triangle).
      INTEGER NPTS
      PARAMETER (NPTS = 3)
      INTEGER PTS_X(NPTS), PTS_Y(NPTS)

      INTEGER LASOPS

!     Initialise the X11 backend.
!     XVINITGRAPHIQUE alone does NOT create the GC on X11 (per prpr/xvtest0.f
!     comments line 17-22: window+pixmap only created by XVOUVRIR -> XVINFO ->
!     XVINIT chain). Calling drawing primitives before XVOUVRIR SIGSEGVs at
!     XChangeGC. XVOUVRIR opens window AND creates GC. For batch use, set
!     MEFISTO_BATCH_X11=1 + MEFISTO_XVSOURIS_AUTOEXIT=1 in the environment
!     before invocation so xvsouris returns immediately.
      CALL XVOUVRIR

!     Begin PostScript capture. handleLasops(1) opens TEMPORAIRE.EPS.
      LASOPS = 1
      CALL XVPOSTSCRIPT(LASOPS)

!     1. Epaisseur change — PS "epais" opcode.
      CALL XVEPAISSEUR(3)

!     2. Typetrait change — PS "typet" opcode.
      CALL XVTYPETRAIT(2)

!     3. Five line segments. Coordinates are constant integers;
!     PsEmitter applies the Y-flip inside line(). Canvas is whatever
!     xvinitgraphique sets — the test compares byte streams produced by
!     the same backend on the same input, so the ypixels value is
!     internally consistent.
      CALL XVTRAIT( 10,  20,  30,  40)
      CALL XVTRAIT( 50,  60,  70,  80)
      CALL XVTRAIT( 90, 100, 110, 120)
      CALL XVTRAIT(130, 140, 150, 160)
      CALL XVTRAIT(170, 180, 190, 200)

!     4. Two filled faces (triangles).
      PTS_X(1) =  10
      PTS_Y(1) =  20
      PTS_X(2) =  30
      PTS_Y(2) =  40
      PTS_X(3) =  50
      PTS_Y(3) =  60
      CALL XVFACE(PTS_X, PTS_Y, NPTS)

      PTS_X(1) =  70
      PTS_Y(1) =  80
      PTS_X(2) =  90
      PTS_Y(2) = 100
      PTS_X(3) = 110
      PTS_Y(3) = 120
      CALL XVFACE(PTS_X, PTS_Y, NPTS)

!     5. Ellipse arc. xvuelc.c emits "el" opcode.
!     Args: x, y, rx, ry, angle_start, angle_extent.
      CALL XVBORDARCELLIPSE(100, 100, 50, 30, 0.0, 360.0)

!     6. Two text strings — exercises the "T" opcode and chargefonte.
!     XVCHARGEFONTE legacy signature is (nofont0, nofont, largpx, hautpx) —
!     4 INTEGERS. The earlier ('Courier', 7, 12, 10, 0) form was wrong arity
!     (Fortran-C string passing adds a hidden length, mismatching the C decl
!     of 4 INTs). nofont0=0 = default font slot; nofont=1 = font #1; the
!     legacy backend resolves the actual font name via xvinfo's namefonts[].
      CALL XVCHARGEFONTE(0, 1, 12, 10)
      CALL XVTEXTE('Hello PS', 8, 50, 50)
      CALL XVTEXTE('MEFISTO', 7, 100, 100)

!     End PostScript capture. handleLasops(0) flushes concat_ and closes
!     fpo_. TEMPORAIRE.EPS is now finalized in cwd.
      LASOPS = 0
      CALL XVPOSTSCRIPT(LASOPS)

      CALL XVFERMER()
      END
