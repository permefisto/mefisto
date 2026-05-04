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
! How to materialize the golden (Plan 06 Task 3 procedure):
!   $ cd /tmp
!   $ cp $MEFISTO/xvue/qt/tests/golden/scene01_driver.f .
!   $ gfortran -I$MEFISTO/incl -c scene01_driver.f
!   $ gfortran scene01_driver.o $MEFISTO/xvue/xvuelc.o \
!         $MEFISTO/xvue/*.o $MEFISTO/util/*.o \
!         -L/usr/X11R6/lib -lX11 -lXt -o scene01_x11
!   $ ./scene01_x11        # produces TEMPORAIRE.EPS in cwd
!   $ cp TEMPORAIRE.EPS $MEFISTO/xvue/qt/tests/golden/scene01.eps
!   $ cd $MEFISTO && git add xvue/qt/tests/golden/scene01.eps
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

!     Initialise an offscreen-equivalent X11 surface so emit calls have
!     valid mempx + ypixels. The xvinitgraphique_ entry sets up the X11
!     pixmap behind the scenes; for this golden producer we run headless
!     on Xvfb (the human procedure above wraps the binary in xvfb-run).
      CALL XVINITGRAPHIQUE()

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
      CALL XVCHARGEFONTE('Courier', 7, 12, 10, 0)
      CALL XVTEXTE('Hello PS', 8, 50, 50)
      CALL XVTEXTE('MEFISTO', 7, 100, 100)

!     End PostScript capture. handleLasops(0) flushes concat_ and closes
!     fpo_. TEMPORAIRE.EPS is now finalized in cwd.
      LASOPS = 0
      CALL XVPOSTSCRIPT(LASOPS)

      CALL XVFERMER()
      END
