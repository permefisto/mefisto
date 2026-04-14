      PROGRAM XVTEST0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TESTER LE CYCLE D'OUVERTURE/FERMETURE + LES PRIMITIVES DE TRACE
C       (Phase 1 SHELL + Phase 2 DRAW-01..09 coverage driver)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : xvue-qt migration team   AVRIL 2026
C2345X7..............................................................012
C
C     Driver minimal pour la validation Phase 1 du backend Qt 6 :
C     chaque appel a XVOUVRIR doit afficher une fenetre blanche
C     "MEFISTO" et chaque XVFERMER doit la detruire sans toucher au
C     QApplication. Deux cycles sont executes pour verifier que la
C     reouverture ne declenche pas l'assertion "QApplication: there
C     can only be one".
C
C     XVOUVRIR est l'entree officielle des deux backends (legacy X11
C     et Qt 6) : elle appelle XTINIT (ouverture display/QApplication)
C     puis XVINIT (ouverture fenetre + pixmap + couleurs + polices).
C     Appeler xvinitgraphique_ directement ne cree la fenetre et le
C     pixmap que dans la backend Qt ; sur le backend X11 la fenetre
C     et le pixmap ne sont crees que par XVINFO appele depuis XVINIT,
C     ce qui imposait autrefois l'usage de XVOUVRIR.
C
C     Entre le premier XVOUVRIR et le premier XVFERMER, une section
C     "draw-coverage" (D-36, Phase 2) exerce chaque primitive
C     DRAW-01..09 et chaque style de pinceau (0/1/2). En Wave 0 les
C     corps des primitives sont encore des stubs warn-once, donc le
C     run produit des messages stdout et exit 0 sans dessiner.
C
      INTEGER*2 PTS(2,8)
      REAL*4    A1, A2
      INTEGER   NOFONT0, NOFONT, NPLACA, NPHACA
      INTEGER   ICOL, IX1, IY1, IX2, IY2
      INTEGER   NLEN, NPXLA, NPXHA, NOPALC, NBCELLS, I
      REAL      PROUGE_CUS(10), PVERT_CUS(10), PBLEU_CUS(10)
      CHARACTER*32 LABEL
      CHARACTER*32 SCENE
C
      PRINT *
      PRINT *,'==========================================='
      PRINT *,'Phase 1+2: cycle open/close + primitives'
      PRINT *,'==========================================='
C
C     ====================================================
C     PHASE 4 COVERAGE — pixmap save/restore round-trip
C     (runs FIRST per 04-RESEARCH Open Question 2, so the
C     base scene is a clean black canvas; skipped entirely
C     when MEFISTO_XVTEST0_SCENE is blank to preserve Phase
C     1/2/3 coverage in the default driver run.)
C     ====================================================
      SCENE = ' '
      CALL GETENV('MEFISTO_XVTEST0_SCENE', SCENE)
      IF (SCENE .NE. ' ') THEN
        PRINT *,'[xvtest0] PHASE 4 scene = ', SCENE
        CALL XVOUVRIR
        IF (SCENE .EQ. 'P4_BG') THEN
C         Background-only: XVINITGRAPHIQUE + EFFACER, nothing else.
          CALL EFFACER
          CALL MEMPXFENETRE
          CALL XVVOIR
          CALL XVFERMER
          STOP
        END IF
C
C       Default base scene for the other branches.
        CALL P4BASE
C
        IF (SCENE .EQ. 'P4_CTRL') THEN
C         Mesh-only control — no save/restore, no overlays.
          CALL MEMPXFENETRE
          CALL XVVOIR
          CALL XVFERMER
          STOP
        ELSE IF (SCENE .EQ. 'P4_SAVERESTORE') THEN
          CALL SAUVEFENETRE
          CALL XVCOULEUR(6)
          CALL XVBORDRECTANGLE(200, 200, 100, 100)
          CALL RESTAUREFENETRE
          CALL MEMPXFENETRE
          CALL XVVOIR
          CALL XVFERMER
          STOP
        ELSE IF (SCENE .EQ. 'P4_MEMPX_SAVERESTORE') THEN
          CALL SAUVEMEMPX
          CALL XVCOULEUR(6)
          CALL XVBORDRECTANGLE(200, 200, 100, 100)
          CALL RESTAUREMEMPX
          CALL MEMPXFENETRE
          CALL XVVOIR
          CALL XVFERMER
          STOP
        ELSE IF (SCENE .EQ. 'P4_EFFACEMEMPX') THEN
          CALL EFFACEMEMPX
          CALL MEMPXFENETRE
          CALL XVVOIR
          CALL XVFERMER
          STOP
        ELSE IF (SCENE .EQ. 'P4_FENETREMEMPX') THEN
C         P4BASE already drew the base. On Qt both FENETREMEMPX and
C         MEMPXFENETRE are intentional no-ops (Phase 4 D-04), so a
C         capture produced by sandwiching the final XVVOIR between
C         them MUST equal the P4_CTRL capture pixel-for-pixel.
C         [Deviation Rule 1] The planner's original action inserted
C         an extra XVTRAIT here; that would have produced a scene
C         that could never equal P4_CTRL. Removing it is the only
C         way to actually prove the no-ops are transparent.
          CALL FENETREMEMPX
          CALL MEMPXFENETRE
          CALL XVVOIR
          CALL XVFERMER
          STOP
        END IF
C       Unknown scene → fall through to legacy path for safety.
        CALL XVFERMER
      END IF
C
      PRINT *,'[xvtest0] premier appel XVOUVRIR'
      CALL XVOUVRIR
C
C     -- Phase 2 draw-coverage section (D-36) --------------------
C     pen style 0 (solid), width default
      CALL XVTYPETRAIT(0)
      CALL XVEPAISSEUR(1)
C     one line (DRAW-02)
      CALL XVTRAIT(100, 100, 300, 100)
C     polyline >= 3 pts (DRAW-02)
      PTS(1,1) = 120
      PTS(2,1) = 150
      PTS(1,2) = 220
      PTS(2,2) = 200
      PTS(1,3) = 320
      PTS(2,3) = 150
      CALL XVTRAITS(3, PTS)
C     filled polygon >= 4 pts (DRAW-03)
      PTS(1,1) = 400
      PTS(2,1) = 100
      PTS(1,2) = 500
      PTS(2,2) = 150
      PTS(1,3) = 480
      PTS(2,3) = 250
      PTS(1,4) = 400
      PTS(2,4) = 220
      CALL XVFACE(4, PTS)
C     rectangles (DRAW-04) — all four symbols
      CALL XVBORDRECTANGLE( 50, 300, 100,  60)
      CALL XVRECTANGLE    (170, 300, 100,  60)
      CALL XVFRECTANGLE   (290, 300, 100,  60)
      CALL XVFBORDRECTANGLE(410, 300, 100, 60)
C     ellipse arcs (DRAW-05) — REAL*4 angles per corrected ABI
      A1 =  0.0
      A2 = 90.0
      CALL XVARCELLIPSE    (600, 150,  40,  40, A1, A2)
      A1 = -45.0
      A2 = 180.0
      CALL XVBORDARCELLIPSE(600, 150,  60,  60, A1, A2)
C     flush the pre-effacer scene and hold for visual inspection:
C     Checks 1 (geometry), 2 (pen styles), 3 (antialiasing), 5 (resize).
C     MEMPXFENETRE copies the off-screen pixmap to fenetre_mef — xvvoir_
C     itself only calls XRaiseWindow+XFlush since 1999 (window-manager
C     blocking bug) and does NOT copy the pixmap. Without MEMPXFENETRE
C     the drawings stay invisible in mempx. (cf. xvue/xvuelc.c:2476.)
      CALL MEMPXFENETRE
      CALL XVVOIR
      CALL SLEEP(15)
C     clear mid-sequence, then draw with new pen style
      CALL EFFACER
      CALL XVTYPETRAIT(1)
      CALL XVTRAIT( 50, 450, 750, 450)
      CALL XVTYPETRAIT(2)
      CALL XVEPAISSEUR(2)
      CALL XVTRAIT( 50, 470, 750, 470)
      CALL MEMPXFENETRE
      CALL XVVOIR
C     -- end draw-coverage section --------------------------------
C
C     -- begin TEXT coverage section (Phase 3, D-24) ---------------
C     Exercises TEXT-01..TEXT-05 Phase 3 entry points.
C
C     (a) Load font index 2 (12pt DejaVu Sans Mono).
      NOFONT0 = 0
      NOFONT  = 2
      CALL XVCHARGEFONTE(NOFONT0, NOFONT, NPLACA, NPHACA)
C
C     (b) Draw 8 labeled colored lines, one per imposed-default.
      DO 100 ICOL = 0, 7
        CALL XVCOULEUR(ICOL)
        IX1 = 50
        IY1 = 100 + ICOL*25
        IX2 = 250
        IY2 = IY1
        CALL XVTRAIT(IX1, IY1, IX2, IY2)
        WRITE(LABEL,'(A,I1)') 'COLOR-', ICOL
        NLEN = 7
        CALL XVTEXTE(LABEL, NLEN, 270, IY1)
  100 CONTINUE
C
C     (c) Exercise XVACTIVERVB with a ramp at indices 0..8.
      NOPALC  = 0
      NBCELLS = 9
      DO 110 I = 1, 9
        PROUGE_CUS(I) = 0.1 * I
        PVERT_CUS (I) = 1.0 - 0.1 * I
        PBLEU_CUS (I) = 0.5
  110 CONTINUE
      CALL XVACTIVERVB(NOPALC, NBCELLS, PROUGE_CUS, PVERT_CUS,
     +                 PBLEU_CUS)
      CALL XVCOULEUR(8)
      CALL XVTRAIT(50, 320, 250, 320)
      CALL XVTEXTE('RVB-CUSTOM', 10, 270, 320)
C
C     (d) Exercise XVNBPIXELTEXTE + bounding box.
      CALL XVCOULEUR(1)
      CALL XVNBPIXELTEXTE('MESURE-XVNBPIXELTEXTE', 21, NPXLA,NPXHA)
      CALL XVTEXTE('MESURE-XVNBPIXELTEXTE', 21, 50, 360)
      CALL XVBORDRECTANGLE(50, 360 - NPXHA, NPXLA, NPXHA)
C     -- end TEXT coverage section ----------------------------------
C
      CALL MEMPXFENETRE
      CALL XVVOIR
C     Hold the window on screen long enough to be visually verified
C     (SHELL-01, SHELL-06). XVINITGRAPHIQUE has already pumped the
C     event loop until the window is exposed, so the window is mapped
C     and painted before SLEEP starts; the X server keeps the Expose
C     content on screen while the Fortran main blocks in SLEEP.
C     Second hold shows the post-effacer state only — verifies
C     Check 4 (effacer cleared the pre-effacer scene) and DRAW-06
C     pen style variety on the dashed lines.
      CALL SLEEP(10)
      PRINT *,'[xvtest0] premier appel XVFERMER'
      CALL XVFERMER
C
      PRINT *,'[xvtest0] second appel XVOUVRIR (reopen)'
      CALL XVOUVRIR
      CALL SLEEP(3)
      PRINT *,'[xvtest0] second appel XVFERMER'
      CALL XVFERMER
C
      PRINT *
      PRINT *,'[xvtest0] OK — cycle open/close/open/close + draws'
      STOP
      END
C
      SUBROUTINE P4BASE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : Phase 4 — base scene for round-trip tests.
C       Draws a 4x4 colored checker (via XVFACE) and a horizontal line.
C       Deterministic: same input → same backing pixels.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER*2 PTS(2,4)
      INTEGER   IR, IC, X0, Y0, W, H, ICOL
      W = 60
      H = 60
      DO 200 IR = 0, 3
        DO 210 IC = 0, 3
          X0 = 100 + IC * W
          Y0 = 100 + IR * H
          ICOL = 1 + MOD(IR + IC, 7)
          CALL XVCOULEUR(ICOL)
          PTS(1,1) = X0
          PTS(2,1) = Y0
          PTS(1,2) = X0 + W
          PTS(2,2) = Y0
          PTS(1,3) = X0 + W
          PTS(2,3) = Y0 + H
          PTS(1,4) = X0
          PTS(2,4) = Y0 + H
          CALL XVFACE(4, PTS)
  210   CONTINUE
  200 CONTINUE
      CALL XVCOULEUR(7)
      CALL XVTRAIT(80, 400, 380, 400)
      RETURN
      END
