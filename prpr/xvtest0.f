      PROGRAM XVTEST0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TESTER LE CYCLE D'OUVERTURE/FERMETURE + LES PRIMITIVES DE TRACE
C       (Phase 1 SHELL + Phase 2 DRAW-01..09 coverage driver)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : xvue-qt migration team   AVRIL 2026
C2345X7..............................................................012
C
C     Driver minimal pour la validation Phase 1 du backend Qt 6 :
C     chaque appel a XVINITGRAPHIQUE doit afficher une fenetre blanche
C     "MEFISTO" de 800x600 pixels, et chaque XVFERMER doit la detruire
C     sans toucher au QApplication. Deux cycles sont executes pour
C     verifier que la reouverture ne declenche pas l'assertion
C     "QApplication: there can only be one".
C
C     Entre le premier XVINITGRAPHIQUE et le premier XVFERMER, une
C     section "draw-coverage" (D-36, Phase 2) exerce chaque primitive
C     DRAW-01..09 et chaque style de pinceau (0/1/2). En Wave 0 les
C     corps des primitives sont encore des stubs warn-once, donc le
C     run produit des messages stdout et exit 0 sans dessiner.
C
      INTEGER*2 PTS(2,8)
      REAL*4    A1, A2
C
      PRINT *
      PRINT *,'==========================================='
      PRINT *,'Phase 1+2: cycle open/close + primitives'
      PRINT *,'==========================================='
C
      PRINT *,'[xvtest0] premier appel XVINITGRAPHIQUE'
      CALL XVINITGRAPHIQUE
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
      CALL XVVOIR
      CALL SLEEP(15)
C     clear mid-sequence, then draw with new pen style
      CALL EFFACER
      CALL XVTYPETRAIT(1)
      CALL XVTRAIT( 50, 450, 750, 450)
      CALL XVTYPETRAIT(2)
      CALL XVEPAISSEUR(2)
      CALL XVTRAIT( 50, 470, 750, 470)
      CALL XVVOIR
C     -- end draw-coverage section --------------------------------
C
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
      PRINT *,'[xvtest0] second appel XVINITGRAPHIQUE (reopen)'
      CALL XVINITGRAPHIQUE
      CALL SLEEP(3)
      PRINT *,'[xvtest0] second appel XVFERMER'
      CALL XVFERMER
C
      PRINT *
      PRINT *,'[xvtest0] OK — cycle open/close/open/close + draws'
      STOP
      END
