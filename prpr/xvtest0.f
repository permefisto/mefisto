      PROGRAM XVTEST0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TESTER LE CYCLE D'OUVERTURE/FERMETURE DE LA FENETRE QT 6 XVUE
C       (Phase 1 — SHELL-01, SHELL-02 : reopen sans assertion singleton)
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
      PRINT *
      PRINT *,'==========================================='
      PRINT *,'Phase 1: test du cycle open/close de XVUE-Qt'
      PRINT *,'==========================================='
C
      PRINT *,'[xvtest0] premier appel XVINITGRAPHIQUE'
      CALL XVINITGRAPHIQUE
C     Hold the window on screen long enough to be visually verified
C     (SHELL-01, SHELL-06). XVINITGRAPHIQUE has already pumped the
C     event loop until the window is exposed, so the window is mapped
C     and painted before SLEEP starts; the X server keeps the Expose
C     content on screen while the Fortran main blocks in SLEEP.
      CALL SLEEP(1)
      PRINT *,'[xvtest0] premier appel XVFERMER'
      CALL XVFERMER
C
      PRINT *,'[xvtest0] second appel XVINITGRAPHIQUE (reopen)'
      CALL XVINITGRAPHIQUE
      CALL SLEEP(1)
      PRINT *,'[xvtest0] second appel XVFERMER'
      CALL XVFERMER
C
      PRINT *
      PRINT *,'[xvtest0] OK — cycle open/close/open/close complet'
      STOP
      END
