      SUBROUTINE NLDATADS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  POUR UN PROBLEME THERMIQUE NON LINEAIRE ou NLSE
C -----  DETRUIRE DANS MCN LES 3 TABLEAUX
C        MNTHDL(NBPOLY,2), MNTHDLn(NBPOLY,2), MNTHDL0(NBPOLY,2)
C        Cf  ./incl/cthet.inc
C        CES TABLEAUX STOCKENT LES VALEURS DES DL DE L'ONDE SUR UN EF
C        AUX TEMPS tn+1, tn, t0
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray  Aout 2011
C23456---------------------------------------------------------------012
      include"./incl/cthet.inc"
C
      IF( MNTHDL .GT. 0 ) THEN
         CALL TNMCDS( 'REEL2', 3 * 2*NBDLMXEF, MNTHDL )
      ENDIF
C
C     REPARTITION
      MNTHDLn = 0
      MNTHDL0 = 0
C
      RETURN
      END
