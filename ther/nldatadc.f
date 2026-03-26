      SUBROUTINE NLDATADC( NBDLMX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  POUR UN PROBLEME THERMIQUE NON LINEAIRE ou NLSE
C -----  DECLARER DANS MCN LES 3 TABLEAUX
C        MNTHDL(NBPOLY,2), MNTHDLn(NBPOLY,2), MNTHDL0(NBPOLY,2)
C        Cf  ./incl/cthet.inc
C        CES TABLEAUX STOCKENT LES VALEURS DES DL DE L'ONDE SUR UN EF
C        AUX TEMPS tn+1, tn, t0
C
C ENTREE :
C --------
C NBDLMX : NOMBRE MAXIMAL DE DL POUR UN EF OU
C          NOMBRE MAXIMAL DE POLYNOMES DE BASE DE L'INTERPOLATION SUR UN EF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Octobre 2013
C23456---------------------------------------------------------------012
      include"./incl/cthet.inc"
C
      IF( NBDLMXEF .NE. NBDLMX ) THEN
         IF( NBDLMXEF .GT. 0 ) THEN
            CALL TNMCDS( 'REEL2', 3 * 2*NBDLMXEF, MNTHDL )
         ENDIF
      ENDIF

      IF( MNTHDL .LE. 0 ) THEN
         CALL TNMCDC( 'REEL2', 3 * 2*NBDLMX, MNTHDL )
         NBDLMXEF = NBDLMX
      ENDIF
C
C     REPARTITION
      MOREE2  = MOTVAR(6)
      MNTHDLn = MNTHDL  + MOREE2 * 2*NBDLMX
      MNTHDL0 = MNTHDLn + MOREE2 * 2*NBDLMX
C
      RETURN
      END
