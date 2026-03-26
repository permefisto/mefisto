      SUBROUTINE TL12P1D( EXPO,  X, NDSM, TEMPEF,   DELTA, INTEGP )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA NORME L1 DES NDSM TEMPERATURES**EXPO
C -----    TRIANGLE LAGRANGE DE DEGRE 1
C
C ENTREES:
C --------
C EXPO   : EXPOSANT DE LA TEMPERATURE DE NORME A CALCULER
C X      : COORDONNEES DES 3 SOMMETS DU TRIANGLE
C TEMPEF : NDSM TEMPERATURES AUX NBPOLY NOEUDS DE L'EF
C
C SORTIES:
C --------
C DELTA  : SURFACE DU TRIANGLE
C INTEGP : Som Omegal  (û**EXPO)(Sl) Delta(Sl)
C         l=1...,3
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET  LJLL UPMC & SAINT PIERRE DU PERRAY   MAI 2009
C23456---------------------------------------------------------------012
      REAL              X(3,2)
      DOUBLE PRECISION  EXPO, TEMPEF(3,NDSM), INTEGP(NDSM)
      DOUBLE PRECISION  X21, Y21, X13, Y13, X32, Y32, DELTA
C
C     SURFACE DE L'EF
      X21 = X(2,1) - X(1,1)
      X32 = X(3,1) - X(2,1)
      X13 = X(1,1) - X(3,1)
C
      Y21 = X(2,2) - X(1,2)
      Y32 = X(3,2) - X(2,2)
      Y13 = X(1,2) - X(3,2)
C
      DELTA = ABS( X13 * Y21 - X21 * Y13 ) / 2D0
C
      DO 30 N=1, NDSM
C
C        NORME L1 DE U**EXPO INTEGREE NUMERIQUEMENT AUX 3 SOMMETS
         INTEGP(N) = INTEGP(N)
     %             + DELTA * ( TEMPEF(1,N)**EXPO
     %                       + TEMPEF(2,N)**EXPO
     %                       + TEMPEF(3,N)**EXPO ) / 3D0
C
 30   CONTINUE
C
      RETURN
      END
