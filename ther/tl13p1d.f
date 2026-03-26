      SUBROUTINE TL13P1D( EXPO,   X, NBVECT, TEMPEF,
     %                    VOLTET, INTEGP )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA NORME L1 DES NBVECT TEMPERATURES**EXPO
C -----    TETRAEDRE LAGRANGE DE DEGRE 1
C
C ENTREES:
C --------
C EXPO   : EXPOSANT DE LA TEMPERATURE DE NORME A CALCULER
C X      : COORDONNEES DES 4 SOMMETS DU TETRAEDRE
C TEMPEF : NBVECT TEMPERATURES AUX NBPOLY NOEUDS DE L'EF
C
C SORTIES:
C --------
C VOLTET : VOLUME DU TETRAEDRE
C INTEGP : Som   Omegal  (û**EXPO)(Sl) Delta(Sl)
C         l=1...,4
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET  LJLL UPMC & SAINT PIERRE DU PERRAY   MAI 2009
C23456---------------------------------------------------------------012
      REAL              X(4,3)
      DOUBLE PRECISION  EXPO, TEMPEF(4,NBVECT), INTEGP(NBVECT), VOLTET
      DOUBLE PRECISION  X12, Y12, Z12, X13, Y13, Z13, X14, Y14, Z14
      DOUBLE PRECISION  D, E, F
C
      X12 = X(1,1) - X(2,1)
      Y12 = X(1,2) - X(2,2)
      Z12 = X(1,3) - X(2,3)
C
      X13 = X(1,1) - X(3,1)
      Y13 = X(1,2) - X(3,2)
      Z13 = X(1,3) - X(3,3)
C
      X14 = X(1,1) - X(4,1)
      Y14 = X(1,2) - X(4,2)
      Z14 = X(1,3) - X(4,3)
C
      D = Y13 * Z14 - Y14 * Z13
      E = Y14 * Z12 - Y12 * Z14
      F = Y12 * Z13 - Y13 * Z12
C
      VOLTET = - ( X12 * D + X13 * E + X14 * F ) / 6D0
C
      DO 30 N=1, NBVECT
C
C        NORME L1 DE U**EXPO INTEGREE NUMERIQUEMENT AUX 4 SOMMETS
         INTEGP(N) = INTEGP(N)
     %             + VOLTET * ( TEMPEF(1,N)**EXPO
     %                        + TEMPEF(2,N)**EXPO
     %                        + TEMPEF(3,N)**EXPO
     %                        + TEMPEF(4,N)**EXPO ) * 0.25D0
C
 30   CONTINUE
C
      RETURN
      END
