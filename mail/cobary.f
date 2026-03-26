         SUBROUTINE COBARY( PXYZ, SXYZ, LAMBDA1, LAMBDA2, LAMBDA3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES COORDONNEES BARYCENTRIQUES D'UN POINT
C -----    DANS UN TRIANGLE NON DEGENERE
C ENTREES :
C ---------
C PXYZ    : COORDONNEES DU POINT
C SXYZ    : COORDONNEES DES SOMMETS DU TRIANGLE
C
C SORTIES :
C ---------
C LAMBDA  : COORDONNEES BARYCENTRIQUES DU POINT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY  ANALYSE NUMERIQUE UPMC PARIS        FEVRIER 1989
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      COMMON / EPSSSS / EPZERO,EPSXYZ
      REAL              PXYZ(2),SXYZ(2,3),LAMBDA1,LAMBDA2,LAMBDA3
C
C     RESOLUTION DU SYSTEME LINEAIRE
C
C     A * LAMBDA(1) + B * LAMBDA(2) = E
C     C * LAMBDA(1) + D * LAMBDA(2) = F
C
      A = SXYZ(1,1) - SXYZ(1,3)
      B = SXYZ(1,2) - SXYZ(1,3)
      C = SXYZ(2,1) - SXYZ(2,3)
      D = SXYZ(2,2) - SXYZ(2,3)
      E =  PXYZ(1)  - SXYZ(1,3)
      F =  PXYZ(2)  - SXYZ(2,3)
C
      DETER = A * D - B * C
      IF ( ABS(DETER) .LT. EPSXYZ ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15) , '(E15.6)') DETER
         KERR(1) = 'COBARY:TRIANGLE DE SURFACE'//
     %              KERR(MXLGER)(1:15)
         DETER = 1.
      END IF
C
      DETER = 1. / DETER
      LAMBDA1 = (   D * E - B * F ) * DETER
      LAMBDA2 = ( - C * E + A * F ) * DETER
      LAMBDA3 = 1. - LAMBDA1 - LAMBDA2
C
      RETURN
      END
