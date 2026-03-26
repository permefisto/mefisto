      SUBROUTINE F3EX4P1BP1( XYZSOM, NONOTE, DELTAT, CoGrPr,
     %                       NBSOM,  PRESS0, VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TETRAEDRE BREZZI FORTIN
C -----    VE = -DELTAT CoGrPr Integrale sur ef  tV(x) grad P(tn+1,x) dx
C          P1+BULLE CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          P1       CONTINU POUR LA PRESSION
C          VE =  DELTAT CoGrPr Integrale sur ef t(d/dxk P1B) Pression dx
C                 avec oubli de la normale sur les faces ...
C               -DELTAT CoGrPr Integrale  VitP1B.n   Pression  dgamma )
C ENTREES:
C --------
C XYZSOM : 3 COORDONNEES DES SOMMETS DE LA TETRAEDRISATION
C NONOTR : NONOTR(I) NO GLOBAL DU I-EME SOMMET DU TETRAEDRE I=1,...,4
C          NONOTR(5) = NBSOMMET + NUMERO DU TETRAEDRE
C DELTAT : PAS DE TEMPS POUR L'INTEGRATION EN TEMPS
C CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
C          = NOMBRE DE DL DE LA PRESSION
C PRESS0 : DL DE LA PRESSION AUX SOMMETS DE LA TETRAEDRISATION
C          INTERPOLATION P1 aux SOMMETS du TETRAEDRE
C
C SORTIES:
C --------
C VE     : VE(i,k) = DL ELEMENTAIRE Integrale tP1Bi CoefP GRADk P(tn) dX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Novembre 2012
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      REAL              XYZSOM(3,NBSOM)
      DOUBLE PRECISION  PRESS0(NBSOM)
      INTEGER           NBSOM
      INTEGER           NONOTE(5), I, J, K, L
      DOUBLE PRECISION  VE(5,3)
      DOUBLE PRECISION  DELTAT, CoGrPr, COEF
      DOUBLE PRECISION  DFM1(3,3), DF(3,3)
      DOUBLE PRECISION  X1, Y1, Z1, S

C     DOUBLE PRECISION  DP1Bla3d(5,3,4)
C     DP1Bla3d(i,k,j)   = integrale dPi/dxk Lambdaj dX
      include"./incl/dp1bla3d.inc"

!     COEFFICIENT DE L'INTEGRALE
      COEF = DELTAT * CoGrPr
C
C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
C     =================================================
      I = NONOTE(1)
      X1 = XYZSOM(1,I)
      Y1 = XYZSOM(2,I)
      Z1 = XYZSOM(3,I)

      J = NONOTE(2)
      DF(1,1) = XYZSOM(1,J) - X1
      DF(1,2) = XYZSOM(2,J) - Y1
      DF(1,3) = XYZSOM(3,J) - Z1

      K = NONOTE(3)
      DF(2,1) = XYZSOM(1,K) - X1
      DF(2,2) = XYZSOM(2,K) - Y1
      DF(2,3) = XYZSOM(3,K) - Z1

      L = NONOTE(4)
      DF(3,1) = XYZSOM(1,L) - X1
      DF(3,2) = XYZSOM(2,L) - Y1
      DF(3,3) = XYZSOM(3,L) - Z1
C
cccC     LE DETERMINANT DE DF
ccc      DELTAe = ABS( DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
ccc     %            + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
ccc     %            + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ))
C     LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
C
C     LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DELTAe
      DFM1(1,1) = ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) ) * COEF
      DFM1(2,1) = ( DF(2,3) * DF(3,1) - DF(3,3) * DF(2,1) ) * COEF
      DFM1(3,1) = ( DF(2,1) * DF(3,2) - DF(3,1) * DF(2,2) ) * COEF
C
      DFM1(1,2) = ( DF(1,3) * DF(3,2) - DF(1,2) * DF(3,3) ) * COEF
      DFM1(2,2) = ( DF(1,1) * DF(3,3) - DF(1,3) * DF(3,1) ) * COEF
      DFM1(3,2) = ( DF(1,2) * DF(3,1) - DF(1,1) * DF(3,2) ) * COEF
C
      DFM1(1,3) = ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ) * COEF
      DFM1(2,3) = ( DF(1,3) * DF(2,1) - DF(2,3) * DF(1,1) ) * COEF
      DFM1(3,3) = ( DF(1,1) * DF(2,2) - DF(2,1) * DF(1,2) ) * COEF
C
C     L'INTEGRALE COEF t(d/dxk P1B)  Pression dX
C     ------------------------------------------
      DO I=1,5
         DO K=1,3
            S = 0D0
            DO J=1,3
               DO L=1,4
                  S = S +DFM1(K,J) *DP1Bla3d(I,J,L) *PRESS0( NONOTE(L) )
               ENDDO
            ENDDO
            VE(I,K) = S
         ENDDO
      ENDDO
C
      RETURN
      END
