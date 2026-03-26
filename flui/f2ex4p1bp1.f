      SUBROUTINE F2EX4P1BP1( XYZSOM, NONOEF, DELTAT, CoGrPr,
     %                       NBSOM,  PRESS0, VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TRIANGLE BREZZI FORTIN
C -----    VEk=  DELTAT CoGrPr Integrale sur ef t(d/dxk P1B) Pression dx
C          (avec oubli de la normale sur les faces ...
C               -DELTAT CoGrPr Integrale  VitP1B.n   Pression  dgamma )
C          P1+BULLE CONTINU POUR LES 2 COMPOSANTES DE LA VITESSE
C          P1       CONTINU POUR LA PRESSION
C
C ENTREES:
C --------
C XYZSOM : 3 COORDONNEES DES SOMMETS DE LA TRIANGULATION
C          ICI LA COMPOSANTE Z=XYZSOM(3,.) N'EST PAS UTILISEE
C NONOEF : NONOEF(I) NO GLOBAL DU I-EME SOMMET DU TRIANGLE I=1,...,3
C          NONOEF(4) = NB SOMMET + NUMERO DU TRIANGLE
C DELTAT : PAS DE TEMPS
C CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
C          = NOMBRE DE DL DE LA PRESSION
C PRESS0 : DL DE LA PRESSION AUX SOMMETS DE LA TRIANGULATION
C
C SORTIES:
C --------
C VE     : DL ELEMENTAIRES DE LA VITESSE A L'INSTANT temps+deltat
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Novembre 2012
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      REAL              XYZSOM(3,NBSOM)
      INTEGER           NBSOM, NONOEF(4), I, J, K, L, NS1, NS2, NS3
      DOUBLE PRECISION  PRESS0(NBSOM), DELTAT, CoGrPr, VE(4,2)
      DOUBLE PRECISION  DFM1(2,2), S, COEF

C     DOUBLE PRECISION  DP1BLa2d(4,2,3)
C     DP1BLa2d(i,k,j)   = integrale dPi/dxk Lambdaj dX
      include"./incl/dp1bla2d.inc"

!     COEFFICIENT DE L'INTEGRALE
      COEF = DELTAT * CoGrPr

!     NUMERO DES 3 SOMMETS DU TRIANGLE NUELEM
      NS1 = NONOEF(1)
      NS2 = NONOEF(2)
      NS3 = NONOEF(3)

!     CALCUL DE LA MATRICE JACOBIENNE INVERSE
      DFM1(1,1) = ( XYZSOM(2,NS3) - XYZSOM(2,NS1) ) * COEF
      DFM1(2,1) = ( XYZSOM(1,NS1) - XYZSOM(1,NS3) ) * COEF

      DFM1(1,2) = ( XYZSOM(2,NS1) - XYZSOM(2,NS2) ) * COEF
      DFM1(2,2) = ( XYZSOM(1,NS2) - XYZSOM(1,NS1) ) * COEF
C
C     L'INTEGRALE COEF t(d/dxk P1B)  Pression dX
C     ------------------------------------------
      DO I=1,4
         DO K=1,2
            S = 0D0
            DO J=1,2
               DO L=1,3
                  S = S +DFM1(K,J) *DP1BLa2d(I,J,L) *PRESS0( NONOEF(L) )
               ENDDO
            ENDDO
            VE(I,K) = S
         ENDDO
      ENDDO
C
      RETURN
      END
