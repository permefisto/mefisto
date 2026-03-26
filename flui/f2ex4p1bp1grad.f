      SUBROUTINE F2EX4P1BP1GRAD( XYZSOM, NONOEF, DELTAT, CoGrPr,
     %                           NBSOM,  PRESP1,  BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL DU SECOND MEMBRE DU TRIANGLE BREZZI FORTIN
C -----  BEk = -DELTAT CoGrPr Integrale sur ef  tV(x) gradk PRESP1(x) dx
C        P1+BULLE CONTINU POUR LES 2 COMPOSANTES DE LA VITESSE k=1,2
C        P1       CONTINU POUR LA PRESSION

C ENTREES:
C --------
C XYZSOM : 3 COORDONNEES DES SOMMETS DE LA TRIANGULATION
C          ICI LA COMPOSANTE Z=XYZSOM(3,.) N'EST PAS UTILISEE
C NONOEF : NONOEF(I) NO GLOBAL DU I-EME SOMMET DU TRIANGLE I=1,...,3
C          NONOEF(4) = NB SOMMETS NBSOM + NUMERO DU TRIANGLE
C DELTAT : PAS DE TEMPS
C CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
C          = NOMBRE DE DL DE LA PRESSION
C PRESP1 : DL DE LA PRESSION AUX SOMMETS DE LA TRIANGULATION P1

C SORTIE :
C --------
C BE     : DL ELEMENTAIRES DE LA VITESSE A L'INSTANT temps+deltat
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Decembre 2010
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      REAL              XYZSOM(3,NBSOM)
      INTEGER           NBSOM, NONOEF(4), I, J, K
      DOUBLE PRECISION  PRESP1(NBSOM), DELTAT, CoGrPr, BE(4,2)
      DOUBLE PRECISION  DFM1DLa(2,3), S, COEF, D

      DOUBLE PRECISION  INTP1B(4)
C     INTP1B(i) = Integrale sur e chapeau P1Bi dx dy
C     INTP1B(1:3)=11/120, INTP1B(4)=9/40
      DATA  INTP1B / 0.0916666666666666667D0, 0.0916666666666666667D0,
     %               0.0916666666666666667D0, 0.225D0 /

C     DFM1 DLambda(2,3)
C     -----------------
      I = NONOEF(1)
      J = NONOEF(2)
      K = NONOEF(3)

      DFM1DLa(1,1) = XYZSOM(2,J) - XYZSOM(2,K)
      DFM1DLa(2,1) = XYZSOM(1,K) - XYZSOM(1,J)

      DFM1DLa(1,2) = XYZSOM(2,K) - XYZSOM(2,I)
      DFM1DLa(2,2) = XYZSOM(1,I) - XYZSOM(1,K)
C
      DFM1DLa(1,3) = XYZSOM(2,I) - XYZSOM(2,J)
      DFM1DLa(2,3) = XYZSOM(1,J) - XYZSOM(1,I)

C     - DELTAT * CoGrPr INTEGRALE tVitesse  Grad Pression dX
C     ------------------------------------------------------
      COEF = - DELTAT * CoGrPr
      DO I=1,4

C        - DELTAT * CoGrPr  Integrale tP1B dX
         D = COEF * INTP1B(I)

         DO K=1,2

            S = 0D0
            DO J=1,3
C              NONOEF(J) NO GLOBAL DU SOMMET J DU TRIANGLE
               S = S + DFM1DLa(K,J) * PRESP1( NONOEF(J) )
            ENDDO

C           LE COEFFICIENT BE(I,K) EST INITIALISE
            BE(I,K) = S * D

         ENDDO

      ENDDO
C
      RETURN
      END
