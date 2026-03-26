      SUBROUTINE F3EX4P1BP1GRAD( XYZSOM, NONOTE, DELTAT, CoGrPr,
     %                           NBSOM,  PRESS, BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TETRAEDRE BREZZI FORTIN
C -----    BE = -DELTAT CoGrPr Integrale sur ef  tV(x) grad PRESS(x)  dx
C          P1+BULLE CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          P1       CONTINU POUR LA PRESSION
C ENTREES:
C --------
C XYZSOM : 3 COORDONNEES DES SOMMETS DE LA TETRAEDRISATION
C NONOTR : NONOTR(I) NO GLOBAL DU I-EME SOMMET DU TETRAEDRE I=1,...,4
C          NONOTR(5) = NBSOMMET + NUMERO DU TETRAEDRE
C DELTAT : PAS DE TEMPS POUR L'INTEGRATION EN TEMPS
C CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
C          = NOMBRE DE DL DE LA PRESSION
C PRESS  : DL DE LA PRESSION AUX SOMMETS DE LA TETRAEDRISATION
C          INTERPOLATION P1 aux SOMMETS du TETRAEDRE

C SORTIE :
C --------
C BE     : BE(i,k) = DL ELEMENTAIRE Integrale tP1Bi CoefP GRADk PRESS dX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Decembre 2010
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      REAL              XYZSOM(3,NBSOM)
      DOUBLE PRECISION  PRESS(NBSOM)
      INTEGER           NBSOM
      INTEGER           NONOTE(5), I, J, K
      DOUBLE PRECISION  BE(5,3)
      DOUBLE PRECISION  DELTAT, CoGrPr, COEF
      DOUBLE PRECISION  DFM1(3,3), DF(3,3), DFM1DLa(3,4)
      DOUBLE PRECISION  X1, Y1, Z1, S, PRJ
      INTRINSIC         ABS

      DOUBLE PRECISION  INTP1B(5)
C     INTP1B(i) = Integrale sur e chapeau P1Bi dx dy dz

C     MISE A ZERO GENERALE EXCEPTE LE PREMIER BLOC DIAGONAL
C     =====================================================
      DO K = 1,3
         DO I = 1,5
            BE(I,K) = 0D0
         ENDDO
      ENDDO

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

      K = NONOTE(4)
      DF(3,1) = XYZSOM(1,K) - X1
      DF(3,2) = XYZSOM(2,K) - Y1
      DF(3,3) = XYZSOM(3,K) - Z1

cccC     LE DETERMINANT DE DF
ccc      DELTAe = ABS( DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
ccc     %            + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
ccc     %            + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ))
C     LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL

C     LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DELTAe
      DFM1(1,1) = ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
      DFM1(2,1) = ( DF(2,3) * DF(3,1) - DF(3,3) * DF(2,1) )
      DFM1(3,1) = ( DF(2,1) * DF(3,2) - DF(3,1) * DF(2,2) )

      DFM1(1,2) = ( DF(1,3) * DF(3,2) - DF(1,2) * DF(3,3) )
      DFM1(2,2) = ( DF(1,1) * DF(3,3) - DF(1,3) * DF(3,1) )
      DFM1(3,2) = ( DF(1,2) * DF(3,1) - DF(1,1) * DF(3,2) )

      DFM1(1,3) = ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) )
      DFM1(2,3) = ( DF(1,3) * DF(2,1) - DF(2,3) * DF(1,1) )
      DFM1(3,3) = ( DF(1,1) * DF(2,2) - DF(2,1) * DF(1,2) )

C     [DFM1] [DLambda]
      DFM1DLa(1,1) = -DFM1(1,1) - DFM1(1,2) - DFM1(1,3)
      DFM1DLa(1,2) =  DFM1(1,1)
      DFM1DLa(1,3) =  DFM1(1,2)
      DFM1DLa(1,4) =  DFM1(1,3)

      DFM1DLa(2,1) = -DFM1(2,1) - DFM1(2,2) - DFM1(2,3)
      DFM1DLa(2,2) =  DFM1(2,1)
      DFM1DLa(2,3) =  DFM1(2,2)
      DFM1DLa(2,4) =  DFM1(2,3)

      DFM1DLa(3,1) = -DFM1(3,1) - DFM1(3,2) - DFM1(3,3)
      DFM1DLa(3,2) =  DFM1(3,1)
      DFM1DLa(3,3) =  DFM1(3,2)
      DFM1DLa(3,4) =  DFM1(3,3)

C     L'INTEGRALE COEF tP1B  Gradk PressionP1 dX
C     ------------------------------------------
      S = 73.D0 / 2520.D0
      INTP1B(1) = S
      INTP1B(2) = S
      INTP1B(3) = S
      INTP1B(4) = S
      INTP1B(5) = 16.D0 / 315.D0

      COEF = - DELTAT * CoGrPr
      DO I=1,5

         S = INTP1B(I) * COEF
         DO J=1,4
C           NONOTE(J) NO GLOBAL DU SOMMET J DU TETRAEDRE P1
            PRJ = PRESS( NONOTE(J) ) * S
            DO K=1,3
               BE(I,K) = BE(I,K) + PRJ * DFM1DLa(K,J)
            ENDDO
         ENDDO

      ENDDO
C
      RETURN
      END
