      SUBROUTINE F3EX4P2P1GRAD( XYZNOE, NONOTE, NONOSO, DELTAT, CoGrPr,
     %                          NBNOVI, NBSOM,  PRESP1,   BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL DU SECOND MEMBRE DU TETRAEDRE TAYLOR YOUNG
C -----  BEk = -DELTAT CoGrPr Integrale sur ef  tV(x) gradk PRESP1(x) dx
C        P2 CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE k=1,3
C        P1 CONTINU POUR LA PRESSION

C ENTREES:
C --------
C XYZNOE : XYZNOE(3,NBNOVI)  3 COORDONNEES DES NOEUDS DU MAILLAGE P2
C NONOTE : NONOTE(I) NO GLOBAL DU I-EME NOEUD DU TETRAEDRE P2 I=1,...,10
C NONOSO : NONOSO(I) = NUMERO DE SOMMET DE 1 A NBSOM DU I-EME NOEUD GLOBAL
C DELTAT : PAS DE TEMPS
C CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C NBNOVI : NOMBRE DE NOEUDS VITESSE P2
C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION P2
C          = NOMBRE DE DL DE LA PRESSION
C PRESP1 : DL DE LA PRESSION P1 AUX SOMMETS DE LA TETRAEDRISATION P2
C          INTERPOLATION P1 aux SOMMETS du TETRAEDRE

C SORTIE :
C --------
C BE     : BE(i,k) = DL ELEMENTAIRE Integrale tP2i COEF GRADk PRESP1 dX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray     Avril 2011
C23456---------------------------------------------------------------012
      IMPLICIT   NONE
      REAL              XYZNOE(3,NBNOVI)
      DOUBLE PRECISION  PRESP1(NBSOM)
      INTEGER           NBSOM, NBNOVI
      INTEGER           NONOTE(10), NONOSO(*), I, J, K, NSJ
      DOUBLE PRECISION  BE(10,3)
      DOUBLE PRECISION  DELTAT, CoGrPr, DELTAe, DFM1(3,3), DF(3,3),
     %                  DFM1DLa(3,4), X1, Y1, Z1, S, COEF

      DOUBLE PRECISION  INTP2(10)
C     INTP2(i) = Integrale sur e chapeau Pi dx dy dz =
C               INTP2(1:4)=-1D0/120D0  INTP2(5:10)=1D0/30D0 
      DATA INTP2/ 4 * -8.3333333333333332D-3,
     %            6 *  3.3333333333333333D-2 /

C     MISE A ZERO de BE
      DO K = 1, 3
         DO I = 1, 10
            BE(I,K) = 0D0
         ENDDO
      ENDDO

C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
C     =================================================
      I = NONOTE(1)
      X1 = XYZNOE(1,I)
      Y1 = XYZNOE(2,I)
      Z1 = XYZNOE(3,I)

      DO J=1,3
         I = NONOTE(J+1)
         DF(J,1) = XYZNOE(1,I) - X1
         DF(J,2) = XYZNOE(2,I) - Y1
         DF(J,3) = XYZNOE(3,I) - Z1
      ENDDO

C     LE DETERMINANT DE DF
      DELTAe = ABS( DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %            + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %            + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ))
C     LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL

C     LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DELTAe
      DFM1(1,1) = ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
      DFM1(2,1) = ( DF(2,3) * DF(3,1) - DF(3,3) * DF(2,1) )
      DFM1(3,1) = ( DF(2,1) * DF(3,2) - DF(3,1) * DF(2,2) )
C
      DFM1(1,2) = ( DF(1,3) * DF(3,2) - DF(1,2) * DF(3,3) )
      DFM1(2,2) = ( DF(1,1) * DF(3,3) - DF(1,3) * DF(3,1) )
      DFM1(3,2) = ( DF(1,2) * DF(3,1) - DF(1,1) * DF(3,2) )
C
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

C     L'INTEGRALE COEF tP2  Gradk PressionP1 dX
C     -----------------------------------------
      COEF = - DELTAT * CoGrPr
      DO I = 1, 10

         S = INTP2(I) * COEF
         DO J = 1, 4
C           NO GLOBAL DU SOMMET J DU TETRAEDRE P2
            NSJ = NONOSO( NONOTE(J) )
            DO K = 1, 3
               BE(I,K) = BE(I,K) + S * DFM1DLa(K,J) * PRESP1(NSJ)
            ENDDO
         ENDDO

      ENDDO

      RETURN
      END
