      SUBROUTINE F3EX4P2P1( CoefPr, NBNOVI, XYZNOE, NONOTE, NONOSO,
     %                      NBSOM,  PRESS0,  VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TETRAEDRE TAYLOR-HOOD
C -----    P2 CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          P1 CONTINU POUR LA PRESSION
C          CoefPr Integrale  Div VitP2       Pression  dX
C       (=-CoefPr Integrale      VitP2  Grad Pression  dX
C          avec oubli de la normale sur les faces ...
C          -CoefPr Integrale    VitP2.n     Pression  dgamma )
C ENTREES:
C --------
C CoefPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C NBNOVI : NOMBRE DE NOEUDS DU MAILLAGE
C XYZNOE : 3 COORDONNEES DES NBNOVI NOEUDS DE LA TETRAEDRISATION
C NONOTE : NONOTE(I) NO GLOBAL DU I-EME NOEUD DU TETRAEDRE I=1,...,10
C NONOSO : NONOSO(I) = NUMERO DE SOMMET DE 1 A NBSOM DU I-EME NOEUD GLOBAL
C
C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
C          = NOMBRE DE DL DE LA PRESSION
C PRESS0 : DL DE LA PRESSION AUX SOMMETS DE LA TETRAEDRISATION
C          INTERPOLATION P1 aux SOMMETS du TETRAEDRE
C
C SORTIE :
C --------
C VE     : VE(i,k) = CoefPr Integrale  dVitP2i/dxk  Press(tn)  dX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Paris &St Pierre du Perray Juin 2011
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           LECTEU, IMPRIM, NUNITE
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NBSOM, NBNOVI, NONOTE(10), NONOSO(NBNOVI)
      INTEGER           I, J, K, L, NSJ
      DOUBLE PRECISION  PRESS0(NBSOM), VE(10,3)
      DOUBLE PRECISION  CoefPr, DFM1(3,3), DF(3,3), DELTAe,
     %                  X1, Y1, Z1, S
C
      include "./incl/p1dp23d.inc"
C     INTEGRALES P1 DP2 SUR LE TETRAEDRE de REFERENCE
C     P1DP23D(i,k,j) = integrale P1i DP2j/dxk dx dy dz sur TETRAEDRE UNITE
C     DOUBLE PRECISION  P1DP23D(4,3,10)
C
C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
C     =================================================
      I = NONOTE(1)
      X1 = XYZNOE(1,I)
      Y1 = XYZNOE(2,I)
      Z1 = XYZNOE(3,I)
C
      DO J=1,3
         I = NONOTE(J+1)
         DF(J,1) = XYZNOE(1,I) - X1
         DF(J,2) = XYZNOE(2,I) - Y1
         DF(J,3) = XYZNOE(3,I) - Z1
      ENDDO
C
C     LE DETERMINANT DE DF
      DELTAe = ABS( DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %            + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %            + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ))
      IF( DELTAe .LE. 0D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*)'F3EX4P2P1: ATTENTION EF de VOLUME*6=',
     %                      DELTAe,' NON PRIS EN COMPTE'
          ELSE
            WRITE(IMPRIM,*)'F3EX4P2P1: ATTENTION FE of VOLUME*6=',
     %                      DELTAe,' is NOT COMPUTED'
          ENDIF
          GOTO 9999
      ENDIF
C     ICI LE TETRAEDRE EST DE VOLUME NON NUL
C
C     LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DELTAe
      DFM1(1,1) = ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) ) * CoefPr
      DFM1(2,1) = ( DF(2,3) * DF(3,1) - DF(3,3) * DF(2,1) ) * CoefPr
      DFM1(3,1) = ( DF(2,1) * DF(3,2) - DF(3,1) * DF(2,2) ) * CoefPr
C
      DFM1(1,2) = ( DF(1,3) * DF(3,2) - DF(1,2) * DF(3,3) ) * CoefPr
      DFM1(2,2) = ( DF(1,1) * DF(3,3) - DF(1,3) * DF(3,1) ) * CoefPr
      DFM1(3,2) = ( DF(1,2) * DF(3,1) - DF(1,1) * DF(3,2) ) * CoefPr
C
      DFM1(1,3) = ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ) * CoefPr
      DFM1(2,3) = ( DF(1,3) * DF(2,1) - DF(2,3) * DF(1,1) ) * CoefPr
      DFM1(3,3) = ( DF(1,1) * DF(2,2) - DF(2,1) * DF(1,2) ) * CoefPr
C
C     L'INTEGRALE  CoefPr (Div Vitesse) Pression dx
C     CoefPr integrale  ( DF-1 ligne1  DP2
C                        +DF-1 ligne2  DP2
C                        +DF-1 ligne3  DP2 ) P1 dX {dl Pression}
C     -----------------------------------------------------------------
      DO I=1,10
C
         DO K=1,3
C
            S = 0D0
            DO J=1,4
C
C              NO GLOBAL DU SOMMET J DU TETRAEDRE
               NSJ = NONOSO( NONOTE(J) )
C
C              P1DP23D(j,l,i) = integrale P1j dP2i/dxl dx dy
               DO L=1,3
                  S = S + DFM1(K,L) * P1DP23D(J,L,I) * PRESS0( NSJ )
               ENDDO
C
            ENDDO
C
C           LE COEFFICIENT VE(I,K)
            VE(I,K) = S
C
         ENDDO
C
      ENDDO
C
 9999 RETURN
      END
