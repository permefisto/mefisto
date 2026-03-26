      SUBROUTINE F3EX4P2P2( XYZEF,  NONOTE, CoefPr,
     %                      NBNOVI, PRESS0,  VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TETRAEDRE TAYLOR-HOOD
C -----    P2 CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          P2 CONTINU POUR LA PRESSION
C          CoefPr Integrale  Div VitP2       Pression  dX
C         (=-CoefPr Integrale    VitP2  Grad Pression  dX
C            avec oubli de la normale sur les faces ...)
C           -CoefPr Integrale    VitP2.n     Pression  dgamma )
C ENTREES:
C --------
C XYZEF  : 3 COORDONNEES DES 10 NOEUDS DU TETRAEDRE
C NONOTE : NONOTE(I) NO GLOBAL DU I-EME NOEUD DU TETRAEDRE I=1,...,10
C
C CoefPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C NBNOVI : NOMBRE DE NOEUDS DE LA TETRAEDRISATION
C          = NOMBRE DE DL DE LA PRESSION
C PRESS0 : DL DE LA PRESSION AUX SOMMETS DE LA TETRAEDRISATION
C          INTERPOLATION P2 aux SOMMETS du TETRAEDRE
C
C SORTIE :
C --------
C VE     : VE(i,k) = CoefPr Integrale  dVitP2i/dxk  Press(tn)  dX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Paris &St Pierre du Perray Juin 2012
C23456---------------------------------------------------------------012
      IMPLICIT   NONE
      REAL              XYZEF(10,3)
      INTEGER           NBNOVI, NONOTE(10), I, J, K, L
      DOUBLE PRECISION  PRESS0(NBNOVI), VE(10,3)
      DOUBLE PRECISION  CoefPr, DFM1(3,3), DF(3,3),
     %                  X1, Y1, Z1, S
C
      include "./incl/p2dp23d.inc"
C     INTEGRALES P2 DP2 SUR LE TETRAEDRE de REFERENCE
C     P2DP23D(i,k,j) = integrale P2i DP2j/dxk dx dy dz sur TETRAEDRE UNITE
C     DOUBLE PRECISION  P2DP23D(10,3,10)
C
C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
C     =================================================
      X1 = XYZEF(1,1)
      Y1 = XYZEF(1,2)
      Z1 = XYZEF(1,3)
      DF(1,1) = XYZEF(2,1) - X1
      DF(1,2) = XYZEF(2,2) - Y1
      DF(1,3) = XYZEF(2,3) - Z1
C
      DF(2,1) = XYZEF(3,1) - X1
      DF(2,2) = XYZEF(3,2) - Y1
      DF(2,3) = XYZEF(3,3) - Z1
C
      DF(3,1) = XYZEF(4,1) - X1
      DF(3,2) = XYZEF(4,2) - Y1
      DF(3,3) = XYZEF(4,3) - Z1
C
cccC  LE DETERMINANT DE DF
ccc   DELTAe = ABS( DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
ccc  %            + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
ccc  %            + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ))
cccC  LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
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
C                        +DF-1 ligne3  DP2 ) P2 dX {dl Pression}
C     -----------------------------------------------------------------
      DO I=1,10
C
         DO K=1,3
C
            S = 0D0
            DO J=1,10
C
C              P2DP23D(j,l,i) = integrale P2j dP2i/dxl dx dy
               DO L=1,3
                  S = S + DFM1(K,L) * P2DP23D(J,L,I) * PRESS0(NONOTE(J))
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
      RETURN
      END
