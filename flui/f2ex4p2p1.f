      SUBROUTINE F2EX4P2P1( CoGrPr, NBNOVI, XYZNOE, NONOTR, NONOSO, 
     %                      NBSOM,  PRESS0,   VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TRIANGLE TAYLOR HOOD
C -----    P2 CONTINU POUR LES 2 COMPOSANTES DE LA VITESSE
C          P1 CONTINU POUR LA PRESSION
C          CoGrPr Integrale  Div VitP2       Pression  dX
C         (=-CoGrPr Integrale    VitP2  Grad Pression  dX
C            avec oubli de la normale sur les faces ...
C           -CoGrPr Integrale    VitP2.n     Pression  dgamma )
C ENTREES:
C --------
C CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C NBNOVI : NOMBRE DE NOEUDS DU MAILLAGE
C XYZNOE : 3 COORDONNEES DES NBNOVI NOEUDS DE LA TRIANGULATION
C NONOTR : NONOTR(I) NO GLOBAL DU I-EME NOEUD DU TRIANGLE I=1,...,6
C NONOSO : NONOSO(I) = NUMERO DE SOMMET DE 1 A NBSOM DU I-EME NOEUD GLOBAL

C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
C          = NOMBRE DE DL DE LA PRESSION
C PRESS0 : DL DE LA PRESSION AUX SOMMETS DE LA TETRAEDRISATION
C
C SORTIE :
C --------
C VE     : VE(i,k) = Integrale t(div P2i) CoGrPr P1(tn) dX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Paris & St Pierre du Perray Mai 2011
C23456---------------------------------------------------------------012
!$    USE OMP_LIB
      IMPLICIT   NONE
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NBSOM, NBNOVI, NONOTR(6), NONOSO(NBNOVI)
      INTEGER           I, J, K, L,NSJ
      DOUBLE PRECISION  PRESS0(NBSOM), CoGrPr, VE(6,2)
      DOUBLE PRECISION  DFM1(2,2), S
C
      include "./incl/p1dp22d.inc"
C     INTEGRALES P1 DP2 SUR LE TRIANGLE de REFERENCE
C     DOUBLE PRECISION  P1DP22D(3,2,6)
C     P1DP22D(i,k,j) = integrale P1i DP2j/dxk dx dy sur TRIANGLE UNITE
C
C     CALCUL DE LA MATRICE JACOBIENNE INVERSE DFM1
C     LES 4 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DELTAe
C     CAR COMPENSE PAR LA MULTIPLICATION PAR DELTA DE L'INTEGRATION
      I = NONOTR(1)
      DFM1(1,1) = ( XYZNOE(2,NONOTR(3)) - XYZNOE(2,I) ) * CoGrPr
      DFM1(2,1) = ( XYZNOE(1,I) - XYZNOE(1,NONOTR(3)) ) * CoGrPr
C
      DFM1(1,2) = ( XYZNOE(2,I) - XYZNOE(2,NONOTR(2)) ) * CoGrPr
      DFM1(2,2) = ( XYZNOE(1,NONOTR(2)) - XYZNOE(1,I) ) * CoGrPr
C
C     L'INTEGRALE  CoGrPr  t(Div Vitesse) Pression dx
C     CoGrPr integrale t( DF-1 ligne1  DP2
C                        +DF-1 ligne2  DP2 ) P1 dX {dl Pression}
C     ----------------------------------------------------------
      DO I=1,6
C
         DO K=1,2
C
            S = 0D0
            DO J=1,3
C
C              NO GLOBAL DU SOMMET J DU TRIANGLE
               NSJ = NONOSO( NONOTR(J) )
C
C              P1DP22D(j,l,i) = integrale P1j dP2i/dxl dx dy
               DO L=1,2
                  S = S + DFM1(K,L) * P1DP22D(J,L,I) * PRESS0( NSJ )
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
