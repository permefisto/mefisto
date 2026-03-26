      SUBROUTINE F3VEP2( CoefM, CoefV, XYZEF,  AE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE VISCOSITE Ae DE LA VITESSE P2
C -----    DU TETRAEDRE TAYLOR HOOD
C          Ae = Integrale CoefM P2 P2 + CoefV DP2 DP2 dX
C          INTEGRATION EN XYZ EXACTE POUR CoefM et CoefV SUPPOSES CONSTANTS
C
C ENTREES:	
C --------
C CoefM  : COEFFICIENT DE LA MATRICE DE MASSE      P2  P2
C CoefV  : COEFFICIENT DE LA MATRICE DE VISCOSITE DP2 DP2
C XYZEF  : 3 COORDONNEES DES 10 NOEUDS DU TETRAEDRE
C
C SORTIE :
C --------
C AE     : MATRICE ELEMENTAIRE STOCKEE SYMETRIQUE (10*11/2 coefficients)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray Avril 2011
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      REAL               XYZEF(10,3)
      DOUBLE PRECISION   CoefM, CoefV,  AE(55)
      DOUBLE PRECISION   DELTAe, DFM1(3,3), DF(3,3), TDFDF(3,3)
      DOUBLE PRECISION   XD, YD, ZD, C1, C2, S
      INTEGER            I, J, K, L, KE
C
C     INTEGRALES P2 P2 SUR LE TETRAEDRE P2 de REFERENCE
C     DOUBLE PRECISION  P2P23D(10,10)
      include"./incl/p2p23d.inc"
C
C     INTEGRALES DP2 DP2 SUR LE TETRAEDRE de REFERENCE
C     DP2DP2(i,k,j) = integrale DP2i DP2j/dxk dx dy dz sur TETRAEDRE UNITE
C     DOUBLE PRECISION   DP2DP23D(3,10,3,10)
      include"./incl/dp2dp23d.inc"
C
C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
      XD = XYZEF(1,1)
      YD = XYZEF(1,2)
      ZD = XYZEF(1,3)
      DF(1,1) = XYZEF(2,1) - XD
      DF(1,2) = XYZEF(2,2) - YD
      DF(1,3) = XYZEF(2,3) - ZD
C
      DF(2,1) = XYZEF(3,1) - XD
      DF(2,2) = XYZEF(3,2) - YD
      DF(2,3) = XYZEF(3,3) - ZD
C
      DF(3,1) = XYZEF(4,1) - XD
      DF(3,2) = XYZEF(4,2) - YD
      DF(3,3) = XYZEF(4,3) - ZD
C
C     [DF]-1
      CALL M33INV( DF, DELTAe, DFM1 )
C
C     t[DF]-1 [DF]-1  EST UNE MATRICE SYMETRIQUE
C     CALCUL DES COEFFICIENTS DE LA MATRICE DE VISCOSITE
C     tDF-1 * DF-1 * DELTAe * CoefV
      C2 = DELTAe * CoefV
      TDFDF(1,1) = ( DFM1(1,1) * DFM1(1,1)
     %             + DFM1(2,1) * DFM1(2,1)
     %             + DFM1(3,1) * DFM1(3,1) ) * C2
C
      TDFDF(2,1) = ( DFM1(1,2) * DFM1(1,1)
     %             + DFM1(2,2) * DFM1(2,1)
     %             + DFM1(3,2) * DFM1(3,1) ) * C2
      TDFDF(1,2) = TDFDF(2,1)
C
      TDFDF(3,1) = ( DFM1(1,3) * DFM1(1,1)
     %             + DFM1(2,3) * DFM1(2,1)
     %             + DFM1(3,3) * DFM1(3,1) ) * C2
      TDFDF(1,3) = TDFDF(3,1)
C
      TDFDF(2,2) = ( DFM1(1,2) * DFM1(1,2)
     %             + DFM1(2,2) * DFM1(2,2)
     %             + DFM1(3,2) * DFM1(3,2) ) * C2
C
      TDFDF(3,2) = ( DFM1(1,3) * DFM1(1,2)
     %             + DFM1(2,3) * DFM1(2,2)
     %             + DFM1(3,3) * DFM1(3,2) ) * C2
      TDFDF(2,3) = TDFDF(3,2)
C
      TDFDF(3,3) = ( DFM1(1,3) * DFM1(1,3)
     %             + DFM1(2,3) * DFM1(2,3)
     %             + DFM1(3,3) * DFM1(3,3) ) * C2
C
C     AFFECTATION DES COEFFICIENTS DE LA MATRICE TRIANGULAIRE INFERIEURE
C     LIGNE PAR LIGNE ET DE GAUCHE A DROITE JUSQU'A LA DIAGONALE
C
C     COEFFICIENT DE Integrale P2(i) P2(j) dX
      C1 = DELTAe * CoefM
      KE = 0
      DO I = 1, 10
         DO J = 1, I
C
C           CoefM Integrale P2i P2j * deltae
            S = P2P23D(I,J) * C1
C
C           CoefV Integrale tdP2i/dxk dP2j/dxl * deltae
            DO K=1,3
               DO L=1,3
                  S = S + DP2DP23D(K,I,L,J) * TDFDF(K,L)
               ENDDO
            ENDDO
C
C           COEFFICIENT Ae(I,J) de la MATRICE CoefM P2 P2 + CoefV DP2 DP2
            KE = KE + 1
            AE(KE) = S
C
         ENDDO
      ENDDO
C
      RETURN
      END
