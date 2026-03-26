      SUBROUTINE FR3P2P2( CoefM, CoefV, XYZEF, TP2P2, TDP2DP2,  AE )
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
C TP2P2  : INTEGRALE SUR L'EF REFERENCE des  pi  pj dX avec p=P2
C TDP2DP2: INTEGRALE SUR L'EF REFERENCE des dpi dpj dX avec p=P2
C
C SORTIE :
C --------
C AE     : MATRICE ELEMENTAIRE STOCKEE SYMETRIQUE (10*11/2 coefficients)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray Avril 2011
C23456---------------------------------------------------------------012
!$    USE OMP_LIB
      IMPLICIT NONE
      include"./incl/langue.inc"
      REAL               XYZEF(10,3)
      DOUBLE PRECISION   CoefM, CoefV, TP2P2(10,10), TDP2DP2(3,10,3,10),
     %                   AE(55)
      DOUBLE PRECISION   DELTAe, DFM1(3,3), DF(3,3), TDFDF(3,3)
      DOUBLE PRECISION   XD, YD, ZD, C1, C2, S
      INTEGER            I, J, K, L, KE
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
      IF( DELTAe .LE. 0D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            print *,'FR3P2P2: VOLUME du EF=',DELTAe,' <=0 !'
            print *,'FR3P2P2: EF NON PRIS EN COMPTE'
         ELSE
            print *,'FR3P2P2: FE VOLUME=',DELTAe,' <=0 !'
            print *,'FR3P2P2: FE NOT COMPUTED'
      ENDIF
         GOTO 9999
      ENDIF
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
            S = TP2P2(I,J) * C1
C
C           CoefV Integrale tdP2i/dxk dP2j/dxl * deltae
            DO K=1,3
               DO L=1,3
                  S = S + TDP2DP2(K,I,L,J) * TDFDF(K,L)
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
 9999 RETURN
      END
