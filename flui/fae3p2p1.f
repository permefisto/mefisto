      SUBROUTINE FAE3P2P2( CoefM, CoefV, NONOEF, NBNOVI, XYZNOE,  AE )
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
C NONOEF : NUMERO GLOBAL DES 10 NOEUDS DU TETRAEDRE
C NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et MILIEUX DES ARETES DES EF
C XYZNOE : XYZNOE(3,NBNOVI) 3 COORDONNEES DES NOEUDS DU MAILLAGE
C
C SORTIE :
C --------
C AE     : MATRICE ELEMENTAIRE STOCKEE SYMETRIQUE (10*11/2 coefficients)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC&Saint Pierre du Perray Octobre 2012
C23456---------------------------------------------------------------012
!$    USE OMP_LIB
      IMPLICIT NONE
      include"./incl/langue.inc"
      include"./incl/p2p23d.inc"
      include"./incl/dp2dp23d.inc"

      INTEGER            NONOEF(10), NBNOVI
      REAL               XYZNOE(3,NBNOVI)
      DOUBLE PRECISION   CoefM, CoefV, AE(55)
      DOUBLE PRECISION   DELTAe, DFM1(3,3), DF(3,3), TDFDF(3,3)
      DOUBLE PRECISION   XD, YD, ZD, C1, C2, S
      INTEGER            I, J, K, L, KE
C
C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
      I = NONOEF(1)
      XD = XYZNOE(1,I)
      YD = XYZNOE(2,I)
      ZD = XYZNOE(3,I)

      J = NONOEF(2)
      DF(1,1) = XYZNOE(1,J) - XD
      DF(1,2) = XYZNOE(2,J) - YD
      DF(1,3) = XYZNOE(3,J) - ZD
C
      K = NONOEF(3)
      DF(2,1) = XYZNOE(1,K) - XD
      DF(2,2) = XYZNOE(2,K) - YD
      DF(2,3) = XYZNOE(3,K) - ZD
C
      L = NONOEF(4)
      DF(3,1) = XYZNOE(1,L) - XD
      DF(3,2) = XYZNOE(2,L) - YD
      DF(3,3) = XYZNOE(3,L) - ZD
C
C     [DF]-1
      CALL M33INV( DF, DELTAe, DFM1 )
C
      IF( DELTAe .LE. 0D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            print *,'FAE3P2P2: VOLUME du EF=',DELTAe,' <=0 !'
            print *,'FAE3P2P2: EF NON PRIS EN COMPTE'
         ELSE
            print *,'FAE3P2P2: FE VOLUME=',DELTAe,' <=0 !'
            print *,'FAE3P2P2: FE NOT COMPUTED'
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
C     P2P23D  : INTEGRALE SUR L'EF REFERENCE des  pi  pj dX avec p=P2
C     DP2DP23D: INTEGRALE SUR L'EF REFERENCE des dpi dpj dX avec p=P2
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
 9999 RETURN
      END
