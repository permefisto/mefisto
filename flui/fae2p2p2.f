      SUBROUTINE FAE2P2P2( Rho, Coef, NONOEF, NBNOVI, XYZNOE,  AE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE VISCOSITE Ae DE LA VITESSE P2
C -----    DU TRIANGLE TAYLOR HOOD
C          Ae = Integrale Rho P2 P2 + Coef DP2 DP2 dX
C          INTEGRATION EN XYZ EXACTE POUR Rho Coef SUPPOSES CONSTANTS
C
C ENTREES:	
C --------
C Rho    : COEFFICIENT DE LA MATRICE DE MASSE      P2  P2
C Coef   : COEFFICIENT DE LA MATRICE DE VISCOSITE DP2 DP2
C NONOEF : NUMERO GLOBAL DES 6 NOEUDS DU TRIANGLE
C NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et MILIEUX DES ARETES DES EF
C XYZNOE : XYZNOE(3,NBNOVI) 3 COORDONNEES DES NOEUDS DU MAILLAGE
C
C SORTIE :
C --------
C AE     : MATRICE ELEMENTAIRE STOCKEE SYMETRIQUE (6*7/2 coefficients)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC&Saint Pierre du Perray Octobre 2012
C23456---------------------------------------------------------------012
!$    USE OMP_LIB
      IMPLICIT NONE
      include"./incl/langue.inc"
      include"./incl/p2p22d.inc"
      include"./incl/dp2dp22d.inc"

      REAL               XYZNOE(3,NBNOVI)
      DOUBLE PRECISION   Rho, Coef,
     %                   AE(21), ABS, C1, D, S,
     %                   DELTAe, X21, X31, Y21, Y31, TDFDF(2,2)
      INTEGER            NONOEF(6), NBNOVI, I, J, K, L, KE
C
C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
      I = NONOEF(1)
      X21 = XYZNOE(1,NONOEF(2)) - XYZNOE(1,I)
      X31 = XYZNOE(1,NONOEF(3)) - XYZNOE(1,I)
C
      Y21 = XYZNOE(2,NONOEF(2)) - XYZNOE(2,I)
      Y31 = XYZNOE(2,NONOEF(3)) - XYZNOE(2,I)
C
C     CALCUL DU JACOBIEN DE Fe AVEC COMPOSANTES POLYNOMES DE DEGRE 1
      DELTAe = ABS( X21*Y31 - X31*Y21 )
C
      IF( DELTAe .LE. 0D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            print *,'FAE2P2P2: SURFACE EF=',DELTAe,' <=0 !'
            print *,'FAE2P2P2: EF NON PRIS EN COMPTE'
         ELSE
            print *,'FAE2P2P2: FE SURFACE=',DELTAe,' <=0 !'
            print *,'FAE2P2P2: FE NOT COMPUTED'
         ENDIF
         GOTO 9999
      ENDIF
C
C     COEFFICIENT DE Integrale P2(i) P2(j) dX
      C1 = DELTAe * Rho
C
C     COEFFICIENTS DE LA MATRICE Coef tDF-1 * DF-1 / DELTAe
      D = Coef / DELTAe
      TDFDF(1,1) = D * (  X31 * X31 + Y31 * Y31 )
      TDFDF(1,2) = D * ( -X21 * X31 - Y21 * Y31 )
      TDFDF(2,1) = TDFDF(1,2)
      TDFDF(2,2) = D * (  X21 * X21 + Y21 * Y21 )
C
C     AFFECTATION DES COEFFICIENTS DE LA MATRICE TRIANGULAIRE INFERIEURE
C     LIGNE PAR LIGNE ET DE GAUCHE A DROITE JUSQU'A LA DIAGONALE
C     P2P22D   : INTEGRALE SUR L'EF REFERENCE des  pi  pj dX avec p=P2
C     DP2DP22D : INTEGRALE SUR L'EF REFERENCE des dpi dpj dX avec p=P2
      KE = 0
      DO I = 1, 6
         DO J = 1, I
C
C           Integrale P2i P2j * DELTAe * Rho
            S = P2P22D(I,J) * C1
C
C           Integrale DP2i DP2j * DELTAe * Coef
            DO K=1,2
               DO L=1,2
                  S = S + DP2DP22D(K,I,L,J) * TDFDF(K,L)
               ENDDO
            ENDDO
C
C           COEFFICIENT Ae(I,J) de la MATRICE Rho P2 P2 + Coef DP2 DP2
            KE = KE + 1
            AE(KE) = S
C
         ENDDO
      ENDDO
C
 9999 RETURN
      END
