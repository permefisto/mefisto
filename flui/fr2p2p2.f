      SUBROUTINE FR2P2P2( Rho, Coef, XYZEF, TP2P2, TDP2DP2,  AE )
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
C XYZEF  : 2 COORDONNEES DES 6 NOEUDS DU TRIANGLE
C TP2P2  : INTEGRALE SUR L'EF REFERENCE des  pi  pj dX avec p=P2
C TDP2DP2: INTEGRALE SUR L'EF REFERENCE des dpi dpj dX avec p=P2
C
C SORTIE :
C --------
C AE     : MATRICE ELEMENTAIRE STOCKEE SYMETRIQUE (6*7/2 coefficients)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray Avril 2011
C23456---------------------------------------------------------------012
!$    USE OMP_LIB
      IMPLICIT NONE
      include"./incl/langue.inc"
      REAL               XYZEF(6,2)
      DOUBLE PRECISION   Rho, Coef, TP2P2(6,6), TDP2DP2(2,6,2,6),
     %                   AE(21), ABS, C1, D, S,
     %                   DELTAe, X21, X31, Y21, Y31, TDFDF(2,2)
      INTEGER            I, J, K, L, KE
C
C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
      X21 = XYZEF(2,1) - XYZEF(1,1)
      X31 = XYZEF(3,1) - XYZEF(1,1)
C
      Y21 = XYZEF(2,2) - XYZEF(1,2)
      Y31 = XYZEF(3,2) - XYZEF(1,2)
C
C     CALCUL DU JACOBIEN DE Fe AVEC COMPOSANTES POLYNOMES DE DEGRE 1
      DELTAe = ABS( X21*Y31 - X31*Y21 )
C
      IF( DELTAe .LE. 0D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            print *,'FR2P2P2: SURFACE EF=',DELTAe,' <=0 !'
            print *,'FR2P2P2: EF NON PRIS EN COMPTE'
         ELSE
            print *,'FR2P2P2: FE SURFACE=',DELTAe,' <=0 !'
            print *,'FR2P2P2: FE NOT COMPUTED'
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
      KE = 0
      DO I = 1, 6
         DO J = 1, I
C
C           Integrale P2i P2j * DELTAe * Rho
            S = TP2P2(I,J) * C1
C
C           Integrale DP2i DP2j * DELTAe * Coef
            DO K=1,2
               DO L=1,2
                  S = S + TDP2DP2(K,I,L,J) * TDFDF(K,L)
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
