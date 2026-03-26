      SUBROUTINE FM3P2P2( Rho, X, TP2P2,  AE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE MASSE Ae DE LA VITESSE P2
C -----    DU TETRAEDRE TAYLOR HOOD   Ae = Integrale P2 P2 DX
C          INTEGRATION EN XYZ EXACTE POUR Rho SUPPOSE CONSTANT
C
C ENTREES:	
C --------
C Rho    : COEFFICIENT DE LA MATRICE DE MASSE P2 P2
C X      : 3 COORDONNEES DES 10 NOEUDS DU TETRAEDRE
C TP2P2  : INTEGRALE SUR L'EF REFERENCE des pi pj dX avec p=P2
C
C SORTIE :
C --------
C AE     : MATRICE ELEMENTAIRE STOCKEE SYMETRIQUE (10*11/2 coefficients)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray Avril 2011
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      REAL               X(10,3)
      DOUBLE PRECISION   Rho, TP2P2(10,10), AE(55)
      DOUBLE PRECISION   ABS, DETM33, X1, Y1, Z1, DELTAe
      INTEGER            I, J, KE
C
C     CALCUL DU JACOBIEN DE Fe AVEC COMPOSANTES POLYNOMES DE DEGRE 1
      X1 = X(1,1)
      Y1 = X(1,2)
      Z1 = X(1,3)
      DELTAe = ABS( DETM33( X(2,1)-X1, X(3,1)-X1, X(4,1)-X1,
     %                      X(2,2)-Y1, X(3,2)-Y1, X(4,2)-Y1,
     %                      X(2,3)-Z1, X(3,3)-Z1, X(4,3)-Z1 ) )
C
C     PRISE EN COMPTE DU COEFFICIENT Rho
      DELTAe = DELTAe * Rho
C
C     AFFECTATION DES COEFFICIENTS DE LA MATRICE TRIANGULAIRE INFERIEURE
C     LIGNE PAR LIGNE ET DE GAUCHE A DROITE JUSQU'A LA DIAGONALE
      KE = 0
      DO I = 1, 10
         DO J = 1, I
C
C           LE COEFFICIENT Ae(I,J) de la MATRICE de MASSE P2 P2
            KE = KE + 1
C
C           Integrale Pi Pj * delta * Rho
            AE( KE ) = TP2P2(I,J) * DELTAe
C
         ENDDO
      ENDDO
C
      RETURN
      END
