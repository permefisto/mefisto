      SUBROUTINE BF3DMV( Rho, XYZEF, NONOEF, MDVG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES COEFFICIENTS ELEMENTAIRES ET LES ASSEMBLER
C -----    DANS LA MATRICE GLOBALE DIAGONALE D'UNE COMPOSANTE DE LA VITESSE
C
C ENTREES:
C --------
C Rho    : DENSITE VOLUMIQUE DE MASSE DU FLUIDE
C XYZEF  : 3 COORDONNEES DES 4 SOMMETS DE L'EF
C NONOEF : NUMERO DES 5 NOEUDS DU TETRAEDRE BREZZI-FORTIN
C
C SORTIE :
C --------
C MDVG   : MATRICE GLOBALE DIAGONALE D'UNE COMPOSANTE DE LA VITESSE
C          Integrale Rho P1B P1B dX avec une formule d'integration
C          numerique aux 4 sommets + Barycentre => Matrice Diagonale
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Decembre 2010
C23456---------------------------------------------------------------012
      REAL               XYZEF(4,3)
      DOUBLE PRECISION   Rho, DELTA, DETM33, MDVG(1:*), S, X1, Y1, Z1
      INTEGER            NONOEF(5)
      INTRINSIC          ABS
C
C     CALCUL DU JACOBIEN DE CE TETRAEDRE BREZZI-FORTIN
      X1 = XYZEF(1,1)
      Y1 = XYZEF(1,2)
      Z1 = XYZEF(1,3)
      DELTA=ABS( DETM33( XYZEF(2,1)-X1, XYZEF(3,1)-X1, XYZEF(4,1)-X1,
     %                   XYZEF(2,2)-Y1, XYZEF(3,2)-Y1, XYZEF(4,2)-Y1,
     %                   XYZEF(2,3)-Z1, XYZEF(3,3)-Z1, XYZEF(4,3)-Z1) )
C
C     CALCUL DE LA MATRICE ELEMENTAIRE DE MASSE AVEC POINTS
C     D'INTEGRATION AUX 4 SOMMETS + BARYCENTRE
C     ASSEMBLAGE DANS LA MATRICE DIAGONALE MDVG(1:NBNOVI)
      S = Rho * DELTA / 120.D0
      DO I=1,4
         NS = NONOEF(I)
         MDVG( NS ) = MDVG(NS ) + S
      ENDDO
C
C     DE MEME AU BARYCENTRE DU TETRAEDRE AVEC LE POIDS CORRESPONDANT
      NS = NONOEF(5)
      MDVG( NS ) = MDVG(NS ) + S * 16.D0
C
      RETURN
      END
