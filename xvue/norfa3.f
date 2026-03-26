      SUBROUTINE NORFA3( P1, P2, P3, COORNO, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DES COMPOSANTES DU VECTEUR NORMAL DE NORME 1 A UNE FACE
C -----   TRIANGULAIRE DEFINIE PAR LES 3 POINTS P1 P2 P3

C PARAMETRES D ENTREE :
C ---------------------
C P1 P2 P3 : LES 3 FOIS 3 COORDONNEES DES SOMMETS DU TRIANGLE

C PARAMETRE RESULTAT :
C --------------------
C COORNO : LES 3 COMPOSANTES DU VECTEUR NORMAL DE NORME 1 A LA FACE
C IERR   : =0 SI VECTEUR NORMAL CORRECTEMENT CALCULE
C          >0 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS             MARS 1987
C ......................................................................
      REAL  P1(3), P2(3), P3(3), COORNO(3), COORSO(3,3)

C     LES 3 POINTS SONT RANGES DANS COORSO POUR APPELER LE SP NORFAC
      DO I=1,3
         COORSO(I,1) = P1(I)
         COORSO(I,2) = P2(I)
         COORSO(I,3) = P3(I)
      ENDDO

C     LE CALCUL DE LA NORMALE
      CALL NORFAC( 3, COORSO, COORNO, IERR )

      RETURN
      END
