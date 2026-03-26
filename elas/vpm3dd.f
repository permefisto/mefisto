      SUBROUTINE  VPM3DD( A, VP, DIRP )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DES VALEURS ET VECTEURS PROPRES D-UNE MATRICE 2X2
C ----- SYMETRIQUE EN DOUBLE PRECISION
C
C ENTREE :
C --------
C A      : LES COEFFICIENTS DE LA MATRICE 3x3 SYMETRIQUE
C
C SORTIES:
C --------
C VP     : LES 3 VALEURS PROPRES
C DIRP   : COMPOSANTES DES 3 VECTEURS PROPRES
C          DIRP(1:3,K) EST LE k-EME VECTEUR PROPRE ORTHONORMALISE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN  PERRONNET  ANALYSE NUMERIQUE UPMC PARIS 6   MAI 1998
C ....................................................................
      DOUBLE PRECISION A(3,3),VP(3),DIRP(3,3)
C
C     APPEL DE LA METHODE DE JACOBI (NUMERICAL RECIPES)
      CALL JACOBI( A, 3, 3, VP, DIRP, NBROTA )
      RETURN
      END
