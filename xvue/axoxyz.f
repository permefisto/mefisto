      SUBROUTINE AXOXYZ( AXYZ, XYZ )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL DES 3 COORDONNEES DANS LE REPERE XYZ
C -----  A PARTIR DES 3 COORDONNEES AXONOMETRIQUES
C
C ENTREE :
C --------
C AXYZ   : LES 3 COORDONNEES DANS LE REPERE AXONOMETRIQUE
C
C SORTIE :
C --------
C XYZ    : LES 3 COORDONNEES DANS LE REPERE INITIAL DES OBJETS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS         MAI 1994
C2345X789............................................................012
      include"./incl/trvari.inc"
      REAL     AXYZ(3), XYZ(3)
C
C     PROTECTION POUR APPEL AVEC AXYZ = XYZ
      X = AXYZ(1)
      Y = AXYZ(2)
      Z = AXYZ(3)
C
      XYZ(1) = AXOMAT(1,1) * X
     %       + AXOMAT(1,2) * Y
     %       + AXOMAT(1,3) * Z + AXOPTV(1)
      XYZ(2) = AXOMAT(2,1) * X
     %       + AXOMAT(2,2) * Y
     %       + AXOMAT(2,3) * Z + AXOPTV(2)
      XYZ(3) = AXOMAT(3,1) * X
     %       + AXOMAT(3,2) * Y
     %       + AXOMAT(3,3) * Z + AXOPTV(3)

      RETURN
      END
