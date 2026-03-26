      SUBROUTINE XYZAXO( XYZ, AXYZ )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL DES 3 COORDONNEES AXONOMETRIQUES
C -----  A PARTIR DES 3 COORDONNEES DU POINT XYZ
C
C ENTREE :
C --------
C XYZ    : LES 3 COORDONNEES DANS LE REPERE INITIAL
C
C SORTIE :
C --------
C AXYZ   : LES 3 COORDONNEES DANS LE REPERE AXONOMETRIQUE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS         MAI 1994
C2345X789............................................................012
      include"./incl/trvari.inc"
      REAL     XYZ(3), AXYZ(3)
C
C     PROTECTION POUR APPEL AVEC XYZ = AXYZ
      X = XYZ(1) - AXOPTV(1)
      Y = XYZ(2) - AXOPTV(2)
      Z = XYZ(3) - AXOPTV(3)
C
      AXYZ(1) = AXOMAT(1,1) * X + AXOMAT(2,1) * Y + AXOMAT(3,1) * Z
      AXYZ(2) = AXOMAT(1,2) * X + AXOMAT(2,2) * Y + AXOMAT(3,2) * Z
      AXYZ(3) = AXOMAT(1,3) * X + AXOMAT(2,3) * Y + AXOMAT(3,3) * Z

      RETURN
      END
