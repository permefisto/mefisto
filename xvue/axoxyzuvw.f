      SUBROUTINE AXOXYZUVW( AXYZUVW, XYZUVW )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES 6 COORDONNEES DANS LE REPERE XYZUVW
C -----    A PARTIR DES 6 COORDONNEES AXONOMETRIQUES
C
C ENTREE :
C --------
C AXYZUVW: LES 6 COORDONNEES DANS LE REPERE AXONOMETRIQUE
C
C SORTIE :
C --------
C XYZUVW : LES 6 COORDONNEES DANS LE REPERE INITIAL DES OBJETS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire JL LIONS UPMC Paris Decembre 2006
C2345X789............................................................012
      include"./incl/trvari.inc"
      REAL     AXYZUVW(6), XYZUVW(6)
C
C     PROTECTION POUR APPEL AVEC AXYZUVW = XYZUVW
      X = AXYZUVW(1)
      Y = AXYZUVW(2)
      Z = AXYZUVW(3)
      U = AXYZUVW(4)
      V = AXYZUVW(5)
      W = AXYZUVW(6)
C
      XYZUVW(1) = AXOMAT(1,1) * X
     %          + AXOMAT(1,2) * Y
     %          + AXOMAT(1,3) * Z
     %          + AXOMAT(1,4) * U
     %          + AXOMAT(1,5) * V
     %          + AXOMAT(1,6) * W + AXOPTV(1)
C
      XYZUVW(2) = AXOMAT(2,1) * X
     %          + AXOMAT(2,2) * Y
     %          + AXOMAT(2,3) * Z
     %          + AXOMAT(2,4) * U
     %          + AXOMAT(2,5) * V
     %          + AXOMAT(2,6) * W + AXOPTV(2)
C
      XYZUVW(3) = AXOMAT(3,1) * X
     %          + AXOMAT(3,2) * Y
     %          + AXOMAT(3,3) * Z
     %          + AXOMAT(3,4) * U
     %          + AXOMAT(3,5) * V
     %          + AXOMAT(3,6) * W + AXOPTV(3)
C
      XYZUVW(4) = AXOMAT(4,1) * X
     %          + AXOMAT(4,2) * Y
     %          + AXOMAT(4,3) * Z
     %          + AXOMAT(4,4) * U
     %          + AXOMAT(4,5) * V
     %          + AXOMAT(4,6) * W + AXOPTV(4)
C
      XYZUVW(5) = AXOMAT(5,1) * X
     %          + AXOMAT(5,2) * Y
     %          + AXOMAT(5,3) * Z
     %          + AXOMAT(5,4) * U
     %          + AXOMAT(5,5) * V
     %          + AXOMAT(5,6) * W + AXOPTV(5)
C
      XYZUVW(6) = AXOMAT(6,1) * X
     %          + AXOMAT(6,2) * Y
     %          + AXOMAT(6,3) * Z
     %          + AXOMAT(6,4) * U
     %          + AXOMAT(6,5) * V
     %          + AXOMAT(6,6) * W + AXOPTV(6)
      END
