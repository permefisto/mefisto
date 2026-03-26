      SUBROUTINE XYZUVWAXO( XYZUVW, AXYZUVW )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES 6 COORDONNEES AXONOMETRIQUES
C -----    A PARTIR DES 6 COORDONNEES DU POINT XYZUVW
C
C ENTREE :
C --------
C XYZUVW : LES 6 COORDONNEES DANS LE REPERE INITIAL
C
C SORTIE :
C --------
C AXYZUVW: LES 6 COORDONNEES DANS LE REPERE AXONOMETRIQUE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire JL LIONS UPMC Paris Decembre 2006
C2345X789............................................................012
      include"./incl/trvari.inc"
      REAL     XYZUVW(6), AXYZUVW(6)
C
C     PROTECTION POUR APPEL AVEC XYZUVW = AXYZUVW
      X = XYZUVW(1) - AXOPTV(1)
      Y = XYZUVW(2) - AXOPTV(2)
      Z = XYZUVW(3) - AXOPTV(3)
      U = XYZUVW(4) - AXOPTV(4)
      V = XYZUVW(5) - AXOPTV(5)
      W = XYZUVW(6) - AXOPTV(6)
C
      AXYZUVW(1) = AXOMAT(1,1) * X + AXOMAT(2,1) * Y + AXOMAT(3,1) * Z
     %           + AXOMAT(4,1) * U + AXOMAT(5,1) * V + AXOMAT(6,1) * W
      AXYZUVW(2) = AXOMAT(1,2) * X + AXOMAT(2,2) * Y + AXOMAT(3,2) * Z
     %           + AXOMAT(4,2) * U + AXOMAT(5,2) * V + AXOMAT(6,2) * W
      AXYZUVW(3) = AXOMAT(1,3) * X + AXOMAT(2,3) * Y + AXOMAT(3,3) * Z
     %           + AXOMAT(4,3) * U + AXOMAT(5,3) * V + AXOMAT(6,3) * W
      AXYZUVW(4) = AXOMAT(1,4) * X + AXOMAT(2,4) * Y + AXOMAT(3,4) * Z
     %           + AXOMAT(4,4) * U + AXOMAT(5,4) * V + AXOMAT(6,4) * W
      AXYZUVW(5) = AXOMAT(1,5) * X + AXOMAT(2,5) * Y + AXOMAT(3,5) * Z
     %           + AXOMAT(4,5) * U + AXOMAT(5,5) * V + AXOMAT(6,5) * W
      AXYZUVW(6) = AXOMAT(1,6) * X + AXOMAT(2,6) * Y + AXOMAT(3,6) * Z
     %           + AXOMAT(4,6) * U + AXOMAT(5,6) * V + AXOMAT(6,6) * W
C
      END
