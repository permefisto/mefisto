       SUBROUTINE SDRES13(QI,VE,NU,VS,N)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INTERPOLATION AUX POINTS DU JOINT  (METHODE DES JOINTS)
C -----
C
C ENTREE :
C ---------
C QI : LA MATRICE D'INTERPOLATION
C NU : LES NUMEROS DES VALEURS A INTERPOLER
C VE : LES VALEURS EN ENTREE
C N  : LE NOMBRE DE POINTS D'INTERPOLATION
C
C SORTIE :
C --------
C VS : LES VALEURS EN SORTIE
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1993
C23456---------------------------------------------------------------012
C
      DOUBLE PRECISION  QI(2,N),VE(N),VS(N)
      DIMENSION NU(2,N)
C
      DO I = 1 , N
         VS(I) = QI(1,I) * VE(NU(1,I)) + QI(2,I) * VE(NU(2,I))
      ENDDO
C
      RETURN
      END
