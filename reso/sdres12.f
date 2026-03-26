       SUBROUTINE SDRES12(Q,SM,P,N)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : MULTIPLICATION PAR LA MATRICE Q  (METHODE DES JOINTS)
C -----
C        Q(1,I) = SOUS-DIAGONALE PRINCIPALE
C        Q(2,I) = DIAGONALE
C        Q(3,I) = SUR-DIAGONALE PRINCIPALE
C
C ENTREE :
C ---------
C Q  : LA MATRICE
C SM : LE VECTEUR
C N  : LE NOMBRE DE VALEURS
C
C SORTIE :
C --------
C P : P = Q * SM
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1993
C23456---------------------------------------------------------------012
C
      DOUBLE PRECISION  Q(3,N),SM(N),P(N)
C
C     LA PREMIERE LIGNE
C     P(1) = Q(2,1) * SM(1) + Q(3,1) * SM(2)
C     LES AUTRES
      DO 1 I  = 2 , N - 1
         P(I) = Q(1,I) * SM(I-1) + Q(2,I) * SM(I) + Q(3,I) * SM(I+1)
 1    CONTINUE
C     LA DERNIERE LIGNE
C     P(N) = Q(1,N) * SM(N-1) + Q(2,N) * SM(N)
C
      RETURN
      END
