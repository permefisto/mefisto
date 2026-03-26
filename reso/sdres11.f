       SUBROUTINE SDRES11(LU,X,R,N)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RESOLUTION DU SYSTEME LINEAIRE LU * X = R  (METHODE DES JOINTS)
C -----
C        LU(1,I) = SOUS-DIAGONALE PRINCIPALE
C        LU(2,I) = DIAGONALE
C        LU(3,I) = SUR-DIAGONALE PRINCIPALE
C
C ENTREES :
C ---------
C LU : LES FACTEURS DE GAUSS
C R  : LE SECOND MEMBRE
C N  : LE NOMBRE D'INCONNUES
C
C SORTIE :
C --------
C X   : LA SOLUTION
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1993
C23456---------------------------------------------------------------012
C
      DOUBLE PRECISION  LU(3,N),X(N),R(N)
C
C     LA DESCENTE
      X(1) = R(1)
      DO 1 I  = 2 , N
         X(I) = R(I) - LU(1,I) * X(I-1)
 1    CONTINUE
C     LA REMONTEE
      X(N) = X(N) / LU(2,N)
      DO 2 I  = N-1 , 1 , -1
         X(I) = ( X(I) - LU(3,I) * X(I+1) ) / LU(2,I)
 2    CONTINUE
C
      RETURN
      END
