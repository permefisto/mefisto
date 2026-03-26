      SUBROUTINE COMPRT( N0, A0, N1, A1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   SUPPRIMER LES VALEURS MULTIPLES DE A0
C -----   ET LES STOCKER DANS A1
C
C ENTREES:
C --------
C N0    : NOMBRE DE VALEURS DANS A0
C A0    : LES N0 VALEURS
C
C SORTIES :
C ---------
C N1    : NOMBRE DE VALEURS DIFFERENTES DANS A0
C A1    : LES N1 VALEURS DIFFERENTES
C
C A0 ET A1 PEUVENT ETRE CONFONDUS A L'APPEL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    FEVRIER 1992
C2345X7..............................................................012
      INTEGER   A0(N0),A1(N0)
C
      N1 = 0
      DO 20 I=1,N0
         K = A0(I)
         DO 10 J=1,N1
            IF( A0(J) .EQ. K ) GOTO 20
 10      CONTINUE
C        K NON RETROUVE DANS LES VALEURS PRECEDENTES
         N1 = N1 + 1
         A1( N1 ) = K
 20   CONTINUE
      END
