      SUBROUTINE PRCOMX( NTDL, NCODSA, MU, A, AMAX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA VALEUR ABSOLUE DU COEFFICIENT DIAGONAL MAXIMAL
C -----
C
C ENTREES:
C --------
C NTDL   : NOMBRE TOTAL DE DL BLOQUES+LIBRES
C NCODSA : CODE DE STOCKAGE DE LA MATRICE A PROFIL
C          1 SYMETRIQUE
C          0 DIAGONALE
C         -1 NON SYMETRIQUE
C MU     : MU(0)=0
C          MU(I)=NUMERO DANS A DU I-EME COEFFICIENT DIAGONAL DE A
C A      : MATRICE PROFIL EN MEMOIRE CENTRALE A COMPRIMER
C
C SORTIE :
C --------
C AMAX   : VALEUR ABSOLUE DU COEFFICIENT MAXIMAL
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS UPMC     OCTOBRE 1998
C....67..............................................................012
      INTEGER           MU(0:NTDL)
      DOUBLE PRECISION  A(*), AMAX
C
      AMAX = 0D0
      IF( NCODSA .EQ. 0 ) THEN
C
C        LA MATRICE EST DIAGONALE
C        ------------------------
         DO 10 I=1,NTDL
            AMAX = MAX( AMAX, ABS( A(I) ) )
 10      CONTINUE
C
      ELSE
C
C        LA MATRICE N'EST PAS DIAGONALE
C        ------------------------------
         DO 20 I=1,NTDL
            AMAX = MAX( AMAX, ABS( A(MU(I)) ) )
 20      CONTINUE
C
      ENDIF
      END
