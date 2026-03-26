      SUBROUTINE CTEPRO( AMAX, NTDL, NCODSA, MU, A )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MULTIPLIER PAR AMAX TOUS LES COEFFICIENTS DE LA MATRICE A
C -----
C
C ENTREES:
C --------
C AMAX   : VALEUR A MULTIPLIER
C NCODSA : CODE DE STOCKAGE DE LA MATRICE A PROFIL
C          1 SYMETRIQUE
C          0 DIAGONALE
C         -1 NON SYMETRIQUE
C MU     : MU(0)=0
C          MU(I)=NUMERO DANS A DU I-EME COEFFICIENT DIAGONAL DE A
C
C MODIFIE:
C --------
C A      : MATRICE PROFIL EN MEMOIRE CENTRALE A COMPRIMER
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS UPMC     OCTOBRE 1998
C....67..............................................................012
      INTEGER           MU(0:NTDL)
      DOUBLE PRECISION  A(*), AMAX
C
      IF( NCODSA .EQ. 0 ) THEN
C
C        LA MATRICE EST DIAGONALE
C        ------------------------
         NB = NTDL
C
      ELSE
C
C        LA MATRICE N'EST PAS DIAGONALE
C        ------------------------------
         NB = MU(NTDL)
      ENDIF
C
      DO 10 I=1,NB
         A(I) = AMAX * A(I)
 10   CONTINUE
      END
