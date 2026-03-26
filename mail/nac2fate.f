      SUBROUTINE NAC2FATE( NF1, NF2, NAR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LE NUMERO LOCAL DE L'ARETE COMMUNE AUX 2 FACES
C -----    NF1 NF2 D'UN MEME TETRAEDRE
C
C ENTREES:
C --------
C NF1,NF2: NUMERO DES 2 FACES DU TETRAEDRE 123 234 341 412
C
C SORTIE :
C --------
C NAR    : NUMERO DE 1 A 6 DE L'ARETE RETROUVEE  12 23 31 41 42 43
C          0 SI ARETE NON RETROUVEE NF1 OU NF2 INCORRECTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   NOVEMBRE 1993
C MODIFS : ALAIN PERRONNET LJLL UPMC et St Pierre du Perray Fevrier 2016
C2345X7..............................................................012
C     SAUVEGARDE DES VALEURS ENTREES
      IF( NF1 .LT. NF2 ) THEN
         NFF1 = NF1
         NFF2 = NF2
      ELSE
         NFF1 = NF2
         NFF2 = NF1
      ENDIF

      GOTO( 10, 20, 30, 40 ), NFF1

C     FACE 1
 10   GOTO( 40, 12, 13, 14 ), NFF2

 12   NAR = 2
      GOTO 50

 13   NAR = 3
      GOTO 50

 14   NAR = 1
      GOTO 50

C     FACE 2
 20   GOTO( 40, 40, 23, 24 ), NFF2

 23   NAR = 6
      GOTO 50

 24   NAR = 5
      GOTO 50

C     FACE 3
 30   GOTO( 40, 40, 40, 34 ), NFF2

 34   NAR = 4
      GOTO 50

C     ARETE NON RETROUVEE NF1 OU NF2 INCORRECTS
 40   NAR = 0

 50   RETURN
      END
