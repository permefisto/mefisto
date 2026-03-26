      SUBROUTINE NUFSOT( NUOT, NOSMMT, NBFACE, NUFACE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETROUVER LES NUMEROS DES FACES COMMUNES AU SOMMET NOSMMT
C -----  D'UN OT DE NUMERO NUOT
C
C ENTREES:
C --------
C NUOT   : NUMERO DE L'OT, >0 SI OCTAEDRE,  <0 SI TETRAEDRE
C NOSMMT : 1 A 4 POUR UN TETRAEDRE ( NUOT<0 )
C          1 A 6 POUR UN OCTAEDRE  ( NUOT>0 )
C
C SORTIES:
C --------
C NBFACE : NOMBRE DE FACES AYANT EN COMMUN LE SOMMET NOSMMT
C NUFACE : NUMEROS DES NBFACE FACES DE SOMMET NOSMMT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1992
C2345X7..............................................................012
      INTEGER  NUFACE(4)
C
      IF( NUOT .GT. 0 ) THEN
C        OCTAEDRE
         NBFACE = 4
         IF( NOSMMT .EQ. 1 ) THEN
             NUFACE(1) = 1
             NUFACE(2) = 2
             NUFACE(3) = 3
             NUFACE(4) = 4
         ELSE IF( NOSMMT .EQ. 6 ) THEN
             NUFACE(1) = 5
             NUFACE(2) = 6
             NUFACE(3) = 7
             NUFACE(4) = 8
         ELSE IF( NOSMMT .EQ. 2 ) THEN
             NUFACE(1) = 4
             NUFACE(2) = 1
             NUFACE(3) = 8
             NUFACE(4) = 5
         ELSE
             NUFACE(1) = NOSMMT - 2
             NUFACE(2) = NOSMMT - 1
             NUFACE(3) = NOSMMT + 2
             NUFACE(4) = NOSMMT + 3
         ENDIF
      ELSE
C        TETRAEDRE
         NBFACE = 3
         IF( NOSMMT .EQ. 1 ) THEN
             NUFACE(1) = 1
             NUFACE(2) = 3
             NUFACE(3) = 4
         ELSE IF( NOSMMT .EQ. 2 ) THEN
             NUFACE(1) = 1
             NUFACE(2) = 2
             NUFACE(3) = 4
         ELSE IF( NOSMMT .EQ. 3 ) THEN
             NUFACE(1) = 1
             NUFACE(2) = 2
             NUFACE(3) = 3
         ELSE
             NUFACE(1) = 2
             NUFACE(2) = 3
             NUFACE(3) = 4
         ENDIF
      ENDIF
      END
