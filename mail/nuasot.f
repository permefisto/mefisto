      SUBROUTINE NUASOT( NUOT, NOSMMT, NBARET, NUARET )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETROUVER LES NUMEROS DES ARETES COMMUNES AU SOMMET NOSMMT
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
C NBARET : NOMBRE DE ARETES AYANT EN COMMUN LE SOMMET NOSMMT
C NUARET : NUMEROS DES NBARET ARETES DE SOMMET NOSMMT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1992
C2345X7..............................................................012
      INTEGER NUARET(4)
C
      IF( NUOT .GT. 0 ) THEN
C        OCTAEDRE
         NBARET = 4
         IF( NOSMMT .EQ. 1 ) THEN
             NUARET(1) = 1
             NUARET(2) = 2
             NUARET(3) = 3
             NUARET(4) = 4
         ELSE IF( NOSMMT .EQ. 6 ) THEN
             NUARET(1) = 9
             NUARET(2) = 10
             NUARET(3) = 11
             NUARET(4) = 12
         ELSE IF( NOSMMT .EQ. 2 ) THEN
             NUARET(1) = 1
             NUARET(2) = 8
             NUARET(3) = 5
             NUARET(4) = 9
         ELSE
             NUARET(1) = NOSMMT - 1
             NUARET(2) = NOSMMT + 2
             NUARET(3) = NOSMMT + 3
             NUARET(4) = NOSMMT + 7
         ENDIF
      ELSE
C        TETRAEDRE
         NBARET = 3
         IF( NOSMMT .EQ. 1 ) THEN
             NUARET(1) = 1
             NUARET(2) = 3
             NUARET(3) = 4
         ELSE IF( NOSMMT .EQ. 2 ) THEN
             NUARET(1) = 1
             NUARET(2) = 2
             NUARET(3) = 5
         ELSE IF( NOSMMT .EQ. 3 ) THEN
             NUARET(1) = 2
             NUARET(2) = 3
             NUARET(3) = 6
         ELSE
             NUARET(1) = 4
             NUARET(2) = 5
             NUARET(3) = 6
         ENDIF
      ENDIF
      END
