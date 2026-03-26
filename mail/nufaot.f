      SUBROUTINE NUFAOT( NUOT, NOARET, NUFAC1, NUFAC2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETROUVER LES 2 NUMEROS DES FACES COMMUNES A L'ARETE NOARET
C -----  D'UN OT DE NUMERO NUOT
C
C ENTREES:
C --------
C NUOT   : NUMERO DE L'OT, >0 SI OCTAEDRE,  <0 SI TETRAEDRE
C NOARET : 1 A 6  POUR UN TETRAEDRE ( NUOT<0 )
C          1 A 12 POUR UN OCTAEDRE  ( NUOT>0 )
C
C SORTIES:
C --------
C NUFAC1, NUFAC2 : LES 2 NUMEROS DES FACES COMMUNES A L'ARETE NOARET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1992
C2345X7..............................................................012
      IF( NUOT .GE. 0 ) THEN
C        OCTAEDRE
         IF( NOARET .LE. 4 ) THEN
            IF( NOARET .EQ. 1 ) THEN
               NUFAC1 = 4
               NUFAC2 = 1
            ELSE
               NUFAC1 = NOARET - 1
               NUFAC2 = NOARET
            ENDIF
         ELSE IF( NOARET .LE. 8 ) THEN
            NUFAC1 = NOARET - 4
            NUFAC2 = NOARET
         ELSE
            IF( NOARET .EQ. 9 ) THEN
               NUFAC1 = 5
               NUFAC2 = 8
            ELSE
               NUFAC1 = NOARET - 4
               NUFAC2 = NOARET - 5
            ENDIF
         ENDIF
      ELSE
C        TETRAEDRE
         IF( NOARET .LE. 3 ) THEN
            NUFAC1 = 1
            IF( NOARET .EQ. 1 ) THEN
               NUFAC2 = 4
            ELSE
               NUFAC2 = NOARET
            ENDIF
         ELSE
            IF( NOARET .EQ. 4 ) THEN
               NUFAC1 = 3
               NUFAC2 = 4
            ELSE IF( NOARET .EQ. 5 ) THEN
               NUFAC1 = 4
               NUFAC2 = 2
            ELSE
               NUFAC1 = 2
               NUFAC2 = 3
            ENDIF
         ENDIF
      ENDIF
      END
