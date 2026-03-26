      SUBROUTINE NU2SATE( NUARET, NS1, NS2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LES 2 NUMEROS DES SOMMETS DE L'ARETE NUARET
C -----    D'UN TETRAEDRE
C
C ENTREES:
C --------
C NUARET : NUMERO ( 1 A 6 ) DE L'ARETE DU TETRAEDRE
C
C SORTIES:
C --------
C NS1, NS2 : NUMERO DES 2 SOMMETS DE L'ARETE NUARET DU TETRAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1992
C2345X7..............................................................012
      IF( NUARET .LE. 3 ) THEN
         NS1 = NUARET
         NS2 = NUARET + 1
         IF( NS2 .EQ. 4 ) NS2 = 1
      ELSE
         NS1 = NUARET - 3
         NS2 = 4
      ENDIF

      RETURN
      END
