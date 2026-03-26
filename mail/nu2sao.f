      SUBROUTINE NU2SAO( NUARET, NS1, NS2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LES 2 NUMEROS DES SOMMETS DE L'ARETE NUARET
C -----    D'UN OCTAEDRE
C
C ENTREES:
C --------
C NUARET : NUMERO ( 1 A 12 ) DE L'ARETE DE L'OCTAEDRE
C
C SORTIES:
C --------
C NS1, NS2 : NUMERO DES 2 SOMMETS DE L'ARETE NUARET DE L'OCTAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1992
C2345X7..............................................................012
      IF( NUARET .LE. 4 ) THEN
         NS1 = 1
         NS2 = NUARET + 1
      ELSE IF( NUARET .LE. 8 ) THEN
         NS1 = NUARET - 3
         NS2 = NUARET - 2
         IF( NS2 .EQ. 6 ) NS2 = 2
      ELSE
         NS1 = NUARET - 7
         NS2 = 6
      ENDIF
      END
