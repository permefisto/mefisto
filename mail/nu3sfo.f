      SUBROUTINE NU3SFO( NUFACE, NS1, NS2, NS3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LES 3 NUMEROS DES SOMMETS DE LA FACE NUFACE
C -----    D'UN OCTAEDRE
C
C ENTREES:
C --------
C NUFACE : NUMERO ( 1 A 8 ) DE LA FACE DE L'OCTAEDRE
C
C SORTIES:
C --------
C NS1,NS2,NS3 : NUMERO DES 3 SOMMETS DE LA FACE NUFACE DE L'OCTAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1992
C2345X7..............................................................012
      GOTO( 1, 2, 3, 4, 5, 6, 7, 8 ), NUFACE
C
 1    NS1 = 1
      NS2 = 2
      NS3 = 3
      GOTO 10
C
 2    NS1 = 1
      NS2 = 3
      NS3 = 4
      GOTO 10
C
 3    NS1 = 1
      NS2 = 4
      NS3 = 5
      GOTO 10
C
 4    NS1 = 1
      NS2 = 5
      NS3 = 2
      GOTO 10
C
 5    NS1 = 6
      NS2 = 3
      NS3 = 2
      GOTO 10
C
 6    NS1 = 6
      NS2 = 4
      NS3 = 3
      GOTO 10
C
 7    NS1 = 6
      NS2 = 5
      NS3 = 4
      GOTO 10
C
 8    NS1 = 6
      NS2 = 2
      NS3 = 5
C
 10   RETURN
      END
