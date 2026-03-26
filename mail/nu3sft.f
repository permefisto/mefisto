      SUBROUTINE NU3SFT( NUFACE, NS1, NS2, NS3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LES 3 NUMEROS DES SOMMETS DE LA FACE NUFACE
C -----    D'UN TETRAEDRE
C
C ENTREES:
C --------
C NUFACE : NUMERO ( 1 A 4 ) DE LA FACE DU TETRAEDRE
C
C SORTIES:
C --------
C NS1,NS2,NS3 : NUMERO DES 3 SOMMETS DE LA FACE NUFACE DU TETRAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1992
C2345X7..............................................................012
      GOTO( 1, 2, 3, 4 ), NUFACE
C
 1    NS1 = 1
      NS2 = 3
      NS3 = 2
      GOTO 10
C
 2    NS1 = 2
      NS2 = 3
      NS3 = 4
      GOTO 10
C
 3    NS1 = 3
      NS2 = 1
      NS3 = 4
      GOTO 10
C
 4    NS1 = 1
      NS2 = 2
      NS3 = 4
C
 10   RETURN
      END
