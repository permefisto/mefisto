      SUBROUTINE NSC2FATE( NF1, NF2, NS1, NS2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   DANS UN TETRAEDRE, RETROUVER LE NUMERO NS1 NS2 ( 1 A 4 )
C -----   DES 2 SOMMETS DE L'ARETE COMMUNE AUX 2 FACES
C         NF1 (1 A 3 ) < NF2 ( 2 A 4 )
C         POUR LES NUMEROTATIONS LOCALES AU TETRAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Janvier 2017
C2345X7..............................................................012
      GOTO( 110, 120, 130 ), NF1

C     FACE 1
 110  GOTO( 110, 112, 113, 114 ), NF2

C     FACE 2
 112  NS1 = 2
      NS2 = 3
      GOTO 900

C     FACE 3
 113  NS1 = 1
      NS2 = 3
      GOTO 900

C     FACE 4
 114  NS1 = 1
      NS2 = 2
      GOTO 900

C     FACE 2
 120  GOTO( 120, 120, 123, 124 ), NF2

C     FACE 3
 123  NS1 = 3
      NS2 = 4
      GOTO 900

C     FACE 4
 124  NS1 = 4
      NS2 = 2
      GOTO 900

C     FACE 3
 130  NS1 = 4
      NS2 = 1

 900  RETURN
      END
