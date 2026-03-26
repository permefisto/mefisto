      SUBROUTINE CORTLT( MXTRIA, NOTRIA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : VERIFIER LA COHERENCE DE TOUS LES TRIANGLES RANGES DANS NOTRIA
C -----
C ENTREES :
C ---------
C MXTRIA : NOMBRE DE TRIANGLES DELARABLES DANS NOTRIA
C NOTRIA : NUMERO DES 3 SOMMETS ET TRIANGLES OPPOSES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS         MAI 1995
C2345X7..............................................................012
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER  NOTRIA(1:6,1:MXTRIA)
C
      DO 10 NT=1,MXTRIA
         CALL CORTRI( NT, NOTRIA )
 10   CONTINUE
      END
