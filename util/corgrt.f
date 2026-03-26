      SUBROUTINE CORGRT( NBT, NOTR, NOTRIA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : VERIFIER LA COHERENCE DES TRIANGLES NOTR(1..NBT) DE NOTRIA
C -----
C ENTREES :
C ---------
C NBT    : NOMBRE DE TRIANGLES
C NOTR   : NUMERO DANS NOTRIA DES NBT TRIANGLES
C NOTRIA : NUMERO DES 3 SOMMETS ET TRIANGLES OPPOSES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS         MAI 1995
C2345X7..............................................................012
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER   NOTR(1:NBT), NOTRIA(1:6,*)
C
      DO 10 N=1,NBT
         NT = NOTR(N)
         CALL CORTRI( NT, NOTRIA )
 10   CONTINUE
      END
