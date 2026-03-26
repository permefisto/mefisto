      SUBROUTINE EMPITR( NB, NOSOM1, MXSOMM, MXPILE, LHPILE, LAPILE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EMPILER DANS LA PILE LAPILE LE NO DES 3 SOMMETS DES
C -----    SOUS-TRIANGLES DU TRIANGLE RECTANGLE UNITE
C
C ENTREES :
C ---------
C NB     : NOMBRE D INTERVALLES DE SUBDIVISION DE CHAQUE ARETE
C NOSOM1 : NO - 1 DU 1-ER SOMMET A GENERER
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS   DANS LA PILE DES SOMMETS
C MXPILE : NOMBRE MAXIMAL DE TRIANGLES DANS LA PILE LAPILE
C
C PARAMETRES MODIFIES :
C ---------------------
C LHPILE : POINTEUR SUR LE SOMMET DE LA PILE LAPILE
C LAPILE : PILE ( 3 , MXPILE ) LA PILE DES 3 NO DES SOMMETS DE TRIANGLES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1990
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      INTEGER           LAPILE(3,MXPILE)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     GENERATION DU NO DES SOMMETS DES SOUS-TRIANGLES
C     ===============================================
      IF ( NOSOM1 + (NB + 1) * (NB + 2) / 2 .GT. MXSOMM ) GOTO 100
      K1 = NOSOM1
C
      DO 20 I=1,NB
C
C        LE TRIANGLE (I+1,1) , (I+1,2) , (I,J)
C        (I,J) => K = I * (I-1) / 2 + J
C        -------------------------------------
         K1 = K1 + 1
         K2 = K1 + I
         K3 = K2 + 1
         K4 = K1 + 1
         CALL EMPIL3 ( LHPILE , MXPILE , LAPILE , K2 , K3 , K1 )
C
         I1 = I - 1
         IF ( I1 .EQ. 0 ) GOTO 20
         DO 10 J=1,I1
C
C           LE TRIANGLE (I,J+1) , (I,J) , (I+1,J+1)
C           ---------------------------------------
            CALL EMPIL3 ( LHPILE , MXPILE , LAPILE , K4 , K1 , K3 )
C
C           LE TRIANGLE (I+1,J+1) , (I+1,J+2) , (I,J+1)
C           -------------------------------------------
            K1 = K4
            K2 = K3
            K4 = K4 + 1
            K3 = K3 + 1
            CALL EMPIL3 ( LHPILE , MXPILE , LAPILE , K2 , K3 , K1 )
   10    CONTINUE
   20 CONTINUE
      RETURN
C
C     PILE SATUREE
C     ============
 100  NBLGRC(NRERR) = 1
      KERR(1) = 'ERREUR EMPITR: PILE SATUREE'
      CALL LEREUR
      END
