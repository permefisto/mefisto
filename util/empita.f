      SUBROUTINE EMPITA( NB, NULGAR, MXPILE, LHPILE, LAPILE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EMPILER DANS LA PILE LAPILE
C -----    LE NO DE LIGNE DES ARETES DES SOUS-TRIANGLES DU TRIANGLE UNITE
C          ORDRE : DU HAUT EN BAS ET DE LA GAUCHE VERS LA DROITE
C ENTREES :
C ---------
C NB     : NOMBRE D INTERVALLES DE SUBDIVISION DE CHAQUE ARETE
C NULGAR : NO DE LIGNE DES ARETES DU TRIANGLE UNITE
C          0 SI ARETE INTERNE
C MXPILE : NOMBRE MAXIMAL DE TRIANGLES DANS LA PILE LAPILE
C
C MODIFIES :
C ----------
C LHPILE : POINTEUR SUR LE SOMMET DE LA PILE LAPILE
C LAPILE : PILE ( 3, MXPILE ) LA PILE DES 3 NO DE LIGNES DES ARETES
C          DES TRIANGLES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1995
C23456---------------------------------------------------------------012
      INTEGER           NULGAR(3), LAPILE(3,MXPILE)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     GENERATION DU NO DES SOMMETS DES SOUS-TRIANGLES
C     ===============================================
C     TRAITEMENT DE L'ARETE 2 DU PREMIER TRIANGLE
      NLC2 = NULGAR(2)
      I1 = 0
      DO 20 I=1,NB
C
C        LE TRIANGLE (I+1,1) , (I+1,2) , (I,J) 1-ER TRIANGLE DE LA LIGNE
C        (I,J) => K = I * (I-1) / 2 + J
C        -------------------------------------
         CALL EMPIL3 ( LHPILE, MXPILE, LAPILE, 0, NLC2, NULGAR(3) )
         NLC2 = 0
C
         I1 = I - 1
         IF ( I1 .EQ. 0 ) GOTO 20
         DO 10 J=1,I1
C
C           LE TRIANGLE (I,J+1) , (I,J) , (I+1,J+1) GAUCHE
C           ---------------------------------------
            CALL EMPIL3 ( LHPILE , MXPILE , LAPILE , 0, 0, 0 )
C
C           LE TRIANGLE (I+1,J+1) , (I+1,J+2) , (I,J+1) DROITE
C           -------------------------------------------
            CALL EMPIL3 ( LHPILE , MXPILE , LAPILE , 0, 0, 0 )
   10    CONTINUE
C
C        RECTIFICATION DU DERNIER TRIANGLE (ARETE 2) DE CHAQUE LIGNE
         LAPILE(2,LHPILE) = NULGAR(2)
   20 CONTINUE
C
C     RECTIFICATION DE L'ARETE 1 DES TRIANGLES DE LA DERNIERE LIGNE
      NLC2 = LHPILE - 2 * I1
      DO 30 I=NLC2,LHPILE,2
         LAPILE(1,I) = NULGAR(1)
 30   CONTINUE
C
      RETURN
      END
