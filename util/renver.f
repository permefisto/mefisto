      SUBROUTINE RENVER( NBMOTS, MTAB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RENVERSER SUR LUI MEME L'ORDRE DE RANGEMENT DU TABLEAU MTAB
C ----- MTAB(1) EST PERMUTE AVEC MTAB(NBMOTS)
C       MTAB(2) EST PERMUTE AVEC MTAB(NBMOTS-1) ...
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS   SEPTEMBRE 1996
C2345X7..............................................................012
      INTEGER  MTAB(NBMOTS)
C
      DO 10 I = 1,NBMOTS/2
         J = NBMOTS + 1 - I
         K       = MTAB(I)
         MTAB(I) = MTAB(J)
         MTAB(J) = K
 10   CONTINUE
      END
