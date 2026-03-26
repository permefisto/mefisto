      SUBROUTINE DFM1R2( DF , DELTA , DF1 )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DU DETERMINANT D UNE MATRICE 2 * 2 ET DE SON INVERSE
C ----- METHODE DIRECTE
C         ( 11  12 )           ( 22  -12 )
C                    **-1  =               /( 11 * 22 - 12 * 21 )
C         ( 21  22 )           (-21   11 )
C
C PARAMETRE D ENTREE :
C --------------------
C DF     : MATRICE 2 * 2
C
C PARAMETRES RESULTATS :
C ----------------------
C DELTA  : DETERMINANT DE LA MATRICE DF
C DF1    : INVERSE DE LA MATRICE DF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1990
C23456---------------------------------------------------------------012
      DOUBLE PRECISION DF(2,2),DELTA,DF1(2,2),S
C
C     LE DETERMINANT
C     --------------
      S     = DF(1,1)
      DELTA = S * DF(2,2) - DF(1,2) * DF(2,1)
C
C     LA MATRICE DF1 INVERSE DE DF
C     ----------------------------
      DF1(1,1) =  DF(2,2) / DELTA
      DF1(2,2) =  S       / DELTA
      DF1(1,2) = -DF(1,2) / DELTA
      DF1(2,1) = -DF(2,1) / DELTA
      END
