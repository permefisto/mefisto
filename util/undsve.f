      SUBROUTINE UNDSVE( N, NB, NVECT, NO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETOUVER LE NUMERO NO DANS NVECT DE N
C -----
C ENTREES :
C ---------
C N      : LE NOMBRE A RETROUVER
C NB     : NOMBRE DE COMPOSANTES DU VECTEUR
C NVECT  : LES VALEURS DU VECTEUR
C
C SORTIES :
C ---------
C NO     : LE NUMERO DE LA PREMIERE COMPOSANTE DE NVECT EGALE A N
C          0 SI N N'EST PAS UNE COMPOSANTE DE NVECT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1993
C23456...............................................................012
      INTEGER    NVECT(1:NB)
C
      DO 10 NO=1,NB
         IF( N .EQ. NVECT(NO) ) RETURN
 10   CONTINUE
C
C     N NON RETROUVE
      NO = 0
      END
