      SUBROUTINE DFM1R3( DF, DFM1, IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DU DETERMINANT D UNE MATRICE 3 * 3 ET DE SON INVERSE
C ----- METHODE DE GAUSS A PIVOT TOTAL
C
C PARAMETRE D ENTREE PUIS MODIFIE :
C ---------------------------------
C DF     : MATRICE 3 * 3
C
C PARAMETRES RESULTATS :
C ----------------------
C DFM1   : INVERSE DE LA MATRICE DF (DIFFERENTE OU NON DE DF)
C IERR   : 0 SI PAS D'ERREUR, 1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1995
C23456---------------------------------------------------------------012
      DOUBLE PRECISION DF (3,3), DFM1(3,3)
      DOUBLE PRECISION DF1(3,6)
C
C     LA MATRICE IDENTITE SECONDS MEMBRES DU SYSTEME LINEAIRE A RESOUDRE
C     ------------------------------------------------------------------
      DO 20 I=1,3
         DO 10 J=1,3
            DF1(I,  J) = DF(I,J)
            DF1(I,3+J) = 0D0
   10    CONTINUE
         DF1(I,3+I) = 1D0
   20 CONTINUE
C
C     INVERSION DE LA MATRICE PAR GAUSS PIVOT TOTAL
C     ---------------------------------------------
      CALL GAUSPT( 3, 3, DF1, IERR)
      IF( IERR .NE. 0 ) RETURN
C
C     LES RESULTATS SONT REPORTES DANS DFM1
      DO 50 I=1,3
         DO 40 J=1,3
            DFM1(I,J) = DF1(I,3+J)
   40    CONTINUE
   50 CONTINUE
C
      RETURN
      END
