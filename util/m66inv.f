      SUBROUTINE M66INV( DF,  DFM1, DETDF, IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE INVERSE DE LA MATRICE DF 6*6
C -----    ET DU DETERMINANT DE DF PAR LA METHODE DE GAUSS A PIVOT TOTAL
C
C ENTREE :
C --------
C DF     : MATRICE 6 * 6
C
C SORTIES:
C --------
C DFM1   : INVERSE DE LA MATRICE DF  (DIFFERENTE OU NON DE DF)
C DETDF  : DETERMINANT DE LA MATRICE DF (=> 0 SI A NON INVERSIBLE)
C IERR   : 0 SI PAS D'ERREUR, 1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University           Juillet 2005
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  DF (6,6), DFM1(6,6), DETDF
      DOUBLE PRECISION  DF1(6,12)
C
C     MATRICE IDENTITE = SECONDS MEMBRES DU SYSTEME LINEAIRE A RESOUDRE
C     -----------------------------------------------------------------
      DO 20 I=1,6
         DO 10 J=1,6
            DF1(I,  J) = DF(I,J)
            DF1(I,6+J) = 0D0
   10    CONTINUE
         DF1(I,6+I) = 1D0
   20 CONTINUE
C
C     INVERSION DE LA MATRICE PAR GAUSS PIVOT TOTAL
C     ---------------------------------------------
      CALL GAUSDE( 6, 6, DF1, DETDF, IERR )
      IF( IERR .NE. 0 ) RETURN
      IF( ABS(DETDF) .LT. 1D-40 )  THEN
C        LA MATRICE EST CONSIDEREE COMME NON INVERSIBLE
         DETDF = 0D0
         RETURN
      ENDIF
C
C     LES RESULTATS SONT REPORTES DANS DFM1
      DO 50 I=1,6
         DO 40 J=1,6
            DFM1(I,J) = DF1(I,6+J)
   40    CONTINUE
   50 CONTINUE
      END
