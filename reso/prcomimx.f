      SUBROUTINE PRCOMIMX( NBCOEFA, A, AMIN, AMAX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA VALEUR MINIMALE et MAXIMALE des COEFFICIENTS
C -----    de la MATRICE A
C
C ENTREES:
C --------
C NBCOEFA: NOMBRE DE COEFFICIENTS DE LA MATRICE
C A      : MATRICE EN MEMOIRE CENTRALE
C
C SORTIE :
C --------
C AMIN   : VALEUR du COEFFICIENT MINIMAL de A
C AMAX   : VALEUR du COEFFICIENT MAXIMAL de A
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Octobre 2010
C....67..............................................................012
      DOUBLE PRECISION  A(NBCOEFA), AMIN, AMAX
C
      AMIN = A(1)
      AMAX = A(1)
C
      DO I=2,NBCOEFA
         IF( A(I) .LT. AMIN ) THEN
            AMIN = A(I)
         ELSE IF( A(I) .GT. AMAX ) THEN
            AMAX = A(I)
         ENDIF
      ENDDO
C
      RETURN
      END
