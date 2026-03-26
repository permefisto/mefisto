      SUBROUTINE LONGITUD( X, Y, R, ALFA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     CALCULER LA LONGITUDE D'UN POINT LE LONG D'UN CERCLE
C ----
C ENTREES:
C --------
C X, Y   : 2 COORDONNEES DU POINT SUR LE CERCLE
C R      : RAYON DU CERCLE  UTILE POUR LA RELATIVITE A ZERO DE X
C
C SORTIE :
C --------
C ALFA   : L'ANGLE LONGITUDE EN RADIANS DE 0 A 2 PI DU POINT X Y
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2011
C2345X7..............................................................012
      DOUBLE PRECISION  PI, ALFAD
      INTRINSIC         ATAN, ABS
C
      PI = ATAN(1D0) * 4D0
C
      IF( ABS(X) .LT. 1E-4*R ) THEN
         IF( Y .GE. 0 ) THEN
            ALFAD = PI/2D0
         ELSE
            ALFAD =-PI/2D0
         ENDIF
      ELSE
         ALFAD = ATAN( Y / X )
      ENDIF
C
      IF( X .LT. 0 ) ALFAD = ALFAD + PI
C
      IF( ALFAD .LT. -1D-4 ) ALFAD = ALFAD + PI * 2D0
C
C     CONTRER LES ERREURS D'ARRONDI
      IF( ALFAD .LT. 0D0 ) ALFAD = 0
C
      ALFA = REAL( ALFAD )
      RETURN
      END
