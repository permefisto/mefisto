      SUBROUTINEXYZLONLAT( X, Y, Z, EPS,  LONGITUDE, LATITUDE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA LONGITUDE ET LA LATITUDE EN RADIANS DU POINT
C -----    DE COORDONNEES X Y Z
C
C ENTREES:
C --------
C X Y Z  : 3 COORDONNEES DU POINT
C EPS    : PRECISION EN X ou Y ou Z AU DESSOUS DE LAQUELLE LA COORDONNEE
C          EST SUPPOSEE NULLE POUR EVITER UNE DIVISION PAR ZERO
C
C SORTIES:
C --------
C LONGITUDE : en RADIANS DU POINT dans [0, 2PI]
C LATITUDE  : en RADIANS DU POINT dans [-PI/2, PI/2]
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Laboratoire J-L LIONS UMPC PARIS Decembre 2006
C.......................................................................
      DOUBLE PRECISION  X, Y, Z, EPS, LONGITUDE, LATITUDE
      DOUBLE PRECISION  ATAN, RXY, PIS2
C
C     PI/2
      PIS2 = ATAN( 1D0 ) * 2D0
C
C     LATITUDE DU SOMMET
      RXY = SQRT( X*X + Y*Y )
      IF( RXY .GT. EPS ) THEN
         LATITUDE = ATAN( Z / RXY )
      ELSE
         IF( Z .GT. EPS ) THEN
            LATITUDE = PIS2
         ELSE IF( Z .LT. -EPS ) THEN
            LATITUDE = -PIS2
         ELSE
            LATITUDE  = 0D0
            LONGITUDE = 0D0
            RETURN
         ENDIF
      ENDIF
C
C     LONGITUDE DU SOMMET
      IF( ABS(X) .GT. EPS ) THEN
         LONGITUDE = ATAN( Y / X )
         IF( X .LT. 0D0 ) THEN
            LONGITUDE = LONGITUDE + PIS2 * 2D0
         ENDIF
      ELSE
         IF( Y .GT. EPS ) THEN
            LONGITUDE = PIS2
         ELSE IF( Y .LT. -EPS ) THEN
            LONGITUDE = -PIS2
         ELSE
            LONGITUDE = 0D0
         ENDIF
      ENDIF
C
      RETURN
      END
