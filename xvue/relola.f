      SUBROUTINE RELOLA( PTV, OEIL, DEGLON, DEGLAT )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA LONGITUDE ET LATITUDE DEFINIE PAR LE VECTEUR
C -----    PTV -> OEIL
C
C ENTREES :
C ---------
C PTV     : LE POINT VISE
C OEIL    : LA POSITION DE L'OEIL
C
C SORTIES :
C ---------
C DEGLON  : DEGRES DE LONGITUDE
C DEGLAT  : DEGRES DE LATITUDE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS        JUIN 1994
C ......................................................................
      REAL  PTV(3), OEIL(3)
C
C     PTV ET OEIL SONT ILS CONFONDUS ?
      CALL XYZIDE( PTV, OEIL, I )
      IF( I .NE. 0 ) THEN
C        OUI : LONGITUDE = LATITUDE = 0
         DEGLON = 0
         DEGLAT = 0
         RETURN
      ENDIF
C
C     LES COMPOSANTES DU VECTEUR PTV->OEIL
      X = OEIL(1) - PTV(1)
      Y = OEIL(2) - PTV(2)
      Z = OEIL(3) - PTV(3)
      D = SQRT( X*X + Y*Y + Z*Z )
      R = 45.0 / ATAN(1.0)
C
C     LA LONGITUDE
      IF( ABS(X) .LE. 1E-4 * D ) THEN
         IF( Y .LT. 0 ) THEN
            DEGLON = -90
         ELSE
            DEGLON =  90
         ENDIF
      ELSE
         DEGLON = ATAN( Y / X ) * R
C        LONGITUDE ENTRE 0 ET 360 DEGRES
         IF( X .LT. 0 .AND. Y .LT. 0 ) DEGLON = DEGLON + 180
         IF( X .LT. 0 .AND. Y .GT. 0 ) DEGLON = DEGLON + 180
         IF( DEGLON .GT. 180 ) DEGLON = DEGLON - 360
      ENDIF
C
C     LA LATITUDE
      X = SQRT( X*X + Y*Y )
      IF( ABS(X) .LE. 1E-4 * D ) THEN
         IF( Z .LT. 0 ) THEN
            DEGLAT = -90
         ELSE
            DEGLAT =  90
         ENDIF
      ELSE
         DEGLAT = ATAN( Z / X ) * R
      ENDIF
      END
