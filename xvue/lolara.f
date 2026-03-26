      SUBROUTINE LOLARA( XC, YC, ZC, XP, YP, ZP,
     %                   RADLON, RADLAT, LEPOLE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA LONGITUDE ET LATITUDE D'UN POINT SUR LA SPHERE
C -----    DE CENTRE (XC,YC,ZC) ET DE RAYON RAYON
C
C ENTREES:
C --------
C RAYON  : RAYON DE LA SPHERE
C XC     : COORDONNEE X DU CENTRE DE LA SPHERE
C YC     : COORDONNEE Y DU CENTRE DE LA SPHERE
C ZC     : COORDONNEE Z DU CENTRE DE LA SPHERE
C XP     : COORDONNEE X DU POINT SUR LA SPHERE
C YP     : COORDONNEE Y DU POINT SUR LA SPHERE
C ZP     : COORDONNEE Z DU POINT SUR LA SPHERE
C
C SORTIES:
C --------
C RADLON : LONGITUDE EN RADIANS DU POINT (XP,YP,ZP) SUR LA SPHERE
C RADLAT : LATITUDE  EN RADIANS DU POINT (XP,YP,ZP) SUR LA SPHERE
C LEPOLE : 1 POLE NORD, -1 POLE SUD, 0 NON POLE DE LA SPHERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A. PERRONNET  ANALYSE NUMERIQUE UPMC  PARIS         JUIN 1996
C234567--------------------------------------------------------------012
      PI = ATAN(1.0) * 4.0
C
C     RAYON  : RAYON DE LA SPHERE
      RAYON = SQRT( (XP-XC)**2+(YP-YC)**2+(ZP-ZC)**2 )
C
      X = (XP - XC) / RAYON
      Y = (YP - YC) / RAYON
      Z = (ZP - ZC) / RAYON
C
      R = SQRT( X * X + Y * Y )
C
C     TRAITEMENT DE LA LATITUDE PROCHE DE PI/2
      IF( R .LT. 1E-4 ) THEN
C        ANGLE JUGE DE PI/2
         RADLAT = PI / 2
         IF( Z .LT. 0 ) THEN
C           POLE SUD
            RADLAT = -RADLAT
            LEPOLE = -1
         ELSE
C           POLE NORD
            LEPOLE = 1
         ENDIF
C        LA LONGITUDE EST FIXEE ARBITRAIREMENT A ZERO
         RADLON = 0
         RETURN
      ELSE
C        ANGLE DIFFERENT
         RADLAT = ATAN( Z / R )
         LEPOLE = 0
      ENDIF
C
C     TRAITEMENT DE LA LONGITUDE PROCHE DE PI/2
      XX = ABS(X)
      IF( XX .LT. 1E-4 ) THEN
C         ANGLE JUGE DE PI/2
          RADLON = PI / 2
          IF( Y .LT. 0 ) RADLON = -RADLON
      ELSE
C        ANGLE DIFFERENT
         RADLON = ATAN( Y / X )
         IF( X .LT. 0 ) THEN
C           ON AJOUTE PI SI X<0 => RESULTAT DANS [-PI/2, 3PI/2]
            RADLON = RADLON + PI
            IF( RADLON .GT. PI ) RADLON = RADLON - 2 * PI
         ENDIF
      ENDIF
C
      RETURN
      END
