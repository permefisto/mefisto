      SUBROUTINE XYZCUSP( LCUouSP, ARNOYAU, NANOYAU,
     %                    CUB3DI,  CUB3RG, NC3CUB, NBSTCH,
     %                    NUSTDF,  XYZDF, XYZCH )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES 3 COORDONNEES DES SOMMETS D'UN MAILLAGE DE 3-CUBES
C -----    FORMES PAR 3 DIFFERENCES FINIES ET NC3CUB COUCHES HOMOTHETIQUES
C          EN PROGRESSION GEOMETRIQUE
C
C ENTREES:
C --------
C LCUouSP: 0 SI CUBE DEMANDE, 1 SI SPHERE DEMANDEE
C ARNOYAU: LONGUEUR D'UNE ARETE DANS UNE DIRECTION DU NOYAU
C NANOYAU: NOMBRE D'ARETES DANS UNE DIRECTION DU NOYAU (DIFFERENCES FINIES)
C CUB3DI : LARGEUR TOTALE DU 3-CUBE ou SPHERE
C CUB3RG : PROGRESSION GEOMETRIQUE DES ARETES DES COUCHES
C NC3CUB : NOMBRE DE COUCHES HOMOTHETIQUES
C NBSTCH : NOMBRE DE SOMMETS D'UNE COUCHE
C
C SORTIES:
C --------
C NUSTDF : NO GLOBAL DES SOMMETS FRONTALIERS POUR LES DF DU NOYAU
C XYZDF  : 3 COORDONNEES DES SOMMETS DIFFERENCES FINIES  DU NOYAU
C XYZCH  : 3 COORDONNEES DES SOMMETS DES COUCHES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Laboratoire J-L LIONS UMPC PARIS Decembre 2006
C.......................................................................
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
C     CONSTRUCTION DU TABLEAU XYZ DU MAILLAGE EN HEXAEDRES
      REAL     XYZDF(3,0:NANOYAU,0:NANOYAU,0:NANOYAU)
      INTEGER  NUSTDF(0:NANOYAU,0:NANOYAU,0:NANOYAU)
      REAL     XYZCH(3,1:NBSTCH,1:NC3CUB)
      DOUBLE PRECISION  DBLE, CUB3DI2, PIS2, EPS, ATAN
      DOUBLE PRECISION  RXY, RAYON, LONGITUDE, LATITUDE, X, Y, Z, ARETEN
      LOGICAL  NBAIMPAIR
      INTEGER  NBAMILIEU
C
C     LE RAYON DU NOYAU EST EGAL A SA LARGEUR/2
      RAYON = ( ARNOYAU * NANOYAU ) * 0.5D0
C
C     CALCUL DES 3 COORDONNEES DES SOMMETS DIFFERENCES FINIES
C     LE CENTRE DU NOYAU EST L'ORIGINE DES AXES XYZ
      X0 = REAL( - RAYON )
C
C     LES SOMMETS DIFFERENCES FINIES DU NOYAU ET DE LA FRONTIERE
      NBSTDF = (NANOYAU+1)**3
      IF( LCUouSP .EQ. 0 ) THEN
C
C        CUBE DEMANDE
C        ------------
         DO 30 K=0,NANOYAU
            DO 20 J=0,NANOYAU
               DO 10 I=0,NANOYAU
C                 LES 3 COORDONNEES DU SOMMET DIFFERENCES FINIES DU NOYAU
                  XYZDF(1,I,J,K) = X0 + I * ARNOYAU
                  XYZDF(2,I,J,K) = X0 + J * ARNOYAU
                  XYZDF(3,I,J,K) = X0 + K * ARNOYAU
C
C                 LE SOMMET EST IL SUR LA FRONTIERE
C                 DES DIFFERENCES FINIES DU NOYAU?
                  IF( I.EQ.0 .OR. I.EQ.NANOYAU .OR.
     %                J.EQ.0 .OR. J.EQ.NANOYAU .OR.
     %                K.EQ.0 .OR. K.EQ.NANOYAU ) THEN
C                     SOMMET SUR LA FRONTIERE
                      NBSTDF = NBSTDF + 1
C                     SON NUMERO GLOBAL
                      NUSTDF(I,J,K) = NBSTDF
                  ENDIF
 10            CONTINUE
 20         CONTINUE
 30      CONTINUE
C
      ELSE
C
C        SPHERE DEMANDEE
C        ---------------
         PIS2   = ATAN( 1D0 ) * 2D0
         ARETEN = DBLE( ARNOYAU )
         EPS    = ARETEN * 1D-5
         NBAIMPAIR = MOD(NANOYAU,2) .EQ. 1
         NBAMILIEU = (NANOYAU+1)/2
C
         DO 80 K=0,NANOYAU
            NZ = ABS( NBAMILIEU - K )
            IF( NBAIMPAIR .AND. (K .GE. NBAMILIEU) ) NZ=NZ+1
C
            DO 70 J=0,NANOYAU
               NY = ABS( NBAMILIEU - J )
               IF( NBAIMPAIR .AND. (J .GE. NBAMILIEU) ) NY=NY+1
C
               DO 60 I=0,NANOYAU
                  NX = ABS( NBAMILIEU - I )
                  IF( NBAIMPAIR .AND. (I .GE. NBAMILIEU) ) NX=NX+1
C
C                 LE MAXIMUM DES 3 NUMEROS EST LE NUMERO DU RAYON
                  NR = MAX( NZ, MAX(NY,NX) )
C                 LE RAYON DU POINT DANS LA SPHERE
                  RAYON = NR * ARNOYAU
C
C                 LES 3 COORDONNEES DU SOMMET DIFFERENCES FINIES DU NOYAU
                  X = X0 + I * ARETEN
                  Y = X0 + J * ARETEN
                  Z = X0 + K * ARETEN
C
C                 LATITUDE DU SOMMET
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
                        GOTO 50
                     ENDIF
                  ENDIF
C
C                 LONGITUDE DU SOMMET
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
C                 LES XYZ DU SOMMET SUR LA SPHERE DE RAYON NR
 50               RXY = RAYON * COS( LATITUDE )
                  XYZDF(1,I,J,K) = REAL( RXY   * COS( LONGITUDE ) )
                  XYZDF(2,I,J,K) = REAL( RXY   * SIN( LONGITUDE ) )
                  XYZDF(3,I,J,K) = REAL( RAYON * SIN( LATITUDE  ) )
C
C                 LE SOMMET EST IL SUR LA FRONTIERE
C                 DES DIFFERENCES FINIES DU NOYAU?
                  IF( I.EQ.0 .OR. I.EQ.NANOYAU .OR.
     %                J.EQ.0 .OR. J.EQ.NANOYAU .OR.
     %                K.EQ.0 .OR. K.EQ.NANOYAU ) THEN
C                     SOMMET SUR LA FRONTIERE
                      NBSTDF = NBSTDF + 1
C                     SON NUMERO GLOBAL
                      NUSTDF(I,J,K) = NBSTDF
                  ENDIF
C
 60            CONTINUE
C
 70         CONTINUE
C
 80      CONTINUE
      ENDIF
C
C     TRAITEMENT DES COUCHES HOMOTHETIQUES A PARTIR DE LA SURFACE DU NOYAU
C     ====================================================================
C     NOMBRE DE SOMMETS FRONTALIERS OU D'UNE COUCHE
      WRITE(IMPRIM,*) 'NB SOMMETS FRONTALIERS DU NOYAU=', NBSTCH
C
C     CALCUL DES 3 COORDONNEES DES SOMMETS DES NC3CUB COUCHES HOMOTHETIQUES
C     => NBSTCH SOMMETS NOUVEAUX PAR COUCHE
      IF( NC3CUB .LE. 0 ) GOTO 9000
      CUB3DI2 = CUB3DI / 2D0
      DEPS  = 0.001*CUB3DI
      NUSDF = 0
      DO 130 K=0,NANOYAU
         DO 120 J=0,NANOYAU
            DO 110 I=0,NANOYAU
C
C              LE SOMMET EST IL SUR LA FRONTIERE DES DIFFERENCES FINIES?
               IF( I.EQ.0 .OR. I.EQ.NANOYAU .OR.
     %             J.EQ.0 .OR. J.EQ.NANOYAU .OR.
     %             K.EQ.0 .OR. K.EQ.NANOYAU ) THEN
C                  SOMMET SUR LA FRONTIERE
                  NUSDF = NUSDF + 1
C
C                 LES NC3CUB-1 COUCHES
                  DO 108 NBC=1,NC3CUB-1
                     AMPLI = CUB3RG**NBC
                     DO 106 NC=1,3
                        R = XYZDF(NC,I,J,K) * AMPLI
                        XYZCH(NC,NUSDF,NBC) = R
 106                 CONTINUE
 108              CONTINUE
C
C                 LA DERNIERE COUCHE POUR EVITER LES ERREURS D'ARRONDI
                  AMPLI = CUB3RG**NC3CUB
                  DO 109 NC=1,3
                     R = XYZDF(NC,I,J,K) * AMPLI
C                    3-CUBE CENTRE AU MIEUX
                     IF( ABS(R-CUB3DI2) .LE. DEPS ) THEN
                        IF( R .GE. 0 ) THEN
                           R = REAL( CUB3DI2 )
                        ELSE
                           R = REAL( -CUB3DI2 )
                        ENDIF
                     ENDIF
                     XYZCH(NC,NUSDF,NBC) = R
 109              CONTINUE
C
               ENDIF
C
 110        CONTINUE
 120     CONTINUE
 130  CONTINUE
C
 9000 WRITE(IMPRIM,*) 'NB SOMMETS des HEXAEDRES=',
     %                (NANOYAU+1)**3+NBSTCH*NC3CUB
      RETURN
      END
