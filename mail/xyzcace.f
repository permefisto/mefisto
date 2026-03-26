      SUBROUTINE XYZCACE( LCaouce, ARNOYAU, NANOYAU,
     %                    CACEDI,  CACERG,  NC2CAR, NBSTCH,
     %                    NUSTDF,  XYZDF,   XYZCH )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES 2 COORDONNEES DES SOMMETS D'UNE QUADRANGULATION
C -----    FORMEE PAR 2 DIFFERENCES FINIES ET NC2CAR COUCHES HOMOTHETIQUES
C          EN PROGRESSION GEOMETRIQUE
C
C ENTREES:
C --------
C LCaouce: 0 SI CARRE DEMANDE, 1 SI CERCLE DEMANDE
C ARNOYAU: LONGUEUR D'UNE ARETE DANS UNE DIRECTION DU NOYAU
C NANOYAU: NOMBRE D'ARETES DANS UNE DIRECTION DU NOYAU (DIFFERENCES FINIES)
C CACEDI : LARGEUR TOTALE DU CARRE ou CERCLE
C CACERG : PROGRESSION GEOMETRIQUE DES ARETES DES COUCHES
C NC2CAR : NOMBRE DE COUCHES HOMOTHETIQUES
C NBSTCH : NOMBRE DE SOMMETS D'UNE COUCHE
C
C SORTIES:
C --------
C NUSTDF : NO GLOBAL DES SOMMETS FRONTALIERS POUR LES DF DU NOYAU
C XYZDF  : 3 COORDONNEES DES SOMMETS DIFFERENCES FINIES  DU NOYAU
C XYZCH  : 3 COORDONNEES DES SOMMETS DES COUCHES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & LJUBLJANA SLOVENIE   Novembre 2011
C.......................................................................
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
C     CONSTRUCTION DU TABLEAU XYZ DU MAILLAGE EN QUADRANGLES
      REAL              XYZDF(3,0:NANOYAU,0:NANOYAU)
      INTEGER           NUSTDF( 0:NANOYAU,0:NANOYAU)
      REAL              XYZCH(3,1:NBSTCH,1:NC2CAR)
      DOUBLE PRECISION  CACEDI2, PIS2, EPS, ATAN
      DOUBLE PRECISION  RAYON, LONGITUDE, X, Y, ARETEN
      INTRINSIC         REAL, DBLE
      LOGICAL           NBAIMPAIR
      INTEGER           NBAMILIEU
C
C     LE RAYON DU NOYAU EST EGAL A SA LARGEUR/2
      RAYON = ( ARNOYAU * NANOYAU ) * 0.5D0
C
C     CALCUL DES 2 COORDONNEES DES SOMMETS DIFFERENCES FINIES
C     LE CENTRE DU NOYAU EST L'ORIGINE DES AXES XYZ
      X0 = REAL( - RAYON )
C
C     LES SOMMETS DIFFERENCES FINIES DU NOYAU ET DE LA FRONTIERE
      NBSTDF = (NANOYAU+1)**2
      IF( LCaouce .EQ. 0 ) THEN
C
C        CARRE DEMANDE
C        -------------
         DO J=0,NANOYAU
            DO I=0,NANOYAU
C
C              LES 3 COORDONNEES DU SOMMET DIFFERENCES FINIES DU NOYAU
               XYZDF(1,I,J) = X0 + I * ARNOYAU
               XYZDF(2,I,J) = X0 + J * ARNOYAU
               XYZDF(3,I,J) = 0.
C
C              LE SOMMET EST IL SUR LA FRONTIERE
C              DES DIFFERENCES FINIES DU NOYAU?
               IF( I.EQ.0 .OR. I.EQ.NANOYAU .OR.
     %             J.EQ.0 .OR. J.EQ.NANOYAU ) THEN
C                 SOMMET SUR LA FRONTIERE
                  NBSTDF = NBSTDF + 1
C                 SON NUMERO GLOBAL
                  NUSTDF(I,J) = NBSTDF
               ENDIF
            ENDDO
         ENDDO
C
      ELSE
C
C        CERCLE DEMANDE
C        --------------
         PIS2   = ATAN( 1D0 ) * 2D0
         ARETEN = DBLE( ARNOYAU )
         EPS    = ARETEN * 1D-5
         NBAIMPAIR = MOD(NANOYAU,2) .EQ. 1
         NBAMILIEU = (NANOYAU+1)/2
C
         DO J=0,NANOYAU
            NY = ABS( NBAMILIEU - J )
            IF( NBAIMPAIR .AND. (J .GE. NBAMILIEU) ) NY=NY+1
C
            DO I=0,NANOYAU
               NX = ABS( NBAMILIEU - I )
               IF( NBAIMPAIR .AND. (I .GE. NBAMILIEU) ) NX=NX+1
C
C              LE MAXIMUM DES 2 NUMEROS EST LE NUMERO DU RAYON
               NR = MAX(NY,NX)
C              LE RAYON DU POINT DANS LE CERCLE
               RAYON = NR * ARNOYAU
C
C              LES 3 COORDONNEES DU SOMMET DIFFERENCES FINIES DU NOYAU
               X = X0 + I * ARETEN
               Y = X0 + J * ARETEN
C
C              LONGITUDE DU SOMMET
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
C              LES XYZ DU SOMMET SUR LE CERCLE DE RAYON NR
               XYZDF(1,I,J) = REAL( RAYON * COS( LONGITUDE ) )
               XYZDF(2,I,J) = REAL( RAYON * SIN( LONGITUDE ) )
               XYZDF(3,I,J) = 0.
C
C              LE SOMMET EST IL SUR LA FRONTIERE
C              DES DIFFERENCES FINIES DU NOYAU?
               IF( I.EQ.0 .OR. I.EQ.NANOYAU   .OR.
     %             J.EQ.0 .OR. J.EQ.NANOYAU ) THEN
C                  SOMMET SUR LA FRONTIERE
                   NBSTDF = NBSTDF + 1
C                  SON NUMERO GLOBAL
                   NUSTDF(I,J) = NBSTDF
                ENDIF
             ENDDO
          ENDDO
      ENDIF
C
C     TRAITEMENT DES COUCHES HOMOTHETIQUES A PARTIR DE LA SURFACE DU NOYAU
C     ====================================================================
C     NOMBRE DE SOMMETS FRONTALIERS OU D'UNE COUCHE
      WRITE(IMPRIM,*) 'NB SOMMETS FRONTALIERS DU NOYAU=', NBSTCH
C
C     CALCUL DES 3 COORDONNEES DES SOMMETS DES NC2CAR COUCHES HOMOTHETIQUES
C     => NBSTCH SOMMETS NOUVEAUX PAR COUCHE
      IF( NC2CAR .LE. 0 ) GOTO 9000
      CACEDI2 = CACEDI / 2D0
      DEPS  = 0.001*CACEDI
      NUSDF = 0
      DO 120 J=0,NANOYAU
         DO 110 I=0,NANOYAU
C
C           LE SOMMET EST IL SUR LA FRONTIERE DES DIFFERENCES FINIES?
            IF( I.EQ.0 .OR. I.EQ.NANOYAU .OR.
     %          J.EQ.0 .OR. J.EQ.NANOYAU ) THEN
C              SOMMET SUR LA FRONTIERE
               NUSDF = NUSDF + 1
C
C              LES NC2CAR-1 COUCHES
               DO 108 NBC=1,NC2CAR-1
                  AMPLI = CACERG**NBC
                  DO 106 NC=1,3
                     R = XYZDF(NC,I,J) * AMPLI
                     XYZCH(NC,NUSDF,NBC) = R
 106              CONTINUE
 108           CONTINUE
C
C              LA DERNIERE COUCHE POUR EVITER LES ERREURS D'ARRONDI
               AMPLI = CACERG**NC2CAR
               DO 109 NC=1,3
                  R = XYZDF(NC,I,J) * AMPLI
C                 3-CARRE CENTRE AU MIEUX
                  IF( ABS(R-CACEDI2) .LE. DEPS ) THEN
                     IF( R .GE. 0 ) THEN
                        R = REAL( CACEDI2 )
                     ELSE
                        R = REAL( -CACEDI2 )
                     ENDIF
                  ENDIF
                  XYZCH(NC,NUSDF,NBC) = R
 109           CONTINUE
C
            ENDIF
C
 110     CONTINUE
 120  CONTINUE
C
 9000 WRITE(IMPRIM,*) 'NB SOMMETS des QUADRANGLES=',
     %                (NANOYAU+1)**2+NBSTCH*NC2CAR
      RETURN
      END
