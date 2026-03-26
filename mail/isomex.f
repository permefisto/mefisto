      SUBROUTINE ISOMEX( MNDFTR, DMATRI, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CONSTRUIRE LA MATRICE D'ISOMETRIE A PARTIR DE LA DEFINITION
C ----- ( ISOMETRIE=TRANSLATION ROTATION DILATATION OU SYMETRIES )
C
C ENTREES :
C ---------
C MNDFTR : ADRESSE MCN DU TABLEAU ~>TRANSFO>>DEFINITION
C
C SORTIES :
C ---------
C DMATRI : LA MATRICE D'ISOMETRIE  ( DOUBLE PRECISION(3,4) )
C IERR   : 0 SI PAS D'ERREUR , >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/ntmnlt.inc"
      include"./incl/a_point__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a_transfo__definition.inc"
C
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (RMCN(1),MCN(1))
      CHARACTER*24      KNOM1
      REAL              XYZ(3,8),COEPLA(4)
      INTEGER           NBPTTY(6)
      DOUBLE PRECISION  A,B,C,D,D1,D2,D3,DD(3),S,AXE1(3,3),AXE2(3,3),
     %                  DMATRI(3,4),DBLE
      EQUIVALENCE      (D1,DD(1)),(D2,DD(2)),(D3,DD(3))
C
      DATA              NBPTTY/0,8,6,3,2,1/
C
C     LE TYPE DE L'ISOMETRIE
      NUTYIS = MCN( MNDFTR + WUTYIS )
      IF( NUTYIS .LT. 0 .OR. NUTYIS .GT. 7 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'isomex: ISOMETRIE NON DEFINIE'
         ELSE
            KERR(1) = 'isomex: ISOMETRY NOT DEFINED'
         ENDIF
         CALL LEREUR
         GOTO 9999
      ENDIF
C
C     MISE A ZERO DU VECTEUR ET IDENTITE DE LA MATRICE D'ISOMETRIE
      DO 26 I=1,4
         DO 24 J=1,3
            DMATRI( J , I ) = 0D0
 24      CONTINUE
 26   CONTINUE
      DO 28 I=1,3
         DMATRI( I , I ) = 1D0
 28   CONTINUE
C
C     INITIALISATION SELON LE TYPE
C     ============================
      IF( NUTYIS .EQ. 0 ) THEN
C        TRANSLATION D'UN VECTEUR DE R3
         DO 70 I=1,3
            DMATRI( I , 4 ) = RMCN( MNDFTR + WRAVEC - 1 + I )
 70      CONTINUE
         GOTO 1900
      ENDIF
C
      IF( NUTYIS .EQ. 1 ) THEN
         IF( MCN(MNDFTR+WUP1AX) .LE. 0 ) THEN
C           PAS DE POINT DE ROTATION
            NROT2D = 0
            NBPOIN = 0
         ELSE IF( MCN(MNDFTR+WUP2AX) .LE. 0 ) THEN
C           1 SEUL POINT DE ROTATION => 2D
            NROT2D = 2
            NBPOIN = 1
         ELSE
C           2 POINTS DE ROTATION => 3D
            NROT2D = 3
            NBPOIN = 2
         ENDIF
C
C        ANGLE DE ROTATION
         ANGLE = RMCN(MNDFTR+WANGLE)
         IF( ANGLE .EQ. 0. ) THEN
C           PAS DE ROTATION
            NROT2D = 0
            NBPOIN = 0
         ENDIF
C
C        TRANSLATION
         IF( RMCN(MNDFTR+WRANSX).NE.0. .OR.
     %       RMCN(MNDFTR+WRANSY).NE.0. .OR.
     %       RMCN(MNDFTR+WRANSZ).NE.0. ) THEN
            NTRANS = 1
         ELSE
            NTRANS = 0
         ENDIF
C
C        DILATATION
         IF( RMCN(MNDFTR+WILATX) .NE. 0.  .OR.
     %       RMCN(MNDFTR+WILATY) .NE. 0.  .OR.
     %       RMCN(MNDFTR+WILATZ) .NE. 0. ) THEN
            NDILAT = 1
         ELSE
            NDILAT = 0
         ENDIF
C
      ELSE
C
C        TYPE D'ISOMETRIE DIFFERENT DE 1
         NTRANS = 0
         NROT2D = 0
         NDILAT = 0
         NBPOIN = NBPTTY( NUTYIS )
      ENDIF
C
C     CALCUL DE LA MATRICE D'ISOMETRIE
C     ================================
C     VERIFICATION ET RECHERCHE DES SOMMETS DES NBPOIN POINTS
      DO 1020 I=1,NBPOIN
C        LE NUMERO DU POINT
         N = MCN( MNDFTR + WUTYIS + I )
         IF( N .LE. 0 ) THEN
C           LE POINT I N'EST PAS INITIALISE
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'isomex: UN DES POINTS N''EST PAS INITIALISE'
            ELSE
               KERR(1) = 'isomex: ONE of POINTS IS NOT INITIALIZED'
            ENDIF
            CALL LEREUR
            GOTO 9999
         ENDIF
C        LES COORDONNEES DU I-EME POINT
         CALL XYZPOI( N, MN, IERR )
         IF( IERR .NE. 0 ) GOTO 9999
         XYZ(1,I) = RMCN( MN )
         XYZ(2,I) = RMCN( MN + 1 )
         XYZ(3,I) = RMCN( MN + 2 )
 1020 CONTINUE
C
C     GENERATION DE LA MATRICE SELON LE TYPE
      IF( NTRANS .NE. 0 .OR. NROT2D .NE. 0 .OR. NDILAT .NE. 0 ) THEN
C
C        LE TYPE 1 DEFORMATION AVEC POINTS ET REELS
C        ==========================================
         IF( NTRANS .NE. 0 ) THEN
C           TRANSLATION A L'AIDE DE 3 REELS
            DO 1030 I=1,3
               DMATRI( I , 4 ) = RMCN( MNDFTR + WRANSX - 1 + I )
 1030       CONTINUE
         ELSE
            DO 1040 I=1,3
               DMATRI( I , 4 ) = 0D0
 1040       CONTINUE
         ENDIF
C
         IF( NROT2D .NE. 0 ) THEN
C           ROTATION . CONVERSION DE L'ANGLE DEGRES EN RADIANS
            D  = RMCN( MNDFTR+WANGLE ) * 3.14159265358979312D0 / 180.D0
            D1 = COS( D )
            D2 = SIN( D )
C
            DMATRI(1,1) =  D1
            DMATRI(2,1) =  D2
            DMATRI(3,1) =  0D0
C
            DMATRI(1,2) = -D2
            DMATRI(2,2) =  D1
            DMATRI(3,2) =  0D0
C
            DMATRI(1,3) =  0D0
            DMATRI(2,3) =  0D0
            DMATRI(3,3) =  1D0
C
            IF( NROT2D .EQ. 2 ) THEN
C
C              ROTATION 2D CALCUL DU VECTEUR
C              -----------------------------
               D1 = 1.D0 - D1
               DMATRI(1,4) = DMATRI(1,4) + D1 * XYZ(1,1) + D2 * XYZ(2,1)
               DMATRI(2,4) = DMATRI(2,4) + D1 * XYZ(2,1) - D2 * XYZ(1,1)
C
            ELSE IF( NROT2D .EQ. 3 ) THEN
C
C              ROTATION 3D AUTOUR D'UN AXE DEFINI PAR 2 POINTS
C              -----------------------------------------------
C              AXE3 DU REPERE 1 LIE AU 1-ER SOLIDE AVANT ROTATION
               DO 1050 I=1,3
                  AXE1(I,3) = XYZ(I,2) - XYZ(I,1)
 1050          CONTINUE
C              NORMALISATION A 1. EN LONGUEUR
               CALL NORME1( 3 , AXE1(1,3), IERR )
               IF( IERR .NE. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'isomex: AXE 1 REDUIT A UN POINT'
                  ELSE
                     KERR(1) = 'isomex: AXIS 1 REDUCED TO ONE POINT'
                  ENDIF
                  CALL LEREUR
                  GOTO 9999
               ENDIF
C              AXE2 DU REPERE 1 LIE AU 1-ER SOLIDE AVANT ROTATION
               IF( AXE1(1,3).EQ.0D0 .AND. AXE1(2,3).EQ.0D0 ) THEN
C                 AXE PARALLELE A OZ
                  AXE1(1,2) = 0D0
                  AXE1(2,2) = 1D0
                  AXE1(3,2) = 0D0
               ELSE
                  AXE1(1,2) = -AXE1(2,3)
                  AXE1(2,2) =  AXE1(1,3)
                  AXE1(3,2) = 0D0
                  CALL NORME1( 3 , AXE1(1,2), IERR )
                  IF( IERR .NE. 0 ) THEN
                     NBLGRC(NRERR) = 1
                     IF( LANGAG .EQ. 0 ) THEN
                        KERR(1) = 'isomex: AXE 2 REDUIT A UN POINT'
                     ELSE
                        KERR(1) = 'isomex: AXIS 2 REDUCED TO ONE POINT'
                     ENDIF
                     CALL LEREUR
                     GOTO 9999
                  ENDIF
               ENDIF
C              AXE1 DU REPERE 1 = AXE2 PRODUIT VECTORIEL AXE3
               CALL PROVEC( AXE1(1,2) , AXE1(1,3) , AXE1(1,1) )
               CALL NORME1( 3 , AXE1(1,1), IERR )
               IF( IERR .NE. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'isomex: AXE 3 REDUIT A UN POINT'
                  ELSE
                     KERR(1) = 'isomex: AXIS 3 REDUCED TO ONE POINT'
                  ENDIF
                  CALL LEREUR
                  GOTO 9999
               ENDIF
C
C              DMATRI = AXE1 * DMATRI * TRANSPOSEE( AXE1 )
               CALL ATB0D( 3 , 3 , 3 , DMATRI , AXE1 , AXE2   )
               CALL AB0D ( 3 , 3 , 3 , AXE1   , AXE2 , DMATRI )
C
C              SI DILATATION : DMATRI = DMATRI * DILATATION
               IF( NDILAT .NE. 0 ) THEN
                  DO 1070 I=1,3
                     D = RMCN( MNDFTR + WILATX - 1 + I )
                     DO 1060 J=1,3
                        DMATRI( J , I ) = D * DMATRI( J , I )
 1060                CONTINUE
 1070             CONTINUE
               ENDIF
C
C              VECTOR = VECTOR + POINT1 - DMATRI * POINT1
               DO 1090 I=1,3
                  D = 0D0
                  DO 1080 J=1,3
                     D = D + DMATRI(I,J) * XYZ(J,1)
 1080             CONTINUE
                  DMATRI(I,4) = DMATRI(I,4) + XYZ(I,1) - D
 1090          CONTINUE
            ENDIF
C
         ELSE
C
C           PAS DE ROTATION MAIS UNE DILATATION
            DMATRI(1,1) = RMCN( MNDFTR + WILATX )
            DMATRI(2,2) = RMCN( MNDFTR + WILATY )
            DMATRI(3,3) = RMCN( MNDFTR + WILATZ )
         ENDIF
      ENDIF
C
C     LES TYPES 2 A 7 D'ISOMETRIES
C     ============================
      GOTO( 1900 , 1100 , 1100 , 1300 , 1400 , 1600, 1700 ) , NUTYIS
C
C     SI DEFORMATION NBPOIN = 8
C     SI DEPLACEMENT NBPOIN = 6
C     LES 3 POINTS DE CHAQUE REPERE SONT-ILS ALIGNES ?
 1100 NBPOI2 = NBPOIN / 2
      K      = 3
      CALL EQPLAN( XYZ , COEPLA , J )
      N      = 0
      IF( J .GT. 0 ) GOTO 1105
      CALL EQPLAN( XYZ(1,1+NBPOI2) , COEPLA , J )
      N      = NBPOI2
      IF( J .LE. 0 ) GOTO 1110
C
C     ERREUR : 3 POINTS COLINEAIRES => REPERE NON DEFINI
 1105 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'isomex: 3 POINTS CONFONDUS OU ALIGNES'
      ELSE
         KERR(1) = 'isomex: 3 POINTS ON A SAME STRAIGHT LINE'
      ENDIF
      CALL LEREUR
      DO 1108 I=1,K
C        DECODAGE DU NOM DU POINT
         CALL NMOBNU( 'POINT' , MCN( MNDFTR+WUTYIS+I+N ) , KNOM1 )
         WRITE(IMPRIM,11105) KNOM1,(XYZ(J,N+I),J=1,3)
11105 FORMAT(' POINT ',A,' :X=',G15.7,' Y=',G15.7,' Z=',G15.7)
 1108 CONTINUE
      GOTO 9999
C
C     LES COORDONNEES SONT TRANSFORMEES EN DOUBLE PRECISION
 1110 DO 1120 I=1,3
         DO 1115 J=1,3
            DMATRI(J,I) = XYZ(J,I)
            AXE1  (J,I) = XYZ(J,I+NBPOI2)
 1115    CONTINUE
 1120 CONTINUE
C
C     LES AXE2 DU REPERE APRES ISOMETRIE
      CALL P3AXE3( AXE1, AXE2, IERR )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'isomex: AXES FINAUX INCORRECTS'
         ELSE
            KERR(1) = 'isomex: FINAL AXES NOT CORRECT'
         ENDIF
         CALL LEREUR
         GOTO 9999
      ENDIF
C
C     LES AXE1 DU REPERE AVANT ISOMETRIE
      CALL P3AXE3( DMATRI, AXE1, IERR )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'isomex: AXES INITIAUX INCORRECTS'
         ELSE
            KERR(1) = 'isomex: INITIAL AXES NOT CORRECT'
         ENDIF
         CALL LEREUR
         GOTO 9999
      ENDIF
C
C     PRISE EN COMPTE EVENTUELLE DE LA DILATATION
      IF( NUTYIS .EQ. 2 ) THEN
         DO 1140 I=1,3
            D1 = 0D0
            D2 = 0D0
            N  = I + 1
            DO 1130 J=1,3
               D1 = D1 + ( XYZ(J,N)        - XYZ(J,1)         ) ** 2
               D2 = D2 + ( XYZ(J,N+NBPOI2) - XYZ(J,1+NBPOI2)  ) ** 2
 1130       CONTINUE
C
C           SI LA DISTANCE D1 EST NULLE LA DILATATION VAUT 1.
            IF( D1 .LE. 0D0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') I
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = KERR(MXLGER)(1:4)//
     %                   '-EME DILATATION ERRONEE FIXEE A 1.'
               ELSE
                  KERR(1) = KERR(MXLGER)(1:4)//
     %                   '-TH DILATATION WITH ERROR FIXED TO 1.'
               ENDIF
               CALL LEREUR
               DMATRI(I,4) = 1D0
            ELSE
C              DILATATION CORRECTE
               DMATRI(I,4) = SQRT( D2 / D1 )
            ENDIF
 1140    CONTINUE
C
C        LE PRODUIT AXE2 * DILATATION
         DO 1160 I=1,3
            DO 1150 J=1,3
               AXE2(J,I) = AXE2(J,I) * DMATRI(I,4)
 1150       CONTINUE
 1160    CONTINUE
      ENDIF
C
C     DMATRI = ( R2 ) * (DILATATION) * TRANSPOSEE( R1 )
      CALL ATB0D( 3 , 3 , 3 , AXE2 , AXE1 , DMATRI )
C
C     VECTOR = ORIGINE2 -DMATRI * ORIGINE1
      DO 1180 I=1,3
         D1 = 0D0
         DO 1170 J=1,3
            D1 = D1 + DMATRI(I,J) * XYZ(J,1)
 1170    CONTINUE
         DMATRI(I,4) = XYZ(I,NBPOI2+1) - D1
 1180 CONTINUE
      GOTO 1900
C
C     SYMETRIE/PLAN
C     =============
C     EQUATION DU PLAN
 1300 CALL EQPLAN( XYZ , COEPLA , J )
      IF( J .GT. 0 ) THEN
C        ERREUR POINTS CONFONDUS OU ALIGNES
         N = 0
         K = 3
         GOTO 1105
      ENDIF
C
C     LA MATRICE ET VECTEUR DE SYMETRIE/PLAN
      D1 = DBLE( COEPLA(1) ) ** 2
      D2 = DBLE( COEPLA(2) ) ** 2
      D3 = DBLE( COEPLA(3) ) ** 2
      D  = -2D0 / ( D1 + D2 + D3 )
C
      DMATRI(1,1) = 1D0 + D * D1
      DMATRI(2,2) = 1D0 + D * D2
      DMATRI(3,3) = 1D0 + D * D3
C
      D1 = D * COEPLA(1)
      D2 = D * COEPLA(2)
      D3 = D * COEPLA(3)
C
      S           = D1 * COEPLA(2)
      DMATRI(1,2) = S
      DMATRI(2,1) = S
      S           = D1 * COEPLA(3)
      DMATRI(1,3) = S
      DMATRI(3,1) = S
      S           = D2 * COEPLA(3)
      DMATRI(2,3) = S
      DMATRI(3,2) = S
C
      D  = COEPLA(4)
      DMATRI(1,4) = D * D1
      DMATRI(2,4) = D * D2
      DMATRI(3,4) = D * D3
      GOTO 1900
C
C     SYMETRIE/DROITE DANS R ** 3
C     ===========================
 1400 D1 = XYZ(1,1)
      D2 = XYZ(2,1)
      D3 = XYZ(3,1)
      A  = XYZ(1,2) - D1
      B  = XYZ(2,2) - D2
      C  = XYZ(3,2) - D3
      D  = A * A + B * B + C * C
C
      IF( D .EQ. 0D0 ) THEN
C        ERREUR : POINTS CONFONDUS
         N = 0
         K = 2
         GOTO 1105
      ENDIF
C
      D = 2D0 / D
      S = D * A
C
      DMATRI(1,1) = S * A - 1D0
      DMATRI(1,2) = S * B
      DMATRI(2,1) = S * B
      DMATRI(1,3) = S * C
      DMATRI(3,1) = S * C
C
      S = D * B
      DMATRI(2,2) = S * B - 1D0
      DMATRI(2,3) = S * C
      DMATRI(3,2) = S * C
      DMATRI(3,3) = D * C * C - 1D0
C
C     VECTOR = -DMATRI * P1 + P1
      DO 1420 I=1,3
C        DD EN EQUIVALENCE AVEC D1,D2,D3
         S = DD(I)
         DO 1410 J=1,3
            S = S - DMATRI(I,J) * DD(J)
 1410    CONTINUE
         DMATRI(I,4) = S
 1420 CONTINUE
      GOTO 1900
C
C     SYMETRIE/POINT   DANS R ** 3
C     ==============
 1600 DO 1620 I=1,3
         DMATRI(I,4) = 2D0 * XYZ(I,1)
         DO 1610 J=1,3
            DMATRI(I,J) = 0D0
 1610    CONTINUE
         DMATRI(I,I) = -1D0
 1620 CONTINUE
      GOTO 1900
C
C     DILATATION  DANS R ** 3
C     ==========
 1700 DO 1720 I=1,3
         DO 1710 J=1,4
            DMATRI(I,J) = 0D0
 1710    CONTINUE
         DMATRI(I,I) = DBLE( RMCN(MNDFTR+WILA7X-1+I) )
 1720 CONTINUE
C
 1900 IERR = 0
      RETURN
C
C     ERREUR
 9999 IERR = 1
      END
