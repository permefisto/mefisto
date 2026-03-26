      SUBROUTINE PIGNON

      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
      include"./incl/a___xyzsommet.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))


C     CADRE DE L'EXTENSION
      XL= 46.5 + 186 + 155 + 186 + 45.
      XH = XL / 2
      COOEXT(1,1) = -10.
      COOEXT(1,2) = 10. +  XL

      COOEXT(2,1) = -10.
      YH = 390.
      COOEXT(2,2) = 10. + YH

      COOEXT(3,1) = 0.
      COOEXT(3,2) = 0.

      LORBITE = 1
      NORBITE = 0
      NOTYEV  = 0

      XOBMIN = COOEXT(1,1)
      XOBMAX = COOEXT(1,2)
      YOBMIN = COOEXT(2,1)
      YOBMAX = COOEXT(2,2)

      CALL ZOOM2D0( NOTYEV )

 10   CALL ZOOM2D1( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 9000

C     TRACE DES AXES
      CALL TRAXE2

C     LE PIGNON
      CALL XVEPAISSEUR( 5 )
      YG = 263.
      YD = 263.5
      CALL TRAIT2D( NCROUG,  0., 0.,  XL, 0. )
      CALL TRAIT2D( NCROUG,  XL, 0.,  XL, YD )
      CALL TRAIT2D( NCROUG,  XL, YD,  XH, YH )
      CALL TRAIT2D( NCROUG,  XH, YH,  0., YG )
      CALL TRAIT2D( NCROUG,  0., YG,  0., 0. )

C     LA BAIE GAUCHE
      YB   = 210.
      XBG1 = 46.5
      XBG2 = XBG1 + 186.
      CALL TRAIT2D( NCROUG, XBG1, 0.,  XBG2, 0. )
      CALL TRAIT2D( NCROUG, XBG2, 0.,  XBG2, YB )
      CALL TRAIT2D( NCROUG, XBG2, YB,  XBG1, YB )
      CALL TRAIT2D( NCROUG, XBG1, YB,  XBG1, 0. )

C     LA BAIE DROITE
      XBD1 = XBG2 + 155.
      XBD2 = XBD1 + 186.
      CALL TRAIT2D( NCROUG, XBD1, 0.,  XBD2, 0. )
      CALL TRAIT2D( NCROUG, XBD2, 0.,  XBD2, YB )
      CALL TRAIT2D( NCROUG, XBD2, YB,  XBD1, YB )
      CALL TRAIT2D( NCROUG, XBD1, YB,  XBD1, 0. )

C     UN ASSEMBLAGE DE 2 CARREAUX
      CALL XVEPAISSEUR( 1 )
      C3 = 30.75
      C1 = 10.25

C     AU DESSUS DE LA DIAGONALE
      X0 = 0.
      X  = X0
      Y0 = 0.
      Y  = Y0
      DO K = 1, 9

         DO L=1,22

            IF( X .GE. -C3 ) THEN

               IF( X .GE. XBG1 .AND.  X .LE. XBG2-C3-C1  .AND.
     %             Y .GE. 0    .AND.  Y .LE. YB-C3 ) GOTO 20

               IF( X .GE. XBD1 .AND.  X .LE. XBD2-C3-C1  .AND.
     %             Y .GE. 0    .AND.  Y .LE. YB-C3 ) GOTO 20

C           LE CARREAU 30x30
            CALL TRAIT2D( NCBLEU,  X,    Y,     X+C3, Y )
            CALL TRAIT2D( NCBLEU,  X+C3, Y,     X+C3, Y+C3 )
            CALL TRAIT2D( NCBLEU,  X+C3, Y+C3,  X,    Y+C3 )
            CALL TRAIT2D( NCBLEU,  X,    Y+C3,  X,    Y )

C           LE CARREAU 10x10
            CALL TRAIT2D( NCBLEU,  X+C3,    Y,     X+C3+C1, Y    )
            CALL TRAIT2D( NCBLEU,  X+C3+C1, Y,     X+C3+C1, Y+C1 )
            CALL TRAIT2D( NCBLEU,  X+C3+C1, Y+C1,  X+C3,    Y+C1 )
            CALL TRAIT2D( NCBLEU,  X+C3,    Y+C1,  X+C3,    Y    )

            ENDIF

C           PASSAGE AU CARREAU SUIVANT
 20         X = X + C3
            Y = Y + C1
         ENDDO

         X = X0 - K * (C3 + C1)
         Y = Y0 + K * (C3 - C1)
      ENDDO


C     AU DESSOUS DE LA DIAGONALE
      X0 = C3 + C1
      X  = X0
      Y0 = -C3 + C1
      Y  = Y0
      DO K = 1, 6

         DO L=1,30

            IF( X .LE. XL .AND. Y .GE. -C3 ) THEN

               IF( X .GE. XBG1 .AND.  X .LE. XBG2-C3-C1  .AND.
     %             Y .GE. 0    .AND.  Y .LE. YB-C3 ) GOTO 30

               IF( X .GE. XBD1 .AND.  X .LE. XBD2-C3-C1  .AND.
     %             Y .GE. 0    .AND.  Y .LE. YB-C3 ) GOTO 30

C           LE CARREAU 30x30
            CALL TRAIT2D( NCBLEU,  X,    Y,     X+C3, Y )
            CALL TRAIT2D( NCBLEU,  X+C3, Y,     X+C3, Y+C3 )
            CALL TRAIT2D( NCBLEU,  X+C3, Y+C3,  X,    Y+C3 )
            CALL TRAIT2D( NCBLEU,  X,    Y+C3,  X,    Y )

C           LE CARREAU 10x10
            CALL TRAIT2D( NCBLEU,  X+C3,    Y,     X+C3+C1, Y    )
            CALL TRAIT2D( NCBLEU,  X+C3+C1, Y,     X+C3+C1, Y+C1 )
            CALL TRAIT2D( NCBLEU,  X+C3+C1, Y+C1,  X+C3,    Y+C1 )
            CALL TRAIT2D( NCBLEU,  X+C3,    Y+C1,  X+C3,    Y    )

            ENDIF

C           PASSAGE AU CARREAU SUIVANT
 30         X = X + C3
            Y = Y + C1
         ENDDO

         X = X0 + K * (C3 + C1)
         Y = Y0 - K * (C3 - C1)
      ENDDO

C     LE TITRE
      CALL TRFINS( 'TRAVERTIN SUR LE PIGNON' )
C
C     REPRISE DE TRANSLATION ZOOM
      IF( LORBITE .GT. 0 ) GOTO 10
C
 9000 RETURN
      END
