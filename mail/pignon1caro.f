      SUBROUTINE PIGNON1CARO

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

      CALL XVEPAISSEUR( 5 )

C     LE PIGNON
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

      CALL XVEPAISSEUR( 1 )

C     UN CARREAU
      C3 = 31.
      C1 = 0.

C     DECALAGE Y POUR ALIGNER CARREAUX ET HAUT DES 2 BAIES VITREES
      YD = NINT(YB / C3) * C3 - YB

      DO K = 0, 12

         Y = -YD + K * C3

         DO 20 L = 0, 20

            X = -C3/2 + L * C3

            IF( X .GE. XBG1 .AND.  X .LE. XBG2-C3  .AND.
     %          Y .LE. YB-C3 ) GOTO 20

            IF( X .GE. XBD1 .AND.  X .LE. XBD2-C3  .AND.
     %          Y .LE. YB-C3 ) GOTO 20

C           LE CARREAU 30x30
            CALL TRAIT2D( NCBLEU,  X,    Y,     X+C3, Y )
            CALL TRAIT2D( NCBLEU,  X+C3, Y,     X+C3, Y+C3 )
            CALL TRAIT2D( NCBLEU,  X+C3, Y+C3,  X,    Y+C3 )
            CALL TRAIT2D( NCBLEU,  X,    Y+C3,  X,    Y )

 20      ENDDO

      ENDDO


C     LE TITRE
      CALL TRFINS( 'CARREAUX 31x31 de TRAVERTIN sur le PIGNON' )
C
C     REPRISE DE TRANSLATION ZOOM
      IF( LORBITE .GT. 0 ) GOTO 10
C
 9000 RETURN
      END
