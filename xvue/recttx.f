      SUBROUTINE RECTTX( NUR, TEXTE, NETEXT, NATEXT )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACE DES LIGNES DE CARACTERES DU RECTANGLE DE NUMERO NUR
C -----           VERSION xvue et PS
C
C ENTREE :
C --------
C NUR    : LE NUMERO DU RECTANGLE A TRACER
C TEXTE  : LES LIGNES DE TEXTE A TRACER
C NETEXT : 0 PAS DE TRACE DES INTERLIGNES
C          1 TRACE DES INTERLIGNES SANS PREOCCUPATION DE NATEXT
C          2 TRACE DES INTERLIGNES
C NATEXT : LE NUMERO ASSOCIE A CHAQUE LIGNE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS        MAI 1994
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/epombr.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      CHARACTER*(*)   TEXTE(*)
      INTEGER         NATEXT(1:*)
      CHARACTER*4     KNO
      CHARACTER*120   KTXT
C
      IF( INTERA .LT. 3 ) RETURN
C
C     SI PAS DE LIGNES DE TEXTE: RETOUR
      IF( NBLGRC(NUR) .LE. 0 ) RETURN
C
C     LES DIMENSIONS DU RECTANGLE
C     LE NOMBRE MAXIMAL DE PIXELS OCCUPES PAR LES LIGNES DU TEXTE
      NBP = LAMXPXTXT( NBLGRC(NUR), TEXTE ) + 8
C     LE NOMBRE MAXIMAL DE CARACTERES NON BLANCS
      MDRECT(NUR) = MXCANB( NBLGRC(NUR), TEXTE )
      IF( NETEXT .EQ. 2 .AND. NUR .EQ. NRMENU ) THEN
C        5 CARACTERES DE PLUS POUR LE NUMERO OU MENU:
         MDRECT(NUR) = MDRECT(NUR) + 5
         NBP = NINT( REAL(NBP) * (MDRECT(NUR)+4.0) / MDRECT(NUR) )
      ENDIF
C     SI PAS DE CARACTERES NON BLANCS DANS TEXTE: RETOUR
      IF( MDRECT(NUR) .LE. 0 ) RETURN
C     LA LARGEUR DU RECTANGLE
      DXRECT(NUR) = MAX( DXRECT(NUR), ECARLR(NUR)*2 + NBP ) + 2
C     LA HAUTEUR DU RECTANGLE
      DYRECT(NUR) = DYLGRC(NUR) * NBLGRC(NUR)
C
C     SI LE MODE POSTSCRIPT EST EN COURS
      IF ( LASOPS .NE. 0 ) THEN
        IF ( LASOPS .EQ. 1 ) THEN
C         MISE EN SUSPEND DU POSTSCRIPT PENDANT L ECRITURE MENU
          LASOPS = -1
        ELSE
          IF ( LASOPS .EQ. 2 ) THEN
C           ORIENTATION VERS LE MENU CORRESPONDANT
            LASOPS = 3 + NUR
            LASOPS = MIN(10,LASOPS)
          ELSE
            LASOPS = 0
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'RECTTX : MAUVAISE VALEUR DE LASOPS'
               KERR(2) = '         ARRET DU DESSIN POSTSCRIPT'
            ELSE
               KERR(1) = 'RECTTX : BAD VALUE of LASOPS'
               KERR(2) = '         STOP of POSTSCRIPT DRAWING'
            ENDIF
            CALL LEREUR
          ENDIF
        ENDIF
        CALL XVPOSTSCRIPT( LASOPS )
      ENDIF
C
C     LE TRACE DU FOND
      CALL RECTTR( NUR )
C
C     LE TRACE DES LIGNES INTERMEDIAIRES
      CALL RECTLG( NUR, NETEXT, NATEXT )
C
C     LA COULEUR DES CARACTERES
      CALL XVCOULEUR( NKRECT( NUR ) )
C
C     LE TRACE DES LIGNES DU TEXTE
      IF( NBLGRC(NUR) .GT. 0 ) THEN
C
C        POSITION EN PIXELS DU TEXTE
         NX = XRECT( NUR ) + ECARLR( NUR )
         NY = YRECT( NUR ) - ECARLR( NUR ) - NBEPO4
C
         DO 10 I=1,NBLGRC( NUR )
C
C           TRACE DE LA LIGNE I  DANS LE RECTANGLE NUR
C           ------------------------------------------
            NY = NY + DYLGRC( NUR )
            IF( NUR .NE. NRMENU ) THEN
C
C              RECTANGLE DIFFERENT D'UN MENU
C              =============================
               IF( NUR .NE. NRLGLU ) THEN
C
C                 RECTANGLE NON MENU NO LIGNE LUE
C                 ===============================
                  KTXT = TEXTE(I)(1:MDRECT(NUR))
C
               ELSE
C
C                 RECTANGLE DES LIGNES LUES
C                 =========================
                  IF( LANGAG .EQ. 0 ) THEN
C                    FRANCAIS
                     KTXT = 'LIGNE LUE:' // TEXTE(I)
                  ELSE
C                    ENGLISH
                     KTXT = 'READ LINE:' // TEXTE(I)
                  ENDIF
               ENDIF
C
               L = NUDCNB( KTXT )
               CALL XVFTEXTE( KTXT(1:L), L, NX, NY+1 )
C
            ELSE
C
C              RECTANGLE MENU
C              ==============
C              LE NUMERO ASSOCIE A LA LIGNE
               KNO = '    '
               IF( NETEXT .EQ. 2 ) THEN
                  IF( NATEXT(I) .EQ. NATEXT(I+1) ) THEN
                     KNO = '    '
                  ELSE IF( NATEXT(I) .EQ. INMENU-1 ) THEN
                     KNO = '@ ; '
                  ELSE IF( NATEXT(I) .EQ. INMENU+1 ) THEN
                     KNO = '? ; '
                  ELSE IF( NATEXT(I) .EQ. INMENU ) THEN
C                    NATEXT NON INITIALISE
                     KNO = '    '
                  ELSE
                     WRITE( KNO(1:2) , '(I2)' , IOSTAT=L ) NATEXT(I)
                     IF( L .EQ. 0 ) THEN
                        KNO(3:3) = ';'
                     ELSE
                        KNO = '    '
                     ENDIF
                  ENDIF
               ENDIF
C
C              LA PREMIERE LIGNE DU MENU EST AUGMENTEE DE 'MENU:'
               IF( I .EQ. 1 ) THEN
                  L = NUDCNB( TEXTE(1) )
                  KTXT = 'MENU: ' // TEXTE(1)(1:L)
               ELSE
                  KTXT = KNO // TEXTE(I)(1:MDRECT(NUR))
               ENDIF
               L = NUDCNB( KTXT )
               CALL XVFTEXTE( KTXT(1:L), L, NX, NY+1 )
            ENDIF
 10      CONTINUE
      ENDIF
C
C     SI LE MODE POSTSCRIPT EST EN COURS
      IF ( LASOPS .NE. 0 ) THEN
C       RETOUR AU MODE POSTSCRIPT NORMAL
        IF ( ABS(LASOPS) .EQ. 1 ) THEN
          LASOPS = 1
        ELSE
          LASOPS = 2
        ENDIF
        CALL XVPOSTSCRIPT( LASOPS )
      ENDIF
C
      RETURN
      END
