      SUBROUTINE LEOPLS( LOPTRA )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   LECTURE DES OPTIONS COMMUNES DES TRACES DES LIGNES et SURFACES
C -----
C
C SORTIE :
C --------
C LOPTRA : L'OPTION DE TRACE
C          <=0 ANNIHILE LA DEMANDE DE TRACE SUITE A UNE ERREUR OU ABANDON
C          >=26 & <=40 PERMET UN TRACE IMMEDIAT
C          SINON RETOUR
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN LJLL UPMC & St PIERRE DU PERRAY DECEMBRE 2011
C ...................................................................012
      include"./incl/langue.inc"
      include"./incl/mxsuaf.inc"
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
      CHARACTER*72   KNMTS
      CHARACTER*10   NMTYOB, KTYOBJ
      CHARACTER*26   KNOM
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NUNITE(27)
C
C     LE TRAITEMENT SELON L'OPTION
      IF( LOPTRA .LT. 26 .OR. LOPTRA .GT. 40 ) GOTO 999
C
      GOTO(                          260, 270, 999, 290, 300,
     %      310, 320, 330, 340, 350, 360, 370, 380, 390, 400), LOPTRA-25
C
C     RE-DEFINITION des XYZ EXTREMES
 260  IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'DEFINIR l''HEXAEDRE de VISION :'
      ELSE
         KERR(1) = 'DEFINE the HEXAHEDRON of VISION :'
      ENDIF
10260 FORMAT(A1,' MIN: ',G13.6,T22,A1,' MAX: ',G13.6)
      WRITE(KERR(2),10260) 'X',COOEXT(1,1),'X',COOEXT(1,2)
      WRITE(KERR(3),10260) 'Y',COOEXT(2,1),'Y',COOEXT(2,2)
      WRITE(KERR(4),10260) 'Z',COOEXT(3,1),'Z',COOEXT(3,2)
      NBLGRC(NRERR) = 4
      CALL LERESU
C
      CALL INVITE( 107 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , COOEXT(1,1) )
      IF( NCVALS .LE. 0 ) GOTO 999
      CALL INVITE( 105 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , COOEXT(1,2) )
      IF( NCVALS .LE. 0 ) GOTO 999
C
      CALL INVITE( 116 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , COOEXT(2,1) )
      IF( NCVALS .LE. 0 ) GOTO 999
      CALL INVITE( 114 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , COOEXT(2,2) )
      IF( NCVALS .LE. 0 ) GOTO 999
C
      CALL INVITE( 124 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , COOEXT(3,1) )
      IF( NCVALS .LE. 0 ) GOTO 999
      CALL INVITE( 123 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , COOEXT(3,2) )
C
C     VISION TOTALE des XYZ EXTREMES
 270  NOTYVI = 0
      GOTO 999
C
C     PLANS SPECIAUX XY YZ XZ ...
 290  CALL LIMTCL( 'vuesplan' , NMTCL1 )
      IF( NMTCL1 .LE.  0 ) GOTO 999
      GOTO( 291, 292, 293, 294, 295, 296 ),NMTCL1
C
C     VUE du DESSUS
 291  AXOLON = 0.
      AXOLAT = 90.
      GOTO 305
C
C     VUE du DESSOUS
 292  AXOLON = 0.
      AXOLAT =-90.
      GOTO 305
C
C     VUE du GAUCHE
 293  AXOLON =-90.
      AXOLAT = 0.
      GOTO 305
C
C     VUE du DROITE
 294  AXOLON = 90.
      AXOLAT = 0.
      GOTO 305
C
C     VUE du FACE
 295  AXOLON = 0.
      AXOLAT = 0.
      GOTO 305
C
C     VUE de DERRIERE
 296  AXOLON = 180.
      AXOLAT = 0.
      GOTO 305
C
C     VISEE en LONGITUDE et LATITUDE (degres)
C     POINT VU LARGEUR HAUTEUR ARRIERE AVANT INCHANGES
 300  CALL INVITE( 27 )
      NCVALS = 0
      CALL LIRRSP( NCVALS , AXOLON )
      IF( NCVALS .LE. 0 ) GOTO 999
      AXOLON = MAX( -180.0 , MIN( 180.0 , AXOLON ) )
C
      CALL INVITE( 26 )
      NCVALS = 0
      CALL LIRRSP( NCVALS , AXOLAT )
      IF( NCVALS .LE. 0 ) GOTO 999
C
C     DEFINITION DE PTV ET OEIL A PARTIR DE AXOLON ET AXOLAT
 305  CALL LONLAT( AXOLON, AXOLAT )
C
C     L'AXONOMETRIE
      CALL MATAXO
      CALL ISOFENETRE( -AXOLAR, AXOLAR, -AXOHAU, AXOHAU )
      NOTYVI = 11
      GOTO 999
C
C     XYZ des POINT VU et OEIL
 310  CALL INVITE( 112 )
      NCVALS = 5
      CALL LIRXYZ( NCVALS , AXOPTV )
      IF( NCVALS .LE. 0 ) GOTO 999
C
      CALL INVITE( 111 )
      NCVALS = 5
      CALL LIRXYZ( NCVALS , AXOEIL )
      IF( NCVALS .LE. 0 ) GOTO 999
C
C     L'AXONOMETRIE
      GOTO 322
C
C     DEMI LARGEUR ET HAUTEUR DE LA SCENE
C     ===================================
 320  CALL INVITE( 25 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , AXOLAR )
      IF( NCVALS .LE. 0 ) GOTO 999
C
      CALL INVITE( 24 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , AXOHAU )
      IF( NCVALS .LE. 0 ) GOTO 999
C
C     L'AXONOMETRIE
 322  CALL MATAXO
      CALL ISOFENETRE( -AXOLAR, AXOLAR, -AXOHAU, AXOHAU )
      NOTYVI = 11
      GOTO 999
C
C     PLAN de DECOUPE ARRIERE AVANT DU POINT VISE
C     ===========================================
 330  CALL INVITE( 121 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , AXOARR )
      IF( NCVALS .LE. 0 ) GOTO 999
C
      CALL INVITE( 122 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , AXOAVA )
      IF( NCVALS .LE. 0 ) GOTO 999
C
C     VERIFICATION
      IF( AXOARR .GE. AXOAVA ) THEN
C        ERREUR : PAS DE TRONCATURE
         AXOARR = 0
         AXOAVA = 0
      ENDIF
      GOTO 999
C
C     LOUPE >1 OU <1 ou GROSSISSEMENT du TRACE
C     ========================================
 340  CALL INVITE( 22 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , GROSSI )
      IF( NCVALS .LE. 0 ) GOTO 999
      IF( GROSSI .LT. 0.01 ) GOTO 330
C
      IF( ABS(GROSSI-1.0) .LT. 0.01 ) GOTO 999
      AXOLAR = AXOLAR / GROSSI
      AXOHAU = AXOHAU / GROSSI
C
C     L'AXONOMETRIE
      CALL MATAXO
      IF( NDIMLI .LE. 2 ) THEN
         CALL ISOFENETRE( AXOPTV(1)-AXOLAR, AXOPTV(1)+AXOLAR,
     %                    AXOPTV(2)-AXOHAU, AXOPTV(2)+AXOHAU )
      ELSE
         CALL ISOFENETRE( -AXOLAR, AXOLAR, -AXOHAU, AXOHAU )
      ENDIF
      NOTYVI = 11
      GOTO 999
C
C     DEGRES DE ROTATION AUTOUR DE L'AXE Z
C     ====================================
C     LA LONGITUDE ET LATITUDE AVANT
 350  CALL RELOLA( AXOPTV, AXOEIL, DEGLON, DEGLAT )
C
 355  CALL INVITE( 5 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , R )
      IF( NCVALS .LE. 0 ) GOTO 999
      IF( R .LT. -360.0 .OR. R .GT. 360.0 ) GOTO 355
      DEGLON = DEGLON + R
 358  IF( DEGLON .GT. 360.0 ) THEN
         DEGLON = DEGLON - 360.0
         GOTO 358
      ENDIF
      IF( DEGLON .GT. 180.0 ) DEGLON = DEGLON - 360.0
C
C     PRISE EN COMPTE DE LA NOUVELLE LONGITUDE ET DE LA LATITUDE ANTERIEURE
      CALL LONLAT( DEGLON, DEGLAT )
C
C     FORMATION DE L'AXONOMETRIE
      CALL MATAXO
      CALL ISOFENETRE( -AXOLAR, AXOLAR, -AXOHAU, AXOHAU )
      NOTYVI = 11
C
C     SAUVEGARDE DU POINT SUR LE FICHIER FRAPPE
      IF( LANGAG .EQ. 0 ) THEN
       WRITE(NFFRAP,*) '{ POINT VU OEIL LARGEUR HAUTEUR ARRIERE AVANT }'
      ELSE
       WRITE(NFFRAP,*) '{ VIEW POINT; EYE; WIDTH; HEIGHT; BACK; AHEAD }'
      ENDIF
      WRITE(NFFRAP,*) (AXOPTV(I),'; ',I=1,3)
      WRITE(NFFRAP,*) (AXOEIL(I),'; ',I=1,3)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(NFFRAP,*) '32; { Demi LARGEUR et HAUTEUR de la SCENE }'
      ELSE
         WRITE(NFFRAP,*) '32; { HALF WIDTH and HEIGHT of the SCENE }'
      ENDIF
      WRITE(NFFRAP,*)  AXOLAR,'; ', AXOHAU,'; '
      WRITE(NFFRAP,*)  AXOARR,'; ', AXOAVA,'; '
      GOTO 999
C
C     CHANGER LA PALETTE DES COULEURS
C     ===============================
 360  CALL INVITE( 83 )
      N      = NOPACL
      NCVALS = 4
      CALL LIRENT( NCVALS , N )
      IF( NCVALS .LE. 0 ) GOTO 999
      CALL PALCDE( N )
      GOTO 999
C
C     DEFINIR LES COULEURS DES ARETES, FACES ... D'UN OBJET
C     =====================================================
 370  CALL INVITE( 23 )
      CALL LIMTCL( 'typ_objt' , N )
      IF( N .LE. 0 ) GOTO 999
      KTYOBJ = NMTYOB( N )
      I      = INDEX( KTYOBJ, ' ' )
      KNMTS  = 'NOM DE ' // KTYOBJ(1:I-2)
      GOTO( 371, 372, 373, 374, 375 ),N
 371  CALL INVITE( 51 )
      GOTO 379
 372  CALL INVITE( 40 )
      GOTO 379
 373  CALL INVITE( 42 )
      GOTO 379
 374  CALL INVITE( 60 )
      GOTO 379
 375  CALL INVITE( 45 )
C
 379  NCVALS = 0
      KNOM   = ' '
      CALL LIRCAR( NCVALS , KNOM )
      IF( NCVALS .LE. 0 ) GOTO 999
      I = INDEX( KNOM , ' ' )
      IF( I .LE. 1 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNOM
         KERR(1) = 'NOM INCORRECT'
         CALL LEREUR
         GOTO 360
      ENDIF
      KNMTS = '~>' // KTYOBJ(1:N-1) // '>' // KNOM(1:I-1) //
     %        '>TRACE'
      CALL MOTSTD( '~>>>TRACE' , KNMTS , N )
      IF( N .NE. 0 ) THEN
         NBLGRC(NRERR)  = 2
         KERR(1)(1:I-1) = KNOM(1:I-1)
         KERR(2)        = 'NOM INCORRECT'
         CALL LEREUR
      ENDIF
      GOTO 360
C
C     RETOUR AUX COULEURS PAR DEFAUT
C     ==============================
 380  CALL PALCDE( 12 )
      CALL COUDEF
C     LA COULEUR DU FOND EST REGENEREE
 381  CALL XVFOND( NCOFON )
      GOTO 999
C
C     COULEUR du FOND
C     ===============
 390  CALL INVITE( 17 )
      CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 999
      ELSE IF( I .EQ. 0 ) THEN
C        LA COULEUR NOIRE
         NCOFON = 0
      ELSE
         NCOFON = N1COEL + I
      ENDIF
      GOTO 381
C
C     QUALITE DU MAILLAGE
C     ===================
 400  IF( LCRITR .EQ. 0 ) THEN
         LCRITR = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE du MAILLAGE avec sa QUALITE'
         ELSE
            KERR(1) = 'DRAWING of the MESH WITH its QUALITY'
         ENDIF
         CALL LERESU
C        LA PALETTE ARC EN CIEL
         CALL PALCDE( 12 )
      ELSE
         LCRITR = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE du MAILLAGE SANS la QUALITE'
         ELSE
            KERR(1) = 'DRAWING of the MESH WITHOUT its QUALITY'
         ENDIF
         CALL LERESU
C        LA PALETTE DES GRIS
         CALL PALCDE( 10 )
      ENDIF
C
 999  RETURN
      END
