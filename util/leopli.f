      SUBROUTINE LEOPLI( LOPTRA )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   LECTURE DES OPTIONS DE TRACE DES ARETES D'UNE SURFACE OU VOLUME
C -----
C
C SORTIE :
C --------
C LOPTRA : L'OPTION DE TRACE
C         <=0 : RETOUR SANS TRACE
C          >0 : TRACE  IMMEDIAT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS        MARS 1991
C ...................................................................012
      PARAMETER  (MXSUAR=16)
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
C     L'OPTION A TRAITER
 1    CALL LIMTCL( 'opt_lign' , LOPTRA )
      IF( LOPTRA .LE. 0  ) RETURN
      IF( LOPTRA .EQ. 49 ) GOTO 490
      IF( LOPTRA .EQ. 50 ) GOTO 5000
      IF( LOPTRA .EQ. 75 ) GOTO 750
      IF( LOPTRA .EQ. 76 ) GOTO 760
      IF( LOPTRA .EQ. 90 ) RETURN
      IF( LOPTRA .LE. 20 ) THEN
         GOTO(  10,   1,  30,   1,  50,   1,  70,  80,  90, 100,
     %         110, 120, 130, 140, 150, 160, 170,   1,   1, 200),LOPTRA
      ELSE IF( 27 .LE. LOPTRA .AND. LOPTRA .LE. 40 ) THEN
         CALL LEOPLS( LOPTRA )
         RETURN
      ENDIF
      GOTO 1
C
C     'bascule NOIR et BLANC ou COULEURS'
 10   IF( NBRCOU .EQ. 0 ) THEN
C
C        EN COULEURS
         NBRCOU = NDCOUL - N1COUL
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)')  NBRCOU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE AVEC ' // KERR(MXLGER)(1:4) // ' COULEURS'
         ELSE
            KERR(1) = 'DRAWING WITH ' // KERR(MXLGER)(1:4) // ' COLORS'
         ENDIF
         CALL LERESU
C
      ELSE
C
C        NOIR ET BLANC
         NBRCOU = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE NOIR ET BLANC'
         ELSE
            KERR(1) = 'BLACK and WHITE DRAWING'
         ENDIF
         CALL LERESU
C        ARETES BLANCHES
         NCOUAR = NCBLAN
C        FOND NOIR DES CARACTERES
         NCOUFO = NCNOIR
C        CARACTERES BLANCS
         NCOLIG = NCBLAN
C
      ENDIF
      GOTO 1
C
C     'COULEUR des ARETES'
 30   CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 9999
      ELSE IF( I .EQ. 0 ) THEN
C        LA COULEUR NOIRE
         NCOUAL = 0
      ELSE
         NCOUAL = N1COEL + I
      ENDIF
      GOTO 1
C
C     'POURCENTAGE de REDUCTION des ARETES'
 50   CALL INVITE( 129 )
      CALL LIRRSP( NCVALS , R )
      IF( NCVALS .LE. 0 ) GOTO 9999
      PREDUA = R
      PREDUA = MIN( 100.0 , PREDUA )
      PREDUA = MAX(   0.0 , PREDUA )
      NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:10),'(G10.2)' ) PREDUA
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = '% de REDUCTION des ARETES= ' // KERR(MXLGER)(1:10)
      ELSE
         KERR(1) = 'EDGE REDUCTION %= ' // KERR(MXLGER)(1:10)
      ENDIF
      CALL LERESU
      GOTO 1
C
C     'COULEUR du NOM de la LIGNE'
 70   CALL INVITE( 18 )
      CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 9999
      ELSE IF( I .EQ. 0 ) THEN
C        LA COULEUR NOIRE
         NCOLIG = 0
      ELSE
         NCOLIG = N1COEL + I
      ENDIF
      GOTO 1
C
C     'bascule TRACE ou NON du NOM de la LIGNE'
 80   IF( LPLIGN .EQ. 0 ) THEN
         LPLIGN = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE avec NOM de la LIGNE'
         ELSE
            KERR(1) = 'DRAWING with the LINE NAMES'
         ENDIF
         CALL LERESU
      ELSE
         LPLIGN = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE sans NOM de la LIGNE'
         ELSE
            KERR(1) = 'DRAWING without the LINE NAMES'
         ENDIF
         CALL LERESU
      ENDIF
      GOTO 1
C
C     'bascule TRACE du TITRE ou NON'
 90   IF( IAVTIT .EQ. 0 ) THEN
         IAVTIT = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE avec TITRE'
         ELSE
            KERR(1) = 'DRAWING with TITLE'
         ENDIF
         CALL LERESU
      ELSE
         IAVTIT = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE sans TITRE'
         ELSE
            KERR(1) = 'DRAWING without TITLE'
         ENDIF
         CALL LERESU
      ENDIF
      GOTO 1
C
C     'bascule TRACE du NO des SOMMETS'
 100  IF( IAVNSO .EQ. 0 ) THEN
         IAVNSO = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE avec NO DES SOMMETS'
         ELSE
            KERR(1) = 'DRAWING with VERTEX NUMBERS'
         ENDIF
         CALL LERESU
      ELSE
         IAVNSO = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE sans NO DES SOMMETS'
         ELSE
            KERR(1) = 'DRAWING without VERTEX NUMBERS'
         ENDIF
         CALL LERESU
      ENDIF
      GOTO 1
C
C     'COULEUR du NO des SOMMETS'
 110  CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 9999
      ELSE IF( I .EQ. 0 ) THEN
C        LA COULEUR NOIRE
         NCONSO = 0
      ELSE
         NCONSO = N1COEL + I
      ENDIF
      GOTO 1
C
C     'bascule TRACE du NO des EF (ARETES)'
 120  IF( IAVNEF .EQ. 0 ) THEN
         IAVNEF = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE avec NO des EF'
         ELSE
            KERR(1) = 'DRAWING with EDGE NUMBERS'
         ENDIF
         CALL LERESU
      ELSE
         IAVNEF = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE sans NO des  EF'
         ELSE
            KERR(1) = 'DRAWING without EDGE NUMBERS'
         ENDIF
         CALL LERESU
      ENDIF
      GOTO 1
C
C     'COULEUR du NO des EF'
 130  CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 9999
      ELSE IF( I .EQ. 0 ) THEN
C        LA COULEUR NOIRE
         NCONEF = 0
      ELSE
         NCONEF = N1COEL + I
      ENDIF
      GOTO 1
C
C     'Nbre SUBDIVISIONS une ARETE P3'
 140  CALL INVITE( 34 )
      CALL LIRENT(  NCVALS , I )
      IF( NCVALS .LE. 0 ) GOTO 9999
C     PROTECTION CONTRE LES MAUVAISES ENTREES
      NBSUAR = MAX( 2, MIN( I, MXSUAR ) )
      NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:2),'(I2)' ) NBSUAR
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'Nombre SUBDIVISIONS d''une ARETE P3='
     %        // KERR(MXLGER)(1:2)
      ELSE
         KERR(1) = 'SUBDIVISION NUMBER of P3 EDGE='
     %        // KERR(MXLGER)(1:2)
      ENDIF
      CALL LERESU
      GOTO 1
C
C     'TRACE ou NON TANGENTES des ARETES'
 150  IF( IAVTGA .EQ. 0 ) THEN
         IAVTGA = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE DES TANGENTES DES ARETES'
         ELSE
            KERR(1) = 'DRAWING with EDGE TANGENTS'
         ENDIF
         CALL LERESU
      ELSE
         IAVTGA = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE sans les TANGENTES des ARETES'
         ELSE
            KERR(1) = 'DRAWING without EDGE TANGENTS'
         ENDIF
         CALL LERESU
      ENDIF
      GOTO 1
C
C     'COULEUR du trace des TANGENTES'
 160  CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 9999
      ELSE IF( I .EQ. 0 ) THEN
C        COULEUR NULLE
         NCOTGA = 0
      ELSE
         NCOTGA = N1COEL + I
      ENDIF
      GOTO 1
C
C     'TYPE de TRAIT des TANGENTES'
 170  CALL LIMTCL( 'typtrait' , I )
      IF( I .LT. 0 ) THEN
         GOTO 9999
      ELSE
C        LE TYPE DU TRAIT 0:CONTINU, 1:TIRET'E, 2:DOUBLE TIRETS
         NTRTGA = I
      ENDIF
      GOTO 1
C
C     'EPAISSEUR du trace des ARETES'
 200  CALL INVITE( 77 )
      NCVALS = 4
      I      = NEPARL
      CALL LIRENT( NCVALS , I )
      IF( NCVALS .LE. 0 ) GOTO 9999
      I      = MAX( 0 , I )
      I      = MIN(10 , I )
      NEPARL = I
      GOTO 1
C
C     PLANS SPECIAUX XY YZ XZ ...
      CALL LIMTCL( 'vuesplan' , NMTCL1 )
      IF( NMTCL1 .LE.  0 ) GOTO 10
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
      CALL INVITE( 27 )
      NCVALS = 0
      CALL LIRRSP( NCVALS , AXOLON )
      IF( NCVALS .LE. 0 ) GOTO 1
      AXOLON = MAX( -180.0 , MIN( 180.0 , AXOLON ) )
C
      CALL INVITE( 26 )
      NCVALS = 0
      CALL LIRRSP( NCVALS , AXOLAT )
      IF( NCVALS .LE. 0 ) GOTO 1
C
C     DEFINITION DE PTV ET OEIL A PARTIR DE AXOLON ET AXOLAT
 305  CALL LONLAT( AXOLON, AXOLAT )
C
C     L'AXONOMETRIE
      CALL MATAXO
      CALL ISOFENETRE( -AXOLAR, AXOLAR, -AXOHAU, AXOHAU )
      NOTYVI = 11
      GOTO 1
C
C     XYZ des POINT VU et OEIL
      CALL INVITE( 112 )
      NCVALS = 5
      CALL LIRXYZ( NCVALS , AXOPTV )
      IF( NCVALS .LE. 0 ) GOTO 1
C
      CALL INVITE( 111 )
      NCVALS = 5
      CALL LIRXYZ( NCVALS , AXOEIL )
      IF( NCVALS .LE. 0 ) GOTO 1
C
C     L'AXONOMETRIE
      CALL MATAXO
      CALL ISOFENETRE( -AXOLAR, AXOLAR, -AXOHAU, AXOHAU )
      NOTYVI = 11
      GOTO 1
C
C     AXES DU REPERE DE L'OBJET
 490  CALL LIMTCL( 'tracaxes', NTRAXE )
      IF( NTRAXE .GT. 0 ) THEN
         CALL TRAXES
      ELSE
         NTRAXE = 0
      ENDIF
      GOTO 1
C
C     MODE DE TRACE en ORBITE ZOOM TRANSLATION
 750  LORBITE=1
      NBLGRC(NRERR) = 1
       IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'TRACE avec TRANSLATION ORBITE ZOOM'
      ELSE
         KERR(1) = 'DRAWING with TRANSLATION ORBIT ZOOM'
      ENDIF
      CALL LERESU
      GOTO 1
C
C     PLUS DE ORBITE TRANSLATION ZOOM
 760  LORBITE=0
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'TRACE sans TRANSLATION ORBITE ZOOM'
      ELSE
         KERR(1) = 'DRAWING without TRANSLATION ORBIT ZOOM'
      ENDIF
      CALL LERESU
      GOTO 1
C
C     EFFACER LE TRACE ACTUEL
 5000 CALL EFFACE
C     PLUS D'ITEMS VISIBLES
      CALL ITEMS0
      GOTO 1
C
 9999 LOPTRA = 0
      END
