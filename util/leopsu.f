      SUBROUTINE LEOPSU( LOPTRA )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   LECTURE DES OPTIONS DE TRACE DES FACES d'UNE SURFACE OU VOLUME
C -----
C
C SORTIE :
C --------
C LOPTRA : L'OPTION DE TRACE
C       <=0 : ANNIHILE LA DEMANDE DE TRACE SUITE A UNE ERREUR OU ABANDON
C        >0 : PERMET UN TRACE IMMEDIAT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS        MARS 1991
C ...................................................................012
      include"./incl/langue.inc"
      include"./incl/mxsuaf.inc"
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
      include"./incl/mecoit.inc"
C
C     LE NUMERO D'UNITE DU CLAVIER, FENETRE DES AFFICHAGES ET
C     LE PARAMETRE DE NIVEAU D'INTERACTIVITE
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
C     =========================================
C     SAISIE DE L'OPTION DE TRACE DE LA SURFACE
C     =========================================
 1    CALL LIMTCL( 'opt_surf' , LOPTRA )
C
      IF( LOPTRA .LE. 0  ) RETURN
      IF( LOPTRA .EQ. 49 ) GOTO 490
      IF( LOPTRA .EQ. 50 ) GOTO 5000
      IF( LOPTRA .EQ. 75 ) GOTO 750
      IF( LOPTRA .EQ. 76 ) GOTO 760
      IF( LOPTRA .EQ. 90 ) RETURN
      IF( LOPTRA .LE. 23 ) THEN
         GOTO(  10,  20,  30,  40,  50,  60,  70,  80,  90, 100,
     %         110, 120, 130, 140, 150, 160, 170,   1,   1, 200,
     %         210, 220, 230,   1,   1) , LOPTRA
      ELSE IF( 26 .LE. LOPTRA .AND. LOPTRA .LE. 40 ) THEN
         CALL LEOPLS( LOPTRA )
         RETURN
      ENDIF
      GOTO 1
C
C     'bascule NOIR et BLANC ou COULEURS'
 10   IF( NBRCOU .EQ. 0 ) THEN
         NBRCOU = NDCOUL - N1COUL
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)')  NBRCOU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE AVEC ' // KERR(MXLGER)(1:4) // ' COULEURS'
         ELSE
            KERR(1) = 'DRAWING WITH ' // KERR(MXLGER)(1:4) // ' COLORS'
         ENDIF
         CALL LERESU
      ELSE
         NBRCOU = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE NOIR ET BLANC'
         ELSE
            KERR(1) = 'BLACK and WHITE DRAWING'
         ENDIF
         CALL LERESU
C        FACES NOIRES
         NCOUFA = NCNOIR
C        ARETES BLANCHES
         NCOUAF = NCBLAN
C        CARACTERES BLANCS
         NCOSUR = NCBLAN
C        FOND NOIR DES CARACTERES
         NCOUFO = NCNOIR
      ENDIF
      GOTO 1
C
C     'bascule TRACE ou NON des ARETES des FACES'
 20   IF( IAVARE .EQ. 0 ) THEN
         IAVARE = 1
         NCOUAF = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE des ARETES'
         ELSE
            KERR(1) = 'DRAWING of EDGES'
         ENDIF
         CALL LERESU
      ELSE
         IAVARE = 0
C        COULEUR INVISIBLE
         NCOUAF = -2
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'PAS DE TRACE DES ARETES'
         ELSE
            KERR(1) = 'NO EDGE DRAWING'
         ENDIF
         CALL LERESU
      ENDIF
      GOTO 1
C
C     'COULEUR des ARETES des FACES'
 30   CALL LIMTCL( 'couleur0' , I )
      IF( I .EQ. -2 ) THEN
         NCOUAF = -2
      ELSEIF( I .LT.  0 ) THEN
         GOTO 9999
      ELSE IF( I .EQ. 0 ) THEN
C        COULEUR NOIRE
         NCOUAF = 0
      ELSE
         NCOUAF = N1COEL + I
      ENDIF
      GOTO 1
C
C     'bascule TRACE ou NON des FACES'
 40   IF( IAVFAC .EQ. 0 ) THEN
         IAVFAC = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE DES FACES'
         ELSE
            KERR(1) = 'DRAWING of FACES'
         ENDIF
         CALL LERESU
      ELSE
         IAVFAC = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'PAS DE TRACE DES FACES'
         ELSE
            KERR(1) = 'NO FACE DRAWING'
         ENDIF
         CALL LERESU
      ENDIF
      GOTO 1
C
C     'POURCENTAGE de REDUCTION des FACES'
 50   CALL INVITE( 130 )
      CALL LIRRSP( NCVALS , R )
      IF( NCVALS .LE. 0 ) GOTO 9999
      PREDUF = R
      PREDUF = MIN( 100.0 , PREDUF )
      PREDUF = MAX(   0.0 , PREDUF )
      NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:10),'(G10.2)' ) PREDUF
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = '% de REDUCTION des FACES= ' // KERR(MXLGER)(1:10)
      ELSE
         KERR(1) = 'FACE REDUCTION %= ' // KERR(MXLGER)(1:10)
      ENDIF
      CALL LERESU
      GOTO 1
C
C     'bascule ELOIGNEMENT ou NON'
 60   IF( IAVELO .EQ. 0 ) THEN
         IAVELO = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE avec ELOIGNEMENT'
         ELSE
            KERR(1) = 'DRAWING with DARKNESS'
        ENDIF
         CALL LERESU
      ELSE
         IAVELO = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE sans ELOIGNEMENT'
         ELSE
            KERR(1) = 'DRAWING without DARKNESS'
        ENDIF
         CALL LERESU
      ENDIF
      GOTO 1
C
C     'COULEUR du NOM des SURFACES'
 70   CALL INVITE( 19 )
      CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 9999
      ELSE IF( I .EQ. 0 ) THEN
C        COULEUR NOIRE
         NCOSUR = 0
      ELSE
         NCOSUR = N1COEL + I
      ENDIF
      GOTO 1
C
C     'bascule TRACE ou NON du NOM de la SURFACE'
 80   IF( LPSURF .EQ. 0 ) THEN
         LPSURF = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE avec NOM de la SURFACE'
         ELSE
            KERR(1) = 'DRAWING with the SURFACE NAMES'
         ENDIF
         CALL LERESU
      ELSE
         LPSURF = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE sans NOM de la SURFACE'
         ELSE
            KERR(1) = 'DRAWING without the SURFACE NAMES'
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
C        COULEUR NOIRE
         NCONSO = 0
      ELSE
         NCONSO = N1COEL + I
      ENDIF
      GOTO 1
C
C     'bascule TRACE du NO des EF'
 120  IF( IAVNEF .EQ. 0 ) THEN
         IAVNEF = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE avec NO des EF'
         ELSE
            KERR(1) = 'DRAWING with FE NUMBERS'
         ENDIF
         CALL LERESU
      ELSE
         IAVNEF = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE sans NO des EF'
         ELSE
            KERR(1) = 'DRAWING without FE NUMBERS'
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
C        COULEUR NOIRE
         NCONEF = 0
      ELSE
         NCONEF = N1COEL + I
      ENDIF
      GOTO 1
C
C     'Max SUBDIVISIONS 1 ARETE des FACES P3'
 140  CALL INVITE( 29 )
      CALL LIRENT(  NCVALS , I )
      IF( NCVALS .LE. 0 ) GOTO 9999
C     PROTECTION CONTRE LES MAUVAISES ENTREES
      NBSUAF = MAX( 2, MIN( I, MXSUAF ) )
      NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:2),'(I2)' ) NBSUAF
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'Nombre SUBDIVISIONS d''une ARETE des FACES P3='
     %        // KERR(MXLGER)(1:2)
      ELSE
         KERR(1) = 'SUBDIVISION NUMBER of EDGE of P3 FACE='
     %        // KERR(MXLGER)(1:2)
      ENDIF
      CALL LERESU
      GOTO 1
C
C     'TRACE ou NON des TANGENTES des FACES'
 150  IF( IAVTGF .EQ. 0 ) THEN
         IAVTGF = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE DES TANGENTES DES FACES'
         ELSE
            KERR(1) = 'DRAWING with FACE TANGENTS'
         ENDIF
         CALL LERESU
      ELSE
         IAVTGF = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE sans les TANGENTES des FACES'
         ELSE
            KERR(1) = 'DRAWING without FACE TANGENTS'
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
         NCOTGF = 0
      ELSE
         NCOTGF = N1COEL + I
      ENDIF
      GOTO 1
C
C     'TYPE de TRAIT des TANGENTES'
 170  CALL LIMTCL( 'typtrait' , I )
      IF( I .LT. 0 ) THEN
         GOTO 9999
      ELSE
C        LE TYPE DU TRAIT 0:CONTINU, 1:TIRET'E, 2:DOUBLE TIRETS
         NTRTGF = I
      ENDIF
      GOTO 1
C
C     ligne EPAISSIE des ARETES DES FACES
 200  CALL INVITE( 77 )
      NCVALS = 4
      I      = NEPARF
      CALL LIRENT( NCVALS , I )
      IF( NCVALS .LE. 0 ) GOTO 9999
      I      = MAX( 0 , I )
      I      = MIN(10 , I )
      NEPARF = I
      GOTO 1
C
C     'TRACE ou NON du VECTEUR NORMAL aux FACES'
 210  IF( IAVNRF .EQ. 0 ) THEN
         IAVNRF = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE du VECTEUR NORMAL AUX FACES'
         ELSE
            KERR(1) = 'DRAWING with FACE NORMAL VECTOR'
         ENDIF
         CALL LERESU
      ELSE
         IAVNRF = 0
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE sans le VECTEUR NORMAL AUX FACES'
         ELSE
            KERR(1) = 'DRAWING without FACE NORMAL VECTOR'
         ENDIF
         CALL LERESU
      ENDIF
      GOTO 1
C
C     'COULEUR du trace du VECTEUR NORMAL aux FACES'
 220  CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 9999
      ELSE IF( I .EQ. 0 ) THEN
C        COULEUR NULLE
         NCONRF = 0
      ELSE
         NCONRF = N1COEL + I
      ENDIF
      GOTO 1
C
C     'TYPE de TRAIT du VECTEUR NORMAL aux FACES'
 230  CALL LIMTCL( 'typtrait' , I )
      IF( I .LT. 0 ) THEN
         GOTO 9999
      ELSE
C        LE TYPE DU TRAIT 0:CONTINU, 1:TIRET'E, 2:DOUBLE TIRETS
         NTRNRF = I
      ENDIF
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
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'TRACE sans TRANSLATION ORBIT ZOOM'
      ELSE
         KERR(1) = 'DRAWING without TRANSLATION ORBIT ZOOM'
      ENDIF
      NBLGRC(NRERR) = 1
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
