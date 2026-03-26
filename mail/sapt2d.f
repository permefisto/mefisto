      SUBROUTINE SAPT2D
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SAISIR A LA SOURIS DES POINTS 2D
C ----- SIMULER SUR LE FICHIER FRAPPE LA FRAPPE DE CES POINTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1990
C2345X7..............................................................012
      IMPLICIT           INTEGER(W)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pilect.inc"
      include"./incl/xyzext.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a_point__definition.inc"
      COMMON / EPSSSS / EPZERO,EPSXYZ
      COMMON / TRACE1 / PTCOUR (3) , MUETR1 (7)
      REAL              RMCN(1)
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      EQUIVALENCE       (MCN(1),RMCN(1))
      COMMON / UNITES /  LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NUNITE(27)
      CHARACTER*24       KNOMPO,KNOMPT,KNOMTR,KNOM,KNOMDP
      CHARACTER*25       KSYMBO
      CHARACTER*4        KENTIE
      CHARACTER*24       KAIMAX,KAIMAY
      REAL               CADRSA(2,3),CADRS0(2,3),P1(3),XYZ(3),XY(2)
      INTEGER            MNSO(3)
      DOUBLE PRECISION   D2D3(3,3)
      LOGICAL            AIMANX,AIMANY
      SAVE               NBPOIN,KNOMPO
      DATA               NBPOIN / 0 /
      DATA               KNOMPO / 'P' /
C
C     VALEURS PAR DEFAUT
C     ==================
C     CADRE DE SAISIE CADRSA(X Y , MIN MAX PAS)
      IF( COOEXT(1,1) .LT. RINFO('GRAND') )  THEN
         ECAMAX = MAX(COOEXT(1,2)-COOEXT(1,1) , COOEXT(2,2)-COOEXT(2,1))
         IF( ECAMAX .LE. 0. ) ECAMAX = 1.
         ECAMAX = ECAMAX * 0.1
C        ABSCISSES
         CADRS0(1,1) = COOEXT(1,1) - ECAMAX
         CADRS0(1,2) = COOEXT(1,2) + ECAMAX
C        ORDONNEES
         CADRS0(2,1) = COOEXT(2,1) - ECAMAX
         CADRS0(2,2) = COOEXT(2,2) + ECAMAX
C        PAS
         NBPASX = 10
         CADRS0(1,3) = ABS( ( CADRS0(1,2) - CADRS0(1,1) ) / NBPASX )
         NBPASY = 10
         CADRS0(2,3) = ABS( ( CADRS0(2,2) - CADRS0(2,1) ) / NBPASY )
      ELSE
C        ABSCISSES
         CADRS0(1,1) =-1.0
         CADRS0(1,2) = 11.0
C        ORDONNEES
         CADRS0(2,1) =-1.0
         CADRS0(2,2) = 11.0
C        PAS
         CADRS0(1,3) = 1.0
         CADRS0(2,3) = 1.0
      ENDIF
C
C     ATTRACTION PAR DEFAUT
      IF( LANGAG .EQ. 0 ) THEN
         KAIMAX = 'NE PAS ATTRACTER EN X'
         KAIMAY = 'NE PAS ATTRACTER EN Y'
         NAIMAX = 21
         NAIMAY = 21
      ELSE
         KAIMAX = 'NO X-ATTRACTION'
         KAIMAY = 'NO Y-ATTRACTION'
         NAIMAX = 15
         NAIMAY = 15
      ENDIF
      AIMANX = .TRUE.
      AIMANY = .TRUE.
C
C     LE NOM GENERIQUE DE LA TRANSFORMATION
      KNOMTR = 'I'
      NUTRAN = 1
C     PAS DE DEFINITION DE PLAN
      LEPLAN = 0
      XCTE   = 0
      YCTE   = 0
      ZCTE   = 0
      LDEPLA = 0
C
C     LECTURE DES DONNEES
 10   CALL LIMTCL( 'saisipt2' , NUMTCL )
      IF( NUMTCL .EQ. -1 ) GOTO 9900
      IF( NUMTCL .EQ. 20 ) GOTO 500
      GOTO( 100 , 200 , 300 , 400 , 450 , 470 , 475 , 480 , 490 ),NUMTCL
C
C     LA DEMANDE DE NOM GENERIQUE DES POINTS
C     ======================================
 100  CALL INVITE( 62 )
      NCVALS = 2
      KNOMPO = 'P                       '
      CALL LIRCAR( NCVALS , KNOMPO )
      NBPOIN = 0
      GOTO 10
C
C     LA DEMANDE DE NOM DE LA TRANSFORMATION
C     ======================================
 200  CALL INVITE( 44 )
      NCVALS = 2
      CALL LIRCAR( NCVALS , KNOMTR )
      GOTO 10
C
C     ABSCISSES MINIMALE MAXIMALE et PAS du CADRE des POINTS a saisir
C     ================================================================
 300  CALL INVITE( 110 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , CADRS0(1,1) )
      IF( NCVALS .LE. 0 ) GOTO 10
      CALL INVITE( 109 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , CADRS0(1,2) )
      IF( NCVALS .LE. 0 ) GOTO 10
      IF( CADRS0(1,2) .LT. CADRS0(1,1) ) THEN
C        PERMUTATION
         X           = CADRS0(1,1)
         CADRS0(1,1) = CADRS0(1,2)
         CADRS0(1,2) = X
      ENDIF
      CALL INVITE( 87 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , CADRS0(1,3) )
      IF( CADRS0(1,3) .GT. CADRS0(1,2)-CADRS0(1,1) ) THEN
         CADRS0(1,3) = ABS( (CADRS0(1,2)-CADRS0(1,1))/10 )
      ENDIF
      GOTO 10
C
C     ORDONNEES MINIMALE MAXIMALE et PAS du CADRE des POINTS a saisir
C     ================================================================
 400  CALL INVITE( 119 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , CADRS0(2,1) )
      IF( NCVALS .LE. 0 ) GOTO 10
      CALL INVITE( 118 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , CADRS0(2,2) )
      IF( NCVALS .LE. 0 ) GOTO 10
      IF( CADRS0(2,2) .LT. CADRS0(2,1) ) THEN
C        PERMUTATION
         X           = CADRS0(2,1)
         CADRS0(2,1) = CADRS0(2,2)
         CADRS0(2,2) = X
      ENDIF
      CALL INVITE( 88 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , CADRS0(2,3) )
      IF( CADRS0(2,3) .GT. CADRS0(2,2)-CADRS0(2,1) ) THEN
         CADRS0(2,3) = ABS( (CADRS0(2,2)-CADRS0(2,1))/10 )
      ENDIF
      GOTO 10
C
C     NOMS DES 3 POINTS DE DEFINITION DES POINTS A SAISIR
C     ===================================================
 450  DO 460 I=1,3
         MNSO(I) = 0
         GOTO(451,452,453),I
 451     CALL INVITE( 52 )
         GOTO 455
 452     CALL INVITE( 53 )
         GOTO 455
 453     CALL INVITE( 54 )
 455     NCVALS = 0
         CALL LIRCAR( NCVALS , KNOM )
         CALL LXLXOU( NTPOIN , KNOM , NTLXPO , MNLXPO )
         IF( NTLXPO .LE. 0 ) THEN
            NBLGRC(1) = 1
            KERR(1)   = 'POINT INCONNU : '//KNOM
            CALL LEREUR
            GOTO 10
         ENDIF
         CALL LXTSOU( NTLXPO , 'XYZSOMMET' , NTSO , MNSO(I) )
         IF( NTSO .LE. 0 ) THEN
            NBLGRC(1) = 2
            KERR(1) = 'POINT: '//KNOM
            KERR(2) = 'SANS SOMMETS'
            CALL LEREUR
            GOTO 10
         ENDIF
 460  CONTINUE
C
C     GENERATION DE LA MATRICE ET LE SECOND MEMBRE DU CHANGEMENT
C     DE REPERE
      CALL DF3D2D( RMCN(MNSO(1)+WYZSOM) ,
     %             RMCN(MNSO(2)+WYZSOM) ,
     %             RMCN(MNSO(3)+WYZSOM) , D2D3 , IERR )
      IF( IERR .GT. 0 ) GOTO 10
C
C     PROTECTION DU POINT 1 EN CAS DE MODIFICATION
      P1(1)  = RMCN(MNSO(1)+WYZSOM)
      P1(2)  = RMCN(MNSO(1)+WYZSOM+1)
      P1(3)  = RMCN(MNSO(1)+WYZSOM+2)
      LEPLAN = 4
      XCTE   = 0
      YCTE   = 0
      ZCTE   = 0
      GOTO 10
C
C     PLAN A X=CTE
C     ============
 470  CALL INVITE( 104 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , XCTE )
      YCTE   = 0
      ZCTE   = 0
      LEPLAN = 1
      GOTO 10
C
C     PLAN A Y=CTE
C     ============
 475  CALL INVITE( 113 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , YCTE )
      XCTE   = 0
      ZCTE   = 0
      LEPLAN = 2
      GOTO 10
C
C     PLAN A Z=CTE
C     ============
 480  CALL INVITE( 120 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , ZCTE )
      XCTE   = 0
      YCTE   = 0
      LEPLAN = 3
      GOTO 10
C
C     PAS DE PLAN DE DEFINITION DES POINTS A SAISIR
C     =============================================
 490  LEPLAN = 0
      XCTE   = 0
      YCTE   = 0
      ZCTE   = 0
      GOTO 10
C
C     =============================
C     SAISIE a la SOURIS des POINTS
C     =============================
C     LE NUMERO DE LA TRANSFORMATION
 500  CALL NUOBNM( 'TRANSFO' , KNOMTR , NUTRAN )
      IF( NUTRAN .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRANSFORMATION INCONNUE'
         ELSE
            KERR(1) = 'UNKNOWN MAPPING'
         ENDIF
         CALL LEREUR
         GOTO 200
      ENDIF
C
C     RECHERCHE DE LA DIMENSION MAXIMALE ECAMAX DE L'OBJET
      IF( CADRS0(1,1) .GT. CADRS0(1,2) ) THEN
         X           = CADRS0(1,1)
         CADRS0(1,1) = CADRS0(1,2)
         CADRS0(1,2) = X
      ENDIF
      ECAMAX = CADRS0(1,2) - CADRS0(1,1)
      IF( ECAMAX .LE. 0. .OR. ECAMAX .GT. 1E30 ) ECAMAX = 10.
C
      IF( CADRS0(2,1) .GT. CADRS0(2,2) ) THEN
         X           = CADRS0(2,1)
         CADRS0(2,1) = CADRS0(2,2)
         CADRS0(2,2) = X
      ENDIF
      ECAMAY = CADRS0(2,2) - CADRS0(2,1)
      IF( ECAMAY .LE. 0. .OR. ECAMAY .GT. 1E30 ) ECAMAY = 10.
C
      RAPECR = (LHPXFE / CYMMPX) / (LAPXFE / CXMMPX)
      IF( ECAMAY .LE. ECAMAX*RAPECR )THEN
         ECAMAY = ECAMAX * RAPECR
      ELSE
         ECAMAX = ECAMAY / RAPECR
      ENDIF
C
C     LE MILIEU
      X = ( CADRS0(1,1) + CADRS0(1,2) ) * 0.5
      CADRSA(1,1) = X - ECAMAX * 0.5
      CADRSA(1,2) = X + ECAMAX * 0.5
      CADRSA(1,3) = CADRS0(1,3)
C
      Y = ( CADRS0(2,1) + CADRS0(2,2) ) * 0.5
      CADRSA(2,1) = Y - ECAMAY * 0.5
      CADRSA(2,2) = Y + ECAMAY * 0.5
      CADRSA(2,3) = CADRS0(2,3)
C
C     LE CADRSA OBJET EN UNITES UTILISATEUR
      CALL ISOFENETRE( CADRSA(1,1), CADRSA(1,2),
     %                 CADRSA(2,1), CADRSA(2,2) )
C
C     MISE A JOUR
      CADRSA(1,1) = XOBMIN
      CADRSA(1,2) = XOBMAX
      CADRSA(2,1) = YOBMIN
      CADRSA(2,2) = YOBMAX
C
      AXOPTV(1) = ( XOBMIN + XOBMAX ) * 0.5
      AXOPTV(2) = ( YOBMIN + YOBMAX ) * 0.5
      AXOLAR = ( XOBMAX - XOBMIN ) * 0.5
      AXOHAU = ( YOBMAX - YOBMIN ) * 0.5
C
C     ==================
C     BOUCLE DES SAISIES
C     ==================
C     L'ECRAN EST EFFACE
 600  CALL EFFACE
C
C     REMISE A ZERO DES ITEMS
      CALL ITEMS0
C
C     LE TRACE EVENTUEL DES VERTICALES
C     --------------------------------
      CALL XVEPAISSEUR( 1 )
      ECAR   = ECAMAX * 0.02
      X      = ( CADRSA(1,2) - CADRSA(1,1) ) / CADRS0(1,3)
      X      = MAX( 2.0 , X )
      NBPASX = NINT( X )
      NX     = NBPASX / 2
      X      = CADRS0(1,1) - NX * CADRS0(1,3)
      XC     = ECAR / 3
      DO 610 I = -NX, 3*NX
         IF( X .LT. CADRSA(1,1) ) GOTO 607
         CALL TRAIT2D( NCGRIS, X, CADRSA(2,1), X, CADRSA(2,2) )
         WRITE(KNOM(1:13),'(G11.3)') X
         CALL TEXTE2D( NCROUG, X-ECAR, CADRSA(2,2)-ECAR+XC, KNOM(1:13) )
         CALL TEXTE2D( NCROUG, X-ECAR, CADRSA(2,1)+ECAR-XC, KNOM(1:13) )
 607     X  = X + CADRS0(1,3)
         XC = -XC
         IF( X .GT. CADRSA(1,2) ) GOTO 615
 610  CONTINUE
C
C     LE TRACE EVENTUEL DES HORIZONTALES
C     ----------------------------------
 615  X      =( CADRSA(2,2) - CADRSA(2,1) ) / CADRS0(2,3)
      X      = MAX( 2.0 , X )
      NBPASY = NINT( X )
      NY     = NBPASY / 2
      X      = CADRS0(2,1) - NY * CADRS0(2,3)
      DO 618 I = -NY , 3*NY
         IF( X .LT. CADRSA(2,1) ) GOTO 617
         CALL TRAIT2D( NCGRIS, CADRSA(1,1), X, CADRSA(1,2) ,X )
         WRITE(KNOM(1:13),'(G11.3)') X
         CALL TEXTE2D( NCVERT, CADRSA(1,2)-ECAR*3.0 , X+ECAR*0.3,
     %                 KNOM(1:13) )
         CALL TEXTE2D( NCVERT, CADRSA(1,1)+ECAR, X+ECAR*0.3,
     %                 KNOM(1:13) )
 617     X = X + CADRS0(2,3)
         IF( X .GT. CADRSA(2,2) ) GOTO 620
 618  CONTINUE
C
C     TRACE DU RECTANGLE DU MENU A DROITE ET EN HAUT DE L'ECRAN
C     =========================================================
 620  NBMXCA = 21
      NBLIGM = 5
      LIGLPX = NPLACA * NBMXCA + 4*ECARLR(NRMENU)
      LIGHPX = NPHACA + 2*ECARLR(NRMENU)
C
C     LE COIN SUPERIEUR DROIT DU MENU DANS LA FENETRE
      MENUX2 = LAPXFE - 60
      MENUY2 = 30
C     LE COIN INFERIEUR GAUCHE DU MENU DANS LA FENETRE
      MENUX1 = MENUX2 - LIGLPX
      MENUY1 = MENUY2 + LIGHPX * NBLIGM
C
C     LE TRACE DU RECTANGLE DU MENU
      CALL XVCOULEUR( NCBLEU )
      CALL XVRECTANGLE( MENUX1, MENUY2, MENUX2-MENUX1, MENUY1-MENUY2 )
      CALL XVEPAISSEUR( 2 )
      CALL XVCOULEUR( NCCYAN )
      CALL XVBORDRECTANGLE( MENUX1, MENUY2,
     %                       MENUX2-MENUX1, MENUY1-MENUY2 )
C
C     LES LIGNES INTERMEDIAIRES
      CALL XVCOULEUR( NCCYAN )
      NY = MENUY1
      DO 625 I=1,NBLIGM-1
         NY = NY - LIGHPX
         CALL XVTRAIT( MENUX1, NY, MENUX2 , NY )
 625  CONTINUE
C
C     LE TEXTE DU MENU
      NX = MENUX1 + 2*ECARLR(NRMENU)
      NY = MENUY1 - ECARLR(NRMENU)
      CALL XVCOULEUR( NCBLAN )
      IF( LANGAG .EQ. 0 ) THEN
         CALL XVTEXTE( 'ABANDONNER', 10, NX, NY )
      ELSE
         CALL XVTEXTE( 'ESCAPE', 6, NX, NY )
      ENDIF
C
      NY = NY - LIGHPX
      CALL XVTEXTE( KAIMAY, NAIMAY, NX, NY )
C
      NY = NY - LIGHPX
      CALL XVTEXTE( KAIMAX, NAIMAX, NX, NY )
C
      NY = NY - LIGHPX
      IF( LANGAG .EQ. 0 ) THEN
         CALL XVTEXTE( 'TUER un POINT', 13, NX, NY )
      ELSE
         CALL XVTEXTE( 'DELETE a POINT', 14, NX, NY )
      ENDIF
C
      NY = NY - LIGHPX
      IF( LANGAG .EQ. 0 ) THEN
         CALL XVTEXTE( 'DEPLACER un POINT', 17, NX, NY )
      ELSE
         CALL XVTEXTE( 'MOVE a POINT', 12, NX, NY )
      ENDIF
C
C     LE TRACE DE X ET Y
C     ------------------
      CALL XVTEXTE('--X-->', 6, LAPXFE-9*NPLACA, LHPXFE-ECARLR(NRMENU))
      CALL XVTEXTE( 'Y', 1, 10*NPLACA, LIGHPX )
C
C     TRACE des POINTS EXISTANTS
C     ==========================
C     LE DEBUT DU CHAINAGE DES POINTS OCCUPES DANS LE LEXIQUE
      CALL TAMSOU( NTPOIN , MNPOIN )
      NUPOIN = MCN( MNPOIN + 5 )
C     NOMBRE D'ENTIERS POUR STOCKER UN NOM DE POINT
      NBENNM = MCN( MNPOIN + 2 )
C
C     LA BOUCLE SUR LES POINTS OCCUPES
 642  IF( NUPOIN .GT. 0 ) THEN
C        ADRESSE MCN DU DEBUT DU POINT DANS LE LEXIQUE
         MNOBJ = MNPOIN + MCN(MNPOIN) * NUPOIN
C        LE LEXIQUE DE CE POINT EXISTE-T-IL ?
         NTOBJ = MCN( MNOBJ + NBENNM + 2 )
         IF( NTOBJ .GT. 0 ) THEN
C           CE POINT EXISTE : RECHERCHE DE SON NOM
            CALL ENTNOM( NBENNM , MCN( MNOBJ ) , KNOMPT )
C           TRACE DE CE POINT DE NOM KNOMPT
            NBC = INDEX( KNOMPT , ' ' ) - 1
            IF( NBC .LE. 0 ) NBC = LEN( KNOMPT )
            CALL LXTSOU( NTOBJ , 'XYZSOMMET' , NTSOM , MNSOM )
            IF( NTSOM .GT. 0 ) THEN
               IF( LEPLAN .EQ. 0 ) THEN
                  XYZ(1) = RMCN(MNSOM+WYZSOM)
                  XYZ(2) = RMCN(MNSOM+WYZSOM+1)
                  XYZ(3) = RMCN(MNSOM+WYZSOM+2)
               ELSE IF( LEPLAN .EQ. 1 ) THEN
                  XYZ(1) = RMCN(MNSOM+WYZSOM+1)
                  XYZ(2) = RMCN(MNSOM+WYZSOM+2)
                  XYZ(3) = RMCN(MNSOM+WYZSOM)-XCTE
               ELSE IF( LEPLAN .EQ. 2 ) THEN
                  XYZ(1) = RMCN(MNSOM+WYZSOM)
                  XYZ(2) = RMCN(MNSOM+WYZSOM+2)
                  XYZ(3) = RMCN(MNSOM+WYZSOM+1)-YCTE
               ELSE IF( LEPLAN .EQ. 3 ) THEN
                  XYZ(1) = RMCN(MNSOM+WYZSOM)
                  XYZ(2) = RMCN(MNSOM+WYZSOM+1)
                  XYZ(3) = RMCN(MNSOM+WYZSOM+2)-ZCTE
               ELSE IF ( LEPLAN .EQ. 4 ) THEN
C                 PASSAGE DANS LE SYSTEME D'AXES AVEC LE PLAN
                  CALL CH3D3D( P1,D2D3,
     %                         RMCN(MNSOM+WYZSOM),XYZ )
               ENDIF
C              TRACE DU POINT S'IL EST DANS LE PLAN
               IF( ABS(XYZ(3)) .LE. EPSXYZ ) THEN
                  CALL ITEMP2( XYZ,KNOMPT(1:NBC),NUPOIN )
               ENDIF
            ENDIF
         ENDIF
C        PASSAGE AU POINT SUIVANT
         NUPOIN = MCN( MNOBJ + NBENNM )
         GOTO 642
      ENDIF
C
C     LE NOMBRE DE CARACTERES NON BLANCS DU NOM DU POINT
      NBCAPO = NUDCNB( KNOMPO )
      NBCAPO = MIN( NBCAPO , 20 )
      KNOMPT = KNOMPO
C
      CALL MEMPXFENETRE
C
C     =======================================
C     SAISIE D'UN POINT PAR CLIC DE LA SOURIS
C     =======================================
 630  CALL SAIPTC( NOTYEV, NX, NY, NOCHAR )
      IF( NOTYEV .EQ. 0 ) GOTO 9900
      IF( NOTYEV .LT. 0 ) GOTO 630
C
C     LE POINT SAISI EST IL DANS LEMENU  ?
      IF( NX .GE. MENUX1 .AND. NX .LE. MENUX2 .AND.
     %    NY .LE. MENUY1 .AND. NY .GE. MENUY2 ) THEN
C
C        OUI  : DANS QUELLE LIGNE ?
         IF( NY .GE. MENUY1 - LIGHPX ) THEN
C           ABANDON DE LA SAISIE
            GOTO 10
C
         ELSE IF( NY .GE. MENUY1-2*LIGHPX ) THEN
C           ATTRACTER OU NON EN Y
            AIMANY = .NOT. AIMANY
            IF( AIMANY ) THEN
               IF( LANGAG .EQ. 0 ) THEN
                  KAIMAY = 'NE PAS ATTRACTER EN Y'
                  NAIMAY = 21
               ELSE
                  KAIMAY = 'NO Y-ATTRACTION'
                  NAIMAY = 15
               ENDIF
            ELSE
               IF( LANGAG .EQ. 0 ) THEN
                  KAIMAY = 'ATTRACTER EN Y'
                  NAIMAY = 14
               ELSE
                  KAIMAY = 'Y-ATTRACTION'
                  NAIMAY = 12
               ENDIF
            ENDIF
            GOTO 600
C
         ELSE IF( NY .GE. MENUY1-3*LIGHPX )THEN
C           ATTRACTER OU NON EN X
            AIMANX = .NOT. AIMANX
            IF( AIMANX ) THEN
               IF( LANGAG .EQ. 0 ) THEN
                  KAIMAX = 'NE PAS ATTRACTER EN X'
                  NAIMAX = 21
               ELSE
                  KAIMAX = 'NO X-ATTRACTION'
                  NAIMAX = 15
               ENDIF
            ELSE
               IF( LANGAG .EQ. 0 ) THEN
                  KAIMAX = 'ATTRACTER EN X'
                  NAIMAX = 14
               ELSE
                  KAIMAX = 'X-ATTRACTION'
                  NAIMAX = 12
               ENDIF
            ENDIF
            GOTO 600
C
         ELSE
C
C           TUER OU DEPLACER UN POINT A CLIQUER
            CALL PTPLPR( 1 , NUPOIN )
            IF( NUPOIN .LE. 0 ) GOTO 600
C           DESTRUCTION DU LEXIQUE DU POINT
            CALL NMOBNU( 'POINT' , NUPOIN , KNOMPT )
            CALL LXLXDS(  NTPOIN , KNOMPT )
C           SAUVEGARDE DE L'OPERATION DESTRUCTION SUR LE FICHIER FRAPPE
            KERR(MXLGER) = KNOMPT//' ; I; 82;'
            CALL SANSDBL( KERR(MXLGER), I )
            WRITE(NFFRAP,*) KERR(MXLGER)(1:I)
            IF( Y .LE. MENUY1-5*LIGHPX ) THEN
C              DEPLACER
               LDEPLA = 1
               KNOMDP = KNOMPT
            ELSE
C              TUER
               LDEPLA = 0
            ENDIF
            GOTO 600
         ENDIF
      ENDIF
C
C     PRISE EN COMPTE OU NON DE L'ATTRACTION SUR LA GRILLE
C     ====================================================
C     COORDONNEES OBJET 2D
      X = XOB2PX( NX )
      Y = YOB2PX( NY )
C
      IF( AIMANX ) THEN
C        OUI: ATTRACTION SUR LA GRILLE EN X
         I = NINT( ( X - CADRS0(1,1) ) / CADRS0(1,3) )
         X = CADRS0(1,1) + I * CADRS0(1,3)
      ENDIF
C
      IF( AIMANY ) THEN
C        OUI: ATTRACTION SUR LA GRILLE EN Y
         I = NINT( ( Y - CADRS0(2,1) ) / CADRS0(2,3) )
         Y = CADRS0(2,1) + I * CADRS0(2,3)
      ENDIF
C
C     RECHERCHE DU NUMERO POUR LE NOM DU POINT
C     ========================================
      IF( LDEPLA .NE. 0 ) THEN
C        DEPLACEMENT DONC NOM CONNU AVEC NUMERO NUPOIN
         KNOMPT = KNOMDP
         GOTO 670
      ENDIF
C     RECHERCHE DU PLUS PETIT ENTIER TEL QUE LE NOM GENERIQUE CONCATENE
C     AVEC CET ENTIER NE SOIT PAS DEJA LE NOM D'UN POINT
      NBPOIN = 0
 650  NBPOIN = NBPOIN + 1
      WRITE( KENTIE, '(I4)' ) NBPOIN
      DO 655 J=4,1,-1
         IF( KENTIE(J:J) .EQ. ' ' ) GOTO 660
 655  CONTINUE
      J = 0
 660  KNOMPT(NBCAPO+1:NBCAPO+4-J) = KENTIE(J+1:4)
C
 670  CALL LXLXOU( NTPOIN , KNOMPT , NTLXPO , MNLXPO )
      IF( NTLXPO .GT. 0 ) GOTO 650
C     ICI LE POINT NBPOIN EST LIBRE
C
C     DECLARATION DU POINT
C     ====================
      CALL LXLXDC( NTPOIN , KNOMPT , 24 , 8 )
      CALL LXLXOU( NTPOIN , KNOMPT , NTLXPO , MNLXPO )
      IF( NTLXPO .LE. 0 ) GOTO 9000
C     LE NUMERO DU POINT DANS LE LEXIQUE DES POINTS
      CALL NUOBNM( 'POINT' , KNOMPT , NUPOIN )
C
C     CONSTRUCTION DU TABLEAU 'DEFINITION'
C     ====================================
      CALL LXTNDC( NTLXPO , 'DEFINITION' , 'MOTS'  , WOORPO + 3 )
      CALL LXTSOU( NTLXPO , 'DEFINITION' ,  NTDFPO , MNDFPO )
      IF( NTDFPO .LE. 0 ) GOTO 9000
C
C     LA TRANSFORMATION
      MCN( MNDFPO + WTYTRP ) = NUTRAN
C
C     LE TYPE DU POINT
      MCN( MNDFPO + WUTYPO ) = 1
C
C     X Y Z
C     PRISE EN COMPTE DU FACTEUR MULTIPLICATIF
      XE = X
      YE = Y
      IF( LEPLAN .EQ. 1 ) THEN
C        XCTE
         RMCN( MNDFPO + WOORPO     ) = XCTE
         RMCN( MNDFPO + WOORPO + 1 ) = XE
         RMCN( MNDFPO + WOORPO + 2 ) = YE
      ELSE IF( LEPLAN .EQ. 2 ) THEN
C        YCTE
         RMCN( MNDFPO + WOORPO     ) = XE
         RMCN( MNDFPO + WOORPO + 1 ) = YCTE
         RMCN( MNDFPO + WOORPO + 2 ) = YE
      ELSE IF( LEPLAN .EQ. 3 ) THEN
C        ZCTE
         RMCN( MNDFPO + WOORPO     ) = XE
         RMCN( MNDFPO + WOORPO + 1 ) = YE
         RMCN( MNDFPO + WOORPO + 2 ) = ZCTE
      ELSE IF( LEPLAN .EQ. 4 ) THEN
C        TRANSFERT DU PLAN DANS R3
         XY(1) = X
         XY(2) = Y
         CALL CH2D3D( P1,D2D3,XY,
     %                RMCN(MNDFPO+WOORPO) )
      ELSE
         RMCN( MNDFPO + WOORPO     ) = XE
         RMCN( MNDFPO + WOORPO + 1 ) = YE
         RMCN( MNDFPO + WOORPO + 2 ) = 0.0
      ENDIF
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNDFPO) )
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNDFPO + MOTVAR(6) ) = NONMTD( '~>POINT>>DEFINITION' )
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET'
C     =================================
      CALL LXTNDC( NTLXPO , 'XYZSOMMET' , 'MOTS'  , WYZSOM + 3 )
      CALL LXTSOU( NTLXPO , 'XYZSOMMET' ,  NTSOPO , MNSOPO )
      IF( NTSOPO .LE. 0 ) GOTO 9000
C
C     LE NOMBRE DE SOMMETS
      MCN( MNSOPO + WNBSOM ) = 1
C
C     X Y Z
      RMCN( MNSOPO + WYZSOM     ) = RMCN( MNDFPO + WOORPO )
      RMCN( MNSOPO + WYZSOM + 1 ) = RMCN( MNDFPO + WOORPO + 1 )
      RMCN( MNSOPO + WYZSOM + 2 ) = RMCN( MNDFPO + WOORPO + 2 )
C
C     MISE A JOUR DES EXTREMES
      COOEXT(1,1) = MIN( COOEXT(1,1) , RMCN( MNDFPO + WOORPO ) )
      COOEXT(1,2) = MAX( COOEXT(1,2) , RMCN( MNDFPO + WOORPO ) )
      COOEXT(2,1) = MIN( COOEXT(2,1) , RMCN( MNDFPO + WOORPO + 1 ) )
      COOEXT(2,2) = MAX( COOEXT(2,2) , RMCN( MNDFPO + WOORPO + 1 ) )
      COOEXT(3,1) = MIN( COOEXT(3,1) , RMCN( MNDFPO + WOORPO + 2 ) )
      COOEXT(3,2) = MAX( COOEXT(3,2) , RMCN( MNDFPO + WOORPO + 2 ) )
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNSOPO) )
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOPO + WNBTGS ) = 0
      MCN( MNSOPO + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     LA TRANSFORMATION DU POINT
C     ==========================
      IF( NUTRAN .GT. 1 ) THEN
          CALL TRPLSV( 1 , NUPOIN , NTSOPO , MNSOPO , IERR )
      ENDIF
C
C     TRACE DU POINT
      KSYMBO = 'x' // KNOMPT
      CALL SYMBOLE2D( NCGRIS, X, Y, KSYMBO )
C
C     LE POINT EST MIS DANS LE TABLEAU MNITEP
      IF( MCN(MNITEP+2) .GE. MCN(MNITEP+1) ) THEN
C        TABLEAU TROP PETIT:LA TAILLE DU TABLEAU EST AUGMENTEE
         CALL ITEMAU( MNITEP )
      ENDIF
C     LE NOMBRE D'ITEMS EST AUGMENTE DE 1
      MCN(MNITEP+2) = MCN(MNITEP+2) + 1
C     L'ADRESSE MCN DE L'ITEM
      MNI  = MNITEP + MCN(MNITEP) * MCN(MNITEP+2)
C     LES COORDONNEES PIXELS DE CET ITEM , SON CODE
      MCN( MNI     ) = NX
      MCN( MNI + 1 ) = NY
C     NUMERO DU POINT DANS SON LEXIQUE
      MCN( MNI + 2 ) = NUPOIN
C
C     RETOUR AU POINT GENERAL
      LDEPLA = 0
C
C     SAUVEGARDE DU POINT SUR LE FICHIER FRAPPE
C     =========================================
      KERR(MXLGER) = KNOMPT//' ; '// KNOMTR // ' ; '//'1;'
      CALL SANSDBL( KERR(MXLGER), I )
      WRITE(NFFRAP,*) KERR(MXLGER)(1:I)
      WRITE(KERR(MXLGER)(1:25) , '(G25.17)' ) RMCN( MNDFPO+WOORPO )
      KERR(MXLGER)(26:NBCAER) = ' ;  {X du POINT} '
      CALL SANSDBL( KERR(MXLGER), I )
      WRITE(NFFRAP,*) KERR(MXLGER)(1:I)
      WRITE(KERR(MXLGER)(1:25) , '(G25.17)' ) RMCN( MNDFPO+WOORPO+1 )
      KERR(MXLGER)(26:NBCAER) = ' ;  {Y du POINT} '
      CALL SANSDBL( KERR(MXLGER), I )
      WRITE(NFFRAP,*) KERR(MXLGER)(1:I)
      WRITE(KERR(MXLGER)(1:25) , '(G25.17)' ) RMCN( MNDFPO+WOORPO+2 )
      KERR(MXLGER)(26:NBCAER) = ' ;  {Z du POINT} '
      CALL SANSDBL( KERR(MXLGER), I )
      WRITE(NFFRAP,*) KERR(MXLGER)(1:I)
      GOTO 600
C
C     ERREUR
 9000 NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'SAPT2D: LE SUPER TABLEAU MCN EST SATURE'
         KERR(2) = 'SORTIR AVEC SAUVEGARDE et RELANCER MAILLER'
      ELSE
         KERR(1) = 'SAPT2D: SUPER ARRAY MCN is SATURATED'
         KERR(2) = 'EXIT WITH SAVE and EXECUTE AGAIN MAILLER'
      ENDIF
      CALL LEREUR
      GOTO 9999
C
C     GENERATION D'UNE REMONTEE VERS LE MENU DE DEPART
 9900 WRITE(NFFRAP,*) '@; '
 9999 RETURN
      END
