      SUBROUTINE TRVMTR31( MISTRE, KNOMOB, NBTYEL, MNELEM,
     %                     MODECO, NCAS,   CONTMN, CONTMX,
     %                     NBST,   MNCRIT, NOINTC,
     %                     NDPGST, MNPOGE, NBPOI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER EN 3D LE CRITERE DES CONTRAINTES DE VON MISES TRESCA
C -----    PAR COULEURS ARC EN CIEL DANS UNE COUPE PLANE DES EF 3D
C       pour  si la i-eme contrainte principale
C       Von MISES = sqrt( (s1-s2)**2 + (s2-s3)**2 + (s3-s1)**2 ) / sqrt(2)
C       TRESCA    = MAX( abs(s1-s2), abs(s2-s3), abs(s3-s1) ) / 2
C
C ENTREES :
C ---------
C MISTRE : 1 POUR CRITERE DE VON MISES
C          2 POUR CRITERE DE TRESCA
C KNOMOB : NOM DE L'OBJET
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS DANS CETTE TOPOLOGIE
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE
C MNELEM : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE
C MODECO : MODE DES VECTEURS SOLUTIONS
C          =1 CE SONT DEPLACEMENTS
C          =2 CE SONT DES MODES PROPRES
C NCAS   : NUMERO DU CAS A TRAITER
C CONTMN : CONTRAINTE MIN
C CONTMX : CONTRAINTE MAX
C NBST   : NOMBRE DE SOMMETS DES TYPES D'EF DU CRITERE
C MNCRIT : ADRESSE MCN DES TABLEAUX CRITERE(NBPIEX,NBELFI)
C NDPGST : CODE TRAITEMENT DES SOMMETS POINTS NOEUDS DU MAILLAGE
C          0 : NOEUDS=POINTS=SOMMETS
C          1 : NOEUDS=POINTS#SOMMETS
C          2 : NOEUDS#POINTS=SOMMETS
C          3 : NOEUDS#POINTS#SOMMETS
C MNPOGE : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES(=NOEUDS) DE L'OBJET
C MNNOEU : ADRESSE MCN DU TABLEAU NOEUDS D'INTERPOLATION DU MAILLAGE
C NBPOI  : NOMBRE TOTAL DE POINTS (=NOEUDS) DU MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Laboratoire J-L. LIONS UPMC PARIS    MAI 2007
C23456---------------------------------------------------------------012
      IMPLICIT INTEGER (W)
      PARAMETER     (LIGCON=0)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___face.inc"
      include"./incl/a___aretefr.inc"
      include"./incl/a___contrainte.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ctemps.inc"
      include"./incl/xyzext.inc"
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     KNOMOB
      CHARACTER*160     KNOM
      REAL              CONTMN, CONTMX
      INTEGER           NBST(4), MNCRIT(4), NOINTC(4)
      REAL              HXMIMX(6,2)
C
C     VARIABLES PERSO
      INTEGER           N
      REAL              NORM,EPS
C     HMMX POUR CALCUL DES MINMAX SUIVANTS DES DIRECTIONS
C          DEFINIES DANS POINTS (1:3,4) : TEMPORAIRE
      REAL              HMMX(2), POINTS(3,4)
C
C     POUR EVITER DES VECTEURS NORMAUX TROP PETITS
      EPS = 1E-8
C
      NBFPLA = 0
      MOTSFA = 0
      MOXYZP = 0
      MOREE2 = MOTVAR(6)
C
      LPEXIST = 0
      MNPLAV = 0
      MNSOLE = 0
      MNPILE = 0
      MNVALS = 0
      MNXYZP = 0
      MNTSFA = 0
      MNBARY = 0
      MNNUFA = 0
      NBF    = 0
      NBPLAN = 1
      MXPLAN = NBPLAN
      NBCOUL = NDCOUL - N1COUL + 1
C
C     PAR DEFAUT SECTION SELON Z
      NOAXE  = 3
C     NBCOOR : NOMBRE DE COORDONNEES DES NOEUDS=POINTS (3 ou 6)
      NBCOOR = 3
C
C     MIN ET MAX DES COORDONNEES DES POINTS DU MAILLAGE DE L'OBJET
      CALL MAJEXT( MNPOGE )
      CALL MIMXPT( NBCOOR, NBPOI, RMCN(MNPOGE+WYZPOI), HXMIMX )
      SECMIN = HXMIMX(NOAXE,1)
      SECMAX = HXMIMX(NOAXE,2)
C
C     PARAMETRES DES TAILLES DES TABLEAUX. A MODIFIER EVENTUELLEMENT
C     MXSOUI : MAXIMUM D INTERVALLES DU SEGMENT UNITE
C              MXSOUI ** 3 SOUS-TETRAEDRES CREES AU PLUS DANS UN TETRAEDRE
      MXSOUI = 4
C
C     CREATION OU REDECOUVERTE DU TMS OBJET>>>FACE
      CALL HACHOB( KNOMOB, 4, NTFAOB, MNFAOB, IERR )
C
C     CREATION DU HACHAGE DES ARETES DES FACES FRONTALIERES DE L'OBJET
      CALL HACHAF( KNOMOB, 0, NTFAOB, MNFAOB,
     %             NTAFOB, MNAFOB, I )
C
C     MXPILE : MAXIMUM DE SOUS-TETRAEDRES GENERES DANS L EF REFERENCE
      MXPILE = 6 * MXSOUI * MXSOUI * MXSOUI
C     NPILE(4,MXPILE) PILE DES 4 SOMMETS DE CHAQUE SOUS-TETRAEDRE
      MOPILE = 4 * MXPILE
      CALL TNMCDC( 'ENTIER', MOPILE, MNPILE )
C
C     SOLEL (MXNOEL) SOLUTION AUX NOEUDS DE L ELEMENT FINI DU MAILLAGE
C     FBASE (MXNOEL) VALEUR DES FONCTIONS DE BASE EN UN POINT
C     COPOE (MXPOEL,NBCOOR) COORDONNEES DES POINTS DE L EF COURANT
      MOSOLE = MXNOEL+MXNOEL+MXPOEL*NBCOOR
      CALL TNMCDC( 'REEL2', MOSOLE, MNSOLE )
      MNFBAS = MNSOLE + MOREE2 * MXNOEL
      MNCOPO = MNFBAS + MOREE2 * MXNOEL
C
C     MXSOMM : MAXIMUM DE SOMMETS DES SOUS-TETRAEDRES DE L EF REFERENCE
      MXSOMM = ( MXSOUI + 1 ) ** 3
C
C     VALST (5,MXSOMM) 3 COORDONNEES AUX SOMMETS DES SOUS-TETRAEDRES
C                      DE L'EF DE REFERENCE DANS R3
C                    + VALEUR DE LA SOLUTION AUX SOMMETS DES SOUS-TETRAEDRE
C                    + VALEUR DE LA COORDONNEE NOAXE SUR L'EF COURANT
      MOVALS = 5 * MXSOMM
      CALL TNMCDC( 'REEL', MOVALS, MNVALS )
C
C     PALETTE ARC EN CIEL
      CALL PALCDE( 11 )
C
C     RE-INITIALISER LA DONNEE DES PLANS
C     ----------------------------------
 80   IF( MOXYZP .GT. 0 ) CALL TNMCDS( 'REEL', MOXYZP, MNXYZP )
      IF( MOTSFA .GT. 0 ) CALL TNMCDS( 'REEL', MOTSFA, MNTSFA )
      MOXYZP = 0
      MOTSFA = 0
      MNXYZP = 0
      MNTSFA = 0
      IF( MNPLAV .GT. 0 ) THEN
C        DESTRUCTION DU TABLEAU DES VALEURS DES PLANS
         CALL TNMCDS( 'REEL', MXPLAN+1, MNPLAV )
      ENDIF
      MNPLAV = 0
      NBFPLA = 0
C     PAR DEFAUT SECTION SELON Z
      NOAXE  = 3
      LPEXIST = 0
C
C     LECTURE DES DONNEES
C     ===================
 100  CALL LIMTCL( 'sectplan', NMTCL )
      IF( NMTCL .LT.  0 ) GOTO 9000
      IF( NMTCL .EQ. 19 ) GOTO  190
      IF( NMTCL .EQ. 80 ) GOTO   80
      IF( NMTCL .EQ. 90 ) THEN
         IF( NBFPLA .LE. 0 ) THEN
            GOTO 200
         ELSE
            GOTO 300
         ENDIF
      ENDIF
      IF( NMTCL .GE. 1 .AND. NMTCL .LE. 3 ) THEN
C
C        COUPE ORTHOGONALE A L'AXE NMTCL
         NOAXE = NMTCL
         SECMIN = HXMIMX(NOAXE,1)
         SECMAX = HXMIMX(NOAXE,2)
         LPEXIST = 0
         GOTO 145
C
      ELSE IF( NMTCL .EQ. 4 ) THEN
C
C        COUPE SELON UN PLAN A DEFINIR
         CALL LIMTCL( 'defplan', NMTCL1 )
         IF( NMTCL1 .LE. 0 ) GOTO 9000
         IF( NMTCL1 .EQ. 1 ) THEN
C           DEFINITION DU PLAN PAR 3 POINTS
            N = 3
         ELSE
C           DEFINITION DU PLAN PAR 1 POINT + 1 VECTEUR NORMAL
            N = 2
         ENDIF
C
         DO I=1,N
C           ENTREE XI
            CALL INVITE( 98 )
            NCVALS = 0
            CALL LIRRSP( NCVALS, POINTS(1,I) )
            IF (NCVALS .LE. 0) GOTO 100
C
C           ENTREE YI
            NCVALS = 0
            CALL INVITE( 99 )
            CALL LIRRSP( NCVALS, POINTS(2,I) )
            IF (NCVALS .LE. 0) GOTO 100
C
C           ENTREE ZI
            NCVALS = 0
            CALL INVITE( 100 )
            CALL LIRRSP( NCVALS, POINTS(3,I) )
            IF (NCVALS .LE. 0) GOTO 100
         ENDDO
C
C        POUR 3 POINTS => POINTS + VECTEUR PAR PRODUIT VECTORIEL
         IF (N .EQ. 3) THEN
            DO I=1,3
               POINTS(I,3) = POINTS(I,3) - POINTS(I,1)
               POINTS(I,4) = POINTS(I,2) - POINTS(I,1)
            ENDDO
            CALL PROVER( POINTS(1,4), POINTS(1,3), POINTS(1,2) )
         ENDIF
C
C        NORMALISATION DU VECTEUR NORMAL
         NORM = SQRT(PROSCR(POINTS(1,2),POINTS(1,2),3))
C
         IF (NORM .LT. EPS) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'TRVMTR31: VECTEUR NORMAL AU PLAN NUL'
               KERR(2) = '          VEUILLEZ RECOMMENCER'
            ELSE
               KERR(1) = 'TRVMTR31: NULL NORMAL VECTOR'
               KERR(2) = '          GIVE IT AGAIN'
            ENDIF
            CALL LEREUR
            GOTO 100
         ENDIF
C
         DO I=1,3
            POINTS(I,2) = POINTS(I,2)/NORM
         ENDDO
C
         CALL MMPSPT( NBCOOR, NBPOI, POINTS, RMCN(MNPOGE+WYZPOI), HMMX )
         SECMIN = HMMX(1)
         SECMAX = HMMX(2)
         NOAXE = 4
         LPEXIST = 0
         GOTO 145
C
      ELSE IF( NMTCL .GT. 14 ) THEN
C        ERREUR
         GOTO 100
      ENDIF
C
      IF( NMTCL .EQ. 5 ) THEN
C
C        COULEUR des ARETES de la FRONTIERE
C        ..................................
         CALL LIMTCL( 'couleur0' , I )
         IF( I .EQ. -1 ) GOTO 100
         IF( I .EQ. -2 ) THEN
C           COULEUR INVISIBLE A NE PAS TRACER
            NCOAFR = -2
         ELSE IF( I .EQ. 0 ) THEN
C           LA COULEUR NOIRE
            NCOAFR = 0
         ELSE
            NCOAFR = N1COEL + I
         ENDIF
         GOTO 100
C
      ENDIF
C
      IF( NMTCL .EQ. 6 ) THEN
C
C        TYPE du TRAIT des ARETES de la FRONTIERE
C        ........................................
         CALL LIMTCL( 'typtrait' , I )
         IF( I .EQ. -1 ) GOTO 100
         NTLAFR = I
         GOTO 100
C
      ENDIF
C
      IF( NMTCL .EQ. 7 ) THEN
C
C        COULEUR des ARETES dans le plan de SECTION
C        ..........................................
         CALL LIMTCL( 'couleur0' , I )
         IF( I .EQ. -1 ) GOTO 100
         IF( I .EQ. -2 ) THEN
C           COULEUR INVISIBLE A NE PAS TRACER
            NCOAPL = -2
         ELSE IF( I .EQ. 0 ) THEN
C           LA COULEUR NOIRE
            NCOAPL = 0
         ELSE
            NCOAPL = N1COEL + I
         ENDIF
         GOTO 100
C
      ENDIF
C
      IF( NMTCL .EQ. 8 ) THEN
C
C        TYPE du TRAIT des ARETES dans le PLAN de SECTION
C        ................................................
         CALL LIMTCL( 'typtrait' , I )
         IF( I .EQ. -1 ) GOTO 100
         NTLAPL = I
         GOTO 100
C
      ENDIF
      GOTO( 110, 130, 140, 150 ) , NMTCL-10
C
C     NOMBRE DE PLANS DE SECTION
C     ...........................
 110  IF( MNPLAV .GT. 0 ) THEN
C        DESTRUCTION DU TABLEAU DES VALEURS DES PLANS
         CALL TNMCDS( 'REEL', MXPLAN+1, MNPLAV )
      ENDIF
      CALL INVITE( 71 )
      NCVALS = 4
      CALL LIRENT( NCVALS, NBPLAN )
      IF( NCVALS .LE. 0 ) GOTO 100
      NBPLAN = ABS( NBPLAN )
C     AU MOINS UN PLAN
      NBPLAN  = MAX( 1, NBPLAN )
      MXPLAN  = NBPLAN
      CALL TNMCDC( 'REEL', MXPLAN+1, MNPLAV )
C     PAR DEFAUT PLANS DISTANTS REGULIEREMENT
C
C     DISTANCE REGULIERE ENTRE LES PLANS
C     ..................................
 130  X = ( SECMAX - SECMIN ) / 100.
      IF( X .EQ. 0 ) X = 1.0
      SECMIN = SECMIN + X
      SECMAX = SECMAX - X
      GOTO 145
C
C     COORDONNEE MINIMALE ET MAXIMALE DES PLANS DE SECTION
C     ....................................................
 140  NCVALS = 0
      CALL INVITE( 11 )
      CALL LIRRSP( NCVALS, SECMIN )
      IF( NCVALS .LE. 0 ) GOTO 100
      CALL INVITE( 10 )
      NCVALS = 0
      CALL LIRRSP( NCVALS, SECMAX )
      IF( NCVALS .LE. 0      ) GOTO 100
      IF( SECMIN .GT. SECMAX ) GOTO 140
C
 145  I = NBPLAN - 1
      IF( I .EQ. 0 ) I = 1
      X = ( SECMAX - SECMIN ) / I
      F = SECMIN
      DO 148 I=1,NBPLAN
         RMCN( MNPLAV - 1 + I ) = F
         F                      = F + X
 148  CONTINUE
      RMCN( MNPLAV + NBPLAN ) = SECMAX
      LPEXIST = 0
      IF( NBPLAN .EQ. 1 ) RMCN(MNPLAV) = (SECMIN + SECMAX) * 0.5
      IF( (NMTCL.GE.1 .AND. NMTCL.LE.4) .OR. NMTCL .EQ. 11 ) GOTO 100
      GOTO 200
C
C     COORDONNEE DES PLANS DE SECTION
C     ................................
 150  KNOM = 'COORDONNEE DU PLAN'
      IF( MNPLAV .EQ. 0 ) THEN
C        INITIALISATIONS PAR DEFAUT
         NBPLAN = ABS( NBPLAN )
         NBPLAN = MAX( 1, NBPLAN )
         MXPLAN = NBPLAN
         CALL TNMCDC( 'REEL', MXPLAN+1, MNPLAV )
      ENDIF
      DO 158 I=1,NBPLAN
         WRITE( KNOM(19:24), '(I6)' ) I
         CALL INVITD( KNOM )
         NCVALS = 0
         CALL LIRRSP( NCVALS, RMCN(MNPLAV-1+I) )
         IF( NCVALS .LE. 0 ) GOTO 100
 158  CONTINUE
      CALL TRIREE( NBPLAN, RMCN(MNPLAV) )
      GOTO 200
C
C     POURCENTAGE de REDUCTION des FACES
C     ..................................
 190  CALL INVITE( 130 )
      CALL LIRRSP( NCVALS , R )
      IF( NCVALS .LE. 0 ) GOTO 100
      PREDUF = R
      PREDUF = MIN( 100.0 , PREDUF )
      PREDUF = MAX(   0.0 , PREDUF )
      NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:10),'(G10.2)' ) PREDUF
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = '% de REDUCTION des FACES= ' // KERR(MXLGER)(1:10)
      ELSE
         KERR(1) = '% of FACE REDUCTION= ' // KERR(MXLGER)(1:10)
      ENDIF
      CALL LERESU
      GOTO 100
C
C     EXECUTION DU TRACE DES FACES INTERSECTION AVEC LES PLANS
C     ========================================================
 200  IF( MNPLAV .EQ. 0 ) THEN
C        INITIALISATIONS PAR DEFAUT
         NBPLAN = ABS( NBPLAN )
         NBPLAN = MAX( 1, NBPLAN )
         MXPLAN = NBPLAN
         CALL TNMCDC( 'REEL', MXPLAN+1, MNPLAV )
C        MIN ET MAX DE LA COORDONNEE NOAXE
         IF (NOAXE .EQ. 4) THEN
C           C'EST DU SURPLUS A OPTIMISER
            CALL MMPSPT(NBCOOR, NBPOI, POINTS, RMCN(MNPOGE+WYZPOI),HMMX)
            SECMIN = HMMX(1)
            SECMAX = HMMX(2)
         ELSE
            SECMIN = HXMIMX(NOAXE,1)
            SECMAX = HXMIMX(NOAXE,2)
         ENDIF
         X = ( SECMAX - SECMIN ) / 100.
         IF( X .EQ. 0 ) X = 1.0
         SECMIN = SECMIN + X
         SECMAX = SECMAX - X
         I = NBPLAN - 1
         IF( I .LE. 0 ) I=1
         X    = ( SECMAX - SECMIN ) / I
         F    = SECMIN
         DO 210 I=1,NBPLAN
            RMCN( MNPLAV - 1 + I ) = F
            F                      = F + X
 210     CONTINUE
         IF( NBPLAN .EQ. 1 ) RMCN(MNPLAV) = (SECMIN + SECMAX) * 0.5
      ENDIF
C
      IF (NOAXE .EQ. 4) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10310)
     %        NCAS, (POINTS(I,1),I=1,3), (POINTS(I,2),I=1,3),
     %        (I,RMCN(MNPLAV-1+I),I=1,NBPLAN)
         ELSE
            WRITE(IMPRIM,20310)
     %        NCAS, (POINTS(I,1),I=1,3), (POINTS(I,2),I=1,3),
     %        (I,RMCN(MNPLAV-1+I),I=1,NBPLAN)
         ENDIF
      ELSE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10210)
     %        NCAS, NOAXE, (I,RMCN(MNPLAV-1+I),I=1,NBPLAN)
         ELSE
            WRITE(IMPRIM,20210)
     %        NCAS, NOAXE, (I,RMCN(MNPLAV-1+I),I=1,NBPLAN)
         ENDIF
      ENDIF
C
10210 FORMAT(/'CAS ',I3,
     %' : VALEURS DE LA COORDONNEE ',I1,' CONSTANTE DES PLANS'/
     % 5(I4,' : ',G13.5))
20210 FORMAT(/'CASE ',I3,
     %' : VALUE of the COORDINATE ',I1,' CONSTANT of PLANES'/
     % 5(I4,' : ',G13.5))
10310 FORMAT(/'CAS ',I3,
     %     ' : VALEURS DE LA COORDONNEE RELATIVES AU PLAN :'/
     %     ' POINT=(', 3(G13.5), ')'/
     %     ' VECTEUR NORMAL=(', 3(G13.5),')'/
     %     5(I3,' : ',G13.5))
20310 FORMAT(/'CASE ',I3,
     %     ' : VALUE of the RELATIVE COORDINATE at PLANE :'/
     %     ' POINT=(', 3(G13.5), ')'/
     %     ' NORMAL VECTOR=(', 3(G13.5),')'/
     %     5(I3,' : ',G13.5))
C
C     CREATION DU TABLEAU DES XYZ DES SOMMETS DANS LES PLANS DE SECTION
      NBCOPS = 4
      MXFPLA  = NBFPLA + 2 * NBPOI * NBPLAN
      MOXYZP1 = NBCOPS * 4 * (NBFPLA+MXFPLA)
      IF( MNXYZP .LE. 0 ) THEN
         CALL TNMCDC( 'REEL', MOXYZP1, MNXYZP )
      ELSE
         CALL TNMCAU( 'REEL', MOXYZP, MOXYZP1, NBCOPS*4*NBFPLA, MNXYZP )
      ENDIF
      MOXYZP = MOXYZP1
      MOTSFA1 = 4 * (NBFPLA+MXFPLA)
      IF( MNTSFA .LE. 0 ) THEN
         CALL TNMCDC( 'REEL', MOTSFA1, MNTSFA )
      ELSE
         CALL TNMCAU( 'REEL', MOTSFA, MOTSFA1, 4*NBFPLA, MNTSFA )
      ENDIF
      MOTSFA = MOTSFA1
C
C     BOUCLE SUR LES DIFFERENTS TYPES D'ELEMENTS FINIS DU MAILLAGE
C     ============================================================
      DO 250 I = 0, NBTYEL-1
C
C        L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF
         MNELE = MCN( MNELEM + I )
C
C        LE NUMERO DU TYPE DES ELEMENTS FINIS
         NUTYEL = MCN( MNELE + WUTYEL )
         IF( NUTYEL .EQ. 30 ) GOTO 250
C        PAS DE TRACE DES 6CUBES, SEULEMENT DES 3Q1C EXTRAITS
C
C        LE NOMBRE DE TELS ELEMENTS
         NBELEM = MCN( MNELE + WBELEM )
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELTYCA( NUTYEL )
C
C        L'ADRESSE MCN DU TABLEAU 'NPEF' POUR CE TYPE D'EF
         MNPGEL = MNELE + WUNDEL
         IF( NDPGST .GE. 2 ) THEN
            MNPGEL = MNPGEL + MCN(MNELE+WBELEM) * MCN(MNELE+WBNDEL)
         ENDIF
C
C        NOMBRE DE SOMMETS POUR LE CRITERE
         NBPIEX = NBST(I+1)
C
C        CALCUL DES FACES APPARTENANT AUX PLANS DE SECTION DE L'OBJET
         CALL TRVMTR32( NOAXE,  POINTS, NBPLAN, RMCN(MNPLAV), NBCOOR,
     %                  NOINTC(I+1), NBPIEX, NBELEM, RMCN(MNCRIT(I+1)),
     %                  NBPOE,  MCN(MNPGEL), NBPOI,  MCN(MNPOGE+WYZPOI),
     %                  MXSOMM, MCN(MNSOLE), MCN(MNCOPO),
     %                  MXPILE, MCN(MNPILE), MCN(MNVALS), MCN(MNFBAS),
     %                  NBCOPS, TRIANG, MXFPLA, NBFPLA,
     %                  RMCN(MNXYZP), RMCN(MNTSFA) )
 250  CONTINUE
C     LES SECTIONS ONT ETE CALCULEES
      LPEXIST = 1
C
      IF( NMTCL .EQ. 90 ) GOTO 300
      GOTO 100
C
C     -------------------------------------------------------------
C     OPTIONS DE LA VISEE POUR VOIR L'OBJET ET LES PLANS DE SECTION
C     -------------------------------------------------------------
 300  IF( LPEXIST .EQ. 0 ) GOTO 200
      IF( NBFPLA  .LE. 0 ) GOTO 100
C
      CALL VISE3D( NMTCL )
      IF( NMTCL .LT. 0 ) GOTO 100
C
C     LE NOMBRE D'ENTIERS PAR ARETE FRONTALIERE
      MOARFR = MCN( MNAFOB + WOARFR )
C     LA MAJORATION DU NOMBRE DES ARETES FRONTALIERES
      MXARFR = MCN( MNAFOB + WXARFR )
C     LE NUMERO DANS LAREFR DE LA PREMIERE ARETE FRONTALIERE
      L1ARFR = MCN( MNAFOB + W1ARFR )
C     LE NOMBRE D'ARETES FRONTALIERES DANS LE CHAINAGE
      NBARFR = MCN( MNAFOB + WBARFR )
C
      WRITE(IMPRIM,*) 'NOMBRE DE FACES DANS LES PLANS=',NBFPLA
      WRITE(IMPRIM,*) 'NOMBRE D''ARETES FRONTALIERES =',NBARFR
C
      IF( MNBARY .GT. 0 ) CALL TNMCDS( 'REEL',   NBF, MNBARY )
      IF( MNNUFA .GT. 0 ) CALL TNMCDS( 'ENTIER', NBF, MNNUFA )
      MNBARY = 0
      MNNUFA = 0
C
      NBF = NBFPLA + NBARFR
      IF( MNBARY .EQ. 0 ) CALL TNMCDC( 'REEL',   NBF, MNBARY )
      IF( MNNUFA .EQ. 0 ) CALL TNMCDC( 'ENTIER', NBF, MNNUFA )
C
C     INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
      IF( LORBITE .NE. 0 ) THEN
         CALL ORBITE0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 300
      ENDIF
C
 400  IF( LORBITE .NE. 0 ) THEN
C        ORBITE OU ZOOM OU TRANSLATION ACTIFS
         CALL ORBITE1( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 300
      ENDIF
C
C     FORMATION DES TABLEAUX NUMERO DES FACES ET DISTANCE AXONOMETRIQUE
C     POUR LES FACES DES PLANS DE SECTION
      CALL TRBFIS( TRIANG, NBFPLA, NBCOPS, 0, RMCN(MNXYZP),
     %             MCN(MNNUFA), RMCN(MNBARY) )
C
C     FORMATION DES TABLEAUX NUMERO DES FACES ET DISTANCE AXONOMETRIQUE
C     POUR LES ARETES FRONTALIERES DE L'OBJET
C     ATTENTION LES ARETES FRONTALIERES DOIVENT OBLIGATOIREMENT ETRE APRES
C               LES FACES DES PLANS DE SECTION
      CALL TRBARF( MOARFR, MXARFR, L1ARFR, MCN(MNAFOB+WAREFR),
     %             MCN(MNPOGE+WBCOOP), RMCN(MNPOGE+WYZPOI),
     %             NBFPLA, MCN(MNNUFA), RMCN(MNBARY))
C
C     LE TRI PAR TAS DE CETTE DISTANCE
C     LA FACE LA PLUS PROCHE EST LA PREMIERE
      CALL TRITRP( NBF, RMCN(MNBARY), MCN(MNNUFA) )
C
C     TRACE EFFECTIF DES ARETES FRONTALIERES ET DES FACES ISOVALEURS
C     --------------------------------------------------------------
      CALL TRAXE3
      CALL T3AFP3( NBARFR, NBFPLA, MCN(MNNUFA),
     %             MOARFR, MXARFR, MCN(MNAFOB+WAREFR),
     %             RMCN(MNPOGE+WYZPOI),
     %             TRIANG, NBCOPS, RMCN(MNXYZP), RMCN(MNTSFA),
     %             NCAS,   NCAS,   NCAS,
     %             CONTMN, CONTMX )
C
C     TRACE DU POINT DE TEMPERATURE MINIMALE ET MAXIMALE
      CALL XVTYPETRAIT( LIGCON )
C
C     EFFACEMENT DE LA LEGENDE SUR POSTSCRIPT
      IF ( LASOPS .NE. 0 ) THEN
        IF ( LASOPS .EQ. 1 ) THEN
          LASOPS = -11
        ELSE
          IF ( LASOPS .EQ. 2 ) THEN
            LASOPS = -12
          ELSE
            LASOPS = 0
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'TRVMTR31: MAUVAISE VALEUR DE LASOPS'
               KERR(2) = '          ARRET DU POSTSCRIPT'
            ELSE
               KERR(1) = 'TRVMTR31: BAD VALUE of LASOPS'
               KERR(2) = '          STOP of POSTSCRIPT'
            ENDIF
            CALL LEREUR
            GOTO 100
          ENDIF
        ENDIF
        CALL XVPOSTSCRIPT(LASOPS)
        LASOPS = - LASOPS
        CALL XVPOSTSCRIPT(LASOPS)
      ENDIF
C
C     LE TRACE DU TITRE FINAL
C     =======================
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'OBJET: ' // KNOMOB
      ELSE
         KNOM = 'OBJECT: ' // KNOMOB
      ENDIF
      I = NUDCNB( KNOM )
      CALL XVCOULEUR( NCGRIS )
      CALL XVTEXTE( KNOM(1:I), I, 50, 30 )
C
C     LE TRACE DE LA LEGENDE : COULEURS => VALEURS
      NBCOUL = NDCOUL - N1COUL
      NCPAS  = NBCOUL / 10
      TPAS   = (CONTMX-CONTMN) / 10
      T      = CONTMN
C
C     TRACE DE 11 VALEURS
      NCOUL = N1COUL
      NX    = LAPXFE - 170
      NY    = LHPXFE - 30
      DO 60 I=0,10
         CALL XVCOULEUR( NCOUL )
         CALL XVRECTANGLE( NX, NY, 30, 10 )
         WRITE( KNOM(1:10), '(G10.3)' ) T
         CALL XVTEXTE( KNOM(1:10), 10, NX+40, NY+10 )
         NCOUL = NCOUL + NCPAS
         T     = T  + TPAS
         NY    = NY - 15
 60   CONTINUE
C
C     RETOUR AU TRACE NORMAL POUR POSTSCRIPT
      IF ( LASOPS .NE. 0 ) THEN
        LASOPS = LASOPS - 10
        CALL XVPOSTSCRIPT( LASOPS )
      ENDIF
C
C     FIN DU TRACE
      IAVTIT = 1
      WRITE( KNOM(1:4), '(I4)' ) NCAS
      WRITE( KNOM(5:19), '(G15.6)' ) TEMPS
      IF( LANGAG .EQ. 0 ) THEN
         IF( MISTRE .EQ. 1 ) THEN
            IF( MODECO .EQ. 1 ) THEN
               CALL TRFINS( 'CRITERE de PLASTICITE de Von MISES: Cas'
     %                     //KNOM(1:4) // ' au temps ' // KNOM(5:19) )
            ELSE
               CALL TRFINS( 'CRITERE de PLASTICITE de Von MISES: '//
     %    'Frequence Propre '// KNOM(1:4) //': ' // KNOM(5:19) // ' Hz')
            ENDIF
         ELSE
            IF( MODECO .EQ. 1 ) THEN
               CALL TRFINS( 'CRITERE de PLASTICITE de TRESCA: Cas'
     %                      //KNOM(1:4) // ' au temps ' // KNOM(5:19) )
            ELSE
               CALL TRFINS( 'CRITERE de PLASTICITE de TRESCA: '//
     %    'Frequence Propre '// KNOM(1:4) //': ' // KNOM(5:19) // ' Hz')
            ENDIF
         ENDIF
      ELSE
         IF( MISTRE .EQ. 1 ) THEN
            IF( MODECO .EQ. 1 ) THEN
               CALL TRFINS( 'Von MISES''s PLASTICITY Criterion: Case'
     %                      //KNOM(1:4) // ' at time ' // KNOM(5:19) )
            ELSE
               CALL TRFINS( 'Von MISES''s PLASTICITY Criterion: '//
     %    'EigenFrequency '// KNOM(1:4) //': ' // KNOM(5:19) // ' Hz')
            ENDIF
         ELSE
            IF( MODECO .EQ. 1 ) THEN
               CALL TRFINS( 'TRESCA''s PLASTICITY Criterion: Case'
     %                     //KNOM(1:4) // ' at time ' // KNOM(5:19) )
            ELSE
               CALL TRFINS( 'TRESCA''s PLASTICITY Criterion: '//
     %    'EigenFrequency '// KNOM(1:4) //': ' // KNOM(5:19) // ' Hz')
            ENDIF
         ENDIF
      ENDIF
C
C     RETOUR POUR UNE NOUVELLE VISEE
      IF( LORBITE .NE. 0 ) GOTO 400
C
      CALL CLICSO
      GOTO 300
C
C     SORTIE DU TRACE DES PLANS DE SECTION
C     ====================================
 9000 IF( MNPLAV .GT. 0 ) THEN
         CALL TNMCDS( 'REEL', MXPLAN+1, MNPLAV )
      ENDIF
      IF( MNBARY .GT. 0 ) CALL TNMCDS( 'REEL',   NBF,    MNBARY )
      IF( MNNUFA .GT. 0 ) CALL TNMCDS( 'ENTIER', NBF,    MNNUFA )
      IF( MNSOLE .GT. 0 ) CALL TNMCDS( 'REEL2' , MOSOLE, MNSOLE )
      IF( MNPILE .GT. 0 ) CALL TNMCDS( 'ENTIER', MOPILE, MNPILE )
      IF( MNVALS .GT. 0 ) CALL TNMCDS( 'REEL'  , MOVALS, MNVALS )
      IF( MNXYZP .GT. 0 ) CALL TNMCDS( 'REEL'  , MOXYZP, MNXYZP )
      IF( MNTSFA .GT. 0 ) CALL TNMCDS( 'REEL'  , MOTSFA, MNTSFA )
C     RETOUR AUX PARAMETRES INITIAUX
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
      RETURN
      END
