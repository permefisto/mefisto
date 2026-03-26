      SUBROUTINE TRPLSE( MODSEC, KNOMOB, MODECO,
     %                   NBTYEL, MNELEM, MNPOGE, NDPGST,
     %                   NCAS0,  NCAS1,  NTDL,   NTYP, TEMPER, dptemp,
     %                   TMIN,   NOEMIN, NCAMIN, TMAX, NOEMAX, NCAMAX,
     %                   TIMES )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DE LA TEMPERATURE PAR ZONES DE COULEURS DANS DES PLANS
C -----    X OU Y OU Z = CONSTANTE  OU PLAN QUELCONQUE
C          SUR UN OBJET 3D ou 6D
C          SOIT DANS LE PLAN, SOIT ORTHOGONALEMENT AU PLAN
C ENTREES:
C --------
C MODSEC : 0 TRACE DANS LES PLANS DE SECTION
C          1 TRACE ORTHOGONALEMENT AUX PLANS DE SECTION (PROFIL)
C NOPROJ : TYPE DE PROJECTION
C          0 FIXE LA COORDONNEE A ZERO
C         -1 PAS DE PROJECTION TRAITEMENT en XYZ NORMAL
C          1 : 'X Y Z 0 0 0'
C          2 : 'X Y 0 U 0 0'
C          3 : 'X 0 0 U V 0'
C          4 : '0 0 0 U V W'
C KNOMOB : NOM DE L'OBJET
C MODECO : MODE DE TRACE DES VECTEURS DU TMS D'ADRESSE MNDEPL
C         = 1 CE SONT DES TEMPERATURES
C         = 2 CE SONT DES MODES PROPRES
C         = 3 CE SONT DES PRESSIONS P1 A PARTIR D'UNE INTERPOLATION
C             SOIT en 2D P2 SOIT P1+BULLE P3 OU en 3D P1 OU en 3D P2
C         = 4 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C         = 5 CAS PARTICULIER DU 1D TRACE DES VECTEURS TEMPERATURES
C         = 6 CE SONT LES NORMES D'UNE VITESSE
C         = 7 FONCTION COURANT D'UN FLUIDE
C         = 8 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C             DU MODULE D'UNE ONDE COMPLEXE
C            (CALCUL DU MODULE DE L'ERREUR COMPLEXE A FAIRE ICI)
C         = 9  CE SONT DES PARTIES REELLES     D'UNE ONDE COMPLEXE
C         = 10 CE SONT DES PARTIES IMAGINAIRES D'UNE ONDE COMPLEXE
C         = 11 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C              DU MODULE DE LA VITESSE D'UN FLUIDE
C         = 12 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C         = 13 FLECHES REPRESENTANT DES VECTEURS VITESSE
C         = 14 ROTATIONNEL DE VITESSES ou TOURBILLONS ou VORTICITES

C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS DANS CETTE TOPOLOGIE
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE
C MNELEM : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE
C MNPOGE : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C MNNOEU : ADRESSE MCN DU TABLEAU NOEUDS D'INTERPOLATION DU MAILLAGE
C NDPGST : CODE TRAITEMENT DES SOMMETS POINTS NOEUDS DU MAILLAGE
C          0 : NOEUDS=POINTS=SOMMETS
C          1 : NOEUDS=POINTS#SOMMETS
C          2 : NOEUDS#POINTS=SOMMETS
C          3 : NOEUDS#POINTS#SOMMETS

C NCAS0  : NUMERO DU PREMIER CAS OU VECTEUR SOLUTION A TRACER
C NCAS1  : NUMERO DU DERNIER CAS OU VECTEUR SOLUTION A TRACER
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE EN TEMPERATURE
C NTYP   : =0 EMPLOI de TEMPER
C          =1 EMPLOI de dptemp
C TEMPER : TABLEAU DES NTDL * NCAS0:NCAS1 TEMPERATURES
C dptemp : TABLEAU DES NCAS0:NCAS1 TABLEAU(NTDL) TEMPERATURES

C TMIN   : TEMPERATURE MINIMALE DES CAS NCAS0:NCAS1
C NOEMIN : NUMERO DU NOEUD OU LA TEMPERATURE EST MINIMALE
C NCAMIN : NO DU CAS NCAS0:NCAS1 OU LA TEMPERATURE EST MINIMALE
C TMAX   : TEMPERATURE MAXIMALE DES CAS NCAS0:NCAS1
C NOEMAX : NUMERO DU NOEUD OU LA TEMPERATURE EST MAXIMALE
C NCAMAX : NO DU CAS NCAS0:NCAS1 OU LA TEMPERATURE EST MAXIMALE
C TIMES  : TEMPS DU CALCUL DES CAS NCAS0:NCAS1 des VECTEURS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C MODIFS : ALAIN PERRONNET TEXAS A & M UNIVERSITY           JUILLET 2005
C MODIFS : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  OCTOBRE 2010
C MODIFS : ALAIN PERRONNET SAINT PIERRE DU PERRAY           FEVRIER 2021
C23456---------------------------------------------------------------012
      PARAMETER ( LIGCON=0, TRIANG=1E25 )
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___face.inc"
      include"./incl/a___aretefr.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ctemps.inc"
      include"./incl/traaxe.inc"
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))

      DOUBLE PRECISION  TEMPER(NTDL,NCAS0:NCAS1)
      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(NCAS0:NCAS1) :: dptemp

      CHARACTER*(*)     KNOMOB
      CHARACTER*120     KNOM
      CHARACTER*24      NOMFGIF

      REAL              TIMES( NCAS0:NCAS1 ),
     %                  HXMIMX(6,2), NORMALE(3), VNORMAL(3),
     %                  NORM, EPS

      REAL, allocatable, dimension(:):: XYZSFP, TEMSFP

C     HMMX POUR CALCUL DES MINMAX SUIVANTS DES DIRECTIONS
C          DEFINIES DANS POINTS (1:3,4) : TEMPORAIRE
      REAL              HMMX(2), POINTS(3,4)

C     POUR EVITER DES VECTEURS NORMAUX TROP PETITS
      EPS = 1E-8

C     SAUVEGARDE DES MIN MAX DE LA SOLUTION
      TMIN0 = TMIN
      TMAX0 = TMAX

      MOREE2 = MOTVAR(6)
      IERXYZSFP = 1
      IERTEMSFP = 1
      NBFPLA = 0
      LPEXIST= 0
      MNPLAV = 0
      MNSOLE = 0
      MNPILE = 0
      MNVALS = 0
      MOVALS = 0
      MOXYZSFP= 0
      MNXYZIP= 0
      MOTSFA = 0
      MNBARY = 0
      MNNUFA = 0
      NBF    = 0
      NBPLAN = 1
      MXPLAN = NBPLAN
C     COULEUR DES ARETES DANS LES PLANS
ccc      NCOPLAN = -2 => NON TRACE
      NCOPLAN = NCGRIS
      NTLPLAN = LIGCON
      NBCOUL = NDCOUL - N1COUL + 1
C
C     PAR DEFAUT SECTION SELON Z
      NOAXE  = 3
C     PAR DEFAUT NORMALE AU PLAN Z
      NORMALE(1) = 0
      NORMALE(2) = 0
      NORMALE(3) = 1.0
C     COEFFICIENT D'AMPLIFICATION
      AMPLI = 1.0
C
C     NBPOI  NOMBRE DE POINTS DU MAILLAGE DE L'OBJET
      NBPOI = MCN(MNPOGE+WNBPOI)
C
C     NBCOOR : NOMBRE DE COORDONNEES DES NOEUDS=POINTS (3 ou 6)
      NBCOOR = MCN( MNPOGE + WBCOOP )
C
C     MIN ET MAX DES COORDONNEES DES POINTS DU MAILLAGE DE L'OBJET
      CALL MAJEXT( MNPOGE )
      CALL MIMXPT( NBCOOR, NBPOI, RMCN(MNPOGE+WYZPOI), HXMIMX )
      SECMIN = HXMIMX(NOAXE,1)
      SECMAX = HXMIMX(NOAXE,2)
C     LONGUEUR DE LA DIAGONALE DANS R3
      DIAGON = SQRT( (HXMIMX(1,2)-HXMIMX(1,1))**2
     %             + (HXMIMX(2,2)-HXMIMX(2,1))**2
     %             + (HXMIMX(3,2)-HXMIMX(3,1))**2 )
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
      IF( MNPILE .LE. 0 ) GOTO 9999
C
C     SOLEL (MXNOEL) SOLUTION AUX NOEUDS DE L ELEMENT FINI DU MAILLAGE
C     FBASE (MXNOEL) VALEUR DES FONCTIONS DE BASE EN UN POINT
C     COPOE (MXPOEL,NBCOOR) COORDONNEES DES POINTS DE L EF COURANT
      MOSOLE = MXNOEL+MXNOEL+MXPOEL*NBCOOR
      CALL TNMCDC( 'REEL2', MOSOLE, MNSOLE )
      IF( MNSOLE .LE. 0 ) GOTO 9999
      MNFBAS = MNSOLE + MOREE2 * MXNOEL
      MNCOPO = MNFBAS + MOREE2 * MXNOEL
C
C     MXSOMM : MAXIMUM DE SOMMETS DES SOUS-TETRAEDRES DE L EF REFERENCE
      MXSOMM = ( MXSOUI + 1 ) ** 3
C
C     VALST (4,MXSOMM) 3 COORDONNEES AUX SOMMETS DES SOUS-TETRAEDRES
C                      DE L'EF DE REFERENCE DANS R3
C                    + VALEUR DE LA COORDONNEE NOAXE SUR L'EF COURANT
      IF( MNVALS .GT. 0 ) CALL TNMCDS( 'REEL', MOVALS, MNVALS )
      MOVALS = 4 * MXSOMM
      CALL TNMCDC( 'REEL', MOVALS, MNVALS )
      IF( MNVALS .LE. 0 ) GOTO 9999
C
      IF( MODECO .NE. 4  .AND. MODECO .NE. 8   .AND.
     %    MODECO .NE. 11 .AND. MODECO .NE. 12 ) THEN
C        PALETTE ARC EN CIEL
         CALL PALCDE( 11 )
      ELSE
C        PALETTE JAUNE ROUGE MARRON NOIR DES ERREURS
         CALL PALCDE( 5 )
      ENDIF
C
C     RE-INITIALISER LA COORDONNEE DES PLANS
C     --------------------------------------
 80   IF( MNPLAV .GT. 0 ) THEN
C        DESTRUCTION DU TABLEAU DES VALEURS DES PLANS
         CALL TNMCDS( 'REEL', MXPLAN+1, MNPLAV )
      ENDIF
      MNPLAV = 0
      NBFPLA = 0
C     PAR DEFAUT SECTION SELON Z
      NOAXE   = 3
      LPEXIST = 0
C
C     LECTURE DES DONNEES DE L'UTILISATEUR
C     ====================================
 100  IF( MODSEC .EQ. 0 ) THEN
C        TRACE DANS LES PLANS DE SECTION
         CALL LIMTCL( 'sectplan', NMTCL )
C
C        NOM DU FICHIER VIDEO  SELON MODECO // 'sepl'
         CALL VIDEONM( MODECO, 'sepl', NOMFGIF )
C
      ELSE
C        TRACE ORTHOGONALEMENT AUX PLANS DE SECTION
         CALL LIMTCL( 'profplan', NMTCL )
C
C        NOM DU FICHIER VIDEO  SELON MODECO // 'prpl'
         CALL VIDEONM( MODECO, 'prpl', NOMFGIF )
C
      ENDIF
C
      IF( NMTCL .LT.  0 ) GOTO 9000
      IF( NMTCL .EQ. 19 ) GOTO  190
      IF( NMTCL .EQ. 20 ) GOTO  195
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
         NOAXE  = NMTCL
         SECMIN = HXMIMX(NOAXE,1)
         SECMAX = HXMIMX(NOAXE,2)
C        NORMALE AU PLAN NMTCL
         NORMALE(1) = 0
         NORMALE(2) = 0
         NORMALE(3) = 0
         NORMALE(NMTCL) = 1.0
         LPEXIST= 0
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
               KERR(1) = 'TRPLSE: VECTEUR NORMAL AU PLAN NUL'
               KERR(2) = '        VEUILLEZ RECOMMENCER'
            ELSE
               KERR(1) = 'TRPLSE: NULL NORMAL VECTOR'
               KERR(2) = '        GIVE IT AGAIN'
            ENDIF
            CALL LEREUR
            GOTO 100
         ENDIF
C
         DO I=1,3
            POINTS(I,2) = POINTS(I,2)/NORM
         ENDDO
C
C        MINIMUM ET MAXIMUM SONT MIS DANS LE TABLEAU HMMX
C        DES COORDONNEES DES POINTS PROJETES SUR LA DROITE DEFINIE
C        PAR UN POINT ET UN VECTEUR
         CALL MMPSPT( NBCOOR, NBPOI, POINTS, RMCN(MNPOGE+WYZPOI), HMMX )
         SECMIN = HMMX(1)
         SECMAX = HMMX(2)
         NOAXE  = 4
         NORMALE(1) = POINTS(1,2)
         NORMALE(2) = POINTS(2,2)
         NORMALE(3) = POINTS(3,2)
         LPEXIST = 0
         GOTO 145
C
      ELSE IF( NMTCL .GT. 15 ) THEN
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
C
      IF( NMTCL .EQ. 9 ) THEN
C
C        COULEUR des ARETES dans les PLANS
C        .................................
         CALL LIMTCL( 'couleur0' , I )
         IF( I .EQ. -1 ) GOTO 100
         IF( I .EQ. -2 ) THEN
C           COULEUR INVISIBLE A NE PAS TRACER
            NCOPLAN = -2
         ELSE IF( I .EQ. 0 ) THEN
C           LA COULEUR NOIRE
            NCOPLAN = 0
         ELSE
            NCOPLAN = N1COEL + I
         ENDIF
         GOTO 100
C
      ENDIF
C
      IF( NMTCL .EQ. 10 ) THEN
C
C        TYPE du TRAIT des ARETES dans les PLANS
C        .......................................
         CALL LIMTCL( 'typtrait' , I )
         IF( I .EQ. -1 ) GOTO 100
         NTLPLAN = I
         GOTO 100
C
      ENDIF
C
      GOTO( 110, 130, 140, 150, 160 ) , NMTCL-10
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
      IF( MNPLAV .LE. 0 ) GOTO 9999
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
         IF( MNPLAV .LE. 0 ) GOTO 9999
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
C     COEFFICIENT d'AMPLIFICATION POUR PROFIL
C     .......................................
 160  CALL INVITD( 'AMPLIFICATION' )
      NCVALS = 0
      CALL LIRRSP( NCVALS, R )
      IF( NCVALS .LE. 0 ) GOTO 100
      IF( R .EQ. 0 ) GOTO 160
      AMPLI = R
      GOTO 80
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
C     SEUILS MINIMUM MAXIMUM de la SOLUTION a TRACER
C     ..............................................
 195  CALL INVITE( 154 )
      CALL LIRRSP( NCVALS , R )
      IF( NCVALS .LE. 0 ) GOTO 100
      TMIN = R
      NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:10),'(G10.2)' ) TMIN
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'Seuil Minimum du TRACE= ' // KERR(MXLGER)(1:10)
      ELSE
         KERR(1) = 'Drawing Minimum TRHRESHOLD= ' // KERR(MXLGER)(1:10)
      ENDIF
      CALL LERESU
C
      CALL INVITE( 155 )
      CALL LIRRSP( NCVALS , R )
      IF( NCVALS .LE. 0 ) GOTO 100
      TMAX = R
C
      IF( TMIN .GT. TMAX ) THEN
         R    = TMIN
         TMIN = TMAX
         TMAX = R
      ENDIF
C
      IF( TMIN .LT. TMIN0 ) TMIN=TMIN0
      IF( TMAX .GT. TMAX0 ) TMAX=TMAX0
C
      NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:10),'(G10.2)' ) TMAX
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'Seuil Maximum du TRACE= ' // KERR(MXLGER)(1:10)
      ELSE
         KERR(1) = 'Drawing Maximum THRESHOLD= ' // KERR(MXLGER)(1:10)
      ENDIF
      CALL LERESU
      GOTO 100
C
C
C     1-ERE ETAPE: DEFINITION DE LA COORDONNEE DES PLANS DE SECTION
C     =============================================================
 200  IF( MNPLAV .EQ. 0 ) THEN
C        INITIALISATIONS PAR DEFAUT
         NBPLAN = ABS( NBPLAN )
         NBPLAN = MAX( 1, NBPLAN )
         MXPLAN = NBPLAN
         CALL TNMCDC( 'REEL', MXPLAN+1, MNPLAV )
         IF( MNPLAV .LE. 0 ) GOTO 9999
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
C     AFFICHAGES DES COORDONNEES DES PLANS
      IF( NOAXE .EQ. 4 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10310)
     %        (POINTS(I,1),I=1,3), (POINTS(I,2),I=1,3),
     %        (I,RMCN(MNPLAV-1+I),I=1,NBPLAN)
         ELSE
            WRITE(IMPRIM,20310)
     %        (POINTS(I,1),I=1,3), (POINTS(I,2),I=1,3),
     %        (I,RMCN(MNPLAV-1+I),I=1,NBPLAN)
         ENDIF
      ELSE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10210)
     %        NOAXE, (I,RMCN(MNPLAV-1+I),I=1,NBPLAN)
         ELSE
            WRITE(IMPRIM,20210)
     %        NOAXE, (I,RMCN(MNPLAV-1+I),I=1,NBPLAN)
         ENDIF
      ENDIF
C
10210 FORMAT(/'VALEURS DE LA COORDONNEE ',I1,' CONSTANTE DES PLANS'/
     % 5(I4,' : ',G13.5))
20210 FORMAT(/'VALUE of the COORDINATE ',I1,' CONSTANT of PLANES'/
     % 5(I4,' : ',G13.5))
10310 FORMAT(/'VALEURS DE LA COORDONNEE RELATIVES AU PLAN :'/
     %     ' POINT=(', 3(G13.5), ')'/
     %     ' VECTEUR NORMAL=(', 3(G13.5),')'/
     %     5(I3,' : ',G13.5))
20310 FORMAT(/'VALUE of the RELATIVE COORDINATE of PLANE :'/
     %     ' POINT=(', 3(G13.5), ')'/
     %     ' NORMAL VECTOR=(', 3(G13.5),')'/
     %     5(I3,' : ',G13.5))
C
C     2-EME ETAPE: CALCUL DES TABLEAUX DES INTERSECTIONS AVEC LES PLANS
C     NBFPLA : NOMBRE TOTAL DE FACES DES PLANS DE SECTION
C     XYZSFP : XYZ DES 4 SOMMETS DES FACES DES PLANS DE SECTION
C     TEMSFP : TEMPERATURE DES 4 SOMMETS AU PLUS DES FACES DES PLANS
C              POUR LES CAS NCAS0 A NCAS1
C     =================================================================
C     SI PLUSIEURS TEMPS, SECTION SELON UNE SEULE FAMILLE DE PLANS
C     PAS DE CUMUL POSSIBLE AVEC LES SECTIONS PRECEDENTES
      IF( NCAS0 .NE. NCAS1 ) NBFPLA=0
      IF( MODSEC .EQ. 0 ) THEN
C        STOCKAGE DE XYZ DES POINTS D'INTERSECTION AVEC LES PLANS
         NBCOPS = 3
      ELSE
C        STOCKAGE DE XYZ DES POINTS D'INTERSECTION AVEC LES PLANS
C        PLUS LA COORDONNEE CTE DU PLAN K POUR LE TRACE DES PROFILS
         NBCOPS = 4
      ENDIF
      NBFPLA = 0

C     TABLEAU DES COORDONNEES DES POINTS D'INTERSECTION DES EF AVEC LES PLANS
C     -----------------------------------------------------------------------
      IF( IERXYZSFP .EQ. 0 ) DEALLOCATE( XYZSFP )
ccc   MXFPLA  = NBFPLA + MAX( 4096, NBPOI ) * NBPLAN
ccc   MXFPLA  = NBFPLA + MIN( MAX( 4096, NBPOI ), 256 * 1024 )
      MXFPLA  = NBPLAN * MAX( 4096, NBPOI ) * 2
      MOXYZSFP = NBCOPS * 4 * MXFPLA
      IF( MOXYZSFP .LT. 0 ) THEN
         MOXYZSFP = 2**29
      ENDIF
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'trplse: DEMANDE  ALLOCATION XYZSFP(',MOXYZSFP,
     %          ') REELS'
         ALLOCATE ( XYZSFP( 1:MOXYZSFP ), STAT=IERXYZSFP )
         IF( IERXYZSFP .NE. 0 ) THEN
           PRINT*,'trplse: ERREUR ALLOCATION XYZSFP(',MOXYZSFP,') REELS'
            IERR = IERXYZSFP
            GOTO 9999
         ENDIF
         PRINT*,'trplse: CORRECTE ALLOCATION XYZSFP(',MOXYZSFP,') REELS'
      ELSE
         PRINT*,'trplse: ALLOCATION DEMAND  XYZSFP(',MOXYZSFP,') REALS'
         ALLOCATE ( XYZSFP( 1:MOXYZSFP ), STAT=IERXYZSFP )
         IF( IERXYZSFP .NE. 0 ) THEN
            PRINT*,'trplse: ALLOCATION ERROR XYZSFP(',MOXYZSFP,') REALS'
            IERR = IERXYZSFP
            GOTO 9999
         ENDIF
         PRINT*,'trplse: CORRECT ALLOCATION XYZSFP(',MOXYZSFP,') REALS'
      ENDIF

C     ALLOCATION du TABLEAU DE LA VALEUR DE LA SOLUTION AUX POINTS
C     DES FACES D'INTERSECTION DES EF AVEC LES PLANS
C     ------------------------------------------------------------
      IF( IERTEMSFP .EQ. 0 ) DEALLOCATE( TEMSFP )
      MOTSFA = 4 * ( NCAS1 - NCAS0 + 1 ) * MXFPLA
      IF( MOTSFA .LT. 0 ) THEN
         MOTSFA = 2**29
      ENDIF
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'trplse: DEMANDE  ALLOCATION TEMSFP(',MOTSFA,') REELS'
         ALLOCATE ( TEMSFP( 1:MOTSFA ), STAT=IERTEMSFP )
         IF( IERTEMSFP .NE. 0 ) GOTO 222
         PRINT*,'trplse: CORRECTE ALLOCATION TEMSFP(',MOTSFA,') REELS'
      ELSE
         PRINT*,'trplse: ALLOCATION DEMAND  TEMSFP(',MOTSFA,') REALS'
         ALLOCATE ( TEMSFP( 1:MOTSFA ), STAT=IERTEMSFP )
         IF( IERTEMSFP .NE. 0 ) GOTO 222
         PRINT*,'trplse: CORRECT ALLOCATION TEMSFP(',MOTSFA,') REALS'
      ENDIF
      GOTO 225

C     TABLEAU TROP GRAND => PAS ASSEZ de MEMOIRE
 222  IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'trplse: ERREUR ALLOCATION TEMSFP(',MOTSFA,') REELS'
         IERR = IERTEMSFP
         GOTO 9999
      ELSE
         PRINT*,'trplse: ALLOCATION ERROR TEMSFP(',MOTSFA,') REALS'
         IERR = IERTEMSFP
         GOTO 9999
      ENDIF


C     MISE A L'ECHELLE SEULEMENT ACTIF POUR LE TRACE DES PROFILS
 225  IF( TMIN .NE. TMAX ) THEN
         RAPPORT = AMPLI / (TMAX-TMIN) * DIAGON * 0.333333
      ELSE
         RAPPORT = 1.0
      ENDIF
      VNORMAL(1) = NORMALE(1) * RAPPORT
      VNORMAL(2) = NORMALE(2) * RAPPORT
      VNORMAL(3) = NORMALE(3) * RAPPORT

C     CONSTRUCTION DU TABLEAU DES FACES D'INTERSECTION AVEC LES EF
C     BOUCLE SUR LES DIFFERENTS TYPES D'ELEMENTS FINIS DU MAILLAGE
C     ------------------------------------------------------------
      DO 250 I = 0, NBTYEL-1
C
C        L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF DE CE TYPE D'EF
         MNELE = MCN( MNELEM + I )
C
C        LE NUMERO DU TYPE DES ELEMENTS FINIS
         NUTYEL = MCN( MNELE + WUTYEL )
         IF( NUTYEL .EQ. 30 ) GOTO 250
C        PAS DE TRACE DES 6CUBES, SEULEMENT DES 3Q1C EXTRAITS
C
C        LE NOMBRE DE TELS ELEMENTS
         NBELEM = MCN(MNELE + WBELEM )
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
C        CALCUL DES FACES APPARTENANT AUX PLANS SECTION DE L'OBJET
         CALL TRPLS3( MODSEC, NOAXE,  POINTS, NBPLAN, RMCN(MNPLAV),
     %                NBCOOR, NCAS0,  NCAS1,  NTDL,
     %                NTYP,   TEMPER, dptemp,
     %                NUTYEL, NBELEM,
     %                NBNOE,  MCN(MNELE+WUNDEL), NBPOE, MCN(MNPGEL),
     %                NBPOI,  MCN(MNPOGE+WYZPOI),
     %                MXSOMM, MCN(MNSOLE), MCN(MNCOPO),
     %                MXPILE, MCN(MNPILE),
     %                MCN(MNVALS), MCN(MNFBAS),
     %                NBCOPS, TRIANG, MXFPLA, NBFPLA,
     %                XYZSFP, TEMSFP )

 250  CONTINUE
C     LES FACES SECTIONS AVEC LES PLANS ONT ETE CALCULEES
      LPEXIST = 1
C
      IF( NMTCL .EQ. 90 ) GOTO 300
      GOTO 100
C
C     3-EME ETAPE: RECUEIL DES OPTIONS DE LA VISEE POUR VOIR
C     LES FACES FRONTALIERES DE L'OBJET ET LES PLANS DE SECTION
C     =========================================================
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
C     TABLEAU POUR LE TRI SELON LA DISTANCE A L'OEIL
      IF( MNBARY .GT. 0 ) CALL TNMCDS( 'REEL',   NBF, MNBARY )
      IF( MNNUFA .GT. 0 ) CALL TNMCDS( 'ENTIER', NBF, MNNUFA )
      MNBARY = 0
      MNNUFA = 0
      NBF    = NBARFR + NBFPLA * ( 1 + MODSEC )
      CALL TNMCDC( 'REEL',   NBF, MNBARY )
      IF( MNBARY .LE. 0 ) GOTO 9999
      CALL TNMCDC( 'ENTIER', NBF, MNNUFA )
      IF( MNNUFA .LE. 0 ) GOTO 9999
C
C     INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
C     -------------------------------------------
      IF( LORBITE .NE. 0 ) THEN
         CALL ORBITE0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 300
      ENDIF
C
C     4-EME ETAPE: FORMATION DES TABLEAUX NUMERO DES FACES ET
C                  DISTANCE AXONOMETRIQUE A L'OEIL
C                  POUR LE BARYCENTRE DES FACES DES PLANS DE SECTION
C     ==============================================================
 400  CALL TRBFIS( TRIANG, NBFPLA, NBCOPS, NOAXE, XYZSFP,
     %             MCN(MNNUFA), RMCN(MNBARY) )
C
C     FORMATION DES TABLEAUX NUMERO DES FACES ET DISTANCE AXONOMETRIQUE
C     POUR LES ARETES FRONTALIERES DE L'OBJET
C     ATTENTION LES ARETES FRONTALIERES DOIVENT OBLIGATOIREMENT ETRE APRES
C               LES FACES DES PLANS DE SECTION
      CALL TRBARF( MOARFR, MXARFR, L1ARFR, MCN(MNAFOB+WAREFR),
     %             MCN(MNPOGE+WBCOOP), RMCN(MNPOGE+WYZPOI),
     %             NBFPLA*(1+MODSEC), MCN(MNNUFA), RMCN(MNBARY))
C
C     LE TRI PAR TAS DE CETTE DISTANCE
C     LA FACE LA PLUS PROCHE EST LA PREMIERE
      CALL TRITRP( NBF, RMCN(MNBARY), MCN(MNNUFA) )
C
C     5-EME ETAPE: TRACE EFFECTIF DES ARETES FRONTALIERES ET
C                  DES FACES D'INTERSECTION AVEC LES PLANS
C     ======================================================
C
C     *****************************************************************
C     BOUCLE SUR LES VECTEURS SOLUTIONS A TRACER
C     *****************************************************************
      DO 1000 NCAS = NCAS0, NCAS1
C
C        TEMPS DE CALCUL DU VECTEUR NCAS
         TEMPS = TIMES( NCAS )
C
C        L'ECRAN PIXELS EST EFFACE (PAS DE SCINTILLEMENT)
         CALL EFFACEMEMPX

C        TRACE DES AXES 3D
         NETAXE = 0
         CALL TRAXE3
C
         IF( MODSEC .EQ. 0 ) THEN
C
C           TRACE DANS DES PLANS DE SECTION
            CALL T3AFP3( NBARFR, NBFPLA, MCN(MNNUFA),
     %                   MOARFR, MXARFR, MCN(MNAFOB+WAREFR),
     %                   RMCN(MNPOGE+WYZPOI),
     %                   TRIANG, NBCOPS, XYZSFP, TEMSFP,
     %                   NCAS0,  NCAS1,  NCAS,
     %                   TMIN,   TMAX )
         ELSE
C
C           TRACE ORTHOGONALEMENT A DES PLANS DE SECTION
            CALL T3AFP4( NBARFR, NBFPLA, MCN(MNNUFA),
     %                   MOARFR, MXARFR, MCN(MNAFOB+WAREFR),
     %                   NBCOOR, RMCN(MNPOGE+WYZPOI),
     %                   TRIANG, NOAXE,  XYZSFP, TEMSFP,
     %                   NCAS0,  NCAS1,  NCAS,
     %                   TMIN,   TMAX,   VNORMAL, NCOPLAN, NTLPLAN )
         ENDIF
C
C        TRACE DU POINT DE TEMPERATURE MINIMALE ET MAXIMALE
         CALL XVTYPETRAIT( LIGCON )
         IF( NCAS .EQ. NCAMIN ) THEN
            MN = MNPOGE + WYZPOI + NBCOOR*NOEMIN -NBCOOR
            CALL SYMBOLE3D( NCNOIR, RMCN(MN), '.Min' )
         ENDIF
         IF( NCAS .EQ. NCAMAX ) THEN
            MN = MNPOGE + WYZPOI + NBCOOR*NOEMAX -NBCOOR
            CALL SYMBOLE3D( NCNOIR, RMCN(MN), '.Max' )
         ENDIF
cccc
cccc  modifs pour trace TAMU sur les atomes
c-------------------------------------------------------------------
cccc     Trace TAMU du repere et de H2+
ccc      call trah2p
c-------------------------------------------------------------------
cccc     l'axe x entre les groupes de protons
ccc      call xvepaisseur( 2 )
ccc      xyz1(1) =-0.5
ccc      xyz1(2) =  0
ccc      xyz1(3) =  0
cccc
ccc      xyz2(1) = 0.5
ccc      xyz2(2) =  0
ccc      xyz2(3) =  0
ccc      call trait3d( ncblan, xyz1, xyz2 )
ccc      call traxe3
cccc
cccc     2 protons
ccc      call symbole3d( ncroug, xyz1, '++' )
cccc
cccc     1 proton
ccc      call symbole3d( ncroug, xyz2, '+' )
C
C        EFFACEMENT DE LA LEGENDE SUR POSTSCRIPT
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
                     KERR(1) = 'TRPLSE: MAUVAISE VALEUR DE LASOPS'
                     KERR(2) = '        ARRET DU POSTSCRIPT'
                  ELSE
                     KERR(1) = 'TRPLSE: BAD VALUE of LASOPS'
                     KERR(2) = '        STOP of POSTSCRIPT'
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
C        6-EME ETAPE: LE TRACE DU TITRE FINAL
C        ====================================
C        LE TRACE DE LA LEGENDE : COULEURS => VALEURS
         CALL  LEGCOULSO( TMIN, TMAX )
C
C        RETOUR AU TRACE NORMAL POUR POSTSCRIPT
         IF( LASOPS .NE. 0 ) THEN
            LASOPS = LASOPS - 10
            CALL XVPOSTSCRIPT(LASOPS)
         ENDIF
C
C        RETOUR AUX PARAMETRES INITIAUX
         CALL XVEPAISSEUR( 1 )
         CALL XVTYPETRAIT( LIGCON )
C
cccC        DEFINITION DU TITRE ET FIN DU TRACE
ccc         CALL LETITR( NOPROJ, MODECO, NCAS, TEMPS, KNOM )
C
C        TRACE DE LA 2-EME LIGNE DU TITRE DU TRACE
         CALL TIT2LG( KNOMOB, MODECO )
C
C        DEFINITION DE LA 3-EME LIGNE DU TITRE
         CALL TIT3LG( MODECO, NCAS, TEMPS, TMIN, TMAX, KNOM )
C
C        TRACE DU TITRE
         CALL TRFINS( KNOM )
C
C        MISE SUR FICHIER NomfgifBoImage.xwd puis NomfgifNoImage.jpg
C        DE LA PIXMAP de la FENETRE X11 ACTUELLE
         CALL VIDEO1( NOMFGIF, NCAS )
C
C        ATTENDRE POUR LIRE LE TRACE
         CALL ATTENDSEC( TEMP2TRAC )
C
C        FIN DE LA BOUCLE SUR LES CAS
 1000 CONTINUE
C
C     CONSTRUIRE le FICHIER VIDEO Nomfic.gif A PARTIR DES FICHIERS
C     CONSTRUITS de NOMS NomfgifNoImag.jpg
C     ------------------------------------------------------------
      CALL VIDEOFIN( NOMFGIF )
C
C     RETOUR POUR UNE NOUVELLE VISEE
C     ------------------------------
      IF( LORBITE .NE. 0 ) THEN
         IF( NCAS0 .EQ. NCAS1 ) THEN
C           ORBITE BOUTON ENFONCE et DEPLACE
            CALL ORBITE1( NOTYEV )
         ELSE
C           ORBITE BOUTON ENFONCE et DEPLACE et RELACHE
            CALL ORBITE3( NOTYEV )
         ENDIF
         IF( NOTYEV .EQ. 0 ) GOTO 300
         GOTO 400
      ELSE
C        POUR LIRE LE TRACE AVANT D'AFFICHER UN MENU
         CALL CLICSO
      ENDIF
      GOTO 300
C
C     7-EME ETAPE: SORTIE DU TRACE DES PLANS DE SECTION
C     =================================================
 9000 IF( IERXYZSFP .EQ. 0 ) THEN
         PRINT*,' DEALLOCATE( XYZSFP )'
         DEALLOCATE( XYZSFP )
         IERXYZSFP = 1
      ENDIF
      IF( IERTEMSFP .EQ. 0 ) THEN
         PRINT*,'DEALLOCATE( TEMSFP )'
         DEALLOCATE( TEMSFP )
         IERTEMSFP = 1
      ENDIF
      IF( MNPLAV .GT. 0 ) CALL TNMCDS( 'REEL',  MXPLAN+1, MNPLAV )
      IF( MNBARY .GT. 0 ) CALL TNMCDS( 'REEL',   NBF,     MNBARY )
      IF( MNNUFA .GT. 0 ) CALL TNMCDS( 'ENTIER', NBF,     MNNUFA )
      IF( MNSOLE .GT. 0 ) CALL TNMCDS( 'REEL2' , MOSOLE,  MNSOLE )
      IF( MNPILE .GT. 0 ) CALL TNMCDS( 'ENTIER', MOPILE,  MNPILE )
      IF( MNVALS .GT. 0 ) CALL TNMCDS( 'REEL'  , MOVALS,  MNVALS )
      RETURN
C
C     PAS ASSEZ DE MOTS DANS MCN
C     ==========================
 9999 NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) ='PAS ASSEZ de MOTS pour le SUPER-TABLEAU MCN'
         KERR(2) ='AUGMENTER ce NOMBRE dans incl.pp.inc'
      ELSE
         KERR(1) ='NO SUFFICIENT WORD NUMBER of SUPER-ARRAY MCN'
         KERR(2) ='AUGMENT IT in incl/pp.inc'
      ENDIF
      CALL LEREUR
      GOTO 9000

      END
