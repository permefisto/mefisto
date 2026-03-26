      SUBROUTINE TRZTXY( NOFOUT, NDIM,   KNOMOB, MODECO,
     %                   NBTYEL, MNELEM, MNPOGE,
     %                   NCAS0,  NCAS1,  NTDL,  NTYP, TEMPER, dptemp,
     %                   TMIN0,NOEMIN0,NCAMIN0, TMAX0,NOEMAX0,NCAMAX0,
     %                   TIMES )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE 3D DE LA SURFACE MAILLEE DE L'OBJET 2D
C -----    EN CONSTRUISANT LES TMS XYZSOMMET ET NSEF DES EF SURFACIQUES
C          ET EN PORTANT EN Z LA TEMPERATURE CALCULEE OU L'ERREUR
C
C ENTREES:
C --------
C NOFOUT : NUMERO DE LA FONCTION TEMPERATURE_EXACTE(t,x,y,z)
C          ou DEPLACEMENT_EXACT(t,x,y,z,nc)
C          ou PARTIE_REELLE_EXACTE(t,x,y,z,nc)
C          ou PARTIE_IMAGINAIRE_EXACTE(t,x,y,z,nc)
C          ou VITESSE_EXACTE(t,x,y,z,nocomp) ou EXACT_VELOCITY(t,x,y,z,nocomp)
C          ou PRESSION_EXACTE(t,x,y,z) ou EXACT_PRESSURE(t,x,y,z)
C          DANS LE LEXIQUE DES FONCTIONS DE L'UTILISATEUR
C NDIM   : DIMENSION DE L'ESPACE DE L'OBJET (2 OU 3)
C KNOMOB : NOM DE L'OBJET
C MODECO : MODE DE TRACE DES VECTEURS
C         =1  CE SONT DES TEMPERATURES
C         =2  CE SONT DES VECTEURS PROPRES
C         =3  CE SONT DES PRESSIONS P1 A PARTIR D'UNE INTERPOLATION
C             SOIT P2 SOIT P1+BULLE P3
C         =4  CE SONT DES ERREURS PONCTUELLES SOLEX-SOLCAL
C             => TRACE DE L'ERREUR ENTRE TEMPERATURE_EXACTE
C                ET LA TEMPERATURE CALCULEE
C         =6  CE SONT DES NORMES DE VITESSE AUX NOEUDS D'UN MAILLAGE
C         =8  CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C             DU MODULE D'UNE ONDE COMPLEXE
C            (CALCUL DU MODULE DE L'ERREUR COMPLEXE A FAIRE ICI)
C         =9  PARTIE REELLE     D'UNE ONDE COMPLEXE NLSE
C         =10 PARTIE IMAGINAIRE D'UNE ONDE COMPLEXE NLSE
C         =11 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C             DU MODULE DE LA VITESSE D'UN FLUIDE
C            (CALCUL DE L'ERREUR(VITESSE EXACTE - CALCULEE) DEJA CALCULE)
C         =12 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C             DE LA PRESSION P1 DANS UN FLUIDE
C            (CALCUL DE L'ERREUR(PRESSION EXACTE-CALCULEE) DEJA CALCULE)
C
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS FINIS DANS CETTE TOPOLOGIE
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE
C MNELEM : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE
C MNPOGE : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C MNNOEU : ADRESSE MCN DU TABLEAU NOEUDS D'INTERPOLATION DU MAILLAGE
C NDPGST : CODE TRAITEMENT DES SOMMETS POINTS NOEUDS DU MAILLAGE
C          0 : NOEUDS=POINTS=SOMMETS
C          1 : NOEUDS=POINTS#SOMMETS
C          2 : NOEUDS#POINTS=SOMMETS
C          3 : NOEUDS#POINTS#SOMMETS
C NCAS0  : NUMERO DU PREMIER CAS A TRAITER
C NCAS1  : NUMERO DU DERNIER CAS A TRAITER
C NCAS0:NCAS1 : NCAS1-NCAS0+1 NOMBRE DE CAS OU VECTEURS TEMPER
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE
C NTYP   : =0 EMPLOI de TEMPER
C         =1  EMPLOI de dptemp
C TEMPER : TABLEAU DES NTDL * NCAS0:NCAS1 TEMPERATURES
C dptemp : TABLEAU DES NCAS0:NCAS1 TABLEAU(NTDL) TEMPERATURES
C TMIN0  : TEMPERATURE MINIMALE
C NOEMIN0: NUMERO DU NOEUD OU LA TEMPERATURE EST MINIMALE
C NCAMIN0: NO DU CAS DE 1 A NCAS0:NCAS1 OU LA TEMPERATURE EST MINIMALE
C TMAX0  : TEMPERATURE MAXIMALE DES NCAS0:NCAS1 CAS
C NOEMAX0: NUMERO DU NOEUD OU LA TEMPERATURE EST MAXIMALE
C NCAMAX0: NO DU CAS DE 1 A NCAS0:NCAS1 OU LA TEMPERATURE EST MAXIMALE
C TIMES  : TEMPS DU CALCUL DES NCAS0:NCAS1 VECTEURS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1998
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray              Mars 2021
C23456---------------------------------------------------------------012
      IMPLICIT    INTEGER(W)
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/inteel.inc"
      include"./incl/ctemps.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
      include"./incl/xvfontes.inc"
      include"./incl/traaxe.inc"
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     KNOMOB
      CHARACTER*120     KNOM
      CHARACTER*4       NOMELE(2)
      CHARACTER*24      NOMFGIF
      INTEGER           NOEUDS(9)
      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(NCAS0:NCAS1) :: dptemp
      DOUBLE PRECISION  TEMPER( NTDL, NCAS0:NCAS1 ),
     %                  PBASE(9),T
      DOUBLE PRECISION  DEC, D, DTX, IDTX, DPARAF(5)
      REAL              TIMES(NCAS0:NCAS1)
      INTEGER           NOFOWE(2)
C
      IERR   = 0
      MNXYZS = 0
      MNNSEF = 0
      MNZTXY = 0
      NBSOM  = 0
      NBEFOB = 0
C
C     SAUVEGARDE DE XYZ MIN MAX ACTUELS
      XMIN0 = COOEXT(1,1)
      XMAX0 = COOEXT(1,2)
      YMIN0 = COOEXT(2,1)
      YMAX0 = COOEXT(2,2)
      ZMIN0 = COOEXT(3,1)
      ZMAX0 = COOEXT(3,2)
C
C     PROTECTION DES 4 PARAMETRES
      NOEMIN = NOEMIN0
      NCAMIN = NCAMIN0
      NOEMAX = NOEMAX0
      NCAMAX = NCAMAX0
      TMIN   = TMIN0
      TMAX   = TMAX0
C
C     NOM DU FICHIER VIDEO SUIVI DU NO DE L'IMAGE AVEC 4 CHIFFRES
C     SELON LA VALEUR DE MODECO
C          =1 CE SONT DES TEMPERATURES
C          =2 CE SONT DES VECTEURS PROPRES
C          =3 CE SONT DES PRESSIONS P1 A PARTIR D'UNE INTERPOLATION
C             SOIT P2 SOIT P1+BULLE P3
C          =4 CE SONT DES ERREURS PONCTUELLES SOLEX-SOLCAL
C             => TRACE DE L'ERREUR ENTRE TEMPERATURE_EXACTE
C                ET LA TEMPERATURE CALCULEE
C          =6 CE SONT DES NORMES DE VITESSE AUX NOEUDS D'UN MAILLAGE
C          =8 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C             DU MODULE D'UNE ONDE COMPLEXE
C            (CALCUL DU MODULE DE L'ERREUR COMPLEXE A FAIRE ICI)
C          =9 CE SONT DES PARTIES REELLES     D'UNE ONDE COMPLEXE
C          =10CE SONT DES PARTIES IMAGINAIRES D'UNE ONDE COMPLEXE
C          =11 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C             DE LA NORME DE LA VITESSE D'UN FLUIDE 2D
C            (CALCUL DE L'ERREUR(VITESSE EXACTE - CALCULEE) A FAIRE ICI)
C          =12 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C             DE LA PRESSION P1 DANS UN FLUIDE
      CALL VIDEONM( MODECO, 'zfxy', NOMFGIF )
C
C     NDIM LA DIMENSION 1 OU 2 OU 3 DE L'ESPACE DES COORDONNEES
C     ATTENTION: CETTE PROGRAMMATION SOUS ENTEND NOEUDS=POINTS ET 2D
      NBPOI = MCN(MNPOGE+WNBPOI)
      CALL DIMCOO( NBPOI, MCN(MNPOGE+WYZPOI), NDIM )
      IF( NDIM .NE. 2 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNOMOB
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) ='ERREUR: OBJET NON EN 2D'
         ELSE
            KERR(2) ='ERROR: NOT a 2D-OBJECT'
         ENDIF
         CALL LEREUR
         IERR = 4
         GOTO 9999
      ENDIF
C     LE NOMBRE DE SOMMETS DES EF DECOUPES EVENTUELLEMENT EN 4
      NBSOM = NBPOI
C
C     CALCUL DU NOMBRE DE TRIANGLES ET QUADRANGLES
C     ============================================
      LDEGRE = 1
      NBTRI  = 0
      NBQUA  = 0
      DO 10 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"xxxx
         MNELE = MCN( MNELEM - 1 + NOTYEL )
C
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        NCOGEL LE CODE GEOMETRIQUE DE L'ELEMENT FINI EST CALCULE
         CALL ELNUCG( NUTYEL, NCOGEL )
         IF( NCOGEL .NE. 3 .AND. NCOGEL .NE. 4 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = KNOMOB
            IF( LANGAG .EQ. 0 ) THEN
               KERR(2) ='TRZTXY: EF NON TRAITE EN DIMENSION 2'
            ELSE
               KERR(2) ='TRZTXY: FE NOT TREATED. ONLY 2D'
            ENDIF
            CALL LEREUR
            IERR = 3
            GOTO 9999
         ENDIF
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        LE NOMBRE D'EF ET DE SOMMETS DU MAILLAGE 2D
         IF( NCOGEL .EQ. 3 ) THEN
            NBTRI = NBTRI + NBELEM
         ELSE
            NBQUA = NBQUA + NBELEM
         ENDIF
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL, NOMELE )
         CALL ELTYCA( NUTYEL )
C
C        LAGRANGE DE DEGRE 1 OU 2?
         IF( NBNOE .GT. 4 ) LDEGRE = 2
C
 10   CONTINUE
C
C     CAS SPECIAL DU BARYCENTRE DES EF DE BREZZI-FORTIN
C     NBPOI = NB SOMMETS + NB BARYCENTRES
      NBEF  = NBTRI + NBQUA
      IF( NUTYEL .EQ. 13 .AND. NBPOI .GT. NBEF ) THEN
C        TRIANGLE BREZZI-FORTIN AVEC VITESSE AU BARYCENTRE
         NBPOI = NBPOI - NBEF
      ENDIF
C
C     CREATION DU TMC XYZSOMMET DE L'OBJET
C     ====================================
      IF( LDEGRE .EQ. 2 .AND. NBQUA .GT. 0 ) THEN
C        LE BARYCENTRE DE CHACUN DES QUAD 2Q2C EST AJOUTE
         NBSOM = NBPOI + NBQUA
      ENDIF
C     LES SOMMETS DU MAILLAGE DE DEGRE 1
      CALL TNMCDC( 'MOTS', WYZSOM+3*NBSOM, MNXYZS )
C
C     CREATION DU TMC Z(T,X,Y) DE L'OBJET
C     ===================================
C     LA TEMPERATURE EN CES SOMMETS POUR LES NCAS0:NCAS1 CAS
      NBVECT = NCAS1-NCAS0+1
      MOZTXY = NBSOM * NBVECT
      CALL TNMCDC( 'REEL', MOZTXY, MNZTXY )
      MNZT = MNZTXY - 1
C
C     LE NOMBRE DE SOMMETS
      MCN( MNXYZS + WNBSOM ) = NBSOM
C
C     LE NOMBRE DE TANGENTES
      MCN( MNXYZS + WNBTGS ) = 0
C
C     ADRESSE MCN DES XYZ DES SOMMETS
      MNS = MNXYZS + WYZSOM
C     ADRESSE MCN DES XYZ DES NOEUDS=POINTS
      MNN = MNPOGE + WYZPOI
C
      XMIN = RINFO('GRAND')
      XMAX =-XMIN
      YMIN = XMIN
      YMAX = XMAX
      IF( MODECO .EQ.  4 .OR. MODECO .EQ. 8 ) THEN
C        POUR CALCULER LE MIN MAX DES ERREURS CALCULEES ICI
         TMIN =  1E28
         TMAX = -1E28
      ENDIF
C
C     CALCUL DU NOMBRE DE PARAMETRES DE LA FONCTION EXACTE
C     TEMPERATURE_EXACTE(t,x,y,z) ou DEPLACEMENT_EXACT(t,x,y,z,nc)
C     ou PARTIE_PARTIE_REELLE_EXACTE(t,x,y,z)
C     ou PARTIE_IMAGINAIRE_EXACTE(t,x,y,z)
C     ou VITESSE_EXACTE(t,x,y,z,nocomp) ou EXACT_VELOCITY(t,x,y,z,nocomp)
C     ou PRESSION_EXACTE(t,x,y,z) ou EXACT_PRESSURE(t,x,y,z)
      NOFOWE(1) = NOFOPREX()
      NOFOWE(2) = NOFOPIEX()
C
      IF( NOFOTEEX() .EQ. NOFOUT ) THEN
C        TEMPERATURE_EXACTE(t,x,y,z)
         NBARGS = 4
      ELSE IF( NOFODEEX() .EQ. NOFOUT ) THEN
C        DEPLACEMENT_EXACT(t,x,y,z,nc)
         NBARGS = 5
      ELSE IF( NOFOWE(1) .EQ. NOFOUT ) THEN
C        PARTIE_PARTIE_REELLE_EXACTE(t,x,y,z)
         NBARGS = 4
      ELSE IF( NOFOWE(2) .EQ. NOFOUT ) THEN
C        PARTIE_IMAGINAIRE_EXACTE(t,x,y,z)
         NBARGS = 4
      ELSE IF( NOFOWE(2) .EQ. NOFOUT ) THEN
C        PARTIE_IMAGINAIRE_EXACTE(t,x,y,z)
         NBARGS = 4
      ELSE IF( NOFOVITE() .EQ. NOFOUT ) THEN
C        VITESSE_EXACTE(t,x,y,z,nocomp)
         NBARGS = 5
      ELSE IF( NOFOPRES() .EQ. NOFOUT ) THEN
C        PRESSION_EXACTE(t,x,y,z)
         NBARGS = 4
      ELSE
C        ERREUR PAS DE TRACE
         IERR = 1
         GOTO 9999
      ENDIF
C
C     TRAITEMENT DE LA TEMPERATURE
      DO 70 NS=1,NBPOI
C
C        L'ABSCISSE
         X = RMCN( MNN )
         IF( X .LT. XMIN ) XMIN = X
         IF( X .GT. XMAX ) XMAX = X
         RMCN( MNS     ) = X
C
C        L'ORDONNEE
         X = RMCN( MNN + 1 )
         IF( X .LT. YMIN ) YMIN = X
         IF( X .GT. YMAX ) YMAX = X
         RMCN( MNS + 1 ) = X
C
C        LA COTE
         RMCN( MNS + 2 ) = 0.0
C
         IF( MODECO .EQ. 4 .OR. MODECO .EQ. 8 ) THEN
C
C           CE SONT DES ERREURS A CALCULER ICI
            DO NCAS = NCAS0, NCAS1
C
C              CALCUL DE LA TEMPERATURE_EXACTE EN CE NOEUD
               TEMPS = TIMES( NCAS )
               DPARAF(1) = TEMPS
               DPARAF(2) = RMCN(MNS)
               DPARAF(3) = RMCN(MNS+1)
C              VALEUR DE Z
               DPARAF(4) = 0D0
C              NUMERO DE LA SEULE COMPOSANTE 1 ICI
               DPARAF(5) = 1D0
C              TEMPERATURE_EXACTE(TEMPS,X,Y,Z) ou
C              DEPLACEMENT_EXACT (TEMPS,X,Y,Z,NOCOMP)
C              PRESSION_EXACTE(t,x,y,z)
               IF( MODECO .EQ. 4 ) THEN
C                 1 SEULE FONCTION UTILISATEUR A CALCULER
                  CALL FONVAL( NOFOUT, NBARGS, DPARAF, NCODEV, DTX )
               ELSE IF( MODECO .EQ. 8 ) THEN
C                 2 FONCTIONS UTILISATEUR A CALCULER ET MODULE A CALCULER
C                 FONCTION PARTIE_REELLE_EXACTE(t,x,y,z)
                  CALL FONVAL( NOFOWE(1), NBARGS, DPARAF, NCODEV,  DTX )
C                 FONCTION PARTIE_IMAGINAIRE_EXACTE(t,x,y,z)
                  CALL FONVAL( NOFOWE(2), NBARGS, DPARAF, NCODEV, IDTX )
C                 LE MODULE
                  DTX = SQRT( DTX**2 + IDTX**2 )
               ENDIF
C
               IF( NCODEV .GT. 0 ) THEN
C                 L'ERREUR NON ABSOLUE
                  IF( NTYP .EQ. 0 ) THEN
                     T = TEMPER( NS, NCAS )
                  ELSE
                     T = dptemp( NCAS )%dptab( NS )
                  ENDIF
                  D = DTX - T
                  IF( D .LT. TMIN ) THEN
                     NOEMIN = NS
                     TMIN   = REAL( D )
                     NCAMIN = NCAS
                  ENDIF
                  IF( D .GT. TMAX ) THEN
                     NOEMAX = NS
                     TMAX   = REAL( D )
                     NCAMAX = NCAS
                  ENDIF
               ENDIF
C
C              LA SOLUTION EXACTE
               MN = MNZT + ( NCAS - NCAS0 ) * NBSOM
               RMCN( MN + NS ) = REAL( DTX )
C
            ENDDO
         ENDIF
C
         MNN = MNN + 3
         MNS = MNS + 3
 70   CONTINUE
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNXYZS) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNXYZS + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     MISE A L'ECHELLE DE Z POUR AVOIR LE MAX DE L'ECART EN X ET Y
      X = TMAX - TMIN
cccc
cccc      AMPLITUDE EN Z DU TRACE => MEME ECHELLE POUR LES 3 TRACES
cccc      X=1.23 pour ther2abc 3.30 pour ther3abc 15.0 pour ther1abc
cccc          POUR LE CAS 2d du pb -delta u -u**3 de  CHEN-ZHOU
cccc
      IF( X .EQ. 0 ) THEN
C        TEMPERATURE CONSTANTE
         X = 1.0
      ENDIF
      DEC = MAX( XMAX-XMIN, YMAX-YMIN ) / 1.85
C     Z EST RAMENE DANS L'INTERVALLE 1/2 LONGUEUR
      ECHELLE = REAL( DEC / X )
C
C     COTE DU MAILLAGE PAR RAPPORT A LA SURFACE Z=TEMPERATURE(X,Y)
      DECALZ  = REAL( DEC * (-0.15D0) )
C
C     CREATION DU TMS NSEF DE L'OBJET
C     ===============================
C     LE NOMBRE D'EF
      NBEFOB = ( NBTRI + NBQUA ) * LDEGRE**2
      CALL TNMCDC( 'MOTS', WUSOEF+4*NBEFOB, MNNSEF )
C
C     LE TYPE DE L'OBJET  ICI = 3 SURFACE!
      MCN( MNNSEF + WUTYOB ) = 3
C     SURFACE NON FERMEE
      MCN( MNNSEF + WUTFMA ) = 0
C     4 SOMMETS PAR EF
      MCN( MNNSEF + WBSOEF ) = 4
C     0 TANGENTE PAR EF
      MCN( MNNSEF + WBTGEF ) = 0
C     NOMBRE D'EF
      MCN( MNNSEF + WBEFOB ) = NBEFOB
C     NOMBRE D'EF A TGS
      MCN( MNNSEF + WBEFTG ) = 0
C     NOMBRE D'EF AVEC POINTEUR
      MCN( MNNSEF + WBEFAP ) = 0
C     MAILLAGE NON STRUCTURE
      MCN( MNNSEF + WUTYMA ) = 0
C
C     ADRESSE MCN DU NUMERO DES 4 SOMMETS DES EF
      MNN = MNNSEF + WUSOEF - 5
      NBP = NBPOI
C
C     ADRESSE MCN DES XYZ DES SOMMETS
      MNS = MNXYZS + WYZSOM - 4
C
      DO 99 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"xxxx
         MNELE = MCN( MNELEM - 1 + NOTYEL )
C
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        NCOGEL LE CODE GEOMETRIQUE DE L'ELEMENT FINI EST CALCULE
         CALL ELNUCG( NUTYEL, NCOGEL )
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL, NOMELE )
         CALL ELTYCA( NUTYEL )
C
         IF( LDEGRE .EQ. 2 .AND. NBNOE .EQ. 8 ) THEN
C           RETROUVER LE NUMERO DES INTERPOLATIONS DE LA TEMPERATURE
C           ET DES COMPOSANTES DE LA TRANSFORMATION:ELT REFERENCE->ELEMENT
            CALL ELINTE( 'THERMIQUE', NUTYEL, NDIMF, NOINTF,
     &                   NBINVA, NUINVA, NUINTI, NBNDIN )
C           LA TEMPERATURE EST ICI LA SEULE INCONNUE VARIATIONNELLE
            NOINTE = NUINTI(1)
C           VALEUR DES POLYNOMES DE BASE AU BARYCENTRE DU CARRE UNITE
            CALL INTERP( NOINTE, 0.5D0, 0.5D0, 0D0, K, PBASE )
         ENDIF
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        LA BOUCLE SUR LES EF DE CE TYPE NUTYEL
         DO 95 NUELEM = 1, NBELEM
C
C           LES NOEUDS DE L'ELEMENT FINI
            CALL EFNOEU( MNELE, NUELEM, NBNDEL, NOEUDS )
C
            IF( LDEGRE .EQ. 1 ) THEN
C
C              1 EF OBJET => 1 EF SURFACE
               MNN = MNN + 4
               DO I=1,NCOGEL
C                 LE NUMERO DE SOMMET DU NOEUD I
                  MCN( MNN + I ) = NOEUDS(I)
               ENDDO
C              LE 4-EME SOMMET D'UN TRIANGLE EST MIS A ZERO
               IF( NCOGEL .EQ. 3 )  MCN( MNN+4 ) = 0
C
            ELSE
C
C              1 EF OBJET => 4 EF SURFACE
               IF( NCOGEL .EQ. 3 ) THEN
C
C                 TRIANGLE 1 4 6
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(1)
                  MCN( MNN + 2 ) = NOEUDS(4)
                  MCN( MNN + 3 ) = NOEUDS(6)
                  MCN( MNN + 4 ) = 0
C                 TRIANGLE 4 2 5
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(4)
                  MCN( MNN + 2 ) = NOEUDS(2)
                  MCN( MNN + 3 ) = NOEUDS(5)
                  MCN( MNN + 4 ) = 0
C                 TRIANGLE 4 5 6
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(4)
                  MCN( MNN + 2 ) = NOEUDS(5)
                  MCN( MNN + 3 ) = NOEUDS(6)
                  MCN( MNN + 4 ) = 0
C                 TRIANGLE 6 5 3
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(6)
                  MCN( MNN + 2 ) = NOEUDS(5)
                  MCN( MNN + 3 ) = NOEUDS(3)
                  MCN( MNN + 4 ) = 0
C
               ELSE
C
C                 CONSTRUCTION DU BARYCENTRE DU QUADRANGLE
                  NBP = NBP + 1
                  NOEUDS(9) = NBP
C                 LES 2 COORDONNEES X Y DU BARYCENTRE
                  MNP = MNS + 3 * NBP
                  DO 84 K=1,2
                     D = 0D0
                     DO 82 L=1,NBPOE
                        D = D + PBASE(L) * RMCN(MNS+3*NOEUDS(L)+K)
 82                  CONTINUE
                     RMCN(MNP+K) = REAL( D )
 84               CONTINUE
C
C                 LA TEMPERATURE AU BARYCENTRE NOEUD=POINT AJOUTE
                  DO NCAS = NCAS0, NCAS1
                     TEMPS = TIMES( NCAS )
                     D = 0D0
                     DO 86 L=1,NBPOE
                        K = NOEUDS(L)
                        IF( NTYP .EQ. 0 ) THEN
                           T = TEMPER( K, NCAS )
                        ELSE
                           T = dptemp( NCAS )%dptab( K )
                        ENDIF
                        D = D + PBASE(L) * T
 86                  CONTINUE
                     MN = MNZT + ( NCAS - NCAS0 ) * NBSOM
                     RMCN(MN+NBP) = REAL( (D-TMIN) * ECHELLE )
C
                     IF( MODECO .EQ. 4 .OR. MODECO .EQ. 8 ) THEN
C
C                       CALCUL DE LA SOLUTION_EXACTE EN CE BARYCENTRE
C                       LA VALEUR DES PARAMETRES DE LA FONCTION
                        DPARAF(1) = TEMPS
                        DPARAF(2) = RMCN(MNP+1)
                        DPARAF(3) = RMCN(MNP+2)
C                       VALEUR DE Z
                        DPARAF(4) = 0D0
C                       NUMERO DE LA SEULE COMPOSANTE 1 A PRIORI
                        DPARAF(5) = 1D0
C                       TEMPERATURE_EXACTE(TEMPS,X,Y,Z) ou
C                       DEPLACEMENT_EXACT (TEMPS,X,Y,Z,NOCOMP)
C                       PRESSION_EXACTE(T,X,Y,Z)
C                       VITESSE_EXACTE(T,X,Y,Z,NOCOMP)
                        IF( NOFOUT .LE. 0 ) THEN
                          WRITE(IMPRIM,*)'NO TEMPERATURE_EXACTE NOFOUT='
     %                                   ,NOFOUT
                          GOTO 9999
                        ENDIF
                        IF( MODECO .EQ. 4 ) THEN
C
C                          1 SEULE FONCTION UTILISATEUR A CALCULER
                           CALL FONVAL( NOFOUT, NBARGS, DPARAF, NCODEV,
     %                                  DTX )
                        ELSE IF( MODECO .EQ. 8 ) THEN
C
C                          2 FONCTIONS UTILISATEUR A CALCULER
C                          ET MODULE A CALCULER
C                          FONCTION PARTIE_REELLE_EXACTE(t,x,y,z)
                           CALL FONVAL( NOFOWE(1), NBARGS, DPARAF,
     %                                  NCODEV,  DTX )
C                          FONCTION PARTIE_IMAGINAIRE_EXACTE(t,x,y,z)
                           CALL FONVAL( NOFOWE(2), NBARGS, DPARAF,
     %                                  NCODEV, IDTX )
C                          LE MODULE
                           DTX = SQRT( DTX**2 + IDTX**2 )
C
                        ENDIF
C
                        IF( NCODEV .GT. 0 ) THEN
C                          L'ERREUR NON ABSOLUE AU BARYCENTRE ET A L'ECHELLE
                           RMCN(MN+NBP) = REAL(((DTX-D)-TMIN)*ECHELLE)
                        ENDIF
C
                     ENDIF
                  ENDDO
C
C                 CALCUL DES 2 COORDONNEES X Y DU BARYCENTRE DU QUADRANGLE 2Q2C
C                 QUADRANGLE 1 5 9 8
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(1)
                  MCN( MNN + 2 ) = NOEUDS(5)
                  MCN( MNN + 3 ) = NOEUDS(9)
                  MCN( MNN + 4 ) = NOEUDS(8)
C                 QUADRANGLE 5 2 6 9
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(5)
                  MCN( MNN + 2 ) = NOEUDS(2)
                  MCN( MNN + 3 ) = NOEUDS(6)
                  MCN( MNN + 4 ) = NOEUDS(9)
C                 QUADRANGLE 8 9 7 4
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(8)
                  MCN( MNN + 2 ) = NOEUDS(9)
                  MCN( MNN + 3 ) = NOEUDS(7)
                  MCN( MNN + 4 ) = NOEUDS(4)
C                 QUADRANGLE 9 6 3 7
                  MNN = MNN + 4
                  MCN( MNN + 1 ) = NOEUDS(9)
                  MCN( MNN + 2 ) = NOEUDS(6)
                  MCN( MNN + 3 ) = NOEUDS(3)
                  MCN( MNN + 4 ) = NOEUDS(7)
C
               ENDIF
C
            ENDIF
C
 95      CONTINUE
 99   CONTINUE
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNNSEF) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNNSEF + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C     SELON LE CAS LA TEMPERATURE OU L'ERREUR EST PORTEE EN Z
C     =======================================================
      ZZZMIN = 1E28
      ZZZMAX =-1E28
      ERRMIN = 1E28
      ERRMAX =-1E28
      IF( MODECO .NE. 4 .AND. MODECO .NE. 8 ) THEN
C        LA SOLUTION EST MISE A L'ECHELLE AUX NOEUDS DE L'OBJET
         DO I=1,NBPOI
            DO NCAS = NCAS0, NCAS1
               MN = MNZT + ( NCAS - NCAS0 ) * NBSOM
               IF( NTYP .EQ. 0 ) THEN
                  T = TEMPER( I, NCAS )
               ELSE
                  T = dptemp( NCAS )%dptab( I )
               ENDIF
               RMCN(MN+I) = REAL( (T-TMIN) * ECHELLE )
               IF( RMCN(MN+I) .LT. ZZZMIN ) ZZZMIN = RMCN(MN+I)
               IF( RMCN(MN+I) .GT. ZZZMAX ) ZZZMAX = RMCN(MN+I)
            ENDDO
         ENDDO
      ELSE
C        ERREUR MISE A L'ECHELLE
         DO I=1,NBPOI
            DO NCAS = NCAS0, NCAS1
               MN = MNZT + ( NCAS - NCAS0 ) * NBSOM
               IF( NTYP .EQ. 0 ) THEN
                  T = TEMPER( I, NCAS )
               ELSE
                  T = dptemp( NCAS )%dptab( I )
               ENDIF
               ERRT = REAL( RMCN(MN+I)-T )
               RMCN(MN+I) = (ERRT-TMIN) * ECHELLE
               ERRMIN = MIN( ERRMIN, ERRT )
               ERRMAX = MAX( ERRMAX, ERRT )
               IF( RMCN(MN+I) .LT. ZZZMIN ) ZZZMIN = RMCN(MN+I)
               IF( RMCN(MN+I) .GT. ZZZMAX ) ZZZMAX = RMCN(MN+I)
            ENDDO
         ENDDO
      ENDIF
C
C     LCRITR = -2 => TRACE T31FCO DE LA SURFACE EN ARC EN CIEL
C                    SELON LA VALEUR DE ZZZMIN et ZZZMAX
      LCRITR = -2
C
C     MISE A JOUR DU TABLEAU COOEXT DU COMMON / XYZEXT /
      COOEXT(1,1) = XMIN
      COOEXT(1,2) = XMAX
      COOEXT(2,1) = YMIN
      COOEXT(2,2) = YMAX
ccc      COOEXT(3,1) = (TMIN-TMIN) * ECHELLE + DECALZ  3/2/2010
      COOEXT(3,1) = (TMIN-TMIN) * ECHELLE
      COOEXT(3,2) = (TMAX-TMIN) * ECHELLE
cccc
cccc      MISE A LA VALEUR INDIQUEE EN Z => AMPLITUDE EN Z
cccc      COOEXT(3,2) = 1.23 * ECHELLE       pour ther2d2abc de CHEN-ZHOU
cccc      COOEXT(3,2) = 3.30 * ECHELLE       pour ther2d3abc de CHEN-ZHOU
cccc      COOEXT(3,2) = 15.0 * ECHELLE       pour ther2d1abc de CHEN-ZHOU
C
C     MISE EN ECRAN POUR LE TRACE DE LA SURFACE=TEMPERATURE
C     SAUVEGARDE AVEC OU SANS CRITERE DE QUALITE
      LCRIT0 = LCRITR
C     IAVFAC : avec ou non trace des faces
      IAVFA0 = IAVFAC
C     PREDUF : pourcentage de reduction des faces
      PREDU0 = PREDUF
C     NEPARF : nombre d'epaisseurs des ARETES des FACES
      NEPAR0 = NEPARF
C     NCOUAF : COULEUR des ARETES des FACES
      NCOUA0 = NCOUAF
C     NCOUFA : COULEUR des FACES
      NCOUF0 = NCOUFA
C     IAVELO : avec ou sans eloignement
      IAVEL0 = IAVELO
      IAVELO = 0
C     LE TYPE DE LA VISEE
      NOTYVI = 0
C
      IF( MODECO .NE. 4 .AND. MODECO .NE. 8  .AND.
     %    MODECO .NE.11 .AND. MODECO .NE.12 ) THEN
C        LA PALETTE ARC EN CIEL POUR TRACER LA TEMPERATURE=Z
         CALL PALCDE( 11 )
      ELSE
C        LA PALETTE ORANGE ROUGE POUR TRACER LES ERREURS
         CALL PALCDE( 5 )
C        IAVELO sans eloignement
         IAVELO = 0
      ENDIF
C
C     LONGITUDE ET LATITUDE POUR VOIR "NATURELLEMENT" LE MAILLAGE
      AXOLON = -70.0
      AXOLAT =  18.0
C     DEFINITION DE PTV ET OEIL A PARTIR DE AXOLON ET AXOLAT
      CALL LONLAT( AXOLON, AXOLAT )
C     L'AXONOMETRIE
      CALL MATAXO
      CALL ISOFENETRE( -AXOLAR, AXOLAR, -AXOHAU, AXOHAU )
      NOTYVI = 0
      CALL VISEE0
C
C     LONGITUDE ET LATITUDE POUR VOIR "NATURELLEMENT" LE MAILLAGE
      AXOLON = -85.0
      AXOLAT =  18.0
C     DEFINITION DE PTV ET OEIL A PARTIR DE AXOLON ET AXOLAT
      CALL LONLAT( AXOLON, AXOLAT )
C     L'AXONOMETRIE
      CALL MATAXO
      CALL ISOFENETRE( -AXOLAR, AXOLAR, -AXOHAU, AXOHAU )
      NOTYVI = 11
C
C     ===========================================
C     TRACE DE LA SURFACE AVEC Z=TEMPERATURE(X,Y)
C     ===========================================
 100  CALL LIMTCL( 'valztxy', NMTCL )
      IF( NMTCL .LE.  0 ) GOTO 9900
      IF( NMTCL .EQ. 90 ) GOTO 300
      GOTO( 300, 120, 130, 140, 150, 160, 170, 180 ), NMTCL
C
C     bascule TRACE ou NON des ARETES des FACES
C     .........................................
 120  IF( IAVARE .EQ. 0 ) THEN
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
            KERR(1) = 'NO DRAWING of EDGES'
         ENDIF
         CALL LERESU
      ENDIF
      GOTO 100
C
C     POURCENTAGE de REDUCTION des FACES
C     ..................................
 130  CALL INVITE( 130 )
      CALL LIRRSP(  NCVALS , R )
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
C     bascule TRACE ou NON des FACES
C     ..............................
 140  IF( IAVFAC .EQ. 0 ) THEN
         IAVFAC = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TRACE des FACES'
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
            KERR(1) = 'NO DRAWING of FACES'
         ENDIF
         CALL LERESU
      ENDIF
      GOTO 100
C
C     COULEUR des ARETES du MAILLAGE
C     ..............................
 150  CALL LIMTCL( 'couleur0' , I )
      IF( I .EQ. -1 ) THEN
         GOTO 100
      ELSE IF( I .EQ. -2 ) THEN
         NCOUAF = -2
      ELSE IF( I .EQ. 0 ) THEN
C        COULEUR NOIRE
         NCOUAF = 0
      ELSE
C        COULEUR RESERVEE
         NCOUAF = N1COEL + I
      ENDIF
      GOTO 100
C
C     TYPE du TRAIT des ARETES du MAILLAGE
C     ....................................
 160  CALL LIMTCL( 'typtrait' , I )
      IF( I .EQ. -1 ) GOTO 100
      NTLAFR = I
      GOTO 100
C
C     COULEUR des ARETES sur la TEMPERATURE EN Z
C     ..........................................
 170  CALL LIMTCL( 'couleur0' , I )
      IF( I .EQ. -1 ) GOTO 100
      IF( I .EQ. -2 ) THEN
C        COULEUR INVISIBLE A NE PAS TRACER
         NCOAPL = -2
      ELSE IF( I .EQ. 0 ) THEN
C        LA COULEUR NOIRE
         NCOAPL = 0
      ELSE
         NCOAPL = N1COEL + I
      ENDIF
      GOTO 100
C
C     TYPE du TRAIT des ARETES des SURFACES ISOTHERMES
C     ................................................
 180  CALL LIMTCL( 'typtrait' , I )
      IF( I .EQ. -1 ) GOTO 100
      NTLAPL = I
      GOTO 100
C
C     EXECUTION DU TRACE DU MAILLAGE ET DE LA TEMPERATURE EN Z
C     ========================================================
C     OPTIONS DE LA VISEE
 300  CALL VISE3D( NMTCL )
      IF( NMTCL .LT. 0 ) GOTO 100
C
C     LE NUMERO DANS LE LEXIQUE DE L'OBJET KNOMOB
      CALL LXNMNO( NTOBJE, KNOMOB, NUOBJT, I )
C
C     TRACE DES AXES
      NTRAXZ = 1
      IF( MODECO .EQ. 4 .OR. MODECO .EQ. 8 ) THEN
         ZMIAXZ = ERRMIN
         ZMXAXZ = ERRMAX
      ELSE
ccc         ZMIAXZ = TMIN0 + DECALZ  3/2/2010
         ZMIAXZ = TMIN0
         ZMXAXZ = TMAX0
      ENDIF
C
C     INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
      IF( LORBITE .NE. 0 ) THEN
         CALL ORBITE0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 300
      ENDIF
C
C     LES DIFFERENTS CAS
 400  DO NCAS = NCAS0, NCAS1
c
cccc        INITIALISATION DU TRACE PS
ccc         CALL XVINITIERPS( 1 )
C
C        L'ECRAN PIXELS EST EFFACE (PAS DE SCINTILLEMENT)
         CALL EFFACEMEMPX

C        TRACE EFFECTIF DES AXES
         NETAXE = 0
C
C        LE TEMPS
         TEMPS = TIMES( NCAS )
C
C        ADRESSE DE LA SOLUTION
         MN = MNZT + ( NCAS - NCAS0 ) * NBSOM
C
C        TRACE DES AXES 3D
         CALL TRAXE3
C
cccC     TRACE FIL DE FER DES ARETES DU MAILLAGE DANS LE PLAN XOY 3/2/2010
cccC     --------------------------------------------------------
cccC     ADRESSE MCN DES XYZ DES SOMMETS
ccc      MNS = MNXYZS + WYZSOM
ccc      DO 310 I=1,NBSOM
cccC        Z EST LE DECALAGE EN Z POUR TRACER LE MAILLAGE DANS UN PLAN // XOY
ccc         RMCN( MNS + 2 ) = DECALZ
ccc         MNS = MNS + 3
ccc 310  CONTINUE
cccC     NEPARF : nombre d'epaisseurs des ARETES des FACES
ccc      NEPARF = 0
ccc      CALL XVEPAISSEUR( NEPARF )
cccC     NCOUAF : COULEUR des ARETES des FACES
cccC     TRACE DE LA QUALITE DES EF
ccc      LCRITR0 = LCRITR
ccc      LCRITR = 1
ccc      CALL T31FCO( KNOMOB, NUOBJT, MNNSEF, MNXYZS )
ccc      LCRITR = LCRITR0
C
C        LE TRACE DE LA SURFACE AVEC Z=TEMPERATURE * ECHELLE POUR NCAS
C        -------------------------------------------------------------
C        LA TEMPERATURE EST IMPOSEE EN Z AVEC UNE
C        MISE A L'ECHELLE DE LA PLUS GRANDE DIFFERENCE DE COORDONNEE
C        ADRESSE MCN DES XYZ DES SOMMETS
         MNS = MNXYZS + WYZSOM
         DO 320 I=1,NBSOM
C           Z DU NOEUD EST LA TEMPERATURE MISE A L'ECHELLE
            RMCN( MNS + 2 ) = RMCN( MN + I )
            MNS = MNS + 3
 320     CONTINUE
C
C        NCOUAF : COULEUR des ARETES des FACES SI TRACE DEMANDE
         NC     = NCOUAF
         NCOUAF = NCOAPL
C        NEPARF : nombre d'epaisseurs des ARETES des FACES
         NEPARF = 0
         CALL XVEPAISSEUR( NEPARF )
         CALL T31FCO( KNOMOB, NUOBJT, MNNSEF, MNXYZS )
         NCOUAF = NC
C
C        LE TRACE DE LA LEGENDE : COULEURS => VALEURS
C        --------------------------------------------
         IF( IAVFAC .GT. 0 ) THEN
C           TRACE SI LES FACES SONT TRACEES
            CALL LEGCOULSO( TMIN, TMAX )
         ENDIF
C
C        TRACE DE LA 2-EME LIGNE DU TITRE DU TRACE
C        -----------------------------------------
         CALL TIT2LG( KNOMOB, MODECO )
C
C        CONSTRUCTION DE LA 3-EME ET DERNIERE LIGNE DU TITRE DU TRACE
C        ------------------------------------------------------------
         CALL TIT3LG( MODECO, NCAS, TEMPS, TMIN, TMAX,  KNOM )
C
C        TRACE DU TITRE DU TRACE
C        -----------------------
         CALL TRFINS( KNOM )
C
C        MISE SUR FICHIER NomfgifNoImage.xwd puis NomfgifNoImage.jpg
C        DE LA PIXMAP de la FENETRE X11 ACTUELLE SI VIDEO DEMANDEE
C        -----------------------------------------------------------
         CALL VIDEO1( NOMFGIF, NCAS )
C
C        ATTENDRE POUR PERMETTRE LA LECTURE DU TRACE
         CALL ATTENDSEC( TEMP2TRAC )
C
C        FIN DE LA BOUCLE SUR LES CAS
      ENDDO
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
C     RESTAURATION DES OPTIONS INITIALES DE TRACE
C     -------------------------------------------
 9900 COOEXT(1,1) = XMIN0
      COOEXT(1,2) = XMAX0
      COOEXT(2,1) = YMIN0
      COOEXT(2,2) = YMAX0
      COOEXT(3,1) = ZMIN0
      COOEXT(3,2) = ZMAX0
      NOTYVI = 0
C     RESTAURATION DES PARAMETRES INITIAUX
      LCRITR = LCRIT0
C     IAVFAC : avec ou non trace des faces
      IAVFAC = IAVFA0
C     PREDUF : pourcentage de reduction des faces
      PREDUF = PREDU0
C     NEPARF : nombre d'epaisseurs des ARETES des FACES
      NEPARF = NEPAR0
C     NCOUAF : COULEUR des ARETES des FACES
      NCOUAF = NCOUA0
C     NCOUFA : COULEUR des FACES
      NCOUFA = NCOUF0
C     IAVELO : avec ou sans eloignement
      IAVELO = IAVEL0
      NTRAXZ = 0
      NDIM   = 2
C
C     DESTRUCTION DES TMC INUTILES
 9999 IF( MNZTXY .GT. 0 ) CALL TNMCDS( 'REEL', MOZTXY, MNZTXY )
      IF( MNXYZS .GT. 0 ) CALL TNMCDS( 'MOTS', WYZSOM+3*NBSOM,  MNXYZS )
      IF( MNNSEF .GT. 0 ) CALL TNMCDS( 'MOTS', WUSOEF+4*NBEFOB, MNNSEF )
C
      RETURN
      END
