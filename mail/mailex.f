      SUBROUTINE MAILEX( NMTOBJ, NOMOBJ,
     %                   NTNSEF, MNNSEF, NTSOMM, MNSOMM, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUIRE LES TMS XYZSOMMET ET NSEF D'UN PLSV
C -----    A PARTIR DE SON TMS DEFINITION
C
C ATTENTION : TRAITEMENT RECURSIF DES UNIONS ET TRANSFORMATIONS ...
C =========== C-A-D DES PLSV DEFINIS PAR D'AUTRES PLSV
C
C ENTREES:
C --------
C NMTOBJ : NOM DU TYPE DE L'OBJET 'POINT' 'LIGNE' ... 'VOLUME'
C NOMOBJ : NOM DU PLSV
C
C SORTIES:
C --------
C NTNSEF : NUMERO DU TMS 'NSEF' POUR LIGNE SURFACE VOLUME
C          0 POUR UN POINT
C MNNSEF : ADRESSE MCN DU TMS 'NSEF' DU LSV
C          0 POUR UN 'POINT'
C NTSOMM : NUMERO DU TMS 'XYZSOMMET' DU PLSV
C MNSOMM : ADRESSE MCN DU TMS 'XYZSOMMET' DU PLSV
C IERR   : 0  SI PAS D'ERREUR A LA GENERATION DU MAILLAGE
C          1  OBJET INCONNU
C          2  TABLEAU 'DEFINITION' INCONNU DU PLSV
C          3  TYPE PLSV INCONNU
C          4  OBJET DE NUMERO INCORRECT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/pp.inc"
      include"./incl/langue.inc"
      include"./incl/majmai.inc"
      include"./incl/gsmenu.inc"
      include"./incl/typnoobj.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
C
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*(*)     NMTOBJ,NOMOBJ
      CHARACTER*8       MMTOBJ,MMTOB
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (RMCN(1),MCN(1))
      CHARACTER*10      NMTYOB,NMTOB
      CHARACTER*24      KNOM
      LOGICAL           AVANT
      INTEGER           NUMOB, NUMOBJ
      EQUIVALENCE      (NUMOB, NUMOBJ)
      INTEGER           NUTYOB
C
C     INITIALISATION
CCC      WRITE(IMPRIM,*) 'entree mailex avec ',nomobj
CCC      WRITE(IMPRIM,*) (MCN(I),I=1,5)
      IERR   = 0
      MAFAIT = 0
      NTNSEF = 0
      MNNSEF = 0
      NTSOMM = 0
      MNSOMM = 0
      MNPILE = 0

C     REMPLISSAGE DES COMMONS de ./incl/typnoobj.inc
      NUMTYPOBJ = 0
      NUMOBJLX  = 0
      KNMTYPOBJ = NMTOBJ
      KNMOBJLX  = NOMOBJ

C     NUMERO DU TYPE D'OBJET A MAILLER
      I = NUDCNB( NMTOBJ )
      MMTOBJ = NMTOBJ(1:I)
      IF( INDEX( MMTOBJ, 'LIGNE' ) .GT. 0 ) MMTOBJ = 'LINE'
      IF( INDEX( MMTOBJ, 'OBJET' ) .GT. 0 ) MMTOBJ = 'OBJECT'
      NUTYOB = NTYOBJ( NMTOBJ )
      IF( NUTYOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'MAILEX: TYPE INCORRECT d''OBJET '//NMTOBJ
         ELSE
            KERR(1) = 'MAILEX: INCORRECT TYPE of OBJECT '//NMTOBJ
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9999
      ENDIF
C
C     NUMERO TS DU LEXIQUE DES OBJETS DE CE TYPE D'OBJET
      NTLX = NTMN( NUTYOB )
C
C     NUMERO TS ET ADRESSE DU LEXIQUE DE CET OBJET
      CALL LXLXOU( NTLX, NOMOBJ, NTLXOB, I )
ccc      print *,'mailex: ntlxob=',ntlxob
      IF( NTLXOB .LE. 0 ) THEN
C        OBJET INCONNU
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'MAILEX: '//MMTOBJ//' INCONNU : '//NOMOBJ
         ELSE
            KERR(1) = 'MAILEX: '//MMTOBJ//' UNKNOWN : '//NOMOBJ
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9999
      ENDIF
C
C     LE NUMERO DE L'OBJET DANS SON LEXIQUE
      CALL LXNMNO( NTLX, NOMOBJ, NUMOBJ, I )
ccc      print *,'mailex: nomobj=',nomobj,' numobj=',numobj,' i=',I

C     REMPLISSAGE DES COMMONS de ./incl/typnoobj.inc
      NUMTYPOBJ = NUTYOB 
      NUMOBJLX  = NUMOBJ
ccc   KNMTYPOBJ = NMTOBJ
ccc   KNMOBJLX  = NOMOBJ

C     DECLARATION D'UNE PILE DES OBJETS A TRAITER
      MXPILE = 1024
      CALL TNMCDC( 'ENTIER',MXPILE,MNPILE )
      LHPILE = 0
C
C     L'OBJET A TRAITER EST EMPILE:
C     SON TYPE NUTYOB
C     PUIS SON NUMERO NUMOBJ DANS LE TYPE
C          (1:POINT 2:LIGNE ... 5:OBJET 6:FONCTIONS 7:TRANSFORMATIONS)
C     PUIS L'ABSENCE D'OBJET QU'IL QUALIFIE (DONT IL EST UN CONSTITUANT)
      MNDFSO = 0
      CALL PILEOB( NUTYOB, NUMOBJ, MNDFSO ,
     %             LHPILE, MXPILE, MNPILE )
ccc      print *,'mailex: lhpile=',lhpile,'  mnpile=',mnpile
C
C     ==================================================================
C     TANT QUE LA PILE EST NON VIDE TRAITER L'OBJET DU HAUT DE LA PILE
C     ==================================================================
 10   IF( IERR   .GT. 0 ) GOTO 9999
      IF( LHPILE .GT. 0 ) THEN
C
C        LE TYPE ET LE NUMERO DE L'OBJET A TRAITER
         I      = MNPILE - 1 + LHPILE
C        L'ADRESSE DU TABLEAU DEFINITION DU SUPER-OBJET
C        DONT L'OBJET COURANT EST UN CONSTITUANT
         MNDFSO = MCN( I     )
C        LE NUMERO DE L'OBJET DANS LE LEXIQUE DE CE TYPE D'OBJETS
         NUMOB  = MCN( I - 1 )
C        LE NUMERO DU TYPE DE L'OBJET
C       (1:POINT 2:LIGNE ... 5:OBJET 6:FONCTIONS 7:TRANSFORMATIONS)
         NUTYO  = MCN( I - 2 )
         NUTYA = ABS( NUTYO )
C        PROTECTION
ccc         print *,'mailex: mnpile=',mnpile,' mndfso=',mndfso,
ccc     %           ' numob=',numob,' nutyo=',nutyo
         IF( NUMOB .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(MXLGER-1)(1:4),'(I4)') NUMOB
            WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYA
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='MAILEX: NO OBJET'//KERR(MXLGER-1)(1:4)//' <=0'
               KERR(2)='de TYPE '//KERR(MXLGER)(1:4)
            ELSE
               KERR(1)='MAILEX: OBJECT NUMBER '//KERR(MXLGER-1)(1:4)//
     %                 ' <=0'
               KERR(2)='of TYPE '//KERR(MXLGER)(1:4)
            ENDIF
            CALL LEREUR
            IERR = 4
            GOTO 9999
         ENDIF
C
C        L'OBJET EST DEPILE
         LHPILE = LHPILE - 3
C
CCCC        SI L'OBJET EST UNE TRANSFORMATION
CCCC        ALORS RIEN N'EST FAIT. ELLE EST RECALCULEE A CHAQUE FOIS
CCC         IF( NUTYA .EQ. 7 ) GOTO 10
C
C        LE NOM CORRESPONDANT AU TYPE DE CET OBJET COURANT
         NMTOB = NMTYOB( NUTYA )
         MMTOB = NMTOB(1:8)
         IF( INDEX( MMTOBJ, 'LIGNE' ) .GT. 0 ) MMTOB = 'LINE'
         IF( INDEX( MMTOBJ, 'OBJET' ) .GT. 0 ) MMTOB = 'OBJECT'
C
C        ADRESSE DU LEXIQUE DES OBJETS DE TYPE NUTYA
         CALL TAMSOU( NTMN(NUTYA), MNLXO )
C        OUVERTURE DU LEXIQUE NTOB DE L'OBJET DE NUMERO NUMOB
         MN = MNLXO + MCN(MNLXO) * NUMOB
C        LE TABLEAU MS LEXIQUE DE CET OBJET
         NBENN = MCN( MNLXO + 2 )
         NTLXO = MCN( MN + NBENN + 2 )
C        LE NOM KNOM DE L'OBJET CONVERTI D'ENTIERS EN CARACTERES
         CALL ENTNOM( NBENN, MCN(MN), KNOM )
         IF( NTLXO .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'MAILEX: '//MMTOB//' INCONNU: '//KNOM
            ELSE
               KERR(1) = 'MAILEX: '//MMTOB//' UNKNOWN: '//KNOM
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 9999
         ENDIF
C
         IF( NUTYA .NE. 6 ) THEN
C           OUVERTURE S'IL EXISTE DU TABLEAU 'DEFINITION' DE CET OBJET
            CALL LXTSOU( NTLXO,'DEFINITION', NTDFO,MNDFO )
            IF( NTDFO .LE. 0 ) THEN
C              LE TABLEAU 'DEFINITION' N'EXISTE PAS
               NBLGRC(NRERR) = 2
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) ='MAILEX: '//MMTOB//' SANS DEFINITION:'
               ELSE
                  KERR(1) ='MAILEX: '//MMTOB//' WITHOUT DEFINITION:'
               ENDIF
               KERR(2) = KNOM
               CALL LEREUR
               IERR = 2
               GOTO 9999
            ENDIF
C
C           LE TYPE DEFINITION DE L'OBJET
            NTYDFO = MCN( MNDFO + WUTYLI )
C
C           LE NUMERO DE LA TRANSFORMATION
            IF( NUTYA .NE. 5 ) THEN
               NUTRAN = MCN( MNDFO + WTYTRL )
            ELSE
               NUTRAN = 1
            ENDIF
         ELSE
C           OUVERTURE S'IL EXISTE DU TABLEAU 'ARBRE' DE CETTE FONCTION
            CALL LXTSOU( NTLXO,'ARBRE', NTDFO,MNDFO )
            IF( NTDFO .LE. 0 ) THEN
C              LE TABLEAU 'ARBRE' N'EXISTE PAS
               NBLGRC(NRERR) = 2
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) ='MAILEX: FONCTION SANS ARBRE'
               ELSE
                  KERR(1) ='MAILEX: FUNCTION WITHOUT ''ARBRE'''
               ENDIF
               KERR(2) = KNOM
               CALL LEREUR
               IERR = 2
               GOTO 9999
            ENDIF
C           LE TYPE DEFINITION DE LA FONCTION
            NTYDFO = 0
C           LE NUMERO DE LA TRANSFORMATION
            NUTRAN = 1
         ENDIF
C
C        SI PAS DE MISE A JOUR AUTOMATIQUE DES MAILLAGES ALORS ALLER EN 500
         IF( MAJAUT .EQ. 0 ) GOTO 500
C
C        SI SECOND TRAITEMENT DE L'OBJET ALORS ALLER EN 500
         IF( NUTYO .LT. 0 ) GOTO 500
C
C***********************************************************************
C
C        PREMIER TRAITEMENT DE L'OBJET COURANT
C        SES COMPOSANTS SONT GENERALEMENT EMPILES
C
C***********************************************************************
C
C        ICI LE MAILLAGE DE L'OBJET A TRAITER N'EXISTE PAS
C        -------------------------------------------------
C        L'OBJET EST EMPILE AVEC UN TYPE OBJET DE SIGNE NEGATIF
         CALL PILEOB( -NUTYA, NUMOB,  MNDFSO,
     %                LHPILE, MXPILE, MNPILE )
C        LA HAUTEUR DANS LA PILE DE L'OBJET COURANT
         MNDFSO = MNDFO
C
C        RECHERCHE DES POINTEURS SUR LES LEXIQUES DU TABLEAU DEFINITION
CCC         IF( NUTYA .EQ. 7 ) THEN
C           TRANSFORMATION : RIEN A FAIRE .
C                            ELLE EST RECALCULEE A CHAQUE FOIS
CCC            GOTO 10
CCC         ENDIF
         CALL XXTSTD( NUTYA,NUMOB,NBPTLX,MNPTLX )
C
C        EMPILAGE DES COUPLES ( NO TYPE,NO OBJET )
         MN = MNPTLX
         DO 20 J = 1, NBPTLX
            IF( MCN(MN) .NE. 7 .OR. MCN(MN+1) .NE. 1 ) THEN
               CALL PILEOB( MCN(MN), MCN(MN+1), MNDFSO,
     %                      LHPILE,  MXPILE,    MNPILE )
            ENDIF
            MN = MN + 2
 20      CONTINUE
C
C        DESTRUCTION DU TABLEAU INTERMEDIAIRE
         IF( NBPTLX .GT. 0 ) THEN
            CALL TNMCDS( 'ENTIER',NBPTLX*2,MNPTLX )
         ENDIF
         GOTO 10
C
C***********************************************************************
C
C        ICI TOUS LES COMPOSANTS DE L'OBJET ONT ETE MIS A JOUR
C        L'OBJET EST DONC A MAILLER SI CE N'EST PAS DEJA FAIT
C
C***********************************************************************
C
C        LA DATE DE DEFINITION DU SUPER-OBJET DE L'OBJET COURANT EST
C        MISE AU MAX( DATE OBJET COURANT,DATE SUPER OBJET )
 500     IF( MNDFSO .GT. 0 ) THEN
            IF( AVANT( MCN(MNDFSO),MCN(MNDFO) ) ) THEN
C              SI LE SUPER-OBJET A UN CONSTITUANT PLUS JEUNE QUE LUI
C              ALORS LA DATE DU SUPER-OBJET DEVIENT LA DATE ACTUELLE
               CALL ECDATE( MCN(MNDFSO) )
            ENDIF
         ENDIF
C
C        OUVERTURE S'IL EXISTE DU TABLEAU 'XYZSOMMET'
         IF( NUTYA .GE. 1 .AND. NUTYA .LE. 5 ) THEN
            GOTO 30
C
         ELSE IF( NUTYA .EQ. 6 .OR. NUTYA .EQ. 7 ) THEN
C           FONCTION OU TRANSFORMATION
C           PAS D'AUTRE TABLEAU QUE CELUI DE DEFINITION
            GOTO 10
C
         ELSE
            NBLGRC(NRERR) = 1
            WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYA
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'MAILEX: TYPE INCONNU: '//KERR(MXLGER)(1:4)
            ELSE
               KERR(1) = 'MAILEX: UNKNOWN TYPE: '//KERR(MXLGER)(1:4)
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 9999
         ENDIF
C
C        OUVERTURE DU TABLEAU 'XYZSOMMET'
 30      CALL LXTSOU( NTLXO,'XYZSOMMET',NTSOM ,MNSOM  )
C
         IF( NTSOM  .GT. 0 ) THEN
C           LE TABLEAU 'XYZSOMMET' EXISTE
            IF( AVANT( MCN(MNDFO),MCN(MNSOM) ) ) THEN
C              LA DEFINITION A PRECEDE LE MAILLAGE
               GOTO 10
            ELSE
C              LA DEFINITION A ETE MODIFIEE DEPUIS LE MAILLAGE
               IF( NUTYA .NE. 4 .OR. NTYDFO .NE. 54 ) THEN
C                 CONSERVATION DU MAILLAGE SI VOLUME OPTION 54
C                 C-A-D LE VOLUME MULTI-MATERIAUX DEVIENT MONO-MATERIAU
C                 DESTRUCTION DE L'ANCIEN MAILLAGE
                  CALL MAILDS( NMTOB, KNOM )
                  NTSOM = 0
                  MNSOM = 0
               ENDIF
            ENDIF
         ENDIF
C
C        LE MAILLAGE DOIT ETRE FAIT
C        --------------------------
         MAFAIT = 1
         IF( NTYDFO .GE. 0  .AND.
     %      (NTYDFO .LT. 50 .OR. NTYDFO .EQ. 53
     %                      .OR. NTYDFO .EQ. 54) ) THEN
            IF( NUTYA .EQ. 1 ) THEN
C
C              POINT
C              ==========================================================
               IF( NTYDFO .EQ. 1 ) THEN
C
C                 POINT = X Y Z
C                 =============
                  CALL POEX01( NTLXO, RMCN(MNDFO),
     %                         NTSOM ,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 2 ) THEN
C
C                 POINT = R TETA Z
C                 ================
                  CALL POEX02( NTLXO, RMCN(MNDFO),
     %                         NTSOM ,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 3 ) THEN
C
C                 POINT = R TETA PHI
C                 ==================
                  CALL POEX03( NTLXO, RMCN(MNDFO),
     %                         NTSOM ,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 35 ) THEN
C
C                 POINT = SOMMETS EXTRAITS D'UNE LIGNE
C                 ====================================
                  CALL POEX35( NTLXO, MCN(MNDFO) ,
     %                         NTSOM ,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 36 ) THEN
C
C                 POINT = SOMMETS EXTRAITS SUR CRITERE D'UN PLSV
C                 ==============================================
                  CALL POEX36( NTLXO, MCN(MNDFO) ,
     %                         NTSOM ,MNSOM,IERR )
C
               ENDIF
C
            ELSE IF( NUTYA .EQ. 2 ) THEN
C
C              LIGNE
C              ==========================================================
               IF( NTYDFO .EQ. 2 ) THEN
C
C                 DROITE
C                 ======
                  CALL LIEX02( NTLXO, MCN(MNDFO),RMCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 3 .OR. NTYDFO .EQ. 4 ) THEN
C
C                 ARC_CERCLE DEFINI PAR 3 POINTS DE R**3
C                 ARC_CERCLE DEFINI PAR 2 POINTS DE R**3 RAYON ET UN POINT PLAN
C                 =============================================================
                  CALL LIEX03( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 5 ) THEN
C
C                 ARC D'ELLIPSE 2D
C                 ================
                  CALL LIEX05( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 6 ) THEN
C
C                 POINTS_UTILISATEUR_DE_LA_LIGNE
C                 ==============================
                  CALL LIEX06( NTLXO, MCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 7 ) THEN
C
C                 POINTS_DE_LA_LIGNE_DEFINIS PAR 3 COORDONNEES
C                 ============================================
                  CALL LIEX07( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 8 ) THEN
C
C                 CERCLE 3D
C                 =========
                  CALL LIEX08( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 9 ) THEN
C
C                 ELLIPSE 2D
C                 ==========
                  CALL LIEX09( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 10 ) THEN

C                 DONNEE DES TMS XYZSOMMET et NSEF
C                 ================================
                  CALL TAMSOU( NTLXO, M )
                  NTSOM  = MCN( MNDFO + WUTSOL )
                  CALL LXTSOU( NTLXO, 'XYZSOMMET', NTSOM, MNSOM )
                  NTNSEO = MCN( MNDFO + WUTSSL )
                  CALL LXTSOU( NTLXO, 'NSEF', NTNSEO, MNNSEO )
                  IERR = 0

               ELSE IF( NTYDFO .EQ. 11 ) THEN
C
C                 B-SPLINE UNIFORME OUVERTE D'INTERPOLATION
C                 =========================================
                  CALL LIEX11( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 21 ) THEN
C
C                 B-SPLINE UNIFORME FERMEE D'INTERPOLATION
C                 ========================================
                  CALL LIEX21( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 12 .OR. NTYDFO .EQ. 13 .OR.
     %                  NTYDFO .EQ. 14 .OR. NTYDFO .EQ. 22 .OR.
     %                  NTYDFO .EQ. 23 ) THEN
C
C                 B-SPLINE UNIFORME OU NON, OUVERTE OU FERMEE,
C                 POLYNOMIALE OU RATIONNELLE
C                 ===================================================
                  CALL LIEX14( NTLXO, MCN(MNDFO),RMCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
              ELSE IF( NTYDFO .EQ. 20 ) THEN
C
C                 LIGNE EXTRUDEE A PARTIR D'UN POLYGONE DEFINI PAR DES POINTS
C                 ===========================================================
                  CALL LIEX20( NUMOB, NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
              ELSE IF( NTYDFO .EQ. 24 ) THEN
C
C                 LIGNE INTERSECTION PLAN TRONC DE CONE
C                 ===================================================
                  CALL LIEX24( NTLXO, MCN(MNDFO),RMCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
              ELSE IF( NTYDFO .EQ. 25 ) THEN
C
C                 LIGNE INTERSECTION PLAN ELLIPSOIDE
C                 ===================================================
                  CALL LIEX25( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 26 ) THEN
C
C                 INTERSECTION DE 2 TRONCS DE CONES
C                 =================================
                  CALL LIEX26( NTLXO,MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ELSE IF( NTYDFO .EQ. 27 ) THEN
C
C                 INTERSECTION DE 2 CYLINDRES D'AXES Y Z
C                 ======================================
                  CALL LIEX27( NTLXO,MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ELSE IF( NTYDFO .EQ. 28 ) THEN
C
C                 LIGNE CONGE
C                 ===================================================
                  CALL LIEX28( NTLXO, MCN(MNDFO),RMCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 30 ) THEN
C
C                 AMELIORATION DE LA QUALITE DES ARETES D'UNE LIGNE
C                 A PARTIR DE LA FONCTION UTILISATEUR TAILLE_IDEALE(x,y,z)
C                 ou EDGE_LENGTH(x,y,z) ou sinon LA TAILLE D'ARETE PAR DEFAUT
C                 ===========================================================
                  CALL LIEX30( NTLXO, MCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 35 ) THEN
C
C                 EXTRACTION D'ARETES  D'UNE LIGNE
C                 ================================
                  CALL LIEX35( NTLXO, MCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 36 ) THEN
C
C                 EXTRACTION D'UNE SOUS-LIGNE D'UNE LIGNE
C                 =======================================
                  CALL LIEX36( NTLXO, MCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 37 ) THEN
C
C                 EXTRACTION D'UNE LIGNE D'UNE SURFACE
C                 ====================================
                  CALL LIEX37( NTLXO, MCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 40 ) THEN
C
C                 RESTRUCTURATION D'UNE LIGNE NON STRUCTUREE
C                 ==========================================
                  CALL LIEX40( NTLXO,MCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ELSE IF( NTYDFO .EQ. 43 ) THEN
C
C                 INVERSION DE L'ORDRE DES SOMMETS D'UNE LIGNE STRUCTUREE
C                 =======================================================
                  CALL LIEX43( NTLXO,  MCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
ccc               ELSE IF( NTYDFO .EQ. 44 ) THEN
cccC
cccC                 RENUMEROTER LES SOMMETS D'UNE LIGNE FERMEE A PARTIR D'UN SO
cccC                 ===========================================================
ccc                  CALL LIEX44( NTLXO,  MCN(MNDFO),
ccc     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ELSE IF( NTYDFO .EQ. 48 ) THEN
C
C                 LIGNE 48 = IMPOSER DES POINTS COMME PLUS PROCHE SOMMET
C                 ======================================================
                  CALL LIEX48( NTLXO,  MCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
CCC               ELSE IF( NTYDFO .EQ. 49 ) THEN
C
C                 DESTRUCTURATION D'UNE LIGNE STRUCTUREE
C                 ======================================
CCCccc               CALL LIEX49( NTLXO,  MCN(MNDFO),
CCCccc     %                      NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
ccc                  CALL STNOST( NTLXO,
ccc     %                         NTNSEO, MNNSEO, IERR )
C
               ENDIF
C
            ELSE IF( NUTYA .EQ. 3 ) THEN
C
C              SURFACE
C              ==========================================================
               IF( NTYDFO .EQ. 1 .OR. NTYDFO .EQ. 2 ) THEN
C
C                 SURFACE 1 = QUADRANGLE ALGEBRIQUE STRUCTURE GORDON,..
C                 SURFACE 2 = QUADRANGLE STRUCTURE  GORDON+WINSLOW
C                 =====================================================
                  CALL SUEX01( NTLXO,MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .GE. 3 .AND. NTYDFO .LE. 4 ) THEN
C
C                 SURFACE 3 A 4 = QUADRANGLE B-SPLINE STRUCTURE D'INTERPOLATION
C                 =============================================================
                  CALL SUEX03( NTLXO,MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 5 ) THEN
C
C                 SURFACE 5 = TRIANGLE BEZIER
C                 ===========================
                  CALL SUEX05( NTLXO,MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 6 ) THEN
C
C                 SURFACE 6 = TRIANGLE STRUCTURE
C                 ==============================
                  CALL SUEX06( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 7 ) THEN
C
C                 SURFACE 7 = QUADRANGLE ALGEBRIQUE STRUCTURE GORDON
C                             SUIVI DE PROJECTIONS SUR CERCLES OU CYLINDRES
C                 =========================================================
                  CALL SUEX01( NTLXO,MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 8 ) THEN
C
C                 SURFACE 8 = QUADRANGULATION d'un CARRE ou CERCLE
C                             AVEC UN NOYAU DE n DIFFERENCES FINIES et
C                             des COUCHES GENEREES PAR HOMOTHETIE
C                             EN PROGRESSION GEOMETRIQUE
C                 ====================================================
                  CALL SUEX08( NTLXO,MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 9 ) THEN
C
C                 SURFACE 9 = TRIANGULATION DE L'ARBRE DES TRIANGLES EQUILATERAU
C                 ==============================================================
                  NUTY = 9
                  CALL SUEX09( NUTY, NTLXO, MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )

               ELSE IF( NTYDFO .EQ. 10 ) THEN

C                 DONNEE DES TMS XYZSOMMET et NSEF
C                 ================================
                  CALL TAMSOU( NTLXO, M )
                  NTSOM  = MCN( MNDFO + WUTSOS )
                  CALL LXTSOU( NTLXO, 'XYZSOMMET', NTSOM, MNSOM )
                  NTNSEO = MCN( MNDFO + WUTSSS )
                  CALL LXTSOU( NTLXO, 'NSEF', NTNSEO, MNNSEO )
                  IERR = 0

               ELSE IF( NTYDFO .EQ. 11 ) THEN

C                 SURFACE 11 = RECTANGLE
C                 ======================
                  CALL SUEX11( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )

               ELSE IF( NTYDFO .EQ. 12 ) THEN
C
C                 SURFACE 12 = TRIANGLE
C                 ======================
                  CALL SUEX12( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 13 ) THEN
C
C                 SURFACE 13 = QUADRANGULATION D'UNE LIGNE CONVEXE FERMEE
C                 =======================================================
                  CALL SUEX13( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 14 ) THEN
C
C                 SURFACE 14 = SURFACE GENEREE PAR COUCHES D'UNE LIGNE
C                 ====================================================
                  CALL SUEX14( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 15 ) THEN
C
C                 SURFACE 15 = SURFACE GENEREE PAR UNE LIGNE EN ROTATION
C                 ======================================================
                  CALL SUEX15( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 16 ) THEN
C
C                 SURFACE 16 = GROSSIR PAR COUCHES NORMALES UNE LIGNE
C                 ===================================================
                  CALL SUEX16( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 17 .OR. NTYDFO .EQ. 18  ) THEN
C
C                 SURFACE 17 18 = TRIANGULATION REGULIERE 1/8,1 SPHERE
C                 ======================================================
                  CALL SUEX18( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 19 ) THEN
C
C                 SURFACE 19 = TRIANGULATION PAR ETOILES DES POINTS DE LA FRONTI
C                 ==============================================================
                  NUTY = 19
                  CALL SUEX19( NUTY, NTLXO, MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
              ELSE IF( NTYDFO .EQ. 20 ) THEN
C
C                 SURFACE EXTRUDEE A PARTIR D'UN POLYGONE DEFINI PAR DES POINTS
C                 =============================================================
                  CALL SUEX20( NUMOB,  NTLXO, MCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
              ELSE IF( NTYDFO .EQ. 21 ) THEN
C
C                 SURFACE FERMEE des 6 FACES d'UNE BOITE
C                 ======================================
                  CALL SUEX21( NUMOB,  NTLXO, MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
              ELSE IF( NTYDFO .EQ. 22 ) THEN
C
C                 SURFACE FERMEE d'UN CYLINDRE ou CONE
C                 ====================================
                  CALL SUEX22( NUMOB,  MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
              ELSE IF( NTYDFO .EQ. 23 ) THEN
C
C                 SURFACE FERMEE d'UN TORE
C                 ========================
                  CALL SUEX23( NUMOB,  NTLXO, MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
              ELSE IF( NTYDFO .EQ. 25 ) THEN
C
C                 SURFACE FERMEE d'UN ARBRE et de ses RACINES
C                 ===========================================
                  CALL SUEX25( NUMOB,  MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ELSE IF( NTYDFO .EQ. 27 ) THEN
C
C                 SURFACE 27 = TRIANGULATION D'UN TRONC DE CONE ou CYLINDRE
C                 =========================================================
                  CALL SUEX27( NTLXO, MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 28 ) THEN
C
C                 SURFACE 28 = SURFACE AVEC CONGE
C                 ===============================
                  CALL SUEX28( NTLXO, MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ELSE IF( NTYDFO .EQ. 29 ) THEN
C
C                 SURFACE 29 = MODIFS INTERACTIVES DU MAILLAGE de SURFACE 2D ou
C                 ==============================================================
                  CALL SUEX29( NUMOB, NTLXO, MCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
ccc               ELSE IF( NTYDFO .EQ. 29 ) THEN
cccC
cccC                 SURFACE 29 = MOLECULE
cccC                 =====================
ccc                  CALL SUEX29( NTLXO, MCN(MNDFO),RMCN(MNDFO),
ccc     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 30 ) THEN
C
C                 SURFACE 30 = AMELIORATION DES ANGLES D'UNE TRIANGULATION
C                 ========================================================
                  CALL SUEX30( NUMOB, NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 31 ) THEN
C
C                 SURFACE 31 = TRIANGULATION D'UNE QUADRANGULATION
C                 ================================================
                  CALL SUEX31( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 32 ) THEN
C
C                 SURFACE 32 = JONCTION POINT ARETES D'UNE LIGNE
C                 ==============================================
                  CALL SUEX32( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 33 ) THEN
C
C                 SURFACE 33 = JONCTION DES ARETES DE PLUSIEURS LIGNES
C                 ====================================================
                  CALL SUEX33( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 34 ) THEN
C
C                 SURFACE 34 = EXTRACTION D'UNE SURFACE DANS UN HEXAEDRE STRUCTU
C                 ==============================================================
                  CALL SUEX34( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 35 ) THEN
C
C                 SURFACE 35 = EXTRACTION D'UNE SURFACE DANS UN QUADRANGLE STRUC
C                 ==============================================================
                  CALL SUEX35( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 36 ) THEN
C
C                 EXTRACTION D'UNE SOUS-SURFACE D'UNE SURFACE
C                 ===========================================
                  CALL SUEX36( NTLXO, MCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 37 ) THEN
C
C                 EXTRACTION D'UNE SURFACE D'UN VOLUME PAR FONCTION CRITERE
C                 =========================================================
                  CALL SUEX37( NTLXO, MCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 38 ) THEN
C
C                 EXTRACTION D'UNE SOUS-SURFACE CLIQUEE D'UNE SURFACE
C                 ===================================================
                  CALL SUEX38( NTLXO, MCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 39 ) THEN
C
C                 SURFACE 39 = DIVISION EN 4 SOUS EF DES EF 2D
C                 ============================================
                  CALL SUEX39( NTLXO,  MCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ELSE IF( NTYDFO .EQ. 40 ) THEN
C
C                 SURFACE EXTRAITE = SURFACE1 - SURFACE2
C                 ======================================
                  CALL SUEX40( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 41 ) THEN
C
C                 EXTRACTION D'UNE SOUS-SURFACE (NORMALE,VECTEUR)>=0
C                 ==================================================
                  CALL SUEX41( NTLXO, MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 42 ) THEN
C
C                 SURFACE RENDUE G1-CONTINUE PAR PROJECTION
C                 DES TANGENTES SUR LE PLAN A DISTANCE MINIMALE
C                 =============================================
                  CALL SUEX42( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 43 ) THEN
C
C                 RESTRUCTURATION D'UN QUADRANGLE NON STRUCTURE
C                 =============================================
                  CALL SUEX43( NTLXO, MCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )

               ELSE IF( NTYDFO .EQ. 45 ) THEN
C
C                 OPERATEURS LOGIQUES R2 SUR DES 2 TRIANGULATIONS PLANES
C                 OPERATEURS LOGIQUES R3 SUR DES 2 TRIANGULATIONS FERMEES DE R3
C                 =============================================================
                  CALL SUEX45( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 46 ) THEN
C
C                  COUDRE DES ARETES SIMPLES POUR FERMER UNE SURFACE MAILLEE
C                 ==========================================================
                  CALL SUEX46( NUMOB,  NTLXO,  MCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ELSE IF( NTYDFO .EQ. 48 ) THEN
C
C                 SURFACE 48 = IMPOSER DES POINTS COMME PLUS PROCHE SOMMET
C                 ========================================================
                  CALL SUEX48( NTLXO,  MCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
ccc               ELSE IF( NTYDFO .EQ. 49 ) THEN
C
C                 DESTRUCTURATION D'UNE SURFACE STRUCTUREE
C                 ========================================
CCCccc               CALL SUEX49( NTLXO, MCN(MNDFO) ,
CCCccc     %                      NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
ccc                  CALL STNOST( NTLXO,
ccc     %                         NTNSEO, MNNSEO, IERR )
C
               ENDIF
C
            ELSE IF( NUTYA .EQ. 4 ) THEN
C
C              VOLUME
C              ==========================================================
               IF( (NTYDFO .EQ. 1) .OR. (NTYDFO.EQ.2) ) THEN
C
C                 VOLUME 1 = HEXAEDRE ALGEBRIQUE STRUCTURE GORDON,..
C                 ===============================================
C                 VOLUME 2 = HEXAEDRE STRUCTURE  GORDON+WINSLOW
C                 ===============================================
                  CALL VOEX01( NTYDFO,NTLXO, MCN(MNDFO),
     %                         RMCN(MNDFO), NTNSEO,MNNSEO,NTSOM ,
     %                         MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 3 ) THEN
C
C                 VOLUME 3 = TETRAEDRE STRUCTURE
C                 ==============================
                  CALL VOEX03( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 4 ) THEN
C
C                 VOLUME 4 = PENTAEDRE STRUCTURE
C                 ==============================
                  CALL VOEX04( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
CCC               ELSE IF( NTYDFO .EQ. 7 ) THEN
CCCC
CCCC                 VOLUME = TETRAEDRISATION FRONTALE
CCCC                 =================================
CCC                  CALL VOEX07( NUMOB,NTLXO,MCN(MNDFO), RMCN(MNDFO),
CCC     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
CCCC
               ELSE IF( NTYDFO .EQ. 8 ) THEN
C
C                 VOLUME = UN VOLUME D'UN VOLUME 9 => ERREUR
C                 ==========================================
                  NTNSEO = 0
                  MNNSEO = 0
                  NTSOM  = 0
                  MNSOM  = 0
                  IERR   = 3
C
               ELSE IF( NTYDFO .EQ. 9 ) THEN
C
C                 VOLUME = TETRAEDRISATION1 DE VOLUMES DE TYPE 8
C                 ==============================================
                  CALL VOEX09( NUMOB ,NTLXO  ,
     %                         MCN(MNDFO), RMCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 10 ) THEN

C                 DONNEE DES TMS XYZSOMMET et NSEF
C                 ================================
                  CALL TAMSOU( NTLXO, M )
                  NTSOM  = MCN( MNDFO + WUTSOV )
                  CALL LXTSOU( NTLXO, 'XYZSOMMET', NTSOM, MNSOM )
                  NTNSEO = MCN( MNDFO + WUTSSV )
                  CALL LXTSOU( NTLXO, 'NSEF', NTNSEO, MNNSEO )
                  IERR = 0

               ELSE IF( NTYDFO .EQ. 11 ) THEN
C
C                 VOLUME 11 = PARALLELIPIPEDE_RECTANGLE
C                 =====================================
                  CALL VOEX11( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 12 ) THEN
C
C                 VOLUME 12 = TETRAEDRE DROIT
C                 ===========================
                  CALL VOEX12( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 13 ) THEN
C
C                 VOLUME 13 = PENTAEDRE DROIT
C                 ===========================
                  CALL VOEX13( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 14 ) THEN
C
C                 VOLUME 14 = VOLUME GENERE PAR L'EXTRUSION D'UNE SURFACE
C                 =======================================================
                  CALL VOEX14( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 15 ) THEN
C
C                 VOLUME 15 = VOLUME GENERE PAR ROTATION D'UNE SURFACE
C                 ====================================================
                  CALL VOEX15( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 16 ) THEN
C
C                 VOLUME 16 = GROSSIR PAR COUCHES NORMALES UNE LIGNE
C                 ===================================================
                  CALL VOEX16( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
CCCC
CCC               ELSE IF( NTYDFO .EQ. 17 .OR. NTYDFO .EQ. 18 ) THEN
CCCC
CCCC                 VOLUME 17 18 = 1/8 OU 1 SPHERE
CCCC                 ===============================
CCC                  CALL VOEX18( NTLXO, MCN(MNDFO),RMCN(MNDFO) ,
CCC     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
CCCC
C
               ELSE IF( NTYDFO .EQ. 19 ) THEN
C
C                 VOLUME = TETRAEDRISATION2 DE VOLUMES DE TYPE 8
C                 ==============================================
                  CALL VOEX19( NUMOB ,NTLXO,
     %                         MCN(MNDFO), RMCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
CC
               ELSE IF( NTYDFO .EQ. 21 ) THEN
C
C                 VOLUME 21 = CYLINDRE
C                 ====================
                  CALL VOEX21( NTLXO, MCN(MNDFO),RMCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 22 ) THEN
C
C                 VOLUME 22 = 1/2 CYLINDRE
C                 ========================
                  CALL VOEX22( NTLXO, MCN(MNDFO),RMCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 23 .OR. NTYDFO .EQ. 24 ) THEN
C
C                 VOLUME 23 = CONE ET 24 = 1/2 CONE
C                 =================================
                  CALL VOEX24( NTLXO,  MCN(MNDFO),RMCN(MNDFO) ,
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ELSE IF( NTYDFO .EQ. 25 ) THEN
C
C                 VOLUME 25 = ARBRE a PARTIR de CERCLES a Z=Cte
C                 =============================================
                  CALL VOEX25( NUMOB,  MCN(MNDFO),RMCN(MNDFO) ,
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ELSE IF( NTYDFO .EQ. 26 ) THEN
C
C                 VOLUME 26 = M a J de la LONGUEUR des ARETES des TETRAEDRES
C                 ==========================================================
                  CALL VOEX26( NUMOB,  MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ELSE IF( NTYDFO .EQ. 30 ) THEN
C
C                 VOLUME 30 = AMELIORATION D'UNE TETRAEDRISATION
C                 ==============================================
                  CALL VOEX30( NUMOB,  NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ELSE IF( NTYDFO .EQ. 31 ) THEN
C
C                 VOLUME 31 = TETRAEDRISATION D'UN MAILLAGE DE
C                             TETRAEDRES-PYRAMIDES-PENTAEDRES-HEXAEDRES
C                 =====================================================
                  CALL VOEX31( NUMOB,  NTLXO,  MCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ELSE IF( NTYDFO .EQ. 33 ) THEN
C
C                 VOLUME 33 = HEXAEDRISATION D'UN CUBE OU SPHERE
C                 ==============================================
                  CALL VOEX33( NTLXO,  MCN(MNDFO), RMCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ELSE IF( NTYDFO .EQ. 35 ) THEN
C
C                 VOLUME 35 = EXTRACTION D'UN SOUS VOLUME D'UN VOLUME STRUCTURE
C                 =============================================================
                  CALL VOEX35( NTLXO, MCN(MNDFO),
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 36 ) THEN
C
C                 EXTRACTION D'UN SOUS-VOLUME D'UN VOLUME
C                 =======================================
                  CALL VOEX36( NTLXO, MCN(MNDFO) ,
     %                         NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
               ELSE IF( NTYDFO .EQ. 39 ) THEN
C
C                 VOLUME 39 = DIVISION EN 8 SOUS EF DES EF 3D
C                 =====================================================
                  CALL VOEX39( NTLXO,  MCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ELSE IF( NTYDFO .EQ. 48 ) THEN
C
C                 VOLUME 48 = IMPOSER DES POINTS COMME PLUS PROCHE SOMMET
C                 =======================================================
                  CALL VOEX48( NTLXO,  MCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
ccc               ELSE IF( NTYDFO .EQ. 49 ) THEN
C
C                 DESTRUCTURATION D'UN VOLUME STRUCTURE
C                 =====================================
CCCccc               CALL VOEX49( NTLXO, MCN(MNDFO) ,
CCCccc     %                      NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
ccc                  CALL STNOST( NTLXO,
ccc     %                         NTNSEO, MNNSEO, IERR )
C
               ELSE IF( NTYDFO .EQ. 53 ) THEN
C
C                 VOLUME 53 = SEPARATION DES MATERIAUX MAILLES
C                 ============================================
                  CALL VOEX53( NUMOB,  MCN(MNDFO),
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C                 ATTENTION: LE NOM DE L'OBJET EST REMPLACE PAR SON NOM.1
                  CALL NMOBNU( 'VOLUME', NUMOB, NOMOBJ )
C                 ICI .1 A ETE AJOUTE AU NOM DU VOLUME
C
               ELSE IF( NTYDFO .EQ. 54 ) THEN
C
C                 VOLUME 54 = SUPPRESSION DES VOLUMES EN UN SEUL VOLUME
C                       LE VOLUME MULTI-MATERIAUX DEVIENT MONO-MATERIAU
C                 =====================================================
                  CALL VOEX54( NUMOB,
     %                         NTNSEO, MNNSEO, NTSOM, MNSOM, IERR )
C
               ENDIF
C
            ENDIF
C
         ELSE IF( NTYDFO .EQ. 50 ) THEN
C
C           TRANSFORMATION DE CET OBJET
C           ===========================
C           LE NUMERO DE L'OBJET INITIAL
            NUOBJ1 = MCN( MNDFO + WULIIN )
C
C           COPIE SIMPLE DE 'XYZSOMMET',EVENTUELLEMENT 'NSEF'
C           DE L'OBJET INITIAL DANS L'OBJET FINAL
            CALL COPOBJ( NMTOB ,NUOBJ1,NUMOB,
     %                   NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C           LA TRANSFORMATION EST ASSUREE PLUS TARD PAR CALL TRPLSV
C           CF DESSOUS
C
         ELSE IF( NTYDFO .EQ. 51 .OR. NTYDFO .EQ. 52 ) THEN
C
C           UNION D'OBJETS ALORS IDENTIFICATION DES SOMMETS
C                                MISE A JOUR DANS LE MAILLAGE
C                                51 => PLSV     C0-CONTINUS
C                                52 => SURFACES G1-CONTINUES
C           =================================================
C           ATTENTION:IL FAUT PASSER MCN(MNDFO+WBLIUN) ET NON SA VALEUR
C                     CAR ELLE PEUT ETRE MODIFIEE SI UN NOM EST DOUBLE
            CALL UNPLSV( NTYDFO, NUTYA, NUMOB,
     %                   MCN(MNDFO+WBLIUN), MCN(MNDFO+WULIUN),
     %                   NTNSEO,MNNSEO ,
     %                   NTSOM ,MNSOM ,NTUNIO,MNUNIO,IERR )
            IF( IERR .NE. 0 ) GOTO 9999
C
         ELSE IF( NTYDFO .EQ. 61 ) THEN
C
C           6D VOLUME 61 = 6-CUBE de type DIFFERENCES FINIES
C           ================================================
            CALL VOEX61( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                   NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
         ELSE IF( NTYDFO .EQ. 62 ) THEN
C
C           6D VOLUME 62 = 6-CUBE PAR HOMOTHETIES D'UN 6-CUBE NOYAU
C           =======================================================
            CALL VOEX62( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                   NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
         ELSE IF( NTYDFO .EQ. 63 ) THEN
C
C           6D VOLUME 63 = 6-CUBE 3 DIFFERENCES FINIES et HOMOTHETIES
C           =========================================================
            CALL VOEX63( NTLXO, MCN(MNDFO),RMCN(MNDFO),
     %                   NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C
         ELSE IF( NTYDFO .EQ. 81 ) THEN
C
C           RENOMMER UN P L S V
C           ===================
            CALL PLSV81( NUTYA,NUMOB ,NTLXO, MCN(MNDFO) ,
     %                   NTNSEO,MNNSEO,NTSOM,MNSOM,IERR )
C           LE NOUVEAU NUMERO DU P L S V EST EN FAIT L'ANCIEN
            NTLXOB = NTLXO
C
         ELSE IF( NTYDFO .EQ. 82 ) THEN
C
C           TUER UN P L S V
C           ===============
            CALL PLSV82( NUTYA,NUMOB ,NTLXO, MCN(MNDFO) ,
     %                   IERR )
C           EN SORTIE NTLXO A ETE MIS A ZERO CAR LE LEXIQUE EST DETRUIT
            GOTO 9999
C
         ELSE
C
C           TYPE NON DEFINI DANS LE TABLEAU DEFINITION
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'MAILEX: '//MMTOB//'TYPE NON DEFINI POUR'
            ELSE
               KERR(1) = 'MAILEX: '//MMTOB//'TYPE NOT DEFINED FOR'
            ENDIF
            KERR(2) =  KNOM
            CALL LEREUR
            IERR = 2
            GOTO 9999
         ENDIF
C
C        LA TRANSFORMATION DES SOMMETS
C        -----------------------------
         IF( NUTRAN .GT. 1 .AND. IERR .EQ. 0 ) THEN
            CALL TRPLSV( NUTYA,NUMOB,NTSOM,MNSOM,IERR )
         ENDIF
         GOTO 10
C
C     ==================================================================
C     FIN DU TANT QUE PILE DES OBJETS NON VIDE  FAIRE ...
C     ==================================================================
      ENDIF
C
      IF( IERR .NE. 0 ) THEN
C        UNE ERREUR EST INTERVENUE
C        DESTRUCTION DES TMS XYZSOMMET, NSEF, DEFINITION S'ILS EXISTENT
         CALL LXTSOU( NTLXOB, 'XYZSOMMET', NTSOMM, MNSOMM )
         IF( NTSOMM .GT. 0 ) THEN
            CALL LXTSDS( NTLXOB, 'XYZSOMMET' )
         ENDIF
         CALL LXTSOU( NTLXOB, 'NSEF', NTNSEF, MNNSEF )
         IF( NTNSEF .GT. 0 ) THEN
            CALL LXTSDS( NTLXOB, 'NSEF' )
         ENDIF
         CALL LXTSOU( NTLXOB, 'DEFINITION', NTDF, MNDF )
         IF( NTDF .GT. 0 ) THEN
            CALL LXTSDS( NTLXOB, 'DEFINITION' )
         ENDIF
         RETURN
      ENDIF
C
C     LE NUMERO ET ADRESSE DU TABLEAU 'NSEF' DE L'OBJET
C     -------------------------------------------------
      IF( NUTYOB .GT. 1 ) THEN
         CALL LXTSOU( NTLXOB,'NSEF',NTNSEF,MNNSEF )
      ELSE
C        POINTS
         NTNSEF = 0
         MNNSEF = 0
      ENDIF
C
C     LE NUMERO ET ADRESSE DU TABLEAU 'XYZSOMMET' DE L'OBJET
C     ------------------------------------------------------
      I = 0
      CALL LXTSOU( NTLXOB,'XYZSOMMET',NTSOMM,MNSOMM )
C
C     MISE A JOUR DU CADRE MAXIMAL DES SOMMETS ACTUELS
      IF( MNSOMM .GT. 0 ) CALL MAJEXT( MNSOMM )
C
C     AFFICHAGE DU NOMBRE DE SOMMETS ET D'EF
      IF( NUTYOB .GE. 2 .AND. NUTYOB .LE. 4 .AND. LORBITE .EQ. 0 ) THEN
         K1 = NUDCNB( MMTOBJ )
         K2 = NUDCNB( NOMOBJ )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,19001) MMTOBJ(1:K1), NOMOBJ(1:K2), NTYDFO,
     %           MCN(MNSOMM+WNBSOM), MCN(MNSOMM+WNBTGS),
     %           MCN(MNNSEF+WBEFOB), MCN(MNNSEF+WBEFTG)
         ELSE
            WRITE(IMPRIM,29001) MMTOBJ(1:K1), NOMOBJ(1:K2), NTYDFO,
     %           MCN(MNSOMM+WNBSOM), MCN(MNSOMM+WNBTGS),
     %           MCN(MNNSEF+WBEFOB), MCN(MNNSEF+WBEFTG)
         ENDIF
      ENDIF
19001 FORMAT(1X,A,' ',A,' de TYPE',I3,' : ',I9,' Sommets,', I9,
     % ' Tangentes,', I9,' Elements Finis,', I9, ' EF avec Tangentes' )
29001 FORMAT(1X,A,' ',A,' of TYPE',I3,' : ',I9,' Vertices,', I9,
     % ' Tangents,', I9,' Finite Elements,', I9, ' FE with Tangents' )
C
      IF( NUTYOB .EQ. 3 .OR. NUTYOB .EQ. 4 ) THEN
C
         IF( LCRITR .NE. 0  .OR. MAFAIT .NE. 0 ) THEN

C           IMPRESSION DE LA QUALITE DU MAILLAGE DE LA SURFACE OU VOLUME
 9001       CALL IMPQUA( NUTYOB, NOMOBJ, MNNSEF, MNSOMM, NBEFMQ, QUAMIN,
     %                   SURVOLEF )

            IF( NUTYOB .EQ. 3 .AND. NBEFMQ .GT. 0 .AND.
     %          MAFAIT .NE. 0 .AND. I .EQ. 0 ) THEN
C
C              TENTATIVE DE TRANSFORMER LES QUADRANGLES AYANT 2 SOMMETS
C              CONSECUTIFS IDENTIQUES EN TRIANGLES
               CALL QDEGET( NUTYOB, NOMOBJ, MNNSEF, MNSOMM, NBEFMQ,
     %                      IERR )
               IF( IERR .GT. 0 ) GOTO 9999
               I = I + 1
C              REIMPRESSION APRES CORRECTION
               GOTO 9001
C
            ENDIF
         ENDIF
      ENDIF
C
C     SI L'OBJET EST UNE LIGNE OU UNE SURFACE AVEC UN TYPE
C     INCONNU DE FERMETURE ALORS SON TYPE EST MIS A JOUR
C     ----------------------------------------------------
      IF( NUTYOB .EQ. 2 .OR. NUTYOB .EQ. 3 ) THEN
C        LIGNE OU SURFACE
         IF( MNNSEF .GT. 0 ) THEN
C           EXISTENCE DU TABLEAU 'NSEF'
C           SON TYPE DE FERMETURE EST FORCE INCONNU => MISE A JOUR
            MCN( MNNSEF + WUTFMA ) = -1
cccC           SANS AFFICHAGE DES ARETES SIMPLES (DANS 1 SEULE FACE)
ccc            CALL OBJFER( NUTYOB, NUMOBJ, 0, I )
C           AVEC AFFICHAGE DES ARETES SIMPLES (DANS 1 SEULE FACE)
            CALL OBJFER( NUTYOB, NUMOBJ, 1, I )
         ENDIF
      ENDIF
C
C     LE TYPE DE L'OBJET TRAITE EST EFFACE
 9999 NUMTYPOBJ = 0
      NUMOBJLX  = 0

      IF( MNPILE .GT. 0 ) CALL TNMCDS( 'ENTIER',MXPILE, MNPILE )

CCC      WRITE(IMPRIM,*) (MCN(I),I=1,5)
CCC      WRITE(IMPRIM,*) 'sortie mailex avec ',nomobj

      RETURN
      END
