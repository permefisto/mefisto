      SUBROUTINE SUEX23( NUSFTO, NTLXSU, LADEFI, RADEFI,
     %                   NTNSTO, MNNSTO, NTSTTO, MNSTTO, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     GENERER LE MAILLAGE DE LA SURFACE FERMEE D'UN TORE
C -----     SOIT AVEC UNE LONGUEUR PAR DEFAUT DES ARETES
C           SOIT AVEC EMPLOI DE LA FONCTION TAILLE_IDEALE(X,Y,Z)
C                DEFINIE PAR L'UTILISATEUR
C
C ENTREES :
C --------
C NUSFTO  : NUMERO DE LA SURFACE DANS LE LEXIQUE DES SURFACES
C NTLXSU  : NUMERO DU TABLEAU TS DU LEXIQUE DES SURFACES
C LADEFI  : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
C RADEFI  : TABLEAU REEL   DE DEFINITION DE LA SURFACE
C           CF '~td/d/a_surface__definition'
C
C SORTIES :
C ---------
C NTNSTO  : NUMERO      DU TMS 'NSEF' DES NUMEROS DES ELEMENTS FINIS
C MNNSTO  : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES ELEMENTS FINIS
C           CF '~td/d/a___nsef'
C NTSTTO  : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSTTO  : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C           CF '~td/d/a___xyzsommet'
C IERR    : 0 SI PAS D'ERREUR
C         > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain LJLL UPMC & St Pierre du Perray Novembre 2011
C2345X7..............................................................012
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/darete.inc"
C
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN( 1 )
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
      INTRINSIC         ATAN, MAX
      CHARACTER*24      NOMSUR
C
      IERR  = 0
      MNSOM = 0
      MXSOM = 0
      MXEF  = 0
      MNEF  = 0
C
C     PARAMETRES DU MAILLAGE
C     ======================
C     variable RAPCTO 'rayon du petit cercle generateur' reel ;
      RAPCTO = RADEFI( WAPCTO )
C     variable RAGDTO 'rayon du grand cercle parcours'   reel ;
      RAGCTO = RADEFI( WAGCTO )
C     variable PTCETO 'nom du point centre du tore'    ^~>POINT ;
      NTCETO = LADEFI( WTCETO )
C     variable PTAXTO 'nom du point sur axe X du tore' ^~>POINT ;
      NTAXTO = LADEFI( WTAXTO )
C     variable PTAYTO 'nom du point sur axe Y du tore' ^~>POINT ;
      NTAYTO = LADEFI( WTAYTO )
C
C     VERIFICATION DES PARAMETRES
C     ===========================
      IF( RAPCTO .LE. 0.0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15),'(G15.7)') RAPCTO
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'RAYON INCORRECT PETIT CERCLE du TORE='
     %             // KERR(MXLGER)(1:15)
         ELSE
            KERR(1) = 'INCORRECT RADIUS of SMALL CIRCLE of TORE'
     %             // KERR(MXLGER)(1:15)
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
      IF( RAGCTO .LE. 0.0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15),'(G15.7)') RAGCTO
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'RAYON INCORRECT GRAND CERCLE du TORE='
     %              // KERR(MXLGER)(1:15)
         ELSE
            KERR(1) = 'INCORRECT RADIUS of GREAT CIRCLE of TORE'
     %             // KERR(MXLGER)(1:15)
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
      IF( RAPCTO .GE. RAGCTO-RAPCTO ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)( 1:15),'(G15.7)') RAPCTO
         WRITE(KERR(MXLGER)(16:30),'(G15.7)') RAGCTO
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TORE avec des RAYONS INCOHERENTS='
     %              // KERR(MXLGER)(1:30)
         ELSE
            KERR(1) = 'TORUS with INCORRECT RADIUS of CIRCLES='
     %             // KERR(MXLGER)(1:30)
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     EXISTE-T-IL UNE FONCTION TAILLE_IDEALE(X,Y,Z)
C                             ou EDGE_LENGTH(X,Y,Z)?
C     ----------------------------------------------
      NOFOTI = NOFOTIEL()
      IF( NOFOTI .LE. 0 ) THEN
C
C        NON: LA LONGUEUR PAR DEFAUT DES ARETES EST ELLE DEFINIE?
C        --------------------------------------------------------
         IF( DARETE .LE. 0D0 ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='LONGUEUR PAR DEFAUT DES ARETES NON INITIALISEE'
               KERR(2)='UTILISER L''OPTION 0 DU MENU DEBUT'
            ELSE
               KERR(1)='DEFAULT LENGTH of EDGES NOT INITIALIZED'
               KERR(2)='USE THE OPTION 0 of THE DEBUT MENU'
            ENDIF
            CALL LEREUR
            IERR = 3
            GOTO 9999
         ENDIF
      ENDIF
C
C     RECUPERATION DES 3 COORDONNEES DU CENTRE DU TORE
C     ================================================
      CALL LXNLOU( NTPOIN, NTCETO, NTLXSC, MN )
      CALL LXTSOU( NTLXSC, 'XYZSOMMET', NTSOSC, MNSOSC )
      MNCTTO = MNSOSC + WYZSOM
C
C     RECUPERATION DES 3 COORDONNEES DU POINT SUR AXE X DU TORE
C     =========================================================
      CALL LXNLOU( NTPOIN, NTAXTO, NTLXSC, MN )
      CALL LXTSOU( NTLXSC, 'XYZSOMMET', NTSOSC, MNSOSC )
      MNAXTO = MNSOSC + WYZSOM
C
C     RECUPERATION DES 3 COORDONNEES DU POINT SUR AXE Y DU TORE
C     =========================================================
      CALL LXNLOU( NTPOIN, NTAYTO, NTLXSC, MN )
      CALL LXTSOU( NTLXSC, 'XYZSOMMET', NTSOSC, MNSOSC )
      MNAYTO = MNSOSC + WYZSOM
C
C     CENTRE et POINT AXE X SONT ILS IDENTIQUES?
      CALL XYZIDE( RMCN(MNCTTO), RMCN(MNAXTO), IDENTQ )
      IF( IDENTQ .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'POINTS CONFONDUS CENTRE et AXE X du TORE'
         ELSE
            KERR(1) = 'IDENTICAL CENTRE and AXIS X POINT of TORUS'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     CENTRE et POINT AXE Y SONT ILS IDENTIQUES?
      CALL XYZIDE( RMCN(MNCTTO), RMCN(MNAYTO), IDENTQ )
      IF( IDENTQ .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'POINTS CONFONDUS CENTRE et AXE Y du TORE'
         ELSE
            KERR(1) = 'CENTRE and AXIS Y POINT of TORUS are IDENTICAL'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     LES POINTS des AXES SONT ILS IDENTIQUES?
      CALL XYZIDE( RMCN(MNAXTO), RMCN(MNAYTO), IDENTQ )
      IF( IDENTQ .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'POINTS CONFONDUS AXE X et Y du TORE'
         ELSE
            KERR(1) = 'AXIS X et Y of TORUS are IDENTICAL'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     LONGUEUR DU PETIT CERCLE ET NOMBRE D'ARETES
      R = REAL( ATAN(1D0) * 8D0 * RAPCTO / DARETE )
C     NOMBRE D'ARETES SUR LE PETIT CERCLE
      NBARPC = NINT( R )
      IF( NBARPC .LE. 2 ) NBARPC=3
C
C     LONGUEUR DU GRAND CERCLE ET NOMBRE D'ARETES
      R = REAL( ATAN(1D0) * 8D0 * RAGCTO / DARETE )
C     NOMBRE D'ARETES SUR LE PETIT CERCLE
      NBARGC = NINT( R )
      IF( NBARGC .LE. 2 ) NBARGC=3
C
C     NOMBRE DE SOMMETS ET EF A PARTIR DE L'ARETE PAR DEFAUT
      IF( NOFOTI .LE. 0 ) THEN
         MXSOM = NBARPC * NBARGC
         MXEF  = MXSOM
      ENDIF
C
C     TRIANGULATION DE LA SURFACE FERMEE DU TORE
C     ==========================================
      IF( NOFOTI .LE. 0 ) THEN
C
C        EMPLOI DE LA TAILLE DES ARETES PAR DEFAUT
C        -----------------------------------------
         CALL TNMCDC( 'REEL',   MXSOM*3, MNSOM )
         CALL TNMCDC( 'ENTIER', MXEF*4,  MNEF  )
         CALL TOREDA( RAPCTO, NBARPC, RAGCTO, NBARGC,
     %                RMCN(MNCTTO), RMCN(MNAXTO), RMCN(MNAYTO),
     %                NBSOM, RMCN(MNSOM), NBEF, MCN(MNEF),
     %                IERR )
      ELSE
C
C        EMPLOI DE LA FONCTION TAILLE_IDEALE(X,Y,Z)
C        ------------------------------------------
         CALL TORETI( NUSFTO, RAPCTO, RAGCTO,
     %                RMCN(MNCTTO), RMCN(MNAXTO), RMCN(MNAYTO),
     %                NTFASU, MNFASU, NTSOFA, MNSOFA, NOMSUR,
     %                IERR )
C        NOMBRE DE SOMMETS DU MAILLAGE
         NBSOM = MCN( MNSOFA + WNBSOM )
C        ADRESSE MCN DES XYZ DES SOMMETS
         MNSOM = MNSOFA + WYZSOM
C        NOMBRE D'EF DU MAILLAGE
         NBEF  = MCN( MNFASU + WBEFOB )
C        ADRESSE MCN DES NO DES SOMMETS DES EF
         MNEF  = MNFASU + WUSOEF
C
      ENDIF
C
      IF( IERR .NE. 0 ) THEN
C        PAS DE TMS XYZSOMMET ET NSEF
         NTNSTO=0
         MNNSTO=0
         NTSTTO=0
         MNSTTO=0
         GOTO 9999
      ENDIF
C
C     GENERATION DES XYZ DES SOMMETS
C     ==============================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET'
      CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'MOTS', WYZSOM+3*NBSOM )
      CALL LXTSOU( NTLXSU, 'XYZSOMMET', NTSTTO, MNSTTO )
C
C     GENERATION DES NUMEROS DES NOEUDS DES EF
C     ========================================
C     CONSTRUCTION DU TABLEAU 'NSEF'
      CALL LXTNDC( NTLXSU, 'NSEF', 'ENTIER', WUSOEF+4*NBEF )
      CALL LXTSOU( NTLXSU, 'NSEF', NTNSTO,   MNNSTO )
C
C     MISE A JOUR DES TABLEAUX
C     ========================
C     MISE A JOUR DU TABLEAU 'XYZSOMMET' DE CETTE SURFACE FERMEE
C     NBSOM 'nombre de sommets'
      MCN( MNSTTO + WNBSOM ) = NBSOM
C     LE NOMBRE DE TANGENTES
      MCN( MNSTTO + WNBTGS ) = 0
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSTTO + WBCOOR ) = 3
C
C     LES 3 COORDONNEES DES NBSOM SOMMETS
      CALL TRTATA( RMCN(MNSOM), RMCN(MNSTTO+WYZSOM), 3*NBSOM )
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSTTO) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSTTO + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     MISE A JOUR DU TABLEAU 'NSEF' DE CETTE SURFACE FERMEE
C     TYPE DE L'OBJET : SURFACE
      MCN( MNNSTO + WUTYOB ) = 3
C     LE TYPE FERME DE FERMETURE DU MAILLAGE
      MCN( MNNSTO + WUTFMA ) = 1
C     PAS DE TANGENTES STOCKEES
      MCN( MNNSTO + WBTGEF ) = 0
      MCN( MNNSTO + WBEFAP ) = 0
      MCN( MNNSTO + WBEFTG ) = 0
C     NUMERO DU TYPE DE MAILLAGE : NON STRUCTURE
      MCN( MNNSTO + WUTYMA ) = 0
C     NBSOEF 'nombre de sommets par sous-objets'
      MCN( MNNSTO + WBSOEF ) = 4
C     NBEFOB 'nombre de sous-objets de l''objet'
      MCN( MNNSTO + WBEFOB ) = NBEF
C
C     LES 4 NUMEROS DES SOMMETS DES NBEF EF
      CALL TRTATA( RMCN(MNEF), RMCN(MNNSTO+WUSOEF), 4*NBEF )
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNNSTO) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNNSTO + MOTVAR(6) ) = NONMTD ( '~>>>NSEF' )
C
      IF( NOFOTI .GT. 0 ) THEN
C        DESTRUCTION DE LA SURFACE DES PARAMETRES
         CALL LXLXOU( NTSURF, NOMSUR, NT1SF, MN1SF )
         IF( MN1SF .GT. 0 ) CALL LXTSDS( NTSURF, NOMSUR )
      ENDIF
C
C     DESTRUCTION DES TABLEAUX INUTILES
 9999 IF( MXSOM .GT. 0 ) CALL TNMCDS( 'REEL',   MXSOM*3, MNSOM )
      IF( MXEF  .GT. 0 ) CALL TNMCDS( 'ENTIER', MXEF*4,  MNEF  )
      RETURN
      END
