      SUBROUTINE NLSEIMPL( KNOMOB, MOREE2, NDIM,  NDPGST, MNXYZN,
     %                     NBTYEL, MNNPEF,
     %                     MNTPOB, MNTAUX,
     %                     NUMIOB, NUMAOB,
     %                     MNDOEL, IESOPO,
     %                     RELMIN, MNTHER, MNTAEL, MNX,
     %                     NORESO, MNLPLI, MNLPCO, NBRDKG,NCODSK,
     %                     MNUG00, RDeltaT,DTSTOC, TPSINI, TPSFIN,
     %                     NBVECT, MXVECT, NTDL,
     %                     NTVECT, MNVECT, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : PB NLSE: CALCULER L'ONDE COMPLEXE ET SES FLUX DANS UN DOMAINE
C ----- 1D ou 2D ou 3D OU AXISYMETRIQUE DE l'EQUATION NON LINEAIRE
C       DE SCHRODINGER (NLSE)
C       i Rho dU(t,X)/dt - Alfa LAPLACIEN U(t,X) + N(U(t,X)**2) U
C                        - i OmegaZ (x dU/dy - y dU/dx) = F
C       SELON  TESTNL=7 METHODE IMPLICITE D'EULER AVEC PAS DE TEMPS CONSTANT
C       MASSE, CONDUCTIVITE, ECHANGE INDEPENDANTS DU TEMPS ET ONDE
C       FORCE, FIXATION, COEFFICIENT DU TERME NON LINEAIRE PEUVENT
C       DEPENDRE DU TEMPS ET DES ITERATIONS DE POINT FIXE SONT NECESSAIRES
C       L'ONDE IMPOSEE (CONDITION DE DIRICHLET) EST TRAITEE PAR PENALISATION
C       AVEC ECHANGE=1/EPSILON et FORCE=ONDE/EPSILON ET 1/EPSILON=PENALI=1D30
C***********************************************************************
C   Integration en temps par EULER a PAS CONSTANT et IMPLICITE MODIFIE
C   n=0  Vn=v0     Wn=w0
C
C(3)m=0  Vn+1m=Vn  Wn+1m=Wn
C
C A L'ETAPE m+1, LE PROBLEME CONSISTE A TROUVER Un+1 SOLUTION DE
C   [-(-Alfa Laplac+N(V0**2+W0**2));  Rho/dt   ]{Vn+1m+1}
C(4)[  Rho/dt;     -Alfa Laplac +N(V0**2+W0**2)]{Wn+1m+1} =
C
C { -Fr(tn+1)  - (N(V0**2+W0**2)-N(Vn+1m**2+Wn+1m**2)) Vn+1m
C   +Rho/dt Wn +  OmegaZ (x d/dy-y d/dx) Wn+1m }
C {  Fi(tn+1)  + (N(V0**2+W0**2)-N(Vn+1m**2+Wn+1m**2)) Wn+1m 
C   +Rho/dt Vn +  OmegaZ (x d/dy-y d/dx) Vn+1m }

C  Si ||Vn+1m+1-Vn+1m||>Eps||Vn+1m+1|| ou ||Wn+1m+1-Wn+1m||>Eps||Wn+1m+1||
C  Alors m=m+1  Aller en (4)
C  Sinon Vn+1=Vn+1m+1  Wn+1=Wn+1m+1
C
C  Si TPSINI+(n+1)dt < TPSFIN Alors n=n+1; Aller en (3)
C  Fin
C
C L'IDEE EST DE CALCULER UNE SEULE FOIS LA MATRICE GLOBALE PROFIL 2n x 2n
C   [-Alfa Laplac+N(V0**2+W0**2); Rho/dt                     ]
C(  [Rho/dt;                     -Alfa Laplac +N(V0**2+W0**2)]
C et de la FACTORISER L D tL de CROUT POUR RESOUDRE LE SYSTEME LINEAIRE
C
C A L'ETAPE 0 IL FAUT RECUPERER UG0 = VG0 + i WG0
C***********************************************************************
C
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET A TRAITER
C MOREE2 : NOMBRE DE MOTS   D'UNE VARIABLE REELLE DOUBLE PRECISION
C NDIM   : DIMENSION DES COORDONNEES DES POINTS ( 1 OU 2 OU 3 )
C NDPGST : CODE D'IDENTIFICATION DES SOMMETS POINTS NOEUDS
C MNXYZN : ADRESSE MCN DU TABLEAU XYZNOEUD DE L'OBJET KNOMOB
C
C NBTYEL : NOMBRE DE TYPES D'EF DU MAILLAGE DE CET OBJET
C MNNPEF : ADRESSE MCN DU TABLEAU DES ADRESSES MCN DES TMS NPEF"TYPE EF
C NBTTEF : NOMBRE TOTAL D'EF DU MAILLAGE
C MNTPOB : ADRESSE MCN DU TABLEAU POINTEUR SUR LES TABLEAUX POBA DES EF
C NBDLMX : NOMBRE MAXIMAL DE DEGRES DE LIBERTE (TEMPERATURES) D'UN EF
C MOAUX  : NOMBRE DE OTS AUXILIAIRES NECESSAIRES AU CALCUL DES EF
C MNTAUX : ADRESSE MCN DU TABLEAU AUXILIAIRE POUR LES EF
C
C NUMIOB : NUMERO MINIMAL DU PLSV DANS LA DEFINITION DE L'OBJET
C NUMAOB : NUMERO MAXIMAL DU PLSV DANS LA DEFINITION DE L'OBJET
C NBOBIN : NOMBRE DE VOLUMES en 3D, SURFACES en 2D DE L'OBJET
C MNOBIN : ADRESSE MCN DU DEBUT DU TABLEAU NUOBIN (PARTIE DE TOPOLOGIE)
C NBOBCL : NOMBRE DE PLS en 3D, PL en 2D DE L'OBJET
C MNOBCL : ADRESSE MCN DU DEBUT DU TABLEAU NUOBCL (PARTIE DE TOPOLOGIE)
C MXTYEL : NOMBRE MAXIMAL DE TYPES D'EF (7)
C MXDOEL : NOMBRE DE MOTS DECLARES DES TABLEAUX D'ADRESSE DANS MNDOEL
C MNDOEL : LES 4 ADRESSES MCN DES TABLEAUX DES ADRESSES DES
C          TABLEAUX DECRIVANT LES DONNEES THERMIQUES DE L'OBJET COMPLET
C
C IEMAST : NOMBRE DE TMS MASSE            DES SV DE L'OBJET RETROUVES
C IECHMA : NOMBRE DE TMS CHALEURMASSIQUE  DES SV DE L'OBJET RETROUVES
C IECOND : NOMBRE DE TMS CONDUCTIVITE     DES SV DE L'OBJET RETROUVES
C IEDILA : NOMBRE DE TMS DILATATION       DES SV DE L'OBJET RETROUVES
C IEVIFL : NOMBRE DE TMS VITESSEFLUIDE    DES SV DE L'OBJET RETROUVES
C IECOET : NOMBRE DE TMS COEFTEMPERATURE  DES SV DE L'OBJET RETROUVES
C IESOIN : NOMBRE DE TMS FORCE "INTERNE"  DES SV DE L'OBJET RETROUVES
C
C IECONT : NOMBRE DE TMS FIXATION            DES PLS DE L'OBJET RETROUVES
C IEECHA : NOMBRE DE TMS ECHANGE             DES PLS DE L'OBJET RETROUVES
C IESOCL : NOMBRE DE TMS FORCE "AUX LIMITES" DES PLS DE L'OBJET RETROUVES
C IESOPO : NOMBRE DE TMS FORCE               DES P   DE L'OBJET RETROUVES
C PENALI : COEFFICIENT DE PENALISATION DES TEMPERATURES FIXEES
C          ICI PENALI VAUT 1D30 POUR LE PRENDRE EN COMPTE
C RELMIN : PLUS PETIT REEL SERVANT DE MARQUEUR DE NON UTILISATION
C
C MNTHER : 128 REELS DOUBLE PRECISION POUR LA MATRICE DE CONDUCTIVITE
C MOTAEL : NOMBRE DE MOTS DECLARES DU TABLEAU DES TABLEAUX ELEMENTAIRES
C MNTAEL : ADRESSE MCN DES TABLEAUX ELEMENTAIRES
C MNX    : ADRESSE MCN DU TABLEAU DES 3 COORDONNEES DES NOEUDS=POINTS DE L'EF
C
C NORESO : CODE RESOLUTION ET STOCKAGE DU SYSTEME LINEAIRE
C          1 FACTORISATION DE CROUT SUR UNE MATRICE PROFIL
C MNLPLP : ADRESSE MCN DES POINTEURS SUR LES COEFFICIENTS DIAGONAUX DE LA MATRIC
C MNLPLI : ADRESSE MCN DES POINTEURS SUR LES COEFFICIENTS DIAGONAUX DE LA MATRIC
C MNLPCO : ADRESSE MCN DU NUMERO DE COLONNE DES COEFFICIENTS DE LA MATRICE MORSE
C NIVEAU : NIVEAU DE FACTORISATION INCOMPLETE DE LA MATRICE DE PRECONDITIONNEMEN
C          A CHOISIR PARMI 0 1 2
C NBRDKG : NOMBRE DE REELS DOUBLE PRECISION DES MATRICES M ET K
C NCODSK : 1 SI MATRICE DE CONDUCTIVITE SYMETRIQUE
C          0 SI MATRICE DE CONDUCTIVITE DIAGONALE
C         -1 SI MATRICE DE CONDUCTIVITE NON SYMETRIQUE
C MNUG00 : ADRESSE MCN DE U0(NBNOMA,2) TEMPERATURE A L'INSTANT INITIAL
C
C RDeltaT: PAS CONSTANT DU TEMPS REEL SIMPLE PRECISION
C DTSTOC : PAS CONSTANT DU TEMPS ENTRE 2 STOCKAGES DU VECTEUR"ONDENLSE
C TPSINI : TEMPS INITIAL DU CALCUL
C TPSFIN : TEMPS FINAL   DU CALCUL
C NBVECT : NUMERO DU DERNIER VECTEUR TEMPERATURE STOCKE
C MXVECT : NOMBRE DE VECTEUR"ONDENLSE A STOCKER
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE DE L'OBJET =
C          2 FOIS LE NOMBRE DE NOEUDS DU MAILLAGE POUR LA PARTIE REELLE
C          ET IMAGINAIRE DE L'ONDE
C
C SORTIES:
C --------
C NTVECT : NUMERO      DU TMS VECTEUR"ONDENLSE DE L'OBJET
C MNVECT : ADRESSE MCN DU TMS VECTEUR"ONDENLSE DE L'OBJET
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
C MODIFS : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray  Aout 2011
C23456---------------------------------------------------------------012
      PARAMETER         (MXTYEL=7)
      PARAMETER         (ITERMX=10)
      PARAMETER         (ITERMXKG=8)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/donthe.inc"
      include"./incl/donele.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__erreurth.inc"
      include"./incl/a___arete.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___conductivite.inc"
      include"./incl/a___dilatation.inc"
      include"./incl/a___source.inc"
      include"./incl/a___contact.inc"
      include"./incl/a___fixation.inc"
      include"./incl/a___force.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___profil.inc"
      include"./incl/a___morse.inc"
      include"./incl/a___fluxpt.inc"
      include"./incl/a___dtemperature.inc"
      include"./incl/a___temperinit.inc"
      include"./incl/msvaau.inc"
      include"./incl/ponoel.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
      include"./incl/xyzext.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      INTEGER           NUMIOB(4), NUMAOB(4), MNDOEL(4)
      CHARACTER*(*)     KNOMOB
      CHARACTER*160     KNOM
      CHARACTER*16      TEXT
C
      DOUBLE PRECISION, allocatable, dimension(:) :: KG
      INTEGER           IERKGALLOC
      INTRINSIC         ALLOCATED
C
      REAL              RDeltaT, DTSTOC, TSTOC, TPSINI, TPSFIN
      DOUBLE PRECISION   DeltaT, DINFO, D, Testm,
     %                  MODUMX0,   MODUMX,
     %                  WCAMIN(2), WCAMAX(2), WEXMIN(2), WEXMAX(2),
     %                  ERRMAX(2), ERRMXR, ERRMXI,
     %                  CPU, CPU0, CPUVG, CPUPKG, CPUFAC, CPUITE
      INTEGER           NOFOWE(2)
      INTRINSIC         REAL
C
      DOUBLE PRECISION  RELMIN, PENALI, TGV, UNSTGV
C     RELMIN : PLUS PETIT REEL SERVANT DE MARQUEUR DE NON UTILISATION
      RELMIN = -1D28
C
C     TEMPS DE CONSTRUCTION DE KG et FACTORISATION DE CROUT
      CPUVG  = 0D0
      CPUPKG = 0D0
      CPUFAC = 0D0
      CPUITE = 0D0
C
Cccc     PENALISATION DE LA CONDITION DE DIRICHLET PAR FOURIER
Cccc     ECHANGE=PENALI=1/EPSILON ET FORCE=TEMPERATURE/EPSILON
ccc      PENALI = 1D30
      PENALI = 0D0
C
C     COEFFICIENT DE PRISE EN COMPTE DES DL FIXES DE DIRICHLET
      TGV = 1D30
C
C     PAS DE TEMPS EN DOUBLE PRECISION
      DeltaT = RDeltaT
C
C     L'ADRESSE DE L'ONDE A L'ITERATION 0 DU POINT FIXE
      MNTHET  = MNUG00
      MNTHET0 = MNUG00
C
C     ADRESSES MCN DES TABLEAUX DE TRAVAIL
      MNBG   = 0
      MNNFNX = 0
      MNVFNX = 0
      MNNDLX = 0
      MNVDLX = 0
      MNERRT = 0
      MNMODU = 0
      MNPIL3 = 0
      MNSOLE = 0
      MOVALS = 0
      MNVALS = 0
      MNFBAS = 0
      MNCOPO = 0
      MNXYZS = 0

C     NBCOOR = NOMBRE DE COORDONNEES DES NOEUDS
      NBCOOR = MCN(MNXYZN+WBCOON )
C
C     NOMBRE DE NOEUDS DU MAILLAGE
      NBNOMA  = NTDL / 2
      NBNOEMA = NBNOMA

      IF( INTERA .GE. 1 ) THEN

C        CONSTRUCTION DES DONNEES POUR LE TRACE DU MODULE DE L'ONDE A CHAQUE m
C        =====================================================================
C        DECLARATION DU TABLEAU DU MODULE AUX NOEUDS DE L'ONDE
         CALL TNMCDC( 'REEL2', NBNOMA, MNMODU )

C        TRACE DU MODULE D'UNE ONDE COMPLEXE AUX NOEUDS DU MAILLAGE
         MODECO = 8
C        LA PALETTE ARC EN CIEL
         CALL PALCDE( 14 )

C        CADRE MAXIMAL DE L'OBJET
         CALL MIMXPT( NBCOOR, NBNOMA, MCN(MNXYZN+WYZNOE), COOEXT )

C        PARAMETRES DE LA VISEE
         NOTYVI  = 0
         LORBITE = 0
         NETAXE  = 0

         IF( NDIM .EQ. 2 ) THEN
C           LA VISEE EN 2D
            CALL VISEE0

C           PARAMETRES DES TABLEAUX POUR LE TRACE DES ZONES DE COULEURS
C           MXSOUI : MAXIMUM D INTERVALLES DU SEGMENT UNITE
C           MXSOUI ** 2 SOUS-TRIANGLES CREES AU PLUS DANS UN TRIANGLE
            MXSOUI = 4
C           MXPIL3 : MAXIMUM DE SOUS-TRIANGLES GENERES DANS L ELEMENT REFERENCE
            MXPIL3 = 2 * MXSOUI * MXSOUI + 16
C           MAXDLE : MAXIMUM DE D.L. D UNE INTERPOLATION (CF SP INTERP)
            MAXDLE = 30
C
C           NPIL3 (3,MXPIL3) PILE DES 3 SOMMETS DE CHAQUE SOUS-TRIANGLE
            MOPIL3 = 3 * MXPIL3
            CALL TNMCDC( 'ENTIER', MOPIL3, MNPIL3 )
C
C           SOLEL (MXNOEL) SOLUTION AUX NOEUDS DE L ELEMENT DU MAILLAGE
C           FBASE (MAXDLE) VALEUR DES FONCTIONS DE BASE EN UN POINT
C           COPOE (MXPOEL,3) COORDONNEES DES POINTS DE L ELEMENT
            MOSOLE = MXNOEL+MAXDLE+MXPOEL*3
            CALL TNMCDC( 'REEL2', MOSOLE, MNSOLE )
            MNFBAS = MNSOLE + MOREE2 * MXNOEL
            MNCOPO = MNFBAS + MOREE2 * MAXDLE
C
C           MXSOMM : MAXIMUM DE SOMMETS DES SOUS-TRIANGLES DE L EF REFERENCE
            MXSOMM = ( MXSOUI + 1 ) ** 2 + 16
C           VALST (3,MXSOMM) COORDONNEES ET VALEUR DE LA SOLUTION
C           XYZSOM(3,MXSOMM) COORDONNEES DES SOMMETS DES SOUS-TRIANGLES
C           AUX SOMMETS DES SOUS-TRIANGLES (DANS R2 OU R3)
            IF( MNVALS .GT. 0 ) CALL TNMCDS( 'REEL', MOVALS, MNVALS )
            MOVALS = 6 * MXSOMM
            CALL TNMCDC( 'REEL', MOVALS, MNVALS )
            MNXYZS = MNVALS + 3 * MXSOMM
         ENDIF
      ENDIF
C
C     DEBUT DU TEMPS DE CALCUL
      CPU0 = DINFO( 'CPU' )
C
C     AFFICHAGE DE L'ONDE AU TEMPS INITIAL UG00(NBNOMA,2) PAR COMPOSANTES
C     -------------------------------------------------------------------
      CALL AFNLSE( 5,   1, MNXYZN, NBNOMA, 1, MCN(MNUG00),
     %             WCAMIN, WCAMAX, NOFOWE, WEXMIN, WEXMAX, ERRMAX )
C
C     CALCUL du Max |UG00(NOEUD)| POUR DETECTER UNE EXPLOSION DE L'ONDE
C     -----------------------------------------------------------------
      MN0 = ( MNUG00 - 1 ) / MOREE2
      MODUMX0 = DMCN(MN0+1)
      DO I=1,NBNOMA
         MN  = MN0 + I
         D = SQRT( DMCN(MN)**2 + DMCN(MN+NBNOMA)**2 )
         IF( D .GT. MODUMX0 ) MODUMX0 = D
      ENDDO
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'Max |Onde(t=',TEMPS,',NOEUD)|=',MODUMX0
      ELSE
         WRITE(IMPRIM,*) 'Max |Wave(t=',TEMPS,',NODE)|=',MODUMX0
      ENDIF
C
C     VECTEUR GLOBAL BG(NTDL) AUXILIAIRE
C     VECTEUR GLOBAL UGtn  DE L'ONDE A L'INSTANT tn      EST DECLARE(NBNOMA,2)
C     VECTEUR GLOBAL UGtm  DE L'ONDE A L'INSTANT tn+1m   EST DECLARE(2,NBNOMA)
C     VECTEUR GLOBAL UGtm1 DE L'ONDE A L'INSTANT tn+1m+1 EST DECLARE(2,NBNOMA)
C     VECTEUR GLOBAL UGtmC DE L'ONDE A L'INSTANT tn+1m   EST DECLARE(NBNOMA,2)
C     ========================================================================
      CALL TNMCDC( 'REEL2', 5*NTDL, MNBG )
      MNUGtm  = MNBG    + MOREE2 * NTDL
      MNUGtm1 = MNUGtm  + MOREE2 * NTDL
      MNUGtmC = MNUGtm1 + MOREE2 * NTDL
      MNUGtn  = MNUGtmC + MOREE2 * NTDL

C     SOLUTION UG00(NBNOMA,2) AU TEMPS t0 => UGtn(NBNOMA,2) AU TEMPS t0
      CALL TRTATD( MCN(MNUG00), MCN(MNUGtn), NTDL )
      MNTHETn = MNUGtn
C     SOLUTION UG00(NBNOMA,2) AU TEMPS t0 => UGtmC(NBNOMA,2) AU TEMPS t0
      CALL TRTATD( MCN(MNUG00), MCN(MNUGtmC), NTDL )
C
C     SORTIE DE THED1T LE VECTEUR UG00 EST DONNE PAR COMPOSANTES UG00(NBNOMA,2)
C     LES DL DE LA PARTIE REELLE PUIS CEUX DE LA PARTIE IMAGINAIRE
C     UG1 DOIT ETRE RANGE PAR NOEUDS EN CHAQUE NOEUD, PARTIE RELLE PUIS IMAGINAIRE
C     UG00(NBNOMA,2) => UGtm(2,NBNOMA)
C     ============================================================================
      CALL DLCPND( 2, NBNOMA, MCN(MNUG00), MCN(MNUGtm) )
C
C     TABLEAU MC: TEMPS + Testm Max |U|
C                       + eventuellement ERREUR REEL et IMAG
C               A CHAQUE PAS DE TEMPS CALCULE et VALEUR de m
C     ------------------------------------------------------
      IF( NOFOWE(1) .GT. 0 .AND. NOFOWE(2) .GT. 0 ) THEN
C        STOKAGE: Testm, Max|U|, ERREUR PR, ERREUR PI
         NBVERR = 2 + 2
      ELSE
C        STOKAGE: Testm, Max|U| seulement
         NBVERR = 2
      ENDIF
C     1+NBVERR POUR STOCKER LE TEMPS,
C     INFO COMPLEMENTAIRE DANS VECTEUR"TESTM_ERREUR
      MOERRT = (1+NBVERR) * MXVECT * 2
      CALL TNMCDC( 'REEL', MOERRT, MNERRT )
C
C     BILAN DE L'INSTANT INITIAL HORS MATRICE GLOBAL KG
C     =================================================
C     LE PROCHAIN TEMPS POUR STOCKER LE VECTEUR"ONDENLSE
      TSTOC  = TPSINI + DTSTOC
      MNTEMP = MNVECT + WECTEU + MOREE2 * NTDL * (NBVECT-1)
      MNTIME = MNVECT + WECTEU + MOREE2 * NTDL * MXVECT - 1
C     LE NOMBRE DE DL D'UN VECTEUR
      MCN( MNVECT + WBCOVE ) = NTDL
C     LE NOMBRE DE VECTEURS STOCKES
      MCN( MNVECT + WBVECT ) = NBVECT
C     AU DEPART PAS DE TEMPS INDIQUES
      MCN( MNVECT + WBCPIN ) = 0
C
C     DECLARATION DE LA MATRICE GLOBALE [KG]
C     ======================================
      IERKGALLOC = 1
      IF( NBRDKG .LE. 0 ) THEN
C         PLACE MEMOIRE INSUFFISANTE
          NBLGRC(NRERR) = 3
          WRITE(KERR(MXLGER-1)(1:25),'(G25.0)') DMAT
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1) ='ERREUR: PLACE MEMOIRE INSUFFISANTE'
             KERR(2) = KERR(MXLGER-1)(1:25) //
     %               ' MOTS NECESSAIRES pour [M] [K] [R]'
             KERR(3) = 'REDUIRE le MAILLAGE'
          ELSE
             KERR(1) ='ERROR: NOT ENOUGH MEMORY'
             KERR(2) = KERR(MXLGER-1)(1:25) //
     %               ' NECESSARY WORDS to store [M] [K] [R]'
             KERR(3) ='REDUCE the MESH'
          ENDIF
          CALL LEREUR
          IF( INTERA .LE. 1 ) CALL ARRET( 100 )
          GOTO 9999
      ENDIF
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10290) NBRDKG, NBRDKG/NTDL
      ELSE
         WRITE(IMPRIM,20290) NBRDKG, NBRDKG/NTDL
      ENDIF
10290 FORMAT('UNE MATRICE PROFIL DE',I15,
     %' REELS DOUBLE PRECISION avec 1/2 LARGEUR DE BANDE MOYENNE =',I9)
20290 FORMAT('ONE SKYLINE MATRIX of',I15,' DOUBLE REALS with an HALF WID
     %TH AVERAGE=',I9)
C
C     ALLOCATION DYNAMIQUE EN FORTRAN 90 DE LA MATRICE PROFIL KG
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) 'ALLOCATION DEMAND  of',NBRDKG,
     %                ' DOUBLE PRECISION of [KG] MATRIX'
      ALLOCATE ( KG(1:NBRDKG), STAT=IERKGALLOC )
      IF( IERKGALLOC .NE. 0 ) THEN
       WRITE(IMPRIM,*)'ALLOCATION ERROR   of',NBRDKG,
     %                ' DOUBLE PRECISION of [KG] MATRIX'
         IERR = IERKGALLOC
         GOTO 9999
      ENDIF
      WRITE(IMPRIM,*) 'ALLOCATION CORRECT of',NBRDKG,
     %                ' DOUBLE PRECISION of [KG] MATRIX'
      WRITE(IMPRIM,*)
C
C     TEMPS CALCUL DES RESERVATIONS DE TABLEAUX
      CPU   = DINFO( 'CPU' )
      CPUVG = CPU - CPU0
C
C     ===============================================================
C     CONSTRUCTION DE LA MATRICE [KG] SYMETRIQUE NON DIAGONALE PROFIL
C     [-(-Alfa Laplac+N(V0**2+W0**2)); Rho/dt                     ]
C     [  Rho/dt;                      -Alfa Laplac +N(V0**2+W0**2)]
C     ===============================================================
      NBPADT = 0
      NBLOCm0= 0
      NBCABG = 0
      NBSMKG = 0
      ITERKG = 0
      NTERKG = 0
C
C     RETOUR ICI EN CAS DE MODIFICATION DU PAS DE TEMPS
C     DEBUT DU TEMPS DE CONSTRUCTION DE LA MATRICE KG[2NBNOMA,2NBNOMA]
 50   CPU0   = DINFO( 'CPU' )
      ID00tm = 1
      CALL NLSEIMKG( DeltaT, PENALI, NDIM,   MNXYZN,
     %               NBTYEL, MNNPEF, NDPGST,
     %               MNTPOB, MNTAUX, NUMIOB, NUMAOB, MNDOEL,
     %               MNTHER, MNTAEL, MNX,
     %               NORESO, NCODSK, NBRDKG, MNLPLI, MNLPCO, KG,
     %               NBPTAF, IERR  )
      IF( IERR .NE. 0 ) GOTO 9999

      call affvect( 'KG PROFIL SANS CL', 20    , KG )
      call afl1ve(  'KG PROFIL SANS CL', NBRDKG, KG )
C
C     PRISE EN COMPTE DES DL FIXES (NUMEROTATION DES DL PAR NOEUDS)
C     ============================
C     CONSTRUCTION DE LA LISTE ET DES VALEURS DES DL FIXES DIRICHLET
      CALL DLFXLI( 1,      NBNOMA, NDIM,   NBTYEL, MNNPEF, NDPGST,
     %             MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %             NBRDLX, MONDLX, MNNDLX, MNVDLX, IERR )
C
C     PRISE EN COMPTE DES DL FIXES SUR LA MATRICE KG
C     ----------------------------------------------
      CALL DLFXMG( NBRDLX, MCN(MNNDLX), TGV, NTDL, MCN(MNLPLI), KG )
C
C     TEMPS DE CALCUL DE CONSTRUCTION DE KG
      CPU    = DINFO( 'CPU' )
      CPUPKG = CPUPKG + CPU - CPU0
      CPU0   = CPU
C
C     FACTORISATION DE CROUT DE [KG]=L D tL
C     =====================================
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10050) NTERKG+1
      ELSE
         WRITE(IMPRIM,20050) NTERKG+1
      ENDIF
10050 FORMAT('Factorisation',I5,' de Crout: [KG]=[L] [D] t[L]')
20050 FORMAT('Crout Factorization',I5,' : [KG]=[L] [D] t[L]')
C
      CALL CRMC1D( MCN(MNLPLI), KG, NTDL, 1E-10, 1, KG, IERR )
C     IERR : 0 SI   AUCUN     PIVOT<EPS
C            1 SI AU MOINS UN PIVOT<EPS
C     ON ENTEND PAR PIVOT (A0(IDIAGONAL)-SA)/A0(IDIAGONAL)

      call affvect( 'CROUT L D tL PROFIL', 20    , KG )
      call afl1ve(  'CROUT L D tL PROFIL', NBRDKG, KG )

      IF( IERR .NE. 0 ) THEN
C        MATRICE NON INVERSIBLE
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: MATRICE NON INVERSIBLE'
            KERR(2) = 'REVOYEZ LES CONDITIONS AUX LIMITES'
         ELSE
            KERR(1) = 'ERROR: NON INVERSIBLE MATRIX'
            KERR(2) = 'SEE AGAIN BOUNDARY CONDITIONS'
         ENDIF
         CALL LEREUR
         IF( INTERA .LE. 1 ) CALL ARRET( 100 )
         IERR = 7
         GOTO 9999
      ENDIF
C
C     LE TEMPS INITIAL EST TEMPSINI
C
C     TEMPS CUMULE DE CALCUL DE LA FACTORISATION DE CROUT DE KG
      CPU    = DINFO( 'CPU' )
      CPUFAC = CPUFAC + CPU - CPU0
      CPU0   = CPU

C     ##################################################################
C     ##                                                              ##
C     ##  LA BOUCLE EN TEMPS AVEC DES PAS DE TEMPS CONSTANTS = DeltaT ##
C     ##                                                              ##
C     ##################################################################
C
C     MNUGtm CONTIENT UG0(2,NBNOMA)  L'ONDE A L'INSTANT INITIAL
C
C                   [ K(-Alfa +N(UG0**2)), M(Rho/DeltaT)       ]
C     KG = L D tL =
C                   [ M(Rho/DeltaT)        K(-Alfa +N(UG0**2)) ]
C
C     ICI LE NOUVEAU TEMPS tn+1 OU SE FAIT LE CALCUL
C     ==============================================
 100  TEMPS = REAL( TEMPS + DeltaT )
      ITERm = 0

C     MODIFICATION DU VECTEUR UGtmC POUR EVITER QUE BG=0 QUAND Un+1m=U00
      IF( ID00tm .EQ. 1 ) THEN
         MN = (MNUGtmC-1) / MOREE2
         DO I=1,NTDL
            DMCN(MN+I) = 1.01D0 * DMCN(MN+I)
         ENDDO
         ID00tm = 0
      ENDIF
C
C     ICI LES ITERATIONS DE POINT FIXE SONT MISES EN OEUVRE
C     =====================================================
C     UNE NOUVELLE ITERATION DE POINT FIXE
 110  ITERm = ITERm + 1
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10100) TEMPS, DeltaT, ITERm, NBCABG
      ELSE
         WRITE(IMPRIM,20100) TEMPS, DeltaT, ITERm, NBCABG
      ENDIF
10100 FORMAT(/'Au TEMPS',G14.6,' PAS de TEMPS=',G15.7,' ITERm=',I3,
     %'  Nb ITERm=',I10,'  -> CALCUL de l''ONDE COMPLEXE'/140('*'))
20100 FORMAT(/'At TIME',G14.6,' TIME STEP=',G15.7,' ITERm=',I3,
     %'  ITERm Nb=',I10,'  -> COMPUTATION of the COMPLEX WAVE'/140('*'))
C
C     CALCUL DU SECOND MEMBRE ELEMENTAIRE BG=Fn+1,m+1(2,NBNOMA) A L'INSTANT tn+1
C     --------------------------------------------------------------------------
C     TESTNL=7: SCHEMA IMPLICITE
C     BG1 =-Fr(tn+1) - (N(V0**2+W0**2)-N(Vn+1m**2+Wn+1m**2)) Vn+1m
C                       +Rho/dt Wn +OmegaZ ( x dWn+1m/dy - y dWn+1m/dx )
C     BG2 = Fi(tn+1) + (N(V0**2+W0**2)-N(Vn+1m**2+Wn+1m**2)) Wn+1m
C                    + Rho/dt Vn + OmegaZ ( x dVn+1m/dy - y dVn+1m/dx )
      MNTHET = MNUGtm
      CALL NLSEIMBG( DeltaT, PENALI, NDIM,   MNXYZN, NBNOMA,
     %               NBTYEL, MNNPEF, NDPGST,
     %               MNTPOB, MNTAUX, NUMIOB, NUMAOB, MNDOEL,
     %               MNX,    MCN(MNUGtn), MCN(MNUGtmC),
     %               MCN(MNBG), MCN(MNBG), IERR  )
      IF( IERR .NE. 0 ) GOTO 9999
C
C     AJOUT DES EVENTUELLES FORCES PONCTUELLES DANS BG(2,NBNOMA)
C     ----------------------------------------------------------
      IF( IESOPO .GT. 0 .AND. NDIM .GT. 1 ) THEN
C
C        CONSTRUCTION NODL(NBNOMA,2) ET VALEUR DES FORCES PONCTUELLES
         CALL THSNFX( NTDL,   NBTYEL, MNNPEF, NDPGST,
     %                MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %                NBFNFX, MONFNX, MNNFNX, MOVFNX, MNVFNX )
C
C        ASSEMBLAGE DES FORCES PONCTUELLES DANS BG(2,NBNOMA)
         CALL ASFONOEUD( 2, NBNOMA, 1,
     %                   NBFNFX, MCN(MNNFNX), MCN(MNVFNX), RELMIN,
     %                   MCN(MNBG) )
C
      ENDIF
C
C     PRISE EN COMPTE DES DL FIXES SUR LE VECTEUR GLOBAL BG(2,NBNOMA)
C     ---------------------------------------------------------------
      CALL DLFXVG( NBRDLX, MCN(MNNDLX), MCN(MNVDLX), TGV, NTDL,
     %             MCN(MNBG) )
C
C     RESOLUTION DU SYSTEME FACTORISE PAR CROUT KG UG1 = L D tL UG1 = BG
C     -------------------------------------------------------------------
      CALL DRCRPR( NTDL, NCODSK, MCN(MNLPLI), KG, MCN(MNBG), 3,
     %             MCN(MNUGtm1), IERR )
C
C     NETTOYAGE DES PETITES VALEURS
      UNSTGV = 1D0 / TGV
      MN1    = ( MNUGtm1 - 1 ) / MOREE2
      DO I=1,NTDL
         D = DMCN(MN1+I)
         IF( ABS(D) .LT. UNSTGV )  THEN
            DMCN(MN1+I)=0D0
         ENDIF
      ENDDO
C
C     CONSTRUCTION A PARTIR DE UGtm1(2,NBNOMA) DE UGtmC(NBNOMA,2)
C     -----------------------------------------------------------
      CALL DLNDCP( 2, NBNOMA, MCN(MNUGtm1), MCN(MNUGtmC) )
C     TERME NON LINEAIRE CALCULE AVEC CETTE NOUVELLE SOLUTION
      MNTHET = MNUGtmC
C
C     AFFICHAGE DE L'ONDE AU TEMPS ET POUR L'ITERATION DU POINT FIXE
C     AFFICHAGE DE L'ERREUR SI SOLUTION EXACTE CONNUE
C     --------------------------------------------------------------
      CALL AFNLSE( 5,   1, MNXYZN, NBNOMA, 1, MCN(MNUGtmC),
     %             WCAMIN, WCAMAX, NOFOWE, WEXMIN, WEXMAX, ERRMAX )
C
C     AFFICHER LES NORMES ABSOLUES ET RELATIVES DE
C     UGtm=UGn+1m(2,NBNOMA), UGtm1=UGn+1m+1(2,NBNOMA)
C     ET CALCULER LA VALEUR POUR TESTER LA FIN DES ITERATIONS m
C     ---------------------------------------------------------
      CALL NLSEITER( 1, TEMPS, ITERm, NBNOMA,
     %            MCN(MNUGtmC), MCN(MNUGtm), MCN(MNUGtmC), MCN(MNUGtm1),
     %            MODUMX, Testm )
C     MODUMX=Max | UG(tn+1,m+1,N) | aux NOEUDS N du MAILLAGE
C     Testm =Som||UGm+1|(N)-|UGm|(N)|/Som|UGm+1|(N) aux NOEUDS N du MAILLAGE

C     ERRMXR = MAX |RP EXACT-RP COMPUTED|(Node)/(MAX RP(Node)-MIN RP(Node))
C     ERRMXI = MAX |IP EXACT-IP COMPUTED|(Node)/(MAX IP(Node)-MIN IP(Node))
C     ---------------------------------------------------------------------
      IF( NOFOWE(1) .GT. 0 .AND. NOFOWE(2) .GT. 0 ) THEN
         ERRMXR = ERRMAX(1) / (WEXMAX(1)-WEXMIN(1))
         ERRMXI = ERRMAX(2) / (WEXMAX(2)-WEXMIN(2))
      ELSE
         ERRMXR = 0D0
         ERRMXI = 0D0
      ENDIF
C
C     STOCKAGE TEMPS + Testm + ModuleMax + ERREUR Max DANS TMC MNERRT
C     ---------------------------------------------------------------
      IF( (1+NBVERR)*NBCABG .GE. MOERRT ) THEN
C        AUGMENTATION DU TABLEAU
         CALL TNMCAU( 'REEL', MOERRT, 2*MOERRT, MOERRT, MNERRT )
         MOERRT = 2 * MOERRT
      ENDIF
C
      MN = MNERRT + (1+NBVERR) * NBCABG
C     LE TEMPS ACTUEL du CALCUL tn+1
      RMCN( MN ) = TEMPS
C     TEST D'ARRET DE L'ITERATION m =  Som| |Um+1|-|Um| | / Som |Um+1|
      RMCN( MN+1 ) = REAL( Testm )
C     Max|U(Node)|
      RMCN( MN+2 ) = REAL( MODUMX )
C
      IF( NBVERR .GT. 2 ) THEN
C        ERREUR RELATIVE PARTIE REELLE
         RMCN( MN+3 ) = REAL( ERRMXR )
C        ERREUR RELATIVE PARTIE IMAGINAIRE
         RMCN( MN+4 ) = REAL( ERRMXI )
      ENDIF
C
C     UN CALCUL ITERm DE PLUS EFFECTUE
      NBCABG = NBCABG + 1
      NBSMKG = NBSMKG + 1
C
C     TRACE DU MODULE DE L'ONDE A L'ITERATION m
C     -----------------------------------------
      IF( INTERA .GE. 1 ) THEN

C        L'ECRAN PIXELS EST EFFACE (PAS DE SCINTILLEMENT)
         CALL EFFACEMEMPX

C        CONSTRUCTION DU MODULE AUX NBNOMA NOEUDS DE L'ONDE
         MN  = ( MNMODU  - 1 ) / MOREE2
         MN1 = ( MNUGtmC - 1 ) / MOREE2
         DO I = 1, NBNOMA
            DMCN(MN+I) = SQRT( DMCN(MN1+I) ** 2
     %                       + DMCN(MN1+NBNOMA+I) ** 2 )
         ENDDO

         IF( NDIM .EQ. 2 ) THEN

C           TRACE DES AXES 2D
            CALL TRAXE2
C
C           BOUCLE SUR LES DIFFERENTS TYPES D'ELEMENTS FINIS DU MAILLAGE
            DO I = 0, NBTYEL-1
C
C              L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF DE CE TYPE D'EF
               MNELE = MCN( MNNPEF + I )
C
C              LE NUMERO DU TYPE DES ELEMENTS FINIS
               NUTYEL = MCN( MNELE + WUTYEL )
C
C              LE NOMBRE DE TELS ELEMENTS
               NBELEM = MCN( MNELE + WBELEM )
C
C              LES CARACTERISTIQUES DE L'ELEMENT FINI
               CALL ELTYCA( NUTYEL )
C
C              L'ADRESSE MCN DU TABLEAU 'NPEF' POUR CE TYPE D'EF
               MNPGEL = MNELE + WUNDEL
               IF( NDPGST .GE. 2 ) THEN
                  MNPGEL = MNPGEL + MCN(MNELE+WBELEM)*MCN(MNELE+WBNDEL)
               ENDIF

C              TRACE 2D EFFECTIF DES ZONES ISO-MODULE DE L'ONDE
               NBNOEU = MCN(MNXYZN+WNBNOE)
               CALL TRZON2( 0.0, REAL(MODUMX), NDIM, 1, 1, 1,
     %                      NBNOMA, MCN(MNMODU),
     %                      NUTYEL, NBELEM,
     %                      NBNOE,  MCN(MNELE+WUNDEL),
     %                      NBPOE,  MCN(MNPGEL),
     %                      NBCOOR, NBNOEU, MCN(MNXYZN+WYZNOE),
     %                      NBNOEU, MCN(MNXYZN+WYZNOE),
     %                      MXSOMM, MCN(MNSOLE), MCN(MNCOPO),
     %                      MXPIL3, MCN(MNPIL3),
     %                      MCN(MNVALS), MCN(MNFBAS), MCN(MNXYZS) )
            ENDDO

C        ELSE 

C           NDIM = 3   A FAIRE

         ENDIF
C
C        LE TRACE DE LA LEGENDE: COULEURS => VALEURS
         CALL LEGCOULSO( 0.0, REAL(MODUMX) )

C        TRACE DE LA 2-EME LIGNE DU TITRE DU TRACE
         CALL TIT2LG( KNOMOB, MODECO )
C
C        DEFINITION DE LA 3-EME LIGNE DU TITRE
         IF( LANGAG .EQ. 0 ) THEN
            KNOM = 'MODULE de l''ONDE: Temps= '
         ELSE
            KNOM = 'WAVE MAGNITUDE: Time= '
         ENDIF

C        LE TEMPS
         WRITE( TEXT, '(G14.6)' ) TEMPS
C        SUPPRESSION DES BLANCS DE DEBUT ET INTERMEDIAIRES
         CALL TEXTSB( TEXT, L )
         I = NUDCNB( KNOM )
         KNOM(I+1:I+L) = TEXT(1:L)

C        LE PAS DE TEMPS
         I = NUDCNB( KNOM )
         IF( LANGAG .EQ. 0 ) THEN
            KNOM(I+1:I+11) = ' Pas Temps='
         ELSE
            KNOM(I+1:I+11) = ' Time Step='
         ENDIF
         WRITE( TEXT, '(G14.6)' ) DeltaT
         CALL TEXTSB( TEXT, L )
         I = NUDCNB( KNOM )
         KNOM(I+1:I+L) = TEXT(1:L)

C        LE NOMBRE DE PAS DE TEMPS CALCULES
         I = NUDCNB( KNOM )
         IF( LANGAG .EQ. 0 ) THEN
            KNOM(I+1:I+14) = ' Nb Pas Temps='
         ELSE
            KNOM(I+1:I+14) = ' Time Step NB='
         ENDIF
         WRITE( TEXT, '(I9)' ) NBPADT
         CALL TEXTSB( TEXT, L )
         I = NUDCNB( KNOM )
         KNOM(I+1:I+L) = TEXT(1:L)

C        LE NOMBRE D'ITERATIONS m a ce TEMPS
         I = NUDCNB( KNOM )
         KNOM(I+1:I+7) = ' ITERm='
         WRITE( TEXT, '(I9)' ) ITERm
         CALL TEXTSB( TEXT, L )
         I = NUDCNB( KNOM )
         KNOM(I+1:I+L) = TEXT(1:L)

C        NBCABG LE NOMBRE TOTAL D'ITERATIONS m 
         I = NUDCNB( KNOM )
         KNOM(I+1:I+13) = ' Total ITERm='
         WRITE( TEXT, '(I9)' ) NBCABG
         CALL TEXTSB( TEXT, L )
         I = NUDCNB( KNOM )
         KNOM(I+1:I+L) = TEXT(1:L)

C        LE MAX DU MODULE EN UN NOEUD
         I = NUDCNB( KNOM )
         KNOM(I+1:I+11) = ' Max|U(X)|='
         WRITE( TEXT, '(G14.6)' ) MODUMX
         CALL TEXTSB( TEXT, L )
         I = NUDCNB( KNOM )
         KNOM(I+1:I+L) = TEXT(1:L)

C        Testm pour la convergence des iterations m
         I = NUDCNB( KNOM )
         KNOM(I+1:I+7) = ' Testm='
         WRITE( TEXT, '(G14.6)' ) Testm
         CALL TEXTSB( TEXT, L )
         I = NUDCNB( KNOM )
         KNOM(I+1:I+L) = TEXT(1:L)

C        RETOUR AUX PARAMETRES INITIAUX
         CALL XVEPAISSEUR( 1 )
         CALL XVTYPETRAIT( LIGCON )
C
C        TRACE DU TITRE
         CALL TRFINS( KNOM )
         print *,KNOM

      ENDIF
C
C     TEMPS DE CALCUL DES ITERATIONS DE POINT FIXE ET EN TEMPS
      CPU    = DINFO( 'CPU' )
      CPUITE = CPUITE + CPU - CPU0
      CPU0   = CPU
C
C     CONTROLE DES RESULTATS DE L'ITERATION m
C     =======================================
      IF( MODUMX .GT. 1000D0 * MODUMX0 ) THEN
C
C        TROP GRAND ACCROISSEMENT DU MAX EN UN NOEUD DE L'ONDE => EXPLOSION
C        ------------------------------------------------------------------
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10131) TEMPS
         ELSE
            WRITE(IMPRIM,20131) TEMPS
         ENDIF
10131    FORMAT(/100('*')/'EXPLOSION DETECTEE au TEMPS', G14.6/100('*'))
20131    FORMAT(/100('*')/'DETECTION of WAVE EXPLOSION at TIME', G14.6/
     %          100('*'))
C        SAUVEGARDE DES RESULTATS ACTUELS
         GOTO 9900
      ENDIF
C
 150  IF( NBCABG .GT. 4 .AND. Testm .GT. 1D-2 ) THEN
C
C        LE TEST m EST TROP GRAND => REDUIRE LE PAS DE TEMPS DeltaT
C        EST IL ENCORE POSSIBLE DE DIVISER LE PAS DE TEMPS PAR 2?
C        ----------------------------------------------------------
         IF( ITERKG .GE. ITERMXKG ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='STOP apres    DIVISIONS du PAS DE TEMPS'
               KERR(2)='Reduire le PAS de TEMPS INITIAL TROP GRAND?'
            ELSE
               KERR(1)='STOP after    DIVISIONS of the TIME STEP'
               KERR(2)='Reduce the INITIAL TIME STEP TOO GREAT?'
            ENDIF
            WRITE(KERR(1)(12:13),'(I2)') ITERMXKG
            CALL LEREUR
            IERR = 29
C           MISE DANS LE VECTEUR"ONDENLSE DES VECTEURS DEJA STOCKES
            GOTO 9900
         ENDIF
C
C        OUI: RECONSTRUCTION DE KG avec un NOUVEAU PAS DE TEMPS/2
C        ========================================================
         ITERKG = ITERKG + 1
         NTERKG = NTERKG + 1
C        Retour aux RESULTATS du dernier TEMPS CALCULE
         TEMPS  = REAL( TEMPS - DeltaT )
C        LE NOUVEAU TEMPS INITIAL EST tn
         TEMPSINI = TEMPS
C        LA SOLUTION INITIALE EST REMPLACEE PAR LA DERNIERE SOLUTION CALCULEE
C        CONSTRUCTION DE UG00(NBNOMA,2) A PARTIR DE UG0(2,NBNOMA)
         CALL DLNDCP( 2, NBNOMA, MCN(MNUGtm), MCN(MNUG00) )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10200) DeltaT, DeltaT/2D0
         ELSE
            WRITE(IMPRIM,20200) DeltaT, DeltaT/2D0
         ENDIF
10200    FORMAT(/'Changement du PAS de TEMPS',G14.6,
     %   ' par le NOUVEAU PAS de TEMPS=',G14.6,
     %   '  + RECONSTRUCTION [KG]')
20200    FORMAT(/'Changement of the TIME STEP=',G14.6,
     %   ' to the NEW TIME STEP=',G14.6,
     %   '  + [KG] is CONSTRUCT AGAIN')
C        NOUVEAU PAS DE TEMPS
         DeltaT = DeltaT / 2D0
         PasTemps = REAL( DeltaT )
         NBSMKG = 0
         GOTO 50
C
      ENDIF
C
      IF( Testm .GT. 5D-4 ) THEN
C
C        NON CONVERGENCE DU POINT FIXE => UNE ITERATION DE PLUS A FAIRE
C        --------------------------------------------------------------
         IF( ITERm .GE. ITERMX ) THEN
            NBLGRC(NRERR) = 3
            WRITE(KERR(5)(1:6),'(I6)') ITERMX
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: PAS de POINT FIXE ATTEINT apres'
               KERR(3) = 'PAS INITIAL de TEMPS TROP GRAND a REDUIRE?'
            ELSE
               KERR(1) = 'ERROR: FIX POINT NOT CONVERGED after'
               KERR(3) = 'INITIAL TIME STEP TOO GREAT to REDUCE?'
            ENDIF
            KERR(2) = KERR(5)(1:6) // ' ITERATIONS'
            CALL LEREUR
            Testm = 1D0
            GOTO 150
ccc            IERR = 29
cccC           MISE DANS LE VECTEUR"ONDENLSE DES VECTEURS DEJA STOCKES
ccc            GOTO 9900
         ENDIF
C
C        LE NOUVEAU VECTEUR ONDE DE L'ITERATION m-1 EST UG1 CALCULE
         MN      = MNUGtm
         MNUGtm  = MNUGtm1
         MNUGtm1 = MN
         GOTO 110
C
      ENDIF
C
C     CONVERGENCE DES ITERATIONS DE POINT FIXE
C     ----------------------------------------
C     UN PAS DE TEMPS CALCULE EN PLUS
      NBPADT = NBPADT + 1
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10140) TEMPS, ITERm
      ELSE
         WRITE(IMPRIM,20140) TEMPS, ITERm
      ENDIF
C
10140 FORMAT('Au TEMPS',G14.6,' ITER m=',I3,'  CONVERGENCE des ITERATION
     %S m')
20140 FORMAT('At TIME', G14.6,' ITER m=',I3,'  m ITERATION CONVERGENCE')

c     LE VECTEUR DES DL DE L'ONDE ACTUELLE EST IL A STOCKER?
C     ------------------------------------------------------
      IF( TEMPS .GE. TSTOC*0.9999 .AND. NBVECT .LT. MXVECT ) THEN
C
C        OUI: STOCKAGE DE L'ONDE BG=UG1(NBNOMA,2) A CET INSTANT TEMPS
         MNTEMP = MNTEMP + NTDL * MOREE2
         CALL TRTATD( MCN(MNUGtmC), MCN(MNTEMP), NTDL )
C        LE NOMBRE DE VECTEURS TEMPERATURE STOCKES
         NBVECT = NBVECT + 1
C        LE TEMPS DE STOCKAGE
         RMCN(MNTIME+NBVECT) = TEMPS
C
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10150) TEMPS,NBVECT,NTDL
         ELSE
            WRITE(IMPRIM,20150) TEMPS,NBVECT,NTDL
         ENDIF
C
10150 FORMAT('Au TEMPS',G14.6,' STOCKAGE DE l''ONDE',I5,
     %' de',I9,' DEGRES de LIBERTE')
20150 FORMAT('At TIME',G14.6,' STORAGE of the WAVE',I5,
     %' of',I9,' DEGREES of FREEDOM')
C
C        LE PROCHAIN TEMPS DE STOCKAGE
         TSTOC = TSTOC + DTSTOC

      ENDIF
C
C     MISE A JOUR DE MNUGtm ET MNUGtm1
C     LA DERNIERE SOLUTION CALCULEE EST MAINTENANT A L'ADRESSE MNUGtm
      MN      = MNUGtm
      MNUGtm  = MNUGtm1
      MNUGtm1 = MN
C
C     ETUDE EN TEMPS TERMINEE?
C     ------------------------
      IF( TEMPS + DeltaT .LT. TPSFIN*1.00001 ) THEN
C
C        NON: La SOLUTION UG0C(NBNOMA,2) AU TEMPS tn+1
C             EST COPIEE dans UGtn(NBNOMA,2) AU NOUVEAU TEMPS tn
         CALL TRTATD( MCN(MNUGtmC), MCN(MNUGtn), NTDL )
C
         IF( (Testm .LT. 5D-4 .AND. ITERm .LE. 1) .OR.
     %        MOD(NBSMKG,100) .EQ. 0 ) THEN
C
C           POUR ACCELERER RECONSTRUCTION DE KG avec un PAS DE TEMPS*1.5
C           et SOLUTION INITIALE REMPLACEE PAR LA DERNIERE SOLUTION CALCULEE
C           ================================================================
            ITERKG = ITERKG - 1
            NTERKG = NTERKG + 1
C           LE NOUVEAU TEMPS INITIAL EST tn
            TEMPSINI = TEMPS
C           SOLUTION INITIALE REMPLACEE PAR LA DERNIERE SOLUTION CALCULEE
C           CONSTRUCTION DE UG00(NBNOMA,2) A PARTIR DE UG0(2,NBNOMA)
            CALL DLNDCP( 2, NBNOMA, MCN(MNUGtm), MCN(MNUG00) )
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10200) DeltaT, DeltaT*1.5D0
            ELSE
               WRITE(IMPRIM,20200) DeltaT, DeltaT*1.5D0
            ENDIF
C           NOUVEAU PAS DE TEMPS
            DeltaT = DeltaT * 1.5D0
            PasTemps = REAL( DeltaT )
            NBSMKG = 0
            GOTO 50
C
         ENDIF
C
         NBLOCm = NBPADT / 100
         IF( NBPADT .LE. 5 .OR. NBLOCm .GT. NBLOCm0 ) THEN
C
C           RECONSTRUCTION DE KG avec U0=(Un+1,m+1) ET MEME PAS DE TEMPS
C           ============================================================
C           LE NOUVEAU TEMPS INITIAL EST tn+1
            TEMPSINI = TEMPS
C           SOLUTION INITIALE REMPLACEE PAR LA DERNIERE SOLUTION CALCULEE
            CALL TRTATD( MCN(MNUGtm1), MCN(MNUG00), NTDL )
C
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10201) TEMPS, DeltaT
            ELSE
               WRITE(IMPRIM,20201) TEMPS, DeltaT
            ENDIF
10201       FORMAT(/'RECONSTRUCTION de [KG] avec U0=Un+1,m+1 au temps=',
     %               G14.6,' avec un PAS de TEMPS=',G14.6)
20201       FORMAT(/'[KG] is CONSTRUCT AGAIN with U0=Un+1,m+1 at time='
     %               G14.6,' and a TIME STEP=',G14.6)
C           PAS DE TEMPS INCHANGE
            NBLOCm0 = NBLOCm
            GOTO 50
C
         ENDIF
C
         GOTO 100

      ENDIF
C
C    ##############################################################
C    ##                                                          ##
C    ##                FIN DE LA BOUCLE EN TEMPS                 ##
C    ##                                                          ##
C    ##############################################################
C
ccc 9900 IF( RMCN(MNTIME+NBVECT) .NE. TEMPS ) THEN
C     STOCKAGE DE L'ONDE BG=UG1(NBNOMA,2) AU DERNIER INSTANT TEMPS
 9900 IF( NBVECT .LT. MXVECT ) THEN
         MNTEMP = MNTEMP + NTDL * MOREE2
C        LE NOMBRE DE VECTEURS TEMPERATURE STOCKES
         NBVECT = NBVECT + 1
      ELSE
C        LE VECTEUR SOLUTION MXVECT
         NBVECT = MXVECT
      ENDIF
      CALL TRTATD( MCN(MNUGtmC), MCN(MNTEMP), NTDL )
C     LE TEMPS DE STOCKAGE
      RMCN(MNTIME+NBVECT) = TEMPS
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10150) TEMPS, NBVECT, NTDL
      ELSE
         WRITE(IMPRIM,20150) TEMPS, NBVECT, NTDL
      ENDIF
C
C     MISE A JOUR DU TMS 'VECTEUR"ONDENLSE'
C     =====================================
      MCN( MNVECT + WBCOVE ) = NTDL
      MCN( MNVECT + WBVECT ) = NBVECT
      MCN( MNVECT + WBCPIN ) = NBVECT
      IF( NBVECT .LT. MXVECT ) THEN
C        LE TMS EST RACOURCI
         L  = MNVECT + WECTEU + NTDL * MXVECT * MOREE2 - 1
         L1 = MNVECT + WECTEU + NTDL * NBVECT * MOREE2 - 1
         DO I=1,NBVECT
            RMCN(L1+I) = RMCN(L+I)
         ENDDO
         CALL TAMSRA( NTVECT, WECTEU+NTDL*NBVECT*MOREE2+NBVECT )
      ENDIF
C     LA DATE
      CALL ECDATE( MCN(MNVECT) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNVECT + MOREE2 ) = NONMTD( '~>>>VECTEUR' )

C     L'ADRESSE -1 DU PREMIER TEMPS STOCKE DERRRIERE LES VECTEURS TEMPERATURE
      MNTIME = MNVECT + WECTEU + NTDL * NBVECT * MOREE2 -1

C     CONSTRUCTION DE VECTEUR"TESTM_ERREUR A PARTIR DE L'EVENTUEL ANCIEN
C     AVEC AJOUT DU TMC MNERRT ACTUEL
C     ==================================================================
      CALL TITSMXERR( KNOMOB, MOREE2, TPSINI, NBVERR, NBCABG, MNERRT )
C
C     DESTRUCTION DES TMC DEVENUS INUTILES
C     ====================================
 9999 IF( IERKGALLOC .EQ. 0 ) DEALLOCATE( KG )
      IF( MNBG   .GT. 0 ) CALL TNMCDS( 'REEL2',  NTDL*5,  MNBG   )
      IF( MNNFNX .GT. 0 ) CALL TNMCDS( 'ENTIER', MONFNX,  MNNFNX )
      IF( MNVFNX .GT. 0 ) CALL TNMCDS( 'REEL2',  MOVFNX,  MNVFNX )
      IF( MNNDLX .GT. 0 ) CALL TNMCDS( 'ENTIER', MONDLX,  MNNDLX )
      IF( MNVDLX .GT. 0 ) CALL TNMCDS( 'REEL2',  MONDLX,  MNVDLX )
      IF( MNLPLI .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL+1,  MNLPLI )
      IF( MNERRT .GT. 0 ) CALL TNMCDS( 'REEL',   MOERRT,  MNERRT )
      IF( MNPIL3 .GT. 0 ) CALL TNMCDS( 'ENTIER', MOPIL3,  MNPIL3 )
      IF( MNSOLE .GT. 0 ) CALL TNMCDS( 'REEL2',  MOSOLE,  MNSOLE )
      IF( MNVALS .GT. 0 ) CALL TNMCDS( 'REEL',   MOVALS,  MNVALS )
      IF( MNMODU .GT. 0 ) CALL TNMCDS( 'REEL2',  NBNOMA,  MNMODU )
C
C     BILAN SUR LA PLACE MEMOIRE CENTRALE OCCUPEE PAR LES MATRICES ...
C     ================================================================
C     CROUT PROFIL
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,19002) NTDL + MOREE2 * NBRDKG
      ELSE
         WRITE(IMPRIM,29002) NTDL + MOREE2 * NBRDKG
      ENDIF
19002 FORMAT(/'NLSEIMPL: STOCKAGE MATRICE PROFIL=',I15,' MOTS'/ )
29002 FORMAT(/'NLSEIMPL: SKYLINE MATRIX STORAGE=',I15,' MEMORY WORDS')
C
C     AFFICHAGE DES TEMPS CALCUL
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,19000) CPUVG, CPUPKG, CPUFAC, CPUITE
      ELSE
         WRITE(IMPRIM,29000) CPUVG, CPUPKG, CPUFAC, CPUITE
      ENDIF
19000 FORMAT(/
     %'CPU DECLARATION DES VECTEURSVG=',F12.2,' SECONDES CPU'/
     %'CPU FORMATION DE LA MATRICE KG=',F12.2,' SECONDES CPU'/
     %'CPU FACTORISATION   MATRICE KG=',F12.2,' SECONDES CPU'/
     %'CPU ITERATIONS TEMPS + PT FIXE=',F12.2,' SECONDES CPU')
29000 FORMAT(/
     %'CPU of VECTOR VG DECLARATION  =',F12.2,' CPU SECONDS'/
     %'CPU of MATRIX KG CONSTRUCTION =',F12.2,' CPU SECONDS'/
     %'CPU of MATRIX KG FACTORIZATION=',F12.2,' CPU SECONDS'/
     %'CPU of TIME+FIX PT ITERATIONS =',F12.2,' CPU SECONDS')
C
      RETURN
      END
