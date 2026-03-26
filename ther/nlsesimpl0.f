      SUBROUTINE NLSESIMPL0( KNOMOB, MOREE2, NDIM,   NDPGST, MNXYZN,
     %                       NBTYEL, MNNPEF, MNTPOB, MNTAUX,
     %                       NUMIOB, NUMAOB, MNDOEL, IESOPO,
     %                       RELMIN, MNTHER, MNTAEL, MNX,
     %                       NORESO, MNLPLI, MNLPCO, NBRDKG, NCODSK,
     %                       MNUG00, RDeltaT,DTSTOC, TPSINI, TPSFIN,
     %                       NBVECT, MXVECT, NTDL,
     %                       NTVECT, MNVECT, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : NLSE: CALCULER L'ONDE COMPLEXE ET SES FLUX DANS UN DOMAINE
C -----       2D ou 3D OU AXISYMETRIQUE DE L'EQUATION NON LINEAIRE
C       COMPLEXE DE SCHRODINGER (NLSE) 
C       i Rho dU(t,X)/dt - Alfa LAPLACIEN U(t,X) + N(U(t,X)**2) U
C                        - i OmegaZ (x dU/dy - y dU/dx) = F
C       SELON TESTNL=6  UNE METHODE SEMI-IMPLICITE D'EULER MODIFIEE
C       AVEC PAS DE TEMPS CONSTANT 
C       MASSE, CONDUCTIVITE, ECHANGE INDEPENDANTS DU TEMPS ET ONDE
C       FORCE, FIXATION, COEFFICIENT DU TERME NON LINEAIRE PEUVENT
C       DEPENDRE DU TEMPS ET DES ITERATIONS DE POINT FIXE SONT NECESSAIRES
C       L'ONDE IMPOSEE (CONDITION DE DIRICHLET) EST TRAITEE PAR PENALISATION
C       AVEC ECHANGE=1/EPSILON et FORCE=ONDE/EPSILON ET 1/EPSILON=PENALI=1D30
C       ET DOIT ETRE AUX MEMES NOEUDS POUR LES PARTIES REELLE ET IMAGINAIRE
C       MAIS AVEC DES POSSIBLES VALEURS DIFFERENTES SI NECESSAIRE
C***********************************************************************
C   n=0  V=v0      W=w0
C
C(1)m=0  Vn+1m=Vn  Wn+1m=Wn
C
C(2)[N(V0**2+W0**2) -Alfa LAPLACIEN] Vn+1m+1 = Fr(tn+1)
C   + Rho/dt (Wn+1m-Wn) - OmegaZ ( x dWn+1m/dy - y dWn+1m/dx )
C   + ( N(V0**2+W0**2)-N(Vn+1m**2+Wn+1m**2) ) Vn+1m         
C
C   [N(V0**2+W0**2) -Alfa LAPLACIEN] Wn+1m+1 = Fi(tn+1)
C   - Rho/dt (Vn+1m-Vn) + OmegaZ ( x dVn+1m/dy - y dVn+1m/dx )
C   + ( N(V0**2-W0**2)-N(Vn+1m**2+Wn+1m**2) ) Wn+1m
C                       
C    Si ||Vn+1m+1-Vn+1m||>Eps||Vn+1m+1|| ou ||Wn+1m+1-Wn+1m||>Eps||Wn+1m+1||
C    Alors m=m+1  Aller en (2)
C    Sinon Vn+1=Vn+1m+1  Wn+1=Wn+1m+1
C
C   Si TPSINI+(n+1)dt < TPSFIN Alors n=n+1; Aller en (1)
C   Fin
C
C Ajout de N(V0**2+W0**2) pour avoir une matrice definie
C meme en cas de condition de NEUMANN nulle sur la frontiere
C
C L'IDEE EST DE CALCULER UNE SEULE FOIS LA MATRICE GLOBALE PROFIL  n x n
C [ N(V0**2+W0**2) -Alfa LAPLACIEN] et de la FACTORISER L D tL CROUT
C  POUR RESOUDRE LE SYSTEME
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
C         (2 GRADIENT CONJUGUE SUR UNE MATRICE MORSE  NON PROGRAMME)
C MNLPLP : ADRESSE MCN DES POINTEURS SUR LES COEFFICIENTS DIAGONAUX DE LA MATRIC
C MNLPLI : ADRESSE MCN DES POINTEURS SUR LES COEFFICIENTS DIAGONAUX DE LA MATRIC
C MNLPCO : ADRESSE MCN DU NUMERO DE COLONNE DES COEFFICIENTS DE LA MATRICE MORSE
C NIVEAU : NIVEAU DE FACTORISATION INCOMPLETE DE LA MATRICE DE PRECONDITIONNEMEN
C          A CHOISIR PARMI 0 1 2
C NBRDKG : NOMBRE DE REELS DOUBLE PRECISION DES MATRICES M ET K
C NCODSK : 1 SI MATRICE DE CONDUCTIVITE SYMETRIQUE
C          0 SI MATRICE DE CONDUCTIVITE DIAGONALE
C         -1 SI MATRICE DE CONDUCTIVITE NON SYMETRIQUE
C MNUG00 : ADRESSE MCN DE L'ONDE A L'INSTANT INITIAL (NBNOMA,2)
C
C RDeltaT: PAS CONSTANT DU TEMPS REEL SIMPLE PRECISION
C DTSTOC : PAS CONSTANT DU TEMPS ENTRE 2 STOCKAGES DU VECTEUR"ONDENLSE
C TPSINI : TEMPS INITIAL DU CALCUL
C TPSFIN : TEMPS FINAL   DU CALCUL
C NBVECT : NUMERO DU DERNIER VECTEUR TEMPERATURE STOCKE
C MXVECT : NOMBRE DE VECTEUR"ONDENLSE A STOCKER
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE DE L'OBJET
C          =2 FOIS LE NOMBRE DE NOEUDS DU MAILLAGE
C          POUR LA PARTIE REELLE ET IMAGINAIRE DE L'ONDE AUX NOEUDS
C
C SORTIES:
C --------
C NTVECT : NUMERO      DU TMS VECTEUR"ONDENLSE DE L'OBJET
C MNVECT : ADRESSE MCN DU TMS VECTEUR"ONDENLSE DE L'OBJET
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray    Aout 2011
C MODIFS : ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Octobre 2013
C23456---------------------------------------------------------------012
      PARAMETER         (MXTYEL=7)
      PARAMETER         (ITERmMX=32)
      PARAMETER         (NBFOKGMX=8)
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
      REAL              RDeltaT, DTSTOC, TPSINI, TPSFIN
      DOUBLE PRECISION   DeltaT,
     %                  MODUMX0, MODUMX, Testm, WCAMIN(2), WCAMAX(2),
     %                  WEXMIN(2), WEXMAX(2), ERRMAX(2), ERRMXR, ERRMXI,
     %                  D, DINFO, CPU0, CPU1, CPUPKG, CPUITE
      INTEGER           NOFOWE(2)

      DOUBLE PRECISION  RELMIN, PENALI, TGV, UNSTGV

      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) 'Depart NLSESIMPL0:'

C     RELMIN : PLUS PETIT REEL SERVANT DE MARQUEUR DE NON UTILISATION
      RELMIN = -1D28
C
C     DEBUT DU TEMPS DE CALCUL
      CPU0   = DINFO( 'CPU' )
      CPUPKG = 0D0
      CPUITE = 0D0
C
C     ATTENTION: FACTORISATION DE CROUT IMPOSEE  [KG] = L D tL
      NORESO = 1
C
Cccc     PENALISATION DE LA CONDITION DE DIRICHLET PAR FOURIER
Cccc     ECHANGE=PENALI=EPSILON ET FORCE=TEMPERATURE*EPSILON
ccc      PENALI = 1D30
C     CONDITION DE DIRICHLET PRISE EN COMPTE DIRECTEMENT
      PENALI = 0D0

C     COEFFICIENT DE PRISE EN COMPTE DES DL FIXES DE DIRICHLET
      TGV    = 1D30
      UNSTGV = 1D0 / TGV
C
C     PAS DE TEMPS EN DOUBLE PRECISION
      DeltaT = RDeltaT
C
C     L'ADRESSE DE L'ONDE A L'ITERATION 0 DU POINT FIXE
      MNTHET  = MNUG00
C     PASSAGE DE UG00 DANS LE common de incl/cthet.inc
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
C     NOMBRE DE NOEUDS DU MAILLAGE dans incl/cthet.inc
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
C     AFFICHAGE DE L'ONDE AU TEMPS INITIAL UG00(NBNOMA,2) PAR COMPOSANTES
C     -------------------------------------------------------------------
      CALL AFNLSE( 5,   1, MNXYZN, NBNOMA, 1, MCN(MNUG00),
     %             WCAMIN, WCAMAX, NOFOWE, WEXMIN, WEXMAX, ERRMAX )
C
C     CALCUL du Max |UG00(NOEUD)| POUR DETECTER PLUS TARD UNE EXPLOSION DE L'ONDE
C     ---------------------------------------------------------------------------
      MN0 = ( MNUG00 - 1 ) / MOREE2
      MODUMX0 = 0D0
      DO I=1,NBNOMA
         MN = MN0 + I
         D  = DMCN(MN)**2 + DMCN(MN+NBNOMA)**2
         IF( D .GT. MODUMX0 ) MODUMX0 = D
      ENDDO
      MODUMX0 = SQRT( MODUMX0 )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'Max |Onde(t=',TEMPS,',NOEUD)|=',MODUMX0
      ELSE
         WRITE(IMPRIM,*) 'Max |Wave(t=',TEMPS,',NODE)|=',MODUMX0
      ENDIF

C     DECLARATION DES VECTEURS GLOBAUX NECESSAIRES
C     VECTEUR GLOBAL BG(NTDL) SECOND MEMBRE GLOBAL
C     VECTEUR GLOBAL UGtn  DE L'ONDE A L'INSTANT tn      EST DECLARE(NBNOMA,2)
C     VECTEUR GLOBAL UGtm  DE L'ONDE A L'INSTANT tn+1m   EST DECLARE(NBNOMA,2)
C     VECTEUR GLOBAL UGtm1 DE L'ONDE A L'INSTANT tn+1m+1 EST DECLARE(NBNOMA,2)
C     ========================================================================
      CALL TNMCDC( 'REEL2', 4*NTDL, MNBG )
      MNUGtn  = MNBG    + MOREE2 * NTDL
      MNUGtm  = MNUGtn  + MOREE2 * NTDL
      MNUGtm1 = MNUGtm  + MOREE2 * NTDL
C     PASSAGE DE L'ADRESSE de UGtn DANS LE common de incl/cthet.inc
      MNTHETn= MNUGtn

C     DECLARATION DE LA MATRICE GLOBALE [KG] NBNOMA x NBNOMA
C     ======================================================
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
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10290) NBRDKG, NBRDKG/NBNOMA
      ELSE
         WRITE(IMPRIM,20290) NBRDKG, NBRDKG/NBNOMA
      ENDIF
10290 FORMAT('UNE MATRICE PROFIL DE',I15,
     %' REELS DOUBLE PRECISION avec 1/2 LARGEUR DE BANDE MOYENNE =',I9)
20290 FORMAT('ONE SKYLINE MATRIX of',I15,' DOUBLE REALS with an HALF WID
     %TH AVERAGE=',I9)
C
C     ALLOCATION DYNAMIQUE EN FORTRAN 90 DE LA MATRICE PROFIL KG
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

C     LE PROCHAIN TEMPS POUR STOCKER LE VECTEUR"ONDENLSE
C     --------------------------------------------------
      TSTOC  = TPSINI + DTSTOC
      MNTEMP = MNVECT + WECTEU + MOREE2 * NTDL * (NBVECT-1)
      MNTIME = MNVECT + WECTEU + MOREE2 * NTDL * MXVECT - 1
C     LE NOMBRE DE DL D'UN VECTEUR
      MCN( MNVECT + WBCOVE ) = NTDL
C     LE NOMBRE DE VECTEURS STOCKES
      MCN( MNVECT + WBVECT ) = NBVECT
C     AU DEPART PAS DE TEMPS INDIQUES
      MCN( MNVECT + WBCPIN ) = 0

C     STOKAGE: Testm, Max |U|, eventuellement ERREUR PR, ERREUR PI
C     ------------------------------------------------------------
      IF( NOFOWE(1) .GT. 0 .AND. NOFOWE(2) .GT. 0 ) THEN
C        STOKAGE: Testm, Max |U|, eventuellement ERREUR PR, ERREUR PI
         NBVERR = 2 + 2
      ELSE
C        STOKAGE: Testm, Max U(N)  seulement
         NBVERR = 2
      ENDIF
C     1+NBVERR POUR STOCKER LE TEMPS en PLUS
C     INFO COMPLEMENTAIRE DANS VECTEUR"TESTM_ERREUR
      MOERRT = (1+NBVERR) * MXVECT
      CALL TNMCDC( 'REEL', MOERRT, MNERRT )

C     TEMPS = TEMPS INITIAL
      RMCN( MNERRT ) = TEMPS
C     TEST D'ARRET DES ITERATIONS m
      RMCN( MNERRT+1 ) = 0
C     Max |Um|
      RMCN( MNERRT+2 ) = REAL( MODUMX0 )
      IF( NBVERR .GT. 2 ) THEN
C        ERREUR RELATIVE PARTIE REELLE
         RMCN( MNERRT+3 ) = 0
C        ERREUR RELATIVE PARTIE IMAGINAIRE
         RMCN( MNERRT+4 ) = 0
      ENDIF

C     TEMPS DE CALCUL DE LA DECLARATION DE KG ET AUTRES VECTEURS
      CPU1   = DINFO( 'CPU' )
      CPUPKG = CPU1 - CPU0
      CPU0   = CPU1

C     NOMBRE TOTAL D'ITERATIONS m EFFECTUEES
      NBCABG = 0

C     NOMBRE DE PAS DE TEMPS CALCULES
      NBPADT  = 0
      NBLOCm0 = 0

C     NOMBRE DE FORMATION DE LA MATRICE KG (-1 POUR LA TOUTE PREMIERE)
      NBFOKG = -1

C     CONSTRUCTION DE LA MATRICE [KG] SYMETRIQUE NON DIAGONALE PROFIL
C     [ K(Alfa,g,N,V0,W0)] = [N(V0**2+W0**2) - Alfa LAPLACIEN]
C     RETOUR ICI EN CAS DE MODIFICATION DU PAS DE TEMPS
C     ===============================================================
C     LE VECTEUR UG00 EST DONNE PAR COMPOSANTES UG00(NBNOMA,2)
C     LES DL DE LA PARTIE REELLE PUIS CEUX DE LA PARTIE IMAGINAIRE
C     MNUG00 EST RANGE DANS MNTHET POUR UTILISER UG00 A PARTIR DE incl/cthet.inc
 50   MNTHET = MNUG00
C     PASSAGE DE UG00 DANS LE common de incl/cthet.inc
      MNTHET0 = MNUG00

C     COPIE DE UG00(NBNOMA,2) DANS UGtm(NBNOMA,2)
      CALL TRTATD( MCN(MNUG00), MCN(MNUGtm), NTDL )
      ID00tm = 1

ccc      IF( NBFOKG .LT. 0 ) THEN
cccC        ANNULATION DU VECTEUR UG00 POUR EVITER QUE BG=0
ccc         CALL  AZEROD( NTDL, MCN(MNUG00) )
ccc      ENDIF

C     FORMATION DE KG
      CALL NLSESIKG0( PENALI, NDIM,   MNXYZN, NBTYEL, MNNPEF, NDPGST,
     %                MNTPOB, MNTAUX, NUMIOB, NUMAOB, MNDOEL,
     %                MNTHER, MNTAEL, MNX,
     %                NORESO, NCODSK, NBRDKG, MNLPLI, MNLPCO, KG,
     %                NBPTAF, IERR  )
      IF( IERR .NE. 0 ) GOTO 9999

      call affvect( 'KG PROFIL SANS CL', 20    , KG )
      call afl1ve(  'KG PROFIL SANS CL', NBRDKG, KG )

C     CONSTRUCTION DE LA LISTE ET DES VALEURS DES DL FIXES DIRICHLET
C     AVEC UNE NUMEROTATION GLOBALE DES DL FX PAR NOEUDS (2,NBNOMA)
      CALL DLFXLI( 1,      NBNOMA, NDIM,   NBTYEL, MNNPEF, NDPGST,
     %             MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %             NBRDLX, MONDLX, MNNDLX, MNVDLX, IERR )
C
C     PRISE EN COMPTE DES DL FIXES SUR LA MATRICE KG NBNOMA x NBNOMA
C     RESTRICTION UN NOEUD DIRICHLET => PARTIES REELLE ET IMAG FIXEES
      CALL DLFXKG( NBRDLX, MCN(MNNDLX), TGV, NBNOMA, MCN(MNLPLI), KG )

C     FACTORISATION DE CROUT DE [KG]=L D tL
      CALL CRMC1D( MCN(MNLPLI), KG, NBNOMA, 1E-10, 1, KG, IERR )
C     IERR : 0 SI   AUCUN     PIVOT<EPS
C            1 SI AU MOINS UN PIVOT<EPS
C     ON ENTEND PAR PIVOT (A0(IDIAGONAL)-SA)/A0(IDIAGONAL)

      call affvect( 'CROUT L D tL PROFIL', 20    , KG )
      call afl1ve(  'CROUT L D tL PROFIL', NBRDKG, KG )

      IF( IERR .NE. 0 ) THEN
C        MATRICE NON INVERSIBLE => ARRET
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

C     SOMME DES TEMPS DE CALCUL DE CONSTRUCTION DE KG
      CPU1   = DINFO( 'CPU' )
      CPUPKG = CPUPKG + CPU1  - CPU0
      CPU0   = CPU1

C     ICI KG = L D tL =  [K(-Alfa +N(UG0**2))] MATRICE GLOBALE NBNOMA x NBNOMA
C     ICI MNUGtm CONTIENT UGtm(NBNOMA,2) L'ONDE A L'INSTANT INITIAL

C     ##################################################################
C     ##                                                              ##
C     ##  LA BOUCLE EN TEMPS AVEC DES PAS DE TEMPS CONSTANTS = DeltaT ##
C     ##                                                              ##
C     ##################################################################
C
C     ICI LE NOUVEAU TEMPS OU SE FAIT LE CALCUL
C     =========================================
 100  TEMPS = REAL( TEMPS + DeltaT )
      ITERm = 0

C     SOLUTION UGtm AU TEMPS tn => UGtn(NBNOMA,2) AU TEMPS tn
      CALL TRTATD( MCN(MNUGtm), MCN(MNUGtn), NTDL )

C     MODIFICATION DU VECTEUR UGtm POUR EVITER QUE BG=0 QUAND Un+1m=U00
      IF( ID00tm .EQ. 1 ) THEN
         MN = (MNUGtm-1)/MOREE2
         DO I=1,NTDL
            DMCN(MN+I) = 1.01D0 * DMCN(MN+I)
         ENDDO
         ID00tm = 0
      ENDIF

C     ICI LES ITERATIONS m DE POINT FIXE SONT MISES EN OEUVRE
C     =======================================================
C     UNE NOUVELLE ITERATION m DE POINT FIXE
 110  ITERm = ITERm + 1
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10100) TEMPS, DeltaT, ITERm, NBCABG
      ELSE
         WRITE(IMPRIM,20100) TEMPS, DeltaT, ITERm, NBCABG
      ENDIF
10100 FORMAT(/'Au TEMPS',G14.6,'  PAS de TEMPS=',G15.7,' ITERm=',I3,
     %'  Nb ITERm total=',I10,'  -> CALCUL de l''ONDE COMPLEXE'/
     % 140('*'))
20100 FORMAT(/'At TIME',G14.6,' TIME STEP=',G15.7,' ITERm=',I3,
     %'  Total Nb ITERm=',I10,'  -> COMPUTATION of the COMPLEX WAVE'/
     % 140('*'))
C
C     GENERATION DU SECOND MEMBRE BG=Fn+1,m+1(NBNOMA,2) AU TEMPS tn+1
C     BG1 = Fr(tn+1)
C         + Rho/dt (Wn+1m-Wn) -OmegaZ ( x dWn+1m/dy - y dWn+1m/dx )
C         + ( N(V0**2+W0**2)-N(Vn+1m**2+Wn+1m**2) ) Vn+1m
C     BG2 = Fi(tn+1)
C         - Rho/dt (Vn+1m-Vn) +OmegaZ ( x dVn+1m/dy - y dVn+1m/dx )
C         + ( N(V0**2-W0**2)-N(Vn+1m**2+Wn+1m**2) ) Wn+1m
C     ---------------------------------------------------------------
      MNTHET = MNUGtm
      call afl1ve( 'UGtn avant NLSEIMBG', NTDL, DMCN((MNUGtn+1)/MOREE2))
      call afl1ve( 'UGtm avant NLSEIMBG', NTDL, DMCN((MNUGtm+1)/MOREE2))
      CALL NLSEIMBG( DeltaT, PENALI, NDIM,   MNXYZN, NBNOMA,
     %               NBTYEL, MNNPEF, NDPGST,
     %               MNTPOB, MNTAUX, NUMIOB, NUMAOB, MNDOEL,
     %               MNX,    MCN(MNUGtn), MCN(MNUGtm),
     %               MCN(MNBG), MCN(MNBG), IERR  )
      call afl1ve( 'BG n+1 m PR apres NLSEIMBG', NBNOMA,
     %              DMCN((MNBG+1)/MOREE2) )
      call afl1ve( 'BG n+1 m PI apres NLSEIMBG', NBNOMA, 
     %                                DMCN((MNBG+1)/MOREE2+NBNOMA) )
      IF( IERR .NE. 0 ) GOTO 9999

C     AJOUT DES EVENTUELLES FORCES PONCTUELLES
C     ----------------------------------------
      IF( IESOPO .GT. 0 .AND. NDIM .GT. 1 ) THEN
C
C        CONSTRUCTION NODL(NBNOMA,2) ET VALEUR DES FORCES PONCTUELLES
         CALL THSNFX( NTDL,   NBTYEL, MNNPEF, NDPGST,
     %                MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %                NBFNFX, MONFNX, MNNFNX, MOVFNX, MNVFNX )
C
C        ASSEMBLAGE DES FORCES PONCTUELLES DANS BG(NBNOMA,2)
         CALL ASFONO( NTDL, 1, NBFNFX,MCN(MNNFNX),MCN(MNVFNX),RELMIN,
     %                MCN(MNBG) )
C
      ENDIF
C
C     PRISE EN COMPTE DES DL FIXES(2,NBNOMA) SUR LE VECTEUR GLOBAL BG(NBNOMA,2)
C     -------------------------------------------------------------------------
      CALL DLFXBG( NBRDLX, MCN(MNNDLX), MCN(MNVDLX), TGV, NBNOMA,
     %             MCN(MNBG) )
C
C     RESOLUTION DU SYSTEME FACTORISE PAR CROUT de KG = L D tL
C     L D tL UGtm1(NBNOMA,1) = BG(NBNOMA,1) => Vtn+1m+1 PARTIE REELLE DE L'ONDE
C     -------------------------------------------------------------------------
      CALL DRCRPR( NBNOMA, NCODSK, MCN(MNLPLI), KG,
     %             MCN(MNBG),    3,
     %             MCN(MNUGtm1), IERR )
      call afl1ve( 'VG n+1 m Partie Reelle', NBNOMA,
     %             DMCN((MNUGtm1+1)/MOREE2) )

C     RESOLUTION DU SYSTEME FACTORISE PAR CROUT de KG = L D tL
C     L D tL UGtm1(NBNOMA,2) = BG(NBNOMA,2) => Wtn+1m+1 PARTIE IMAGINAIRE DE L'ONDE
C     -----------------------------------------------------------------------------
      CALL DRCRPR( NBNOMA, NCODSK, MCN(MNLPLI), KG,
     %             MCN(MNBG   +MOREE2*NBNOMA), 3,
     %             MCN(MNUGtm1+MOREE2*NBNOMA), IERR )
      call afl1ve( 'WG n+1 m Partie Imaginaire',
     %             NBNOMA, DMCN((MNUGtm1+1)/MOREE2+NBNOMA) )

C     NETTOYAGE DES PETITES VALEURS DE UGtm1 DUES A TGV
      MN = ( MNUGtm1 - 1 ) / MOREE2
      DO I=1,NTDL
         D = ABS( DMCN(MN+I) )
         IF( D .LT. UNSTGV ) THEN
            DMCN(MN+I)=0D0
         ENDIF
      ENDDO

C     AFFICHAGE DE L'ONDE UGtm1(NBNOMA,2) AU TEMPS tn+1m+1
C     AFFICHAGE DU MODULE DE L'ERREUR SI SOLUTION EXACTE CONNUE
C     ---------------------------------------------------------
      CALL AFNLSE( 5,   1, MNXYZN, NBNOMA, 1, MCN(MNUGtm1),
     %             WCAMIN, WCAMAX, NOFOWE, WEXMIN, WEXMAX, ERRMAX )

C     AFFICHER LES NORMES ABSOLUES ET RELATIVES DE 
C     UGn+1m(NBNOMA,2), UGn+1,m+1(NBNOMA,2)
C     ET CALCULER LA VALEUR Testm POUR TESTER LA FIN DES ITERATIONS m
C     ---------------------------------------------------------------
      CALL NLSEITER( 0, TEMPS, ITERm, NBNOMA,
     %             MCN(MNUGtm), MCN(MNUGtm), MCN(MNUGtm1), MCN(MNUGtm1),
     %             MODUMX, Testm )
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
C     TEMPS
      RMCN( MN ) = TEMPS
C     TEST D'ARRET DES ITERATIONS m
      RMCN( MN+1 ) = REAL( Testm )
C     Max Um
      RMCN( MN+2 ) = REAL( MODUMX )

      IF( NBVERR .GT. 2 ) THEN
C        ERREUR RELATIVE PARTIE REELLE
         RMCN( MN+3 ) = REAL( ERRMXR )
C        ERREUR RELATIVE PARTIE IMAGINAIRE
         RMCN( MN+4 ) = REAL( ERRMXI )
      ENDIF

C     UN CALCUL DE UGn+1 m+1 DE PLUS EFFECTUE
      NBCABG = NBCABG + 1
C
C     TRACE DU MODULE DE L'ONDE A L'ITERATION m
C     -----------------------------------------
      IF( INTERA .GE. 1 ) THEN

C        L'ECRAN PIXELS EST EFFACE (PAS DE SCINTILLEMENT)
         CALL EFFACEMEMPX

C        CONSTRUCTION DU MODULE AUX NBNOMA NOEUDS DE L'ONDE
         MN  = ( MNMODU  - 1 ) / MOREE2
         MN1 = ( MNUGtm1 - 1 ) / MOREE2
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
         KNOM(I+1:I+14) = ' Time Step NB='
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
C     TEMPS DE CALCUL DE L'ITERATION DE POINT FIXE
      CPU1   = DINFO( 'CPU' )
      CPUITE = CPUITE + CPU1 - CPU0
      CPU0   = CPU1
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
10131    FORMAT(/100('&')/'EXPLOSION DETECTEE au TEMPS', G14.6/100('&'))
20131    FORMAT(/100('&')/'DETECTION of WAVE EXPLOSION at TIME', G14.6/
     %          100('&'))
C        SAUVEGARDE DES RESULTATS ACTUELS UGtm1
         GOTO 9900
      ENDIF
C
      IF( MOD(NBCABG,1000) .EQ. 0 .AND. NBFOKG .LE. 1 ) THEN
C
C        TOUTES LES 1000 ITERATIONS SANS MODIFICATION DE KG
C        LA MATRICE KG EST REMISE A JOUR POUR SON TERME NON LINEAIRE
C        C-A-D Vtn+1 Wtn+1 => VO et W0 et DeltaT = DeltaT * 1.5D0
C        -----------------------------------------------------------
C        Retour aux RESULTATS du dernier TEMPS CALCULE
         TEMPS = REAL( TEMPS - DeltaT )
C        LE NOUVEAU TEMPS INITIAL EST tn
         TEMPSINI = TEMPS
C        LA SOLUTION INITIALE UG00 EST REMPLACEE PAR LA DERNIERE
C        SOLUTION CALCULEE C'EST A DIRE UGtm1(NBNOMA,2)
         CALL TRTATD( MCN(MNUGtm1), MCN(MNUG00), NTDL )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10200) DeltaT, DeltaT*1.5D0
         ELSE
            WRITE(IMPRIM,20200) DeltaT, DeltaT*1.5D0
         ENDIF
C        NOUVEAU PAS DE TEMPS
         DeltaT   = DeltaT * 1.5D0
         PasTemps = REAL( DeltaT )
         NBFOKG   = 0
         GOTO 50
C
      ENDIF
C
 150  IF( ITERm .GT. 4 .AND. Testm .GT. 1D-2 ) THEN
C
C        LE TEST m EST TROP GRAND => REDUIRE LE PAS DE TEMPS DeltaT
C        EST IL ENCORE POSSIBLE DE DIVISER LE PAS DE TEMPS PAR 2?
C        ----------------------------------------------------------
         IF( NBFOKG .GE. NBFOKGMX ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='STOP apres    DIVISIONS du PAS DE TEMPS'
               KERR(2)='Reduire le PAS de TEMPS INITIAL TROP GRAND?'
            ELSE
               KERR(1)='STOP after    DIVISIONS of the TIME STEP'
               KERR(2)='Reduce the INITIAL TIME STEP TOO GREAT?'
            ENDIF
            WRITE(KERR(1)(12:13),'(I2)') NBFOKGMX
            CALL LEREUR
            IERR = 29
C           MISE DANS LE VECTEUR"ONDENLSE DES VECTEURS DEJA STOCKES
            GOTO 9900
         ENDIF
C
C        OUI: RECONSTRUCTION DE KG avec un NOUVEAU PAS DE TEMPS/2
C        ========================================================
         NBFOKG = NBFOKG + 1
C        Retour aux RESULTATS du dernier TEMPS CALCULE
         TEMPS = REAL( TEMPS - DeltaT )
C        LE NOUVEAU TEMPS INITIAL EST tn
         TEMPSINI = TEMPS
C        LA SOLUTION INITIALE UG00 EST REMPLACEE PAR LA DERNIERE
C        SOLUTION CALCULEE C'EST A DIRE UGtn(NBNOMA,2)
         CALL TRTATD( MCN(MNUGtn), MCN(MNUG00), NTDL )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10200) DeltaT, DeltaT/2D0
         ELSE
            WRITE(IMPRIM,20200) DeltaT, DeltaT/2D0
         ENDIF
10200    FORMAT(/'Changement du PAS de TEMPS',G14.6,
     %   ' par le NOUVEAU PAS de TEMPS=',G14.6,
     %   '  + RECONSTRUCTION [KG]')
20200    FORMAT(/'Changement of TIME STEP=',G14.6,
     %   ' to the NEW TIME STEP=',G14.6,
     %   '  + [KG] is CONSTRUCT AGAIN')
C        NOUVEAU PAS DE TEMPS
         DeltaT = DeltaT / 2D0
         PasTemps = REAL( DeltaT )
         GOTO 50
C
      ENDIF
C
      IF( Testm .GT. 1D-3 ) THEN
C
C        NON CONVERGENCE DU POINT FIXE => UNE ITERATION DE PLUS A FAIRE
C        --------------------------------------------------------------
         IF( ITERm .GE. ITERmMX ) THEN
            NBLGRC(NRERR) = 3
            WRITE(KERR(5)(1:6),'(I6)') ITERmMX
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: POINT FIXE NON ATTEINT apres'
               KERR(3) = 'PAS de TEMPS TROP GRAND a REDUIRE?'
            ELSE
               KERR(1) = 'ERROR: FIX POINT NOT CONVERGED after'
               KERR(3) = 'TIME STEP TOO GREAT to REDUCE?'
            ENDIF
            KERR(2) = KERR(5)(1:6) // ' ITERATIONS'
            CALL LEREUR
            Testm = 1D0
            GOTO 150
ccc            IERR = 29
cccC           MISE DE UGtm1 DANS LE VECTEUR"ONDENLSE DES VECTEURS DEJA STOCKES
ccc            GOTO 9900
         ENDIF
C
C        LE NOUVEAU VECTEUR ONDE DE L'ITERATION m-1 EST UGtm1 CALCULE
         MN      = MNUGtm
         MNUGtm  = MNUGtm1
         MNUGtm1 = MN
         GOTO 110
      ENDIF
C
C     CONVERGENCE DES ITERATIONS DE POINT FIXE
C     ----------------------------------------
      NBPADT = NBPADT + 1
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10140) TEMPS, ITERm
      ELSE
         WRITE(IMPRIM,20140) TEMPS, ITERm
      ENDIF
C
10140 FORMAT('Au TEMPS',G14.6,' ITER m=',I3,' CONVERGENCE')
20140 FORMAT('At TIME', G14.6,' ITER m=',I3,' CONVERGENCE')
C
c     LE VECTEUR DES DL DE L'ONDE ACTUELLE EST IL A STOCKER?
C     ------------------------------------------------------
      IF( TEMPS .GE. TSTOC*0.9999 .AND. NBVECT .LT. MXVECT ) THEN
C
C        OUI: STOCKAGE DE L'ONDE UGtm1(NBNOMA,2) A CET INSTANT TEMPS
         MNTEMP = MNTEMP + NTDL * MOREE2
         CALL TRTATD( MCN(MNUGtm1), MCN(MNTEMP), NTDL )
C        LE NOMBRE DE VECTEURS TEMPERATURE STOCKES
         NBVECT = NBVECT + 1
C        LE TEMPS DE STOCKAGE
         RMCN(MNTIME+NBVECT) = TEMPS
C
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10150) TEMPS, NBVECT, NTDL, TSTOC, DTSTOC
         ELSE
            WRITE(IMPRIM,20150) TEMPS, NBVECT, NTDL, TSTOC, DTSTOC
         ENDIF
C
10150 FORMAT('Au TEMPS',G14.6,' STOCKAGE DE l''ONDE',I5,
     %' de',I9,' DEGRES de LIBERTE TSTOC=',G14.6,' dtSTOC=',G14.6/)
20150 FORMAT('At TIME',G14.6,' STORAGE of the WAVE',I5,
     %' of',I9,' DEGREES of FREEDOM TSTOC=',G14.6,' dtSTOC=',G14.6/)
C
C        LE PROCHAIN TEMPS DE STOCKAGE
         TSTOC = TSTOC + DTSTOC

      ENDIF
C
C     ETUDE EN TEMPS TERMINEE?
C     ------------------------
      IF( TEMPS + RDeltaT .LT. TPSFIN*1.00001 ) THEN
C
         IF( Testm .LT. 5D-4 .AND. ITERm .LE. 1 ) THEN
C
C           POUR ACCELERER RECONSTRUCTION DE KG avec un PAS DE TEMPS*1.5
C           ============================================================
            NBFOKG = MAX( NBFOKG-1, 0 )
C           LE NOUVEAU TEMPS INITIAL EST tn
            TEMPSINI = TEMPS
C           SOLUTION INITIALE REMPLACEE PAR LA DERNIERE SOLUTION CALCULEE
            CALL TRTATD( MCN(MNUGtm1), MCN(MNUGtn), NTDL )
            CALL TRTATD( MCN(MNUGtm1), MCN(MNUG00), NTDL )
C
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10200) DeltaT, DeltaT*1.5D0
            ELSE
               WRITE(IMPRIM,20200) DeltaT, DeltaT*1.5D0
            ENDIF
C           NOUVEAU PAS DE TEMPS
            DeltaT = DeltaT * 1.5D0
            PasTemps = REAL( DeltaT )
            GOTO 50
C
         ENDIF
C
         NBLOCm = NBPADT / 2
         IF( NBPADT .LT. 10 .OR. NBLOCm .GT. NBLOCm0 ) THEN
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
C        MISE A JOUR: MNUGtm1 DEVIENT MNUGtm
C        LA DERNIERE SOLUTION CALCULEE EST MAINTENANT A L'ADRESSE MNUGtm
         MN      = MNUGtm
         MNUGtm  = MNUGtm1
         MNUGtm1 = MN
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
C     STOCKAGE DE L'ONDE UGtm1(NBNOMA,2) A CET INSTANT TEMPS
 9900 IF( NBVECT .LT. MXVECT ) THEN
C        LE VECTEUR SOLUTION SUIVANT
         MNTEMP = MNTEMP + NTDL * MOREE2
C        LE NOMBRE DE VECTEURS TEMPERATURE STOCKES
         NBVECT = NBVECT + 1
      ELSE
C        LE VECTEUR SOLUTION MXVECT
         NBVECT = MXVECT
      ENDIF
      CALL TRTATD( MCN(MNUGtm1), MCN(MNTEMP), NTDL )
C     LE TEMPS DE STOCKAGE
      RMCN(MNTIME+NBVECT) = TEMPS
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10150) TEMPS, NBVECT, NTDL, TSTOC, DTSTOC
      ELSE
         WRITE(IMPRIM,20150) TEMPS, NBVECT, NTDL, TSTOC, DTSTOC
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
C
C     LES TEMPERATURES A L'INSTANT INITIAL SONT A L'ADRESSE MNTEMP
      MNTEMP = MNVECT + WECTEU
C     L'ADRESSE -1 DU PREMIER TEMPS STOCKE DERRRIERE LES VECTEURS TEMPERATURE
      MNTIME = MNTEMP + NTDL * NBVECT * MOREE2 -1

C     CONSTRUCTION DE VECTEUR"TESTM_ERREUR A PARTIR DE L'EVENTUEL ANCIEN
C     AVEC AJOUT DU TMC MNERRT ACTUEL
C     ==================================================================
      CALL TITSMXERR( KNOMOB, MOREE2, TPSINI, NBVERR, NBCABG, MNERRT )
C
C     DESTRUCTION DES TMC DEVENUS INUTILES
C     ====================================
 9999 IF( IERKGALLOC .EQ. 0 ) DEALLOCATE( KG )
      IF( MNBG   .GT. 0 ) CALL TNMCDS( 'REEL2',  NTDL*4,  MNBG   )
      IF( MNNFNX .GT. 0 ) CALL TNMCDS( 'ENTIER', MONFNX,  MNNFNX )
      IF( MNVFNX .GT. 0 ) CALL TNMCDS( 'REEL2',  MOVFNX,  MNVFNX )
      IF( MNNDLX .GT. 0 ) CALL TNMCDS( 'ENTIER', MONDLX,  MNNDLX )
      IF( MNVDLX .GT. 0 ) CALL TNMCDS( 'REEL2',  MONDLX,  MNVDLX )
      IF( MNLPLI .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL/2+1,MNLPLI )
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
19002 FORMAT(/'NLSESIMPL0: STOCKAGE MATRICE PROFIL =',I15,' MOTS'/ )
29002 FORMAT(/'NLSESIMPL0: SKYLINE MATRIX STORAGE=',I15,' MEMORY WORDS')
C
C     AFFICHAGE DES TEMPS CALCUL
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,19000) CPUPKG, CPUITE, CPUPKG+CPUITE
      ELSE
         WRITE(IMPRIM,29000) CPUPKG, CPUITE, CPUPKG+CPUITE
      ENDIF
19000 FORMAT(/
     %'CPU TRAITEMENT DES MATRICES MG KG =',F12.2,' SECONDES CPU'/
     %'CPU ITERATIONS TEMPS + PT FIXE    =',F12.2,' SECONDES CPU'/
     %'CPU TOTAL DE CALCUL DE LA SOLUTION=',F12.2,' SECONDES CPU')

29000 FORMAT(/
     %'CPU of MG & KG MATRIX CONSTRUCTION=',F12.2,' CPU SECONDS'/
     %'CPU of TIME + FIX PT ITERATIONS   =',F12.2,' CPU SECONDS'/
     %'CPU TOTAL SOLUTION COMPUTATION    =',F12.2,' CPU SECONDS')

      RETURN
      END
