      SUBROUTINE NLSESIMPL( KNOMOB, MOREE2, NDIM,   NDPGST, MNXYZN,
     %                      NBTYEL, MNNPEF, MNTPOB, MNTAUX,
     %                      NUMIOB, NUMAOB, MNDOEL, IESOPO,
     %                      RELMIN, MNTHER, MNTAEL, MNX,
     %                      NORESO, MNLPLI, MNLPCO, NBRDKG, NCODSK,
     %                      MNUG00, RDeltaT,DTSTOC, TPSINI, TPSFIN,
     %                      NBVECT, MXVECT, NTDL,
     %                      NTVECT, MNVECT, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : NLSE: CALCULER L'ONDE COMPLEXE ET SES FLUX DANS UN DOMAINE
C ----- 2D ou 3D DE L'EQUATION NON LINEAIRE COMPLEXE DE SCHRODINGER
C       i Rho dU(t,X)/dt - Alfa LAPLACIEN U(t,X) + N(U(t,X)**2) U(t,X)
C                        -i OmegaZ (x d/dy - y d/dx) U(t,X) = F(t,X)
C       SELON TESTNL=6  UNE METHODE SEMI-IMPLICITE D'EULER MODIFIEE
C       AVEC UN PAS DE TEMPS INITIAL QUI S'ADAPTE 
C       MASSE, CONDUCTIVITE, ECHANGE INDEPENDANTS DU TEMPS ET ONDE
C       FORCE, FIXATION, COEFFICIENT DU TERME NON LINEAIRE PEUVENT
C       DEPENDRE DU TEMPS ET DES ITERATIONS DE POINT FIXE SONT NECESSAIRES
C       L'ONDE IMPOSEE (CONDITION DE DIRICHLET) EST TRAITEE PAR PENALISATION
C       AVEC ECHANGE=1/EPSILON et FORCE=ONDE/EPSILON ET 1/EPSILON=PENALI=1D30
C       ET DOIT ETRE AUX MEMES NOEUDS POUR LES PARTIES REELLE ET IMAGINAIRE
C       MAIS AVEC DES POSSIBLES VALEURS DIFFERENTES SI NECESSAIRE
C***********************************************************************
C VERSION 2:
C    n=0  V=v0      W=w0
C
C(1) m=0  Vn+1m=Vn  Wn+1m=Wn
C
C(2) [Rho/dt -Alfa LAPLACIEN +N(V0**2+W0**2)] Vn+1m+1 = 
C     Fr(tn+1) +Fi(tn+1) +Rho/dt (Vn +Wn+1m-Wn) +Alfa LAPLACIEN Wn+1m
C     +OmegaZ ( x d/dy - y d/dx )(Vn+1m-Wn+1m)
C     -N(Vn+1m**2+Wn+1m**2)(Vn+1m+Wn+1m) + N(V0**2+W0**2) Vn+1m
C
C    [Rho/dt -Alfa LAPLACIEN +N(V0**2+W0**2)] Wn+1m+1 =
C    -Fr(tn+1) +Fi(tn+1) +Rho/dt (-Vn+1m+Vn +Wn) -Alfa LAPLACIEN Vn+1m
C    +OmegaZ ( x d/dy - y d/dx )(Vn+1m+Wn+1m)
C    +N(Vn+1m**2+Wn+1m**2)(Vn+1m-Wn+1m) + N(V0**2+W0**2) Wn+1m
C                       
C    Si ||Vn+1m+1-Vn+1m||>Eps||Vn+1m+1|| ou ||Wn+1m+1-Wn+1m||>Eps||Wn+1m+1||
C    Alors m=m+1  Aller en (2)
C    Sinon Vn+1=Vn+1m+1  Wn+1=Wn+1m+1
C
C   Si TPSINI+(n+1)dt < TPSFIN Alors n=n+1; Aller en (1)
C   Fin
C
C L'IDEE EST DE CALCULER LA MATRICE GLOBALE PROFIL  n x n
C [Rho/dt -Alfa LAPLACIEN +N(V0**2+W0**2)] et de la FACTORISER L D tL CROUT
C POUR RESOUDRE LE SYSTEME
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
C          ICI PENALI VAUT 1D30 POUR LE PRENDRE EN COMPTE ou 0D0
C RELMIN : PLUS PETIT REEL SERVANT DE MARQUEUR DE NON UTILISATION
C
C MNTHER : 128 REELS DOUBLE PRECISION POUR LA MATRICE DE CONDUCTIVITE
C MOTAEL : NOMBRE DE MOTS DECLARES DU TABLEAU DES TABLEAUX ELEMENTAIRES
C MNTAEL : ADRESSE MCN DES TABLEAUX ELEMENTAIRES
C MNX    : ADRESSE MCN DU TABLEAU DES 3 COORDONNEES DES NOEUDS=POINTS DE L'EF
C
C NORESO : CODE RESOLUTION ET STOCKAGE DU SYSTEME LINEAIRE
C          1 FACTORISATION DE CROUT SUR UNE MATRICE PROFIL KG nxn
C          2 GRADIENT CONJUGUE SIMPLE STOCKAGE MORSE DE KG
C MNLPLP : ADRESSE MCN DES POINTEURS SUR LES COEFFICIENTS DIAGONAUX DE LA MATRIC
C MNLPLI : ADRESSE MCN DES POINTEURS SUR LES COEFFICIENTS DIAGONAUX DE LA MATRIC
C MNLPCO : ADRESSE MCN DU NUMERO DE COLONNE DES COEFFICIENTS DE LA MATRICE MORSE
C NIVEAU : NIVEAU DE FACTORISATION INCOMPLETE DE LA MATRICE DE PRECONDITIONNEMEN
C          A CHOISIR PARMI 0 1 2
C NBRDKG : NOMBRE DE REELS DOUBLE PRECISION DE LA  MATRICE KG
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
C AUTEUR : ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray    Mars 2014
C23456---------------------------------------------------------------012
      PARAMETER         (MXTYEL=7)
      PARAMETER         (ITERMX=32)
      PARAMETER         (MXREDT=8)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/donthe.inc"
      include"./incl/donele.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__erreurth.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___arete.inc"
      include"./incl/a___aretefr.inc"
      include"./incl/a___face.inc"
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
      include"./incl/nctyef.inc"
      include"./incl/trvari.inc"
      include"./incl/traaxe.inc"
      include"./incl/xyzext.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      INTEGER           NUMIOB(4), NUMAOB(4), MNDOEL(4)
      CHARACTER*(*)     KNOMOB
      CHARACTER*160     KNOM
      CHARACTER*16      TEXT
C
      DOUBLE PRECISION, allocatable, dimension(:) :: KG
      INTEGER           IERKGALLOC
      INTRINSIC         ALLOCATED

      REAL              RDeltaT, DTSTOC, TPSINI, TPSFIN, TSTOC
      DOUBLE PRECISION   DeltaT,
     %                  MODUMX0, MODUMX, Testm, WCAMIN(2), WCAMAX(2),
     %                  WEXMIN(2), WEXMAX(2), ERRMAX(2), ERRMXR, ERRMXI,
     %                  D, DINFO, CPU0, CPU1, CPUPKG, CPUITE
      INTEGER           NOFOWE(2), KGCUET(2)
      INTRINSIC         SQRT, REAL, DBLE

      DOUBLE PRECISION  RELMIN, PENALI, TGV, UNSTGV, STGV
C     RELMIN : PLUS PETIT REEL SERVANT DE MARQUEUR DE NON UTILISATION
      RELMIN = -1D28

      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) 'Depart NLSESIMPL:'
C
C     DEBUT DU TEMPS DE CALCUL
      CPU0   = DINFO( 'CPU' )
      CPUPKG = 0D0
      CPUITE = 0D0
C
Cccc     PENALISATION DE LA CONDITION DE DIRICHLET PAR FOURIER
Cccc     ECHANGE=PENALI=EPSILON ET FORCE=TEMPERATURE*EPSILON
ccc      PENALI = 1D30
C     CONDITION DE DIRICHLET PRISE EN COMPTE DIRECTEMENT
      PENALI = 0D0
C
C     MODE DE PRISE EN COMPTE DES CONDITIONS AUX LIMITES DE DIRICHLET
C     SUR LES SYSTEMES LINEAIRES avec LA MATRICE KG
C     ---------------------------------------------------------------
      TGV = 1D30
      IF( NORESO .EQ. 1 ) THEN
C        MATRICE PROFIL
         STGV = TGV
      ELSE IF( NORESO .EQ. 2 ) THEN
C        MATRICE MORSE
         STGV = 1D0
      ELSE
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='ERREUR: METHODE DE RESOLUTION Ax=b NON PROGRAMMEE'
         ELSE
            KERR(1)='ERROR: METHOD of Ax=b SOLVER NOT PROGRAMMED'
         ENDIF
         CALL LEREUR
         CALL ARRET( 100 )
         RETURN
      ENDIF

C     SEUIL DE MISE A ZERO DES COEFFICIENTS DE L'ONDE
      UNSTGV = 1000D0 / TGV
C
C     PAS DE TEMPS EN DOUBLE PRECISION
      DeltaT = RDeltaT

C     TEMPS DU STOCKAGE DE L'ONDE CALCULEE
      IF( LANGAG .EQ. 0 ) THEN
         print*,'NLSESIMPL: temps=',TEMPS,
     %       ' temps initial=',TPSINI,
     %       ' temps fin=',TPSFIN,
     %       ' Pas temps stockage=',DTSTOC
      ELSE
         print*,'NLSESIMPL: time=',TEMPS,
     %       ' initial time=',TPSINI,
     %       ' final time=',TPSFIN,
     %       ' storage time step=',DTSTOC
      ENDIF
C
C     ADRESSES MCN DES TABLEAUX DE TRAVAIL
      MNAUX1 = 0
      MNAUX2 = 0
      MNAUX3 = 0
      MNBG   = 0
      MNNFNX = 0
      MNVFNX = 0
      MNNDLX = 0
      MNVDLX = 0
      MNERRT = 0
      MNNDPFX= 0
      MNMODU = 0
      MNPIL3 = 0
      MNSOLE = 0
      MOVALS = 0
      MNVALS = 0
      MNNOPO = 0
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

      IF( NORESO .EQ. 1 ) THEN

C        MATRICE PROFIL POUR CROUT
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10290) NBRDKG, NBRDKG/NBNOMA
         ELSE
            WRITE(IMPRIM,20290) NBRDKG, NBRDKG/NBNOMA
         ENDIF
10290 FORMAT('UNE MATRICE PROFIL DE',I15,
     %' REELS DOUBLE PRECISION avec 1/2 LARGEUR DE BANDE MOYENNE =',I9)
20290 FORMAT('ONE SKYLINE MATRIX of',I15,' DOUBLE REALS with an HALF WID
     %TH AVERAGE=',I9)

      ELSE

C        MATRICE MORSE POUR GRADIENT CONJUGUE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10292) NBRDKG, NBRDKG/NBNOMA
         ELSE
            WRITE(IMPRIM,20292) NBRDKG, NBRDKG/NBNOMA
         ENDIF
10292 FORMAT('UNE MATRICE MORSE DE',I15,
     %' REELS DOUBLE PRECISION avec un nombre de COEFFICIENTS par LIGNE'
     %,I9)
20292 FORMAT('ONE MORSE MATRIX of',I15,' DOUBLE REALS WITH A LINE COEFFI
     %CIENT MEAN NUMBER=',I9)

      ENDIF

C     ALLOCATION DYNAMIQUE EN FORTRAN 90 DE LA MATRICE KG
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) 'ALLOCATION DEMAND  of',NBRDKG,
     %                ' DOUBLE PRECISION of [KG] MATRIX'
      ALLOCATE ( KG(1:NBRDKG), STAT=IERKGALLOC )
      IF( IERKGALLOC .NE. 0 ) THEN
         WRITE(IMPRIM,*)'ALLOCATION ERROR   of',NBRDKG,
     %                  ' DOUBLE PRECISION of [KG] MATRIX'
         IERR = IERKGALLOC
         GOTO 9999
      ENDIF
      WRITE(IMPRIM,*) 'ALLOCATION CORRECT of',NBRDKG,
     %                ' DOUBLE PRECISION of [KG] MATRIX'
      WRITE(IMPRIM,*)
C
C     VECTEUR GLOBAL BG(NTDL) AUXILIAIRE                 EST DECLARE(NBNOMA,2)
C     VECTEUR GLOBAL UGtn  DE L'ONDE A L'INSTANT tn      EST DECLARE(NBNOMA,2)
C     VECTEUR GLOBAL UGtm  DE L'ONDE A L'INSTANT tn+1m   EST DECLARE(NBNOMA,2)
C     VECTEUR GLOBAL UGtm1 DE L'ONDE A L'INSTANT tn+1m+1 EST DECLARE(NBNOMA,2)
C     ========================================================================
      CALL TNMCDC( 'REEL2', 4*NTDL, MNBG )
      MNUGtn  = MNBG    + MOREE2 * NTDL
      MNUGtm  = MNUGtn  + MOREE2 * NTDL
      MNUGtm1 = MNUGtm  + MOREE2 * NTDL
C
C     AFFICHAGE DE L'ONDE AU TEMPS INITIAL UG00(NBNOMA,2) PAR COMPOSANTES
C     -------------------------------------------------------------------
      CALL AFNLSE( 2,   1, MNXYZN, NBNOMA, 1, MCN(MNUG00),
     %             WCAMIN, WCAMAX, NOFOWE, WEXMIN, WEXMAX, ERRMAX )
C
C     CALCUL du Max |UG00(NOEUD)| POUR DETECTER UNE EXPLOSION DE L'ONDE
C     -----------------------------------------------------------------
ccc      Max( Max-Min Partie Reelle et imaginaire )
ccc      MODUMX0 = MAX( ABS(WCAMAX(1)-WCAMIN(1)),
ccc     %               ABS(WCAMAX(2)-WCAMIN(2)) )
      MN0 = ( MNUG00 - 1 ) / MOREE2
      MODUMX0 = 0D0
      DO I=1,NBNOMA
         MN = MN0 + I
         D  = DMCN(MN)**2 + DMCN(MN+NBNOMA)**2
         IF( D .GT. MODUMX0 ) MODUMX0 = D
      ENDDO
      MODUMX0 = SQRT( MODUMX0 )
      IF( LANGAG .EQ. 0 ) THEN
      WRITE(IMPRIM,*)'NLSESIMPL: Max |Onde(t=',TEMPS,',NOEUD)|=',MODUMX0
      ELSE
       WRITE(IMPRIM,*)'NLSESIMPL: Max |Wave(t=',TEMPS,',NODE)|=',MODUMX0
      ENDIF

C     BILAN DE L'INSTANT INITIAL
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

C     TABLEAU MC: TEMPS + Testm Max |U|
C                       + eventuellement ERREUR REEL et IMAG
C               A CHAQUE PAS DE TEMPS CALCULE et VALEUR de m
C     ------------------------------------------------------
      IF( NOFOWE(1) .GT. 0 .AND. NOFOWE(2) .GT. 0 ) THEN
C        STOKAGE: Temps, Testm, Max |U|, ERREUR PR, ERREUR PI
         NBVERR = 1 + 2 + 2
      ELSE
C        STOKAGE: Temps, Testm, Max U(N)  seulement
         NBVERR = 1 + 2
      ENDIF
C     LE TEMPS en PLUS car INFO COMPLEMENTAIRE DANS VECTEUR"TESTM_ERREUR
      MOERRT = NBVERR * MXVECT * 2
      CALL TNMCDC( 'REEL', MOERRT, MNERRT )
C     TEMPS
      RMCN( MNERRT ) = TEMPS
C     TEST D'ARRET DES ITERATIONS m
      RMCN( MNERRT+1 ) = 0
C     Max Um
      RMCN( MNERRT+2 ) = REAL( MODUMX0 )
      IF( NBVERR .GT. 3 ) THEN
C        ERREUR RELATIVE PARTIE REELLE
         RMCN( MNERRT+3 ) = 0
C        ERREUR RELATIVE PARTIE IMAGINAIRE
         RMCN( MNERRT+4 ) = 0
      ENDIF
C
C     GESTION DES ADRESSES DES TABLEAUX DE L'ONDE AUX DIFFERENTS TEMPS
C     ----------------------------------------------------------------
C     PASSAGE DE UG00 DANS LE common de incl/cthet.inc
      MNTHET0 = MNUG00
C     PASSAGE DE UGtn DANS LE common de incl/cthet.inc
      MNTHETn= MNUGtn

      IF( NORESO .GE. 2 ) THEN
C
C        KG MORSE : LES 3 TABLEAUX AUXILIAIRES du SOUS PROGRAMME Gcaxbk
C        MNAUX1, MNAUX2, MNAUX3 : ADRESSES DE TABLEAUX AUXILIAIRES MCN
C        MNAUX1 SERT AUSSI DE TABLEAU POUR CALCULER LA NORME L2
         LOAUX = NBNOMA * 3
         CALL TNMCDC( 'REEL2', LOAUX, MNAUX1 )
         IF( MNAUX1 .LE. 0 ) GOTO 9999
         MNAUX2 = MNAUX1 + NBNOMA * MOREE2
         MNAUX3 = MNAUX2 + NBNOMA * MOREE2

C        LE NOMBRE D'ITERATIONS DU GC MIS A JOUR APRES la m-CONVERGENCE
         KGCUET(1) = NBNOMA
         KGCUET(2) = NBNOMA

      ENDIF

      NBPADT  = 0
      NBLOCm0 = 0
      NBCABG  = 0
      NBCAKG  = 0
      NBREDT  = 0
C
C     TEMPS DE CALCUL DE DECLARATION DE KG
      CPU1   = DINFO( 'CPU' )
      CPUPKG = CPU1 - CPU0
      CPU0   = CPU1
C
C     ===============================================================
C     CONSTRUCTION DE LA MATRICE [KG] SYMETRIQUE NON DIAGONALE
C     [KG] = [Rho/dt - LAPLACIEN/2 + N(V0**2+W0**2)]
C     RETOUR ICI EN CAS DE MODIFICATION DU PAS DE TEMPS
C     ===============================================================
C     SORTI DE THED1T LE VECTEUR UG00 EST DONNE PAR COMPOSANTES UG00(NBNOMA,2)
C     LES DL DE LA PARTIE REELLE PUIS CEUX DE LA PARTIE IMAGINAIRE
C     MNUG00 EST RANGE DANS MNTHET POUR UTILISER UG00 A PARTIR DE incl/cthet.inc
 50   MNTHET = MNUG00
C     COPIE DE UG00(NBNOMA,2) DANS UGtm(NBNOMA,2)
      CALL TRTATD( MCN(MNUG00), MCN(MNUGtm), NTDL )
      CALL NLSESIKG( PENALI, NDIM,   MNXYZN,
     %               NBTYEL, MNNPEF, NDPGST,
     %               MNTPOB, MNTAUX, NUMIOB, NUMAOB, MNDOEL,
     %               MNTHER, MNTAEL, MNX,
     %               NORESO, NCODSK, NBRDKG, MNLPLI, MNLPCO, KG,
     %               NBPTAF, IERR  )
      IF( IERR .NE. 0 ) GOTO 9999
      NBCAKG = NBCAKG + 1

      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10050) NBCAKG
10050    FORMAT('NOMBRE DE CALCUL DE LA MATRICE [KG] =',I3)
      ELSE
         WRITE(IMPRIM,20050) NBCAKG
20050    FORMAT('[KG] COMPUTATION NUMBER=',I3)
      ENDIF

      call affvect( '[KG] SANS CL', 20,     KG )
      call afl1ve(  '[KG] SANS CL', NBRDKG, KG )

C     CONSTRUCTION DE LA LISTE ET DES VALEURS DES DL FIXES DIRICHLET
C     AVEC UNE NUMEROTATION GLOBALE DES DL FX PAR NOEUDS (2,NBNOMA)
C     --------------------------------------------------------------
      CALL DLFXLI( 1,      NBNOMA, NDIM,   NBTYEL, MNNPEF, NDPGST,
     %             MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %             NBRDLX, MONDLX, MNNDLX, MNVDLX, IERR )
C
C     PRISE EN COMPTE DES DL FIXES SUR LA MATRICE KG NBNOMA x NBNOMA
C     RESTRICTION UN NOEUD DIRICHLET => PARTIES REELLE ET IMAG FIXEES
C     ---------------------------------------------------------------
      CALL DLFXKG( NBRDLX, MCN(MNNDLX), STGV, NBNOMA, MCN(MNLPLI), KG )

      IF( NORESO .GE. 2 ) THEN
C
C        CONSTRUCTION DU TABLEAU des DL FIXES 1 OU NON 0 POUR GC
         CALL TNMCDC( 'ENTIER', NTDL, MNNDPFX )
C        A PRIORI TOUS LES DL SONT LIBRES => 0
         CALL AZEROI( NTDL, MCN(MNNDPFX) )

C        MISE AU NUMERO DU DL FIXE DANS LES TABLEAUX DES DL FIXES
         DO I = 1, NBRDLX
C           LE NO GLOBAL DU DL FIXE: NDL=2*NOEUD-1 ou 2*NOEUD
            NDL = MCN(MNNDLX - 1 + I )
C
C           NUMERO DU NOEUD DU DL
            NOE = (NDL+1) / 2
C
C           NUMERO DE LA PARTIE REELLE 1 ou IMAGINAIRE 2
            IF( 2*NOE .EQ. NDL ) THEN
C              PARTIE IMAGINAIRE => AU DELA DE NBNOMA
               NDL = NBNOMA + NOE
            ELSE
C              PARTIE REELLE
               NDL = NOE
            ENDIF

C           MISE A I DU DL FIXE NDL
            MCN( MNNDPFX-1 + NDL ) = I
         ENDDO
C
      ENDIF

      call affvect( '[KG] APRES CL', 20,     KG )
      call afl1ve(  '[KG] APRES CL', NBRDKG, KG )

      IF( NORESO .EQ. 1 ) THEN
C
C        FACTORISATION DE CROUT DE [KG]=L D tL
C        =====================================
         CALL CRMC1D( MCN(MNLPLI), KG, NBNOMA, 1E-10, 1, KG, IERR )
C        IERR : 0 SI   AUCUN     PIVOT<EPS
C               1 SI AU MOINS UN PIVOT<EPS
C        ON ENTEND PAR PIVOT (A0(IDIAGONAL)-SA)/A0(IDIAGONAL)
         IF( IERR .NE. 0 ) THEN
C           MATRICE NON INVERSIBLE
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

      ENDIF
C
C     SOMME DES TEMPS DE CALCUL DE CONSTRUCTION DE LA MATRICE KG
      CPU1   = DINFO( 'CPU' )
      CPUPKG = CPUPKG + CPU1  - CPU0
      CPU0   = CPU1
C
C     ICI MNUGtm CONTIENT UG00(NBNOMA,2) L'ONDE A L'INSTANT INITIAL
C
C     KG = [K(Rho,dt,Alfa,g,N(V0**2+W0**2)] MATRICE GLOBALE NBNOMA x NBNOMA
C
C     ###################################################################
C     ##                                                               ##
C     ##  LA BOUCLE EN TEMPS AVEC DES PAS DE TEMPS CONSTANTS = DeltaT  ##
C     ##                                                               ##
C     ###################################################################
C
C     ICI LE NOUVEAU TEMPS OU SE FAIT LE CALCUL
C     =========================================
 100  TEMPS = REAL( TEMPS + DeltaT )
      ITERm = 0

C     SOLUTION UGtm AU TEMPS tn => UGtn(NBNOMA,2) AU TEMPS tn
      CALL TRTATD( MCN(MNUGtm), MCN(MNUGtn), NTDL )
C
C     ICI LES ITERATIONS m DE POINT FIXE SONT MISES EN OEUVRE
C     =======================================================
C     UNE NOUVELLE ITERATION DE POINT FIXE
 110  ITERm = ITERm + 1
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10100) TEMPS, DeltaT, NBPADT, ITERm, NBCABG+1
      ELSE
         WRITE(IMPRIM,20100) TEMPS, DeltaT, NBPADT, ITERm, NBCABG+1
      ENDIF
10100 FORMAT(/'Au TEMPS',G14.6,'  PAS de TEMPS=',G15.7,
     %' Nb PAS TEMPS=',I10,
     %' ITERm=',I3,
     %'  Total ITERm=',I10,' -> CALCUL de l''ONDE COMPLEXE NLSE'/
     %145('*'))
20100 FORMAT(/'At TIME',G14.6,' TIME STEP=',G15.7,
     %' TIME STEP Nb=',I10,
     %' ITERm=',I3,
     %'  Total ITERm=',I10,' -> COMPUTATION of the NLSE COMPLEX WAVE'/
     %145('*'))
C
C     GENERATION DU SECOND MEMBRE PARTIE RELLE BG=Fn+1,m+1 AU TEMPS tn+1
C     REAL PART=
C     Fr(tn+1) +Fi(tn+1) + Rho/dt (Vn +Wn+1m-Wn) +Alfa LAPLACIEN Wn+1m
C     +OmegaZ ( x d/dy - y d/dx )(Vn+1m-Wn+1m)
C     -N(Vn+1m**2+Wn+1m**2)(Vn+1m+Wn+1m) + N(V0**2+W0**2) Vn+1m
C
C     IMAG PART:
C    -Fr(tn+1) +Fi(tn+1) + Rho/dt (-Vn+1m+Vn +Wn) -Alfa LAPLACIEN Vn+1m
C    +OmegaZ ( x d/dy - y d/dx )(Vn+1m+Wn+1m)
C    +N(Vn+1m**2+Wn+1m**2)(Vn+1m-Wn+1m) + N(V0**2+W0**2) Wn+1m
C     ==================================================================
      MNTHET = MNUGtm
      CALL NLSESIBG( DeltaT, PENALI, NDIM,   MNXYZN, NBNOMA,
     %               NBTYEL, MNNPEF, NDPGST,
     %               MNTPOB, MNTAUX, NUMIOB, NUMAOB, MNDOEL,
     %               MNX,    MCN(MNUGtn), MCN(MNUGtm),
     %               MCN(MNBG),  IERR  )
      call afl1ve( 'ITER m BG SANS CL', NTDL, DMCN((MNBG+1)/MOREE2) )
      IF( IERR .NE. 0 ) GOTO 9999
C
C     AJOUT DES EVENTUELLES FORCES PONCTUELLES
C     ----------------------------------------
      IF( IESOPO .GT. 0 .AND. NDIM .GT. 1 ) THEN
C
C        CONSTRUCTION NO(NBNOMA,2) ET VALEUR DES FORCES PONCTUELLES
         CALL THSNFX( NTDL,   NBTYEL, MNNPEF, NDPGST,
     %                MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %                NBFNFX, MONFNX, MNNFNX, MOVFNX, MNVFNX )
C
C        ASSEMBLAGE DES FORCES PONCTUELLES(NBNOMA,2)
         CALL ASFONO( NTDL, 1, NBFNFX,MCN(MNNFNX),MCN(MNVFNX),RELMIN,
     %                MCN(MNBG) )
C
      ENDIF
C
C     PRISE EN COMPTE DES DL FIXES(2,NBNOMA) SUR LE VECTEUR GLOBAL BG(NBNOMA,2)
C     -------------------------------------------------------------------------
      CALL DLFXBG( NBRDLX, MCN(MNNDLX), MCN(MNVDLX), STGV, NBNOMA,
     %             MCN(MNBG) )

      M = MOREE2 * NBNOMA
      IF( NORESO .EQ. 1 ) THEN

C        RESOLUTION DU SYSTEME FACTORISE PAR CROUT de KG = L D tL
C        L D tL UGtm1(NBNOMA,1) = BG(*,1) => Vtn+1m+1 PARTIE REELLE DE L'ONDE
C        ====================================================================
         CALL DRCRPR( NBNOMA, NCODSK, MCN(MNLPLI), KG,  MCN(MNBG), 3,
     %                MCN(MNUGtm1), IERR )

C        RESOLUTION DU SYSTEME FACTORISE PAR CROUT de KG = L D tL
C        L D tL UGtm1(NBNOMA,2) = BG(*,2) => Wtn+1m+1 PARTIE IMAGINAIRE DE L'ONDE
C        ========================================================================
         CALL DRCRPR( NBNOMA, NCODSK, MCN(MNLPLI), KG,  MCN(MNBG+M), 3,
     %                MCN(MNUGtm1+M), IERR )

      ELSE

C        RESOLUTION DU SYSTEME PAR GRADIENT CONJUGUE SIMPLE
C        KG UGtm1(NBNOMA,1) = BG(*,1) => Vtn+1m+1 PARTIE REELLE DE L'ONDE
C        ================================================================
         CALL GCAxbk( NBNOMA,       MCN(MNNDPFX),
     %                MCN(MNLPLI),  MCN(MNLPCO),  KG,
     %                MCN(MNBG),    MCN(MNUGtm),
     %                MCN(MNAUX1),  MCN(MNAUX2),  MCN(MNAUX3),
     %                MCN(MNUGtm1), KGCUET(1) )

C        RESOLUTION DU SYSTEME PAR GRADIENT CONJUGUE SIMPLE
C        KG UGtm1(NBNOMA,2) = BG(*,2) => Vtn+1m+1 PARTIE IMAGINAIRE DE L'ONDE
C        ====================================================================
         CALL GCAxbk( NBNOMA,         MCN(MNNDPFX+NBNOMA),
     %                MCN(MNLPLI),    MCN(MNLPCO),   KG,
     %                MCN(MNBG+M),    MCN(MNUGtm+M),
     %                MCN(MNAUX1),    MCN(MNAUX2),   MCN(MNAUX3),
     %                MCN(MNUGtm1+M), KGCUET(2) )

      ENDIF

C     NETTOYAGE DES PETITES VALEURS DE UGtm1
      MN1 = ( MNUGtm1 - 1 ) / MOREE2
      DO I=1,NTDL
         D = ABS( DMCN(MN1+I) )
         IF( D .LT. UNSTGV ) THEN
            DMCN(MN1+I)=0D0
         ENDIF
      ENDDO

C     UN CALCUL DE PLUS EFFECTUE
      NBCABG = NBCABG + 1
C
C     AFFICHAGE DE L'ONDE UGtm1(NBNOMA,2) AU TEMPS tn+1m+1
C     AFFICHAGE DU MODULE DE L'ERREUR SI SOLUTION EXACTE CONNUE
C     ---------------------------------------------------------
      CALL AFNLSE( 2,   1, MNXYZN, NBNOMA, 1, MCN(MNUGtm1),
     %             WCAMIN, WCAMAX, NOFOWE, WEXMIN, WEXMAX, ERRMAX )
C     ERRMXR = MAX |RP EXACT-RP COMPUTED|(Node)/(MAX RP(Node)-MIN RP(Node))
C     ERRMXI = MAX |IP EXACT-IP COMPUTED|(Node)/(MAX IP(Node)-MIN IP(Node))
      IF( NOFOWE(1) .GT. 0 .AND. NOFOWE(2) .GT. 0 ) THEN
         ERRMXR = ERRMAX(1) / (WEXMAX(1)-WEXMIN(1))
         ERRMXI = ERRMAX(2) / (WEXMAX(2)-WEXMIN(2))
      ELSE
         ERRMXR = 0D0
         ERRMXI = 0D0
      ENDIF
C
C     AFFICHER LES NORMES ABSOLUES ET RELATIVES DE UGtm, UGtm1(2,NBNOMA)
C     ET CALCULER LA VALEUR Testm POUR TESTER LA FIN DES ITERATIONS m
C     ------------------------------------------------------------------
      CALL NLSEITER( 0, TEMPS, ITERm, NBNOMA,
     %             MCN(MNUGtm), MCN(MNUGtm), MCN(MNUGtm1), MCN(MNUGtm1),
     %             MODUMX, Testm )
C     MODUMX=Max | UG(tn+1,m+1,N) | aux NOEUDS N du MAILLAGE
C     Testm =Som||UGm+1|(N)-|UGm|(N)|/Som|UGm+1|(N) aux NOEUDS N du MAILLAGE

C     STOCKAGE TEMPS + Testm + ModuleMax + ERREUR Max DANS TMC MNERRT
C     ---------------------------------------------------------------
      IF( NBVERR*NBCABG .GE. MOERRT ) THEN
C        AUGMENTATION DU TABLEAU
         CALL TNMCAU( 'REEL', MOERRT, 2*MOERRT, MOERRT, MNERRT )
         MOERRT = 2 * MOERRT
      ENDIF

      MN = MNERRT + NBVERR * NBCABG
C     TEMPS
      RMCN( MN ) = TEMPS
C     TEST D'ARRET DES ITERATIONS m
      RMCN( MN+1 ) = REAL( Testm )
C     Max Um
      RMCN( MN+2 ) = REAL( MODUMX )

      IF( NBVERR .GT. 3 ) THEN
C        ERREUR RELATIVE sur la PARTIE REELLE
         RMCN( MN+3 ) = REAL( ERRMXR )
C        ERREUR RELATIVE sur la PARTIE IMAGINAIRE
         RMCN( MN+4 ) = REAL( ERRMXI )
      ENDIF
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

C           PAS DE TRACE DES ARETES DU MAILLAGE
            NCOAFR = -2
            NCOAPL = -2
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
               IF( MODUMX .LE. 0D0 ) MODUMX = 1D0
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

C     TEMPS DE CALCUL DE L'ITERATION m DE POINT FIXE
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
            WRITE(IMPRIM,10131) TEMPS, MODUMX
         ELSE
            WRITE(IMPRIM,20131) TEMPS, MODUMX
         ENDIF
10131    FORMAT(/100('*')/'EXPLOSION DETECTEE',
     %                    ' Max |Onde(t=',G14.6,',X)|=',G15.7/100('*'))
20131    FORMAT(/100('*')/'DETECTION of WAVE EXPLOSION'
     %                    ' Max |Wave(t=',G14.6,',X)|=',G15.7/100('*'))
C        SAUVEGARDE DES RESULTATS ACTUELS UGtm1 et SORTIE
         GOTO 9900
      ENDIF
C
 150  IF( ITERm .GE. 8 .AND. Testm .GT. 5D-2 ) THEN
C
C        LE TEST m EST TROP GRAND => REDUIRE LE PAS DE TEMPS DeltaT
C        S'IL EST ENCORE POSSIBLE DE DIVISER LE PAS DE TEMPS PAR 2
C        ----------------------------------------------------------
         IF( NBREDT .GE. MXREDT ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='STOP apres    DIVISIONS du PAS DE TEMPS'
               KERR(2)='Reduire le PAS de TEMPS INITIAL TROP GRAND?'
            ELSE
               KERR(1)='STOP after    DIVISIONS of the TIME STEP'
               KERR(2)='Reduce the INITIAL TIME STEP TOO GREAT?'
            ENDIF
            WRITE(KERR(1)(12:13),'(I2)') MXREDT
            CALL LEREUR
            IERR = 29
C           MISE DANS LE VECTEUR"ONDENLSE DES VECTEURS DEJA STOCKES
            GOTO 9900
         ENDIF
C
C        OUI: RECONSTRUCTION DE KG avec un NOUVEAU PAS DE TEMPS/2
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
     %   '  et RECONSTRUCTION [KG]')
20200    FORMAT(/'TIME STEP=',G14.6,
     %   ' is CHANGED to the NEW TIME STEP=',G14.6,
     %   '  and [KG] is CONSTRUCT AGAIN')
C        NOUVEAU PAS DE TEMPS REDUIT
         NBREDT = NBREDT + 1
         DeltaT = DeltaT / 2D0
         PasTemps = REAL( DeltaT )
         GOTO 50
C
      ENDIF

      IF( Testm .GT. 5D-4 ) THEN
C
C        NON CONVERGENCE DU POINT FIXE => UNE ITERATION DE PLUS A FAIRE
C        --------------------------------------------------------------
         IF( ITERm .GE. ITERMX ) THEN
            NBLGRC(NRERR) = 3
            WRITE(KERR(5)(1:6),'(I6)') ITERMX
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
10140 FORMAT('Au TEMPS',G14.6,' ITER m=',I3,'  CONVERGENCE des ITERATION
     %S m')
20140 FORMAT('At TIME', G14.6,' ITER m=',I3,'  m ITERATION CONVERGENCE')
C
c     LE VECTEUR DES DL DE L'ONDE ACTUELLE EST IL A STOCKER?
C     ------------------------------------------------------
      IF( TEMPS .GE. TSTOC*0.9999 .AND. NBVECT .LT. MXVECT ) THEN

C        OUI: STOCKAGE DE L'ONDE UGtm1(NBNOMA,2) A CET INSTANT TEMPS
         MNTEMP = MNTEMP + NTDL * MOREE2
         CALL TRTATD( MCN(MNUGtm1), MCN(MNTEMP), NTDL )
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
10150 FORMAT('Au TEMPS',G14.6,' STOCKAGE de l''ONDE',I5,
     %' de',I9,' DEGRES de LIBERTE')
20150 FORMAT('At TIME',G14.6,' STORAGE of the WAVE',I5,
     %' of',I9,' DEGREES of FREEDOM')
C
C        LE PROCHAIN TEMPS DE STOCKAGE
         IF( LANGAG .EQ. 0 ) THEN
            print*,'nlsesimpl: temps=',TEMPS,
     %       ' temps0=',TPSINI,
     %       ' temps fin=',TPSFIN,
     %       ' Temps STOCKAGE=',TSTOC,
     %       ' Pas Temps STOCKAGE=',DTSTOC
         ELSE
            print*,'nlsesimpl: Time=',TEMPS,
     %       ' INITIAL Time=',TPSINI,
     %       ' FINAL Time=',TPSFIN,
     %       ' Storage Time=',TSTOC,
     %       ' Storage Time Step=',DTSTOC
         ENDIF

C        PROCHAIN TEMPS DE STOCKAGE DE L'ONDE
         TSTOC = TSTOC + DTSTOC

      ENDIF
C
C     ETUDE EN TEMPS TERMINEE?
C     ------------------------
      IF( TEMPS + RDeltaT .LT. TPSFIN*1.00001 ) THEN

C        NON: UN PAS DE TEMPS DE PLUS EST CALCULE

         IF( Testm .LT. 5D-4 .AND. ITERm .LE. 1 ) THEN
C
C           POUR ACCELERER: RECONSTRUCTION DE KG avec un PAS DE TEMPS*1.5
C           =============================================================
C           LE NOUVEAU TEMPS INITIAL EST tn
            TEMPSINI = TEMPS
C           SOLUTION INITIALE REMPLACEE PAR LA DERNIERE SOLUTION CALCULEE
            CALL TRTATD( MCN(MNUGtm1), MCN(MNUG00), NTDL )
C
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10200) DeltaT, DeltaT*1.5D0
            ELSE
               WRITE(IMPRIM,20200) DeltaT, DeltaT*1.5D0
            ENDIF
C           NOUVEAU PAS DE TEMPS
            NBREDT = NBREDT - 1
            DeltaT = DeltaT * 1.5D0
            PasTemps = REAL( DeltaT )
            GOTO 50
C
         ENDIF
C
         NBLOCm = NBPADT / 50
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
C     STOCKAGE DE L'ONDE UGtm1(NBNOMA,2) AU DERNIER INSTANT TEMPS
 9900 IF( NBVECT .LT. MXVECT ) THEN
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
C
C     LES TEMPERATURES A L'INSTANT INITIAL SONT A L'ADRESSE MNTEMP
      MNTEMP = MNVECT + WECTEU
C     L'ADRESSE -1 DU PREMIER TEMPS STOCKE DERRRIERE LES VECTEURS TEMPERATURE
      MNTIME = MNTEMP + NTDL * NBVECT * MOREE2 -1

C     CONSTRUCTION DE VECTEUR"TESTM_ERREUR A PARTIR DE L'EVENTUEL ANCIEN
C     AVEC AJOUT DU TMC MNERRT ACTUEL
C     ==================================================================
      CALL TITSMXERR( KNOMOB, MOREE2, TPSINI, NBVERR-1, NBCABG, MNERRT )
C
C     DESTRUCTION DES TMC DEVENUS INUTILES
C     ====================================
 9999 IF( IERKGALLOC .EQ. 0 ) DEALLOCATE( KG )
      IF( MNAUX1 .GT. 0 ) CALL TNMCDS( 'REEL2',  LOAUX,   MNAUX1 )
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
      IF( NORESO .EQ. 1 ) THEN
C        CROUT PROFIL
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,19002) NTDL*4 + NBNOMA + MOREE2 * NBRDKG * 2
         ELSE
            WRITE(IMPRIM,29002) NTDL*4 + NBNOMA + MOREE2 * NBRDKG * 2
         ENDIF
19002 FORMAT(/'NLSESIMPL: STOCKAGE MATRICE PROFIL + VECTEURS =',I15,
     %' MOTS'/ )
29002 FORMAT(/'NLSESIMPL: SKYLINE MATRIX + VECTORS STORAGE =',I15,
     %' MEMORY WORDS')
      ELSE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,19003) NTDL*5 + NBNOMA*3 + MOREE2 * NBRDKG * 2
         ELSE
            WRITE(IMPRIM,29003) NTDL*5 + NBNOMA*3 + MOREE2 * NBRDKG * 2
         ENDIF
19003 FORMAT(/'NLSESIMPL: STOCKAGE MATRICE MORSE + VECTEURS=',I15,
     %' MOTS MEMOIRE'/ )
29003 FORMAT(/'NLSESIMPL: CONDENSED MATRIX + VECTORS STORAGE=',I15,
     %' MEMORY WORDS')
      ENDIF
C
C     AFFICHAGE DES TEMPS CALCUL
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,19000) CPUPKG, CPUITE, CPUPKG+CPUITE
      ELSE
         WRITE(IMPRIM,29000) CPUPKG, CPUITE, CPUPKG+CPUITE
      ENDIF
19000 FORMAT(/
     %'CPU TRAITEMENT DE LA MATRICE KG   =',F12.2,' SECONDES CPU'/
     %'CPU ITERATIONS TEMPS + PT FIXE    =',F12.2,' SECONDES CPU'/
     %'CPU TOTAL DE CALCUL DE LA SOLUTION=',F12.2,' SECONDES CPU')

29000 FORMAT(/
     %'CPU of KG MATRICE CONSTRUCTION =',F12.2,' CPU SECONDS'/
     %'CPU of TIME + FIX PT ITERATIONS=',F12.2,' CPU SECONDS'/
     %'CPU TOTAL SOLUTION COMPUTATION =',F12.2,' CPU SECONDS')

      RETURN
      END
