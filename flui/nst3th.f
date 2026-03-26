      SUBROUTINE NST3TH( KNOMOB, NAVSTO, NTLXOB, MOREE2,RELMIN, KNOMFIC,
     %                   NDIM,   MNXYZN, NBSOM,  NBNOVI,
     %                   NBTYEL, MNNPEF, NUTYEL, MNTPOB, NBDL1EF, MOAUX,
     %                   NUMIOB, NUMAOB, MNDOEL, MNDTEL, NDPGST, IEBLPR,
     %                   NBTEFX, MONTEFX, MNNTEFX, MNVTEFX,
     %                   NORESO, LP2LIGN, LP2COLO, NBCMVG,
     %                   DT,     DTSTOC,  TPSINI, TPSFIN, MNTIMES,
     %                   NBPAST, NOVVIPR, NTDLVP, NTDLTE,
     %                   NDDLNO, VXYZPN,  VITMAX, VITMOY, TEMPER0,
     %                   N1VPMIMX,MXPAST, VPMIMX,
     %                   DFABG,  DFACTO,  DVITPR, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA VITESSE+PRESSION+TEMPERATURE DANS UN FLUIDE 2D ou 3D
C ----- DU PROBLEME DE NAVIER-STOKES + THERMIQUE EN INSTATIONNAIRE
C       SELON L'APPROXIMATION DE BOUSSINESQ

C       https://fr.wikipedia.org/wiki/Convection
Cf      https://en.wikipedia.org/wiki/Boussinesq_approximation_(buoyancy)
Cf      C. BERNARDI, B. METIVET, B.PERNAUD-THOMAS
C       Couplage des equtions de Navier-Stokes et de la chaleur:
C       Le modele et son approximation par elements finis
C       RAIRO Modelisation Mathematique et Analyse NUmerique
C       tome 29,no7(1995), p. 871-921

C       Les equations traitees avec u la vitesse, T la temperature et
C       la pression statique locale en deux termes : P = Ph + Pd
C       ou Pd est la pression dynamique due au mouvement (convection) du fluide
C          Ph designe la pression hydrostatique (-Rho g Z) (zmax-z)
C             soit encore grad Ph = -Rho0 g Z
C          Z  le vecteur de la direction Z (opposee a la pesanteur)
C          FGth les sources internes de chaleur
C          FGU  les forces  internes au fluide

C       Rho Cp( dT/dt+(u.D)T ) - div( Conduc grad T ) = FGth
C       Rho   ( du/dt+(u.D)u ) - div( Mhu grad u ) + grad Pd =
C                              - CoBous Rho (T-T0) g (-Z) + FGu
C       div u = 0

C**************************************************************************

C       NAVSTO=3: Discretisation avec la methode des caracteristiques retro
C       A L'ETAPE n+1, LE PROBLEME CONSISTE A TROUVER
C                      {u(tn+1), p(tn+1), T(tn+1)=te(tn+1)} SOLUTION de

C      1. Calcul de la temperature {te(tn+1,x)} :
C      ( Rho Cp - dt Conduc Laplacien ) {te(tn+1,x)} =
C      soit
C       = Rho Cp {te(X(tn;tn+1,x)} + dt {FGte(tn+1)})
C       ou X(tn;tn+1,x) est la position de la molecule du fluide a tn
C            qui sera en x a tn+1 (methode des caracteristiques en tn+1)
C      soit
C       = Rho Cp {te(tn,x) - Rho Cp dt {Som  Ui(tn,x) . d/dxi te(tn,x)}
c                                     i=1,NDIM   (MOINS IMPLICITE car en tn)
C        + dt {FGte(tn+1)})

C      Condition aux limites sur la temperature te
C       te(tn+1,X) = teG(tn+1,X)  avec X sur GammaD
C       - Conduc d te(tn+1,x)/dn = CoefEch ( te(tn+1,x) - teExt(tn+1,x) )
C       condition de Fourier-Robin avec le coefficient d'echange CoefEch
C       teExt(t,x) la temperature exterieure a l'instant t et en x
C                  sur la frontiere GammaF
C     ...................................................................

C      2. Calcul des NDIM composantes de la vitesse intermediaire {U*i(tn+1,x)}
C      ( Rho - dt Mhu Laplacien ) {U*i(tn+1,x)} =
C            Rho {Ui(tn,X(tn;tn+1,x))}
C           - dt CoGrPr gradi p(tn,x)
C           - dt Rho CoBous g (-Z) {te(tn+1)-te(t0)}
C           + dt {FGui(tn+1)

C       ou X(tn;tn+1,x) est la position de la molecule du fluide a tn
C                       qui sera en x a tn+1 (methode des caracteristiques)
C          g l'acceleration de la pesanteur est ici prise en compte
C            dans la direction opposee a l'axe Z car c'est la force
C            qui fait apparaitre la CONVECTION NATURELLE du fluide
C          CoBous est le coefficient de DILATATION THERMIQUE
C                     du fluide (THERMAL EXPANSION COEFFICIENT)
C          Cf https://www.comsol.fr/multiphysics/boussinesq-approximation

C      Condition aux limites sur U*i(tn+1)
C       U*i(tn+1,X) = Ui(tn+1,X)        avec X sur GammaU
C       U*i(tn+1,X) = Ui(X(tn;tn+1,x))  avec X sur GammaU de sortie de fluide
C       -dt Mhu dU*i(tn+1,X)/dn = 0     avec X sur Frontiere-GammaU 
C                               pour i=1,...,NDIM
C     .......................................................................

C      3. Calcul de la difference des pressions {P(tn+1)-P(tn)}
C     -dt CoGrPr Laplacien {P(tn+1)-P(tn)} =
C       - (Rho -dt Mhu Laplacien) (Div{U*(tn+1)} -Integrale Div U*(tn+1) /Vol) }
C     Condition aux limites
C       Pn+1,m+1(Gamma) = P(tn+1,Gamma) pour Gamma sur la frontiere GammaP
C      -dt CoGrPr dPn+1,m+1/dn = -n . (Rho -dt Mhu Laplacien) {U*(tn+1)}
C     Si dt Mhu est petit devant Rho, le terme dt Mhu Laplacien est neglige
C     .......................................................................

C      4. Calcul des NDIM composantes de la difference des vitesses {U(tn+1)-U*}
C     ( Rho - dt Mhu Laplacien ) {U(tn+1)-U*}i =
C                                 - dt CoGrPr d/dxi (P(tn+1)-P(tn))
C       Condition aux limites
C       Ui(tn+1,X)-U* = 0           avec X sur GammaU
C       -dt Mhu dUi(tn+1,X)/dn = 0  avec X sur Frontiere-Gamma 
C                              pour i=1,...,NDIM
C     .......................................................................

C      5. Calcul des NDIM composantes de la vitesse {U(tn+1)}
C      U(tn+1,Noeuds) = U(tn,Noeuds) + (U(tn+1) - U(tn))(Noeuds)
C     .......................................................................

C      6. Calcul de la pression {P(tn+1)}
C      P(tn+1,Sommets) = P(tn,Sommets) + ( P(tn+1,Sommets)-P(tn,Sommets) )
C     .......................................................................


C       A L'ETAPE t0 initiale IL FAUT RECUPERER U(t0) P(t0) Te(t0)
C       puis, les iterations en temps peuvent demarrer
C     .......................................................................

C       QUELQUES RESTRICTIONS du PROGRAMME nst3th:

C       LE PAS de TEMPS EST CONSTANT
C       La direction de la pesanteur est en 3d selon l'axe -Z  (en 2d -Y)
C       LA METHODE DES CARACTERISTIQUES RETROGRADES PERMET
C       LE CALCUL des TERMES de TRANSPORT
C       LE MAILLAGE est FORME de TRIANGLES ou TETRAEDRES P2 en TEMPERATURE
C       de TAYLOR-HOOD P2 EN VITESSE et P1 EN PRESSION
C
C       DENSITE DE MASSE Rho, VISCOSITE Mhu SONT INDEPENDANTS DU TEMPS
C       DENSITE DE MASSE Rho, VISCOSITE DYNAMIQUE Mhu,
C       COEFFICIENT de la PRESSION CoGrPr sont RETROUVES AU BARYCENTRE
C       DE L'EF 1

C       BLOCAGE DE LA TEMPERATURE IMPOSEE
C       FLUX THERMIQUE SUR LA FRONTIERE
C       FORCE EXTERNE  SUR LA FRONTIERE
C       BLOCAGE DE LA VITESSE  IMPOSEE
C       BLOCAGE DE LA PRESSION IMPOSEE   PEUVENT DEPENDRE DU TEMPS
C       mais TOUTES LES NDIM COMPOSANTES DE LA VITESSE EN UN NOEUD
C       DOIVENT ETRE FIXEES mais avec DES VALEURS NON FORCEMENT EGALES
C       et QUI PEUVENT DEPENDRE DU TEMPS

C***************************************************************************

C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET A TRAITER
C NAVSTO : =3 RESOLUTION du PROBLEME DE NAVIER-STOKES+THERMIQUE BOUSSINESQ
C NTLXOB : NUMERO DU TMS DU LEXIQUE DE L'OBJET KNOMOB
C MOREE2 : NOMBRE DE MOTS D'UNE VARIABLE REELLE DOUBLE PRECISION
C RELMIN : PLUS PETIT REEL SERVANT DE MARQUEUR DE NON UTILISATION
C KNOMFIC : NOM DU FICHIER SUPPORT DU VECTEUR VITESSE+PRESSION

C NDIM   : DIMENSION DES COORDONNEES DES POINTS ( 2 OU 3 )
C MNXYZN : ADRESSE MCN DU TABLEAU XYZNOEUD DE L'OBJET KNOMOB
C          XYZ des SOMMETS+BARYCENTRES    DES EF BREZZI-FORTIN
C          XYZ des SOMMETS+MILIEUX ARETES DES EF TAYLOR-HOOD
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C          SUPPORTS DES DL DE LA PRESSION
C NBNOVI : NOMBRE DE NOEUDS DU MAILLAGE (SOMMETS ET MILIEUX DES ARETES)
C          SUPPORTS DES DL DE LA TEMPERATURE et DES VITESSES

C NBTYEL : NOMBRE DE TYPES D'EF DU MAILLAGE DE CET OBJET
C MNNPEF : ADRESSE MCN DU TABLEAU DES ADRESSES MCN DES TMS NPEF"TYPE EF
C NUTYEL : NUMERO TYPE D'EF TAYLOR-HOOD 15 (2d) ou 20 (3d)
C MNTPOB : ADRESSE MCN DU TABLEAU POINTEUR SUR LES TABLEAUX POBA DES EF
C NBDL1EF: NOMBRE DE DL D'UN EF DU FLUIDE
C MOAUX  : NOMBRE MAXIMAL DE MOTS POUR LES TABLEAUX AUXILIAIRES DES EF

C NUMIOB : NUMERO MINIMAL DU PLSV DANS LA DEFINITION DE L'OBJET
C NUMAOB : NUMERO MAXIMAL DU PLSV DANS LA DEFINITION DE L'OBJET
C MNDOEL : LES 4 ADRESSES MCN DES TABLEAUX DES ADRESSES DES
C          TABLEAUX DECRIVANT LES DONNEES FLUIDE DU FLUIDE DE L'OBJET COMPLET
C MNDTEL : LES 4 ADRESSES MCN DES TABLEAUX DES ADRESSES DES
C          TABLEAUX DECRIVANT LES DONNEES THERMIQUES DU FLUIDE DE L'OBJET

C NDPGST : CODE TRAITEMENT DES SOMMETS POINTS NOEUDS DU MAILLAGE
C          1: noeuds = points # sommets les tms XYZSOMMET XYZNOEUD existent
C IEBLPR : NOMBRE DE TMS BLPRESSION DES PLS DE L'OBJET TROUVE

C NBTEFX : NOMBRE DE NOEUDS DE TEMPERATURE FIXEE calcule dans tempert0.f
C MONTEFX: NOMBRE DE MOTS DECLARES DU TABLEAU MC NO DES TEMPERATURE FIXEES
C MNNTEFX: ADRESSE MCN DU TABLEAU MC DES NUMEROS DES TEMPERATURES FIXEES
C          =0 SINON
C MNVTEFX: ADRESSE MCN DU TABLEAU MC DES VALEURS DES TEMPERATURES FIXEES
C          =0 SINON

C NORESO : CODE RESOLUTION DES SYSTEMES LINEAIRES AVEC LES MATRICES VG TG et PG
C          1 FACTORISATION COMPLETE DE CROUT ET MATRICE PROFIL VG TG et PG
C          2 GRADIENT CONJUGUE SIMPLE avec STOCKAGE MORSE DE VG TG PG

C LP2LIGN: TABLEAU des POINTEURS SUR LES COEFFICIENTS DIAGONAUX
C          DE LA MATRICE MORSE VG et TG
C LP2COLO: NO DES COLONNES DES COEFFICIENTS STOCKES DE CHAQUE LIGNE
C          DE LA MATRICE MORSE VG et TG NON DIAGONALES
C NBCMVG : NOMBRE DE COEFFICIENTS STOCKES DE LA MATRICE GLOBALE VG et TG

C DT     : PAS CONSTANT DU TEMPS
C DTSTOC : PAS CONSTANT DU TEMPS ENTRE 2 STOCKAGES DU VECTEUR"VITESSEPRESSION
C TPSINI : TEMPS INITIAL DU CALCUL
C TPSFIN : TEMPS FINAL   DU CALCUL
C MNTIMES: ADRESSE MCN DES TEMPS DES VECTEURS VITESSE+PRESSION

C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTE DU FLUIDE
C          DL VITESSES (P2 INTERPOLATION) et PRESSION (P1 INTERPOLATION)
C NTDLTE : NOMBRE TOTAL DE DEGRES DE LIBERTE TEMPERATURE (NOEUDS P2)
C          =NBNOVI

C NDDLNO : NUMERO DU DERNIER DL DE CHAQUE NOEUD VITESSE (0:NBNOVI)
C          A UTILISER POUR RETROUVER LES DL PRESSION DES VECTEURS VXYZP
C VXYZPN : TABLEAU VITESSEXYZ+PRESSION AUX NOEUDS A L'INSTANT INITIAL
C TEMPER0: TABLEAU TEMPERATURE AUX  NBNOVI NOEUDS A L'INSTANT INITIAL

C N1VPMIMX:=11 NOMBRE DE VALEURS STOCKES DANS VPMIMX POUR CHAQUE PAS DE TEMPS
C MXPAST : NOMBRE MAXIMAL DE PAS DE TEMPS CALCULABLE DANS VPMIMX

C MODIFIES:
C --------
C NBPAST : NOMBRE DE PAS DE TEMPS D'INTEGRATION des EQUATIONS de NAVIER-STOKES
C NOVVIPR: NUMERO DU DERNIER VECTEUR VITESSEPRESSION STOCKE
C VITMAX : NORME MAXIMALE INITIALE PUIS FINALE DE LA VITESSE AUX NOEUDS
C VITMOY : NORME MOYENNE  INITIALE PUIS FINALE DE LA VITESSE AUX NOEUDS

C SORTIES:
C --------
C VXYZPN : TABLEAU VITESSEXYZ+PRESSION AUX NOEUDS A L'INSTANT FINAL
C NBPAST : NOMBRE DE PAS DE TEMPS CALCULES ET DE RESULTATS STOCKES DANS VPMIMX
C VPMIMX : (TEMPS, VitesseMoyenne, VitesseMax,
C           PressionMoyenne, Pression Max-Min,
C           FLUX-, FLUX+, No du Pas TEMPS,
C           TemperatureMoyenne, TemperatureMin, TemperatureMax )
C           a chaque temps calcule

C DFABG  : TEMPS CPU DES FORMATION      DES MATRICES GLOBALES
C DFACTO : TEMPS CPU DES FACTORISATIONS DES MATRICES GLOBALES
C DVITPR : TEMPS CPU DES ITERATIONS DE CALCUL DE LA VITESSE PRESSION TEMPERATURE

C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY          Janvier 2023
C23456---------------------------------------------------------------012
C$    use OMP_LIB

C     VALEUR A METTRE SUR LA DIAGONALE DE LA MATRICE PROFIL
C     POUR TRAITER LES DL FIXES
      DOUBLE PRECISION  TGV
      PARAMETER        (TGV = 1D30)

C     ACCELERATION DE LA PESANTEUR TERRESTRE en M/S/S
      DOUBLE PRECISION   G
      PARAMETER        ( G = 9.8D0 )

      DOUBLE PRECISION   D2PI
C     NOMBRE 2 x PI DANS UNE VARIABLE REELLE DOUBLE PRECISION
      PARAMETER        ( D2PI = ATAN(1D0) * 8D0 )

      include"./incl/lu.inc"
      include"./incl/xvfontes.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/nmproj.inc"
      include"./incl/donflu.inc"
      include"./incl/donthe.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a___arete.inc"
      include"./incl/a___face.inc"
      include"./incl/a___aretefr.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___profil.inc"
      include"./incl/a___morse.inc"
      include"./incl/a___temperinit.inc"
      include"./incl/msvaau.inc"
      include"./incl/ponoel.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
      include"./incl/nctyef.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/traaxe.inc"
      include"./incl/xyzext.inc"
      include"./incl/p2p22d.inc"
      include"./incl/p2p23d.inc"
      include"./incl/threads.inc"
      include"./incl/pp.inc"

      COMMON             MCN(MOTMCN)
      REAL              RMCN(MOTMCN)
      DOUBLE PRECISION  DMCN(MOTMCN/2)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)

      CHARACTER*(*)     KNOMOB, KNOMFIC
      CHARACTER*4       NOMELE(2)
      INTEGER           NUMIOB(4), NUMAOB(4), MNDOEL(4), MNDTEL(4)
      INTEGER           NOOBVC(1), NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NDDLNO(0:*), LP2LIGN(0:NBNOVI), LP2COLO(NBCMVG)
      INTEGER           NOSOEF(4), NONOEF(10)
      INTEGER           KGCTEMP, KGCUET(3), KGCPRE, KGCUMU(3)
      REAL              DT, DTSTOC, TPSINI, TPSFIN, XYZEF(30)

      DOUBLE PRECISION  VPMIMX(N1VPMIMX,0:MXPAST)
      DOUBLE PRECISION  RELMIN, VXYZPN(NTDLVP), TEMPER0(NBNOVI)
      DOUBLE PRECISION  Rho, Mhu,   CoGrPr,
     %                  Cp,  RhoCp, Conduc,
     %                  CoBOUS, CoBOUSfl, CoBOUSth,
     %                  VOLUME, INTDIVV,INTPRES
      DOUBLE PRECISION  VITMIN, VITMAX0,VITMAX,  VITMAXn,
     %                  VITMOY0,VITMOY, VITMOYn, VITFXMAX,
     %                  PREMIN, PREMAX, PREMOY,
     %                  TEMMIN, TEMMAX, TEMMOY, TEMMINt0, TEMMAXt0,
     %                  TEMMINFX, TEMMAXFX, TE, TEMMINECR, TEMMAXECR,
     %                  FLUNEG, FLUPOS, SEMINVIT, SEMAXVIT
      DOUBLE PRECISION  VITECR
ccc      DOUBLE PRECISION  dtMhu
      DOUBLE PRECISION  S, STGV, DELTAT, Pi, PENALI
      DOUBLE PRECISION  DINFO, DFABG, DFACTO, DVITPR, DMOTSTO
      DOUBLE PRECISION  DATE00, DATE0, DATE, SECONDES, TIMEMOY1DT
      DOUBLE PRECISION  t_cpu_0, t_cpu_1, t_cpu_it0, t_cpu_it1
      integer           t0, t1, nbclockps
      DOUBLE PRECISION  tt0

      INTRINSIC         DBLE, ABS

      DOUBLE PRECISION, allocatable, dimension(:,:) :: VXVYVZPR
      DOUBLE PRECISION, allocatable, dimension(:)   :: PG, VG, TG, V3P2,
     %                  PVAUX, Pr1Pr, DAUX1, DAUX2, DAUX3, VITC
 
      INTEGER,          allocatable, dimension(:) :: NO1EFN,  NONOSO,
     %                  NODLTEFX, NODLPRFX, NODLVIFX, LP1LIGN, LP1COLO

      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ), dimension(:), allocatable :: temp

      PRINT *
      PRINT *,'Start of nst3th IMPLICIT HEAT NAVIER-STOKES SOLVER on the
     % OBJECT ',KNOMOB

C     LA DATE EN SECONDES DEPUIS LE 1/1/70 MINUIT
      call secondes1970( DATE00 )
      call system_clock(count=t0, count_rate=nbclockps)

C     TEMPS CPU donne par le systeme
      call cpu_time(t_cpu_0)

C     Le TEMPS CPU donne par le systeme a la fin de l'iteration 0
      call cpu_time(t_cpu_it0)

C///////////////////////////////////////////////////////////////////////
C$OMP PARALLEL SHARED( tt0 )
      tt0 = OMP_GET_WTIME()
      PRINT*,'nst3th.f: OMP_GET_NUM_THREADS=',NBTHREADS
C$OMP END PARALLEL
C///////////////////////////////////////////////////////////////////////
      TIMEMOY1DT = 0D0
      Pi = ATAN( 1D0 ) * 4D0

      IERNO1EFN  = 1
      IERNONOSO  = 1
      IERNODLTEFX  = 1
      IERNODLPRFX  = 1
      IERNODLVIFX  = 1
      IERPGALLOC = 1
      IERVGALLOC = 1
      IERTGALLOC = 1
      IERVXVYVZPR= 1
      IERV3P2    = 1
      IERVITC    = 1
      IERDAUX1   = 1
      IERDAUX2   = 1
      IERDAUX3   = 1
      IERPr1Pr   = 1
      IERPVAUX   = 1
      IERLP1LIGN = 1
      IERLP1COLO = 1

      NBPRFX  = 0
      MNNPRFX = 0
      MNVPRFX = 0

      NBVCFX  = 0
      NBVIFX  = 0
      MNNVIFX = 0
      MNVVIFX = 0

      MNFVSF  = 0
      MOSFOB  = 0
      NBCMPG  = 0

      MNVVEF  = 0
      MNSFEF  = 0
      MNLAEF  = 0
      MNPSEF  = 0

      MOARET  = 0
      MXARET  = 0
      MNLARE  = 0
      MOFACE  = 0
      MXFACE  = 0
      MNLFAC  = 0

      MNTAUX = 0
      MNTAEL = 0
      MNTHER = 0
      MNNODL = 0
      MONODL = 0
      MNX    = 0
      MNBG   = 0
      MNBGD  = 0
      MNUTE0 = 0
      NBCHTROL= 0

C     ADRESSE MCN DE TABLEAUX DANS cthet.inc
ccc      MNTHET0= 0
      MNTHETn= 0
      MNTHET = 0

C     DECLARATION de temp POUR LE TRACE DES TEMPERATURES
      IALtemp = 1
      allocate( temp( 1:1 ), STAT=IALtemp )
      IF(IALtemp .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'ERREUR en ALLOCATION de temp(1:1)'
         ELSE
            PRINT*,'ALLOCATION ERROR of temp(1:1)'
         ENDIF
      ELSE
         NULLIFY( temp(1)%dptab )
      ENDIF

C     LE PAS DE TEMPS EN DOUBLE PRECISION
      DELTAT = DT

C     MODE DE PRISE EN COMPTE DES CONDITIONS AUX LIMITES DE DIRICHLET
C     SUR LES SYSTEMES LINEAIRES avec LES MATRICES VG PG TG
      IF( NORESO .EQ. 1 ) THEN
C        MATRICE PROFIL
         STGV = TGV
      ELSE
C        MATRICE MORSE
         STGV = 1D0
      ENDIF

C     PARAMETRES POUR LA FACTORISATION DE CROUT A=L D tL
C     NENTRE : =0 RETOUR AU PROGRAMME APPELANT SI ABS(PIVOT)<EPS
C              =1 LES CALCULS SE POURSUIVENT SAUF SI PIVOT=0
      NENTRE  = 1
C     SEUIL DES FACTORISATIONS CORRECTES L D tL DES MATRICES PROFIL
      EPSCROUT = 0.0

      NOFONT0 = NOFONT

C     LA VITESSE MAXIMALE ET MOYENNE INITIALE
      VITMAX0 = VITMAX
      VITMOY0 = VITMOY

C     RECUPERATION DU MAILLAGE EN ELEMENTS FINIS P2 DE TAYLOR-HOOD
C     ============================================================
C     MNELE : ADRESSE DU TABLEAU NPEF"TYPE EF (TRIA 2P2C ou TETR 3P2C)
      MNELE  = MCN( MNNPEF )

C     NOMBRE D'ELEMENTS FINIS DE CE TYPE
      NBELEM = MCN( MNELE + WBELEM )

C     LE NUMERO DU TYPE DE L'ELEMENT FINI (TRIA 2P2C ou TETR 3P2C)
      NUTYEL = MCN( MNELE + WUTYEL )

C     MNPGEL ADRESSE MCN DES NUMEROS NOEUDS ET POINTS GEOMETRIQUES DES EF
      MNPGEL = MNELE + WUNDEL

C     LES CARACTERISTIQUES DE L'ELEMENT FINI
      CALL ELNUNM( NUTYEL, NOMELE )
      CALL ELTYCA( NUTYEL )
      IF( NUTYEL .EQ. 15 ) THEN

C        TRIANGLE DE TAYLOR-HOOD
         NDIM   = 2
C        NOMBRE DE NO DE NOEUDS D'UN EF DANS NPEF
         NBNOEF = 6
C        NOMBRE DE SOMMETS DE L'EF
         NBPSEF = 3
C        NOMBRE D'ARETES D'UN EF
         NBLAEF = 3
C        NOMBRE DE FACES D'UN EF
         NBSFEF = 1
C        NOMBRE DE VOLUME D'UN EF
         NBVVEF = 0

      ELSE IF( NUTYEL .EQ. 20 ) THEN

C        TETRAEDRE DE TAYLOR-HOOD
         NDIM   = 3
C        NOMBRE DE NO DE NOEUDS D'UN EF DANS NPEF
         NBNOEF = 10
C        NOMBRE DE SOMMETS DE L'EF
         NBPSEF = 4
C        NOMBRE D'ARETES D'UN EF
         NBLAEF = 6
C        NOMBRE DE FACES D'UN EF
         NBSFEF = 4
C        NOMBRE DE VOLUME D'UN EF
         NBVVEF = 1

      ELSE

C        MAILLAGE D'EF NON EF DE TAYLOR-HOOD
         IERR = 1
         IF( LANGAG .EQ. 0 ) THEN
            PRINT *,'nst3th: MAILLAGE SANS ELEMENT FINI P2 TAYLOR-HOOD'
         ELSE
            PRINT *,'nst3th: MESH WITHOUT P2 TAYLOR-HOOD FINITE ELEMENT'
         ENDIF
         GOTO 9999

      ENDIF

C     NOMBRE TOTAL DE DL DE LA VITESSE SUR LE MAILLAGE
      NTDLVI = NDIM * NBNOVI

C     NBSOM EST LE NOMBRE DE DL DE LA PRESSION EGAL AU NOMBRE DE SOMMETS
C     NBSOM = NBSOM

C     NOMBRE TOTAL DE DL VITESSE  (P2) AUX NOEUDS  DU MAILLAGE
C                      + PRESSION (P1) AUX SOMMETS DU MAILLAGE
      NTDLVP = NTDLVI + NBSOM

C     NOMBRE TOTAL DE DL TEMPERATURE (P2) AUX NOEUDS DU MAILLAGE
      NTDLTE = NBNOVI

      IF( NDIM .EQ. 2 ) THEN

C        TRIANGLES: RECUPERATION DU TMS DES ARETES DE L'OBJET 2D
c        CONSTRUCTION IMPOSEE DU TMS ARETE
C        POUR REMONTER LES CARACTERISTIQUES
C        -------------------------------------------------------
         CALL LXTSOU( NTLXOB, 'ARETE', NTAROB, MNAROB )
         IF( NTAROB .GT. 0 ) CALL LXTSDS( NTLXOB, 'ARETE' )
C
C        LE TABLEAU N'EXISTE PAS => IL EST CREE
C        CALCUL PAR HACHAGE DES ARETES DE L'OBJET A PARTIR DE TOPO+NPEF"...
C        NECESSAIRE POUR CONNAITRE LES EF ADJACENTS PAR ARETE
         CALL HAC2AF( KNOMOB, 3, NTAROB, MNAROB, IERR )
         IF( IERR .NE. 0 .OR. NTAROB .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'OBJET SANS TABLEAU DES ARETES'
            ELSE
               KERR(1) = 'OBJECT WITHOUT THE ARRAY OF EDGES'
            ENDIF
            CALL LEREUR
            GOTO 9999
         ENDIF
C        LE NOMBRE D'ENTIERS PAR ARETE
         MOARET = MCN( MNAROB + WOARET )
C        LA MAJORATION DU NOMBRE D'ARETES
         MXARET = MCN( MNAROB + WXARET )
C        LE NOMBRE D'ARETES FRONTALIERES NON SUR LIGNES UTILISATEUR
C        NBARFB = MCN( MNAROB + WBARFB )
C        LE NOMBRE D'ARETES INTERFACES   NON SUR LIGNES UTILISATEUR
C        NBARIN = MCN( MNAROB + WBARIN )
C        LE NUMERO MINIMAL DE LIGNE DE L'OBJET
         NUMILF = MCN( MNAROB + WUMILF )
C        LE NUMERO MAXIMAL DE LIGNE DE L'OBJET
         NUMXLF = MCN( MNAROB + WUMXLF )
C        LE NUMERO DE LA PREMIERE ARETE FRONTALIERE NON SUR LIGNES DE L'OBJET
C        L1ARFB = MCN( MNAROB + W1ARFB )
C        LE NUMERO DE LA PREMIERE ARETE INTERFACE NON SUR LIGNES DE L'OBJET
C        L1ARIN = MCN( MNAROB + W1ARIN )

C        ADRESSE MCN DU 1-ER MOT DU TABLEAU LARETE
         MNLARE = MNAROB + W1LGFR + NUMXLF - NUMILF + 1
C        LARETE : TABLEAU DES ARETES DU MAILLAGE
C        LARETE(1,I)= NO DU 1-ER  SOMMET DE L'ARETE (0 SI PAS D'ARETE)
C        LARETE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C        LARETE(3,I)= CHAINAGE HACHAGE SUR ARETE SUIVANTE
C        LARETE(4,I)= NUMERO DU 1-ER TRIANGLE CONTENANT CETTE ARETE
C                   + NO TYPE EF (1 a NBTYEF) * 100 000 000
C                     0 SI PAS DE 1-ER  TRIANGLE
C        LARETE(5,I)= NUMERO DU 2-EME TRIANGLE CONTENANT CETTE ARETE
C                   + NO TYPE EF (1 a NBTYEF) * 100 000 000
C                     0 SI PAS DE 2-EME TRIANGLE
C        LARETE(6,I)= NUMERO DANS LARETE DE L'ARETE SUIVANTE
C                 SOIT DANS LE CHAINAGE D'UNE LIGNE J ENTRE NUMILF ET NUMXLF
C                 SOIT DANS LE CHAINAGE DES ARETES FRONTALIERES
C                 0 SI C'EST LA DERNIERE

      ELSE

C        TETRAEDRES: RECUPERATION DU TMS DES FACES DE L'OBJET 3D
C        CONSTRUCTION IMPOSEE DU TMS FACE
C        POUR REMONTER LES CARACTERISTIQUES
C        -------------------------------------------------------
         CALL LXTSOU( NTLXOB, 'FACE', NTFAOB, MNFAOB )
         IF( NTFAOB .GT. 0 ) CALL LXTSDS( NTLXOB, 'FACE' )

C        LE TABLEAU N'EXISTE PAS => IL EST CREE
C        CALCUL PAR HACHAGE DES FACES DE L'OBJET A PARTIR DE TOPO+NPEF"...
C        NECESSAIRE POUR CONNAITRE LES EF ADJACENTS PAR FACE
         CALL HACHOB( KNOMOB, 4, NTFAOB, MNFAOB, IERR )
         IF( IERR .NE. 0 .OR. NTFAOB .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'OBJET SANS TABLEAU DES FACES'
            ELSE
               KERR(1) = 'OBJECT WITHOUT THE ARRAY OF FACES'
            ENDIF
            CALL LEREUR
            GOTO 9999
         ENDIF
C        LE NOMBRE D'ENTIERS PAR FACE
         MOFACE = MCN( MNFAOB + WOFACE )
C        LA MAJORATION DU NOMBRE DE FACES DU TABLEAU LFACES
         MXFACE = MCN( MNFAOB + WXFACE )
C        LE NUMERO DE LA PREMIERE FACE FRONTALIERE (DANS UN SEUL EF)
         L1FAFR = MCN( MNFAOB + W1FAFR )
C        LE NOMBRE DE FACES FRONTALIERES
         NBFAFR = MCN( MNFAOB + WBFAFR )

C        ADRESSE MCN DU 1-ER MOT DU TABLEAU LFACES
         MNLFAC= MNFAOB + WFACES

C        POUR LE TRACE DES FLECHES
C        CREATION DU HACHAGE DES ARETES DES FACES FRONTALIERES DE L'OBJET
         CALL HACHAF( KNOMOB, 0, NTFAOB, MNFAOB,
     %                NTAFOB, MNAFOB, I )

C        LE NOMBRE D'ENTIERS PAR ARETE FRONTALIERE
         MOARFR = MCN( MNAFOB + WOARFR )
C        LA MAJORATION DU NOMBRE DES ARETES FRONTALIERES
         MXARFR = MCN( MNAFOB + WXARFR )
C        LE NUMERO DANS LAREFR DE LA PREMIERE ARETE FRONTALIERE
         L1ARFR = MCN( MNAFOB + W1ARFR )

C        ATTENTION: MODIFICATION TEMPORAIRE DU TABLEAU LFACES
C                   NECESSAIRE POUR IDENTIFIER VITE UNE FACE FRONTALIERE
C        LES FACES FRONTALIERES RECUPERENT ZERO COMME SECOND EF
C        AU LIEU DU CHAINAGE DES FACES FRONTALIERES DANS LE TABLEAU LFACES
C        LFACES(1,I)= NO DU 1-ER  SOMMET DE LA FACE
C        LFACES(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C        LFACES(3,I)= NO DU 3-EME SOMMET DE LA FACE
C        LFACES(4,I)= NO DU 4-EME SOMMET DE LA FACE
C                     0 SI TRIANGLE
C        LFACES(5,I)= CHAINAGE HACHAGE SUR FACE SUIVANTE
C        LFACES(6,I)= NUMERO DU 1-ER  CUBE CONTENANT CETTE FACE
C                     0 SI PAS DE 1-ER  CUBE
C        LFACES(7,I)= NUMERO DU 2-EME CUBE CONTENANT CETTE FACE
C                     ou CHAINAGE SUR LA FACE FRONTALIERE SUIVANTE
C        Cette valeur est remise a zero dans ce sous programme
C        puis remise au lien sur la face frontaliere suivante
C        LFACES(8,I)= NUMERO DE FACE A TANGENTES
         CALL MAZFAF( MOFACE, MXFACE, MCN(MNLFAC), L1FAFR )

C        DEFINITION DE LA VISEE AXONOMETRIQUE
C        ------------------------------------
C        POINT VU LE CENTRE DE L'HEXAEDRE ENGLOBANT
         AXOPTV(1) = ( COOEXT(1,1) + COOEXT(1,2) ) * 0.5
         AXOPTV(2) = ( COOEXT(2,1) + COOEXT(2,2) ) * 0.5
         AXOPTV(3) = ( COOEXT(3,1) + COOEXT(3,2) ) * 0.5
C        POSITION DE L'OEIL
         AXOEIL(1) = AXOPTV(1)
         AXOEIL(2) = AXOPTV(2) - 0.52 * ( COOEXT(2,2) - COOEXT(2,1) )
         AXOEIL(3) = AXOPTV(3)
C        LARGEUR/2 et HAUTEUR/2 de la FENETRE
         AXOLAR = ( COOEXT(1,2) - COOEXT(1,1) ) * 0.5
         AXOHAU = ( COOEXT(3,2) - COOEXT(3,1) ) * 0.5
C        PAS DE PLAN ARRIERE ET AVANT
         AXOARR = 0
         AXOAVA = 0
         IF( INTERA .GE. 1 )
     %   CALL AXONOMETRIE(AXOPTV, AXOEIL, AXOLAR, AXOHAU, AXOARR,AXOAVA)

      ENDIF

C     A PRIORI: CONSTRUCTION DES TABLEAUX
C     en 3D: NUVVEF(NBELEM), NUSFEF(NBSFEF,NBELEM), 
C            NULAEF(NBLAEF,NBELEM), NUPSEF(NBPSEF,NBELEM)
C     en 2D: NUSFEF(NBELEM), NULAEF(NBLAEF,NBELEM), NUPSEF(NBPSEF,NBELEM)
C     -------------------------------------------------------------------
C     MAIS COMME LE TABLEAU NUPSEF=MCN(MNPSEF) N'EST PAS ENSUITE UTILISE
C     LE TABLEAU MCN(MNPSEF) N'EST PAS CONSTRUIT
      NBPSEF = 0

      IF( NDIM .EQ. 3 ) THEN
C        EN NDIM=3 LE TABLEAU NULAEF=MCN(MNLAEF) N'EST PAS ENSUITE UTILISE
C        LE TABLEAU MCN(MNPSEF) N'EST PAS CONSTRUIT
         NBLAEF = 0
      ENDIF

      CALL EFVFAS( MNELE, NBVVEF, NBSFEF, NBLAEF, NBPSEF,
     %                    MNVVEF, MNSFEF, MNLAEF, MNPSEF, IERR )
      IF( IERR .NE. 0 ) GOTO 9993

C     ADRESSE A 1 POUR EVITER LES PROBLEMES D'APPEL de MCN(0)
      IF( NBVVEF .EQ. 0 ) MNVVEF = 1
      IF( NBSFEF .EQ. 0 ) MNSFEF = 1
      IF( NBLAEF .EQ. 0 ) MNLAEF = 1
      IF( NBPSEF .EQ. 0 ) MNPSEF = 1


C     RECUPERATION DES DONNEES PHYSIQUES et THERMIQUES A PARTIR DE L'EF 1
C     SUPPOSEES CONSTANTES DANS LE FLUIDE SUR LE PREMIER ELEMENT FINI
C     DU MAILLAGE: Rho, Mhu, CoGrPr, CoBOUS
C     ===================================================================
C     NO DES NOEUDS DE L'ELEMENT FINI 1
      NUELEM = 1
      CALL EFNOEU( MNELE, NUELEM, NBNOEF, NONOEF )

C     NO DES POINTS LIGNES SURFACES VOLUMES DES SOMMETS ARETES FACES VOLUME
      CALL EFPLSV( MNELE , NUELEM,
     %             NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %             NOOBVC, NOOBSF, NOOBLA, NOOBPS, IERR )

C     INTERPOLATION DU FLUIDE avec F: e chapeau->e  P1 ndim
C     COORDONNEES XYZEF(NBNOEF,3) DES SOMMETS=POINTS DE L'EF
      CALL EFXYZP( NDIM, MNXYZN, NBELEM, NUELEM, MNPGEL, NBNOEF,
     %             XYZEF )


C     RECHERCHE PARMI LES DONNEES du FLUIDE de
C     la DENSITE DE MASSE Rho, la VISCOSITE DYNAMIQUE Mhu,
C     le COEFFICIENT DE LA PRESSION CoGrPr AU BARYCENTRE DE L'EF 1
C     CONSIDERES comme CONSTANT pour tous les EF
C     ------------------------------------------------------------
      IF( NDIM .EQ. 2 ) THEN
         CALL REDOFL( NBNOEF, NDIM, XYZEF,
     %                NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)),
     %                Rho, Mhu, CoGrPr, CoBOUSfl )
      ELSE
         CALL REDOFL( NBNOEF, NDIM, XYZEF,
     %                NOOBVC, NUMIOB(4), NUMAOB(4), MCN(MNDOEL(4)),
     %                Rho, Mhu, CoGrPr, CoBOUSfl )
      ENDIF

      IF( LANGAG .EQ. 0 ) THEN
         PRINT 10010, Rho, Mhu, CoGrPr, CoBOUSfl
      ELSE
         PRINT 20010, Rho, Mhu, CoGrPr, CoBOUSfl
      ENDIF
10010 FORMAT(/' MASSE VOLUMIQUE Rho et VISCOSITE DYNAMIQUE Mhu NE DEPEND
     %ENT NI DU TEMPS NI DE LA VITESSE'/
     %' Rho=',g13.6,'  Mhu=',g13.6,'  CoGrPression=',g13.6,
     %' Coef Fluide BOUSSINESQ=',g13.6)
20010 FORMAT(/' Rho DENSITY of MASS and Mhu DYNAMIC VISCOSITY ARE INDEPE
     %NDENT of TIME and VELOCITY'/
     %' Rho=',g13.6,'  Mhu=',g13.6,'  CoGrPressure=',g13.6,
     %' Coef Fluide BOUSSINESQ=',g13.6)


C     RECHERCHE PARMI LES DONNEES THERMIQUES du FLUIDE de
C     LA DENSITE DE MASSE Rho, (La MEME QUE CI-DESSUS)
C     LA CAPACITE THERMIQUE MASSIQUE Cp,
C     LA CONDUCTIVITE THERMIQUE Conduc,
C     LE COEFFICIENT DE PROPORTIONNALITE DE LA DIFFERENCE DE
C     TEMPERATURE A LA VARIATION DE LA MASSE VOLUMIQUE CoBOUS
C     de l'APPROXIMATION de BOUSSINESQ Rho0-Rho=Rho CoBOUS (Temp-Temp0)
C     -----------------------------------------------------------------
      IF( NDIM .EQ. 2 ) THEN
         CALL REDOTH( NBNOEF, NDIM, XYZEF,
     %                NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDTEL(3)),
     %                Rho, Cp, Conduc, CoBOUSth )
      ELSE
         CALL REDOTH( NBNOEF, NDIM, XYZEF,
     %                NOOBVC, NUMIOB(4), NUMAOB(4), MCN(MNDTEL(4)),
     %                Rho, Cp, Conduc, CoBOUSth )
      ENDIF

      IF( LANGAG .EQ. 0 ) THEN
         PRINT 10015, Rho, Cp, Conduc, CoBOUSth
      ELSE
         PRINT 20015, Rho, Cp, Conduc, CoBOUSth
      ENDIF
10015 FORMAT(/' MASSE VOLUMIQUE Rho, CAPACITE THERMIQUE et CONDUCTIVITE 
     %NE DEPENDENT NI DU TEMPS NI DE LA VITESSE'/
     %' Rho=',g13.6,'  Cp=',g13.6,'  Conduc=',g13.6)
20015 FORMAT(/' Rho DENSITY of MASS, HEAT CAPACITY and CONDUCTIVITY ARE
     %INDEPENDENT of TIME and VELOCITY'/
     %' Rho=',g13.6,'  Cp=',g13.6,'  Conduc=',g13.6,
     %' Coef thermal BOUSSINESQ=',g13.6)

C     LE COEFFICIENT FINAL de BOUSSINESQ EST CELUI donne en FLUIDE ou THERMIQUE
      IF( CoBOUSfl .LT. 0D0 .OR. CoBOUSth .LT. 0D0 ) THEN
C        Si NEGATIF, C'est le MIN
         CoBOUS = MIN( CoBOUSfl, CoBOUSth )
      ELSE
C        Si POSITI C'est le MAX
         CoBOUS = MAX( CoBOUSfl, CoBOUSth ) 
      ENDIF

      PRINT 10016, CoBOUS
10016 FORMAT(' FINAL Coefficient BOUSSINESQ=',g13.6)


C     CONSTRUCTION DU TABLEAU NONOSO: NO NOEUD P2 => NO SOMMET P1
C     -----------------------------------------------------------
      ALLOCATE ( NONOSO(1:NBNOVI), STAT=IERNONOSO )
      IF( IERNONOSO .NE. 0 ) GOTO 9993
      CALL AZEROI( NBNOVI,  NONOSO )

C     NOMBRE TOTAL NBSOM DES SOMMETS DES EF TAYLOR-HOOD
      NDIM1 = NDIM + 1
      DO NUELEM = 1, NBELEM
         CALL EFNOEU( MNELE, NUELEM, NBNOEF, NONOEF )
         DO I = 1, NDIM1
            NONOSO( NONOEF(I) ) = 1
         ENDDO
      ENDDO

      NBSOM = 0
      DO I = 1, NBNOVI
         IF( NONOSO( I ) .GT. 0 ) THEN
            NBSOM = NBSOM + 1
            NONOSO( I ) = NBSOM
         ENDIF
      ENDDO

C     CONSTRUCTION DU TABLEAU NO1EFN(NBNOVI): No NOEUD P2 => No 1EF LE CONTENANT
C     --------------------------------------------------------------------------
      ALLOCATE ( NO1EFN(1:NBNOVI), STAT=IERNO1EFN )
      IF( IERNO1EFN .NE. 0 ) GOTO 9993
      CALL CO1EFN( NBNOEF, NBELEM, MCN(MNELE+WUNDEL), NBNOVI,  NO1EFN )

C     TEMPS CPU DE FABRICATION DES VECTEURS GLOBAUX
      DFABG = DINFO( 'DELTA CPU' )

C     ================================================================
C     CONSTRUCTION DE LA MATRICE GLOBALE DE PRESSION [PG(CoGrPr)]
C     DeltaT CoGrPr Integrale grad P1 grad P1 dX avec INTERPOLATION P1
C     DeltaT CoGrPr dPression/dn = 0 sur GammaP
C     ================================================================
C     MATRICE PRESSION SYMETRIQUE NON DIAGONALE
      NCODSP = 1

C     CONSTRUCTION DU POINTEUR SUR LA DIAGONALE DE PG (INTERPOLATION P1)
      ALLOCATE( LP1LIGN(0:NBSOM), STAT=IERLP1LIGN )
      IF( IERLP1LIGN .NE. 0 ) GOTO 9993

      IF( NORESO .EQ. 1 ) THEN

C        CROUT PROFIL POUR LA PRESSION
         IF( LANGAG .EQ. 0 ) THEN
            PRINT 10001
         ELSE
            PRINT 20001
         ENDIF

C        MAX POUR CALCULER LE MIN DES VOISINS D'UN SOMMET
         LP1LIGN( 0 ) = 0
         DO N = 1, NBSOM
            LP1LIGN( N ) = NBSOM
         ENDDO

         DO NUELEM = 1, NBELEM
C           LES NUMEROS DES NOEUDS DE L'EF NUELEM
            CALL EFNOEU( MNELE, NUELEM, NBNOEF, NONOEF )
C           RECHERCHE DU NUMERO MINIMAL DES NDIM1 SOMMETS DE l'EF
            NOSOEF(1) = NONOSO( NONOEF(1) )
            NOMIN     = NOSOEF(1)
            DO I = 2, NDIM1
C              NUMERO DU SOMMET I
               NOSOEF(I) = NONOSO( NONOEF(I) )
C              NUMERO MINIMAL DES SOMMETS
               NOMIN = MIN( NOMIN, NOSOEF(I) )
            ENDDO
            DO I = 1, NDIM1
C              NUMERO DU SOMMET I
               N = NOSOEF(I)
C              NUMERO MINIMUM DES VOISINS DU SOMMET N
               LP1LIGN( N ) = MIN( NOMIN, LP1LIGN( N ) )
            ENDDO
         ENDDO

C        POINTEUR SUR LE COEFFICIENT DIAGONAL DE LA MATRICE PRESSION
         LP1LIGN( 0 ) = 0
         DO N = 1, NBSOM
C           NOMBRE DE COEFFICIENTS ENTRE LE MIN ET DIAGONAL
            LP1LIGN( N ) = LP1LIGN( N-1 ) + N - LP1LIGN( N ) + 1
         ENDDO
C        NOMBRE DE COEFFICIENTS DE LA MATRICE PG PROFIL DE PRESSION
         NBCMPG = LP1LIGN( NBSOM )

      ELSE IF( NORESO .EQ. 2 ) THEN

C        GRADIENT CONJUGUE SUR MATRICE MORSE DE LA PRESSION
         IF( LANGAG .EQ. 0 ) THEN
            PRINT 10002
         ELSE
            PRINT 20002
         ENDIF

C        A PARTIR DES POINTEURS DE LA MATRICE MORSE VITESSE  P2 PAR NOEUDS
C        CALCUL   DES POINTEURS DE LA MATRICE MORSE PRESSION P1 PAR SOMMETS
         CALL PRGCLP2P1( NBSOM,   NBNOVI,  NONOSO,
     %                   NCODSP,  LP2LIGN, LP2COLO,
     %                   LP1LIGN, IERR )
         IF( IERR .NE. 0 ) GOTO 9993

C        NOMBRE DE COEFFICIENTS DE LA MATRICE MORSE PRESSION AUX SOMMETS
         NBCMPG = LP1LIGN( NBSOM )
C        LE POINTEUR SUR LES COLONNES DES SOMMETS P1
         ALLOCATE( LP1COLO(1:NBCMPG), STAT=IERLP1COLO )
         IF( IERLP1COLO .NE. 0 ) GOTO 9993

         CALL PRGCCP2P1( NBNOVI,  NONOSO,  NBSOM,
     %                   NCODSP,  LP2LIGN, LP2COLO, NBCMPG,
     %                   LP1LIGN, LP1COLO, IERR )
         IF( IERR .NE. 0 ) GOTO 9993

      ELSE

C        ERREUR SUR NORESO. METHODE DE RESOLUTION Ax=b INCONNUE
         IF( LANGAG .EQ. 0 ) THEN
            PRINT 10004, NORESO
         ELSE
            PRINT 20004, NORESO
         ENDIF
         IERR = 1
         GOTO 9999

      ENDIF

10001 FORMAT('CONSTRUCTION de la MATRICE PROFIL [PG] de PRESSION'/
     %' RESOLUTION par FACTORISATION COMPLETE PG=L D tL de CROUT')
20001 FORMAT(' CONSTRUCTION of the PRESSURE SKYLINE MATRIX [PG]'/
     %' SOLUTION by COMPLETE CROUT FACTORIZATION PG=L D Lt')

10002 FORMAT('CONSTRUCTION de la MATRICE MORSE [PG] de PRESSION'/
     %'RESOLUTION par GRADIENT CONJUGUE')
20002 FORMAT(' CONSTRUCTION of the PRESSURE CONDENSED MATRIX [PG]'/
     %'SOLUTION by CONJUGATE GRADIENT')

10004 FORMAT('CODE DE RESOLUTION INCONNU  NORESO=',I5/)
20004 FORMAT('UNKNOWN SOLUTION CODE  NORESO=',I5/)

C     MATRICE DE PRESSION SYMETRIQUE NON DIAGONALE
C     NBCMPG : NOMBRE DE REELS DOUBLE PRECISION DE LA MATRICE PG P1
C     ALLOCATION de la MATRICE PG en LANGAGE FORTRAN 90
      PRINT*
      PRINT*, 'ALLOCATION DEMAND  of',NBCMPG,
     %        ' DOUBLE PRECISION of [PG] MATRIX'
      ALLOCATE ( PG(1:NBCMPG), STAT=IERPGALLOC )
      IF( IERPGALLOC .NE. 0 ) THEN
         PRINT*,'ALLOCATION ERROR   of',NBCMPG,
     %          ' DOUBLE PRECISION of [PG] MATRIX'
         IERR = IERPGALLOC
         GOTO 9993
      ENDIF
      PRINT*, 'ALLOCATION CORRECT of',NBCMPG,
     %        ' DOUBLE PRECISION of [PG] MATRIX'


C     CONSTRUCTION DES TABLEAUX DES NO ET VALEURS DES PRESSIONS FIXEES
C     ----------------------------------------------------------------
      IF( IEBLPR .GT. 0 ) THEN
         CALL PRESFXST( RELMIN, NBSOM,   NDIM,   MNXYZN, NONOSO,
     %                  NBTYEL, MNNPEF,  NUMIOB, MNDOEL,
     %                  NBPRFX, MNNPRFX, MNVPRFX )
      ELSE
         NBPRFX = 0
      ENDIF

C     LE COEFFICIENT DEVANT LA MATRICE DE PRESSION PG
      S = DELTAT * CoGrPr

      IF( NORESO .EQ. 1 ) THEN

C        ASSEMBLAGE PROFIL: dt CoGrPr ( GRAD p, GRAD q ) INTERPOLATION P1
C        DE LA MATRICE PG GLOBALE DE LA PRESSION
C        ----------------------------------------------------------------
         IF( NUTYEL .EQ. 15 ) THEN
C           TRIANGLE 2D TAYLOR-HOOD
            CALL P12DAGPR( S,      NBSOM,  MCN(MNXYZN+WYZNOE),
     %                     NBNOEF, NBELEM, MCN(MNELE+WUNDEL),
     %                     NONOSO, NBCMPG, LP1LIGN, PG )
         ELSE
C           TETRAEDRE TAYLOR-HOOD
            CALL P13DAGPR( S,      NBSOM,  MCN(MNXYZN+WYZNOE),
     %                     NBNOEF, NBELEM, MCN(MNELE+WUNDEL),
     %                     NONOSO, NBCMPG, LP1LIGN, PG )
         ENDIF

         call affvect( 'PG PROFIL ASSEMBLE avant CL', 5,      PG )
         call afl1ve(  'PG PROFIL ASSEMBLE avant CL', NBCMPG, PG )

cccC        SUPPRESSION DES COEFFICIENTS D'ERREURS D'ARRONDIS
cccC        et MULTIPLICATION PAR  dt CoGrPr
ccc         CALL PRARMU( 1D-4, S, NBCMPG, PG )
ccc         call affvect( 'PG PROFIL ASSEMBLE corrige',5,  PG )

C        PRISE EN COMPTE SUR PG DES PRESSIONS FIXEES
         DO I = 1, NBPRFX
C           LE NO GLOBAL DU DL FIXE
            NDL = MCN( MNNPRFX - 1 + I )
C           LE NO DU COEFFICIENT DIAGONAL NDL
            NDIAG = LP1LIGN( NDL )
C           MODIFICATION DU COEFFICIENT DIAGONAL
C           CETTE VALEUR PEUT SERVIR DE MARQUEUR D'UN DL FIXE POUR PG
            PG( NDIAG ) = STGV
         ENDDO
         PRINT *,' PG avec STGV=',STGV,' Nb PRESSIONS FIXEES=',NBPRFX
         call affvect( 'PG PROFIL apres CONDITIONS aux LIMITES', 5,
     %                  PG )
         call afl1ve(  'PG PROFIL apres CONDITIONS aux LIMITES', NBCMPG,
     %                  PG )

C        FACTORISATION COMPLETE DE CROUT  PG = L D tL
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*, 'DEBUT de FACTORISATION CROUT PG=L D tL Nb lignes=',
     %               NBSOM
         ELSE
            PRINT*, 'STARTING of CROUT FACTORIZATION PG=L D tL Lines Num
     %ber=',NBSOM
         ENDIF
         CALL CRMC1D( LP1LIGN, PG, NBSOM, EPSCROUT, NENTRE,  PG, IERR )
         call affvect( 'PG=L D tL PROFIL FACTORISEE', 5, PG )
         call afl1ve(  'PG=L D tL PROFIL FACTORISEE', NBCMPG, PG )

      ELSE IF( NORESO .EQ. 2 ) THEN

C        ASSEMBLAGE MORSE dt CoGrPr ( GRAD p, GRAD q ) INTERPOLATION P1
C        DE LA MATRICE GLOBALE PG DE LA PRESSION
C        RESOLUTION PAR GRADIENT CONJUGUE PRECONDITIONNE
C        --------------------------------------------------------------
         IF( NUTYEL .EQ. 15 ) THEN
C           TRIANGLE 2D  PRESSION P1
            CALL P12DAGGC( S,      NBSOM,   MCN(MNXYZN+WYZNOE),
     %                     NBNOEF, NBELEM,  MCN(MNELE+WUNDEL), NONOSO,
     %                     NBCMPG, LP1LIGN, LP1COLO, PG )
         ELSE
C           TETRAEDRE 3D  PRESSION P1
            CALL P13DAGGC( S,      NBSOM,   MCN(MNXYZN+WYZNOE),
     %                     NBNOEF, NBELEM,  MCN(MNELE+WUNDEL), NONOSO,
     %                     NBCMPG, LP1LIGN, LP1COLO, PG )
         ENDIF
         call affvect( 'PG MORSE avant CL', 5,      PG )
         call afl1ve(  'PG MORSE avant CL', NBCMPG, PG )

C        CONSTRUCTION DU TABLEAU des DL PRESSION FIXES 1 OU NON 0
         ALLOCATE ( NODLPRFX(1:NBSOM), STAT=IERNODLPRFX )
         IF( IERNODLPRFX .NE. 0 ) GOTO 9993

C        A PRIORI TOUS LES DL SONT LIBRES => 0
         CALL AZEROI( NBSOM, NODLPRFX )

C        MISE AU NUMERO DU DL FIXE DANS LES TABLEAUX
C        DES DL FIXES DE LA PRESSION
         DO I = 1, NBPRFX
C           LE NO GLOBAL DU DL PRESSION FIXEE
            NDL = MCN( MNNPRFX - 1 + I )
C           MISE A I DU DL PRESSION FIXEE NDL
            NODLPRFX( NDL ) = I
C           MODIFICATION DU COEFFICIENT DIAGONAL
C           LE NO DU COEFFICIENT DIAGONAL NDL
            NDIAG = LP1LIGN( NDL )
C           CETTE VALEUR PEUT SERVIR DE MARQUEUR D'UN DL FIXE POUR PG
            PG( NDIAG ) = 1D0
         ENDDO
         call affvect( 'PG MORSE apres CL', 5,      PG )
         call afl1ve(  'PG MORSE apres CL', NBCMPG, PG )

C        LES 3 TABLEAUX AUXILIAIRES du SOUS PROGRAMME GCAxbk
         ALLOCATE ( DAUX1(1:NBNOVI), STAT=IERDAUX1 )
         IF( IERDAUX1 .NE. 0 ) THEN
            PRINT*,'ALLOCATION ERROR   of',NBNOVI,
     %             ' DOUBLE PRECISION REALS of DAUX1'
            IERR = IERDAUX1
            GOTO 9993
         ENDIF

         ALLOCATE ( DAUX2(1:NBNOVI), STAT=IERDAUX2 )
         IF( IERDAUX2 .NE. 0 ) THEN
            PRINT*,'ALLOCATION ERROR   of',NBNOVI,
     %             ' DOUBLE PRECISION REALS of DAUX2'
            IERR = IERDAUX2
            GOTO 9993
         ENDIF

         ALLOCATE ( DAUX3(1:NBNOVI), STAT=IERDAUX3 )
         IF( IERDAUX3 .NE. 0 ) THEN
            PRINT*,'ALLOCATION ERROR   of',NBNOVI,
     %             ' DOUBLE PRECISION REALS of DAUX3'
            IERR = IERDAUX3
            GOTO 9993
         ENDIF

         ALLOCATE ( PVAUX(1:NBNOVI), STAT=IERPVAUX )
         IF( IERPVAUX .NE. 0 ) THEN
            PRINT*,'ALLOCATION ERROR   of',NBNOVI,
     %             ' DOUBLE PRECISION REALS of PVAUX'
            IERR = IERPVAUX
            GOTO 9993
         ENDIF

      ENDIF

C     LE TABLEAU AUXILIAIRE de P(tn+1)-P(tn)
      ALLOCATE ( Pr1Pr(1:NBSOM), STAT=IERPr1Pr )
      IF( IERPr1Pr .NE. 0 ) THEN
         PRINT*,'ALLOCATION ERROR   of',NBSOM,
     %          ' DOUBLE PRECISION REALS of Pr1Pr'
         IERR = IERPr1Pr
         GOTO 9993
      ELSE
         CALL AZEROD( NBSOM, Pr1Pr )
      ENDIF

C     TEMPS CALCUL FORMATION DE LA MATRICE GLOBALE PG ET FACTORISATION
      DFACTO = DINFO( 'DELTA CPU' )
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'TEMPS FORMATION FACTORISATION MATRICE PG=',DFACTO
      ELSE
        PRINT*,'PG MATRIX FORMATION FACTORIZATION TIME=',DFACTO
      ENDIF

C     ===========================================================
C     DECLARATION DE LA MATRICE VG D'UNE COMPOSANTE DE LA VITESSE
C     MATRICE SYMETRIQUE NON DIAGONALE => NCODSV = 1
C     Integrale Rho P2 P2 + deltat Mhu grad P2 grad P2 dX
C     avec INTEGRATION EXACTE SUR L'ELEMENT FINI
C     ===========================================================
      NCODSV = 1
C     NBCMVG : NOMBRE DE REELS DOUBLE PRECISION DE LA MATRICE VG
C     VERSION ALLOCATION de la MATRICE VG en LANGAGE FORTRAN 90
      PRINT*
      PRINT*, 'ALLOCATION DEMAND  of',NBCMVG,
     %        ' DOUBLE PRECISION of [VG] MATRIX'
      ALLOCATE ( VG(1:NBCMVG), STAT=IERVGALLOC )
      IF( IERVGALLOC .NE. 0 ) THEN
         PRINT*,'ALLOCATION ERROR   of',NBCMVG,
     %          ' DOUBLE PRECISION of [VG] MATRIX'
         IERR = IERVGALLOC
         GOTO 9993
      ENDIF
      PRINT*, 'ALLOCATION CORRECT of',NBCMVG,
     %        ' DOUBLE PRECISION of [VG] MATRIX'

C     CONSTRUCTION ET ASSEMBLAGE DE LA MATRICE GLOBALE VG
C     ---------------------------------------------------
      CALL AGVITTH( Rho,    DELTAT*Mhu,  MCN(MNXYZN+WYZNOE),
     %              NUTYEL, NBNOEF, NBELEM, MCN(MNELE+WUNDEL),
     %              NORESO, NBNOVI, NBCMVG, LP2LIGN, LP2COLO, VG )
      call affvect( 'CONSTRUCTION VG avant CL', 5,      VG )
      call afl1ve(  'CONSTRUCTION VG avant CL', NBCMVG, VG )

C     LISTE DES DL VITESSES FIXEES EN TOUT TEMPS
C     ATTENTION: CETTE PROGRAMMATION SUPPOSE QUE TOUS LES DL EN UN NOEUD
C                SONT FIXES ET LES NOEUDS FIXES RESTENT LES MEMES POUR TOUT TEMPS
C                MAIS PAS FORCEMENT LA VITESSE A FIXER QUI ELLE DEPEND   du TEMPS
C    ----------------------------------------------------------------------------
C     LISTE DES DL VITESSES FIXEES EN TOUT TEMPS
      CALL VITEFX( RELMIN, NTDLVI, NDIM,   MNXYZN,
     %             NBTYEL, MNNPEF, NUMIOB, MNDOEL,
     %             NBVCFX, NBVIFX, MNNVIFX, MNVVIFX, VITFXMAX )
      VITMAX0 = MAX( VITMAX0, VITFXMAX )

      IF( NORESO .EQ. 2 ) THEN
C        GRADIENT CONJUGUE SIMPLE SUR VG MATRICE MORSE
C        CONSTRUCTION DU TABLEAU des DL 1VITESSE FIXE I OU NON 0
         ALLOCATE ( NODLVIFX(1:NBNOVI), STAT=IERNODLVIFX )
         IF( IERNODLVIFX .NE. 0 ) GOTO 9993

C        A PRIORI TOUS LES DL SONT LIBRES => 0
         CALL AZEROI( NBNOVI, NODLVIFX )

C        MISE AU NUMERO DU DL FIXE DANS LES TABLEAUX VIFX
C        DES DL FIXES DE LA 1-ERE COMPOSANTE DE LA VITESSE
         DO I=1, NBVIFX/NDIM
C           LE NO GLOBAL DU DL 1VITESSE FIXEE
            NDL = MCN( MNNVIFX - 1 + I )
C           MISE A I DU DL 1 VITESSE FIXEE NDL
            NODLVIFX( NDL ) = I
         ENDDO
      ENDIF

C     *******************************************************************
C     RESTRICTION: TOUTES LES NDIM COMPOSANTES DE LA VITESSE EN UN NOEUD
C                  DOIVENT ETRE FIXEES DE VALEURS POUVANT VARIER EN TEMPS
C     *******************************************************************
C     PRISE EN COMPTE DE LA VITESSE FIXEE SUR LA MATRICE VG
      DO I = 1, NBVIFX/NDIM
C        LE NO GLOBAL DU DL FIXE DE LA PREMIERE COMPOSANTE DE LA VITESSE
         NDL = MCN( MNNVIFX - 1 + I )
C        LE NO DU COEFFICIENT DIAGONAL NDL
         NDIAG = LP2LIGN( NDL )
C        MODIFICATION DU COEFFICIENT DIAGONAL
C        CETTE VALEUR PEUT SERVIR DE MARQUEUR D'UN DL FIXE POUR VG
         VG( NDIAG ) = STGV
      ENDDO
      PRINT *,' VG avec STGV=',STGV,' Nb VITESSES FIXEES=',NBVIFX/NDIM
ccc      call affvect( 'VG apres CONDITIONS aux LIMITES', 5,      VG )
ccc      call afl1ve(  'VG apres CONDITIONS aux LIMITES', NBCMVG, VG )

      IF( NORESO .EQ. 1 ) THEN

C        FACTORISATION COMPLETE DE CROUT  VG = L D tL
C        --------------------------------------------
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'DEBUT DE FACTORISATION CROUT VG=L D tL  Nb lignes=',
     %              NBNOVI
         ELSE
            PRINT*,'STARTING of CROUT FACTORIZATION VG=L D tL Lines Numb
     %er=',NBNOVI
         ENDIF
         CALL CRMC1D( LP2LIGN, VG, NBNOVI, EPSCROUT, NENTRE,  VG, IERR )
ccc         call affvect( 'VG=L D tL PROFIL FACTORISEE', 5,      VG )
ccc         call afl1ve(  'VG=L D tL PROFIL FACTORISEE', NBCMVG, VG )

      ENDIF

C     TEMPS CALCUL DE FORMATION DES MATRICES GLOBALES PG VG ET FACTORISATIONS
      DFACTO = DFACTO + DINFO( 'DELTA CPU' )
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'TEMPS FORMATION FACTORISATION MATRICE PG + VG=',DFACTO
      ELSE
        PRINT*,'PG + VG MATRIX FORMATION FACTORIZATION TIME=',DFACTO
      ENDIF

C     ===========================================================
C     DECLARATION DE LA MATRICE TG GLOBALE DE LA TEMPERATURE
C     MATRICE TG SYMETRIQUE NON DIAGONALE => NCODST = 1
C     Integrale Rho Cp P2 P2 + deltat Conduc grad P2 grad P2 dX
C     avec INTEGRATION EXACTE SUR L'ELEMENT FINI
C     ===========================================================
      NCODST = 1
C     NBCMTG : NOMBRE DE REELS DOUBLE PRECISION DE LA MATRICE TG
C              IDENTIQUE A CELUI DE LA MATRICE VG
      NBCMTG = NBCMVG
C     VERSION ALLOCATION de la MATRICE TG en LANGAGE FORTRAN 90
      PRINT*
      PRINT*, 'ALLOCATION DEMAND  of',NBCMTG,
     %        ' DOUBLE PRECISION of [TG] MATRIX'
      ALLOCATE ( TG(1:NBCMTG), STAT=IERTGALLOC )
      IF( IERTGALLOC .NE. 0 ) THEN
         PRINT*,'ALLOCATION ERROR   of',NBCMTG,
     %          ' DOUBLE PRECISION of [TG] MATRIX'
         IERR = IERTGALLOC
         GOTO 9993
      ENDIF
      PRINT*, 'ALLOCATION CORRECT of',NBCMTG,
     %        ' DOUBLE PRECISION of [TG] MATRIX'

C     CONSTRUCTION ET ASSEMBLAGE DE LA MATRICE GLOBALE TG
C     ---------------------------------------------------
      CALL AGVITTH( Rho*Cp, DELTAT*Conduc,  MCN(MNXYZN+WYZNOE),
     %              NUTYEL, NBNOEF, NBELEM, MCN(MNELE+WUNDEL),
     %              NORESO, NBNOVI, NBCMTG, LP2LIGN, LP2COLO, TG )
      call affvect( 'CONSTRUCTION TG avant CL', 5,      TG )
      call afl1ve(  'CONSTRUCTION TG avant CL', NBCMTG, TG )

C     LISTE DES NUMEROS ET VALEURS DES TEMPERATURES FIXEES AU TEMPS
C     -------------------------------------------------------------
      CALL THDLFX(      1, NTDLTE,  NDIM,
     %             NBTYEL, MNNPEF,  NDPGST,
     %             MNXYZN, NUMIOB,  MNDTEL,  RELMIN,
     %             NBTEFX, MONTEFX, MNNTEFX, MNVTEFX, IERR )

C     LISTE EVENTUELLE DES DL TEMPERATURES FIXEES POUR TOUS LES TEMPS
C    ----------------------------------------------------------------
      IF( NORESO .EQ. 2 ) THEN
C        GRADIENT CONJUGUE SIMPLE SUR TG MATRICE MORSE
C        CONSTRUCTION DU TABLEAU des DL 1TEMPERATURE FIXE I OU NON 0
         ALLOCATE ( NODLTEFX(1:NTDLTE), STAT=IERNODLTEFX )
         IF( IERNODLTEFX .NE. 0 ) GOTO 9993

C        A PRIORI TOUS LES DL SONT LIBRES => 0
         CALL AZEROI( NTDLTE, NODLTEFX )

C        MISE AU NUMERO DU DL FIXE DANS LES TABLEAUX TEFX
         DO I=1, NBTEFX
C           LE NO GLOBAL DU DL TEMPERATURE FIXEE
            NDL = MCN( MNNTEFX - 1 + I )
C           MISE A I DU DL TEMPERATURE FIXEE NDL
            NODLTEFX( NDL ) = I
         ENDDO
      ENDIF

C     PRISE EN COMPTE DES TEMPERATURES FIXEES SUR LA MATRICE TG
C     ---------------------------------------------------------
      DO I = 1, NBTEFX
C        LE NO GLOBAL DU DL FIXE DE LA PREMIERE COMPOSANTE DE LA VITESSE
         NDL = MCN( MNNTEFX - 1 + I )
C        LE NO DU COEFFICIENT DIAGONAL NDL
         NDIAG = LP2LIGN( NDL )
C        MODIFICATION DU COEFFICIENT DIAGONAL
C        CETTE VALEUR PEUT SERVIR DE MARQUEUR D'UN DL FIXE POUR TG
         TG( NDIAG ) = STGV
      ENDDO
      PRINT *,' TG avec STGV=',STGV,' NB TEMPERATURES FIXEES=',NBTEFX
      call affvect( 'TG apres CONDITIONS aux LIMITES', 5,      TG )
      call afl1ve(  'TG apres CONDITIONS aux LIMITES', NBCMTG, TG )

      IF( NORESO .EQ. 1 ) THEN

C        FACTORISATION COMPLETE DE CROUT  TG = L D tL
C        --------------------------------------------
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'DEBUT DE FACTORISATION CROUT TG=L D tL  Nb lignes=',
     %              NBNOVI
         ELSE
            PRINT*,'STARTING of CROUT FACTORIZATION TG=L D tL Lines Numb
     %er=',NBNOVI
         ENDIF
         CALL CRMC1D( LP2LIGN, TG, NBNOVI, EPSCROUT, NENTRE,  TG, IERR )
         call affvect( 'TG=L D tL PROFIL FACTORISEE', 5,      TG )
         call afl1ve(  'TG=L D tL PROFIL FACTORISEE', NBCMVG, TG )

      ENDIF

C     TEMPS CALCUL DE FORMATION DES 3 MATRICES GLOBALES ET FACTORISATIONS
      DFACTO = DFACTO + DINFO( 'DELTA CPU' )
      IF( LANGAG .EQ. 0 ) THEN
        PRINT*,'TEMPS FORMATION FACTORISATION MATRICES PG VG TG=',DFACTO
      ELSE
        PRINT*,'PG VG TG MATRICES FORMATION FACTORIATION TIME=',DFACTO
      ENDIF


C     DECLARATIONS DES VECTEURS GLOBAUX V3P0 et V3P1(NTDLVP)
C     ------------------------------------------------------
C     DECLARATION DU VECTEUR V3P0 SOLUTION AU TEMPS 0 EN 4 VECTEURS VX VY VZ PR
C     DECLARATION DU VECTEUR V3P1 SOLUTION AU TEMPS 1 EN 4 VECTEURS VX VY VZ PR
      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'nst3th: DEMANDE  ALLOCATION VXVYVZPR(',NTDLVP,
     %          ', 0:1) DREELS'
         ALLOCATE ( VXVYVZPR( 1:NTDLVP, 0:1 ), STAT=IERVXVYVZPR )
         IF( IERVXVYVZPR .NE. 0 ) THEN
            PRINT*,'nst3th: ERREUR ALLOCATION VXVYVZPR(',NTDLVP,
     %             ', 0:1) DREELS'
            IERR = IERVXVYVZPR
            GOTO 9993
         ENDIF
         PRINT*,'nst3th: CORRECTE ALLOCATION VXVYVZPR(',NTDLVP,
     %          ', 0:1) DREELS'
      ELSE
         PRINT*,'nst3th: ALLOCATION DEMAND  VXVYVZPR(',NTDLVP,
     %          ', 0:1) DREALS'
         ALLOCATE ( VXVYVZPR( 1:NTDLVP, 0:1 ), STAT=IERVXVYVZPR )
         IF( IERVXVYVZPR .NE. 0 ) THEN
            PRINT*,'nst3th: ALLOCATION ERROR VXVYVZPR(',NTDLVP,
     %             ', 0:1) DREALS'
            IERR = IERVXVYVZPR
            GOTO 9993
         ENDIF
         PRINT*,'nst3th: CORRECT ALLOCATION VXVYVZPR(',NTDLVP,
     %          ', 0:1) DREALS'
      ENDIF

C     NUMERO DES VECTEURS VITESSE+PRESSION DANS VXVYVZPR
C     POUR EVITER DES PERMUTATIONS COUTEUSES...
      NVP0 = 0
      NVP1 = 1

C     NO DU PREMIER VITX VITY VITZ PRESSION DANS VXVYVZPR
      N1X = 1
      N1Y = N1X + NBNOVI
      N1Z = N1Y + NBNOVI
      N1P = N1X + NBNOVI * NDIM
      N0P = N1P - 1

      CALL VXVYVZPR1( NDIM,   NUTYEL, NDDLNO, 
     %                NTDLVP, 1,      VXYZPN, NONOSO,
     %                NBNOVI, VXVYVZPR(N1X,NVP0), VXVYVZPR(N1Y,NVP0),
     %                        VXVYVZPR(N1Z,NVP0),
     %                NBSOM,  VXVYVZPR(N1P,NVP0) )

C     TABLEAU AUXILIAIRE V3P2(NTDLVI)=V3P2(NBNOVI,NDIM)
      ALLOCATE ( V3P2( 1:NTDLVI ), STAT=IERV3P2 )
      IF( IERV3P2 .NE. 0 ) THEN
         PRINT*,'nst3th: ERREUR ALLOCATION V3P2(',NTDLVI,') DREELS'
         IERR = IERV3P2
         GOTO 9993
      ENDIF

C     TABLEAU AUXILIAIRE VITC DE LA VITESSE CONVECTEE AUX NOEUDS DU MAILLAGE
      IF( NBVCFX .GT. 0 ) THEN
         MOVITC = NTDLVI
      ELSE
         MOVITC = 1
      ENDIF
      ALLOCATE ( VITC( 1:MOVITC ), STAT=IERVITC )
      IF( IERVITC .NE. 0 ) THEN
         PRINT*,'nst3th: ERREUR ALLOCATION VITC(',MOVITC,') DREELS'
         IERR = IERVITC
         GOTO 9993
      ELSE
         PRINT*,'nst3th: ALLOCATION CORRECTE VITC(',MOVITC,') DREELS'
      ENDIF
      CALL AZEROD( MOVITC, VITC )

cccC     L'ADRESSE DE LA VITESSEPRESSION A L'ITERATION 0 DU POINT FIXE
cccC     NON UTILISE ICI PUISQUE LA NON LINEARITE EST PRISE EN COMPTE
cccC     PAR L'INTEGRATION RETROGRADE EN TEMPS LE LONG DES CARACTERISTIQUES
ccc      MNVITI = 0

      IF( NDIM .EQ. 3 ) THEN
C        FLUX DE LA VITESSE SUR LES SURFACES DE L'OBJET POUR UN TEMPS DONNE
         MOSFOB = NUMAOB(3) - NUMIOB(3) + 1
      ELSE
C        FLUX DE LA VITESSE SUR LES ARETES DE L'OBJET POUR UN TEMPS DONNE
         MOSFOB = NUMAOB(2) - NUMIOB(2) + 1
      ENDIF
      CALL TNMCDC( 'REEL2', MOSFOB, MNFVSF )

C     BILAN DE L'INSTANT INITIAL
C     LE PROCHAIN TEMPS POUR STOCKER LE VECTEUR VITESSEPRESSION
      TSTOC = TPSINI + DTSTOC
C
C     AFFICHAGE DU VECTEUR"VITESSEPRESSION AU TEMPS INITIAL
C     -----------------------------------------------------
      CALL AFVIPRTE( NUTYEL, NDIM,   NBNOVI, MNXYZN, NDDLNO,
     %               TEMPS,  NTDLVP, MIN(2,NBNOVI),  VXYZPN, TEMPER0,
     %               VITMIN, NOVMIN, VITMAX0,NOVMAX, VITMOY0,
     %               PREMIN, NOPMIN, PREMAX, NOPMAX, PREMOY,
     %               TEMMIN, NOTMIN, TEMMAX, NOTMAX, TEMMOY )
      TEMMINt0 = TEMMIN
      TEMMAXt0 = TEMMAX

cccC     VITESSE D'ECRETAGE INITIALE
ccc      IF( NBPAST .GT. 0 ) THEN
cccC        LA VITESSE RECUPEREE A DEJA SUBI DES ECRETAGES
ccc         VITECR = VITMAX0
ccc      ELSE
ccc         VITECR = 10 * VITMOY0 + VITMAX0
ccc      ENDIF

      IF( NDIM .EQ. 3 ) THEN

C        FLUX DE LA VITESSE A TRAVERS LES FACES DES SURFACES DE L'OBJET
C        --------------------------------------------------------------
         CALL FLUXVITSF( NBNOVI, MCN(MNXYZN+WYZNOE),
     %                   NUTYEL, NBNOEF, NBELEM, MCN(MNELE+WUNDEL),
     %                   NDDLNO, NBSFEF, MCN(MNSFEF),
     %                   VXYZPN, NUMIOB(3), NUMAOB(3),
     %                   FLUNEG, FLUPOS, MCN(MNFVSF) )
      ELSE

C        FLUX DE LA VITESSE A TRAVERS LES ARETES DES LIGNES DE L'OBJET
C        -------------------------------------------------------------
         CALL FLUXVITLG( NBNOVI, MCN(MNXYZN+WYZNOE),
     %                   NUTYEL, NBNOEF, NBELEM, MCN(MNELE+WUNDEL),
     %                   NDDLNO, NBLAEF, MCN(MNLAEF),
     %                   VXYZPN, NUMIOB(2), NUMAOB(2),
     %                   FLUNEG, FLUPOS, MCN(MNFVSF) )
      ENDIF

C     STOCKAGE AU TEMPS INITIAL DES MIN-MAX DE VITESSE et PRESSION
      VPMIMX(1,NBPAST) = TEMPS
      VPMIMX(2,NBPAST) = VITMOY0
      VPMIMX(3,NBPAST) = VITMAX0
      VPMIMX(4,NBPAST) = PREMOY
      VPMIMX(5,NBPAST) = PREMAX - PREMIN
      VPMIMX(6,NBPAST) = FLUNEG
      VPMIMX(7,NBPAST) = FLUPOS
      VPMIMX(8,NBPAST) = NBPAST
      VPMIMX(9,NBPAST) = TEMMOY
      VPMIMX(10,NBPAST) = TEMMIN
      VPMIMX(11,NBPAST) = TEMMAX


cccC     INITIALISATION DU VECTEUR des TEMPERATURES INITIALES
cccC     A PARTIR DES DONNEES SUR L'OBJET
cccC     ====================================================
cccC     ENTREE DE LA TEMPERATURE INITIALE DE L'OBJET:
cccC     CONSTRUCTION DU TMS '~>OBJET>NomObjet>TEMPERINIT'
ccc      L1 = NUDCNB(KNOMOB)
ccc      NOMTS  = '~>OBJET>' // KNOMOB(1:L1) // '>TEMPERINIT'
ccc      KNOMTD = '~>>>TEMPERINIT'
ccc      L1 = NUDCNB(NOMTS)
ccc      CALL MOTSTD( KNOMTD, NOMTS(1:L1), IERR )
ccc      IF( IERR .NE. 0 ) THEN
ccc         NBLGRC(NRERR) = 1
ccc         IF( LANGAG .EQ. 0 ) THEN
ccc            KERR(1)='IMPOSSIBLE RECUPERER TEMPERATURE INITIALE de l''OBJET'
ccc         ELSE
ccc       KERR(1)='IMPOSSIBLE to RECOVER INITIAL TEMPERATURE of the OBJECT'
ccc         ENDIF
ccc         CALL LEREUR
ccc         IERR = 1
ccc         RETURN
ccc      ENDIF

C     ADRESSAGE DES TABLEAUX AUXILIAIRES ET ELEMENTAIRES DES EF P2
C     DEJA RETROUVES DANS fluideNS.f
C     ============================================================
      CALL TNMCDC( 'REEL2', MOAUX, MNTAUX )
C
C     MOTAEL = 2 MATRICES ELEMENTAIRES (DE CAPACITE + DE CONDUCTIVITE)
C              ET MXVECT SECOND MEMBRES ELEMENTAIRES
      MXVECT = 1
      MOTAEL = NBDL1EF * (NBDL1EF+1) + NBDL1EF * (1+MXVECT)
      CALL TNMCDC( 'REEL2', MOTAEL, MNTAEL )

C     MNTHER : CONDUCTIVITE AUX POINTS D'INTEGRATION NUMERIQUE DE L'EF
      CALL TNMCDC( 'REEL2', 128, MNTHER )
C
C     LE TABLEAU DU NUMERO DES DEGRES DE LIBERTE D'UN ELEMENT FINI P2
      MONODL = NBDL1EF * 2
      CALL TNMCDC( 'ENTIER', MONODL, MNNODL )
C
C     LE TABLEAU DES 3 COORDONNEES DES NBDL1EF POINTS DE L'EF MAXIMAL
      MNX = 0
      CALL TNMCDC( 'REEL', NBDL1EF*3, MNX )
C     LE TABLEAU DES 3 COORDONNEES DES NBDL1EF POINTS DE L'EF MAXIMAL

C     LE VECTEUR GLOBAL SECOND MEMBRE DES TEMPERATURES
      NBJEUX = 1
      IF( MNBG .LE. 0 ) THEN
         CALL TNMCDC( 'REEL2', NTDLTE, MNBG )
         MNBGD = ( MNBG+1 ) / MOREE2
         CALL AZEROD( NTDLTE, DMCN(MNBGD) )
      ENDIF

cccC     INITIALISATION DU VECTEUR DMCN(MNTHET0) DES TEMPERATURES
cccC     AUX NOEUDS DU MAILLAGE ET A L'INSTANT TEMPS INITIAL
cccC     ========================================================
cccC     MNTHET0 EST LE VECTEUR GLOBAL DES TEMPERATURES INITIALES
ccc      CALL TEMINI( KNOMOB, NTLXOB, MOREE2, NTDLTE, TPSINI, IETEIN,
ccc     %             NBTYEL, MNNPEF, NDPGST,
ccc     %             MNXYZN, NUMIOB, MNDTEL,
ccc     %             NTVECT, MNVECT, NBVTEMP, MNTHET0, IERR )
ccc      IF( IERR .NE. 0 ) GOTO 9999
ccc      MNTHET0D = ( MNTHET0 + 1 ) / MOREE2
cccC     COPIE DE LA TEMPERATURE INITIALE DANS TEMPERATURE t0
ccc      CALL TRTATD( TEMPER0, DMCN(MNTHET0D), NBNOVI )

C     LE VECTEUR GLOBAL DES TEMPERATURES tn EST DECLARE
      IF( MNTHETn .LE. 0 ) THEN
         CALL TNMCDC( 'REEL2', NTDLTE, MNTHETn )
      ENDIF
      MNTHETnD = ( MNTHETn + 1 ) / MOREE2
C     COPIE DE LA TEMPERATURE INITIALE DANS TEMPERATURE tn
      CALL TRTATD( TEMPER0, DMCN(MNTHETnD), NTDLTE )

C     LE VECTEUR GLOBAL DES TEMPERATURES tn+1 EST DECLARE
      IF( MNTHET .LE. 0 ) THEN
         CALL TNMCDC( 'REEL2', NTDLTE, MNTHET )
      ENDIF
      MNTHETD = ( MNTHET + 1 ) / MOREE2
C     COPIE DE LA TEMPERATURE INITIALE DANS TEMPERATURE tn+1
      CALL TRTATD( TEMPER0, DMCN(MNTHETD),  NTDLTE )

C     TEMPS CPU DE FABRICATION DES VECTEURS GLOBAUX
      DFABG = DFABG + DINFO( 'DELTA CPU' )

C     LA DATE EN SECONDES DEPUIS LE 1/1/70 MINUIT
      CALL SECONDES1970( DATE0 )
      TIMEMOY1DT = 0D0

C     LE NOMBRE D'ITERATIONS DU GC MIS A JOUR APRES CONVERGENCE
      KGCTEMP   = NBNOVI
      KGCUET(1) = NBNOVI
      KGCUET(2) = NBNOVI
      KGCUET(3) = NBNOVI
      KGCPRE    = NBSOM
      KGCUMU(1) = NBNOVI
      KGCUMU(2) = NBNOVI
      KGCUMU(3) = NBNOVI

C     VITESSE MOY et MAX AU TEMPS PRECEDENT. ICI TEMPS INITIAL
      VITMOYn = VITMOY0
      VITMAXn = VITMAX0

C     DONNEES DU CALCUL DU SECOND MEMBRE DES TEMPERATURES
      IEMG   = 0
      IEKG   = 0
      IEBG   = 1
      NBRDTG = 0
      NORESOO= 0
      MNLPLI = 0
      MNLPCO = 0
      PENALI = 0D0
      NCOM = 0
      NCOK = 0

C     LES DONNEES POUR LE TRACE DES TEMPERATURES ET DES FLECHES VITESSES
      CALL VISEE0
      NOTYVI  = 0
      LORBITE = 0
      NETAXE  = 0
      CMFLEC  = 2.5

C     PAS DE TRACE des ARETES
      IAVARE = 0
C     PAS DE TRACE DES FACES
      IAVFAC = 0
C     PAS DE TRACE DES ARETES DES FACES AVEC LA COULEUR
      NCOUAF = -2
C     PAS DE TRACE DES ARETES DES FACES FRONTIERE
      NCOAFR = -2

C     POURCENTAGE DE REDUCTION DES FACES
      PREDUF = 30.

      IF( NDIM .EQ. 3 ) THEN

C        LA VISEE EN 3D:
C        LONGITUDE et LATITUDE
ccc            CALL LONLAT( -80.0, 10. )
ccc            CALL LONLAT( -82.0, 8. )
ccc            CALL LONLAT( -110., 20. )
         CALL LONLAT( -85., 3. )

C        LOUPE GROSSISSANTE
         GROSSI = 0.75
         AXOLAR = AXOLAR / GROSSI
         AXOHAU = AXOHAU / GROSSI

      ENDIF


C     ##################################################################
C     ##                                                              ##
C     ##    LA BOUCLE EN TEMPS AVEC UN PAS DE TEMPS CONSTANT = DT     ##
C     ##                                                              ##
C     ##################################################################

C     NUMERO du PAS DE TEMPS DT A CALCULER
 100  NBPAST = NBPAST + 1

C     LES DONNEES SONT au TEMPS tn=TEMPS et
C     les CALCULS SUIVANTS pour le TEMPS tn+1 = tn + DT
C     ............................................................
      tn    = TEMPS
      TEMPS = TEMPS + DT
      tn1   = TEMPS

      IF( LANGAG .EQ. 0 ) THEN
         PRINT 10101, NBPAST, MXPAST-1, TEMPS
      ELSE
         PRINT 20101, NBPAST, MXPAST-1, TEMPS
      ENDIF
10101 FORMAT(/'nst3th: Au Pas',I7,'/',I7,' Temps ',G14.6,
     %        '  Calcul de la VITESSE+PRESSION+TEMPERATURE ',80('='))
20101 FORMAT(/'nst3th: At Step',I7,'/',I7,' Time ',G14.6,
     %        '  VELOCITY PRESSURE TEMPERATURE computation ',80('='))

ccc      PRINT*,'sin(Pi*(t/2.5))=',sin( Pi*(temps/2.5) ),
ccc     %      ' sin(Pi*(t/2.5 + 2./3))=',sin( Pi*(temps/2.5 + 2./3.) ),
ccc     %      ' sin(Pi*(t/2.5 + 4./3))=',sin( Pi*(temps/2.5 + 4./3.) )

ccc      call affvect( 'V3P0', 5,  VXVYVZPR(1,  NVP0) )
ccc      call affvect( 'PRES0',5,  VXVYVZPR(N1P,NVP0) )

      IF( NBVCFX .GT. 0 ) THEN

C        CALCUL DES DONNEES VITESSES0 CONVECTEES U(tn,X(tn;tn+1,NOEUDS))
C        ---------------------------------------------------------------
         CALL CONVECTH( tn, tn1, NDIM,  NBNOVI, MCN(MNXYZN+WYZNOE),
     %                  NBNOEF, NBELEM, MCN(MNELE+WUNDEL), NO1EFN,
     %                  MOARET, MXARET, MNLARE,
     %                  MOFACE, MXFACE, MCN(MNLFAC),
     %                  NDDLNO, VXYZPN, VITMAXn,
     %                  VITC,   IERR )

      ENDIF

C     =================================================================
C     1. LE CALCUL DES TEMPERATURES te(tn+1,x) AUX NOEUDS P2 DU FLUIDE
C     ( Rho Cp - dt Conduc Laplacien ) {te(tn+1,x)} = [TG] {te(tn+1,x)}
C      =Rho Cp {te(X(tn;tn+1,x)) + dt {FGte(tn+1,x)})
C     =================================================================

C     CALCUL DU SECOND MEMBRE GLOBAL DE L'EQUATION A L'INSTANT TEMPS tn+1
C     CONSTRUCTION DU SECOND MEMBRE GLOBAL DU SYSTEME POUR LA TEMPERATURE

C     1.1 CONTRIBUTION a BG des SOURCES de CHALEUR INTERNES et EXTERNES
C     dt (Integrale tP2 Source(tn+1) dX + Integrale SOURCE sur les FRONTIERES)
C     ------------------------------------------------------------------------
C     ATTENTION: Ne PAS DECLARER de TRANSPORT de la CHALEUR par la VITESSE
C                car CET ASPECT EST TRAITE ENSUITE PAR LA METHODE
C                DES CARACTERISTIQUES RETROGRADES
      CALL THEMKB( NBJEUX, IEMG,   IEKG,   IEBG,  PENALI,
     %             D2PI,   NDIM,   NTDLTE, VXVYVZPR(N1X,NVP0),
     %             NBTYEL, MNNPEF, NDPGST,
     %             MNTPOB, MXPOBA, MNTAUX,
     %             MNXYZN, NUMIOB, NUMAOB, MNDTEL,
     %             MNTHER, MNTAEL, MNX,    MNNODL,
     %             NORESOO,MNLPLI, MNLPCO,
     %             NBRDTG, TG,     NBRDTG, TG,   MNBG,
     %             NCOM,   NCOK,   NBPTAF, IERR )

      DO K=0,NTDLTE-1
         DMCN(MNBGD+K) = DMCN(MNBGD+K) * DELTAT
      ENDDO
ccc      call affvect('Fin THEMKB   BG:', 5,      DMCN(MNBGD) )
      call afl1ve( 'Fin THEMKB   BG:', NTDLTE, DMCN(MNBGD) )


C     1.2 CONTRIBUTION a BG du TRANSPORT de la TEMPERATURE par
C         la VITESSE U(tn) et la METHODE DES CARACTERISTIQUES RETROGRADES
C       + Integrale tP2 Rho Cp te(tn,X(tn;tn+1,x)) dx 
C     -------------------------------------------------------------------
C     RhoCp: DENSITE DE MASSE*CHALEUR SPECIFIQUE=CAPACITE THERMIQUE
      RhoCp = Rho*Cp
      CALL BGTHBOUS( tn, tn1,RhoCp,NBNOVI,NDDLNO,RMCN(MNXYZN+WYZNOE),
     %               NUTYEL, NBNOEF, NBELEM, MCN(MNELE+WUNDEL),
     %               MOARET, MXARET, MNLARE,
     %               MOFACE, MXFACE, MCN(MNLFAC),
     %               VXYZPN, VITMAXn, DMCN(MNTHETnD),
     %               DMCN(MNBGD), NBCHTROL, IERR )
ccc      call affvect('Fin BGTHBOUS BG:', 5,      DMCN(MNBGD) )
      call afl1ve( 'Fin BGTHBOUS BG:', NTDLTE, DMCN(MNBGD) )

      IF( IERR .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'nst3th: Sortie de bgthbous avec IERR=',IERR
         ELSE
            PRINT*,'nst3th: Problem exit bgthbous with IERR=',IERR
         ENDIF
         GOTO 9999
      ENDIF

      IF( IERR .NE. 0 ) GOTO 9999

C     1.3 PRISE EN COMPTE DES CONDITIONS AUX LIMITES DE LA TEMPERATURE
C     ----------------------------------------------------------------  
C     LISTE DES NUMEROS ET VALEURS DES TEMPERATURES FIXEES AU TEMPS tn+1
      CALL THDLFX( 1,      NTDLTE, NDIM,
     %             NBTYEL, MNNPEF, NDPGST,
     %             MNXYZN, NUMIOB, MNDTEL, RELMIN,
     %             NBTEFX, MONTEFX, MNNTEFX, MNVTEFX, IERR )

C     PRISE EN COMPTE DES CONDITIONS AUX LIMITES de TEMPERATURE SUR MNBG
      CALL BLDLFX( NTDLTE, 1, NBTEFX, MCN(MNNTEFX), MCN(MNVTEFX),
     %             STGV, DMCN(MNBGD) )

ccc      call affvect('Fin THDLFX BG:', 5,      DMCN(MNBGD) )
ccc      call afl1ve( 'Fin THDLFX BG:', NTDLTE, DMCN(MNBGD) )
cccC     CARRE DE LA NORME DE BG INITIAL
Cccc     DOUBLE PRECISION PROSCD
ccc      print *,'nst3th: (bg,bg)=',
ccc     %         PROSCD( DMCN(MNBGD), DMCN(MNBGD), NTDLTE )

C     LES TEMPERATURES MINIMALES et MAXIMALES FIXEES au TEMPS tn+1
      TEMMINFX =  1D100
      TEMMAXFX = -1D100
      MNT0 = ( MNVTEFX - 1 ) / MOREE2
      DO I=1, NBTEFX
C        LE NO GLOBAL DU DL TEMPERATURE FIXE
         NDL = MCN( MNNTEFX - 1 + I )
C        VALEUR FIXEE DE LA TEMPERATURE(tn+1,ndl)
         TE = DMCN( MNT0 + I )
         IF( TE .LT. TEMMINFX ) THEN
            TEMMINFX = TE
         ENDIF
         IF( TE .GT. TEMMAXFX ) THEN
            TEMMAXFX = TE
         ENDIF
      ENDDO

C     1.4 RESOLUTION DU SYSTEME LINEAIRE -> Te(tn+1) = DMCN(MNTHETD)
C     --------------------------------------------------------------
      IF( NORESO .EQ. 1 ) THEN

C        DESCENTE et REMONTEE DU PROFIL DE LA MATRICE TEMPERATURE TG
C        L D tL COMPLETE
         CALL DRCRPR( NTDLTE, NCODST, LP2LIGN, TG, DMCN(MNBGD), 3,
     %                DMCN(MNTHETD), IERR )

ccc         call affvect('Fin drcrpr: TEMPERATURE n',5,   DMCN(MNTHETD))
         call afl1ve( 'Fin drcrpr: TEMPERATURE n', NTDLTE,DMCN(MNTHETD))

      ELSE

C        RESOLUTION PAR GC SIMPLE DE LA TEMPERATURE
C        LE VECTEUR INITIAL DU GC EST LE DERNIER VECTEUR TEMPERATURE CALCULE
C        EN SORTIE DE GCAXBK LA SOLUTION EST DANS TEMPER tn+1 = DMCN(MNTHETD)
         CALL GCAXBK( NTDLTE, NODLTEFX, LP2LIGN, LP2COLO, TG,
     %                DMCN(MNBGD), DMCN(MNTHETnD),
     %                DAUX1, DAUX2, DAUX3,
     %                DMCN(MNTHETD), KGCTEMP )

ccc         call affvect('Fin GcAxbk: TEMPERATURE n',5,   DMCN(MNTHETD))
         call afl1ve( 'Fin GcAxbk: TEMPERATURE n', NTDLTE,DMCN(MNTHETD))

C        VERIFICATION DE LA CONVERGENCE DU GC POUR LA TEMPERATURE
         IF( IERR .NE. 0 ) THEN
C           IL N'Y A PAS EU CONVERGENCE DU GC
C           LA METHODE DU GC NE CONVERGE PAS => ABANDON DES CALCULS
            WRITE(KERR(MXLGER)(1:10),'(I10)') NTDLTE
            NBLGRC(NRERR) = 3
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'NON CONVERGENCE apres '//KERR(MXLGER)(1:10)
     %                   // ' ITERATIONS du GC'
               KERR(2) = 'ABANDON DE LA METHODE DU GRADIENT CONJUGUE'
               KERR(3) = 'ESSAYER UNE AUTRE METHODE DE RESOLUTION'
            ELSE
               KERR(1) = 'NO CONVERGENCE after '// KERR(MXLGER)(1:10)
     %                   // ' ITERATIONS of CONJUGATE GRADIENT'
               KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
               KERR(3) = 'TRY an OTHER METHOD'
            ENDIF
            CALL LEREUR
            IERR = 9
            GOTO 9999
         ENDIF

      ENDIF

C     ECRETAGE DE LA TEMPERATURE tn+1 A LA TEMPERATURE MIN MAX PRECEDENTE
C     -------------------------------------------------------------------
      TEMMINECR = MIN( TEMMINt0, TEMMINFX )
ccc      TEMMAXECR = MAX( 0.95D0 * TEMMAX, 1.05D0 * TEMMAX, TEMMAXFX )
      TEMMAXECR = MAX( TEMMAXt0, TEMMAXFX )

      NBTECRMIN = 0
      NBTECRMAX = 0
      MNT = MNTHETD - 1
      DO K=1,NTDLTE
         MNT = MNT + 1
         TE  = DMCN(MNT)
         IF( TE .LT. TEMMINECR ) THEN
            NBTECRMIN = NBTECRMIN + 1
C           ECRETAGE AU MINIMUM DE LA TEMPERATURE 
            DMCN(MNT) = TEMMINECR
         ENDIF
         IF( TE .GT. TEMMAXECR ) THEN
            NBTECRMAX = NBTECRMAX + 1
C           ECRETAGE AU MAXIMUM DE LA TEMPERATURE 
            DMCN(MNT) = TEMMAXECR
         ENDIF
      ENDDO
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'Au TEMPS',TEMPS,': ',NBTECRMIN,
     %          ' TEMPERATURES ECRETEES au MIN',TEMMINECR,
     %          ' et',NBTECRMAX,' au MAX',TEMMAXECR
      ELSE
         PRINT*,'At TIME',TEMPS,': ',NBTECRMIN,
     %          ' ECRETED TEMPERATURES at MIN',TEMMINECR,
     %          ' and',NBTECRMAX,' at MAX',TEMMAXECR
      ENDIF

C     BLOCAGE FINAL DES TEMPERATURES FIXEES AU TEMPS tn+1
C     ---------------------------------------------------
      MNT0 = ( MNVTEFX - 1 ) / MOREE2
      MNT  = MNTHETD - 1
      DO I=1, NBTEFX
C        LE NO GLOBAL DU DL TEMPERATURE FIXE
         NDL = MCN( MNNTEFX - 1 + I )
C        VALEUR FIXEE DE LA TEMPERATURE(tn+1,ndl)
         DMCN( MNT+NDL ) = DMCN( MNT0 + I )
      ENDDO


C     TRACE DES ZONES DE COULEURS des TEMPERATURES tn+1 en 2D et ISOTHERMES en 3D
C     ---------------------------------------------------------------------------
      IF( INTERA .GE. 1 ) THEN

C        MODE GRAPHIQUE AVEC X11: LA MEMOIRE PIXELS EST EFFACEE
         CALL EFFACEMEMPX

         NCAS   = 1
         NCAMIN = 1
         NCAMAX = 1
         TMIN = REAL( TEMMINECR )
         TMAX = REAL( TEMMAXECR )

         IF( NDIM .EQ. 2 ) THEN

ccc         PRINT *,'trzont 2d Temps=',TEMPS,'  ---------->'
C           TRACE DES ZONES DE COULEURS DANS LE MAILLAGE 2D
            CALL TRZONT( 0,      NDIM,    KNOMOB, -11,
     %                   NBTYEL, MNNPEF,  MNXYZN, MNXYZN, NDPGST,
     %                   NCAS,   NCAS,    NTDLTE,
     %                   0,      DMCN(MNTHETD),  temp,
     %                   TMIN,   NOTMIN, NCAMIN, TMAX, NOTMAX, NCAMAX,
     %                   TEMPS )

         ELSE

ccc         PRINT *,'trisot 3d Temps=',TEMPS,'  ---------->'
C           TRACE DES SURFACES ISOTHERMES
            CALL TRISOT( NDIM,   KNOMOB, -11,
     %                   NBTYEL, MNNPEF, MNXYZN, NDPGST,
     %                   NCAS,   NCAS,   NTDLTE,
     %                   0,      DMCN(MNTHETD),  temp,
     %                   TMIN,   NOTMIN, NCAMIN, TMAX, NOTMAX, NCAMAX,
     %                   TEMPS )

         ENDIF

cccC        TEMPS UTILISATEUR EN SECONDES EN ATTENTE AVANT DE CONTINUER
ccc         CALL ATTENDSEC( 1D0 )
      ENDIF


C     ================================================================
C     2. EF TAYLOR-HOOD P2 EN VITESSE U* et P1 EN PRESSION
C        VG U*(tn+1) = ( Rho - dt Mhu Laplacien ) U*(tn+1) =
C      = Rho U(tn,X(tn;tn+1,x) - dt CoGrPr GRAD P(tn)
C      - dt Rho g CoBOUS (Temper(tn+1)-Temper(t0)) axe Z  (ou Y en 2d)
C      + dt Force(tn+1)
C     ================================================================

C     CALCUL du SECOND MEMBRE DU SYSTEME LINEAIRE:
C     Integrale Rho tP2 Ui(tn,X(tn;tn+1,bl) dX
C   - Integrale tP2 dt CoGrPr GRAD P(tn) dX
C   + Integrale dt tP2 Force(tn+1) dX
C   - Integrale Rho g tP2 CoBOUS (Temper(tn+1)-Temper(t0)) dX axe -Z (ou -Y 2d)
C     Le dernier terme BOUSSINESQ PREND EN COMPTE la TEMPERATURE
C     ------------------------------------------------------------------------
      CALL BGUIBOUS( tn, tn1, Rho, G, CoGrPr, CoBOUS,
     %               NDIM,   NBSOM,  NBNOVI, NTDLVI,
     %               NDDLNO, NONOSO, RMCN(MNXYZN+WYZNOE),
     %               NUTYEL, NBNOEF, NBELEM, MCN(MNELE+WUNDEL),
     %               NBVVEF, NBSFEF, NBLAEF,
     %               MCN(MNVVEF), MCN(MNSFEF), MCN(MNLAEF),
     %               NUMIOB(2), NUMAOB(2), MCN(MNDOEL(2)),
     %               NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)),
     %               NUMIOB(4), NUMAOB(4), MCN(MNDOEL(4)),
     %               MOARET, MXARET, MNLARE,
     %               MOFACE, MXFACE, MCN(MNLFAC),
     %               VXYZPN, VITMAXn, VXVYVZPR(N1P,NVP0),
     %               TEMPER0,DMCN(MNTHETD),
     %               P2P22D, P2P23D,
     %               V3P2,   NBCHTROL, IERR )
      IF( IERR .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'nst3th: Sortie PB de bguibous avec IERR=',IERR
         ELSE
            PRINT*,'nst3th: Problem exit bguibous with IERR=',IERR
         ENDIF
         GOTO 9999
      ENDIF

C     CONSTRUCTION DE LA LISTE DES DL VITESSES FIXEES AU TEMPS tn+1
C     -------------------------------------------------------------
      CALL VITEFX( RELMIN, NTDLVI, NDIM,   MNXYZN,
     %             NBTYEL, MNNPEF, NUMIOB, MNDOEL,
     %             NBVCFX, NBVIFX, MNNVIFX, MNVVIFX, VITFXMAX )
      VITMAXn = MAX( VITMAXn, VITFXMAX )

C     BLOCAGE DES COMPOSANTES FIXEES DE LA VITESSE AU TEMPS tn+1 SUR V3P2
C     -------------------------------------------------------------------
      CALL VITEFXBG( NBVIFX, MCN(MNNVIFX), MCN(MNVVIFX), STGV,
     %               VITC,   NDIM, NBNOVI, V3P2 )

C     RESOLUTION des NDIM COMPOSANTES de la VITESSE U*i: i=1,...,NDIM
C     Integrale  Rho tP2 P2 + dt Mhu tGrad P2 Grad P2 dx U*i(tn+1)
C   = Integrale  Rho tP2 Ui(tn,X(tn;tn+1,bl) dX
C   - Integrale  dt CoGrPr tP2 Gradi P(tn) dX
C   - Integrale Rho g tP2 CoBOUS (Temper(tn+1)-Temper(t0)) dX axe Z(ou Y en 2d)
C   + Integrale  dt        tP2 Forcei(tn)  dX     => U*(tn+1)
C     ---------------------------------------------------------------------
      MV = 1
      IF( NORESO .EQ. 1 ) THEN

C        DESCENTE et REMONTEE DU PROFIL DE LA MATRICE VITESSE VG
C        L D tL COMPLETE DES NDIM COMPOSANTES DE LA VITESSE U*(tn+1)
         DO K = 1, NDIM
            CALL DRCRPR( NBNOVI, NCODSV, LP2LIGN, VG, V3P2(MV), 3,
     %                   VXVYVZPR(MV,NVP1), IERR )

ccc         print *,'Composante ',K,' de U*  -------------------------'
ccc        call affvect('Apres DRCRPR: U*i',5,       VXVYVZPR(MV,NVP1) )
           call afl1ve( 'Apres DRCRPR: U*i', NBNOVI, VXVYVZPR(MV,NVP1) )
            MV = MV + NBNOVI
         ENDDO

      ELSE

C        RESOLUTION PAR GC DES NDIM COMPOSANTES DE LA VITESSE U*(tn+1)
C        LE VECTEUR INITIAL DU GC EST LE DERNIER VECTEUR CALCULE
C        EN SORTIE DE GCAXBK LA SOLUTION EST DANS VXVYVZPR
         DO K = 1, NDIM
            CALL GCAXBK( NBNOVI, NODLVIFX, LP2LIGN, LP2COLO, VG,
     %                   V3P2(MV), VXVYVZPR(MV,NVP0),
     %                   DAUX1, DAUX2, DAUX3,
     %                   VXVYVZPR(MV,NVP1),   KGCUET(K) )

ccc         print *,'Composante ',K,' de U*  ----------------------'
ccc         call affvect('Apres GCAxbk: U*i',5,      VXVYVZPR(MV,NVP1) )
            call afl1ve( 'Apres GCAxbk: U*i', NBNOVI,VXVYVZPR(MV,NVP1) )
            MV = MV + NBNOVI
         ENDDO

      ENDIF


C     ==========================================================================
C     3. CALCUL DE LA DIFFERENCE DE PRESSION P(tn+1)-P(tn) par la RESOLUTION de
C     -dt CoGrPr LAPLACIEN(P(tn+1)-P(tn))=-Rho (Div U* -(Integrale Div U*)/Vol))
C     -dt CoGrPr      dP(tn+1)-P(tn))/dn = 0
C     Ici le coefficient dt Mhu est SUPPOSE TRES PETIT
C     L'integrale Div U* /Vol IMPOSE MODEREMENT, EN MOYENNE VOLUMIQUE,  Div U*=0
C     ==========================================================================
C     CALCUL DU SECOND MEMBRE  BG = - Rho (Div U* -(Integrale Div U*)/Volume))
      IF( NUTYEL .EQ. 15 ) THEN

C        TRIANGLE TAYLOR-HOOD
C        CALCUL DE Integrale Div VITXY dX et Integrale dX sur le MAILLAGE
         CALL MOYDIVTH2D( NBNOVI, RMCN(MNXYZN+WYZNOE),
     %                    NBELEM, MCN(MNELE+WUNDEL), VXVYVZPR(1,NVP1),
     %                    VOLUME, INTDIVV )
C        Integrale (-Rho) * tP1 ( Div VITXYZ - Integrale Div VITXYZ / VOLUME )
         CALL BGDIVTH2( NBSOM,  RMCN(MNXYZN+WYZNOE),
     %                  NBNOEF, NBELEM,  MCN(MNELE+WUNDEL), NONOSO,
     %                  -Rho,   NBNOVI,  VXVYVZPR(  1,NVP1),
     %                  VOLUME, INTDIVV, VXVYVZPR(N1P,NVP1) )

      ELSE

C        TETRAEDRE TAYLOR-HOOD
C        CALCUL DE Integrale Div VITXYZ dX et Integrale dX sur le MAILLAGE
         CALL MOYDIVTH3D( NBNOVI, RMCN(MNXYZN+WYZNOE),
     %                    NBELEM, MCN(MNELE+WUNDEL), VXVYVZPR(1,NVP1),
     %                    VOLUME, INTDIVV )
C        Integrale (-Rho) * tP1 ( Div VITXYZ - Integrale Div VITXYZ / VOLUME )
         CALL BGDIVTH3( NBSOM,  RMCN(MNXYZN+WYZNOE),
     %                  NBNOEF, NBELEM, MCN(MNELE+WUNDEL), NONOSO,
     %                  -Rho,   NBNOVI,  VXVYVZPR(  1,NVP1),
     %                  VOLUME, INTDIVV, VXVYVZPR(N1P,NVP1) )

      ENDIF


cccC     =========================================================================
cccC     3. CALCUL DE LA DIFFERENCE DE PRESSION P(tn+1)-P(tn) PAR LA RESOLUTION de
cccC     -dt CoGrPr LAPLACIEN(P(tn+1)-P(tn)) = (Rho-dt Mhu Laplacien) (-Div U*)
cccC     -dt CoGrPr       dP(tn+1)-P(tn))/dn = dt Mhu      Laplacien    n . U*
cccC     =========================================================================
ccc      dtMhu = DELTAT * Mhu
ccc      IF( NUTYEL .EQ. 15 ) THEN
cccC        TRIANGLE TAYLOR-HOOD
ccc         CALL BGDIVLATH2( Rho,    dtMhu,  VXVYVZPR(1,NVP1),
ccc     %                    NBNOVI, MCN(MNXYZN+WYZNOE),
ccc     %                    NONOSO, NBELEM, MCN(MNELE+WUNDEL),
ccc     %                    MCN(MNLAEF), NBSOM,    VXVYVZPR(N1P,NVP1) )
ccc      ELSE
cccC        TETRAEDRE TAYLOR-HOOD
ccc         CALL BGDIVLATH3( Rho,    dtMhu,  VXVYVZPR(1,NVP1),
ccc     %                    NBNOVI, MCN(MNXYZN+WYZNOE),
ccc     %                    NONOSO, NBELEM, MCN(MNELE+WUNDEL),
ccc     %                    MCN(MNSFEF), NBSOM,    VXVYVZPR(N1P,NVP1) )
ccc      ENDIF

ccc      call affvect( 'Sd Membre BGDIVTH',5,      VXVYVZPR(N1P,NVP1) )
ccc      call afl1ve(  'Sd Membre BGDIVTH', NBSOM, VXVYVZPR(N1P,NVP1) )


C     CONSTRUCTION DES TABLEAUX NO ET VALEURS DES PRESSIONS FIXEES a tn+1
C     -------------------------------------------------------------------
      IF( IEBLPR .GT. 0 ) THEN
         CALL PRESFXST( RELMIN, NBSOM,   NDIM,   MNXYZN, NONOSO,
     %                  NBTYEL, MNNPEF,  NUMIOB, MNDOEL,
     %                  NBPRFX, MNNPRFX, MNVPRFX )
      ELSE
         NBPRFX = 0
      ENDIF

      IF( NORESO .EQ. 1 ) THEN

C        PRISE EN COMPTE DES DEGRES DE LIBERTE PRESSION FIXES P(tn+1)-P(tn)
C        ------------------------------------------------------------------
         MNV = ( MNVPRFX - 1 ) / MOREE2
         DO I=1, NBPRFX
C           LE NO SOMMET GLOBAL DU DL FIXE
            NDL = MCN( MNNPRFX - 1 + I )
C           VALEUR FIXEE DE LA PRESSION  P(tn+1)-P(tn)
            VXVYVZPR(N0P+NDL,NVP1) = STGV *
     %                           ( DMCN(MNV+I)-VXVYVZPR(N0P+NDL,NVP0) )
         ENDDO

C        RESOLUTION de P(tn+1)-P(tn) = [L D tL]-1 (-(Rho-dt Mhu Laplacien)Div U*)
C        PAR DESCENTE REMONTEE DE LA MATRICE FACTORISEE PG = L D tL PROFIL
C        ------------------------------------------------------------------------
         CALL DRCRPR( NBSOM, NCODSP, LP1LIGN, PG, VXVYVZPR(N1P,NVP1), 3,
     %                Pr1Pr, IERR )

ccc      call affvect('Apres DRCRPR P(tn+1)-P(tn)=', 5,     Pr1Pr )
         call afl1ve( 'Apres DRCRPR P(tn+1)-P(tn)=', NBSOM, Pr1Pr )

      ELSE

C        COPIE DE PRESSION1 DANS PVAUX POUR LE GC
         CALL TRTATD( VXVYVZPR(N1P,NVP1), PVAUX, NBSOM )

C        PRISE EN COMPTE DES DEGRES DE LIBERTE PRESSION FIXES P(tn+1)-P(tn)
C        ------------------------------------------------------------------
         MNV = ( MNVPRFX - 1 ) / MOREE2
         DO I=1, NBPRFX
C           LE NO SOMMET GLOBAL DU DL FIXE
            NDL = MCN( MNNPRFX - 1 + I )
C           VALEUR FIXEE DE LA PRESSION  P(tn+1)-P(tn)
            PVAUX(NDL) = STGV * ( DMCN(MNV+I) - VXVYVZPR(N0P+NDL,NVP0) )
         ENDDO

C        --------------------------------------------------------------------
C        RESOLUTION de P(tn+1)-P(tn) = [PG]-1 (-(Rho-dt Mhu Laplacien)Div U*)
C        PAR GRADIENT CONJUGUE AVEC MATRICE SYMETRIQUE M
C        --------------------------------------------------------------------
C        RESOLUTION PAR GC SIMPLE
C        LE VECTEUR INITIAL DU GC EST LE DERNIER VECTEUR PRESSION Pr1Pr CALCULE
         CALL GCAxbk( NBSOM, NODLPRFX, LP1LIGN, LP1COLO, PG,
     %                PVAUX, Pr1Pr,
     %                DAUX1, DAUX2, DAUX3,
     %                Pr1Pr, KGCPRE )

ccc      call affvect('Apres GCAxbk P(tn+1)-P(tn)=', 5,     Pr1Pr )
         call afl1ve( 'Apres GCAxbk P(tn+1)-P(tn)=', NBSOM, Pr1Pr )

      ENDIF


C     =========================================================================
C     4. (Rho-dt Mhu Laplacien)(ui(tn+1)-u*i)=- dt CoGrPr gradi(p(tn+1)- p(tn))
C     V3P2i = - dt CoGrPr gradi (p(tn+1)- p(tn))  i=1,...,3
C     =========================================================================
      CALL BGGRATH( DELTAT, CoGrPr, NDIM,  NBSOM, NBNOVI, NTDLVI,
     %              NONOSO, MCN(MNXYZN+WYZNOE),
     %              NUTYEL, NBNOEF, NBELEM, MCN(MNELE+WUNDEL),
     %              Pr1Pr,  V3P2 )
ccc      call affvect( 'BGGRATH avant CL VITESSE',5,           V3P2 )
ccc      call afl1ve(  'BGGRATH avant CL VITESSE', NDIM*NBNOVI,V3P2 )

C     BLOCAGE DES DL VITESSES  U(tn+1)-U*  FIXEES A ZERO
C     CAR U*=U(tn+1) sur la FRONTIERE DIRICHLET DU PROBLEME
C     -----------------------------------------------------
      DO I=1, NBVIFX
C        LE NO GLOBAL DU DL FIXE
         NDL = MCN( MNNVIFX - 1 + I )
C        VALEUR FIXEE DE LA VITESSE  ( Ui(tn+1)-Ui* )=0
         V3P2( NDL ) = 0D0
      ENDDO

C     RESOLUTION des NDIM COMPOSANTES DE LA VITESSE Ui(tn+1)-U*i: i=1,...,NDIM
C     Integrale Rho tP2 P2 + dt Mhu tGrad P2 Grad P2 dx (Ui(tn+1)-U*i)
C  = -Integrale dt CoGrPr tP2 Grad (P(tn+1)-P(tn)) dX   =>  Ui(tn+1)-U*i
C     ------------------------------------------------------------------------
      MV = 1
      IF( NORESO .EQ. 1 ) THEN

C        DESCENTE et REMONTEE DU PROFIL DE LA MATRICE VITESSE
C        L D tL COMPLETE DES NDIM COMPOSANTES DE LA VITESSE U*(tn+1)
         DO K = 1, NDIM
            CALL DRCRPR( NBNOVI, NCODSV, LP2LIGN, VG, V3P2(MV), 3,
     %                   V3P2(MV), IERR )
ccc         call affvect('Apres DRCRPR: (U(tn+1)-U*)i',    5, V3P2(MV) )
            call afl1ve( 'Apres DRCRPR: (U(tn+1)-U*)i',NBSOM, V3P2(MV) )
            MV = MV + NBNOVI
         ENDDO

      ELSE

C        RESOLUTION PAR GC DES NDIM COMPOSANTES DE LA VITESSE U(tn+1)-U*
         DO K = 1, NDIM

C           COPIE DE V3P2 DANS PVAUX POUR LE PROTEGER
C           CAR B DOIT ETRE DIFFERENT DE X DANS LE GC
            CALL TRTATD( V3P2(MV), PVAUX, NBNOVI )

C           MISE A ZERO DU VECTEUR INITIAL DU GRADIENT CONJUGUE
            CALL AZEROD( NBNOVI, V3P2(MV) )

C           RESOLUTION PAR GRADIENT CONJUGUE
            CALL GCAXBK( NBNOVI, NODLVIFX, LP2LIGN, LP2COLO, VG,
     %                   PVAUX,  V3P2(MV),
     %                   DAUX1,  DAUX2,  DAUX3,
     %                   V3P2(MV),  KGCUMU(K) )

cccC           VERIFICATION DE LA CONVERGENCE DU GC
ccc            IF( IERR .NE. 0 ) THEN
cccC              IL N'Y A PAS EU CONVERGENCE DU GC
cccC              LA METHODE DU GC NE CONVERGE PAS => ABANDON DES CALCULS
ccc               WRITE(KERR(MXLGER)(1:10),'(I10)') NBNOVI
ccc               NBLGRC(NRERR) = 3
ccc               IF( LANGAG .EQ. 0 ) THEN
ccc                  KERR(1) = 'NON CONVERGENCE apres '//KERR(MXLGER)(1:10)
ccc     %                   // ' ITERATIONS DU GRADIENT CONJUGUE'
ccc                  KERR(2) = 'ABANDON DE LA METHODE DU GRADIENT CONJUGUE'
ccc                  KERR(3) = 'ESSAYER UNE AUTRE METHODE DE RESOLUTION'
ccc               ELSE
ccc                  KERR(1) = 'NO CONVERGENCE after '// KERR(MXLGER)(1:10)
ccc     %                   // ' ITERATIONS of CONJUGATE GRADIENT'
ccc                  KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
ccc                  KERR(3) = 'TRY an OTHER METHOD'
ccc               ENDIF
ccc               CALL LEREUR
ccc               IERR = 9
ccc               GOTO 199
ccc            ENDIF

ccc            call affvect('Apres GCAxbk (U(tn+1)-U*)i=', 5,V3P2(MV))
            call afl1ve( 'Apres GCAxbk (U(tn+1)-U*)i=',NBNOVI,V3P2(MV))
            MV = MV + NBNOVI

         ENDDO

      ENDIF
      PRINT*


C     ============================================================
C     5. U(tn+1) = U* + (U(tn+1)-U*)  => U(tn+1)=VXVYVZPR(MV,NVP1)
C     ============================================================
      CALL SOM2VED( NTDLVI, VXVYVZPR(1,NVP1), V3P2, VXVYVZPR(1,NVP1) )
      MV = 1
      DO K = 1, NDIM
         PRINT*,' COMPOSANTE',K,' U(tn+1) = U* + (U(tn+1)-U*) ---------'
ccc         call affvect( 'U(tn+1) = U* + (U(tn+1)-U*)',5,
ccc     %                 VXVYVZPR(MV,NVP1))
         call afl1ve('U(tn+1) = U* + (U(tn+1)-U*)',NBNOVI,
     %             VXVYVZPR(MV,NVP1))
         MV = MV + NBNOVI
      ENDDO

C     BLOCAGE DES DL VITESSES U(tn+1) FIXEES
C     --------------------------------------
      MNV0 = ( MNVVIFX - 1 ) / MOREE2      
      DO I=1, NBVIFX
C        LE NO GLOBAL DU DL FIXE
         NDL = MCN( MNNVIFX - 1 + I )
C        VALEUR FIXEE DE LA VITESSE  U(tn+1,ndl)
         VXVYVZPR(NDL,NVP1) = DMCN( MNV0 + I )
      ENDDO
      PRINT*


C     ================================================================
C     6. CALCUL DE P(tn+1) = P(tn) + (P(tn+1)-P(tn))  => P(tn+1)=PRES1
C     ================================================================
      CALL SOM2VED( NBSOM, VXVYVZPR(N1P,NVP0), Pr1Pr,
     %                     VXVYVZPR(N1P,NVP1) )
ccc      call affvect( 'P(tn+1)=P(tn) +(P(tn+1)-P(tn))', 5,
ccc     %               VXVYVZPR(N1P,NVP1))
      call afl1ve(  'P(tn+1)=P(tn) +(P(tn+1)-P(tn))',
     %              NBSOM, VXVYVZPR(N1P,NVP1) )

      IF ( NBPRFX .GT. 0 ) THEN

C        PRISE EN COMPTE DES DEGRES DE LIBERTE PRESSION FIXES
C        ----------------------------------------------------
         MNP0 = ( MNVPRFX - 1 ) / MOREE2
         DO I=1, NBPRFX
C           LE NO GLOBAL DU DL FIXE
            NDL = MCN( MNNPRFX - 1 + I )
C           VALEUR FIXEE DE LA PRESSION
            VXVYVZPR( N0P+NDL, NVP1 ) = DMCN( MNP0 + I )
         ENDDO

      ELSE

C        CALCUL DE LA MOYENNE VOLUMIQUE DE LA PRESSION
C        POUR LUI IMPOSER UNE MOYENNE VOLUMIQUE NULLE SUR LE MAILLAGE
C        ------------------------------------------------------------
         IF( NUTYEL .EQ. 15 ) THEN
C           TRIANGLE TAYLOR-HOOD
            CALL MOYP12D( NBNOVI, RMCN(MNXYZN+WYZNOE),
     %                    NBELEM, NBNOEF, MCN(MNELE+WUNDEL),NONOSO,
     %                    NBSOM,  VXVYVZPR(N1P,NVP1), VOLUME, INTPRES )
         ELSE
C           TETRAEDRE TAYLOR-HOOD
            CALL MOYP13D( NBNOVI, RMCN(MNXYZN+WYZNOE),
     %                    NBELEM, NBNOEF, MCN(MNELE+WUNDEL),NONOSO,
     %                    NBSOM,  VXVYVZPR(N1P,NVP1), VOLUME, INTPRES )
         ENDIF

C        TRANSLATION POUR IMPOSER UNE MOYENNE VOLUMIQUE NULLE DE LA PRESSION
         INTPRES = INTPRES / VOLUME
         DO NN = 1, NBSOM
            VXVYVZPR(N0P+NN,NVP1) = VXVYVZPR(N0P+NN,NVP1) - INTPRES
         ENDDO

cccC        PAS DE CONDITION AUX LIMITES SUR LA PRESSION
cccC        => ELLE EST DEFINIE A UNE CONSTANTE PRES
cccC        RECHERCHE DE LA PLUS BASSE PRESSION POUR SERVIR DE TRANSLATION
cccC        ELLE EST SOUSTRAITE DES AUTRES PRESSIONS ET DEVIENT ZERO
cccC        --------------------------------------------------------------
ccc         PREMIN = 1D100
ccc         DO N=1,NBSOM
ccc            IF( VXVYVZPR(N0P+N,NVP1) .LT. PREMIN ) THEN
ccc               PREMIN = VXVYVZPR(N0P+N,NVP1)
ccc            ENDIF
ccc         ENDDO
cccC        TRANSLATION DE PREMIN A ZERO
ccc         DO N=1,NBSOM
ccc            VXVYVZPR(N0P+N,NVP1) = VXVYVZPR(N0P+N,NVP1) - PREMIN
ccc         ENDDO
ccc         IF( LANGAG .EQ. 0 ) THEN
ccc            PRINT*, 'Au TEMPS',TEMPS,' PRESSION MIN',PREMIN,
ccc     %                      ' est TRANSLATEE a ZERO'
ccc         ELSE
ccc            PRINT*, 'At TIME',TEMPS,' PRESSURE MIN',PREMIN,
ccc     %                      ' is TRANSLATED to ZERO'
ccc         ENDIF

      ENDIF

ccc      call affvect( 'PRESSION tn+1 FINAL',     5, VXVYVZPR(N1P,NVP1) )
ccc      call afl1ve(  'PRESSION tn+1 FINAL', NBSOM, VXVYVZPR(N1P,NVP1) )


C     BILAN: CALCUL DES |VITESSES tn+1| MIN ET MAX ET MOYENNE
C     -------------------------------------------------------
      CALL MAXVIT( NDIM, NBNOVI, 1,
     %    VXVYVZPR(N1X,NVP1), VXVYVZPR(N1Y,NVP1), VXVYVZPR(N1Z,NVP1),
     %             NOVMIN, NCVMIN, VITMIN,
     %             NOVMAX, NCVMAX, VITMAX, VITMOY )

ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         PRINT 10050, TEMPS, VITMOY, VITMAX, VITMIN, NBPAST
ccc      ELSE
ccc         PRINT 20050, TEMPS, VITMOY, VITMAX, VITMIN, NBPAST
ccc      ENDIF
ccc10050 FORMAT(' Au Temps ', G14.7,
ccc     %' |VITESSE|Moyenne=',G14.7,
ccc     %' |VITESSE|Max=',    G14.7,
ccc     %' |VITESSE|Min=',    G14.7, ' Pas de temps',I7)
ccc20050 FORMAT(' At Time ',  G14.7,
ccc     %' |VELOCITY| Mean=', G14.7,
ccc     %' |VELOCITY| Max=',  G14.7,
ccc     %' |VELOCITY| Min=',  G14.7, ' Time Step',I7)


C     =========================================================================
C     PAS D'ECRETAGE ............... supprime le  13/4/2021 et remis 1/10/2022
C                                 et supprime ici
      IF( NBSOM .LT. 0 ) THEN   
         IF( NBPAST .GT. 50 ) THEN
C           ECRETAGE DE LA VITESSE AU TEMPS tn+1 A LA VITESSE VITECR
C           apres 50 PAS DE TEMPS AU DEMARRAGE
C           --------------------------------------------------------
C           MISE A JOUR DU SEUIL D'ECRETAGE TOUTES LES 200 ITERATIONS
            IF( MOD( NBPAST, 200 ) .EQ. 0 ) THEN
ccc               VITECR = MIN( VITECR, VITMOY + VITMAX )
               VITECR = VITMAX
            ENDIF

cccC        ECRETAGE SI VITESSE MAX > 5% DE PLUS DE LA VALEUR MAXIMALE PRECEDENTE
ccc         IF( VITMAX .GT. 1.05D0 * VITMAXn ) THEN   14/4/2021
cccccc            VITECR = 1.02D0 * VITMAXn
ccc            VITECR = 1.01D0 * VITMAXn
ccc             GOTO 55
ccc         ENDIF

C        ECRETAGE SI  VITESSE MAX SUPERIEURE A LA VITESSE D'ECRETAGE
         IF( VITMAX .GT. VITECR ) THEN
            CALL VITECRET( NDIM, NBNOVI,
     %                     VXVYVZPR(N1X,NVP1), VXVYVZPR(N1Y,NVP1),
     %                     VXVYVZPR(N1Z,NVP1), VITECR )
         ENDIF
         ENDIF
      ENDIF
C     =========================================================================


C     PASSAGE DE VITX,VITY,VITZ,PRESSION PAR COMPOSANTES DE NBNOVI 
C                VITX,VITY,VITZ,PRESSION PAR NOEUDS
C     ===============================================================
      CALL VXYZPRE( NDIM, NUTYEL, NDDLNO, 1, NBNOVI,
     %        VXVYVZPR(N1X,NVP1), VXVYVZPR(N1Y,NVP1),VXVYVZPR(N1Z,NVP1),
     %              NBSOM, NONOSO, VXVYVZPR(N1P,NVP1), NTDLVP,
     %              VXYZPN )

C     AFFICHAGE DU VECTEUR"VITESSEPRESSION+TEMPERATURE AU TEMPS apres ECRETAGE
C     ------------------------------------------------------------------------
      CALL AFVIPRTE( NUTYEL, NDIM,   NBNOVI, MNXYZN, NDDLNO,
     %               TEMPS,  NTDLVP, MIN(1,NBNOVI),
     %               VXYZPN, DMCN(MNTHETD),
     %               VITMIN, NOVMIN, VITMAX, NOVMAX, VITMOY,
     %               PREMIN, NOPMIN, PREMAX, NOPMAX, PREMOY,
     %               TEMMIN, NOTMIN, TEMMAX, NOTMAX, TEMMOY )

      IF( NDIM .EQ. 3 ) THEN

C        FLUX DE LA VITESSE A TRAVERS LES FACES DES SURFACES DE L'OBJET
C        --------------------------------------------------------------
         CALL FLUXVITSF( NBNOVI, MCN(MNXYZN+WYZNOE),
     %                   NUTYEL, NBNOEF, NBELEM, MCN(MNELE+WUNDEL),
     %                   NDDLNO, NBSFEF, MCN(MNSFEF),
     %                   VXYZPN, NUMIOB(3), NUMAOB(3),
     %                   FLUNEG, FLUPOS, MCN(MNFVSF) )
      ELSE

C        FLUX DE LA VITESSE A TRAVERS LES ARETES DES LIGNES DE L'OBJET
C        -------------------------------------------------------------
         CALL FLUXVITLG( NBNOVI, MCN(MNXYZN+WYZNOE),
     %                   NUTYEL, NBNOEF, NBELEM, MCN(MNELE+WUNDEL),
     %                   NDDLNO, NBLAEF, MCN(MNLAEF),
     %                   VXYZPN, NUMIOB(2), NUMAOB(2),
     %                   FLUNEG, FLUPOS, MCN(MNFVSF) )
      ENDIF

C     STOCKAGE DES RESULTATS DU PAS DE TEMPS
      VPMIMX( 1,NBPAST) = TEMPS
      VPMIMX( 2,NBPAST) = VITMOY
      VPMIMX( 3,NBPAST) = VITMAX
      VPMIMX( 4,NBPAST) = PREMOY
      VPMIMX( 5,NBPAST) = PREMAX - PREMIN
      VPMIMX( 6,NBPAST) = FLUNEG
      VPMIMX( 7,NBPAST) = FLUPOS
      VPMIMX( 8,NBPAST) = NBPAST
      VPMIMX( 9,NBPAST) = TEMMOY
      VPMIMX(10,NBPAST) = TEMMIN
      VPMIMX(11,NBPAST) = TEMMAX


C     STOCKAGE OU PAS DU VECTEUR VITESSE+PRESSION+TEMPERATURE?
C     --------------------------------------------------------
      IF( TEMPS .GE. TSTOC*0.9999 ) THEN

C        LE NOMBRE DE VECTEURS VITESSE+PRESSION+TEMPERATURE STOCKES
         NOVVIPR = NOVVIPR + 1

C        STOCKAGE DE LA VITESSE+PRESSION+TEMPERATURE NOVVIPR
C        A CET INSTANT TEMPS SUR UN FICHIER DU REPERTOIRE DU PROJET
         CALL ECFIVIPRTE( KNOMOB, TEMPS,  NOVVIPR, NAVSTO, NBPAST,
     %                    NDIM,   NBNOVI, NBSOM,  NTDLTE,
     %                    NTDLVP, VXYZPN, NTDLTE, DMCN(MNTHETD),
     %                    KNOMFIC, IERR )

C        LE TEMPS DU STOCKAGE
         RMCN(MNTIMES-1+NOVVIPR) = TEMPS

         NBK = NUDCNB( NMPROJ )
         IF( LANGAG .EQ. 0 ) THEN
            PRINT 10150, TEMPS, NTDLVP+NTDLTE, NOVVIPR, NMPROJ(1:NBK)
         ELSE
            PRINT 20150, TEMPS, NTDLVP+NTDLTE, NOVVIPR, NMPROJ(1:NBK)
         ENDIF

10150 FORMAT('nst3th: Au TEMPS',G14.6,' STOCKAGE de',I11,
     %       ' D.L. VITESSE+PRESSION+TEMPERATURE sur FICHIER',I5,
     %       ' du REPERTOIRE du PROJET ',A)
20150 FORMAT('nst3th: At TIME',G14.6,' STORAGE of',I11,
     %       ' VELOCITY+PRESSURE+TEMPERATURE DoF on FILE',I5,
     %       ' of PROJECT DIRECTORY ',A)

C        LE PROCHAIN TEMPS DE STOCKAGE
         TSTOC = TSTOC + DTSTOC

      ENDIF

C     TRACE DES VITESSES EN TOUS LES NOEUDS AVEC SEUIL MIN A TRACER
C     =============================================================
      IF( INTERA .GE. 1 .AND. VITMAX .GT. 0D0 ) THEN

ccc      PRINT *,'trvite Temps=',TEMPS,'  ---------->'
C        MODE GRAPHIQUE AVEC X11: LA MEMOIRE PIXELS EST EFFACEE
         CALL EFFACEMEMPX

cccC        TEMPS UTILISATEUR EN SECONDES D'ATTENTE apres LE TRACE 
cccC        POUR VOIR CE TRACE si UN AUTRE SUIT IMMEDIATEMENT
ccc         TEMP2TRAC = 1.D0

C        TRACE DES FLECHES VITESSES
         CMFLEC = 2.5
         CMVITE = REAL( CMFLEC / VITMAX )

         IF( NDIM .EQ. 2 ) THEN

C           LA VISEE EN 2D:
C           INTERVALLE DE TRACE DES FLECHES VITESSES
            SEMINVIT = 0
            SEMAXVIT = VITMAX
            CALL TRVITE2D( NBNOVI, MCN(MNXYZN+WYZNOE),
     %                     MOARET, MXARET, MCN(MNLARE),
     %                     SEMINVIT, SEMAXVIT, CMVITE, 1,
     %                     VXVYVZPR(N1X,NVP1), VXVYVZPR(N1Y,NVP1),
     %                     VITMAX )

         ELSE

C           LA VISEE EN 3D:
C           INTERVALLE DE TRACE DES FLECHES VITESSES
            SEMINVIT = MIN( 2D0*VITMOY, 0.5D0*VITMAX )
            SEMAXVIT = VITMAX
            CALL TRVITE3D( NBNOVI, MCN(MNXYZN+WYZNOE),
     %                     MOARFR, MXARFR, L1ARFR, MCN(MNAFOB+WAREFR),
     %                     SEMINVIT, SEMAXVIT, CMVITE, 1,
     %                     VXVYVZPR(N1X,NVP1), VXVYVZPR(N1Y,NVP1),
     %                     VXVYVZPR(N1Z,NVP1), VITMAX )

         ENDIF

C        LE TRACE DU TITRE FINAL
         CALL LEGVIT6( KNOMOB, 1, VITMOY, VITMIN, VITMAX, CMVITE )

cccC        TEMPS UTILISATEUR EN SECONDES EN ATTENTE AVANT DE CONTINUER
ccc         CALL ATTENDSEC( 1D0 )

      ENDIF

      IF( NORESO .EQ. 2 ) THEN
         NBITERGC = KGCTEMP + KGCUET(1) + KGCUET(2) + KGCUET(3)
     %            + KGCPRE  + KGCUMU(1) + KGCUMU(2) + KGCUMU(3)
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'Au TEMPS',TEMPS,' Nb ITERATIONS des GC=',NBITERGC
         ELSE
            PRINT*,'At TIME',TEMPS,' CG ITERATIONS NUMBER=',NBITERGC
         ENDIF
      ENDIF

C     TEMPS CPU donne par le systeme a la fin de l'iteration n+1
      call cpu_time(t_cpu_it1)
C     LA DATE EN SECONDES DEPUIS LE 1/1/70 MINUIT
      CALL  SECONDES1970( DATE )
      SECONDES = DATE - DATE0
      TIMEMOY1DT = TIMEMOY1DT + SECONDES
      IF( LANGAG .EQ. 0 ) THEN
         PRINT 10105, t_cpu_it1-t_cpu_it0, SECONDES
10105    FORMAT(' nst3th: Pour ce pas de temps: ',g12.3,
     %' secondes CPU en ',g12.3,' secondes reelles')
      ELSE
         PRINT 20105, t_cpu_it1-t_cpu_it0, SECONDES
20105    FORMAT(' nst3th: For this time step: ',g12.3,
     %' CPU seconds in ',g12.3,' real seconds')
      ENDIF
      t_cpu_it0 = t_cpu_it1
      DATE0     = DATE


C     PASSAGE AU PAS DE TEMPS SUIVANT tn+1 des TABLEAUX
C     *************************************************
C     PERMUTATION DE V3P0 ET V3P1
C     ===========================
      N    = NVP0
      NVP0 = NVP1
      NVP1 = N

C     PERMUTATION DE TEMPERtn et TEMPERtn+1
C     =====================================
      N        = MNTHETn
      MNTHETn  = MNTHET
      MNTHET   = N
      MNTHETD  = ( MNTHET  + 1 ) / MOREE2
      MNTHETnD = ( MNTHETn + 1 ) / MOREE2

      VITMOYn = VITMOY
      VITMAXn = VITMAX

      IF( TEMPS + DT .LT. TPSFIN*1.0001 ) GOTO 100
C    ##############################################################
C    ##                                                          ##
C    ##                FIN DE LA BOUCLE EN TEMPS                 ##
C    ##                                                          ##
C    ##############################################################


C     STOCKAGE OU PAS DU VECTEUR VITESSE+PRESSION A CE TEMPS?
C     -------------------------------------------------------
      IF( RMCN(MNTIMES-1+NOVVIPR) .NE. TEMPS ) THEN

C        LE NOMBRE DE VECTEURS VITESSEPRESSION STOCKES
         NOVVIPR = NOVVIPR + 1

C        STOCKAGE DE LA VITESSEPRESSION NOVVIPR A CET INSTANT TEMPS
C        SUR FICHIER DU REPERTOIRE DU PROJET
         CALL ECFIVIPRTE( KNOMOB, TEMPS,  NOVVIPR, NAVSTO, NBPAST,
     %                    NDIM,   NBNOVI, NBSOM,   NTDLTE,
     %                    NTDLVP, VXYZPN, NTDLTE,  DMCN(MNTHETnD),
     %                    KNOMFIC, IERR )

C        LE TEMPS DU STOCKAGE
         RMCN(MNTIMES-1+NOVVIPR) = TEMPS

         IF( LANGAG .EQ. 0 ) THEN
            PRINT 10150, TEMPS, NTDLVP, NOVVIPR, NMPROJ
         ELSE
            PRINT 20150, TEMPS, NTDLVP, NOVVIPR, NMPROJ
         ENDIF

      ENDIF


C     BILAN SUR LES REMONTEES DES CARACTERISTIQUES TROP LONGUES
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*, NBCHTROL,
     %    ' REMONTEES TROP LONGUES de la CARACTERISTIQUE'
         IF( NBCHTROL .GT. 0 ) THEN
            PRINT*, '=> DIMINUER le PAS de TEMPS'
         ENDIF
      ELSE
         PRINT*, NBCHTROL,
     %      ' TOO LONG COMES BACK of CHARACTERISTICS'
         IF( NBCHTROL .GT. 0 ) THEN
            PRINT*, '=> REDUCE the TIME STEP'
         ENDIF
      ENDIF
      GOTO 9999


C     PLACE MEMOIRE INSUFFISANTE
 9993 NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'ERREUR: PLACE MEMOIRE INSUFFISANTE'
         KERR(2) = 'REDUIRE LE MAILLAGE'
      ELSE
         KERR(1) = 'ERROR: NOT ENOUGH LARGE MEMORY'
         KERR(2) = 'REDUCE the MESH'
      ENDIF
      CALL LEREUR
      IERR = 20


C     DESTRUCTION DES TMC DEVENUS INUTILES
C     ====================================
C     FREE MEMORY USED BY ARRAY VG and PG is DEALLOCATED
 9999 IF( IERVITC .EQ. 0 ) then
         PRINT*,'nst3th: VITC est DESALLOUE'
         DEALLOCATE( VITC )
      ENDIF

      IF( IERVXVYVZPR  .EQ. 0 ) then
         PRINT*,'nst3th: VXVYVZPR est DESALLOUE'
         DEALLOCATE( VXVYVZPR )
      ENDIF

      IF( IERPr1Pr  .EQ. 0 ) then
         PRINT*,'nst3th: Pr1Pr est DESALLOUE'
         DEALLOCATE( Pr1Pr )
      ENDIF

      IF( IERPVAUX  .EQ. 0 ) then
         PRINT*,'nst3th: PVAUX est DESALLOUE'
         DEALLOCATE( PVAUX )
      ENDIF

      IF( IERLP1COLO  .EQ. 0 ) then
         PRINT*,'nst3th: LP1COLO est DESALLOUE'
         DEALLOCATE( LP1COLO )
      ENDIF

      IF( IERLP1LIGN  .EQ. 0 ) then
         PRINT*,'nst3th: LP1LIGN est DESALLOUE'
         DEALLOCATE( LP1LIGN )
      ENDIF

      IF( IERDAUX3  .EQ. 0 ) then
         PRINT*,'nst3th: DAUX3 est DESALLOUE'
         DEALLOCATE( DAUX3 )
      ENDIF

      IF( IERDAUX2  .EQ. 0 ) then
         PRINT*,'nst3th: DAUX2 est DESALLOUE'
         DEALLOCATE( DAUX2 )
      ENDIF

      IF( IERDAUX1  .EQ. 0 ) then
         PRINT*,'nst3th: DAUX1 est DESALLOUE'
         DEALLOCATE( DAUX1 )
      ENDIF

      IF( IERV3P2  .EQ. 0 ) then
         PRINT*,'nst3th: V3P2 est DESALLOUE'
         DEALLOCATE( V3P2 )
      ENDIF

      IF( IERNODLTEFX  .EQ. 0 ) then
         PRINT*,'nst3th: NODLTEFX est DESALLOUE'
         DEALLOCATE( NODLTEFX )
      ENDIF

      IF( IERNODLVIFX  .EQ. 0 ) then
         PRINT*,'nst3th: NODLVIFX est DESALLOUE'
         DEALLOCATE( NODLVIFX )
      ENDIF

      IF( IERNODLPRFX  .EQ. 0 ) then
         PRINT*,'nst3th: NODLPRFX est DESALLOUE'
         DEALLOCATE( NODLPRFX )
      ENDIF

      IF( IERVGALLOC  .EQ. 0 ) then
         PRINT*,'nst3th: VG est DESALLOUEE'
         DEALLOCATE( VG )
      ENDIF
      IF( IERPGALLOC  .EQ. 0 ) then
         PRINT*,'nst3th: PG est DESALLOUEE'
         DEALLOCATE( PG )
      ENDIF
      
      IF( MNTHET .NE. 0 ) CALL TNMCDS( 'REEL2',  NTDLTE, MNTHET )
      IF( MNTHETn.NE. 0 ) CALL TNMCDS( 'REEL2',  NTDLTE, MNTHETn )
ccc      IF( MNTHET0.NE. 0 ) CALL TNMCDS( 'REEL2',  NTDLTE, MNTHET0 )
      IF( MNBG   .GT. 0 ) CALL TNMCDS( 'REEL2',  NTDLTE, MNBG )

      IF( NBVVEF .GT. 0 ) CALL TNMCDS( 'ENTIER', NBELEM*NBVVEF, MNVVEF )
      IF( NBSFEF .GT. 0 ) CALL TNMCDS( 'ENTIER', NBELEM*NBSFEF, MNSFEF )
      IF( NBLAEF .GT. 0 ) CALL TNMCDS( 'ENTIER', NBELEM*NBLAEF, MNLAEF )
      IF( NBPSEF .GT. 0 ) CALL TNMCDS( 'ENTIER', NBELEM*NBPSEF, MNPSEF )

      IF( NBTEFX .GT. 0 ) CALL TNMCDS( 'ENTIER', NBTEFX, MNNTEFX )
      IF( NBTEFX .GT. 0 ) CALL TNMCDS( 'REEL2',  NBTEFX, MNVTEFX )
      IF( NBPRFX .GT. 0 ) CALL TNMCDS( 'ENTIER', NBPRFX, MNNPRFX )
      IF( NBPRFX .GT. 0 ) CALL TNMCDS( 'REEL2',  NBPRFX, MNVPRFX )
      IF( NBVIFX .GT. 0 ) CALL TNMCDS( 'ENTIER', NBVIFX, MNNVIFX )
      IF( NBVIFX .GT. 0 ) CALL TNMCDS( 'REEL2',  NBVIFX, MNVVIFX )

      IF( MNFVSF .GT. 0 ) CALL TNMCDS( 'REEL2',  MOSFOB, MNFVSF  )

C     LES TABLEAUX AUXILIAIRES POUR LES ELEMENTS FINIS P2
      IF( MNTAUX .GT. 0 ) CALL TNMCDS( 'REEL2',  MOAUX,    MNTAUX )
      IF( MNTAEL .GT. 0 ) CALL TNMCDS( 'REEL2',  MOTAEL,   MNTAEL )
      IF( MNTHER .GT. 0 ) CALL TNMCDS( 'REEL2',  128,      MNTHER )
      IF( MNNODL .GT. 0 ) CALL TNMCDS( 'ENTIER', MONODL,   MNNODL )
      IF( MNX    .GT. 0 ) CALL TNMCDS( 'REEL',   NBDL1EF*3,MNX )

C     BILAN SUR LA PLACE MEMOIRE CENTRALE OCCUPEE
C     ===========================================
      DMOTSTO = 901.D0+NBNOVI*2+NBTEFX+NBPRFX+NBVIFX+NBSOM*3
     %        + MOREE2*( NBSOM+NBNOVI+NTDLVI+NTDLVP*3+NTDLTE*3
     %           +NBTEFX+NBPRFX+NBVIFX+ MOSFOB+NBCMVG+NBCMPG+NBCMTG )
     %        + NOVVIPR + MOREE2 *(NTDLVP*NOVVIPR+MOAUX)

      IF( NORESO .GE. 2 ) THEN
C        GRADIENT CONJUGUE STOCKAGE MORSE
         DMOTSTO = DMOTSTO + NBNOVI+1 +NBCMPG+NBCMVG+NBCMTG
      ENDIF
      IF( NORESO .EQ. 1 ) THEN

C        FACTORISATION COMPLETE CROUT STOCKAGE PROFIL
         IF( LANGAG .EQ. 0 ) THEN
            PRINT 19002, DMOTSTO
         ELSE
            PRINT 29002, DMOTSTO
         ENDIF
19002 FORMAT(/'STOCKAGE PROFIL des MATRICES et VECTEURS=',G25.17,
     %' MOTS MEMOIRE'/)
29002 FORMAT(/'SKYLINE STORAGE of MATRICES and VECTORS=', G25.17,
     %' MEMORY WORDS'/)

      ENDIF

C     TOTAL DES MOTS NECESSAIRES A LA RESOLUTION DU PROBLEME
C     ------------------------------------------------------
      IF( LANGAG .EQ. 0 ) THEN
         PRINT 19005, DMOTSTO
      ELSE
         PRINT 29005, DMOTSTO
      ENDIF
19005 FORMAT(' STOCKAGE TOTAL des MATRICES et VECTEURS=',G25.17,
     %' MOTS MEMOIRE')
29005 FORMAT(' TOTAL STORAGE of MATRICES and VECTORS=', G25.17,
     %' MEMORY WORDS')

C     CHAINAGE DES FACES FRONTALIERES C-A-D APPARTENANT A UN SEUL EF
C     A PARTIR DU ZERO EN POSITION 7 DE LFACES
C     --------------------------------------------------------------
      IF( MNLFAC .GT. 0 ) THEN
         CALL MACFAF( MOFACE, MXFACE, MCN(MNLFAC),  L1FAFR, NBFAFR )
C        MISE A JOUR DU NUMERO DE LA PREMIERE FACE FRONTALIERE (DANS UN SEUL EF)
         MCN( MNFAOB + W1FAFR ) = L1FAFR
      ENDIF

C     REMISE EN ETAT DE LA FONTE COURANTE
      IF( INTERA .GT. 0 ) CALL CHARGEFONTE( NOFONT0 )

      IF( IERNONOSO  .EQ. 0 ) then
         PRINT*,'nst3th: NONOSO est DESALLOUE'
         DEALLOCATE( NONOSO)
      ENDIF

      IF( IERNO1EFN  .EQ. 0 ) then
         PRINT*,'nst3th: NO1EFN est DESALLOUE'
         DEALLOCATE( NO1EFN)
      ENDIF

C     TEMPS CALCUL DES VITESSES PRESSIONS
C     ===================================
      DVITPR = DINFO( 'DELTA CPU' )

C     TEMPS REEL D'EXECUTION DE nst3th
      CALL SECONDES1970( DATE )
      SECONDES = DATE - DATE00
      IF( LANGAG .EQ. 0 ) THEN
         PRINT 19000, SECONDES
19000    FORMAT(' nst3th: Secondes REELLES execution=',F10.0)
      ELSE
         PRINT 29000, SECONDES
29000    FORMAT(' nst3th: Elapsed REAL seconds of execution=',F10.0)
      ENDIF

C     AUTRE CALCUL DU TEMPS DES THREADS
      call system_clock( count=t1, count_rate=nbclockps )
      S = real( t1-t0, kind=8 ) / real( nbclockps, kind=8 )
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'nst3th: Nombre des secondes d''EXECUTION=',S,
     %          ' sur',nbthreads,' threads'
      ELSE
         PRINT*,'nst3th: Thread EXECUTION Seconds=',S,
     %          ' on',nbthreads,' threads'
      ENDIF

C$OMP PARALLEL SHARED( tt0 )
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'nst3th: Numero Thread=',OMP_GET_THREAD_NUM(),
     %                  ' OMP seconds=', OMP_GET_WTIME()-tt0
      ELSE
         PRINT*,'nst3th: Thread Number=',OMP_GET_THREAD_NUM(),
     %                  ' OMP seconds=', OMP_GET_WTIME()-tt0
      ENDIF
C$OMP END PARALLEL

      call cpu_time(t_cpu_1)
      PRINT*,'nst3th: CPU_time=',t_cpu_1-t_cpu_0

      IF( NBPAST .GT. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'nst3th: Temps moyen EXECUTION 1 pas de temps=',
     %              TIMEMOY1DT/NBPAST,' pour',NBPAST,' pas de temps'
         ELSE
          PRINT*,'nst3th: Mean EXECUTION Time of 1 time step=',
     %            TIMEMOY1DT/NBPAST,' for',NBPAST,' time steps'
         ENDIF
      ENDIF

C     TEMPS UTILISATEUR EN ATTENTE EN SECONDES apres UN TRACE
      TEMP2TRAC = 0.D0

      RETURN
      END
