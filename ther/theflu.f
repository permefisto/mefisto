      SUBROUTINE THEFLU( KNOMOB, NTLXOB, MNTOPO, NOAXIS, D2PI,
     &                   NDIM,   MOREE2, NDSM,   NTDL,
     &                   NBTYEL, MNNPEF, NDPGST, MNTPOB,
     &                   MNTAUX, MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     &                   MOTAEL, MNTAEL, MNX,    MNVECT )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL LES FLUX AUX POINTS D'INTEGRATION NUMERIQUE DES FACES
C -----    DES ELEMENTS FINIS LAGRANGE ISOPARAMETRIQUES 2D 3D AXISYMETRIQUES
C          CONSTRUCTION DES TMS FLUXPT"NMTYEL DE L'OBJET
C
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET DE FLUX A CALCULER
C NTLXOB : NUMERO DU TMS DU LEXIQUE DE L'OBJET
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE DE L'OBJET
C NOAXIS : 1 SI PROBLEME AXI-SYMETRIQUE (=> NDIM=2 R=x et Z=y et 0=z)
C          0 SINON
C D2PI   : 2 PI REEL DOUBLE PRECISION
C NDIM   : ESPACE 1 OU 2 OU 3 DE L'OBJET (2 EN AXISYMETRIQUE)
C MOREE2 : 2 SI UN DOUBLE PRECISION OCCUPE 2 MOTS, 1 SINON
C NDSM   : NOMBRE DE CAS DE CHARGE OU VECTEUR"TEMPERATURE OU
C          VECTEUR"VALEURPROPRE
C          CE NOMBRE PEUT ETRE >1 PAR EXEMPLE EN INSTATIONNAIRE
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE EN TEMPERATURE
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX NPEF"TYPE EF
C NDPGST : CODE DE STOCKAGE DES SOMMETS POINTS NOEUDS
C MNTPOB : ADRESSE MCN DES TABLEAUX POLYNOMES DE BASE DES TYPES D'EF
C MNTAUX : ADRESSE MCN DES TABLEAUX AUXILIAIRES
C MNXYZP : ADRESSE MCN DE TMS XYZSOMMET DE L'OBJET
C NUMIOB : NUMERO MINIMAL DES OBJETS
C NUMAOB : NUMERO MAXIMAL DES OBJETS
C MNDOEL : ADRESSE MCN DES DONNEES DE L'OBJET
C MNTHER : ADRESSE MCN DU TABLEAU TENSEUR DE CONDUCTIVITE,...
C MOTAEL : NOMBRE DE REEL2 DE LA DECLARATION DU TABLEAU TAEL
C          LA MATRICE ELEMENTAIRE ET LES NDSM SECONDS MEMBRES
C          MOTAEL = NBDLMX * (NBDLMX+1) / 2 + NBDLMX * NDSM
C MNTAEL : ADRESSE MCN DES TABLEAUX ELEMENTAIRES POUR LES CONTRAINTES
C MNX    : ADRESSE MCN DES COORDONNEES DES POINTS DE L'EF COURANT
C MNVECT : ADRESSE MCN DU TMS VECTEUR"TEMPERATURE OU VECTEUR"VALEURPROPRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1998
C MODIFS: ALAIN PERRONNET  TIMS NTU TAIPEI TAIWAN           Octobre 2009
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ponoel.inc"
      include"./incl/gsmenu.inc"
      include"./incl/donthe.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___arete.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___fluxpt.inc"
      include"./incl/a___fluxfr.inc"
      include"./incl/a___dtemperature.inc"
      include"./incl/cthet.inc"
      include"./incl/ctemps.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
C
      INTEGER           NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NUMIOB(4), NUMAOB(4),  MNDOEL(4)
      CHARACTER*4       CHARX, NOMELE(2)
      CHARACTER*(*)     KNOMOB
      CHARACTER*160     KNOM
      DOUBLE PRECISION  DELTA, D2PI
C
C***************************************************************************
C***************************************************************************
C     CALCUL DES FLUX SUR CHAQUE FRONTIERE DE TYPE MAX DE L'OBJET
C    (1D => 2 SOMMETS, 2D=> 3 ou 4 ARETES, 3D=>4 ou 5 ou 6 FACES)
C***************************************************************************
C***************************************************************************
C
      IF( NDSM .LE. 0 ) RETURN
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10545)
      ELSE
         WRITE(IMPRIM,20545)
      ENDIF
10545 FORMAT(/' Les FLUX NORMAUX de CHALEUR aux FRONTIERES de l''OBJET'
     %,/1X,80(1H=))
20545 FORMAT(/' The NORMAL HEAT FLUXES on the OBJECT BOUNDARIES',
     %/1X,80(1H=))
C
      MNFLTO = 0
      MNFLPT = 0
      MNFMIX = 0
      MNPOLQ = 0
      MNPOLA = 0
      MNPOL  = 0
      MNPOIQ = 0
      MNPOID = 0
      MNPOIA = 0
      MNPDEL = 0
      MNF2   = 0
      MNDPOQ = 0
      MNDPOL = 0
      MNDPOA = 0
      MNDPA1 = 0
      MNDPA2 = 0
      MNDFM3 = 0
      MNDFM2 = 0
      MNDDPO = 0
      MNDDP  = 0
C
C     LE NOMBRE D'OBJETS INTERNES ET FRONTALIERS
      NBOBIN = MCN( MNTOPO + WBOBIN )
      NBOBCL = MCN( MNTOPO + WBOBCL )
C     L'ADRESSE MCN DU DEBUT DU TABLEAU NUOBCL
      MNOBIN = MNTOPO + WMTYEL + NBTYEL
      MNOBCL = MNOBIN + MOTVAR(13) * NBOBIN
C
C     MNTEMP : ADRESSE MCN DU 1-ER MOT DU PREMIER VECTEUR DU
C              VECTEUR"TEMPERATURE OU VECTEUR"VALEURPROPRE
      MNTEMP = MNVECT + WECTEU
      MNTHET = MNTEMP
C
C     LES TEMPS ONT ILS ETE STOCKES?
      IF( MCN( MNVECT + WBCPIN ) .GT. 0 ) THEN
C        OUI: LE TEMPS INITIAL EST CELUI DU VECTEUR TEMPERATURE
         NBTEMP = MCN( MNVECT + WBVECT )
         NTDL   = MCN( MNVECT + WBCOVE )
         MNTIME = MNVECT + WECTEU + NTDL * NBTEMP * MOREE2 -1
         TEMPS  = RMCN( MNTIME + NBTEMP )
      ENDIF
C
C     TYPE DES OBJETS SUR LA FRONTIERE
C    ( POINT=>1 EN 1D, LIGNE=>2 EN 2D, SURFACE=>3 EN 3D )
      NTYF = NDIM
C
C     LE TABLEAU DES FLUXPT PAR POINT D'UN EF  (MAX NBPOLY * NPIA )
      MOFLPT = 270
C     NBPOLY<=27 ET NPIA<=9
      CALL TNMCDC( 'REEL2', MOFLPT, MNFLPT )
C
C     LE TABLEAU DES FLUXFR PAR OBJET FRONTIERE ET POUR CHAQUE TEMPS CACULE
      NBPLSF = ( NUMAOB(NTYF)-NUMIOB(NTYF)+1 )
      MOFLTO = NDSM * NBPLSF
      CALL TNMCDC( 'REEL2', MOFLTO, MNFLTO )
      CALL AZEROD( MOFLTO, MCN(MNFLTO) )
C
C     LE TABLEAU ELEMENTAIRE DES FLUX SUR LES FRONTIERES DE L'EF COURANT
      I = NDSM * 4
      IF( I .GT. MOTAEL ) THEN
C        LE TABLEAU TAEL EST TROP PETIT
C        IL EST DETRUIT ET REDECLARE
         CALL TNMCDS( 'REEL2', MOTAEL, MNTAEL )
         MOTAEL = I
         CALL TNMCDC( 'REEL2', MOTAEL, MNTAEL )
      ENDIF
C
      IF( NDIM .EQ. 2 ) THEN
C
C        CALCUL PAR HACHAGE DES ARETES DE L'OBJET A PARTIR DE TOPO+NPEF"...
C        NECESSAIRE AU CALCUL DU SAUT DE LA DERIVEE NORMALE SUR CHAQUE ARETE
         CALL HAC2AF( KNOMOB, 3, NTARSU, MNARSU, IERR )
         IF( IERR .NE. 0 ) GOTO 9999
C
C        LE NOMBRE D'ENTIERS PAR ARETE
         MOARET = MCN( MNARSU + WOARET )
C        LA MAJORATION DU NOMBRE D'ARETES
         MXARET = MCN( MNARSU + WXARET )
C        LE NOMBRE D'ARETES FRONTALIERES
         NBARFB = MCN( MNARSU + WBARFB )
C        LE NOMBRE D'ARETES INTERFACES
         NBARIN = MCN( MNARSU + WBARIN )
C        LE NUMERO MINIMAL DE LIGNE DE L'OBJET
         NUMILF = MCN( MNARSU + WUMILF )
C        LE NUMERO MAXIMAL DE LIGNE DE L'OBJET
         NUMXLF = MCN( MNARSU + WUMXLF )
C        LE NUMERO DE LA PREMIERE ARETE FRONTALIERE
         L1ARFB = MCN( MNARSU + W1ARFB )
C        LE NUMERO DE LA PREMIERE ARETE INTERFACE
         L1ARIN = MCN( MNARSU + W1ARIN )
C        ADRESSE MCN DU 1-ER MOT DU TABLEAU LARETE
         MNARET = MNARSU + W1LGFR + NUMXLF - NUMILF + 1
C
      ENDIF
C
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
C     ----------------------------------------
      NUTTEF = 0
      DO 200 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"TYPE_EF
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES ELEMENTS
         MNNDEL = MNELE + WUNDEL
         MNPGEL = MNNDEL
         IF( NDPGST .GE. 2 ) THEN
            MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
         ENDIF
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL, NOMELE )
         CALL ELTYCA( NUTYEL )
C
C        ----------------------------------------------
C        LES ADRESSAGES POUR RECUPERER LES INFORMATIONS
C        EN FONCTION DU TYPE DE L'ELEMENT FINI
C        ----------------------------------------------
         GOTO(  51,  51,  51,  51,  50,  50,  50,  50,  50,  50,
     &          50,  50,  63,  50,  51,  51,  50,  51,  63,  51,
     &          51,  51,  51,  51,  50,  50,  50,  51,  51,  50,
     &          51,  51,  51,  50 ) ,NUTYEL
C
C        ERREUR
C        ------
   50    NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
           KERR(1)='ERREUR: TYPE d''ELEMENT FINI '//NOMELE(1)//NOMELE(2)
           KERR(2)='NON PROGRAMME'
         ELSE
           KERR(1) ='ERROR: FINITE ELEMENT TYPE '//NOMELE(1)//NOMELE(2)
           KERR(2) ='NOT PROGRAMMED'
         ENDIF
         CALL LEREUR
         GOTO 9999
C
C        ELEMENTS AXISYMETRIQUES ET NON AXISYMETRIQUES 2D LAGRANGE
C        ELEMENTS 3D LAGRANGE ISOPARAMETRIQUES
C        *********************************************************
C        L'ELEMENT FINI : 2D ARETE DE REFERENCE, 3D SURFACE DE REFERENCE
   51    L      = MNTPOB + (NOTYEL-1) * MXPOBA
         IF( NDIM .EQ. 1 ) GOTO 55
C
         IA     = MCN( L  )
C        DIMENSION DE L ESPACE
C        NDIMA  = MCN( IA )
C        NOMBRE DE POLYNOMES D INTERPOLATION
         NBPOLA = MCN( IA + 1 )
C        NOMBRE DE POINTS D INTEGRATION NUMERIQUE
         NPIA   = MCN( IA + 2 )
C        ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
         MNPOIA = IA + 8
C        ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
         MNPOLA = MNPOIA + MCN(IA + 3) * MOREE2 * NPIA
C        ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D I
         MNDPOA = MNPOLA + MCN(IA + 4) * MOREE2 * NBPOLA * NPIA
C
         IF( NBNSOM .EQ. 5 .OR. NBNSOM .EQ. 6 ) THEN
C           EF AVEC 2 TYPES DE FACES : PAR EXEMPLE PYRAMIDE et PENTAEDRE
C                                      POLA => TRIANGLE
C                                      POLQ => QUADRANGLE
C
C           LES POLYNOMES DE L'EF DE DIMENSION - 1 POUR LE SECOND TYPE DE
C           FACE=QUADRANGLE POSITIONNE EN 4-EME POSITION (CF SP ETTAEL)
C           ...........................................................
            IA     = MCN( L + 3  )
C           DIMENSION DE L ESPACE
C           NDIMQ  = MCN( IA )
C           NOMBRE DE POLYNOMES D INTERPOLATION
            NBPOLQ = MCN( IA + 1 )
C           NOMBRE DE POINTS D INTEGRATION NUMERIQUE
            NPIQ   = MCN( IA + 2 )
C           ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
            MNPOIQ = IA + 8
C           ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
            MNPOLQ = MNPOIQ + MCN(IA + 3) * MOREE2 * NPIQ
C           ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D'INTEGRAT
            MNDPOQ = MNPOLQ + MCN(IA + 4) * MOREE2 * NBPOLQ * NPIQ
C
         ELSE
C
C           POLYNOMES DES FACES DE TYPE 2 = CEUX DE LA FACE DE TYPE 1
C           NDIMQ  = MCN( IA )
C           NOMBRE DE POLYNOMES D INTERPOLATION
            NBPOLQ = NBPOLA
C           NOMBRE DE POINTS D INTEGRATION NUMERIQUE
            NPIQ   = NPIA
C           ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
            MNPOIQ = MNPOIA
C           ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
            MNPOLQ = MNPOLA
C           ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D'INTEGRAT
            MNDPOQ = MNDPOA
C
         ENDIF
C
C        LES VALEURS DES POLYNOMES DE L'EF DE DIMENSION TOTALE
C        EN 2D SURFACE DE REFERENCE, EN 3D VOLUME DE REFERENCE
C        ......................................................
C        RECHERCHE DU TABLEAU DE POBA ET PARTAGE EN P ET DP
         L      = L + 1
 55      IA     = MCN( L )
C        DIMENSION DE L ESPACE
C        NDIM   = MCN( IA )
C        NOMBRE DE POLYNOMES D INTERPOLATION
         NBPOLY = MCN( IA + 1 )
C        NOMBRE DE POINTS D INTEGRATION NUMERIQUE
         NPI    = MCN( IA + 2 )
C        ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
         MNPOID = IA + 8
C        ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
         MNPOL  = MNPOID + MCN(IA + 3) * MOREE2 * NPI
C        ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS
C        D'INTEGRATION DE L'ELEMENT FINI
         MNDPOL = MNPOL  + MCN(IA + 4) * MOREE2 * NBPOLY * NPI
C        ADRESSE DES VALEURS DES DERIVEES SECONDES DES POLYNOMES AUX POINTS
C        D'INTEGRATION DE L'ELEMENT FINI
         MNDDPO = MNDPOL + MCN(IA + 5) * MOREE2 * NDIM * NBPOLY * NPI
C
C        LES DERIVEES DES POLYNOMES AUX POINTS D'INTEGRATION
C        DES ARETES DE L'ELEMENT DE REFERENCE  DPAR
         L = L + 1
C        LE DECALAGE POUR ATTEINDRE CHAQUE ARETE EN 2D, CHAQUE FACE EN 3D
         MNDPA1 = MCN(L)
C        LE DEBUT DU TABLEAU DPAR
         MNDPA2 = MNDPA1 + 6
C
C        LES TABLEAUX AUXILIAIRES
         NDIM2 = NDIM * ( NDIM + 1 ) / 2
         MNF1   = MNTAUX
         MNF2   = MNF1   + MOREE2 * NPI
         MNF3   = MNF2   + MOREE2 * NPI
         MNPDEL = MNF1   + MOREE2 * NPI * NDIM
         MNDP   = MNPDEL + MOREE2 * NPI
         MNDFM1 = MNDP   + MOREE2 * NPI * NDIM  * NBPOLY
C
         MNDFM2 = MNDFM1 + MOREE2 * NPI * NDIM  * NDIM
         MNDFM3 = MNDFM2 + MOREE2 * NPI * NDIM2 * NDIM
         MNDDP  = MNDFM3 + MOREE2 * NPI * NDIM2 * NDIM2
C        ADRESSE DERNIER MOT = MNDDP+ MOREE2 * NPI * NDIM2 * NBPOLY
C        MOTAUX = MOREE2 * NPI * ( NDIM + 1 + NDIM * (NDIM+NBPOLY)
C               + NDIM*(NDIM+1)/2 * ( NDIM + NDIM*(NDIM+1)/2 + NBPOLY ) )
         GOTO 90
C
C
C        TRIANGLE 2P1D  ET  TETRAEDRE 3P1D
C        =================================
   63    NBPOLY = NDIM + 1
         NPI    = 1
         NPIQ   = 1
         NPIA   = 1
         MNF1   = MNTAUX
         MNF2   = MNF1 + MOREE2
         MNF3   = MNF2 + MOREE2
         MNDP   = MNF3 + MOREE2
         MNDFM1 = MNDP + MOREE2 * NDIM * NBPOLY
C        AU TOTAL = MNDFM1 + MOREE2 * NDIM * NDIM
C
C        LA DECLARATION DU TABLEAU 'FLUXPT"NMTYEL'
C        -----------------------------------------
   90    KNOM = 'FLUXPT"' // CHARX( MCN(MNTOPO+WMTYEL+NOTYEL-1) )
         CALL LXTSOU( NTLXOB, KNOM, NTFLUX, MNFLUX )
         IF( NTFLUX .GT. 0 ) THEN
C           LE TABLEAU EST DETRUIT POUR ETRE REDECLARE
            CALL LXTSDS( NTLXOB, KNOM )
         ENDIF
C        DECLARATION OUVERTURE DU TABLEAU 'FLUX'
         IF( NDIM .EQ. 3 ) THEN
            NBFAFX = NFACE
C           NBPNFX=NOMBRE DE POINTS OU LES FLUX NORMAUX SONT CALCULES
            IF( NFACE .EQ. 5 ) THEN
               NBPNFX = 3 * NPIQ + 2 * NPIA
            ELSE
               NBPNFX = NFACE * NPIA
            ENDIF
         ELSE IF( NDIM .EQ. 2 ) THEN
            NBFAFX = NARET
            NBPNFX = NARET * NPIA
         ELSE IF( NDIM .EQ. 1 ) THEN
            NBFAFX = 2
            NBPNFX = 2
         ENDIF
         I = MOREE2 * NBPNFX * NBELEM * NDSM
         L = WLUXNP + I + 6 * NBPNFX * NBELEM
         CALL LXTNDC( NTLXOB, KNOM, 'MOTS', L )
         CALL LXTSOU( NTLXOB, KNOM, NTFLUX, MNFLUX )
         MNFLNP = MNFLUX + WLUXNP
         MNCOPN = MNFLNP + I
C
C        LA DECLARATION DU TABLEAU 'DTEMPERATURE"NMTYEL'
C        -----------------------------------------------
         KNOM = 'DTEMPERATURE"' // CHARX( MCN(MNTOPO+WMTYEL+NOTYEL-1) )
         CALL LXTSOU( NTLXOB, KNOM, NTDTEM, MNDTEM )
         IF( NTDTEM .GT. 0 ) THEN
C           LE TABLEAU EST DETRUIT POUR ETRE REDECLARE
            CALL LXTSDS( NTLXOB, KNOM )
         ENDIF
         NBM = MOREE2 * NBELEM * NPI * NDIM * NDSM
         L   = WTEMPE + NBM + NBELEM * NPI * NDIM
         CALL LXTNDC( NTLXOB, KNOM, 'MOTS', L )
         CALL LXTSOU( NTLXOB, KNOM, NTDTEM, MNDTEM )
C        LES COORDONNEES DES POINTS D'INTEGRATION
         MNCOPI = MNDTEM + WTEMPE + NBM
C
C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NUTYEL
C        ==================================================
         DO 190 NUELEM = 1, NBELEM
C           LE NUMERO GLOBAL DE CET EF
            NUTTEF = NUTTEF + 1
C
C           LES NOEUDS DE L'ELEMENT
C           -----------------------
CCC            CALL EFNOEU( MNELE, NUELEM, NBNDEL, MCN(MNNODL) )
C
C           LES POINTS GEOMETRIQUES DE L'ELEMENT
C           ------------------------------------
CCC            CALL EFPOGE( MNELE, NUELEM, NBPGEF, MCN(MNPOEF) )
C
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
C           ----------------------------------------
            CALL EFPLSV( MNELE , NUELEM,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                   NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C           LES COORDONNEES DES POINTS DE L'EF COURANT
C           ------------------------------------------
            DO  100 I=1,NBPOE
C              LE NUMERO DU I-EME POINT DE L'ELEMENT NUELEM
               N    = MCN( MNPGEL-1 + NUELEM + NBELEM*(I-1) )
               MNCE = MNXYZP + WYZPOI + (N-1) * 3
               L    = MNX - 1 + I
               RMCN( L ) = RMCN(MNCE)
               IF( NDIM .GE. 2 ) RMCN(L+NBPOE      ) = RMCN(MNCE+1)
               IF( NDIM .EQ. 3 ) RMCN(L+NBPOE+NBPOE) = RMCN(MNCE+2)
  100       CONTINUE
C
C           --------------------------------------------
C           CALCUL DES TABLEAUX AUXILIAIRES ET DES FLUX
C           ELEMENTAIRES SELON LE TYPE DE L'ELEMENT FINI
C           --------------------------------------------
            GOTO( 102, 102, 102, 102,  50,  50,  50,  50,  50,  50,
     &             50,  50, 113,  50, 102, 102,  50, 102, 119, 111,
     &            111, 111, 111, 111,  50,  50,  50, 101, 102,  50,
     &            111, 111, 101,  50),NUTYEL
C
C           ***********
C           1D LAGRANGE
C           ***********
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN 1D
 101        IF( RMCN(MNX) .GT. RMCN(MNX+1) ) THEN
C              PERMUTATION DES 2 SOMMETS POUR AVOIR UN EF DE MESURE>0
               XS = RMCN(MNX)
               RMCN(MNX) = RMCN(MNX+1)
               RMCN(MNX+1 ) = XS
C
               MNCE = MNPGEL-1 + NUELEM
               NST = MCN( MNCE )
               MCN( MNCE ) = MCN( MNCE + NBELEM )
               MCN( MNCE + NBELEM ) = NST
C
               NST       = NOOBPS(1)
               NOOBPS(1) = NOOBPS(2)
               NOOBPS(2) = NST
            ENDIF
            CALL E11LAG( NBPOLY, NPI, MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX) ,  MCN(MNF1),
     %                   MCN(MNPDEL), MCN(MNDP), DELTA )
C
C           CALCUL DU FLUX DE CHAQUE SOMMET ET DU FLUX NORMAL AUX POINTS
            CALL TF1LAG( RMCN(MNX), NDSM, NBELEM, NUELEM,
     %                   NTDL, NBPOLY, MCN(MNNDEL), NBPNFX,
     %                   NOOBPS, NUMIOB(1), NUMAOB(1),
     %                   NOOBLA, NUMIOB(2), NUMAOB(2), MCN(MNDOEL(2)),
     %                   MCN(MNDPA1), MCN(MNDPA2), MCN(MNTEMP),
     %                   MCN(MNFLPT),
     %                   MCN(MNCOPN), MCN(MNFLNP), MCN(MNFLTO) )
            GOTO 180
C
C           ********************************************
C           2D OU AXISYMETRIQUE LAGRANGE ISOPARAMETRIQUE
C           ********************************************
  102       CALL E32LAG( NBPOLY, NPI, MCN(MNPOID),
     %                  MCN(MNPOL), MCN(MNDPOL), MCN(MNDDPO),
     %                  RMCN(MNX),
     %                  MCN(MNF1) ,  MCN(MNF2),
     %                  MCN(MNPDEL), MCN(MNDP),   MCN(MNDFM1),
     %                  MCN(MNDDP),  MCN(MNDFM2), MCN(MNDFM3) )
C           CALCUL DU FLUX DE CHAQUE ARETE ET DU FLUX NORMAUX AUX ARETES
            CALL TF2LAG( D2PI, NOAXIS, RMCN(MNX),
     %                   NDSM, NBELEM, NUELEM,
     %                   NTDL, NBPOLY, MCN(MNNDEL), NBPNFX,
     %                   NBPOLA, NPIA,
     %                   MCN(MNPOIA), MCN(MNPOLA), MCN(MNDPOA),
     %                   NARET, NOOBLA, NUMIOB(2), NUMAOB(2),
     %                   NOOBSF(1), NUMIOB(3), NUMAOB(3),MCN(MNDOEL(3)),
     %                   MCN(MNDPA1), MCN(MNDPA2), MCN(MNTEMP),
     %                   MCN(MNFLPT), MCN(MNCOPN),
     %                   MCN(MNFLNP), MCN(MNFLTO) )
            GOTO 180
C
C           ***************************
C           3D LAGRANGE ISOPARAMETRIQUE
C           ***************************
  111       CALL E13LAG ( NBPOLY, NPI, MCN(MNPOID),
     %                   MCN(MNPOL), MCN(MNDPOL),
     %                   RMCN(MNX) , MCN(MNF1),
     %                   MCN(MNPDEL), MCN(MNDP), MCN(MNDFM1) )
C           CALCUL DU FLUX DE CHAQUE FACE ET DU FLUX NORMAUX AUX FACES
            CALL TF3LAG( RMCN(MNX), NDSM, NUTYEL, NUELEM, NBELEM,
     %                   NTDL, MCN(MNNDEL),  NBPNFX,
     %                   NBPOLA, NPIA,  MCN(MNPOIA),
     %                   MCN(MNPOLA),   MCN(MNDPOA),
     %                   NBPOLQ, NPIQ,  MCN(MNPOIQ),
     %                   MCN(MNPOLQ),   MCN(MNDPOQ),
     %                   NFACE, NOOBSF, NUMIOB(3), NUMAOB(3), NBPOLY,
     %                   NOOBVC, NUMIOB(4), NUMAOB(4), MCN(MNDOEL(4)),
     %                   MCN(MNDPA1), MCN(MNDPA2), MCN(MNTEMP),
     %                   MCN(MNFLPT), MCN(MNCOPN),
     %                   MCN(MNFLNP), MCN(MNFLTO) )
            GOTO 180
C
C           ****************************
C           TRIANGLE 2P1D LAGRANGE DROIT
C           ****************************
  113       CALL E12P1D( RMCN(MNX), MCN(MNF1), MCN(MNF2), MCN(MNDP) )
            CALL TF2P1D( RMCN(MNX), NDSM, NBELEM, NUELEM,
     %                   NTDL, MCN(MNNDEL),
     %                   NOOBLA, NUMIOB(2), NUMAOB(2),
     %                   NOOBSF(1), NUMIOB(3), NUMAOB(3),MCN(MNDOEL(3)),
     %                   MCN(MNTEMP), MCN(MNFLPT),
     %                   MCN(MNCOPN), MCN(MNFLNP), MCN(MNFLTO) )
            GOTO 180
C
C           *****************************
C           TETRAEDRE 3P1D LAGRANGE DROIT
C           *****************************
C           CALCUL DU FLUX SUR CHAQUE FACE ET FLUX NORMAL AU BARYCENTRE DES FACE
  119       CALL E13P1D( RMCN(MNX), MCN(MNF1),
     %                   DELTA, MCN(MNDFM1), MCN(MNDP) )
            CALL TF3P1D( RMCN(MNX), NDSM, NUELEM, NBELEM, MCN(MNDP),
     %                   NTDL, MCN(MNNDEL),
     %                   NOOBSF, NUMIOB(3), NUMAOB(3),
     %                   NOOBVC, NUMIOB(4), NUMAOB(4), MCN(MNDOEL(4)),
     %                   MCN(MNTEMP), MCN(MNFLPT),
     %                   MCN(MNCOPN), MCN(MNFLNP), MCN(MNFLTO) )
C           GOTO  180
C
C           -----------------------------------------------------------------
C           CALCUL DES COORDONNEES DES POINTS D'INTEGRATION NUMERIQUE DE L'EF
C           OU EST CALCULE LE GRADIENT DE LA TEMPERATURE
C           -----------------------------------------------------------------
  180       CALL CPINEL( NBELEM, NUELEM, NPI, NDIM,
     %                   MCN(MNF1), MCN(MNF2), MCN(MNF3), RMCN(MNCOPI) )
C
C           CALCUL GRADIENT DE LA TEMPERATURE AUX POINTS D'INTEGRATION DE L'EF
            CALL DERTEM( NBELEM, NUELEM, NDSM,
     %                   NDIM, NBPOLY, NPI, MCN(MNDP),
     %                   NTDL, MCN(MNNDEL),
     %                   MCN(MNTEMP), MCN(MNDTEM+WTEMPE) )
C
  190    CONTINUE
C        FIN DE LA BOUCLE SUR LES EF DE CE TYPE D'EF
C
C        MISE A JOUR DU TMS 'FLUXPT"NMTYEL'
C        ==================================
         MCN( MNFLUX + WBCAFX ) = NDSM
         MCN( MNFLUX + WDIMFX ) = NDIM
         MCN( MNFLUX + WUTYFX ) = NUTYEL
         MCN( MNFLUX + WBELFX ) = NBELEM
         MCN( MNFLUX + WBPNFX ) = NBPNFX
C        LA DATE
         CALL ECDATE( MCN(MNFLUX) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNFLUX + MOREE2 ) = NONMTD( '~>>>FLUXPT' )
C
C        MISE A JOUR DU TMS 'DTEMPERATURE"NMTYEL'
C        ========================================
         MCN( MNDTEM + WBJECD ) = NDSM
         MCN( MNDTEM + WDIMED ) = NDIM
         MCN( MNDTEM + WUTYED ) = NUTYEL
         MCN( MNDTEM + WBELFD ) = NBELEM
         MCN( MNDTEM + WBPEFD ) = NPI
C        LA DATE
         CALL ECDATE( MCN(MNDTEM) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNDTEM + MOREE2 ) = NONMTD( '~>>>DTEMPERATURE' )
C
  200 CONTINUE
C     FIN DE LA BOUCLE SUR LES TYPES D'EF
C
C     CREATION DU TMS FLUXFR des FLUX NORMAUX aux PLS FRONTIERES
C     ==========================================================
      NBFMIX = NUMAOB(NTYF) - NUMIOB(NTYF) + 1
      CALL TNMCDC( 'ENTIER', NBFMIX, MNFMIX )
      CALL  AZEROI( NBFMIX, MCN(MNFMIX) )
C     CALCUL DU NOMBRE DE PLS DES FRONTIERES DE L'OBJET
      NBPLSF = 0
      DO 220 N = NUMIOB(NTYF), NUMAOB(NTYF)
C        N EST IL UN PLS DE L'OBJET?
         MNNO = MNOBCL
         DO 210 I=1,NBOBCL
            IF( MCN(MNNO) .EQ. NDIM .AND. MCN(MNNO+1) .EQ. N ) GOTO 215
C           TYPE ET NUMERO CORRESPONDENT
C           SI LE POINT EN 1D, LA LIGNE EN 2D, LA SURFACE EN 3D
C           EST AUX LIMITES DE L'OBJET ALORS LE FLUX A ETE CALCULE
            MNNO = MNNO + 2
 210     CONTINUE
C        PLS NON DANS L'OBJET TRAITE
         GOTO 220
C        LE NO DU PLS N EST STOCKE
 215     MCN( MNFMIX + N - NUMIOB(NTYF) ) = N
         NBPLSF = NBPLSF + 1
 220  CONTINUE
C
      CALL LXTSOU( NTLXOB, 'FLUXFR', NTFLFR, MNFLFR )
      IF( NTFLFR .GT. 0 ) THEN
C        LE TABLEAU EST DETRUIT POUR ETRE REDECLARE
         CALL LXTSDS( NTLXOB, 'FLUXFR' )
      ENDIF
C
      MOFLFR = MOREE2 * NBPLSF * NDSM
      MOTS   = WLUXFR + MOFLFR + NBPLSF
      CALL LXTNDC( NTLXOB, 'FLUXFR', 'MOTS', MOTS )
      CALL LXTSOU( NTLXOB, 'FLUXFR', NTFLFR, MNFLFR )
      MCN( MNFLFR + WBCAFF ) = NDSM
      MCN( MNFLFR + WTPLSF ) = NDIM
      MCN( MNFLFR + WBPLSF ) = NBPLSF
C     COPIE DES FLUX AUX FRONTIERES DE L'OBJET
      MNFL = ( MNFLTO - 1 ) / MOREE2
      MNFX = ( MNFLFR + WLUXFR - 1 ) / MOREE2
      MNNO = MNFLFR + WLUXFR + MOFLFR
      NB   = 0
      DO N = 1, NBFMIX
         NUPLS = MCN(MNFMIX-1+N)
         IF( NUPLS .NE. 0 ) THEN
C           NO DU PLSF DANS SON LX
            MCN(MNNO+NB) = NUPLS
            NB = NB + 1
C           COPIE DES FLUX DE CE PLSF
            M1 = MNFL + N
            M2 = MNFX + NB
            DO K=1,NDSM
               DMCN(M2) = DMCN(M1)
               M1 = M1 + NBFMIX
               M2 = M2 + NBPLSF
            ENDDO
         ENDIF
      ENDDO
C     LA DATE
      CALL ECDATE( MCN(MNFLFR) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFLFR + MOREE2 ) = NONMTD( '~>>>FLUXFR' )
C
C     AFFICHAGE DES FLUXFR DE CHALEUR SUR LES OBJETS AUX LIMITES
C     ==========================================================
      CALL AFFLUXFR( KNOMOB )
C
C     DESTRUCTION DES TABLEAUX TEMPORAIRES
C     ====================================
 9999 IF( MNFLTO .GT. 0 ) CALL TNMCDS( 'REEL2',  MOFLTO, MNFLTO )
      IF( MNFLPT .GT. 0 ) CALL TNMCDS( 'REEL2',  MOFLPT, MNFLPT )
      IF( MNFMIX .GT. 0 ) CALL TNMCDC( 'ENTIER', NBFMIX, MNFMIX )
      RETURN
      END
