      SUBROUTINE TFLUIDE( KNOMOB, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES VITESSES, PRESSIONS, TOURBILLONS, PARTICULES
C -----    D'UN OBJET FLUIDE INCOMPRESSIBLE REGI PAR LES EQUATIONS
C          DE STOKES ou NAVIER STOKESDE, FLUIDE de NOM KNOMOB APRES
C          CALCUL de SES VITESSES et PRESSIONS AUX NOEUDS du MAILLAGE
C          a RECUPERER sur des FICHIERS du REPERTOIRE du PROJET

C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET CONTENANT LE FLUIDE

C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 2000
C MODIFS : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  OCTOBRE 2010
C MODIFS : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  JANVIER 2011
C MODIFS : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY    AVRIL 2012
C MODIFS : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY     MARS 2013
C MODIFS : ALAIN PERRONNET VEULETTES & ST PIERRE DU PERRAY     AOUT 2020
C MODIFS : ALAIN PERRONNET Saint PIERRE du PERRAY           FEVRIER 2021
C MODIFS : ALAIN PERRONNET Saint PIERRE du PERRAY               Mai 2021
C MODIFS : ALAIN PERRONNET Saint PIERRE du PERRAY           Fevrier 2022
C23456---------------------------------------------------------------012
      PARAMETER  (MXTYEL=7, LIGCON=0, MXVPFILE=2048)
      include"./incl/langue.inc"
      include"./incl/nmproj.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/msvaau.inc"
      include"./incl/trvari.inc"
      include"./incl/xvfontes.inc"
      include"./incl/mecoit.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___arete.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/xyzext.inc"
      include"./incl/ponoel.inc"
      include"./incl/donflu.inc"

      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      INTEGER           MCN
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(MOTMCN)
      DOUBLE PRECISION DMCN(MOTMCN/2)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      COMMON / MSSFTA / MSSF(28),NTADAM

      CHARACTER*(*)     KNOMOB
      CHARACTER*4       NOMELE(2)
      CHARACTER*12      NMSOLU
      CHARACTER*24      KNMSURF
      CHARACTER*128     NMVPFILE(MXVPFILE)
      CHARACTER*128     KNOMFIC
      INTEGER           NUMIOB(4), NUMAOB(4), MNDOEL(4), MXDOEL(4)
      INTEGER           NONOEF(10), VouW

      DOUBLE PRECISION  VITMIN, VITMAX, VITMOY, VITMOD, VITMX,
     %                  VITMI1, VITMA1, VITMO1,
     %                  PREMIN, PREMAX, PREMOY,
     %                  PREMI1, PREMA1, PREMO1,
     %                  TEMMIN, TEMMAX, TEMMOY,
     %                  PSIMIN, PSIMAX, ROTMIN, ROTMAX, ROTMOY,
     %                  SEMINVIT, SEMAXVIT,
     %                  DINFO,  DCPU, D, V, P, D2PI

      REAL              VIMIN,  VIMAX,  PRMIN,  PRMAX,  PSMIN, PSMAX,
     %                  VitErMIN, VitErMAX, PreErMIN, PreErMAX,
     %                  TEMIN, TEMAX,
     %                  ROMIN, ROMAX
      LOGICAL           AVANT

      EXTERNAL          ETTAEL

      DOUBLE PRECISION, allocatable, dimension(:)    :: PG,VXYZPN,TEMPER
      DOUBLE PRECISION, allocatable, dimension(:,:)  :: VPAUX
      DOUBLE PRECISION, allocatable, dimension(:,:,:):: ROTV

      INTRINSIC         ALLOCATED

C     DECLARATION BIDON POUR LE PASSAGE des types avec pointers...
      DOUBLE PRECISION VITXYZ(1:3, 1:2, 1:4 ), PRESSION(1:2,1:2)

      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ), dimension(:), allocatable :: vitx, vity, vitz,
     %                                                vitm, pres, temp

      IF( INTERA .LE. 0 ) THEN
C        DEMANDE DE TRACE EN MODE BATCH => ARRET DE MEFISTO
         CALL ARRET( 100 )
      ENDIF

      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10001)
      ELSE
         WRITE(IMPRIM,20001)
      ENDIF
10001 FORMAT(/120('-')/'TRACE de la VITESSE la PRESSION et les ERREURS'
     %       /120('-') )
20001 FORMAT(/120('-')/'DRAWINGS of VELOCITIES PRESSURES and ERRORS'
     %       /120('-') )

C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)

C     2 * Pi
      D2PI   = ATAN(1D0) * 8D0

C     INITIALISATIONS POUR LA REMANENCE DES VALEURS
      NDIM1  = 1
      VouW   = 0

C     INDICE=1 DE NON ALLOCATION des TABLEAUX DYNAMIQUES
      IALvitx    = 1
      IALvity    = 1
      IALvitz    = 1
      IALvitm    = 1
      IALpres    = 1
      IALtemp    = 1
      IALPG      = 1
      IALVXYZPN  = 1
      IALTEMPER  = 1
      IALVPAUX   = 1
      IALROTV    = 1

      NBSOMT   = 0
      NTDLVP   = 0
      NBVPFILE = 0
      NBCOMP   = 1
      NDPGSTP1 = 0
      MOXYZNP1 = 0
      MNXYZNP1 = 0
      MNNPEFP1 = 0
      MNELEP1  = 0
      MOELP1 = 0
      MNAGD  = 0
      NBVVIPR= 0
      MNNPEF = 0
      MNELE  = 0
      TEMPS  = 0.0
      NOPT   = 1
      CMFLEC = 2.5
      CMVITE = 0.0
      IERR   = 0
      MODECO = 0
      MNVECT = 0
      MNTIMES= 0
      MNTAUX = 0
      MOTAUX = 0
      MNTAEL = 0
      MOTAEL = 0
      MNX    = 0
      NBDLMX = 0
      MNNDDL = 0
      MNPSI  = 0
      MNSTFR = 0
      MNPTDG = 0
      MNLPCO = 0
      MONOSO = 0
      MNNOSO = 0
      MOSONO = 0
      MNSONO = 0
      MNWITE = 0
      MNVITE = 0
      MNFVSF = 0
      NBSFTR = 0
      MNNOSFTR = 0
      MNXYSFTR = 0
      MNEFSFTR = 0
      SEMINVIT = 0D0
      SEMAXVIT = 1D100
      MOAUX1 = 0
      MNAUX1 = 0
      NOCOMP = 1

      MNBPPAR  = 0
      MNDISTEM = 0
      MNNOITEM = 0

C     LE CAS A VISUALISER PAR DEFAUT
      NCAS0  = 1
      NCAST0 = 1
      NCAS1  = 1
      NCAST1 = 1
      NCAS   = 1

C     PAS DE TABLEAU DES XYZ VXYZ RAYON TEMPS DE PARTICULES
      NBPART = 0
      MNPART = 0

C     OUVERTURE DE L'OBJET FLUIDE
C     ---------------------------
      CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'tfluide ERREUR: OBJET INCONNU ' // KNOMOB
         ELSE
            KERR(1) = 'tfluide ERROR: UNKNOWN OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9997
      ENDIF

C     OUVERTURE DU TABLEAU DEFINITION DE L'OBJET FLUIDE
C     -------------------------------------------------
      CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )
C
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='tfluide ERREUR: DEFINITION INCONNUE de l''OBJET '
     %               // KNOMOB
         ELSE
            KERR(1)='tfluide ERROR: UNKNOWN DEFINITION of the OBJECT '
     %               // KNOMOB
         ENDIF
         CALL LEREUR
         IERR = 2
         GOTO 9997
      ENDIF

C     RECUPERATION DU MAILLAGE ELEMENTS FINIS DE L'OBJET FLUIDE
C     ---------------------------------------------------------
C     RECHERCHE DES TABLEAUX XYZSOMMET XYZNOEUD XYZPOINT NPEF" DE L'OBJET
      CALL MIMAOB( 1,      NTLXOB, MXDOFL, NTTOPO, MNTOPO,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %             MXTYEL, NBTYEL, MTNPEF, MNNPEF,
     %             NUMIOB, NUMAOB,
     %             NDPGST, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             MXDOEL, MNDOEL, IERR )
      IF( IERR .NE. 0 ) GOTO 9997
C     ICI NUMIOB CONTIENT LES 4 NUMEROS MINIMA DES PLSV
C         NUMAOB          LES 4 NUMEROS MAXIMA DES PLSV
C     MNDOEL LES 4 ADRESSES MCN DES TABLEAUX DES ADRESSES DES
C     TABLEAUX DECRIVANT LA PHYSIQUE DE L'OBJET COMPLET

C     LES POINTS SONT IDENTIQUES AUX NOEUDS   3/9/2010
      MNXYZP = MNXYZN
C
C     NOMBRE DE NOEUDS SUPPORT DE LA VITESSE DU MAILLAGE
      NBNOEU = MCN(MNXYZN+WNBNOE)

C     CALCULER DU MINIMUM ET MAXIMUM DES 3-ERES COORDONNEES
C     DES NOEUDS DANS LE TABLEAU COOEXT
      CALL MIMXPT( 3, NBNOEU, RMCN(MNXYZN+WYZNOE), COOEXT )

C     NDIM DIMENSION EFFECTIVE DE L'ESPACE DES COORDONNEES
      CALL DIMCOO( NBNOEU, MCN(MNXYZN+WYZNOE), NDIM )
      NDIMLI = NDIM
      NDIM1  = NDIM + 1

C     PAR DEFAUT TRACE du MODULE DE LA VITESSE EN POSITION NDIM+1 DES VITESSES
      NOCOMP = NDIM1
C
C     MNELE : ADRESSE DU TABLEAU NPEF"TYPE EF (P1D ou P2C)
      MNELE = MCN( MNNPEF )
C
C     NOMBRE D'ELEMENTS FINIS DE CE TYPE
      NBELEM = MCN( MNELE + WBELEM )
C
C     LE NUMERO DU TYPE DE L'ELEMENT FINI
C     (13:TRIA 2P1D ou 15:TRIA 2P2C ou 19:TETR 3P1D ou 20:TETR 3P2C)
      NUTYEL = MCN( MNELE + WUTYEL )
C
C     NDPGST: CODE TRAITEMENT DES XYZ DES SOMMETS POINTS NOEUDS
C              0 : NOEUDS=POINTS=SOMMETS
C              1 : NOEUDS=POINTS#SOMMETS
C              2 : NOEUDS#POINTS=SOMMETS
C              3 : NOEUDS#POINTS#SOMMETS
      NDPGST = MCN( MNTOPO + WDPGST )
C
C     LES CARACTERISTIQUES DE L'ELEMENT FINI
      CALL ELNUNM( NUTYEL, NOMELE )
C
C     LES CARACTERISTIQUES DE L'ELEMENT FINI
      CALL ELTYCA( NUTYEL )
C
      IF( NUTYEL .NE. 13 .AND. NUTYEL .NE. 19 .AND.
     %    NUTYEL .NE. 15 .AND. NUTYEL .NE. 20 ) THEN
          NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: ELEMENT FINI FUIDE INCONNU ' // NOMELE(1)
     %             // NOMELE(2)
         ELSE
            KERR(1)= 'ERROR: UNKNOWN FLUID FINITE ELEMENT ' // NOMELE(1)
     %             // NOMELE(2)
         ENDIF
         CALL LEREUR
         IERR = 4
         GOTO 9997
      ENDIF
C
      IF( NUTYEL .EQ. 13 .OR. NUTYEL .EQ. 19 ) THEN
C
C        EF BREZZI-FORTIN P1+BULLE P1:  AJOUT DES DL VITESSE AUX BARYCENTRES
C        -------------------------------------------------------------------
C        NOMBRE DE SOMMETS DU MAILLAGE
         IF( NBNOEU .GT. NBELEM ) THEN
            NBSOMT = NBNOEU - NBELEM
         ELSE
            NBSOMT = NBNOEU
         ENDIF
C        NBSOMT = NOMBRE DE NOEUDS-SOMMETS SUPPORT DE LA PRESSION
C        NOMBRE DE NOEUDS SUPPORTANT LA VITESSE LES BARYCENTRES EN PLUS
         NBNOVI = NBSOMT + NBELEM
C        INTERPOLATION DE LA VITESSE EN P1+BULLE
         INTERPV = 3

C        TABLEAU NONOSO INEXISTANT
         MNNOSO = 1
C        TABLEAU NOSONO INEXISTANT
         MNSONO = 1
C
C        AJOUT DES XYZ DU BARYCENTRE DES EF BREZZI-FORTIN COMME NOEUDS
C        DANS LE TMS XYZNOEUD
C       (MAIS PAS DE NO DU BARYCENTRE DANS NPEF".P1D)
         CALL XYZNOEBF( NTLXOB, MNNPEF, NTXYZN, MNXYZN )

C        EF BREZZI-FORTIN  => 1 SEUL MAILLAGE P1 DES XYZ DES SOMMETS
C                             1 SEUL MAILLAGE P1 NPEF
         MNXYZNP1 = MNXYZN
         MOELP1   = 0
         MNELEP1  = MNELE
         MNNPEFP1 = MNNPEF

C        ICI 
C
      ELSE IF( NUTYEL .EQ. 15 .OR. NUTYEL .EQ. 20 ) THEN
C
C        EF TAYLOR-HOOD P2: CONSTRUCTION DU MAILLAGE P1 A PARTIR de P2
C        -------------------------------------------------------------
C        CONSTRUCTION DU TABLEAU NOSO: NO NOEUD  => NO SOMMET
C        ....................................................
         IF( MONOSO .GT. 0 ) CALL TNMCDS('ENTIER', MONOSO, MNNOSO )
         MONOSO = NBNOEU
         CALL TNMCDC( 'ENTIER', MONOSO, MNNOSO )
         IF( MNNOSO .LE. 0 ) GOTO 9900
         CALL AZEROI( NBNOEU,  MCN(MNNOSO) )
         MNNOS1 = MNNOSO-1
C        NDIM1=NDIM+1 = NOMBRE DE SOMMETS D'UN EF (TRIANGLE ou TETRAEDRE)
         DO NUELEM = 1, NBELEM
            CALL EFNOEU( MNELE, NUELEM, NBNOEF, NONOEF )
            DO I = 1, NDIM1
C              TEMOIN DE SOMMET
               MCN( MNNOS1+NONOEF(I) ) = 1
            ENDDO
         ENDDO
C        NUMEROTATION CROISSANTE DES SOMMETS
         NBSOMT = 0
         DO I = 1, NBNOEU
            IF( MCN(MNNOS1+I ) .GT. 0 ) THEN
               NBSOMT = NBSOMT + 1
               MCN( MNNOS1 + I ) = NBSOMT
            ENDIF
         ENDDO

C        CONSTRUCTION DU TABLEAU SONO: NO SOMMET => NO NOEUD
C        ...................................................
         IF( MOSONO .GT. 0 ) CALL TNMCDS('ENTIER', MOSONO, MNSONO )
         MOSONO = NBSOMT
         CALL TNMCDC( 'ENTIER', MOSONO, MNSONO )
         IF( MNSONO .LE. 0 ) GOTO 9900
         CALL AZEROI( NBSOMT,  MCN(MNSONO) )
         MNSON1 = MNSONO - 1
C        NUMEROTATION CROISSANTE DES SOMMETS
         NBSOMT = 0
         DO I = 1, NBNOEU
            IF( MCN(MNNOS1+I ) .GT. 0 ) THEN
               NBSOMT = NBSOMT + 1
               MCN( MNSON1 + NBSOMT ) = I
            ENDIF
         ENDDO

C        INTERPOLATION P2 DE LA VITESSE
         INTERPV = 2

C        NBSOMT = NOMBRE DE NOEUDS-SOMMETS SUPPORT DE LA PRESSION
C        NOMBRE DE NOEUDS P2 SUPPORTANT LA VITESSE
         NBNOVI = NBNOEU

C        CONSTRUCTION TABLEAU NPEF"P1 NO DES NDIM+1 SOMMETS DE CHAQUE EF
C        POUR LES TRACES DES SOLUTIONS P1: PRESSION, ROTV, ...
         CALL TNMCDC( 'ENTIER', 1, MNNPEFP1 )
         MOELP1  = WUNDEL + NBELEM * NDIM1
         CALL TNMCDC( 'ENTIER', MOELP1, MNELEP1 )
C        MNELEP1 : ADRESSE DU TABLEAU NPEF"TYPE EF (P1D)
         MCN( MNNPEFP1 ) = MNELEP1
C
C        variable NUTYEL 'Numero du type des elements finis' nP1D
         IF( NDIM .EQ. 2 ) THEN
C           'TRIA','2P1D'
            NUTYP1 = 13
         ELSE
C           'TETR','3P1D'
            NUTYP1 = 19
         ENDIF
         MCN( MNELEP1 + WUTYEL ) = NUTYP1
C        variable NBELEM 'Nombre d''EF de ce TYPE'    
         MCN( MNELEP1 + WBELEM ) = NBELEM
C        variable NBELAP 'Nombre d''EF avec POINTEUR sur EF a TG'  
         MCN( MNELEP1 + WBELAP ) = 0
C        variable NBELTG 'Nombre d''EF a tg'          
         MCN( MNELEP1 + WBELTG ) = 0
C        variable NBNDEL 'Nombre de noeuds d''un EF'  
         MCN( MNELEP1 + WBNDEL ) = NDIM1
C        variable NBPGEL 'Nombre de points d''un EF (0 si points=noeuds)' 
         MCN( MNELEP1 + WBPGEL ) = 0
C        variable NBTGEL 'Nombre de tangentes d''un EF'
         MCN( MNELEP1 + WBTGEL ) = 0
C        variable NBPSEL 'Nombre de SOMMETS d''EF sur un POINT utilisateur'  
         MCN( MNELEP1 + WBPSEL ) = 0
C        variable MOPSEL 'Nombre de variables des tableaux NLPSEL et NEPSEL' 
         MCN( MNELEP1 + WOPSEL ) = 0
C        variable NBLAEL 'Nombre d''ARETES d''EF sur une LIGNE utilisateur'  
         MCN( MNELEP1 + WBLAEL ) = 0
C        variable MOLAEL 'Nombre de variables des tableaux NLLAEL et NELAEL' 
         MCN( MNELEP1 + WOLAEL ) = 0
C        variable NBSFEL 'Nombre de FACES d''EF sur une SURFACE utilisateur' 
         MCN( MNELEP1 + WBSFEL ) = 0
C        variable MOSFEL 'Nombre de variables des tableaux NLSFEL et NESFEL' 
         MCN( MNELEP1 + WOSFEL ) = 0
C        variable NBVCEL 'Nombre d''EF dans un VOLUME utilisateur'  
         MCN( MNELEP1 + WBVCEL ) = 0
C        variable MOVCEL 'Nombre de variables des tableaux NLVCEL et NEVCEL'  
         MCN( MNELEP1 + WOVCEL ) = 0
C        tableau NUNDEL(1..NBELEM,1..NBNDEL) 'Numero des NOEUDS de chaque EF' 
C        ici le NUMERO EST CELUI DU SOMMET DANS LA NUMEROTATION GLOBALE des SOMMETS
         MNP1 = MNELEP1 + WUNDEL -1
         MNP2 = MNELE   + WUNDEL -1
         DO K = 0, NDIM
            DO N = 1, NBELEM
C              LE NUMERO GLOBAL DE SOMMET DU NOEUD K+1 DE L'EF N
               MCN( MNP1 + N ) = MCN( MNNOS1 + MCN( MNP2 + N ) )
            ENDDO
            MNP1 = MNP1 + NBELEM
            MNP2 = MNP2 + NBELEM
         ENDDO
C        LA DATE
         CALL ECDATE( MCN(MNELEP1) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNELEP1 + MOREE2 ) = NONMTD( '~>>>NPEF' )

C        LE TABLEAU XYZSOMMET DES SOMMETS DU MAILLAGE
         MOXYZNP1 = WYZPOI + 3 * NBSOMT
         CALL TNMCDC( 'ENTIER', MOXYZNP1, MNXYZNP1 )
C        variable NBNOE 'Nombre de noeuds' P1
         MCN( MNXYZNP1 + WNBNOE ) = NBSOMT
C        variable NBTGN 'Nombre de tangentes' 
         MCN( MNXYZNP1 + WNBTGN ) = 0
C        variable NBCOON 'Nombre coordonnees d''un point'
         MCN( MNXYZNP1 + WBCOON ) = 3
         MNNOS1 = MNNOSO - 1
         DO N = 1, NBNOEU
C           NS NUMERO DE SOMMET DU NOEUD N
            NS = MCN( MNNOS1 + N )
            IF( NS .GT. 0 ) THEN
               MNP1 = MNXYZNP1+ WYZNOE + 3 * NS - 4
               MNP2 = MNXYZN  + WYZNOE + 3 * N  - 4
               DO K = 1, 3
                  RMCN(MNP1+K) = RMCN(MNP2+K)
               ENDDO
            ENDIF
         ENDDO
C        LA DATE
         CALL ECDATE( MCN(MNXYZNP1) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNXYZNP1 + MOREE2 ) = NONMTD( '~>>>XYZPOINT' )

      ENDIF

C     INTERPOLATION P1 DE LA PRESSION
      INTERPP = 1

C     AFFICHAGE DES DONNEES
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10002) INTERPV, NBNOVI, INTERPP, NBSOMT
      ELSE
         WRITE(IMPRIM,20002) INTERPV, NBNOVI, INTERPP, NBSOMT
      ENDIF

10002 FORMAT('TYPE d''INTERPOLATION de la VITESSE', T36,I9/
     %       'NOMBRE de NOEUDS de la VITESSE',      T36,I9/
     %       'TYPE d''INTERPOLATION de la PRESSION',T36,I9/
     %       'NOMBRE de NOEUDS de la PRESSION',     T36,I9/)

20002 FORMAT('VELOCITY INTERPOLATION TYPE NUMBER',T35,I9/
     %       'VELOCITY NODE NUMBER',              T35,I9/
     %       'PRESSURE INTERPOLATION TYPE NUMBER',T35,I9/
     %       'PRESSURE NODE NUMBER',T35,I9/)

C     RECONSTRUCTION DU TABLEAU NDDLNO(0:NBNOEU)
C     POINTEUR SUR LE DERNIER DL DE CHAQUE NOEUD VITESSE DU FLUIDE
C     ------------------------------------------------------------
      MNNDDL = 0
      CALL TNMCDC( 'ENTIER', 1+NBNOVI, MNNDDL )
      IF( MNNDDL .LE. 0 ) GOTO 9900
      CALL AZEROI( 1+NBNOVI, MCN(MNNDDL) )
      CALL PTDLFL( MNELE,  NDIM, NBELEM, NBNOEU, NUTYEL,
     %             NBSOMT, NBNOVI, MCN(MNNDDL), IERR )
      IF( IERR .NE. 0 ) GOTO 9997

C     RECUPERATION OBLIGATOIRE des TEMPS a partir du NOM des FICHIERS des
C     VECTEURS VITESSE+PRESSION SAUVEGARDES SUR DISQUE dans le
C     REPERTOIRE du PROJET
C     -------------------------------------------------------------------
      CALL LINMFIVPT( KNOMOB,  MXVPFILE,
     %                MNTIMES, NMVPFILE, NBVPFILE, NCAST0, NCAST1 )
C     LE TEMPS DE CHAQUE FICHIER VECTEUR VITESSE+PRESSION
C     RMCN(MNTIMES:MNTIMES-1+NBVPFILE)=TIMES(NCAST0:NCAST1)
      IF( NBVPFILE .LE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
             PRINT*,'tfluide: AUCUN VECTEUR VITESSE+PRESSION dans le REP
     %ERTOIRE du PROJET ',NMPROJ
         ELSE
             PRINT*,'tfluide: NO VELOCITY+PRESSURE VECTOR in PROJECT LEX
     %ICON ',NMPROJ
         ENDIF
         GOTO 9997
      ENDIF

C     LE CAS A VISUALISER PAR DEFAUT
      NCAS0 = NCAST0
      NCAS1 = NCAST1
      NCAS  = NCAST1

C     TEMPS DU DERNIER CAS TIMES(NCAST1) LUS DES FICHIERS
      TEMPS = RMCN( MNTIMES-1+NBVPFILE )


C     ==================================================================
C     LECTURE DES DONNEES DE TRACES DE L'UTILISATEUR
C     ==================================================================

C     LECTURE du Nombre de SURFACES a TRACER avec LES VITESSES-PRESSIONS
C     ------------------------------------------------------------------
      IF( NDIM .EQ. 3 ) THEN
         CALL INVITE( 37 )
         CALL LIRENT( NCVALS, NBSFTR )
         IF( NCVALS .EQ. -1 ) GOTO 9997
      ELSE
         NBSFTR = 0
      ENDIF

      IF( NBSFTR .LE. 0 ) THEN
         NBSFTR = 0
         GOTO 15
      ENDIF

C     OUVERTURE DU LEXIQUE DES SURFACES
      CALL LXLXOU( NTADAM, 'SURFACE', NTSURF, MNSURF )

C     LECTURE DU NOM DES NBSFTR A TRACER AVEC LES VITESSES-PRESSIONS
C     --------------------------------------------------------------
C     RESERVATION DU TABLEAU NOSFTR
      CALL TNMCDC( 'ENTIER', NBSFTR, MNNOSFTR )
      CALL TNMCDC( 'ENTIER', NBSFTR, MNXYSFTR )
      CALL TNMCDC( 'ENTIER', NBSFTR, MNEFSFTR )
      IF( MNEFSFTR .LE. 0 ) GOTO 9997

C     NBDOBJ = NOMBRE DE PLSV de L'OBJET
      NBDOBJ = MCN( MNDFOB + WBDOBJ )
      MNTYOB = MNDFOB + WTYOBJ - 2

      DO K = 1, NBSFTR

C        INVITE POUR ENTRER LE NOM DE LA SURFACE
 10      CALL INVITE( 43 )
         CALL OBJENU( 'SURFACE', NUSURFTR )

         IF( NUSURFTR .LT. 0 ) THEN
C           '@' DONNE
            GOTO 9997
         ENDIF

         IF( NUSURFTR .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='tfluide: NOM INCONNU de SURFACE. A REDONNER'
            ELSE
               KERR(1)='tfluide: UNKNOWN SURFACE NAME. GIVE AGAIN'
            ENDIF
            CALL LEREUR
            GOTO 10
         ENDIF

C        VERIFICATION: CETTE SURFACE EST ELLE UNE SURFACE DE L'OBJET?
         N = 0
         DO J = 1, NBDOBJ
            MN = MNTYOB + 2 * J
            IF( MCN(MN) .EQ. 3 ) THEN
C              LE TYPE EST CELUI D'UNE SURFACE
               N = N + 1
               IF( MCN(MN+1) .EQ. NUSURFTR ) GOTO 12
            ENDIF
         ENDDO

C        SURFACE NON RETROUVEE DANS L'OBJET
         NBLGRC(NRERR) = 2 + N
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'tfluide: NOM DE SURFACE NON RECONNUE PARMI'
            KERR(2) = 'LES NOMS DE SURFACE DE L''OBJET'
         ELSE
            KERR(1) = 'tfluide: SURFACE NAME NOT KNOWN AMONG'
            KERR(2) = 'SURFACE NAMES of the OBJECT'
         ENDIF
         N = 2
         DO J = 1, NBDOBJ
            MN = MNTYOB + 2 * J
            IF( MCN(MN) .EQ. 3 ) THEN
C              LE TYPE EST CELUI D'UNE SURFACE
               N = N + 1
C              LE NOM DE LA SURFACE
               CALL NMOBNU( 'SURFACE', MCN(MN+1), KERR(N) )
            ENDIF
         ENDDO
         CALL LEREUR
         GOTO 10

C        SURFACE DEJA DONNEE?
 12      DO J = 1, K-1
            IF( MCN( MNNOSFTR-1+J ) .EQ. NUSURFTR ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'NOM DE SURFACE DEJA DONNE'
               ELSE
                  KERR(1) = 'SURFACE NAME ALREADY GIVEN'
               ENDIF
               CALL LEREUR
               GOTO 10
            ENDIF
         ENDDO

C        SURFACE RETROUVEE DANS L'OBJET
         MCN( MNNOSFTR - 1 + K ) = NUSURFTR

C        LE NOM DE LA SURFACE A TRACER
         CALL NMOBNU( 'SURFACE', NUSURFTR, KNMSURF )
         CALL SANSBL( KNMSURF, NBK )

C        OUVERTURE DE LA SURFACE NUSURFTR
         CALL LXNLOU( NTSURF, NUSURFTR, NTLXSFTR, MNLXSFTR )
         IF( NTLXSFTR .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'LA SURFACE '//KNMSURF(1:NBK)//' EST INCONNUE'
            ELSE
               KERR(1) = 'SURFACE '//KNMSURF(1:NBK)//' IS UNKNOWN'
            ENDIF
            CALL LEREUR
            GOTO 10
         ENDIF

C        RECUPERATION DES COORDONNEES XYZ DES SOMMETS DE LA SURFACE K
         CALL LXTSOU( NTLXSFTR, 'XYZSOMMET', NTXYZSU, MNXYZSU )
         IF( NTXYZSU .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='SURFACE '//KNMSURF(1:NBK)//' SANS XYZSOMMET'
            ELSE
               KERR(1)='SURFACE '//KNMSURF(1:NBK)//' WITHOUT XYZSOMMET'
            ENDIF
            CALL LEREUR
            GOTO 10
         ENDIF
         MCN( MNXYSFTR - 1 + K ) = MNXYZSU

C        RECUPERATION DES NO DES SOMMETS DES EF DE LA SURFACE K
         CALL LXTSOU( NTLXSFTR, 'NSEF', NTNSEFSU, MNNSEFSU )
         IF( NTNSEFSU .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='SURFACE '//KNMSURF(1:NBK)//' SANS NSEF'
            ELSE
               KERR(1)='SURFACE '//KNMSURF(1:NBK)//' WITHOUT NSEF'
            ENDIF
            CALL LEREUR
            GOTO 10
         ENDIF
         MCN( MNEFSFTR - 1 + K ) = MNNSEFSU

C        LE NOM DE LA SURFACE A TRACER EST AFFICHE
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'tfluide: La surface ',KNMSURF(1:NBK),' SERA TRACEE'
         ELSE
            PRINT*,'tfluide: ',KNMSURF(1:NBK),' SURFACE WILL BE DRAWN'
         ENDIF

      ENDDO

C     VALEUR PAR DEFAUT DES PARAMETRES DE TRACES ( COULEURS, ...)
C     AVANT RETOUCHE POSSIBLE A L'AIDE DES MENUS DE TRACES
C     -----------------------------------------------------------
      LCRITR = 0
C     TRACE DES FACES
      IAVFAC = 1
C     TRACE DES ARETES
      IAVARE = 1
C     COULEUR PAR DEFAUT DES ARETES DES FACES FRONTIERE
      NCOAFR = NCGRIC
C     COULEUR PAR DEFAUT DES ARETES DES FACES
      NCOUAF = NCGRIM


C     LECTURE DU CHOIX DES TEMPS DES FICHIERS VITESSES+PRESSIONS A TRACER
C     DEFINITION des VECTEURS VITESSES PRESSIONS ERREURS a RETROUVER
C     SOIT PAR NUMERO INITIAL et FINAL DES FICHIERS SUPPORT VECTEURS V+P
C     SOIT PAR TEMPS  INITIAL et FINAL DES FICHIERS SUPPORT VECTEURS V+P
C     =====================================================================
 15   CALL LIMTCL( 'nbvitpre', NMTCL )
      IF( NMTCL .LE. 0  ) GOTO 9997

      IF( NMTCL .EQ. 1 ) THEN

C        Numeros INITIAL et FINAL des VECTEURS SOLUTIONS a TRACER
C        --------------------------------------------------------
         CALL INVITE( 158 )
         CALL LIRENT( NCVALS, NCAS0 )
         IF( NCVALS .EQ. -1 ) GOTO 15
         IF( NCAS0 .LT. NCAST0 ) NCAS0 = NCAST0

         CALL INVITE( 159 )
         CALL LIRENT( NCVALS, NCAS1 )
         IF( NCVALS .EQ. -1 ) GOTO 15

         IF( NCAS0 .GT. NCAS1 ) THEN
C           PERMUTATION DES CAS
            N     = NCAS0
            NCAS0 = NCAS1
            NCAS1 = N
         ENDIF

         IF( NCAS0 .LT. NCAST0 ) NCAS0 = NCAST0
         IF( NCAS1 .GT. NCAST1 ) NCAS1 = NCAST1

      ELSE

C        TEMPS INITIAL et FINAL des VECTEURS SOLUTIONS a TRACER
C        ------------------------------------------------------
         CALL INVITE( 95 )
         CALL LIRRSP( NCVALS, TEMPS0 )
         IF( NCVALS .EQ. -1 ) GOTO 15

         CALL INVITE( 94 )
         CALL LIRRSP( NCVALS, TEMPS1 )
         IF( NCVALS .EQ. -1 ) GOTO 15

         IF( TEMPS0 .GT. TEMPS1 ) THEN
C           PERMUTATION DES TEMPS
            T      = TEMPS0
            TEMPS0 = TEMPS1
            TEMPS1 = T
         ENDIF

         IF( TEMPS0 .LT. RMCN(MNTIMES) ) TEMPS0 = RMCN(MNTIMES)
         IF( TEMPS1 .GT. RMCN(MNTIMES-1+NBVPFILE) )
     %       TEMPS1 = RMCN(MNTIMES-1+NBVPFILE)

C        RECHERCHE DU NO NCAS0 DU TEMPS LE PLUS PROCHE DE TEMPS0
         TPROCH = RINFO( 'GRAND' )
         DO I=1,NBVPFILE
            TTT = ABS( RMCN(MNTIMES-1+I) - TEMPS0 )
            IF( TTT .LT. TPROCH ) THEN
               NCAS0  = NCAST0-1+I
               TPROCH = TTT
            ENDIF
         ENDDO
         NCAS0 = MAX( NCAS0, NCAST0 )

C        RECHERCHE DU NO NCAS1 DU TEMPS LE PLUS PROCHE DE TEMPS1
         TPROCH = RINFO( 'GRAND' )
         DO I=1,NBVPFILE
            TTT = ABS( RMCN(MNTIMES-1+I) - TEMPS1 )
            IF( TTT .LT. TPROCH ) THEN
               NCAS1  = NCAST0-1+I
               TPROCH = TTT
            ENDIF
         ENDDO
         NCAS1 = MIN( NCAS1, NCAST1 )

      ENDIF

C     NOMBRE DEMANDE DE VECTEURS VITESSE+PRESSION A TRAITER
      NBVVIPR = NCAS1 - NCAS0 + 1

C     RECUPERATION DES VECTEURS VITESSE+PRESSION NCAS0 a NCAS1 SUR FICHIERS
C     =====================================================================
      TEMPS0 = RMCN( MNTIMES + NCAS0 - NCAST0 )
      TEMPS1 = RMCN( MNTIMES + NCAS1 - NCAST0 )

      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'tfluide: Recuperation des VECTEURS VITESSE+PRESSION',
     %           NCAS0,' a',NCAS1,' des TEMPS',TEMPS0,' a',TEMPS1
      ELSE
         PRINT*,'tfluide: Restore VELOCITY+PRESSURE VECTORS',
     %           NCAS0,' to',NCAS1,' of TIMES',TEMPS0,' to',TEMPS1
      ENDIF
      IF( NBVVIPR .LE. 0 ) GOTO 15

C     OUVERTURE ET LECTURE PARTIELLE DU FICHIER NCAS0 POUR RECUPERER
C     LA VALEUR DE NTDLVP NOMBRE DE DL D'UN VECTEUR VITESSE+PRESSION
C     LA VALEUR DE NTDLTE NOMBRE DE DL DE LA TEMPERATURE du FLUIDE
      TEMPS  = RMCN( MNTIMES + NCAS0 - NCAST0 )
      NTDLVP = -1
      NTDLTE = 0
      CALL LIFIVIPRTE( KNOMOB,  TEMPS,   NCAS0,  NAVSTO, NBPASDT,
     %                 NDIM,    NBNOVI,  NBSOMT, NBNOTE,
     %                 NTDLVP,  DMCN(1), NTDLTE, DMCN(1),
     %                 KNOMFIC, IERR )
      IF( IERR .NE. 0 ) THEN
         GOTO 9997
      ENDIF

C     ALLOCATION DYNAMIQUE d'un VECTEUR TEMPERATURE
C     ---------------------------------------------
      IF( NTDLTE .EQ. 0 ) THEN
         NTDLTE0 = 1
      ELSE
         NTDLTE0 = NTDLTE
      ENDIF
      IF( IALTEMPER .EQ. 0 ) DEALLOCATE( TEMPER )
      ALLOCATE( TEMPER(1:NTDLTE0), STAT=IALTEMPER )

      IF( IALTEMPER .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'ERREUR en ALLOCATION de',NTDLTE0,
     %             ' DOUBLE PRECISION de TEMPER(',NTDLTE0,')'
         ELSE
            PRINT*, 'ALLOCATION ERROR of',NTDLTE0,
     %              ' DOUBLE PRECISION of TEMPER(',NTDLTE0,')'
         ENDIF
         IERR = IALTEMPER
         GOTO 9997
      ELSE
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'CORRECTE ALLOCATION de',NTDLTE0,
     %             ' DOUBLE PRECISION de TEMPER(',NTDLTE0,')'
         ELSE
            PRINT*, 'ALLOCATION CORRECT of',NTDLTE0,
     %              ' DOUBLE PRECISION of TEMPER(',NTDLTE0,')'
         ENDIF
      ENDIF

C     ALLOCATION DYNAMIQUE d'un VECTEUR VITESSEXYZ+PRESSION
C     -----------------------------------------------------
      IF( IALVXYZPN .EQ. 0 ) DEALLOCATE( VXYZPN )
      ALLOCATE( VXYZPN(1:NTDLVP), STAT=IALVXYZPN )

      IF( IALVXYZPN .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'ERREUR en ALLOCATION de',NTDLVP,
     %             ' DOUBLE PRECISION de VXYZPN(',NTDLVP,')'
         ELSE
            PRINT*, 'ALLOCATION ERROR of',NTDLVP,
     %              ' DOUBLE PRECISION of VXYZPN(',NTDLVP,')'
         ENDIF
         IERR = IALVXYZPN
         GOTO 9997
      ELSE
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'CORRECTE ALLOCATION de',NTDLVP,
     %             ' DOUBLE PRECISION de VXYZPN(',NTDLVP,')'
         ELSE
            PRINT*, 'ALLOCATION CORRECT of',NTDLVP,
     %              ' DOUBLE PRECISION of VXYZPN(',NTDLVP,')'
         ENDIF
      ENDIF

C     ALLOCATION DYNAMIQUE du TABLEAU des NBVVIPR POINTEURS sur les VECTEURS
C     VX VY VZ VMOD PRESSION TEMPERATURE des NBVVIPR TEMPS A LIRE SUR FICHIERS
C     ------------------------------------------------------------------------
      allocate( vitx( NCAS0:NCAS1 ), STAT=IALvitx )
      IF(IALvitx .EQ. 0 ) GOTO 13
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'ERREUR en ALLOCATION de vitx(',NCAS0,':',
     %           NCAS1,')'
      ELSE
         PRINT*,'ALLOCATION ERROR of vitx(',NCAS0,':',
     %           NCAS1,')'
      ENDIF
      IERR = 5
      GOTO 9997

 13   allocate( vity( NCAS0:NCAS1 ), STAT=IALvity )
      IF(IALvity .EQ. 0 ) GOTO 14
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'ERREUR en ALLOCATION de vity(',NCAS0,':',
     %           NCAS1,')'
      ELSE
         PRINT*,'ALLOCATION ERROR of vity(',NCAS0,':',
     %           NCAS1,')'
      ENDIF
      IERR = 5
      GOTO 9997

 14   allocate( vitz( NCAS0:NCAS1 ), STAT=IALvitz )
      IF(IALvitz .EQ. 0 ) GOTO 16
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'ERREUR en ALLOCATION de vitz(',NCAS0,':',
     %           NCAS1,')'
      ELSE
         PRINT*,'ALLOCATION ERROR of vitz(',NCAS0,':',
     %           NCAS1,')'
      ENDIF
      IERR = 5
      GOTO 9997

 16   allocate( vitm( NCAS0:NCAS1 ), STAT=IALvitm )
      IF(IALvitm .EQ. 0 ) GOTO 17
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'ERREUR en ALLOCATION de vitm(',NCAS0,':',
     %           NCAS1,')'
      ELSE
         PRINT*,'ALLOCATION ERROR of vitm(',NCAS0,':',
     %           NCAS1,')'
      ENDIF
      IERR = 5
      GOTO 9997

 17   allocate( pres( NCAS0:NCAS1 ), STAT=IALpres )
      IF(IALpres .EQ. 0 ) GOTO 18
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'ERREUR en ALLOCATION de pres(',NCAS0,':',
     %           NCAS1,')'
      ELSE
         PRINT*,'ALLOCATION ERROR of pres(',NCAS0,':',
     %           NCAS1,')'
      ENDIF
      IERR = 5
      GOTO 9997

 18   IF( NTDLTE .GT. 0 ) THEN
         allocate( temp( NCAS0:NCAS1 ), STAT=IALtemp )
         IF(IALtemp .EQ. 0 ) GOTO 19
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'ERREUR en ALLOCATION de temp(',NCAS0,':',
     %              NCAS1,')'
         ELSE
            PRINT*,'ALLOCATION ERROR of temp(',NCAS0,':',
     %              NCAS1,')'
         ENDIF
         IERR = 5
         GOTO 9997
      ENDIF

 19   IF( LANGAG .EQ. 0 ) THEN
        PRINT*,'CORRECTE ALLOCATION de vitx y z m pres temp(',NCAS0,':',
     %          NCAS1,')'
      ELSE
        PRINT*, 'CORRECT ALLOCATION of vitx y z m pres temp(',NCAS0,':',
     %           NCAS1,')'
      ENDIF

C     POUR EVITER LES DESAGREMENTS LORS D'UNE DEALLOCATION INTEMPESTIVE
      do n = ncas0, ncas1
         IF( NTDLTE .GT. 0 ) NULLIFY( temp( n )%dptab )
         NULLIFY( pres( n )%dptab )
         NULLIFY( vitm( n )%dptab )
         NULLIFY( vitz( n )%dptab )
         NULLIFY( vity( n )%dptab )
         NULLIFY( vitx( n )%dptab )
       enddo

C     INITIALISATION DES TABLEAUX
C     vitx(NCAS0:NCAS1)%dptab(NBNOVI)
C     vity(NCAS0:NCAS1)%dptab(NBNOVI)
C     vitz(NCAS0:NCAS1)%dptab(NBNOVI)
C     vitm(NCAS0:NCAS1)%dptab(NBNOVI)
C     pres(NCAS0:NCAS1)%dptab(NBSOMT)
C     A PARTIR des VECTEURS NCAS0 a NCAS1 VITESSEPRESSION LUS
C     sur les FICHIERS du REPERTOIRE du PROJET et
C     CALCUL de la VITESSE MIN MAX et MOYENNE des NBVVIPR VECTEURS
C     ============================================================
      NBVVIPRLU= 0
      NOEVMIN  = 0
      NCAVMIN  = 0
      NOEVMAX  = 0
      NCAVMAX  = 0

      VITMIN = 1D111
      VITMAX =-1D111
      VITMOY = 0D0

      NOEPMIN = 0
      NCAPMIN = 0
      NOEPMAX = 0
      NCAPMAX = 0

      PREMIN = 1D111
      PREMAX =-1D111
      PREMOY = 0D0

      TEMMIN = 1D111
      TEMMAX =-1D111
      TEMMOY = 0D0

      IF( NUTYEL .EQ. 15 .OR. NUTYEL .EQ. 20 ) THEN

C        TAYLOR-HOOD
C        ===========
         DO 22 NCAS = NCAS0, NCAS1, 1

C           LECTURE DU VECTEUR VITESSE-PRESSION AU TEMPS NOVVIPR PAR NOEUDS
C           FICHIER DU REPERTOIRE PROJET DANS LE TABLEAU MNVXYZP0
            TEMPS = RMCN( MNTIMES + NCAS - NCAST0 )

            CALL LIFIVIPRTE( KNOMOB,  TEMPS,  NCAS,   NAVSTO, NBPASDT,
     %                       NDIM,    NBNOVI, NBSOMT, NBNOTE,
     %                       NTDLVP,  VXYZPN, NTDLTE, TEMPER,
     %                       KNOMFIC, IERR )
            IF( IERR .NE. 0 ) GOTO 9997

C           NUMERO DU VECTEUR LU DANS VX VY VZ PR
            NBVVIPRLU = NBVVIPRLU + 1

            allocate( vitx( NCAS )%dptab(1:NBNOVI),STAT=IALLOC)
            IF( IALLOC .NE. 0 ) GOTO 9990

            allocate( vity( NCAS )%dptab(1:NBNOVI),STAT=IALLOC)
            IF( IALLOC .NE. 0 ) GOTO 9990

            allocate( vitz( NCAS )%dptab(1:NBNOVI),STAT=IALLOC)
            IF( IALLOC .NE. 0 ) GOTO 9990

            allocate( vitm( NCAS )%dptab(1:NBNOVI),STAT=IALLOC)
            IF( IALLOC .NE. 0 ) GOTO 9990

            allocate( pres( NCAS )%dptab(1:NBSOMT),STAT=IALLOC)
            IF( IALLOC .NE. 0 ) GOTO 9990

            IF( NTDLTE .GT. 0 ) THEN
               allocate( temp( NCAS )%dptab(1:NBNOVI),STAT=IALLOC)
               IF( IALLOC .NE. 0 ) GOTO 9990
            ENDIF

C           NTDLVP LE NOMBRE DE DEGRES DE LIBERTE VITESSES + PRESSIONS
            DO I = 1, NBNOVI
C              LE NUMERO DU DL AVANT LE NOEUD I
               NODGLI = MCN(MNNDDL+I-1)

C              LE NOMBRE DE DL DU NOEUD I
               NDLP = MCN(MNNDDL+I) - NODGLI

C              COMPOSANTE 1 DE LA VITESSE
               V = VXYZPN(NODGLI+1)
               vitx( NCAS )%dptab( I ) = V
               VITMOD = V * V

C              COMPOSANTE 2 DE LA VITESSE
               V = VXYZPN(NODGLI+2)
               vity( NCAS )%dptab( I ) = V
               VITMOD = VITMOD + V * V

               IF( NDIM .EQ. 3 ) THEN
C                 COMPOSANTE 3 DE LA VITESSE
                  V = VXYZPN(NODGLI+3)
                  vitz(NCAS)%dptab( I ) = V
                  VITMOD = VITMOD + V * V
               ENDIF

C              MODULE DE LA VITESSE EN COMPOSANTE NDIM1
               VITMOD = SQRT( VITMOD )
               vitm( NCAS )%dptab(I) = VITMOD

               VITMOY = VITMOY + VITMOD
               IF( VITMOD .LT. VITMIN ) THEN
                  VITMIN  = VITMOD
                  NOEVMIN = I
                  NCAVMIN = NCAS
               ENDIF
               IF( VITMOD .GT. VITMAX ) THEN
                  VITMAX  = VITMOD
                  NOEVMAX = I
                  NCAVMAX = NCAS
               ENDIF

               IF( NDLP .EQ. NDIM1 ) THEN
C                 LA PRESSION EN CE NOEUD=SOMMET
                  V = VXYZPN(NODGLI+NDIM1)
C                 NO DE SOMMET DU NOEUD I
                  NST = MCN( MNNOSO-1+I )
                  pres(NCAS)%dptab(NST) = V
                  PREMOY = PREMOY + V
                  IF( V .LT. PREMIN ) THEN
                     PREMIN  = V
                     NOEPMIN = I
                     NCAPMIN = NCAS
                  ENDIF
                  IF( V .GT. PREMAX ) THEN
                     PREMAX  = V
                     NOEPMAX = I
                     NCAPMAX = NCAS
                  ENDIF
               ENDIF

C              TEMPERATURE AU NOEUD I
               IF( NTDLTE .GT. 0 ) THEN
                  V = TEMPER( I )
                  temp( NCAS )%dptab( I ) = V

                  TEMMOY = TEMMOY + V
                  IF( V .LT. TEMMIN ) THEN
                     TEMMIN  = V
                     NOETMIN = I
                     NCATMIN = NCAS
                  ENDIF
                  IF( V .GT. TEMMAX ) THEN
                     TEMMAX  = V
                     NOETMAX = I
                     NCATMAX = NCAS
                  ENDIF
               ENDIF

            ENDDO

 22      ENDDO

      ELSE IF( NUTYEL .EQ. 13 .OR. NUTYEL .EQ. 19 ) THEN

C        EF TRIANGLE ou TETRAEDRE BREZZI-FORTIN
C        --------------------------------------
C        NBNOVI = SOMMETS + NOEUDS BARYCENTRE DE CHAQUE EF EN VITESSE
C        ATTENTION: IMPLICITEMENT LES SOMMETS SONT NUMEROTES DE 1 A NBSOMT
C                   LES AUTRES NOEUDS BARYCENTRE SONT NUMEROTES NBSOMT+NO EF

         DO 24 NCAS = NCAS0, NCAS1, 1

C           TRANSFERT DE LA VITESSE-PRESSION AU TEMPS NOVVIPR PAR NOEUDS
C           FICHIER DU REPERTOIRE PROJET DANS LE TABLEAU MNVXYZP0
            TEMPS = RMCN( MNTIMES + NCAS - NCAST0 )

            CALL LIFIVIPRTE( KNOMOB, TEMPS,  NCAS,   NAVSTO, NBPASDT,
     %                       NDIM,   NBNOVI, NBSOMT, NBNOTE,
     %                       NTDLVP, VXYZPN, NTDLTE, TEMPER,
     %                       KNOMFIC, IERR )
            IF( IERR .NE. 0 ) GOTO 9997

C           NUMERO DU VECTEUR LU DANS VX VY VZ PR
            NBVVIPRLU = NBVVIPRLU + 1

            allocate( vitx( NCAS )%dptab(1:NBNOVI),STAT=IALLOC)
            IF( IALLOC .NE. 0 ) GOTO 9990

            allocate( vity( NCAS )%dptab(1:NBNOVI),STAT=IALLOC)
            IF( IALLOC .NE. 0 ) GOTO 9990

            allocate( vitz( NCAS )%dptab(1:NBNOVI),STAT=IALLOC)
            IF( IALLOC .NE. 0 ) GOTO 9990

            allocate( vitm( NCAS )%dptab(1:NBNOVI),STAT=IALLOC)
            IF( IALLOC .NE. 0 ) GOTO 9990

            allocate( pres( NCAS )%dptab(1:NBSOMT),STAT=IALLOC)
            IF( IALLOC .NE. 0 ) GOTO 9990

            IF( NTDLTE .GT. 0 ) THEN
               allocate( temp( NCAS )%dptab(1:NBSOMT),STAT=IALLOC)
               IF( IALLOC .NE. 0 ) GOTO 9990
            ENDIF

C           NBNOVI LE NOMBRE DE NOEUDS SUPPORTS des VITESSES + PRESSIONS
            DO I = 1, NBNOVI

C              LE NUMERO DU DL AVANT LE NOEUD I
               NODGLI = MCN(MNNDDL+I-1)

C              LE NOMBRE DE DL PRESSION AU NOEUD I
               NDLP = MCN(MNNDDL+I) - NODGLI - NDIM

C              COMPOSANTE 1 DE LA VITESSE
               V = VXYZPN(NODGLI+1)
               vitx( NCAS )%dptab( I ) = V
               VITMOD = V * V

C              COMPOSANTE 2 DE LA VITESSE
               V = VXYZPN(NODGLI+2)
               vity( NCAS )%dptab( I ) = V
               VITMOD = VITMOD + V * V

               IF( NDIM .EQ. 3 ) THEN
C                 COMPOSANTE 3 DE LA VITESSE
                  V = VXYZPN(NODGLI+3)
                  vitz(NCAS)%dptab( I ) = V
                  VITMOD = VITMOD + V * V
               ENDIF

C              MODULE DE LA VITESSE EN COMPOSANTE NDIM1
               VITMOD = SQRT( VITMOD )
               vitm( NCAS )%dptab( I ) = VITMOD

               VITMOY = VITMOY + VITMOD
               IF( VITMOD .LT. VITMIN ) THEN
                  VITMIN  = VITMOD
                  NOEVMIN = I
                  NCAVMIN = NCAS
               ENDIF
               IF( VITMOD .GT. VITMAX ) THEN
                  VITMAX  = VITMOD
                  NOEVMAX = I
                  NCAVMAX = NCAS
               ENDIF

               IF( NDLP .GT. 0 ) THEN

C                 LA PRESSION EN CE SOMMET DE NO AVANT LES BARYCENTRES
                  P = VXYZPN(NODGLI+NDIM1)
                  pres(NCAS)%dptab( I ) = P
                  PREMOY = PREMOY + P
                  IF( P .LT. PREMIN ) THEN
                     PREMIN  = P
                     NOEPMIN = I
                     NCAPMIN = NCAS
                  ENDIF
                  IF( P .GT. PREMAX ) THEN
                     PREMAX  = P
                     NOEPMAX = I
                     NCAPMAX = NCAS
                  ENDIF

C                 TEMPERATURE AU SOMMET I
                  IF( NTDLTE .GT. 0 ) THEN
                     V = TEMPER(I)
                     temp( NCAS )%dptab( I ) = V
                     TEMMOY = TEMMOY + V
                     IF( V .LT. TEMMIN ) THEN
                        TEMMIN  = V
                        NOETMIN = I
                        NCATMIN = NCAS
                     ENDIF
                     IF( V .GT. TEMMAX ) THEN
                        TEMMAX  = V
                        NOETMAX = I
                        NCATMAX = NCAS
                     ENDIF

                  ENDIF

               ENDIF
            ENDDO

 24      ENDDO

      ENDIF

C     LES VITESSES
C     CALCUL de la VITESSE et PRESSION MOYENNE des NBVVIPRLU VECTEURS
      VITMOY = VITMOY / ( NBNOVI * NBVVIPRLU )
      PREMOY = PREMOY / ( NBSOMT * NBVVIPRLU )
      IF( NUTYEL .EQ. 13 .OR. NUTYEL .EQ. 19 ) THEN
         NBT = NBSOMT
      ELSE
         NBT = NBNOVI
      ENDIF
      TEMMOY = TEMMOY / ( NBT * NBVVIPRLU )

      IF( NBVVIPRLU .NE. NBVVIPR ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'tfluide: ERREUR',NBVVIPRLU,
     %             ' VECTEURS LUS au lieu de',NBVVIPR,' VECTEURS'
         ELSE
            PRINT*,'tfluide: ERROR',NBVVIPRLU,
     %             ' READ VECTORS IN PLACE of',NBVVIPR,' VECTORS'
         ENDIF
         NBVVIPR = NBVVIPRLU
         IERR   = 1
      ENDIF

      VITMX = VITMAX

      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10100) NBNOVI, VITMOY
      ELSE
         WRITE(IMPRIM,20100) NBNOVI, VITMOY
      ENDIF
10100 FORMAT(' VITESSE MOYENNE en TEMPS et des',I9,' NOEUDS=',G13.6)
20100 FORMAT(' VELOCITY MEAN in TIME and',I9,' NODES=',G13.6)


ccc      IF( NOEVMIN .LE. NBNOEU ) THEN
      MN = MNXYZN + WYZNOE + 3 * NOEVMIN - 3
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10101) VITMIN,NOEVMIN,(RMCN(MN+K),K=0,2),NCAVMIN,
     %                       RMCN(MNTIMES-NCAST0+NCAVMIN)
      ELSE
         WRITE(IMPRIM,20101) VITMIN,NOEVMIN,(RMCN(MN+K),K=0,2),NCAVMIN,
     %                       RMCN(MNTIMES-NCAST0+NCAVMIN)
      ENDIF
ccc      ENDIF
10101 FORMAT(' VITESSE MINIMALE =',G13.6,' au NOEUD',I9,
     %' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' Cas',I5,' Temps',G14.6)
20101 FORMAT(' MINIMUM VELOCITY =',G13.6,' at Node',I9,
     %' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' Case',I5,' Time',G14.6)


ccc      IF( NOEVMAX .LE. NBNOEU ) THEN
      MN = MNXYZN + WYZNOE + 3 * NOEVMAX - 3
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10104) VITMAX,NOEVMAX,(RMCN(MN+K),K=0,2),NCAVMAX,
     %                       RMCN(MNTIMES-NCAST0+NCAVMAX)
      ELSE
         WRITE(IMPRIM,20104) VITMAX,NOEVMAX,(RMCN(MN+K),K=0,2),NCAVMAX,
     %                       RMCN(MNTIMES-NCAST0+NCAVMAX)
      ENDIF
ccc      ENDIF
10104 FORMAT(' VITESSE MAXIMALE =',G13.6,' au NOEUD',I9,
     %' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' Cas',I5,' Temps',G14.6)
20104 FORMAT(' MAXIMUM VELOCITY =',G13.6,' at Node',I9,
     %' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' Case',I5,' Time',G14.6)

C     MODIFICATION FAIBLE EN CAS DE VITMIN=VITMAX POUR LA VISEE
      IF( ABS(VITMAX-VITMIN) .LT. 1D-5*(ABS(VITMIN)+ABS(VITMAX)) ) THEN
         D = ABS(VITMAX) / 1000
         IF( D .EQ. 0D0 ) D = 1D0
         VITMAX = VITMAX + D
         VITMIN = VITMIN - D
      ENDIF

C     LES PRESSIONS
C     LE SOMMET DE PRESSION MINIMALE EST NOEPMIN
      MN = MNXYZN + WYZNOE - 3 + 3 * NOEPMIN
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,12010) NBSOMT,PREMOY,PREMIN,NOEPMIN,
     %                       (RMCN(MN+K),K=0,2),NCAPMIN,
     %                       RMCN(MNTIMES-NCAST0+NCAPMIN)
      ELSE
         WRITE(IMPRIM,22010) NBSOMT,PREMOY,PREMIN,NOEPMIN,
     %                       (RMCN(MN+K),K=0,2),NCAPMIN,
     %                       RMCN(MNTIMES-NCAST0+NCAPMIN)
      ENDIF
12010 FORMAT(/' PRESSION MOYENNE en TEMPS et aux',I9,' SOMMETS=',G13.6/
     %        ' PRESSION MINIMALE =',G13.6,
     %' au NOEUD',I7,' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' CAS',I5,' Temps',G14.6)
22010 FORMAT(/' MEAN    PRESSURE in TIME and at',I9,' VERTICES=',G13.6/
     %        ' MINIMUM PRESSURE =',G13.6,
     %' at Node',I7,' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' CASE',I5,' TIME',G14.6)
C
C     LE SOMMET DE PRESSION MAXIMALE EST NOEPMAX
      MN = MNXYZN + WYZNOE - 3 + 3 * NOEPMAX
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,12001) PREMAX,NOEPMAX,
     %                       (RMCN(MN+K),K=0,2),NCAPMAX,
     %                       RMCN(MNTIMES-NCAST0+NCAPMAX)
      ELSE
         WRITE(IMPRIM,22001) PREMAX,NOEPMAX,
     %                       (RMCN(MN+K),K=0,2),NCAPMAX,
     %                       RMCN(MNTIMES-NCAST0+NCAPMAX)
      ENDIF
12001 FORMAT(' PRESSION MAXIMALE =',G13.6,
     %' AU NOEUD',I7,' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' CAS',I5,' Temps',G14.6)
22001 FORMAT(' MAXIMUM PRESSURE =',G13.6,
     %' at Node',I7,' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' CASE',I5,' TIME',G14.6)
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,12002) PREMAX-PREMIN
      ELSE
         WRITE(IMPRIM,22002) PREMAX-PREMIN
      ENDIF
12002 FORMAT(' PRESSION MAXIMALE - MINIMALE =',G13.6 )
22002 FORMAT(' MAXIMUM - MINIMUM PRESSURE =',G13.6 )

C     MODIFICATION FAIBLE EN CAS DE PREMIN=PREMAX
      IF( ABS(PREMAX-PREMIN) .LT. 1D-5*ABS(PREMAX) ) THEN
         D = ABS(PREMAX) / 1000
         IF( D .EQ. 0D0 ) D = 1D0
         PREMAX = PREMAX + D
         PREMIN = PREMIN - D
      ENDIF

C     LES TEMPERATURES
      IF( NTDLTE .GT. 0 ) THEN

C        LE SOMMET DE TEMPERATURE MINIMALE EST NOETMIN
         MN = MNXYZN + WYZNOE - 3 + 3 * NOETMIN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,12011) NBSOMT,TEMMOY,TEMMIN,NOETMIN,
     %                         (RMCN(MN+K),K=0,2),NCATMIN,
     %                          RMCN(MNTIMES-NCAST0+NCATMIN)
         ELSE
            WRITE(IMPRIM,22011) NBSOMT,TEMMOY,TEMMIN,NOETMIN,
     %                         (RMCN(MN+K),K=0,2),NCATMIN,
     %                          RMCN(MNTIMES-NCAST0+NCATMIN)
         ENDIF
12011 FORMAT(/' TEMPERATURE MOYENNE en TEMPS et aux',I9,
     %' NOEUDS=',G13.6/' TEMPERATURE MINIMALE =',G13.6,
     %' au NOEUD',I7,' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' CAS',I5,' Temps',G14.6)
22011 FORMAT(/' MEAN    TEMPERATURE in TIME and at',
     %I9,' NODES=',G13.6/
     %        ' MINIMUM TEMPERATURE =',G13.6,
     %' at Node',I7,' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' CASE',I5,' TIME',G14.6)

C        LE SOMMET DE TEMPERATURE MAXIMALE EST NOETMAX
         MN = MNXYZN + WYZNOE - 3 + 3 * NOETMAX
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,12021) TEMMAX,NOETMAX,
     %                       (RMCN(MN+K),K=0,2),NCATMAX,
     %                       RMCN(MNTIMES-NCAST0+NCATMAX)
         ELSE
            WRITE(IMPRIM,22021) TEMMAX,NOETMAX,
     %                       (RMCN(MN+K),K=0,2),NCATMAX,
     %                       RMCN(MNTIMES-NCAST0+NCATMAX)
      ENDIF
12021 FORMAT(' TEMPERATURE MAXIMALE =',G13.6,
     %' AU NOEUD',I7,' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' CAS',I5,' Temps',G14.6)
22021 FORMAT(' MAXIMUM TEMPERATURE =',G13.6,
     %' at Node',I7,' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' CASE',I5,' TIME',G14.6)

C        MODIFICATION FAIBLE EN CAS DE TEMMIN=TEMMAX
         IF( ABS(TEMMAX-TEMMIN) .LT. 1D-5*ABS(TEMMAX) ) THEN
            D = ABS(TEMMAX) / 1000
            IF( D .EQ. 0D0 ) D = 1D0
            TEMMAX = TEMMAX + D
            TEMMIN = TEMMIN - D
         ENDIF

      ENDIF


C     ALLOCATION DYNAMIQUE D'UN TABLEAU VPAUX AUXILIAIRE ENCOMBRANT
C     -------------------------------------------------------------
C     NOFOVI : NUMERO DE LA FONCTION VITESSE_EXACTE(t,x,y,z,nc)
      NOFOVI = NOFOVITE()

C     NOFOPR : NUMERO DE LA FONCTION PRESSION_EXACTE(t,x,y,z)
      NOFOPR = NOFOPRES()

      IF( NOFOVI .GT. 0 .OR. NOFOPR .GT. 0 ) THEN
         IF( IALVPAUX .NE. 0 ) THEN
            MOVP = NBNOVI * NBVVIPR
            PRINT*
ccc            PRINT*, 'ALLOCATION DEMAND  of',MOVP,
ccc     %              ' DOUBLE PRECISION of VPAUX(',NBNOVI,',',NBVVIPR,')'

            ALLOCATE ( VPAUX(1:NBNOVI,NCAS0:NCAS1), STAT=IALVPAUX )
            IF( IALVPAUX .NE. 0 ) THEN
               IF( LANGAG .EQ. 0 ) THEN
                  PRINT*,'ERREUR en ALLOCATION de',MOVP,
     %                   ' DOUBLE PRECISION de VPAUX(',NBNOVI,',',
     %                    NCAS0,':',NCAS1,')'
               ELSE
                  PRINT*,'ALLOCATION ERROR   of',MOVP,
     %                   ' DOUBLE PRECISION of VPAUX(',NBNOVI,',',
     %                    NCAS0,':',NCAS1,')'
               ENDIF
               IALR = IALVPAUX
               GOTO 9997
            ENDIF

            IF( LANGAG .EQ. 0 ) THEN
               PRINT*, 'CORRECTE ALLOCATION de',MOVP,
     %                 ' DOUBLE PRECISION de VPAUX(',NBNOVI,',',
     %                   NCAS0,':',NCAS1,')'
            ELSE
               PRINT*, 'ALLOCATION CORRECT of',MOVP,
     %                 ' DOUBLE PRECISION of VPAUX(',NBNOVI,',',
     %                   NCAS0,':',NCAS1,')'
            ENDIF
         ENDIF
      ENDIF
      PRINT*


C     ******************************************************************
C     ******************************************************************
C     TRACES de la VITESSE ou de la PRESSION ou TEMPERATURE ou AFFICHAGE
C     ******************************************************************
C     ******************************************************************
 55   CALL RECTEF( NRHIST )
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO

      IF( NDIM .LE. 2 ) THEN
         CALL LIMTCL( 'vitepr2d', NMTCL )
      ELSE
         CALL LIMTCL( 'vitepr3d', NMTCL )
      ENDIF

      IF( NMTCL .LE. 0  ) GOTO  9997
      IF( NMTCL .EQ. 90 ) GOTO 19000
      IF( NMTCL .EQ. 91 ) GOTO 19100
      IF( NMTCL .EQ. 92 ) GOTO 19200
      GOTO ( 2000, 200, 3000, 4000, 5000, 6000, 7000, 55, 9000, 10000,
     %       55, 12000, 13000 ),NMTCL
      GOTO 55


C     ******************************************************************
C     TRACE DE LA PRESSION PORTEE PAR LES SOMMETS EN INTERPOLATION P1
C     ******************************************************************
 2000 PRMIN = REAL( PREMIN )
      PRMAX = REAL( PREMAX )
C
C     OPTIONS DU TRACE DE LA PRESSION NCAS0 A NCAS1
C     PROVENANCE PRESSION D'UN FLUIDE
      MODECO = 3
C
C     MODIFICATION DU NUMERO DE NOEUDS DES SOMMETS DANS LE TABLEAU
C     DES FACES POUR LE TETRAEDRE DE TAYLOR-HOOD
      IF( NUTYEL .EQ. 20 ) THEN
         CALL ARETFRP2P1( KNOMOB, NBNOEU, MCN(MNNOSO), IERR )
         IF( IERR .NE. 0 ) GOTO 55
      ENDIF
C
C     LES OPTIONS DE TRACE DE LA PRESSION
C     ******************************************************************
 2005 IF( NDIM .LE. 2 ) THEN
         CALL LIMTCL( 'tracpr2d', NMTCL0 )
      ELSE
         CALL LIMTCL( 'tracpr3d', NMTCL0 )
      ENDIF

      IF( NMTCL0 .LE. 0 ) THEN
C        RETOUR AU NUMERO DE NOEUDS DES SOMMETS DANS LE TABLEAU DES FACES
C        POUR LE TETRAEDRE DE TAYLOR-HOOD
         IF( NUTYEL .EQ. 20 ) THEN
            CALL ARETFRP1P2( KNOMOB, NBSOMT, MCN(MNSONO) )
         ENDIF
         GOTO 55
      ENDIF

      GOTO( 2010, 2020, 2030, 2040, 2050,
     %      2060, 2070, 2080, 2090       ) ,NMTCL0
      GOTO 2005
C
C     TRACE DES ISO-PRESSIONS (LIGNES en 2D, SURFACES en 3D)
C     ======================================================
 2010 CALL TRISOT( NDIM,   KNOMOB,   MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1, NDPGSTP1,
     %             NCAS0,  NCAS1,    NBSOMT,
     %             1,      PRESSION(1,NCAS0), pres,
     %             PRMIN,  NOEPMIN, NCAPMIN, PRMAX,  NOEPMAX, NCAPMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 2005
C
C     TRACE DES ZONES DE COULEURS ISOTHERMES
C     ======================================
 2020 CALL TRZONT( 0,      NDIM,    KNOMOB,   MODECO,
     %             NBTYEL, MNNPEFP1,MNXYZNP1, MNXYZNP1, NDPGSTP1,
     %             NCAS0,  NCAS1,   NBSOMT,
     %             1,      PRESSION, pres,
     %             PRMIN,  NOEPMIN, NCAPMIN, PRMAX,  NOEPMAX, NCAPMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 2005
C
C     TRACE DES ZONES DE COULEURS ISOTHERMES PAR SECTIONS X ou Y ou Z=CTE
C     ===================================================================
 2030  IF( NDIM .EQ. 2 ) GOTO 2020
C      LA PRESSION EN 3D
       CALL TRPLSE( 0,      KNOMOB,   MODECO,
     %              NBTYEL, MNNPEFP1, MNXYZNP1, NDPGSTP1,
     %              NCAS0,  NCAS1,    NBSOMT, 
     %              1,      PRESSION, pres,
     %              PRMIN,  NOEPMIN, NCAPMIN, PRMAX,  NOEPMAX, NCAPMAX,
     %              RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 2005
C
C     TRACE DES PROFILS DE COULEURS PAR SECTIONS X ou Y ou Z=CTE
C     ==========================================================
 2040 IF( NDIM .EQ. 2 ) GOTO 2090
C      LA PRESSION EN 3D
       CALL TRPLSE( 1,      KNOMOB,   MODECO,
     %              NBTYEL, MNNPEFP1, MNXYZNP1, NDPGSTP1,
     %              NCAS0,  NCAS1,    NBSOMT,
     %              1,      PRESSION, pres,
     %              PRMIN,  NOEPMIN, NCAPMIN, PRMAX,  NOEPMAX, NCAPMAX,
     %              RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 2005
C
C     TRACE DE LA PRESSION LE LONG D'UNE DROITE DEFINIE PAR 2 POINTS
C     ==============================================================
 2050 IF( NDIM .EQ. 2 ) GOTO 2020
C     LA PRESSION EN 3D
      CALL TRLLDR( 0,      NDIM,     KNOMOB, MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1, NDPGSTP1,
     %             NCAS0,  NCAS1,    NBSOMT,
     %             1,      PRESSION, pres,
     %             PRMIN,  NOEPMIN,  NCAPMIN, PRMAX,  NOEPMAX, NCAPMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 2005
C
C     TRACE DE LA PRESSION SUR UNE OU PLUSIEURS DES SURFACES DE L'OBJET
C     =================================================================
 2060 IF( NDIM .EQ. 2 ) GOTO 2020
C     LA PRESSION P1 EN 3D
      CALL TRSO1SO( INTERPP, 0, MODECO, KNOMOB, NTLXOB, MNDFOB,
     %              NBSOMT, NCAS0, NCAS1,
     %              1,      PRESSION, pres,
     %              PRMIN,  PRMAX, RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 2005
C
C     TRACE DE LA PRESSION SUR TOUTES LES SURFACES FRONTALIERES
C     =========================================================
 2070 IF( NDIM .EQ. 2 ) GOTO 2020
C     LA PRESSION P1 EN 3D
      CALL TRSOSF( 0,      KNOMOB, MODECO, MNXYZN,
     %             NBSOMT, NCAS0,  NCAS1,
     %             1,      PRESSION, pres,
     %             PRMIN,  PRMAX,  RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 2005
C
C     CALCUL DE L'INTEGRALE DE LA PRESSION SUR DES SURFACES DE L'OBJET
C     ================================================================
 2080 IF( NDIM .EQ. 2 ) GOTO 2020
      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,
     %   'INTEGRALE DE LA PRESSION SUR LES SURFACES DE L''OBJET'
         PRINT*,
     %   '----------------------------------------------------'
      ELSE
         PRINT*,
     %   'PRESSURE INTEGRAL ON OBJECT SURFACES'
         PRINT*,
     %   '------------------------------------'
      ENDIF
      IF( LANGAG .EQ. 0 ) THEN
         NMSOLU = 'PRESSION '
      ELSE
         NMSOLU = 'PRESSURE '
      ENDIF
      CALL INSOSF( 0,      INTERPP, NTLXOB, MNDFOB,
     %             NCAS0,  NCAS1,   pres, NMSOLU,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 2005
C
C     TRACE EN 2D SURFACE(X,Y,PRESSION(X,Y))
C     ======================================
 2090 IF( NDIM .NE. 2 ) GOTO 2005
      CALL TRZTXY(      0, NDIM,     KNOMOB, MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1,
     %             NCAS0,  NCAS1,    NBSOMT,
     %             1,      PRESSION, pres,
     %             PRMIN,  NOEPMIN, NCAPMIN, PRMAX,  NOEPMAX, NCAPMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 2005


C     ********************************************************
C     TRACE DE LA VITESSE SELON des FLECHES
C     ********************************************************
 200  VITMX = VITMAX

C     AUCUN ITEM VISIBLE
      CALL ITEMS0

C     SI LA SOURIS EST DISPONIBLE EN INTERACTIF ELLE EST ACTIVEE
C     POUR MODIFIER LE CADRE DES TRACES
      IF( INTERA .GE. 2 ) LORBITE=1

C     LA PREPARATION DU TRACE: VISEE PAR DEFAUT
      CALL VISEE0
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO

C     LECTURE DES OPTIONS DE TRACE DES VITESSES
C     =========================================
 202  CALL LIMTCL( 'trflevit', NMTCL )
      IF( NMTCL .LE.  0 ) GOTO 55
      IF( NMTCL .EQ. 50 ) GOTO 225
      IF( NMTCL .EQ. 90 ) GOTO 300
      GOTO( 202, 212, 213, 214, 215, 216, 217, 218, 219 ), NMTCL

C     FLECHE MAXIMALE TRACEE EN CM
C     ----------------------------
 212  NOPT   = 1
      NCVALS = 0
      CALL INVITE( 66 )
      CALL LIRRSP(NCVALS,CMFLEC)
      IF( NCVALS .EQ. -1 ) GOTO 200
      IF( CMFLEC .LT. 0. ) CMFLEC = -CMFLEC
      VITMX = VITMAX
      GOTO 202

C     1CM  VAUT EN NORME DE LA VITESSE
C     --------------------------------
 213  NOPT   = 2
      NCVALS = 0
      CALL INVITE( 132 )
      CALL LIRRDP(NCVALS,VITMX)
      IF( NCVALS .EQ. -1 ) GOTO 200
      IF( VITMX .LT. 0. ) VITMX = -VITMX
      IF( VITMX .EQ. 0. ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VITESSE MAX NULLE pour 1CM. A MODIFIER'
         ELSE
            KERR(1) = 'For 1CM MAXIMUM VELOCITY NULL. MODIFY it'
         ENDIF
         CALL LEREUR
         GOTO 202
      ENDIF
      GOTO 202

C     ECHELLE PRECEDENTE
C     ------------------
 214  IF( CMVITE .LE. 0. ) THEN
         WRITE(KERR(MXLGER)(1:13),'(g13.6)') CMVITE
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: 1CM/VITESSE='// KERR(MXLGER)(1:13)
            KERR(2) = 'ECHELLE INTERDITE :'
     %           // ' LONGUEUR FLECHE MAXIMALE FORCEE A 1 CM'
         ELSE
            KERR(1) = 'ERROR: 1CM/VELOCITY='// KERR(MXLGER)(1:13)
            KERR(2) = 'VALUE FORBIDDEN :'
     %           // ' MAXIMUM ARROW LENGTH FORCED at 1 CM'
         ENDIF
         CALL LEREUR
         NOPT   = 1
         CMFLEC = 2.5
      ELSE
         NOPT   = 3
      ENDIF
      GOTO 202
C
C     COULEUR des ARETES du MAILLAGE
C     ------------------------------
 215  CALL LIMTCL( 'couleur0' , I )
      IF( I .EQ. -1 ) THEN
         GOTO 202
      ELSE IF( I .EQ. -2 ) THEN
         NCOUAF = -2
      ELSE IF( I .EQ. 0 ) THEN
C        COULEUR NOIRE
         NCOUAF = 0
      ELSE
C        COULEUR RESERVEE
         NCOUAF = N1COEL + I
      ENDIF
      GOTO 202
C
C     TYPE du TRAIT des ARETES du MAILLAGE
C     ------------------------------------
 216  CALL LIMTCL( 'typtrait' , I )
      IF( I .EQ. -1 ) GOTO 202
      NTLAFR = I
      GOTO 202
C
C     COULEUR des FLECHES
C     -------------------
 217  CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 202
      ELSE IF( I .EQ. 0 ) THEN
C        COULEUR NOIRE
         NCOUFL = 0
      ELSE
         NCOUFL = N1COEL + I
      ENDIF
      GOTO 202
C
C     NOMBRE D'EPAISSEURS DE TRAIT D'UNE FLECHE
C     -----------------------------------------
 218  NCVALS = 0
      CALL INVITE( 78 )
      CALL LIRENT( NCVALS, NEPFLE )
      IF( NCVALS .EQ. -1 ) GOTO 202
      IF( NEPFLE .LT. 0 ) NEPFLE=0
      IF( NEPFLE .GT. 5 ) NEPFLE=5
      GOTO 202
C
C     PAS DU NOMBRE DE FLECHES DE VITESSE A TRACER (1/5, 1/2, ... )
C     -------------------------------------------------------------
 219  NCVALS = 0
      CALL INVITE( 156 )
      CALL LIRENT( NCVALS, NPAFLE )
      IF( NCVALS .EQ. -1 ) GOTO 202
      IF( NPAFLE .LE. 0     ) NPAFLE=1
      IF( NPAFLE .GT. 10000 ) NPAFLE=10000
      GOTO 202
C
C     EFFACER LE TRACE ACTUEL
C     -----------------------
 225  CALL EFFACE
C     PLUS D'ITEMS VISIBLES
      CALL ITEMS0
      CALL TRAXES
      GOTO 202
C
C     OPTION DE VISEE DU TRACE DES FLECHES DES VITESSES
C     =================================================
 300  IF( NDIM  .EQ. 2 ) THEN
         CALL VISEE(NMTCL)
         IF( NMTCL .LT. 0 ) GOTO 202
      ELSE
C        TRACE TOUTES LES VITESSES OU CELLES DES EF SECTIONNES PAR UN PLAN
         CALL LIMTCL( 'vitetype', NMTCLV )
         IF( NMTCLV .LE. 0 ) GOTO 202
      ENDIF
C
C     EXECUTION DU TRACE DES VITESSES
      IF( NOPT .NE. 1 .OR. CMFLEC .GT. 0. ) GOTO 310
C
C     ERREUR
 305  WRITE(KERR(MXLGER-3)(1:13),'(I13)')   NOPT
      WRITE(KERR(MXLGER-2)(1:13),'(g13.6)') CMFLEC
      WRITE(KERR(MXLGER-1)(1:13),'(g13.6)') VITMAX
      WRITE(KERR(MXLGER)(1:13),'(g13.6)')   CMVITE
      NBLGRC(NRERR) = 4
      KERR(1) = ' ERREUR :  OPTION            '
     %        // KERR(MXLGER-3)(1:13)
      KERR(2) = ' FLECHE MAXIMALE EN CM       '
     %        // KERR(MXLGER-2)(1:13)
      KERR(3) = ' VITESSE D''UN CM DE TRACE'
     %        // KERR(MXLGER-1)(1:13)
      KERR(4) = ' CM / VITESSE             '
     %        // KERR(MXLGER)(1:13)
      CALL LEREUR
      GOTO 202
C
 310  IF( NOPT .EQ. 2 .AND. VITMX  .LE. 0. ) GOTO 305
      IF( NOPT .EQ. 3 .AND. CMVITE .LE. 0. ) GOTO 305
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10310) NCAS0, NCAS1
      ELSE
         WRITE(IMPRIM,20310) NCAS0, NCAS1
      ENDIF
10310 FORMAT(' VECTEURS VITESSE VISUALISE =',I7,' A',I7/)
20310 FORMAT(' DRAWN VELOCITY VECTOR NUMBERS =',I7,' to',I7/)

C     MISE A JOUR DE L'ECHELLE
      IF( NOPT .EQ. 1 ) THEN
         CMVITE = REAL( CMFLEC / VITMX )
      ELSE IF( NOPT .EQ. 2 ) THEN
         CMVITE = REAL( 1D0 / VITMX )
      ENDIF

C     TRACE EFFECTIF DES FLECHES DES VITESSES
C     =======================================
      IF( NDIM .EQ. 2 ) THEN

         IF( LANGAG .EQ. 0 ) THEN
            PRINT*, 'Les TEMPS des VECTEURS VITESSES'
         ELSE
            PRINT*, 'The TIMES of VELOCITY VECTORS'
         ENDIF
         WRITE(IMPRIM,10320) (K,RMCN(MNTIMES-NCAST0+K),K=1,NBVVIPR)
10320    FORMAT(5(I6,':',G13.6,'   '))
C
         CALL TRVI2D( KNOMOB, NBNOVI, MCN(MNXYZN+WYZNOE), MNNPEF,
     %                NCAS0,  NCAS1,  vitx,  vity,
     %                VITMIN, VITMOY, VITMAX,
     %                RMCN(MNTIMES-NCAST0+NCAS0), CMVITE )
         GOTO 300

      ELSE

         IF( NMTCLV .EQ. 1 ) THEN

C           TRACE DES VITESSES DES EF SECTIONNES PAR UN PLAN POUR LE NCAS1
            TEMPS = RMCN(MNTIMES-NCAST0+NCAS1)
            CALL TRVI3DPLEF(KNOMOB,  NBSOMT, NBNOVI, MCN(MNXYZN+WYZNOE),
     %                      NBTYEL,  MNNPEF, CMVITE,
     %                      NCAS0,   NCAS1,  NCAS1,
     %                      vitx,    vity,   vitz,
     %                      RMCN(MNTIMES-NCAST0+NCAS0) )
            GOTO 300

         ELSE IF( NMTCLV .EQ. 2 ) THEN

C           TRACE DES VITESSES EN TOUS LES POINTS
            CALL TRVI3D( KNOMOB, NBNOVI, MCN(MNXYZN+WYZNOE),
     %                   0D0,    VITMAX, CMVITE,
     %                   NCAS0,  NCAS1,  vitx,  vity,  vitz,
     %                   VITMOY, VITMIN, VITMAX,
     %                   RMCN(MNTIMES-NCAST0+NCAS0) )
            GOTO 300

         ELSE IF( NMTCLV .EQ. 3 ) THEN

C           Entree du SEUIL MINIMAL de la NORME d''une VITESSE a tracer
            CALL INVITE( 150 )
            NCVALS = 0
            CALL LIRRDP( NCVALS, SEMINVIT )
            IF( NCVALS .LT. 0 ) THEN
C              ABANDON DE LA LECTURE DES DONNEES
               SEMINVIT = 0D0
               GOTO 300
            ENDIF

C           Entree du SEUIL MAXIMAL de la NORME d''une VITESSE a tracer
            CALL INVITE( 153 )
            NCVALS = 0
            CALL LIRRDP( NCVALS, SEMAXVIT )
            IF( NCVALS .LT. 0 ) THEN
C              ABANDON DE LA LECTURE DES DONNEES
               SEMAXVIT = 1D100
               GOTO 300
            ENDIF

            IF( SEMINVIT .LT. 0D0    ) SEMINVIT = 0D0
            IF( SEMAXVIT .GT. VITMAX ) SEMAXVIT = VITMAX
            IF( SEMINVIT .GT. SEMAXVIT ) THEN
               D        = SEMINVIT
               SEMINVIT = SEMAXVIT
               SEMAXVIT = D
            ENDIF

            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10330) SEMINVIT, SEMAXVIT
            ELSE
               WRITE(IMPRIM,20330) SEMINVIT, SEMAXVIT
            ENDIF
10330 FORMAT('SEUIL MINIMAL de la NORME d''une VITESSE a tracer',G13.6/
     %       'SEUIL MAXIMAL de la NORME d''une VITESSE a tracer',G13.6)
20330 FORMAT('MINIMUM of the VELOCITY MAGNITUDE to draw',G13.6/
     %       'MAXIMUM of the VELOCITY MAGNITUDE to draw',G13.6)
            GOTO 300

         ELSE IF( NMTCLV .EQ. 4 ) THEN

C           TRACE DES VITESSES EN TOUS LES POINTS AVEC SEUIL MIN A TRACER
            IF( SEMINVIT . LT. 0D0    ) SEMINVIT = 0D0
            IF( SEMAXVIT . GT. VITMAX ) SEMAXVIT = VITMAX
            CALL TRVI3D( KNOMOB, NBNOVI, MCN(MNXYZN+WYZNOE),
     %                   SEMINVIT, SEMAXVIT, CMVITE,
     %                   NCAS0,    NCAS1, vitx, vity, vitz,
     %                   VITMOY, VITMIN, VITMAX,
     %                   RMCN(MNTIMES-NCAST0+NCAS0) )
            GOTO 300
         ENDIF

      ENDIF
      GOTO 300


C     *****************************************************************
C     TRACE D'UNE COMPOSANTE NOCOMP ou du MODULE DE LA VITESSE
C     *****************************************************************
C     MIN MAX de la COMPOSANTE NOCOMP de la Vitesse aux NOEUDS
 3000 VIMIN = REAL( VITMIN )
      VIMAX = REAL( VITMAX )

C     TRACE D'UNE COMPOSANTE NOCOMP de la VITESSE des TEMPS NCAS0 a NCAS1
      MODECO = 6

C     PAR DEFAUT TRACE DU MODULE DE LA VITESSE
      NOCOMP = NDIM1

 350  CALL LIMTCL( 'tracmovi', NMTCL0 )
      IF( NMTCL0 .LT. 0 ) THEN
         GOTO 55
      ENDIF

      IF( NMTCL0 .EQ. 0 ) THEN
C        CHOIX du No de la COMPOSANTE ou MODULE de la VITESSE
         CALL LIMTCL( 'noxyzmod', NOCOMP )
         IF( NOCOMP .LE. 0 ) THEN
            GOTO 55
         ENDIF
         IF( NOCOMP .EQ. 4 ) THEN
C           NECESSAIRE POUR NDIM=2
            NOCOMP = NDIM1
         ENDIF
C        MODULE IMPOSE DE LA VITESSE
         NOCOMP = NDIM1
         GOTO 350
      ENDIF

      GOTO( 410, 420, 430, 440, 450, 460, 350, 350, 490, 55 ),NMTCL0


C     TRACE DES ISO-MODULE VITESSES (LIGNES en 2D, SURFACES en 3D)
C     ============================================================
 410  CALL TRISOT( NDIM,   KNOMOB,  MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             1, VITXYZ(1,1,NOCOMP), vitm,
     %             VIMIN,  NOEVMIN, NCAVMIN, VIMAX,  NOEVMAX, NCAVMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 350

C     TRACE DES ZONES DE COULEURS ISOMODULE VITESSE
C     =============================================
 420  CALL TRZONT( 0,      NDIM,    KNOMOB, MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN, MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             1, VITXYZ(1,1,NOCOMP), vitm,
     %             VIMIN,  NOEVMIN, NCAVMIN, VIMAX, NOEVMAX, NCAVMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 350

C     TRACE DES ZONES DE COULEURS ISOMODULES PAR SECTIONS X ou Y ou Z=CTE
C     ===================================================================
 430  IF( NDIM .EQ. 2 ) GOTO 420
C     LE MODULE DES VITESSES EN 3D
      CALL TRPLSE( 0,      KNOMOB,  MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             1, VITXYZ(1,1,NOCOMP), vitm,
     %             VIMIN,  NOEVMIN, NCAVMIN, VIMAX, NOEVMAX, NCAVMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 350

C     TRACE DES PROFILS DE COULEURS MODULE VITESSE PAR SECTIONS X ou Y ou Z=CTE
C     =========================================================================
 440  IF( NDIM .EQ. 2 ) GOTO 490
C     LE MODULE DES VITESSES EN 3D
      CALL TRPLSE( 1,      KNOMOB,  MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             1, VITXYZ(1,1,NOCOMP), vitm,
     %             VIMIN,  NOEVMIN, NCAVMIN, VIMAX,  NOEVMAX, NCAVMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 350

C     TRACE DU MODULE VITESSE LE LONG D'UNE DROITE DEFINIE PAR 2 POINTS
C     =================================================================
 450  IF( NDIM .EQ. 2 ) GOTO 420
C     LE MODULE VITESSE EN 3D
      CALL TRLLDR( 0,      NDIM,    KNOMOB, MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             1, VITXYZ(1,1,NOCOMP), vitm,
     %             VIMIN,  NOEVMIN, NCAVMIN, VIMAX,  NOEVMAX, NCAVMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 350

C     TRACE DU MODULE DE LA VITESSE SUR DES SURFACES DE L'OBJET
C     =========================================================
 460  IF( NDIM .EQ. 2 ) GOTO 420
C     LE MODULE DE LA VITESSE P2 POUR TAYLOR-YOUNG et
C     P1+BULLE POUR BREZZI-FORTIN
      CALL TRSO1SO( INTERPV,  0,  MODECO, KNOMOB, NTLXOB, MNDFOB,
     %              NBNOVI, NCAS0, NCAS1,
     %              1, VITXYZ(1,1,NOCOMP), vitm,
     %              VIMIN,  VIMAX, RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 350

C     TRACE EN 2D SURFACE(X,Y,MODULE VITESSE(X,Y))
C     ============================================
 490  IF( NDIM .NE. 2 ) GOTO 350
      CALL TRZTXY(      0, NDIM,    KNOMOB, MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             1, VITXYZ(1,1,NOCOMP), vitm,
     %             VIMIN,  NOEVMIN, NCAVMIN, VIMAX, NOEVMAX, NCAVMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 350


C     **********************************************************************
C     TRACE de |VITESSE|Moyenne  |VITESSE|Maximum  PRESSION Max-Min=f(Temps)
C     EN TOUS LES TEMPS ou LES VECTEURS VITESSE-PRESSION ont ete CALCULES
C     **********************************************************************
4000  CALL TRMXVIPR( KNOMOB, IERR, DCPU )
C     DEFINITION DE LA VISEE POUR LA TOTALITE DE L'OBJET
      CALL VISEE0
      GOTO 55


C     *****************************************************************
C     TRACES de la TEMPERATURE du FLUIDE BOUSSINESQ
C     *****************************************************************
 5000 IF( NTDLTE .LE. 0 ) GOTO 55
C     MIN MAX de la TEMPERATURE
      TEMIN = REAL( TEMMIN )
      TEMAX = REAL( TEMMAX )

C     TRACE de la TEMPERATURE des TEMPS NCAS0 a NCAS1
      MODECO = 1

 5001 IF( NDIM .EQ. 2 ) THEN
         CALL LIMTCL( 'tractem2', NMTCL0 )
      ELSE
         CALL LIMTCL( 'tractem3', NMTCL0 )
      ENDIF
      IF( NMTCL0 .LT. 0 ) THEN
         GOTO 55
      ENDIF

      GOTO( 5010, 5020, 5030, 5040, 5050, 5001, 5070, 5001, 5001),NMTCL0

C     TRACE DES ISOTHERMES par des LIGNES en 2D ou des SURFACES en 3D
C     ===============================================================
 5010 CALL TRISOT( NDIM,   KNOMOB,  MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             1,      TEMPER,  temp,
     %             TEMIN,  NOETMIN, NCATMIN, TEMAX,  NOETMAX, NCATMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 5001

C     TRACE DES ZONES DE COULEURS DES TEMPERATURES
C     =============================================
 5020 CALL TRZONT( 0,      NDIM,    KNOMOB, MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN, MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             1,      TEMPER,  temp,
     %             TEMIN,  NOETMIN, NCATMIN, TEMAX, NOETMAX, NCATMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 5001

C     TRACE DES ZONES DE COULEURS TEMPERATURES PAR SECTIONS X ou Y ou Z=CTE
C     =====================================================================
 5030 IF( NDIM .EQ. 2 ) GOTO 5020
      CALL TRPLSE( 0,      KNOMOB,  MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             1,      TEMPER,  temp,
     %             TEMIN,  NOETMIN, NCATMIN, TEMAX, NOETMAX, NCATMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 5001

C     TRACE DES PROFILS DE COULEURS TEMPERATURES PAR SECTIONS X ou Y ou Z=CTE
C     =======================================================================
 5040 IF( NDIM .EQ. 2 ) GOTO 5020
      CALL TRPLSE( 1,      KNOMOB,  MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             1,      TEMPER,  temp,
     %             TEMIN,  NOETMIN, NCATMIN, TEMAX, NOETMAX, NCATMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 5001

C     TRACE DU MODULE VITESSE LE LONG D'UNE DROITE DEFINIE PAR 2 POINTS
C     =================================================================
 5050 IF( NDIM .EQ. 2 ) GOTO 5020
      CALL TRLLDR( 0,      NDIM,    KNOMOB, MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             1,      TEMPER,  temp,
     %             TEMIN,  NOETMIN, NCATMIN, TEMAX, NOETMAX, NCATMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 5001

C     TRACE EN 2D SURFACE(X,Y,MODULE VITESSE(X,Y))
C     ============================================
 5070  IF( NDIM .NE. 2 ) GOTO 5001
      CALL TRZTXY(      0, NDIM,    KNOMOB, MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             1,      TEMPER,  temp,
     %             TEMIN,  NOETMIN, NCATMIN, TEMAX, NOETMAX, NCATMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 5001


C     *****************************************************************
C     TRACES des FLUX de la TEMPERATURE du FLUIDE BOUSSINESQ
C     *****************************************************************
 6000 IF( NTDLTE .LE. 0 ) GOTO 55

C     CONSTRUCTION DU TMS VECTEUR"TEMPERATURE DU CAS NCAS1
      NDSM = 1
      CALL LXTSOU( NTLXOB, 'VECTEUR"TEMPERATURE', NTVECT, MNVECT )
      IF( NTVECT .GT. 0 ) THEN
C        LE VECTEUR DECLARE EST DETRUIT
         CALL LXTSDS( NTLXOB, 'VECTEUR"TEMPERATURE' )
      ENDIF

C     RECUPERATION DES TABLEAUX POBA NECESSAIRES A LA
C     CONSTRUCTION DES TABLEAUX ELEMENTAIRES DES EF TAYLOR-HOOD P2
      CALL TAPOBA( NBTYEL, MNNPEF, ETTAEL,
     %             MNTPOB, NBDLMX, MOTAUX, NBTTEF, NOAXIS, I, IERR )
      IF( IERR   .NE. 0 ) GOTO 55
      IF( MOTAUX .GT. 0 .AND. MNTAUX .LE. 0 ) THEN
         CALL TNMCDC( 'REEL2', MOTAUX, MNTAUX )
         IF( MNTAUX .EQ. 0 ) GOTO 55
      ENDIF

C     MOTAEL: NOMBRE DE REEL2 DE LA DECLARATION DU TABLEAU TAEL
      MOTAEL = NBDLMX * (NBDLMX+1) / 2 + NBDLMX * NDSM
      IF( MNTAEL .LE. 0 ) THEN
         CALL TNMCDC( 'REEL2', MOTAEL, MNTAEL )
         IF( MNTAEL .EQ. 0 ) GOTO 6099
      ENDIF

C     LE TABLEAU DES NBCOOR COORDONNEES DES NBDLMX POINTS DE L'EF MAXIMAL
      IF( MNX .LE. 0 ) THEN
         CALL TNMCDC( 'REEL', NBDLMX*3, MNX )
      ENDIF

C     RETROUVER LES ADRESSES MCN DES DONNEES THERMIQUES
C     DES   SV "OBJETS INTERNES"    DE L'OBJET
C     DES PLS  "OBJETS AUX LIMITES" DE L'OBJET
C     RETRIEVE THE THERMIC DATA of SURFACES or VOLUMES
C              THE BOUNDARY CONDITIONS on POINTS LINES (SURFACES)
      NBJEUX = 1
      CALL THEDON( NUMIOB, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             NBJEUX, MNDOEL,
     %             IEMAST, IECHMA, IECOND, IEDILA, IEVIFL, IECOET,
     %             IESOIN, IECONT, IEECHA, IESOCL, IESOPO,
     %             IETEIN, IEVIIN, IEVIANT,IECOBO,
     %             IERR )
      IF( IECOND .NE. NBOBIN ) THEN
          NBLGRC(NRERR) = 3
          IF( LANGAG .EQ. 0 ) THEN
             IF( NDIM .EQ. 1 ) THEN
                KERR(1) = 'UNE LIGNE SANS CONDUCTIVITE'
                KERR(3) = 'DONNER UNE CONDUCTIVITE A CETTE LIGNE'
             ELSE IF( NDIM .EQ. 2 ) THEN
                KERR(1) = 'UNE SURFACE SANS CONDUCTIVITE'
                KERR(3) = 'DONNER UNE CONDUCTIVITE A CETTE SURFACE'
             ELSE
                KERR(1) = 'UN VOLUME SANS CONDUCTIVITE'
                KERR(3) = 'DONNER UNE CONDUCTIVITE A CE VOLUME'
             ENDIF
             KERR(2) = '=> PROBLEME THERMIQUE SANS SOLUTION'
             KERR(3) = 'DONNER UNE CONDUCTIVITE A CETTE SURFACE'
          ELSE
             IF( NDIM .EQ. 1 ) THEN
                KERR(1) = 'A LINE WITHOUT CONDUCTIVITY'
                KERR(3) = 'GIVE A CONDUCTIVITY at this LINE'
             ELSE IF( NDIM .EQ. 2 ) THEN
                KERR(1) = 'A SURFACE WITHOUT CONDUCTIVITY'
                KERR(3) = 'GIVE A CONDUCTIVITY at this SURFACE'
             ELSE
                KERR(1) = 'A VOLUME WITHOUT CONDUCTIVITY'
                KERR(3) = 'GIVE A CONDUCTIVITY at this VOLUME'
             ENDIF
             KERR(2) = '=> HEAT PROBLEM WITHOUT SOLUTION'
          ENDIF
          CALL LEREUR
          IERR = 11
      ENDIF
      IF( IERR .NE. 0 ) GOTO 6099

C     DECLARATION DU TMS 'VECTEUR"TEMPERATURE'
      L = WECTEU + NTDLTE * MOREE2 + NDSM
      CALL LXTNDC( NTLXOB, 'VECTEUR"TEMPERATURE', 'MOTS', L )
      CALL LXTSOU( NTLXOB, 'VECTEUR"TEMPERATURE', NTVECT, MNVECT )
      IF( MNVECT .LE. 0 ) GOTO 6099

C     MISE A JOUR DU TMS 'VECTEUR"TEMPERATURE'
      MCN( MNVECT + WBCOVE ) = NTDLTE
      MCN( MNVECT + WBVECT ) = NDSM
      MCN( MNVECT + WBCPIN ) = NDSM

C     COPIE DES NBDLTE TEMPERATURES NCAS1
      MNTHDL  = MNVECT + WECTEU
      MNTHDLD = ( MNTHDL + 1 ) / MOREE2
      CALL TRTATD( temp(NCAS1)%dptab(1), DMCN(MNTHDLD), NTDLTE )

C     L'ADRESSE DU TEMPS DERRIERE LA TEMPERATURE DU CAS NCAS1
      TEMPS = RMCN( MNTIMES - NCAST0 + NCAS1 )
      RMCN( MNVECT + WECTEU + NTDLTE * MOREE2 ) = TEMPS

C     LA DATE
      CALL ECDATE( MCN(MNVECT) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNVECT + MOREE2 ) = NONMTD( '~>>>VECTEUR' )

C     CALCUL LES FLUX AUX POINTS D'INTEGRATION NUMERIQUE DES FACES
C     DES ELEMENTS FINIS LAGRANGE ISOPARAMETRIQUES 2D 3D AXISYMETRIQUES
C     CONSTRUCTION DES TMS FLUXPT"NMTYEL DE L'OBJET
      CALL THEFLU( KNOMOB, NTLXOB, MNTOPO, NOAXIS, D2PI,
     &             NDIM,   MOREE2, NDSM,   NTDLTE,
     &             NBTYEL, MNNPEF, NDPGST, MNTPOB,
     &             MNTAUX, MNXYZN, NUMIOB, NUMAOB, MNDOEL,
     &             MOTAEL, MNTAEL, MNX,    MNVECT )

C     TRACE DES TEMPERATURES NCAS1
      CALL TRTHER( KNOMOB, 1, IERR ) 

 6099 IF( MNX    .GT. 0 ) CALL TNMCDS( 'REEL',  NBDLMX*3, MNX )
      IF( MNTAEL .GT. 0 ) CALL TNMCDS( 'REEL2', MOTAEL, MNTAEL )
      IF( MNTAUX .GT. 0 ) CALL TNMCDS( 'REEL2', MOTAUX, MNTAUX )
      IF( NTVECT .GT. 0 ) THEN
C        LES TMS VECTEUR FLUXPT DTEMPERATURE SONT DETRUITS
         CALL LXTSDS( NTLXOB, 'VECTEUR"TEMPERATURE' )
         CALL LXTSDS( NTLXOB, 'FLUXPT"2P2C' )
         CALL LXTSDS( NTLXOB, 'DTEMPERATURE"2P2C' )
      ENDIF
      GOTO 55

C     *****************************************************************
C     PARCOURS DES PARTICULES DURANT L'INTERVALLE DE TEMPS[NCAS0:NCAS1]
C     *****************************************************************
C     MNPART: ADRESSE MCN DU TABLEAU DES XYZ + XYZVIT + RAYON + TEMPS0
C             des NBPART PARTICULES
C             en 2D: (Z0, VZ0, RAYON, NE SONT PAS ICI UTILISES)
 7000 IF( NCAS0 .EQ. NCAS1 ) GOTO 55
      CALL PARTICULE( KNOMOB, NTLXOB, NDIM, MNNPEF, MNXYZP,
     %                NBSOMT, NBNOVI, MCN(MNNOSO),
     %                RMCN(MNTIMES-NCAST0+NCAS0), NCAS0, NCAS1,
     %                vitx, vity, vitz, pres, VITMAX, VITMOY,
     %                NBSFTR, MNNOSFTR, MNXYSFTR, MNEFSFTR, IERR )
      GOTO 55


C     ******************************************************************
C     en 2D SEULEMENT :  CALCUL de la FONCTION COURANT PSI 
C     DEFINIE PAR (Vx,Vy) = (dPsi/dy, -dPsi/dx) et
C     - Laplacien Psi = -d2Psi/dx2 -d2Psi/dy2 = dVy/dx - dVx/dy
C     dPsi/dn = - Vy nx + Vx ny   sur la FRONTIERE du DOMAINE
C     pour TRACER les LIGNES de COURANT de l'ECOULEMENT
C     ******************************************************************
 9000 IF( NDIM .NE. 2 ) GOTO 55
      CALL LXTSOU( NTLXOB, 'VECTEUR"COURANT', NTCOUR, MNCOUR )
      IF( NTCOUR .GT. 0 ) THEN
         IF( AVANT( MCN(MNVECT), MCN(MNCOUR) ) ) THEN
C           VITESSES-PRESSIONS ANTERIEURES AUX FONCTIONS COURANTES
C           DEJA CALCULEES ET RECUPEREES SANS LES RECALCULER
            GOTO 6040
         ENDIF
         CALL LXTSDS( NTLXOB, 'VECTEUR"COURANT' )
      ENDIF

C     CALCUL DE LA FONCTION COURANT AUX NBVVIPR TEMPS
C     TEMPS CPU DE CALCUL DES FONCTIONS COURANT
      DCPU = DINFO( 'CPU' )

      IF( LANGAG .EQ. 0 ) THEN
         PRINT*, 'CALCUL des',NBVVIPR,' FONCTIONS COURANT'
      ELSE
         PRINT*, 'COMPUTATION of',NBVVIPR,' STREAM FUNCTIONS'
      ENDIF

C     CONSTRUCTION DES NBVVIPR VECTEURS DE LA FONCTION COURANT PSI
      L = WECTEU + NBSOMT * NBVVIPR * MOREE2 + NBVVIPR
      CALL LXTNDC( NTLXOB, 'VECTEUR"COURANT', 'REEL2', L )
      CALL LXTSOU( NTLXOB, 'VECTEUR"COURANT', NTCOUR, MNCOUR )
      IF( NTCOUR .LE. 0 ) GOTO 55

      MNPSI = MNCOUR + WECTEU

C     DESALLOCATION DE LA MATRICE MORSE SI ELLE EXISTE
      IF( IALPG .EQ. 0 ) THEN
         DEALLOCATE( PG )
         IALPG = 1
      ENDIF

C     DESALLOCATION DES TABLEAUX DU GC
      IF( MNAUX1 .GT. 0 ) CALL TNMCDS( 'REEL2',  MOAUX1, MNAUX1 )
      IF( MNSTFR .GT. 0 ) CALL TNMCDS( 'ENTIER', NBSOMT, MNSTFR )
      IF( MNPTDG .GT. 0 ) CALL TNMCDS( 'ENTIER', NBSOMT, MNPTDG )
      IF( MNLPCO .GT. 0 ) CALL TNMCDS( 'ENTIER', NBCOPG, MNLPCO )

C     TRIANGLES: CONSTRUCTION IMPOSEE du tms arete de l'OBJET 2D
      CALL LXTSOU( NTLXOB, 'ARETE', NTAROB, MNAROB )
      IF( NTAROB .LE. 0 ) THEN
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
            GOTO 9997
         ENDIF
      ENDIF

C     LE NOMBRE D'ENTIERS PAR ARETE
      MOARET = MCN( MNAROB + WOARET )
C     LA MAJORATION DU NOMBRE D'ARETES
      MXARET = MCN( MNAROB + WXARET )
C     LE NOMBRE D'ARETES FRONTALIERES NON SUR LIGNES UTILISATEUR
C     NBARFB = MCN( MNAROB + WBARFB )
C     LE NOMBRE D'ARETES INTERFACES   NON SUR LIGNES UTILISATEUR
C     NBARIN = MCN( MNAROB + WBARIN )
C     LE NUMERO MINIMAL DE LIGNE DE L'OBJET
      NUMILF = MCN( MNAROB + WUMILF )
C     LE NUMERO MAXIMAL DE LIGNE DE L'OBJET
      NUMXLF = MCN( MNAROB + WUMXLF )
C     LE NUMERO DE LA PREMIERE ARETE FRONTALIERE NON SUR LIGNES DE L'OBJET
C     L1ARFB = MCN( MNAROB + W1ARFB )
C     LE NUMERO DE LA PREMIALE ARETE INTERFACE NON SUR LIGNES DE L'OBJET
C     L1ARIN = MCN( MNAROB + W1ARIN )

C     ADRESSE MCN DU 1-ER MOT DU TABLEAU LARETE
      MNLARE = MNAROB + W1LGFR + NUMXLF - NUMILF + 1
C     LARETE : TABLEAU DES ARETES DU MAILLAGE
C     LARETE(1,I)= NO DU 1-ER  SOMMET DE L'ARETE (0 SI PAS D'ARETE)
C     LARETE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C     LARETE(3,I)= CHAINAGE HACHAGE SUR ARETE SUIVANTE
C     LARETE(4,I)= NUMERO DU 1-ER TRIANGLE CONTENANT CETTE ARETE
C                + NO TYPE EF (1 a NBTYEF) * 100 000 000
C                  0 SI PAS DE 1-ER  TRIANGLE
C     LARETE(5,I)= NUMERO DU 2-EME TRIANGLE CONTENANT CETTE ARETE
C                + NO TYPE EF (1 a NBTYEF) * 100 000 000
C                  0 SI PAS DE 2-EME TRIANGLE
C     LARETE(6,I)= NUMERO DANS LARETE DE L'ARETE SUIVANTE
C                 SOIT DANS LE CHAINAGE D'UNE LIGNE J ENTRE NUMILF ET NUMXLF
C                 SOIT DANS LE CHAINAGE DES ARETES FRONTALIERES
C                 0 SI C'EST LA DERNIERE

C     CONSTRUCTION DES TABLEAUX DE POINTEURS SUR LA MATRICE MORSE
C     DE L'OPERATEUR LAPLACIEN
C
      IF( NUTYEL .EQ. 13 .OR. NUTYEL .EQ. 19 ) THEN

C        EF BREZZI-FORTIN
C        NOMBRE DE NOEUDS=SOMMETS DE L'EF
         NBNOEF = NDIM + 1
C        MOMENTANEMENT LE NOMBRE DE NOEUDS = NB SOMMETS
         NBN = MCN( MNXYZN + WNBNOE )
         MCN( MNXYZN + WNBNOE ) = NBSOMT
         CALL PRGCMC( MNTOPO, MCN(MNNPEF), MNXYZN, 1, 1,
     %                MNPTDG, MNLPCO, IERR )
C        RESTAURATION DE LA VALEUR INITIALE
         MCN( MNXYZN + WNBNOE ) = NBN

      ELSE

C        EF TAYLOR-HOOD
C        CALCUL DES POINTEURS DE LA MATRICE MORSE P1 PAR SOMMETS
C        ALORS QUE LE MAILLAGE EST P2
         CALL PRGCMCST( NBSOMT, MCN(MNNOSO), MNTOPO, MCN(MNNPEF),MNXYZN,
     %                  1, 1, MNPTDG, MNLPCO, IERR )

      ENDIF

C     VERSION ALLOCATION de la MATRICE PG en LANGAGE FORTRAN 90
      NBCOPG = MCN( MNPTDG + NBSOMT )
      PRINT*, 'ALLOCATION DEMAND  of',NBCOPG,
     %                ' DOUBLE PRECISION of [PG] MATRIX'
      ALLOCATE( PG(1:NBCOPG), STAT=IALPG )
      IF( IALPG .NE. 0 ) THEN
       PRINT*,'ALLOCATION ERROR   of',NBCOPG,
     %                ' DOUBLE PRECISION of [PG] MATRIX'
         IALR = IALPG
         GOTO 9997
      ENDIF
      PRINT*, 'ALLOCATION CORRECT of',NBCOPG,
     %                ' DOUBLE PRECISION of [PG] MATRIX'
      PRINT*

C     CONSTRUCTION DE LA MATRICE MORSE AVEC EF P1 DE L'OPERATEUR LAPLACIEN
      CALL LAPLACA( 1D0,    NDIM,   NBSOMT,  MCN(MNXYZN+WYZNOE),
     %              NBNOEF, NBELEM, MCN(MNELE+WUNDEL),
     %              NBNOEU, MCN(MNNOSO),
     %              NBCOPG, MCN(MNPTDG), MCN(MNLPCO), PG )

C     AFFICHAGE PARTIEL DE LA MATRICE PG
      call affvect( 'PG MORSE PSI SANS CL', 20, PG )

C     CONSTRUCTION DU TABLEAU des DL FIXES
C     (1 SI SOMMET FRONTALIER, 0 SI INTERNE
      CALL TNMCDC( 'ENTIER', NBSOMT, MNSTFR )
      IF( MNSTFR .LE. 0 ) GOTO 9900

cccC     VERSION CONDITION AUX LIMITES de DIRICHLET HOMOGENE
cccC     CONSTRUCTION D'UN INDICATEUR DE TOUS LES SOMMETS FRONTALIERS
cccC     NOSTFR(N)=0 SI SOMMET N NON FRONTALIER
cccC              =1 SI SOMMET N     FRONTALIER
ccc      CALL NUSTFR( NBSOMT, MOARET, MXARET, MCN(MNLARE),
ccc     %             NBSTFR, MCN(MNSTFR) )
cccC
ccc      DO N=1,NBSOMT
ccc         IF( MCN(MNSTFR-1+N) .NE. 0 ) THEN
cccC           SOMMET FRONTALIER
ccc            PG( MCN(MNPTDG + N ) ) = 1D0
ccc         ENDIF
ccc      ENDDO

C     A PRIORI TOUS LES DL SONT LIBRES => 0
      CALL AZEROI( NBSOMT, MCN(MNSTFR) )

C     VERSION CL NEUMANN: +Integrale de  tLambda ( -v nx + u ny ) dGamma
C     + 1 DL FIXE ARBITRAIREMENT (LE PREMIER DL GLOBAL)
C     MISE A 1 DU DL FIXE
      MCN( MNSTFR ) = 1
C     MODIFICATION DU COEFFICIENT DIAGONAL
C     LE NO DU COEFFICIENT DIAGONAL 1 EST LE PREMIER COEFFICIENT STOCKE
C     CETTE VALEUR POURRAIT SERVIR DE MARQUEUR D'UN DL FIXE POUR PG
      PG( 1 ) = 1D0

cccC     AFFICHAGE PARTIEL DE LA MATRICE PG
ccc      call affvect( 'PG MORSE PSI AVEC CL NEUMANN', 20, PG )

C     LES 3 TABLEAUX AUXILIAIRES du SOUS PROGRAMME Gcaxb + BG
C     MNAUX1, MNAUX2, MNAUX3 : ADRESSES DE TABLEAUX AUXILIAIRES
      MOAUX1 = NBSOMT * 4
      CALL TNMCDC( 'REEL2', MOAUX1, MNAUX1 )
      IF( MNAUX1 .LE. 0 ) GOTO 9900
      MNAUX2 = MNAUX1 + NBSOMT * MOREE2
      MNAUX3 = MNAUX2 + NBSOMT * MOREE2
      MNBG   = MNAUX3 + NBSOMT * MOREE2

C     CONSTRUCTION DU SECOND MEMBRE DU SYSTEME = dv/dx-du/dy
      MNPS0 = MNPSI
      MNPS  = MNPSI

C     VECTEUR INITIAL NUL POUR LE PREMIER GRADIENT CONJUGUE
      CALL AZEROD( NBSOMT, MCN(MNPS0) )

      DO NCAS = NCAS0, NCAS1

C        MISE A ZERO DU SECOND MEMBRE
         CALL AZEROD( NBSOMT, MCN(MNBG) )
C
C        CALCUL DU SECOND MEMBRE DU SYSTEME: -LAPLACIEN PSI = dv/dx-du/dy
C        C-A-D  Integrale e  tLambda ( dv/dx-du/dy ) dX
         CALL LAPLACB(NBSOMT, MCN(MNXYZN+WYZNOE),
     %                NBNOEF, NBELEM, MCN(MNELE+WUNDEL), MCN(MNNOSO),
     %                NBNOVI,  NCAS0,  NCAS1, NCAS, vitx, vity,
     %                MCN(MNBG) )

C        CONDITION de NEUMANN
C        + Integrale de  tLambda ( -v nx + u ny ) dGamma
         CALL LAPLACF(NBSOMT, MCN(MNXYZN+WYZNOE), NBELEM, MCN(MNNOSO),
     %                MOARET, MXARET, MCN(MNLARE),
     %                NBNOVI, NCAS0,  NCAS1, NCAS, vitx, vity,
     %                MCN(MNBG) )

C        MISE A ZERO DU DL 1 POUR ASSURER L'INVERSIBILITE DU SYSTEME
         DMCN( (MNBG+1)/MOREE2 ) = 0D0

cccC        CONDITION de DIRICHLET HOMOGENE POUR PSI FONCTION COURANT
ccc         MN  = (MNBG -1)/MOREE2
ccc         DO N=1,NBSOMT
ccc            IF( MCN(MNSTFR-1+N) .NE. 0 ) THEN
cccC              SOMMET FRONTALIER => DIRICHLET HOMOGENE SUR LE SECOND MEMBRE
ccc               DMCN( MN + N ) = 0D0
ccc            ENDIF
ccc         ENDDO
C
C        RESOLUTION PAR GC SIMPLE SANS PRECONDITIONNEMENT
C        LE VECTEUR INITIAL DU GC EST LE DERNIER VECTEUR CALCULE
ccc         print *
ccc         print *,'A COMPARER a Psi exact CAS',1,' au TEMPS 0.00'
ccc         call affvect( 'Avant Gcaxb BG=', 30, MCN(MNBG) )
         CALL GCAXB( NBSOMT,      MCN(MNSTFR),
     %               MCN(MNPTDG), MCN(MNLPCO),  PG,
     %               MCN(MNBG),   MCN(MNPS0),
     %               MCN(MNAUX1), MCN(MNAUX2),  MCN(MNAUX3),
     %               MCN(MNPS),   IERR )

ccc         call affvect( 'Apres Gcaxb PSI=', 30, MCN(MNPS) )
cccC        PSI = PSI - min PSI
ccc         call tramin( NBSOMT,  MCN(MNPS) )
C
C        PASSAGE AU CAS SUIVANT
         MNPS0 = MNPS
         MNPS  = MNPS + MOREE2 * NBSOMT

      ENDDO

ccc      print *
ccc      print *,'A COMPARER a Psi exact CAS',1,' au TEMPS',
ccc     %         RMCN(MNTIMES-NCAST0+NCAS0)
ccc      call affvect( 'Apres Gcaxb Fonction COURANT PSI=',
ccc     %               500, MCN(MNPSI) )
C
C     MISE A JOUR DU TMS 'VECTEUR"COURANT'
      MCN( MNCOUR + WBCOVE ) = NBSOMT
      MCN( MNCOUR + WBVECT ) = NBVVIPR
      MCN( MNCOUR + WBCPIN ) = NBVVIPR
C     LA DATE
      CALL ECDATE( MCN(MNCOUR) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNCOUR + MOREE2 ) = NONMTD( '~>>>VECTEUR' )
C     COPIE DES TEMPS DES VECTEURS STOCKES DERRRIERE LES NBVVIPR PSI
      L = WECTEU + NBSOMT * NBVVIPR * MOREE2
      CALL TRTATA( MCN(MNTIMES), MCN(MNCOUR+L), NBVVIPR )

C     RENDRE LA PLACE MEMOIRE OCCUPEE INUTILEMENT
      IF( IALPG .EQ. 0 ) THEN
         DEALLOCATE( PG )
         IALPG = 1
      ENDIF
      IF( MNPTDG .GT. 0 ) CALL TNMCDS( 'ENTIER', NBSOMT, MNPTDG )
      IF( MNLPCO .GT. 0 ) CALL TNMCDS( 'ENTIER', NBCOPG, MNLPCO )
C
C     TEMPS CPU DE CALCUL DES FONCTIONS COURANT
      DCPU = DINFO( 'DELTA CPU' )
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*, 'TEMPS CPU CALCUL DES ',NBVVIPR,
     % ' FONCTIONS COURANT=',DCPU,' SECONDES'
      ELSE
         PRINT*, 'CPU TIME of',NBVVIPR,' STREAM FUNCTIONS=',
     %   DCPU,' SECONDS'
      ENDIF

C     CALCUL DES MIN ET MAX DES NBVVIPR VECTEURS PSI FONCTION COURANT en 2D
C     =====================================================================
 6040 NBSOMT = MCN( MNCOUR + WBCOVE )
      NBVVIPR = MCN( MNCOUR + WBVECT )
      MDVECT = (MNCOUR+WECTEU+1)/MOREE2
      CALL MXVECT( NBSOMT, NBVVIPR, DMCN(MDVECT),
     %             PSIMIN, NOPSMIN, NCPMIN, PSIMAX, NOPSMAX, NCPMAX )
      IF( LANGAG .EQ. 0 ) THEN
       WRITE(IMPRIM,16040) PSIMIN, NOPSMIN, NCPMIN,
     %                     RMCN(MNTIMES-NCAST0+NCPMIN),
     %                     PSIMAX, NOPSMAX, NCPMAX,
     %                     RMCN(MNTIMES-NCAST0+NCPMAX)
      ELSE
       WRITE(IMPRIM,26040) PSIMIN, NOPSMIN, NCPMIN,
     %                     RMCN(MNTIMES-NCAST0+NCPMIN),
     %                     PSIMAX, NOPSMAX, NCPMAX,
     %                     RMCN(MNTIMES-NCAST0+NCPMAX)
      ENDIF
16040 FORMAT(' FONCTION COURANT MINIMALE =',G13.6,' au NOEUD',I9,
     %' Cas',I5,' Temps',G14.6/
     %       ' FONCTION COURANT MAXIMALE =',G13.6,' au NOEUD',I9,
     %' Cas',I5,' Temps',G14.6)
26040 FORMAT(' MINIMUM STREAM FUNCTION=',G13.6,' at Node',I9,
     %' Case',I5,' Time',G14.6/
     %       ' MAXIMUM STREAM FUNCTION=',G13.6,' at Node',I9,
     %' Case',I5,' Time',G14.6)

C     OPTIONS DU TRACE DE LA FONCTION COURANT 2D des CAS NCAS0 A NCAS1
C     ================================================================
      MDVECT = (MNCOUR+1+WECTEU)/MOREE2
      MODECO = 7
      PSMIN  = REAL( PSIMIN )
      PSMAX  = REAL( PSIMAX )

 6050 CALL LIMTCL( 'tracoura', NMTCL0 )
      IF( NMTCL0 .LE. 0 ) THEN
         IF( MNSTFR .GT. 0 ) CALL TNMCDS( 'ENTIER', NBSOMT, MNSTFR )
         IF( MNAUX1 .GT. 0 ) CALL TNMCDS( 'REEL2',  MOAUX1, MNAUX1 )
         GOTO 55
      ENDIF
      GOTO( 6100, 6200, 6200, 6200, 6900, 6050, 6050, 6050, 6900),NMTCL0

C     TRACE DES LIGNES ISO-FONCTION COURANT
C     =====================================
 6100 CALL TRISOT( NDIM,   KNOMOB,   MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1, NDPGSTP1,
     %             NCAS0,  NCAS1,    NBSOMT,
     %             0,      DMCN(MDVECT), vitm,
     %             PSMIN,  NOPSMIN,  NCPMIN, PSMAX,  NOPSMAX, NCPMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 6050

C     TRACE DES ZONES DE COULEURS FONCTION COURANT
C     ============================================
 6200 CALL TRZONT( 0,      NDIM,     KNOMOB, MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1, MNXYZNP1, NDPGSTP1,
     %             NCAS0,  NCAS1,    NBSOMT,
     %             0,      DMCN(MDVECT), vitm,
     %             PSMIN,  NOPSMIN,  NCPMIN, PSMAX,  NOPSMAX, NCPMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 6050

C     TRACE EN 2D SURFACE(X,Y,FONCTION COURANT(X,Y))
C     ==============================================
 6900 CALL TRZTXY( 0,      NDIM,     KNOMOB, MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1,
     %             NCAS0,  NCAS1,    NBSOMT,
     %             0,      DMCN(MDVECT), vitm,
     %             PSMIN,  NOPSMIN,  NCPMIN, PSMAX,  NOPSMAX, NCPMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 6050


C     ******************************************************************
C     CALCUL du VECTEUR VORTICITE =  Rot( Vitesse )
C     ROT V: en 3d: RotX=dw/dy-dv/dz   VECTEUR A 3 COMPOSANTES
C                   RotY=du/dz-dw/dx
C                   RotZ=dv/dx-du/dy et
C                   MODULE=SQRT(RotX**2+RotY**2+RotZ**2)
C            en 2d: dv/dx-du/dy EST UN SCALAIRE
C     ******************************************************************
C     ALLOCATION DYNAMIQUE D'UN TABLEAU AUXILIAIRE ROTV ENCOMBRANT
10000 IF( IALROTV .NE. 0 ) THEN

         IF( NDIM .EQ. 3 ) THEN
C           NOMBRE DE COMPOSANTES STOCKEES DANS ROTV
            NBCOMP = 4
         ELSE
            NBCOMP = 1
         ENDIF
         MOROT = NBSOMT * NBVVIPR * NBCOMP
         PRINT*
         PRINT*, 'ALLOCATION DEMAND  of',MOROT,
     %   ' DOUBLE PRECISION of ROTV(',NBSOMT,',',NBVVIPR,',',NBCOMP,
     %   ') VECTORS'
         ALLOCATE ( ROTV(1:NBSOMT,NCAS0:NCAS1,1:NBCOMP), STAT=IALROTV )
         IF( IALROTV .NE. 0 ) THEN
            PRINT*, 'ALLOCATION ERROR   of',MOROT,
     %       ' DOUBLE PRECISION of ROTV(',NBSOMT,',',NBVVIPR,',',NBCOMP,
     %       ') VECTORS'
            GOTO 55
         ENDIF
         PRINT*, 'ALLOCATION CORRECT of',MOROT,
     %       ' DOUBLE PRECISION of ROTV(',NBSOMT,',',NBVVIPR,',',NBCOMP,
     %       ') VECTORS'
         PRINT*

      ENDIF

C     RESERVATION DE LA MATRICE DIAGONALE DE PROJECTION P1
      IF( MNAGD .EQ. 0 ) CALL TNMCDC( 'REEL2', NBSOMT, MNAGD )

      IF( NUTYEL .EQ. 13 ) THEN
C
C        TRIANGLE Brezzi Fortin   Calcul ROTV(*,*) 2D
         CALL VORTBF2( NBNOVI, MCN(MNXYZN+WYZNOE),
     %                 NBELEM, MCN(MNELE+WUNDEL),
     %                 NCAS0,  NCAS1, vitx, vity,
     %                 NBSOMT, MCN(MNAGD), ROTV )

      ELSE IF( NUTYEL .EQ. 15 ) THEN
C
C        TETRAEDRE Taylor-Hood  Calcul ROTV(*,*)
         CALL VORTTH2( NBNOVI, MCN(MNXYZN+WYZNOE),
     %                 NBELEM, MCN(MNELE+WUNDEL), MCN(MNNOSO),
     %                 NCAS0,  NCAS1, vitx, vity,
     %                 NBSOMT, MCN(MNAGD), ROTV )

      ELSE IF( NUTYEL .EQ. 19 ) THEN
C
C        TRIANGLE Brezzi Fortin   Calcul ROTV(*,*,1:4)
         CALL VORTBF3( NBNOVI, MCN(MNXYZN+WYZNOE),
     %                 NBELEM, MCN(MNELE+WUNDEL),
     %                 NCAS0,  NCAS1, vitx, vity, vitz,
     %                 NBSOMT, MCN(MNAGD),  ROTV )

      ELSE IF( NUTYEL .EQ. 20 ) THEN
C
C        TETRAEDRE Taylor-Hood  Calcul ROTV(*,*,1:4)
         CALL VORTTH3( NBNOVI, MCN(MNXYZN+WYZNOE),
     %                 NBELEM, MCN(MNELE+WUNDEL), MCN(MNNOSO),
     %                 NCAS0,  NCAS1, vitx, vity, vitz,
     %                 NBSOMT, MCN(MNAGD),  ROTV )

      ENDIF

C     MODIFICATION DU NUMERO DE NOEUDS DES SOMMETS DANS LE TABLEAU
C     DES FACES POUR LE TETRAEDRE DE TAYLOR-HOOD
      IF( NUTYEL .EQ. 20 ) THEN
         CALL ARETFRP2P1( KNOMOB, NBNOEU, MCN(MNNOSO), IERR )
         IF( IERR .NE. 0 ) GOTO 55
      ENDIF

C     PAR DEFAUT TRACE DU MODULE DU ROTATIONNEL DE LA VITESSE
      NOCOMP = NBCOMP

C     MIN MAX de ROT Vitesse aux SOMMETS
C     ----------------------------------
 7002 IF( NDIM .LE. 2 ) NOCOMP=1
      CALL MAXCOMP( NDIM,   NOCOMP,  NBNOVI, MCN(MNXYZN+WYZNOE),
     %              NBSOMT, NBVVIPR, ROTV,
     %              COOEXT, NOROMIN, NRAMIN, ROTMIN,
     %                      NOROMAX, NRAMAX, ROTMAX, ROTMOY )
C
      ROMIN = REAL( ROTMIN )
      ROMAX = REAL( ROTMAX )

C     OPTIONS DU TRACE de la COMPOSANTE NOCOMP de ROT Vitesse
C     des TEMPS NCAS0 a NCAS1
      MODECO = 14

C     LES OPTIONS DE TRACE DU ROTATIONNEL DES VITESSES
C     *****************************************************************
 7007 IF( NDIM .LE. 2 ) THEN
         CALL LIMTCL( 'trarot2d', NMTCL0 )
      ELSE
         CALL LIMTCL( 'trarot3d', NMTCL0 )
      ENDIF

      IF( NMTCL0 .LT. 0 ) THEN
C        RETOUR AU NUMERO DE NOEUDS DES SOMMETS DANS LE TABLEAU DES FACES
C        POUR LE TETRAEDRE DE TAYLOR-HOOD
         IF( NUTYEL .EQ. 20 ) THEN
            CALL ARETFRP1P2( KNOMOB, NBSOMT, MCN(MNSONO) )
         ENDIF
         GOTO 55
      ENDIF

      IF( NMTCL0 .EQ. 0 ) THEN
C        CHOIX du No de la COMPOSANTE ou MODULE du RotV
         CALL LIMTCL( 'noxyzmod', NOCOMP )
         IF( NOCOMP .LE. 0 ) THEN
            GOTO 55
         ENDIF
         IF( NDIM .LE. 2 ) THEN
C           NECESSAIRE POUR NDIM=2
            NOCOMP = 1
         ENDIF
         GOTO 7007
      ENDIF

      IF( NDIM .EQ. 3 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'Trace 3D de la COMPOSANTE',NOCOMP,' de Rot Vitesse'
         ELSE

            PRINT*,'3D Drawing of the VORTICITY COMPONENT',NOCOMP
         ENDIF
      ENDIF

      GOTO( 7010, 7020, 7030, 7040, 7050, 7060, 55, 55, 7090, 55,
     %      7110, 7120, 7130, 7140 ), NMTCL0
      GOTO 7007
C
C     TRACE DES ISO-ROTATIONNEL VITESSES (LIGNES en 2d, SURFACES en 3d)
C     =================================================================
 7010 CALL TRISOT( NDIM,   KNOMOB,   MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1, NDPGSTP1,
     %             NCAS0,  NCAS1,    NBSOMT,
     %             0,      ROTV(1,NCAS0,NOCOMP),   vitm,
     %             ROMIN,  NOROMIN,  NRAMIN,   ROMAX, NOROMAX, NRAMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 7007
C
C     TRACE DES ZONES DE COULEURS ISO-ROTATIONNEL VITESSE
C     ===================================================
 7020 CALL TRZONT( 0,      NDIM,     KNOMOB,  MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1,MNXYZNP1, NDPGSTP1,
     %             NCAS0,  NCAS1,    NBSOMT,
     %             0,      ROTV(1,NCAS0,NOCOMP),  vitm,
     %             ROMIN,  NOROMIN,  NRAMIN,  ROMAX, NOROMAX, NRAMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 7007
C
C     TRACE DES ZONES DE COULEURS ISO-ROTV PAR SECTIONS X ou Y ou Z=CTE
C     =================================================================
 7030 IF( NDIM .EQ. 2 ) GOTO 7020
      CALL TRPLSE( 0,      KNOMOB,   MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1, NDPGSTP1,
     %             NCAS0,  NCAS1,   NBSOMT,
     %             0,      ROTV(1,NCAS0,NOCOMP), vitm,
     %             ROMIN,  NOROMIN, NRAMIN, ROMAX, NOROMAX, NRAMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 7007
C
C     TRACE PROFILS DE COULEURS Rot VITESSE PAR SECTIONS X ou Y ou Z=CTE
C     ==================================================================
 7040 IF( NDIM .EQ. 2 ) GOTO 7090
      CALL TRPLSE( 1,      KNOMOB,   MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1, NDPGSTP1,
     %             NCAS0,  NCAS1,   NBSOMT,
     %             0,      ROTV(1,NCAS0,NOCOMP), vitm,
     %             ROMIN,  NOROMIN, NRAMIN, ROMAX, NOROMAX, NRAMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 7007
C
C     TRACE DE ROT VITESSE LE LONG D'UNE DROITE 3D DEFINIE PAR 2 POINTS
C     =================================================================
 7050 IF( NDIM .EQ. 2 ) GOTO 7090
      CALL TRLLDR( 0,      NDIM,     KNOMOB, MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1, NDPGSTP1,
     %             NCAS0,  NCAS1,    NBSOMT,
     %             0,      ROTV(1,NCAS0,NOCOMP), vitm,
     %             ROMIN,  NOROMIN,  NRAMIN, ROMAX, NOROMAX, NRAMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 7007
C
C     TRACE DU MODULE DU ROT DE LA VITESSE SUR DES SURFACES DE L'OBJET 3D
C     ===================================================================
C     LA COMPOSANTE NOCOMP du ROTATIONNEL DE  LA VITESSE
C     P2 POUR TAYLOR-YOUNG ou P1+BULLE POUR BREZZI-FORTIN
 7060 IF( NDIM .EQ. 2 ) GOTO 7090
      CALL TRSO1SO( INTERPP, 0,  MODECO, KNOMOB, NTLXOB, MNDFOB,
     %              NBSOMT,  NCAS0, NCAS1,
     %              0,       ROTV(1,NCAS0,NOCOMP), vitm,
     %              ROMIN,   ROMAX, RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 7007
C
C     TRACE EN 2D SURFACE(X,Y,ROTV(X,Y))
C     ==================================
 7090 IF( NDIM .NE. 2 ) GOTO 7007
      CALL TRZTXY(      0, NDIM,     KNOMOB, MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1,
     %             NCAS0,  NCAS1,   NBSOMT,
     %             0,      ROTV,    vitm,
     %             ROMIN,  NOROMIN, NRAMIN, ROMAX, NOROMAX, NRAMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 7007
C
C     TRACE DE LA COMPOSANTE X DU ROTATIONNEL DE LA VITESSE
C     =====================================================
 7110 NOCOMP = 1
      GOTO 7002
C
C     TRACE DE LA COMPOSANTE Y DU ROTATIONNEL DE LA VITESSE
C     =====================================================
 7120 IF( NDIM .GT. 2 ) THEN
         NOCOMP = 2
      ELSE
         NOCOMP = 1
      ENDIF
      GOTO 7002

C     TRACE DE LA COMPOSANTE Z DU ROTATIONNEL DE LA VITESSE
C     =====================================================
 7130 IF( NDIM .GT. 2 ) THEN
         NOCOMP = 3
      ELSE
         NOCOMP = 1
      ENDIF
      GOTO 7002

C     TRACE DU MODULE DU ROTATIONNEL DE LA VITESSE
C     ============================================
 7140 IF( NDIM .GT. 2 ) THEN
         NOCOMP = 4
      ELSE
         NOCOMP = 1
      ENDIF
      GOTO 7002


C     CALCUL DE L'INTEGRALE DE LA PRESSION d'INTERPOLATION P1 SUR TOUTES
C     LES SURFACES DE L'OBJET POUR TOUS LES VECTEURS STOCKES 1 A NBVVIPR
C     et AUSSI CALCUL DU COUPLE EXERCE SUR TOUTES LES SURFACES
C     ==================================================================
12000 IF( NDIM .NE. 3 ) GOTO 55
      PRINT*
      PRINT*,'PRESSURE INTEGRAL OVER ALL 3D OBJECT SURFACES:'
      WRITE(IMPRIM,19991)
19991 FORMAT(160('='))
      IF( LANGAG .EQ. 0 ) THEN
         NMSOLU = 'PRESSION '
      ELSE
         NMSOLU = 'PRESSURE '
      ENDIF
      CALL INSOSF( 0,     INTERPP, NTLXOB,  MNDFOB,
     %             NCAS0, NCAS1, pres, NMSOLU,
     %             RMCN(MNTIMES-NCAST0+NCAS0)  )
      GOTO 55


C     FLUX NORMAL DE LA VITESSE A TRAVERS LES FACES DES SURFACES 3D DE L'OBJET
C     ========================================================================
13000 IF( NDIM .NE. 3 ) GOTO 55
      MOSFOB = 2 * ( NUMAOB(3) - NUMIOB(3) + 1 ) * NBVVIPR
      CALL TNMCDC( 'REEL2', MOSFOB, MNFVSF )
      print *,'tfluide: fluvitasf avec NBVVIPR=',NBVVIPR,
     %        ' NBNOVI=',NBNOVI,' mnele=',mnele
      CALL FLUVITASF( MNXYZN, MNELE, MCN(MNNDDL),
     %                NCAS0,  NCAS1, vitx, vity, vitz,
     %                NUMIOB(3), NUMAOB(3), RMCN(MNTIMES-NCAST0+NCAS0),
     %                MCN(MNFVSF) )
      CALL TNMCDS( 'REEL2', MOSFOB, MNFVSF )
      GOTO 55


C     Affichage des VITESSES PRESSIONS ERREURS des cas NCAS0 A NCAS1
C     ==============================================================
19000 CALL AFVIPRER( KNOMOB, KNOMFIC, NDIM, NUTYEL, MNXYZN, MCN(MNNDDL),
     %               NBNOVI, NBSOMT, NTDLVP, NTDLTE,
     %               NCAS0,  NCAS1, NBSOMT/2+1,  NBSOMT/2+10,
     %               RMCN(MNTIMES-NCAST0+NCAS0), VXYZPN, TEMPER,
     %               NOFOVI, NOFOPR,
     %               VITMI1, VITMA1, VITMO1, PREMI1, PREMA1, PREMO1 )
      GOTO 55


C     ******************************************************************
C     TRACE DES ERREURS SUR LE MODULE DE LA VITESSE
C     ******************************************************************
19100 NOFOVI = NOFOVITE()
      IF( NOFOVI .LE. 0 ) GOTO 55

C     CALCUL DES ERREURS SUR LA VITESSE AUX NOEUDS ET AU TEMPS NCAS0:NCAS1
C     --------------------------------------------------------------------
      CALL ERREURVIT( KNOMOB, KNOMFIC, NDIM, NUTYEL, MNELE,
     %                MNXYZN, MCN(MNNDDL), NBNOVI, NTDLVP, NTDLTE,
     %                NCAS0,  NCAS1, RMCN(MNTIMES-NCAST0+NCAS0),
     %                VXYZPN, TEMPER,
     %                NOFOVI, VPAUX, VitErMOY,
     %                VitErMIN, NOVitErMIN, NCVitErMIN,
     %                VitErMAX, NOVitErMAX, NCVitErMAX )

ccc      IF( NOEVMIN .LE. NBNOEU ) THEN
      MN = MNXYZN + WYZNOE + 3 * NOEVMIN - 3
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,19104)VitErMIN,NOVitErMIN,(RMCN(MN+K),K=0,2),
     %                      NCVitErMIN,
     %                      RMCN(MNTIMES-NCAST0+NCVitErMIN)
      ELSE
         WRITE(IMPRIM,29104)VitErMIN,NOVitErMIN,(RMCN(MN+K),K=0,2),
     %                      NCVitErMIN,
     %                      RMCN(MNTIMES-NCAST0+NCVitErMIN)
      ENDIF
ccc      ENDIF
19104 FORMAT(' |VITESSE| ERREUR MINIMALE=',G13.6,' au NOEUD',I9,
     %' X=',G14.6,' Y=',G14.6,' Z=',G14.6,
     %' Cas',I5,' Temps',G14.6)
29104 FORMAT(' MINIMUM VELOCITY ERROR=',G13.6,' at Node',I9,
     %' X=',G14.6,' Y=',G14.6,' Z=',G14.6,
     %' Case',I5, ' Time',G14.6)

ccc      IF( NOEVMAX .LE. NBNOEU ) THEN
      MN = MNXYZN + WYZNOE + 3 * NOEVMAX - 3
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,19105)VitErMAX,NOVitErMAX,(RMCN(MN+K),K=0,2),
     %                      NCVitErMAX,
     %                      RMCN(MNTIMES-NCAST0+NCVitErMAX)
      ELSE
         WRITE(IMPRIM,29105)VitErMAX,NOVitErMAX,(RMCN(MN+K),K=0,2),
     %                      NCVitErMAX,
     %                      RMCN(MNTIMES-NCAST0+NCVitErMAX)
      ENDIF
ccc      ENDIF
19105 FORMAT(' |VITESSE| ERREUR MAXIMALE=',G13.6,' au NOEUD',I9,
     %' X=',G14.6,' Y=',G14.6,' Z=',G14.6,
     %' Cas',I5,' Temps',G14.6)
29105 FORMAT(' MAXIMUM |VELOCITY| ERROR=',G13.6,' at Node',I9,
     %' X=',G14.6,' Y=',G14.6,' Z=',G14.6,
     %' Case',I5, ' Time',G14.6)

C     OPTIONS DU TRACE DE L'ERREUR SUR LA |VITESSE| NCAS0 A NCAS1
C     PROVENANCE ERREUR SUR LE LE MODULE DE LA VITESSE D'UN FLUIDE
      MODECO = 11

C
C     OPTIONS DE TRACE DE L'ERREUR SUR LE MODULE DE LA |VITESSE|
C     **********************************************************
 9105 CALL LIMTCL( 'traerrvi', NMTCL0 )
      IF( NMTCL0 .LE. 0 ) GOTO 55
      GOTO( 9110, 9120, 9130, 9140, 9150, 9160, 9105, 9105, 9191),NMTCL0
C
C     TRACE DES ISO-MODULE VITESSES (LIGNES en 2D, SURFACES en 3D)
C     ============================================================
 9110 CALL TRISOT( NDIM,   KNOMOB, MODECO,
     %             NBTYEL, MNNPEF, MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,  NBNOVI,
     %             0,      VPAUX,  vitm,
     %             VitErMIN,NOVitErMIN,NCVitErMIN,
     %             VitErMAX,NOVitErMAX,NCVitErMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 9105
C
C     TRACE DES ZONES DE COULEURS ISOMODULE VITESSE
C     =============================================
 9120  CALL TRZONT( 0,     NDIM,    KNOMOB, MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN, MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             0,      VPAUX,   vitm,
     %             VitErMIN,NOVitErMIN,NCVitErMIN,
     %             VitErMAX,NOVitErMAX,NCVitErMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 9105
C
C     TRACE DES ZONES DE COULEURS ISOMODULES PAR SECTIONS X ou Y ou Z=CTE
C     ===================================================================
 9130 IF( NDIM .EQ. 2 ) GOTO 9120
C     LE MODULE DES VITESSES EN 3D
      CALL TRPLSE( 0,      KNOMOB,  MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             0,      VPAUX,   vitm,
     %             VitErMIN,NOVitErMIN,NCVitErMIN,
     %             VitErMAX,NOVitErMAX,NCVitErMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 9105
C
C     TRACE DES PROFILS DE COULEURS MODULE VITESSE PAR SECTIONS X ou Y ou Z=CTE
C     =========================================================================
 9140 IF( NDIM .EQ. 2 ) GOTO 9191
C     LE MODULE DES VITESSES EN 3D
      CALL TRPLSE( 1,      KNOMOB,  MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN,  NDPGST,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             0,      VPAUX,   vitm,
     %             VitErMIN,NOVitErMIN,NCVitErMIN,
     %             VitErMAX,NOVitErMAX,NCVitErMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 9105
C
C     TRACE DU MODULE VITESSE LE LONG D'UNE DROITE DEFINIE PAR 2 POINTS
C     =================================================================
 9150 IF( NDIM .EQ. 2 ) GOTO 9120
C     LE MODULE VITESSE EN 3D
      CALL TRLLDR( 0,      NDIM,    KNOMOB, MODECO,
     %             NBTYEL, MNNPEF,  MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,   NBNOVI,
     %             0,      VPAUX,   vitm,
     %             VitErMIN,NOVitErMIN,NCVitErMIN,
     %             VitErMAX,NOVitErMAX,NCVitErMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 9105
C
C     TRACE DU MODULE DE LA VITESSE SUR DES SURFACES DE L'OBJET
C     =========================================================
 9160 IF( NDIM .EQ. 2 ) GOTO 9120
C     LE MODULE DE LA VITESSE P2 POUR TAYLOR-YOUNG et
C     P1+BULLE POUR BREZZI-FORTIN
      CALL TRSO1SO( INTERPV, 0,  MODECO, KNOMOB, NTLXOB, MNDFOB,
     %              NBNOVI, NCAS0,  NCAS1,
     %              0,      VPAUX,  vitm,
     %              VitErMIN, VitErMAX, RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 9105

C     TRACE EN 2D SURFACE(X,Y,MODULE ERREUR VITESSE(X,Y))
C     ===================================================
 9191 IF( NDIM .NE. 2 ) GOTO 9105
      CALL TRZTXY( 0,      NDIM,   KNOMOB, MODECO,
     %             NBTYEL, MNNPEF, MNXYZN,
     %             NCAS0,  NCAS1,  NBNOVI,
     %             0,      VPAUX,  vitm,
     %             VitErMIN,NOVitErMIN,NCVitErMIN,
     %             VitErMAX,NOVitErMAX,NCVitErMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 9105


C     ******************************************************************
C     TRACE DE L'ERREUR SUR LA PRESSION DANS TOUT LE MAILLAGE
C     ******************************************************************
19200 NOFOPR = NOFOPRES()
      IF( NOFOPR .LE. 0 ) GOTO 55

C     CALCUL DES ERREURS SUR LA PRESSION AUX SOMMETS ET AU TEMPS NCAS0:NCAS1
C     ----------------------------------------------------------------------
      CALL ERREURPRES(KNOMOB, KNOMFIC, NDIM, NUTYEL, MNXYZNP1,
     %                MCN(MNNDDL), MCN(MNNOSO),
     %                NBSOMT, NBNOVI, NTDLVP, NTDLTE,
     %                NCAS0,  NCAS1, RMCN(MNTIMES-NCAST0+NCAS0),
     %                VXYZPN, TEMPER,
     %                NOFOPR, VPAUX, PreErMOY,
     %                PreErMIN, NOPreErMIN, NCPreErMIN,
     %                PreErMAX, NOPreErMAX, NCPreErMAX )
C
C     LE SOMMET DE PRESSION MINIMALE EST NOPreErMIN
      MN = MNXYZN + WYZNOE - 3 + 3 * NOPreErMIN
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,19204) PreErMIN,NOPreErMIN,
     %                      (RMCN(MN+K),K=0,2),NCPreErMIN,
     %                       RMCN(MNTIMES-NCAST0+NCPreErMIN)
      ELSE
         WRITE(IMPRIM,29204) PreErMIN,NOPreErMIN,
     %                      (RMCN(MN+K),K=0,2),NCPreErMIN,
     %                       RMCN(MNTIMES-NCAST0+NCPreErMIN)
      ENDIF

19204 FORMAT(/' ERREUR PRESSION MINIMALE=',G13.6,
     %' au NOEUD',I7,' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' Cas',I5,' Temps',G14.6)
29204 FORMAT(/' MINIMUM PRESSURE ERROR=',G13.6,
     %' at Node',I7,' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' CASE',I5,' TIME',G14.6)
C
C     LE SOMMET DE PRESSION MAXIMALE EST NOPreErMAX
      MN = MNXYZN + WYZNOE - 3 + 3 * NOPreErMAX
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,19205) PreErMAX,NOPreErMAX,
     %                      (RMCN(MN+K),K=0,2),NCPreErMAX,
     %                       RMCN(MNTIMES-NCAST0+NCPreErMAX)
      ELSE
         WRITE(IMPRIM,29205) PreErMAX,NOPreErMAX,
     %                      (RMCN(MN+K),K=0,2),NCPreErMAX,
     %                       RMCN(MNTIMES-NCAST0+NCPreErMAX)
      ENDIF
19205 FORMAT(' ERREUR PRESSION MAXIMALE=',G13.6,
     %' AU NOEUD',I7,' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' Cas',I5,' Temps',G14.6)
29205 FORMAT(' MAXIMUM PRESSURE ERROR=',G13.6,
     %' at Node',I7,' X=',g13.6,' Y=',g13.6,' Z=',g13.6,
     %' CASE',I5,' TIME',G14.6)
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,19206) PreErMAX-PreErMIN
      ELSE
         WRITE(IMPRIM,29206) PreErMAX-PreErMIN
      ENDIF
19206 FORMAT(' ERREUR PRESSION MAXIMALE - MINIMALE =',G13.6 )
29206 FORMAT(' MAXIMUM - MINIMUM PRESSURE ERROR =',G13.6)
C
C     OPTIONS DU TRACE DE LA PRESSION NCAS0 A NCAS1
      PreErMIN = REAL( PreErMIN )
      PreErMAX = REAL( PreErMAX )
C
C     PROVENANCE ERREUR SUR LA PRESSION D'UN FLUIDE
      MODECO = 12

C     MODIFICATION DU NUMERO DE NOEUDS DES SOMMETS DANS LE TABLEAU
C     DES FACES POUR LE TETRAEDRE DE TAYLOR-HOOD
      IF( NUTYEL .EQ. 20 ) THEN
         CALL ARETFRP2P1( KNOMOB, NBNOEU, MCN(MNNOSO), IERR )
         IF( IERR .NE. 0 ) GOTO 55
      ENDIF

C     LES OPTIONS DE TRACE DE L'ERREUR SUR LA PRESSION
C     ************************************************
 9208 CALL LIMTCL( 'traerrpr', NMTCL0 )

      IF( NMTCL0 .LE. 0 ) THEN
C        RETOUR AU NUMERO DE NOEUDS DES SOMMETS DANS LE TABLEAU
C        DES FACES POUR LE TETRAEDRE DE TAYLOR-HOOD
         IF( NUTYEL .EQ. 20 ) THEN
            CALL ARETFRP1P2( KNOMOB, NBSOMT, MCN(MNSONO) )
         ENDIF
         GOTO 55
      ENDIF

      GOTO( 9210, 9292, 9230, 9240, 9250, 9260, 9270, 55, 9290 ) ,NMTCL0
      GOTO 9208

C     TRACE DES ISO-ERREUR PRESSIONS (LIGNES en 2D, SURFACES en 3D)
C     =============================================================
 9210 CALL TRISOT( NDIM,   KNOMOB,   MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1,  NDPGSTP1,
     %             NCAS0,  NCAS1,   NBSOMT,
     %             0,      VPAUX,   vitm,
     %             PreErMIN,NOPreErMIN,NCPreErMIN,
     %             PreErMAX,NOPreErMAX,NCPreErMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 9208
C
C     TRACE DES ZONES DE COULEURS ISO-ERREUR PRESSIONS
C     ================================================
 9292 CALL TRZONT( 0,      NDIM,     KNOMOB, MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1, MNXYZNP1, NDPGSTP1,
     %             NCAS0,  NCAS1,   NBSOMT,
     %             0,      VPAUX,   vitm,
     %             PreErMIN,NOPreErMIN,NCPreErMIN,
     %             PreErMAX,NOPreErMAX,NCPreErMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 9208
C
C     TRACE DES ZONES DE COULEURS ISO-ERREUR PRESSIONS PAR SECTIONS XouYouZ=CTE
C     =========================================================================
 9230  IF( NDIM .EQ. 2 ) GOTO 9208
C      LA ERREUR PRESSION EN 3D
       CALL TRPLSE( 0,      KNOMOB,   MODECO,
     %              NBTYEL, MNNPEFP1, MNXYZNP1, NDPGSTP1,
     %              NCAS0,  NCAS1,   NBSOMT,
     %              0,      VPAUX,   vitm,
     %              PreErMIN,NOPreErMIN,NCPreErMIN,
     %              PreErMAX,NOPreErMAX,NCPreErMAX,
     %              RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 9208
C
C     TRACE DES PROFILS ERREUR PRESSIONS DE COULEURS PAR SECTIONS XouYou Z=CTE
C     ========================================================================
 9240 IF( NDIM .EQ. 2 ) GOTO 9290
C      L'ERREUR PRESSION EN 3D
       CALL TRPLSE( 1,      KNOMOB,   MODECO,
     %              NBTYEL, MNNPEFP1, MNXYZNP1, NDPGSTP1,
     %              NCAS0,  NCAS1,    NBSOMT,
     %              0,      VPAUX,    vitm,
     %              PreErMIN,NOPreErMIN,NCPreErMIN,
     %              PreErMAX,NOPreErMAX,NCPreErMAX,
     %              RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 9208
C
C     TRACE DE L'ERREUR PRESSION LE LONG D'UNE DROITE DEFINIE PAR 2 POINTS
C     ====================================================================
 9250 IF( NDIM .EQ. 2 ) GOTO 9208
C     L'ERREUR PRESSION EN 3D
      CALL TRLLDR( 0,      NDIM,     KNOMOB, MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1, NDPGSTP1,
     %             NCAS0,  NCAS1,    NBSOMT,
     %              0,     VPAUX,    vitm,
     %             PreErMIN,NOPreErMIN,NCPreErMIN,
     %             PreErMAX,NOPreErMAX,NCPreErMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 9208
C
C     TRACE DE L'ERREUR PRESSION SUR UNE OU PLUSIEURS DES SURFACES DE L'OBJET
C     =======================================================================
 9260 IF( NDIM .EQ. 2 ) GOTO 9208
C     LA ERREUR PRESSION P1 EN 3D
      CALL TRSO1SO( INTERPP, 0, MODECO, KNOMOB, NTLXOB, MNDFOB,
     %              NBSOMT,  NCAS0, NCAS1,
     %              0,       VPAUX, vitm,
     %              PreErMIN,  PreErMAX, RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 9208
C
C     TRACE DE L'ERREUR PRESSION SUR TOUTES LES SURFACES FRONTALIERES
C     ===============================================================
 9270 IF( NDIM .EQ. 2 ) GOTO 9208
C     LA ERREUR PRESSION P1 EN 3D
      CALL TRSOSF( 0,       KNOMOB, MODECO, MNXYZN,
     %             NBSOMT,  NCAS0,  NCAS1,
     %             0,       VPAUX,  vitm,
     %             PreErMIN,  PreErMAX, RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 9208

C     TRACE EN 2D SURFACE(X,Y,ERREUR PRESSION(t,X,Y))
C     ===============================================
 9290 IF( NDIM .NE. 2 ) GOTO 9208
      CALL TRZTXY(      0, NDIM,     KNOMOB, MODECO,
     %             NBTYEL, MNNPEFP1, MNXYZNP1,
     %             NCAS0,  NCAS1,    NBSOMT,
     %             0,      VPAUX,    vitm,
     %             PreErMIN,NOPreErMIN,NCPreErMIN,
     %             PreErMAX,NOPreErMAX,NCPreErMAX,
     %             RMCN(MNTIMES-NCAST0+NCAS0) )
      GOTO 9208


cccC     ********************************************************************
cccC     CALCUL et TRACE de ||VITESSE||L2=f(Temps) et ||PRESSION||L2=f(Temps)
cccC     AUX SEULS TEMPS ou LES VECTEURS VITESSE-PRESSION ont ete STOCKES
cccC     ********************************************************************
ccc19400 IF( NDIM .LE. 3 ) GOTO 55
cccC     a mettre a jour avec vitx vity vitz vitm
ccc      CALL FL2NORMES( KNOMOB, IERR, DCPU )
ccc      GOTO 55


C     PAS ASSEZ DE MOTS DANS MCN
C     ==========================
 9900 NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) ='PAS ASSEZ de MOTS pour le SUPER-TABLEAU MCN'
         KERR(2) ='AUGMENTER ce NOMBRE dans incl.pp.inc'
      ELSE
         KERR(1) ='NO SUFFICIENT WORD NUMBER of SUPER-ARRAY MCN'
         KERR(2) ='AUGMENT IT in incl/pp.inc'
      ENDIF
      CALL LEREUR
      GOTO 9997


C     ERREUR D'ALLOCATION vitx1t vity1t vitz1t vitm1t ou pres1t de vxyzmp(NCAS)
C     =========================================================================
 9990 IF( LANGAG .EQ. 0 ) THEN
         PRINT*, 'ERREUR en ALLOCATION de vxyzmp(',NCAS,
     %           ')%vitx1t vity1t vitz1t vitm1t ou pres1t'
      ELSE
         PRINT*, 'ALLOCATION ERROR of vxyzmp(',NCAS,
     %           ')%vitx1t vity1t vitz1t vitm1t ou pres1t'
      ENDIF
      IERR = IALLOC


C     FIN DE L'EXECUTION
C     ==================
 9997 IF( NTDLTE .GT. 0 ) THEN
         if( allocated(temp) ) THEN
            do n = ncas1, ncas0, -1
               deallocate( temp( n )%dptab )
            enddo
            deallocate( temp )
         endif
      ENDIF

      if( allocated(pres) ) THEN
         do n = ncas1, ncas0, -1
            deallocate( pres( n )%dptab )
         enddo
         deallocate( pres )
      endif
      if( IALvitm .EQ. 0 ) THEN
         do n = ncas1, ncas0, -1
            deallocate( vitm( n )%dptab )
         enddo
         deallocate( vitm )
      endif

      if( IALvitz .EQ. 0 ) THEN
         do n = ncas1, ncas0, -1
            deallocate( vitz( n )%dptab )
         enddo
         deallocate( vitz )
      endif

      if( IALvity .EQ. 0 ) THEN
         do n = ncas1, ncas0, -1
            deallocate( vity( n )%dptab )
         enddo
         deallocate( vity )
      endif

      if( IALvitx .EQ. 0 ) THEN
         do n = ncas1, ncas0, -1
            deallocate( vitx( n )%dptab )
         enddo
         deallocate( vitx )
      endif

      IF( IALROTV   .EQ. 0 ) DEALLOCATE( ROTV   )
      IF( IALVPAUX  .EQ. 0 ) DEALLOCATE( VPAUX  )
      IF( IALPG     .EQ. 0 ) DEALLOCATE( PG     )
      IF( IALVXYZPN .EQ. 0 ) DEALLOCATE( VXYZPN )
      IF( IALTEMPER .EQ. 0 ) DEALLOCATE( TEMPER )

C     DESTRUCTION DU TABLEAU DES DONNEES DES PARTICULES
      IF( NBPART .GT. 0 .AND. MNPART .GT. 0 ) THEN
         CALL TNMCDS( 'REEL', 8*NBPART, MNPART )
         NBPART = 0
      ENDIF

      IF( MNEFSFTR.GT. 0 ) CALL TNMCDS( 'ENTIER', NBSFTR,   MNEFSFTR)
      IF( MNXYSFTR.GT. 0 ) CALL TNMCDS( 'ENTIER', NBSFTR,   MNXYSFTR)
      IF( MNNOSFTR.GT. 0 ) CALL TNMCDS( 'ENTIER', NBSFTR,   MNNOSFTR)
      IF( MNAGD   .GT. 0 ) CALL TNMCDS( 'REEL2',  NBSOMT,   MNAGD   )
      IF( MONOSO  .GT. 0 ) CALL TNMCDS( 'ENTIER', MONOSO,   MNNOSO  )
      IF( MOSONO  .GT. 0 ) CALL TNMCDS( 'ENTIER', MOSONO,   MNSONO  )
      IF( MNPTDG  .GT. 0 ) CALL TNMCDS( 'ENTIER', NBSOMT,   MNPTDG  )
      IF( MNLPCO  .GT. 0 ) CALL TNMCDS( 'ENTIER', NBCOPG,   MNLPCO  )
      IF( MNSTFR  .GT. 0 ) CALL TNMCDS( 'ENTIER', NBSOMT,   MNSTFR  )
      IF( MNAUX1  .GT. 0 ) CALL TNMCDS( 'REEL2',  MOAUX1,   MNAUX1  )
      IF( MNNPEF  .GT. 0 ) CALL TNMCDS( 'ENTIER', 2*MXTYEL, MNNPEF  )
      IF( MNTAUX  .GT. 0 ) CALL TNMCDS( 'REEL2',  MOTAUX,   MNTAUX  )
      IF( MNNDDL  .GT. 0 ) CALL TNMCDS( 'ENTIER', 1+NBNOVI, MNNDDL  )
      IF( MNTIMES .GT. 0 ) CALL TNMCDS( 'REEL',   MXVPFILE, MNTIMES )
      IF( MOXYZNP1.GT. 0 ) CALL TNMCDS( 'ENTIER', MOXYZNP1, MNXYZNP1)
      IF( MOELP1  .GT. 0 ) CALL TNMCDS( 'ENTIER', MOELP1,   MNELEP1 )
      IF( MOELP1  .GT. 0 ) CALL TNMCDS( 'ENTIER', 1,        MNNPEFP1)

C     RETOUR AUX PARAMETRES INITIAUX
      CALL XVEPAISSEUR( 0 )
      CALL XVTYPETRAIT( LIGCON )

      RETURN
      END
