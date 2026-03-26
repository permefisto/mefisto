      SUBROUTINE THEVVP( KNOMOB, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES VALEURS ET VECTEURS PROPRES DE L'OPERATEUR
C -----    DE LA CHALEUR LINEAIRE POUR DES ELEMENTS FINIS LAGRANGE DE
C          DEGRE 1 OU 2 EN 2D OU 3D OU 6D
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET DE VALEURS ET VECTEURS PROPRES A CALCULER
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1998
C MODIFS : ALAIN PERRONNET TEXAS A & M UNIVERSITY           JUILLET 2003
C MODIFS : ALAIN PERRONNET TEXAS A & M UNIVERSITY           JUILLET 2005
C MODIFS : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR  FEVRIER 2011
C23456---------------------------------------------------------------012
      DOUBLE PRECISION   PENALI
      PARAMETER         (PENALI=0D0)
      PARAMETER         (MXTYEL=7)
      PARAMETER         (MOPAGE=512)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/donthe.inc"
      include"./incl/donele.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___contact.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___profil.inc"
      include"./incl/a___fluxpt.inc"
      include"./incl/a___tableau1r.inc"
      include"./incl/msvaau.inc"
      include"./incl/ponoel.inc"
      include"./incl/cthet.inc"
      include"./incl/ctemps.inc"
      include"./incl/cnonlin.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      EXTERNAL          ETTAEL
      DOUBLE PRECISION  RELMIN, D2PI, EIGMIN, EIGMAX
      DOUBLE PRECISION  DINFO,  DCPU, DPREP, DMOPR, DMAT
      INTEGER           NUMIOB(4), NUMAOB(4), MNDOEL(4), MXDOEL(4)
C
      DOUBLE PRECISION, allocatable, dimension(:) :: MG
      DOUBLE PRECISION, allocatable, dimension(:) :: KG
      INTEGER           IERMGALLOC, IERKGALLOC
      INTRINSIC         ALLOCATED
C
      CHARACTER*(*)     KNOMOB
      DATA              RELMIN/-1D28/
C     TABLEAU NON UTILISE (VITESSE D'UN FLUIDE NON ICI CALCULE)
      DOUBLE PRECISION  VITEGt(5,3)
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000) KNOMOB
      ELSE
         WRITE(IMPRIM,20000) KNOMOB
      ENDIF
10000 FORMAT(/1X,80('=')/
     %' RESOLUTION des VALEURS et VECTEURS PROPRES de l''OBJET: ',A/
     %1X,80('='))
20000 FORMAT(/1X,80('=')/
     %' SOLUTION of EIGENVALUES and EIGENVECTORS of OBJECT: ',A/
     %1X,80('='))
C
C     QUELQUES INITIALISATIONS
      NBJEUX = 1
      TESTNL = 0
      NBCOOR = 0
      DCPU   = DINFO( 'CPU' )
      DPREP  = 0D0
      DMODR  = 0D0
C     TEMPS POUR LA THERMIQUE
      TEMPS  = 0.0
      IERR   = 0
C     2 PI
      D2PI   = ATAN( 1D0 ) * 8D0
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
C
C     PROTECTION DES ADRESSES POUR EVITER DES PROBLEMES LORS
C     DE LA DESTRUCTION DES TABLEAUX
      MNNPEF = 0
      MNTPOB = 0
      MNTAUX = 0
      MNTAEL = 0
      MNNODL = 0
      MNTHER = 0
      MNX    = 0
      MNMUMG = 0
      MNMUKG = 0
      DO I=1,4
         MNDOEL(I) = 0
         MXDOEL(I) = 0
      ENDDO
      NBTYEL = 0
      MOAUX  = 0
      MOTAEL = 0
      NBDLMX = 0
      NTDL   = 0
      MOFLTO = 0
      MOFLPT = 0
      MONDLX = 0
      MNNDLX = 0
      MNVDLX = 0
      NDSM   = 1
      IERMGALLOC = 1
      IERKGALLOC = 1
C
C     AFFICHAGE ET VERIFICATION DU NOM_DE_L'OBJET
C     ===========================================
C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
      CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET INCONNU ' // KNOMOB
         ELSE
            KERR(1) = 'ERROR: UNKNOWN OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         CALL LXIM( NTOBJE )
         IERR = 1
         RETURN
      ENDIF
C
C     RECHERCHE DU TABLEAU DEFINITION
      CALL LXTSOU( NTLXOB,'DEFINITION', NTDFOB, MNDFOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='ERREUR: DEFINITION INCONNUE de l''OBJET ' //KNOMOB
         ELSE
            KERR(1) ='ERROR: UNKNOWN DEFINITION for the OBJECT '//KNOMOB
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     ============================================
C     CALCUL DES MATRICES DE MASSE ET CONDUCTIVITE
C     ============================================
C     RECHERCHE DES TABLEAUX XYZSOMMET XYZNOEUD XYZPOINT NPEF"
C     ASSOCIES A L'OBJET
      CALL MIMAOB( 1,      NTLXOB, MXDOTH, NTTOPO, MNTOPO,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %             MXTYEL, NBTYEL, MTNPEF, MNNPEF,
     %             NUMIOB, NUMAOB,
     %             NDPGST, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             MXDOEL, MNDOEL, IERR )
      IF( IERR .NE. 0 ) GOTO 9999
C     ICI NUMIOB CONTIENT LES 4 NUMEROS MINIMA DES OBJETS
C         NUMAOB          LES 4 NUMEROS MAXIMA DES OBJETS
C     MNDOEL LES 4 ADRESSES MCN DES TABLEAUX DES ADRESSES DES
C     TABLEAUX DECRIVANT LA THERMIQUE DE L'OBJET COMPLET
C
C     OUVERTURE DES TABLEAUX DES DONNEES THERMIQUES DES PLSV DE L'OBJET
C     =================================================================
      CALL THEDON( NUMIOB, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             NBJEUX, MNDOEL,
     %             IEMAST, IECHMA, IECOND, IEDILA, IEVIFL, IECOET,
     %             IESOIN, IECONT, IEECHA, IESOCL, IESOPO,
     %             IETEIN, IEVIIN, IEVIANT,IECOBO,
     %             IERR )
C
CCCC     BILAN DES DONNEES THERMIQUES
CCC      IF( IECONT .EQ. 0 ) THEN
CCC         NBLGRC(NRERR) = 2
CCC         KERR(1) = 'ERREUR AUCUN CONTACT SUR L''OBJET '// KNOMOB
CCC         KERR(2) = 'AJOUTER DES CONTACTS'
CCC         CALL LEREUR
CCC
CCC      ENDIF
      IF( IERR .GT. 0 ) THEN
         IERR = 4
         GOTO 9999
      ENDIF
C
C     INITIALISATIONS DE VARIABLES ET AFFICHAGES
C     ==========================================
C     NBCOOR : NOMBRE DE COORDONNEES DES NOEUDS=POINTS (3 ou 6)
      NBCOOR = MCN( MNXYZP + WBCOOP )
C
C     LE NOMBRE DE POINTS DU MAILLAGE DE L'OBJET
      NBPOI = MCN( MNXYZP + WNBPOI )
C     NDIM LA DIMENSION 2 OU 3 OU 6 DE L'ESPACE DES COORDONNEES
      CALL DIMCOO( NBPOI, MCN(MNXYZP+WYZPOI), NDIM )
C     PARTICULARITE DES 6-CUBES (3Q1C NON A TRAITER EN THERMIQUE)
      IF( NBCOOR .EQ. 6 ) THEN
         NBTYEL = 1
         NDIM   = NBCOOR
      ENDIF
C
C     LE NOMBRE DE NOEUDS DU MAILLAGE DE L'OBJET=NOMBRE TOTAL DE DEGRES DE LIBER
      NTDL = MCN( MNXYZN + WNBNOE )
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10210) NDIM,NTDL
      ELSE
         WRITE(IMPRIM,20210) NDIM,NTDL
      ENDIF
10210 FORMAT(/' DIMENSION 2 OU 3 OU 6 DE L''ESPACE',T42,'=',I6/
     %' NOMBRE DE COMPOSANTES DE CHAQUE VECTEUR',   T42,'=',I6)
20210 FORMAT(/' SPACE DIMENSION (2 or 3 or 6)',T40,'=',I6/
     %' NUMBER of COMPONENTS of each VECTOR',  T40,'=',I6)
C
C     RECUPERATION DES TABLEAUX POBA NECESSAIRES A LA
C     CONSTRUCTION DES TABLEAUX ELEMENTAIRES
C     ===============================================
      CALL TAPOBA( NBTYEL, MNNPEF, ETTAEL,
     %             MNTPOB, NBDLMX, MOAUX, NBTTEF, NOAXIS, NCODSM, IERR )
      IF( IERR .NE. 0 ) GOTO 9999
C
C     ADRESSAGE DES TABLEAUX AUXILIAIRES ET ELEMENTAIRES
C     ===================================================
      CALL TNMCDC( 'REEL2', MOAUX, MNTAUX )
C
C     LES 2 MATRICES ELEMENTAIRES ET LES NDSM VECTEURS ELEMENTAIRES
      MOTAEL = NBDLMX * (NBDLMX+1) + NBDLMX * NDSM
      CALL TNMCDC( 'REEL2', MOTAEL, MNTAEL )
C
C     LE NUMERO DES DEGRES DE LIBERTE GLOBAUX DES DL D'UN EF
      CALL TNMCDC( 'ENTIER', NBDLMX, MNNODL )
C
C     LE TENSEUR DE CONDUCTIVITE ou COEFFICIENT TEMPERATURE(Pts INTEGRATION)
      CALL TNMCDC( 'REEL2', 128, MNTHER )
C
C     LE TABLEAU DES NBCOOR COORDONNEES DES NBDLMX NOEUDS D'UN ELEMENT FINI
      CALL TNMCDC( 'REEL', NBDLMX*NBCOOR, MNX )
C
C     PREPARATION DE LA RESOLUTION STOCKAGE PROFIL DE K ET M
C     ============================ =========================
C     CALCUL DU PROFIL DES MATRICES ICI SYMETRIQUES NON DIAGONALES
      NCODSK = 1
      CALL TNMCDC( 'ENTIER', 1+NTDL, MNMUKG )
      CALL PRPRMC( MNTOPO, MCN(MNNPEF), MNXYZN, 1, NCODSK,
     %             MCN(MNMUKG), IERR )
      IF( IERR .GT. 0 ) GOTO 9999
C
C     COPIE DU TABLEAU POINTEUR (ENSUITE PEUT ETRE MODIFIE PAR LES CL)
      CALL TNMCDC( 'ENTIER', 1+NTDL, MNMUMG )
      CALL TRTATA( MCN(MNMUKG), MCN(MNMUMG), 1+NTDL )
C
C     DECLARATION DES 2 MATRICES PROFIL SYMETRIQUE
C    ( +1 DANS CALVVP POUR LA MATRICE DE PROTECTION APRES CL => PLUS PETITE)
      NBRDMG = MCN( MNMUMG + NTDL )
      NBRDKG = MCN( MNMUKG + NTDL )
      DMAT   = NBRDKG * 3D0
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
10290 FORMAT(' 3 MATRICES PROFIL CHACUNE DE',I15,
     %' REELS DOUBLE PRECISION'/
     %' 1/2 LARGEUR DE BANDE MOYENNE =',I9)
20290 FORMAT(' 3 SKYLINE MATRICES EACH of',I15,' DOUBLE REALS'/
     %' HALF WIDTH AVERAGE=',I9)
C
C     ALLOCATION DYNAMIQUE EN FORTRAN 90 DES MATRICES PROFIL MG et KG
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) 'ALLOCATION DEMAND  of',NBRDKG,
     %                ' DOUBLE PRECISION of the [MG] and [KG] MATRICES'
      ALLOCATE ( MG(1:NBRDMG), STAT=IERMGALLOC )
      IF( IERMGALLOC .NE. 0 ) THEN
       WRITE(IMPRIM,*)'ALLOCATION ERROR   of',NBRDKG,
     %                ' DOUBLE PRECISION of the [MG] MATRIX'
         IERR = IERMGALLOC
         GOTO 9999
      ENDIF
      ALLOCATE ( KG(1:NBRDKG), STAT=IERKGALLOC )
      IF( IERKGALLOC .NE. 0 ) THEN
       WRITE(IMPRIM,*)'ALLOCATION ERROR   of',NBRDKG,
     %                ' DOUBLE PRECISION of the [KG] MATRIX'
         IERR = IERKGALLOC
         GOTO 9999
      ENDIF
      WRITE(IMPRIM,*) 'ALLOCATION CORRECT of',NBRDKG,
     %                ' DOUBLE PRECISION of the [MG] and [KG] MATRICES'
      WRITE(IMPRIM,*)
C
C     KG & MG SYMETRIQUES NON DIAGONALES PROFIL
      NCODSK = 1
      NCODSM = 1
C
C     ICI: PAS DE SECOND MEMBRE GLOBAL
      MNBG = 0
C
C     CALCUL DES MATRICES DE CAPACITE ET CONDUCTIVITE SUPPOSEES
C     INDEPENDANTES DU TEMPS ET DE LA TEMPERATURE
C     ---------------------------------------------------------
      IEMG = 1
      IEKG = 1
      IEBG = 0
      CALL THEMKB( NBJEUX, IEMG,   IEKG,   IEBG,   PENALI,
     &             D2PI,   NDIM,   NTDL,   VITEGt,
     &             NBTYEL, MNNPEF, NDPGST,
     &             MNTPOB, MXPOBA, MNTAUX,
     &             MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     &             MNTHER, MNTAEL, MNX,    MNNODL,
     &             1,      MNMUKG, 0,
     &             NBRDMG, MG,     NBRDKG, KG,   MNBG,
     &             NCODSM, NCODSK, NBPTAF, IERR )
      IF( IERR .NE. 0 ) GOTO 9999
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)
     %      'MATRICES [M] et [K] CONSTANTES CALCULEES'
      ELSE
         WRITE(IMPRIM,*)
     %      '[M] and [K] CONSTANT MATRICES COMPUTED'
      ENDIF
C
C     CONSTRUCTION DES TABLEAUX DU NUMERO ET VALEUR DES DL FIXES
C     ----------------------------------------------------------
      CALL THDLFX( 0,      NTDL,   NDIM,
     &             NBTYEL, MNNPEF, NDPGST,
     &             MNXYZN, NUMIOB, MNDOEL, RELMIN,
     &             NBDLFX, MONDLX, MNNDLX, MNVDLX, IERR )
      IF( IERR .NE. 0 ) GOTO 9999
C
C     TEMPS CALCUL DE K ET M
      DPREP = DINFO( 'DELTA CPU' )
C
C     CALCUL DES VALEURS EIGV ET VECTEURS PROPRES VP TELS QUE
C     ( K + EIGV M ) V = 0
C     =======================================================
      METVVP = 0
      CALL CALVVP( 'THERMIQUE', NTLXOB, NTDL,
     %             METVVP, EIGMIN, EIGMAX, NBROOT,
     %             NBDLFX, MNNDLX, MNVDLX,
     %             NCODSK, MNMUKG, KG,
     %             NCODSM, MNMUMG, MG,
     %             NTVVPR, MNVVPR, IERR )
C     LE NOMBRE DE VALEURS PROPRES
      NBVALP = MCN(MNVVPR+WBVECT)
C
C     DEALLOCATION DES MATRICES GLOBALES POUR REDONNER DE LA PLACE
C     ============================================================
      IF(MNMUMG .GT. 0) CALL TNMCDS( 'ENTIER', NTDL+1, MNMUMG )
      IF(MNMUKG .GT. 0) CALL TNMCDS( 'ENTIER', NTDL+1, MNMUKG )
C     FREE MEMORY USED BY ARRAY KG and MG
      IF( IERKGALLOC .EQ. 0 ) THEN
         DEALLOCATE( KG )
         IERKGALLOC = 1
      ENDIF
      IF( IERMGALLOC .EQ. 0 ) THEN
         DEALLOCATE( MG )
         IERMGALLOC = 1
      ENDIF
C
C     COUT CALCUL DES VALEURS ET VECTEURS PROPRES
C     ===========================================
      DMOPR = DINFO( 'DELTA CPU' )
      IF( IERR .NE. 0 ) GOTO 9999
C
C **************************************************************************
C --------------------------------------------------------------------------
C **************************************************************************
C
C     CALCUL DES FLUX EN CHAQUE POINT D INTEGRATION DES FACES DE
C     CHAQUE ELEMENT FINI DE CHAQUE TYPE D'EF DE L'OBJET
C     ==========================================================
      IF( NBCOOR .LE. 3 ) THEN
         CALL THEFLU( KNOMOB, NTLXOB, MNTOPO, NOAXIS, D2PI,
     %                NDIM,   MOREE2, NBVALP, NTDL,
     %                NBTYEL, MNNPEF, NDPGST, MNTPOB,
     %                MNTAUX, MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     %                MOTAEL, MNTAEL, MNX,    MNVVPR )
      ENDIF
C
C     DESTRUCTION DES TABLEAUX TEMPORAIRES
C     ====================================
C     FREE MEMORY USED BY ARRAYS KG & MG
 9999 IF( IERKGALLOC .EQ. 0 ) THEN
         DEALLOCATE( KG )
         IERKGALLOC = 1
      ENDIF
      IF( IERMGALLOC .EQ. 0 ) THEN
         DEALLOCATE( MG )
         IERMGALLOC = 1
      ENDIF
C
      IF( MNNPEF .GT. 0 ) CALL TNMCDS( 'ENTIER',2*MXTYEL, MNNPEF )
      DO I=1,4
         IF( MNDOEL(I) .GT. 0 .AND. MXDOEL(I) .GT. 0 ) THEN
            CALL TNMCDS( 'ENTIER', MXDOEL(I), MNDOEL(I) )
         ENDIF
      ENDDO
      IF( MNTPOB .GT. 0 ) CALL TNMCDS( 'ENTIER', NBTYEL*MXPOBA,MNTPOB )
      IF( MNTAUX .GT. 0 ) CALL TNMCDS( 'REEL2' , MOAUX , MNTAUX )
      IF( MNTAEL .GT. 0 ) CALL TNMCDS( 'REEL2' , MOTAEL, MNTAEL )
      IF( MNNODL .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLMX, MNNODL )
      IF( MNTHER .GT. 0 ) CALL TNMCDS( 'REEL2' , 128   , MNTHER )
      IF( MNX    .GT. 0 ) CALL TNMCDS( 'REEL',   NBDLMX*NBCOOR, MNX )
      IF( MNNDLX .GT. 0 ) CALL TNMCDS( 'ENTIER', MONDLX, MNNDLX )
      IF( MNVDLX .GT. 0 ) CALL TNMCDS( 'REEL2',  MONDLX, MNVDLX )
      IF( MNMUMG .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL+1, MNMUMG )
      IF( MNMUKG .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL+1, MNMUKG )
C
C     GESTION DES ERREURS
C     ===================
      IF( IERR .EQ. 7 ) THEN
C        RETOUR SI MATRICE NON INVERSIBLE
         IERR = 0
         RETURN
      ENDIF
      IF( IERR .NE. 0 ) RETURN
C
C     COUT CALCUL DES FLUX DE TEMPERATURE
      DCPU = DINFO( 'DELTA CPU' )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,12001) DPREP,DMOPR,DCPU,DPREP+DMOPR+DCPU
      ELSE
         WRITE(IMPRIM,22001) DPREP,DMOPR,DCPU,DPREP+DMOPR+DCPU
      ENDIF
12001 FORMAT(/
     %' TEMPS CALCUL DE LA PREPARATION   =',F12.2,' SECONDES CPU'/
     %' TEMPS CALCUL DES VALEURS PROPRES =',F12.2,' SECONDES CPU'/
     %' TEMPS CALCUL DES FLUX            =',F12.2,' SECONDES CPU'/
     %' TEMPS CALCUL TOTAL               =',F12.2,' SECONDES CPU')
22001 FORMAT(/
     %' PREPARATION TIME =',F12.2,' CPU SECONDS'/
     %' EIGENVALUES TIME =',F12.2,' CPU SECONDS'/
     %' NORMAL FLUX TIME =',F12.2,' CPU SECONDS'/
     %' TOTAL       TIME =',F12.2,' CPU SECONDS')
C
      RETURN
      END
