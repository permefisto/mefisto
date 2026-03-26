      SUBROUTINE THENERGI( KNOMOB, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER L'ENERGIE CINETIQUE
C -----    ( CONDUC gradient temperature, gradient temperature )
C          PAR PRODUIT DE [K] avec {TEMPERATURES}
C          ET L'ENERGIE POTENTIELLE
C          Integrale temperature**2 (ROT C ) dX
C          PAR PRODUIT DE [M] avec {TEMPERATURES}
C
C          PAS DE CONDITIONS AUX LIMITES PRISES SUR [M] et [K]
C
C          POUR DES ELEMENTS FINIS LAGRANGE DE DEGRE 1 OU 2 EN 2D OU 3D
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET A CALCULER
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS         MAI 2006
C MODIFS: ALAIN PERRONNET  LJLL UPMC & SAINT PIERRE DU PERRAY  JUIN 2008
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
      CHARACTER*(*)     KNOMOB
      DOUBLE PRECISION  RELMIN,D2PI
      DOUBLE PRECISION  DINFO,DPREP,DMOPR
      INTEGER           NUMIOB(4),NUMAOB(4),MNDOEL(4),MXDOEL(4)
C
      DOUBLE PRECISION, allocatable, dimension(:) :: MG
      DOUBLE PRECISION, allocatable, dimension(:) :: KG
      INTEGER           IERMGALLOC, IERKGALLOC
      INTRINSIC         ALLOCATED
      DATA              RELMIN/-1D28/
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000) KNOMOB
      ELSE
         WRITE(IMPRIM,20000) KNOMOB
      ENDIF
10000 FORMAT(/90('=')/
     %'Calcul des ENERGIES CONNAISSANT les VECTEURS SOLUTIONS de l''OBJE
     %T: ',A/90('='))
20000 FORMAT(/87('=')/
     %'Computation of ENERGIES with KNOWN SOLUTION VECTORS of OBJECT: '
     %,A/87('='))
C
C     QUELQUES INITIALISATIONS
      TESTNL = 0
      DPREP = DINFO( 'CPU' )
      DMODR = 0D0
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
      MNMUKG = 0
      MNMUMG = 0
      DO 2 I=1,4
         MNDOEL(I) = 0
         MXDOEL(I) = 0
 2    CONTINUE
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
      NBCOOR = 0
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
C     RECHERCHE DES TEMPERATURES DE L'OBJET
      CALL  LXTSOU( NTLXOB, 'VECTEUR"TEMPERATURE', NTTEMP, MNTEMP )
      IF( NTTEMP .LE. 0 ) THEN
         CALL LXTSOU( NTLXOB, 'VECTEUR"VALEURPROPRE', NTTEMP, MNTEMP )
         IF( NTTEMP .LE. 0 ) THEN
C
C           RECHERCHE DES DEPLACEMENTS DE L'ONDE
            CALL  LXLXOU( NTLXOB, 'VECTEUR"DEPLACT', NTTEMP, MNTEMP )
            IF( NTTEMP .LE. 0 ) THEN
C              ERREUR PAS DE VECTEUR A VISUALISER
               NBLGRC(NRERR) = 2
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1)='ERREUR: OBJET '// KNOMOB
                  KERR(2)='SANS TEMPERATURE, DEPLACT OU VALEURS PROPRES'
               ELSE
                  KERR(1)='ERROR: OBJECT '// KNOMOB
                  KERR(2)='WITHOUT TEMPERATURE, DISPLACT or EIGENVALUES'
               ENDIF
               CALL LEREUR
               GOTO 9998
            ELSE
C              MODE DE TRACE DU DEPLACEMENT DE L'ONDE
               MODECO = 1
            ENDIF
         ELSE
C           MODE DE TRACE DES VALEURS ET VECTEURS PROPRES
            MODECO = 2
         ENDIF
      ELSE
C        MODE DE TRACE DES TEMPERATURES
         MODECO = 1
      ENDIF
      IF( NTTEMP .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='ERREUR: OBJET '// KNOMOB
            KERR(2)='SANS TEMPERATURE, DEPLACT OU VALEURS PROPRES'
         ELSE
            KERR(1)='ERROR: OBJECT '// KNOMOB
            KERR(2)='WITHOUT TEMPERATURE, DISPLACT or EIGENVALUES'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9998
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
C     NOMBRE (3 ou 6) DES COORDONNEES DES NOEUDS
      NBCOOR = MCN( MNXYZN + WBCOON )
C
C     OUVERTURE DES TABLEAUX DES DONNEES THERMIQUES DES PLSV DE L'OBJET
C     =================================================================
      CALL THEDON( NUMIOB, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             1,      MNDOEL,
     %             IEMAST, IECHMA, IECOND, IEDILA, IEVIFL, IECOET,
     %             IESOIN, IECONT, IEECHA, IESOCL, IESOPO,
     %             IETEIN, IEVIIN, IEVIANT,IECOBO,
     %             IERR )
C
C     BILAN DES DONNEES THERMIQUES
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
C     LE NOMBRE DE POINTS DU MAILLAGE DE L'OBJET
      NBPOI = MCN( MNXYZP + WNBPOI )
C     NDIM LA DIMENSION 2 OU 3 DE L'ESPACE DES COORDONNEES
      CALL DIMCOO( NBPOI, MCN(MNXYZP+WYZPOI), NDIM )
C     LE NOMBRE DE NOEUDS DU MAILLAGE DE L'OBJET=NOMBRE TOTAL DE DL
      NTDL = MCN( MNXYZN + WNBNOE )
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10210) NDIM, NTDL
      ELSE
         WRITE(IMPRIM,20210) NDIM, NTDL
      ENDIF
10210 FORMAT(/' DIMENSION 2 OU 3 OU 6 DE L''ESPACE',T42,'=',I6/
     %' NOMBRE DE COMPOSANTES DE CHAQUE VECTEUR',   T42,'=',I6/)
20210 FORMAT(/' SPACE DIMENSION (2 or 3 or 6)',T40,'=',I6/
     %' NUMBER of COMPONENTS of each VECTOR',  T40,'=',I6/)
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
C     LE TENSEUR DE CONDUCTIVITE,...
      CALL TNMCDC( 'REEL2', 128,MNTHER )
C
C     LE TABLEAU DES 3 COORDONNEES DES NBDLMX NOEUDS D'UN ELEMENT FINI
      CALL TNMCDC( 'REEL', NBDLMX*NBCOOR, MNX )
C
C     PREPARATION DE LA RESOLUTION STOCKAGE PROFIL DE K ET M
C     ============================ =========================
C     CALCUL DU PROFIL DES MATRICES ICI SYMETRIQUES NON DIAGONALES
      NCODSK = 1
      CALL TNMCDC( 'ENTIER', NTDL+1, MNMUKG )
      CALL PRPRMC( MNTOPO, MCN(MNNPEF), MNXYZN, 1, NCODSK,
     %             MCN(MNMUKG), IERR )
      IF( IERR .NE. 0 ) GOTO 9998
C
C     NOMBRE DE DOUBLE PRECISION DE CHACUNE DES 2 MATRICES PROFIL SYMETRIQUES
      NBRDKG = MCN( MNMUKG + NTDL )
C
C     DECLARATION DE LA MATRICE GLOBALE DE CAPACITE [MG] ET CONDUCTIVITE [KG]
C     =======================================================================
      IF( NBRDKG .LE. 0 ) THEN
C         PLACE MEMOIRE INSUFFISANTE
          NBLGRC(NRERR) = 3
          WRITE(KERR(MXLGER-1)(1:25),'(G25.0)') DMAT
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1) ='ERREUR: PLACE MEMOIRE INSUFFISANTE'
             KERR(2) = KERR(MXLGER-1)(1:25) //
     %               ' MOTS NECESSAIRES pour [M] [K]'
             KERR(3) = 'REDUIRE le MAILLAGE'
          ELSE
             KERR(1) ='ERROR: NOT ENOUGH MEMORY'
             KERR(2) = KERR(MXLGER-1)(1:25) //
     %               ' NECESSARY WORDS to store [M] [K]'
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
10290 FORMAT(' 2 MATRICES PROFIL CHACUNE DE',I15,
     %' REELS DOUBLE PRECISION'/
     %' 1/2 LARGEUR DE BANDE MOYENNE =',I9)
20290 FORMAT(' 2 SKYLINE MATRICES EACH of',I15,' DOUBLE REALS'/
     %' HALF WIDTH AVERAGE=',I9)
C
C     ALLOCATION DYNAMIQUE EN FORTRAN 90 DES MATRICES PROFIL MG et KG
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) 'ALLOCATION DEMAND  of',NBRDKG,
     %                ' DOUBLE PRECISION of the [MG] and [KG] MATRICES'
      ALLOCATE ( MG(1:NBRDKG), STAT=IERMGALLOC )
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
      CALL THEMKB( 1,      IEMG,   IEKG,   IEBG,   PENALI,
     &             D2PI,   NDIM,   NTDL,
     &             NBTYEL, MNNPEF, NDPGST,
     &             MNTPOB, MXPOBA, MNTAUX,
     &             MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     &             MNTHER, MNTAEL, MNX,    MNNODL,
     &             1,      MNMUKG, 0,
     &             NBRDKG, MG,     NBRDKG, KG,     MNBG,
     &             NCODSM, NCODSK, NBPTAF, IERR )
C
C     COPIE DE MU DE KG DANS MU DE MG
      CALL TNMCDC( 'ENTIER', 1+NTDL, MNMUMG )
      CALL TRTATA( MCN(MNMUKG), MCN(MNMUMG), 1+NTDL )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)
     %   'CONSTANTES MATRICES [M] & [K] CALCULEES'
      ELSE
         WRITE(IMPRIM,*)
     %   'CONSTANT MATRICES [M] & [K] COMPUTED'
      ENDIF
C
C     TEMPS CALCUL DE K ET M
      DPREP = DINFO( 'DELTA CPU' )
C
C     NOMBRE DE DEGRES DE LIBERTE
      NTDL0 = MCN( MNTEMP + WBCOVE )
      IF( NTDL0 .NE. NTDL ) THEN
         NBLGRC(NRERR) = 4
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET ' // KNOMOB
            KERR(2) = 'NOMBRE DES NOEUDS DU MAILLAGE NON EGAL AU'
            KERR(3) = 'NOMBRE DES COMPOSANTES DES VECTEURS'
            KERR(4) = 'RESOUDRE A NOUVEAU LE PROBLEME'
         ELSE
            KERR(1) = 'ERROR: OBJECT ' // KNOMOB
            KERR(2) = 'NUMBER of MESH NODES NOT EQUAL to'
            KERR(3) = 'NUMBER of SOLUTION VECTOR COMPONENTS'
            KERR(4) = 'SOLVE AGAIN the PROBLEM'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9998
      ENDIF
C
C     NOMBRE TOTAL DE CAS
      NDSM = MCN( MNTEMP + WBVECT )
      IF( NDSM .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: PAS DE VECTEUR STOCKE'
         ELSE
            KERR(1) = 'ERROR: NO VECTOR STORED'
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9999
      ENDIF
C
C     L'ADRESSE DES TEMPS SI PB INSTATIONNAIRE
      MNTIME = MNTEMP + WECTEU + NTDL * NDSM * MOREE2 - 1
C
C     SI PB INSTATIONNAIRE LE DERNIER VECTEUR EST CHOISI
      IF( MCN(MNTEMP + WBCPIN) .GT. 0 ) THEN
C        LE TEMPS FINAL
         NCAS  = NDSM
         TEMPS = RMCN(MNTIME+NDSM)
      ENDIF
C
C     PRODUIT DES MATRICES K et M avec les TEMPERATURES CALCULEES
      MNENER0 = 0
      CALL TNMCDC( 'REEL2', NDSM * 2, MNENER0 )
      MNENER1 = MNENER0 + NDSM*MOREE2
      MNAUX = 0
      CALL TNMCDC( 'REEL2', NTDL, MNAUX )
C
C     ( KG * VECTEUR, VECTEUR ) => ENERGIE CINETIQUE
      CALL ENERGY( NTDL, NDSM,
     %             NCODSK, MCN(MNMUKG), KG,
     %             MCN(MNTEMP+WECTEU),  MCN(MNAUX),
     %             MCN(MNENER0) )
C
C     ( MG * VECTEUR, VECTEUR ) => ENERGIE POTENTIELLE
      CALL ENERGY( NTDL,   NDSM,
     %             NCODSM, MCN(MNMUMG), MG,
     %             MCN(MNTEMP+WECTEU),  MCN(MNAUX),
     %             MCN(MNENER1) )
C
C     AFFICHAGE DES ENERGIES CINETIQUE ET POTENTIELLE
      WRITE(IMPRIM,*)
      MN0 = ( MNENER0 - 1 ) / 2
      MN1 = ( MNENER1 - 1 ) / 2
      DO 600 I=1,NDSM
         IF( MCN(MNTEMP + WBCPIN) .GT. 0 ) THEN
            VAP = RMCN(MNTIME+I)
         ELSE
            VAP = 0
         ENDIF
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10600) I, VAP, DMCN(MN0+I), DMCN(MN1+I),
     %                          DMCN(MN0+I)/DMCN(MN1+I)
         ELSE
            WRITE(IMPRIM,20600) I, VAP, DMCN(MN0+I), DMCN(MN1+I),
     %                          DMCN(MN1+I)/DMCN(MN0+I)
         ENDIF
 600  CONTINUE
10600 FORMAT('VALEUR PROPRE',I5,' =',G14.6,
     %     '  ENERGIE CINETIQUE KE=',G14.6,
     %     '  ENERGIE POTENTIELLE PE=',G14.6,
     %     '  PE/KE=',G14.6)
20600 FORMAT('EIGENVALUE',I5,' =',G14.6,
     %     '  KINEMATIC ENERGY KE=',G14.6,
     %     '  POTENTIAL ENERGY PE=',G14.6,
     %     '  PE/KE=',G14.6)
C
cccC     POUR GOONG CHEN cas NIVEAUX D'ENERGIES DE L'HYDROGEN
cccC     AMELIORATION PAR LE VIRIAL THEOREME Cf Extrapolation_rev.pdf
cccC     ------------------------------------------------------------
ccc      NU(2) = 1D0 / SQRT(5D0) pour Hydrogen
ccc      NU(1) = 2D0 * NU(2)
cccC     POUR -1/2 Delta P + r2/2 P = E P
ccc      NU(1) = 1D0 / SQRT(2D0)
ccc      NU(2) =-NU(1)
ccc      MN0 = ( MNENER0 - 1 ) / 2
ccc      MN1 = ( MNENER1 - 1 ) / 2
ccc      WRITE(IMPRIM,*)
ccc      WRITE(IMPRIM,*)
ccc     %'Computation -1/2 Delta P + r2/2 P = E P  Energy levels and improv
ccc     %ements'
cccccc     %'Computation of Hydrogen energy levels and improvements'
ccc      DO 640 I=1,NDSM
cccC        I-th COMPUTED EIGENVALUE
ccc         VAP = RMCN(MNTIME+I)
cccC        I-th KINEMATIC ENERGY
ccc         KE  = DMCN(MN0+I)
cccC        I-th POTENTIAL ENERGY
ccc         PE  = DMCN(MN1+I)
cccC        NOTATION of the ARTICLE
ccc         XH(1) = KE
ccc         XH(2) = PE
ccc         EHE = NU(1) * XH(1) + NU(2) * XH(2)
ccc         XHE(1) = XH(1) - EHE * NU(1)
ccc         XHE(2) = XH(2) - EHE * NU(2)
ccc         WRITE(IMPRIM,10640) I, VAP, KE, PE, PE/KE,
ccc     %                   I, XHE(1)+XHE(2), XHE(1), XHE(2), XHE(2)/XHE(1)
ccc 640  CONTINUE
ccc10640 FORMAT(
ccc     %'Eh',I2,'=',G14.6,' (KEh=',G14.6,',PEh=',G14.6,') PEh/KEh=',G14.6/
ccc     %'Ee',I2,'=',G14.6,' (KEe=',G14.6,',PEe=',G14.6,') PEe/KEe=',G14.6/
ccc     %)
cccC
ccc
ccc      DO 601 I=1,NDSM
ccc         IF( MCN(MNTEMP + WBCPIN) .GT. 0 ) THEN
ccc            VAP = RMCN(MNTIME+I)
ccc         ELSE
ccc            VAP = 0
ccc         ENDIF
ccc         WRITE(IMPRIM,*) '    %',VAP,', ', DMCN(MN0+I),', ',DMCN(MN1+I)
ccc 601  CONTINUE
ccc
cccC
cccC     POUR GOONG CHEN cas HeH++: VERIFICATION KE + PE - VAP + 2/R = 0
cccC     ---------------------------------------------------------------
ccc      WRITE(IMPRIM,*)
ccc      MN0 = ( MNENER0 - 1 ) / 2
ccc      MN1 = ( MNENER1 - 1 ) / 2
ccc      DO 630 I=1,NDSM
ccc         IF( MCN(MNTEMP + WBCPIN) .GT. 0 ) THEN
ccc            VAP = RMCN(MNTIME+I)
ccc         ELSE
ccc            VAP = 0
ccc         ENDIF
ccc         IF( LANGAG .EQ. 0 ) THEN
ccc            WRITE(IMPRIM,10630) I, VAP,
ccc     %                          DMCN(MN0+I)+DMCN(MN1+I)-VAP+2D0,
ccc     %                          DMCN(MN0+I)/DMCN(MN1+I)
ccc         ELSE
ccc            WRITE(IMPRIM,20630) I, VAP,
ccc     %                          DMCN(MN0+I)+DMCN(MN1+I)-VAP+2D0,
ccc     %                          DMCN(MN0+I)/DMCN(MN1+I)
ccc         ENDIF
ccc 630  CONTINUE
ccc10630 FORMAT('VALEUR PROPRE',I5,' =',G15.7,
ccc     %       '  KE+PE-VP+2/R=',G15.7,'   KE/PE=',G15.7)
ccc20630 FORMAT('EIGENVALUE',I5,' =',G15.7,
ccc     %       '  KE+PE-VP+2/R=',G15.7,'   KE/PE=',G15.7)
cccC
ccc
ccc
cccC     POUR GOONG CHEN cas 1E-6 * r**4 POTENTIAL: CALCUL PE/KE=1/2
cccC     -----------------------------------------------------------
ccc      WRITE(IMPRIM,*)
ccc      MN0 = ( MNENER0 - 1 ) / 2
ccc      MN1 = ( MNENER1 - 1 ) / 2
ccc      DO 630 I=1,NDSM
ccc         IF( MCN(MNTEMP + WBCPIN) .GT. 0 ) THEN
ccc            VAP = RMCN(MNTIME+I)
ccc         ELSE
ccc            VAP = 0
ccc         ENDIF
ccc         E = VAP
ccc         KE = DMCN(MN0+I)
ccc         PE = DMCN(MN1+I)
ccc         XI1 = KE + ( -KE + 2D0*PE ) / 5D0
ccc         XI2 = PE - ( -KE + 2D0*PE ) * 2D0 / 5D0
ccc         EI  = XI1 + XI2
ccc         WRITE(IMPRIM,10630) I, E, KE, PE, PE/KE, XI1, XI2, EI
ccc 630  CONTINUE
ccc10630 FORMAT('Eh',I2,'=',F7.5,'  Xh=(',F7.5,',',F7.5,')  PEh/KEh=',F7.5,
ccc     %'  XhI=(',F7.5,',',F7.5,')  EhI=',F7.5 )
ccc
ccc
cccC     POUR GOONG CHEN cas MORSE POTENTIAL: VERIFICATION KE + PE - VAP = 0
cccC     -------------------------------------------------------------------
ccc         WRITE(IMPRIM,10630) I, VAP, DMCN(MN0+I), DMCN(MN1+I),
ccc     %                       DMCN(MN0+I)+DMCN(MN1+I)-VAP
ccc 630  CONTINUE
ccc10630 FORMAT('E',I3,'=',G15.7,'  KE=',G15.7,'  PE=',G15.7,
ccc     %       '  KE+PE-E=',G15.7)
ccc
      CALL TNMCDS( 'REEL2', NTDL, MNAUX )
      CALL TNMCDS( 'REEL2', NDSM * 2, MNENER0 )
      IF( IERR .NE. 0 ) GOTO 9998
C
C     FERMETURE DES MATRICES POUR REDONNER DE LA PLACE EN MC
C     ======================================================
 9998 IF( IERKGALLOC .EQ. 0 ) DEALLOCATE( KG )
      IF( IERMGALLOC .EQ. 0 ) DEALLOCATE( MG )
      IF(MNMUMG .GT. 0) CALL TNMCDS( 'ENTIER', NTDL+1, MNMUMG )
      IF(MNMUKG .GT. 0) CALL TNMCDS( 'ENTIER', NTDL+1, MNMUKG )
      IF(IERR   .NE. 0) GOTO 9999
C
C **************************************************************************
C --------------------------------------------------------------------------
C **************************************************************************
C
C     DESTRUCTION DES TABLEAUX TEMPORAIRES
C     ====================================
9999  IF( MNNPEF .GT. 0 ) CALL TNMCDS( 'ENTIER',2*MXTYEL, MNNPEF )
      DO 11000 I=1,4
         IF( MNDOEL(I) .GT. 0 .AND. MXDOEL(I) .GT. 0 ) THEN
            CALL TNMCDS( 'ENTIER', MXDOEL(I), MNDOEL(I) )
         ENDIF
11000 CONTINUE
      IF( MNTPOB .GT. 0 ) CALL TNMCDS( 'ENTIER', NBTYEL*MXPOBA,MNTPOB )
      IF( MNTAUX .GT. 0 ) CALL TNMCDS( 'REEL2' , MOAUX , MNTAUX )
      IF( MNTAEL .GT. 0 ) CALL TNMCDS( 'REEL2' , MOTAEL, MNTAEL )
      IF( MNNODL .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLMX, MNNODL )
      IF( MNTHER .GT. 0 ) CALL TNMCDS( 'REEL2' , 128   , MNTHER )
      IF( MNX    .GT. 0 ) CALL TNMCDS( 'REEL',   NBDLMX*NBCOOR, MNX  )
      IF( MNNDLX .GT. 0 ) CALL TNMCDS( 'ENTIER', MONDLX, MNNDLX )
      IF( MNVDLX .GT. 0 ) CALL TNMCDS( 'REEL2',  MONDLX, MNVDLX )
C
C     GESTION DES ERREURS
C     ===================
      IF( IERR .EQ. 7 ) THEN
C        RETOUR SI MATRICE NON INVERSIBLE
         IERR = 0
         RETURN
      ENDIF
C
C     COUT CALCUL DES ENERGIES
C     ========================
      DMOPR = DINFO( 'DELTA CPU' )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,12001) DPREP,DMOPR,DPREP+DMOPR
      ELSE
         WRITE(IMPRIM,22001) DPREP,DMOPR,DPREP+DMOPR
      ENDIF
12001 FORMAT(/
     %' TEMPS CALCUL DE LA PREPARATION =',D15.6,' SECONDES CPU'/
     %' TEMPS CALCUL DES ENERGIES      =',D15.6,' SECONDES CPU'/
     %' TEMPS CALCUL TOTAL             =',D15.6,' SECONDES CPU')
22001 FORMAT(/
     %' PREPARATION TIME =',F12.2,' CPU SECONDS'/
     %' ENERGIES    TIME =',F12.2,' CPU SECONDS'/
     %' TOTAL       TIME =',F12.2,' CPU SECONDS')
C
      RETURN
      END
