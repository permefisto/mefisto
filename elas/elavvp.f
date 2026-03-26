      SUBROUTINE ELAVVP( KNOMOB, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES MODES PROPRES (VALEURS ET VECTEURS PROPRES)
C -----    EN ELASTICITE LINEAIRE
C          POUR DES ELEMENTS FINIS LAGRANGE DE DEGRE 1 OU 2 EN 2D OU 3D
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET DE MODES PROPRES A CALCULER
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET      ANALYSE NUMERIQUE UPMC PARIS    AOUT 1998
C23456---------------------------------------------------------------012
      PARAMETER         (MXTYEL=7)
      PARAMETER         (MOPAGE=512)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/donela.inc"
      include"./incl/donele.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___fixation.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___profil.inc"
      include"./incl/a___contrainte.inc"
      include"./incl/a___tableau1r.inc"
      include"./incl/msvaau.inc"
      include"./incl/homdir.inc"
      include"./incl/ponoel.inc"
      include"./incl/ctemps.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / ELFINI / NELFI(512)
C
      DOUBLE PRECISION, allocatable, dimension(:) :: MG
      DOUBLE PRECISION, allocatable, dimension(:) :: KG
      INTEGER           IERMGALLOC, IERKGALLOC
      INTRINSIC         ALLOCATED
CC
      INTEGER           NOMTAB(5)
      DOUBLE PRECISION  RELMIN,D2PI,PENALI,CONMIN,CONMAX,EIGMIN,EIGMAX
      DOUBLE PRECISION  DINFO,DCPU,DPREP,DMOPR,DMAT
      INTEGER           NUMIOB(4),NUMAOB(4),MNDOEL(4),MXDOET(4)
      INTRINSIC         SQRT, REAL
C
      CHARACTER*10      NMTYOB,KNM
      CHARACTER*(*)     KNOMOB
      CHARACTER*160     KNOM
      DATA              RELMIN/-1D28/
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000) KNOMOB
      ELSE
         WRITE(IMPRIM,20000) KNOMOB
      ENDIF
10000 FORMAT(/1X,80('=')/
     %' RESOLUTION des MODES PROPRES de l''OBJET: ',A/1X,80('='))
20000 FORMAT(/1X,80('=')/
     %' EIGENSOLUTIONS (FREQUENCIES and MODE SHAPE VECTORS) of OBJECT: '
     %,A/1X,80('='))
C
C     QUELQUES INITIALISATIONS
      DCPU   = DINFO( 'CPU' )
C     TEMPS POUR L'ELASTICITE
      TEMPS  = 0.0
      DPREP  = 0D0
      DMOPR  = 0D0
      DCPU   = 0D0
      IERR   = 0
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
      MNIP   = 0
      MNELAS = 0
      MNFORC = 0
      MNX    = 0
      MNMUMG = 0
      MNMUKG = 0
      DO 2 I=1,4
         MNDOEL(I) = 0
         MXDOET(I) = 0
 2    CONTINUE
      NBTYEL = 0
      MOAUX  = 0
      MOTAEL = 0
      NDSM   = 1
      NBDLMX = 0
      NTDL   = 0
      MOFLTO = 0
      MOFLPT = 0
      MONDLX = 0
      MNNDLX = 0
      MNVDLX = 0
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
      PENALI = 0D0
C
C     =======================================
C     CALCUL DES MATRICES DE MASSE ET RAIDEUR
C     =======================================
C
C     ADRESSAGE DES ADRESSES DES TABLEAUX NPEF" DE CET OBJET
      CALL TNMCDC( 'ENTIER', 2*MXTYEL, MNNPEF )
      MNTELE = MNNPEF + MXTYEL
C
C     RECHERCHE DES TABLEAUX XYZSOMMET XYZNOEUD XYZPOINT DE L'OBJET
C     -------------------------------------------------------------
      CALL NDPGEL( NTLXOB, NTTOPO, MNTOPO,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %             NBTYEL, MCN(MNTELE), MCN(MNNPEF), IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     LE TYPE DE MAILLAGE ET LES OBJETS INTERNES ET AUX LIMITES
      NDPGST = MCN( MNTOPO + WDPGST )
      NBOBIN = MCN( MNTOPO + WBOBIN )
      NBOBCL = MCN( MNTOPO + WBOBCL )
C
C     RECHERCHE DU MIN ET MAX DES NUMEROS DES OBJETS IMPLIQUES
C     ========================================================
C     DANS L'OBJET ( POINTS, LIGNES, SURFACES, VOLUMES )
C     LES 4 MINIMA D'ABORD, LES 4 MAXIMA ENSUITE
      J = IINFO( 'GRAND' )
      K = - J
      DO 110 I=1,4
         NUMIOB(I) = J
         NUMAOB(I) = K
 110  CONTINUE
C
C     L'ADRESSE MCN DU DEBUT DU TABLEAU NUOBIN
      MNOBIN = MNTOPO + WMTYEL + NBTYEL
      MN     = MNOBIN - 2
      DO 120 I=1,NBOBIN
C        LE TYPE DE L'OBJET INTERNE
         MN   = MN + 2
         NYOB = MCN( MN )
         NUOB = MCN( MN + 1 )
C        LE MINIMUM DES NUMEROS D'OBJETS
         NUMIOB( NYOB ) = MIN( NUMIOB(NYOB), NUOB )
C        LE MAXIMUM DES NUMEROS D'OBJETS
         NUMAOB( NYOB ) = MAX( NUMAOB(NYOB), NUOB )
 120  CONTINUE
C
C     L'ADRESSE MCN DU DEBUT DU TABLEAU NUOBCL
      MNOBCL = MNOBIN + MOTVAR(13) * NBOBIN
      MN     = MNOBCL - 2
      DO 130 I=1,NBOBCL
C        LE TYPE DE L'OBJET AUX LIMITES
         MN   = MN + 2
         NYOB = MCN( MN )
         NUOB = MCN( MN + 1 )
C        LE MINIMUM DES NUMEROS D'OBJETS
         NUMIOB( NYOB ) = MIN( NUMIOB(NYOB), NUOB )
C        LE MAXIMUM DES NUMEROS D'OBJETS
         NUMAOB( NYOB ) = MAX( NUMAOB(NYOB), NUOB )
 130  CONTINUE
C
C     DECLARATION DES TABLEAUX DES DONNEES DES OBJETS
      DO 140 I=1,4
C        NOMBRE DE VARIABLES DES TABLEAUX
         IF( NUMIOB(I) .EQ. J ) THEN
C           OBJET NON REPERTORIE
            NUMAOB(I) = 0
            MNDOEL(I) = 0
         ELSE
            MN = NUMAOB(I) - NUMIOB(I) + 1
            MXDOET(I) = MN * MXDOEL
            CALL TNMCDC( 'ENTIER', MXDOET(I), MNDOEL(I) )
            CALL AZEROI( MXDOET(I), MCN( MNDOEL(I) ) )
         ENDIF
 140  CONTINUE
C     ICI NUMIOB CONTIENT LES 4 NUMEROS MINIMA DES OBJETS
C         NUMAOB          LES 4 NUMEROS MAXIMA DES OBJETS
C     MNDOEL LES 4 ADRESSES MCN DES TABLEAUX DES ADRESSES DES
C     TABLEAUX DECRIVANT L'ELASTICITE DE L'OBJET COMPLET
C
C     OUVERTURE DES TABLEAUX DES DONNEES DES OBJETS
C     *********************************************
C
C     BOUCLE SUR LES OBJETS "INTERNES"
C     ================================
      IERR   = 0
      NERR   = 0
      MN     = MNOBIN - 2
      DO 200 I=1,NBOBIN
C        LE TYPE DE L'OBJET
         MN   = MN + 2
         NYOB = MCN( MN )
         NUOB = MCN( MN + 1 )
C
         IF( NYOB .LE. 2 ) THEN
C           POINT OU LIGNE OBJET INTERNE => SA DIMENSION EST <=2 => ERREUR
            KNM = NMTYOB( NYOB )
            CALL NMOBNU( KNM, NUOB, KNOM )
            NBLGRC(NRERR) = 2
            N = NUDCNB( KNOM )
            KERR(1) = KNM // ' : ' // KNOM(1:N)
            IF( LANGAG .EQ. 0 ) THEN
               KERR(2) = 'NE PEUT ETRE INTERNE'
            ELSE
               KERR(2) = 'CAN''T BE INTERNAL'
            ENDIF
            CALL LEREUR
            NERR = NERR + 1
            GOTO 200
         ENDIF
C
C        OUVERTURE DE L'OBJET
         CALL LXNLOU( NTMN(NYOB), NUOB, NTOB, MNOB )
         IF( NTOB .LE. 0 ) THEN
            KNM = NMTYOB( NYOB )
            CALL NMOBNU( KNM, NUOB, KNOM )
            NBLGRC(NRERR) = 2
            N = NUDCNB( KNOM )
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: ' // KNOM(1:N)
               KERR(2) = KNM // ' INTERNE INCONNU'
            ELSE
               KERR(1) = 'ERROR: ' // KNOM(1:N)
               KERR(2) = KNM // ' UNKNOWN INTERNAL'
            ENDIF
            CALL LEREUR
            NERR = NERR + 1
            GOTO 200
         ENDIF
C
C        L'ADRESSE DU DEBUT DES TABLEAUX DES DONNEES ELASTICITE DE L'OBJET
         MN1 = MNDOEL( NYOB )
         MN1 = MN1 + MXDOEL * ( NUOB - NUMIOB(NYOB) ) - 1
C
C        OUVERTURE DU TABLEAU YOUNG ( ET POISSON )
         CALL LXTSOU( NTOB, 'YOUNG', NT, MCN(MN1+LPYOUN) )
         IF( MCN(MN1+LPYOUN) .LE. 0 ) THEN
            KNM = NMTYOB( NYOB )
            CALL NMOBNU( KNM, NUOB, KNOM )
            NBLGRC(NRERR) = 2
            N = NUDCNB( KNOM )
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR OBJET '// KNOM(1:N)
               KERR(2) = 'SANS DONNEE de YOUNG et POISSON'
            ELSE
               KERR(1) = 'ERROR OBJECT '// KNOM(1:N)
               KERR(2) = 'WITHOUT INPUT of YOUNG and POISSON'
            ENDIF
            CALL LEREUR
            IERR = IERR + 1
            GOTO 9999
         ENDIF
C
C        OUVERTURE DU TABLEAU MASSE
         CALL LXTSOU( NTOB, 'MASSE', NT, MCN(MN1+LPMASS) )
         IF( MCN(MN1+LPMASS) .LE. 0 ) THEN
            KNM = NMTYOB( NYOB )
            CALL NMOBNU( KNM, NUOB, KNOM )
            NBLGRC(NRERR) = 2
            N = NUDCNB( KNOM )
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR OBJET '// KNOM(1:N)
               KERR(2) = 'SANS DONNEE DE MASSE VOLUMIQUE'
            ELSE
               KERR(1) = 'ERROR OBJECT '// KNOM(1:N)
               KERR(2) = 'WITHOUT INPUT of VOLUMINAL MASS'
            ENDIF
            CALL LEREUR
            IERR = IERR + 1
            GOTO 9999
         ENDIF
C
 200  CONTINUE
C
C     BOUCLE SUR LES OBJETS AUX LIMITES DE L'OBJET
C     ============================================
      MN = MNOBCL - 2
      KDFIXA = 0
      DO 210 I=1,NBOBCL
C        LE TYPE DE L'OBJET
         MN   = MN + 2
         NYOB = MCN( MN )
         NUOB = MCN( MN + 1 )
C        OUVERTURE DE L'OBJET
         CALL LXNLOU( NTMN(NYOB), NUOB, NTOB, MNOB )
         IF( NTOB .LE. 0 ) THEN
            KNM = NMTYOB( NYOB )
            CALL NMOBNU( KNM, NUOB, KNOM )
            NBLGRC(NRERR) = 3
            N = NUDCNB( KNOM )
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: AUX LIMITES '
               KERR(2) = KNOM(1:N)
               KERR(3) = KNM // ' INCONNU'
            ELSE
               KERR(1) = 'ERROR: At BOUNDARY '
               KERR(2) = KNOM(1:N)
               KERR(3) = KNM // ' UNKNOWN'
            ENDIF
            CALL LEREUR
            NERR = NERR + 1
            GOTO 210
         ENDIF
C
C        L'ADRESSE DU DEBUT DES TABLEAUX DES DONNEES ELASTICITE DE L'OBJET
         MN1 = MNDOEL( NYOB )
         MN1 = MN1 + MXDOEL * ( NUOB - NUMIOB(NYOB) ) - 1
C
C        OUVERTURE DU TABLEAU FIXATION
         CALL LXTSOU( NTOB, 'FIXATION', NT, MN2 )
         MCN(MN1+LPFIXA) = MN2
         IF( MN2 .GT. 0 ) KDFIXA = 1
C
 210  CONTINUE
C
C     BILAN DES DONNEES
C     =================
      IF( KDFIXA .EQ. 0 ) THEN
         NBLGRC(NRERR) = 2
         N = NUDCNB( KNOMOB )
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='ERREUR: AUCUNE FIXATION SUR L''OBJET '//KNOMOB(1:N)
            KERR(2)='AJOUTER DES FIXATIONS'
         ELSE
            KERR(1)='ERROR: NO FIXATION of OBJECT '// KNOMOB(1:N)
            KERR(2)='ADD FIXATIONS'
         ENDIF
         CALL LEREUR
         NERR = NERR + 1
      ENDIF
      IF( NERR .GT. 0 ) THEN
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
C     LE NOMBRE DE NOEUDS DU MAILLAGE DE L'OBJET
      NBNOEU = MCN( MNXYZN + WNBNOE )
C     LE NOMBRE TOTAL DE DEGRES DE LIBERTE
      NTDL   = NDIM * NBNOEU
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10210) NDIM,NBNOEU,NTDL
      ELSE
         WRITE(IMPRIM,20210) NDIM,NBNOEU,NTDL
      ENDIF
10210 FORMAT(' DIMENSION 2 OU 3 DE L''ESPACE',T39,'=',I6/
     %' NOMBRE DE NOEUDS'                    ,T39,'=',I6/
     %' NOMBRE DE COMPOSANTES DE CHAQUE MODE',T39,'=',I6)
20210 FORMAT(' SPACE DIMENSION (2 or 3)'  ,T37,'=',I6/
     %' NUMBER of NODES'                  ,T37,'=',I6/
     %' NUMBER of COMPONENTS of each MODE',T37,'=',I6)
C
C     LES TABLEAUX POBA NECESSAIRES AUX TABLEAUX ELEMENTAIRES
C     =======================================================
C     OUVERTURE DU FICHIER DIRECT POBA SUPPORT DE TABLEAUX ELEMENTAIRES
      CALL TRUNIT( NFPOBA )
      KNOM = HOMDIR // '/pp/pxyz'
      N    = NUDCNB( KNOM )
      OPEN( UNIT=NFPOBA, ERR=9900, STATUS='OLD',
     %      FILE=KNOM(1:N), ACCESS='DIRECT', FORM='UNFORMATTED',
     %      RECL=MOPAGE*NBCHMO )
C     LECTURE DE LA PREMIERE PAGE
      READ (UNIT=NFPOBA,REC=1) NELFI
C
      CALL TNMCDC( 'ENTIER', NBTYEL*MXPOBA, MNTPOB )
      MOAUX  = 1
      NOAXIS = 0
      DO 250 NOTYEL=1,NBTYEL
C        L'ADRESSE DU TABLEAU ELEMENTS
         MNELE  = MCN( MNNPEF - 1 + NOTYEL )
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
         IF (NUTYEL.LE.4) THEN
            NOAXIS=1
         ENDIF
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELTYCA( NUTYEL )
C        LES CARACTERISTIQUES DES TABLEAUX DE L'ELEMENT FINI
         CALL EETAEL(NUTYEL,NDSM,NBDL,NCODEM,NTPOBA,NOMTAB,MOTAUX)
C
C        NBDL      : NOMBRE DE DL DE L ELEMENT
C        NTPOBA    : NOMBRE DE TABLEAUX A LIRE SUR LE FICHIER POBA
C        NOMTAB(J) : NOM DU J-EME TABLEAU DE POBA A LIRE 1<=J<=NTPOBA
C        MOTAUX    : NOMBRE DE MOTS AUXILIAIRES NECESSAIRES A L ELEMENT
C
C        DECLARATION ET LECTURE DES TABLEAUX DE POBA
         IF( NTPOBA .GT. 0 ) THEN
            DO 240 J=1,NTPOBA
               CALL FINDEL(NOMTAB(J),MCN(MNTPOB-1+(NOTYEL-1)*MXPOBA+J),
     %                     NFPOBA,MOPAGE)
 240        CONTINUE
         ENDIF
C
C        BILAN DES TABLEAUX AUXILIAIRES ET ELEMENTAIRES
         IF( MOAUX  .LT. MOTAUX ) MOAUX  = MOTAUX
         IF( NBDLMX .LT. NBDL   ) NBDLMX = NBDL
 250  CONTINUE
C
C     FERMETURE DU FICHIER POBA
      CLOSE( UNIT=NFPOBA )
C
C     ADRESSAGE DES TABLEAUX AUXILIAIRES ET ELEMENTAIRES
C     ===================================================
      CALL TNMCDC( 'REEL2', MOAUX, MNTAUX )
C     LA MATRICE ELEMENTAIRE ET LES VECTEURS NECESSAIRES
      MOTAEL = NBDLMX * (NBDLMX+1) / 2 + NBDLMX * NDSM
      CALL TNMCDC( 'REEL2',  MOTAEL, MNTAEL )
      CALL TNMCDC( 'ENTIER', NBDLMX, MNNODL )
      CALL TNMCDC( 'ENTIER', NBDLMX, MNIP )
      CALL TNMCDC( 'REEL2',  27,     MNELAS )
C     6*NDSM A CAUSE DES CONTRAINTES INITIALES
      MOFORC = 6 * NDSM
      CALL TNMCDC( 'REEL2',  MOFORC, MNFORC )
C     LE TABLEAU DU NUMERO DES NOEUDS D'UN ELEMENT FINI
      CALL TNMCDC( 'ENTIER', NBDLMX*NDIM, MNX )
C
C     PREPARATION DE LA RESOLUTION STOCKAGE PROFIL DE K ET M
C     ============================ =========================
C     CALCUL DU PROFIL DES MATRICES
C     (LES MATRICES PROFIL SONT ICI SYMETRIQUES)
      NCODSK = 1
      NCODSM = 1
      CALL TNMCDC( 'ENTIER', 1+NTDL, MNMUKG )
      CALL PRPRMC( MNTOPO, MCN(MNNPEF), MNXYZN, NDIM, NCODSK,
     %             MCN(MNMUKG), IERR )
      IF( IERR .GT. 0 ) GOTO 9999
C
C     COPIE DU TABLEAU POINTEUR (ENSUITE PEUT ETRE MODIFIE PAR LES CL)
      CALL TNMCDC( 'ENTIER', 1+NTDL, MNMUMG )
      CALL TRTATA( MCN(MNMUKG), MCN(MNMUMG), 1+NTDL )
C
C     DECLARATION DES 2 MATRICES PROFIL SYMETRIQUE
C                    +1 MATRICE DE PROTECTION DANS CALVVP
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
     %               ' MOTS NECESSAIRES pour MG + KG + RG'
             KERR(3) = 'REDUIRE le MAILLAGE'
          ELSE
             KERR(1) ='ERROR: NOT ENOUGH MEMORY'
             KERR(2) = KERR(MXLGER-1)(1:25) //
     %               ' NECESSARY WORDS to store MG + KG + RG'
             KERR(3) ='REDUCE the MESH'
          ENDIF
          CALL LEREUR
          IF( INTERA .LE. 0 ) CALL ARRET( 100 )
          GOTO 9999
      ENDIF
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10290) NBRDKG, NBRDKG/NTDL
      ELSE
         WRITE(IMPRIM,20290) NBRDKG, NBRDKG/NTDL
      ENDIF
10290 FORMAT('3 MATRICES PROFIL CHACUNE DE',I15,
     %' REELS DOUBLE PRECISION'/
     %'1/2 LARGEUR DE BANDE MOYENNE =',I9)
20290 FORMAT('3 SKYLINE MATRICES EACH of',I15,' DOUBLE PRECISION REALS'/
     %'HALF WIDTH AVERAGE=',I9)
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
C     ADRESSAGE DE KG & MG SYMETRIQUES NON DIAGONALES PROFIL
      NCODSK = 1
      NCODSM = 1
C
C     CALCUL DES MATRICES DE MASSE ET RAIDEUR SUPPOSES INDEPENDANTES
C     DU TEMPS ET DES DEPLACEMENTS
C     --------------------------------------------------------------
      CALL ELAMKB( 1,      1,      0,
     %             PENALI, D2PI,   NDIM,
     %             NBTYEL, MNNPEF, NDPGST,
     %             MNTPOB, MXPOBA, MNTAUX,
     %             MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     %             MNELAS, MNFORC, MNTAEL, MNX,    MNIP,   MNNODL,
     %             0,      0,      1,      MNMUKG, 0,
     %             NBRDMG, MG,     NBRDKG, KG,   NTDL, 0,
     %             NCODSM, NCODSK, NPIMAX )
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
      CALL ELDLFX( NTDL,   NDIM,   NBTYEL, MNNPEF, NDPGST,
     %             MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %             NBRDLX, MONDLX, MNNDLX, MNVDLX,
     %             IERR )
      IF( IERR .NE. 0 ) GOTO 9999
C
C     TEMPS CALCUL DE K ET M
      DPREP = DINFO( 'DELTA CPU' )
C
C     CALCUL DES VALEURS EIGV ET VECTEURS PROPRES VP TELS QUE
C     ( K + EIGV M ) V = 0
C     =======================================================
      METVVP = 0
      CALL CALVVP( 'ELASTICITE', NTLXOB, NTDL,
     %             METVVP, EIGMIN,EIGMAX, NBROOT,
     %             NBRDLX, MNNDLX, MNVDLX,
     %             NCODSK, MNMUKG, KG,
     %             NCODSM, MNMUMG, MG,
     %             NTVVPR, MNVVPR, IERR )
      IF( IERR .NE. 0 ) GOTO 9998
C
C     PASSAGE DES VALEURS PROPRES OMEGA**2 AUX FREQUENCES EN Hz
C     PAR LA FORMULE:  FREQUENCE = SQRT( OMEGA**2 ) / 2 Pi
C     =========================================================
C     LE NOMBRE DE VALEURS PROPRES ET DE DL
      NBROOT = MCN(MNVVPR+WBVECT)
      MN     = MNVVPR + WECTEU + MOREE2*NTDL*NBROOT
      DO 30 I=1,NBROOT
C        LA VALEUR PROPRE I = OMEGA AU CARRE >=0
         S = RMCN(MN)
         IF( S .LT. 0 ) THEN
            NBLGRC(NRERR) = 3
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'VALEUR PROPRE      NEGATIVE'
               WRITE(KERR(1)(15:19),'(I5)') I
               KERR(2) = 'VALEUR PROPRE='
               WRITE(KERR(2)(15:29),'(G15.6)') S
               KERR(3) = 'VALEUR IMPOSEE A ZERO'
            ELSE
               KERR(1) = 'EIGENVALUE       NEGATIVE'
               WRITE(KERR(1)(13:17),'(I5)') I
               KERR(2) = 'EIGENVALUE='
               WRITE(KERR(2)(12:26),'(G15.6)') S
               KERR(3) = 'VALUE IMPOSED NULL'
            ENDIF
            CALL LEREUR
            S = 0
         ENDIF
         RMCN(MN) = REAL( SQRT( S ) / D2PI )
         MN = MN + 1
 30   CONTINUE
C
C
C     FERMETURE DE LA MATRICE POUR REDONNER DE LA PLACE EN MC
C     =======================================================
 9998 IF(MNMUKG .GT. 0) CALL TNMCDS( 'ENTIER', NTDL+1, MNMUKG )
      IF(MNMUMG .GT. 0) CALL TNMCDS( 'ENTIER', NTDL+1, MNMUMG )
C     FREE MEMORY USED BY ARRAY VG
      IF( IERKGALLOC .EQ. 0 ) THEN
         DEALLOCATE( KG )
         IERKGALLOC = 1
      ENDIF
      IF( IERMGALLOC .EQ. 0 ) THEN
         DEALLOCATE( MG )
         IERMGALLOC = 1
      ENDIF
      IF(IERR   .NE. 0) GOTO 9999
C
C     COUT CALCUL DES MODES PROPRES
C     =============================
      DMOPR = DINFO( 'DELTA CPU' )
C
C **************************************************************************
C --------------------------------------------------------------------------
C **************************************************************************
C
C     CALCUL DES CONTRAINTES EN CHAQUE POINT D INTEGRATION DE
C     CHAQUE ELEMENT FINI DE CHAQUE TYPE D'EF DE L'OBJET
C     =======================================================
      CALL ELASTR( NTLXOB, MNTOPO, NOAXIS, NDIM,   MOREE2,
     %             NPIMAX, NBROOT, NBROOT, NTDL,
     %             NBTYEL, MNNPEF, NDPGST, MNTPOB,
     %             MNTAUX, MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     %             MNELAS, MOTAEL, MNTAEL, MNX,
     %             MNVVPR, 1,      0,
     %             CONMIN, CONMAX  )
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
      IF(MNMUMG .GT. 0) CALL TNMCDS( 'ENTIER', NTDL+1, MNMUMG )
      IF(MNMUKG .GT. 0) CALL TNMCDS( 'ENTIER', NTDL+1, MNMUKG )
C
      IF( MNNPEF .GT. 0 ) CALL TNMCDS( 'ENTIER',2*MXTYEL, MNNPEF )
      DO 11000 I=1,4
         IF( MNDOEL(I) .GT. 0 .AND. MXDOET(I) .GT. 0 ) THEN
            CALL TNMCDS( 'ENTIER', MXDOET(I), MNDOEL(I) )
         ENDIF
11000 CONTINUE
      IF( MNTPOB .GT. 0 ) CALL TNMCDS( 'ENTIER', NBTYEL*MXPOBA,MNTPOB )
      IF( MNTAUX .GT. 0 ) CALL TNMCDS( 'REEL2' , MOAUX , MNTAUX )
      IF( MNTAEL .GT. 0 ) CALL TNMCDS( 'REEL2' , MOTAEL, MNTAEL )
      IF( MNNODL .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLMX, MNNODL )
      IF( MNIP   .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLMX, MNIP )
      IF( MNELAS .GT. 0 ) CALL TNMCDS( 'REEL2' , 27    , MNELAS )
      IF( MNFORC .GT. 0 ) CALL TNMCDS( 'REEL2' , MOFORC, MNFORC )
      IF( MNX    .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLMX*NDIM, MNX )
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
C     COUT CALCUL DES CONTRAINTES
      DCPU = DINFO( 'DELTA CPU' )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,12001) DPREP,DMOPR,DCPU,DPREP+DMOPR+DCPU
      ELSE
         WRITE(IMPRIM,22001) DPREP,DMOPR,DCPU,DPREP+DMOPR+DCPU
      ENDIF
12001 FORMAT(/
     %' TEMPS CALCUL DE LA PREPARATION =',F12.2,' SECONDES CPU'/
     %' TEMPS CALCUL DES MODES PROPRES =',F12.2,' SECONDES CPU'/
     %' TEMPS CALCUL DES CONTRAINTES   =',F12.2,' SECONDES CPU'/
     %' TEMPS CALCUL TOTAL             =',F12.2,' SECONDES CPU')
22001 FORMAT(/
     %' PREPARATION   TIME =',F12.2,' CPU SECONDS'/
     %' EIGENVALUES   TIME =',F12.2,' CPU SECONDS'/
     %' MAIN STRESSES TIME =',F12.2,' CPU SECONDS'/
     %' TOTAL         TIME =',F12.2,' CPU SECONDS')
      RETURN
C
C     ERREUR A L'OUVERTURE DU FICHIER POBA
 9900 NBLGRC(NRERR) = 2
      N = NUDCNB( KNOM )
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'IMPOSSIBLE D''OUVRIR LE FICHIER POBA'
         KERR(2) = KNOM(1:N)
      ELSE
         KERR(1) = 'IMPOSSIBLE to OPEN the FILE POBA'
         KERR(2) = KNOM(1:N)
      ENDIF
      CALL LEREUR

      RETURN
      END
