      SUBROUTINE THENORM( KNOMOB, EXPO, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA NORME L1 DES SOLUTIONS U
C -----
C          normp = Som       Som    Omegal  (û**EXPO)(bl) Delta(bl)
C                  e dans E l=1...,L
C          Vol   = Som       Som    Omegal Delta(bl)
C                  e dans E l=1...,L
C
C          POUR DES ELEMENTS FINIS LAGRANGE DE DEGRE 1 OU 2 EN 2D OU 3D
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET A CALCULER
C EXPO   : PUISSANCE DE LA SOLUTION DE NORME L1 A CALCULER
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET  LJLL UPMC & SAINT PIERRE DU PERRAY   MAI 2009
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
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
      EXTERNAL          ETTAEL
      DOUBLE PRECISION  RELMIN
      INTEGER           NUMIOB(4),NUMAOB(4),MNDOEL(4),MXDOEL(4)
C
      DOUBLE PRECISION  EXPO
      CHARACTER*(*)     KNOMOB
      DATA              RELMIN/-1D28/
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000) EXPO, KNOMOB
      ELSE
         WRITE(IMPRIM,20000) EXPO, KNOMOB
      ENDIF
10000 FORMAT(/100('=')/
     %'Calcul de la NORME L1 des VECTEURS SOLUTIONS **',D15.8,
     %' de l''OBJET',A/100('='))
20000 FORMAT(/100('=')/
     %'Computation of the L1-NORM of SOLUTION VECTORS **',D15.8,
     %' of the OBJECT: ',A/100('='))
CC
C     QUELQUES INITIALISATIONS
      TESTNL = 0
C     TEMPS POUR LA THERMIQUE
      TEMPS  = 0.0
      IERR   = 0
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
      MNX    = 0
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
      NDSM   = 1
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
               GOTO 9999
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
         GOTO 9999
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
C     INITIALISATIONS DE VARIABLES ET AFFICHAGES
C     ==========================================
C     LE NOMBRE DE POINTS DU MAILLAGE DE L'OBJET
      NBPOI = MCN( MNXYZP + WNBPOI )
C
C     NOMBRE (3 ou 6) DES COORDONNEES DES NOEUDS
      NBCOOR = MCN( MNXYZN + WBCOON )
C
C     NDIM LA DIMENSION 2 OU 3 DE L'ESPACE DES COORDONNEES
      CALL DIMCOO( NBPOI, MCN(MNXYZP+WYZPOI), NDIM )
C
C     LE NOMBRE DE NOEUDS DU MAILLAGE DE L'OBJET=NOMBRE DE DEGRES DE LIBERTE
      NTDL = MCN( MNXYZN + WNBNOE )
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10210) NDIM,NDSM,NTDL
      ELSE
         WRITE(IMPRIM,20210) NDIM,NDSM,NTDL
      ENDIF
10210 FORMAT(/' DIMENSION 2 OU 3 OU 6 DE L''ESPACE',T42,'=',I6/
     %        ' NOMBRE DE VECTEURS SOLUTIONS',T42,'=',I6/
     %        ' NOMBRE DE COMPOSANTES DE CHAQUE VECTEUR',   T42,'=',I6/)
20210 FORMAT(/' SPACE DIMENSION (2 or 3 or 6)',T40,'=',I6/
     %        ' NUMBER of SOLUTION VECTORS',   T40,'=',I6/
     %        ' NUMBER of COMPONENTS of each VECTOR',T40,'=',I6/)
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
C     LE TABLEAU DES 3 COORDONNEES DES NBDLMX NOEUDS D'UN ELEMENT FINI
      CALL TNMCDC( 'REEL', NBDLMX*NBCOOR, MNX )
C
C     CALCUL DE LA NORME L1 DES NDSM TEMPERATURES**EXPO
C     =================================================
      CALL THNORM( EXPO,   NDIM,   NTDL,   NDSM,   MCN(MNTEMP+WECTEU),
     %             NBCOOR, MNX,    MNNODL, NBTYEL, MNNPEF, NDPGST,
     %             MNTPOB, MXPOBA, MNTAUX, MNXYZP,
     %             IERR  )
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
      IF( MNX    .GT. 0 ) CALL TNMCDS( 'REEL',   NBDLMX*NBCOOR, MNX  )
C
      RETURN
      END
