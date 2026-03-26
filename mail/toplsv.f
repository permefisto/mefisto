      SUBROUTINE TOPLSV( NOMOBJ, NTPLSV, NMPLSV, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER L'OBJET NOMOBJ + un PLSV pour introduire des DONNEES
C------
C ENTREES:
C --------
C NOMOBJ : NOM DE L'OBJET AVEC INTERPOLATION A TRACER
C NTPLSV : 1:POINT, 2:LIGNE, 3:SURFACE, 4:VOLUME
C NMPLSV : NOM DU PLSV    A TRACER AVEC L'OBJET
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR, >0 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 2002
C.......................................................................
      IMPLICIT          INTEGER (W)
      PARAMETER        (MXTYEL=7, MXPILE=128)
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      COMMON / MSSFTA / MSSF(28),NTADAM
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (RMCN(1),MCN(1))
      INTEGER           NTPLSV
      CHARACTER*(*)     NOMOBJ, NMPLSV
      CHARACTER*24      KNOMOB
      CHARACTER*10      NMTYOB
C
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/xyzext.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___aretefr.inc"
      include"./incl/ponoel.inc"
C
C     L'OBJET EXISTE-T-IL ?
C     =====================
      IERR = 0
      CALL LXLXOU( NTOBJE, NOMOBJ, NTLXOB, MNLXOB )
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR: OBJET INCONNU ' // NOMOBJ
         CALL LEREUR
         IERR   = 1
         KNOMOB = NOMOBJ
         GOTO 9999
      ENDIF
C
C     L'OBJET EST IL AVEC OU SANS INTERPOLATION ?
C     ===========================================
C     LE TABLEAU TOPOLOGIE DE L'OBJET
      CALL LXTSOU( NTLXOB, 'TOPOLOGIE', NTTOPO, MNTOPO )
      IF( NTTOPO .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR: OBJET SANS INTERPOLATION ' // NOMOBJ
         CALL LEREUR
         IERR   = 2
         RETURN
      ENDIF
C
C     LE PLSV EXISTE-T-IL ?
C     =====================
      IF( NTPLSV .LE. 0 .OR. NTPLSV .GT. 4 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR: TYPE <0 ou >4 INCORRECT DE PLSV'
         CALL LEREUR
         IERR   = 3
         RETURN
      ENDIF
      CALL LXLXOU( NTMN(NTPLSV), NMPLSV, NTLXPL, MNLXPL )
      IF( NTLXPL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR: PLSV INCONNU ' // NMPLSV
         CALL LEREUR
         IERR   = 4
         RETURN
      ENDIF
C     NUMERO DU PLSV DANS SON LEXIQUE
      CALL NUOBNM( NMTYOB( NTPLSV ), NMPLSV, NUPLSV )
C
C     OBJET AVEC INTERPOLATION
C     ========================
C     RECHERCHE DES TABLEAUX SOMMETS NOEUDS POINTS ASSOCIES A L'OBJET
C     ADRESSAGE DES ADRESSES DES TABLEAUX ELEMENTS DE CET OBJET
      MNELEM = 0
      CALL TNMCDC( 'ENTIER', 2*MXTYEL, MNELEM )
      MNTELE = MNELEM + MXTYEL
      CALL NDPGEL( NTLXOB, NTTOPO, MNTOPO,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %             NBTYEL, MCN(MNTELE), MCN(MNELEM), IERR )
C     NTTOPO : NUMERO      DU TMS 'TOPOLOGIE' DE L'OBJET
C     MNTOPO : ADRESSE MCN DU TMS 'TOPOLOGIE' DE L'OBJET
C     NTXYZP : NUMERO      DU TMS 'XYZPOINT'    DE L'OBJET
C     MNXYZP : ADRESSE MCN DU TMS 'XYZPOINT'    DE L'OBJET
C     NTXYZN : NUMERO      DU TMS 'XYZNOEUD'    DE L'OBJET
C     MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD'    DE L'OBJET
C     NBTYEL : NOMBRE DE TYPES D'ELEMENTS FINIS DU MAILLAGE
C     NTELEM : NUMERO      DU TMS DES NBTYEL TYPES D'ELEMENTS
C     MNELEM : ADRESSE MCN DU TMS DES NBTYEL TYPES D'ELEMENTS
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR: INTERPOLATION INCORRECTE OBJET ' // NOMOBJ
         CALL LEREUR
         IERR   = 5
         RETURN
      ENDIF
C
C     NDPGST : CODE TRAITEMENT DES XYZ DES SOMMETS POINTS NOEUDS DU MAILLAGE
C              0 : NOEUDS=POINTS=SOMMETS
C              1 : NOEUDS=POINTS#SOMMETS
C              2 : NOEUDS#POINTS=SOMMETS
C              3 : NOEUDS#POINTS#SOMMETS
      NDPGST = MCN( MNTOPO + WDPGST )
C
C     NDIM DIMENSION EFFECTIVE DE L'ESPACE DES COORDONNEES
      NBCOOR = MCN(MNXYZP+WBCOOP)
      NBPOMA = MCN(MNXYZP+WNBPOI)
      IF( NBCOOR .EQ. 6 ) THEN
         NDIM = 3
      ELSE
         CALL DIMCOO( NBPOMA, MCN(MNXYZP+WYZPOI), NDIM )
      ENDIF
      NDIMLI = NDIM
      IF( NDIM .LE. 2 ) RETURN
C
C     OBJET 3D AVEC INTERPOLATION
C     ===========================
C
C     LES ARETES FRONTALIERES DES VOLUMES DE L'OBJET
C     ----------------------------------------------
C     CREATION OU REDECOUVERTE DU TMS OBJET>>>FACE
      CALL HACHOB( NOMOBJ, 4, NTFAOB, MNFAOB, IERR )
      IF( IERR .GT. 0 .OR. NTFAOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'ERREUR TOPLSV: OBJET ' // NOMOBJ
         KERR(2) = 'CALCUL IMPOSSIBLE DE SES FACES'
         CALL LEREUR
         RETURN
      ENDIF
C
C     CREATION DU HACHAGE DES ARETES DES FACES FRONTALIERES DE L'OBJET
      CALL HACHAF( NOMOBJ, 0, NTFAOB, MNFAOB, NTAFOB, MNAFOB, I )
C
C     LE NOMBRE D'ENTIERS PAR ARETE FRONTALIERE
      MOARFR = MCN( MNAFOB + WOARFR )
C     LA MAJORATION DU NOMBRE DES ARETES FRONTALIERES
      MXARFR = MCN( MNAFOB + WXARFR )
C     LE NUMERO DANS LAREFR DE LA PREMIERE ARETE FRONTALIERE
      L1ARFR = MCN( MNAFOB + W1ARFR )
C     LE NOMBRE D'ARETES FRONTALIERES DANS LE CHAINAGE
      NBARFR = MCN( MNAFOB + WBARFR )
C
C     BOUCLE SUR LES DIFFERENTS TYPES D'ELEMENTS FINIS DU MAILLAGE
C     POUR CALCULER LE NOMBRE DE FACES, ARETES, SOMMETS APPARTENANT A
C     DES SURFACES, LIGNES ET POINTS NOMMES PAR L'UTILISATEUR
C     ===============================================================
      NBFA = 0
      NBAR = 0
      NBSO = 0
C
      DO 110 I = 0, NBTYEL-1
C
C        L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF DE CE TYPE D'EF
         MNELE = MCN( MNELEM + I )
C
C        LES SOMMETS IDENTIQUES A DES POINTS UTILISATEUR
         NBSO = NBSO + MCN( MNELE + WBPSEL )
C
C        LES ARETES APPARTENANT A UNE LIGNE
         NBAR = NBAR + MCN( MNELE + WBLAEL )
C
C        LES FACES APPARTENANT A UNE SURFACE
         NBFA = NBFA + MCN( MNELE + WBSFEL )
C
 110  CONTINUE
C
C     CREATION DU TABLEAU DES NUMEROS DANS XYZPOINT DES SOMMETS
C     DES FACES ARETES SOMMETS APPARTENANT A UNE SURFACE LIGNE POINT
C     --------------------------------------------------------------
      MNFASU = 0
      MNARLI = 0
      MNSOPO = 0
      IF( NBFA .GT. 0 ) THEN
         MNFASU = 0
         CALL TNMCDC( 'ENTIER', 5*NBFA, MNFASU )
      ENDIF
      MNF = MNFASU - 1
C
      IF( NBAR .GT. 0 ) THEN
         MNARLI = 0
         CALL TNMCDC( 'ENTIER', 3*NBAR, MNARLI )
      ENDIF
      MNA = MNARLI - 1
C
      IF( NBSO .GT. 0 ) THEN
         MNSOPO = 0
         CALL TNMCDC( 'ENTIER', 2*NBSO, MNSOPO )
      ENDIF
      MNS = MNSOPO - 1
C
C     INITIALISATION EVENTUELLE DES 3 TABLEAUX
      NBSP = 0
      NBAL = 0
      NBFS = 0
      DO 200 I = 0, NBTYEL-1
C
C        L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF DE CE TYPE D'EF
         MNELE = MCN( MNELEM + I )
C
C        LE NUMERO DU TYPE DES ELEMENTS FINIS
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        LE NOMBRE DE TELS ELEMENTS
         NBELEM = MCN( MNELE + WBELEM )
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELTYCA( NUTYEL )
C
C        L'ADRESSE MCN DU TABLEAU 'NPEF' POUR CE TYPE D'EF
         NBNDEL = MCN( MNELE + WBNDEL )
         NBPGEL = MCN( MNELE + WBPGEL )
C        L'ADRESSE DE NUNDEL PUIS NUPGEL
         MNPGEL = MNELE + WUNDEL
         IF( NDPGST .LT. 2 ) THEN
C           NOEUDS=POINTS
            NBPGEL = NBNDEL
         ELSE
C           NOEUDS#POINTS
            NBPGEL = MCN( MNELE + WBPGEL )
            MNPGEL = MNPGEL + NBELEM * NBNDEL
         ENDIF
C
C        ATTENTION : RECUL D'UN MOT POUR EVITER LA SOUSTRACTION DANS LES BOUCLES
         MNPGEL = MNPGEL - 1
C
C        LES SOMMETS = POINTS UTILISATEUR
         NBPSEL = MCN( MNELE + WBPSEL )
         MOPSEL = MCN( MNELE + WOPSEL )
C        ADRESSE JUSTE AVANT LE TABLEAU NLPSEL
         MNPSEL = MNPGEL + NBELEM * NBPGEL + NBPSEL
         IF( NTPLSV .EQ. 1 ) THEN
            DO 120 J = 1,NBPSEL
C              LE NUMERO DU POINT
               NUM = MCN( MNPSEL - NBPSEL + J )
               IF( NUM .EQ. NUPLSV ) THEN
C                 LE NUMERO LOCAL DU SOMMET DANS L'EF
                  N   = MCN( MNPSEL + J )
C                 LE NUMERO DE L'EF DE CE SOMMET
                  NEF = MCN( MNPSEL + MOPSEL + J )
C                 LE NUMERO XYZPOI DE CE SOMMET
                  NS  = MCN( MNPGEL + NEF + NBELEM * ( N - 1 ) )
C                 CE NUMERO EST STOCKE DANS LE TABLEAU DES SOMMETS=POINTS
                  MCN( MNS + 1 ) = NS
                  MCN( MNS + 2 ) = NUM
                  MNS  = MNS + 2
                  NBSP = NBSP + 1
               ENDIF
 120        CONTINUE
         ENDIF
C
C        LES ARETES APPARTENANT A UNE LIGNE
         NBLAEL = MCN( MNELE + WBLAEL )
         MOLAEL = MCN( MNELE + WOLAEL )
C        ADRESSE JUSTE AVANT LE TABLEAU NLLAEL
         MNLAEL = MNPSEL + MOPSEL + MOPSEL + NBLAEL
         IF( NTPLSV .EQ. 2 ) THEN
         DO 130 J = 1,NBLAEL
C           LE NUMERO DE LA LIGNE
            NUM = MCN( MNLAEL - NBLAEL + J )
            IF( NUM .EQ. NUPLSV ) THEN
C              LE NUMERO LOCAL DE L'ARETE DANS L'EF
               N   = MCN( MNLAEL + J )
C              LE NUMERO DE L'EF DE CETTE ARETE
               NEF = MCN( MNLAEL + MOLAEL + J )
C              LE NUMERO XYZPOI DE SES 2 SOMMETS
               DO 125 K=1,2
C                 LE NUMERO LOCAL DU SOMMET DANS L'EF
                  NS = NOSOAR( K, N )
                  NS = MCN( MNPGEL + NEF + NBELEM * ( NS - 1 ) )
C                 CE NUMERO EST STOCKE DANS LE TABLEAU DES SOMMETS=POINTS
                  MCN( MNA + K ) = NS
 125           CONTINUE
               MCN( MNA + 3 ) = NUM
               MNA  = MNA + 3
               NBAL = NBAL + 1
            ENDIF
 130     CONTINUE
         ENDIF
C
C        LES FACES APPARTENANT A UNE SURFACE
         NBSFEL = MCN( MNELE + WBSFEL )
         MOSFEL = MCN( MNELE + WOSFEL )
C        ADRESSE JUSTE AVANT LE TABLEAU NLSFEL
         MNSFEL = MNLAEL + MOLAEL + MOLAEL + NBSFEL
         IF( NTPLSV .EQ. 3 ) THEN
            DO 140 J = 1,NBSFEL
C              LE NUMERO DE LA SURFACE
               NUM = MCN( MNSFEL - NBSFEL + J )
               IF( NUM .EQ. NUPLSV ) THEN
C                 LE NUMERO LOCAL DE LA FACE DANS L'EF
                  N   = MCN( MNSFEL + J )
C                 LE NUMERO DE L'EF DE CETTE FACE
                  NEF = MCN( MNSFEL + MOSFEL + J )
C                 LE NOMBRE DE SOMMETS DE CETTE FACE
                  NBSTFA = NBSOFA(N)
C                 TEMOIN DE TRIANGLE, ECRASE SI CETTE FACE A 4 SOMMETS
                  MCN( MNF + 4 ) = 0
C                 LE NUMERO XYZPOI DE SES NBSTFA SOMMETS
                  DO 135 K=1,NBSTFA
C                    LE NUMERO LOCAL DU SOMMET DANS L'EF
                     NS = NOSOFA( K, N )
C                    LE NUMERO XYZPOI DE CE SOMMET
                     NS = MCN( MNPGEL + NEF + NBELEM * ( NS - 1 ) )
C                    CE NUMERO EST STOCKE DANS LE TABLEAU DES SOMMETS=POINTS
                     MCN( MNF + K ) = NS
 135              CONTINUE
                  MCN( MNF + 5 ) = NUM
                  MNF  = MNF + 5
                  NBFS = NBFS + 1
               ENDIF
 140        CONTINUE
         ENDIF
 200  CONTINUE
C
C     LE TABLEAU DES Z AXONOMETRIQUES DES BARYCENTRES DES FACES ARETES ET SOMMET
      MNBARY = 0
      MNNUBA = 0
      NBBARY = NBSP + NBAL + NBFS + NBARFR
      CALL TNMCDC( 'REEL',   NBBARY, MNBARY )
      CALL TNMCDC( 'ENTIER', NBBARY, MNNUBA )
C
C     FORMATION DES TABLEAUX Z-AXONOMETRIQUES ET NUMERO DES BARYCENTRES AVANT TR
C     --------------------------------------------------------------------------
      CALL TBAFAS( NBCOOR, NBPOMA, RMCN(MNXYZP+WYZPOI),
     %             NBSP,   MCN(MNSOPO),
     %             NBAL,   MCN(MNARLI),
     %             NBFS,   MCN(MNFASU),
     %             MOARFR, MXARFR, L1ARFR, MCN(MNAFOB+WAREFR),
     %             NBBARY, MCN(MNNUBA), RMCN(MNBARY),
     %             NUMXPO, NUMXLI, NUMXSU )
C
C     LE TRI PAR TAS DE CETTE COTE AXONOMETRIQUE
      CALL TRITRP( NBBARY, RMCN(MNBARY), MCN(MNNUBA) )
C
C     TRACE EFFECTIF DES ARETES FRONTALIERES SURFACES LIGNES ET POINTS
C     ----------------------------------------------------------------
      MNPLS = 0
      MXPLS = NUMXPO + NUMXLI + NUMXSU
      IF( MXPLS .GT. 0 ) THEN
         CALL TNMCDC( 'ENTIER', MXPLS, MNPLS )
         CALL AZEROI( MXPLS, MCN(MNPLS) )
      ELSE
C        ADRESSE BIDON NON UTILISEE ENSUITE
         MNPLS = MNNUBA
      ENDIF
      CALL T3AFAS( NBCOOR, RMCN(MNXYZP+WYZPOI),
     %             NBSP,   MCN(MNSOPO),
     %             NBAL,   MCN(MNARLI),
     %             NBFS,   MCN(MNFASU),
     %             MOARFR, MXARFR, MCN(MNAFOB+WAREFR),
     %             NBBARY, MCN(MNNUBA),
     %             NUMXPO, NUMXLI, NUMXSU, MCN(MNPLS) )
      CALL TRFINS( NMPLSV )
C
C     DESTRUCTION DES TABLEAUX DEVENUS INUTILES
      IF( MXPLS .GT. 0 ) CALL TNMCDS( 'ENTIER', MXPLS,  MNPLS  )
      IF( NBSO  .GT. 0 ) CALL TNMCDS( 'ENTIER', 2*NBSO, MNSOPO )
      IF( NBAR  .GT. 0 ) CALL TNMCDS( 'ENTIER', 3*NBAR, MNARLI )
      IF( NBFA  .GT. 0 ) CALL TNMCDS( 'ENTIER', 5*NBFA, MNFASU )
      IF( MNBARY .GT. 0 )CALL TNMCDS( 'REEL',   NBBARY, MNBARY )
      IF( MNNUBA .GT. 0 )CALL TNMCDS( 'ENTIER', NBBARY, MNNUBA )
      CALL TNMCDS( 'ENTIER', 2*MXTYEL, MNELEM )
      RETURN
C
C     OBJET INCONNU => ERREUR
C     =======================
 9999 NBLGRC(NRERR) = 2
      KERR(1) = 'OBJET INCONNU:' // KNOMOB
      KERR(2) = 'A choisir parmi :'
      CALL LEREUR
      WRITE(IMPRIM,19999) 'OBJET',KNOMOB
19999 FORMAT(1X,A,1X,A,' NON RETROUVE PARMI')
C     OUVERTURE DES OBJETS
      CALL LXLXOU( NTADAM, 'OBJET', NTOBJE, MNOBJE )
      CALL LXIM0( MNOBJE )
      RETURN
      END
