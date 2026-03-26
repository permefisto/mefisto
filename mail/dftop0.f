      SUBROUTINE DFTOP0( KNOMOB, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DEFINIR LES RECOUVREMENTS ENTRE POINTS LIGNES SURFACES VOLUMES
C ----- 6-CUBES  D'UN OBJET  EN TERMES DE SOMMETS SEULEMENT
C       PRENDRE EN COMPTE L'INTERPOLATION CHOISIE EN CREEANT
C       LES NOEUDS ET LES POINTS A PARTIR DES SOMMETS DES EF
C       RENUMEROTER EVENTUELLEMENT LES NOEUDS DE L'OBJET

C ENTREE:
C -------
C KNOMOB: NOM DE L'OBJET DE TOPOLOGIE A CONSTRUIRE

C SORTIE:
C -------
C IERR  :=0 PAS D'ERREUR
C        >0 SI ERREUR RENCONTREE(=3 =>LE LEXIQUE DE L'OBJET EST DETRUIT)
C        -1 SI ABANDON DEMANDE
C         2 SI TAILLE DE TABLEAU TROP GRANDE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE      UPMC PARIS AVRIL 1989
C MODIF  : ALAIN PERRONNET  ANALYSE NUMERIQUE LJLL UPMC PARIS MARS  2004
C MODIF  : ALAIN PERRONNET  TEXAS A & M UNIVERSITY            JULY  2005
C MODIF  : ALAIN PERRONNET  Laboratoire J.L. LIONS UPMC Paris Mars  2007
C MODIFS : ALAIN PERRONNET  Saint PIERRE du PERRAY            Mars  2021
C2345X7..............................................................012
C     RUSE POUR CALCULER 2**31-1 PLUS GRAND ENTIER STOCKABLE SUR 32 BITS
C                                        -2 pour la FIN de BOUCLE
      PARAMETER         ( MAXENTIER=2**30-2+2**30 )
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___union.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a_objet__nouvnoeud.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     KNOMOB
      INTEGER           NBELEM(9), MOELEM(9), MNELEM(9), LIENEL(9),
     %                  NUDREL(9), NTNPEF(9), MNNPEF(9), NMNPEF(9)
      INTEGER           NBELTG(9), MOELTG(9), MNELTG(9), NUELTG(9)
      LOGICAL           AVANT
      CHARACTER*10      NMTYOB,KNOMTY
      CHARACTER*24      KNOM
      INTEGER, allocatable, dimension(:) :: LPTDVOI, LISTVOI


C     MOELEM = NOMBRE DE SOMMETS + LIEN + 1EF + 1EF A TG
      DATA              MOELEM/ 4, 5, 6, 7, 7, 9, 11, 65, 8 /
      DATA              LIENEL/ 2, 3, 4, 5, 5, 7,  9,  0, 6 /
C     MOELTG = NOMBRE DE TANGENTES SELON LE TYPE DE L'EF
      DATA              MOELTG/ 1, 2, 6, 8, 12, 18, 24, 0, 15 /

      IERR   = 0
      PRINT*,'dftop0: 2**31-2=MAXENTIER=',MAXENTIER

C     TABLEAU LISTVOI NON ALLOUE
      IERALISTVOI = 1
      IERALPTDVOI = 1
      NTTOPO = 0
      MNOBPR = 0
      MNMNSO = 0
      MNOBIN = 0
      MNRENU = 0

      CALL AZEROI( 9, MNELEM )
      CALL AZEROI( 9, MNNPEF )
      CALL AZEROI( 9, NTNPEF )

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
         IERR = 1
         GOTO 9999
      ENDIF

C     TRACE EVENTUEL DU NOM DE L'OBJET DANS L'HISTORIQUE
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO

C     RECHERCHE DU TABLEAU DEFINITION
      CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )

C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: DEFINITION INCONNUE OBJET ' // KNOMOB
         ELSE
            KERR(1) = 'ERROR: UNKNOWN DEFINITION for OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         IERR = 2
         GOTO 9999
      ENDIF

      L = NUDCNB( KNOMOB )
      WRITE(IMPRIM,10000) KNOMOB(1:L)
10000 FORMAT(/' CREATION TMS: ~>OBJET>',A,'>TOPOLOGIE',/,1X,79('='))
C
C     L'OBJET DOIT-IL ETRE STRUCTURE EN SOUS-DOMAINES ou JOINTS?
C     ==========================================================
cccC     TRAITEMENT CLASSIQUE FORCE  ...... 3 aout 2006
ccc      IF( MCN( MNDFOB + WDOUNO ) .GT. 0 ) MCN( MNDFOB + WDOUNO ) = 0
ccc
      NDOUNO = MCN( MNDFOB + WDOUNO )
      IF( NDOUNO .EQ. 1 ) THEN
C        TRAITEMENT PAR SOUS-DOMAINES
         GOTO 500
      ELSE IF( NDOUNO .EQ. 2 ) THEN
C        TRAITEMENT PAR JOINTS
         GOTO 510
      ENDIF

C     ==================================================================
C     TRAITEMENT CLASSIQUE AVEC OU SANS RENUMEROTATION DES NOEUDS
C     ==================================================================
C     LE TABLEAU TOPOLOGIE EST RECHERCHE
      CALL LXTSOU( NTLXOB, 'TOPOLOGIE', NTTOPO, MNTOPO )
      IF( NTTOPO .GT. 0 ) THEN
C        SI LA DEFINITION A ETE MODIFIEE ALORS DESTRUCTION
C        DU TABLEAU TOPOLOGIE DE CET OBJET ET CALCUL
         IF( AVANT( MCN(MNTOPO), MCN(MNDFOB) ) ) THEN
C           DESTRUCTION DES TABLEAUX XYZNOEUD XYZPOINT NPEF TOPOLOGIE
            CALL DSNPET( NTLXOB )
            NTTOPO = 0
            MNTOPO = 0
            NBLGRC(NRERR) = 3
            KERR(1) =  KNOMOB
            IF( LANGAG .EQ. 0 ) THEN
               KERR(2) = 'L''OBJET A ETE MODIFIE'
               KERR(3) = 'SON INTERPOLATION EST RECALCULEE'
            ELSE
               KERR(2) = 'The OBJECT HAS BEEN MODIFIED'
               KERR(3) = 'Its INTERPOLATION IS COMPUTED AGAIN'
            ENDIF
            CALL LERESU
         ELSE
            NBLGRC(NRERR) = 4
            KERR(1) =  KNOMOB
            IF( LANGAG .EQ. 0 ) THEN
               KERR(2) = 'L''OBJET N''A PAS ETE MODIFIE'
               KERR(3) = 'SON INTERPOLATION EXISTE DEJA'
               KERR(4) = 'EST ELLE A RECALCULER?'
            ELSE
               KERR(2) = 'The OBJECT HAS NOT BEEN MODIFIED'
               KERR(3) = 'Its INTERPOLATION ALREADY EXISTS'
               KERR(4) = 'WILL YOU COMPUTE IT AGAIN?'
            ENDIF
            CALL LERESU
            CALL LIMTCL( 'non_oui', N )
            IF( N .EQ. 0 ) GOTO 9999
C           DESTRUCTION DES TABLEAUX XYZNOEUD XYZPOINT NPEF TOPOLOGIE
            CALL DSNPET( NTLXOB )
         ENDIF
      ENDIF

C     LES OBJETS DE L'OBJET SONT EXPRIMES EN TERMES DE PLSV SANS OBJETS
C     EXCLUSIVEMENT C'EST A DIRE DE PLSV 'PREMIERS'
C     ATTENTION: LES TYPES SONT TRIES DANS L'ORDRE V, S, L, P
C     =================================================================
      CALL OBJPRE( KNOMOB, NBOBPR, MOOBPR, MNOBPR, IERR )
      IF( IERR .NE. 0 ) GOTO 9990

C     UNION DES TABLEAUX XYZSOMMET DES OBJETS PREMIERS DE L'OBJET
C     IDENTIFICATION DES SOMMETS   ET CREATION D'UNE SEULE NUMEROTATION
C     IDENTIFICATION DES TANGENTES ET CREATION D'UNE SEULE NUMEROTATION
C     =================================================================
      CALL NUOBNM( 'OBJET', KNOMOB, NUOBJE )
      CALL UNIOBJ( NUOBJE, NBOBPR, MCN(MNOBPR),
     %             NTSOMM, MNSOMM, NBSOM,
     %             NTUNIO, MNUNIO, IERR   )
      IF( IERR .NE. 0 ) GOTO 9990

C     RECHERCHE DE LA VRAIE DIMENSION DE L'ESPACE DES COORDONNEES
      NBCOOR = MCN(MNSOMM+WBCOOR)
      IF( NBCOOR .LE. 3 ) THEN
         CALL DIMCOO( NBSOM, MCN(MNSOMM+WYZSOM), NDIMEN )
      ELSE
         NDIMEN = NBCOOR
      ENDIF

C     OUVERTURE DES TABLEAUX NSEF DE CES NBOBPR OBJETS
      CALL TNMCDC( 'ENTIER', NBOBPR, MNMNSO )

C     BOUCLE SUR LES NO SOMMET DES OBJETS PREMIERS
      DO 20 N = 0, NBOBPR-1
C        LE NUMERO DU TYPE DE L'OBJET PREMIER
         MNOB   = MNOBPR + N + N
         NUTYOB = MCN( MNOB )
C        LE NUMERO DE L'OBJET DANS CE TYPE
         NUOBJE = MCN( MNOB + 1 )

C        LE LEXIQUE DE CE TYPE D'OBJETS
         NTLX = NTMN( NUTYOB )

C        LE LEXIQUE DE CET OBJET
         CALL LXNLOU( NTLX, NUOBJE, NTLXO, NT )

C        LE TABLEAU NSEF DE CET OBJET EXISTE-T-IL ?
         CALL LXTSOU( NTLXO, 'NSEF', NT, MCN(MNMNSO+N) )
         IF( NT .LE. 0 ) THEN
            IF( NUTYOB .EQ. 1 ) THEN
C              L'OBJET EST UN POINT : SON TABLEAU NSEF EST CREE
C              AFIN DE RENDRE POINTS LIGNES, ... HOMOGENES
               CALL LXTSOU( NTLXO, 'XYZSOMMET', NT, L )
C              LE NOMBRE DE SOMMETS DE CE POINT
               J = MCN( L + WNBSOM )
               CALL LXTNDC( NTLXO, 'NSEF', 'MOTS', WBARNS+1 )
               CALL LXTSOU( NTLXO, 'NSEF', NT, L )
               MCN( MNMNSO + N ) = L
C              LE TYPE DE L'OBJET = POINT
               MCN( L + WUTYOB ) = 1
C              POINT => INCONNU NI OUVERT NI FERME
               MCN( L + WUTFMA ) = -1
C              LE NOMBRE DE SOMMETS ET TANGENTES PAR NOEUD
               MCN( L + WBSOEF ) = 1
               MCN( L + WBTGEF ) = 0
C              LE NOMBRE DE NOEUDSOMMETS
               MCN( L + WBEFOB ) = J
               MCN( L + WBEFTG ) = 0
               MCN( L + WBEFAP ) = 0
C              LE TYPE DE MAILLAGE : STRUCTURE
               MCN( L + WUTYMA ) = 1
C              LE NOMBRE DE NOEUDSOMMETS
               MCN( L + WBARNS ) = J
C              LA DATE
               CALL ECDATE( MCN(L) )
C              LE NUMERO DU TABLEAU DESCRIPTEUR
               MCN( L + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
            ELSE
C              LE NOM DE L'OBJET
               KNOMTY = NMTYOB( NUTYOB )
               CALL NMOBNU( KNOMTY, NUOBJE, KNOM )
               NBLGRC(NRERR) = 2
               KERR(1) = KNOM
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(2) = 'OBJET SANS NSEF => A MAILLER AVANT'
               ELSE
                  KERR(2) = 'OBJECT WITHOUT NSEF => MESH IT BEFORE'
               ENDIF
               CALL LEREUR
               IERR = 1
               GOTO 20
            ENDIF
         ENDIF
 20   ENDDO
      IF( IERR .NE. 0 ) GOTO 9990

C     CALCUL DU NOMBRE DE NOEUDSOMMETS,ARETES,TRIANGLES,QUADRANGLES,
C     TETRAEDRES,PENTAEDRES,HEXAEDRES 6-CUBES ET
C     DES EF A TG DE CHAQUE TYPE
C     ==============================================================
      CALL DFTOP1( NBOBPR, MCN(MNMNSO), NBELEM, NBELTG, IERR )
      IF( IERR .NE. 0 ) GOTO 9990
      IF( NBELEM(9) .GT. 0 ) THEN
         NBLGRC(NRERR) = 3
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'OBJET ' // KNOM
            KERR(2) = 'OBJET AVEC PYRAMIDES NON PROGRAMME'
            KERR(3)='LES DECOUPER EN 2 TETRAEDRES OPTION 31 des VOLUMES'
         ELSE
            KERR(1) = 'OBJECT ' // KNOM
            KERR(2) = 'OBJECT WITH PYRAMIDS NOT PROGRAMED'
            KERR(3) = 'CUT THEM INTO 2 TETRAHEDRA BY VOLUME OPTION 31'
         ENDIF
         CALL LEREUR
         IERR = 2
         GOTO 9990
      ENDIF

C     DECLARATION SANS INITIALISATION DES TABLEAUX DE HACHAGE DES
C     NOEUDSOMMETS, SEGMENTS , TRIANGLES, QUADRANGLES,
C     TETRAEDRES , PENTAEDRES, HEXAEDRES, 6-CUBES
C     DECLARATION DES TABLEAUX DES NUMEROS DES TG DES EF A TG
C     ===========================================================
      CALL DFTOP2( MOELEM, NBELEM, MNELEM,
     %             MOELTG, NBELTG, MNELTG )

C     DECLARATION DES TABLEAUX NUOBIN, NUOBCL
      CALL TNMCDC( 'ENTIER', NBOBPR * 2, MNOBIN )
      CALL AZEROI( NBOBPR * 2, MCN(MNOBIN) )
      MNOBCL = MNOBIN + NBOBPR

C     LES SOMMETS DES POINTS, LES ARETES DES LIGNES,
C     LES TRIANGLES DES SURFACES, LES QUADRANGLES DES SURFACES,
C     LES TETRAEDRES DES VOLUMES, LES PENTAEDRES DES VOLUMES,
C     LES HEXAEDRES DES VOLUMES PREMIERS, LES 6-CUBES
C     SONT LISTES PAR HACHAGE
C     SANS RECOUPEMENT ENTRE OBJETS PREMIERS DE TYPES DIFFERENTS
C     ==========================================================
      MNDESO = MNUNIO + WDSCOU
      MNNUSO = MNDESO + 1 + NBOBPR
      MNDETG = MNNUSO + MCN(MNDESO+NBOBPR)
      MNNUTG = MNDETG + 1 + NBOBPR
      IF( NBELEM(8) .EQ. 0 ) THEN
C        CONSTRUCTION DES HACHAGES DES DIFFERENTS EF
         CALL DFTOP3( MOELEM, NBELEM, MNELEM,
     %                MOELTG, NBELTG, MNELTG, NUELTG,
     %                LIENEL, NUDREL,
     %                NBOBPR, MCN(MNOBPR), MCN(MNMNSO),
     %                MCN(MNDESO), MCN(MNNUSO),
     %                MCN(MNDETG), MCN(MNNUTG),
     %                NDIMEN, MCN(MNSOMM+WYZSOM),
     %                IERR )
         IF( IERR .NE. 0 ) GOTO 9990
      ELSE
C        CONSTRUCTION DIRECTE DU TMS NPEF"6Q1C DES 6-CUBES
         CALL DFT36C( NTLXOB, NBELEM, MNELEM,
     %                NBOBPR, MCN(MNMNSO),
     %                MCN(MNDESO), MCN(MNNUSO),
     %                NTNPEF, MNNPEF, NMNPEF,
     %                MCN(MNOBIN), IERR )
      ENDIF

C     LECTURE DU TYPE DE L'INTERPOLATION A EFFECTUER
C     ==============================================
      CALL LIMTCL( 'interpol', NOINTE )
C     NOINTE = 1 : AXISYMETRIQUE de DEGRE 1
C              2 : AXISYMETRIQUE de DEGRE 2
C              3 : LAGRANGE de DEGRE 1
C              4 : LAGRANGE de DEGRE 2
C             -1 : ABANDON DEMANDE
      IF( NOINTE .EQ. -1 ) THEN
C        ABANDON DEMANDE
         IERR = NOINTE
         RETURN
      ENDIF
      IF( NBELEM(8) .GT. 0 ) THEN
C        6-CUBES => 6Q1C interpolation LAGRANGE de DEGRE 1
         NOINTE = 3
         GOTO 40
      ENDIF

      IF( NOINTE .LE. 2 ) THEN
C
C        VERIFICATION QUE LA COORDONNEE R GARDE TOUJOURS LE MEME SIGNE
         MN = MNSOMM+WYZSOM
C        RECHERCHE DU PREMIER RAYON NON NUL
         DO I=1,NBSOM
            IF( RMCN(MN) .NE. 0.0 ) GOTO 126
            MN = MN + 3
         ENDDO
 126     R = RMCN(MN)
C        EXISTE T IL UN RAYON DE SIGNE DIFFERENT?
         DO J=I+1,NBSOM
            MN = MN + 3
            IF( R * RMCN(MN) .LT. 0.0 ) THEN
C              2 RAYONS DE SIGNES DIFFERENTS
               NBLGRC(NRERR) = 4
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'EN AXISYMETRIE LES RAYONS R'
                  KERR(2) = 'DOIVENT TOUS AVOIR LE MEME SIGNE'
                  WRITE(KERR(6)( 1:15),'(G15.7)') R
                  WRITE(KERR(6)(16:30),'(G15.7)') RMCN(MN)
                  KERR(3) = 'ICI ' // KERR(6)(1:15) // ' ET '
     %                             // KERR(6)(16:30)
                  KERR(4) = 'CORRIGER LE SIGNE DES RAYONS R'
               ELSE
                  KERR(1) = 'FOR AXISYMMETRY PROBLEM the RADII R'
                  KERR(2) = 'MUST HAVE ALL THE SAME SIGN'
                  WRITE(KERR(6)( 1:15),'(G15.7)') R
                  WRITE(KERR(6)(16:30),'(G15.7)') RMCN(MN)
                  KERR(3) = 'HERE ' // KERR(6)(1:15) // ' and '
     %                              // KERR(6)(16:30)
                  KERR(4) = 'CORRECT the SIGN of RADII R'
               ENDIF
               CALL LEREUR
               IERR = 5
               GOTO 9990
            ENDIF
         ENDDO
      ENDIF

C     FORMATION( DES TABLEAUX 'NPEF"TYPE_EF' PAR TYPE D'EF

C     EN 3D: UN POUR LES HEXAEDRES, UN POUR LES PENTAEDRES, UN POUR TETRAEDRES
C     EN 2D: UN POUR LES QUADRANGLES, UN POUR LES TRIANGLES

C     ATTENTION: APRES TRAITEMENT DE TOUS LES VOLUMES, LES EF DE DIMENSION
C     INFERIEURE (QUADRANGLES, ..., SOMMETS) ONT DU TOUS ETRE RETROUVES
C     COMME DES LIMITES (FACES, ARETES, SOMMETS) DES EF DES VOLUMES
C     DANS LE CAS CONTRAIRE, IL Y A ERREUR DE DONNEE CAR AUCUNE CONDITION AUX
C     LIMITES NE POURRA ETRE DEFINIE SUR CES LAISSES POUR COMPTE
C     ==========================================================================
      CALL DFTOP4( NTLXOB, MOELEM, NBELEM, MNELEM,
     %             MOELTG, NBELTG, MNELTG, NUELTG, LIENEL,
     %             NOINTE, NBOBPR, MCN(MNOBPR), MCN(MNSOMM+WYZSOM),
     %             NTNPEF, MNNPEF, NMNPEF, MCN(MNOBIN), MCN(MNOBCL),
     %             IERR   )
      IF( IERR .NE. 0 ) GOTO 9990

C     LE NOMBRE DE TABLEAUX 'NPEF'
 40   NBTYEL = 0
      DO I=1,9
         IF( NTNPEF(I) .GT. 0 ) NBTYEL = NBTYEL + 1
      ENDDO

C     CALCUL DE NBOBIN NBOBCL
      NBOBIN = 0
      NBOBCL = 0
      DO I=0,NBOBPR-1
         IF( MCN(MNOBIN+I) .GT. 0 ) NBOBIN = NBOBIN + 1
         IF( MCN(MNOBCL+I) .GT. 0 ) NBOBCL = NBOBCL + 1
         IF( MCN(MNOBIN+I) .GT. 0 .AND.
     %       MCN(MNOBCL+I) .GT. 0 ) THEN
             J      = MNOBPR + I + I
             KNOMTY = NMTYOB( MCN(J) )
             CALL NMOBNU( KNOMTY, MCN(J+1), KNOM )
             NBLGRC(NRERR) = 2
             IF( LANGAG .EQ. 0 ) THEN
                KERR(1) = 'ERREUR: ' // KNOMTY // ' '  // KNOM
                KERR(2) = 'A LA FOIS INTERNE ET AUX LIMITES'
             ELSE
                KERR(1) = 'ERROR: ' // KNOMTY // ' '  // KNOM
                KERR(2) = 'INTERNAL and ON THE BOUNDARY'
             ENDIF
             CALL LEREUR
             IERR = 1
         ENDIF
         IF( MCN(MNOBIN+I) .LE. 0 .AND. MCN(MNOBCL+I) .LE. 0 ) THEN
             J      = MNOBPR + I + I
             KNOMTY = NMTYOB( MCN(J) )
             CALL NMOBNU( KNOMTY, MCN(J+1), KNOM )
             NBLGRC(NRERR) = 3
             IF( LANGAG .EQ. 0 ) THEN
                KERR(1) = 'ERREUR: ' // KNOMTY // ' ' // KNOM
                KERR(2) = 'NI MATERIAU INTERNE NI SUR LA FRONTIERE'
                KERR(3) = 'EN FAIT NON RETROUVE DANS LE MAILLAGE'
             ELSE
                KERR(1)= 'ERROR: ' // KNOMTY // ' ' // KNOM
                KERR(2)= 'NEITHER INTERNAL MATERIAL NOR ON THE BOUNDARY'
                KERR(3)= 'IN FACT NOT RETRIEVED IN THE MESH'
             ENDIF
             CALL LEREUR
             IERR = 1
         ENDIF
      ENDDO

      IF( NBOBIN .LE. 0 ) THEN
         NBLGRC(NRERR) = 3
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET ' // KNOMOB
            KERR(2) = 'SANS OBJET INTERNE c-a-d MATERIAU'
            KERR(3) = '=> PAS DE CARACTERISTIQUES PHYSIQUES IMPOSABLES'
         ELSE
            KERR(1) = 'ERROR: OBJECT ' // KNOMOB
            KERR(2) = 'WITHOUT MATERIAL'
            KERR(3) = '=> ZERO PHYSICAL CHARACTERISTICS to IMPOSE'
         ENDIF
         CALL LEREUR
         IERR = 2
      ENDIF

      IF( NBOBCL .LE. 0 .AND. NBELEM(8) .EQ. 0 ) THEN
         NBLGRC(NRERR) = 3
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET ' // KNOMOB
            KERR(2) = 'SANS POINT ou LIGNE (ou SURFACE) AUX LIMITES'
            KERR(3) = '=> PAS DE CONDITIONS AUX LIMITES IMPOSABLES'
         ELSE
            KERR(1) = 'ERROR: OBJECT ' // KNOMOB
          KERR(2) = 'WITHOUT POINT or LINE (or SURFACE) at the BOUNDARY'
            KERR(3) = '=> NO BOUNDARY CONDITIONS to IMPOSE'
         ENDIF
         CALL LEREUR
         IERR = 3
      ENDIF
      IF( IERR .NE. 0 ) GOTO 9990

C     FORMATION DU TABLEAU 'TOPOLOGIE'
C     ================================
      L = MOTVAR( NUMTYP('TYPEOBJET') )
      J = WMTYEL + NBTYEL + (NBOBIN+NBOBCL) * L
      CALL LXTNDC( NTLXOB, 'TOPOLOGIE', 'MOTS', J )
      CALL LXTSOU( NTLXOB, 'TOPOLOGIE', NTTOPO, MNTOPO )

C     INITIALISATION DU TABLEAU TOPOLOGIE
C     LE NOMBRE DES OBJETS INTERNES ET AUX LIMITES
      MCN( MNTOPO + WBOBIN ) = NBOBIN
      MCN( MNTOPO + WBOBCL ) = NBOBCL
C     LE NOMBRE DE TYPE D'ELEMENTS FINIS
      MCN( MNTOPO + WBTYEL ) = NBTYEL
C     POUR L'INSTANT: NOEUDS = POINTS = SOMMETS
      MCN( MNTOPO + WDPGST ) = 0
C     LES NOMS DES TYPES D'ELEMENTS FINIS
      J = MNTOPO + WMTYEL
      DO I=1,9
         IF( NTNPEF(I) .GT. 0 ) THEN
            MCN( J ) = NMNPEF( I )
            J = J + 1
         ENDIF
      ENDDO

      MNOBI = MNTOPO + WMTYEL + NBTYEL
      MNOBC = MNOBI  + L * NBOBIN
      DO I=0,NBOBPR - 1

C        LE RECENSEMENT DES OBJETS INTERNES
         IF( MCN(MNOBIN+I) .GT. 0 ) THEN
C           LE TYPE DE L'OBJET
            MCN( MNOBI  ) = MCN( MNOBPR+I+I)
C           LE NUMERO DE L'OBJET
            MCN( MNOBI+1) = MCN( MNOBPR+I+I+1)
            MNOBI = MNOBI + L
         ENDIF

C        LE RECENSEMENT DES OBJETS AUX LIMITES
         IF( MCN(MNOBCL+I) .GT. 0 ) THEN
C           LE TYPE DE L'OBJET
            MCN( MNOBC  ) = MCN( MNOBPR+I+I)
C           LE NUMERO DE L'OBJET
            MCN( MNOBC+1) = MCN( MNOBPR+I+I+1)
            MNOBC = MNOBC + L
         ENDIF
      ENDDO

C     LA DATE
      CALL ECDATE( MCN(MNTOPO) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNTOPO + MOTVAR(6) ) = NONMTD( '~>OBJET>>TOPOLOGIE' )

C     DESTRUCTION DES TABLEAUX DE L'OBJET DEVENUS INUTILES
      IF( MNOBPR .GT. 0 ) CALL TNMCDS( 'ENTIER', MOOBPR, MNOBPR )
      IF( MNMNSO .GT. 0 ) CALL TNMCDS( 'ENTIER', NBOBPR, MNMNSO )
C     DESTRUCTION DES TABLEAUX DE HACHAGE DEVENUS INUTILES
      DO I=1,9
         IF( MNELEM(I) .GT. 0 ) THEN
            CALL TNMCDS( 'ENTIER', MOELEM(I)*NBELEM(I), MNELEM(I) )
            MNELEM(I) = 0
         ENDIF
         IF( MNELTG(I) .GT. 0 ) THEN
            CALL TNMCDS( 'ENTIER', MOELTG(I)*NBELTG(I), MNELTG(I) )
            MNELTG(I) = 0
         ENDIF
      ENDDO
      IF( MNOBIN .GT. 0 ) CALL TNMCDS( 'ENTIER',NBOBPR*2,MNOBIN )

C     AJOUT DU NUMERO DES NOEUDS NON SOMMETS
C     AJOUT DU NUMERO DES POINTS S'ILS SONT DIFFERENTS DES NOEUDS
C     AJOUT DES TABLEAUX NOEUDS ET POINTS GEOMETRIQUES
C     ===========================================================
      IF( IERR .EQ. 0 ) THEN
         CALL DFTOP5( NTLXOB, NDIMEN, NOINTE,
     %                NBELEM, NTNPEF,
     %                NBSOM,  MNSOMM,
     %                NTXYZN, MNXYZN, NTXYZP, MNXYZP,
     %                NDPGST, IERR   )
      ENDIF
      IF( IERR .GT. 0 ) GOTO 9990
C     MISE A JOUR DANS TOPOLOGIE DU CODE DE TRAITEMENT
C     DES TABLEAUX NOEUDS POINTS GEOMETRIQUES SOMMETS
      IF( NDPGST .GT. 0 ) THEN
         MCN( MNTOPO + WDPGST ) = NDPGST
      ENDIF

C     RENUMEROTATION DES NOEUDS DES ELEMENTS FINIS
C     INITIALISATION,NUMEROS DE TMS ET DES ADRESSES MCN
C     =================================================
      CALL NDPGEL( NTLXOB, NTTOPO, MNTOPO,
     &             NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     &             NBTYEL, NTNPEF, MNNPEF, IERR )
      IF( IERR.NE.0 ) GOTO 9990

C     TABLEAU DES NOEUDS RENUMEROTES
      CALL LXTSOU( NTLXOB, 'NOUVNOEUD', NTNONO, MNNONO )
      IF( NTNONO .GT. 0 ) THEN
C        LE TABLEAU EXISTANT EST DETRUIT
         CALL LXTSDS( NTLXOB, 'NOUVNOEUD' )
      ENDIF
C     NOMBRE TOTAL DE NOEUDS
      NBNOE = MCN( MNXYZN + WNBNOE )
C     UN MOT DE PLUS EST NECESSAIRE DANS GIBB9 LIGNES 145 et 150
C     POUR LE TABLEAU DE RENUMEROTATION DE TAILLE NBNOE
C     MAIS DONT LA VARIABLE NBNOE+1 EST UTILISEE (INDICE 0)!
      CALL LXTNDC( NTLXOB, 'NOUVNOEUD', 'MOTS', WONOEU+NBNOE+1 )
      CALL LXTSOU( NTLXOB, 'NOUVNOEUD', NTNONO, MNNONO )
C     ADRESSE DE DEBUT DES NOUVEAUX NUMEROS DES NOEUDS
      MNRENU = MNNONO + WONOEU

ccc      IF(NDIMEN .GT. 1 .AND. NDOUNO .GE. 0 .AND. NBELEM(8) .LE. 0) THEN
      IF( NDOUNO .GE. 0 .AND. NBELEM(8) .LE. 0 ) THEN

C        CONSTRUCTION DE LA LISTE DES NOEUDS VOISINS DES NOEUDS
C        ALLOCATION DU TABLEAU LPTDVOI DES POINTEURS SUR LE DERNIER VOISIN
C        =================================================================
         NBNOE1 = 1 + NBNOE
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'dftop0: DEMANDE  ALLOCATION LPTDVOI(',NBNOE1,
     %             ') ENTIERS'
            ALLOCATE ( LPTDVOI( 0:NBNOE ), STAT=IERALPTDVOI )
            IF( IERALPTDVOI .NE. 0 ) THEN
               PRINT*,'dftop0: ERREUR ALLOCATION LPTDVOI(',NBNOE1,
     %                ') ENTIERS'
               IERR = IERALPTDVOI
               GOTO 9990
            ENDIF
            PRINT*,'dftop0: CORRECTE ALLOCATION LPTDVOI(',NBNOE1,
     %             ') ENTIERS'
         ELSE
            PRINT*,'dftop0: ALLOCATION DEMAND  LPTDVOI(',NBNOE1,
     %             ') INTEGERS'
            ALLOCATE ( LPTDVOI( 0:NBNOE ), STAT=IERALPTDVOI )
            IF( IERALPTDVOI .NE. 0 ) THEN
               PRINT*,'dftop0: ALLOCATION ERROR LPTDVOI(',NBNOE1,
     %                ') INTEGERS'
               IERR = IERALPTDVOI
               GOTO 9990
            ENDIF
            PRINT*,'dftop0: CORRECT ALLOCATION LPTDVOI(',NBNOE1,
     %             ') INTEGERS'
         ENDIF

C        CREER LA LISTE DU NOMBRE D'ELEMENTS FINIS DE CHAQUE NOEUD
         CALL NBEFNOEUD( MNTOPO, MNNPEF, NBNOE,
     %                   LPTDVOI(1), MXNOE1EF, IERR )

C        NOINTE=1 : AXISYMETRIQUE DEGRE 1 <-> NOEUDS=SOMMETS
C               2 : AXISYMETRIQUE DEGRE 2 <-> NOEUDS=SOMMETS+MILIEUX
C               3 : LAGRANGE de DEGRE 1   <-> NOEUDS=SOMMETS
C               4 : LAGRANGE de DEGRE 2   <-> NOEUDS=SOMMETS+MILIEUX

C        ESTIMATION DU NOMBRE MOYEN DE NOEUDS VOISINS DANS CHAQUE EF
         NBNOVO = MXNOE1EF/2

C        ALLOCATION DU TABLEAU LISTVOI DES NOEUDS VOISINS
         NBAGMX = 0

 60      NBAGMX = NBAGMX + 1
         IF( NBAGMX .GT. 6 ) THEN
C           TROP D'AUGMENTATIONS DE LA VALEUR DE NBNOVO
            NBLGRC(NRERR) = 3
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='dftop0: TROP d''AUGMENTATIONS de NBNOVO'
               KERR(2)='POUR DECLARER LE TABLEAU DES VOISINS DES NOEUDS'
               KERR(3)='REDUIRE le MAILLAGE'
            ELSE
               KERR(1)='dftop0: TOO AUGMENTATIONS of NBNOVO'
               KERR(2)='to DECLARE the ARRAY of NODE NEIGHBORS of NODES'
               KERR(3)='REDUCE the MESH'
            ENDIF
            CALL LEREUR
            IERR = 2
            GOTO 9990
         ELSE
            NBNOVO = NBNOVO + 1
            PRINT*,'dftop0: NB VOISINS 1 EF NBNOVO=',NBNOVO,
     %             'MAX NOEUDS 1EF=',MXNOE1EF
         ENDIF

         IF( NBAGMX .GT. 1 ) THEN
C           RECREER LA LISTE DU NOMBRE D'ELEMENTS FINIS DE CHAQUE NOEUD
            CALL NBEFNOEUD( MNTOPO, MNNPEF, NBNOE,
     %                      LPTDVOI(1), MXNOE1EF, IERR )
         ENDIF

C        ESTIMATION DU NOMBRE DE VOISINS DE CHAQUE NOEUD A PARTIR DU
C        NOMBRE D'EF DE CHAQUE NOEUD ET LE NOMBRE MOYEN DE VOISINS
         LPTDVOI(0) = 0
         DO K = 1, NBNOE
C           NOMBRE D'EF DU NOEUD K
            NBEF = LPTDVOI(K)
C           NOMBRE MOYEN DE VOISINS PAR EF
            IF( NBEF .LE. 6 ) THEN
               NBNV = MXNOE1EF - 1
            ELSE
               NBNV = NBNOVO
            ENDIF
C           ESTIATION DU NOMBRE DE NOEUDS VOISINS DU NOEUD K
            LPTDVOI(K) = LPTDVOI(K-1) + NBNV * NBEF
         ENDDO
         MXLISTV = LPTDVOI(NBNOE)

         IF( MXLISTV .LE. 0 ) THEN
C           MXLISTV > 2**31-1
            NBLGRC(NRERR) = 3
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='dftop0: SATURATION DU MAX DES ENTIERS'
               KERR(2)='POUR DECLARER LE TABLEAU DES VOISINS DES NOEUDS'
               KERR(3)='REDUIRE le MAILLAGE'
            ELSE
               KERR(1)='dftop0: SATURATION of INTEGER MAXIMUM'
               KERR(2)='to DECLARE the ARRAY of NODE NEIGHBORS of NODES'
               KERR(3)='REDUCE the MESH'
            ENDIF
            CALL LEREUR
            IERR = 2
            GOTO 9990
         ENDIF

C        DECLARATION DU TABLEAU LISTVOI(MXLISTV) DES NOEUDS VOISINS
C        SI TAILLE<MAXENTIER=2**31-1
         IF( MXLISTV .LE. 0 .OR. MXLISTV .GE. MAXENTIER ) THEN
C           MXLISTV EST TROP GRAND POUR ETRE DECLARE
            PRINT*,'dftop0: declaration de LISTVOI(',MXLISTV,')'
            NBLGRC(NRERR) = 3
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: TAILLE DE LA LISTE DES VOISINS avec'
               KERR(2) = 'DEPASSEMENT DU MAX DES ENTIERS 2**31-1'
               KERR(3) = 'REDUIRE LA TAILLE DU MAILLAGE'
            ELSE
               KERR(1) = 'ERROR: NEIGHBORS LIST SIZE'
               KERR(2) = 'OVERFLOW of INTEGER MAX 2**31-1'
               KERR(3) = 'REDUCE THE SIZE OF THE MESH'
            ENDIF
            CALL LERESU
            IERR = 2
            GOTO 9990
         ENDIF

         IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'DEMANDE  ALLOCATION LISTVOI(',MXLISTV,') ENTIERS'
         ALLOCATE ( LISTVOI( 1:MXLISTV ), STAT=IERALISTVOI )
         IF( IERALISTVOI .NE. 0 ) THEN
            PRINT*,'ERREUR ALLOCATION LISTVOI(',MXLISTV,') ENTIERS'
            IERR = IERALISTVOI
            GOTO 9990
         ENDIF
         PRINT*,'CORRECTE ALLOCATION LISTVOI(',MXLISTV,') ENTIERS'
         ELSE
         PRINT*,'ALLOCATION DEMAND  LISTVOI(',MXLISTV,') INTEGERS'
         ALLOCATE ( LISTVOI( 1:MXLISTV ), STAT=IERALISTVOI )
         IF( IERALISTVOI .NE. 0 ) THEN
            PRINT*,'ALLOCATION ERROR LISTVOI(',MXLISTV,') INTEGERS'
            IERR = IERALISTVOI
            GOTO 9990
         ENDIF
         PRINT*,'CORRECT ALLOCATION LISTVOI(',MXLISTV,') INTEGERS'
         ENDIF

         CALL LISTNOVO( MNTOPO, MNNPEF, NBNOE, MXLISTV,
     %                  LPTDVOI(0), LISTVOI, MXVOIS, IERR )

         IF( IERR .EQ. 5 ) THEN
C           REALLOCATION DE LISTVOI DEMANDEE
            IF( IERALISTVOI .EQ. 0 ) then
               IF( LANGAG .EQ. 0 ) THEN
                  PRINT*,'LISTVOI est DESALLOUE'
               ELSE
                  PRINT*,'LISTVOI is DESALLOCATED'
               ENDIF
               DEALLOCATE( LISTVOI )
               IERALISTVOI = 1
            ENDIF
            GOTO 60
         ENDIF

         IF( IERR .NE. 0 ) GOTO 9990

C        RENUMEROTATION DES NOEUDS DES EF SELON L'ALGORITHME DE GIBBS
C        ============================================================
         CALL GIBBS( NBNOE, LPTDVOI(0), LISTVOI, MNRENU )

C        RENUMEROTATION EFFECTIVE DES NOEUDS DES EF
         CALL GIBBA( MNTOPO, MNNPEF, MNRENU, IERR )

C        RENUMEROTATION EFFECTIVE DES XYZ DES NOEUDS DES EF
         IF( IERR .EQ. 0 ) CALL GIBBB( MNXYZN, MNRENU )

      ELSE

C        PAS DE RENUMEROTATION DES NOEUDS POUR LES
C        6-CUBES MAILLAGE STRUCTURE DONC NON NECESSAIRE!
C        ===============================================
         DO I = 1, NBNOE+1
            MCN( MNRENU-1 + I ) = I
         ENDDO

      ENDIF

      IF( IERALISTVOI .EQ. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'LISTVOI est DESALLOUE'
         ELSE
            PRINT*,'LISTVOI is DESALLOCATED'
         ENDIF
         DEALLOCATE( LISTVOI )
         IERALISTVOI = 1
      ENDIF

      IF( IERALPTDVOI .EQ. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'LPTDVOI est DESALLOUE'
         ELSE
            PRINT*,'LPTDVOI is DESALLOCATED'
         ENDIF
         DEALLOCATE( LPTDVOI )
         IERALPTDVOI = 1
      ENDIF

C     SAUVEGARDE DE LA RENUMEROTATION DANS LE TMS NOUVNOEUD
C     =====================================================
      MCN( MNNONO + WBNONO ) = NBNOE
C     LA DATE
      CALL ECDATE( MCN(MNNONO) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNNONO + MOTVAR(6) ) = NONMTD( '~>OBJET>>NOUVNOEUD' )

C     AFFICHAGES
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'NOMBRE de SOMMETS       =',NBSOM
         WRITE(IMPRIM,*) 'NOMBRE de NOEUDS        =',NBNOE
         WRITE(IMPRIM,*) 'NOMBRE de SOMMETS-POINTS=',NBELEM(1)
         WRITE(IMPRIM,*) 'NOMBRE de SEGMENTS      =',NBELEM(2)
         WRITE(IMPRIM,*) 'NOMBRE de TRIANGLES     =',NBELEM(3)
         WRITE(IMPRIM,*) 'NOMBRE de QUADRANGLES   =',NBELEM(4)
         WRITE(IMPRIM,*) 'NOMBRE de TETRAEDRES    =',NBELEM(5)
         WRITE(IMPRIM,*) 'NOMBRE de PYRAMIDES     =',NBELEM(9)
         WRITE(IMPRIM,*) 'NOMBRE de PENTAEDRES    =',NBELEM(6)
         WRITE(IMPRIM,*) 'NOMBRE d''HEXAEDRES      =',NBELEM(7)
         WRITE(IMPRIM,*) 'NOMBRE de 6-CUBES       =',NBELEM(8)
      ELSE
         WRITE(IMPRIM,*) 'NUMBER of VERTICES       =',NBSOM
         WRITE(IMPRIM,*) 'NUMBER of NODES          =',NBNOE
         WRITE(IMPRIM,*) 'NUMBER of VERTICES-POINTS=',NBELEM(1)
         WRITE(IMPRIM,*) 'NUMBER of SEGMENTS       =',NBELEM(2)
         WRITE(IMPRIM,*) 'NUMBER of TRIANGLES      =',NBELEM(3)
         WRITE(IMPRIM,*) 'NUMBER of QUADRANGLES    =',NBELEM(4)
         WRITE(IMPRIM,*) 'NUMBER of TETRAHEDRA     =',NBELEM(5)
         WRITE(IMPRIM,*) 'NUMBER of PYRAMIDS       =',NBELEM(9)
         WRITE(IMPRIM,*) 'NUMBER of PENTAHEDRA     =',NBELEM(6)
         WRITE(IMPRIM,*) 'NUMBER of HEXAHEDRA      =',NBELEM(7)
         WRITE(IMPRIM,*) 'NUMBER of 6-CUBES        =',NBELEM(8)
      ENDIF

C     FERMETURE DE L'OBJET
      CALL LXLXFE( NTOBJE, KNOMOB )
      IF( IERR .NE. 0 ) GOTO 9990
      GOTO 800

C     ==========================================================
C             STRUCTURATION DE L'OBJET EN SOUS-DOMAINES
C     ==========================================================
C     METHODE DES SOUS-DOMAINES
 500  CALL SOUDOM( KNOMOB, NTLXOB, MNDFOB, IERR )
      IF (IERR .NE. 0) THEN
         GOTO 9990
      ELSE
         GOTO 800
      ENDIF

C     ==========================================================
C                   METHODE DES JOINTS
C     ==========================================================
 510  CALL JOINTS( KNOMOB, NTLXOB, MNDFOB, IERR )
      IF (IERR .NE. 0) GOTO 9990

C     ==========================================================
C     AFFICHAGE DE FIN DU CALCUL DE L'INTERPOLATION DE CET OBJET
C     ==========================================================
 800  L = NUDCNB( KNOMOB )
      IF( LANGAG .EQ. 0 ) THEN
         KERR(2) = 'dftop0: FIN de CREATION du TMS: ~>OBJET>'
     %              // KNOMOB(1:L)
      ELSE
         KERR(2) = 'dftop0: CREATION END of TMS: ~>OBJET>'
     %              // KNOMOB(1:L)
      ENDIF
      L = NUDCNB( KERR(2) )
      KERR(1) = KERR(2)(1:L) // '>TOPOLOGIE'
      WRITE(IMPRIM,*) KERR(1)
      WRITE(IMPRIM,*)
      GOTO 9999

C     ERREUR DETECTEE  SUPPRESSION DU LEXIQUE DE L'OBJET
C     ===============
 9990 NBLGRC(NRERR) = 3
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'OBJET de NOM ' // KNOMOB
         KERR(2) = 'TOPOLOGIE IMPOSSIBLE A CONSTRUIRE. REDEFINIR'
         KERR(3) = 'COMPLETEMENT CET OBJET A PARTIR DE L''OPTION 5;'
      ELSE
         KERR(1) = 'OBJECT NAME ' // KNOMOB
         KERR(2) = 'IMPOSSIBLE TO CONSTRUCT THE TOPOLOGY'
         KERR(3) = 'DEFINE AGAIN THE OBJECT (OPTION 5;)'
      ENDIF
      CALL LEREUR
      CALL LXLXDS( NTOBJE, KNOMOB )
      IERR = 3

C     SORTIE
C     ======
 9999 IF( MNOBIN .GT. 0 ) CALL TNMCDS( 'ENTIER', NBOBPR*2, MNOBIN )
      IF( MNMNSO .GT. 0 ) CALL TNMCDS( 'ENTIER', NBOBPR,   MNMNSO )

      IF( IERALISTVOI .EQ. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'LISTVOI est DESALLOUE'
         ELSE
            PRINT*,'LISTVOI is DESALLOCATED'
         ENDIF
         DEALLOCATE( LISTVOI )
         IERALISTVOI = 1
      ENDIF

      IF( IERALPTDVOI .EQ. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'LPTDVOI est DESALLOUE'
         ELSE
            PRINT*,'LPTDVOI is DESALLOCATED'
         ENDIF
         DEALLOCATE( LPTDVOI )
         IERALPTDVOI = 1
      ENDIF

      RETURN
      END


