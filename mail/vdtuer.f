      SUBROUTINE VDTUER( KNMVOL, MXSOMM, PTXYZD, NBVOPA,
     %                   MXTETR, N1TEVI, N1TETS, NOTETR,
     %                   MXFACO, LEFACO,
     %                   MNF1VO, CHATVO, NBTETR, NUDTETR, NUDSOMM,IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DESTRUCTION DES TETRAEDRES EXTERIEURS AU DOMAINE
C -----    CHAINAGE DES TETRAEDRES DE CHACUN DES VOLUMES DU DOMAINE

C ENTREES:
C --------
C KNMVOL : NOM DU VOLUME A TETRAEDRISER
C MXSOMM : NOMBRE MAXIMUM DE POINTS DECLARABLES DANS PTXYZD
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C NBVOPA : NOMBRE DE VOLUMES DE LA PARTITION
C MXTETR : NOMBRE MAXIMUM DE TETRAEDRES DECLARABLES DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412

C MXFACO : NOMBRE MAXIMAL DE FACES DECLARABLES DANS LEFACO
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          1:   =0 POUR UNE FACE VIDE
C          123: NO (DANS PTXYZD) DU SOMMET 1 < SOMMET 2 < SOMMET 3
C          45:  NO (DANS NUVOPA 0 SINON) DU VOLUME1 , VOLUME2 DE LA FACE
C          678: NO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3

C          9: ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C             => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C             LEFACO(9,*) -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)

C          10: HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C              LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C              NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C              SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C              NF = LEFACO( 9, NF )  ...

C          11: >0  NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE,
C              =0  SINON

C SORTIES :
C ---------
C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NOTETR
C N1TETS : N1TETS(NS) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET NS
C NUDTETR: NU NOTETR DU DERNIER TETRAEDRE ACTIF
C NUDSOMM: NU PTXYZD DU DERNIER SOMMET ACTIF

C MNF1VO : F1VO(1,I) = NUMERO DU 1-ER TETRAEDRE DE CHAQUE VOLUME
C          F1VO(2,I) = NUMERO LOCAL DE LA FACE
C          F1VO(3,I) = NUMERO DU VOLUME DE 1 A NBVOPA
C          F1VO TABLEAU ENTIER( 1:3 , 1:NBVOPA )
C CHATVO : AU DEBUT D'EXECUTION DE VDTUER
C          >0  NUMERO DU TETRAEDRE SUIVANT DANS LE MEME VOLUME
C          =0  DERNIER TETRAEDRE DU VOLUME
C          =-1 TETRAEDRE JAMAIS CHAINE (EXTERIEUR) OU SUPPRIME DU CHAINAGE

C          A LA FIN DE L'EXECUTION DE VDTUER
C          >0  NUMERO DU VOLUME DE 1 A NBVOPA DE CHAQUE TETRAEDRE
C          =0  TETRAEDRE INACTIF

C NBTETR : NOMBRE DE TETRAEDRES ACTUELS DES NBVOPA VOLUMES
C IERR   : =0 SI PAS D'ERREUR
C          =5 SI AUCUNE FACE LEFACO RETROUVEE
C          >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        JUIN 1991
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY     MARS 2016
C MODIFS : ALAIN PERRONNET Saint PIERRE du PERRAY         Septembre 2018
C2345X7..............................................................012
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/langue.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0
C     STOPTE = FAUX ==> PAS D'ARRET APRES LE TRACE DE CHAQUE TETRAEDRE
C              VRAI ==> DEMANDE D'UN CARACTERE POUR REDEMARRER
C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C              VRAI ==> TRACE DES TETRAEDRES DE L'ETOILE

      INTEGER           LEFACO(11,0:MXFACO),N1TETS(*),NOTETR(8,*),
     %                  CHATVO(MXTETR)
      DOUBLE PRECISION  PTXYZD(4,*)
      CHARACTER*24      KNMVOL

ccc      WRITE(IMPRIM,10000)
ccc10000 FORMAT(/' ENTREE vdtuer: SUPPRESSION DES TETRAEDRES EXTERIEURS')
ccc      CALL TRFALEFACO( MXFACO, LEFACO, PTXYZD )
ccc      do 1 i=1,mxfaco
ccc         if( lefaco(1,i) .le. 0 ) GOTO 1
ccc         print *,'lefaco(',i,')=',(lefaco(k,i),k=1,11)
ccc         if( lefaco(4,i).le.0 .and. lefaco(5,i).le.0 ) then
ccc            print *,'PB avec lefaco(',i,')'
ccc         endif
ccc 1    continue

      PRINT *
      PRINT *,'vdtuer: Debut du RETRAIT des TETRAEDRES EXTERIEURS'

C     AU DEBUT F1VO(1,I) = NUMERO 1 FACE APPARTENANT AU SEUL VOLUME I
C     A LA FIN             NUMERO 1-ER TETRAEDRE DE CHAQUE VOLUME
      L = 3*NBVOPA
      CALL TNMCDC( 'ENTIER', L, MNF1VO )
      CALL AZEROI( L, MCN(MNF1VO) )

C     PIFA PILE DES FACES AFIN DE PARCOURIR LES TETRAEDRES D'UN VOLUME
      MXPIFA = MAX( 4096, MXTETR/10 )
      CALL TNMCDC( 'ENTIER', MXPIFA, MNPIFA )

C     RECHERCHE A PARTIR DU POINT 1 EXTERIEUR D'UNE FACE APPARTENANT
C     A UN VOLUME AVEC DESTRUCTION DES TETRAEDRES EXTERIEURS
C     ==============================================================
      NF = 2
C     2 ET NON 1 CAR SOMMET PLUS PROCHE DE L'HEXAEDRE ENGLOBANT

 10   NT1 = N1TETS( NF )
      IF( NT1 .LE. 0 ) THEN
         NF = NF + 1
         GOTO 10
      ENDIF

C     LES FACES DU TETRAEDRE EXTERIEUR NT1 SONT EMPILEES
C     LE NUMERO NT1 EST RETIRE DES FACES DE NT1
      NBPIFA = 0
      DO NF1=1,4

C        LA FACE APPARTIENT ELLE A UN VOLUME ?
         CALL NULEFT( NF1, NT1, NOTETR, MXFACO, LEFACO, NOFACO )
         IF( NOFACO .GT. 0 ) THEN
C           LA FACE APPARTIENT A UN VOLUME
            GOTO 100
         ENDIF

C        LA FACE N'APPARTIENT PAS A UN VOLUME. ELLE EST EMPILEE
         IF( NBPIFA+2 .GT. MXPIFA ) THEN
            WRITE(IMPRIM,*) 'vdtuer: AUGMENTER MXPIFA=',MXPIFA
            IERR = 5
            RETURN
         ENDIF

C        LA FACE NF1 DE NT1 EST EMPILEE
         MCN( MNPIFA + NBPIFA     ) = NT1
         MCN( MNPIFA + NBPIFA + 1 ) = NF1
         NBPIFA = NBPIFA + 2

      ENDDO

C     TANT QUE LA PILE EST NON VIDE
C        DEPILER LA FACE
C        NT2 L'AUTRE TETRAEDRE DE LA FACE
C        EMPILER SES FACES N'APPARTENANT PAS A UN VOLUME
C           SINON ALLER EN 100
C     ==================================================
 40   IF( NBPIFA .GT. 0 ) THEN

C        LA FACE EST DEPILEE
         NBPIFA = NBPIFA - 2

C        LE TETRAEDRE
         MN  = MNPIFA + NBPIFA
         NT1 = MCN( MN )

C        LE NUMERO LOCAL DE LA FACE DANS LE TETRAEDRE
         NF1 = MCN( MN + 1 )

C        LE TETRAEDRE OPPOSE A NT1 PAR CETTE FACE
         NT2 = NOTETR( 4+NF1, NT1 )

C        LA FACE NF1 DE NT1 EST MARQUEE TRAITEE
         NOTETR( 4+NF1, NT1 ) = -ABS( NT2 )

         IF( NT2 .LE. 0 ) GOTO 40

C        BOUCLE SUR LES 4 FACES DU TETRAEDRE NT2
         DO 50 NF2=1,4

C           LA FACE NF2 N'APPARTIENT PAS A UN VOLUME
            NTOP = NOTETR( 4+NF2, NT2 )
            IF( NTOP .LE. 0 ) THEN
C               LA FACE EST DEJA TRAITEE DANS LE TETRAEDRE OPPOSE
               GOTO 50
            ENDIF

C           RECHERCHE DE LA FACE NF2 DANS NT2 QUI EST LA FACE NF1 DE NT1
            IF( NTOP .EQ. NT1 ) THEN
C              LA FACE NF2 DE NT2 EST MARQUEE TRAITEE
               NOTETR( 4+NF2, NT2 ) = -NT1
               GOTO 50
            ENDIF

C           LA FACE NF2 APPARTIENT ELLE A UN VOLUME ?
            CALL NULEFT( NF2, NT2, NOTETR, MXFACO, LEFACO, NOFACO )
            IF( NOFACO .GT. 0 ) THEN
C              LA FACE APPARTIENT A UN VOLUME
               NT1 = NT2
               NF1 = NF2
               GOTO 100
            ENDIF

C           LA FACE NF2 N'APPARTIENT PAS A UN VOLUME
C           LE TETRAEDRE OPPOSE N'EST PAS TRAITE => EMPILE
            IF( NBPIFA+2 .GT. MXPIFA ) THEN

C              LA PILE TROP COURTE EST AUGMENTEE
               IF( MXPIFA .GE. MXTETR ) THEN
C                 PROBLEME: PAS DE FACE LEFACO RETROUVEE
C                 IL Y A UN TROU DANS LA FRONTIERE DE L'OBJET
                  WRITE(IMPRIM,*)
     %            'vdtuer: PB PAS DE FACE LEFACO RETROUVEE'
                  WRITE(IMPRIM,*) 'ESSAYER UNE AUTRE VALEUR DE ARETGR'
                  TRACTE0 = TRACTE
                  TRACTE = .TRUE.
                  CALL TRLEFACO( KNMVOL, MXFACO, LEFACO, PTXYZD )
                  TRACTE = TRACTE0
                  IERR = 5
                  RETURN
               ENDIF
               CALL TNMCAU( 'ENTIER', MXPIFA, MXPIFA*2,
     %                       NBPIFA , MNPIFA )
               MXPIFA = MXPIFA*2

            ENDIF

C           LE TETRAEDRE ET NUMERO LOCAL DE FACE SONT EMPILES
            MCN( MNPIFA + NBPIFA     ) = NT2
            MCN( MNPIFA + NBPIFA + 1 ) = NF2
            NBPIFA = NBPIFA + 2

 50      ENDDO
         GOTO 40
C        FIN TANT QUE
      ENDIF

      WRITE(IMPRIM,*) 'vdtuer ERREUR: 0 FACE APPARTIENT A 1 VOLUME'
      IERR = 1
      RETURN

C     ----------------------------------------------------------------
C     LE TETRAEDRE NT1 EST EXTERNE . SA FACE NF1 EST EN CONTACT
C     AVEC UN MATERIAU. LE TETRAEDRE OPPOSE SERT DE DEPART AU PARCOURS
C     DES TETRAEDRES INTERIEURS AUX VOLUMES
C     ----------------------------------------------------------------
 100  NOVOLU = LEFACO(4,NOFACO)
      IF( NOVOLU .LE. 0 ) THEN
         NOVOLU = LEFACO(5,NOFACO)
      ENDIF
C
C     LE 1-ER VOLUME TRAITE EST NOVOLU
      NBVOPT = 1
C     LE TETRAEDRE INTERNE OPPOSE A NT1
      NT2 = NOTETR( 4+NF1, NT1 )
C     LE NUMERO DE FACE DANS NT2
      DO NF2=1,4
         IF( NOTETR(4+NF2,NT2) .EQ. NT1 ) GOTO 108
      ENDDO
C     LE TETRAEDRE NT2 EST STOCKE
 108  MCN( MNF1VO     ) = NT2
C     LE NUMERO LOCAL DE LA FACE
      MCN( MNF1VO + 1 ) = NF2
C     SUIVIE DU NUMERO DE SON VOLUME
      MCN( MNF1VO + 2 ) = NOVOLU
C
C     LE NOMBRE TOTAL DE TETRAEDRES
      NBTETR = 0
C
C     BOUCLE SUR LES VOLUMES DE LA PARTITION DANS UN ORDRE QUELCONQUE
C     ===============================================================
      DO 400 NOVOL = 1 , NBVOPA

C        RECHERCHE D'UN VOLUME NON TRAITE
         MN = MNF1VO
         DO I=1,NBVOPA
            IF( MCN( MN ) .GT. 0 ) GOTO 120
            MN = MN + 3
         ENDDO
         WRITE(IMPRIM,*) 'ERREUR vdtuer: VOLUME',NOVOL,' SANS FACE'
         IERR = 1
         RETURN
C
C        RECHERCHE DE LA FACE DU 1-ER TETRAEDRE A L'INTERIEUR DU VOLUME
 120     NT1    = MCN( MN )
         NF1    = MCN( MN + 1 )
C        SON VOLUME
         NOVOLU = MCN( MN + 2 )
C
C        LE CHAINAGE SUR LE 1-ER TETRAEDRE EST EFFECTUE DANS F1VO
C        LE SIGNE - POUR MONTRER LE TRAITEMENT DE CE VOLUME
         MCN( MN ) = -NT1
C
C        LA PILE DES FACES DU VOLUME EST INITIALISEE
         NBPIFA = 2
         MCN( MNPIFA     ) = NT1
         MCN( MNPIFA + 1 ) = NF1
C
C        LE TETRAEDRE PRECEDENT NT1 N'EXISTE PAS
         NT0 = 0
C
C        TANT QUE LA PILE DES FACES EST NON VIDE
C        =======================================
C           DEPILER LA FACE
C           CHERCHER LE TETRAEDRE NT1 AU DELA
C           CHAINER LE TETRAEDRE NT1
C           EMPILER SES FACES NON FRONTALIERES DU VOLUME
C                    ET SUPPRIMER NT1 POUR CETTE FACE
C           OU AJOUTER LA FACE A F1VO POUR LE VOLUME OPPOSE
C
 200     IF( NBPIFA .GT. 0 ) THEN
            NBPIFA = NBPIFA - 2
C           LE TETRAEDRE ET SA FACE EST DEPILE
            NT1 = MCN( MNPIFA + NBPIFA )
            IF( NT1 .LE. 0 ) GOTO 200
            DO NF1 = 1,4
C              LE TETRAEDRE EST IL TRAITE ?
               IF( NOTETR(4+NF1,NT1) .LE. 0 ) GOTO 200
C              NON
            ENDDO
C
C           LA FACE ET LE TETRAEDRE SONT MARQUES
            NF1 = MCN( MNPIFA + NBPIFA + 1 )
            NOTETR(4+NF1,NT1) = -ABS( NOTETR(4+NF1,NT1) )

            DO 260 NF1=1,4
C
C              LE TETRAEDRE OPPOSE A LA FACE NF1 DE NT1
               NT2 = NOTETR(4+NF1,NT1)
               IF( NT2 .GT. 0 ) THEN
C                 FACE NON TRAITEE. LA FACE EST MARQUEE
                  NOTETR(4+NF1,NT1) = -ABS( NT2 )
C
C                 LA FACE NF1 APPARTIENT ELLE A LEFACO?
                  CALL NULEFT( NF1, NT1, NOTETR, MXFACO, LEFACO, NOFACO)
                  IF( NOFACO .EQ. 0 ) THEN
C
C                    LA FACE N'EST PAS FRONTALIERE. ELLE EST EMPILEE
C                    -----------------------------
                     IF( NBPIFA+2 .GT. MXPIFA ) THEN
C                       LA PILE TROP COURTE EST AUGMENTEE
                        CALL TNMCAU( 'ENTIER', MXPIFA, MXPIFA+1024,
     %                                NBPIFA , MNPIFA )
                        MXPIFA = MXPIFA + 1024
                     ENDIF
C                    RECHERCHE DE LA FACE DANS NT2
                     DO NF2=1,4
                        IF( ABS(NOTETR(4+NF2,NT2)) .EQ. NT1 ) GOTO 220
                     ENDDO
                     print *,'vdtuer: Anomalie boucle 210 avec'
                     print *,'TETRAEDRE NT1',NT1,(NOTETR(K,NT1),K=1,8)
                     print *,'TETRAEDRE NT2',NT2,(NOTETR(K,NT2),K=1,8)

 220                 IF( NOTETR(4+NF2,NT2) .GT. 0 ) THEN
C                       FACE NON TRAITEE => EMPILEE
                        MCN( MNPIFA + NBPIFA     ) = NT2
                        MCN( MNPIFA + NBPIFA + 1 ) = NF2
                        NBPIFA = NBPIFA + 2
                     ENDIF
C
                  ELSE
C
C                    LA FACE EST FRONTALIERE ACTIVE
C                    ------------------------------
C                    EXISTE-T-IL UN AUTRE VOLUME DIFFERENT DE NOVOLU
                     NV1 = LEFACO(4,NOFACO)
                     NV2 = LEFACO(5,NOFACO)
                     IF( NV1 .NE. NOVOLU .AND. NV2 .NE. NOVOLU ) THEN
C                       ERREUR NOVOLU N'EST PAS UN VOLUME DE CETTE FACE
                        WRITE(IMPRIM,10230) NOVOLU,NV1,NV2,NF2,NT2
10230 FORMAT(' ERREUR vdtuer: VOLUME',I5,' NON EGAL A',2I5/
     %' FACE',I1,' TETRAEDRE',I8/)
                     ENDIF
C
C                    LE VOLUME NOVOL2 AU DELA DE CELUI TRAITE
                     NOVOL2 = NV1
                     IF( NOVOL2 .EQ. NOVOLU ) NOVOL2 = NV2
                     IF( NOVOL2 .EQ. 0 ) GOTO 260
C
C                    CE VOLUME NOVOL2 A T IL DEJA ETE RECENSE ?
                     MN1 = MNF1VO + 2
                     DO J=1,NBVOPT
                        IF( MCN(MN1) .EQ. NOVOL2 ) GOTO 260
                        MN1 = MN1 + 3
                     ENDDO
C
C                    LE VOLUME EST AJOUTE
                     MCN( MN1     ) = NOVOL2
                     MCN( MN1 - 1 ) = NF2
                     MCN( MN1 - 2 ) = NT2
                     NBVOPT = NBVOPT + 1
                  ENDIF
               ENDIF
 260        ENDDO
C
C           LE TETRAEDRE NT1 EST CHAINE A CEUX DU VOLUME NOVOLU
C           COMME ETANT LE SUIVANT DE NT0 TETRAEDRE PRECEDENT DE NT1
            IF( NT0 .GT. 0 ) THEN
                CHATVO( NT0 ) = NT1
            ENDIF
            NT0 = NT1

            GOTO 200
C           FIN TANT QUE SUR LES FACES DE LA PILE
         ENDIF
C
C        FIN DU CHAINAGE DES TETRAEDRES DU VOLUME NOVOLU
         CHATVO( NT0 ) = 0

 400  ENDDO

C     LE CHAINAGE SUR LES TETRAEDRES REPREND UN SIGNE POSITIF
C     -------------------------------------------------------
      MN = MNF1VO
      DO NOVOL=1,NBVOPA
         IF( MCN( MN ) .LT. 0 ) THEN
            MCN(MN) = -MCN(MN)
         ELSE
            WRITE(IMPRIM,*) 'ERREUR vdtuer:VOLUME',NOVOL,' NON CHAINE'
         ENDIF
         MN = MN + 3
      ENDDO

C     DANS NOTETR, LE NUMERO DE FACE OPPOSEE NEGATIF REDEVIENT POSITIF
C     REMISE A ZERO DES TETRAEDRES OPPOSES EXTERIEURS
C     COMPTAGE DES TETRAEDRES DES VOLUMES ET DE LA PARTITION
C     ----------------------------------------------------------------
      NBTETR = 0
      MN = MNF1VO - 3
      DO NOVOL=1,NBVOPA

C        LE 1-ER TETRAEDRE DU VOLUME
         NBTETV = 0
         MN  = MN + 3
         NT0 = MCN(MN)
C
 415     IF( NT0 .GT. 0 ) THEN
            NBTETV = NBTETV + 1
            DO I=5,8
               NTOP = ABS( NOTETR(I,NT0) )
               IF( NTOP .GT. 0 .AND. CHATVO(NTOP) .LT. 0 ) THEN
C                 NTOP EST UN TETRAEDRE EXTERIEUR AUX NBVOPA VOLUMES
                  NOTETR(I,NT0) = 0
               ELSE
                  NOTETR(I,NT0) = NTOP
               ENDIF
            ENDDO
C           LE TETRAEDRE SUIVANT
            NT0 = CHATVO( NT0 )
            GOTO 415
         ENDIF

         IF( LANGAG .EQ. 0 ) THEN
           WRITE(IMPRIM,*)'vdtuer:',NBTETV,' TETRAEDRES dans le Volume',
     %                              NOVOLU
         ELSE
            WRITE(IMPRIM,*)'vdtuer:',NBTETV,' TETRAEDRA in Volume',
     %                               NOVOLU
         ENDIF

C        NOMBRE TOTAL DE TETRAEDRES
         NBTETR = NBTETR + NBTETV

      ENDDO

      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)'vdtuer:',NBTETR,' TETRAEDRES TOTAL de TOUS les 
     %VOLUMES'
      ELSE
         WRITE(IMPRIM,*)'vdtuer:',NBTETR,' TETRAHEDRA in ALL VOLUMES'
      ENDIF

C     REMISE A ZERO DES TETRAEDRES EXTERIEURS ET N1TETS
C     -------------------------------------------------
      DO N = 1, MXSOMM
         N1TETS( N ) = 0
      ENDDO
      NUDSOMM = 0
      NUDTETR = 0
      N1TEVI  = 0
      DO N = MXTETR, 1, -1
         IF( CHATVO( N ) .LT. 0 ) THEN
C           N EST UN TETRAEDRE EXTERIEUR -> VIDE
            NOTETR( 1, N ) = 0
C           MISE A JOUR DU 1-ER TETRAEDRE VIDE
            NOTETR( 5, N ) = N1TEVI
            N1TEVI = N
         ELSE
            IF( NUDTETR .EQ. 0 ) NUDTETR = N
            DO K=1,4
               NS = NOTETR( K, N )
               N1TETS( NS ) = N
               NUDSOMM = MAX( NUDSOMM, NS )
            ENDDO
         ENDIF
      ENDDO


C     CHATVO CHAINAGE SUR LE TETRAEDRE SUIVANT DE CHACUN DES NBVOPA VOLUMES
C     DEVIENT LE NUMERO DE VOLUME DE 1 A NBVOPA DE CHAQUE TETRAEDRE ACTIF
C     ---------------------------------------------------------------------
      DO NOVOL = 1, NBVOPA

C        RECHERCHE D'UN VOLUME NON TRAITE
         MN = MNF1VO
         DO I = 1, NBVOPA
            IF( MCN( MN ) .GT. 0 ) GOTO 600
            MN = MN + 3
         ENDDO
         WRITE(IMPRIM,*) 'ERREUR vdtuer: VOLUME',NOVOL,' SANS TETRAEDRE'
         IERR = 1
         RETURN
C
C        RECHERCHE DU 1-ER TETRAEDRE A L'INTERIEUR DU VOLUME
 600     NT0    = MCN( MN )
         NF0    = MCN( MN + 1 )
C        SON VOLUME
         NOVOLU = MCN( MN + 2 )

C        MARQUE DE TRAITEMENT DE CE VOLUME
         MCN( MN ) = -NT0

 610     IF( NT0 .GT. 0 ) THEN
            NT1 = CHATVO( NT0 )
            CHATVO( NT0 ) = NOVOLU
            NT0 = NT1
            GOTO 610
         ENDIF

      ENDDO

C     NUMERO DE VOLUME DES TETRAEDRES INACTIFS
      DO NT0 = 1, MXTETR
         IF( NOTETR(1,NT0) .EQ. 0 ) THEN
C           LE VOLUME D'UN TETRAEDRE INACTIF EST -1
            CHATVO( NT0 ) = -1
         ENDIF
      ENDDO

C     LE CHAINAGE SUR LES TETRAEDRES REPREND UN SIGNE POSITIF
C     -------------------------------------------------------
      MN = MNF1VO
      DO NOVOL = 1, NBVOPA
         IF( MCN( MN ) .LT. 0 ) THEN
            MCN(MN) = -MCN(MN)
         ELSE
            WRITE(IMPRIM,*) 'vdtuer: ERREUR VOLUME',NOVOL,' NON CHAINE'
         ENDIF
         MN = MN + 3
      ENDDO

C     DESTRUCTION DU TABLEAU INUTILE ENSUITE
      CALL TNMCDS( 'ENTIER', MXPIFA, MNPIFA )

      RETURN
      END
