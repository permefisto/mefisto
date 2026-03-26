      SUBROUTINE AMQUTETRA( QUTEMN, ANTE2P, NBITER,
     %                      NBCOOR, NBSOMM, MXSOMM, XYZSOM, NPSOFR,
     %                      NBSOTE, MXTETR, NBTETR, NSTETR,
     %                      NBDM,   NUDMEF,
     %                      VOLUMT, QUALIT, QUALII,
     %                      NO1TSO, MXTESO, NOTESO,
     %                      MXFAET, NFETOI, VOETOI, QUETOI,
     %                      NOSOET, MXTEXA, NOTEXA,
     %                      DCPUTO, IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AMELIORER LA QUALITE D'UNE TETRAEDRISATION
C -----
C ENTREES:
C --------
C QUTEMN : QUALITE MINIMALE AU DESSOUS DE LAQUELLE UN TETRAEDRE DOIT ETRE
C          AMELIORE
C ANTE2P : ANGLE AU DESSOUS DUQUEL 2 FACES ADJACENTES SONT CONSIDEREES
C          COPLANAIRES ET PERMISSION D'ECHANGER LES DIAGONALES
C NBITER : NOMBRE MAXIMAL D'ITERATIONS SOUHAITEES
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS DECLARABLES DANS XYZSOM ET NPSOFR
C XYZSOM : COORDONNEES X Y Z DES NBSOMM SOMMETS DES TETRAEDRES
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NSTETR(>3)
C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES DECLARABLES DANS NSTETR
C NBDM   : 0 SI 1 MATERIAU=VOLUME, SINON NOMBRE DE MATERIAUX DU VOLUME
C MXTESO : NOMBRE MAXIMAL DE NUMERO DE TETRAEDRES DES SOMMETS
C MXFAET : NOMBRE MAXIMAL DE FACES DECLARABLES DANS NFETOI

C TABLEAUX AUXILIAIRES :
C ----------------------
C NFETOI(5,MXFAET)  DES ENTIERS
C NOSOET(MXTETR)    DES ENTIERS
C NOTEXA(MXTEXA)    DES ENTIERS

C MODIFIES :
C ----------
C NBSOMM : NOMBRE ACTUEL DE SOMMETS DE LA TETRAEDRISATION
C NBTETR : NOMBRE DE TETRAEDRES
C NPSOFR : NUMERO 0 SI SOMMET INTERNE
C                 1 SI SOMMET SUR LA FRONTIERE
C                 2 SI SOMMET SUR L'INTERFACE ENTRE 2 MATERIAUX
C                -1 SI SOMMET SUPPRIME LORS DE L'AMELIORATION
C NSTETR : NUMERO DES NBSOTE SOMMETS DE CHAQUE TETRAEDRE DU MAILLAGE
C NUDMEF : NUMERO DE MATERIAU DE CHAQUE TETRAEDRE DU MAILLAGE
C          ATTENTION: CE TABLEAU EXISTE SEULEMENT SI NBDM>0
C DCPUTO : TEMPS CPU TOTAL POUR L'AMELIORATION

C SORTIES:
C --------
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DANS NOTESO DE CHAQUE SOMMET
C NOTESO : NOTESO(1,*) NUMERO DANS NSTETR DU TETRAEDRE
C          NOTESO(2,*) NUMERO DANS NOTESO DU TETRAEDRE SUIVANT
C                      0 SI C'EST LE DERNIER
C QUALIT : QUALITE DES TETRAEDRES DE LA TETRAEDRISATION
C QUALII : PROTECTION DE LA QUALITE DES TETRAEDRES AVANT TRI
C VOLUMT : VOLUME  DES TETRAEDRES DE LA TETRAEDRISATION
C IERR   : =0 SI PAS D'ERREUR
C          >0 EN CAS DE SATURATION D'UN TABLEAU
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   Decembre 1991
C MODIFS : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   Novembre 1993
C MODIFS : ALAIN PERRONNET LJLL UPMC et St Pierre du Perray    Juin 2008
C MODIFS : ALAIN PERRONNET LJLL UPMC et St Pierre du Perray    Aout 2012
C MODIFS : ALAIN PERRONNET LJLL UPMC et St Pierre du Perray Fevrier 2016
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))

      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE,TRACTE0

      INTEGER           NPSOFR(MXSOMM),
     %                  NSTETR(NBSOTE,MXTETR),
     %                  NUDMEF(MXTETR),
     %                  NO1TSO(MXSOMM),
     %                  NOTESO(2,MXTESO) ,
     %                  NFETOI(5,MXFAET),
     %                  NOSOET(MXTETR),
     %                  NOTEXA(MXTEXA)

      REAL              XYZSOM(3,MXSOMM),
     %                  VOLUMT(MXTETR),
     %                  QUALIT(MXTETR),
     %                  QUALII(MXTETR),
     %                  VOETOI(MXFAET),
     %                  QUETOI(MXFAET)

      DOUBLE PRECISION  DINFO, DCPU, DCPUTO
      INTEGER, allocatable, dimension(:,:) :: LFACES

      IERR = 0
      IERALFACES = 1
      TRACTE = .FALSE.
      TRACTE0 = TRACTE

C     A PARTIR DES NBTETR TETRAEDRES FORMATION DE LA LISTE NOTESO
C     DES TETRAEDRES DE CHAQUE SOMMET AVEC SON POINTEUR NO1TSO
C     ===========================================================
      CALL TETSOM( MXSOMM, NBSOTE, MXTETR, NBTETR, NSTETR,
     %             N1TEVI, N1TESO, NO1TSO, MXTESO, NOTESO )

C     DECLARATION DU TABLEAU LFACES DES FACES DES TETRAEDRES
C     POUR DETECTER LES FACES APPARTENANT A 3 TETRAEDRES.
C     A SUPPRIMER ENSUITE
C     ======================================================
      MXFACE = 4 * NBTETR
      MOFACE = 7
      L = MOFACE * MXFACE
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'DEMANDE ALLOCATION de',L,' ENTIERS pour LFACES'
         ALLOCATE ( LFACES( 1:MOFACE, 1:MXFACE ), STAT=IERALFACES )
         IF( IERALFACES .NE. 0 ) THEN
            PRINT*,'ERREUR ALLOCATION de',L,' ENTIERS pour LFACES'
            IERR = IERALFACES
            GOTO 9900
         ENDIF
         PRINT*,'CORRECTE ALLOCATION de LFACES(1:',MOFACE,
     %          ',1:',MXFACE,') ENTIERS'
      ELSE
         PRINT*,'ALLOCATION DEMAND of',L,' INTEGER of LFACES'
         ALLOCATE ( LFACES( 1:MOFACE, 1:MXFACE ), STAT=IERALFACES )
         IF( IERALFACES .NE. 0 ) THEN
            PRINT*,'ALLOCATION ERROR  of',L,' INTEGER of LFACES'
            IERR = IERALFACES
            GOTO 9900
         ENDIF
         PRINT*,'CORRECT ALLOCATION of LFACES(1:',MOFACE,
     %          ',1:',MXFACE,') INTEGER'
      ENDIF

C     DETECTION DES FACES APPARTENANT A 3 TETRAEDRES
      CALL CRHAFAVO( NBSOTE, MXTETR, NSTETR, NUDTETR,
     %               MOFACE, MXFACE, LFACES, NBFAPB, IERR )
      IF( NBFAPB .GT. 0 ) THEN
         PRINT *
         PRINT *,'amqutetra 0):',NBFAPB,' FACES DE 3 TETRAEDRES'
         TRACTE = .TRUE.
      ENDIF

C     TRACE FILAIRE DU MAILLAGE, ARETES ROUGES SI TETRA DE QUALITE<0.15
C     VOLUME ET QUALITE DU MAILLAGE INITIAL
      CALL TRNSTETR( NBSOTE, MXTETR, NSTETR, NBCOOR, XYZSOM,
     %               NBTETR, VOLET00, QUAET00, QUAMY00 )
      TRACTE = TRACTE0

      PRINT *
      PRINT *,'amqutetra step 0 : Debut ITERATION=  0 : Maillage VOLET00
     %=',VOLET00,' QUAET00=',QUAET00,' QUAMY00=',QUAMY00,
     % ' NBFAPB=',NBFAPB,
     %' ==============================================================='

C     INITIALISATION A VIDE DE TOUTES LES FACES D'UNE ETOILE
C     ------------------------------------------------------
      CALL ETOIVIDE( MXFAET, N1FEOC, N1FEVI, NFETOI )

C     ==================================================================
C     LES ITERATIONS D'AMELIORATION
C     ==================================================================
      GRAND = RINFO( 'GRAND' )

C     NUMERO MAXIMUM DES TETRAEDRES ACTUELS
      MXTETA = NBTETR
11000 FORMAT(' amqutetra: ITERATION',I3,' AMELIORATION de la QUALITE de
     % la TETRAEDRISATION:'/150('-'))
21000 FORMAT(' amqutetra: STEP',I3,' TETRAHEDRIZATION QUALITY IMPROVEMEN
     %T'/150('-'))

      DO 1000 ITER=1,NBITER

C        ***************************************************************
C        1-ERE PARTIE : RESORBER LES TETRAEDRES DE TRES MAUVAISE QUALITE
C                       BOUCLE SUR LES TETRAEDRES DE MAUVAISE QUALITE
C        ***************************************************************

C        INITIALISATION DE LA QUALITE ET VOLUME DES TETRAEDRES
         IF( LANGAG .EQ. 0 ) THEN
            print 11000, ITER
         ELSE
            print 21000, ITER
         ENDIF

         CALL QUMESH(' amqutetra: Qualite TETRAEDRES', ITER, GRAND,
     %                XYZSOM, NBSOTE, MXTETR, NSTETR,
     %                VOLUMT, QUALIT, QUTEMN,
     %                NBTU,   VOLT,   QUAMY0, QUAET0,
     %                NBQINF, QUALII, NOSOET )

C        VOLUME MOYEN D'UN TETRAEDRE
         VOLMOYTE = VOLT / NBTU

C        TRI CROISSANT SELON LA QUALITE DES NBQINF MAUVAIS TETRAEDRES
         CALL TRITRP( NBQINF, QUALII, NOSOET )

C        SAUVEGARDE DU NUMERO DES PLUS MAUVAIS TETRAEDRES
C        DANS UN PLUS PETIT TABLEAU NOTEXA
         MXTET = MIN( MXTEXA, NBQINF )

ccc            print *,'amqutetra step 1 : MXTEXA=',mxtexa,
ccc     %              ' nbqinf=',NBQINF,' mxtet=',mxtet
ccc         IF( LANGAG .EQ. 0 ) THEN
ccc         print *,'amqutetra step 1 : NOMBRE MAUVAIS TETRAEDRES TRAITES='
ccc     %           ,MXTET
ccc         ELSE
ccc          print *,'amqutetra step 1 : NUMBER of BAD TREATED TETRAHEDRA='
ccc     %           ,MXTET
ccc         ENDIF

         CALL TRTATA( NOSOET, NOTEXA, MXTET )

C        BOUCLE SUR LES MXTET PLUS MAUVAIS TETRAEDRES
C        --------------------------------------------
         DO NTQF=1,MXTET

C           LE NUMERO NSTETR DU TETRAEDRE NTQF DE NOTEXA
            NT = NOTEXA( NTQF )
            IF( NSTETR(1,NT) .GT. 0 ) THEN

C              AMELIORATION DE LA QUALITE DE CE TETRAEDRE NT DE NSTETR
               IF( NTQF .LE. 3 ) THEN
                 print*,'amqutetra: PLUS MAUVAIS Tetraedre',NTQF,':',NT,
     %                  ' Volume=',VOLUMT(NT),' Qualit=',QUALIT(NT)
                  DO kk=1,4
                     NS = NSTETR(KK,NT)
                     print*,'amqutetra: St',NS,' X=',XYZSOM(1,NS),
     %                      ' Y=',XYZSOM(2,NS),' Z=',XYZSOM(3,NS)
                  ENDDO
               ENDIF

               CALL QUA1TE( NT,  QUTEMN, GRAND,  VOLMOYTE,
     %                   MOFACE, MXFACE, LFACES, NBFAPB,
     %                   NBSOMM, MXSOMM, XYZSOM, NPSOFR, NBDM,   NUDMEF,
     %                   NBSOTE, MXTETR, MXTETA, N1TEVI, NSTETR,
     %                   VOLUMT, QUALIT, NO1TSO, MXTESO, N1TESO, NOTESO,
     %                   MXFAET, N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                   NOSOET, VOLET0, QUAET0, VOLET1, QUAET1, IERR  )

ccc               IF( NBFAPB .GT. 0 ) THEN

ccc                  PRINT *
ccc                  PRINT *,'amqutetra: 1 NTE=',NT,
ccc     %                     NBFAPB,' FACES DE 3 TETRAEDRES'
ccc                  TRACTE = .TRUE.

cccC                 TRACE FILAIRE DU MAILLAGE, ARETES ROUGES SI TETRA DE QUALITE<0.15
cccC                 VOLUME ET QUALITE DU MAILLAGE FINAL
ccc                  CALL TRNSTETR( NBSOTE, MXTETR, NSTETR, NBCOOR, XYZSOM,
ccc     %                           NBTETR, VOLET1, QUAET1, QUAMY1 )

ccc                  PRINT*,'amqutetra step 1 : Tetra',NT,
ccc     %                   ' VOLET0=',VOLET0,' QUAET0=',QUAET0
ccc                  PRINT*,'amqutetra step 1 : Tetra',NT,
ccc     %                   ' VOLET1=',VOLET1,' QUAET1=',QUAET1

ccc               ENDIF

            ENDIF

         ENDDO

C        DETECTION FINALE DES FACES APPARTENANT A 3 TETRAEDRES
         CALL CRHAFAVO( NBSOTE, MXTETR, NSTETR, NUDTETR,
     %                  MOFACE, MXFACE, LFACES, NBFAPB, IERR )
         IF( NBFAPB .GT. 0 ) THEN
            PRINT *,'amqutetra: 1',NBFAPB,' FACES DE 3 TETRAEDRES'
            TRACTE = .TRUE.
         ENDIF

C        TRACE FILAIRE DU MAILLAGE, ARETES ROUGES SI TETRA DE QUALITE<0.15
C        VOLUME ET QUALITE DU MAILLAGE FINAL
         CALL TRNSTETR( NBSOTE, MXTETR, NSTETR, NBCOOR, XYZSOM,
     %                  NBTETR, VOLET1, QUAET1, QUAMY1 )
         TRACTE = TRACTE0

         DCPU   = DINFO('DELTA CPU')
         DCPUTO = DCPUTO + DCPU

         print*
         IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'amqutetra step 1 : Qualite des MAUVAIS TETRAEDRES -----
     %-----------------------------------------------------------------'
         ELSE
         PRINT*,'amqutetra step 1 : TETRAHEDRA BAD QUALITY  ------------
     %-----------------------------------------------------------------'
         ENDIF

         PRINT*,'amqutetra step 1 ITER=',ITER,
     %                     ' VOLET00=',VOLET00,' QUAET00=',QUAET00,
     %                     ' QUAMY00=',QUAMY00
         PRINT*,'amqutetra step 1 ITER=',ITER,
     %                     ' VOLET 1=',VOLET1,' QUAET 1=',QUAET1,
     %                     ' QUAMY 1=',QUAMY1,
     %          ' VOL Diff=',ABS(VOLET1-VOLET00)/VOLET00*100,' %'

         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*)'FIN BOUCLE sur MAUVAIS TETRAEDRES',
     %                      NINT(100D0*DCPU)/100.0,' SECONDES'
         ELSE
            WRITE(IMPRIM,*)'LOOP END on BAD TETRAHEDRA',
     %                      NINT(100D0*DCPU)/100.0,' SECONDS'
         ENDIF


C        *******************************************************************
C        2-EME PARTIE : BOUCLE SUR LES SOMMETS
C                       DEPLACEMENT OU SUPPRESSION DE SOMMETS
C        *******************************************************************
C        INITIALISATION DE LA QUALITE ET VOLUME DES TETRAEDRES
         CALL QUMESH(' amqutetra: Qualite SOMMETS', ITER, GRAND, XYZSOM,
     %                NBSOTE, MXTETR, NSTETR,
     %                VOLUMT, QUALIT, QUTEMN,
     %                NBTU,   VOLT,   QUAMY1, QUAET1,
     %                NBQINF, QUALII, NOSOET )
C
C        GENERATION DU TABLEAU DES MAUVAIS SOMMETS
C        =========================================
         NBQINF = 0
         DO 150 NT=1,MXTETR
            DO 130 J=1,4
C              LE SOMMET J DU TETRAEDRE NT
               NS = NSTETR(J,NT)
               IF( NS .EQ. 0 ) GOTO 150
               IF( NPSOFR( NS ) .EQ. -1 ) GOTO 130
               IF( NO1TSO( NS ) .EQ.  0 ) GOTO 150
               IF( NO1TSO( NS ) .LT.  0 ) THEN
C                 NS EST DEJA STOCKE DANS NOSOET OU A DEJA ETE ANALYSE
                  GOTO 130
               ENDIF
C              NS N'EST PAS ENCORE STOCKE DANS NOSOET. IL EST MARQUE
               NO1TSO( NS ) = -NO1TSO( NS )

C              CALCUL DU VOLUME ET QUALITE DES TETRAEDRES DE SOMMET NS
               QUALI0 = GRAND
C              POSITION DANS NOTESO DU 1-ER TETRAEDRE DE SOMMET NS
               NDT = ABS( NO1TSO( NS ) )
C              TANT QU'IL EXISTE UN TETRAEDRE DE SOMMET NS FAIRE
 120           IF( NDT .GT. 0 ) THEN
C                 LE NUMERO DU TETRAEDRE DANS NSTETR
                  NTQF    = NOTESO(1,NDT)
                  QUALI0 = MIN( QUALI0 , QUALIT(NTQF) )
C                 LE TETRAEDRE SUIVANT
                  NDT = NOTESO(2,NDT)
                  GOTO 120
               ENDIF

C              STOCKAGE DU SOMMET S'IL EST INFERIEUR A LA QUALITE SOUHAITEE
               IF( QUALI0 .LT. QUTEMN ) THEN
                  NBQINF = NBQINF + 1
                  NOSOET( NBQINF ) = NS
C                 LA QUALITE DU SOMMET
                  QUALII( NBQINF ) = QUALI0
               ENDIF
 130        CONTINUE
 150     CONTINUE

C        TRI CROISSANT SELON LA MAUVAISE QUALITE DES SOMMETS
         CALL TRITRP( NBQINF, QUALII, NOSOET )

C        AFFICHAGE DE LA QUALITE MINIMALE DES SOMMETS POUR CETTE ITERATION
C        ET DES COORDONNEES DU PLUS MAUVAIS SOMMET
         NS = NOSOET(1)
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,11150) ITER, NBQINF, QUALII(1),
     %                          NS, (XYZSOM(J,NS),J=1,3)
         ELSE
            WRITE(IMPRIM,21150) ITER, NBQINF, QUALII(1),
     %                          NS, (XYZSOM(J,NS),J=1,3)
         ENDIF
11150 FORMAT( ' ITERATION',I3,
     %        ' NOMBRE de MAUVAIS SOMMETS=',I8,
     %        ' de QUALITE MINIMALE=',G15.6,
     %        ' de PIRE SOMMET',I8,' XYZ=',3G14.6)
21150 FORMAT( ' ITERATION',I3,
     %        ' NUMBER of BAD VERTICES=',I8,
     %        ' of MINIMUM QUALITY=',G15.6,
     %        ' of WORST VERTEX',I8,' XYZ=',3G14.6)

C        SAUVEGARDE DU NUMERO DES PLUS MAUVAIS SOMMETS
         MXTET = MIN( NBSOMM, MXTEXA, NBQINF )
ccc         print *,'amqutetra step 3 : NBSOMM=',nbsomm,' MXTEXA=',mxtexa,
ccc     %           ' nbqinf=',NBQINF,' mxtet=',mxtet
         CALL TRTATA( NOSOET, NOTEXA, MXTET )
C
C        REMISE EN ETAT DU TABLEAU NO1TSO
         DO NS=1,MXSOMM
            IF( NO1TSO(NS) .LT. 0 ) NO1TSO(NS) = -NO1TSO(NS)
         ENDDO

C        BOUCLE SUR LES PLUS MAUVAIS SOMMETS
C        -----------------------------------
ccc         IF( LANGAG .EQ. 0 ) THEN
ccc          print *,'amqutetra step 3 : NOMBRE MAUVAIS SOMMETS TRAITES='
ccc     %           ,MXTET
ccc         ELSE
ccc          print *,'amqutetra step 3 : NUMBER of BAD TREATED VERTICES='
ccc     %           ,MXTET
ccc         ENDIF

         DO NTQF=1,MXTET

C           LE NUMERO DU SOMMET DE QUALITE A AMELIORER
            NS = NOTEXA( NTQF )

C           AMELIORATION DE LA QUALITE DU SOMMET NS
            CALL QUA1ST( QUTEMN, ANTE2P, GRAND,  NS,
     %                   MXSOMM, XYZSOM, NPSOFR, NBDM,   NUDMEF,
     %                   NBSOTE, MXTETR, N1TEVI, NSTETR, MXTETA,
     %                   VOLUMT, QUALIT,
     %                   NO1TSO, MXTESO, N1TESO, NOTESO,
     %                   MXFAET, N1FEOC, N1FEVI, NFETOI, NOSOET,
     %                   IERR  )

ccc               PRINT*,'amqutetra step 2 : Sommet',NS,
ccc     %                ' VOLET1=',VOLET1,' QUAET1=',QUAET1
ccc               PRINT*,'amqutetra step 2:  Sommet',NS,
ccc     %                ' VOLET2=',VOLET2,' QUAET2=',QUAET2
cccC              POUR MARQUER LA FIN DU TRAITEMENT DU TETRAEDRE NT
ccc               print *

         ENDDO

C        DETECTION DES FACES APPARTENANT A 3 TETRAEDRES
         CALL CRHAFAVO( NBSOTE, MXTETR, NSTETR, NUDTETR,
     %                  MOFACE, MXFACE, LFACES, NBFAPB, IERR )
         IF( NBFAPB .GT. 0 ) THEN
            PRINT *,'amqutetra: 2',NBFAPB,' FACES DE 3 TETRAEDRES'
            TRACTE = .TRUE.
         ENDIF

C        TRACE FILAIRE DU MAILLAGE, ARETES ROUGES SI TETRA DE QUALITE<0.15
C        VOLUME ET QUALITE DU MAILLAGE FINAL
         CALL TRNSTETR( NBSOTE, MXTETR, NSTETR, NBCOOR, XYZSOM,
     %                  NBTETR, VOLET2, QUAET2, QUAMY2 )
         TRACTE = TRACTE0

C        LE TEMPS CPU D'EXECUTION EN SECONDES
         DCPU   = DINFO('DELTA CPU')
         DCPUTO = DCPUTO + DCPU

         print*
         IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'amqutetra step 2 : Qualite des MAUVAIS SOMMETS  -------
     %-----------------------------------------------------------------'
         ELSE
         PRINT*,'amqutetra step 2 : VERTICES BAD QUALITY  --------------
     %-----------------------------------------------------------------'
         ENDIF
         PRINT*,'amqutetra step 2 :',
     %                     ' VOLET1=',VOLET1,' QUAET1=',QUAET1,
     %                     ' QUAMY1=',QUAMY1
         PRINT*,'amqutetra step 2 :',
     %                     ' VOLET2=',VOLET2,' QUAET2=',QUAET2,
     %                     ' QUAMY2=',QUAMY2,
     %          ' VOL Diff=',ABS(VOLET2-VOLET1)/VOLET1*100,' %'

         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*)'FIN BOUCLE sur MAUVAIS SOMMETS =',
     %                      NINT(100D0*DCPU)/100.0,' SECONDES'
         ELSE
            WRITE(IMPRIM,*)'LOOP END on BAD VERTICES=',
     %                      NINT(100D0*DCPU)/100.0,' SECONDS'
         ENDIF


C        *******************************************************************
C        3-EME PARTIE : BOUCLE SUR LES SOMMETS INTERNES A UN SEUL MATERIAU
C                       REMPLACEMENT EVENTUEL PAR LE BARYCENTRE DES SOMMETS
C                       DE LEURS TETRAEDRES
C        *******************************************************************

C        BOUCLE DE BARYCENTRAGE DES SOMMETS INTERNES A UN MEME MATERIAU
C        --------------------------------------------------------------
         NBBARF = 0
         DO NS = 1, NBSOMM

            IF( NPSOFR( NS ) .EQ. 0 ) THEN
C              NS EST UN SOMMET INTERNE APPARTENANT A UN SEUL MATERIAU
               CALL QUA1SB( GRAND,  NS,     XYZSOM, NPSOFR,
     %                      NBSOTE, NSTETR, VOLUMT, QUALIT,
     %                      NO1TSO, NOTESO, NBBARF )
            ENDIF

         ENDDO

C        DETECTION DES FACES APPARTENANT A 3 TETRAEDRES
         CALL CRHAFAVO( NBSOTE, MXTETR, NSTETR, NUDTETR,
     %                  MOFACE, MXFACE, LFACES, NBFAPB, IERR )
         IF( NBFAPB .GT. 0 ) THEN
            PRINT *,'amqutetra: 3',NBFAPB,' FACES DE 3 TETRAEDRES'
            TRACTE = .TRUE.
         ENDIF

C        TRACE FINAL FILAIRE DU MAILLAGE, TRACE DES ARETES ROUGES
C        DES TETRAEDRES DE QUALITE<0.1
C        --------------------------------------------------------
         CALL TRNSTETR( NBSOTE, MXTETR, NSTETR, NBCOOR, XYZSOM,
     %                  NBTETR, VOLET3, QUAET3, QUAMY3 )
         TRACTE = TRACTE0

C        LE TEMPS CPU D'EXECUTION EN SECONDES
         DCPU   = DINFO('DELTA CPU')
         DCPUTO = DCPUTO + DCPU

         print*
         IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'amqutetra step 3 : Qualite des BARYCENTRES  -----------
     %-----------------------------------------------------------------'
         PRINT*,'amqutetra step 3 : NOMBRE d''ECHANGES SOMMET->BARYCENTR
     %E FAITS=',NBBARF
         ELSE
         PRINT*,'amqutetra step 3 : BARYCENTERS QUALITY  ---------------
     %-----------------------------------------------------------------'
         PRINT*,'amqutetra step 3 : NUMBER OF EXCHANGES VERTEX->BARYCENT
     %ER=',NBBARF
         ENDIF

         PRINT*,'amqutetra step 3 :',
     %                     ' VOLET2=',VOLET2,' QUAET1=',QUAET2,
     %                     ' QUAMY2=',QUAMY2
         PRINT*,'amqutetra step 3 :',
     %                     ' VOLET3=',VOLET3,' QUAET3=',QUAET3,
     %                     ' QUAMY3=',QUAMY3,
     %          ' VOL Diff=',ABS(VOLET3-VOLET2)/VOLET2*100,' %'

         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*)'FIN BOUCLE sur BARYCENTRAGES=',
     %                      NINT(100D0*DCPU)/100.0,' SECONDES'
         ELSE
           WRITE(IMPRIM,*)'LOOP END on VERTICES BARYCENTER=',
     %                     NINT(100D0*DCPU)/100.0,' SECONDS'
         ENDIF

         print*

C        ARRET DES ITERATIONS SI LA QUALITE EST STATIONNAIRE
C        ---------------------------------------------------
         IF( ITER .GE. 5 ) THEN
            IF( QUAET3 .GE. QUTEMN ) GOTO 9900
            IF( ABS(QUAMY3-QUAMY0) .LE. 0.005*QUAMY00 ) GOTO 9900
         ENDIF

 1000 CONTINUE

C     TRACE FINAL FILAIRE DU MAILLAGE, TRACE DES ARETES ROUGES
C     DES TETRAEDRES DE QUALITE<0.1
C     --------------------------------------------------------
 9900 PRINT*
      TRACTE = .TRUE.
      CALL TRNSTETR( NBSOTE, MXTETR, NSTETR, NBCOOR, XYZSOM,
     %               NBTETR, VOLET3, QUAET3, QUAMY3 )
      TRACTE = TRACTE0

C     LISTE DES TETRAEDRES DE QUALITE <QUTEMI=0.01
C     --------------------------------------------
      QUTEMI = 0.01
      CALL QUMESH(' Fin amqutetra: Qualite TETRAEDRES', ITER, GRAND,
     %              XYZSOM, NBSOTE, MXTETR, NSTETR,
     %              VOLUMT, QUALIT, QUTEMI,
     %              NBTU,   VOLT,   QUAMY0, QUAET0,
     %              NBQINF, QUALII, NOSOET )

cccC     VOLUME MOYEN D'UN TETRAEDRE
ccc      VOLMOYTE = VOLT / NBTU

C     TRI CROISSANT SELON LA QUALITE DES NBQINF MAUVAIS TETRAEDRES
      CALL TRITRP( NBQINF, QUALII, NOSOET )

C     AFFICHAGE DES PLUS MAUVAIS TETRAEDRES DE QUALITE <0.01
C     ET DESTRUCTION DE CES TETRAEDRES
C     ------------------------------------------------------
      DO NTQF=1,NBQINF

C        LE NUMERO NSTETR DU TETRAEDRE NTQF DE NOSOET
         NT = NOSOET( NTQF )
         IF( NSTETR(1,NT) .GT. 0 ) THEN
            PRINT*
            PRINT*,'amqutetra: PLUS MAUVAIS Tetraedre',NTQF,':',NT,
     %             ' Volume=',VOLUMT(NT),' Qualit=',QUALIT(NT)
            DO kk=1,4
               NS = NSTETR(KK,NT)
               PRINT*,'amqutetra: St',NS,' X=',XYZSOM(1,NS),
     %                ' Y=',XYZSOM(2,NS),' Z=',XYZSOM(3,NS)
            ENDDO
         ENDIF

C        SUPPRESSION DU TETRAEDRE NT DANS NSTETR
         CALL DS1TET( NT,     NBSOTE, N1TEVI, NSTETR,
     %                NO1TSO, N1TESO, NOTESO )
         PRINT*,'amqutetra: le TETRAEDRE NT=',NT,' est DETRUIT'
         QUALII( NT ) = GRAND
         VOLUMT( NT ) = 0
         IF( NBDM .GT. 0 ) NUDMEF( NT ) = 0
C        SUPPRESSION DU NO DE MATERIAU

      ENDDO

C     NETTOYAGE
      IF( IERALFACES .EQ. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'LFACES est DESALLOUE'
         ELSE
            PRINT*,'LFACES is DESALLOCATED'
         ENDIF
         DEALLOCATE( LFACES )
         IERALFACES = 1
      ENDIF

C     LE TEMPS FINAL D'EXECUTION EN SECONDES
      DCPU   = DINFO( 'DELTA CPU' )
      DCPUTO = DCPUTO + DCPU
      PRINT*,'amqutetra: sortie avec IERR=',IERR,' en',DCPU,' secondes'
      TRACTE = TRACTE0

      RETURN
      END
