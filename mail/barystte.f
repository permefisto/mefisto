      SUBROUTINE BARYSTTE( QUAMED, NBSOMM,  PTXYZD, NPSOFR,
     %                     NBTETR, NUDTETR, MXTETR, NOTETR, N1TETS,
     %                     MXTEET, NOTEET )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    BARYCENTRAGE DES SOMMETS INTERNES DE LA TETRAEDRISATION
C -----    DONT AUCUN SOMMET VOISIN N'EST FRONTALIER
C          AVEC RECHERCHE D'UNE QUALITE MINIMALE DES TETRAEDRES
C          DE LEUR ETOILE

C ENTREES:
C --------
C QUAMED : QUALITE D'UN TETRAEDRE AU DESSOUS DUQUEL LA QUALITE
C          D'UN TETRAEDRE EST MEDIOCRE ET LEUR NOMBRE CALCULE
C NBSOMM : NOMBRE DES SOMMETS DE LA TETRAEDRISATION
C NPSOFR : NUMERO DE POSITION PAR RAPPORT A LA FRONTIERE DE CHAQUE POINT
C          LE SIGNE DEVIENT NEGATIF SI LE SOMMET EST DEPLACE
C          =  0 SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  3 SI LE POINT EST IMPOSE PAR L'UTILISATEUR
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C                    LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -1 SI LE POINT EST SOMMET RECONNU TROP PROCHE
C          = -3 SI LE POINT A ETE REFUSE DANS LA TETRAEDRISATION
C          = -4 SI LE POINT EST SOMMET D'OT NON TROP PROCHE
C          = -NPSOFR(I) SI POINT I DEPLACE SUR LA SURFACE
C             DANS LEUR SURFACE FERMEE INITIALE OU NO DE POINT INTERNE
C   ou     =0 SI SOMMET INTERNE
C          =1 SI SOMMET FRONTALIER
C          =2 SI SOMMET SUR INTERFACE ENTRE 2 MATERIAUX
C          =3 SI SOMMET IMPOSE PAR L'UTILISATEUR LORS DE LA DEFINITION
C                       DU VOLUME

C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES DECLARABLES
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412

C ENTREE ET SORTIE :
C ------------------
C NBTETR : NOMBRE DE TETRAEDRES ACTIFS DE NOTETR
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C N1TETS : NUMERO DANS NOTETR D'UN TETRAEDRE DE CHAQUE SOMMET

C TABLEAUX AUXILIAIRES :
C ----------------------
C MXTEET : NOMBRE MAXIMAL DE NUMERO NOTETR DE TETRAEDRES DANS NOTEET
C NOTEET : NUMERO NOTETR  DE TETRAEDRES (DE L'ETOILE D'UN SOMMET)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Saint PIERRE du PERRAY           Novembre 2018
C23456...............................................................012
      include"./incl/darete.inc"
      DOUBLE PRECISION  PTXYZD(4,NBSOMM)
      DOUBLE PRECISION  VTE, VOLTET, VTOTGR, VMOYGR,
     %                  XYZDIN(4), XYZBABA(4), XYZBAR(4), VOLTID, VNTE,
     %                  OMEGA, OMEGA1, D, D1, D2, D3, D4, VOLUMT, VMOY0
      DOUBLE PRECISION  ARMIN, ARMAX, SURFTR(4)
      INTEGER           NOTETR(8,MXTETR), NPSOFR(NBSOMM),
     %                  N1TETS(NBSOMM), NOTEET(MXTEET)

      NUDTETR = 0

      N1SOM = NBSOMM
      N2SOM = 1
      NPAS  =-1

      DO 200 ITER = 1, 2

C     BARYCENTRAGE PONDERE DES SOMMETS INTERIEURS
C     -------------------------------------------
      PRINT*,'barystte: ITER=',ITER,' BARYCENTRAGE DES SOMMETS INTERIEUR
     %S',N1SOM,' A',N2SOM,' PAR PAS de',NPAS
      DO 100 NS = N1SOM, N2SOM, NPAS

         NPSOFRNS = NPSOFR( NS )
         IF( ( NPSOFRNS .GE. -4 .AND. NPSOFRNS .LE. 0 ) .AND.
     %         N1TETS( NS ) .GT. 0 ) THEN

C           RECHERCHE DES NBTEET TETRAEDRES DE SOMMET NS
            CALL TETR1S( NS,     N1TETS, NOTETR,
     %                   NBTEET, MXTEET, NOTEET, IERR )

            IF( IERR .EQ. 2 ) THEN
C              CORRECTION: LE POINT NS N'EST PLUS UN SOMMET DES TETRAEDRES
               N1TETS( NS ) = 0
               GOTO 100
            ENDIF

            IF( NBTEET .LE. 0 ) THEN
               GOTO 100
            ENDIF

C           CALCUL DE XYZBABA BARYCENTRE DES NBTEET BARYCENTRES
C           ---------------------------------------------------
            DO K=1,4
C              PROTECTION DES XYZD DU POINT NS
               XYZBABA( K ) = 0D0
               XYZDIN(  K ) = PTXYZD( K, NS )
            ENDDO

            VMOY0  = 0D0
            QMIN0  = 2.
            QMOY0  = 0.
            NBBARY = 0
            DO 40 NT = 1, NBTEET

C              LE NT-EME TETRAEDRE DU SOMMET NS
               NTE = NOTEET( NT )

               IF( NTE .LE. 0 ) GOTO 40
               IF( NOTETR( 1, NTE ) .LE. 0 ) GOTO 40

               DO K = 1, 4
                  NST = NOTETR( K, NTE )
                  NPSOFRNST = NPSOFR( NST )
                  IF( NPSOFRNST .LT. -4 .OR. NPSOFRNST .GT. 0 ) GOTO 100
               ENDDO

C              VOLUME et QUALITE DU TETRAEDRE NT
               CALL QUATETD( PTXYZD( 1, NOTETR( 1, NTE ) ),
     %                       PTXYZD( 1, NOTETR( 2, NTE ) ),
     %                       PTXYZD( 1, NOTETR( 3, NTE ) ),
     %                       PTXYZD( 1, NOTETR( 4, NTE ) ),
     %                       ARMIN, ARMAX, SURFTR, VNTE, QNTE )
               VMOY0 = VMOY0 + VNTE
               QMOY0 = QMOY0 + QNTE
               IF( QNTE .LT. QMIN0 ) QMIN0 = QNTE

               IF( NTE .GT. NUDTETR ) NUDTETR = NTE

C              LE VOLUME du TETRAEDRE IDEAL NTE POUR LA TAILLE SOUHAITEE
C              DE SES ARETES EN CHACUN DE SES 4 SOMMETS
               VOLTID = ( ( ( PTXYZD(4,NOTETR(1,NTE))
     %                      + PTXYZD(4,NOTETR(2,NTE))
     %                      + PTXYZD(4,NOTETR(3,NTE))
     %                      + PTXYZD(4,NOTETR(4,NTE)) ) / 4 ) **3 ) / 6

C              LE VOLUME DU TETRAEDRE NTE
               VTE = VOLTET( PTXYZD(1,NOTETR(1,NTE)),
     %                       PTXYZD(1,NOTETR(2,NTE)),
     %                       PTXYZD(1,NOTETR(3,NTE)),
     %                       PTXYZD(1,NOTETR(4,NTE)) )

               IF( VTE .LT. VOLTID*1D-4 ) THEN
                  GOTO 40
               ENDIF

C              LE BARYCENTRE de NTE AVEC LES POIDS 1/TailleIdeale(St)
               D1 = 1D0 / PTXYZD( 4, NOTETR(1,NTE) )
               D2 = 1D0 / PTXYZD( 4, NOTETR(2,NTE) )
               D3 = 1D0 / PTXYZD( 4, NOTETR(3,NTE) )
               D4 = 1D0 / PTXYZD( 4, NOTETR(4,NTE) )
               D  = 1D0 / ( D1 + D2 + D3 + D4 )
               DO K=1,4
                  XYZBAR( K ) = ( PTXYZD( K, NOTETR(1,NTE) ) * D1
     %                          + PTXYZD( K, NOTETR(2,NTE) ) * D2
     %                          + PTXYZD( K, NOTETR(3,NTE) ) * D3
     %                          + PTXYZD( K, NOTETR(4,NTE) ) * D4 ) * D
               ENDDO

C              AJOUT DU BARYCENTRE DU TETRAEDRE NTE
               NBBARY = NBBARY + 1
               DO K=1,4
C                 BARYCENTRE DES BARYCENTRES
                  XYZBABA( K ) = XYZBABA( K ) + XYZBAR( K )
               ENDDO

 40         ENDDO

            VMOY0 = VMOY0 / NBTEET
            QMOY0 = QMOY0 / NBTEET
            IF( QMOY0 .LE. 0 ) THEN
               print*,'barystte: Barycentrage de NS=',NS,
     %                ' QMIN0=',QMIN0,' QMOY0=',QMOY0,
     %                ' VMOY0=',VMOY0,' => ABANDON du BARYCENTRAGE'
               GOTO 100
            ENDIF

C           XYZD DU BARYCENTRE DES BARYCENTRES
            IF( NBBARY .EQ. 0 ) THEN
               GOTO 100
            ENDIF
            D = 1D0 / DBLE( NBBARY )
            DO K=1,4
               XYZBABA( K ) = XYZBABA( K ) * D
            ENDDO

C           RECHERCHE DU DEPLACEMENT OPTIMAL DU SOMMET NS SUR NS-BARYCENTRE
C           ---------------------------------------------------------------
            NOMEGAMX = -1
            QMINGRMX = -2.0
            DO NOME = -4, 8, ITER

               OMEGA  = 0.25D0 * NOME
               OMEGA1 = 1D0 - OMEGA

C              DEPLACEMENT DE NS ENTRE NS-BARYCENTRE DES BARYCENTRES
               DO K=1,4
                  PTXYZD(K,NS) = OMEGA1 * XYZDIN(K) + OMEGA * XYZBABA(K)
               ENDDO

C              APRES DEPLACEMENT DU SOMMET NS
C              QUALITES ET VOLUMES DU GROUPE DES NBTEET TETRAEDRES
               CALL QUALGRTE( PTXYZD, MXTETR, NOTETR, NBTEET, NOTEET,
     %                        QUAMED, NBTMED,
     %                        QMINGR, QMOYGR, VTOTGR, VMOYGR )

C              LE GROUPE DE TETRAEDRES DE SOMMET NS EST IL DE
C              MEILLEURE QUALITE?
               IF( QMINGR .GT. QMINGRMX ) THEN
C                 LE DEPLACEMENT DE NS APPORTE UNE MEILLEURE QUALITE
                  QMINGRMX = QMINGR
                  NOMEGAMX = NOME
               ENDIF

            ENDDO

            IF( QMINGRMX .LE. 0 ) THEN

C              LE SOMMET NS RESTE A SA PLACE
C              -----------------------------
               DO K=1,4
                  PTXYZD( K, NS ) = XYZDIN( K )
               ENDDO
               GOTO 100

            ENDIF

C           LE SOMMET NS EST DEPLACE AU POINT OPTIMAL AVEC NOMEGAMX
C           -------------------------------------------------------
            OMEGA  = 0.2D0 * NOMEGAMX
            OMEGA1 = 1D0 - OMEGA

ccc            PRINT*,'barystte: No Sommet',NS,' No OMEGA=',NOMEGAMX,
ccc     %      ' MAX QUALITE MIN de l''ETOILE=',QMINGRMX

C           DEPLACEMENT DE NS ENTRE NS-BARYCENTRE DES BARYCENTRES
            DO K=1,4
               PTXYZD(K,NS) = OMEGA1 * XYZDIN(K) + OMEGA * XYZBABA(K)
            ENDDO

C           CALCUL de la TAILLE_IDEALE du SOMMET NS DEPLACE
            CALL TAILIDEA( NOFOTI, PTXYZD(1,NS), NCODEV, PTXYZD(4,NS) )

C           APRES DEPLACEMENT DU SOMMET NS
C           QUALITES ET VOLUMES DU GROUPE DES NBTEET TETRAEDRES
            CALL QUALGRTE( PTXYZD, MXTETR, NOTETR, NBTEET, NOTEET,
     %                     QUAMED, NBTMED,
     %                     QMINGR, QMOYGR, VTOTGR, VMOYGR )

            IF( QMINGR .LE. 0 ) THEN
               print*,'barystte: Attention: Barycentrage de NS=',ns,
     %                ' QMINGR=',QMINGR,' QMOYGR=',QMOYGR,
     %                ' VMOYGR=',VMOYGR
            ENDIF

         ENDIF

 100  ENDDO

C     QUALITE DES TETRAEDRES APRES ITER ITERATIONS DE BARYCENTRAGE
C     ------------------------------------------------------------
ccc      PRINT*,'barystte: ITER=',ITER,' Fin du BARYCENTRAGE'
      CALL QUALTETR( PTXYZD, MXTETR,  NOTETR,
     %               NBTETR, NUDTETR, QUAMIN, QUAMOY, VOLUMT )
  
C     CHANGEMENT DE SENS DE PARCOURS DES SOMMETS INTERNES
      K     = N1SOM
      N1SOM = N2SOM
      N2SOM = K
C     MULTIPLICATION PAR ITER POUR UN ASPECT ALEATOIRE DU BARYCENTRAGE
      NPAS  =-NPAS * ITER

 200  ENDDO

      RETURN
      END
