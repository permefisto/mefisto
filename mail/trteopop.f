      SUBROUTINE TRTEOPOP( NBSOM, XYZSOM, NBTETR, NSTETR, NTE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACE DES ARETES DU TETRAEDRE NTE, DE SES 4 TETRAEDRES VOISINS,
C -----  ET DE LEURS TETRAEDRES VOISINS

C ENTREES:
C --------
C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
C XYZSOM : X  Y  Z DES NBSOM SOMMETS
C NBTETR : NOMBRE DE TETRAEDRES
C NSTETR : NO DES 4 SOMMETS DE CHAQUE TETRAEDRE et 0 0 0 0
C NTE    : NO NSTETR DU TETRAEDRE INITIAL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY         Novembre 2021
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"

      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
C     STOPTE = FAUX ==> PAS D'ARRET APRES LE TRACE DE CHAQUE TETRAEDRE
C              VRAI ==> DEMANDE D'UN CARACTERE POUR REDEMARRER
C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C              VRAI ==> TRACE DES TETRAEDRES DE L'ETOILE

      REAL              XYZSOM(3,NBSOM)
      INTEGER           NSTETR(8,NBTETR)

      CHARACTER*80      KTITRE

      INTEGER           NOTEVO(1+4+4*4), NOTEVO1(4), NOTEVO2(4,4)
      EQUIVALENCE      (NOTEVO( 2),NOTEVO1(1)),
     %                 (NOTEVO( 6),NOTEVO2(1,1)),
     %                 (NOTEVO(10),NOTEVO2(1,2)),
     %                 (NOTEVO(14),NOTEVO2(1,3)),
     %                 (NOTEVO(18),NOTEVO2(1,4))

      REAL              ARMIN, ARMAX, SURFTR(4), VOLUTE
      REAL              R, XYZ(3)

      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE / 1,3,2,  2,3,4, 3,1,4, 4,1,2 /

      IF( .NOT. TRACTE ) RETURN

C     AFFICHAGE DU VOLUME ET DE LA QUALITE DU TETRAEDRE NTE
      CALL QUATET( XYZSOM(1,NSTETR(1,NTE)),
     %             XYZSOM(1,NSTETR(2,NTE)),
     %             XYZSOM(1,NSTETR(3,NTE)),
     %             XYZSOM(1,NSTETR(4,NTE)),
     %             ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
      PRINT*,'trteopop: Trace du TETRAEDRE',NTE,
     %       ' : St',(NSTETR(k,NTE),k=1,4),
     %       ' V=',VOLUTE,' Q=',QUALTE
      DO K=1,4
         NS = NSTETR(K,NTE)
         PRINT*,' Sommet',NS,' X=',XYZSOM(1,NS),
     %          ' Y=',XYZSOM(2,NS),' Z=',XYZSOM(1,NS)
      ENDDO

C     CONSTRUCTION DU NO DES TETRAEDRES VOISINS ET VOISINS DES VOISINS
      NOTEVO(1) = NTE
      CALL TEOPFA( NBTETR, NSTETR, NTE, NOTEVO1 )

      DO K=1,4
         IF( NOTEVO1( K ) .EQ. 0 ) THEN
            PRINT*,' Face',K,' St',(NSTETR(NOSOFATE(M,K),NTE),M=1,3),
     %             ' SANS TETRAEDRE OPPOSE'
         ENDIF
      ENDDO

C     CONSTRUCTION DU NO DES TETRAEDRES VOISINS DES VOISINS
      DO K=1,4
         IF( NOTEVO1(K) .EQ. 0 ) THEN
            DO M=1,4
               NOTEVO2(M,K) = 0
            ENDDO
         ELSE
            CALL TEOPFA( NBTETR, NSTETR, NOTEVO1(K), NOTEVO2(1,K) )
         ENDIF
      ENDDO

C     CADRE COOEXT RESTREINT AUX TETRAEDRES A TRACER
      DO L=1,3
C        LE MINIMUM
         COOEXT(L,1) =  1E25
C        LE MAXIMUM
         COOEXT(L,2) = -1E25
      ENDDO

C     CALCUL DU MIN et MAX des COORDONNEES DU TETRAEDRE NTE
      DO N=1,4
C        SOMMET N DU TETRAEDRE NTE
         NS = NSTETR(N,NTE)
         DO L=1,3
            R = XYZSOM(L,NS)
C           LE MINIMUM
            COOEXT(L,1) = MIN( COOEXT(L,1), R )
C           LE MAXIMUM
            COOEXT(L,2) = MAX( COOEXT(L,2), R )
         ENDDO
      ENDDO

C     CALCUL DU MIN et MAX des COORDONNEES DES TETRAEDRES VV
      DO M = 1, 4
         DO K=1,4
            NT = NOTEVO2(M,K)
            IF( NT .GT. 0 ) THEN
               DO N=1,4
C                 SOMMET N DU TETRAEDRE
                  NS = NSTETR(N,NT)
                  DO L=1,3
                     R = XYZSOM(L,NS)
C                    LE MINIMUM
                     COOEXT(L,1) = MIN( COOEXT(L,1), R )
C                    LE MAXIMUM
                     COOEXT(L,2) = MAX( COOEXT(L,2), R )
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

C     PARAMETRES DE VISEE
      AXOLAR = 0
      DO L=1,3
         AXOPTV(L) = ( COOEXT(L,1) + COOEXT(L,2) ) * 0.5
         AXOEIL(L) = COOEXT(L,2)
         AXOLAR = MAX( AXOLAR, COOEXT(L,2) - COOEXT(L,1) )
      ENDDO
      AXOLAR = AXOLAR * 0.5
      AXOHAU = AXOLAR * 0.75
C     PAS DE PLAN ARRIERE ET AVANT
      AXOARR = 0
      AXOAVA = 0
      CALL AXONOMETRIE( AXOPTV, AXOEIL, AXOLAR, AXOHAU, AXOARR, AXOAVA )

      DISMOY = AXOLAR/30

C     TRACE EN MODE 3 BOUTONS POUR DEPLACEMENT ROTATIONS ZOOM
      PREDU0  = PREDUF
      PREDUF  = 10.0
      LORBITE = 1

      KTITRE = 'TRACE du TETRAEDRE               et ses VOISINS'
      WRITE(KTITRE(20:31),'(I12)') NTE
      CALL SANSDBL( KTITRE, M )

      IF( LORBITE .EQ. 0 ) GOTO 20

C     INITIALISATION DE L'ORBITE
C     ==========================
      CALL ORBITE0( NOTYEV )
      GOTO 20

C     TRACE SELON L'ORBITE OU ZOOM OU TRANSLATION ACTIFS
C     ==================================================
 10   CALL ORBITE1( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 9000

C     TRACE DES AXES 3D
 20   CALL TRAXE3

      DO 50 M = 21, 1, -1

C        LE TETRAEDRE A TRACER
         NT = NOTEVO( M )
         IF( NT .EQ. 0 ) GOTO 50

C        LE TRACE DES 6 ARETES DU TETRAEDRE NSTETR(*,NT)
         IF( M .GE. 6 ) THEN
C           TRACE D'UN TETRAEDRE VOISIN D'UN VOISIN DE NTE
            CALL XVEPAISSEUR( 0 )
            NCA = NCVERT
            GOTO 30
         ENDIF

         IF( M .GT. 1 ) THEN
C           TRACE D'UN TETRAEDRE VOISIN DE NTE PAR UNE FACE
            CALL XVEPAISSEUR( 2 )
            NCA = NCBLEU
         ELSE
C           TRACE DE NTE
            CALL XVEPAISSEUR( 4 )
            NCA = NCROUG
         ENDIF

C        TRACE DES 6 ARETES DU TETRAEDRE NT DE COULEUR NCA
 30      CALL TRTETRAR( NCA, NSTETR(1,NT), XYZSOM )

C        LE TRACE DU NUMERO DES SOMMETS DES TETRAEDRES DE NSTETR
         DO I = 1, 4
            NS = NSTETR(I,NT)
            XYZ(1) = XYZSOM(1,NS)
            XYZ(2) = XYZSOM(2,NS)
            XYZ(3) = XYZSOM(3,NS)
            CALL ENTIER3D( NCNOIR, XYZ, NS )
         ENDDO

 50   ENDDO

C     TITRE ET TRACE EFFECTIF
      CALL TRFINS( KTITRE )

C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 10

 9000 PREDUF = PREDU0

      RETURN
      END
