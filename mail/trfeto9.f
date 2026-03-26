      SUBROUTINE TRFETO9( KTITRE, PTXYZD, NTE, NOTETR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DES ARETES DU TETRAEDRE NTE, DE SES 4 TETRAEDRES VOISINS,
C -----    ET DE LEURS TETRAEDRES VOISINS

C ENTREES:
C --------
C KTITRE : TITRE DU TRACE COMPLETE PAR LE NOMBRE DE TETRAEDRES
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C NTE    : NO NOTETR DU TETRAEDRE INITIAL
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY          Octobre 2018
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

      DOUBLE PRECISION  PTXYZD(4,*),
     %                  ARMIN, ARMAX, SURFTR(4), VOLUTE
      INTEGER           NOTETR(8,*)
      INTEGER           NOTEVO(1+4+4*4), NOSOTR(3)
      CHARACTER*(*)     KTITRE
      REAL              R, XYZ(3)

      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE / 1,3,2,  2,3,4, 3,1,4, 4,1,2 /

      IF( .NOT. TRACTE ) RETURN

C     CONSTRUCTION DU NO DES TETRAEDRES VOISINS ET VOISINS DES VOISINS
      NBTETR = 1
      NOTEVO( 1 ) = NTE

      DO K=1,4
         NTEOP = NOTETR( 4+K, NTE )
         IF( NTEOP .GT. 0 .AND. NOTETR(1,NTEOP) .GT. 0 ) THEN
            NBTETR = NBTETR + 1
            NOTEVO( NBTETR ) = NTEOP
         ENDIF
      ENDDO
      NBTETR1 = NBTETR

      DO M=2,NBTETR1
         NTEOP = NOTEVO( M )
         DO 5 K=1,4
            NTEOPOP = NOTETR( 4+K, NTEOP )
            IF( NTEOPOP .GT. 0 .AND. NOTETR(1,NTEOPOP) .GT. 0 ) THEN
C              TETRAEDRE DEJA RECENSE?
               DO L=1,NBTETR
                  IF( NTEOPOP .EQ. NOTEVO(L) ) GOTO 5
               ENDDO
C              NON IL EST AJOUTE
               NBTETR = NBTETR + 1
               NOTEVO( NBTETR ) = NTEOPOP
            ENDIF
 5       ENDDO
      ENDDO
      NBTETR2 = NBTETR


C     a supprimer ensuite ...
      DO M = 1, NBTETR2
         NT = NOTEVO( M )
            CALL QUATETD( PTXYZD(1,NOTETR(1,NT)),
     %                    PTXYZD(1,NOTETR(2,NT)),
     %                    PTXYZD(1,NOTETR(3,NT)),
     %                    PTXYZD(1,NOTETR(4,NT)),
     %                    ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
         PRINT*,'trfeto9: NOTETR(',NT,'): St',(NOTETR(k,NT),k=1,4),
     %          ' Tetra Op',(NOTETR(k,NT),k=5,8),
     %          ' V=',VOLUTE,' Q=',QUALTE
      ENDDO


C     CADRE COOEXT RESTREINT AUX TETRAEDRES A TRACER
      DO L=1,3
C        LE MINIMUM
         COOEXT(L,1) =  1E25
C        LE MAXIMUM
         COOEXT(L,2) = -1E25
      ENDDO

C     CALCUL DU MIN et MAX des COORDONNEES DES NBTETR2 TETRAEDRES
      DO M = 1, NBTETR2
         NT = NOTEVO( M )
         DO K=1,4
C           SOMMET K DU TETRAEDRE
            NS = NOTETR(K,NT)
            DO L=1,3
               R = REAL( PTXYZD(L,NS) )
C              LE MINIMUM
               COOEXT(L,1) = MIN( COOEXT(L,1), R )
C              LE MAXIMUM
               COOEXT(L,2) = MAX( COOEXT(L,2), R )
            ENDDO
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

      DO 50 M = NBTETR2, 1, -1
         NT = NOTEVO( M )

C        TRACE DU TETRAEDRE NOTETR( NT )
C        LE TRACE DES 6 ARETES DU TETRAEDRE NOTETR(*,NT)
         IF( M .EQ. 1 ) THEN
            CALL XVEPAISSEUR( 5 )
            NCA = NCROUG
         ELSE IF( M .LE. NBTETR1 ) THEN
            CALL XVEPAISSEUR( 3 )
            NCA = NCTURQ
         ELSE
            CALL XVEPAISSEUR( 1 )
            NCA = NCGRIM
         ENDIF

C        LE TETRAEDRE NT A T IL DES FACES SUR LA FRONTIERE?
         NBFAFR = 0
         DO NFE=1,4
            IF( NOTETR(4+NFE,NT) .LE. 0 ) THEN
               NBFAFR = NBFAFR + 1
            ENDIF
         ENDDO

         IF( NBFAFR .GT. 0 ) THEN

C           PRESENCE DE FACES FRONTIERE: TRACE DES 4 FACES DU TETRAEDRE NT
            DO NFE=1,4

C              NUMERO DES 3 SOMMETS DE LA FACE NFE DU TETRAEDRE NT
               NOSOTR(1) = NOTETR( NOSOFATE(1,NFE), NT )
               NOSOTR(2) = NOTETR( NOSOFATE(3,NFE), NT )
               NOSOTR(3) = NOTETR( NOSOFATE(2,NFE), NT )

C              TRACE DU TRIANGLE NOSOTR
               NTOP = NOTETR(4+NFE,NT)
               IF( NTOP .EQ. 0 ) THEN
C                 FACE FRONTIERE
                  NCF = NCORAN
               ELSE IF( NTOP .LT. 0 ) THEN
C                 FACE INCONNUE
                  NCF = NCMAGE
               ELSE
C                 SEULES LES 3 ARETES SONT TRACEES
                  NCF = -1
               ENDIF
               CALL TRFATR( NCF, NCGRIS, NOSOTR, PTXYZD )

            ENDDO

         ENDIF

C        TRACE DES 6 ARETES DU TETRAEDRE NT DE COULEUR NCA
         CALL TRTETRA( NCA, NOTETR(1,NT), PTXYZD )

C        LE TRACE DU NUMERO DES SOMMETS DES TETRAEDRES DE NOTETR
         DO I = 1, 4
            NS = NOTETR(I,NT)
            XYZ(1) = REAL( PTXYZD(1,NS) )
            XYZ(2) = REAL( PTXYZD(2,NS) )
            XYZ(3) = REAL( PTXYZD(3,NS) )
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
