      SUBROUTINE TRNOTETD( MXTETR, NOTETR, PTXYZD, NBTETR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DES ARETES DE TOUS LES TETRAEDRES DEFINIS DANS NOTETR
C -----
C ENTREES:
C --------
C MXTETR : NOMBRE MAXIMAL DECLARABLE DE FACES DU CONTOUR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE

C SORTIE :
C --------
C NBTETR : NOMBRE ACTUEL DE TETRAEDRES RANGES DANS NOTETR ET TRACES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC ET ST PIERRE DU PERRAY   AVRIL 2008
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

      DOUBLE PRECISION  PTXYZD(4,*)
      INTEGER           NOTETR(8,MXTETR)
      DOUBLE PRECISION  ARMIN, ARMAX, SURFTR(4), VOLUME
      REAL              XYZ(3)

      IF( .NOT. TRACTE ) THEN
C         CALCUL SEUL DU NOMBRE DE TETRAEDRES ACTUELS
          CALL NBDTETRA( MXTETR, NOTETR, NBTETR, NDTETRA )
          RETURN
       ENDIF

C     CADRE COOEXT RESTREINT AUX MXTETR ETRAEDRES
      DO L=1,3
C        LE MINIMUM
         COOEXT(L,1) =  1E25
C        LE MAXIMUM
         COOEXT(L,2) = -1E25
      ENDDO

      DO NT=1,MXTETR
         IF( NOTETR(1,NT) .GT. 0 ) THEN
            DO K=1,4
C              NUMERO DU SOMMET K DU TETRAEDRE NT
               NS = NOTETR(K,NT)
               DO L=1,3
                  XYZ(L) = REAL( PTXYZD(L,NS) )
C                 LE MINIMUM
                  COOEXT(L,1) = MIN( COOEXT(L,1), XYZ(L) )
C                 LE MAXIMUM
                  COOEXT(L,2) = MAX( COOEXT(L,2), XYZ(L) )
               ENDDO
            ENDDO
         ENDIF
      ENDDO

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
      CALL TRAXE3

C     LE TRACE DES FACES TRIANGULAIRES DU TABLEAU NOTETR
C     CALCUL DU NOMBRE DE TETRAEDRES
 20   NBTETR = 0
      DO 50 NT=1,MXTETR

         IF( NOTETR(1,NT) .GT. 0 ) THEN

C           CALCUL DE LA QUALITE DU TETRAEDRE NT
            CALL QUATETD( PTXYZD(1,NOTETR(1,NT)),
     %                    PTXYZD(1,NOTETR(2,NT)),
     %                    PTXYZD(1,NOTETR(3,NT)),
     %                    PTXYZD(1,NOTETR(4,NT)),
     %                    ARMIN, ARMAX, SURFTR, VOLUME, QUALIT )

C           COULEUR ET EPAISSEUR DES ARETES SELON LA QUALITE
            IF( QUALIT .LT. 0.01 ) THEN
               NC = NCROUG
               CALL XVEPAISSEUR( 3 )
            ELSE
               NC = NCGRIS
               CALL XVEPAISSEUR( 0 )
            ENDIF

C           TRACE DU TETRAEDRE
            CALL TRTETRA( NC, NOTETR(1,NT), PTXYZD )
            NBTETR = NBTETR + 1
         ENDIF

 50   CONTINUE

C     TITRE ET TRACE EFFECTIF
      CALL TRFINS( 'Les ARETES des TETRAEDRES en ROUGE si QUALITE<0.01')

C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 10

 9000 RETURN
      END
