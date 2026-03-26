      SUBROUTINE TRACEETOILE( KTITRE, PTXYZD, NOTETR, NUPOIN,
     %                        NBTETR, NUTETR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ARETES DE LA LISTE NUTETR DES TETRAEDRES DE NOTETR
C -----
C ENTREES:
C --------
C KTITRE : TITRE DU TRACE
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C NUPOIN : >0 NUMERO DE PTXYZD DU POINT A TRACER
C          =0 SI PAS DE TRACE DEMANDE
C NBTETR : NOMBRE DE TETRAEDRES DE LA LISTE DES TETRAEDRES A TRACER
C NUTETR : NUMERO DANS NOTETR DES TETRAEDRES A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC PARIS - VEULETTES SUR MER AOUT 2014
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

      CHARACTER*(*)     KTITRE
      DOUBLE PRECISION  PTXYZD(4,*)
      INTEGER           NOTETR(8,*), NOSOTR(3)
      INTEGER           NUTETR(NBTETR)
      REAL              XYZ(3)
      INTEGER           NOSOFATE(3,4)
C     NO DES SOMMETS DES FACES POUR QUE VU DU SOMMET MANQUANT
C     LES SOMMETS SOIENT VUS DANS LE SENS DIRECT
      DATA              NOSOFATE/ 1,2,3, 2,4,3, 3,4,1, 4,2,1 /


      IF( .NOT. TRACTE .OR. NBTETR .LE. 0 ) RETURN

C     CADRE COOEXT RESTREINT AUX NBTETR ETRAEDRES
      DO L=1,3
C        LE MINIMUM
         COOEXT(L,1) =  1E25
C        LE MAXIMUM
         COOEXT(L,2) = -1E25
      ENDDO

      DO N=1,NBTETR
         NT = ABS( NUTETR(N) )
         IF( NT .GT. 0 .AND. NOTETR(1,NT) .GT. 0 ) THEN
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

      IF( NUPOIN .GT. 0 ) THEN
         DO L=1,3
            XYZ(L) = REAL( PTXYZD(L,NUPOIN) )
C           LE MINIMUM
            COOEXT(L,1) = MIN( COOEXT(L,1), XYZ(L) )
C           LE MAXIMUM
            COOEXT(L,2) = MAX( COOEXT(L,2), XYZ(L) )
         ENDDO
      ENDIF

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
      PREDUF  = 20.0
      LORBITE = 1
C     TRACE QUALITE DES EF => PAS D'OMBRAGE ET NON UTILISATION DE LA PALETTE
      LCRITR  = 1

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

      DO N=1,NBTETR
         NT = ABS( NUTETR(N) )
         IF( NT .GT. 0 .AND. NOTETR(1,NT) .GT. 0 ) THEN

cccC           LE TRACE DES 6 ARETES DU TETRAEDRE OPPOSE A UNE FACE DE NT
ccc            NC     = NCVERT
ccc            PREDUF = 20.0
ccc            DO K=1,4
ccc               NTEOP = NOTETR( 4+K, NT )
ccc               IF( NTEOP .GT. 0 .AND. NOTETR(1,NTEOP).GT.0 ) THEN
ccc                  CALL TRTETRA( NC, NOTETR(1,NTEOP), PTXYZD )
ccc               ENDIF
ccc            ENDDO

C           LE TRACE DES 6 ARETES DU TETRAEDRE NOTETR(*,NT)
            PREDUF = 10.0
            NCA    = NCBLEU
ccc            IF( N .EQ. NBTETR ) NCA=NCORAN
            CALL TRTETRA( NCA, NOTETR(1,NT), PTXYZD )

            DO NF=1,4
               NTOP = NOTETR( 4+NF, NT )
               IF( NTOP .LE. 0 ) THEN
                  IF( NTOP .EQ. 0 ) THEN
C                    FACE FRONTIERE
                     NCF = NCROSE
                     NCA = NCORAN
                  ELSE
C                    FACE DE TETRAEDRE OPPOSE INCONNU
                     NCF = NCTURQ
                     NCA = NCVERT
                  ENDIF
C                 NF FACE FRONTIERE EST TRACEE
                  DO K=1,3
                     NOSOTR(K) = NOTETR( NOSOFATE(K,NF), NT )
                  ENDDO
                  CALL TRFATR( NCF, NCA, NOSOTR, PTXYZD )
               ENDIF
            ENDDO

C           LE TRACE DU NUMERO DES SOMMETS DU TETRAEDRE NT
            DO I = 1, 4
               NS = NOTETR(I,NT)
               XYZ(1) = REAL( PTXYZD(1,NS) )
               XYZ(2) = REAL( PTXYZD(2,NS) )
               XYZ(3) = REAL( PTXYZD(3,NS) )
               CALL ENTIER3D( NCNOIR, XYZ, NS )
            ENDDO

         ENDIF
      ENDDO

C     LE TRACE DU POINT NUPOIN DE PXYZD
      IF( NUPOIN .GT. 0 ) THEN
         XYZ(1) = REAL( PTXYZD(1,NUPOIN) )
         XYZ(2) = REAL( PTXYZD(2,NUPOIN) )
         XYZ(3) = REAL( PTXYZD(3,NUPOIN) )
ccc         CALL SYMBOLE3D( NCGRIS, XYZ, 'O' )
         CALL SYMBOLE3D( NCROUG, XYZ, '+' )
         CALL ENTIER3D(  NCROUG, XYZ, NUPOIN )
      ENDIF

C     TITRE ET TRACE EFFECTIF
      CALL TRFINS( KTITRE )

C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 10

 9000 PREDUF = PREDU0

      RETURN
      END
