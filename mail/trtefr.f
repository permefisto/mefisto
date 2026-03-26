      SUBROUTINE TRTEFR( KTITRE,  NBSOM1, NBSOM2, PTXYZD,
     %                   NUDTETR, NOTETR, NTARMX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DES ARETES DES FACES FRONTIERE DES TETRAEDRES
C -----    ET DU TETRAEDRE NTARMX D'ARETE MAXIMALE

C ENTREES:
C --------
C KTITRE : TITRE DU TRACE COMPLETE PAR LE NOMBRE DE TETRAEDRES
C NBSOM1 : PREMIER NO DE SOMMET A TRACER
C NBSOM2 : DERNIER NO DE SOMMET A TRACER
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C NUDTETR: NUMERO NOTETR DU DERNIER TETRAEDRE ACTIF
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C NTARMX : TETRAEDRE D'ARETE MAXIMALE
C          =0 PAS DE TRACE DE CE TETRAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY Novembre 2018
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
C
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
C     STOPTE = FAUX ==> PAS D'ARRET APRES LE TRACE DE CHAQUE TETRAEDRE
C              VRAI ==> DEMANDE D'UN CARACTERE POUR REDEMARRER
C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C              VRAI ==> TRACE DES TETRAEDRES DE L'ETOILE
C
      DOUBLE PRECISION  PTXYZD(4,*)
      INTEGER           NOTETR(8,*)
      CHARACTER*(*)     KTITRE
      REAL              R, XYZ(3), XYZ2(3)
C     LES SOMMETS SONT VUS A PARTIR DU BARYCENTRE DANS LE SENS DIRECT
      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE/ 1,2,3, 2,4,3, 3,4,1, 4,2,1 /


      IF( .NOT. TRACTE ) RETURN

C     CADRE COOEXT RESTREINT AUX FACES A TRACER
      DO L=1,3
C        LE MINIMUM
         COOEXT(L,1) =  1E25
C        LE MAXIMUM
         COOEXT(L,2) = -1E25
      ENDDO

C     CALCUL DU MIN MAX DES XYZ DES TETRAEDRES
      DO NTE = 1, NUDTETR
         IF( NOTETR(1,NTE) .GT. 0 ) THEN
            DO K=1,4
C              SOMMET L DU TETRAEDRE
               NS = NOTETR(K,NTE)
               DO L=1,3
                  R = REAL( PTXYZD(L,NS) )
C                 LE MINIMUM
                  COOEXT(L,1) = MIN( COOEXT(L,1), R )
C                 LE MAXIMUM
                  COOEXT(L,2) = MAX( COOEXT(L,2), R )
               ENDDO
            ENDDO
         ENDIF

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
      PREDUF  = 0.0
      LORBITE = 1

      IF( LORBITE .EQ. 0 ) GOTO 20
C
C     INITIALISATION DE L'ORBITE
C     ==========================
      CALL ORBITE0( NOTYEV )
      GOTO 20
C
C     TRACE SELON L'ORBITE OU ZOOM OU TRANSLATION ACTIFS
C     ==================================================
 10   CALL ORBITE1( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 9000
C
C     TRACE DES AXES 3D
 20   CALL TRAXE3
C
C     LE TRACE DES TETRAEDRES FRONTIERE
      CALL XVEPAISSEUR( 1 )
      DO 50 NTE = 1, NUDTETR

C        NUMERO NOTETR DU TETRAEDRE A TRACER
         IF( NOTETR(1,NTE) .GT. 0 ) THEN

C           LE TETRAEDRE A T IL UN SOMMET SUPERIEUR A NBSOM1?
            DO N=1,4
               NS = NOTETR( N, NTE )
               IF( NS .GE. NBSOM1 ) THEN
                  CALL XVEPAISSEUR(1)
                  NCA = NCTURQ
                  CALL TRTETRA( NCA, NOTETR(1,NTE), PTXYZD )
               ENDIF
            ENDDO

C           LE TETRAEDRE A T IL UNE FACE FRONTIERE?
            DO NF=1,4
               NTOP = NOTETR( 4+NF, NTE )
               IF( NTOP .LE. 0 ) THEN

                  IF( NTOP .EQ. 0 ) THEN
C                    FACE FRONTIERE
                     CALL XVEPAISSEUR( 1 )
                     NCA = NCGRIM
                  ELSE
C                    FACE DE TETRAEDRE OPPOSE INCONNU
                     CALL XVEPAISSEUR( 3 )
                     NCA = NCROUG
                  ENDIF

C                 TRACE DE LA FACE NF DE NTE
                  NS1 = NOTETR( NOSOFATE(3,NF), NTE )
                  XYZ(1) = REAL( PTXYZD(1,NS1) )
                  XYZ(2) = REAL( PTXYZD(2,NS1) )
                  XYZ(3) = REAL( PTXYZD(3,NS1) )

                  DO K=1,3
                     NS2 = NOTETR( NOSOFATE(K,NF), NTE )
                     XYZ2(1) = REAL( PTXYZD(1,NS2) )
                     XYZ2(2) = REAL( PTXYZD(2,NS2) )
                     XYZ2(3) = REAL( PTXYZD(3,NS2) )
                     CALL TRAIT3D(  NCA,    XYZ,  XYZ2 )
                     CALL ENTIER3D( NCNOIR, XYZ2, NS2  )

                     NS1    = NS2
                     XYZ(1) = XYZ2(1)
                     XYZ(2) = XYZ2(2)
                     XYZ(3) = XYZ2(3)
                  ENDDO
               ENDIF
            ENDDO

         ENDIF

 50   ENDDO

C     TRACE DU TETRAEDRE NTARMX
      IF( NTARMX .GT. 0 .AND. NOTETR(1,NTARMX) .GT. 0 ) THEN
         CALL XVEPAISSEUR( 3 )
         CALL TRTETRA( NCORAN, NOTETR(1,NTARMX), PTXYZD )
      ENDIF

C     LE TRACE DU NUMERO DES SOMMETS NBSOM1 A NBSOM2
      DO NS = NBSOM1, NBSOM2
         XYZ(1) = REAL( PTXYZD(1,NS) )
         XYZ(2) = REAL( PTXYZD(2,NS) )
         XYZ(3) = REAL( PTXYZD(3,NS) )
         CALL SYMBOLE3D(NCCYAN, XYZ, '*' )
         CALL ENTIER3D( NCBLEU, XYZ,  NS )
      ENDDO

C     TITRE ET TRACE EFFECTIF
      CALL TRFINS( KTITRE )
C
C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 10

 9000 PREDUF = PREDU0
      CALL XVEPAISSEUR( 1 )
      RETURN
      END
