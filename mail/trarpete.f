      SUBROUTINE TRARPETE( MXFACO, LEFACO, PTXYZD, N1TETS, NOTETR,
     %                     MXTE1S, NOTE1S, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DES ARETES DES FACES FRONTIERE et INTERFACE
C -----    DU TABLEAU LEFACO et DES TETRAEDRES DE LEURS 3 SOMMETS

C ENTREES:
C --------
C MXFACO : NOMBRE MAXIMAL DECLARABLE DE FACES DU CONTOUR
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1 , VOLUME2 DE LA FACE
C          NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C
C          ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C          => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C          LEFACO(9,*)  -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C          HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C          LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C          NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C          SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C          NF  = LEFACO( 9, NF )  ...
C          LEFACO(11,.) = NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE, 0 SINON
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C N1TETS : N1TETS(I) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET I
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE DE NUMEROTATION
C          1: 123      2: 234      3: 341      4: 412
C MXTE1S : NOMBRE D'ENTIERS DU TABLEAU NOTE1S
C NOTE1S : TABLEAU AUXILIAIRE

C SORTIES:
C ---------
C IERR   : = 0 PAS D'ERREUR DETECTEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray         Septembre 2019
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
C
      DOUBLE PRECISION  PTXYZD(4,*)
      INTEGER           LEFACO(1:11,0:MXFACO)
      INTEGER           NOTETR(8,*),
     %                  N1TETS(*),
     %                  NOTE1S( MXTE1S )
      REAL              XYZ(3)
      CHARACTER*80      KTITRE

      IF( .NOT. TRACTE ) RETURN

C     TRACE SELON L'ORBITE OU ZOOM OU TRANSLATION ACTIFS
C     ==================================================
C     LE TRACE DES 3 ARETES DES FACES TRIANGULAIRES DU TABLEAU LEFACO
      DO 50 NF=1,MXFACO

         IF( LEFACO(1,NF) .GT. 0 .AND. LEFACO(11,NF) .LE. 0 ) THEN

C           FACE FRONTIERE APPARTENANT A AUCUN TETRAEDRE
            NBTE3S = 0
            DO NS = 1, 3

C              NO DU SOMMET NS DE LA FACE NF PERDUE
               NOST = LEFACO( NS, NF )

C              LES TETRAEDRES DE SOMMET NOST SONT AJOUTES
               CALL TETR1S( NOST,   N1TETS, NOTETR,
     %                      NBTE1S, MXTE1S, NOTE1S(NBTE3S+1), IERR )
               NBTE3S = NBTE3S + NBTE1S

            ENDDO

C           ELIMINATION DES TETRAEDRES REPETES
            NBT = 0
            DO 3 N = 1, NBTE3S
               NTE = NOTE1S( N )
               DO L = 1, NBT
                  IF( NTE .EQ. NOTE1S( L ) ) THEN
C                    NTE EXISTE DEJA DANS NOTE1S
                     GOTO 3
                  ENDIF
               ENDDO
               NBT = NBT + 1
               NOTE1S( NBT ) = NTE
 3          ENDDO
            NBTE3S = NBT
            IF( NBTE3S .LE. 0 ) GOTO 50

C           CADRE COOEXT RESTREINT AUX NBTE3S TETRAEDRES A TRACER
            DO L=1,3
C              LE MINIMUM
               COOEXT(L,1) =  1E25
C              LE MAXIMUM
               COOEXT(L,2) = -1E25
            ENDDO

            DO 5 N = 1, NBTE3S
               NT = ABS( NOTE1S(N) )
               DO K = 1, 4
                  NS = NOTETR( K, NT )
                  IF( NS .LE. 0 ) GOTO 5
                  DO L=1,3
                     R = REAL( PTXYZD(L,NS) )
C                    LE MINIMUM
                     COOEXT(L,1) = MIN( COOEXT(L,1), R )
C                    LE MAXIMUM
                     COOEXT(L,2) = MAX( COOEXT(L,2), R )
                  ENDDO
               ENDDO
 5          ENDDO

            AXOLAR = 0
            DO L=1,3
               AXOPTV(L) = ( COOEXT(L,1) + COOEXT(L,2) ) * 0.5
               AXOEIL(L) = COOEXT(L,2)
               AXOLAR = MAX( AXOLAR, COOEXT(L,2) - COOEXT(L,1) )
            ENDDO
            AXOLAR = AXOLAR * 0.5
            AXOHAU = AXOLAR * 0.75
C           PAS DE PLAN ARRIERE ET AVANT
            AXOARR = 0
            AXOAVA = 0
            CALL AXONOMETRIE( AXOPTV, AXOEIL, AXOLAR, AXOHAU,
     %                        AXOARR, AXOAVA ) 
C
C           TRACE EN MODE 3 BOUTONS POUR DEPLACEMENT ROTATIONS ZOOM
            PREDUF = 20.0

            LORBITE = 1
            IF( LORBITE .EQ. 0 ) GOTO 20
C
C           INITIALISATION DE L'ORBITE
C           ==========================
            CALL ORBITE0( NOTYEV )
            GOTO 20

C           TRACE SELON L'ORBITE OU ZOOM OU TRANSLATION ACTIFS
C            ==================================================
 10         CALL ORBITE1( NOTYEV )
            IF( NOTYEV .EQ. 0 ) GOTO 50

C           TRACE DES AXES 3D
 20         CALL TRAXE3

            DO N=1,NBTE3S
               NTE = NOTE1S( N )
C              TRACE DES 6 ARETES DU TETRAEDRE NTE
               CALL TRTETRA( NCGRIS, NOTETR(1,NTE), PTXYZD )
            ENDDO

C           TRACE DES 3 ARETES DE LA FACE PERDUE
            CALL TRARTR( NCROUG, LEFACO(1,NF), PTXYZD )

C           TRACE DU NO DES 3 SOMMETS DE LA FACE PERDUE
            DO N=1,3
               NS = LEFACO(N,NF)
               DO K=1,3
                  XYZ(K) = REAL( PTXYZD(K,NS) )
               ENDDO
               CALL ENTIER3D( NCORAN, XYZ, NS )
            ENDDO

C           TITRE ET TRACE EFFECTIF
            KTITRE='ARETES de la FACE PERDUE         et des         TETR
     %AEDRES de ses 3 SOMMETS'
            WRITE(KTITRE(26:32),'(I7)') NF
            WRITE(KTITRE(41:47),'(I7)') NBTE3S
            CALL TRFINS( KTITRE )

C           REPRISE DE L'ORBITE
C           ===================
            IF( LORBITE .NE. 0 ) GOTO 10

         ENDIF
 50   ENDDO

      RETURN
      END
