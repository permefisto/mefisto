      SUBROUTINE TRTRIAN( NOMSP, XYZSOM,
     %                    M1TRIA, MXTRPV, NBTRPV, NOTRPV, NOTRIA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DES TRIANGLES NOTRPV DE NOTRIA, LEURS VOISINS et
C -----    LEURS VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT

C ENTREES:
C --------
C NOMSP  : NOM DE LA SUBROUTINE D'APPEL de trtrian.f
C XYZSOM : X  Y  Z DES SOMMETS DU MAILLAGE
C M1TRIA : NOMBRE DE MOTS POUR CHAQUE TRIANGLE NOTRIA
C          =6 ALORS NOTRIA: NS1 NS2 NS3 + NTROP1 NTOP2 NTROP3
C          =4 ALORS NOTRIA: NS1 NS2 NS3 + 0
C NBTRPV : NOMBRE INITIAL DE TRIANGLES DU TABLEAU NOTRPV
C NOTRPV : LISTE DES NBTRPV NUMEROS DANS NOTRIA DES TRIANGLES A TRACER
C          ( 1:NBTRPV ) TRIANGLE INITIAL
C          (NBTRPV+1:NBT2) TRIANGLES OPPOSES
C          (NBT2+1:NBT3) TRIANGLES OPPOSES OPPOSES
C NOTRIA : TABLEAU DU NUMERO DES 3 SOMMETS ET (0 ou 3 TRIANGLES ADJACENTS)
C          M1TRIA=6 NOTRIA: NS1 NS2 NS3 + NTROP1 NTOP2 NTROP3
C                =4 NOTRIA: NS1 NS2 NS3 + 0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY  OCTOBRE 2015
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
      CHARACTER*(*)     NOMSP
      CHARACTER*96      KTITRE
      INTEGER           NOTRIA(M1TRIA,*), NOTRPV(MXTRPV)
      REAL              XYZSOM(3,*), X(3), Y(3), Z(3), XYZ(3)

      IF( .NOT. TRACTE  ) RETURN

C     CADRE DES NBTRPV TRIANGLES
C     ==========================
      IF( NBTRPV .LE. 0 ) RETURN

C     RECHERCHE DU MIN ET MAX DES NBTRPV TRIANGLES
      RMAX   = RINFO( 'GRAND' )
C     MINIMUM
      COOEXT(1,1) = RMAX
      COOEXT(2,1) = RMAX
      COOEXT(3,1) = RMAX
C     MAXIMUM
      COOEXT(1,2) = -RMAX
      COOEXT(2,2) = -RMAX
      COOEXT(3,2) = -RMAX

C     RETRAIT DES VALEURS NEGATIVES DU TABLEAU NOTRPV
      NBT1 = 0
      DO I=1,NBTRPV
         NT = NOTRPV(I)
         IF( NT .GT. 0 .AND. NOTRIA(1,NT) .GT. 0 ) THEN
            NBT1 = NBT1 + 1
            NOTRPV( NBT1 ) = NT
         ENDIF
      ENDDO

ccc      mise a jour INTERDITE car APPEL POSSIBLE avec NBTRPV=Constante 1 
ccc      NBTRPV = NBT1

      IF( M1TRIA .LT. 6 ) THEN
C        PAS DE NO DES TRIANGLES VOISINS
         NBT2 = NBT1
         NBT3 = NBT1
         GOTO 5
      ENDIF

C     AJOUT DES TRIANGLES VOISINS DES NBT1 TRIANGLES INITIAUX
      NBT2 = NBT1
      DO I=1,NBT1
         NT = NOTRPV(I)
         IF( NT .GT. 0 .AND. NOTRIA(1,NT) .GT. 0 ) THEN
         DO 1 J=1,3
            NTOP = NOTRIA( 3+J, NT )
            IF( NTOP .GT. 0 .AND. NOTRIA(1,NTOP) .GT. 0 ) THEN
C              NTOP EST IL DEJA STOCKE?
               DO K=1,NBT2
                  IF( NTOP .EQ. NOTRPV( K ) ) GOTO 1
               ENDDO
               IF( NBT2 .GE. MXTRPV ) THEN
                  NBT3 = NBT2
                  GOTO 5
               ENDIF
               NBT2 = NBT2 + 1
               NOTRPV( NBT2 ) = NTOP
            ENDIF
 1       ENDDO
         ENDIF
      ENDDO

C     AJOUT DES TRIANGLES VOISINS DES VOISINS DES NBT1 TRIANGLES INITIAUX
      NBT3 = NBT2
      DO I=NBT1+1,NBT2
         NT = NOTRPV(I)
         DO 2 J=1,3
            NTOP = NOTRIA( 3+J, NT )
            IF( NTOP .GT. 0 .AND. NOTRIA(1,NTOP) .GT. 0 ) THEN
C              NTOP EST IL DEJA STOCKE?
               DO K=1,NBT3
                  IF( NTOP .EQ. NOTRPV( K ) ) GOTO 2
               ENDDO
               IF( NBT3 .GE. MXTRPV ) THEN
                  GOTO 5
               ENDIF
               NBT3 = NBT3 + 1
               NOTRPV( NBT3 ) = NTOP
            ENDIF
 2       ENDDO
      ENDDO

C     CALCUL DU MIN MAX DES XYZ DES SOMMETS DES NBT3 TRIANGLES
 5    DO I=1,NBT3
         NT = NOTRPV(I)
         IF( NT .GT. 0 .AND. NOTRIA(1,NT) .GT. 0 ) THEN
         DO J=1,3
C           LE J-EME SOMMET DE LA FACE NT
            NS = NOTRIA( J, NT )
            DO K=1,3
               R = XYZSOM(K,NS)
               COOEXT(K,1) = MIN( COOEXT(K,1), R )
               COOEXT(K,2) = MAX( COOEXT(K,2), R )
            ENDDO
         ENDDO
         ENDIF
      ENDDO

C     ELARGISSEMENT DU CADRE
      DO K=1,3
         R = ( COOEXT(K,2) - COOEXT(K,1) ) * 1.5
         COOEXT(K,1) = COOEXT(K,1) - R
         COOEXT(K,2) = COOEXT(K,2) + R
      ENDDO

C     LA MEMOIRE PIXELS EST EFFACEE
      CALL EFFACEMEMPX

C     PARAMETRES DE VISEE
      AXOLAR = 0
      DO L=1,3
         AXOPTV(L) = ( COOEXT(L,1) + COOEXT(L,2) ) * 0.5
         AXOEIL(L) = AXOPTV(L) + COOEXT(L,2)
         AXOLAR = MAX( AXOLAR, COOEXT(L,2) - COOEXT(L,1) )
      ENDDO

      AXOLAR = AXOLAR * 0.25
      AXOHAU = AXOLAR * 0.75
C     PAS DE PLAN ARRIERE ET AVANT
      AXOARR = 0
      AXOAVA = 0
      CALL AXONOMETRIE( AXOPTV, AXOEIL, AXOLAR, AXOHAU, AXOARR, AXOAVA )

C     TRACE EN MODE 3 BOUTONS POUR DEPLACEMENT ROTATIONS ZOOM
      PREDU0  = PREDUF
      PREDUF  = 12.0
      LORBITE = 1
      KTITRE = NOMSP //
     %': les TRIANGLES, LEURS ADJACENTS et LEURS ADJACENTS d''ADJACENTS'

      IAVNSO0 = IAVNSO
      IAVNSO  = 1

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

C     TRACE DES TRIANGLES ET DE LEURS TRIANGLES ADJACENTS
C     ET DE LEURS TRIANGLES ADJACENT D'ADJACENTS
C     DANS L'ORDRE INVERSE POUR MIEUX VOIR LES INITIAUX
      DO 30 K=NBT3,1,-1

C        TRACE DU TRIANGLE ET SES ADJACENTS D'ADJACENTS
         NT = NOTRPV( K )
         IF( NT .LE. 0 ) GOTO 30
         IF( NOTRIA(1,NT) .LE. 0 ) GOTO 30

         IF( K .LE. NBT1 ) THEN
C           TRIANGLE INITIAL
            NCF = NCCYAN
            NCA = NCBLEU
            GOTO 25
         ELSE IF( K .LE. NBT2 ) THEN
C           TRIANGLE OPPOSE AU TRIANGLE INITIAL
            NCF = NCORAN
            NCA = NCROUG
            GOTO 25
         ELSE
C           TRIANGLE OPPOSE AUX TRIANGLES OPPOSES
            NCF = NCGRIM
            NCA = NCGRIS
            GOTO 25
         ENDIF

C        TRACE DU TRIANGLE NT
 25      DO I=1,3
            NS = NOTRIA( I, NT )
            X(I) = XYZSOM(1,NS)
            Y(I) = XYZSOM(2,NS)
            Z(I) = XYZSOM(3,NS)
         ENDDO
         CALL FAP13D( NCF, NCA, PREDUF, 3, X, Y, Z )

C        TRACE DU NUMERO DES 3 SOMMETS DU TRIANGLE
         DO I=1,3
            NS = NOTRIA( I, NT )
            CALL ENTIER3D( NCNOIR, XYZSOM(1,NS), NS )
         ENDDO

C        TRACE  DU NUMERO NT DU TRIANGLE EN SON BARYCENTRE
         DO I=1,3
            XYZ( I ) = ( XYZSOM( I, NOTRIA(1,NT) )
     %                 + XYZSOM( I, NOTRIA(2,NT) )
     %                 + XYZSOM( I, NOTRIA(3,NT) ) ) / 3
         ENDDO
         CALL ENTIER3D( NCMAGE, XYZ, NT )

 30   ENDDO

      CALL TRFINS( KTITRE )

C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 10

 9000 PREDUF = PREDU0
      IAVNSO = IAVNSO0

      RETURN
      END
