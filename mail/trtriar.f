      SUBROUTINE TRTRIAR( NOMSP, XYZSOM, MXTRAR, L1ARET, L2ARET, MNARET,
     %                    NBTRIA, NOTRIA,
     %                    MXTRPV, NBTRPV, NBT2, NBT3, NOTRPV )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUCTION DU TABLEAU LARETE DES ARETES DU MAILLAGE
C -----    A PARTIR DU TABLEAU NOTRIA(4,NBTRIA) et
C          TRACE DES TRIANGLES NOTRPV DE NOTRIA, LEURS VOISINS et
C          LEURS VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT

C ENTREES:
C --------
C NOMSP  : NOM DE LA SUBROUTINE D'APPEL de trtriar.f
C XYZSOM : X  Y  Z DES SOMMETS DU MAILLAGE
C MXTRAR : NOMBRE MAXIMAL DE TRIANGLES ADJACENTS A UNE ARETE
C L1ARET : NOMBRE DE MOTS PAR ARETE DU TABLEAU NARET
C L2ARET : NOMBRE DE ARETES DU TABLEAU NARET
C MNARET : ADRESSE MCN DU TABLEAU LARETE DES ARETES DU MAILLAGE NOTRIA
C          TABLEAU DES ARETES DU MAILLAGE
C          LARETE(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          LARETE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LARETE(3,I)= CHAINAGE HACHAGE SUR ARETE SUIVANTE
C          LARETE(4,I)= NUMERO DU 1-ER TRIANGLE CONTENANT CETTE ARETE
C                       0 SI PAS DE 1-ER  TRIANGLE
C          LARETE(5,I)= NUMERO DU 2-EME TRIANGLE CONTENANT CETTE ARETE
C                       0 SI PAS DE 2-EME TRIANGLE
C NBTRIA : NOMBRE DE TRIANGLES DU TABLEAU NOTRIA
C NOTRIA : TABLEAU DU NUMERO DES 3 SOMMETS ET 0
C          SANS NUMERO NOTRIA des 3 TRIANGLES ADJACENTS
C MXTRPV : NOMBRE MAXIMAL DE NUMERO NOTRIA DE TRIANGLES DANS NOTRPV
C NBTRPV : NOMBRE INITIAL DE TRIANGLES DU TABLEAU NOTRPV

C SORTIES:
C --------
C NBTRPV : NOMBRE DE TRIANGLES DE NOTRPV INITIAUX (LES 0 RETIRES)
C NBT2   : NOMBRE DE TRIANGLES DE NOTRPV INITIAUX + ADJACENTS
C NBT3   : NOMBRE DE TRIANGLES DE NOTRPV INITIAUX + ADJACENTS + ADJACENTS
C NOTRPV : LISTE DES NBTRPV NUMEROS DANS NOTRIA DES TRIANGLES A TRACER
C          ( 1:NBTRPV ) TRIANGLE INITIAL
C          (NBTRPV+1:NBT2) TRIANGLES OPPOSES
C          (NBT2+1:NBT3) TRIANGLES OPPOSES OPPOSES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint Pierre du Perray             Mars 2020
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)

      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE
      CHARACTER*(*)     NOMSP
      CHARACTER*96      KTITRE
      INTEGER           NOTRIA(4,NBTRIA),
     %                  NOTRPV(MXTRPV), NOSOAR(2)
      REAL              XYZSOM(3,*), X(3), Y(3), Z(3), XYZ(3)

C     EN CAS DE SORTIE PREMATUREE...
      NBT2 = NBTRPV
      NBT3 = NBTRPV
      IF( .NOT. TRACTE  .OR. NBTRPV .LE. 0 ) RETURN

C     CONSTRUCTION DU TABLEAU LARETE DES ARETES DU MAILLAGE A PARTIR DE NOTRIA
C     ------------------------------------------------------------------------
C     EVENTUELLE DESTRUCTION DU TABLEAU DES ARETES DU MAILLAGE
      IF( MNARET .GT. 0 ) CALL TNMCDS( 'ENTIER', L1ARET*L2ARET, MNARET )

C     CONSTRUCTION DU TABLEAU DES ARETES DE LA TRIANGULATION
      CALL GEARSU( 4,      NBTRIA, NOTRIA, MXTRAR,
     %             L1ARET, L2ARET, MNARET, IERR )
C     ARETE(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C     ARETE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C     ARETE(3,I)= CHAINAGE HACHAGE SUR ARETE SUIVANTE
C     ARETE(4,I)= NUMERO DU 1-ER TRIANGLE CONTENANT CETTE ARETE
C                 0 SI PAS DE 1-ER  TRIANGLE
C     ARETE(5,I)= NUMERO DU 2-EME TRIANGLE CONTENANT CETTE ARETE
C                 0 SI PAS DE 2-EME TRIANGLE
C     ARETE(3+MXTRAR,I)= NO DU MXTRAR-EME TRIANGLE CONTENANT CETTE ARETE
C                        0 SI PAS DE MXTRAR-EME TRIANGLE
      LIBREF = L2ARET

C     CADRE DES NBTRPV TRIANGLES
C     --------------------------
C     RECHERCHE DU MIN ET MAX XYZ DES NBTRPV TRIANGLES
      RMAX   = RINFO( 'GRAND' )
C     MINIMUM XYZ
      COOEXT(1,1) = RMAX
      COOEXT(2,1) = RMAX
      COOEXT(3,1) = RMAX
C     MAXIMUM XYZ
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
      NBTRPV = NBT1

C     AJOUT DES TRIANGLES VOISINS DES NBT1 TRIANGLES INITIAUX
      NBT2 = NBT1
      DO I = 1, NBT1

         NT = NOTRPV(I)
         IF( NT .GT. 0 .AND. NOTRIA(1,NT) .GT. 0 ) THEN

            DO 1 J1=1,3

C              LES 2 SOMMETS DE L'ARETE J1 DU TRIANGLE NT
               NS1 = NOTRIA( J1, NT )
               IF( J1 .EQ. 3 ) THEN
                  J2 = 1
               ELSE
                  J2 = J1 + 1
               ENDIF
               NS2 = NOTRIA( J2, NT )

C              RECHERCHE DU TRIANGLE NTOP ADJACENT A NT PAR L'ARETE NS1-NS2
               IF( NS1 .LT. NS2 ) THEN
                  NOSOAR(1) = NS1
                  NOSOAR(2) = NS2
               ELSE
                  NOSOAR(1) = NS2
                  NOSOAR(2) = NS1
               ENDIF

C              NAR LE NUMERO DE L'ARETE NS1 NS2 DANS LE TABLEAU LARETE
               CALL HACHAR(2, NOSOAR, L1ARET,L2ARET,MCN(MNARET), 3, NAR)
C              NAR =0 SI LE TABLEAU NOSOAR N'A PAS ETE RETROUVE
C                  >0 SI LE TABLEAU NOSOAR   A     ETE RETROUVE

               IF( NAR .LE. 0 ) THEN
C                 PAS D'ARETE NOSOAR => PAS DE TRIANGLE ADJACENT AU TRIANGLE NT
                  NTOP = 0
                  GOTO 1
               ENDIF

C              QUEL EST LE NOMBRE DE TRIANGLES CONTENANT CETTE ARETE?
               MNA = MNARET + L1ARET * NAR - L1ARET -1
               DO J = 4, L1ARET
                  NTOP = ABS( MCN(MNA+J) )
                  IF( NTOP .GT. 0 .AND. NOTRIA(1,NTOP) .GT. 0 ) THEN
C                    NTOP EST IL DEJA STOCKE?
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
               ENDDO

 1          ENDDO

         ENDIF

      ENDDO

C     AJOUT DES TRIANGLES VOISINS DES TRIANGLES NBT1+1 A NBT2 C-A-D
C     DES VOISINS DE VOISINS DES NBT1 TRIANGLES INITIAUX
      NBT3 = NBT2
      DO I = NBT1+1, NBT2

         NT = NOTRPV(I)
         DO 3 J1=1,3

C           LES 2 SOMMETS DE L'ARETE J1 DU TRIANGLE NT
            NS1 = NOTRIA( J1, NT )
            IF( J1 .EQ. 3 ) THEN
               J2 = 1
            ELSE
               J2 = J1 + 1
            ENDIF
            NS2 = NOTRIA( J2, NT )

C           RECHERCHE DU TRIANGLE NTOP ADJACENT A NT PAR L'ARETE NS1-NS2
            IF( NS1 .LT. NS2 ) THEN
               NOSOAR(1) = NS1
               NOSOAR(2) = NS2
            ELSE
               NOSOAR(1) = NS2
               NOSOAR(2) = NS1
            ENDIF

C           NAR LE NUMERO DE L'ARETE NS1 NS2 DANS LE TABLEAU LARETE
            CALL HACHAR( 2, NOSOAR, L1ARET, L2ARET, MCN(MNARET), 3, NAR)
C           NAR =0 SI LE TABLEAU NOSOAR N'A PAS ETE RETROUVE
C               >0 SI LE TABLEAU NOSOAR   A     ETE RETROUVE

            IF( NAR .LE. 0 ) THEN
C              PAS D'ARETE NOSOAR =>
C              PAS DE TRIANGLE ADJACENT AU TRIANGLE NT
               NTOP = 0
               GOTO 3
            ENDIF

C           QUEL EST LE NOMBRE DE TRIANGLES CONTENANT CETTE ARETE?
            MNA = MNARET + L1ARET * NAR - L1ARET -1
            DO J = 4, L1ARET
               NTOP = ABS( MCN(MNA+J) )
               IF( NTOP .GT. 0 .AND. NOTRIA(1,NTOP) .GT. 0 ) THEN
C                 NTOP EST IL DEJA STOCKE?
                  DO K=1,NBT3
                     IF( NTOP .EQ. NOTRPV( K ) ) GOTO 3
                  ENDDO
                  IF( NBT3 .GE. MXTRPV ) THEN
                     GOTO 5
                  ENDIF
                  NBT3 = NBT3 + 1
                  NOTRPV( NBT3 ) = NTOP    
               ENDIF
            ENDDO

 3       ENDDO

      ENDDO


C     CALCUL DU MIN MAX DES XYZ DES SOMMETS DES NBT3 TRIANGLES
 5    DO I=1,NBT3
         NT = NOTRPV(I)
         IF( NT .GT. 0 .AND. NOTRIA(1,NT) .GT. 0 ) THEN
            DO J=1,3
C              LE J-EME SOMMET DE LA FACE NT
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
         R = ( COOEXT(K,2) - COOEXT(K,1) ) * 0.75
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
      PREDUF  = 15.0
      LORBITE = 1
      KTITRE  = NOMSP //
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
