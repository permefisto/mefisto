      SUBROUTINE TRARFASU( NDIMES, NBSOM,  XYZSOM, NOSOEF,
     %                     L1ARFA, L2ARFA, NARFA,  NBARXF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER EN VERT      TOUTES LES ARETES APPARTENANT A  1 FACE
C -----           EN TURQUOISE TOUTES LES FACES ADJACENTES AUX ARETES 1F
C                 EN NOIR      TOUTES LES ARETES APPARTENANT A  2 FACES
C                 EN ROUGE     TOUTES LES ARETES APPARTENANT A >2 FACES
C                 EN ORANGE    TOUTES LES FACES ADJACENTES AUX ARETES >2F
C          D'UNE SURFACE 2D ou 3D
C
C ENTREES :
C ---------
C NDIMES : ESPACE DE TRAVAIL 2 ou 3
C NBSOM  : NOMBRE DE POINTS DE XYZSOM
C XYZSOM : LES 3 COORDONNEES DES POINTS DE L'OBJET
C NOSOEF : LES 4 NUMEROS DES SOMMETS DES FACES DE LA SURFACE
C L1ARFA : NOMBRE DE MOTS POUR UNE ARETE DU TABLEAU NARFA
C L2ARFA : NOMBRE MAXIMUM D'ARETES DU TABLEAU NARFA
C NARFA  : TABLEAU NARFA DES ARETES DU MAILLAGE
C          NARFA(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          NARFA(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          NARFA(3,I)= CHAINAGE HACHAGE SUR L'ARETE SUIVANTE
C          NARFA(4:L2ARFA,I)= NO NOSOEF DE LA FACE CONTENANT L'ARETE I
C          SI NB FACES>L2ARFA-3 NARFA(4:L1ARFA)=-NO FACE
C NBARXF : NBARXF(n) NOMBRE D'ARETES APPARTENANT A n FACES n=1,2,3
C                    POUR n=3 APPARTENANT a >2 FACES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Novembre 2015
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/traaxe.inc"
      include"./incl/mecoit.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      REAL              XYZSOM(3,NBSOM), XYZFACE(3,4),
     %                  XFACE(4),YFACE(4), ZFACE(4)
      INTEGER           NOSOEF(4,*), NARFA(L1ARFA,L2ARFA), NBARXF(3)
      CHARACTER*120     KTITRE
C
C     EFFACER LA MEMOIRE DE TRACE
      CALL EFFACEMEMPX
C
      IF( INTERA .GE. 3 ) THEN
         LORBITE = 1
      ELSE
         LORBITE = 0
      ENDIF
C
      IF( NDIMES .NE. 2 ) GOTO 300


C     SURFACE 2D: INITIALISATION DU ZOOM DEPLACEMENT
C     ==============================================
      CALL ZOOM2D0( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 9000
      GOTO 210
C
C     ZOOM OU TRANSLATION ACTIFS
 200  CALL ZOOM2D1( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 9000
C
C     TRACE EFFECTIF DES AXES
 210  NETAXE = 0
      CALL TRAXE2
C
C     TRACE DES ARETES DES FACES DE LA SURFACE
C     ----------------------------------------
      NBARXF(1) = 0
      NBARXF(2) = 0
      NBARXF(3) = 0

      NBAR = 0
      DO 250 K = 1, L2ARFA

            IF( NARFA(1,K) .GT. 0 ) THEN

               NBAR = NBAR + 1

               IF( NARFA(4,K) .EQ. 0 ) GOTO 250
               IF( NARFA(5,K) .EQ. 0 ) THEN

C                 ARETE APPARTENANT A 1 SEULE FACE
C                 -> FACE TURQUOISE avec ARETE VERTE
C                 ----------------------------------
                  NBARXF(1) = NBARXF(1) + 1
                  CALL XVEPAISSEUR( 1 )

C                 TRACE DES NBS ARETES DE LA FACE NF
                  NF = ABS( NARFA(4,K) )

                  IF( NOSOEF(4,NF) .EQ. 0 ) THEN
                     NBS = 3
                  ELSE
                     NBS = 4
                  ENDIF

                  DO L=1,NBS
                     XFACE(L) = XYZSOM( 1, NOSOEF(L,NF) )
                     YFACE(L) = XYZSOM( 2, NOSOEF(L,NF) )
                  ENDDO
                  CALL FACE2D( NCTURQ, NCCYAN, NBS, XFACE, YFACE )

C                 LE NUMERO DES 2 SOMMETS DE L'ARETE SIMPLE  K
                  CALL XVEPAISSEUR( 5 )
                  NS1 = NARFA(1,K)
                  NS2 = NARFA(2,K)
                  CALL TRAIT2D( NCVERT, XYZSOM(1,NS1), XYZSOM(2,NS1),
     %                                  XYZSOM(1,NS2), XYZSOM(2,NS2))

               ELSE IF( NARFA(6,K) .EQ. 0 ) THEN

C                 ARETE APPARTENANT A 2 FACES -> NOIR
C                 -----------------------------------
                  NBARXF(2) = NBARXF(2) + 1
                  CALL XVEPAISSEUR( 0 )
C                 LE NUMERO DES 2 SOMMETS DE L'ARETE K
                  NS1 = NARFA(1,K)
                  NS2 = NARFA(2,K)
                  CALL TRAIT2D( NCNOIR, XYZSOM(1,NS1), XYZSOM(2,NS1),
     %                                  XYZSOM(1,NS2), XYZSOM(2,NS2))

               ELSE

C                 ARETE APPARTENANT A >2 FACES ROSES avec l'ARETE ROUGE
C                 -----------------------------------------------------
                  NBARXF(3) = NBARXF(3) + 1
C                 TRACE DES FACES ADJACENTES DE L'ARETE DANS >2 FACES
                  CALL XVEPAISSEUR( 1 )
                  DO M = 4, L1ARFA
                     NF = ABS( NARFA(M,K) )
                     IF( NF .EQ. 0 ) GOTO 250

                     IF( NOSOEF(4,NF) .EQ. 0 ) THEN
                        NBS = 3
                     ELSE
                        NBS = 4
                     ENDIF

C                    TRACE EN ROSE DE LA FACE AVEC UNE ARETE TRIPLE ROUGE
                     DO L=1,NBS
                        XFACE(L) = XYZSOM( 1, NOSOEF(L,NF) )
                        YFACE(L) = XYZSOM( 2, NOSOEF(L,NF) )
                     ENDDO
                     CALL FACE2D( NCROSE, NCORAN, NBS, XFACE, YFACE )

                  ENDDO

C                 TRACE PLUS EPAIS DE L'ARETE K APPARTENANT A >2 FACES
                  CALL XVEPAISSEUR( 7 )
C                 LE NUMERO DES 2 SOMMETS DE L'ARETE K APPARTENANT A >2 FACES
                  NS1 = NARFA(1,K)
                  NS2 = NARFA(2,K)
                  CALL TRAIT2D( NCROUG, XYZSOM(1,NS1), XYZSOM(2,NS1),
     %                                  XYZSOM(1,NS2), XYZSOM(2,NS2) )

               ENDIF

            ENDIF

 250  ENDDO

      GOTO 500


C     SURFACE 3D: INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
C     =======================================================
 300  CALL ORBITE0( NOTYEV )
      GOTO 320
C
C     ORBITE OU ZOOM OU TRANSLATION ACTIFS
 310  CALL ORBITE1( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 9000
C
C     TRACE DES AXES
 320  NETAXE = 0
      CALL TRAXE3
C
C     TRACE DES ARETES DES FACES DE LA SURFACE
C     ----------------------------------------
      NBARXF(1) = 0
      NBARXF(2) = 0
      NBARXF(3) = 0
      NBAR = 0

      DO 350 K = 1, L2ARFA

         IF( NARFA(1,K) .GT. 0 ) THEN

            NBAR = NBAR + 1

            IF( NARFA(4,K) .EQ. 0 ) GOTO 350
            IF( NARFA(5,K) .EQ. 0 ) THEN

C              TRACE EN TURQUOISE DE LA FACE AVEC UNE ARETE SIMPLE VERTE
C              ---------------------------------------------------------
               NBARXF(1) = NBARXF(1) + 1
               CALL XVEPAISSEUR( 1 )

               NF = ABS( NARFA(4,K) )
               IF( NOSOEF(4,NF) .EQ. 0 ) THEN
                  NBS = 3
               ELSE
                  NBS = 4
               ENDIF

               DO M=1,NBS
                  DO N = 1,3
                     XYZFACE( N, M ) = XYZSOM( N, NOSOEF(M,NF) )
                  ENDDO
               ENDDO
               CALL FACE3D( NCTURQ, NCCYAN, NBS, XYZFACE )

C              LE NUMERO DES 2 SOMMETS DE L'ARETE K
               NS1 = NARFA(1,K)
               NS2 = NARFA(2,K)
               CALL XVEPAISSEUR( 5 ) 
               CALL TRAIT3D( NCVERT, XYZSOM(1,NS1), XYZSOM(1,NS2) )
               CALL ENTIER3D(NCBLEU, XYZSOM(1,NS1), NS1)
               CALL ENTIER3D(NCBLEU, XYZSOM(1,NS2), NS2)

            ELSE IF( NARFA(6,K) .EQ. 0 ) THEN

C              ARETE APPARTENANT A 2 FACES -> ARETES en NOIR
C              ---------------------------------------------
               NBARXF(2) = NBARXF(2) + 1
               CALL XVEPAISSEUR( 0 )

cccC              TRACE DE L'ARETE K DE NARFA de SOMMETS NS1-NS2
ccc               NS1 = NARFA(1,K)
ccc               NS2 = NARFA(2,K)
ccc               CALL TRAIT3D( NCNOIR, XYZSOM(1,NS1), XYZSOM(1,NS2) )

C              TRACE DES 2 FACES DE L'ARETE K
               NCF = -1
               NF  = ABS( NARFA(4,K) )
               IF( NOSOEF(4,NF) .EQ. 0 ) THEN
                  NBS = 3
               ELSE
                  NBS = 4
               ENDIF
               DO L=1,NBS
                  XFACE(L) = XYZSOM( 1, NOSOEF(L,NF) )
                  YFACE(L) = XYZSOM( 2, NOSOEF(L,NF) )
                  ZFACE(L) = XYZSOM( 3, NOSOEF(L,NF) )
               ENDDO
               CALL FAP13D( NCF, NCNOIR, PREDUF, NBS, XFACE,YFACE,ZFACE)

               NF  = ABS( NARFA(5,K) )
               IF( NOSOEF(4,NF) .EQ. 0 ) THEN
                  NBS = 3
               ELSE
                  NBS = 4
               ENDIF
               DO L=1,NBS
                  XFACE(L) = XYZSOM( 1, NOSOEF(L,NF) )
                  YFACE(L) = XYZSOM( 2, NOSOEF(L,NF) )
                  ZFACE(L) = XYZSOM( 3, NOSOEF(L,NF) )
               ENDDO
               CALL FAP13D( NCF, NCNOIR, PREDUF, NBS, XFACE,YFACE,ZFACE)

            ELSE

C              TRACE EN ROSE DE LA FACE AVEC UNE ARETE TRIPLE ROUGE
C              ----------------------------------------------------
               NBARXF(3) = NBARXF(3) + 1
               CALL XVEPAISSEUR( 1 )

C              TRACE EN ROUGE DES FACES ADJACENTES PAR CETTE ARETE K
               DO L = 4, L1ARFA

                  NF = ABS( NARFA(L,K) )
                  IF( NF .EQ. 0 ) GOTO 340
C                 TRACE DES 3 ARETES DE LA FACE NF
                  IF( NOSOEF(4,NF) .EQ. 0 ) THEN
                     NBS = 3
                  ELSE
                     NBS = 4
                  ENDIF

                  DO M=1,NBS
                     DO N = 1,3
                        XYZFACE(N,M) = XYZSOM( N, NOSOEF(M,NF) )
                     ENDDO
                  ENDDO
                  CALL FACE3D( NCROSE, NCORAN, NBS, XYZFACE )

               ENDDO

C              LE NUMERO DES 2 SOMMETS DE L'ARETE K -> ROUGE
 340           CALL XVEPAISSEUR( 7 )
               NS1 = NARFA(1,K)
               NS2 = NARFA(2,K)
               CALL TRAIT3D( NCROUG, XYZSOM(1,NS1), XYZSOM(1,NS2) )
               CALL ENTIER3D(NCBLEU, XYZSOM(1,NS1), NS1)
               CALL ENTIER3D(NCBLEU, XYZSOM(1,NS2), NS2)
 
            ENDIF

         ENDIF

 350  ENDDO

C     TRACE DU TITRE
C     ==============
 500  CALL XVEPAISSEUR( 0 )
      IF( LANGAG .EQ. 0 ) THEN
         KTITRE='          ARETES VERTES dans 1 FACE &          ARETES N
     %OIRES dans 2 FACES &           ARETES ROUGES dans >2 FACES'
         WRITE(KTITRE(39:47),'(I9)') NBARXF(2)
         WRITE(KTITRE(77:85),'(I9)') NBARXF(3)
      ELSE
         KTITRE='          GREEN EDGES in 1 FACE &           BLACK EDGES
     % in 2 FACES &            RED EDGES in >2 FACES'
         WRITE(KTITRE(35:43),'(I9)') NBARXF(2)
         WRITE(KTITRE(70:78),'(I9)') NBARXF(3)
      ENDIF
      WRITE(KTITRE(1:9),'(I9)') NBARXF(1)
      CALL SANSDBL( KTITRE, L )
      CALL TRFINS(  KTITRE(1:L) )


C     REPRISE DU TRACE SELON LORBITE
C     ==============================
      IF( LORBITE .GT. 0 ) THEN
         IF( NDIMES .EQ. 2 ) THEN
            GOTO 200
         ELSE
            GOTO 310
         ENDIF
      ENDIF

 9000 CALL XVEPAISSEUR( 0 )

      IF( NBARXF(1) .GT. 0 .OR. NBARXF(3) .GT. 0 ) THEN
         WRITE(IMPRIM,*) 'trarfasu:',NBAR,' ARETES'
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) NBARXF(1),' ARETES APPARTENANT A  1 FACE'
            WRITE(IMPRIM,*) NBARXF(2),' ARETES APPARTENANT A  2 FACES'
            WRITE(IMPRIM,*) NBARXF(3),' ARETES APPARTENANT A >2 FACES'
         ELSE
            WRITE(IMPRIM,*) NBARXF(1),' EDGES in  1 FACE'
            WRITE(IMPRIM,*) NBARXF(2),' EDGES in  2 FACES'
            WRITE(IMPRIM,*) NBARXF(3),' EDGES in >2 FACES'
         ENDIF
      ENDIF

      RETURN
      END
