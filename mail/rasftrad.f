      SUBROUTINE RASFTRAD( COSMAXPL, NBSOM, XYZSOM, NBTRIA, NOTRIA,
     %                     RAP2TRMI, SFMOYTR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DETECTION et MODIFICATION DES COUPLES DE TRIANGLES ADJACENTS
C -----    PAR UNE ARETE et DE SURFACES TRES DIFFERENTES

C ENTREES:
C --------
C COSMAXPL: COSINUS DE L'ANGLE ENTRE LES 2 PLANS AU DESSOUS DUQUEL
C           LES 2 TRIANGLES SONT JUGES NON COPLANAIRES
C NBSOM  :  NOMBRE DE SOMMETS DE LA TRIANGULATION
C NBTRIA :  NOMBRE DE TRIANGLES INITIAUX

C MODIFIE:
C --------
C XYZSOM : X Y Z LES 3 COORDONNEES DES NBSOM SOMMETS DE LA TRIANGULATION
C NOTRIA : NUMERO XYZSOM DES 3 SOMMETS et NUMERO NOTRIA DES 3 TRIANGLES
C          ADJACENTS
C SORTIE :
C --------
C RAP2TRMI: RAPPORT MINIMUM des SURFACES de 2 TRIANGLES ADJACENTS
C SFMOYTR : SURFACE MOYENNE DES NBTRIA TRIANGLES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET  VEULETTES SUR MER                Fevrier 2020
C....................................................................012
      PARAMETER        (MXRATR=1024, MXTRST=64)
      include"./incl/langue.inc"
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0

      REAL            XYZSOM(3,NBSOM)
      INTEGER         NOTRIA(6,NBTRIA)
      INTEGER         LIRATR(MXRATR), NOTRST(MXTRST)
      REAL            STR, SURTRR, S, RAP2TR, RAP2TRMI, RAPARTR(3),
     %                XYZBAR(3), SFMOYTR

      TRACTE0  = TRACTE
      NBRATR   = 0
      NB2DIA   = 0
      RAP2TRMI = 1E27
      SFMOYTR  = 0.
      NBTR     = 0

      DO 100 NTR = 1, NBTRIA

         IF( NOTRIA(1,NTR) .EQ. 0 ) GOTO 100
C        LE TRIANGLE NTR EST ACTIF

C        COMPARAISON DES SURFACES DE NTR et des TRIANGLES ADJACENTS
         STR = SURTRR( XYZSOM( 1, NOTRIA(1,NTR) ),
     %                 XYZSOM( 1, NOTRIA(2,NTR) ),
     %                 XYZSOM( 1, NOTRIA(3,NTR) ) )
         NBRATR0 = NBRATR

         NBTR    = NBTR + 1
         SFMOYTR = SFMOYTR + STR

         DO 10 NA = 1, 3

C           LE TRIANGLE ADJACENT PAR L'ARETE NA DU TRIANGLE NTR
            NTRADJ = NOTRIA( 3+NA, NTR )
            IF( NTRADJ .LE. 0 ) GOTO 10

C           COMPARAISON DES SURFACES DES TRIANGLES NTR et NTRADJ
            S = SURTRR( XYZSOM( 1, NOTRIA(1,NTRADJ) ),
     %                  XYZSOM( 1, NOTRIA(2,NTRADJ) ),
     %                  XYZSOM( 1, NOTRIA(3,NTRADJ) ) )

C           LE RAPPORT DES SURFACES ENTRE 0 et 1
            IF( STR .LE. S ) THEN
               RAP2TR = STR / S
            ELSE
C              NON PRIS EN COMPTE POUR NE PAS TRAITER 2 FOIS LE MEME CAS
C              DU TRIANGLE ET DE SON TRIANGLE ADJACENT PAR CETTE ARETE NA
ccc               RAP2TR = S / STR
               GOTO 10
            ENDIF

            RAPARTR( NA ) = RAP2TR
            RAP2TRMI = MIN( RAP2TRMI, RAP2TR )

            IF( RAP2TR .LT. 0.0667 ) THEN

C              UNE ARETE DE RAPPORT DE SURFACES TROP PETIT
               NBRATR = NBRATR + 1
               PRINT*
               IF( LANGAG .EQ. 0 ) THEN
                  PRINT*,'rasftrad: Le RAPPORT',RAP2TR,
     %                   ' des SURFACES des TRIANGLES',NTR,NTRADJ,
     %                   ' leur TETRAEDRISATION PEUT POSER PROBLEME'
               ELSE
                  PRINT*,'rasftrad: The SURFACE RAPPORT',RAP2TR,
     %                   ' of TRIANGLES',NTR,NTRADJ,
     %                   ' MAY BE a PROBLEM TO TETRAHEDRIZE'
               ENDIF
               PRINT*,'rasftrad: TRIANGLE',NTR,' Surface=',STR,
     %                ' St:',(NOTRIA(K,NTR),K=1,3)
               PRINT*,'rasftrad: TRIANGLE',NTRADJ,' Surface=',S,
     %                ' St:',(NOTRIA(K,NTRADJ),K=1,3)

            ENDIF

 10      ENDDO

         IF( NBRATR0 .LT. NBRATR ) THEN
            PRINT*,'rasftrad: Les 3 RAPPORTS de SURFACES=',RAPARTR

C           TRACE DU TRIANGLE NTR DE NOTRIA, SES VOISINS et
C           SES VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT
            LIRATR( NBRATR ) = NTR
            NBT = 1
            CALL TRTRIAN( 'rasftrad 0', XYZSOM, 6,
     %                     MXRATR-NBRATR, NBT, LIRATR(NBRATR),
     %                     NOTRIA )

            IF( NBRATR0+3 .EQ. NBRATR ) THEN

C              LES 3 ARETES DE NTR PRESENTENT UN RAPPORT FAIBLE
C              LE TRIANGLE NTR EST REMPLACE PAR SON BARYCENTRE
C              ------------------------------------------------
C              CALCUL DU BARYCENTRE DU TRIANGLE NTR
               DO K=1,3
                  XYZBAR( K ) = 0D0
               ENDDO
               DO N=1,3
                  NS = NOTRIA( N, NTR )
                  DO K=1,3
                     XYZBAR( K ) = XYZBAR( K ) + XYZSOM( K, NS )
                  ENDDO
               ENDDO
               DO K=1,3
                  XYZBAR( K ) = XYZBAR( K ) / 3
               ENDDO

C              IDENTIFICATION DES 3 SOMMETS DE NTR EN
C              LE PREMIER SOMMET DU TRIANGLE NTR
               NS1 = NOTRIA( 1, NTR )

               DO K=1,3
                  XYZSOM( K , NS1 ) = XYZBAR( K )
               ENDDO

               DO N=2,3
                  NS = NOTRIA( N, NTR )
                  DO 30 NT=1,NBTRIA
                     IF( NT .NE. NTR ) THEN
                        DO K=1,3
                           IF( NOTRIA(K,NT) .EQ. NS ) THEN
                               NOTRIA(K,NT) = NS1
                               GOTO 30
                           ENDIF
                        ENDDO
                     ENDIF
 30               ENDDO
               ENDDO

C              DESTRUCTION DES 3 TRIANGLES ADJACENTS AU TRIANGLE NTR
               DO N=1,3
                  NTADJ = NOTRIA( 3+N, NTR )
C                 RECHERCHE D'UN TRIANGLE ADJACENT NON NTR
                  DO K=1,3
                     NTADAD = NOTRIA( K, NTADJ )
                     IF( NTADAD .NE. NTR ) THEN
C                       POUR TRACER UN TRIANGLE PROCHE NON DETRUIT
                        LIRATR( NBRATR ) = NTADAD
                        GOTO 40
                     ENDIF
                  ENDDO
 40               DO K=1,6
                     NOTRIA( K, NTADJ ) = 0
                  ENDDO
               ENDDO

C              DESTRUCTION DU TRIANGLE NTR
               DO K=1,6
                  NOTRIA( K, NTR ) = 0
               ENDDO

               PRINT*,'rasftrad: LE TRIANGLE',NTR,' St:',
     %                (NOTRIA(K,NTR),K=1,3),
     %                ' EST REMPLACE PAR SON BARYCENTRE NS=',NS

C              TRACE DU TRIANGLE NTR DE NOTRIA, SES VOISINS et
C              SES VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT
               NBT = 1
               CALL TRTRIAN( 'rasftrad 1', XYZSOM, 6,
     %                        MXRATR-NBRATR, NBT, LIRATR(NBRATR),
     %                        NOTRIA )

C              LES NBTRST NO NOTRIA DES TRIANGLES DE SOMMET NS
               CALL NOTRI1ST( NS,6,3,NBTRIA,NOTRIA,MXTRST,NBTRST,NOTRST)

C              TENTATIVE D'ECHANGE DES DIAGONALES DES ARETES D'ADJACENCE
               NB2DIA0 = NB2DIA
               DO NT=1,NBTRST
                  NTRIA = NOTRST( NT )
                  DO K = 1, 3
                     CALL EC2DIA( COSMAXPL, NTRIA,  K,  NOTRIA,
     %                            NBSOM,    XYZSOM, NTNEW )
                     IF( NTNEW .GT. 0 ) THEN
                        NB2DIA = NB2DIA + 1
                     ENDIF
                  ENDDO
               ENDDO

               IF( NB2DIA0 .LT. NB2DIA ) THEN
C                 TRACE DU TRIANGLE NTR DE NOTRIA, SES VOISINS et
C                 SES VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT
                  NBT = 1
                  CALL TRTRIAN( 'rasftrad 2', XYZSOM, 6,
     %                           MXRATR-NBRATR, NBT, LIRATR(NBRATR),
     %                           NOTRIA )
               ENDIF

               GOTO 100

            ENDIF


            IF( NBRATR0 .LT. NBRATR ) THEN

C              AU MOINS UNE ARETE DE NTR PRESENTE UN RAPPORT FAIBLE
C              ----------------------------------------------------
               DO NN = 1, 3

C                 LE SOMMET NN DU TRIANGLE NTR
                  NS = NOTRIA( NN, NTR )

C                 LES NBTRST NO NOTRIA DES TRIANGLES DE SOMMET NS
                CALL NOTRI1ST(NS,6,3,NBTRIA,NOTRIA,MXTRST,NBTRST,NOTRST)


                  IF( NBTRST .EQ. 4 ) THEN

C                    LE SOMMET NS DE NTR APPARTIENT A 4 TRIANGLES
C                    ALORS IL EST DEPLACE AU BARYCENTRE DE CES 4 TRIANGLES
C                    -----------------------------------------------------
C                    CALCUL DU BARYCENTRE DU QUADRANGLE
C                    ENGLOBANT LES NBTRST TRIANGLES (DONT NTR)
                     DO K=1,3
                        XYZBAR(K) = 0D0
                     ENDDO
                     NBS = 0
                     DO NT=1,NBTRST
                        NTRIA = NOTRST( NT )
                        DO I=1,3
                           NS2 = NOTRIA( I, NTRIA )
                           IF( NS2 .NE. NS ) THEN
                              NBS = NBS + 1
                              DO K=1,3
                                 XYZBAR(K) = XYZBAR(K) + XYZSOM(K,NS2)
                              ENDDO
                           ENDIF
                        ENDDO
                     ENDDO
C                    LE SOMMET NS EST DEPLACE AU BARYCENTRE
C                    DES NBTRST TRIANGLES
                     DO K=1,3
                        XYZSOM(K,NS) = XYZBAR(K) / NBS
                     ENDDO

                     PRINT*,'rasftrad: LE SOMMET',NS,
     %              ' EST DEPLACE AU BARYCENTRE DES',NBTRST,' TRIANGLES'
     %              ,(NOTRST(K),K=1,NBTRST)

C                    TRACE DU TRIANGLE NTR DE NOTRIA, SES VOISINS et
C                    SES VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT
                     NBT = 1
                     CALL TRTRIAN( 'rasftrad 3', XYZSOM, 6,
     %                             MXRATR-NBRATR, NBT, LIRATR(NBRATR),
     %                             NOTRIA )

C                    TENTATIVE D'ECHANGE DES DIAGONALES DES ARETES D'ADJACENCE
                     NB2DIA0 = NB2DIA
                     DO NT=1,NBTRST
                        NTRIA = NOTRST( NT )
                        DO K = 1, 3
                           CALL EC2DIA( COSMAXPL, NTRIA,  K,  NOTRIA,
     %                                  NBSOM,    XYZSOM, NTNEW )
                           IF( NTNEW .GT. 0 ) THEN
                              NB2DIA = NB2DIA + 1
                           ENDIF
                        ENDDO
                     ENDDO

                     IF( NB2DIA0 .LT. NB2DIA ) THEN
C                       TRACE DU TRIANGLE NTR DE NOTRIA, SES VOISINS et
C                       SES VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT
                        NBT = 1
                        CALL TRTRIAN( 'rasftrad 4', XYZSOM, 6,
     %                               MXRATR-NBRATR, NBT, LIRATR(NBRATR),
     %                                 NOTRIA )
                     ENDIF

                     GOTO 100

                  ENDIF


                  IF( NBTRST .EQ. 3 ) THEN

C                    LES 3 TRIANGLES DE SOMMET NS SONT REMPLACES
C                    PAR LE TRIANGLE ENGLOBANT DE 'CENTRE' NS
C                    -------------------------------------------
C                    LE TRIANGLE 1 DE SOMMET NS
                     NTRIA1 = NOTRST( 1 )

C                    NO LOCAL DU SOMMET NS DANS NTRIA1
                     DO N=1,3
                        IF( NOTRIA(N,NTRIA1) .EQ. NS ) GOTO 50
                     ENDDO
C                    L'ARETE OPPOSEE AU SOMMET NS DANS NTRIA1
 50                  IF( N .EQ. 3 ) THEN
                        N1 = 1
                     ELSE
                        N1 = N + 1
                     ENDIF
                     IF( N1 .EQ. 3 ) THEN
                        N2 = 1
                     ELSE
                        N2 = N1 + 1
                     ENDIF

C                    LES 2 SOMMETS DE L'ARETE N1 DE NTRIA1 C-A-D
C                    LES 2 PREMIERS SOMMETS DU TRIANGLE ENGLOBANT
                     NS1 = NOTRIA( N1, NTRIA1 )
                     NS2 = NOTRIA( N2, NTRIA1 )

C                    LE TRIANGLE OPPOSE A L'ARETE N1 NS1-NS2 DE NTRIA1
                     NTROP1 = NOTRIA( 3+N1, NTRIA1 )

C                    RECHERCHE DE L'AUTRE TRIANGLE DE SOMMET NS2
                     DO NT2 = 2, 3
                        NTRIA2 = NOTRST( NT2 )
                        DO K2=1,3
                           IF( NOTRIA(K2,NTRIA2) .EQ. NS2 ) GOTO 60
                        ENDDO
                     ENDDO

C                    NS3 LE 3-EME SOMMET DU TRIANGLE ENGLOBANT
 60                  DO N3=1,3
                        NS3 = NOTRIA( N3, NTRIA2 )
                        IF( NS3 .NE. NS2 .AND. NS3 .NE. NS ) GOTO 65
                     ENDDO

C                    RECHERCHE DE L'ARETE NS2-NS3 EXTERIEURE DE NTRIA2
 65                  IF( K2 .EQ. 3 ) THEN
                        K3 = 1
                     ELSE
                        K3 = K2 + 1
                     ENDIF
                     IF( NOTRIA( K3, NTRIA2 ) .NE. NS3 ) THEN
C                        NS2-NS3 EST L'ARETE PRECEDENTE
                        IF( K2 .EQ. 1 ) THEN
                           K2 = 3
                        ELSE
                           K2 = K2 - 1
                        ENDIF
                     ENDIF
         
C                    LE TRIANGLE OPPOSE A L'ARETE K2 DE NTRIA2
                     NTROP2 = NOTRIA( 3+K2, NTRIA2 )

C                    LE 3-EME TRIANGLE DE SOMMET NS
                     IF( NT2 .EQ. 2 ) THEN
                        NTRIA3 = NOTRST( 3 )
                     ELSE
                        NTRIA3 = NOTRST( 2 )
                     ENDIF

C                    N3 NO DE L'ARETE NS3-NS1 DANS NTRIA3?
                     DO N3=1,3
                        IF( NOTRIA( N3, NTRIA3 ) .EQ. NS3 ) GOTO 70
                     ENDDO

 70                  IF( N3 .EQ. 3 ) THEN
                        K3 = 1
                     ELSE
                        K3 = N3 + 1
                     ENDIF

                     IF( NOTRIA( K3, NTRIA3 ) .NE. NS1 ) THEN
C                       NS3-NS1 EST L'ARETE PRECEDENTE
                        IF( N3 .EQ. 1 ) THEN
                           N3 = 3
                        ELSE
                           N3 = N3 - 1
                        ENDIF
                     ENDIF
         
C                    LE TRIANGLE OPPOSE A L'ARETE N3 DE NTRIA3
                     NTROP3 = NOTRIA( 3+N3, NTRIA3 )

C                    DESTRUCTION DU TRIANGLE NTRIA2
                     DO K=1,6
                        NOTRIA( K, NTRIA2 ) = 0
                     ENDDO

C                    DESTRUCTION DU TRIANGLE NTRIA3
                     DO K=1,6
                        NOTRIA( K, NTRIA3 ) = 0
                     ENDDO

C                    LE TRIANGLE NTRIA1 DEVIENT LE TRIANGLE ENGLOBANT
C                    SES 3 SOMMETS
                     NOTRIA( 1, NTRIA1 ) = NS1
                     NOTRIA( 2, NTRIA1 ) = NS2
                     NOTRIA( 3, NTRIA1 ) = NS3

C                    SES 3 TRIANGLES ADJACENTS
                     NOTRIA( 4, NTRIA1 ) = NTROP1
                     NOTRIA( 5, NTRIA1 ) = NTROP2
                     NOTRIA( 6, NTRIA1 ) = NTROP3

C                    SON LISTAGE POUR LES TRACES
                     LIRATR( NBRATR ) = NTRIA1

                     PRINT*,'rasftrad: LES 3 TRIANGLES du SOMMET',NS,
     %              ' SONT REMPLACES PAR LE TRIANGLE LES ENGLOBANT'

C                    TRACE DU TRIANGLE NTRIA1 DE NOTRIA, SES VOISINS et
C                    SES VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT
                     NBT = 1
                     CALL TRTRIAN( 'rasftrad 5', XYZSOM, 6,
     %                             MXRATR-NBRATR, NBT, LIRATR(NBRATR),
     %                             NOTRIA )

C                    TENTATIVE D'ECHANGE DES DIAGONALES DES ARETES D'ADJACENCE
                     NB2DIA0 = NB2DIA
                     DO K = 1, 3
                        CALL EC2DIA( COSMAXPL, NTRIA1, K,  NOTRIA,
     %                               NBSOM,    XYZSOM, NTNEW )
                        IF( NTNEW .GT. 0 ) THEN
                           NB2DIA = NB2DIA + 1
                        ENDIF
                     ENDDO

                     IF( NB2DIA0 .LT. NB2DIA ) THEN
C                       TRACE DU TRIANGLE NTRIA1 DE NOTRIA, SES VOISINS et
C                       SES VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT
                        NBT = 1
                        CALL TRTRIAN( 'rasftrad 6', XYZSOM, 6,
     %                               MXRATR-NBRATR, NBT, LIRATR(NBRATR),
     %                                 NOTRIA )
                     ENDIF

                     GOTO 100

                  ENDIF

               ENDDO

            ENDIF

         ENDIF

 100  ENDDO

C     SURFACE MOYENNE DES NBTR TRIANGLES ACTIFS
      SFMOYTR = SFMOYTR / NBTR

C     AFFICHAGE DU MINIMUM DES RAPPORTS SUR LES ARETES DES TRIANGLES
      IF( NBRATR .GT. 0 ) THEN

         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'rasftrad: NOMBRE de RAPPORTS des SURFACES de 2 TRIAN
     %GLES ADJACENTS < 0.05 =',NBRATR,' RAPPORT MINIMUM=',RAP2TRMI
            PRINT*,'rasftrad: NOMBRE d''ECHANGES de DIAGONALES=',NB2DIA
            PRINT*,'rasftrad: SURFACE MOYENNE d''UN TRIANGLE=',SFMOYTR
         ELSE
            PRINT*,'rasftrad: 2 ADJACENT TRIANGLE SURFACES RATIO<0.05 NU
     %MBER=',NBRATR,' MINIMUM RATIO=',RAP2TRMI
            PRINT*,'rasftrad: DIAGONAL EXCHANGE NUMBER=',NB2DIA
            PRINT*,'rasftrad: TRIANGLE MEAN  SURFACE=',SFMOYTR
         ENDIF

C        TRACE DE NBRATR TRIANGLES DE NOTRIA, LEURS VOISINS et
C        LEURS VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT
         tracte = .true.
         CALL TRTRIAN( 'rasftrad 7', XYZSOM, 6, MXRATR, NBRATR, LIRATR,
     %                  NOTRIA )
         tracte = tracte0

      ENDIF

      RETURN
      END
