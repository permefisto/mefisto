      SUBROUTINE A1LGARTI( MXSOM,  NBSOM,  XYZSOM,
     %                     L1ARET, L2ARET, LARETE,
     %                     MXTRIA, NBTRIA, NSTRIA,  IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AVEC FONCTION TAILLE_IDEALE(x,y,z) ou DARETE VALEUR PAR DEFAUT
C -----    de include"./incl/darete.inc"
C          AJOUTER DES POINTS SUR CHAQUE ARETE DONT LES SOMMETS
C          NE VERIFIENT PAS LA TAILLE_IDEALE et
C          SOUS TRIANGULER LES TRIANGLES DONT AU MOINS UNE ARETE A
C          UN ou DES POINTS AJOUTES

C ENTREES:
C --------
C MXSOM  : NOMBRE MAXIMAL DE SOMMETS DE LA TRIANGULATION
C NBSOM  : NOMBRE DE SOMMETS DE LA TRIANGULATION
C XYZSOM : X Y Z LES 3 COORDONNEES DES SOMMETS DE LA TRIANGULATION
C L1ARET : NOMBRE DE MOTS PAR ARETE DU TABLEAU NARET
C L2ARET : NOMBRE DE ARETES DU TABLEAU NARET
C MXTRIA : NOMBRE MAXIMAL DE TRIANGLES DECLARABLES DANS NSTRIA
C NBTRIA : NOMBRE DE TRIANGLES DE LA TRIANGULATION

C MODIFIES :
C ----------
C LARETE : TABLEAU DES ARETES DU MAILLAGE
C          LARETE(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          LARETE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LARETE(3,I)= CHAINAGE HACHAGE SUR ARETE SUIVANTE
C          LARETE(4,I)= NUMERO DU 1-ER TRIANGLE CONTENANT CETTE ARETE
C                       0 SI PAS DE 1-ER  TRIANGLE
C          LARETE(5,I)= NUMERO DU 2-EME TRIANGLE CONTENANT CETTE ARETE
C                       0 SI PAS DE 2-EME TRIANGLE
C          LARETE(6,I)= NOMBRE DE POINTS AJOUTES SUR L'ARETE I
C          LARETE(7,I)= NUMERO XYZSOM DU DERNIER POINT AJOUTE SUR I
C NBSOM  : NOMBRE DE SOMMETS DE LA TRIANGULATION DANS XYZSOM
C NSTRIA : NUMERO XYZSOM DES 3 SOMMETS ET 0 POUR CHACUN DES NBTRIA TRIANGLES
C NBTRIA : NOMBRE DE TRIANGLES DE LA TRIANGULATION

C SORTIES:
C --------
C IERR   : 0 SI PAS D'ERREUR DETECTEE
C          1 SI TAILLE_IDEALE NON CALCULABLE EN UN SOMMET
C          2 SI TABLEAU XYZSOM SATURE (AUGMENTER MXSOM)
C          3 SI TABLEAU NSTRIA SATURE (AUGMENTER MXTRIA)
C          4 SI UNE ARETE N'EST PAS RETROUVEE DANS LARETE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC PARIS & VEULETTES sur MER Avril 2017
C....................................................................012
      include"./incl/darete.inc"
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0
      PARAMETER        (MXPTAR=1024, MXTRAT=16192)
      INTEGER           LARETE(L1ARET,L2ARET), NSTRIA(4,MXTRIA),
     %                  NOSOAR(2),NOAR(3), NOPTAR(0:MXPTAR),
     %                  NDPTNTR(3), LITRAT(MXTRAT)
      REAL              XYZSOM(3,MXSOM)
      DOUBLE PRECISION  LGARET, S, C, D, DD, XYZ(3,2), TI(2)

C     -------------------------------------------------------------------------
C     CONSTRUCTION DES POINTS SUR LES ARETES EN RESPECTANT TAILLE_IDEALE(x,y,z)
C     -------------------------------------------------------------------------
      TRACTE0 = TRACTE
      LIBREF  = L2ARET
      NBSOM0  = NBSOM
      NBTRIA0 = NBTRIA
      N3 = 0
 
      DO 50 NAR = 1, L2ARET

C        PAS DE POINT AJOUTE SUR L'ARETE NAR
         LARETE( 6, NAR ) = 0
         LARETE( 7, NAR ) = 0

         IF( LARETE(1,NAR) .EQ. 0 ) GOTO 50

C        XYZ DES 2 SOMMETS DE L'ARETE NAR DE LARETE
         DO N = 1, 2

            NS = LARETE( N, NAR )
            DO K = 1, 3
               XYZ( K, N ) = XYZSOM( K, NS )
            ENDDO

C           CALCUL DE TAILLE_IDEALE(X,Y,Z) AU SOMMET NS
            CALL TAILIDEA( NOFOTI, XYZ(1,N),  NCODEV, TI(N) )

         ENDDO

C        LONGUEUR DE L'ARETE
         LGARET = SQRT( ( XYZ(1,2) - XYZ(1,1) ) ** 2
     %                + ( XYZ(2,2) - XYZ(2,1) ) ** 2
     %                + ( XYZ(3,2) - XYZ(3,1) ) ** 2 )

         IF( NCODEV .EQ. 0 ) THEN
            TI(1) = LGARET
            TI(2) = LGARET
         ENDIF

C        LE NUMERO XYZSOM DES 2 SOMMETS NS1<NS2 DE L'ARETE NAR
         NS1 = LARETE( 1, NAR )
         NS2 = LARETE( 2, NAR )

C        LE MIN DES TAILLE_IDEALE() AUX 2 SOMMETS DE L'ARETE NAR
C        EST ELLE INFERIEURE A LA LONGUEUR ACTUELLE DE L'ARETE?
         NBS = 0

         IF( LGARET .GE. 1.366D0*MIN( TI(1), TI(2) ) ) THEN

C           L'ARETE NAR EST TROP LONGUE PAR RAPPORT A LA TAILLE SOUHAITEE
C           AJOUT DE POINTS SUR L'ARETE

            IF( TI(1) .GT. 1.366D0 * TI(2) ) THEN
C              INVERSION DU SENS D'AJOUT DES POINTS SUR L'ARETE
               NSENS = -1
               NS  = NS1
               NS1 = NS2
               NS2 = NS
               DO K=1,3
                  XYZ( K, 1 ) = XYZSOM( K, NS1 )
                  XYZ( K, 2 ) = XYZSOM( K, NS2 )
               ENDDO
               D     = TI(1)
               TI(1) = TI(2)
               TI(2) = D
            ELSE
               NSENS = 1
            ENDIF

            D = LGARET
            NOPTAR( NBS ) = NS1

C           AJOUT D'UN NOUVEAU POINT A DISTANCE TI(1) DE XYZ(1) ?
 10         IF( D .GE. MIN( TI(1), TI(2) )*1.366D0 ) THEN

C              OUI: AJOUT D'UN POINT A DISTANCE TI(1) DE XYZ(1)
               C  = TI(1) / D
               DD = 0D0
               DO K=1,3
                  S  = XYZ(K,1) + ( XYZ(K,2) - XYZ(K,1) ) * C
                  DD = DD + ( S - XYZ(K,1) ) ** 2
                  XYZ(K,1) = S
               ENDDO
               DD = SQRT( DD )

C              LONGUEUR RESTANTE A TRAITER DE L'ARETE
               D = D - DD

               IF( NBSOM .GE. MXSOM  ) GOTO 9992
               IF( NBS   .GE. MXPTAR ) GOTO 20
               NBSOM  = NBSOM + 1
               DO K=1,3
                  XYZSOM(K,NBSOM) = REAL( XYZ(K,1) )
               ENDDO
C              TAILLE IDEALE AU POINT NBSOM
               CALL TAILIDEA( NOFOTI, XYZ(1,1), NCODEV, TI(1) )

C              POINT AJOUTE DANS NOPTAR
               NBS = NBS + 1
               NOPTAR( NBS ) = NBSOM

               GOTO 10

            ENDIF

C           NS2 EST LE DERNIER POINT DE L'ARETE NAR
 20         NOPTAR( NBS + 1 ) = NS2

C           NOMBRE DE POINTS SUR L'ARETE NAR (EXTREMITES NON COMPRISES)
            LARETE( 6, NAR ) = NBS

C           NUMERO DU DERNIER POINT AJOUTE SUR L'ARETE NAR
C           AVEC LE SIGNE DE NSENS SENS D'AJOUT DES POINTS SUR L'ARETE
            LARETE( 7, NAR ) = NBSOM * NSENS

C           BARYCENTRAGE PONDERE DES POINTS AJOUTES POUR HOMOGENEISATION
            IF( NBS .GE. 2 ) THEN
               DO N = NBS, 1, -1
                  NS  = NOPTAR( N   )
                  NS1 = NOPTAR( N-1 )
                  NS2 = NOPTAR( N+1 )
                  DO K=1,3
                     XYZSOM(K,NS) = 0.7 *   XYZSOM(K,NS)
     %                            + 0.3 * ( XYZSOM(K,NS1)
     %                                    + XYZSOM(K,NS2) ) /2
                  ENDDO
               ENDDO
            ENDIF

         ENDIF

 50   ENDDO

      PRINT*,'a1lgarti:',NBSOM,' SOMMETS   dont',NBSOM-NBSOM0,
     %       ' SOMMETS AJOUTES sur  les ARETES'

      IF( NBSOM .EQ. NBSOM0 ) GOTO 1010
      NBSOM0 = NBSOM

C     --------------------------------------------------------------
C     CONSTRUCTION DES SOUS TRIANGLES DE TOUT TRIANGLE DONT AU MOINS 
C     SUR UNE DE CES 3 ARETES A ETE AJOUTE UN POINT
C     --------------------------------------------------------------
      NBTRIA0 = NBTRIA
      NBTRAT0 = 0
      NBTRAT  = 0
      LIBREF  = L2ARET
      DO 1000 NTR = 1, NBTRIA0

         IF( NSTRIA(1,NTR) .EQ. 0 ) GOTO 1000

         DO N1 = 1, 3
            IF( N1 .LT. 3 ) THEN
               N2 = N1 + 1
            ELSE
               N2 = 1
            ENDIF

            NS1 = NSTRIA( N1, NTR )
            NS2 = NSTRIA( N2, NTR )

            IF( NS1 .LT. NS2 ) THEN
               NOSOAR(1) = NS1
               NOSOAR(2) = NS2
            ELSE
               NOSOAR(1) = NS2
               NOSOAR(2) = NS1
            ENDIF

C           NAR LE NUMERO DE L'ARETE NS1 NS2 DANS LARETE
            CALL HACHAG( 2, NOSOAR, L1ARET, L2ARET, LARETE, 3,
     %                   LIBREF, NAR )
C           NAR =0 SI SATURATION DU TABLEAU LARETE
C               >0 SI LE TABLEAU NOSOAR A ETE RETROUVE
C               <0 SI LE TABLEAU NOSOAR A ETE AJOUTE

            IF( NAR .LE. 0 ) THEN
C             ARETE AJOUTEE OU NON TROUVEE => PROBLEME
              PRINT*,'a1lgarti: ARETE',NS1,NS2,' du TRIANGLE(',NTR,')=',
     %       ')=',(NSTRIA(L,NTR),L=1,3),' NON ARETE de la TRIANGULATION'
               IERR = 6
               GOTO 9999
            ENDIF

            NOAR( N1 ) = NAR

         ENDDO


C        TRIANGULATION DU TRIANGLE NTR EN FONCTION DU NOMBRE DE POINTS
C        AJOUTES SUR SES 3 ARETES
C        =============================================================
         NBPAT = LARETE(6,NOAR(1)) +LARETE(6,NOAR(2)) +LARETE(6,NOAR(3))
         IF( NBPAT .LE. 0 ) GOTO 1000

C        NOMBRE DE TRIANGLES A TRACER
         NBTRI0 = NBTRIA

C        NTR A T IL 2 ARETES SANS POINT AJOUTE?
C        NOMBRE D'ARETES DE NTR SANS POINT AJOUTE
         IF( NBPAT .GE. 2 ) THEN
            N0 = 0
            DO N1 = 1, 3
               IF( LARETE( 6, NOAR(N1) ) .EQ. 0 ) N0 = N0 + 1
            ENDDO
            IF( N0 .EQ. 2 ) GOTO 800
         ENDIF

         IF( NBPAT .EQ. 1 ) THEN

C           1 SEUL POINT AJOUTE SUR L'UNE DES 3 ARETES DU TRIANGLE NTR
c           ==========================================================
C           RECHERCHE DE L'ARETE AVEC 1 POINT AJOUTE
            DO N1 = 1, 3
               IF( LARETE( 6, NOAR(N1) ) .EQ. 1 ) GOTO 101
            ENDDO

C           PERMUTATION POUR AMENER L'ARETE N1 EN ARETE 1 DE NTR
 101        IF( N1 .EQ. 2 ) THEN
               I                = NSTRIA( 1, NTR )
               NSTRIA( 1, NTR ) = NSTRIA( 2, NTR )
               NSTRIA( 2, NTR ) = NSTRIA( 3, NTR )
               NSTRIA( 3, NTR ) = I
               I         = NOAR( 1 )
               NOAR( 1 ) = NOAR( 2 )
               NOAR( 2 ) = NOAR( 3 )
               NOAR( 3 ) = I
            ELSE IF( N1 .EQ. 3 ) THEN
               I                = NSTRIA( 1, NTR )
               NSTRIA( 1, NTR ) = NSTRIA( 3, NTR )
               NSTRIA( 3, NTR ) = NSTRIA( 2, NTR )
               NSTRIA( 2, NTR ) = I
               I         = NOAR( 1 )
               NOAR( 1 ) = NOAR( 3 )
               NOAR( 3 ) = NOAR( 2 )
               NOAR( 2 ) = I
            ENDIF

C           LE POINT AJOUTE NS4 EST ENTRE NS1 ET NS2 PREMIERE ARETE DE NTR
            NS1 = NSTRIA( 1, NTR )
            NS2 = NSTRIA( 2, NTR )
            NS3 = NSTRIA( 3, NTR )
            NS4 = ABS( LARETE( 7, NOAR(1) ) )

C           MODIFICATION DU TRIANGLE NTR NS2->NS4
            NSTRIA( 2, NTR ) = NS4

C           AJOUT DU NOUVEAU TRIANGLE NBTRIA: NS4 NS2 NS3
            IF( NBTRIA .GE. MXTRIA ) GOTO 9993
            NBTRIA = NBTRIA + 1
            NSTRIA( 1, NBTRIA ) = NS4
            NSTRIA( 2, NBTRIA ) = NS2
            NSTRIA( 3, NBTRIA ) = NS3

            GOTO 900


         ELSE IF( NBPAT .EQ. 2 ) THEN

C           2 POINTS AJOUTES SUR LES 3 ARETES DU TRIANGLE NTR
C           =================================================
            N0 = 0
            DO N1 = 1, 3
               IF( LARETE( 6, NOAR(N1) ) .EQ. 0 ) THEN
                  N0 = N0 + 1
C                 EVENTUEL DERNIER NUMERO ARETE A ZERO POINT0 POINT
                  N2 = N1
               ELSE
C                 EVENTUEL NUMERO ARETE A 2 POINTS
                  N3 = N1
               ENDIF
            ENDDO

            IF( N0 .EQ. 1 ) THEN

C              1 POINT AJOUTE SUR 2 DES 3 ARETES DU TRIANGLE NTR
c              -------------------------------------------------
C              L'ARETE AVEC 0 POINT AJOUTE EST N2

C              PERMUTATION POUR AMENER L'ARETE N2 EN ARETE 1 DE NTR
               IF( N2 .EQ. 2 ) THEN
                  I                = NSTRIA( 1, NTR )
                  NSTRIA( 1, NTR ) = NSTRIA( 2, NTR )
                  NSTRIA( 2, NTR ) = NSTRIA( 3, NTR )
                  NSTRIA( 3, NTR ) = I
                  I         = NOAR( 1 )
                  NOAR( 1 ) = NOAR( 2 )
                  NOAR( 2 ) = NOAR( 3 )
                  NOAR( 3 ) = I

               ELSE IF( N2 .EQ. 3 ) THEN
                  I                = NSTRIA( 1, NTR )
                  NSTRIA( 1, NTR ) = NSTRIA( 3, NTR )
                  NSTRIA( 3, NTR ) = NSTRIA( 2, NTR )
                  NSTRIA( 2, NTR ) = I

                  I         = NOAR( 1 )
                  NOAR( 1 ) = NOAR( 3 )
                  NOAR( 3 ) = NOAR( 2 )
                  NOAR( 2 ) = I

               ENDIF

               NS1 = NSTRIA( 1, NTR )
               NS2 = NSTRIA( 2, NTR )
               NS3 = NSTRIA( 3, NTR )

C              LE POINT AJOUTE NS4 EST ENTRE NS2 ET NS3 SECONDE ARETE
               NS4 = ABS( LARETE( 7, NOAR(2) ) )

C              LE POINT AJOUTE NS5 EST ENTRE NS3 ET NS1 TROISIEME ARETE
               NS5 = ABS( LARETE( 7, NOAR(3) ) )

C              MODIFICATION DU TRIANGLE NTR NS3->NS5
               NSTRIA( 3, NTR ) = NS5

C              AJOUT DU NOUVEAU TRIANGLE NBTRIA: NS2 NS4 NS5
               IF( NBTRIA .GE. MXTRIA ) GOTO 9993
               NBTRIA = NBTRIA + 1
               NSTRIA( 1, NBTRIA ) = NS2
               NSTRIA( 2, NBTRIA ) = NS4
               NSTRIA( 3, NBTRIA ) = NS5

C              AJOUT DU NOUVEAU TRIANGLE NBTRIA: NS4 NS3 NS5
               IF( NBTRIA .GE. MXTRIA ) GOTO 9993
               NBTRIA = NBTRIA + 1
               NSTRIA( 1, NBTRIA ) = NS4
               NSTRIA( 2, NBTRIA ) = NS3
               NSTRIA( 3, NBTRIA ) = NS5

               GOTO 900

            ELSE

C              2 POINTS AJOUTES SUR 1 DES 3 ARETES DU TRIANGLE NTR
c              ---------------------------------------------------
C              L'ARETE AVEC 2 POINTS AJOUTES EST N3

C              PERMUTATION POUR AMENER L'ARETE N3 EN ARETE 1 DE NTR
               IF( N3 .EQ. 2 ) THEN
                  I                = NSTRIA( 1, NTR )
                  NSTRIA( 1, NTR ) = NSTRIA( 2, NTR )
                  NSTRIA( 2, NTR ) = NSTRIA( 3, NTR )
                  NSTRIA( 3, NTR ) = I
                  I         = NOAR( 1 )
                  NOAR( 1 ) = NOAR( 2 )
                  NOAR( 2 ) = NOAR( 3 )
                  NOAR( 3 ) = I

               ELSE IF( N3 .EQ. 3 ) THEN
                  I                = NSTRIA( 1, NTR )
                  NSTRIA( 1, NTR ) = NSTRIA( 3, NTR )
                  NSTRIA( 3, NTR ) = NSTRIA( 2, NTR )
                  NSTRIA( 2, NTR ) = I
                  I         = NOAR( 1 )
                  NOAR( 1 ) = NOAR( 3 )
                  NOAR( 3 ) = NOAR( 2 )
                  NOAR( 2 ) = I

               ENDIF

               NS1 = NSTRIA( 1, NTR )
               NS2 = NSTRIA( 2, NTR )
               NS3 = NSTRIA( 3, NTR )

C              LE DERNIER POINT AJOUTE NS5 EST ENTRE NS1 ET NS2 1-ERE ARETE
C              LE PREMIER POINT AJOUTE NS4 EST ENTRE NS1 ET NS2 1-ERE ARETE
               NS = LARETE( 7, NOAR(1) )
               IF( NS .LT. 0 ) THEN
                  NSENS = -1
                  NS    = -NS
               ELSE
                  NSENS = 1
               ENDIF
               IF( NS1 .LT. NS2 ) THEN
                  IF( NSENS .GT. 0 ) THEN
                     NS5 = NS
                     NS4 = NS5 - 1
                  ELSE
                     NS4 = NS
                     NS5 = NS4 - 1
                  ENDIF
               ELSE
                  IF( NSENS .GT. 0 ) THEN
                     NS4 = NS
                     NS5 = NS4 - 1
                  ELSE
                     NS5 = NS
                     NS4 = NS5 - 1
                  ENDIF
               ENDIF

C              MODIFICATION DU TRIANGLE NTR NS2->NS4 NS3->NS5
               NSTRIA( 2, NTR ) = NS4
               NSTRIA( 3, NTR ) = NS3

C              AJOUT DU NOUVEAU TRIANGLE NBTRIA
               IF( NBTRIA .GE. MXTRIA ) GOTO 9993
               NBTRIA = NBTRIA + 1
               NSTRIA( 1, NBTRIA ) = NS4
               NSTRIA( 2, NBTRIA ) = NS5
               NSTRIA( 3, NBTRIA ) = NS3

C              AJOUT DU NOUVEAU TRIANGLE NBTRIA
               IF( NBTRIA .GE. MXTRIA ) GOTO 9993
               NBTRIA = NBTRIA + 1
               NSTRIA( 1, NBTRIA ) = NS5
               NSTRIA( 2, NBTRIA ) = NS2
               NSTRIA( 3, NBTRIA ) = NS3

               GOTO 900

            ENDIF


         ELSE IF( NBPAT .EQ. 3 ) THEN

C           3 POINTS AJOUTES SUR LES 3 ARETES DU TRIANGLE NOTR
C           ==================================================
            N0 = 0
            DO N1 = 1, 3
               IF( LARETE( 6, NOAR(N1) ) .EQ. 0 ) THEN
                  N0 = N0 + 1
C                 EVENTUEL DERNIER NUMERO ARETE A ZERO POINT AJOUTE
                  N2 = N1
               ELSE
C                 EVENTUEL NUMERO ARETE A 3 POINTS AJOUTES
                  N3 = N1
               ENDIF
            ENDDO

            IF( N0 .EQ. 0 ) THEN

C              1 POINT AJOUTE SUR CHACUNE DES 3 ARETES DU TRIANGLE NTR
C              -------------------------------------------------------
               NS1 = NSTRIA( 1, NTR )
               NS2 = NSTRIA( 2, NTR )
               NS3 = NSTRIA( 3, NTR )

C              NS4 POINT AJOUTE ENTRE NS1 ET NS2 PREMIERE ARETE DE NTR
               NS4 = ABS( LARETE( 7, NOAR(1) ) )

C              NS5 POINT AJOUTE ENTRE NS2 ET NS3 SECONDE ARETE DE NTR
               NS5 = ABS( LARETE( 7, NOAR(2) ) )

C              NS6 POINT AJOUTE ENTRE NS2 ET NS3 3-EME ARETE DE NTR
               NS6 = ABS( LARETE( 7, NOAR(3) ) )

C              MODIFICATION DU TRIANGLE NTR
               NSTRIA( 2, NTR ) = NS4
               NSTRIA( 3, NTR ) = NS6

C              AJOUT DU NOUVEAU TRIANGLE NBTRIA
               IF( NBTRIA .GE. MXTRIA ) GOTO 9993
               NBTRIA = NBTRIA + 1
               NSTRIA( 1, NBTRIA ) = NS6
               NSTRIA( 2, NBTRIA ) = NS4
               NSTRIA( 3, NBTRIA ) = NS5

C              AJOUT DU NOUVEAU TRIANGLE NBTRIA
               IF( NBTRIA .GE. MXTRIA ) GOTO 9993
               NBTRIA = NBTRIA + 1
               NSTRIA( 1, NBTRIA ) = NS4
               NSTRIA( 2, NBTRIA ) = NS2
               NSTRIA( 3, NBTRIA ) = NS5

C              AJOUT DU NOUVEAU TRIANGLE NBTRIA
               IF( NBTRIA .GE. MXTRIA ) GOTO 9993
               NBTRIA = NBTRIA + 1
               NSTRIA( 1, NBTRIA ) = NS6
               NSTRIA( 2, NBTRIA ) = NS5
               NSTRIA( 3, NBTRIA ) = NS3

               GOTO 900


            ELSE IF( N0 .EQ. 2 ) THEN

C              3 POINTS AJOUTES SUR UNE ARETE DU TRIANGLE NTR
C              ----------------------------------------------
C              PERMUTATION POUR AMENER L'ARETE N3 EN ARETE 1 DE NTR
               IF( N3 .EQ. 2 ) THEN
                  I                = NSTRIA( 1, NTR )
                  NSTRIA( 1, NTR ) = NSTRIA( 2, NTR )
                  NSTRIA( 2, NTR ) = NSTRIA( 3, NTR )
                  NSTRIA( 3, NTR ) = I
                  I         = NOAR( 1 )
                  NOAR( 1 ) = NOAR( 2 )
                  NOAR( 2 ) = NOAR( 3 )
                  NOAR( 3 ) = I

               ELSE IF( N3 .EQ. 3 ) THEN
                  I                = NSTRIA( 1, NTR )
                  NSTRIA( 1, NTR ) = NSTRIA( 3, NTR )
                  NSTRIA( 3, NTR ) = NSTRIA( 2, NTR )
                  NSTRIA( 2, NTR ) = I
                  I         = NOAR( 1 )
                  NOAR( 1 ) = NOAR( 3 )
                  NOAR( 3 ) = NOAR( 2 )
                  NOAR( 2 ) = I

               ENDIF

               NS1 = NSTRIA( 1, NTR )
               NS2 = NSTRIA( 2, NTR )
               NS3 = NSTRIA( 3, NTR )

C              LE DERNIER POINT AJOUTE NS6 EST ENTRE NS1 ET NS2 1-ERE ARETE
C              LE PREMIER POINT AJOUTE NS4 EST ENTRE NS1 ET NS2 1-ERE ARETE
               NS = LARETE( 7, NOAR(1) )
               IF( NS .LT. 0 ) THEN
                  NSENS = -1
                  NS    = -NS
               ELSE
                  NSENS = 1
               ENDIF
               IF( NS1 .LT. NS2 ) THEN
                  IF( NSENS .GT. 0 ) THEN
                     NS6 = NS
                     NS5 = NS6 - 1
                     NS4 = NS5 - 1
                  ELSE
                     NS4 = NS
                     NS5 = NS4 - 1
                     NS6 = NS5 - 1
                  ENDIF
               ELSE
                  IF( NSENS .GT. 0 ) THEN
                     NS4 = NS
                     NS5 = NS4 - 1
                     NS6 = NS5 - 1
                  ELSE
                     NS6 = NS
                     NS5 = NS6 - 1
                     NS4 = NS5 - 1
                  ENDIF
               ENDIF

C              CONSTRUCTION DU BARYCENTRE DU TRIANGLE NTR
               IF( NBSOM .GE. MXSOM ) GOTO 9992
               NBSOM = NBSOM + 1
               DO K=1,3
                  XYZSOM(K,NBSOM) = ( XYZSOM(K,NS1)
     %                              + XYZSOM(K,NS2)
     %                              + XYZSOM(K,NS3) ) / 3
               ENDDO

C              MODIFICATION DU TRIANGLE NTR
               NSTRIA( 2, NTR ) = NS4
               NSTRIA( 3, NTR ) = NBSOM

C              AJOUT DU NOUVEAU TRIANGLE NBTRIA
               IF( NBTRIA .GE. MXTRIA ) GOTO 9993
               NBTRIA = NBTRIA + 1
               NSTRIA( 1, NBTRIA ) = NS4
               NSTRIA( 2, NBTRIA ) = NS5
               NSTRIA( 3, NBTRIA ) = NBSOM

C              AJOUT DU NOUVEAU TRIANGLE NBTRIA
               IF( NBTRIA .GE. MXTRIA ) GOTO 9993
               NBTRIA = NBTRIA + 1
               NSTRIA( 1, NBTRIA ) = NS5
               NSTRIA( 2, NBTRIA ) = NS6
               NSTRIA( 3, NBTRIA ) = NBSOM

C              AJOUT DU NOUVEAU TRIANGLE NBTRIA
               IF( NBTRIA .GE. MXTRIA ) GOTO 9993
               NBTRIA = NBTRIA + 1
               NSTRIA( 1, NBTRIA ) = NS6
               NSTRIA( 2, NBTRIA ) = NS2
               NSTRIA( 3, NBTRIA ) = NBSOM

C              AJOUT DU NOUVEAU TRIANGLE NBTRIA
               IF( NBTRIA .GE. MXTRIA ) GOTO 9993
               NBTRIA = NBTRIA + 1
               NSTRIA( 1, NBTRIA ) = NS2
               NSTRIA( 2, NBTRIA ) = NS3
               NSTRIA( 3, NBTRIA ) = NBSOM

C              AJOUT DU NOUVEAU TRIANGLE NBTRIA
               IF( NBTRIA .GE. MXTRIA ) GOTO 9993
               NBTRIA = NBTRIA + 1
               NSTRIA( 1, NBTRIA ) = NS3
               NSTRIA( 2, NBTRIA ) = NS1
               NSTRIA( 3, NBTRIA ) = NBSOM

               GOTO 900

            ELSE

C              3 POINTS AJOUTES
C              RECHERCHE DES 2 ARETES, UNE AVEC 1 POINT, L'AUTRE AVEC 2 POINTS
C              ---------------------------------------------------------------
C              L'ARETE N2 A ZERO POINT DEVIENT LA PREMIERE ARETE DE NTR
               IF( N2 .EQ. 2 ) THEN
                  I                = NSTRIA( 1, NTR )
                  NSTRIA( 1, NTR ) = NSTRIA( 2, NTR )
                  NSTRIA( 2, NTR ) = NSTRIA( 3, NTR )
                  NSTRIA( 3, NTR ) = I
                  I         = NOAR( 1 )
                  NOAR( 1 ) = NOAR( 2 )
                  NOAR( 2 ) = NOAR( 3 )
                  NOAR( 3 ) = I

               ELSE IF( N2 .EQ. 3 ) THEN
                  I                = NSTRIA( 1, NTR )
                  NSTRIA( 1, NTR ) = NSTRIA( 3, NTR )
                  NSTRIA( 3, NTR ) = NSTRIA( 2, NTR )
                  NSTRIA( 2, NTR ) = I
                  I         = NOAR( 1 )
                  NOAR( 1 ) = NOAR( 3 )
                  NOAR( 3 ) = NOAR( 2 )
                  NOAR( 2 ) = I

               ENDIF

               NS1 = NSTRIA( 1, NTR )
               NS2 = NSTRIA( 2, NTR )
               NS3 = NSTRIA( 3, NTR )

C              RECHERCHE DE L'ARETE A 1 POINT AJOUTE
               IF( LARETE( 6, NOAR(2) ) .EQ. 1 ) THEN

C                 L'ARETE 2 DE NTR A 1 POINT  AJOUTE
C                 L'ARETE 3 DE NTR A 2 POINTS AJOUTES
C                 -----------------------------------
                  NS4 = ABS( LARETE( 7, NOAR(2) ) )

C                 LE DERNIER POINT AJOUTE NS6 EST ENTRE NS3 ET NS1 3-EME ARETE
C                 LE PREMIER POINT AJOUTE NS5 EST ENTRE NS3 ET NS1 3-EME ARETE
                  NS = LARETE( 7, NOAR(3) )
                  IF( NS .LT. 0 ) THEN
                     NSENS = -1
                     NS    = -NS
                  ELSE
                     NSENS = 1
                  ENDIF
                  IF( NS3 .LT. NS1 ) THEN
                     IF( NSENS .GT. 0 ) THEN
                        NS6 = NS
                        NS5 = NS6 - 1
                     ELSE
                        NS5 = NS
                        NS6 = NS5 - 1
                     ENDIF
                  ELSE
                     IF( NSENS .GT. 0 ) THEN
                        NS5 = NS
                        NS6 = NS5 - 1
                     ELSE
                        NS6 = NS
                        NS5 = NS6 - 1
                     ENDIF
                  ENDIF

C                 MODIFICATION DU TRIANGLE NTR
                  NSTRIA( 3, NTR ) = NS6

C                 AJOUT DU NOUVEAU TRIANGLE NBTRIA
                  IF( NBTRIA .GE. MXTRIA ) GOTO 9993
                  NBTRIA = NBTRIA + 1
                  NSTRIA( 1, NBTRIA ) = NS6
                  NSTRIA( 2, NBTRIA ) = NS2
                  NSTRIA( 3, NBTRIA ) = NS4

C                 AJOUT DU NOUVEAU TRIANGLE NBTRIA
                  IF( NBTRIA .GE. MXTRIA ) GOTO 9993
                  NBTRIA = NBTRIA + 1
                  NSTRIA( 1, NBTRIA ) = NS6
                  NSTRIA( 2, NBTRIA ) = NS4
                  NSTRIA( 3, NBTRIA ) = NS5

C                 AJOUT DU NOUVEAU TRIANGLE NBTRIA
                  IF( NBTRIA .GE. MXTRIA ) GOTO 9993
                  NBTRIA = NBTRIA + 1
                  NSTRIA( 1, NBTRIA ) = NS5
                  NSTRIA( 2, NBTRIA ) = NS4
                  NSTRIA( 3, NBTRIA ) = NS3

               ELSE

C                 L'ARETE 3 DE NTR A 1 POINT  AJOUTE
C                 L'ARETE 2 DE NTR A 2 POINTS AJOUTES
C                 -----------------------------------
                  NS4 = ABS( LARETE( 7, NOAR(3) ) )

C                 LE DERNIER POINT AJOUTE NS6 EST ENTRE NS2 ET NS3 2-EME ARETE
C                 LE PREMIER POINT AJOUTE NS5 EST ENTRE NS2 ET NS3 2-EME ARETE
                  NS = LARETE( 7, NOAR(2) )
                  IF( NS .LT. 0 ) THEN
                     NSENS = -1
                     NS    = -NS
                  ELSE
                     NSENS = 1
                  ENDIF
                  IF( NS2 .LT. NS3 ) THEN
                     IF( NSENS .GT. 0 ) THEN
                        NS6 = NS
                        NS5 = NS6 - 1
                     ELSE
                        NS5 = NS
                        NS6 = NS5 - 1
                     ENDIF
                  ELSE
                     IF( NSENS .GT. 0 ) THEN
                        NS5 = NS
                        NS6 = NS5 - 1
                     ELSE
                        NS6 = NS
                        NS5 = NS6 - 1
                     ENDIF
                  ENDIF

C                 MODIFICATION DU TRIANGLE NTR
                  NSTRIA( 3, NTR ) = NS5

C                 AJOUT DU NOUVEAU TRIANGLE NBTRIA
                  IF( NBTRIA .GE. MXTRIA ) GOTO 9993
                  NBTRIA = NBTRIA + 1
                  NSTRIA( 1, NBTRIA ) = NS1
                  NSTRIA( 2, NBTRIA ) = NS5
                  NSTRIA( 3, NBTRIA ) = NS4

C                 AJOUT DU NOUVEAU TRIANGLE NBTRIA
                  IF( NBTRIA .GE. MXTRIA ) GOTO 9993
                  NBTRIA = NBTRIA + 1
                  NSTRIA( 1, NBTRIA ) = NS4
                  NSTRIA( 2, NBTRIA ) = NS5
                  NSTRIA( 3, NBTRIA ) = NS6

C                 AJOUT DU NOUVEAU TRIANGLE NBTRIA
                  IF( NBTRIA .GE. MXTRIA ) GOTO 9993
                  NBTRIA = NBTRIA + 1
                  NSTRIA( 1, NBTRIA ) = NS4
                  NSTRIA( 2, NBTRIA ) = NS6
                  NSTRIA( 3, NBTRIA ) = NS3

               ENDIF

               GOTO 900

            ENDIF

         ENDIF

C        PLUS DE 3 POINTS AJOUTES SUR LES ARETES DU TRIANGLE NTR
C        LE BARYCENTRE EST JOINT AUX 3 SOMMETS ET TOUS LES POINTS AJOUTES
C        CAS AUSSI DE 2 ARETES DE NTR SANS POINT AJOUTE
C        ================================================================
C        RANGEMENT DANS LE SENS DIRECT DU NUMERO DES SOMMETS ET POINTS AJOUTES
 800     NBS = 0
         DO N1 = 1, 3

C           LE PREMIER POINT EST LE PREMIER SOMMET DE L'ARETE N1
            NBS = NBS + 1
            NS1 = NSTRIA( N1, NTR )
            NOPTAR( NBS ) = NS1
            IF( N1 .EQ. 3 ) THEN
               N = 1
            ELSE
               N = N1+1
            ENDIF
            NS2 = NSTRIA( N, NTR )

C           LE DERNIER POINT AJOUTE DE L'ARETE N1
            N3 = LARETE( 7, NOAR(N1) )
            IF( N3 .LT. 0 ) THEN
               NSENS = -1
               N3    = -N3
            ELSE
               NSENS = 1
            ENDIF

C           LE NUMERO - 1 DU PREMIER POINT DE L'ARETE N1
            NBPAJ = LARETE( 6, NOAR(N1) )
            N2    = N3 - NBPAJ

C           L'ORDRE DOIT SUIVRE LE SENS DIRECT DES ARETES DU TRIANGLE
            IF( NS1 .LT. NS2 ) THEN
               IF( NSENS .GT. 0 ) THEN
                  DO K = 1, NBPAJ
                     NBS = NBS + 1
                     NOPTAR( NBS ) = N2 + K
                  ENDDO
               ELSE
                  DO K = 1, NBPAJ
                     NBS = NBS + 1
                     NOPTAR( NBS ) = N3 - K + 1
                  ENDDO
               ENDIF
            ELSE
               IF( NSENS .GT. 0 ) THEN
                  DO K = 1, NBPAJ
                     NBS = NBS + 1
                     NOPTAR( NBS ) = N3 - K + 1
                  ENDDO
               ELSE
                  DO K = 1, NBPAJ
                     NBS = NBS + 1
                     NOPTAR( NBS ) = N2 + K
                  ENDDO
               ENDIF
            ENDIF

         ENDDO

C        FERMETURE DU NUMERO DES POINTS AJOUTES SUR LES 3 ARETES DE NTR
         NOPTAR( NBS+1 ) = NSTRIA( 1, NTR )

C        NUMERO NOPTAR DU DERNIER POINT SUR L'ARETE I DU TRIANGLE NTR
         N = 0
         DO N1=1,3
            N = N + 1 + LARETE( 6, NOAR(N1) )
            NDPTNTR( N1 ) = N
         ENDDO

C        TRIANGULATION DU TRIANGLE NTR CONSIDERE COMME UN
C        CONTOUR FERME PAR LES SOMMETS DE NOPTAR(1:NBS+1)
C        ------------------------------------------------
         CALL TRIA1TRIA( NOFOTI, MXSOM,  NBSOM,  XYZSOM,
     %                   NTR,    MXTRIA, NBTRIA, NSTRIA,
     %                   NOPTAR, NDPTNTR,  IERR )
         IF( IERR .NE. 0 ) GOTO 9999


C        TOUS LES EF CREES SONT DES TRIANGLES DE 4-EME SOMMET NUL
C        --------------------------------------------------------
 900     NSTRIA( 4, NTR ) = 0
         DO K = NBTRI0+1, NBTRIA
            NSTRIA( 4, K ) = 0
         ENDDO

C        TRACE DES SOUS-TRIANGLES DU TRIANGLE NTR
         IF( NBTRAT+1+NBTRIA-NBTRI0 .LE. MXTRAT ) THEN
            NBTRAT0 = NBTRAT
            NBTRAT  = NBTRAT + 1
            LITRAT( NBTRAT ) = NTR
            DO K = NBTRI0+1, NBTRIA
               NBTRAT = NBTRAT + 1
               LITRAT( NBTRAT ) = K
            ENDDO
C           TRACE DES TRIANGLES NBTRAT0+1 A NBTRAT
            CALL TRTRIAN( 'a1lgarti', XYZSOM, 4, MXTRAT-NBTRAT0,
     %                     NBTRAT-NBTRAT0, LITRAT(NBTRAT0+1), NSTRIA )
         ENDIF

 1000 ENDDO

 1010 PRINT*,'a1lgarti:',NBSOM, ' SOMMETS   dont',NBSOM-NBSOM0, 
     %       ' SOMMETS AJOUTES dans les TRIANGLES'
      PRINT*,'a1lgarti:',NBTRIA,' TRIANGLES dont',NBTRIA-NBTRIA0,
     %' TRIANGLES AJOUTES'

C     TRACE DU TOTAL DES NBTRAT SOUS-TRIANGLES
      CALL TRTRIAN( 'a1lgarti', XYZSOM, 4,MXTRAT,NBTRAT,LITRAT, NSTRIA )

C     PAS D'ERREUR DETECTEE
      IERR = 0
      GOTO 9999


cccC     PROBLEME DE CALCUL DE LA FONCTION TAILLE_IDEALE(x,y,z)
ccc 9991 PRINT *,'a1lgarti: CALCUL INCORRECT de TAILLE_IDEALE() AU SOMMET',
ccc     %         (XYZ(L,N),L=1,3),' TAILLE_IDEALE=',TI
ccc      IERR = 1
ccc      GOTO 9999

C     TABLEAU XYZSOM SATURE
 9992 IERR = 2
      PRINT *,'a1lgarti: TABLEAU XYZSOM SATURE MXSOM=',MXSOM
      GOTO 9995

C     TABLEAU NSTRIA SATURE
 9993 IERR = 3
      PRINT*,'a1lgarti: TABLEAU NSTRIA SATURE MXTRIA=',MXTRIA
      
C     AVERTISSEMENT
 9995 PRINT *,'a1lgarti: REVOIR la TAILLE d''ARETE PAR DEFAUT',DARETE,
     %' ou la fonction TAILLE_IDEALE()'

 9999 TRACTE = TRACTE0
      RETURN
      END
