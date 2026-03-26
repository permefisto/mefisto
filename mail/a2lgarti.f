      SUBROUTINE A2LGARTI( NOFOTI, MXSOM,  NBSOM,  XYZSOM,
     %                     L1ARET, L2ARET, LARETE,
     %                     MXTRIA, NBTRIA, NOTRIA,  IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TOUTE ARETE TROP LONGUE PAR RAPPORT AUX VALEURS DE 
C -----    de la FONCTION TAILLE_IDEALE ou DARETE VALEUR PAR DEFAUT
C          AUX EXTREMITES IL EST AJOUTE UN POINT MILIEU JOINT AUX
C          2 SOMMETS OPPOSES DES TRIANGLES QUI LUI SONT ADJACENTS

C ENTREES:
C --------
C NOFOTI : NUMERO DE LA FONCTION 'TAILLE_IDEALE' DES ARETES SINON 0
C MXSOM  : NOMBRE MAXIMAL DE SOMMETS DE LA TRIANGULATION
C NBSOM  : NOMBRE DE SOMMETS DE LA TRIANGULATION
C XYZSOM : X Y Z LES 3 COORDONNEES DES SOMMETS DE LA TRIANGULATION
C L1ARET : NOMBRE DE MOTS PAR ARETE DU TABLEAU NARET
C L2ARET : NOMBRE DE ARETES DU TABLEAU NARET
C MXTRIA : NOMBRE MAXIMAL DE TRIANGLES DECLARABLES DANS NOTRIA
C NBTRIA : NOMBRE DE TRIANGLES DE LA TRIANGULATION

C MODIFIES :
C ----------
C NBSOM  : NOMBRE DE SOMMETS DE LA TRIANGULATION DANS XYZSOM
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
C NOTRIA : NUMERO XYZSOM DES 3 SOMMETS ET 0 POUR CHACUN DES NBTRIA TRIANGLES
C NBTRIA : NOMBRE DE TRIANGLES DE LA TRIANGULATION

C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR DETECTEE
C          1 SI TAILLE_IDEALE NON CALCULABLE EN UN SOMMET
C          2 SI TABLEAU XYZSOM SATURE (AUGMENTER MXSOM)
C          3 SI TABLEAU NOTRIA SATURE (AUGMENTER MXTRIA)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC PARIS & VEULETTES SUR MER Avril 2017
C....................................................................012
      REAL              XYZSOM(3,MXSOM), LGARET, LGAMAX
      INTEGER           NOTRIA(4,NBTRIA)
      INTEGER           LARETE(L1ARET,L2ARET), NOSOAR(2)
      DOUBLE PRECISION  D, XYZ(3,2), TI(2)

C  ---------------------------------------------------------------------
C     AJOUTER UN MILIEU A LA PLUS LONGUE ARETE, JOINT AU SOMMET OPPOSE
C     DE TOUS LES TRIANGLES ADJACENTS PAR CETTE ARETE (1 ou 2 TRIANGLES)
C  ---------------------------------------------------------------------
      LIBREF  = L2ARET
      NBSOM0  = NBSOM
      NBTRIA0 = NBTRIA

      DO 10 NTR = 1, NBTRIA

C        CALCUL DE LA LONGUEUR DES 3 ARETES DU TRIANGLE NTR
         NMAX   = 0
         LGAMAX = 0.0
         NS1    = NOTRIA(3,NTR)
         DO N = 1, 3
            NS2 = NOTRIA( N, NTR )
            LGARET = ( XYZSOM(1,NS2) - XYZSOM(1,NS1) ) ** 2
     %             + ( XYZSOM(2,NS2) - XYZSOM(2,NS1) ) ** 2
     %             + ( XYZSOM(3,NS2) - XYZSOM(3,NS1) ) ** 2
            IF( LGARET .GT. LGAMAX ) THEN
               NMAX   = N
               LGAMAX = LGARET
            ENDIF
C           PASSAGE A L'ARETE SUIVANTE
            NS1 = NS2
         ENDDO
         LGAMAX = SQRT(LGAMAX)

C        L'ARETE DE LONGUEUR MAX EST NMAX DE SOMMETS NS1AM NS2AM
C        CALCUL DE LA TAILLE_IDEALE AUX 2 SOMMETS
         NS2AM = NOTRIA( NMAX, NTR )
         IF( NMAX .EQ. 1 ) THEN
            NS1AM = 3
         ELSE
            NS1AM = NMAX - 1
         ENDIF
         IF( NS1AM .EQ. 1 ) THEN
            NS3AM = 3
         ELSE
            NS3AM = NS1AM - 1
         ENDIF
         NS1AM = NOTRIA( NS1AM, NTR )
         NS3AM = NOTRIA( NS3AM, NTR )

         DO K = 1, 3
            XYZ( K, 1 ) = XYZSOM( K, NS1AM )
            XYZ( K, 2 ) = XYZSOM( K, NS2AM )
         ENDDO

C        CALCUL DE TAILLE_IDEALE(X,Y,Z) AUX 2 SOMMETS
         CALL TAILIDEA( NOFOTI, XYZ(1,1),  NCODEV, TI(1) )
         IF( NCODEV .EQ. 0 ) THEN
            N = 1
            GOTO 9991
         ENDIF

         CALL TAILIDEA( NOFOTI, XYZ(1,2),  NCODEV, TI(2) )
         IF( NCODEV .EQ. 0 ) THEN
            N = 2
            GOTO 9991
         ENDIF

         IF( LGAMAX .GE. 1.366D0*MIN( TI(1), TI(2) ) ) THEN

C           SI CETTE TROP GRANDE ARETE EST COMMUNE A 2 TRIANGLES
C           ALORS DECOUPAGE EN 2 DE CHAQUE TRIANGLE EN JOIGNANT
C                 LE MILIEU DE L'ARETE AU SOMMET OPPOSE
            IF( NS1AM .LT. NS2AM ) THEN
               NOSOAR(1) = NS1AM
               NOSOAR(2) = NS2AM
            ELSE
               NOSOAR(1) = NS2AM
               NOSOAR(2) = NS1AM
            ENDIF

C           NAR LE NUMERO DE L'ARETE NS1 NS2 DANS LARETE
            CALL HACHAG( 2, NOSOAR, L1ARET, L2ARET, LARETE, 3,
     %                   LIBREF, NAR )

            IF( NAR .LE. 0 ) THEN
C              ARETE AJOUTEE OU NON TROUVEE => PAS D'AJOUT
ccc              PRINT*,'a2lgarti: ARETE',NS1,NS2,' du TRIANGLE(',NTR,')=',
ccc     %        (NOTRIA(L,NTR),L=1,3),' NON ARETE de la TRIANGULATION'
               GOTO 10
            ENDIF

C           NTR0 EST LE TRIANGLE OPPOSE AU TRIANGLE NTR PAR L'ARETE NAR
            IF( ABS( LARETE(4,NAR) ) .EQ. NTR ) THEN
               NTR0 = ABS( LARETE(5,NAR) )
            ELSE IF(  ABS( LARETE(5,NAR) ) .EQ. NTR ) THEN
               NTR0 = ABS( LARETE(4,NAR) )
            ELSE
C              TRIANGLE FRONTALIER PAS D'AJOUT
               NTR0 = 0
               GOTO 10
            ENDIF

C           VERIFICATION QUE NS1AM et NS2AM SONT 2 SOMMETS DE NTR0
            DO N = 1, 3
               IF( NOTRIA( N, NTR0 ) .EQ. NS1AM ) THEN
                  DO N1 = 1, 3
                     IF( NOTRIA( N1, NTR0 ) .EQ. NS2AM ) GOTO 3
                  ENDDO
C                 NTR0 A POUR SOMMET NS1AM MAIS PAS NS2AM
                  GOTO 10
               ENDIF
            ENDDO
C           NTR0 N'A PAS POUR SOMMET NS1AM
            GOTO 10

C           AJOUT EFFECTIF DU MILIEU DE L'ARETE NAR DE NTR ET NTR0
C           ------------------------------------------------------
C           RECHERCHE DU SOMMET NON NS1AM NON NS2AM DANS NTR0
 3          DO N = 1, 3
               NS4 = NOTRIA( N, NTR0 )
               IF( NS4 .NE. NS1AM .AND. NS4 .NE. NS2AM ) GOTO 5
            ENDDO

C           CONSTRUCTION DU "MILIEU" DE L'ARETE AU PRORATA DE LA
C           TAILLE IDEALE AUX 2 SOMMETS DE L'ARETE MAX
 5          IF( NBSOM .GE. MXSOM ) GOTO 9992
            NBSOM = NBSOM + 1

C           POIDS = PERMUTE DE LA TAILLE_IDEALE
            D = TI(1) + TI(2)
            DO K=1,3
               XYZSOM(K,NBSOM)= REAL( ( TI(2) * XYZ(K,1)
     %                                + TI(1) * XYZ(K,2) ) / D )
            ENDDO

C           LE TRIANGLE NTR EST DECOUPE EN 2: NTR + NBTRIA+1
            NOTRIA( 1, NTR ) = NS1AM
            NOTRIA( 2, NTR ) = NBSOM
            NOTRIA( 3, NTR ) = NS3AM
            NOTRIA( 4, NTR ) = 0

C           AJOUT DU NOUVEAU TRIANGLE NBTRIA
            IF( NBTRIA .GE. MXTRIA ) GOTO 9993
            NBTRIA = NBTRIA + 1
            NOTRIA( 1, NBTRIA ) = NS3AM
            NOTRIA( 2, NBTRIA ) = NBSOM
            NOTRIA( 3, NBTRIA ) = NS2AM
            NOTRIA( 4, NBTRIA ) = 0

C           LE TRIANGLE NTR0 EST DECOUPE EN 2: NTR0 + NBTRIA+1
            NOTRIA( 1, NTR0 ) = NS1AM
            NOTRIA( 2, NTR0 ) = NS4
            NOTRIA( 3, NTR0 ) = NBSOM
            NOTRIA( 4, NTR0 ) = 0

C           AJOUT DU NOUVEAU TRIANGLE NBTRIA
            IF( NBTRIA .GE. MXTRIA ) GOTO 9993
            NBTRIA = NBTRIA + 1
            NOTRIA( 1, NBTRIA ) = NBSOM
            NOTRIA( 2, NBTRIA ) = NS4
            NOTRIA( 3, NBTRIA ) = NS2AM
            NOTRIA( 4, NBTRIA ) = 0

         ENDIF

 10   ENDDO

      PRINT*,'a2lgarti:',NBSOM, ' SOMMETS  '
     %      ,NBSOM -NBSOM0 ,' SOMMETS   AJOUTES'
      PRINT*,'a2lgarti:',NBTRIA,' TRIANGLES'
     %      ,NBTRIA-NBTRIA0,' TRIANGLES AJOUTES'
      PRINT*

      GOTO 9999

C     PROBLEME
 9991 PRINT *,'a2lgarti: CALCUL INCORRECT de TAILLE_IDEALE() AU
     % SOMMET', (XYZ(L,N),L=1,3),' TAILLE_IDEALE=',TI
      IERR = 1
      GOTO 9999

C     TABLEAU XYZSOM SATURE
 9992 PRINT *,'a2lgarti: TABLEAU XYZSOM SATURE MXSOM=',MXSOM
      IERR = 2
      GOTO 9999

C     TABLEAU NOTRIA SATURE
 9993 PRINT*,'a2lgarti: TABLEAU NOTRIA SATURE MXTRIA=',MXTRIA
      IERR = 3

 9999 RETURN
      END
