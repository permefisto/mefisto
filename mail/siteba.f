      SUBROUTINE SITEBA( IAR0,   NOTET0, XYZSOM,
     %                   NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                   N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                   NBTEDS, MXTETR, NOTEDS,
     %                   QUALIT, VOLET0, QUAET0, VOLET1, QUAET1,
     %                   P,      NESSAI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SIMULATION DE LA TETRAEDRISATION ETOILEE PAR LE BARYCENTRE
C -----    DU TETRAEDRE NOTET0 A PARTIR D'UNE ARETE DIAGONALE

C ENTREES:
C --------
C IAR0   : NUMERO DANS LE TETRAEDRE NOTET0 DE L'ARETE-DIAGONALE A TRAITER
C NOTET0 : NUMERO DANS NSTETR DU TETRAEDRE CONTENANT LE POINT P
C XYZSOM : LES 3 COORDONNEES DES SOMMETS
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NSTETR(>3)
C NSTETR : NUMERO DES 4 SOMMETS DE CHAQUE TETRAEDRE
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DANS NOTESO DE CHAQUE SOMMET
C NOTESO : NOTESO(1,*) NUMERO DANS NSTETR DU TETRAEDRE
C          NOTESO(2,*) NUMERO DANS NOTESO DU TETRAEDRE SUIVANT
C                      0 SI C'EST LE DERNIER

C MODIFIES :
C ----------
C N1FEVI : NUMERO DANS NFETOI DE LA PREMIERE FACE VIDE
C N1FEOC : NUMERO DANS NFETOI DE LA PREMIERE FACE OCCUPEE
C NFETOI : NUMERO DES 3 SOMMETS, TETRAEDRE ET SUIVANTE DES FACES
C          VUES UNE FOIS DE L'ETOILE CHAINEES N1FEOC PUIS NFETOI(5,*)
C VOETOI : VOLUME  DU TETRAEDRE INITIAL ASSOCIE A LA FACE
C QUETOI : QUALITE DU TETRAEDRE INITIAL ASSOCIE A LA FACE
C QUALIT : QUALITE DES TETRAEDRES ACTUELS

C SORTIES:
C --------
C NBTEDS : NOMBRE DE TETRAEDRES A DETRUIRE
C NOTEDS : NUMERO DANS NSTETR DES TETRAEDRES A DETRUIRE
C QUAET0 : QUALITE INITIALE DE LA TETRAEDRISATION ETOILEE PAR P
C QUAET1 : QUALITE FINALE   DE LA TETRAEDRISATION ETOILEE PAR P
C          0 INDIQUE UNE IMPOSSIBLITE D'OBTENIR UNE ETOILE POUR P
C P      : 3 COORDONNEES DU POINT BARYCENTRE
C NESSAI : =0 PAS    DE    METHODE AUGMENTANT LA QUALITE DE L'ETOILE
C          >0 NUMERO DE LA METHODE AUGMENTANT LA QUALITE DE L'ETOILE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY  Fevrier 2016
C....................................................................012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
C     STOPTE = FAUX ==> PAS D'ARRET APRES LE TRACE DE CHAQUE TETRAEDRE
C              VRAI ==> DEMANDE D'UN CARACTERE POUR REDEMARRER
C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C              VRAI ==> TRACE DES TETRAEDRES

      INTEGER           NSTETR(NBSOTE,*),
     %                  NO1TSO(*),
     %                  NOTESO(2,*),
     %                  NFETOI(5,*),
     %                  NOTEDS(MXTETR)

      REAL              P(3), XYZSOM(3,*),
     %                  VOETOI(*), QUETOI(*), QUALIT(*),
     %                  ARMIN, ARMAX, SURFTR(4), VTE, VOLET0, VOLET1

ccc      CHARACTER*80      KTITRE

C                       NO SOMMETS D'UNE ARETE D'UN TETRAEDRE
      INTEGER           NOSOARTE(2,6)
      DATA              NOSOARTE / 1,2, 2,3, 3,1, 4,1, 4,2, 4,3 /
C                       NO ARETE OPPOSEE DANS UN TETRAEDRE
      INTEGER           NOAROPTE(6)
      DATA              NOAROPTE/  6, 4, 5, 2, 3, 1 /

      GRAND  = RINFO( 'GRAND' )
      QUAET0 = 0
      QUAET1 = 0
      NESSAI = 3
      IAR    = IAR0

C     LES 2 SOMMETS EXTREMITES DE CETTE ARETE IAR
C     -------------------------------------------
 3    NS1 = NSTETR( NOSOARTE(1,IAR), NOTET0 )
      IF( NS1 .LE. 0 ) GOTO 9000

      NS2 = NSTETR( NOSOARTE(2,IAR), NOTET0 )
      IF( NS2 .LE. 0 ) GOTO 9000

C     RECHERCHE DES TETRAEDRES D'ARETE NS1-NS2
C     ----------------------------------------
      CALL NOTEAR( NS1,    NS2,
     %             NBSOTE, NSTETR, NO1TSO, NOTESO,
     %             NBTED1, NOTEDS )
      IF( NBTED1 .GT. MXTETR ) GOTO 9000

C     RECHERCHE DE LA SECONDE DIAGONALE
      IAR2 = NOAROPTE( IAR )

C     LES 2 SOMMETS EXTREMITES DE CETTE ARETE IAR2
C     --------------------------------------------
      NS3 = NSTETR( NOSOARTE(1,IAR2), NOTET0 )
      IF( NS3 .LE. 0 ) GOTO 9000

      NS4 = NSTETR( NOSOARTE(2,IAR2), NOTET0 )
      IF( NS4 .LE. 0 ) GOTO 9000

C     RECHERCHE DES TETRAEDRES D'ARETE NS3-NS4
      CALL NOTEAR( NS3,    NS4,
     %             NBSOTE, NSTETR, NO1TSO, NOTESO,
     %             NBTED2, NOTEDS(NBTED1+1) )
      IF( NBTED1+NBTED2 .GT. MXTETR ) GOTO 9000

C     SUPPRESSION DES DOUBLONS DANS NOTEDS
      NBTEDS = NBTED1
      DO 18 J=NBTED1+1, NBTED1+NBTED2
         NT = NOTEDS( J )
         DO I=1,NBTEDS
            IF( NOTEDS( I ) .EQ. NT ) GOTO 18
         ENDDO
         NBTEDS = NBTEDS + 1
         NOTEDS( NBTEDS ) = NT
 18   ENDDO

C     CALCUL DU VOLUME ET QUALITE DE CES NBTEDS TETRAEDRES
C     ----------------------------------------------------
 20   QUAET0 = GRAND
      VOLET0 = 0
      DO J=1,NBTEDS

         NT = NOTEDS(J)
         CALL QUATET( XYZSOM(1, NSTETR(1,NT) ),
     %                XYZSOM(1, NSTETR(2,NT) ),
     %                XYZSOM(1, NSTETR(3,NT) ),
     %                XYZSOM(1, NSTETR(4,NT) ),
     %                ARMIN, ARMAX, SURFTR, VTE, QTE )

         QUALIT(NT) = QTE
         QUAET0 = MIN( QUAET0, QTE )

         IF( VTE .LT. 0.0 ) THEN
            PRINT *,'siteba: PB nstetr(',nt,')=',(nstetr(k,nt),k=1,4),
     %              ' VOLUME0=',VTE,'<0   QUALITE0=',QTE
            VTE = -VTE
         ENDIF
         VOLET0 = VOLET0 + VTE

      ENDDO

      IF( NESSAI .EQ. 3 .OR. NESSAI .EQ. 5 ) THEN

C        XYZ DU BARYCENTRE DU TETRAEDRE NOTET0
C        -------------------------------------
         DO I=1,3
            P(I) = ( XYZSOM( I, NSTETR(1,NOTET0) )
     %             + XYZSOM( I, NSTETR(2,NOTET0) )
     %             + XYZSOM( I, NSTETR(3,NOTET0) )
     %             + XYZSOM( I, NSTETR(4,NOTET0) ) ) / 4
         ENDDO

      ELSE

C        XYZ DU BARYCENTRE DES BARYCENTRES DES NBTED1 TETRAEDRES NOTEDS
C        --------------------------------------------------------------
         DO I=1,3
            P(I) = 0
         ENDDO

         DO J=1,NBTED1
            NT = NOTEDS(J)
            DO I=1,3
               P(I) = P(I) + ( XYZSOM( I, NSTETR(1,NT) )
     %                       + XYZSOM( I, NSTETR(2,NT) )
     %                       + XYZSOM( I, NSTETR(3,NT) )
     %                       + XYZSOM( I, NSTETR(4,NT) ) ) / 4
            ENDDO
         ENDDO

         DO I=1,3
            P(I) = P(I) / NBTED1
         ENDDO

      ENDIF

C     REINITIALISATION A VIDE DES FACES DE NFETOI
      CALL ETOIOCVI( N1FEOC, N1FEVI, NFETOI )

C     L'ETOILE DES FACES SIMPLES DES TETRAEDRES D'ARETE LES 2 DIAGONALES
C     ------------------------------------------------------------------
      DO J=1,NBTEDS

C        TETRAEDRE J AJOUTE A L'ETOILE
         NT = NOTEDS( J )

         DO I=1,4
C           SI    ( LA FACE I DU TETRAEDRE NOTET0 N'APPARTIENT PAS
C                   AUX FACES DE L'ETOILE NFETOI )
C           ALORS ELLE EST AJOUTEE A L'ETOILE DANS NFETOI
C           SINON ELLE EST EMPILEE DANS NPILE POUR ETRE DETRUITE ENSUITE
C                 ELLE EST SUPPRIMEE DE L'ETOILE NFETOI
            CALL AJFACE( 1, NT, I, NBSOTE, NSTETR,
     %                   N1FEOC, N1FEVI, NFETOI,
     %                   NF  )
            IF( NF .LT. 0 ) THEN
C              SATURATION DES FACES DE L'ETOILE
               GOTO 9000
            ENDIF
         ENDDO

      ENDDO

C     ICI, L'ETOILE EST COMPLETE. CALCUL DE SON VOLUME ET SA QUALITE
C     --------------------------------------------------------------
      VOLET1 = 0
      QUAET1 = GRAND
      NF     = N1FEOC

 30   IF( NF .GT. 0 ) THEN

C        CALCUL DU VOLUME ET QUALITE DU TETRAEDRE DU A LA FACE+P
         CALL QUATET( XYZSOM(1,NFETOI(1,NF)),
     %                XYZSOM(1,NFETOI(2,NF)),
     %                XYZSOM(1,NFETOI(3,NF)),
     %                P,
     %                ARMIN, ARMAX, SURFTR, VTE, QTE )

         IF( VTE .LT. 0 ) THEN

ccc            PRINT*,'siteba: POINT=',P,' NESSAI=',NESSAI
ccc            PRINT*,'siteba: NFETOI(',NF,')=',(NFETOI(K,NF),K=1,5)
ccc            PRINT*,'siteba: VOLUME(',NF,')=',VTE, '<0',
ccc     %                   ' QUALITE(',NF,')=',QTE
 
ccc            TRACTE = .TRUE.
ccc            KTITRE = 'SITEBA: TETRAEDRE FACE+POINT DE VOLUME<0'
ccc            CALL TRETOPT( KTITRE, XYZSOM, P, N1FEOC, NFETOI )

cccC           ESSAI DE REORIENTATION DE LA FACE NFETOI(NF)
ccc            NF1          = NFETOI(2,NF)
ccc            NFETOI(2,NF) = NFETOI(3,NF)
ccc            NFETOI(3,NF) = NF1
ccc            GOTO 35

            GOTO 50

         ENDIF

C        VOLUME DE L'ETOILE DE CENTRE P
         VOETOI(NF) = VTE
         VOLET1 = VOLET1 + ABS( VTE )

C        QUALITE DE L'ETOILE DE CENTRE P
         QUETOI(NF) = QTE
         IF( QTE .LT. QUAET1 ) THEN
            QUAET1 = QTE
         ENDIF

C        PASSAGE A LA FACE SUIVANTE
         NF = NFETOI(5,NF)
         GOTO 30

      ENDIF

C     BILAN SUR LES VOLUMES ET QUALITES DE L'ETOILE AVANT ET APRES
C     ------------------------------------------------------------
      IF( ABS( VOLET1-VOLET0 ) .LE. 0.001 * VOLET0 .AND.
     %    QUAET1 .NE. GRAND .AND. QUAET1 .GT. QUAET0 ) THEN
         GOTO 9999
      ENDIF

 50   IF( NESSAI .EQ. 3 ) THEN

C        ESSAI 4 AVEC LE BARYCENTRE DES NBDES1 TETRAEDRES
C        S'ENROULANT AUTOUR DE LA PREMIERE DIAGONALE
         NESSAI = 4
         NBTEDS = NBTED1
         GOTO 20

      ELSE IF( NESSAI .EQ. 4 ) THEN

C        ESSAI 5 AVEC L'AUTRE DIAGONALE
         NESSAI = 5
         IAR    = IAR2
         GOTO 3

      ELSE IF( NESSAI .EQ. 5 ) THEN

C        ESSAI 6 AVEC L'AUTRE DIAGONALE ET
C        LE BARYCENTRE DES NBDES1 TETRAEDRES
C        S'ENROULANT AUTOUR DE LA PREMIERE DIAGONALE
         NESSAI = 6
         NBTEDS = NBTED1
         GOTO 20

      ELSE

         GOTO 9000

      ENDIF

C     PROBLEME RENCONTRE OU VOLUMES DIFFERENTS OU QUALITE INFERIEURE
C     --------------------------------------------------------------
 9000 NBTEDS = 0
      QUAET1 = 0
      NESSAI = 0

 9999 RETURN
      END
