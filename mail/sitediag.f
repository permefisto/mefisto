      SUBROUTINE SITEDIAG( NAR0,   NOTET0, XYZSOM,
     %                     NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                     NBTEDS, MXTEDS, NOTEDS,
     %                     QUALIT, VOLET0, QUAET0, VOLET1, QUAET1,
     %                     NESSAI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SIMULATION DE LA TETRAEDRISATION D'UNE ETOILE LORS DE
C -----    L'ECHANGE DES DIAGONALES DU TETRAEDRE NOTET0

C ENTREES:
C --------
C NAR0   : NUMERO DANS LE TETRAEDRE NOTET0 DE L'ARETE-DIAGONALE A TRAITER
C NOTET0 : NUMERO DANS NSTETR DU TETRAEDRE CONTENANT LE POINT P
C XYZSOM : LES 3 COORDONNEES DES SOMMETS
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NSTETR(>3)
C NSTETR : NUMERO DES 4 SOMMETS DE CHAQUE TETRAEDRE
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DANS NOTESO DE CHAQUE SOMMET
C NOTESO : NOTESO(1,*) NUMERO DANS NSTETR DU TETRAEDRE
C          NOTESO(2,*) NUMERO DANS NOTESO DU TETRAEDRE SUIVANT
C                      0 SI C'EST LE DERNIER

C SORTIES:
C --------
C NBTEDS : NOMBRE DE TETRAEDRES A DETRUIRE
C NOTEDS : NUMERO DANS NSTETR DES TETRAEDRES A DETRUIRE
C QUALIT : TABLEAU DE LA QUALITE DE CHAQUE TETRAEDRE
C QUAET0 : QUALITE INITIALE DE LA TETRAEDRISATION ETOILEE PAR P
C QUAET1 : QUALITE FINALE   DE LA TETRAEDRISATION ETOILEE PAR P
C          0 INDIQUE UNE IMPOSSIBLITE D'OBTENIR UNE ETOILE POUR P
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
     %                  NOTEDS(MXTEDS),
     %                  NOSOTE(4)
      REAL              QUALIT(*),
     %                  XYZSOM(3,*),
     %                  ARMIN, ARMAX, SURFTR(4),
     %                  QTE0, QTE, QUAET0, QUAET1,
     %                  VTE0, VTE, VOLET0, VOLET1

C                       NO SOMMETS D'UNE ARETE D'UN TETRAEDRE
      INTEGER           NOSOARTE(2,6)
      DATA              NOSOARTE / 1,2, 2,3, 3,1, 4,1, 4,2, 4,3 /
C                       NO ARETE OPPOSEE DANS UN TETRAEDRE
      INTEGER           NOAROPTE(6)
      DATA              NOAROPTE/  6, 4, 5, 2, 3, 1 /

      QUAET0 = 0
      QUAET1 = 0
      NESSAI = 1
      NAR    = NAR0

C     LES 2 SOMMETS EXTREMITES DE CETTE ARETE NAR DIAGONALE 1
C     -------------------------------------------------------
 1    NS1 = NSTETR( NOSOARTE(1,NAR), NOTET0 )
      IF( NS1 .LE. 0 ) GOTO 9000

      NS2 = NSTETR( NOSOARTE(2,NAR), NOTET0 )
      IF( NS2 .LE. 0 ) GOTO 9000

C     RECHERCHE DES TETRAEDRES D'ARETE NS1-NS2
C     ----------------------------------------
      CALL NOTEAR( NS1,    NS2,
     %             NBSOTE, NSTETR, NO1TSO, NOTESO,
     %             NBTEDS, NOTEDS )
      IF( NBTEDS .GT. MXTEDS ) GOTO 9000

C     RECHERCHE DE LA SECONDE DIAGONALE
      NAR2 = NOAROPTE( NAR )

C     LES 2 SOMMETS EXTREMITES DE CETTE ARETE NAR2
C     --------------------------------------------
      NS3 = NSTETR( NOSOARTE(1,NAR2), NOTET0 )
      IF( NS3 .LE. 0 ) GOTO 9000

      NS4 = NSTETR( NOSOARTE(2,NAR2), NOTET0 )
      IF( NS4 .LE. 0 ) GOTO 9000

      IF( NESSAI .LE. 2 ) THEN

C        ESSAI DE TETRAEDRISER LES NBTEDS TETRAEDRES EN REMPLACANT
C        LA PREMIERE DIAGONALE PAR LA SECONDE
C        ---------------------------------------------------------
         QUAET0 = 2.222
         QUAET1 = 2.222
         VOLET0 = 0
         VOLET1 = 0

         DO J=1,NBTEDS

            NT = NOTEDS(J)
            CALL QUATET( XYZSOM( 1, NSTETR(1,NT) ),
     %                   XYZSOM( 1, NSTETR(2,NT) ),
     %                   XYZSOM( 1, NSTETR(3,NT) ),
     %                   XYZSOM( 1, NSTETR(4,NT) ),
     %                   ARMIN, ARMAX, SURFTR, VTE0, QTE0 )

            QUALIT( NT ) = QTE0

C           VOLUME INITIAL ET QUALITE INITIALE DES TETRAEDRES DE L'ETOILE
            VOLET0 = VOLET0 + ABS( VTE0 )
            QUAET0 = MIN( QUAET0, QTE0 )

            IF( NT .NE. NOTET0 ) THEN

C              TETRAEDRE DE L'ESSAI EN ECHANGEANT LES SOMMETS DES DIAGONALES
               DO I=1,4
                  NOSOTE(I) = NSTETR( I, NT )
               ENDDO

               DO I=1,4
                  IF( NOSOTE(I) .EQ. NS1 ) GOTO 5
               ENDDO
               GOTO 50

 5             DO K=1,4
                  IF( NOSOTE(K) .EQ. NS2 ) GOTO 10
               ENDDO
               GOTO 50

C              ECHANGE DE LA 1-ERE DIAGONALE EN LA SECONDE
C              SI LE SOMMET ECHANGE N'EST PAS DEJA SOMMET DE NOSOTE
 10            DO 20 L=1,4
                  IF( L .EQ. I .OR. L .EQ. K ) GOTO 20
                  IF( NOSOTE(L) .EQ. NS3 ) GOTO 30
 20            ENDDO
               NOSOTE(I) = NS3

 30            DO 40 L=1,4
                  IF( L .EQ. I .OR. L .EQ. K ) GOTO 40
                  IF( NOSOTE(L) .EQ. NS4 ) GOTO 50
 40            ENDDO   
               NOSOTE(K) = NS4

 50            CALL QUATET( XYZSOM( 1, NOSOTE(1) ),
     %                      XYZSOM( 1, NOSOTE(2) ),
     %                      XYZSOM( 1, NOSOTE(3) ),
     %                      XYZSOM( 1, NOSOTE(4) ),
     %                      ARMIN, ARMAX, SURFTR, VTE, QTE )


              print*,'sitediag: nstetr=',(nstetr(kk,nt),kk=1,4),
     %               ' V0=',vte0,' Q0=',qte0
              print*,'sitediag: nosote=', nosote,' V1=',vte,' Q1=',qte


               IF( VTE .LT. 0.0 ) THEN
                  PRINT *,'sitediag: NOSOTE=',(NOSOTE(k),k=1,4),
     %              ' VOLUME=',VTE,' QUALITE=',QTE
                  K         = NOSOTE(2)
                  NOSOTE(2) = NOSOTE(3)
                  NOSOTE(3) = K
                  GOTO 50
               ENDIF

C              VOLUME ET QUALITE DES TETRAEDRES DE L'ETOILE ESSAYEE
               VOLET1 = VOLET1 + ABS( VTE )
               QUAET1 = MIN( QUAET1, QTE )

            ENDIF

         ENDDO


         PRINT *,'sitediag: VOLET0=',VOLET0,' VOLET1=',VOLET1,
     %           ' DIFFERENCE', ABS(VOLET1-VOLET0)/VOLET0*100,' %'
         PRINT *,'sitediag: QUAET0=',QUAET0,' QUAET1=',QUAET1
 
         IF( QUAET1 .NE. 2.222 .AND. QUAET1-QUAET0 .GT. 0.001 ) THEN
            GOTO 9999
         ENDIF

         IF( NESSAI .EQ. 1 ) THEN
C           ESSAI 2: ECHANGE DES DIAGONALES A PARTIR DE L'AUTRE DIAGONALE
            NESSAI = 2
            NAR = NAR2
            GOTO 1
         ENDIF

      ENDIF

C     PROBLEME RENCONTRE OU QUALITE INFERIEURE
C     ----------------------------------------
 9000 NBTEDS = 0
      QUAET1 = 0
      NESSAI = 0

 9999 RETURN
      END
