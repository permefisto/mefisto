      SUBROUTINE SITEPT( P,      NOTET0, XYZSOM,
     %                   NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                   N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                   NBTEDS, MXTETR, NOTEDS, VOLUMT, QUALIT,
     %                   VOLET0, QUAET0, VOLET1, QUAET1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SIMULATION DE LA TETRAEDRISATION ETOILEE PAR LE POINT P
C -----    APPARTENANT AU TETRAEDRE NOTET0
C
C ENTREES:
C --------
C P      : POINT ETOILANT
C NOTET0 : NUMERO DANS NSTETR DU TETRAEDRE CONTENANT LE POINT P
C XYZSOM : LES 3 COORDONNEES DES SOMMETS
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NSTETR(>3)
C NSTETR : NUMERO DES 4 SOMMETS DE CHAQUE TETRAEDRE
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DANS NOTESO DE CHAQUE SOMMET
C NOTESO : NOTESO(1,*) NUMERO DANS NSTETR DU TETRAEDRE
C          NOTESO(2,*) NUMERO DANS NOTESO DU TETRAEDRE SUIVANT
C                      0 SI C'EST LE DERNIER
C
C MODIFIES :
C ----------
C N1FEVI : NUMERO DANS NFETOI DE LA PREMIERE FACE VIDE
C N1FEOC : NUMERO DANS NFETOI DE LA PREMIERE FACE OCCUPEE
C NFETOI : NUMERO DES 3 SOMMETS, TETRAEDRE ET SUIVANTE DES FACES
C          VUES UNE FOIS DE L'ETOILE CHAINEES N1FEOC PUIS NFETOI(5,*)
C VOETOI : VOLUME   DU TETRAEDRE INITIAL ASSOCIE A LA FACE
C QUETOI : QUALITE  DU TETRAEDRE INITIAL ASSOCIE A LA FACE
C VOLUMT : VOLUME  DES TETRAEDRES ACTUELS
C QUALIT : QUALITE DES TETRAEDRES ACTUELS
C
C SORTIES:
C --------
C NBTEDS : NOMBRE DE TETRAEDRES A DETRUIRE
C NOTEDS : NUMERO DANS NSTETR DES TETRAEDRES A DETRUIRE
C VOLET0 : VOLUME  INITIAL  DE LA TETRAEDRISATION ETOILEE PAR P
C QUAET0 : QUALITE INITIALE DE LA TETRAEDRISATION ETOILEE PAR P
C VOLET1 : VOLUME  FINAL    DE LA TETRAEDRISATION ETOILEE PAR P
C QUAET1 : QUALITE FINALE   DE LA TETRAEDRISATION ETOILEE PAR P
C          =0 INDIQUE UNE IMPOSSIBLITE D'OBTENIR UNE ETOILE POUR P
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   DECEMBRE 1991
C....................................................................012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           NSTETR(NBSOTE,*),
     %                  NO1TSO(*),
     %                  NOTESO(2,*),
     %                  NFETOI(5,*),
     %                  NOTEDS(MXTETR)
      REAL              P(3), XYZSOM(3,*),
     %                  VOETOI(*), QUETOI(*), QUALIT(*), VOLUMT(*)

      REAL              ARMIN, ARMAX, SURFTR(4),
     %                  VOLET0, QUAET0, VOLET1, QUAET1

      INTEGER           NOSOFATETE(3,4)
C     NO DES SOMMETS DES FACES POUR QUE VU DU SOMMET MANQUANT
C     LES SOMMETS SOIENT VUS DANS LE SENS DIRECT OU ENCORE
C     NORMALE AUX FACES DIRIGEES VERS L'INTERIEUR DU TETRAEDRE
      DATA              NOSOFATETE/ 1,2,3, 2,4,3, 3,4,1, 4,2,1 /

      VOLET1 = 0.0
      GRAND  = RINFO( 'GRAND' )
      QUAET1 = GRAND

C     UN TETRAEDRE A SUPPRIMER
      NBTEDS = 1
      NOTEDS( 1 ) = NOTET0

C     VOLUME ET QUALITE INITIALE DE L'ETOILE
      VOLET0 = ABS( VOLUMT( NOTET0 ) )
      QUAET0 = QUALIT( NOTET0 )
      IF( QUAET0 .GT. 0.5 ) THEN
C        QUALITE DIFFICILE A AMELIORER -> TETRAEDRE NON MODIFIE
         QUAET1 = QUAET0
         VOLET1 = VOLET0
         GOTO 9999
      ENDIF

C     REINITIALISATION A VIDE DES FACES DE NFETOI
C     -------------------------------------------
      CALL ETOIOCVI( N1FEOC, N1FEVI, NFETOI )

C     AJOUT DES 4 FACES DU TETRAEDRE NOTET0 A L'ETOILE
C     ------------------------------------------------
      DO I=1,4
C        SI    ( LA FACE I DU TETRAEDRE NOTET0 N'APPARTIENT PAS
C                AUX FACES DE L'ETOILE NFETOI )
C        ALORS ELLE EST AJOUTEE    A L'ETOILE NFETOI
C        SINON ELLE EST SUPPRIMEE DE L'ETOILE NFETOI
         CALL AJFACE( 1, NOTET0, I, NBSOTE, NSTETR,
     %                N1FEOC, N1FEVI, NFETOI,
     %                NF  )
         IF( NF .LT. 0 ) THEN
C           SATURATION DES FACES DE L'ETOILE
            GOTO 9000
         ELSE IF( NF .GT. 0 ) THEN
C           CALCUL DU VOLUME ET QUALITE DU TETRAEDRE DU A LA FACE
            CALL QUATET( XYZSOM(1,ABS(NFETOI(1,NF))),
     %                   XYZSOM(1,NFETOI(2,NF)),
     %                   XYZSOM(1,NFETOI(3,NF)),
     %                   P,
     %                   ARMIN, ARMAX, SURFTR, VOETOI(NF), QUETOI(NF) )
         ENDIF
      ENDDO

C     RECENSEMENT DES TETRAEDRES DE L'ETOILE DU POINT N SELON UNE
C     RECHERCHE PAR LES FACES RECENSEES DES TETRAEDRES ADJACENTS
C     -----------------------------------------------------------
 40   NF1 = N1FEOC

 45   IF( NF1 .GT. 0 ) THEN

C        LE NO DU 1-ER SOMMET DE LA FACE NF1
         NS1 = NFETOI(1,NF1)
         IF( NS1 .LT. 0 ) THEN

C           FACE DEJA TRAITEE . PASSAGE A LA SUIVANTE
            NF1 = NFETOI(5,NF1)
            GOTO 45

         ELSE

C           LA FACE ET SON TETRAEDRE SONT ANALYSES ENSUITE
C           LE TETRAEDRE NOTET1 DE L'AUTRE COTE DE LA FACE NF1
            CALL NOTEFO( NFETOI(1,NF1), NFETOI(2,NF1), NFETOI(3,NF1),
     %                   NFETOI(4,NF1), NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                   J, NOTET1 )

C           LE TEMOIN DE RECHERCHE EFFECTUEE
            NFETOI(1,NF1) = -NS1

            IF( NOTET1 .LE. 0 ) THEN
C              PAS DE TETRAEDRE OPPOSE : LA FACE EST DONC TRAITEE
C              LA FACE SUIVANTE
               NF1 = NFETOI(5,NF1)
               GOTO 45
            ENDIF

C           VOLUME ET QUALITE DU TETRAEDRE OPPOSE NOTET1
            CALL QUATET( XYZSOM(1,NSTETR(1,NOTET1)),
     %                   XYZSOM(1,NSTETR(2,NOTET1)),
     %                   XYZSOM(1,NSTETR(3,NOTET1)),
     %                   XYZSOM(1,NSTETR(4,NOTET1)),
     %                   ARMIN, ARMAX, SURFTR,
     %                   VOLUMT(NOTET1), QUALIT(NOTET1) )
            V0 = ABS( VOETOI(NF1) ) + ABS( VOLUMT(NOTET1) )

C           CALCUL DU VOLUME ET QUALITE DES 3 FACES DU TETRAEDRE OPPOSE
            VET = 0.0
            QET = GRAND
            DO I=1,4
               IF( I .NE. J ) THEN

C                 FACE DIFFERENTE DE CELLE COMMUNE
                  CALL QUATET( XYZSOM(1,NSTETR(NOSOFATETE(1,I),NOTET1)),
     %                         XYZSOM(1,NSTETR(NOSOFATETE(2,I),NOTET1)),
     %                         XYZSOM(1,NSTETR(NOSOFATETE(3,I),NOTET1)),
     %                         P,
     %                         ARMIN, ARMAX, SURFTR, VF, QF )

C                 VOLUME DE L'ETOILE
                  VET = VET + ABS( VF )

C                 QUALITE DE L'ETOILE
                  QET = MIN( QET, QF )

               ENDIF
            ENDDO

C           COMPARAISON DES VOLUMES DES ETOILES
            IF( ABS(V0-VET) .LT. 0.001 * V0 ) THEN

C              VOLUMES EGAUX . QUALITE MEILLEURE ?
               IF( MIN( QUETOI(NF1), QUALIT(NOTET1) ) .LT. QET ) THEN
C
C                 AVEC CE TETRAEDRE => MEILLEURE QUALITE
C                 UN TETRAEDRE DE PLUS A DETRUIRE
                  IF( NBTEDS .GE. MXTETR ) GOTO 9000

                  NBTEDS = NBTEDS + 1
                  NOTEDS( NBTEDS ) = NOTET1

C                 VOLUME INITIAL
                  VOLET0 = VOLET0 + ABS( VOLUMT( NOTET1 ) )

C                 QUALITE INITIALE
                  QUAET0 = MIN( QUAET0, QUALIT( NOTET1 ) )

C                 AJOUT OU RETRAIT DES 4 FACES DU TETRAEDRE AJOUTE A L'ETOILE
                  DO I=1,4
C                    SI    ( LA FACE I DU TETRAEDRE NOTET1 N'APPARTIENT PAS
C                            AUX FACES DE L'ETOILE NFETOI )
C                    ALORS ELLE EST AJOUTEE    A L'ETOILE NFETOI
C                    SINON ELLE EST SUPPRIMEE DE L'ETOILE NFETOI
                     CALL AJFACE( 1, NOTET1, I, NBSOTE, NSTETR,
     %                            N1FEOC, N1FEVI, NFETOI,
     %                            NF  )
                     IF( NF .LT. 0 ) THEN
C                       SATURATION DES FACES DE L'ETOILE
                        GOTO 9000
                     ELSE IF( NF .GT. 0 ) THEN
C                       CALCUL DU VOLUME ET QUALITE DU TETRAEDRE DU A LA FACE
                        CALL QUATET( XYZSOM(1,ABS(NFETOI(1,NF))),
     %                               XYZSOM(1,NFETOI(2,NF)),
     %                               XYZSOM(1,NFETOI(3,NF)),
     %                               P, ARMIN, ARMAX, SURFTR,
     %                               VOETOI(NF), QUETOI(NF) )
                     ENDIF
                  ENDDO

               ENDIF

            ENDIF

C           LA FACE SUIVANTE
            GOTO 40
         ENDIF
      ENDIF

C     ICI, L'ETOILE EST COMPLETE et CALCUL DU VOLUME ET DE SA QUALITE
C     ---------------------------------------------------------------
      VOLET1 = 0.0
      QUAET1 = GRAND
      NF     = N1FEOC

 90   IF( NF .GT. 0 ) THEN

C        CALCUL DU VOLUME ET QUALITE DU TETRAEDRE DU A LA FACE NF
         CALL QUATET( XYZSOM(1,ABS(NFETOI(1,NF))),
     %                XYZSOM(1,NFETOI(2,NF)),
     %                XYZSOM(1,NFETOI(3,NF)),
     %                P, ARMIN, ARMAX, SURFTR,
     %                VOETOI( NF ), QUETOI( NF ) )

C        VOLUME DE L'ETOILE
         VOLET1 = VOLET1 + ABS( VOETOI( NF ) )

C        QUALITE DE L'ETOILE
         IF( QUETOI(NF) .LT. QUAET1 ) THEN
            QUAET1 = QUETOI( NF )
         ENDIF

C        PASSAGE A LA FACE SUIVANTE
         NF = NFETOI(5,NF)
         GOTO 90

      ENDIF
      GOTO 9999

C     PROBLEME RENCONTRE
 9000 NBTEDS = 0
      QUAET1 = 0.0
      VOLET1 = 0.0
C
 9999 RETURN
      END
