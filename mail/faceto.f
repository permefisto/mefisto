      SUBROUTINE FACETO( NS,     NBSOTE, NOSOTE, NO1TSO, NOTESO,
     %                   N1FEVI, N1FEOC, NFETOI, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CHAINER DANS NFETOI LES FACES DE L'ETOILE DU POINT NS
C -----   VUES UNE SEULE FOIS DANS LES TETRAEDRES DE L'ETOILE DE NS
C
C ENTREES:
C --------
C NS     : NUMERO DANS XYZSOM DU POINT D'ETOILE A TRAITER
C XYZSOM : COORDONNEES X Y Z DES NBSOMM SOMMETS DES TETRAEDRES
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NOSOTE(>3)
C NOSOTE : NUMERO DES 4 SOMMETS DE CHAQUE TETRAEDRE
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DANS NOTESO DE CHAQUE SOMMET
C NOTESO : NOTESO(1,*) NUMERO DANS NOSOTE DU TETRAEDRE
C          NOTESO(2,*) NUMERO DANS NOTESO DU TETRAEDRE SUIVANT
C                      0 SI C'EST LE DERNIER
C
C MODIFIES :
C ----------
C N1FEVI : NUMERO DANS NFETOI DE LA PREMIERE FACE VIDE
C N1FEOC : NUMERO DANS NFETOI DE LA PREMIERE FACE OCCUPEE
C NFETOI : LES FACES VUES UNE FOIS DE L'ETOILE CHAINEES PAR
C          N1FEOC PUIS NFETOI(5,*)
C IERR   : =0 SI PAS D'ERREUR
C          >0 SINON ( SATURATION DES FACES DE L'ETOILE)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   DECEMBRE 1991
C....................................................................012
      include"./incl/gsmenu.inc"
      INTEGER           NOSOTE(NBSOTE,*),
     %                  NO1TSO(*),
     %                  NOTESO(2,*),
     %                  NFETOI(5,*)
      INTEGER           NOSOFATE(1:3,1:4),NSF(1:3),NSFT(1:3),NSFE(1:3)
C     NO DES SOMMETS DES FACES POUR QUE VU DU SOMMET MANQUANT
C     LES SOMMETS SOIENT VUS DANS LE SENS DIRECT
      DATA              NOSOFATE/ 1,2,3, 2,4,3, 3,4,1, 4,2,1 /
C
C     REINITIALISATION A VIDE DES FACES ETOILANT LE POINT N-1
C     -------------------------------------------------------
      NF1 = 0
      NF2 = N1FEOC
 10   IF( NF2 .GT. 0 ) THEN
         NF1 = NF2
         NF2 = NFETOI(5,NF2)
         GOTO 10
      ENDIF
C     FIN DU CHAINAGE DES FACES
      IF( NF1 .GT. 0 ) THEN
C        LA DERNIERE FACE OCCUPEE EST CHAINEE SUR LA 1-ERE VIDE
         NFETOI(5,NF1) = N1FEVI
C        LA 1-ERE OCCUPEE DEVIENT LA 1-ERE VIDE
         N1FEVI = N1FEOC
      ENDIF
C
C     IL N'Y A PAS DE FACES DANS L'ETOILE
      N1FEOC = 0
C
C     PARCOURS DES TETRAEDRES DU SOMMET NS
C     ------------------------------------
C     POSITION DANS NOTESO DU 1-ER TETRAEDRE DE SOMMET NS
      NDT = NO1TSO( NS )
C
C     TANT QU'IL EXISTE UN TETRAEDRE DE SOMMET NS FAIRE
 20   IF( NDT .GT. 0 ) THEN
C
C        LE NUMERO DU TETRAEDRE DANS NOSOTE
         NT = NOTESO(1,NDT)
C
C        AJOUT OU RETRAIT DES 4 FACES DU TETRAEDRE A L'ETOILE
         DO 90 NF=1,4
C           SI    ( LA FACE I DU TETRAEDRE NOTET1 N'APPARTIENT PAS
C                   AUX FACES DE L'ETOILE NFETOI )
C           ALORS ELLE EST AJOUTEE A L'ETOILE DANS NFETOI
C           SINON ELLE EST EMPILEE DANS NPILE POUR ETRE DETRUITE ENSUITE
C                 ELLE EST SUPPRIMEE DE L'ETOILE NFETOI
C
C           LES 3 SOMMETS DE LA FACE NF DU TETRAEDRE VUS SENS DIRECT
            NSF(1) = NOSOTE( NOSOFATE(1,NF), NT )
            NSF(2) = NOSOTE( NOSOFATE(2,NF), NT )
            NSF(3) = NOSOTE( NOSOFATE(3,NF), NT )
C
C           LES 3 SOMMETS PAR ORDRE CROISSANT
            CALL TRI3NO( NSF, NSFT )
C
C           LE DEBUT DU CHAINAGE DES FACES DE L'ETOILE
C           ET BOUCLE SUR LES FACES ACTUELLES DE L'ETOILE
            NF1 = 0
            NF2 = N1FEOC
C
 30         IF( NF2 .GT. 0 ) THEN
C              TRI CROISSANT DES SOMMETS DE LA FACE
C              VOIR PLUS TARD SI STOCKAGE CROISSANT + SENS
               CALL TRI3NO( NFETOI(1,NF2), NSFE )
               IF( NSFE(1) .EQ. NSFT(1) ) THEN
                  IF( NSFE(2) .EQ. NSFT(2) ) THEN
                     IF( NSFE(3) .EQ. NSFT(3) ) THEN
C
C                       LA FACE NF2 EST DEJA DANS L'ETOILE
C                       DESTRUCTION DE LA FACE DE L'ETOILE
                        IF( NF1 .NE. 0 ) THEN
C                          LA FACE NF2,PRECEDEE DE NF1,N'EST PAS LA
C                          PREMIERE DE L'ETOILE
                           NFETOI(5,NF1) = NFETOI(5,NF2)
C                          LA FACE NF2 DEVIENT LA PREMIERE VIDE DE L'ETOILE
                           NFETOI(5,NF2) = N1FEVI
                           N1FEVI = NF2
                        ELSE
C                          LA FACE NF2=N1FEOC EST LA 1-ERE DE L'ETOILE
                           NF1    = NFETOI(5,N1FEOC)
                           NFETOI(5,N1FEOC) = N1FEVI
                           N1FEVI = N1FEOC
                           N1FEOC = NF1
                        ENDIF
                        GOTO 90
                     ENDIF
                  ENDIF
               ENDIF
C              LA FACE EST DIFFERENTE.PASSAGE A LA SUIVANTE
               NF1 = NF2
               NF2 = NFETOI(5,NF2)
               GOTO 30
            ENDIF
C
C           LA FACE UNIQUE EST AJOUTEE A L'ETOILE AU DEBUT DES FACES OCCUPEES
            IF( N1FEVI .LE. 0 ) THEN
C              SATURATION DES FACES DE L'ETOILE
               NBLGRC(NRERR) = 1
               KERR(1) = 'FACETO: SATURATION DES FACES DE L''ETOILE'
               CALL LEREUR
               IERR = 1
               RETURN
            ENDIF
            NF1           = N1FEVI
            N1FEVI        = NFETOI(5,N1FEVI)
C           LES 3 SOMMETS DANS LE SENS DIRECT
            NFETOI(1,NF1) = NSF(1)
            NFETOI(2,NF1) = NSF(2)
            NFETOI(3,NF1) = NSF(3)
C           LA FACE SUIVANTE EST LA PREMIERE AUPARAVANT
            NFETOI(5,NF1) = N1FEOC
            N1FEOC        = NF1
 90      CONTINUE
C
C        LE TETRAEDRE SUIVANT
         NDT = NOTESO(2,NDT)
         GOTO 20
      ENDIF
      END
