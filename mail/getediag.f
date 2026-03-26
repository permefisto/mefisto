      SUBROUTINE GETEDIAG( NESSAI, NBSOMM, XYZSOM, NBDM,   NUDMEF,
     %                     NAR0,   NOTET0, NBSOTE, N1TEVI, NSTETR,
     %                     NO1TSO, N1TESO, NOTESO,
     %                     MXTEDS, NBTEDS, NOTEDS,
     %                     MXTETA, VOLUMT, QUALIT,
     %                     VOLET0, QUAET0, VOLET1, QUAET1, IERR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERATION DE LA TETRAEDRISATION A PARTIR DES TETRAEDRES
C -----    AVEC SUPPRESSION D'UNE DES 2 DIAGONALES DU TETRAEDRE NT
C          A EXECUTER EN CAS DE SUCCES DE sitediag.f

C ENTREES:
C --------
C NESSAI : 1 ou 2 NUMERO DE LA DIAGONALE A ELIMINER

C MXSOMM : NOMBRE MAXIMAL DE SOMMETS DECLARABLES DANS XYZSOM
C XYZSOM : LES 3 COORDONNEES DES SOMMETS
C NBDM   : 0 SI 1 MATERIAU=VOLUME, SINON NOMBRE DE MATERIAUX DU VOLUME
C NUDMEF : NUMERO DE MATERIAU DE CHAQUE TETRAEDRE DU MAILLAGE
C          ATTENTION: CE TABLEAU EXISTE SEULEMENT SI NBDM>0
C NAR0   : NUMERO ARETE DANS NOTET0 DE LA 1-ERE DIAGONALE
C NOTET0 : NUMERO NSTETR DU TETRAEDRE DE DIAGONALE A SUPPRIMER
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NSTETR(>3)
C NSTETR : NUMERO DES 4 SOMMETS DE CHAQUE TETRAEDRE

C NO1TSO : NUMERO DU 1-ER TETRAEDRE DANS NOTESO DE CHAQUE SOMMET
C NOTESO : NOTESO(1,*) NUMERO DANS NSTETR DU TETRAEDRE
C          NOTESO(2,*) NUMERO DANS NOTESO DU TETRAEDRE SUIVANT
C                      0 SI C'EST LE DERNIER
C MODIFIES :
C ----------
C NBTEDS : NOMBRE DE TETRAEDRES A DETRUIRE PUIS CREES
C NOTEDS : NUMERO DANS NSTETR DES TETRAEDRES A DETRUIRE PUIS CREES
C NBSOMM : NOMBRE DE SOMMETS AVANT ET APRES
C MXTETA : NUMERO MAXIMAL D'UN TETRAEDRE ACTUEL
C VOLUMT : VOLUME  DES TETRAEDRES DE L'ETOILE
C QUALIT : QUALITE DES TETRAEDRES DE L'ETOILE
C VOLET1 : VOLUME TOTAL DE L'ETOILE CREEE
C QUAET1: QUALITE MINIMALE DES TETRAEDRES DE L'ETOILE CREEE
C IERR   : =0 SI PAS D'ERREUR
C          >0 SI SATURATION D'UN TABLEAU  => ERREUR DE TETRAEDRISATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY  Fevrier 2016
C....................................................................012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      REAL              XYZSOM(3,*),
     %                  VOLUMT(*),
     %                  QUALIT(*)
      INTEGER           NSTETR(NBSOTE,MXTETA),
     %                  NUDMEF(*),
     %                  NO1TSO(*),
     %                  NOTESO(2,*),
     %                  NOTEDS(MXTEDS),
     %                  NOSOTE(4)
      REAL              ARMIN, ARMAX, SURFTR(4), VOLET1, QUAET1

C                       NO SOMMETS D'UNE ARETE D'UN TETRAEDRE
      INTEGER           NOSOARTE(2,6)
      DATA              NOSOARTE / 1,2, 2,3, 3,1, 4,1, 4,2, 4,3 /
C                       NO ARETE OPPOSEE DANS UN TETRAEDRE
      INTEGER           NOAROPTE(6)
      DATA              NOAROPTE/  6, 4, 5, 2, 3, 1 /

      GRAND  = RINFO( 'GRAND' )
      QUAET0 = GRAND
      QUAET1 = GRAND
      VOLET0 = 0
      VOLET1 = 0
      NBNEWT0 = MXTEDS/2
      NBNEWT1 = NBNEWT0

C     TOUS LES TETRAEDRES NOTEDS DOIVENT ETRE DANS LE MEME
C     MATERIAU NUMATE
      NUMATE = 0
      IF( NBDM .GT. 0 ) THEN
         NUMATE = NUDMEF( NOTEDS(1) )
         DO NF = 2, NBTEDS
            NT = NOTEDS(NF)
            IF( NUMATE .NE. NUDMEF( NOTEDS(NF) ) ) THEN
               IERR = 2
               GOTO 9999
            ENDIF
         ENDDO
      ENDIF

      IF( NESSAI .EQ. 2 ) THEN
C         TRAITEMENT DE LA SECONDE DIAGONALE
          NAR1 = NOAROPTE( NAR0 )
       ELSE
C         TRAITEMENT DE LA PREMIERE DIAGONALE
          NAR1 = NAR0
       ENDIF

C     LES 2 SOMMETS EXTREMITES DE CETTE ARETE NAR1 DIAGONALE
      NS1 = NSTETR( NOSOARTE(1,NAR1), NOTET0 )
      IF( NS1 .LE. 0 ) GOTO 9000

      NS2 = NSTETR( NOSOARTE(2,NAR1), NOTET0 )
      IF( NS2 .LE. 0 ) GOTO 9000

C     RECHERCHE DES TETRAEDRES D'ARETE NS1-NS2
      CALL NOTEAR( NS1,    NS2,
     %             NBSOTE, NSTETR, NO1TSO, NOTESO,
     %             NBTEDS, NOTEDS )
      IF( NBTEDS .GT. MXTEDS ) GOTO 9000

C     LA SECONDE DIAGONALE
      NAR2 = NOAROPTE( NAR1 )

C     LES 2 SOMMETS EXTREMITES DE CETTE ARETE NAR2
      NS3 = NSTETR( NOSOARTE(1,NAR2), NOTET0 )
      IF( NS3 .LE. 0 ) GOTO 9000

      NS4 = NSTETR( NOSOARTE(2,NAR2), NOTET0 )
      IF( NS4 .LE. 0 ) GOTO 9000

C     TETRAEDRISER LES NBTEDS TETRAEDRES EN REMPLACANT
C     LA PREMIERE DIAGONALE PAR LA SECONDE
C     ------------------------------------------------
      DO J=1,NBTEDS

         NT = NOTEDS(J)
         CALL QUATET( XYZSOM( 1, NSTETR(1,NT) ),
     %                XYZSOM( 1, NSTETR(2,NT) ),
     %                XYZSOM( 1, NSTETR(3,NT) ),
     %                XYZSOM( 1, NSTETR(4,NT) ),
     %                ARMIN, ARMAX, SURFTR, VTE0, QTE0 )

         QUALIT( NT ) = QTE0

C        VOLUME INITIAL ET QUALITE INITIALE DES TETRAEDRES DE L'ETOILE
         VOLET0 = VOLET0 + ABS( VTE0 )
         QUAET0 = MIN( QUAET0, QTE0 )

         IF( NT .NE. NOTET0 ) THEN

C           TETRAEDRE EN ECHANGEANT LES SOMMETS DES DIAGONALES
            DO I=1,4
               NOSOTE(I) = NSTETR( I, NT )
            ENDDO

            DO I=1,4
               IF( NOSOTE(I) .EQ. NS1 ) GOTO 5
            ENDDO
            GOTO 50

 5          DO K=1,4
               IF( NOSOTE(K) .EQ. NS2 ) GOTO 10
            ENDDO
            GOTO 50

C           ECHANGE DE LA 1-ERE DIAGONALE EN LA SECONDE
C           SI LE SOMMET ECHANGE N'EST PAS DEJA SOMMET DE NOSOTE
 10         DO 20 L=1,4
               IF( L .EQ. I .OR. L .EQ. K ) GOTO 20
               IF( NOSOTE(L) .EQ. NS3 ) GOTO 30
 20         ENDDO
            NOSOTE(I) = NS3

 30         DO 40 L=1,4
               IF( L .EQ. I .OR. L .EQ. K ) GOTO 40
               IF( NOSOTE(L) .EQ. NS4 ) GOTO 50
 40         ENDDO   
            NOSOTE(K) = NS4

 50         CALL QUATET( XYZSOM( 1, NOSOTE(1) ),
     %                   XYZSOM( 1, NOSOTE(2) ),
     %                   XYZSOM( 1, NOSOTE(3) ),
     %                   XYZSOM( 1, NOSOTE(4) ),
     %                   ARMIN, ARMAX, SURFTR, VTE, QTE )

            print*,'getediag: nstetr=',(nstetr(kk,nt),kk=1,4),
     %             ' V0=',vte0,' Q0=',qte0
            print*,'getediag: nosote=', nosote,' V1=',vte,' Q1=',qte

            IF( VTE .LT. 0.0 ) THEN
               PRINT *,'getediag: NOSOTE=',(NOSOTE(k),k=1,4),
     %                  ' V1=',VTE,' Q1=',QTE,' S2<->S3'
               K         = NOSOTE(2)
               NOSOTE(2) = NOSOTE(3)
               NOSOTE(3) = K
               GOTO 50
            ENDIF

C           VOLUME ET QUALITE DES TETRAEDRES DE L'ETOILE ESSAYEE
            VOLET1 = VOLET1 + ABS( VTE )
            QUAET1 = MIN( QUAET1, QTE )

C           RECHERCHE D'UN TETRAEDRE VIDE
            IF( N1TEVI .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
         KERR(1)='getediag: SATURATION du TABLEAU NSTETR des TETRAEDRES'
               CALL LEREUR
               IERR = 1
               GOTO 9999
            ENDIF
C
C           LE PREMIER TETRAEDRE VIDE
            NT = N1TEVI
C           EST MIS A JOUR
            N1TEVI = NSTETR( 2, N1TEVI )

C           UN TETRAEDRE CREE DE PLUS DANS NOTEDS
            NBNEWT1 = NBNEWT1 + 1
            NOTEDS( NBNEWT1 ) = NT

      call verifns( 'getediag:', nbsomm, nbsote, nstetr, no1tso, noteso,
     %               nbtensv, ierr )

C           FORMATION DU TETRAEDRE DE SOMMET NBSOMM ET DE FACE OPPOSEE NF
            DO L=1,4
               NSTETR(L,NT) = NOSOTE(L)
            ENDDO

            IF( NBDM .GT. 0 ) NUDMEF( NT ) = NUMATE

C           AJOUT DE CE TETRAEDRE DANS LE CHAINAGE NOTESO DES TETRAEDRES
C           DES 4 SOMMETS DU TETRAEDRE NT
            CALL CHTESO( NT, NBSOTE, NSTETR, NO1TSO, N1TESO, NOTESO,
     %                   IERR )
            IF( IERR .NE. 0 ) GOTO 9999

C           QUALITE ET VOLUME DU TETRAEDRE
            CALL QUATET( XYZSOM(1,NSTETR(1,NT)),
     %                   XYZSOM(1,NSTETR(2,NT)),
     %                   XYZSOM(1,NSTETR(3,NT)),
     %                   XYZSOM(1,NSTETR(4,NT)),
     %                   ARMIN, ARMAX, SURFTR, VOLUMT(NT), QUALIT(NT) )

ccc            VOLET1 = VOLET1 + ABS( VOLUMT(NT) )
ccc            QUAET1 = MIN( QUAET1, QUALIT(NT) )

         ENDIF

      ENDDO

C     DESTRUCTION DES NBTEDS TETRAEDRES DE L'ETOILE
      DO L = 1, NBTEDS
         NT = NOTEDS( L )
         CALL DS1TET( NT, NBSOTE, N1TEVI, NSTETR,
     %                NO1TSO, N1TESO, NOTESO )
         QUALIT( NT ) = GRAND
         IF( NBDM .GT. 0 ) NUDMEF( NT ) = 0
      ENDDO

C     MISE AU DEBUT DE NOTEDS DES NBNEWT1-NBNEWT0 TETRAEDRES
      NBTEDS = 0
      DO L = NBNEWT0+1, NBNEWT1
         NBTEDS = NBTEDS + 1
         NOTEDS( NBTEDS ) = NOTEDS( L )
      ENDDO

      PRINT*,'getediag: VOLET0=',VOLET0,' QUAET0=',QUAET0
      PRINT*,'getediag: VOLET1=',VOLET1,' QUAET1=',QUAET1

      GOTO 9999

C     PROBLEME RENCONTRE OU QUALITE INFERIEURE
C     ----------------------------------------
 9000 NBTEDS = 0
      QUAET1 = 0
      VOLET1 = 0
      IERR = 3

 9999 RETURN
      END
