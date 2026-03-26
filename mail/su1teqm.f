      SUBROUTINE SU1TEQM( NTE,    N1TEVI, NUDTETR, NOTETR,
     %                    N1TETS, PTXYZD, NBTECF,  NOTECF,
     %                    NBTNEW, NOTNEW )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   TENTATIVE DE SUPPRIMER CE TETRAEDRE NTE (DE MAUVAISE QUALITE)
C -----   A PARTIR D'UN COUPLE DE TETRAEDRES OPPOSES ET ADJACENTS

C ENTREES:
C --------
C NTE    : NUMERO NOTETR DU TETRAEDRE A TRAITER
C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NOTETR
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C N1TETS : N1TETS(I) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET I
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE

C MODIFIE:
C --------
C NBTECF : NOMBRE DE NO NOTETR DE TETRAEDRES RANGES DANS NOTECF
C NOTECF : NO NOTETR DES NBTECF TETRAEDRES
C          SI MODIFICATION ALORS LES 3 TETRAEDRES SONT SUPPRIMES DE NOTECF
C          ET LES 2 TETRAEDRES CREES SONT AJOUTES

C SORTIES:
C --------
C NBTNEW : =0 PAS DE TETRAEDRE MODIFIE OU SUPPRIME
C          >0 NOMBRE DE NOUVEAUX TETRAEDRES
C NOTNEW : NUMERO NOTETR DES 2 TETRAEDRES CREES SI NBTNEW>0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET VEULETTES & St PIERRE du PERRAY  Janvier 2018
C2345X7..............................................................012
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0

      DOUBLE PRECISION  PTXYZD(4,*)
      INTEGER           NOTETR(8,*), NOTECF(NBTECF), N1TETS(*)
      CHARACTER*80      KTITRE

      INTEGER           NTEAR(2), NOST(2), NOTNEW(2), NVOLTE(1)
      EQUIVALENCE      (NOST(1),  NS1),   (NOST(2),  NS2),
     %                 (NTEAR(1), NTEA1), (NTEAR(2), NTEA2)

      DOUBLE PRECISION  ARMIN,  ARMAX, SURFTR(4), V, V1, V2
      REAL              Q5TE(5)

      INTEGER           NOSOARTE(2,6)
      DATA              NOSOARTE / 1,2, 2,3, 3,1, 4,1, 4,2, 4,3 /
C     NO DES 2 SOMMETS EXTREMITES DES 6 ARETES D'UN TETRAEDRE

      INTEGER           NOFAARTE(2,6)
      DATA              NOFAARTE / 1,4,  1,2,  1,3,  3,4,  2,4,  2,3 /
C     NO DES 2 FACES DU TETRAEDRE ADJACENTES A CHAQUE ARETE DU TETRAEDRE

      TRACTE0 = TRACTE

C     VOLUME ET QUALITE DU TETRAEDRE NTE
      CALL QUATETD( PTXYZD(1,NOTETR(1,NTE)),
     %              PTXYZD(1,NOTETR(2,NTE)),
     %              PTXYZD(1,NOTETR(3,NTE)),
     %              PTXYZD(1,NOTETR(4,NTE)),
     %              ARMIN, ARMAX, SURFTR, V, Q0 )

      PRINT*,'su1teqm: SUPPRESSION? du TETRAEDRE',NTE,' :',
     %       (NOTETR(N,NTE),N=1,8),' NBTECF0=',NBTECF,' V0=',V,' Q0=',Q0
      DO N=1,4
         NTEOP = NOTETR( 4+N, NTE )
         IF( NTEOP .GT. 0 ) THEN
C           VOLUME ET QUALITE DU TETRAEDRE NTEOP
            CALL QUATETD( PTXYZD(1,NOTETR(1,NTEOP)),
     %                    PTXYZD(1,NOTETR(2,NTEOP)),
     %                    PTXYZD(1,NOTETR(3,NTEOP)),
     %                    PTXYZD(1,NOTETR(4,NTEOP)),
     %                    ARMIN, ARMAX, SURFTR, V1, Q1 )
            PRINT*,'su1teqm: TETRAEDRE OPPOSE',N,'(',NTEOP,'):',
     %             (NOTETR(M,NTEOP),M=1,8),' V=',V1,' Q=',Q1
         ENDIF
      ENDDO

C     NOMBRE DE TETRAEDRES CREES
      NBTNEW = 0

      DO 100 NA = 1, 6

         CALL CH3T2T( 1,      MXFACO, LEFACO,  0, 0, NVOLTE,
     %                NTE,    NA,     PTXYZD,
     %                N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                NTEOLD, NTENEW, Q5TE,   IERR )


C         DOIT REMPLACER CE QUI EST DESSOUS


         DO N = 1, 2

C           NO DES TETRAEDRES OPPOSES AUX 2 FACES D'ARETE COMMUNE NA
            NTEAR( N ) = NOTETR( 4 + NOFAARTE(N,NA) , NTE )
            IF( NTEAR( N ) .LE. 0 ) GOTO 100
            IF( NOTETR(1,NTEAR(N)) .LE. 0 ) GOTO 100

C           NO DES 2 SOMMETS NS1 NS2 DE L'ARETE NA DE NTE
            NOST( N ) = NOTETR( NOSOARTE(N,NA), NTE )

         ENDDO

C        LES 2 TETRAEDRES OPPOSES ONT ILS UN SOMMET COMMUN DIFFERENT DE NS1 NS2?
         DO K=1,4
            NS3 = NOTETR(K,NTEA1)
            IF( NS3 .NE. NS1 .AND. NS3 .NE. NS2 ) THEN

C              NS3 SOMMET K DE NTEA1 EST IL UN SOMMET DE NTEA2?
               DO 60 KK=1,4
                  NS4 = NOTETR(KK,NTEA2)
                  IF( NS3 .EQ. NS4 ) THEN

C                    OUI: CHOIX DU SOMMET NS3 N'APPARTENANT PAS A NTE
                     DO L=1,4
                        IF( NS3 .EQ. NOTETR(L,NTE) ) THEN
                           GOTO 60
                        ENDIF
                     ENDDO

C                    NS3 EST SOMMET DES 2 TETRAEDRES OPPOSES ET
C                      N'EST PAS UN SOMMET DE NTE
C                    NS1 NS2 NS3 EST LA FACE COMMUNE AUX 2 TETRAEDRES OPPOSES

C                    NS4 SOMMET DE NTEA1 NON NS1 NS2 NS3
                     DO L=1,4
                        NS4 = NOTETR(L,NTEA1)
                        IF( NS4.NE.NS1 .AND. NS4.NE.NS2 .AND.
     %                      NS4.NE.NS3 ) THEN
                           GOTO 10
                        ENDIF
                     ENDDO

C                    NS5 SOMMET DE NTEA2 NON NS1 NS2 NS3
 10                  DO LL=1,4
                        NS5 = NOTETR(LL,NTEA2)
                        IF( NS5.NE.NS1 .AND. NS5.NE.NS2 .AND.
     %                      NS5.NE.NS3 ) THEN
                           GOTO 20
                        ENDIF
                     ENDDO

C                    L'ARETE NS1-NS2 DES 3 TETRAEDRES NTE, NTEA1, NTEA2
C                    EST REMPLACEE PAR LA FACE COMMUNE NS4 NS3 NS5




C                    DECLARATION DES 2 NOUVEAUX TETRAEDRES
 20                  DO N=1,2
                        IF( N1TEVI .LE. 0 ) THEN
C                          SATURATION DU TABLEAU NOTETR DES TETRAEDRES
                         PRINT*,'su1teqm: SATURATION de NOTETR N1TEVI=',
     %                           N1TEVI
                           GOTO 9999
                        ENDIF
                        NOTNEW(N) = N1TEVI
                        NUDTETR = MAX( NUDTETR, N1TEVI )
                        N1TEVI = NOTETR(5,N1TEVI)
                     ENDDO

C                    REMPLISSAGE DU NUMERO PTXYZD DU 1-ER TETRAEDRE CREE
                     NTEN1 = NOTNEW(1)
                     NOTETR( 1, NTEN1 ) = NS4
 22                  NOTETR( 2, NTEN1 ) = NS3
                     NOTETR( 3, NTEN1 ) = NS5
                     NOTETR( 4, NTEN1 ) = NS1
C                    VOLUME ET QUALITE DU TETRAEDRE NTEN1
                     CALL QUATETD( PTXYZD(1,NOTETR(1,NTEN1)),
     %                             PTXYZD(1,NOTETR(2,NTEN1)),
     %                             PTXYZD(1,NOTETR(3,NTEN1)),
     %                             PTXYZD(1,NOTETR(4,NTEN1)),
     %                             ARMIN, ARMAX, SURFTR, V1, Q1 )
                     IF( Q1 .LE. 0. ) THEN
C                       PERMUTATION POUR OBTENIR UN VOLUME POSITIF
                        L   = NS3
                        NS3 = NS5
                        NS5 = L
                        GOTO 22
                     ENDIF

C                    REMPLISSAGE DU NUMERO PTXYZD DU 2-ND TETRAEDRE CREE
                     NTEN2 = NOTNEW(2)
                     NOTETR( 1, NTEN2 ) = NS4
                     NOTETR( 2, NTEN2 ) = NS5
                     NOTETR( 3, NTEN2 ) = NS3
                     NOTETR( 4, NTEN2 ) = NS2
C 24                 VOLUME ET QUALITE DU TETRAEDRE NTEN2
                     CALL QUATETD( PTXYZD(1,NOTETR(1,NTEN2)),
     %                             PTXYZD(1,NOTETR(2,NTEN2)),
     %                             PTXYZD(1,NOTETR(3,NTEN2)),
     %                             PTXYZD(1,NOTETR(4,NTEN2)),
     %                             ARMIN, ARMAX, SURFTR, V2, Q2 )

ccc                     IF( Q2 .LE. 0. ) THEN
cccC                       NON CAR LA FACE NS4-NS3-NS5 DOIT ETRE
cccC                       PARCOURUE DANS L'AUTRE SENS
cccC                       PERMUTATION POUR OBTENIR UN VOLUME POSITIF
ccc                        NOTETR( 2, NTEN2 ) = NS3
ccc                        NOTETR( 3, NTEN2 ) = NS5
ccc                        GOTO 24
ccc                     ENDIF

                     IF( Q1 .LT. Q0 .OR. Q2 .LT. Q0 ) THEN
C                       ABANDON DE LA SUPPRESSION DE CETTE ARETE NA DE NTE
                        NOTETR(1,NTEN2) = 0
                        NOTETR(5,NTEN2) = N1TEVI
                        N1TEVI = NTEN2
                        NOTETR(1,NTEN1) = 0
                        NOTETR(5,NTEN1) = N1TEVI
                        N1TEVI = NTEN1
                        GOTO 100
                     ENDIF

C                    CHANGEMENT ACCEPTE
C                    REMPLISSAGE DU NUMERO NOTETR DES TETRAEDRES OPPOSES
                     NOTETR( 5, NTEN1 ) = NTEN2
                     NOTETR( 6, NTEN1 ) = -1
                     NOTETR( 7, NTEN1 ) = -1
                     NOTETR( 8, NTEN1 ) = -1

                     NOTETR( 5, NTEN2 ) = NTEN1
                     NOTETR( 6, NTEN2 ) = -1
                     NOTETR( 7, NTEN2 ) = -1
                     NOTETR( 8, NTEN2 ) = -1

                     N1TETS( NS1 ) = NTEN1
                     N1TETS( NS2 ) = NTEN2
                     N1TETS( NS3 ) = NTEN2
                     N1TETS( NS4 ) = NTEN2
                     N1TETS( NS5 ) = NTEN2

C                    CONSTRUCTION DU TABLEAU DES NO DE TETRAEDRES OPPOSES
C                    AUX FACES DE NTE NTEA1 NTEA2 POUR MISE A JOUR
C                    DANS LES TETRAEDRES OPPOSES DE NTEN1 NTEN2
                     NBTOP = NBTECF
                     DO N=1,2
                        DO L=1,4
                           NTOP = NOTETR( 4+L, NTEAR(N) )
                           IF(NTOP .GT. 0 .AND. NTOP .NE. NTE .AND.
     %                        NTOP .NE. NTEA1 .AND. NTOP .NE. NTEA2)THEN
                              NBTOP = NBTOP + 1
                              NOTECF(NBTOP) = NTOP
                           ENDIF
                        ENDDO
                     ENDDO

                     DO L=1,4
                        NTOP = NOTETR( 4+L, NTE )
                        IF( NTOP .GT. 0 .AND. NTOP .NE. NTE .AND.
     %                      NTOP .NE. NTEA1 .AND.
     %                      NTOP .NE. NTEA2 ) THEN
                           NBTOP = NBTOP + 1
                           NOTECF(NBTOP) = NTOP
                        ENDIF
                     ENDDO
C                    NOMBRE DE TETRAEDRES OPPOSES AJOUTES
                     NBTEOP = NBTOP - NBTECF

C                    NTE EST RENDU TETRAEDRE VIDE
                     NOTETR(1,NTE) = 0
                     NOTETR(5,NTE) = N1TEVI
                     N1TEVI = NTE

C                    NTEA1 EST RENDU TETRAEDRE VIDE
                     NOTETR(1,NTEA1) = 0
                     NOTETR(5,NTEA1) = N1TEVI
                     N1TEVI = NTEA1

C                    NTEA2 EST RENDU TETRAEDRE VIDE
                     NOTETR(1,NTEA2) = 0
                     NOTETR(5,NTEA2) = N1TEVI
                     N1TEVI = NTEA2

C                    SUPPRESSION DE NTE NTEAR(1:2) DE NOTECF
C                    AJOUT DE NOTNEW(1:2)
                     DO L=1,NBTECF
                        IF( NOTECF(L) .EQ. NTE ) THEN
                           NOTECF(L) = NTEN1
                           GOTO 27
                        ENDIF
                     ENDDO

 27                  NN = NTEN2
                     DO LL=1,NBTECF
                        IF( NOTECF(LL) .EQ. NTEA1 ) THEN
                           NOTECF(LL) = NN
                           NN = 0
                           GOTO 28
                        ENDIF
                     ENDDO

 28                  DO LLL=1,NBTECF
                        IF( NOTECF(LLL) .EQ. NTEA2 ) THEN
                           NOTECF(LLL) = NN
                           GOTO 30
                        ENDIF
                     ENDDO

C                    COMPRESSION DE NOTECF
 30                  L = 0
                     DO LL=1,NBTOP
                        N = NOTECF( LL )
                        IF( N .GT. 0 ) THEN
                           IF( NOTETR(1,N) .GT. 0 ) THEN
                              L = L + 1
                              NOTECF( L ) = N
C                             MISE A JOUR DE N1TETS POUR LES TETRAEDRES NOTECF
                              DO LLL=1,4
                                 NS = NOTETR( LLL, N )
                                 N1TETS( NS ) = N
                              ENDDO
                           ENDIF
                        ENDIF
                     ENDDO
                     NBTOP  = L
                     NBTECF = L - NBTEOP

C                    MISE A JOUR DES TETRAEDRES OPPOSES
                     CALL MJOPTE( NBTOP,  NOTECF,  N1TETS,
     %                            NOTETR, NUDTETR, N1TEVI,
     %                            PTXYZD, NBFANR )
C                    NBFANR : NOMBRE DE FACES DE TETRAEDRE OPPOSE NON RETROUVE 
C                    LE NUMERO DE TETRAEDRE OPPOSE INCONNU EST MIS A -1
                     PRINT*,'su1teqm: LES TETRAEDRES',NOTNEW,
     %                      ' REMPLACENT LES TETRAEDRES',NTE,NTEAR,
     %                      ' NB FACES FRONTIERE ou NON OPPOSEE=',NBFANR
     %                     ,' NBTECF1=',NBTECF
                     NBTNEW = 2

                     PRINT*,'su1teqm: TETRAEDRE REMPLACANT(',NTEN1,'):',
     %                      (NOTETR(M,NTEN1),M=1,8),' V=',V1,' Q=',Q1
                     PRINT*,'su1teqm: TETRAEDRE REMPLACANT(',NTEN2,'):',
     %                      (NOTETR(M,NTEN2),M=1,8),' V=',V2,' Q=',Q2

C                    TRACE DE NOTECF(1:NBTECF) ET LEURS TETRAEDRES OPPOSES
ccc                     TRACTE = .TRUE.
      KTITRE='       TETRAEDRES apres SUPPRESSION du TETRAEDRE MEDIOCRE 
     %      '
                     WRITE(KTITRE(1:6),  '(I6)') NBTECF
                     WRITE(KTITRE(59:64),'(I6)') NTE
                     CALL TRFETO13(KTITRE, PTXYZD, NBTECF,NOTECF,NOTETR)

                     GOTO 9999

                  ENDIF

 60            ENDDO

            ENDIF
         ENDDO

 100  ENDDO

 9999 TRACTE = TRACTE0
      RETURN
      END
