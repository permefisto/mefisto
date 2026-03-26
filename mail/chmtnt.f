      SUBROUTINE CHMTNT( INFACO, MXFACO, LEFACO, IVOLTE, NVOLTE,
     %                   PTXYZD, NTE0,   NAR,
     %                   N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                   MXTE1A, NBTE1A, NOTE1A,
     %                   MXT2M4, NBT2M4, NOT2M4, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EVENTUEL CHANGEMENT DE m>2 TETRAEDRES AYANT UNE ARETE COMMUNE
C -----    EN 2m-4 TETRAEDRES => L'ARETE NAR EST ECHANGEE EN m-2 FACES
C          SI VOLUMES mT et 2m-4T EGAUX et CHANGEMENT SI et SEULEMENT SI
C          MIN( QUALITE mT INITIAUX ) < MIN( QUALITE DES 2m-4T )

C ENTREES:
C --------
C INFACO : =1 PAS DE mT->2m-4T SI L'UNE DES 3 FACES QUI DISPARAITRAIENT
C             DANS mT->2m-4T EST DANS LEFACO et CHOIX DE LA TRIANGULATION
C             DU POLYGONE POUR OBTENIR LE PLUS DE FACES LEFACO
C          =0 SI PAS DE CONTROLE SUR LA FACE COMMUNE ET SON APPARTENANCE
C             A LEFACO  (LE TABLEAU LEFACO PEUT NE PAS EXISTER)
C MXFACO : NOMBRE MAXIMAL DECLARABLE DE FACES DU CONTOUR
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1 , VOLUME2 DE LA FACE
C          NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C          ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C          => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C          LEFACO(9,*)  -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C          HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C          LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C          NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C          SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C          NF  = LEFACO( 9, NF )  ...
C          11: >0  NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE,
C              =0  SINON
cccC          12: = NO FACEOC DE 1 A NBFACES D'OC

C IVOLTE : 0 PAS       DE TABLEAU NVOLTE A L'APPEL de chmtnt
C          1 EXISTENCE DU TABLEAU NVOLTE A L'APPEL de chmtnt
C NVOLTE : NUMERO DU VOLUME (1 A NBVOPA) DE CHAQUE TETRAEDRE
C          -1 SI VOLUME INCONNU

C NTE0   : NUMERO NOTETR DU PREMIER TETRAEDRE D'ARETE NAR
C NAR    : NUMERO DE 1 A 6 DE L'ARETE COMMUNE DES 3 TETRAEDRES DANS NTE0
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C MXTE1A : NOMBRE MAXIMAL DECLARE D'ENTIERS DU TABLEAU NOTE1A
C          >=NOMBRE MAXIMAL DE TETRAEDRES CONTENANT UN SOMMET
C MXT2M4 : >=20 NOMBRE MAXIMAL DECLARE D'ENTIERS DU TABLEAU NOT2M4 POUR m<=8

C MODIFIES:
C ---------
C N1TETS : N1TETS(I) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET I
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C N1TEVI : NUMERO DU 1 PREMIER TETRAEDRE VIDE DANS LE TABLEAU NOTETR
C          LE CHAINAGE DES TETRAEDRES VIDES SE FAIT SUR NOTETR(5,.)
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE

C SORTIES:
C --------
C NBTE1A : m LE NOMBRE DES TETRAEDRES AUTOUR DE L'ARETE NAR
C NOTE1A : NUMERO NOTETR DES m TETRAEDRES AUTOUR DE L'ARETE NAR
C          DETRUITS SI L'ECHANGE mT -> 2m-4T EST REALISE
C NBT2M4 : 2m-4 LE NOMBRE DE TETRAEDRES REMPLACANT LES NBTE1A TETRAEDRES
C NOT2M4 : NUMERO NOTETR DES 2m-4 TETRAEDRES QUI REMPLACENT NOTE1A
C IERR   : 0 SI PAS D'ERREUR ET mT -> 2m-4T REALISE
C          1 SI LE NOMBRE DE TETRAEDRES D'ARETE COMMUNE NAR EST <=2 ou >10
C            => PAS DE TRANSFORMATION mT=>2m-4T
C          2 SI SATURATION DU TABLEAU NOTETR
C            => PAS DE TRANSFORMATION mT=>2m-4T
C          3 SI LE VOLUME 2m-4T et mT SONT DIFFERENTS
C            OU SI UN DES 3 TETRAEDRES EST DE VOLUME NUL OU NEGATIF
C            OU SI MIN( QUALITE mT INITIAUX ) >= MIN( QUALITE DES 2m-4T )
C            => PAS DE TRANSFORMATION mT=>2m-4T
C          4 VOLUME DES TETRAEDRES DE L'ARETE TROP PETIT DEVANT
C            LEUR VOLUME IDEAL
C          5 SI L'UNE DES 3 FACES APPELEES A DISPARAITRE APPARTIENT
C            A LEFACO ET PAS DE mT->2m-4T
C          6 SI NOTETR(1,NTE0)<=0
C          7 SI ERREUR DANS L'EXECUTION DE TETR1A->TETR1S
C          8 SI SATURATION DU TABLEAU NOT2M4
C          9 SI NTE0 NON RETROUVE DANS NOTE1A(NBTE1A)
C         10 SI NBTE1A>MXTE1A APRES APPEL DE tetr1a
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET  Saint PIERRE DU PERRAY           Octobre 2018
C....................................................................012
      include"./incl/gsmenu.inc"
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
      LOGICAL                          TRACTE0
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  PTXYZD(1:4,1:*), VOLTET, V, V0, V1,
     %                  VOLID, VOLIDT, QD, QDMIN, QDTMAX
      DOUBLE PRECISION  ARMIN, ARMAX, SURFTR(4)
      REAL              Q
      INTEGER           LEFACO(1:11,0:MXFACO)
      INTEGER           NOTETR(8,*), N1TETS(*), NOT2M4(MXT2M4)
      INTEGER           NOSOTR(3), NOTE1A(MXTE1A), NVOLTE(*),
     %                  NOSOCF(32), NOSOCF1(32), NOTE1A1(32)
      CHARACTER*128     KTITRE

C     NO DES 2 SOMMETS EXTREMITES DES ARETES DU TETRAEDRE
      INTEGER           NOSOARTE(2,6)
      DATA              NOSOARTE / 1,2, 2,3, 3,1, 4,1, 4,2, 4,3 /
C     NO DES 3 SOMMETS DES 4 FACES DU TETRAEDRE
      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     NO DES 2 FACES DU TETRAEDRE D'ARETE COMMUNE K
      INTEGER           NOFAARTE(2,6)
      DATA              NOFAARTE / 1,4,  1,2,  1,3,  3,4,  2,4,  2,3 /

      IF( NOTETR(1,NTE0) .LE. 0 ) THEN
         IERR = 6
         RETURN
      ENDIF

      IERR    = 0
      TRACTE0 = TRACTE
      NBT2M4  = 0

C     NUMERO PTXYZD DES 2 SOMMETS DE L'ARETE NAR
      NS1 = NOTETR( NOSOARTE( 1, NAR ), NTE0 )
      NS2 = NOTETR( NOSOARTE( 2, NAR ), NTE0 )

C     L'ARETE NS1-NS2 PERFORE DES FACES COMMUNE A m>2 TETRAEDRES
C     m TETRAEDRES => 2m-4 TETRAEDRES D'ARETE COMMUNE NS1-NS2
      CALL TETR1A( NS1,    NS2,    N1TETS, NOTETR,
     %             NBTE1A, MXTE1A, NOTE1A, IERR )

      IF( NBTE1A .GE. MXTE1A ) THEN
C        TRACE DES TETRAEDRES DU DEBORDEMENT DU TABLEAU NOTE1A
         TRACTE0 = TRACTE
         TRACTE = .TRUE.
         KTITRE ='chmtnt: DEBORDEMENT DU TABLEAU NOTE1A POURQUOI ???'
         CALL TRFETO15( KTITRE, PTXYZD, NBTE1A, NOTE1A, NOTETR,
     %                  NS1, NS2 )
         TRACTE = TRACTE0
      ENDIF

      IF( IERR .NE. 0 ) THEN
         IERR = 7
         GOTO 9999
      ENDIF

      IF( NBTE1A .LE. 2 .OR. NBTE1A .GT. 10 ) THEN
C        m>10 PEUT ETRE UNE SITUATION TROP COMPLIQUEE
C             POUR ETRE AMELIOREE AVEC UN COUT MODERE
ccc         print*,'chmtnt: NBTE1A=',NBTE1A,' => IERR=1 en SORTIE'
         IERR = 1
         GOTO 9999
      ENDIF

ccc      IF( NBTE1A .GE. 10 ) THEN
ccc         TRACTE0 = TRACTE
ccc         TRACTE  = .TRUE.
ccc         KTITRE='chmtnt1:TETRAEDRE               NBTE1A=       NS1=     
ccc     %      NS2=          '
ccc         WRITE( KTITRE(19:28),'(I10)' ) NTE0
ccc         WRITE( KTITRE(40:42),'(I3)'  ) NBTE1A
ccc         WRITE( KTITRE(51:60),'(I10)' ) NS1
ccc         WRITE( KTITRE(66:75),'(I10)' ) NS2
ccc         CALL SANSDBL( KTITRE, NBC )
ccc         CALL TRFETO15( KTITRE(1:NBC),  PTXYZD,
ccc     %                  NBTE1A, NOTE1A, NOTETR, NS1, NS2 )
ccc         TRACTE = TRACTE0
ccc      ENDIF

C     SI LE NO DE VOLUME DES TETRAEDRES (TABLEAU NVOLTE) EST PRESENT
C     ALORS TOUS LES TETRAEDRES DOIVENT ETRE DANS UN MEME VOLUME
C     --------------------------------------------------------------
      NOVOLU = 0
      IF( IVOLTE .NE. 0 ) THEN
         NOVOLU = NVOLTE( NOTE1A(1) )
         DO M = 2, NBTE1A
            IF( NVOLTE( NOTE1A(M) ) .NE. NOVOLU ) THEN
               IERR = 5
               GOTO 9999
            ENDIF
         ENDDO
      ENDIF

C     RECHERCHE DE NTE0 DANS NOTE1A
      DO N = 1, NBTE1A
         IF( NTE0 .EQ. NOTE1A( N ) ) GOTO 5
      ENDDO
      IERR = 9
      GOTO 9999

C     NTE0 EST MIS EN PREMIERE POSITION DE NOTE1A
 5    IF( N .NE. 1 ) THEN
         MM          = NOTE1A( 1 )
         NOTE1A( 1 ) = NOTE1A( N )
         NOTE1A( N ) = MM
      ENDIF

C     TRAITEMENT DE LA FACE 1 D'ARETE NAR DE NTE0
C     -------------------------------------------
C     RECHERCHE DU SOMMET NS3 DE LA FACE 1 D'ARETE NAR DE SOMMETS NS1 NS2
      NFAC1 = NOFAARTE( 1, NAR )
      DO N = 1,3
         NS3 = NOTETR( NOSOFATE( N, NFAC1 ), NTE0 )
         IF( NS3 .NE. NS1 .AND. NS3 .NE. NS2 ) GOTO 10
      ENDDO
 10   NOSOCF(1) = NS3

C     NTEOP TETRAEDRE OPPOSE A LA FACE NFAC1 DE NTE0
      NTEOP = NOTETR( 4+NFAC1, NTE0 )
      IF( NTEOP .LE. 0 ) THEN
C        LES NBTE1A N'ENTOURENT PAS L'ARETE QUI N'EST PAS SUPPRIMABLE
         IERR = 10
         GOTO 9999
      ENDIF

C     PROTECTION DES FACES DE LEFACO QUI DISPARAITRAIENT DANS mT->2m-4T
C     FACE NFAC1 DE NTE0: NS1-NS2-NS3
C     LA FACE DE SOMMETS NOSOTR(1:3) EXISTE-T-ELLE DEJA DANS LEFACO?
C     HACHAGE DANS LEFACO AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
      IF( INFACO .NE. 0 ) THEN
         NOSOTR(1) = NS1
         NOSOTR(2) = NS2
         NOSOTR(3) = NS3
         CALL NULETR( NOSOTR, MXFACO, LEFACO,  NFLEFA )
         IF( NFLEFA .GT. 0 ) THEN
C           OUI: FACE DANS LEFACO -> POUR LA CONSERVER PAS D'ECHANGE
            IERR = 5
            GOTO 9999
         ENDIF
      ENDIF

C     BOUCLE SUR LES TETRAEDRES QUI ENROULENT L'ARETE NAR DE NTE0
C     EXISTE IL UNE RAISON DE NE PAS ELIMINER L'ARETE?
C     -----------------------------------------------------------
      DO 50 M = 1, NBTE1A-1

         NTECF = NOTE1A( M )

C        RECHERCHE DU SOMMET NON NS1 NS2 NOSOCF(M) DE NTECF
         DO N = 1,4
            NS3 = NOTETR( N, NTECF )
            IF( NS3 .NE. NS1 .AND. NS3 .NE. NS2 .AND.
     %          NS3 .NE. NOSOCF(M) ) GOTO 20
         ENDDO
C        CAS IMPOSSIBLE...
         IERR = 11
         GOTO 9999

 20      M1 = M + 1
         NOSOCF( M1 ) = NS3

C        RECHERCHE DE LA FACE NS1 NS2 NOSOCF(M+1)=NS3
         NOSOTR(1) = NS1
         NOSOTR(2) = NS2
         NOSOTR(3) = NS3

         IF( M .EQ. 1 ) THEN
C           FACE NFAC2 DE NTE0
            NFAC2 = NOFAARTE( 2, NAR )
         ELSE
            CALL NUFATRTE( NOSOTR, NOTETR(1,NTECF), NFAC2 )
            IF( NFAC2 .EQ. 0 ) THEN
C              NOSOTR N'EST PAS UNE FACE DE NTECF
               IERR = 12
               GOTO 9999
            ENDIF
         ENDIF

C        NTEOP TETRAEDRE OPPOSE A LA FACE NFAC2 DE NTECF
         NTEOP = NOTETR( 4+NFAC2, NTECF )
         IF( NTEOP .LE. 0 ) THEN
C           LES NBTE1A N'ENTOURENT PAS L'ARETE QUI N'EST PAS SUPPRIMABLE
            IERR = 13
            GOTO 9999
         ENDIF

         IF( INFACO .NE. 0 ) THEN
            CALL NULETR( NOSOTR, MXFACO, LEFACO,  NFLEFA )
            IF( NFLEFA .GT. 0 ) THEN
C              OUI: FACE DANS LEFACO -> POUR LA CONSERVER PAS D'ECHANGE
               IERR = 5
               GOTO 9999
            ENDIF
         ENDIF

C        RECHERCHE DANS NOTE1A DE NTEOP ET PERMUTATION EVENTUELLE
         DO K = M1, NBTE1A
            IF( NOTE1A( K ) .EQ. NTEOP ) THEN
               IF( K .NE. M1 ) THEN
C                 PERMUTATION DU NO DES TETRAEDRES M1 ET K DE NOTE1A
                  MM           = NOTE1A( M1 )
                  NOTE1A( M1 ) = NOTE1A( K )
                  NOTE1A( K  ) = MM
               ENDIF
               GOTO 50
            ENDIF
         ENDDO

 50   ENDDO

C     LES TETRAEDRES NOTE1A S'ENROULENT ILS TOTALEMENT AUTOUR DE NS1-NS2?
C     -------------------------------------------------------------------
C     VERIFICATION NOSOCF(1) EST IL SOMMET DE NOTE1A( NBTE1A )?
      NT = NOTE1A( NBTE1A )
      DO N = 1, 4
         IF( NOTETR( N, NT ) .EQ. NOSOCF( 1 ) ) GOTO 60
      ENDDO

C     NOTE1A(NBTE1A) NE CONTIENT PAS NOSOCF(1). PAS D'ECHANGE
      IERR = 14
      GOTO 9999

C     RECHERCHE DE LA MEILLEURE TRIANGULATION DU POLYGONE NOSOCF
C     ----------------------------------------------------------
 60   IF( NBTE1A .EQ. 3 ) GOTO 70

C     QUADRANGLE => 2 TRIANGULATIONS
C     PENTAGONE  => 5 TRIANGULATIONS ESSAYEES
C     HEXAGONE   => 6 TRIANGULATIONS ESSAYEES (IL EN EXISTE D'AUTRES...)
C     DECAGONE   =>10 TRIANGULATIONS ESSAYEES (IL EN EXISTE D'AUTRES...)
C     NOMBRE DE TRIANGULATIONS A EXAMINER POUR TROUVER LE
C     MAX DES MIN DE LA QUALITE DES TRIANGULATIONS
      IF( NBTE1A .EQ. 4 ) THEN
         NBTR = 2
      ELSE
         NBTR = NBTE1A
      ENDIF

      QDTMAX = -2D0
      NTRMAX = 0
      NBFAMX = 0
      NTFAMX = 0

      DO K = 1, NBTR

C        NOSOCF EST RENDU PERIODIQUE
         NOSOCF( NBTE1A+K ) = NOSOCF( K )
         QDMIN  = 2D0
         NBFACO = 0

         NS0 = NOSOCF( K )
         DO N = 1, NBTE1A-2

C           QD LA QUALITE DU TRIANGLE NOSOTR
            NOSOTR(1) = NS0
            NOSOTR(2) = NOSOCF(K+N  )
            NOSOTR(3) = NOSOCF(K+N+1)
            CALL QUATRID3S( PTXYZD( 1, NOSOTR(1) ),
     %                      PTXYZD( 1, NOSOTR(2) ),
     %                      PTXYZD( 1, NOSOTR(3) ), QD )
            QDMIN = MIN( QDMIN, QD )

            IF( INFACO .NE. 0 ) THEN
               CALL NULETR( NOSOTR, MXFACO, LEFACO,  NFLEFA )
               IF( NFLEFA .GT. 0 ) THEN
C                 OUI: FACE DANS LEFACO
                  NBFACO = NBFACO + 1
               ENDIF
            ENDIF

         ENDDO

         IF( QDMIN .GT. QDTMAX ) THEN
            QDTMAX = QDMIN
            NTRMAX = K
         ENDIF

         IF( INFACO .NE. 0 ) THEN
            IF( NBFACO .GT. NBFAMX ) THEN
               NBFAMX = NBFACO
               NTFAMX = K
            ENDIF
         ENDIF

      ENDDO

      IF( INFACO .NE. 0 ) THEN
         IF( NBFAMX .GT. 0 ) THEN
C           LA CONSERVATION DES FACES LEFACO L'EMPORTE SUR LA MAXIMISATION
C           DE LA QUALITE DES TRIANGLES ET TETRAEDRES CREES
            NTRMAX = NTFAMX
         ENDIF
      ENDIF

      IF( NTRMAX .GT. 1 ) THEN

C        PERMUTATION CIRCULAIRE DU TABLEAU NOSOCF et NOTE1A
         DO K = 1, NBTE1A
            N = NTRMAX-1+K
            NOSOCF1( K ) = NOSOCF( N )
            IF( N .GT. NBTE1A ) N = N - NBTE1A
            NOTE1A1( K ) = NOTE1A( N )
         ENDDO
         DO K = 1, NBTE1A
            NOSOCF( K ) = NOSOCF1( K )
            NOTE1A( K ) = NOTE1A1( K )
         ENDDO

         NOSOCF( NBTE1A+1 ) = NOSOCF( 1 )

      ENDIF

C     LES NBTE1A TETRAEDRES FORMENT ILS UN ENSEMBLE ETOILE?
C     COMPARAISON DES VOLUMES DES TETRAEDRES AVEC ET SANS L'ARETE NAR DE NTE0
C     -----------------------------------------------------------------------
 70   V0     = 0D0
      VOLIDT = 0D0
      QMIN0  = 2.
      DO N = 1, NBTE1A

         NT = NOTE1A( N )
         CALL QUATETD( PTXYZD(1,NOTETR(1,NT)),
     %                 PTXYZD(1,NOTETR(2,NT)),
     %                 PTXYZD(1,NOTETR(3,NT)),
     %                 PTXYZD(1,NOTETR(4,NT)),
     %                 ARMIN, ARMAX, SURFTR, V, Q )
         V0 = V0 + ABS( V )
         QMIN0 = MIN( QMIN0, Q )

C        LE VOLUME du TETRAEDRE IDEAL NT POUR LA TAILLE SOUHAITEE
C        DE SES ARETES EN CHACUN DE SES 4 SOMMETS
C        ( A DIVISER PAR 6 POUR AVOIR DES TETRAEDRES AU LIEU DE CUBES)
         VOLID = ( ( PTXYZD(4,NOTETR(1,NT))
     %             + PTXYZD(4,NOTETR(2,NT))
     %             + PTXYZD(4,NOTETR(3,NT))
     %             + PTXYZD(4,NOTETR(4,NT)) ) / 4 ) **3
         VOLIDT = VOLIDT + VOLID

      ENDDO

C     LES NBTE1A TETRAEDRES ONT ILS UN TRES FAIBLE VOLUME PAR RAPPORT
C     A LEUR VOLUME IDEAL ?
      IF( V0 .LT. 1D-4 * VOLIDT ) THEN
C        OUI: PAS DE MODIFICATION DES TETRAEDRES D'ARETE NAR DE NTE0
C             POUR NE PAS EMPIRER LES TETRAEDRES PLATS
         IERR = 4
C        RETOUR SANS MODIFICATION DES TETRAEDRES
         GOTO 9999
      ENDIF


C     VOLUME >0 ou <0 DES TETRAEDRES NOUVEAUX SANS L'ARETE NAR NS1-NS2?
      V = VOLTET( PTXYZD( 1, NOSOCF(1) ),
     %            PTXYZD( 1, NOSOCF(2) ),
     %            PTXYZD( 1, NOSOCF(3) ),
     %            PTXYZD( 1, NS1       ) )
      IF( V .LT. 0D0 ) THEN
C        PERMUTATION DE NS1 et NS2 POUR OBTENIR DES VOLUMES POSITIFS
ccc         print*,'chmtnt: ECHANGE NS1=',NS1,' et NS2=',NS2
         MM  = NS1
         NS1 = NS2
         NS2 = MM
      ENDIF

C     TRIANGULATION FAITE A PARTIR DU 1-ER SOMMET DU CF NOSOCF
      V1    = 0D0
      QMIN1 = 2.

      NS0 = NOSOCF( 1 )
      NBT = NBTE1A
      DO N = 1, NBTE1A-2
         NBT = NBT + 1
         CALL QUATETD( PTXYZD( 1, NS0         ),
     %                 PTXYZD( 1, NOSOCF(N+1) ),
     %                 PTXYZD( 1, NOSOCF(N+2) ),
     %                 PTXYZD( 1, NS1         ),
     %                 ARMIN, ARMAX, SURFTR, V, Q )
         V1 = V1 + ABS( V )
         QMIN1 = MIN( QMIN1, Q )

         NBT = NBT + 1
         CALL QUATETD( PTXYZD( 1, NS0         ),
     %                 PTXYZD( 1, NOSOCF(N+2) ),
     %                 PTXYZD( 1, NOSOCF(N+1) ),
     %                 PTXYZD( 1, NS2         ),
     %                 ARMIN, ARMAX, SURFTR, V, Q )
         V1 = V1 + ABS( V )
         QMIN1 = MIN( QMIN1, Q )

      ENDDO


C     COMPARAISON DU VOLUME DES mT et 2m-4T
C     et du MAXIMUM des QUALITES MINIMALES mT ou 2m-4T
C     ------------------------------------------------
      IF( ABS(V0-V1) .GT. V0*1D-8  .OR.  QMIN0 .GE. QMIN1 ) THEN

CCC      WRITE(IMPRIM,10050) V0, V1, QMIN0, QMIN1
CCC10050 FORMAT(' chmtnt: VOLUME mT AVANT',G25.16,' NON EGAL A 2m-4T',
ccc     %G25.16,' ou QMIN0=',F8.5,'>= QMIN1=',F8.5)

         IERR = 3
C        RETOUR SANS MODIFICATION DES TETRAEDRES
         GOTO 9999

      ENDIF


C     m Tetra -> 2m-4 Tetra
C     =====================
      NBT2M4 = 2 * NBTE1A - 4

ccc      IF( NBTE1A .GE. 7 ) THEN
ccc         PRINT 10051, NBTE1A, NBT2M4, V0, V1, QMIN0, QMIN1
ccc10051 FORMAT(' chmtnt: m=',I3,'T->',I3,'T. VOLUME AVANT',G25.16,
ccc     %       ' APRES',G25.16,' QMIN0=',F8.5,'< QMIN1=',F8.5)
ccc      ENDIF

C     RECHERCHE DE 2m-4 TETRAEDRES VIDES POUR LES CREER
C     -------------------------------------------------
      IF( NBT2M4 .GT. MXT2M4 ) THEN
C        SATURATION DU TABLEAU NOT2M4
         PRINT*,'chmtnt: SATURATION du TABLEAU NOT2M4 MXT2M4=',MXT2M4
         IERR = 8
         GOTO 9999
      ENDIF

C     RESERVATION DES NBT2MA SOUS-TETRAEDRES
      DO N = 1, NBT2M4

         IF( N1TEVI .LE. 0 ) THEN
C           SATURATION DES TETRAEDRES NOTETR
            NBLGRC(NRERR) = 1
            KERR(1) = 'SATURATION DES TETRAEDRES N1TEVI='
            L = NUDCNB( KERR(1) )
            WRITE(KERR(1)(L+1:L+10),'(I10)') N1TEVI
            CALL LEREUR
            IERR = 2
            GOTO 9999
         ENDIF
         NOT2M4(N) = N1TEVI

         IF( IVOLTE .NE. 0 ) THEN
C           NUMERO DE VOLUME DES SOUS-TETRAEDRES
            NVOLTE( N1TEVI ) = NOVOLU
         ENDIF

C        MISE A JOUR DU DERNIER TETRAEDRE OCCUPE
         NUDTETR = MAX( NUDTETR, N1TEVI )

C        MISE A JOUR DU 1-ER TETRAEDRE VIDE
         N1TEVI  = NOTETR( 5, N1TEVI )

      ENDDO

C     INITIALISATION DES NBT2M4 SOUS-TETRAEDRES
C     -----------------------------------------
C     LES 2 SOUS-TETRAEDRES QUI PRECEDENT LES 2 PREMIERS SOUS-TETRAEDRES

C     LA FACE 4: NS1 NS0 NS3  DU PREMIER NTE1 A UN TETRAEDRE OPPOSE
C                A UNE FACE D'UN TETRAEDRE NOTE1A
      NTE1 = NOT2M4( 1 )
      NS0  = NOSOCF( 1 )
      NS3  = NOSOCF( 2 )

      NOSOTR(1) = NS1
      NOSOTR(2) = NS0
      NOSOTR(3) = NS3
      CALL TRI3NO( NOSOTR, NOSOTR )
      DO NN=1,NBTE1A
         NT1A = NOTE1A(NN)
         CALL NO1F1T( NOSOTR, NOTETR(1,NT1A), NF1A )
         IF( NF1A .GT. 0 ) THEN
C           LE TETRAEDRE OPPOSE
            NTE1PR = NOTETR( 4+NF1A, NT1A )
            IF( NTE1PR .GT. 0 .AND. NOTETR(1,NTE1PR) .GT. 0 ) THEN
               CALL NO1F1T( NOSOTR, NOTETR(1,NTE1PR), NFOP )
               IF( NFOP .GT. 0 ) THEN
                  NOTETR( 4+NFOP, NTE1PR ) = NTE1
               ENDIF
            ELSE
               NTE1PR = 0
            ENDIF
            GOTO 71
         ENDIF
      ENDDO
      print*,'chmtnt: Pb face',nosotr,' NON RETROUVEE dans NOTE1A'
      NTE1PR = 0

C     LA FACE 3: NS3 NS2 NS0  DU PREMIER NTE2 A UN TETRAEDRE OPPOSE
C                A UNE FACE D'UN TETRAEDRE NOTE1A
 71   NTE2 = NOT2M4( 2 )
      NOSOTR(1) = NS2
      NOSOTR(2) = NS0
      NOSOTR(3) = NS3
      CALL TRI3NO( NOSOTR, NOSOTR )
      DO NN=1,NBTE1A
         NT2A = NOTE1A(NN)
         CALL NO1F1T( NOSOTR, NOTETR(1,NT2A), NF2A )
         IF( NF2A .GT. 0 ) THEN
C           LE TETRAEDRE OPPOSE
            NTE2PR = NOTETR( 4+NF2A, NT2A )
            IF( NTE2PR .GT. 0 .AND. NOTETR(1,NTE2PR) .GT. 0 ) THEN
               CALL NO1F1T( NOSOTR, NOTETR(1,NTE2PR), NFOP )
               IF( NFOP .GT. 0 ) THEN
                  NOTETR( 4+NFOP, NTE2PR ) = NTE2
               ENDIF
            ELSE
               NTE2PR = 0
            ENDIF
            GOTO 72
         ENDIF
      ENDDO
      print*,'chmtnt: PB face',nosotr,' NON RETROUVEE dans NOTE1A'
      NTE2PR = 0


C     NOMBRE DE TRIANGLES (QUI REMPLACENT NS1-NS2) A TRAITER
C     CHAQUE TRIANGLE GENERE AVEC NS1 NS2 2 TETRAEDRES
C     ======================================================
 72   NBTR  = NBTE1A - 2

      NB2TE = 0
      DO N = 1, NBTR

C        LES 2 SOMMETS DE L'ARETE OPPOSE A L'ARETE NS1-NS2
         NS3 = NOSOCF(N+1)
         NS4 = NOSOCF(N+2)

C        LES 2 SOUS-TETRAEDRES DE FACE COMMUNE NS0-NS3-NS4 A CREER
         NB2TE = NB2TE + 1
         NTE1  = NOT2M4( NB2TE )

         NB2TE = NB2TE + 1
         NTE2  = NOT2M4( NB2TE )

C        LES SOMMETS DES 2 SOUS-TETRAEDRES
C        ---------------------------------
         NOTETR(1,NTE1) = NS0
         NOTETR(2,NTE1) = NS3
         NOTETR(3,NTE1) = NS4
         NOTETR(4,NTE1) = NS1

         NOTETR(1,NTE2) = NS0
         NOTETR(2,NTE2) = NS4
         NOTETR(3,NTE2) = NS3
         NOTETR(4,NTE2) = NS2

C        LES COUPLES DE TETRAEDRES OPPOSES PAR UNE FACE
C        ----------------------------------------------
C        LES 2 SOUS-TETRAEDRES SUIVANTS DE NTE1 et NTE2
         IF( N .LT. NBTR ) THEN

C           LES 2 TETRAEDRES NTE1 NTE2 SUIVANTS CEUX NTE1 NTE2 A INITIALISER
            NTE1SU = NOT2M4( 2*N+1 )
            NTE2SU = NOT2M4( 2*N+2 )

         ELSE

C           LA FACE NS4 NS1 NS0 DU DERNIER NTE1 A UN TETRAEDRE OPPOSE
C           A UNE FACE D'UN TETRAEDRE NOTE1A
            NOSOTR(1) = NS4
            NOSOTR(2) = NS1
            NOSOTR(3) = NS0
            CALL TRI3NO( NOSOTR, NOSOTR )
            DO NN=1,NBTE1A
               NT1A = NOTE1A(NN)
               CALL NO1F1T( NOSOTR, NOTETR(1,NT1A), NF1A )
               IF( NF1A .GT. 0 ) THEN
C                 LE TETRAEDRE OPPOSE
                  NTE1SU = NOTETR( 4+NF1A, NT1A )
                  IF( NTE1SU .GT. 0 .AND. NOTETR(1,NTE1SU) .GT. 0 ) THEN
                     CALL NO1F1T( NOSOTR, NOTETR(1,NTE1SU), NFOP )
                     IF( NFOP .GT. 0 ) THEN
                        NOTETR( 4+NFOP, NTE1SU ) = NTE1
                     ENDIF
                  ELSE
                     NTE1SU = 0
                  ENDIF
                  GOTO 73
               ENDIF
            ENDDO
            print*,'chmtnt: PB face',nosotr,' NON RETROUVEE dans NOTE1A'
            NTE1SU = 0


C           LA FACE NS2 NS0 NS4 DU DERNIER NTE2 A UN TETRAEDRE OPPOSE
C           A UNE FACE D'UN TETRAEDRE NOTE1A
 73         NOSOTR(1) = NS2
            NOSOTR(2) = NS0
            NOSOTR(3) = NS4
            CALL TRI3NO( NOSOTR, NOSOTR )
            DO NN=1,NBTE1A
               NT2A = NOTE1A(NN)
               CALL NO1F1T( NOSOTR, NOTETR(1,NT2A), NF2A )
               IF( NF2A .GT. 0 ) THEN
C                 LE TETRAEDRE OPPOSE
                  NTE2SU = NOTETR( 4+NF2A, NT2A )
                  IF( NTE2SU .GT. 0 .AND. NOTETR(1,NTE2SU) .GT. 0 ) THEN
                     CALL NO1F1T( NOSOTR, NOTETR(1,NTE2SU), NFOP )
                     IF( NFOP .GT. 0 ) THEN
                        NOTETR( 4+NFOP, NTE2SU ) = NTE2
                     ENDIF
                  ELSE
                     NTE2SU = 0
                  ENDIF
                  GOTO 74
               ENDIF
            ENDDO
            print*,'chmtnt: PB face',nosotr,' NON RETROUVEE dans NOTE1A'
            NTE2SU = 0

         ENDIF

C        LE TETRAEDRE OPPOSE AUX 4 FACES DU SOUS-TETRAEDRE NTE1
C        ------------------------------------------------------
C        LA FACE 1 DE NTE1 EST LA FACE 1 DE NTE2
 74      NOTETR( 5, NTE1 ) = NTE2

C        LA FACE 2 DE NTE1 A UN TETRAEDRE OPPOSE A LA FACE D'UN TETRAEDRE NOTE1A
         NOSOTR(1) = NS3
         NOSOTR(2) = NS4
         NOSOTR(3) = NS1
         CALL TRI3NO( NOSOTR, NOSOTR )
         DO NN=1,NBTE1A
            NT1A = NOTE1A(NN)
            CALL NO1F1T( NOSOTR, NOTETR(1,NT1A), NF1A )
            IF( NF1A .GT. 0 ) THEN
C              LE TETRAEDRE OPPOSE
               NTOP = NOTETR( 4+NF1A, NT1A )
               NOTETR( 6, NTE1 ) = NTOP
               IF( NTOP .GT. 0 .AND. NOTETR(1,NTOP) .GT. 0 ) THEN
                  CALL NO1F1T( NOSOTR, NOTETR(1,NTOP), NFOP )
                  IF( NFOP .GT. 0 ) THEN
                     NOTETR( 4+NFOP, NTOP ) = NTE1
                  ENDIF
               ENDIF
               GOTO 75
            ENDIF
         ENDDO

C        LA FACE 3 DE NTE1: NS4 NS1 NS0 EST LA FACE 3 DU NTE1 SUIVANT
 75      NOTETR( 7, NTE1 ) = NTE1SU

C        LA FACE 4 DE NTE1: NS1 NS0 NS3  EST LA FACE 3 DU PRECEDENT NTE1
         NOTETR( 8, NTE1 ) = NTE1PR

C        LE TETRAEDRE OPPOSE AUX 4 FACES DU SOUS-TETRAEDRE NTE2
C        ------------------------------------------------------
C        LA FACE 1 DE NTE2 EST LA FACE 1 DE NTE1
         NOTETR( 5, NTE2 ) = NTE1

C        LA FACE 2 DE NTE2 A UN TETRAEDRE OPPOSE A LA FACE D'UN TETRAEDRE NOTE1A
         NOSOTR(1) = NS4
         NOSOTR(2) = NS3
         NOSOTR(3) = NS2
         CALL TRI3NO( NOSOTR, NOSOTR )
         DO NN=1,NBTE1A
            NT2A = NOTE1A(NN)
            CALL NO1F1T( NOSOTR, NOTETR(1,NT2A), NF2A )
            IF( NF2A .GT. 0 ) THEN
C              LE TETRAEDRE OPPOSE
               NTOP = NOTETR( 4+NF2A, NT2A )
               NOTETR( 6, NTE2 ) = NTOP
               IF( NTOP .GT. 0 .AND. NOTETR(1,NTOP) .GT. 0 ) THEN
                  CALL NO1F1T( NOSOTR, NOTETR(1,NTOP), NFOP )
                  IF( NFOP .GT. 0 ) THEN
                     NOTETR( 4+NFOP, NTOP ) = NTE2
                  ENDIF
               ENDIF
               GOTO 80
            ENDIF
         ENDDO

C        LA FACE 3 DE NTE2: NS3 NS2 NS0  EST LA FACE 4 DU PRECEDENT NTE2
 80      NOTETR( 7, NTE2 ) = NTE2PR

C        LA FACE 4 DE NTE2: NS2 NS0 NS4 EST LA FACE 3 DU NTE2 SUIVANT
         NOTETR( 8, NTE2 ) = NTE2SU

C        LE TETRAEDRE OPPOSE A LA FACE 4 DU FUTUR NTE1
         NTE1PR = NTE1
C        LE TETRAEDRE OPPOSE A LA FACE 3 DU FUTUR NTE2
         NTE2PR = NTE2

C        NO 1 TETRAEDRE DES SOMMETS NS3 et NS4
         N1TETS( NS3 ) = NTE1
         N1TETS( NS4 ) = NTE2

      ENDDO

C     MISE A JOUR D'UN TETRAEDRE CONTENANT CHAQUE SOMMET
      N1TETS( NS1 ) = NTE1
      N1TETS( NS2 ) = NTE2
      N1TETS( NS0 ) = NTE2

cccC     POUR VOIR LES PLUS NOMBREUX SOUS-TETRAEDRES CREES
ccc      IF( NBT2M4 .GE. 14 ) THEN
ccc         TRACTE0 = TRACTE
ccc         TRACTE  = .TRUE.
ccc         KTITRE='chmtnt: TETRAEDRE               NBT2M4=       NS1=     
ccc     %      NS2=            QMIN0=                QMIN1=        '
ccc         WRITE( KTITRE(19:28),'(I10)' ) NOT2M4(1)
ccc         WRITE( KTITRE(40:42),'(I3)'  ) NBT2M4
ccc         WRITE( KTITRE(51:60),'(I10)' ) NS1
ccc         WRITE( KTITRE(66:75),'(I10)' ) NS2
ccc         WRITE( KTITRE(84:97),  '(G14.6)' ) QMIN0
ccc         WRITE( KTITRE(112:125),'(G14.6)' ) QMIN1
ccc         CALL SANSDBL( KTITRE, NBC )
ccc         PRINT*, KTITRE(1:NBC)
ccc         CALL TRFETO15( KTITRE(1:NBC),  PTXYZD,
ccc     %                  NBT2M4, NOT2M4, NOTETR, NS1, NS2 )
ccc         TRACTE = TRACTE0
ccc      ENDIF

C     DESTRUCTION DES NBTE1A TETRAEDRES NOTE1A
C     QUI REDEVIENNENT VIDES DANS NOTETR
C     ----------------------------------------
      DO N = 1, NBTE1A

C        LE TETRAEDRE A DETRUIRE
         NT = NOTE1A( N )

C        SUPPRIMER EVENTUELLEMENT LES 4 FACES DU TETRAEDRE NT
C        DANS LE TABLEAU LEFACO
         CALL SUTELEFA( NT, NOTETR, INFACO, MXFACO, LEFACO )

C        REMISE A ZERO ET CHAINAGE TETRAEDRE VIDE
         NOTETR( 1, NT ) = 0
         NOTETR( 5, NT ) = N1TEVI
         N1TEVI = NT

         IF( IVOLTE .NE. 0 ) THEN
C           NUMERO INCONNU DE VOLUME DU TETRAEDRE INACTIF
            NVOLTE( NT ) = -1
         ENDIF

      ENDDO


C     CONTROLE DES TETRAEDRES CREES
C     RECHERCHE DE TETRAEDRE CREE AVEC UN TETRAEDRE OPPOSE TRIPLE
C     ET SON TRAITEMENT PAR SUPPRESSION DES TETRAEDRES
C     -----------------------------------------------------------
      DO 90 N=1,NBT2M4

         NTE1 = NOT2M4( N )
         IF( NOTETR(1,NTE1) .EQ. 0 ) GOTO 90

C        RECHERCHE D'UN TETRAEDRE OPPOSE TRIPLE
         CALL TE3TOP( NTE1, NOTETR, NTE1OP3, NTE1OP4 )
         IF( NTE1OP3 .LT. 0 ) GOTO 90
         IF( NTE1OP4 .LE. 0 ) GOTO 90

C        NTE1OP3 A T IL AUSSI UN TETRAEDRE OPPOSE TRIPLE?
         NTE2 = NTE1OP3
         CALL TE3TOP( NTE2, NOTETR, NTE2OP3, NTE2OP4 )
         IF( NTE2OP3 .LT. 0 ) GOTO 90
         IF( NTE2OP4 .LE. 0 ) GOTO 90

ccc         tracte=.true.
         ktitre ='chmtnt: traitement1 des tetraedres opposes triples'
         print*,ktitre
         print*,'chmtnt: NTE0=',NTE0,' NS1=',NS1,' NS2=',NS2,
     %          ' NOSOCF=',(NOSOCF(m),m=1,nbte1a)

         print*,'chmtnt: tableau note1a'
         do k = 1, nbte1a
            ntt = note1a( k )
            print*,'chmtnt: notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)
         enddo

         print*,'chmtnt: tableau not2m4'
         do k=1,nbt2m4
            ntt=not2m4(k)
            print*,'chmtnt: notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)
         enddo
         print*

       ntt=nte1
       print*,'chmtnt nte1   : notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)
       ntt=nte1op3
       print*,'chmtnt nte1op3: notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)
       ntt=nte1op4
       print*,'chmtnt nte1op4: notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)

       ntt=nte2
       print*,'chmtnt nte2   : notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)
       ntt=nte2op3   
       print*,'chmtnt nte2op3: notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)
       ntt=nte2op4
       print*,'chmtnt nte2op4: notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)
       print*

         nbtt = nbt2m4 + 1
         not2m4(nbtt) = NTE1OP3
         nbtt = nbtt + 1
         not2m4(nbtt) = NTE1OP4
         nbtt = nbtt + 1
         not2m4(nbtt) = NTE2OP3
         nbtt = nbtt + 1
         not2m4(nbtt) = NTE2OP4
         call trfeto15( ktitre, ptxyzd, nbtt, not2m4, notetr,
     %                  ns1, ns2 )

C        ECHANGE DES TETRAEDRES OPPOSES UNIQUES DE NTE1 NTE2
C        ET SUPPRESSION DE NTE1OP3 et NTE2OP3
         DO K=5,8
C           RECHERCHE DE LA FACE DE NTE1OP4 DE TETRAEDRE OPPOSE NTE1
            IF( NOTETR(K,NTE1OP4) .EQ. NTE1 ) THEN
                NOTETR(K,NTE1OP4) = NTE2OP4
C               SUPPRESSION DE NTE1OP3
                NOTETR( 1, NTE1OP3 ) = 0
                NOTETR( 5, NTE1OP3 ) = N1TEVI
                N1TEVI = NTE1OP3
                IF( IVOLTE .NE. 0 ) THEN
C                  NUMERO INCONNU DE VOLUME DU TETRAEDRE INACTIF
                   NVOLTE( NTE1OP3 ) = -1
                ENDIF
                DO M=1,4
                   N1TETS( NOTETR(M,NTE1OP3) ) = -1
                ENDDO
                DO M=1,4
                   N1TETS( NOTETR(M,NTE1OP4) ) = NTE1OP4
                ENDDO
                GOTO 85
            ENDIF
         ENDDO

 85      DO L=5,8
C           RECHERCHE DE LA FACE DE NTE2OP4 DE TETRAEDRE OPPOSE NTE2
            IF( NOTETR(L,NTE2OP4) .EQ. NTE2 ) THEN
                NOTETR(L,NTE2OP4) = NTE1OP4
C               SUPPRESSION DE NTE2OP3
                NOTETR( 1, NTE2OP3 ) = 0
                NOTETR( 5, NTE2OP3 ) = N1TEVI
                N1TEVI = NTE2OP3
                IF( IVOLTE .NE. 0 ) THEN
C                  NUMERO INCONNU DE VOLUME DU TETRAEDRE INACTIF
                   NVOLTE( NTE2OP3 ) = -1
                ENDIF
                DO M=1,4
                   N1TETS( NOTETR(M,NTE2OP3) ) = -1
                ENDDO
                DO M=1,4
                   N1TETS( NOTETR(M,NTE2OP4) ) = NTE2OP4
                ENDDO
                GOTO 88
            ENDIF
         ENDDO

C        AFFICHAGE D'INFORMATION A SUPPRIMER ENSUITE
 88      print*

         ktitre ='chmtnt: traitement2 des tetraedres opposes triples'
         print*,ktitre
         print*,'chmtnt: NTE0=',NTE0,' NS1=',NS1,' NS2=',NS2,
     %          ' NOSOCF=',(NOSOCF(m),m=1,nbte1a)

         print*,'chmtnt: tableau note1a'
         do k = 1, nbte1a
            ntt = note1a( k )
            print*,'chmtnt: notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)
         enddo

         print*,'chmtnt: tableau not2m4'
         do k=1,nbt2m4
            ntt=not2m4(k)
            print*,'chmtnt: notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)
         enddo
         print*

       ntt=nte1
       print*,'chmtnt nte1   : notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)
       ntt=nte1op3
       print*,'chmtnt nte1op3: notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)
       ntt=nte1op4
       print*,'chmtnt nte1op4: notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)

       ntt=nte2
       print*,'chmtnt nte2   : notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)
       ntt=nte2op3   
       print*,'chmtnt nte2op3: notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)
       ntt=nte2op4
       print*,'chmtnt nte2op4: notetr(',ntt,')=',(notetr(kk,ntt),kk=1,8)
       print*

ccc         tracte=.true.
         ktitre ='chmtnt: traitement2 des tetraedres opposes triples'
         call trfeto15( ktitre, ptxyzd, nbt2m4, not2m4, notetr,
     %                  ns1, ns2 )

 90   ENDDO

      DO 95 N=1,NBT2M4

C        LE TETRAEDRE CREE N
         NTE1 = NOT2M4( N )
         IF( NOTETR(1,NTE1) .EQ. 0 ) GOTO 95

C        AJOUTER EVENTUELLEMENT LES 4 FACES DU TETRAEDRE NTE1
C        DANS LE TABLEAU LEFACO
         CALL AJTELEFA( NTE1, NOTETR, INFACO, MXFACO, LEFACO )

 95   ENDDO

C     SORTIE SANS ERREUR
      IERR = 0

 9999 TRACTE = TRACTE0
      RETURN
      END
