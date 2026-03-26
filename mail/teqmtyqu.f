      SUBROUTINE TEQMTYQU( PTXYZD,  NTEQM,  N1TEVI, NOTETR, N1TETS,
     %                     NUDTETR, IVOLTE, NVOLTE,
     %                     MXTE1A,  NBTE1A, NOTE1A, MODIFT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     TRAITEMENT DU TETRAEDRE NTEQM DE QUALITE MEDIOCRE
C ----     DE TYPE QUADRANGLE PAR DECOMPOSITION DE TETRAEDRES
C          PAR SUPPRESSION D'UNE DIAGONALE AUTOUR DE
C          LAQUELLE S'ENROULE 4, 5, ... TETRAEDRES

C ENTREES:
C --------
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NTEQM  : NO NOTETR DU TETRAEDRE DE QUALITE MEDIOCRE A TRAITER
C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C N1TETS : N1TETS(NS) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET NS
C NUDTETR: NO NOTETR DU DERNIER TETRAEDRE ACTIF
C IVOLTE : =0 TABLEAU NVOLTE NON PRESENT
C          =1 TABLEAU NVOLTE     PRESENT
C MXTE1A : NOMBRE MAXIMAL DE TETRAEDRES DECLARABLES DANS NOTE1A

C SORTIES:
C --------
C NVOLTE : NUMERO DE VOLUME (1 a NBVOPA) DE CHAQUE TETRAEDRE
C          -1 SI VOLUME INCONNU (Exemple: les TETRAEDRES VIDES DE NOTETR)
C NBTE1A : NOMBRE DE TETRAEDRES MODIFIES DANS NOTE1A
C          =0 SI PAS DE MODIFICATION DES TETRAEDRES
C NOTE1A : NO NOTETR DES NBTE1A TETRAEDRES MODIFIES
C MODIFT : =0 PAS DE MODIFICATION DES TETRAEDRES
C          =1        MODIFICATION DES TETRAEDRES. DECOMPOSITION FAITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET                St PIERRE du PERRAY   Mai 2018
C2345X7..............................................................012
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
      LOGICAL                          TRACTE0
      DOUBLE PRECISION  PTXYZD(4,*)
      DOUBLE PRECISION  ARMIN, ARMAX, SURFTR(4), VTEQM, V, VOLTET,
     %                  VOLET1, VOLET2, VOLTID, TID,
     %                  V1253, V1254, V3451, V3452,
     %                  V34,   V16,   V52,
     %                  V1632, V1624, V1645, V1635,
     %                  V5231, V5214, V5263, V5246
      REAL              QDMIN(2)
      INTEGER           NOTETR(8,*), N1TETS(*), NOTE1A(MXTE1A),
     %                  NVOLTE(*), NOSOTR(3), NOTEQM(8)
      CHARACTER*100     KTITRE
      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

      TRACTE0 = TRACTE
      MODIFT  = 0
      NBTE1A  = 0
      IERR    = 0


ccc      if( nteqm .eq. 3171 ) then
ccc         tracte = .true.
ccc         nt = nteqm
ccc         print*,'teqmtyqu: nteqm=',nteqm,
ccc     %          'notetr(',nt,')=',(notetr(kk,nt),kk=1,8)
ccc         print*
ccc      endif

      IF( NTEQM .LE. 0 ) GOTO 9990
      IF( NOTETR(1,NTEQM) .EQ. 0 ) GOTO 9990

C     VOLUME ET QUALITE DU TETRAEDRE NTEQM
      CALL QUATETD( PTXYZD( 1, NOTETR(1,NTEQM) ),
     %              PTXYZD( 1, NOTETR(2,NTEQM) ),
     %              PTXYZD( 1, NOTETR(3,NTEQM) ),
     %              PTXYZD( 1, NOTETR(4,NTEQM) ),
     %              ARMIN, ARMAX, SURFTR, VTEQM, QTEQM )

ccc      PRINT*,'teqmtyqu: TRAITER le TETRAEDRE',NTEQM,' de Volume=',VTEQM,
ccc     %       ' Qualite=',QTEQM,' St:',(NOTETR(M,NTEQM),M=1,8)


C     QUELLE EST LA NATURE DE LA QUALITE MEDIOCRE DU TETRAEDRE NTEQM?
C     EST IL PLAT DE TYPE QUADRANGLE ou TRIANGLE avec BARYCENTRE?
C     ---------------------------------------------------------------
      CALL TEQTPLAT( PTXYZD,NOTETR(1,NTEQM), NTYPQT,NST1,NST2,NST3,NST4)
C     NTYPQT: =1 TETRAEDRE QUASI-PLAT DE TYPE QUADRANGLE AVEC 2 DIAGONALES
C                DE SOMMETS NST1-NST2 et NST3-NST4 (de 1 a 4) DANS NTEQM
C             =2 TETRAEDRE QUASI-PLAT DE TYPE TRIANGLE AVEC UN SOMMET NST1
C                DE FACE OPPOSEE NST1+1 S'Y PROJETANT EN SON INTERIEUR
C                (NST2=NST3=NST4=0)

      IF( NTYPQT .NE. 1 ) GOTO 9990

C     TRACE DU TETRAEDRE NTEQM ET DES TETRAEDRES OPPOSES A SES 4 FACES
C     ----------------------------------------------------------------
ccc      TRACTE = .TRUE.
      KTITRE='teqmtyqu: TETRAEDRE          de QUALITE=           + 4 TET
     %RAEDRES OPPOSES A TRAITER'
      WRITE( KTITRE(21:28),'(I8)'  ) NTEQM
      WRITE( KTITRE(41:49),'(F9.6)') QTEQM
      CALL SANSDBL( KTITRE, L )
      CALL TRFETO9( KTITRE(1:L), PTXYZD, NTEQM, NOTETR )

C     ================================================================
C     NTEQM EST UN TETRAEDRE QUASI-PLAT DE TYPE QUADRANGLE AVEC 2
C     DIAGONALES DE SOMMETS NST1-NST2 et NST3-NST4 (de 1 a 4) DE NTEQM
C     ================================================================

C     SAUVEGARDE DU NUMERO PTXYZD DES 4 SOMMETS DU TETRAEDRE NTEQM
C     ET DU NUMERO NOTETR DE SES 4 TETRAEDRES OPPOSES
      DO L=1,8
         NOTEQM( L ) = NOTETR( L, NTEQM )
      ENDDO

C     RECHERCHE DE TOUS LES TETRAEDRES D'ARETE LA DIAGONALE NST1-NST2
      MXTEDI = MXTE1A / 2
      NS1 = NOTETR( NST1, NTEQM )
      NS2 = NOTETR( NST2, NTEQM )
      CALL TETR1A( NS1,    NS2,    N1TETS, NOTETR,
     %             NBTE12, MXTEDI, NOTE1A, IERR )
      IF( IERR .NE. 0 ) GOTO 9990

C     RECHERCHE DE TOUS LES TETRAEDRES D'ARETE LA DIAGONALE NST3-NST4
      NS3 = NOTETR( NST3, NTEQM )
      NS4 = NOTETR( NST4, NTEQM )
      CALL TETR1A( NS3,    NS4,    N1TETS,           NOTETR,
     %             NBTE34, MXTEDI, NOTE1A(MXTEDI+1), IERR )
      IF( IERR .NE. 0 ) GOTO 9990

      NBTE1A = MIN( NBTE12, NBTE34 )

      IF( NBTE1A .GE. 5 ) GOTO 9990
      GOTO( 9990, 9990, 3000, 4000 ),NBTE1A


C     ==========================================================
C     UNE DIAGONALE EST UNE ARETE COMMUNE A 3 TETRAEDRES DONT
C     LE TETRAEDRE DE QUALITE MEDIOCRE A TRAITER
C     RECHERCHE DU COUPLE DE FACES DU TETRAEDRE NTEQM TELLES QUE
C     LEUR TETRAEDRE OPPOSE AIENT UN SOMMET OPPOSE COMMUN NS5
C     ==========================================================
 3000 DO 50 NF1 = 1, 3

C        NUMERO DES 3 SOMMETS DE LA FACE NF1 DE NTEQM
         DO J=1,3
            NOSOTR( J ) = NOTETR( NOSOFATE(J,NF1), NTEQM )
         ENDDO

C        NO DU TETRAEDRE OPPOSE A LA FACE NF1
         NTEOP1 = NOTETR( 4+NF1, NTEQM )
         IF( NTEOP1 .LE. 0 ) GOTO 50

C        NS5 NO DU SOMMET DU TETRAEDRE NTEOP1 NON DANS LA FACE NF1
         DO 10 I1=1,4

            NS5 = NOTETR( I1, NTEOP1 )
            DO J=1,3
               IF( NS5 .EQ. NOSOTR(J) ) GOTO 10
            ENDDO
C           NS5 SOMMET DE NTEOP1 N'EST PAS UN SOMMET DE LA
C           FACE NF1 DE NTEQM
            GOTO 20

 10      ENDDO
         PRINT*,'teqmtyqu: ANOMALIE DE SOMMETS DANS NOTETR(',NTEQM,')=',
     %         (NOTETR(L,NTEQM),L=1,8)
         GOTO 9990

C        RECHERCHE D'UN SECOND TETRAEDRE OPPOSE A NTEQM
C        AYANT LE SOMMET NS5 DE NTEOP1
 20      DO NF2 = NF1+1, 4

C           NO DU TETRAEDRE OPPOSE A LA FACE NF2
            NTEOP2 = NOTETR( 4+NF2, NTEQM )
            IF( NTEOP2 .GT. 0 ) THEN
               DO I2=1,4
                  NS6 = NOTETR( I2, NTEOP2 )
                  IF( NS6 .EQ. NS5 ) GOTO 40
               ENDDO
            ENDIF

         ENDDO
         GOTO 50

C        NTEOP1 NTEOP2 ONT UN SOMMET COMMUN NS5 N'APPARTENANT PAS A NTEQM
C        RECHERCHE D'UN EVENTUEL 3-EME TETRAEDRE OPPOSE A NTEQM
C        AYANT LE SOMMET NS5 DE NTEOP1 ET NTEOP2
 40      DO NF3 = NF2+1, 4

C           NO DU TETRAEDRE OPPOSE A LA FACE NF3
            NTEOP3 = NOTETR( 4+NF3, NTEQM )
            IF( NTEOP3 .GT. 0 ) THEN
               DO I3=1,4
                  NS6 = NOTETR( I3, NTEOP3 )
                  IF( NS6 .EQ. NS5 ) THEN

C                    NTEOP1 NTEOP2 NTEOP3 ONT UN SOMMET COMMUN NS1
C                    N'APPARTENANT PAS A NTEQM
C                    FORMANT UN TETRAEDRE NTEQM LES ENGLOBANT ET
C                    LE SOMMET CENTRAL NS1 VA DISPARAITRE DE LA
C                    TETRAEDRISATION ACTUELLE
C                    ICI CE CAS N'EST PAS TRAITE
C                    ---------------------------------------------
                     GOTO 9990

                  ENDIF
               ENDDO
            ENDIF

         ENDDO

C        IL Y A 2 TETRAEDRES OPPOSES A 2 FACES AYANT UN MEME SOMMET NS5
C        IL N'Y A PAS 3 TETRAEDRES OPPOSES A 3 FACES DE NTEQM
C        AYANT UN MEME SOMMET NS5  C-A-D UN SOMMET DE NTEQM 
C        N'EST PAS INTERNE A UN GRAND TETRAEDRE FORME DE
C        3 TETRAEDRES OPPOSES
C        ----------------------------------------------------------------
C        NS5 LE SOMMET I1 DE NTEOP1 EST LE SOMMET I2 DE NTEOP2
C        NS1 NS2 LES 2 SOMMETS DE L'ARETE COMMUNE DES FACES NF1 NF2
C        NS3 NS4 LES 2 SOMMETS DE L'AUTRE DIAGONALE 
C        ESSAI D'ECHANGER NT12345 ( NTEQM NTOP1 NTOP2 ) PAR
C        NT3451 NT3452 DE FACE COMMUNE NS345
C        ----------------------------------------------------------------
         NBTE1A = 0
         IF( NTEOP1 .GT. 0 ) THEN
            NBTE1A = NBTE1A + 1
            NOTE1A( NBTE1A ) = NTEOP1
         ENDIF
         IF( NTEOP2 .GT. 0 ) THEN
            NBTE1A = NBTE1A + 1
            NOTE1A( NBTE1A ) = NTEOP2
         ENDIF
         NBTE1A = NBTE1A + 1
         NOTE1A( NBTE1A ) = NTEQM

cccC        AFFICHAGE DES NBTE1A TETRAEDRES
ccc         DO K=1,NBTE1A
ccc            NTE = NOTE1A( K )
ccc            CALL QUATETD( PTXYZD( 1, NOTETR(1,NTE) ),
ccc     %                    PTXYZD( 1, NOTETR(2,NTE) ),
ccc     %                    PTXYZD( 1, NOTETR(3,NTE) ),
ccc     %                    PTXYZD( 1, NOTETR(4,NTE) ),
ccc     %                    ARMIN, ARMAX, SURFTR, V, Q )
ccc         PRINT *,'teqmtyqu 100:',K,' old tetra',NTE,' St:',
ccc     %          (NOTETR(I,NTE),I=1,8),' Q=',Q,' V=',V
ccc         ENDDO

C        LE TETRAEDRE NTEQM QUASI-PLAT EST DE TYPE QUADRANGLE
C        CALCUL DES 2 SOMMETS NS1 NS2 DE L'ARETE COMMUNE AUX 2 FACES NF1 NF2
C               DES 2 AUTRES SOMMETS NS3 NS4 DU TETRAEDRE
C               DES 2 AUTRES FACES   NF3 NF4 DU TETRAEDRE
         CALL NS4C2FATE( NF1, NF2, NST1, NST2, NST3, NST4, NF3, NF4 )

         NS1 = NOTEQM( NST1 )
         NS2 = NOTEQM( NST2 )
         NS3 = NOTEQM( NST3 )
         NS4 = NOTEQM( NST4 )

C        COMPARAISON DU VOLUME DES 2 CAS
         V1253 = VOLTET( PTXYZD( 1, NOTETR(1,NTEOP1) ),
     %                   PTXYZD( 1, NOTETR(2,NTEOP1) ),
     %                   PTXYZD( 1, NOTETR(3,NTEOP1) ),
     %                   PTXYZD( 1, NOTETR(4,NTEOP1) ) )

         V1254 = VOLTET( PTXYZD( 1, NOTETR(1,NTEOP2) ),
     %                   PTXYZD( 1, NOTETR(2,NTEOP2) ),
     %                   PTXYZD( 1, NOTETR(3,NTEOP2) ),
     %                   PTXYZD( 1, NOTETR(4,NTEOP2) ) )

         V3451 = VOLTET( PTXYZD( 1, NS3 ),
     %                   PTXYZD( 1, NS5 ),
     %                   PTXYZD( 1, NS4 ),
     %                   PTXYZD( 1, NS1 ) )

         IF( V3451 .LT. 0D0 ) THEN
C           PERMUTATION POUR OBTENIR V3451>0
            N    = NST3
            NST3 = NST4
            NST4 = N

            N   = NS3
            NS3 = NS4
            NS4 = N
            V3451 = - V3451
         ENDIF

         V3452 = VOLTET( PTXYZD( 1, NS3 ),
     %                   PTXYZD( 1, NS4 ),
     %                   PTXYZD( 1, NS5 ),
     %                   PTXYZD( 1, NS2 ) )

         VOLET1 = VTEQM + V1253 + V1254
         VOLET2 = V3451 + ABS( V3452 )

C        TAILLE SOUHAITEE MOYENNE AUTOUR DES 5 SOMMETS
         TID = ( PTXYZD(4,NS1) + PTXYZD(4,NS2) + PTXYZD(4,NS3)
     %         + PTXYZD(4,NS4) + PTXYZD(4,NS5) ) / 5D0

C         VOLUME du TETRAEDRE IDEAL PAR RAPPORT A LA TAILLE
C         SOUHAITEE DES ARETES
          VOLTID = ( TID * TID * TID ) / 6D0

ccc      IF( ABS(VOLET1-VOLET2) .GT. VOLMOY*0.02D0 ) THEN
ccc      IF( ABS(VOLET1-VOLET2) .GT. VOLTID*1D-3 ) THEN
         IF( ABS(VOLET1-VOLET2) .GT. VOLTID*1D-4 ) THEN
C           VOLUMES TROP DIFFERENTS
            GOTO 50
         ENDIF

C        VISUALISATION DES TETRAEDRES NTEOP1 NTEOP2 NTEQM
         KTITRE='teqmtyqu: 3 TETRAEDRES VONT ETRE REMPLACES PAR 2 TETRAE
     %DRES'
         CALL SANSDBL( KTITRE, L )
         CALL TRAFNBTE( KTITRE(1:L), PTXYZD, NBTE1A, NOTE1A, NOTETR )

         IF( IVOLTE .NE. 0 ) THEN
            NOVOLU = NVOLTE( NTEOP1 )
            IF( NVOLTE( NTEOP2 ) .NE. NOVOLU ) THEN
C              VOLUMES DIFFERENTS
               GOTO 50
            ENDIF
         ENDIF

C        CREATION DES 2 TETRAEDRES 3451 et 3452
         IF( N1TEVI .LE. 0 ) GOTO 9050
         NT3451  = N1TEVI
         NUDTETR = MAX( NUDTETR, N1TEVI )
         IF( IVOLTE .NE. 0 ) NVOLTE( N1TEVI ) = NOVOLU
         N1TEVI  = NOTETR( 5, N1TEVI )

         IF( N1TEVI .LE. 0 ) GOTO 9050
         NT3452  = N1TEVI
         NUDTETR = MAX( NUDTETR, N1TEVI )
         IF( IVOLTE .NE. 0 ) NVOLTE( N1TEVI ) = NOVOLU
         N1TEVI  = NOTETR( 5, N1TEVI )

C        LES 3 TETRAEDRES DANS LES 2 TETRAEDRES 3451 3452
C        NT3451 DEVIENT LE TETRAEDRE DE SOMMETS 3541
         NOTETR( 1, NT3451 ) = NS3
         NOTETR( 2, NT3451 ) = NS5
         NOTETR( 3, NT3451 ) = NS4
         NOTETR( 4, NT3451 ) = NS1

C        NT3452 DEVIENT LE TETRAEDRE DE SOMMETS 3452
         NOTETR( 1, NT3452 ) = NS3
         NOTETR( 2, NT3452 ) = NS4
         NOTETR( 3, NT3452 ) = NS5
         NOTETR( 4, NT3452 ) = NS2

C        TETRAEDRES OPPOSES A NT3451 et NT3452 PAR LA FACE 345 
         NOTETR( 5, NT3451 ) = NT3452
         NOTETR( 5, NT3452 ) = NT3451

C        STOCKAGE POUR UTILISER UNE BOUCLE
         NOTE1A( 4 ) = NT3451
         NOTE1A( 5 ) = NT3452

         DO M=4,5

C           NT3451 ou NT3452
            NT = NOTE1A( M )

C           FACE 2 A 4 DE NT
            DO K=2,4

C              FACE K DE NT
               NOSOTR( 1 ) = NOTETR( NOSOFATE(1,K), NT )
               NOSOTR( 2 ) = NOTETR( NOSOFATE(2,K), NT )
               NOSOTR( 3 ) = NOTETR( NOSOFATE(3,K), NT )

               NTEOPP = -1
               DO I=1,NBTE1A
                  NTE = NOTE1A( I )
                  CALL NUFATRTE( NOSOTR, NOTETR(1,NTE), NF )
                  IF( NF .GT. 0 ) THEN
                     NTEOPP = NOTETR( 4+NF, NTE )
                     IF( NTEOPP .GT. 0 ) THEN
                        CALL NUFATRTE( NOSOTR, NOTETR(1,NTEOPP), NFF )
                        IF( NFF .GT. 0 ) THEN
                           NOTETR( 4+NFF, NTEOPP ) = NT
                        ELSE
C                         NE DOIT PAS ARRIVER....
                          PRINT*,'teqmtyqu: ANOMALIE NOSOTR=',
     %                           NOSOTR,' NOTETR(',NTEOPP,')=',
     %                          (NOTETR(L,NTEOPP),L=1,8)
                        ENDIF
                        GOTO 45
                     ENDIF
                  ENDIF
               ENDDO

C              LE TETRAEDRE OPPOSE A LA FACE K DE NT
 45            NOTETR( 4+K, NT ) = NTEOPP

            ENDDO

         ENDDO

C        MISE A JOUR DE N1TETS
         N1TETS( NS1 ) = NT3451
         N1TETS( NS2 ) = NT3452
         N1TETS( NS3 ) = NT3451
         N1TETS( NS4 ) = NT3451
         N1TETS( NS5 ) = NT3451

C        SUPPRESSION DES NBTE1A TETRAEDRES DU TABLEAU NOTETR
         DO K = 1, NBTE1A
C           DESTRUCTION DE NTE DU TABLEAU NOTETR
            NTE = NOTE1A( K )
            NOTETR(1,NTE) = 0
            NOTETR(5,NTE) = N1TEVI
            N1TEVI = NTE
            IF( IVOLTE .NE. 0 ) NVOLTE( NTE ) = -1
         ENDDO

C        POUR LA VISUALISATION DES 2 TETRAEDRES CREES
         NBTE1A1 = 2
         NBTE1A  = 2
         NOTE1A( 1 ) = NT3451
         NOTE1A( 2 ) = NT3452

         GOTO 9000

C        FIN DE BOUCLE SUR NF1 LES FACES DE NTEQM
 50   ENDDO

C     FIN 3000 DU CAS 1 DIAGONALE ARETE DE 3 TETRAEDRES
C     =================================================


C     =======================================================================
C     UNE DIAGONALE EST UNE ARETE COMMUNE A 4 TETRAEDRES DONT
C     LE TETRAEDRE DE QUALITE MEDIOCRE A TRAITER
C     ESSAI DE TROUVER UNE DIAGONALE AUTOUR DE LAQUELLE S'ENROULENT >=4
C     TETRAEDRES ET LA FAIRE DISPARAITREE EN FAISANT 4t->4t ou PLUS COMPLIQUE
C     =======================================================================
 4000 DO 100 NOPASS = 1, 2

         QDMIN( NOPASS ) = -2.0

         IF( NOPASS .EQ. 2 ) THEN
C           ESSAI DE FAIRE DISPARAITRE LA DIAGONALE NS1-NS2
C           PERMUTATION DES 2 DIAGONALES POUR RECOMMENCER
C           AVEC LA PREMIERE DIAGONALE A FAIRE DISPARAITRE
            K    = NST3
            NST3 = NST1
            NST1 = K

            K    = NST4
            NST4 = NST2
            NST2 = K
ccc      ELSE
cccC        ESSAI DE FAIRE DISPARAITRE LA DIAGONALE NS3-NS4
         ENDIF

C        LES SOMMETS DE LA DIAGONALE NST1-NST2 A CONSERVER
         NS1 = NOTETR( NST1, NTEQM )
         NS2 = NOTETR( NST2, NTEQM )

C        LES SOMMETS DE L'AUTRE DIAGONALE NST3-NST4 A DISPARAITRE
         NS3 = NOTETR( NST3, NTEQM )
         NS4 = NOTETR( NST4, NTEQM )

C        RECHERCHE DES 2 FACES NF1 NF2 D'ARETE LA DIAGONALE NST3-NST4
         CALL NFA2STE( NST3, NST4, NF1, NF2 )

C        NUMERO NOTETR DES TETRAEDRES OPPOSES A NF1 et NF2
         NTEOP1 = NOTETR( 4+NF1, NTEQM )
         NTEOP2 = NOTETR( 4+NF2, NTEQM )
         IF( NTEOP1 .LE. 0 .OR. NTEOP2 .LE. 0  ) GOTO 100

C        RECHERCHE DU SOMMET NS5 DE NTEOP1 NON SOMMET DE NTEQM
         DO K=1,4
            NS5 = NOTETR( K, NTEOP1 )
            IF( NS5 .NE. NS1 .AND. NS5 .NE. NS2 .AND.
     %          NS5 .NE. NS3 .AND. NS5 .NE. NS4 ) GOTO 4
         ENDDO

C        RECHERCHE DU SOMMET NS6 DE NTEOP2 NON SOMMET DE NTEQM
 4       DO K=1,4
            NS6 = NOTETR( K, NTEOP2 )
            IF( NS6 .NE. NS1 .AND. NS6 .NE. NS2 .AND.
     %          NS6 .NE. NS3 .AND. NS6 .NE. NS4 ) GOTO 8
         ENDDO

C        POUR DECIDER DU SENS DE LA DIAGONALE NS3-NS4
 8       CALL QUATETD( PTXYZD(1,NOTETR(1,NTEOP1)),
     %                 PTXYZD(1,NOTETR(2,NTEOP1)),
     %                 PTXYZD(1,NOTETR(3,NTEOP1)),
     %                 PTXYZD(1,NOTETR(4,NTEOP1)),
     %                 ARMIN, ARMAX, SURFTR, V, Q )

         CALL QUATETD( PTXYZD(1,NS1),
     %                 PTXYZD(1,NS3),
     %                 PTXYZD(1,NS4),
     %                 PTXYZD(1,NS5),
     %                 ARMIN, ARMAX, SURFTR, V34, Q )

      IF( V * V34 .LT. 0D0 ) THEN
C        PERMUTATION NS3 NS4 et NS5 NS6
         K    = NST3
         NST3 = NST4
         NST4 = K

         K   = NS3
         NS3 = NS4
         NS4 = K

         K   = NS5
         NS5 = NS6
         NS6 = K

      ENDIF

C     RECHERCHE DE TOUS LES TETRAEDRES D'ARETE DIAGONALE NST3-NST4
      CALL TETR1A( NS3,    NS4,    N1TETS, NOTETR,
     %             NBTE1A, MXTE1A, NOTE1A, IERR )

      IF( NBTE1A .NE. 4 .OR. IERR .NE. 0 ) THEN
         GOTO 100
      ENDIF

C     CALCUL DE LA QUALITE DES NBTE1A TETRAEDRES D'ARETE NS3-NS4
C     ----------------------------------------------------------
      V34    = 0D0
      Q34MIN = 2.
      DO K=1,NBTE1A
         NTEK = NOTE1A( K )
         CALL QUATETD( PTXYZD(1,NOTETR(1,NTEK)),
     %                 PTXYZD(1,NOTETR(2,NTEK)),
     %                 PTXYZD(1,NOTETR(3,NTEK)),
     %                 PTXYZD(1,NOTETR(4,NTEK)),
     %                 ARMIN, ARMAX, SURFTR, V, Q )
ccc         PRINT*,'teqmtyqu: TETRAEDRE',NTEK,' INITIAL   Volume=',V,
ccc     %       ' Qualite=',Q,' St:',(NOTETR(M,NTEK),M=1,8)
         V34 = V34 + V
         Q34MIN = MIN( Q34MIN, Q )
      ENDDO

C     CALCUL DE LA QUALITE DES 4 TETRAEDRES D'ARETE NS1-NS6
C     -----------------------------------------------------
      CALL QUATETD( PTXYZD(1,NS1), PTXYZD(1,NS6),
     %              PTXYZD(1,NS3), PTXYZD(1,NS2),
     %              ARMIN, ARMAX, SURFTR, V1632, Q1632 )

      CALL QUATETD( PTXYZD(1,NS1), PTXYZD(1,NS6),
     %              PTXYZD(1,NS2), PTXYZD(1,NS4),
     %              ARMIN, ARMAX, SURFTR, V1624, Q1624 )

      CALL QUATETD( PTXYZD(1,NS1), PTXYZD(1,NS6),
     %              PTXYZD(1,NS4), PTXYZD(1,NS5),
     %              ARMIN, ARMAX, SURFTR, V1645, Q1645 )

      CALL QUATETD( PTXYZD(1,NS1), PTXYZD(1,NS6),
     %              PTXYZD(1,NS3), PTXYZD(1,NS5),
     %              ARMIN, ARMAX, SURFTR, V1635, Q1635 )

      V16    = ABS(V1632) + ABS(V1624) + ABS(V1645) + ABS(V1635)
      Q16MIN = MIN( Q1632, Q1624, Q1645, Q1635 )


C     CALCUL DE LA QUALITE DES 4 TETRAEDRES D'ARETE NS5-NS2
C     -----------------------------------------------------
      CALL QUATETD( PTXYZD(1,NS5), PTXYZD(1,NS2),
     %              PTXYZD(1,NS3), PTXYZD(1,NS1),
     %              ARMIN, ARMAX, SURFTR, V5231, Q5231 )

      CALL QUATETD( PTXYZD(1,NS5), PTXYZD(1,NS2),
     %              PTXYZD(1,NS1), PTXYZD(1,NS4),
     %              ARMIN, ARMAX, SURFTR, V5214, Q5214 )

      CALL QUATETD( PTXYZD(1,NS5), PTXYZD(1,NS2),
     %              PTXYZD(1,NS6), PTXYZD(1,NS3),
     %              ARMIN, ARMAX, SURFTR, V5263, Q5263 )

      CALL QUATETD( PTXYZD(1,NS5), PTXYZD(1,NS2),
     %              PTXYZD(1,NS4), PTXYZD(1,NS6),
     %              ARMIN, ARMAX, SURFTR, V5246, Q5246 )

      V52    = ABS(V5231) + ABS(V5214) + ABS(V5263) + ABS(V5246)
      Q52MIN = MIN( Q5231, Q5214, Q5263, Q5246 )


ccc      if( nteqm .eq. 3171 ) then
ccc      PRINT*,'teqmtyqu: NST1=',NST1,' NST2=',NST2,
ccc     %                ' NST3=',NST3,' NST4=',NST4
ccc      PRINT*,'teqmtyqu: NS1=',NS1,' NS2=',NS2,
ccc     %                ' NS3=',NS3,' NS4=',NS4,
ccc     %                ' NS5=',NS5,' NS6=',NS6
ccc      PRINT*,'teqmtyqu: V34=',V34,' Q34MIN=',Q34MIN
ccc      PRINT*,'teqmtyqu: V16=',V16,' Q16MIN=',Q16MIN
ccc      PRINT*,'teqmtyqu: V52=',V52,' Q52MIN=',Q52MIN
ccc      print*
ccc      endif

      QDMIN( NOPASS ) = MIN( Q34MIN, Q16MIN, Q52MIN )

      IF( NOPASS .EQ. 2 ) THEN
         IF( QDMIN(2) .LE. QDMIN(1) ) THEN
C           C'EST POUR EVITER D'ECRASER LES VALEURS OPTIMALES
            QDMIN( 2 ) = -2.
            GOTO 100
         ENDIF
      ENDIF

      IF( Q34MIN .GE. Q16MIN .AND. Q34MIN .GE. Q52MIN ) THEN
C        PAS D'AMELIORATION => PASSAGE A LA SECONDE DIAGONALE
         QDMIN( NOPASS ) = -2.
         GOTO 100
      ENDIF

C     MODIFICATION PERMISE?
      IF( Q16MIN .LE. 0 .AND. Q52MIN .LE. 0 ) THEN
C        NON: TROP DE TETRAEDRES PLATS...
         PRINT*,'teqmtyqu: ABANDON avec Q34MIN=',Q34MIN,
     %          ' Q16MIN=',Q16MIN,' Q52MIN=',Q52MIN
         GOTO 9990
      ENDIF

      IF( Q16MIN .GE. Q34MIN .AND. Q16MIN .GE. Q52MIN ) THEN
C        CHOIX DES TETRAEDRES D'ARETE 16
         IF( ABS(V34-V16) .GT. V16*1D-4 ) THEN
C           NON: VOLUMES DIFFERENTS
            PRINT*,'teqmtyqu: ABANDON avec V34=',V34,' V16=',V16,
     %             ' Q34MIN=',Q34MIN,' Q16MIN=',Q16MIN
            GOTO 9990
         ENDIF
C        CHOIX DES TETRAEDRES D'ARETE 16
         GOTO 16
      ELSE
         IF( ABS(V34-V52) .GT. V52*1D-4 ) THEN
            PRINT*,'teqmtyqu: ABANDON avec V34=',V34,' V52=',V52,
     %             ' Q34MIN=',Q34MIN,' Q52MIN=',Q52MIN
            GOTO 9990
         ENDIF
C        CHOIX DES TETRAEDRES D'ARETE 52
         GOTO 52
      ENDIF

C     REMPLACEMENT DES 4 TETRAEDRES D'ARETE 34
C              PAR LES 4 TETRAEDRES D'ARETE 16
C     ----------------------------------------
 16   NTE1 = NOTE1A( 1 )
      NOTETR(1,NTE1) = NS1
      NOTETR(2,NTE1) = NS6
      NOTETR(3,NTE1) = NS3
      NOTETR(4,NTE1) = NS2

      NTE2 = NOTE1A( 2 )
      NOTETR(1,NTE2) = NS1
      NOTETR(2,NTE2) = NS6
      NOTETR(3,NTE2) = NS2
      NOTETR(4,NTE2) = NS4

      NTE3 = NOTE1A( 3 )
      NOTETR(1,NTE3) = NS1
      NOTETR(2,NTE3) = NS6
      NOTETR(3,NTE3) = NS4
      NOTETR(4,NTE3) = NS5

      NTE4 = NOTE1A( 4 )
      NOTETR(1,NTE4) = NS1
      NOTETR(2,NTE4) = NS6
      NOTETR(3,NTE4) = NS3
      NOTETR(4,NTE4) = NS5

      N1TETS( NS1 ) = NTE1
      N1TETS( NS2 ) = NTE1
      N1TETS( NS3 ) = NTE1
      N1TETS( NS4 ) = NTE2
      N1TETS( NS5 ) = NTE3
      N1TETS( NS6 ) = NTE1

      GOTO 105

C     REMPLACEMENT DES 4 TETRAEDRES D'ARETE 34
C              PAR LES 4 TETRAEDRES D'ARETE 52
C     ----------------------------------------
 52   NTE1 = NOTE1A( 1 )
      NOTETR(1,NTE1) = NS5
      NOTETR(2,NTE1) = NS2
      NOTETR(3,NTE1) = NS3
      NOTETR(4,NTE1) = NS1

      NTE2 = NOTE1A( 2 )
      NOTETR(1,NTE2) = NS5
      NOTETR(2,NTE2) = NS2
      NOTETR(3,NTE2) = NS1
      NOTETR(4,NTE2) = NS4

      NTE3 = NOTE1A( 3 )
      NOTETR(1,NTE3) = NS5
      NOTETR(2,NTE3) = NS2
      NOTETR(3,NTE3) = NS6
      NOTETR(4,NTE3) = NS3

      NTE4 = NOTE1A( 4 )
      NOTETR(1,NTE4) = NS5
      NOTETR(2,NTE4) = NS2
      NOTETR(3,NTE4) = NS4
      NOTETR(4,NTE4) = NS6

      N1TETS( NS1 ) = NTE1
      N1TETS( NS2 ) = NTE1
      N1TETS( NS3 ) = NTE1
      N1TETS( NS4 ) = NTE2
      N1TETS( NS5 ) = NTE1
      N1TETS( NS6 ) = NTE3

      GOTO 105

 100  ENDDO

      IF( QDMIN(1) .EQ. -2. .AND. QDMIN(2) .EQ. -2. ) GOTO 9990


C     AJOUT DANS NOTE1A DES TETRAEDRES OPPOSES AUX NBTE1A TETRAEDRES
C     --------------------------------------------------------------
 105  NBTE1A1 = NBTE1A
      DO 110 K = 1, NBTE1A

C        LE K-EME TETRAEDRE D'ARETE NS1-NS2
         NTEK = NOTE1A( K )

         IF( NTEK .GT. 0 .AND. NOTETR(1,NTEK) .GT. 0 ) THEN

C           VOLUME ET QUALITE DU TETRAEDRE NTEK
            CALL QUATETD( PTXYZD( 1, NOTETR(1,NTEK) ),
     %                    PTXYZD( 1, NOTETR(2,NTEK) ),
     %                    PTXYZD( 1, NOTETR(3,NTEK) ),
     %                    PTXYZD( 1, NOTETR(4,NTEK) ),
     %                    ARMIN, ARMAX, SURFTR, V, Q )

ccc         if( nteqm .eq. 3171 ) then
ccc            PRINT*,'teqmtyqu: TETRAEDRE',NTEK,'   FINAL   Volume=',V,
ccc     %             ' Qualite=',Q,' St:',(NOTETR(M,NTEK),M=1,8)
ccc         endif

C           STOCKAGE NOTE1A DES TETRAEDRES OPPOSES AUX 4 FACES DE NTEK
C           POUR LA MISE A JOUR APRES LES MODIFICATIONS
            DO 108 NF=1,4

C              NO TETRAEDRE OPPOSE A LA FACE NF DE NTEK
               CALL NOFAOP( NF, NTEK, NOTETR, NFOP, NTEOP )

               IF( NTEOP .GT. 0 .AND. NOTETR(1,NTEOP) .GT. 0 ) THEN

C                 NTEK A ETE MODIFIE
C                 LE TETRAEDRE OPPOSE A LA FACE NFOP DE NTEOP DEVIENT INCONNU
                  NOTETR( 4+NFOP, NTEOP ) = -1
C                 LE TETRAEDRE OPPOSE A LA FACE NF   DE NTEK DEVIENT INCONNU
                  NOTETR( 4+NF, NTEK ) = -1

C                 NTEOP EST AJOUTE A NOTE1A SI CE N'EST DEJA FAIT
                  DO L=1,NBTE1A1
                     IF( NOTE1A(L) .EQ. NTEOP ) GOTO 108
                  ENDDO
                  IF( NBTE1A1 .GE. MXTE1A ) THEN
                     PRINT*,'teqmtyqu: AUGMENTER MXTE1A=',MXTE1A
                     GOTO 9990
                  ENDIF
ccc                  PRINT*,'teqmtyqu: tetraedre(',NTEOP,'):',
ccc     %                   (NOTETR(L,NTEOP),L=1,8),
ccc     %                   ' OPPOSE est AJOUTE a NOTE1A'
                  NBTE1A1 = NBTE1A1 + 1
                  NOTE1A( NBTE1A1 ) = NTEOP
               ENDIF

 108        ENDDO

         ENDIF

 110  ENDDO


ccc      if( nteqm .eq. 3171 ) then
ccc         do k=1, NBTE1A1
ccc            nt = NOTE1A( K )
ccc            print*,'Teqmtyqu: NTEQM=',NTEQM,
ccc     %          'NOTETR(',NT,')=',(NOTETR(kk,NT),kk=1,8)
ccc         enddo
ccc         tracte = .true.
ccc      KTITRE='teqmtyqu:           TETRAEDRES AVEC TETRAEDRES OPPOSES DOU
ccc     %BLES'
ccc         WRITE(KTITRE(9:15),'(I7)') NBTETRA
ccc         CALL SANSDBL( KTITRE, NBC )
ccc         CALL TRAFNBTE( KTITRE(1:NBC), PTXYZD, NBTE1A1, NOTE1A, NOTETR )
ccc         print*
ccc      endif

C     COMPLETION DES CHAINAGES DES TETRAEDRES OPPOSES PAR LES FACES
C     -------------------------------------------------------------
 9000 CALL MJOPTE( NBTE1A1, NOTE1A, N1TETS, NOTETR, NUDTETR,
     %             N1TEVI,  PTXYZD, NBFANR )

ccc      if( nteqm .eq. 3171 ) then
ccc         nt = nteqm
ccc         print*,'TEQMtyqU: nteqm=',nteqm,
ccc     %          'notetr(',nt,')=',(notetr(kk,nt),kk=1,8)

ccc         do k=1, NBTE1A1
ccc            nt = NOTE1A( K )
ccc            print*,'TEQMtyqu: NTEQM=',NTEQM,
ccc     %             'NOTETR(',NT,')=',(NOTETR(kk,NT),kk=1,8)
ccc         enddo
ccc         print*
ccc      endif


      IF( NBFANR .GT. 0 ) THEN
C        AU MOINS UNE FACE N'A PAS DE TETRAEDRE OPPOSE
         PRINT*,'teqmtyqu: NTEQM=',NTEQM,' Probleme ',NBFANR,
     %          ' FACES DE TETRAEDRE OPPOSE NON RETROUVE'
      ENDIF

      KTITRE='teqmtyqu: 4 TETRAEDRES SANS la DIAGONALE                 '
      WRITE( KTITRE(42:49),'(I8)') NS3
      WRITE( KTITRE(51:58),'(I8)') NS4
      CALL SANSDBL( KTITRE, L )
      CALL SANSDBL( KTITRE, L )
      CALL TRAFNBTE( KTITRE(1:L), PTXYZD, NBTE1A, NOTE1A, NOTETR )

C     FIN 4000 DU CAS 1 DIAGONALE ARETE DE 4 TETRAEDRES
C     =================================================



C     AU MOINS UN TETRAEDRE NOTE1A A ETE CREE ou MODIFIE
C     ==================================================
C     AFFICHAGE DES NBTE1A TETRAEDRES CREES
      QminGR = 2.
      DO K = 1, NBTE1A
         NTE = NOTE1A( K )
C        VOLUME ET QUALITE DU TETRAEDRE NTE
         CALL QUATETD( PTXYZD( 1, NOTETR(1,NTE) ),
     %                 PTXYZD( 1, NOTETR(2,NTE) ),
     %                 PTXYZD( 1, NOTETR(3,NTE) ),
     %                 PTXYZD( 1, NOTETR(4,NTE) ),
     %                 ARMIN, ARMAX, SURFTR, V, Q )
         QminGR = MIN( QminGR, Q )
      ENDDO

C     VISUALISATION FINALE DES NBTE1A TETRAEDRES MODIFIES
ccc      TRACTE = .TRUE.
      KTITRE='teqmtyqu:                 NOUVEAUX TETRAEDRES DE  QUALITE 
     %MIN            '
      WRITE(KTITRE(11:17),'(I7)'  ) NBTE1A
      WRITE(KTITRE(63:71),'(F9.6)') QminGR
      CALL SANSDBL( KTITRE, L )
      CALL TRAFNBTE( KTITRE(1:L), PTXYZD, NBTE1A,  NOTE1A, NOTETR )
ccc      CALL TRAFNBTE( KTITRE(1:L), PTXYZD, NBTE1A1, NOTE1A, NOTETR )

C     DECOMPOSITION REUSSIE DES TETRAEDRES D'AXE NS1-NS2
      MODIFT = 1
      GOTO 9999


C     TAILLE INSUFFISANTE DU TABLEAU NOTETR
 9050 PRINT *,'teqmtyqu: SATURATION DU TABLEAU NOTETR: N1TEVI=',N1TEVI


C     SORTIE SANS MODIFICATION DU TETRAEDRE ET DE SES VOISINS
C     -------------------------------------------------------
 9990 NBTE1A = 0
      MODIFT = 0

C     BILAN FINAL
 9999 IF( MODIFT .NE. 0 ) THEN
         PRINT *,'teqmtyqu: TETRAEDRE',NTEQM,
     %           ' NOTETR:',(NOTEQM(kk),kk=1,4),' V=',VTEQM,' Q=',QTEQM,
     %           ' MODIFIE pour QminGR=',QminGR
      ENDIF

      TRACTE = TRACTE0
      RETURN
      END
