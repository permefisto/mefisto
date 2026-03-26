      SUBROUTINE TEQMTYTR( QTEMED, PTXYZD, NTEQM,  N1TEVI, NOTETR,
     %                     N1TETS, NUDTETR,
     %                     INFACO, MXFACO, LEFACO, IVOLTE, NVOLTE,
     %                     MXTE1A, NBTE1A, NOTE1A, MODIFT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     TRAITEMENT DU TETRAEDRE NTEQM DE QUALITE MEDIOCRE
C ----     SI DE TYPE TETRAEDRE ECRASE EN UN TRIANGLE ET PAR
C          DECOMPOSITION DE TETRAEDRES EN 2 TETRAEDRES AUTOUR DE LA
C          PLUS LONGUE ARETE DE NTEQM 'CONTENANT' PRESQUE 3 SOMMETS
C          (I.E. UNE DES 4 FACES EST DE SURFACE TRES FAIBLE) 
C          ou 2T->3T AVEC SON TETRAEDRE OPPOSE A SA PLUS GRANDE FACE

C ENTREES:
C --------
C QTEMED : QUALITE MEDIOCRE DES TETRAEDRES AU DESSOUS DE LAQUELLE
C          UNE AMELIORATION DES TETRAEDRES EST DEMANDEE
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
C INFACO : =1 PAS DE 2T->3T SI LA FACE COMMUNE AUX 2T APPARTIENT A LEFACO
C          =0 SI PAS DE CONTROLE SUR LA FACE COMMUNE ET SON APPARTENANCE
C             A LEFACO
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

C MXTE1A : NOMBRE MAXIMAL DE TETRAEDRES MODIFIABLES
C IVOLTE : =0 TABLEAU NVOLTE NON PRESENT
C          =1 TABLEAU NVOLTE     PRESENT

C SORTIES:
C --------
C NVOLTE : NUMERO DE VOLUME (1 a NBVOPA) DE CHAQUE TETRAEDRE
C          -1 SI VOLUME INCONNU (Exemple: les TETRAEDRES VIDES DE NOTETR)
C NUmDTETR: NO NOTETR DU DERNIER TETRAEDRE ACTIF
C NBTE1A : NOMBRE DE SOUS-TETRAEDRES CREES DANS NOTE1A
C          =0 NE SIGNIFIE PAS QU'AUCUN TETRAEDRE N'AIT ETE MODIFIE
C NOTE1A : NO NOTETR DES NBTE1A SOUS-TETRAEDRES CREES
C MODIFT : =0 PAS DE MODIFICATION DES TETRAEDRES NTEQM ...
C          =1        MODIFICATION DES TETRAEDRES. DECOMPOSITION FAITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET                St PIERRE du PERRAY   Mai 2018
C2345X7..............................................................012
      PARAMETER        ( MXTS1S2=1024 )
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
      LOGICAL                          TRACTE0
      DOUBLE PRECISION  PTXYZD(4,*)
      INTEGER           NOTETR(8,*), N1TETS(*), NOTE1A(MXTE1A),
     %                  LEFACO(1:11,0:MXFACO),
     %                  NOTS1S2(MXTS1S2), NVOLTE(*)

      CHARACTER*128     KTITRE
      INTEGER           NOSOTR(3), NTEN23(3), NOTEQM(8)
ccc      INTEGER           NOSOTR2(3)
      DOUBLE PRECISION  ARMIN, ARMAX, SURFTR(4)
      DOUBLE PRECISION  VOLTET, V, VT0, VT12, VTEQM, VOLTID,
     %                  V1, V2, DMAX, D, S, ST, SFMIN, SFMAX
      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

      MODIFT  = 0
      NBTS1S2 = 0
      NBTE1A  = 0
      QminGR  =-2.

ccc      IF( INFACO .NE. 0 ) RETURN

      TRACTE0 = TRACTE

      IF( NTEQM .LE. 0 ) GOTO 9990
      IF( NOTETR(1,NTEQM) .EQ. 0 ) GOTO 9990

      NOVOLU = 0
      IF( IVOLTE .NE. 0 ) THEN
         NOVOLU = NVOLTE( NTEQM )
      ELSE
         NOVOLU = 0
      ENDIF

C     QUELLE EST LA NATURE DE LA QUALITE MEDIOCRE DU TETRAEDRE NTEQM?
C     UNE DES 4 FACES DU TETRAEDRE EST ELLE DE SURFACE TRES FAIBLE?
C     I.E. LE TETRAEDRE NTEQM A T IL 3 SOMMETS SUR UNE MEME ARETE?
C     ---------------------------------------------------------------
C     RECHERCHE DE LA SURFACE>=0 DES 4 FACES DU TETRAEDRE NTEQM
      CALL QUATETD( PTXYZD( 1, NOTETR(1,NTEQM) ),
     %              PTXYZD( 1, NOTETR(2,NTEQM) ),
     %              PTXYZD( 1, NOTETR(3,NTEQM) ),
     %              PTXYZD( 1, NOTETR(4,NTEQM) ),
     %              ARMIN, ARMAX, SURFTR, VTEQM, QTEQM )

      IF( QTEQM .GE. QTEMED ) GOTO 9999

C     SAUVEGARDE DU TETRAEDRE NTEQM POUR AFFICHAGE FINAL DE SA MODIFICATION
      DO L = 1, 8
         NOTEQM( L ) = NOTETR( L, NTEQM )
      ENDDO

C     SOMME DES 4 SURFACES TRIANGULAIRES
      ST = SURFTR(1) + SURFTR(2) + SURFTR(3) + SURFTR(4)

C     RECHERCHE DE LA FACE DE SURFACE MINIMALE ET MAXIMALE
      NFMIN = 0
      SFMIN = 1D100
      NFMAX = 0
      SFMAX =-1D100
      DO NF=1,4
         S = SURFTR(NF)
         IF( S .LT. SFMIN ) THEN
            SFMIN = S
            NFMIN = NF
         ENDIF
         IF( S .GT. SFMAX ) THEN
            SFMAX = S
            NFMAX = NF
         ENDIF
      ENDDO

C     =================================================================
C     LE TETRAEDRE NTEQM AVEC LE TETRAEDRE OPPOSE PAR SA PLUS GRANDE
C     FACE NFMAX PEUT IL ETRE REMPLACE PAR 3 TETRAEDRES D'ARETE COMMUNE
C     LES 2 SOMMETS DES 2 TETRAEDRES N'APPARTENANT PAS A LA FACE MAX?
C     =================================================================
      CALL CH2T3T( INFACO, MXFACO, LEFACO, 1, IVOLTE, NVOLTE,
     %             NTEQM,  NFMAX,  NTOP,
     %             PTXYZD, N1TETS, NOTETR, N1TEVI, NUDTETR,
     %             NTEN23, IERR )
      IF( IERR .EQ. 0 ) THEN
C        DECOMPOSITION REUSSIE 2 TETRAEDRES NTEQM-NTOP -> 3T NTEN23
         MODIFT = 1
         QminGR = 2.
         DO K = 1, 3
C           LE K-EME TETRAEDRE CREE
            NTEK = NTEN23( K )
C           VOLUME ET QUALITE DU TETRAEDRE NTEK
            CALL QUATETD( PTXYZD(1,NOTETR(1,NTEK)),
     %                    PTXYZD(1,NOTETR(2,NTEK)),
     %                    PTXYZD(1,NOTETR(3,NTEK)),
     %                    PTXYZD(1,NOTETR(4,NTEK)),
     %                    ARMIN, ARMAX, SURFTR, V, Q )
            IF( Q .LT. QminGR ) THEN
               QminGR = Q
            ENDIF
         ENDDO
         GOTO 9999
      ENDIF


C     LA FACE NFMIN EST ELLE SUFFISAMMENT PETITE?
C     -------------------------------------------
ccc      IF( SFMIN .LE. (ST-SFMAX) * 0.05D0 ) THEN
      DO NF=1,4
         IF( NF .NE. NFMIN ) THEN
            IF( SFMIN .GE. SURFTR(NF) * 0.05D0 ) GOTO 9990
         ENDIF
      ENDDO

C        =====================================================
C        OUI: LE TETRAEDRE S'ECRASE EN UN TRIANGLE ET
C        LA FACE DE SURFACE MINIMALE EST ASSEZ FAIBLE POUR QUE
C        LE TETRAEDRE AIT UN SOMMET PROCHE D'UNE DE SES ARETES
C        =====================================================

C        LE NO PTXYZD DES 3 SOMMETS DE LA FACE NFMIN DE NTEQM
         DO L=1,3
            NOSOTR(L) = NOTETR( NOSOFATE(L,NFMIN), NTEQM )
         ENDDO

C        LONGUEUR DES 3 ARETES ET NO DE LA PLUS LONGUE DE LA FACE NFMIN
         LMAX = 0
         DMAX = 0D0
         DO L=1,3
            IF( L .EQ. 3 ) THEN
               LL = 1
            ELSE
               LL = L+1
            ENDIF
            D = SQRT(( PTXYZD(1,NOSOTR(LL)) - PTXYZD(1,NOSOTR(L)) )**2
     %              +( PTXYZD(2,NOSOTR(LL)) - PTXYZD(2,NOSOTR(L)) )**2
     %              +( PTXYZD(3,NOSOTR(LL)) - PTXYZD(3,NOSOTR(L)) )**2 )
            IF( D .GT. DMAX ) THEN
               DMAX = D
               LMAX = L
C              LES 2 SOMMETS DE L'ARETE LA PLUS LONGUE DE NFMIN
               NS1 = NOSOTR( L  )
               NS2 = NOSOTR( LL )
            ENDIF
         ENDDO

C        RECHERCHE DE TOUS LES TETRAEDRES D'ARETE NS1-NS2 LA PLUS LONGUE
         CALL TETR1A( NS1,     NS2,     N1TETS,  NOTETR,
     %                NBTS1S2, MXTS1S2, NOTS1S2, IERR )
         IF( NBTS1S2 .LE. 1 .OR. IERR .NE. 0 ) GOTO 9990


ccc         IF( NBTS1S2 .GT. 32 ) THEN
C           TROP DE TETRAEDRES AUTOUR DE L'ARETE NS1-NS2 => A TRACER
ccc            tracte = .true.
            KTITRE='teqmtytr: ARETE                    AVEC           TE
     %TRAEDRES'
            WRITE(KTITRE(17:24),'(I8)') NS1
            WRITE(KTITRE(26:33),'(I8)') NS2
            WRITE(KTITRE(41:45),'(I5)') NBTS1S2
            CALL SANSDBL( KTITRE, NBC )
            CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NS1,
     %                        NBTS1S2, NOTS1S2 )
ccc         ENDIF


C        LE 3-EME SOMMET DE LA FACE NFMIN DE SURFACE MINIMALE
         IF( LMAX .EQ. 1 ) THEN
            NS3 = NOSOTR( 3 )
         ELSE
            NS3 = NOSOTR( LMAX - 1 )
         ENDIF

C        VOLUME du TETRAEDRE IDEAL DE SOMMET NS3
         VOLTID = ( PTXYZD(4,NS3) **3 ) / 6

C        L'ARETE NS1-NS2 EST DECOUPEE EN 2 DEMI ARETES DE "MILIEU" NS3
C        CHAQUE TETRAEDRE AYANT CETTE ARETE NS1-NS2 EST DECOUPE
C        EN 2 TETRAEDRES A PARTIR DU SOMMET NS3 INTERNE A L'ARETE

C        SIMULATION POUR SAVOIR SI LA QUALITE SE DEGRADE OU S'AMELIORE
C        =============================================================
         NBTEK  = 0
         VT0    = 0D0
         VT12   = 0D0
         QMIN0  = 2.
         QminGR = 2.
         NS5 = 0
         DO 45 K = 1, NBTS1S2

C           LE K-EME TETRAEDRE D'ARETE NS1-NS2
            NTEK = NOTS1S2( K )

            IF( NTEK .LE. 0 ) GOTO 9990
C           TABLEAU DES TETRAEDRES D'ARETE NS1 NS2 DEFAILLANT
            IF( NOTETR(1,NTEK) .EQ. 0 ) GOTO 9990

C           VOLUME ET QUALITE DU TETRAEDRE NTEK
            CALL QUATETD( PTXYZD(1,NOTETR(1,NTEK)),
     %                    PTXYZD(1,NOTETR(2,NTEK)),
     %                    PTXYZD(1,NOTETR(3,NTEK)),
     %                    PTXYZD(1,NOTETR(4,NTEK)),
     %                    ARMIN, ARMAX, SURFTR, V, Q )
            VT0 = VT0 + V
            IF( Q .LT. QMIN0 ) THEN
               QMIN0 = Q
            ENDIF

C           NS4 NS5 LES 2 AUTRES SOMMETS DE NTEK
            DO L=1,4
               NS4 = NOTETR( L, NTEK )
               IF( NS4 .NE. NS1 .AND. NS4 .NE. NS2 ) GOTO 30
            ENDDO

 30         DO LL = L+1, 4
               NS5 = NOTETR( LL, NTEK )
               IF( NS5 .NE. NS1 .AND. NS5 .NE. NS2 ) GOTO 40
            ENDDO

C           NTEK NE DOIT PAS ETRE UN TETRAEDRE DE FACE NS1 NS2 NS3
 40         IF( NS4 .EQ. NS3 .OR. NS5 .EQ. NS3 ) GOTO 45 

C           UN TETRAEDRE NTEK DE PLUS A DECOUPER EN 2
            NBTEK = NBTEK + 1

C           POUR DETERMINER L'ORIENTATION DE NS4-NS5
C           VOLUME DU TETRAEDRE NS1 NS2 NS4 NS5 = NTEK
            V2 = VOLTET( PTXYZD(1,NS1), PTXYZD(1,NS2), 
     %                   PTXYZD(1,NS4), PTXYZD(1,NS5) )
            IF( V2 * V .LT. 0D0 ) THEN
               L   = NS4
               NS4 = NS5
               NS5 = L
            ENDIF
C           L'ORIENTATION DE NS1 NS2 NS4 NS5 EST CELLE DE NTEK

C           QUELLE QUALITE AURAIENT LES TETRAEDRES SUBDIVISES?
            CALL QUATETD( PTXYZD(1,NS3),
     %                    PTXYZD(1,NS5),
     %                    PTXYZD(1,NS4),
     %                    PTXYZD(1,NS1),
     %                    ARMIN, ARMAX, SURFTR, V1, Q1 )

            IF( Q1 .LT. QminGR ) THEN
               QminGR = Q1
            ENDIF

            CALL QUATETD( PTXYZD(1,NS3),
     %                    PTXYZD(1,NS4),
     %                    PTXYZD(1,NS5),
     %                    PTXYZD(1,NS2),
     %                    ARMIN, ARMAX, SURFTR, V2, Q2 )

            IF( Q2 .LT. QminGR ) THEN
               QminGR = Q2
            ENDIF
            VT12 = VT12 + V1 + V2

 45      ENDDO

C        SI LES VOLUMES SONT DIFFERENTS, PAS DE MODIFICATION
C        ---------------------------------------------------
         IF( ABS(VT0-VT12) .GE. VT0 * 1D-4 ) GOTO 9990

C        COMPARAISON DES QUALITES AVANT et APRES
C        ---------------------------------------
         IF( QminGR .LE. 0 .OR. QminGR .LE. QMIN0 ) GOTO 9990


C        ICI LE MIN DES QUALITES DES DEMI-TETRAEDRES EST AMELIORE
C        ET LES VOLUMES SONT EGAUX =>
C        DECOMPOSITION EFFECTIVE DES TETRAEDRES D'ARETE NS1-NS2
C        et SUPPRESSION DES 2 TETRAEDRES DE FACE NOSOTR=NS1-NS2-NS3
C        ==========================================================

ccc         tracte = .true.
ccc      KTITRE='teqmtytr: TETRAEDRE          de QUALITE=                 A
ccc     %RETE                 NS3=       NBTS1S2=         et TETRA OPPOSES'
ccc         WRITE( KTITRE(21:28),'(I8)'   ) NTEQM
ccc         WRITE( KTITRE(41:55),'(G15.6)') QTEQM
ccc         WRITE( KTITRE(64:69),'(I6)'   ) NS1
ccc         WRITE( KTITRE(71:76),'(I6)'   ) NS2
ccc         WRITE( KTITRE(84:89),'(I6)'   ) NS3
ccc         WRITE( KTITRE(101:106),'(I6)' ) NBTS1S2
ccc         CALL SANSDBL( KTITRE, L )
ccc         print*
ccc         PRINT*, KTITRE(1:L)
ccc         print*,'teqmtytr: NOTETR(',nteqm,')=',(NOTETR(kk,nteqm),kk=1,8)
ccc         CALL TRFETO13( KTITRE(1:L), PTXYZD, NBTS1S2, NOTS1S2, NOTETR )

C        RECENSEMENT DES TETRAEDRES OPPOSES DES TETRAEDRES D'ARETE NS1-NS2
         NBTE1AOP = 0
         DO K = 1, NBTS1S2

C           LE K-EME TETRAEDRE D'ARETE NS1-NS2
            NTEK = NOTS1S2( K )

            IF( NTEK .LE. 0 ) GOTO 9990
C           TABLEAU DES TETRAEDRES D'ARETE NS1 NS2 DEFAILLANT
            IF( NOTETR(1,NTEK) .EQ. 0 ) GOTO 9990

            DO 50 NF=1,4
               NTEOP = NOTETR( 4+NF, NTEK )
               IF( NTEOP .GT. 0 .AND. NOTETR(1,NTEOP) .GT. 0 ) THEN
C                 UN TETRAEDRE OPPOSE.

C                 EST IL UN DES TETRAEDRES A SUPPRIMER?
                  DO M=1,NBTS1S2
                     IF( NTEOP .EQ. NOTS1S2(M) ) GOTO 50
                  ENDDO

C                 NON: EST IL UN DES TETRAEDRES DEJA RECENSE?
                  DO M=1,NBTE1AOP
                     IF( NTEOP .EQ. NOTE1A(M) ) GOTO 50
                  ENDDO

C                 NTEOP EST AJOUTE A LA LISTE DES TETRAEDRES
C                 OPPOSES POUR ETRE PRIS EN COMPTE PAR mjopte
                  NBTE1AOP = NBTE1AOP + 1
                  NOTE1A( NBTE1AOP ) = NTEOP
ccc                  PRINT*,'teqmtytr: +tetra OPPOSE NOTETR(',NTEOP,')=',
ccc     %                   (NOTETR(kk,NTEOP),kk=1,8)

                  CALL NOFAOP( NF, NTEK, NOTETR, NFOP, NTEOP )
                  IF( NFOP .GT. 0 ) THEN
C                    LE TETRAEDRE OPPOSE A LA FACE NTEOP EST INCONNU
                     NOTETR( 4+NFOP, NTEOP ) = -1
                  ENDIF

               ENDIF
 50         ENDDO

         ENDDO


         IF( NBTEK .LE. 0 ) THEN

C           TOUS LES NBTS1S2 TETRAEDRES ONT LA FACE NS1-NS2-NS3
C           SI LEUR VOLUME EST FAIBLE ALORS ILS SONT SUPPRIMES
C           (POUR EVITER UNE PERTE DE VOLUME)
C           ---------------------------------------------------
ccc            IF( VT0 .LT. VOLTID * 0.01D0 ) THEN  2019/03/15
            IF( VT0 .LT. VOLTID * 1D-4 ) THEN

C              RECENSEMENT DES FACES FRONTIERE
               DO K = 1, NBTS1S2
C                 LE K-EME TETRAEDRE D'ARETE NS1-NS2
                  NTEK = NOTS1S2( K )
                  PRINT*,'teqmtytr: DETRUIT NOTETR(',NTEK,')=',
     %                   (NOTETR(kk,NTEK),kk=1,8),' FRONTIERE?...'
                  DO NF=1,4
                     NTEOP = NOTETR( 4+NF, NTEK )
                     IF( NTEOP .EQ. 0 ) THEN
C                       NF EST UNE FACE FRONTIERE
                        DO NF2=1,4
                           IF( NF2 .NE. NF ) THEN
                              CALL NOFAOP( NF2, NTEK, NOTETR,
     %                                     NFTOP2, NTEOP2 )
                              IF( NFTOP2 .GT. 0 ) THEN
C                                LA FACE NFTOP2 DE NTEOP2
C                                DEVIENT FRONTIERE
                                 NOTETR( 4+NFTOP2, NTEOP2 ) = 0
                              ENDIF
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO

               GOTO 80

            ENDIF

         ENDIF

C        LE TABLEAU NOTE1A EST IL SUFFISAMMENT GRAND POUR AJOUTER
C        LES NBTEK  2 DEMI-TETRAEDRES DECOUPES
C        --------------------------------------------------------
         IF( NBTE1AOP + 2*NBTEK .GT. MXTE1A ) THEN
            PRINT*,'teqmtytr: AUGMENTER MXTE1A=',MXTE1A
C           ABANDON
            GOTO 9990
         ENDIF

         NBTE1A = NBTE1AOP
         DO 70 K = 1, NBTS1S2

C           LE K-EME TETRAEDRE D'ARETE NS1-NS2
            NTEK = NOTS1S2( K )

C           NS4 NS5 LES 2 AUTRES SOMMETS DE NTEK
            DO L=1,4
               NS4 = NOTETR( L, NTEK )
               IF( NS4 .NE. NS1 .AND. NS4 .NE. NS2 ) GOTO 60
            ENDDO

 60         DO LL = L+1, 4
               NS5 = NOTETR( LL, NTEK )
               IF( NS5 .NE. NS1 .AND. NS5 .NE. NS2 ) GOTO 62
            ENDDO

C           NTEK NE DOIT PAS ETRE UN TETRAEDRE DE FACE NS1 NS2 NS3
 62         IF( NS4 .EQ. NS3 .OR. NS5 .EQ. NS3 ) GOTO 70

C           VOLUME ET QUALITE DU TETRAEDRE NTEK
            CALL QUATETD( PTXYZD(1,NOTETR(1,NTEK)),
     %                    PTXYZD(1,NOTETR(2,NTEK)),
     %                    PTXYZD(1,NOTETR(3,NTEK)),
     %                    PTXYZD(1,NOTETR(4,NTEK)),
     %                    ARMIN, ARMAX, SURFTR, V, Q )

C           POUR DETERMINER L'ORIENTATION DE NS4-NS5
C           VOLUME DU TETRAEDRE NS1 NS2 NS4 NS5 = NTEK
            V2 = VOLTET( PTXYZD(1,NS1), PTXYZD(1,NS2), 
     %                   PTXYZD(1,NS4), PTXYZD(1,NS5) )
            IF( V2 * V .LT. 0D0 ) THEN
               L   = NS4
               NS4 = NS5
               NS5 = L
            ENDIF
C           L'ORIENTATION DE NS1 NS2 NS4 NS5 EST CELLE DE NTEK

C           UNE DES FACES DE NTEK EST ELLE LA FACE NOSOTR?
            CALL NUFATRTE( NOSOTR, NOTETR(1,NTEK), NF123 )
C           OUI SI NF123>0, NON SI NF123=0

            IF( NF123 .EQ. 0 ) THEN

C              LE TETRAEDRE NTEK N'A PAS LE SOMMET NS3
C              CE TETRAEDRE NTEK    NS1 NS2 NS4 NS5 EST DECOUPE
C              EN 2 SOUS-TETRAEDRES NS3 NS5 NS4 NS1
C                                   NS3 NS4 NS5 NS2
C              ------------------------------------------------

C              DECLARATION DU PREMIER SOUS TETRAEDRE NTEK1 DE NTEK
               NTEK1   = N1TEVI
               NUDTETR = MAX( NUDTETR, NTEK1 )
               IF( IVOLTE .NE. 0 ) NVOLTE( NTEK1 ) = NOVOLU
               N1TEVI  = NOTETR( 5, N1TEVI )

 64            NOTETR( 1, NTEK1 ) = NS3
               NOTETR( 2, NTEK1 ) = NS5
               NOTETR( 3, NTEK1 ) = NS4
               NOTETR( 4, NTEK1 ) = NS1

C              VERIFICATION: VOLUME DU TETRAEDRE NTEK1
               V1 = VOLTET( PTXYZD(1,NOTETR(1,NTEK1)),
     %                      PTXYZD(1,NOTETR(2,NTEK1)),
     %                      PTXYZD(1,NOTETR(3,NTEK1)),
     %                      PTXYZD(1,NOTETR(4,NTEK1)) )
               IF( V1 .LT. 0D0 ) THEN
                  print*,'teqmtytr: PB V1=',V1,' NE DEVRAIT PAS ETRE <0'
                  M   = NS4
                  NS4 = NS5
                  NS5 = M
                  GOTO 64
               ENDIF

C              MISE A JOUR DU TABLEAU SOMMET -> 1 TETRAEDRE
               N1TETS( NS1 ) = NTEK1
ccc               N1TETS( NS3 ) = NTEK1
ccc               N1TETS( NS4 ) = NTEK1
ccc               N1TETS( NS5 ) = NTEK1

C              LES TETRAEDRES OPPOSES AUX 4 FACES DE NTEK1
               NOTETR( 6, NTEK1 ) = -1
               NOTETR( 7, NTEK1 ) = -1
               NOTETR( 8, NTEK1 ) = -1


C              DECLARATION DU SECOND SOUS TETRAEDRE NTEK2 DE NTEK
               NTEK2   = N1TEVI
               NUDTETR = MAX( NUDTETR, NTEK2 )
               IF( IVOLTE .NE. 0 ) NVOLTE( NTEK2 ) = NOVOLU
               N1TEVI  = NOTETR( 5, N1TEVI )

 66            NOTETR( 1, NTEK2 ) = NS3
               NOTETR( 2, NTEK2 ) = NS4
               NOTETR( 3, NTEK2 ) = NS5
               NOTETR( 4, NTEK2 ) = NS2

C              VERIFICATION: VOLUME DU TETRAEDRE NTEK2
               V2 = VOLTET( PTXYZD(1,NOTETR(1,NTEK2)),
     %                      PTXYZD(1,NOTETR(2,NTEK2)),
     %                      PTXYZD(1,NOTETR(3,NTEK2)),
     %                      PTXYZD(1,NOTETR(4,NTEK2)) )
               IF( V2 .LT. 0D0 ) THEN
                 print*,'teqmtytr: PB V2=',V2,' NE DEVRAIT PAS ETRE <0!'
                  M   = NS4
                  NS4 = NS5
                  NS5 = M
                  GOTO 66
               ENDIF

C              MISE A JOUR DU TABLEAU SOMMET -> 1 TETRAEDRE
               N1TETS( NS2 ) = NTEK2
               N1TETS( NS3 ) = NTEK2
               N1TETS( NS4 ) = NTEK2
               N1TETS( NS5 ) = NTEK2

C              LES TETRAEDRES OPPOSES AUX 4 FACES DE NTEK2
C              NTEK1 et NTEK2 ONT LA FACE NS5 NS4 NS3 COMMUNE
               NOTETR( 5, NTEK2 ) = NTEK1
               NOTETR( 5, NTEK1 ) = NTEK2

               NOTETR( 6, NTEK2 ) = -1
               NOTETR( 7, NTEK2 ) = -1
               NOTETR( 8, NTEK2 ) = -1

cccC              LE TETRAEDRE OPPOSE A LA FACE NS1 NS4 NS5 DE NTEK1
ccc               NOSOTR2( 1 ) = NS1
ccc               NOSOTR2( 2 ) = NS4
ccc               NOSOTR2( 3 ) = NS5
ccc               CALL NUFATRTE( NOSOTR2, NOTETR(1,NTEK), NF145 )
ccc               NOTETR( 6, NTEK1 ) = NOTETR( 4+NF145, NTEK )

cccC              LE TETRAEDRE OPPOSE A LA FACE NS2 NS4 NS5 DE NTEK2
ccc               NOSOTR2( 1 ) = NS2
ccc               NOSOTR2( 2 ) = NS4
ccc               NOSOTR2( 3 ) = NS5
ccc               CALL NUFATRTE( NOSOTR2, NOTETR(1,NTEK), NF245 )
ccc               NOTETR( 6, NTEK2 ) = NOTETR( 4+NF245, NTEK )


cccC              LE TETRAEDRE OPPOSE A LA FACE NS1 NS2 NS4 DE NTEK
ccc               NOSOTR2( 1 ) = NS1
ccc               NOSOTR2( 2 ) = NS2
ccc               NOSOTR2( 3 ) = NS4
ccc               CALL NUFATRTE( NOSOTR2, NOTETR(1,NTEK), NF124 )
ccc               NTOP124 = NOTETR( 4+NF124, NTEK )
ccc               IF( NTOP124 .EQ. 0 ) THEN
cccC                 NF124 EST UNE FACE FRONTIERE QUI ENTRAINE LA
cccC                 FACE 134 DE NTEK1 FRONTIERE
ccc                  NOTETR( 7, NTEK1 ) = 0
cccC                 FACE 234 DE NTEK2 FRONTIERE
ccc                  NOTETR( 8, NTEK2 ) = 0
ccc               ENDIF

cccC              LE TETRAEDRE OPPOSE A LA FACE NS1 NS2 NS5 DE NTEK
ccc               NOSOTR2( 1 ) = NS1
ccc               NOSOTR2( 2 ) = NS2
ccc               NOSOTR2( 3 ) = NS5
ccc               CALL NUFATRTE( NOSOTR2, NOTETR(1,NTEK), NF125 )
ccc               NTOP125 = NOTETR( 4+NF125, NTEK )
ccc               IF( NTOP125 .EQ. 0 ) THEN
cccC                 NF125 EST UNE FACE FRONTIERE QUI ENTRAINE LA
cccC                 FACE 135 DE NTEK1 FRONTIERE
ccc                  NOTETR( 8, NTEK1 ) = 0
cccC                 FACE 235 DE NTEK2 FRONTIERE
ccc                  NOTETR( 7, NTEK2 ) = 0
ccc               ENDIF

C              NTEK1 EST AJOUTE A NOTE1A (MXTE1A EST SUPERIEUR A LA DEMANDE)
               NBTE1A = NBTE1A + 1
               NOTE1A( NBTE1A ) = NTEK1

C              NTEK2 EST AJOUTE A NOTE1A
               NBTE1A = NBTE1A + 1
               NOTE1A( NBTE1A ) = NTEK2

            ENDIF

 70      ENDDO

C        DESTRUCTION DES NBTS1S2 TETRAEDRES NOTS1S2 DANS NOTETR
C        RISQUE DE PERDRE LES FAIBLES VOLUMES DES 2 TETRAEDRES
C        DE FACE NS1 NS2 NS3 S'ILS ONT AU MOINS UNE FACE FRONTIERE
C        ---------------------------------------------------------
 80      DO K = 1, NBTS1S2
            NTEK = NOTS1S2( K )
            NOTETR( 1, NTEK ) = 0
            NOTETR( 5, NTEK ) = N1TEVI
            N1TEVI = NTEK
            IF( IVOLTE .NE. 0 ) NVOLTE( NTEK ) = -1
         ENDDO

C        COMPLETION DES CHAINAGES DES TETRAEDRES OPPOSES PAR LES FACES
C        -------------------------------------------------------------
         CALL MJOPTE( NBTE1A, NOTE1A, N1TETS, NOTETR, NUDTETR,
     %                N1TEVI, PTXYZD, NBFANR )

         IF( NBFANR .GT. 0 ) THEN
C           AU MOINS UNE FACE N'A PAS DE TETRAEDRE OPPOSE
            PRINT*,'teqmtytr: PROBLEME ',NBFANR,
     %             ' FACES DE TETRAEDRE OPPOSE NON RETROUVE'
            PRINT*,'teqmtytr: FACES FRONTIERE?...'
            PRINT*
         ENDIF

C        COMPRESSION DE NOTE1A POUR GARDER LES SOUS-TETRAEDRES GENERES
         DO K=NBTE1AOP+1, NBTE1A
            NOTE1A( K-NBTE1AOP ) = NOTE1A( K )
         ENDDO
         NBTE1A = NBTE1A - NBTE1AOP

ccc         print*,'teqmtytr: LISTE DES TETRAEDRES CREES DE SOMMET',NS3
ccc         do k=1,nbte1a
ccc            nte = note1a( k )
ccc            print*,'teqmtytr: notetr(',nte,')=',(notetr(kk,nte),kk=1,8)
ccc         enddo
cccC        TRACE DES NBTE1A SOUS-TETRAEDRES GENERES
ccc         KTITRE='teqmtytr:            SOUS-TETRAEDRES QminGR=           
ccc     %      REMPLACENT               TETRAEDRES QMIN0=               '
ccc         WRITE( KTITRE(11:20),  '(I10)  ') NBTE1A
ccc         WRITE( KTITRE(46:60),  '(G15.6)') QminGR
ccc         WRITE( KTITRE(73:82),  '(I10)'  ) NBTS1S2
ccc         WRITE( KTITRE(104:118),'(G15.6)') QMIN
ccc         CALL SANSDBL( KTITRE, L )
ccc         PRINT*, KTITRE(1:L)
ccc         CALL TRFETO13( KTITRE(1:L), PTXYZD, NBTE1A, NOTE1A, NOTETR )

C        DECOMPOSITION REUSSIE DES TETRAEDRES NOTS1S2 D'ARETE NS1-NS2
         MODIFT = 1
         GOTO 9999

ccc      ENDIF


C     SORTIE SANS MODIFICATION DU TETRAEDRE ET DE SES VOISINS
C     -------------------------------------------------------
 9990 NBTE1A = 0
      MODIFT = 0


 9999 TRACTE = TRACTE0

ccc      IF( QTEQM .LT. QTEMED ) THEN
ccc         IF( MODIFT .NE. 0 ) THEN
ccc            PRINT*,'teqmtytr: TETRAEDRE',NTEQM,' Volume=',VTEQM,
ccc     %             ' Qualite Initiale=',QTEQM,' a ete AMELIOREE'
ccc         ELSE
ccc            PRINT*,'teqmtytr: TETRAEDRE',NTEQM,' Volume=',VTEQM,
ccc     %             ' Qualite=',QTEQM,' St:',(NOTETR(M,NTEQM),M=1,8),
ccc     %             ' SANS AMELIORATION de sa QUALITE'
ccc         ENDIF
ccc      ENDIF

      IF( MODIFT .GT. 0 ) THEN
         PRINT *,'teqmtytr: TETRAEDRE',NTEQM,' NOTETR:',NOTEQM,
     %' Vol=',VTEQM,' Qual=',QTEQM,' A ETE MODIFIE pour QminGR=',QminGR
      ENDIF

      RETURN
      END
