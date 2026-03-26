      SUBROUTINE TEQMTYTQ( QTEAME, NTEQM,  PTXYZD, NPSOFR,
     %                     INFACO, MXFACO, LEFACO, IVOLTE, NVOLTE,
     %                     MXTETR, N1TEVI, N1TETS, NOTETR, NUDTETR,
     %                     MXTNEW, NBTNEW, NOTNEW, MODIFT, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AMELIORER LA QUALITE DES TETRAEDRES AUTOUR DU TETRAEDRE NTEQM
C -----    DE MEDIOCRE QUALITE EN CHERCHANT DES TETRAEDRES OPPOSES AUX
C          FACES ET AYANT 1 SOMMET COMMUN
C          3 TETRAEDRES OPPOSES A NTEQM AVEC UN SOMMET NST5 COMMUN et
C            UN SOMMET NST1 DE TYPE BARYCENTRE DE NTEQM (4t->1t)
C     ou   2 TETRAEDRES OPPOSES A NTEQM AVEC UN SOMMET NST5 COMMUN
C            (3t->2t)

C ENTREES:
C --------
C QTEAME : QUALITE DES TETRAEDRES AU DESSOUS DE LAQUELLE UNE
C          AMELIORATION DES TETRAEDRES EST DEMANDEE
C NTEQM  : NUMERO NOTETR DU TETRAEDRE DE QUALITE MEDIOCRE A TRAITER
C MXTETR : NOMBRE DE TETRAEDRES DECLARABLES DANS NOTETR
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NPSOFR : NUMERO DE POSITION PAR RAPPORT A LA FRONTIERE DE CHAQUE POINT
C          LE SIGNE DEVIENT NEGATIF SI LE SOMMET EST DEPLACE
C          =  0 SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C                    LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -1 SI LE POINT EST SOMMET RECONNU TROP PROCHE
C          = -4 SI LE POINT EST SOMMET NON TROP PROCHE
C          = -3 SI LE POINT EST SOMMET REFUSE DANS LA TETRAEDRISATION
C          = -NPSOFR(I) SI POINT I DEPLACE SUR LA SURFACE
C             DANS LEUR SURFACE FERMEE INITIALE OU NO DE POINT INTERNE

C INFACO : = 0 PAS DE TABLEAU LEFACO NI CONTROLE DE FACE FRONTIERE
C          = 1 EXISTENCE DU TABLEAU LEFACO ET CONTROLE DE FACE FRONTIERE
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
C IVOLTE : =0 TABLEAU NVOLTE NON PRESENT (No VOLUME DE CHAQUE TETRAEDRE)
C          =1 TABLEAU NVOLTE     PRESENT

C MODIFIES:
C ---------
C NVOLTE : NUMERO DE VOLUME (1 a NBVOPA) DE CHAQUE TETRAEDRE
C          -1 SI VOLUME INCONNU (Exemple: les TETRAEDRES VIDES DE NOTETR)
C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NOTETR
C N1TETS : N1TETS(NS) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET NS
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C NUDTETR: NUMERO NOTETR DU DERNIER TETRAEDRE ACTIF DANS NOTETR

C SORTIES:
C --------
C MODIFT : = 0 AUCUNE MODIFICATION DE TETRAEDRES
C          = 1 DES TETRAEDRES ONT ETE MODIFIES
C IERR   : = 0 PAS D'ERREUR
C          = 1 SATURATION DES SOMMETS    du TABLEAU PTXYZD
C          = 3 SATURATION DES TETRAEDRES du TABLEAU NOTETR
C          =11 TETRAEDRE NTEQM AVEC 0 POUR PREMIER SOMMET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint PIERRE du PERRAY          Fevrier  2017
C MODIFS : ALAIN PERRONNET Saint PIERRE du PERRAY          Novembre 2018
C2345X7..............................................................012
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
      LOGICAL                          TRACTE0

ccc      CHARACTER*100     KTITRE
      DOUBLE PRECISION  PTXYZD(4,*)
      DOUBLE PRECISION  ARMIN, ARMAX, SURFTR(4)
      DOUBLE PRECISION  VOLET1, VOLET2, VTE, VTEQM, V,
     %                  V1253, V1254, V3451, V3452, VOLTID, TID

      INTEGER           NOTETR(8,MXTETR), N1TETS(*),NPSOFR(*),
     %                  LEFACO(1:11,0:MXFACO), NVOLTE(MXTETR),
     %                  NOTNEW(MXTNEW)
      INTEGER           NOSOTR(3), NOTEQM(8), NOSOTR2(3), NOTESU(5)

      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

      INTEGER           NOSOARTE(2,6)
      DATA              NOSOARTE / 1,2, 2,3, 3,1, 4,1, 4,2, 4,3 /

      INTEGER           NOSOARTEOP(2,6)
      DATA              NOSOARTEOP/ 4,3, 4,1, 4,2, 2,3, 3,1, 1,2 /

      INTEGER           NOAROP(6)
      DATA              NOAROP/ 6,   4,   5,   2,   3,   1 /

      TRACTE0 = TRACTE
      NBTNEW  = 0
      NBPASS  = 0
      MODIFT  = 0
      IERR    = 0

C     CARACTERISTIQUES DU TETRAEDRE NTEQM DE QUALITE MEDIOCRE
C     -------------------------------------------------------
      IF( NTEQM .LE. 0 ) THEN
         GOTO 9900
      ENDIF

      IF( NOTETR(1,NTEQM) .LE. 0 ) THEN
         PRINT*,'teqmtytq: PROBLEME Tetraedre VIDE',NTEQM,
     %          ' de St:',(NOTETR(I,NTEQM),I=1,8)
         PRINT*,'Abandon de l''AMELIORATION du TETRAEDRE NTEQM=',NTEQM
         GOTO 9900
      ENDIF

ccc      if( nteqm.eq.28733 ) then
ccc         print*,'TEQMTYTQ: nteqm=',nteqm,
ccc     %          'notetr(',nteqm,')=',(notetr(kk,nteqm),kk=1,8)
ccc         nt = 3
ccc         print*,'TEQMTYTQ: nteqm=',nteqm,
ccc     %          'notetr(',nt,')=',(notetr(kk,nt),kk=1,8)
ccc         nt = 12110
ccc         print*,'TEQMTYTQ: nteqm=',nteqm,
ccc     %          'notetr(',nt,')=',(notetr(kk,nt),kk=1,8)
ccc         print*
ccc         tracte = .true.
ccc      endif

C     VOLUME et QUALITE DU TETRAEDRE NTEQM
      CALL QUATETD( PTXYZD( 1, NOTETR(1,NTEQM) ),
     %              PTXYZD( 1, NOTETR(2,NTEQM) ),
     %              PTXYZD( 1, NOTETR(3,NTEQM) ),
     %              PTXYZD( 1, NOTETR(4,NTEQM) ),
     %              ARMIN, ARMAX, SURFTR, VTEQM, QTEQM )

ccc      PRINT *,'teqmtytq: Debut QTEAME=',QTEAME,
ccc     %' Tetraedre',NTEQM,' de St:',(NOTETR(I,NTEQM),I=1,8),
ccc     %' Q=',QTEQM,' V=',VTEQM

      IF( QTEQM .GE. QTEAME ) THEN
ccc         PRINT*,'teqmtytq: TETRAEDRE',NTEQM,' DE QUALITE',
ccc     %          QTEQM,' de VOLUME',VTEQM,' RETOUR SANS TRAITEMENT'
         GOTO 9900
      ENDIF

C     NUMERO DE VOLUME DU TETRAEDRE NTEQM
      NOVOLU = 0
      IF( IVOLTE .EQ. 0 ) THEN
         NOVOLU = -1
      ELSE
         NOVOLU= NVOLTE( NTEQM )
      ENDIF

cccC     TRACE DU TETRAEDRE NTEQM ET DES 4 TETRAEDRES OPPOSES
ccc      KTITRE='teqmtytq:               TETRAEDRE          de QUALITE=      
ccc     %    + 4 TETRAEDRES OPPOSES'
ccc      WRITE( KTITRE(33:40),'(I8)') NTEQM
ccc      WRITE( KTITRE(53:58),'(F6.3)') QTEQM
ccc      CALL SANSDBL( KTITRE, NBC )
ccc      PRINT*,KTITRE(1:NBC)
ccc      CALL TRFETO9( KTITRE(1:NBC), PTXYZD, NTEQM, NOTETR )


C     RECHERCHE DE CONFIGURATIONS SIMPLES DES TETRAEDRES OPPOSES
C     AUX 4 FACES DE NTEQM, NOTAMMENT DE SOMMETS COMMUNS
C     ----------------------------------------------------------
C     SAUVEGARDE DU NUMERO PTXYZD DES 4 SOMMETS DU TETRAEDRE NTEQM
C     ET DU NUMERO NOTETR DE SES 4 TETRAEDRES OPPOSES
      DO L = 1, 8
         NOTEQM( L ) = NOTETR( L, NTEQM )
      ENDDO

C     RECHERCHE D'UN COUPLE DE FACES DU TETRAEDRE NTEQM TELLES QUE
C     LEUR TETRAEDRE OPPOSE AIENT UN SOMMET OPPOSE COMMUN NST5
C     ------------------------------------------------------------
      DO 500 NF1 = 1, 3

C        NO DU TETRAEDRE OPPOSE A LA FACE NF1
         NTEOP1 = NOTEQM( 4+NF1 )
         IF( NTEOP1 .LE. 0 ) GOTO 500
         IF( NOTETR(1,NTEOP1) .EQ. 0 ) GOTO 500

C        NUMERO DES 3 SOMMETS DE LA FACE NF1 DE NTEQM
         DO J=1,3
            NOSOTR( J ) = NOTEQM( NOSOFATE(J,NF1) )
         ENDDO

         IF( INFACO .NE. 0 ) THEN
            CALL NULETR( NOSOTR, MXFACO, LEFACO,  NFLEFA )
            IF( NFLEFA .GT. 0 ) THEN
C              OUI: FACE DANS LEFACO -> POUR LA CONSERVER PAS D'ECHANGE
               GOTO 500
            ENDIF
         ENDIF

         IF( IVOLTE .NE. 0 ) THEN
            IF( NVOLTE( NTEOP1) .NE. NOVOLU ) THEN
C              LES 2 TETRAEDRES NTEQM et NTEOP1 NE SONT PAS DANS LE MEME VOLUME
               GOTO 500
            ENDIF
         ENDIF

C        NST5 NO DU SOMMET DU TETRAEDRE NTEOP1 NON DANS LA FACE NF1
         DO 15 I1=1,4

            NST5 = NOTETR( I1, NTEOP1 )
            DO J=1,3
               IF( NST5 .EQ. NOSOTR(J) ) GOTO 15
            ENDDO
C           NST5 SOMMET DE NTEOP1 N'EST PAS UN SOMMET DE LA
C           FACE NF1 DE NTEQM
            GOTO 20

 15      ENDDO
         PRINT*,'teqmtytq: ANOMALIE DE SOMMETS DANS NOTETR(',NTEQM,')=',
     %           NOTEQM
         GOTO 9900

C        RECHERCHE D'UN SECOND TETRAEDRE NTEOP2 OPPOSE A NTEQM
C        AYANT LE SOMMET NST5 DE NTEOP1
 20      DO 22 NF2 = NF1+1, 4

C           NO DU TETRAEDRE OPPOSE A LA FACE NF2
            NTEOP2 = NOTEQM( 4+NF2 )
            IF( NTEOP2 .GT. 0 .AND. NTEOP2 .NE. NTEOP1 ) THEN
               IF( NOTETR(1,NTEOP2) .EQ. 0 ) GOTO 22
               DO I2=1,4
                  NST6 = NOTETR( I2, NTEOP2 )
                  IF( NST6 .EQ. NST5 ) THEN

                     IF( IVOLTE .NE. 0 ) THEN
                        IF( NVOLTE( NTEOP2) .NE. NOVOLU ) THEN
C                          LES 3 TETRAEDRES NE SONT PAS DANS LE MEME VOLUME
                           GOTO 22
                        ENDIF
                     ENDIF

                     GOTO 40

                  ENDIF
               ENDDO
            ENDIF

 22      ENDDO
         GOTO 500


C        NTEOP1 NTEOP2 ONT UN SOMMET COMMUN NST5 N'APPARTENANT PAS A NTEQM
C        RECHERCHE D'UN EVENTUEL 3-EME TETRAEDRE OPPOSE A NTEQM
C        AYANT LE SOMMET NST5 DE NTEOP1 ET NTEOP2
C        -----------------------------------------------------------------
 40      DO 44 NF3 = NF2+1, 4

C           NO DU TETRAEDRE OPPOSE A LA FACE NF3
            NTEOP3 = NOTEQM( 4+NF3 )
            IF( NTEOP3 .GT. 0 .AND. NTEOP2 .NE. NTEOP1 .AND.
     %          NTEOP2 .NE. NTEOP2 ) THEN
               IF( NOTETR(1,NTEOP3) .EQ. 0 ) GOTO 44
               DO I3=1,4
                  NST6 = NOTETR( I3, NTEOP3 )
                  IF( NST6 .EQ. NST5 ) THEN

                     IF( IVOLTE .NE. 0 ) THEN
                        IF( NVOLTE( NTEOP3) .NE. NOVOLU ) THEN
C                          LES 4 TETRAEDRES NE SONT PAS DANS LE MEME VOLUME
                           GOTO 44
                        ENDIF
                     ENDIF

                     GOTO 60

                  ENDIF
               ENDDO
            ENDIF

 44      ENDDO

C        ICI: IL Y A 2 TETRAEDRES OPPOSES A 2 FACES AYANT UN MEME SOMMET NST5
C             IL N'Y A PAS 3 TETRAEDRES OPPOSES A 3 FACES DE NTEQM
C             AYANT UN MEME SOMMET NST5  C-A-D UN SOMMET DE NTEQM 
C             N'EST PAS INTERNE A UN GRAND TETRAEDRE FORME DE 3 TETRAEDRES OPPOSES
         GOTO 100


C        =============================================================================
C        NTEOP1 NTEOP2 NTEOP3 ONT UN SOMMET COMMUN NST1 N'APPARTENANT PAS A NTEQM
C        FORMANT UN TETRAEDRE NTEQM LES ENGLOBANT ET LE SOMMET CENTRAL NST1
C        S'IL EST SUR LA FRONTIERE, DISPARAITRAIT DE LA TETRAEDRISATION
C        ET DANS CE CAS LA TRANSFORMATION 4 TETRAEDRES -> 1 TETRAEDRE EST ABANDONNEE
C        PAR CONTRE SI LE SOMMET NST1 EST INTERNE (NPSOFR(NST1)=0) ALORS
C        LA TRANSFORMATION EST FAITE ET LE SOMMET NST1 DISPARAIT DE LA TETRAEDRISATION
C        =============================================================================
C        RECHERCHE DE L'ARETE COMMUNE A NF1 NF2
 60      CALL NSC2FATE( NF1, NF2, NS1, NS2 )

C        RECHERCHE DE L'ARETE COMMUNE A NF2 NF3
         CALL NSC2FATE( NF2, NF3, NS3, NS4 )

C        NS1 DEVIENT LE SOMMET COMMUN AUX 3 FACES NF1 NF2 NF3
C        NS1-NS2 ARETE COMMUNE A NF1 NF2
C        NS1-NS3 ARETE COMMUNE A NF2 NF3
         IF( NS3 .EQ. NS1 ) THEN
            NS3 = NS4
         ELSE IF( NS3 .EQ. NS2 ) THEN
            NS2 = NS1
            NS1 = NS3
            NS3 = NS4
         ELSE IF( NS4 .EQ. NS2 ) THEN
            NS2 = NS1
            NS1 = NS4
         ENDIF

C        NS4 SOMMET DE L'ARETE COMMUNE A NF1 NF3 et NON NS1
C        C'EST LE SOMMET NS4 NON NS1 2 3 DU TETRAEDRE NTEQM
         DO NS4 = 1, 4
            IF( NS4.NE.NS1 .AND. NS4.NE.NS2 .AND. NS4.NE.NS3 ) GOTO 62
         ENDDO

 62      NST1 = NOTEQM( NS1 )

C        NST1 EST INTERNE AU TETRAEDRE NST2 NST4 NST3 NST5
C        --------------------------------------------------
         IF( NPSOFR( NST1 ) .GT. 0 ) THEN
C           LE SOMMET NST1 EST SUR LA FRONTIERE
C           IL N'EST PAS RETIRE DE LA TETRAEDRISATION CAR SA RETETRAEDRISATION
C           POSE DES PROBLEMES SUPPLEMENTAIRES
            GOTO 500
         ENDIF

C        NST1 SOMMET INTERNE AU VOLUME, NON SUR LA FRONTIERE
C        VA DISPARAITRE DE LA TETRAEDRISATION LORS DE LA
C        CREATION DU TETRAEDRE NST 2 4 3 5 ET DE LA
C        DESTRUCTION DES TETRAEDRES NTEQM, NTEOP1, NTEOP2, NTEOP3
C        4 TETRAEDRES -> 1 TETRAEDRE
C        REMARQUE: SI IVOLTE=1 LES 4 TETRAEDRES ONT ETE TESTES
C                              COMME ETANT DANS UN MEME VOLUME
C                  SI INFACO=1 LES FACES DE SOMMET NST1 NE SONT
C                              PAS SUR LA FRONTIERE
C        =========================================================
         NST2 = NOTEQM( NS2 )
         NST3 = NOTEQM( NS3 )
         NST4 = NOTEQM( NS4 )

         IF(NST5 .EQ. NST2 .OR. NST5 .EQ. NST3 .OR. NST5 .EQ. NST4) THEN
            GOTO 500
         ENDIF

C        AFFICHAGE DES TETRAEDRES A SUPPRIMER
         print*,'teqmtytq: 4T->1T NST1=',NST1,' DISPARAIT'
         print*,'teqmtytq: 4T->1T Old NOTETR(',NTEQM,')=',
     %          (NOTEQM(kk),kk=1,8),' Q=',QTEQM,' REMPLACE'
         print*,'teqmtytq: 4T->1T Old NOTETR(',NTEOP1,')=',
     %          (NOTETR(kk,NTEOP1),kk=1,8),' SUPPRIME'
         print*,'teqmtytq: 4T->1T Old NOTETR(',NTEOP2,')=',
     %          (NOTETR(kk,NTEOP2),kk=1,8),' SUPPRIME'
         print*,'teqmtytq: 4T->1T Old NOTETR(',NTEOP3,')=',
     %          (NOTETR(kk,NTEOP3),kk=1,8),' SUPPRIME'
ccc         tracte = .true.

C        PROTECTION DE BOUCLE. V2345 2 FOIS NEGATIF... OUI CELA SE PRODUIT!
         NPASS = 0

C        VOLUME et QUALITE DU TETRAEDRE 2345
 66      CALL QUATETD( PTXYZD( 1, NST2 ),
     %                 PTXYZD( 1, NST3 ),
     %                 PTXYZD( 1, NST4 ),
     %                 PTXYZD( 1, NST5 ),
     %                 ARMIN, ARMAX, SURFTR, VTE, QTE )
         IF( VTE .GE. 0D0 ) GOTO 68

         IF( NPASS .EQ. 0 ) THEN
C           PERMUTATION POUR OBTENIR UN VOLUME POSISTIF
            N    = NST3
            NST3 = NST4
            NST4 = N
            NPASS= 1
            V    = VTE
            GOTO 66
         ELSE
C           VTE 2 FOIS NEGATIF -> ABANDON
            PRINT*,'teqmtytq: V2345 2 FOIS NEGATIF',V,VTE,' ->ABANDON'
            GOTO 500
         ENDIF

C        LE TETRAEDRE NTEQM DEVIENT LE TETRAEDRE NST 2 3 4 5
 68      NTE = NTEQM
         NOTETR( 1, NTE ) = NST2
         NOTETR( 2, NTE ) = NST3
         NOTETR( 3, NTE ) = NST4
         NOTETR( 4, NTE ) = NST5

C        VOLUME et QUALITE DU TETRAEDRE NTE
         CALL QUATETD( PTXYZD( 1, NOTETR(1,NTE) ),
     %                 PTXYZD( 1, NOTETR(2,NTE) ),
     %                 PTXYZD( 1, NOTETR(3,NTE) ),
     %                 PTXYZD( 1, NOTETR(4,NTE) ),
     %                 ARMIN, ARMAX, SURFTR, VTE, QTE )
         QMIN = QTE

C        LE SOMMET NST1 N'APPARTIENT PLUS A UN TETRAEDRE
         N1TETS( NST1 ) = 0

         IF( IVOLTE .NE. 0 ) NVOLTE( NTE ) = NOVOLU

C        RECHERCHE DU TETRAEDRE OPPOSE AUX 4 FACES DE NTE
         NBTESU = 0
         DO 70 NF = 1, 4

C           LES 3 SOMMETS DE LA FACE NF DE NTE
            NOSOTR( 1 ) = NOTETR( NOSOFATE(1,NF), NTE )
            NOSOTR( 2 ) = NOTETR( NOSOFATE(2,NF), NTE )
            NOSOTR( 3 ) = NOTETR( NOSOFATE(3,NF), NTE )

C           CETTE FACE NOSOTR EST ELLE UNE FACE DE NTEQM?
            CALL NUFATRTE( NOSOTR, NOTEQM, NFTE )
            IF( NFTE .GT. 0 ) THEN

C              LA FACE NF DE NTE A POUR TETRAEDRE OPPOSE NTEOP
               NTEOP = NOTEQM( 4+NFTE )

C              LA FACE NF0 DE NTEQM EST LA FACE NFTE DE NTE
               NOTETR( 4+NF, NTE ) = NTEOP

               IF( NTEOP .GT. 0 ) THEN
                  CALL NUFATRTE( NOSOTR, NOTETR(1,NTEOP), NFOP )
                  IF( NFOP .GT. 0 ) THEN
                     NOTETR( 4+NFOP, NTEOP ) = NTE
                  ENDIF
               ENDIF

               GOTO 70

            ENDIF

C           LA FACE NF DE NTE N'EST PAS UNE FACE DE NTEQM
C           MAIS C'EST UNE FACE DE L'UN DES 3 TETRAEDRES OPPOSES A NTEQM
            DO NF0 = 1, 4

C              LA FACE NF0 DE NTEQM A POUR TETRAEDRE OPPOSE NTEOP
               NTEOP = NOTEQM( 4+NF0 )

               IF( NTEOP .GT. 0 .AND. NOTETR(1,NTEOP) .GT. 0 ) THEN
                  CALL NUFATRTE( NOSOTR, NOTETR(1,NTEOP), NFOP )
                  IF( NFOP .GT. 0 ) THEN
C                    LA FACE NF DE NTE EST LA FACE NFOP DE NTEOP
                     NOTETR( 4+NF  , NTE   ) = NTEOP
                     NOTETR( 4+NFOP, NTEOP ) = NTE
C                    NTEOP UN TETRAEDRE A DETRUIRE ENSUITE
                     NBTESU = NBTESU + 1
                     NOTESU( NBTESU ) = NTEOP
                     GOTO 70
                  ENDIF
               ENDIF

            ENDDO

 70      ENDDO

         DO K = 1, NBTESU
C           DESTRUCTION NOTETR DE NTEOP
            NTEOP = NOTESU( K )
            NOTETR( 1, NTEOP ) = 0
            NOTETR( 5, NTEOP ) = N1TEVI
            N1TEVI = NTEOP
            IF( IVOLTE .NE. 0 ) NVOLTE( NTEOP ) = -1
         ENDDO

C        MISE A JOUR DE N1TETS
         DO K = 1, 4
            N1TETS( NOTETR(K,NTE) ) = NTE
         ENDDO
         N1TETS( NST1 ) = 0

C        RAPPEL NTE = NTEQM MODIFIE
         NBTNEW  = 1
         NOTNEW( 1 ) = NTE

         print*,'teqmtytq: 4T->1T New NOTETR(',NTE,')=',
     %          (NOTETR(kk,NTE),kk=1,8),' Q=',QTE,' AJOUTE'
         print*

         GOTO 9000


C        ============================================================
C        NST5 LE SOMMET I1 DE NTEOP1 EST LE SOMMET I2 DE NTEOP2 ET
C        IL N'APPARTIENT PAS AU TETRAEDRE NTEQM
C        NST1 NST2 LES 2 SOMMETS DE L'ARETE COMMUNE DES FACES NF1 NF2
C        NST3 NST4 LES 2 SOMMETS DE L'AUTRE DIAGONALE 
C        ESSAI D'ECHANGER NT12345 ( NTEQM NTOP1 NTOP2 ) PAR
C        (NT3451 NT3452) DE FACE COMMUNE NS345
C        3 TETRAEDRES -> 2 TETRAEDRES
C        REMARQUE: SI IVOLTE=1 LES 3 TETRAEDRES ONT ETE TESTES
C                              COMME ETANT DANS UN MEME VOLUME
C                  SI INFACO=1 LA FACE COMMUNE A NTEOP1 ET NTEOP2
C                              SI ELLE EST FRONTIERE DOIT ETRE GARDEE
C                              PAS SUR LA FRONTIERE
C        ============================================================
C        CALCUL DES 2 SOMMETS NS1 NS2 DE L'ARETE COMMUNE AUX 2 FACES NF1 NF2
C        DES 2 AUTRES SOMMETS NS3 NS4 DU TETRAEDRE
C        DES 2 AUTRES FACES   NF3 NF4 DU TETRAEDRE
 100     CALL NS4C2FATE( NF1, NF2, NS1, NS2, NS3, NS4, NF3, NF4 )

         NST1 = NOTEQM( NS1 )
         NST2 = NOTEQM( NS2 )
         NST3 = NOTEQM( NS3 )
         NST4 = NOTEQM( NS4 )

         IF( INFACO .NE. 0 ) THEN
C           LES 3 SOMMETS DE LA FACE COMMUNE A NTEOP1 et NTEOP2
            NOSOTR2(1) = NST1
            NOSOTR2(2) = NST2
            NOSOTR2(3) = NST5
            CALL NULETR( NOSOTR2, MXFACO, LEFACO,  NFLEFA )
            IF( NFLEFA .GT. 0 ) THEN
C              OUI: FACE DANS LEFACO -> POUR LA CONSERVER PAS D'ECHANGE
               GOTO 500
            ENDIF
         ENDIF

C        3 TETRAEDRES NOTESU(1:3) -> 2 TETRAEDRES
         NBTESU = 0
         IF( NTEOP1 .GT. 0 ) THEN
            NBTESU = NBTESU + 1
            NOTESU( NBTESU ) = NTEOP1
         ENDIF
         IF( NTEOP2 .GT. 0 ) THEN
            NBTESU = NBTESU + 1
            NOTESU( NBTESU ) = NTEOP2
         ENDIF
         NBTESU = NBTESU + 1
         NOTESU( NBTESU ) = NTEQM

C        LE TETRAEDRE NTEQM NTEOP1 NTEOP2 SONT ILS REMPLACABLES
C        PAR S3451 S3452 ?
C        ------------------------------------------------------
C        COMPARAISON DU VOLUME DES 2 CAS
         CALL QUATETD( PTXYZD( 1, NOTETR(1,NTEOP1) ),
     %                 PTXYZD( 1, NOTETR(2,NTEOP1) ),
     %                 PTXYZD( 1, NOTETR(3,NTEOP1) ),
     %                 PTXYZD( 1, NOTETR(4,NTEOP1) ),
     %                 ARMIN, ARMAX, SURFTR, V1253, Q1253 )

         CALL QUATETD( PTXYZD( 1, NOTETR(1,NTEOP2) ),
     %                 PTXYZD( 1, NOTETR(2,NTEOP2) ),
     %                 PTXYZD( 1, NOTETR(3,NTEOP2) ),
     %                 PTXYZD( 1, NOTETR(4,NTEOP2) ),
     %                 ARMIN, ARMAX, SURFTR, V1254, Q1254 )

         VOLET1 = VTEQM + V1253 + V1254
         IF( VOLET1 .LE. 0D0 ) GOTO 500
C        3 TETRAEDRES PLATS-> ABANDON

C        PROTECTION DE BOUCLE. V3451 2 FOIS NEGATIF... OUI CELA SE PRODUIT!
         NPASS = 0

 110     CALL QUATETD( PTXYZD( 1, NST3 ),
     %                 PTXYZD( 1, NST5 ),
     %                 PTXYZD( 1, NST4 ),
     %                 PTXYZD( 1, NST1 ),
     %                 ARMIN, ARMAX, SURFTR, V3451, Q3451 )
         IF( V3451 .GE. 0D0 ) GOTO 111

         IF( NPASS .EQ. 0 ) THEN
C           PERMUTATION POUR OBTENIR UN VOLUME POSISTIF
            N   = NS3
            NS3 = NS4
            NS4 = N
            N    = NST3
            NST3 = NST4
            NST4 = N
            NPASS= 1
            V    = V3451
            GOTO 110
         ELSE
C           V3451 2 FOIS NEGATIF -> ABANDON
            PRINT*,'teqmtytq: V3451 2 FOIS NEGATIF',V,V3451,' ->ABANDON'
            NBTESU = 0
            NBTNEW = 0
            GOTO 500
         ENDIF

 111     CALL QUATETD( PTXYZD( 1, NST3 ),
     %                 PTXYZD( 1, NST4 ),
     %                 PTXYZD( 1, NST5 ),
     %                 PTXYZD( 1, NST2 ),
     %                 ARMIN, ARMAX, SURFTR, V3452, Q3452 )

ccc         VOLET2 = V3451 + ABS( V3452 )
         VOLET2 = V3451 + V3452

         IF( VOLET2 .LE. 0D0 ) THEN
C           VOLUME TROP FAIBLE
            NBTESU = 0
            NBTNEW = 0
            GOTO 500
         ENDIF

C        TAILLE ARETE SOUHAITEE MOYENNE AUTOUR DES 5 SOMMETS
         TID = ( PTXYZD(4,NST1) + PTXYZD(4,NST2) + PTXYZD(4,NST3)
     %         + PTXYZD(4,NST4) + PTXYZD(4,NST5) ) / 5D0

C        VOLUME du TETRAEDRE IDEAL PAR RAPPORT A LA TAILLE
C        SOUHAITEE DES ARETES
         VOLTID = ( TID * TID * TID ) / 6D0

ccc      IF( ABS(VOLET1-VOLET2) .GT. VOLTID*1D-3 ) THEN
         IF( ABS(VOLET1-VOLET2) .GT. VOLTID*1D-4 ) THEN
C           VOLUMES TROP DIFFERENTS
            NBTESU = 0
            NBTNEW = 0
            GOTO 500
         ENDIF

C        LA QUALITE MINIMALE DES 2 TETRAEDRES 3451 3452
         QMIN = MIN( Q3451, Q3452 )

         IF( QMIN .LT. MIN( QTEQM, Q1253, Q1254 ) ) THEN
C           QUALITE TROP FAIBLE DES TETRAEDRES 3451 3452
            NBTESU = 0
            NBTNEW = 0
            GOTO 500
         ENDIF

C        AFFICHAGE DES TETRAEDRES A SUPPRIMER
         print*,'teqmtytq: 3T->2T Old NOTETR(',NTEQM,')=',
     %          (NOTETR(kk,NTEQM),kk=1,8),' Q=',QTEQM,' SUPPRIME'
         print*,'teqmtytq: 3T->2T Old NOTETR(',NTEOP1,')=',
     %          (NOTETR(kk,NTEOP1),kk=1,8),' Q=',Q1253,' SUPPRIME'
         print*,'teqmtytq: 3T->2T Old NOTETR(',NTEOP2,')=',
     %          (NOTETR(kk,NTEOP2),kk=1,8),' Q=',Q1254,' SUPPRIME'


C        CREATION DES 2 TETRAEDRES 3451 et 3452
C        --------------------------------------
         IF( N1TEVI .LE. 0 ) GOTO 9050
         NT3451  = N1TEVI
         N1TEVI  = NOTETR( 5, N1TEVI )
         NUDTETR = MAX( NUDTETR, NT3451 )
         IF( IVOLTE .NE. 0 ) NVOLTE( NT3451 ) = NOVOLU

         IF( N1TEVI .LE. 0 ) GOTO 9050
         NT3452  = N1TEVI
         N1TEVI  = NOTETR( 5, N1TEVI )
         NUDTETR = MAX( NUDTETR, NT3452 )
         IF( IVOLTE .NE. 0 ) NVOLTE( NT3452 ) = NOVOLU

C        LES 3 TETRAEDRES DANS LES 2 TETRAEDRES 3451 3452
C        NT3451 DEVIENT LE TETRAEDRE DE SOMMETS 3541
         NOTETR( 1, NT3451 ) = NST3
         NOTETR( 2, NT3451 ) = NST5
         NOTETR( 3, NT3451 ) = NST4
         NOTETR( 4, NT3451 ) = NST1

C        NT3452 DEVIENT LE TETRAEDRE DE SOMMETS 3452
         NOTETR( 1, NT3452 ) = NST3
         NOTETR( 2, NT3452 ) = NST4
         NOTETR( 3, NT3452 ) = NST5
         NOTETR( 4, NT3452 ) = NST2

C        TETRAEDRES OPPOSES A NT3451 et NT3452 PAR LA FACE 345 
         NOTETR( 5, NT3451 ) = NT3452
         NOTETR( 5, NT3452 ) = NT3451

C        STOCKAGE DES 2 TETRAEDRES CREES POUR UTILISER UNE BOUCLE
         NOTESU( 4 ) = NT3451
         NOTESU( 5 ) = NT3452

         DO M=4,5

C           NT3451 ou NT3452 CREE
            NT = NOTESU( M )

C           FACE 2 A 4 DE NT
            DO 150 K=2,4

C              LES 3 SOMMETS DE LA FACE K DE NT
               NOSOTR( 1 ) = NOTETR( NOSOFATE(1,K), NT )
               NOSOTR( 2 ) = NOTETR( NOSOFATE(2,K), NT )
               NOSOTR( 3 ) = NOTETR( NOSOFATE(3,K), NT )

               DO I = 1, NBTESU

C                 LE TETREDRE I A DISPARAITRE
                  NTE = NOTESU( I )
                  CALL NUFATRTE( NOSOTR, NOTETR(1,NTE), NF )
                  IF( NF .GT. 0 ) THEN
C                    LA FACE NF DE NTE EST LA FACE K DE NT

C                    LE TETRAEDRE OPPOSE A LA FACE NF DE NTE A SUPPRIMER
C                    DEVIENT LE TETRAEDRE OPPOSE A LA FACE K DE NT
                     NTEOPP = NOTETR( 4+NF, NTE )
                     IF( NTEOPP.GT.0 .AND. NOTETR(1,NTEOPP).GT. 0 ) THEN
                        CALL NUFATRTE( NOSOTR, NOTETR(1,NTEOPP), NFOPP )
C                       LA FACE NFOPP DE NTEOPP EST LA FACE K DE NT

                        IF( NFOPP .GT. 0 ) THEN
                           NOTETR( 4+NFOPP, NTEOPP ) = NT
                           NOTETR( 4+K,     NT     ) = NTEOPP
                        ELSE
C                          NE DOIT PAS ARRIVER....
                           PRINT*,'teqmtytq: ANOMALIE NOSOTR=',NOSOTR,
     %                            ' NOTETR(',NTEOPP,')=',
     %                            (NOTETR(L,NTEOPP),L=1,8)
                        ENDIF
                        GOTO 150

                     ELSE
C                       LA FACE K DE NT EST FRONTIERE
                        NOTETR( 4+K, NT ) = 0
                        GOTO 150
                     ENDIF

                  ENDIF

               ENDDO

               PRINT*,'teqmtytq: ANOMALIE FACE',K,' de NOTETR(',NT,')=',
     %               (NOTETR(L,NT),L=1,8),' DE TETRAEDRE OPPOSE INCONNU'

 150        ENDDO

         ENDDO

C        MISE A JOUR DE N1TETS
         N1TETS( NST1 ) = NT3451
         N1TETS( NST2 ) = NT3452
         N1TETS( NST3 ) = NT3451
         N1TETS( NST4 ) = NT3451
         N1TETS( NST5 ) = NT3451

C        SUPPRESSION DES NBTESU TETRAEDRES DU TABLEAU NOTETR
         DO K = 1, NBTESU
C           DESTRUCTION DE NTE DU TABLEAU NOTETR
            NTE = NOTESU( K )
            NOTETR(1,NTE) = 0
            NOTETR(5,NTE) = N1TEVI
            IF( IVOLTE .NE. 0 ) NVOLTE( NTE ) = -1
            N1TEVI = NTE
         ENDDO

C        POUR LA VISUALISATION DES 2 TETRAEDRES CREES
         NBTNEW  = 2
         NOTNEW( 1 ) = NT3451
         NOTNEW( 2 ) = NT3452

         print*,'teqmtytq: 3T->2T New NOTETR(',NT3451,')=',
     %          (NOTETR(kk,NT3451),kk=1,8),' Q=',Q3451,' CREE'
         print*,'teqmtytq: 3T->2T New NOTETR(',NT3452,')=',
     %          (NOTETR(kk,NT3452),kk=1,8),' Q=',Q3452,' CREE'

         GOTO 9000

C     FIN DE BOUCLE SUR LES FACES DE NTEQM
 500  ENDDO

    
C     =================================================================
C     IL N'EXISTE PAS 2 TETRAEDRES OPPOSES A 2 FACES DU TETRAEDRE NTEQM
C     AYANT UN MEME SOMMET NST5
C     LES 4 TETRAEDRES OPPOSES SONT SANS FACE COMMUNE
C     =================================================================
      GOTO 9900


C     AU MOINS UN TETRAEDRE NOTNEW A ETE CREE ou MODIFIE
C     --------------------------------------------------
 9000 MODIFT = 1
      IERR   = 0
      GOTO 9999


C     TAILLE INSUFFISANTE DU TABLEAU NOTETR
C     -------------------------------------
 9050 PRINT *,'teqmtytq: SATURATION DU TABLEAU NOTETR NTEQM=',NTEQM,
     %        ' NST5=',NST5
      IERR   = 3
      MODIFT = 0
      NBTNEW = 0
      GOTO 9999


C     PAS DE MODIFICATION DE TETRAEDRES
C     ---------------------------------
 9900 IERR   = 0
      MODIFT = 0
      NBTESU = 0
      NBTNEW = 0

C     POUR GARDER LE TRACE DES TETRAEDRES EN CAS D'ABANDON
ccc      PRINT *,'teqmtytq: Nst5=',Nst5,' NOTETR(',NTEQM,')=',
ccc     %        (NOTETR(I,NTEQM),I=1,8),' Q=',QTEQM,' V=',VTEQM,
ccc     %        ' PAS d''AMELIORATION de la QUALITE'
ccc      KTITRE='teqmtytq: Point                   TETRAEDRES PAS d''AMELIORA
ccc     %TION de la QUALITE'
ccc      WRITE(KTITRE(15:22),'(I7)') Nst5
ccc      WRITE(KTITRE(24:31),'(I7)') NBTESU
ccc      CALL TRACEETOILE( KTITRE, PTXYZD, NOTETR, Nst5, NBTESU, NOTESU )


C     BILAN FINAL
 9999 IF( MODIFT .NE. 0 ) THEN
         PRINT *,'teqmtytq: TETRAEDRE',NTEQM,' Q=',QTEQM,' V=',VTEQM,
     %           ' NOTETR:',NOTEQM,' A ETE MODIFIE pour QMIN=',QMIN

cccC     VISUALISATION FINALE DES NBTESU TETRAEDRES MODIFIES
ccc      TRACTE = .TRUE.
ccc      KTITRE='teqmtytq: NST5=                 NOUVEAUX TETRAEDRES AUTOUR du
ccc     % TETRAEDRE de QUALITE<       et FIN'
ccc      WRITE(KTITRE(14:20),'(I7)') NST5
ccc      WRITE(KTITRE(22:28),'(I7)') NBTNEW
ccc      WRITE(KTITRE(81:86),'(F6.3)') QTEAME
ccc      CALL SANSDBL( KTITRE, NBC )
ccc      PRINT*,KTITRE(1:NBC)
ccc      CALL TRACEETOILE( KTITRE, PTXYZD, NOTETR, NST5, NBTNEW, NOTNEW )

      ENDIF

ccc      if( nteqm.eq.28733 ) then
ccc         print*,'teqmtytq: nteqm=',nteqm,
ccc     %          'notetr(',nteqm,')=',(notetr(kk,nteqm),kk=1,8)
ccc         nt = 3
ccc         print*,'teqmtytq: nteqm=',nteqm,
ccc     %          'notetr(',nt,')=',(notetr(kk,nt),kk=1,8)
ccc         nt = 12110
ccc         print*,'teqmtytq: nteqm=',nteqm,
ccc     %          'notetr(',nt,')=',(notetr(kk,nt),kk=1,8)
ccc         print*
ccc         tracte = .true.
ccc      endif

      TRACTE = TRACTE0
      RETURN
      END
