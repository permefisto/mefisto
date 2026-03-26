      SUBROUTINE FAPETORE( QTEAME, NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                     MXFAPE, NBFAPE, NOFAPE,
     %                     MXFACO, LEFACO, N1FASC,
     %                     MXTETR, N1TEVI, N1TETS, NOTETR, NUDTETR,
     %                     MXTRCF, NOTRCF, NOSTCF,
     %                     MXARCF, N1ARCF, NOARCF,
     %                     MXFETO, NFETOI, NOSTIS,
     %                     NBABAN, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECUPERATION DE FACES PERDUES DE LA FRONTIERE
C -----    PAR L'ETOILE DE TOUS LES SOMMETS DU CF ET SOMMETS ISOLES
C          ET TETRAEDRISATION DU TORE DES FACES DE L'ETOILE

C ENTREES:
C --------
C QTEAME : QUALITE AU DESSOUS DE LAQUELLE UN TETRAEDRE N'EST PAS CREE ou
C          QUALITE DES TETRAEDRES AU DESSOUS DE LAQUELLE UNE
C          AMELIORATION DE LA QUALITE DES TETRAEDRES EST DEMANDEE
C NBSOMM : NOMBRE DE SOMMETS INITIAUX
C MXSOMM : NOMBRE DE SOMMETS DECLARABLES DANS PTXYZD
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NPSOFR : =  0 SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C              LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -4 SI LE POINT EST SOMMET NON TROP PROCHE PT OU FACE
C          = -1 SI LE POINT EST SOMMET TROP PROCHE PT OU FACE
C          = -3 SI LE POINT EST SOMMET REFUSE DANS LA TETRAEDRISATION
C          = -NPSOFR(I) SI POINT DEPLACE SUR LA SURFACE
C            DANS LEUR SURFACE FERMEE INITIALE ou NO DE POINT INTERNE
C MXFAPE : NOMBRE MAXIMAL DE FACES PERDUES DE LEFACO
C NBFAPE : NOMBRE DE FACES PERDUES DE LEFACO
C NOFAPE : NUMERO DES NBFAPE FACES PERDUES

C MXFACO : MAX DE FACES DECLARABLES DANS LEFACO
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          1:   =0 POUR UNE FACE VIDE
C          123: NO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          45:  NO (DANS NUVOPA 0 SINON) DU VOLUME1 , VOLUME2 DE LA FACE
C          678: NO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3

C          9: ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C             => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C             LEFACO(9,*) -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)

C          10: HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C              LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C              NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C              SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C              NF = LEFACO( 9, NF )  ...

C          11: >0  NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE,
C              =0  SINON

cccC          12: = NO FACEOC DE 1 A NBFACES D'OC
C N1FASC : N1FASC(I) NUMERO D'UN TRIANGLE DE LEFACO DE SOMMET I

C MXTETR : NOMBRE DE TETRAEDRES DECLARABLES DANS NOTETR
C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NOTETR
C N1TETS : N1TETS(NS) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET NS
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412

C MXTRCF : NOMBRE MAXIMAL DECLARABLE D'ARETES OU TRIANGLES
C          OU SOMMETS DANS L'ETOILE
C NOTRCF : TABLEAU DU NUMERO DANS LEFACO DES TRIANGLES PERDUS DU CF
C NOSTCF : TABLEAU DU NUMERO DANS PTXYZD DES SOMMETS DU CF
C MXARCF : MAXIMUM D'ARETES DECLARABLES DANS N1ARCF et NOARCF
C N1ARCF : TABLEAU (0:MXARCF) AUXILIAIRE
C          POINTE SUR LE DEBUT DES ARETES DE CHAQUE LIGNE FERMEE DU CF
C          0 POUR LES PLACES VIDES
C NOARCF : 1:NUMERO DU SOMMET DE L'ARETE DE LA LIGNE DU CONTOUR FERME
C          2:NUMERO DANS NOARCF DE L'ARETE SUIVANTE DU CF
C          3:NUMERO DANS LEFACO DU TRIANGLE ADJACENT OPPOSE A L'ARETE

C MXFETO : NOMBRE MAXIMAL DE FACES DECLARABLES DANS LE CHAINAGE NFETOI
C NFETOI : AU DEBUT VERSION1  FACES DE L'ETOILE
C          1: NUMERO DU TETRAEDRE DANS NOTETR AYANT CETTE FACE
C          2: NUMERO LOCAL AU TETRAEDRE DE LA FACE DE L'ETOILE
C             UN SIGNE NEGATIF INDIQUE UN TRAITEMENT EFFECTUE
C          3: NON UTILISE ICI
C          4: NUMERO DE CETTE FACE DANS LEFACO, 0 SI PAS DANS LEFACO
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES
C          ENSUITE VERSION2
C          1: NUMERO DU TETRAEDRE DANS NOTETR OPPOSE A CETTE FACE
C          2: NUMERO PTXYZD DU SOMMET 1 DE LA FACE
C          3: NUMERO PTXYZD DU SOMMET 2 DE LA FACE
C          4: NUMERO PTXYZD DU SOMMET 3 DE LA FACE
C             S1S2xS1S3 EST DIRIGE VERS L'INTERIEUR DE L'ETOILE
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES
C NOSTIS : NUMERO DES SOMMETS ISOLES N'APPARTENANT PAS AU CONTOUR FERME

C SORTIES :
C ---------
C NBSOMM  : NOMBRE DE SOMMETS FINAUX DANS PTXYZD
C NBABAN  : NOMBRE D'ABANDONS DE TRAITEMENT DES FACES PERDUES
C IERR    : 0 SI PAS D'ERREUR
C          1 ETOILE AVEC MOINS DE 4 FACES
C          2 SATURATION D'UN TABLEAU
C          3 3 SOMMETS NON DANS UN TETRAEDRE
C          9 VOLUME DE L'ETOILE AVANT et APRES DIFFERENTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint PIERRE du PERRAY           Janvier 2018
C2345X7..............................................................012
      PARAMETER        (MXTECR=400000,
     %                  MXETOI=4096, MXPILE=4096, MXSSET=1024 )
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))

      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE
      LOGICAL           TRACTE0

      DOUBLE PRECISION  PTXYZD(4,MXSOMM)
      INTEGER           NOTETR(8,MXTETR), N1TETS(MXSOMM)

      INTEGER           LEFACO(11,0:MXFACO),
     %                  NPSOFR(1:MXSOMM),
     %                  N1FASC(1:MXSOMM),
     %                  NOFAPE(1:MXFAPE),
     %                  N1ARCF(0:MXARCF),
     %                  NOARCF(1:3,1:MXARCF),
     %                  NAETOI(4,MXETOI),
     %                  NOTRCF(MXTRCF),
     %                  NOSTCF(MXTRCF),
     %                  NOSTIS(MXSOMM)

      INTEGER           NOTECR(MXTECR),   NOTEETPR(MXTECR),
     %                  NFETOI(5,MXFETO),
     %                  NO0FAR(3,MXETOI), N1SSET(MXSSET),
     %                  NOSTFAET(MXETOI)

      CHARACTER*200     KTITRE
      DOUBLE PRECISION  VOLET0, VOLET1, VOLUTE, VMOYEN,
     %                  ARMIN, ARMAX, SURFTR(4)
      REAL              QUALTE

      INTEGER           NOSOTR(3), NOSOTR1(3), NOSOTR2(3),
     %                  NOSOTE1(4), NOSOTE2(4), NSTEMX1(4), NSTEMX2(4)

      INTEGER           NOSOARTE(2,6), NOSOAROP(2,6), NOSOFATE(3,4)
      DATA              NOSOARTE / 1,2,  2,3,  3,1,  1,4,  2,4,  3,4 /
      DATA              NOSOAROP / 2,3,  1,4,  2,4,  2,3,  3,1,  1,2 /
      DATA              NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C     VRAI ==> TRACE DES TETRAEDRES DE L'ETOILE
      TRACTE = .FALSE.
      TRACTE0= .FALSE.

C     TRACE EVENTUEL DES FACES PERDUES DANS LEFACO
      CALL TRFAPE( NBFAPE, NOFAPE, MXFACO, LEFACO, PTXYZD )
C
      MNTEFA = 0
      MXTEFA = 1024
      NBABAN = 0

cccC     TRI SELON LES SURFACES DECROISSANTES DES FACES PERDUES
ccc      CALL TRIFAP( NBFAPE, NOFAPE, LEFACO, PTXYZD )


C     ========================================================
C     METHODE: TETRAEDRISER LES FACES PERDUES JUSQU'AUX 
C              FACES DE L'ETOILE EN TORE DES TETRAEDRES DE
C              SOMMETS DU CF ET ISOLES
C     ========================================================
C     ==========================================================
C     BOUCLE SUR LES FACES LEFACO PERDUES DANS LA TETREDRISATION
C     ==========================================================
      NFP = 0

 1    IERR = 0
      NFP  = NFP + 1
      IF( NFP .LE. NBFAPE ) THEN

C        NFLPER NUMERO DE LA FACE PERDUE DANS LEFACO A TRAITER
         NFLPER = NOFAPE( NFP )
         IF( NFLPER .LE. 0 ) GOTO 1

C        LA FACE NFLPER EST ELLE ENCORE PERDUE?
         IF( LEFACO(11,NFLPER) .GT. 0 ) GOTO 1

         IF( NFP .GT. 0 ) THEN
            TRACTE = .TRUE.
         ENDIF

C        PROTECTION DU NOMBRE ACTUEL DE SOMMETS PTXYZD
         IERR    = 0
         NBSOMM0 = NBSOMM
         NBSSET = 0
         NBTEETPR = 0
         NBTECR = 0
         NBPTIN = 0
         NBTRCF = 0
         NBSTIS = 0
         NOPASS = 0
         PRINT *


C        1) FORMATION DU PREMIER CONTOUR FERME FORME DES
C           TRIANGLES LEFACO PERDUS ADJACENTS A PARTIR DE NFLPER
C        -------------------------------------------------------
         CALL VDR1CF( -1.0,   NFLPER, NBSOMM, PTXYZD, N1TETS, NOTETR,
     %                NBFAPE, NOFAPE, MXFACO, LEFACO, N1FASC, NO0FAR,
     %                MXETOI, NAETOI,
     %                MXTRCF, NBTRCF, NOTRCF,
     %                MXARCF, NBCF,   N1ARCF, NOARCF,
     %                MXTRCF, NBSTCF, NOSTCF,
     %                MXSOMM, NBSTIS, NOSTIS,
     %                IERR )

         PRINT *
         PRINT *, 'fapetore: FACE PERDUE',NFP,' /',NBFAPE,
     %            ' LEFACO(',NFLPER,')=',(LEFACO(K,NFLPER),K=1,11),
     %            ' NBTRCF=',NBTRCF,' IERR=',IERR
         PRINT 10002
10002    FORMAT(240('='))

         IF( NBTRCF .LE. 0 .OR. IERR .NE. 0 ) GOTO 9900

C        2) LISTER LES TETRAEDRES DES SOMMETS DU CF OU ISOLES
C        ----------------------------------------------------
         NBTEETPR = 0
         DO NN = 1, NBSTCF+NBSTIS

C           LE SOMMET DES NBTRCF TRIANGLES DU CF
            IF( NN .LE. NBSTCF ) THEN
               NS = NOSTCF( NN )
            ELSE
               NS = NOSTIS( NN-NBSTCF )
            ENDIF

            CALL TETR1S( NS,    N1TETS,          NOTETR,
     %                   NBTNS, MXTECR-NBTEETPR, NOTEETPR(NBTEETPR+1),
     %                   IERR )
            NBTEETPR = NBTEETPR + NBTNS
C           NBTECPR: EN ENTREE NOMBRE DE TETRAEDRES DEJA RANGES DANS NOTEETPR
C                    EN SORTIE AJOUT DES TETRAEDRES AYANT NS COMME SOMMET
C                    = 0 SI SATURATION DU TABLEAU NOTEETPR
C                    =-1 SI UN SOMMET NS N'EST PAS UN SOMMET DE N1TETS(NS)

C           UNICITE DES NO DE TETRAEDRES DU TABLEAU NOTEETPR
            CALL UNITABL( NOTEETPR, NBTEETPR )

            IF( IERR .NE. 0 ) THEN
            print*,'fapetore: SATURATION du TABLEAU NOTEETPR NBTEETPR=',
     %                 NBTEETPR,' MXTEXF=',MXTECR,' IERR=',IERR
               GOTO 9900
            ENDIF

         ENDDO

C        VERIFICATION DES TETRAEDRES OPPOSES AUX FACES DES
C        TETRAEDRES DE SOMMET DU CF OU ISOLES
         CALL VEOPTE( NBTEETPR, NOTEETPR, NOTETR, PTXYZD, NBFANR )
         IF( NBFANR .GT. 0 ) THEN
            PRINT *,'fapetore: 2) NOMBRE INCORRECT DE TETRAEDRES OPPOSES
     %=',NBFANR
            GOTO 9900
         ENDIF

C        CONSTRUCTION NFETOI DE L'ETOILE DES TRIANGLES FACES SIMPLES
C        DES TETRAEDRES ENVELOPPANTS. LES FACES VUES 2 FOIS SONT ELIMINEES
         CALL CRFETOI1( NBTEETPR, NOTEETPR, NOTETR,
     %                  MXFETO, N1FEOC, N1FEVI, NFETOI )
ccc         CALL AFNFETOI( N1FEOC, NFETOI, NOTETR )

C        PASSAGE DU TABLEAU NFETOI DE L'ETOILE de la Version1 en la Version2
         CALL V12NFETOI( N1FEOC, NFETOI, NOTETR, NBFETO )

C        TRACE DES FACES SIMPLES DE L'ETOILE DU CF NFETOI
         CALL TRFETO2( PTXYZD, NOTETR, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR )

         IF( N1FEVI .EQ. -1 ) GOTO 9900
         IF( N1FEOC .LE.  0 ) GOTO 9900


C        3) RECHERCHE D'UN SOMMET QUI JOINT AUX FACES PERDUES
C           GENERENT DES TETRAEDRES DE VOLUME POSITIF
C        ----------------------------------------------------
C        CONSTRUCTION DU TABLEAU DU NO DES SOMMETS DES FACES NFETOI
         NBSTFAET = 0

         NF1 = N1FEOC
 4       IF( NF1 .GT. 0 ) THEN

            DO 5 N=1,3
               NS = ABS( NFETOI(1+N,NF1) )
               DO K=1,NBSTFAET
                  IF( NS .EQ. NOSTFAET(K) ) GOTO 5
               ENDDO
C              LE SOMMET EST AJOUTE
               NBSTFAET = NBSTFAET + 1
               NOSTFAET( NBSTFAET ) = NS
 5          ENDDO

C           LA FACE SUIVANTE DE L'ETOILE A ANALYSER
            NF1 = NFETOI(5,NF1)
            GOTO 4

         ENDIF

C        RECHERCHE DU SOMMET DE L'ETOILE QUI JOINT AUX FACES PERDUES
C        DONNE DES TETRAEDRES DE VOLUME POSITIF
         QUAMXMI = -2.0
         NSTMXMI = 0
         DO 8 N=1,NBSTFAET
            NS = NOSTFAET( N )
            QUAMIST = 2.0
            DO K = 1, NBTRCF
               NTR = NOTRCF( K )
               IF( NTR .GT. 0 ) THEN
C                 FACE DE LEFACO
C                 VOLUME ET QUALITE DU TETRAEDRE LEFACO(NTR)+NS
                  CALL QUATETD( PTXYZD(1,LEFACO(1,NTR)),
     %                          PTXYZD(1,LEFACO(2,NTR)),
     %                          PTXYZD(1,LEFACO(3,NTR)),
     %                          PTXYZD(1,NS),
     %                          ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
                  IF( QUALTE .LE. 0.0 ) GOTO 8
                  IF( QUALTE .LT. QUAMIST ) QUAMIST = QUALTE
               ENDIF
            ENDDO
            IF( QUAMIST .GT. QUAMXMI ) THEN
               QUAMXMI = QUAMIST
               NSTMXMI = NS
            ENDIF
 8       ENDDO

         IF( QUAMXMI .LE. 0.0 ) THEN
C           AUCUN SOMMET
            IERR = 9
            PRINT*,'fapetore:JONCTION SOMMET->FACES PERDUES A AMELIORER'
            print*,'fapetore: FINALEMENT ABANDON de la FACE PERDUE',NFP,
     %             ' avec NBTRCF=',NBTRCF,' et',
     %              NBSOMM-NBSOMM0,' SOMMETS AJOUTES. IERR=',IERR
            NBABAN = NBABAN + 1
            GOTO 1
         ENDIF

C        LE SOMMET NSTMXMI JOINT AUX NBTRCF TRIANGLES
C        DONNE NBTRCF TETRAEDRES DE VOLUME POSITIF
C        CREATION DE CES NBTRCF TETRAEDRES
         VOLET1 = 0D0
         NBTECR = 0
         DO K = 1, NBTRCF
            NTR = NOTRCF( K )
            IF( NTR .GT. 0 ) THEN
C              FACE DE LEFACO
               CALL AJTEET( LEFACO(1,NTR), LEFACO(2,NTR),
     %                      LEFACO(3,NTR), NSTMXMI,
     %                      -1, -1,  -1, -1,
     %                      NBSOMM,  PTXYZD, QTEAME,
     %                      MXTETR,  N1TEVI, NOTETR, 
     %                      NUDTETR, N1TETS,
     %                      MXTECR,  NBTECR, NOTECR, VOLET1,
     %                      VOLUTE,  QUALTE, IERR )
               IF( IERR .NE. 0 ) THEN
                  GOTO 9900
               ENDIF
            ENDIF
         ENDDO

C        RECHERCHE DES TETRAEDRES OPPOSES PAR LES FACES PARMI
C        LES SEULS NBTECR TETRAEDRES CREES
         DO NTC = 1, NBTECR-1

C           LE TETRAEDRE CREE
            NTE = NOTECR( NTC )

C           LES 4 FACES DU TETRAEDRE NTE SONT RECHERCHEES
C           DANS LES AUTRES TETRAEDRES DE LA SOUS-ETOILE
            DO 10 L=1,4

C              LA FACE L DU TETRAEDRE NTE
               IF( NOTETR(4+L,NTE) .LE. 0 ) THEN
C                 FACE DE TETRAEDRE OPPOSE INCONNU
C                 LES 3 SOMMETS DE LA FACE L DE NTE
                  DO K=1,3
                     NOSOTR2(K) = NOTETR( NOSOFATE(K,L), NTE )
                  ENDDO

C                 PARCOURS DES TETRAEDRES AU DELA DE NTC
                  DO NTC1 = NTC+1, NBTECR
                     NTE1 = NOTECR( NTC1 )
C                    LA FACE L DE NTE NOSOTR2 EST ELLE UNE FACE DE NTE1
                     CALL NUFATRTE( NOSOTR2, NOTETR(1,NTE1), LL )
                     IF( LL .GT. 0 ) THEN
C                       LA FACE L  DU TETRAEDRE NTE  EST 
C                       LA FACE LL DU TETRAEDRE NTE1
                        NOTETR(4+L ,NTE ) = NTE1
                        NOTETR(4+LL,NTE1) = NTE
                        GOTO 10
                     ENDIF
                  ENDDO

               ENDIF

 10         ENDDO

         ENDDO

C        ICI LES FACES DES TETRAEDRES 1 A NBTECR DE TETRAEDRE OPPOSE
C        -1 SONT PERIPHERIQUES

         KTITRE='fapetore: Tetraedres A PARTIR des        FACES DU CF'
         WRITE(KTITRE(35:39),'(I5)') NBTRCF
         CALL SANSDBL( KTITRE, L )
         CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECR, NOTECR, NOTETR )


C        4) RECHERCHE DE LA FACE PERIPHERIQUE DE SOMMET NSTMXMI
C        ET DE LA FACE DE L'ETOILE OFFRANT LA CREATION DE
C        2 TETRAEDRES DE BONNE QUALITE POUR CONSTRUIRE
C        UN CYLINDRE DE 3 TETRAEDRES MENANT D'UN TETRAEDRE
C        D'UNE FACE DU CF A LA FACE DE L'ETOILE CE QUI
C        TRANSFORME LE "TORE" (A UN SOMMET PRES) EN UNE "PATATE"
C        TROUEE D'UN CYLINDRE CONNEXE
C        -------------------------------------------------------
         NFETMX = 0
         NTETMX = 0
         NFCFMX = 0
         QMXMX  = -2.0
         DO NTC = 1, NBTECR

C           LE TETRAEDRE CREE
            NTE = NOTECR( NTC )

            DO 30 L=1,4

C              LA FACE L DU TETRAEDRE NTE
               IF( NOTETR(4+L,NTE) .LT. 0 ) THEN

C                 FACE DE TETRAEDRE OPPOSE INCONNU DONC PERIPHERIQUE
C                 LES 3 SOMMETS DE LA FACE L DE NTE
                  NSINCF = 0
                  DO K=1,3
                     NS = NOTETR( NOSOFATE(K,L), NTE )
                     IF( NS .EQ. NSTMXMI )  NSINCF=K
                     NOSOTR2(K) = NOTETR( NOSOFATE(K,L), NTE )
                  ENDDO
                  IF( NSINCF .EQ. 0 ) GOTO 30

C                 LE SOMMET NSTMXMI EST LE PREMIER DE NOSOTR
                  GOTO( 11, 12, 13 ), NSINCF

 11               NOSOTR(1) = NOSOTR2(1)
                  NOSOTR(2) = NOSOTR2(2)
                  NOSOTR(3) = NOSOTR2(3)
                  GOTO 14

 12               NOSOTR(1) = NOSOTR2(2)
                  NOSOTR(2) = NOSOTR2(3)
                  NOSOTR(3) = NOSOTR2(1)
                  GOTO 14

 13               NOSOTR(1) = NOSOTR2(3)
                  NOSOTR(2) = NOSOTR2(1)
                  NOSOTR(3) = NOSOTR2(2)

C                 NOSOTR2 FACE PERIPHERIQUE DE NORMALE
C                 VERS L'EXTERIEUR DU TETRAEDRE
C                 DONC VERS L'INTERIEUR DU DOMAINE A MAILLER

C                 BOUCLE SUR LES FACES DE L'ETOILE HORS
C                 LES TETRAEDRES DE FACES DU CF DE SOMMET NSTMXMI
 14               NF1 = N1FEOC
 18               IF( NF1 .GT. 0 ) THEN

C                    NO DES 3 SOMMETS DE LA FACE NF1
C                    NORMALE VERS L'INTERIEUR DE L'ETOILE A MAILLER
                     NSINET = 0
                     DO K=1,3
                        NS = NFETOI( 1+K, NF1 )
                        IF( NS .EQ. NSTMXMI ) NSINET=K
                        NOSOTR2(K) = NS
                     ENDDO
                     IF( NSINET .EQ. 0 ) GOTO 26

C                    NF1 EST UNE FACE DE L'ETOILE DE SOMMET NSTMXMI
C                    LE SOMMET NSTMXMI EST LE PREMIER DE NOSOTR2
                     GOTO( 21, 22, 23 ), NSINET

 21                  NOSOTR1(1) = NOSOTR2(1)
                     NOSOTR1(2) = NOSOTR2(2)
                     NOSOTR1(3) = NOSOTR2(3)
                     GOTO 24

 22                  NOSOTR1(1) = NOSOTR2(2)
                     NOSOTR1(2) = NOSOTR2(3)
                     NOSOTR1(3) = NOSOTR2(1)
                     GOTO 24

 23                  NOSOTR1(1) = NOSOTR2(3)
                     NOSOTR1(2) = NOSOTR2(1)
                     NOSOTR1(3) = NOSOTR2(2)

C                    CONSTRUCTION DES 2 TETRAEDRES DE LA FACE
C                    NOSOTR PERIPHERIQUE A LA FACE NOSOTR1 DE L'ETOILE
C                    CONFIGURATION 1 DES 2 TETRAEDRES
 24                  NOSOTE1(1) = NOSOTR(1)
                     NOSOTE1(2) = NOSOTR(2)
                     NOSOTE1(3) = NOSOTR(3)
                     NOSOTE1(4) = NOSOTR1(3)

                     NOSOTE2(1) = NOSOTR1(1)
                     NOSOTE2(2) = NOSOTR1(2)
                     NOSOTE2(3) = NOSOTR1(3)
                     NOSOTE2(4) = NOSOTR(3)

C                    VOLUME ET QUALITE DU TETRAEDRE NOSOTE1
                     CALL QUATETD( PTXYZD(1,NOSOTE1(1)),
     %                             PTXYZD(1,NOSOTE1(2)),
     %                             PTXYZD(1,NOSOTE1(3)),
     %                             PTXYZD(1,NOSOTE1(4)),
     %                             ARMIN, ARMAX, SURFTR, VOLUTE, Q1 )
 
                     CALL QUATETD( PTXYZD(1,NOSOTE2(1)),
     %                             PTXYZD(1,NOSOTE2(2)),
     %                             PTXYZD(1,NOSOTE2(3)),
     %                             PTXYZD(1,NOSOTE2(4)),
     %                             ARMIN, ARMAX, SURFTR, VOLUTE, Q2 )

                     Q12 = MIN( Q1, Q2 )

C                    CONFIGURATION 2 DES 2 TETRAEDRES
                     NOSOTE1(4) = NOSOTR1(2)
                     NOSOTE2(4) = NOSOTR(2)

C                    VOLUME ET QUALITE DU TETRAEDRE NOSOTE1
                     CALL QUATETD( PTXYZD(1,NOSOTE1(1)),
     %                             PTXYZD(1,NOSOTE1(2)),
     %                             PTXYZD(1,NOSOTE1(3)),
     %                             PTXYZD(1,NOSOTE1(4)),
     %                             ARMIN, ARMAX, SURFTR, VOLUTE, Q3 )
 
                     CALL QUATETD( PTXYZD(1,NOSOTE2(1)),
     %                             PTXYZD(1,NOSOTE2(2)),
     %                             PTXYZD(1,NOSOTE2(3)),
     %                             PTXYZD(1,NOSOTE2(4)),
     %                             ARMIN, ARMAX, SURFTR, VOLUTE, Q4 )

                     Q34 = MIN( Q3, Q4 )

                     IF( Q12 .LE. 0 .AND. Q34 .LE. 0 ) GOTO 26

C                    MEILLEURE CONFIGURATION?
                     IF( Q12 .GE. Q34 ) THEN
                        QMX = Q12
                        NOSOTE1(4) = NOSOTR1(3)
                        NOSOTE2(4) = NOSOTR(3)
                     ELSE
                        QMX = Q34
                     ENDIF

                     IF( QMX .GT. QMXMX ) THEN
C                       STOCKAGE DES 2 MEILLEURS TETRAEDRES
                        DO M=1,4
                           NSTEMX1(M) = NOSOTE1(M)
                           NSTEMX2(M) = NOSOTE2(M)
                        ENDDO
                        QMXMX  = QMX
                        NFETMX = NF1
                        NTETMX = NTE
                        NFCFMX = L
                     ENDIF

C                    LA FACE SUIVANTE DE L'ETOILE A ANALYSER
 26                  NF1 = NFETOI( 5, NF1 )
                     GOTO 18

                  ENDIF

               ENDIF

 30         ENDDO

         ENDDO

cccC                    LE TETRAEDRE NOSOTR+NST INTERSECTE T IL UN
cccC                    DES TETRAEDRES CREES AVANT?
ccc                     DO K=1,NBTECR
ccc                        NT = NOTECR( K )
ccc                        CALL INTETTET( PTXYZD(1,NOSOTR(1)),
ccc     %                                 PTXYZD(1,NOSOTR(2)),
ccc     %                                 PTXYZD(1,NOSOTR(3)),
ccc     %                                 PTXYZD(1,NST),
ccc     %                                 PTXYZD(1,NOTETR(1,NT)),
ccc     %                                 PTXYZD(1,NOTETR(2,NT)),
ccc     %                                 PTXYZD(1,NOTETR(3,NT)),
ccc     %                                 PTXYZD(1,NOTETR(4,NT)),
ccc     %                                 LINTER )
ccc                        IF( LINTER .GT. 0 ) THEN
cccC                          LE TETRAEDRE NOSOTR+NST INTERSECTE NT
ccc                           GOTO 16
ccc                        ENDIF
ccc                     ENDDO

C        NUMERO NOTETR DU TETRAEDRE OPPOSE A LA FACE 1 DE NSTEMX2
         NTEOP1 = NOTETR( 4+NFCFMX, NTETMX )
         CALL AJTEET( NSTEMX1(1), NSTEMX1(2), NSTEMX1(3), NSTEMX1(4),
     %                NTEOP1,     -1,         -1,         -1,
     %                NBSOMM,  PTXYZD, QTEAME,
     %                MXTETR,  N1TEVI, NOTETR, 
     %                NUDTETR, N1TETS,
     %                MXTECR,  NBTECR, NOTECR, VOLET1,
     %                VOLUTE,  QUALTE, IERR )
         IF( IERR .NE. 0 ) THEN
            GOTO 9900
         ENDIF
C        LE NOUVEAU TETRAEDRE CREE AVEC UNE FACE D'UN TETRAEDRE DU CF
         NTECF = NOTECR( NBTECR )

C        NUMERO NOTETR DU TETRAEDRE OPPOSE A LA FACE 1 DE NSTEMX2
         NTEOP2 = NFETOI( 1, NFETMX )
         CALL AJTEET( NSTEMX2(1), NSTEMX2(2), NSTEMX2(3), NSTEMX2(4),
     %                NTEOP2,     -1,         -1,         -1,
     %                NBSOMM,  PTXYZD, QTEAME,
     %                MXTETR,  N1TEVI, NOTETR, 
     %                NUDTETR, N1TETS,
     %                MXTECR,  NBTECR, NOTECR, VOLET1,
     %                VOLUTE,  QUALTE, IERR )
         IF( IERR .NE. 0 ) THEN
            GOTO 9900
         ENDIF
C        LE NOUVEAU TETRAEDRE CREE AVEC UNE FACE DE L'ETOILE
         NTEET = NOTECR( NBTECR )

         KTITRE='fapetore: 2 TETRAEDRES AVEC 1 FACE TETAEDRE CF OPPOSE e
     %t 1 FACE de l''ETOILE'
         CALL SANSDBL( KTITRE, L )
         PRINT*,KTITRE(1:L)
         PRINT*,'NOTETR1(',NTECF,')=',(NOTETR(L,NTECF),L=1,8)
         PRINT*,'NOTETR2(',NTECF,')=',(NOTETR(L,NTEET),L=1,8)

         CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECR, NOTECR, NOTETR )

C        RECHERCHE DES TETRAEDRES OPPOSES PAR LES FACES PARMI
C        LES NBTECR TETRAEDRES CREES
         DO NTC = 1, NBTECR-1

C           LE TETRAEDRE CREE
            NTE = NOTECR( NTC )

C           LES 4 FACES DU TETRAEDRE NTE SONT RECHERCHEES
C           DANS LES AUTRES TETRAEDRES DE LA SOUS-ETOILE
            DO 33 L=1,4

C              LA FACE L DU TETRAEDRE NTE
               NTEOP = NOTETR(4+L,NTE)
               IF( NTEOP .LE. 0 ) THEN

C                 FACE DE TETRAEDRE OPPOSE INCONNU
C                 LES 3 SOMMETS DE LA FACE L DE NTE
                  DO K=1,3
                     NOSOTR2(K) = NOTETR( NOSOFATE(K,L), NTE )
                  ENDDO

C                 PARCOURS DES TETRAEDRES AU DELA DE NTC DE NOTECR
                  DO NTC1 = NTC+1, NBTECR
                     NTE1 = NOTECR( NTC1 )
C                    LA FACE L DE NTE NOSOTR2 EST ELLE UNE FACE DE NTE1
                     CALL NUFATRTE( NOSOTR2, NOTETR(1,NTE1), LL )
                     IF( LL .GT. 0 ) THEN
C                       LA FACE L  DU TETRAEDRE NTE  EST 
C                       LA FACE LL DU TETRAEDRE NTE1
                        NOTETR(4+L ,NTE ) = NTE1
                        NOTETR(4+LL,NTE1) = NTE
                        GOTO 33
                     ENDIF
                  ENDDO

               ELSE

C                 A SUPPRIMER APRES MISE AU POINT?
C                 VERIFICATION DE L'OPPOSITION DES 2 TETRAEDRES
C                 LA FACE L DE NTE NOSOTR2 EST ELLE UNE FACE DE NTEOP?
                  DO K=1,3
                     NOSOTR2(K) = NOTETR( NOSOFATE(K,L), NTE )
                  ENDDO
                  CALL NUFATRTE( NOSOTR2, NOTETR(1,NTEOP), LL )
                  IF( LL .GT. 0 ) THEN
C                    LA FACE L  DU TETRAEDRE NTE  EST 
C                    LA FACE LL DU TETRAEDRE NTEOP
                     NOTETR(4+L ,NTE  ) = NTEOP
                     NOTETR(4+LL,NTEOP) = NTE
                  ELSE
C                    LA FACE L DU TETRAEDRE NTE N'EST PAS UNE FACE DE NTEOP
                     PRINT*,'fapetore: LA FACE ',L,' du TETRAEDRE',NTE,
     %                      ' N''est PAS une FACE du TETRAEDRE',NTEOP
                     PRINT*,'NOTETR(',NTE,  ')=',(NOTETR(J,NTE)  ,J=1,8)
                     PRINT*,'NOTETR(',NTEOP,')=',(NOTETR(J,NTEOP),J=1,8)
                     GOTO 33
                  ENDIF

               ENDIF

 33         ENDDO

         ENDDO

         KTITRE='fapetore: TETRAEDRES CF + TETRAEDRE CF OPPOSE + TETRAED
     %RE FACE ETOILE'
         CALL SANSDBL( KTITRE, L )
         CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECR, NOTECR, NOTETR )

C        DANS L'ETOILE RETRAIT OU AJOUT DES FACES DES NBTECR TETRAEDRES CREES
ccc         CALL AFNFETOI( N1FEOC, NFETOI, NOTETR )
         NBFASU = 0
         NBFAAJ = 0
         DO NTC = 1, NBTECR

C           LE TETRAEDRE CREE
            NTE = NOTECR( NTC )
            PRINT*,'NOTETR(',NTE,')=',(NOTETR(L,NTE),L=1,8)

            DO J=1,4
C              LA FACE J DU TETRAEDRE NTE EST DANS NFETOI
C              SOIT RETIREE, SOIT AJOUTEE EN FIN DE CHAINAGE
               CALL AJFAET2( NTE,    J,      NOTETR,
     %                       N1FEOC, N1FEVI, NFETOI,
     %                       NFS ,   NBFASU, NBFAAJ )
C              NFS: >0 NUMERO DANS NFETOI DE LA FACE AJOUTEE
C                   =0 SUPPRESSION DE LA FACE N'APPARTENANT PAS A LEFACO
C                   <0 -NUMERO LEFACO DE LA FACE
               IF( N1FEVI .EQ. -1 ) THEN
                  PRINT*,'fapetore: TABLEAU NFETOI SATURE'
                  GOTO 9900
               ENDIF

ccc               IF( NFS .GT. 0 ) THEN
cccC                 INVERSION DE LA NORMALE A LA FACE AJOUTEE POUR
cccC                 DESIGNER L'EXTERIEUR DES TETRAEDRES NSTMXMI-FACE PERDUE
ccc                  N = NFETOI( 3, NFS )
ccc                  NFETOI(3,NFS) =  NFETOI(4,NFS)
ccc                  NFETOI(4,NFS) =  N
ccc               ENDIF

            ENDDO

         ENDDO

C        INVERSION DU CHAINAGE DES FACES DE L'ETOILE POUR AVOIR EN PREMIER
C        LES NOUVELLES FACES CREES ET FACILITER LEUR TETRAEDRISATION
         CALL INVNFETOI( N1FEOC, NFETOI )

         KTITRE='fapetore: Tetraedres St                 Faces Perdues'
         WRITE(KTITRE(25:31),'(I7)') NSTMXMI
         WRITE(KTITRE(33:39),'(I7)') NBTRCF
         CALL SANSDBL( KTITRE, L )
         CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECR, NOTECR, NOTETR )


C        5) L'ETOILE A UNE FORME D'UNE PATATE AVEC UN PUITS DE
C           NBTRCF + 2 TETRAEDRES
C           RECENSEMENT DES SOUS ETOILES DE L'ETOILE N1FEOC
C        -----------------------------------------------------
Cccc        CREATION DES SOUS-ETOILES DE L'ETOILE
ccc         CALL SOUSETO( PTXYZD, NOTETR, N1FEOC, NFETOI,
ccc     %                 MXSSET, NBSSET, N1SSET, IERR )

C        CREATION DE L'UNIQUE SOUS-ETOILE=ETOILE
C        DE TYPE PATATE CONTENANT UN PRISME DE SOMMET SUR ELLE
         NBSSET = 1
         N1SSET( 1 ) = N1FEOC

         IF( IERR .NE. 0 ) THEN
            PRINT*,'fapetore: SORTIE souseto AVEC IERR='
     %            ,IERR,' NBSSET=',NBSSET
            GOTO 9900
         ENDIF

C        CALCUL DU VOLUME DES TETRAEDRES DU TABLEAU NOTEETPR
C        C-A-D DES TETRAEDRES DE L'ETOILE AVANT INTEGRATION
C        DES FACES PERDUES
         VOLET0 = 0D0
         QUAMIN = 2.0
         QUAMOY = 0.0

C        N1TETS SERA MIS A JOUR AVEC LES NOUVEAUX TETRAEDRES
ccc         print*,'fapetore: LISTE des TETRAEDRES A RE-TETRAEDRISER'

         DO K = 1, NBTEETPR

C           NUMERO NOTETR DU K-EME TETRAEDRE DE L'ETOILE
            NTE = NOTEETPR( K )

            IF( NTE .GT. 0 ) THEN
C              SON VOLUME ET SA QUALITE
               CALL QUATETD( PTXYZD(1,NOTETR(1,NTE)),
     %                       PTXYZD(1,NOTETR(2,NTE)),
     %                       PTXYZD(1,NOTETR(3,NTE)),
     %                       PTXYZD(1,NOTETR(4,NTE)),
     %                       ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
ccc            IF( VOLUTE .LE. 0D0 ) THEN
ccc               print*,'fapetore: TETRAEDRE INITIAL',NTE,
ccc     %                ' St:',(NOTETR(L,NTE),L=1,4),
ccc     %                ' de VOLUME=',VOLUTE,' Qualite=',QUALTE
ccc            ENDIF
               VOLET0 = VOLET0 + ABS( VOLUTE )
               QUAMOY = QUAMOY + QUALTE
               IF( QUALTE .LT. QUAMIN ) THEN
                  QUAMIN = QUALTE
               ENDIF
            ENDIF

         ENDDO

C        VOLUME MOYEN DES TETRAEDRES A RE-TETRAEDRISER
         VMOYEN = VOLET0 / NBTEETPR
         QUAMOY = QUAMOY / NBTEETPR
ccc         print *,'fapetore: FACE PERDUE',NFP,' NBTEETPR=',NBTEETPR,
ccc     %           ' VOLET0=',VOLET0,' VMOYEN0=',VMOYEN,
ccc     %           ' QualMOY0=',QUAMOY,' QualMIN0=',QUAMIN

C        TRACE DES NBTEETPR TETRAEDRES A RE-TETRAEDRISER
         KTITRE = '      TETRAEDRES INITIAUX A RE-TETRAEDRISER'
         WRITE(KTITRE(1:5),'(I5)') NBTEETPR
         CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTEETPR, NOTEETPR, NOTETR )


C        6) TETRAEDRISATION DES FACES DE L'ETOILE EN TORE
C           C-A-D DES NBSSET SOUS-ETOILES SANS OU AVEC DES POINTS
C           INTERNES AJOUTES
C        --------------------------------------------------------
C        VOLET1 VOLUME APRES TETRAEDRISATION DE L'ETOILE
         VOLET1 = 0D0
         CALL TETRTORE( QTEAME, KTITRE,
     %                  MXSOMM, NBSOMM0, NBSOMM, PTXYZD,
     %                  NPSOFR, VMOYEN,
     %                  MXTETR, N1TEVI, NOTETR, NUDTETR, N1TETS,
     %                  NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                  N1FEVI, MXFETO, NFETOI,
     %                  MXSSET, NBSSET, N1SSET, NBTEETPR,
     %                  MXTECR, NBTECR, NOTECR, VOLET1, IERR )


C        RE CALCUL DU VOLUME DES TETRAEDRES CREES pour VERIFICATION
         VOLET1 = 0D0
         QUAMIN = 2.0
         QUAMOY = 0.0
ccc         print*,'fapetore: LISTE des TETRAEDRES DE LA RE-TETRAEDRISATION'
         VOLET1 = 0D0
         NBT = 0
         DO K = 1, NBTECR

C           NUMERO NOTETR DU K-EME TETRAEDRE DE L'ETOILE
            NTE = NOTECR( K )

C           PROTECTION EN CAS DE PROBLEME LORS DE LA TETRAEDRISATION
C           DE L'ETOILE POUR REVENIR EN ARRIERE
            IF( NTE .GT. 0 ) THEN
               IF( NOTETR(1,NTE) .LE. 0 ) THEN
                  PRINT*,'fapetore: TETRAEDRE INACTIF NOTETR(',NTE,')=',
     %                   (NOTETR(L,NTE),L=1,8)
                  GOTO 9900
               ENDIF
               NBT = NBT + 1
               NOTECR( NBT ) = NTE
C              SON VOLUME ET SA QUALITE
               CALL QUATETD( PTXYZD(1,NOTETR(1,NTE)),
     %                       PTXYZD(1,NOTETR(2,NTE)),
     %                       PTXYZD(1,NOTETR(3,NTE)),
     %                       PTXYZD(1,NOTETR(4,NTE)),
     %                       ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
ccc            IF( VOLUTE .LE. 0D0 ) THEN
ccc               print*,'fapetore: TETRAEDRE FINAL  ',NTE,
ccc     %                ' St:',(NOTETR(L,NTE),L=1,4),
ccc     %                ' de VOLUME=',VOLUTE,' Qualite=',QUALTE
ccc            ENDIF
               VOLET1 = VOLET1 + ABS( VOLUTE )
               QUAMOY = QUAMOY + QUALTE
               IF( QUALTE .LT. QUAMIN ) THEN
                  QUAMIN = QUALTE
               ENDIF
            ENDIF

         ENDDO

C        VOLUME MOYEN DES NBT TETRAEDRES APRES LA RE-TETRAEDRISATION
         NBTECR = NBT
         VMOYEN = VOLET1 / NBT
         QUAMOY = QUAMOY / NBT
ccc         print *,'fapetore: FACE PERDUE',NFP,' NBTECR=',NBTECR,
ccc     %           ' VOLET1=',VOLET1,' VMOYEN1=',VMOYEN,
ccc     %           ' QualMOY1=',QUAMOY,' QualMIN1=',QUAMIN

         IF( IERR .NE. 0 ) THEN
            GOTO 9900
         ENDIF

C        CONTROLE DU NOMBRE NBTECR DE TETRAEDRES CREES DANS L'ETOILE
         IF( NBTECR .GE. MXTECR ) THEN
            print*,'???????????????????????????????????????????????????'
           print*,'PB fapetore: TROP DE TETRAEDRES CREES     ',NBTECR
           print*,'PB fapetore: NOMBRE DE TETRAEDRES INITIAUX',NBTEETPR
            print*,'BOUCLE INFINIE DE CREATION???'
            print*,'???????????????????????????????????????????????????'
            KTITRE='PB:       TETRAEDRES EN BOUCLE INFINIE DE CREATION?'
            GOTO 9900
         ENDIF


C        7) CONTROLE DU VOLUME DE L'ETOILE AVANT et APRES RETETRAEDRISATION
C        ------------------------------------------------------------------
         print *,'fapetore: NFP=',NFP,' VOLET0=',VOLET0,
     %         ' pour',NBTEETPR,' TETRAEDRES INITIAUX'
         print *,'fapetore: NFP=',NFP,' VOLET1=',VOLET1,
     %         ' pour',NBTECR,' TETRAEDRES FINAUX'

         IF( ABS( VOLET1 - VOLET0 ) .GT. 1D-5 * VOLET0 ) THEN

          print*,'?????????????????????????????????????????????????????'
       print*,'PROBLEME fapetore: VOLUMES TROP DIFFERENTS VOLET0=',
     %              VOLET0,'   VOLET1=',VOLET1
       print*,'PROBLEME fapetore: VOLUMES TROP DIFFERENTS VOLET1=',
     %              VOLET1     
         print*,'?????????????????????????????????????????????????????'
          KTITRE='PB:       TETRAEDRES DE VOLUME DIFFERENT DES INITIAUX 
     %V0=               V1=           '
            WRITE(KTITRE(58:71),'(G14.8)') VOLET0
            WRITE(KTITRE(76:89),'(G14.8)') VOLET1
            GOTO 9900

         ENDIF


C        8) COMPLETION DES CHAINAGES DES NBTECR TETRAEDRES CREES
C           DANS L'ETOILE
C        -------------------------------------------------------
         CALL MJOPTE( NBTECR, NOTECR, N1TETS, NOTETR, MXTETR,
     %                N1TEVI, PTXYZD, NBFANR )
         IF( NBFANR .LE. 0 ) GOTO 240

C        DESSOUS PEUT ETRE A SUPPRIMER
         LHPILE = NBTECR
 202     IF( LHPILE .GT. 0 ) THEN

C           LE HAUT DE LA PILE DES TETRAEDRES CREES DANS L'ETOILE
            NTE = ABS( NOTECR( LHPILE ) )

C           LE TETRAEDRE EST DEPILE
            LHPILE = LHPILE - 1

            IF( NTE .LE. 0 ) GOTO 202

            DO 220 J=1,4

C              FACE J DU TETRAEDRE NTE et NUMERO DE SES 3 SOMMETS
               DO K=1,3
                  NOSOTR(K) = NOTETR( NOSOFATE(K,J), NTE )
               ENDDO
C              TRI CROISSANT DES 3 SOMMETS DE LA FACE J DE NTE
               CALL TRI3NO( NOSOTR, NOSOTR )

C              LA FACE J DE NTE EST ELLE UNE FACE PERDUE DANS LEFACO?
               DO 205 M = 1, NBTRCF

C                 LA FACE M PERDUE DANS LEFACO
                  NTR = NOTRCF( M )
                  IF( NTR .LE. 0 ) GOTO 205

C                 LES 3 SOMMETS DE LA FACE NTR DE LEFACO EN TRI CROISSANT
                  CALL TRI3NO( LEFACO(1,NTR), NOSOTR2 )

                  IF( NOSOTR(1) .EQ. NOSOTR2(1) .AND.
     %                NOSOTR(2) .EQ. NOSOTR2(2) .AND.
     %                NOSOTR(3) .EQ. NOSOTR2(3) ) THEN
C                     LA FACE NTR=NOTRCF(K) EST LA FACE J DU TETRAEDRE NTE
C                     CE QUI DONNE UNE FACE NON PERDUE DANS LEFACO
                      LEFACO(11,NTR) = NTE
C                     MARQUAGE DE LA FACE DE NOTRCF POUR ACCELERER
                      NOTRCF(M) = -NTR
                      GOTO 208
                  ENDIF

 205           ENDDO

C              QUEL EST LE TETRAEDRE NTE1 OPPOSE A LA FACE J DE NTE
C              DE SOMMETS NS1 NS2 NS3 = NOSOTR(1:3)
C              POUR LE NO DES SOMMETS DES FACES
C              NOSOFATE/1,3,2, 2,3,4, 3,1,4, 4,1,2/
 208           IF( NOTETR(J+4,NTE) .GT. 0 ) GOTO 220

C              RECHERCHE DE CETTE FACE NOSOTR PARMI LES FACES
C              DES TETRAEDRES CREES DANS L'ETOILE
               DO 210 I=LHPILE,1,-1

                  NTE1 = ABS( NOTECR( I ) )
C                 REMARQUE: NTE1 EST DIFFERENT DE NTE
                  IF( NTE1 .LE. 0 ) GOTO 210

C                 LA FACE J DE NTE NOSOTR EST ELLE UNE FACE JJ DE NTE1
                  CALL NUFATRTE( NOSOTR, NOTETR(1,NTE1), JJ )

                  IF( JJ .GT. 0 ) THEN

C                    LA FACE J DE NTE EST LA FACE JJ DE NTE1
                     NOTETR( 4+J , NTE  ) = NTE1
                     NOTETR( 4+JJ, NTE1 ) = NTE
                     GOTO 220

                  ENDIF

 210           ENDDO

               PRINT *, 'fapetore: FACE ',NOSOTR, 'du TETRAEDRE',
     %                  (NOTETR(M,NTE),M=1,8),' SANS TETRAEDRE OPPOSE'

 220        ENDDO

C           RETOUR EN HAUT DE PILE
            GOTO 202
         ENDIF

C        VERIFICATION DE LA COMPLETION DES TETRAEDRES OPPOSES
 240     CALL MJOPTE( NBTECR, NOTECR, N1TETS, NOTETR, NUDTETR,
     %                N1TEVI, PTXYZD, NBFANR )
         PRINT*, 'fapetore: sortie mjopte avec NBFANR=',NBFANR
ccc         CALL VEOPTE( NBTECR, NOTECR, NOTETR, PTXYZD, NBFANR )
         IF( NBFANR .GT. 0 ) THEN
            IERR = NBFANR
            GOTO 9900
         ENDIF

C        DESTRUCTION DES NBTEETPR TETRAEDRES NOTEETPR DE PROTECTION
         DO K=NBTEETPR,1,-1
            NTE = ABS( NOTEETPR( K ) )
            NOTETR(1,NTE) = 0
            NOTETR(5,NTE) = N1TEVI
            N1TEVI = NTE
         ENDDO


C        9) SUPPRESSION DANS NOFAPE DES FACES PERDUES RETROUVEES OU A EVITER
C        -------------------------------------------------------------------
         DO 250 K=1,NBTRCF

C           LA FACE RETROUVEE EST EFFACEE DANS NOFAPE
            NTR = ABS( NOTRCF( K ) )

C           LE SIGNE + LUI EST RENDU
            NOTRCF( K ) = NTR

            DO J = 1, NBFAPE
               IF( NTR .EQ. ABS( NOFAPE(J) ) ) THEN
                  NOFAPE( J ) = -NTR
                  GOTO 250
               ENDIF
            ENDDO

 250     ENDDO


C        10) QUALITE et TRACE DES TETRAEDRES REMPLISSANT L'ETOILE OU ACTUELS
C           MISE A JOUR DU NO DE TETRAEDRE DES FACES DE LEFACO ET
C           DES FACES PERDUES
C        -------------------------------------------------------------------
         IF( IERR .EQ. 0 .AND. NBTECR .GT. 0 ) THEN

            QUAMOY = 0.0
            QUAMIN = 2.0
            DO K = 1, NBTECR

C              NUMERO NOTETR DU TETRAEDRE CREE
               NTE = ABS( NOTECR( K ) )

C              SON VOLUME ET SA QUALITE
               CALL QUATETD( PTXYZD(1,NOTETR(1,NTE)),
     %                       PTXYZD(1,NOTETR(2,NTE)),
     %                       PTXYZD(1,NOTETR(3,NTE)),
     %                       PTXYZD(1,NOTETR(4,NTE)),
     %                       ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
               QUAMOY = QUAMOY + QUALTE
               IF( QUALTE .LT. QUAMIN ) THEN
                  QUAMIN = QUALTE
               ENDIF

C              MISE A JOUR DU NO DE TETRAEDRE DES FACES DE LEFACO N1FASC
               DO J=1,4
                  CALL NULEFT( J, NTE, NOTETR, MXFACO, LEFACO,  NOFACO )
                  IF( NOFACO .GT. 0 ) THEN
C                    LE TETRAEDRE NTE CONTIENT LA FACE NOFACO DE LEFACO
                     LEFACO( 11, NOFACO ) = NTE
C                    UN NUMERO DE FACE LEFACO AUX 3 SOMMETS
                     DO N=1,3
                        N1FASC( LEFACO(N,NOFACO) ) = NOFACO
                     ENDDO
                  ENDIF
               ENDDO

ccc               print 109900,NFP,K,NTE,(NOTETR(KK,NTE),KK=1,8)
ccc            ENDDO
ccc109900       FORMAT('FIN fapetore: FACE PERDUE',I5,' TETRAEDRE CREE',I4,
ccc     %             ' NOTETR(',I7,')=',4I7,3X,4I7 )

            ENDDO
            QUAMOY = QUAMOY / NBTECR

            print *,'fapetore: FACE PERDUE',NFP,
     %              ' avec NBTRCF=',NBTRCF,
     %              ' a DONNE',NBTECR,' tetra de QUALITE MOYENNE',
     %              QUAMOY,' MINIMALE',QUAMIN,' avec',NBSOMM-NBSOMM0,
     %              ' SOMMETS AJOUTES pour',
     %              NBTEETPR,' tetra initiaux'

C           BILAN SUR LES FACES PERDUES DU CF RETROUVEES
            NBFPR = 0
            DO 965 K = 1, NBTRCF

C              NO FACE PERDUE DU CF DANS LEFACO
               NTR = ABS( NOTRCF( K ) )
               IF( NTR .LE. 0 ) GOTO 965
            
C              LE TETRAEDRE CONTENANT LA FACE
               NTE = LEFACO(11,NTR)
               IF( NTE .GT. 0 ) THEN
                  IF( NOTETR(1,NTE) .GT. 0 ) THEN

C                    LA FACE NTR DE LEFACO EST RETROUVEE
                     NBFPR = NBFPR + 1
                     PRINT *,'fapetore: FACE',NTR,
     %                '     RETROUVEE dans LEFACO=',
     %                    (LEFACO(L,NTR),L=1,11)
                     GOTO 965

                  ENDIF

               ENDIF

C              LA FACE NTR DE LEFACO N'EST PAS RETROUVEE
               LEFACO(11,NTR) = 0
               PRINT *,'fapetore: Face',NTR,
     %                 ' NON RETROUVEE dans LEFACO=',
     %                 (LEFACO(L,NTR),L=1,11)

 965        ENDDO

C           TRACE DES NBTECR TETRAEDRES CREES POUR CETTE FACE PERDUE NFP
            KTITRE='      TETRAEDRES CREES dans ETOILE FINALE QUALITE MO
     %Y=        QUALITE MIN=      '
            WRITE(KTITRE(1:5),    '(I5)') NBTECR
            WRITE(KTITRE(55:61),'(F7.5)') QUAMOY
            WRITE(KTITRE(75:81),'(F7.5)') QUAMIN
            CALL TRFETO6( KTITRE, PTXYZD,
     %                    NBTRCF, NOTRCF, LEFACO,
     %                    NBTECR, NOTECR, NOTETR )

         ENDIF


C        PASSAGE A LA FACE LEFACO PERDUE SUIVANTE
C        ----------------------------------------
         GOTO 1


C        ===========================================================
C        AFFICHAGE DU PROBLEME RENCONTRE & ABANDON DE LA FACE PERDUE
C        ===========================================================
 9900    PRINT*
         PRINT 19000
19000    FORMAT(150('!'))
         PRINT*,'fapetore: FINALEMENT ABANDON de la FACE PERDUE',NFP,
     %  ' avec NBTRCF=',NBTRCF,' et',
     %   NBSOMM-NBSOMM0,' SOMMETS AJOUTES  IERR=',IERR
         PRINT*,'fapetore: LEFACO(',NFLPER,')=',
     %          (LEFACO(L,NFLPER),L=1,11)
         PRINT 19000
         PRINT*
         NBABAN = NBABAN + 1

C        TRACE DES TETRAEDRES INITIAUX et CREES AVANT ABANDON
         TRACTE0 = TRACTE
         TRACTE = .TRUE.
         KTITRE='fapetore:           TETRAEDRES INITIAUX suite ABANDON'
         WRITE(KTITRE(10:16),'(I7)') NBTEETPR
         CALL TRFETO8( KTITRE,   PTXYZD,
     %                 NBTRCF,   NOTRCF,  LEFACO, NO0FAR,
     %                 NBTEETPR, NOTEETPR,
     %                 NBTECR,   NOTECR,  NOTETR )
C        ABANDON DU TRAITEMENT DE LA FACE PERDUE
         TRACTE = TRACTE0

         IF( NBTEETPR .GT. 0 ) THEN

C           DETRUIRE  LES NBTECR   TETRAEDRES NOTECR,
C           RESTAURER LES NBTEETPR TETRAEDRES NOTEETPR
C           METTRE A JOUR N1TEVI, NOTETR, N1TETS
            CALL RSTETO( PTXYZD, NBTEETPR, NOTEETPR, NBTECR,  NOTECR,
     %                   MXTETR, N1TEVI,   NOTETR,   NUDTETR, N1TETS,
     %                   NBFANR )

C        MISE A JOUR DE LEFACO POUR LES TETRAEDRES A DETRUIRE
         DO K=1,NBTECR
            NTE = NOTECR( K )
C           SUPPRIMER EVENTUELLEMENT LES 4 FACES DU TETRAEDRE NTE
C           DANS LE TABLEAU LEFACO
            CALL SUTELEFA( NTE, NOTETR, INFACO, MXFACO, LEFACO )
         ENDDO

C        MISE A JOUR DE LEFACO POUR LES TETRAEDRES A RESTAURER
         DO K=1,NBTEETPR
            NTE = NOTEETPR( K )
C           AJOUTER EVENTUELLEMENT LES 4 FACES DU TETRAEDRE NTE
C           DANS LE TABLEAU LEFACO
            CALL AJTELEFA( NTE, NOTETR, INFACO, MXFACO, LEFACO )
         ENDDO

C           PLUS AUCUN TETRAEDRE CREE APRES leur DESTRUCTION
            NBTECR = 0
            IF( NBFANR .GT. 0 ) THEN
               PRINT*, 'fapetore: PB RESTAURATION TETRA INITIAUX avec mj
     %opte => NBFANR=',NBFANR,' ?'
            ENDIF

         ENDIF

C        LE NOMBRE DE TETRAEDRES CREES EST REMIS A ZERO
         NBTECR = 0

         DO 9909 K = 1, NBTRCF

C           LA FACE PERDUE ET ABANDONNEE EST MARQUEE DANS NOFAPE
            NTR = ABS( NOTRCF( K ) )

C           LE SIGNE + LUI EST RENDU
            NOTRCF( K ) = NTR

            DO J = 1, NBFAPE
               IF( NTR .EQ. ABS( NOFAPE(J) ) ) THEN
                  NOFAPE( J ) = -NTR
                  GOTO 9909
               ENDIF
            ENDDO

 9909    ENDDO

C        RECUPERATION DES POINTS AJOUTES POUR RETROUVER CETTE FACE PERDUE
         NBSOMM = NBSOMM0
         GOTO 1

      ENDIF


C     SORTIE FINALE de fapetore
C     =========================
C     MISE A JOUR DU TABLEAU N1TETS
      CALL AZEROI( MXSOMM, N1TETS )
      DO NTE = 1, MXTETR
         IF( NOTETR(1,NTE) .GT. 0 ) THEN
            DO K=1,4
               N1TETS( NOTETR(K,NTE) ) = NTE
            ENDDO
         ENDIF
      ENDDO

      PRINT*
      PRINT*,'fapetore: SORTIE avec NBSOMM=',NBSOMM,' et ',NBABAN,
     %' ABANDONS de FACES PERDUES sur',NBFAPE,
     %' FACES INITIALEMENT PERDUES'
      PRINT*

      TRACTE = .FALSE.
      RETURN
      END
