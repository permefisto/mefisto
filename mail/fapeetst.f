      SUBROUTINE FAPEETST( QUAMINEX, QTEAME,
     %                     NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                     MXFAPE, NBFAPE, NOFAPE,
     %                     INFACO, MXFACO, LEFACO, N1FASC,
     %                     MXTETR, N1TEVI, N1TETS,NOTETR, NUDTETR,
     %                     MXTRCF, NOTRCF, NOSTCF,
     %                     MXARCF, N1ARCF, NOARCF,
     %                     MXFETO, NFETOI, NOSTIS,
     %                     IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECUPERATION DE FACES PERDUES DE LA FRONTIERE
C -----    PAR L'ETOILE DE TOUS LES SOMMETS DU CF ET SOMMETS ISOLES

C ENTREES:
C --------
C QUAMINEX: QUALITE MINIMALE EXIGEE POUR CREER UN TETRAEDRE
C QTEAME : QUALITE AU DESSOUS DE LAQUELLE UN TETRAEDRE N'EST PAS CREE ou
C          QUALITE DES TETRAEDRES AU DESSOUS DE LAQUELLE UNE
C          AMELIORATION DE LA QUALITE DES TETRAEDRES EST DEMANDEE
C NOFOTI : NUMERO DE LA FONCTION UTILISATEUR 'TAILLE_IDEALE', 0 SINON
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

C INFACO : = 0 PAS DE TABLEAU LEFACO NI DE CONSERVATION DES
C              FACES FRONTIERE ( NO DE VOLUME CONNU PAR NVOLTE )
C          = 1 EXISTENCE DU TABLEAU LEFACO ET CONSERVATION DES
C              FACES DE LA FRONTIERE
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
C IERR    : 0 SI PAS D'ERREUR
C          1 ETOILE AVEC MOINS DE 4 FACES
C          2 SATURATION D'UN TABLEAU
C          3 3 SOMMETS NON DANS UN TETRAEDRE
C          9 VOLUME DE L'ETOILE AVANT et APRES DIFFERENTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & VEULETTES SUR MER   Novembre 2016
C2345X7..............................................................012
      PARAMETER        (MXTECF=400000,
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

      INTEGER           NOTECF(MXTECF),   NOTECFPR(MXTECF),
     %                  NFETOI(5,MXFETO), LAPILE(2,MXPILE),
     %                  NO0FAR(3,MXETOI), N1SSET(MXSSET)

      CHARACTER*200     KTITRE
      DOUBLE PRECISION  COS2FA, COSMIN, VMOYEN,
     %                  VOLET0, VOLET1, VOLUTE,
     %                  ARMIN, ARMAX, SURFTR(4)
      REAL              QUALTE
      INTEGER           NOSOAR(2,6), NOSOAROP(2,6), NOSOFA(3,4),
     %                  NOSOTR(3), NOSOTR2(3),
     %                  NFA(3), NOFAARTE(2,6)
      DATA              NOSOAR   / 1,2,  2,3,  3,1,  1,4,  2,4,  3,4 /
      DATA              NOSOAROP / 2,3,  1,4,  2,4,  2,3,  3,1,  1,2 /
      DATA              NOSOFA   / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE
C     NUMERO DES 2 FACES DU TETRAEDRE AYANT UNE ARETE COMMUNE K
      DATA              NOFAARTE / 1,4,  1,2,  1,3,  3,4,  2,4,  2,3 /

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
C     METHODE: FORMER 2 DEMI-ETOILES AU DESSUS ET AU DESSOUS
C              DES TETRAEDRES DE SOMMETS LES SOMMETS ET ISOLES
C              DU CF PUIS LES TETRAEDRISER
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

C        PROTECTION DU NOMBRE ACTUEL DE SOMMETS PTXYZD
         IERR    = 0
         NBSOMM0 = NBSOMM
         NBSSET = 0
         NBTECFPR = 0
         NBTECF0  = 0
         NBTECF = 0
         NBTECR = 0
         NBPTIN = 0
         NBTRCF0= 0
         NBTRCF = 0
         NBSTIS = 0
         NOPASS = 0
         PRINT *


C        2.0) FORMATION DU PREMIER CONTOUR FERME FORME DES
C             TRIANGLES LEFACO PERDUS ADJACENTS A PARTIR DE NFLPER
C        ---------------------------------------------------------
         CALL VDR1CF( -1.0,   NFLPER, NBSOMM, PTXYZD, N1TETS, NOTETR,
     %                NBFAPE, NOFAPE, MXFACO, LEFACO, N1FASC, NO0FAR,
     %                MXETOI, NAETOI,
     %                MXTRCF, NBTRCF, NOTRCF,
     %                MXARCF, NBCF,   N1ARCF, NOARCF,
     %                MXTRCF, NBSTCF, NOSTCF,
     %                MXSOMM, NBSTIS, NOSTIS,
     %                IERR )
         NBTRCF0 = NBTRCF

         PRINT *
         PRINT *, 'fapeetst: FACE PERDUE',NFP,' /',NBFAPE,
     %            ' LEFACO(',NFLPER,')=',(LEFACO(K,NFLPER),K=1,11),
     %            ' NBTRCF=',NBTRCF,' IERR=',IERR
         PRINT 10002
10002    FORMAT(240('='))

         IF( NBTRCF .LE. 0 .OR. IERR .NE. 0 ) GOTO 9900

C        2.1) LISTER LES TETRAEDRES DES SOMMETS DU CF OU ISOLES
C        ------------------------------------------------------
         NBTECF = 0
         DO NN = 1, NBSTCF+NBSTIS

C           LE SOMMET DES NBTRCF TRIANGLES DU CF
            IF( NN .LE. NBSTCF ) THEN
               NS = NOSTCF( NN )
            ELSE
               NS = NOSTIS( NN-NBSTCF )
            ENDIF

            CALL TETR1S( NS,    N1TETS,        NOTETR,
     %                   NBTNS, MXTECF-NBTECF, NOTECF(NBTECF+1), IERR )
            NBTECF = NBTECF + NBTNS
C           NBTECF : EN ENTREE NOMBRE DE TETRAEDRES DEJA RANGES DANS NOTECF
C                    EN SORTIE AJOUT DES TETRAEDRES AYANT NS COMME SOMMET
C                    = 0 SI SATURATION DU TABLEAU NOTECF
C                    =-1 SI UN SOMMET NS N'EST PAS UN SOMMET DE N1TETS(NS)

C           UNICITE DES NO DE TETRAEDRES DU TABLEAU NOTECF
            CALL UNITABL( NOTECF, NBTECF )

            IF( IERR .NE. 0 ) THEN
               print*,'fapeetst: SATURATION du TABLEAU NOTECF NBTECF=',
     %                 NBTECF,' MXTEXF=',MXTECF,' IErr=',IERR
               GOTO 9900
            ENDIF

         ENDDO

C        VERIFICATION DES TETRAEDRES OPPOSES AUX FACES DES
C        TETRAEDRES DE SOMMET DU CF OU ISOLES
         CALL VEOPTE( NBTECF, NOTECF, NOTETR, PTXYZD, NBFANR )
         IF( NBFANR .GT. 0 ) THEN
            PRINT *,'fapeetst: 2.1) NOMBRE INCORRECT DE TETRAEDRES OPPOS
     %ES=',NBFANR
            GOTO 9900
         ENDIF

C        VERIFICATION QUE TOUTE ARETE DU CF APPARTIENT A 2 ET
C        SEULEMENT 2 FACES SIMPLES DES TETRAEDRES DE L'ETOILE ACTUELLE
C        -------------------------------------------------------------
C        CONSTRUCTION NFETOI DES TRIANGLES FACES SIMPLES DES TETRAEDRES
C        ENVELOPPANTS. LES FACES VUES 2 FOIS SONT ELIMINEES
         CALL CRFETOI1( NBTECF, NOTECF, NOTETR,
     %                  MXFETO, N1FEOC, N1FEVI, NFETOI )

C        TRACE DES FACES SIMPLES DE L'ETOILE DU CF NFETOI VERSION 1
         CALL TRFETO2( PTXYZD,  NOTETR, N1FEOC, NFETOI,
     %                 NBTRCF0, NOTRCF, LEFACO, NO0FAR )

         IF( N1FEVI .EQ. -1 ) GOTO 9900
         IF( N1FEOC .LE.  0 ) GOTO 9900


C        2.2) AJOUT DES TETRAEDRES OPPOSES COMMUNS A 2 FACES SIMPLES
C             DE L'ETOILE SI AUCUNE DE SES ARETES EST ARETE DU CF
C        -----------------------------------------------------------
 10      NBT0 = NBTECF
         NF1  = N1FEOC
 12      IF( NF1 .GT. 0 ) THEN

C           NUMERO NOTETR DU TETRAEDRE DE LA FACE SIMPLE
            NTE = NFETOI(1,NF1)
C           NUMERO DE LA FACE SIMPLE DANS LE TETRAEDRE NTE
            I = ABS( NFETOI(2,NF1) )
C           TETRAEDRE OPPOSE PAR LA FACE I AU TETRAEDRE NTE
            NTEOP = NOTETR(4+I,NTE)
            IF( NTEOP .LE. 0 ) GOTO 18
            IF( NOTETR(1,NTEOP) .EQ. 0 ) GOTO 18

C           EXISTE T IL UNE AUTRE FACE DE L'ETOILE AYANT CE MEME
C           TETRAEDRE OPPOSE?
            NF2 = NFETOI(5,NF1)
 15         IF( NF2 .GT. 0 ) THEN
               NTE2 = NFETOI(1,NF2)
               I2   = ABS( NFETOI(2,NF2) )
C              TETRAEDRE OPPOSE PAR LA FACE I2 AU TETRAEDRE NTE2
               NTEOP2= NOTETR(4+I2,NTE2)

               IF( NTEOP2 .EQ. NTEOP ) THEN

                  DO J=1,6
C                    LES 2 SOMMETS DE L'ARETE J DU TETRAEDRE NTEOP
                     NS1 = NOTETR( NOSOAR(1,J), NTEOP )
                     NS2 = NOTETR( NOSOAR(2,J), NTEOP )
C                    NS1-NS2 EST ELLE UNE ARETE DU CF?
                     CALL ARARCF( NS1, NS2, NBCF, N1ARCF, NOARCF,
     %                            NFLCF, NAVANT, NA1 )
                     IF( NA1 .NE. 0 ) GOTO 17
C                     OUI => LE TETRAEDRE NTEOP N'EST PAS AJOUTE
C                     POUR EVITER LE TOUR AUTOUR D'UNE ARETE CF
                  ENDDO

C                 AJOUT DU TETRAEDRE NTEOP A L'ETOILE
                  IF( NBTECF .GE. MXTECF ) THEN
                 print*,'fapeetst: SATURATION NOTECF. AUGMENTER MXTECF='
     %                  ,MXTECF
C                    ABANDON
                     GOTO 9900
                  ENDIF
                  NBTECF = NBTECF + 1
                  NOTECF( NBTECF ) = NTEOP
                  print*,'fapeetst: 2.2) ajout',NBTECF,
     %                   ' du tetraedre NOTETR(',nteop,')=',
     %                   (notetr(kk,nteop),kk=1,8)

                  DO J=1,4
C                    LA FACE J DU TETRAEDRE NTEOP EST
C                    SOIT RETIREE, SOIT AJOUTEE EN FIN DE CHAINAGE
                     CALL AJFAET1( NTEOP, J, NOTETR,
     %                             N1FEOC, N1FEVI, NFETOI, NFS )
                     IF( N1FEVI .EQ. -1 ) GOTO 9900
                  ENDDO
                  GOTO 10

               ENDIF

C              LA FACE SUIVANTE DE L'ETOILE DE MEME TETRAEDRE A RETROUVER
 17            NF2 = NFETOI(5,NF2)
               GOTO 15
            ENDIF

C           LA FACE SUIVANTE DE L'ETOILE A ANALYSER
 18         NF1 = NFETOI(5,NF1)
            GOTO 12

         ENDIF

C        TRACE DES FACES SIMPLES DE L'ETOILE DU CF NFETOI VERSION 1
         CALL TRFETO2( PTXYZD,  NOTETR, N1FEOC, NFETOI,
     %                 NBTRCF0, NOTRCF, LEFACO, NO0FAR )

         IF( NBT0 .LT. NBTECF ) THEN
C           AU MOINS UN TETRAEDRE A ETE AJOUTE => NOUVEL ESSAI
            GOTO 10
         ENDIF


C        2.3) RETRAIT DES TETRAEDRES DE L'ETOILE AYANT AU MOINS
C             3 FACES SIMPLES DE L'ETOILE ET AU PLUS UN SOMMET DU CF
C          ou 2 FACES SIMPLES DE L'ETOILE ET ZERO SOMMET DU CF
C        -----------------------------------------------------------
 30      NBTESU = 0
         DO 35 N = 1, NBTECF

C           LE TETRAEDRE N DE L'ETOILE
            NTE = NOTECF( N )
            IF( NTE .LE. 0 ) GOTO 35

C           NOMBRE DE FACES DE L'ETOILE DE CE TETRAEDRE
C           PARCOURS DES FACES SIMPLES DE L'ETOILE EN VERSION 1
            NBFEDT = 0
            NF1    = N1FEOC
 31        IF( NF1 .GT. 0 ) THEN
C              NUMERO NOTETR DU TETRAEDRE DE LA FACE SIMPLE
               IF( NTE .EQ. ABS( NFETOI(1,NF1) ) ) NBFEDT = NBFEDT + 1
C              LA FACE SUIVANTE DE L'ETOILE
               NF1 = NFETOI(5,NF1)
               GOTO 31
            ENDIF

C           NOMBRE DE SOMMETS DE NTE SUR LE CF
            NBSTFR=0
            DO 33 K=1,4
               NS = NOTETR(K,NTE)
               DO L=1,NBSTCF
                  IF( NS .EQ. NOSTCF(L) ) THEN
                     NBSTFR = NBSTFR + 1
                     GOTO 33
                  ENDIF
               ENDDO
 33         ENDDO

ccc            PRINT *,'fapeetst: 2.3) TETRAEDRE',NTE,' St:',
ccc     %       (NOTETR(L,NTE),L=1,8),' avec',NBFEDT,'FacesEtoile+',
ccc     %        NBSTFR,'St du CF'

            IF( NBFEDT .GE. 3 .AND. NBSTFR .LE. 1  .OR.
     %          NBFEDT .EQ. 2 .AND. NBSTFR .LE. 0 ) THEN

C               LE TETRAEDRE NTE EST SUPPRIME DE L'ETOILE
                NBTESU = NBTESU + 1
                PRINT *,'fapeetst: 2.3) TETRAEDRE',NTE,' SUPPRIME',
     %                  ' St:', (NOTETR(L,NTE),L=1,8),' avec',NBFEDT,
     %                  ' FacesEtoile+',NBSTFR,'St du CF'
                NOTECF( N ) = -NTE

C               CONSTRUCTION NFETOI VERSION 1 DES TRIANGLES FACES SIMPLES
C               DES TETRAEDRES. LES FACES VUES 2 FOIS SONT ELIMINEES
                CALL CRFETOI1( NBTECF, NOTECF, NOTETR,
     %                         MXFETO, N1FEOC, N1FEVI, NFETOI )

C               TRACE DES FACES SIMPLES DE L'ETOILE FINALE
C               DU CF NFETOI VERSION 1
                CALL TRFETO2( PTXYZD, NOTETR, N1FEOC, NFETOI,
     %                        NBTRCF0, NOTRCF, LEFACO, NO0FAR )

             ENDIF

 35       ENDDO

C        COMPRESSION DU TABLEAU NOTECF EN ELIMINANT LES NUMEROS <=0
         NBTECF0 = NBTECF
         CALL COMPENTP( NBTECF0, NOTECF,  NBTECF )
         IF( NBTECF .LE. 0 ) GOTO 9900

C        CONTINUER LA SUPPRESSION DES TETRAEDRES EXTERIEURS?
         IF( NBTESU .GT. 0 ) GOTO 30


C        2.4) SUPPRESSION D'UN TETRAEDRE D'ANGLE MAXIMAL AVEC LES
C             TRIANGLES DU CF ET DES TETRAEDRES OPPOSES AUX 4 FACES
C             POUR PERMETTRE LA FORMATION DE 2 DEMI-ETOILES FERMEES
C             APRES L'AJOUT DES TRIANGLES DU CF
C        ------------------------------------------------------------------
C        NOARCF : 1:NUMERO DU SOMMET DE L'ARETE DE LA LIGNE DU CONTOUR FERME
C                 2:NUMERO DANS NOARCF DE L'ARETE SUIVANTE DU CF
C                 3:NUMERO DANS LEFACO DU TRIANGLE ADJACENT OPPOSE A L'ARETE
         NBTESU = 0
         DO NCF = 1, NBCF

C           ARETE INITIALE DU CF NCF
            NA0  = 0
            NA00 = N1ARCF( NCF )
            NA1  = NA00

C           ARETE NA1-NA2 DU CF
 72         NA2 = NOARCF(2,NA1)
C           LES 2 SOMMETS DE L'ARETE NA1 DU CF
            NS1 = NOARCF(1,NA1)
            NS2 = NOARCF(1,NA2)
            NS3 = 0

C           RETROUVER LE TRIANGLE NTR DU CF D'ARETE NS1-NS2
            DO K = 1, NBTRCF
               NTR = NOTRCF(K)
               IF( NTR .GT. 0 ) THEN
C                 FACE DE LEFACO
                  DO L=1,3
                     NAS1 = LEFACO(L,NTR)
                     IF( L .EQ. 3 ) THEN
                        LL = 1
                     ELSE
                        LL = L+1
                     ENDIF
                     NAS2 = LEFACO(LL,NTR)
                     IF( NAS1 .EQ. NS1 .AND. NAS2 .EQ. NS2 .OR.
     %                   NAS2 .EQ. NS1 .AND. NAS1 .EQ. NS2 ) THEN
C                          L'ARETE L DU TRIANGLE NTR EST L'ARETE NA1 DU CF
C                          NS3 NUMERO DU SOMMET DU TRIANGLE NTR NON NS1 ou NS2
                        IF( LL .EQ. 3 ) THEN
                           LLL = 1
                        ELSE
                           LLL = LL+1
                        ENDIF
C                       NS3 NUMERO DU SOMMET DU TRIANGLE NTR NON NS1 et NON NS2
                        NS3 = LEFACO(LLL,NTR)
                        GOTO 74
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO

C           RECHERCHE DES TETRAEDRES NOTECF DE L'ETOILE D'ARETE NS1-NS2
 74         NBTEAR = 0
            COSMIN = 2D0
            NFAMIN = 0
            NTEMIN = 0
            NMIN   = 0

            DO 78 NN = 1, NBTECF

C              LE TETRAEDRE NN DE NOTECF L'ETOILE
               NTE = NOTECF( NN )
               IF( NTE .LE. 0 ) GOTO 78

               DO NA=1,6

C                 LES 2 SOMMETS DE L'ARETE NA DU TETRAEDRE NTE
                  NAS1 = NOTETR( NOSOAR(1,NA), NTE )
                  NAS2 = NOTETR( NOSOAR(2,NA), NTE )

                  IF( NAS1 .EQ. NS1 .AND. NAS2 .EQ. NS2 .OR.
     %                NAS2 .EQ. NS1 .AND. NAS1 .EQ. NS2 ) THEN

C                    L'ARETE NA DU TETRAEDRE NTE EST L'ARETE NS1-NS2 DU CF
C                    ARETE L DE LA FACE LEFACO NTR ou TRIANGLE DU CF
C                    NOFAARTE(1:2,NA) NO DES 2 FACES DE NTE AYANT CETTE ARETE NA
                     NBTEAR = NBTEAR + 1
                     NFA(1) = NOFAARTE(1,NA)
                     NFA(2) = NOFAARTE(2,NA)

                     DO MM=1,2

C                       NS4 NUMERO DU SOMMET DE LA FACE NFA(MM) DE NTE
C                       DIFFERENT DE NAS1 ET NAS2
                        DO KK=1,3
                           NS4 = NOTETR( NOSOFA(KK,NFA(MM)), NTE )
                           IF( NS4 .NE. NAS1 .AND. NS4 .NE. NAS2)THEN
                              GOTO 76
                           ENDIF
                        ENDDO

C                       COSINUS DE L'ANGLE ENTRE LES NORMALES AUX 2 FACES
C                       CELLES DU CF NAS1-NAS2-NS3 et du TETRAEDRE
C                       NAS1-NAS2-NS4 ORIENTEES SELON LE SENS NAS1-NAS2
 76                     CALL COS2TD( PTXYZD(1,NAS1),
     %                               PTXYZD(1,NAS2),
     %                               PTXYZD(1,NS3),
     %                               PTXYZD(1,NAS1),
     %                               PTXYZD(1,NAS2),
     %                               PTXYZD(1,NS4),
     %                               COS2FA, IERR1, IERR2 )
                        IF( IERR1.EQ.0 .AND. IERR2.EQ.0 .AND.
     %                      COS2FA .LT. COSMIN ) THEN
                           COSMIN = COS2FA
                           NMIN   = NN
C                          NUMERO NOTETR DU TETRAEDRE
                           NTEMIN = NTE
C                          NUMERO LOCAL DE SA FACE
                           NFAMIN = NFA(MM)
C                          AUTRE SOMMET QUE NS1-NS2
                           NS4MIN = NS4
                        ENDIF

                     ENDDO
                     GOTO 78


                  ENDIF

               ENDDO

 78         ENDDO

            IF( NMIN .NE. 0 ) THEN

C              IL EXISTE UN TETRAEDRE NTEMIN DE FACE NFAMIN DE COSMIN
C              IL EST SUPPRIME SAUF S'IL A PLUS DE 3 SOMMETS SUR LE CF
               NBSCF = 0
               DO 79 M=1,4

c                 NUMERO DU SOMMET M DE NTEMIN
                  NS = NOTETR(M,NTEMIN)

                  DO NN = 1, NBSTCF
C                    NUMERO DU SOMMET NN DU CF
                     NSC = NOSTCF( NN )
                     IF( NS .EQ. NSC ) THEN
                        NBSCF = NBSCF + 1
                        GOTO 79
                     ENDIF
                  ENDDO

                  DO NN = 1, NBSTIS
C                    NUMERO DU SOMMET NN ISOLE
                     NSC = NOSTIS( NN )
                     IF( NS .EQ. NSC ) THEN
                        NBSCF = NBSCF + 1
                        GOTO 79
                     ENDIF
                  ENDDO

 79            ENDDO

               IF( NBSCF .LE. 2 ) THEN

C                 NTEMIN EST SUPPRIME DE L'ETOILE NOTECF
                  NBTESU = NBTESU + 1
                  NOTECF( NMIN ) = - NOTECF( NMIN )

                  PRINT*,'2.4) RETRAIT NOTECF du TETRAEDRE MIN',
     %                    NTEMIN,' St:',(NOTETR(MM,NTEMIN),MM=1,8),
     %                   ' NB St CF=',NBSCF,' COSMIN=',COSMIN

               ELSE

                  PRINT*,'2.4) GARDE   NOTECF du TETRAEDRE MIN',
     %                    NTEMIN,' St:',(NOTETR(MM,NTEMIN),MM=1,8),
     %                   ' NB St CF=',NBSCF,' COSMIN=',COSMIN

               ENDIF

C              LE TETRAEDRE OPPOSE AUX 4 FACES DU TETRAEDRE NTEMIN
C              INTERSECTE T IL UNE FACE DU CF?
               DO 82 M=1,4
                  NTEOP = NOTETR( 4+M, NTEMIN )
                  IF( NTEOP .GT. 0 ) THEN
C                    RECHERCHE DANS NOTECF DE NTEOP
                     DO KK=1,NBTECF
                        NTE = ABS( NOTECF(KK) )
                        IF( NTE .EQ. NTEOP ) THEN
                           IF( NOTECF(KK) .GT. 0 ) THEN

C                             CE TETRAEDRE INTERSECTE IL UNE FACE DU CF?
                              DO MM=1,NBTRCF
                                 NTR = NOTRCF( MM )
                                 IF( NTR .GT. 0 ) THEN
C                                   FACE DE LEFACO NTR
C                                   INTERSECTION TRIANGLE-TETRAEDRE?
                            CALL INTRITET( PTXYZD( 1, LEFACO(1,NTR) ),
     %                                     PTXYZD( 1, LEFACO(2,NTR) ),
     %                                     PTXYZD( 1, LEFACO(3,NTR) ),
     %                                     PTXYZD( 1, NOTETR(1,NTEOP) ),
     %                                     PTXYZD( 1, NOTETR(2,NTEOP) ),
     %                                     PTXYZD( 1, NOTETR(3,NTEOP) ),
     %                                     PTXYZD( 1, NOTETR(4,NTEOP) ),
     %                                     LINTER )
                                    IF( LINTER  .NE. 0 ) THEN
C                                       OUI: UN POINT D'INTERSECTION
C                                       =>ABANDON DE CE TETRAEDRE
                                       GOTO 82
                                    ENDIF
                                 ENDIF
                              ENDDO

C                             NTEOP N'INTERSECTE PAS UNE DES FACES DU CF
C                             NOMBRE DE SOMMETS DE NTEOP SUR LE CF
                              NONCONS = 0
                              LM0 = -1
                              NBSTFR=0
                              DO 81 KM=1,4
                                 NS = NOTETR(KM,NTE)
                                 DO LM=1,NBSTCF
                                    IF( NS .EQ. NOSTCF(LM) ) THEN
C                                      LE SOMMET KM DE NTE EST
C                                      LE SOMMET LM DU CF
                                       IF( LM .EQ. LM0+1  .OR.
     %                                     LM .EQ. LM0-1  .OR.
     %                                   ( LM .EQ. NBSTCF .AND.
     %                                     LM0.EQ. 1) ) THEN
C                                         LES 2 SOMMETS DU CF SONT
C                                         CONSECUTIFS
                                          NONCONS = 0
                                       ELSE
C                                         LES 2 SOMMETS DU CF NE SONT PAS
C                                         CONSECUTIFS
                                          NONCONS = 1
                                       ENDIF
                                       NBSTFR = NBSTFR + 1
                                       LM0 = LM
                                       GOTO 81
                                    ENDIF
                                 ENDDO
 81                           ENDDO

C                             NTEOP EST SUPPRIME DE L'ETOILE NOTECF
C                             S'IL A SOIT UN ou ZERO SOMMET SUR LE CF
C                                    SOIT 2 SOMMETS CONSECUTIFS DU CF
                              IF( NBSTFR .LE. 1 .OR.
     %                          ( NBSTFR .EQ.2 .AND. NONCONS.EQ.0)) THEN
                                 NBTESU = NBTESU + 1
                                 NOTECF(KK) = - NOTECF(KK)
                              PRINT*,'2.4) retrait NOTECF du tetraedre',
     %                            NTEOP,' st:',(NOTETR(MM,NTEOP),MM=1,8)
                                 GOTO 82
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
 82            ENDDO

            ENDIF

C           PASSAGE A L'ARETE SUIVANTE DU CF
            IF( NA2 .NE. NA00 ) THEN
               NA0 = NA1
               NA1 = NA2
               GOTO 72
            ENDIF

         ENDDO

C        COMPRESSION DU TABLEAU NOTECF EN ELIMINANT LES NUMEROS <=0
         NBTECF0 = NBTECF
         CALL COMPENTP( NBTECF0, NOTECF,  NBTECF )

         PRINT*,'2.4) RETRAIT pour ANGLE MAXIMAL de',NBTESU,
     %          ' TETRAEDRES sur',NBTECF0,' TETRAEDRES de l''ETOILE'


C        2.5) VERIFICATION QUE TOUTE ARETE DU CF APPARTIENT A 2 ET SEULEMENT
C             2 FACES SIMPLES DES TETRAEDRES DE L'ETOILE ACTUELLE ET QUE
C             LEUR ORIENTATION EST BONNE (DANS UN  SENS ET OPPOSE)
C        -------------------------------------------------------------------
C        CONSTRUCTION NFETOI DES TRIANGLES FACES SIMPLES DES
C        NBTECF TETRAEDRES DE L'ETOILE
         CALL CRFETOI1( NBTECF, NOTECF, NOTETR,
     %                  MXFETO, N1FEOC, N1FEVI, NFETOI )
         IF( N1FEVI .EQ. -1 ) GOTO 9900

C        TRACE DES FACES SIMPLES DE L'ETOILE DU CF NFETOI VERSION 1
         CALL TRFETO2( PTXYZD,  NOTETR, N1FEOC, NFETOI,
     %                 NBTRCF0, NOTRCF, LEFACO, NO0FAR )

C        2.6) DECOUPAGE DES FACES DE L'ETOILE EN 2 DEMI-ETOILES REPARTIES
C             DES 2 COTES DES FACES PERDUES QUI CONSTITUENT LE CF
C             AJOUT DES NBTRCF TRIANGLES DU CF POUR FERMER LA DEMI-ETOILE 2
C             AJOUT DES NBTRCF TRIANGLES DU CF POUR FERMER LA DEMI-ETOILE 1
C        ------------------------------------------------------------------
C        RECHERCHE D'UNE ARETE DE DEPART DU CF COMME ARETE DE
C        L'UN DES TRIANGLES DE L'ETOILE
         K = 0
         NF1 = N1FEOC
 151     IF( NF1 .GT. 0 ) THEN

C              NUMERO NOTETR DU TETRAEDRE DE LA FACE SIMPLE
               NTE = NFETOI(1,NF1)
C              NUMERO DE LA FACE SIMPLE DANS LE TETRAEDRE NTE
               I = ABS( NFETOI(2,NF1) )

C              LES 3 SOMMETS DE LA FACE SIMPLE NF1 
C              AVEC UNE NORMALE VERS L'INTERIEUR DE L'ETOILE
               NOSOTR(1) = NOTETR( NOSOFA(1,I), NTE )
               NOSOTR(2) = NOTETR( NOSOFA(3,I), NTE )
               NOSOTR(3) = NOTETR( NOSOFA(2,I), NTE )

C              LES 3 ARETES DE LA FACE TRIANGULAIRE DE L'ETOILE
C              LE NUMERO DE L'ARETE 1:1-2, 2:2-3 3:3-1 DE TRAVERSEE
               DO K=1,3

C                 NS1 NS2 SOMMETS DE L'ARETE K DE LA FACE I DU TETRAEDRE NTE
                  NS1 = NOSOTR(K)
                  IF( K .EQ. 3 ) THEN
                     NS2 = 1
                  ELSE
                     NS2 = K+1
                  ENDIF
                  NS2 = NOSOTR(NS2)

C                 NS1-NS2 EST ELLE UNE ARETE DU CF?
                  CALL ARARCF( NS1, NS2, NBCF, N1ARCF, NOARCF,
     %                         NFLCF, NAVANT, NA1 )
                  IF( NA1 .NE. 0 ) THEN
C                    L'ARETE NA1 DE NOARCF EST L'ARETE K DE LA FACE I
C                    DU TETRAEDRE NTE
C                    LE NUMERO DES 2 SOMMETS DONNANT LE SENS NS1ACF->NS2ACF
C                    DE L'ARETE D'UN TRIANGLE BIEN ORIENTE DE L'ETOILE
C                    DES FACES DES TETRAEDRES
                     NS1ACF = NS1
                     NS2ACF = NS2
                     GOTO 152
                  ENDIF

               ENDDO

C              LA FACE SUIVANTE DE L'ETOILE
               NF1 = NFETOI(5,NF1)
               GOTO 151

         ENDIF

C           L'ARETE K DE LA FACE NF1 DE NFETOI CONTIENT L'ARETE NA1 DU CF
C           CONSTRUCTION DE LA SECONDE DEMI-ETOILE DE DEBUT N2FEOC
C           EN RETIRANT DE L'ETOILE1 TOUS LES TRIANGLES DE L'ETOILE
C           LIES PAR UNE ARETE ADJACENTE A CETTE FACE NF1
C           -------------------------------------------------------------
 152        N2FEOC = 0
C           LA PILE DES ARETES DES FACES DE L'ETOILE A TRAITER
            LHPILE = 1
C           LA FACE NFETOI INITIALE
            LAPILE( 1, 1 ) = NF1
C           LE NUMERO DE L'ARETE 1:1-2, 2:2-3 3:3-1 DE TRAVERSEE
            LAPILE( 2, 1 ) = K

C           PARCOURS DES FACES ADJACENTES A NF2 SANS TRAVERSER LE CF
C           TANT QUE LA PILE DES FACES EST NON VIDE, TRAITER LE HAUT DE PILE
 153        IF( LHPILE .GT. 0 ) THEN

C              LA FACE NFETOI A TRAITER
               NF2 = LAPILE( 1, LHPILE )
               IF( NF2 .LE. 0 ) THEN
                 PRINT*,'fapeetst: PB 153 NF2=LAPILE(1,',LHPILE,')=',NF2
C                 LA FACE EST DEPILEE
                  LHPILE = LHPILE - 1
                  GOTO 153
               ENDIF
C              LE NUMERO DE L'ARETE 1:1-2, 2:2-3 3:3-1 DE TRAVERSEE
               NC  = LAPILE( 2, LHPILE )

C              LA FACE EST DEPILEE
               LHPILE = LHPILE - 1

C              NF2 EST RETIREE DE L'ETOILE1
               IF( NF2 .EQ. N1FEOC ) THEN

C                 LA FACE NF2 EST LA PREMIERE DE L'ETOILE RESTANTE
                  N1FEOC = NFETOI(5,N1FEOC)
                  IF( N1FEOC .LE. 0 ) THEN

C                    LES FACES DE L'ETOILE ONT TOUTES ETE RETROUVEES
C                    DANS LA PREMIERE DEMI-ETOILE. PAS DE SECONDE!
C                    => PROBLEME
                     PRINT*,'fapeetst: N1FEOC=',N1FEOC,
     %               ' PAS de SECONDE DEMI-ETOILE pour NF2=',NF2,
     %               ' TROP de TETRAEDRES RETIRES DE L''ETOILE INITIALE'
ccc                     call afetoi( n2feoc, nfetoi )

C                    TRACE DES FACES SIMPLES DE L'ETOILE DU CF NFETOI VERSION 1
                     CALL TRFETO2( PTXYZD,  NOTETR, N2FEOC, NFETOI,
     %                             NBTRCF0, NOTRCF, LEFACO, NO0FAR )

                     GOTO 9900

                  ENDIF

               ELSE

C                 RECHERCHE DU PRECEDENT NF0 DE NF2 DANS ETOILE1
                  NF0 = 0
                  NF1 = N1FEOC
 154              IF( NF1 .NE. NF2 ) THEN
C                    LA FACE SUIVANTE DE L'ETOILE
                     NF0 = NF1
                     NF1 = NFETOI(5,NF1)
                     IF( NF1 .GT. 0 ) THEN
C                       IL EXISTE UNE FACE NF1
                        GOTO 154
                     ELSE
C                       NF2 N'EST PAS UNE FACE DE ETOILE1
                        GOTO 153
                     ENDIF
                  ENDIF
c                 LE SUIVANT DE NF0 EST LE SUIVANT DE NF2
                  NFETOI(5,NF0) = NFETOI(5,NF2)

               ENDIF

C              LA FACE NF2 EST AJOUTEE A LA DEMI-ETOILE2
C              EN PREMIERE POSITION
               NFETOI(5,NF2) = N2FEOC
               N2FEOC = NF2

C              TRAITEMENT DES 2 AUTRES ARETES DE LA FACE NF2
C              QUI NE FAIT PLUS PARTIE DE L'ETOILE INITIALE
C              NUMERO NOTETR DU TETRAEDRE DE LA FACE SIMPLE
               NTE = NFETOI(1,NF2)
C              NUMERO DE LA FACE SIMPLE DANS LE TETRAEDRE NTE
               I = ABS( NFETOI(2,NF2) )

C              LES 3 SOMMETS DE LA FACE I
               NOSOTR(1) = NOTETR( NOSOFA(1,I), NTE )
               NOSOTR(2) = NOTETR( NOSOFA(3,I), NTE )
               NOSOTR(3) = NOTETR( NOSOFA(2,I), NTE )
ccc               print *,'fapeetst: lhpile=',lhpile+1,' nf2=',nf2,' nc=',nc,
ccc     %                 ' i=',i,' nosotr=',nosotr

               DO 156 L=1,3

C                 L'ARETE L DE LA FACE NF2
                  IF( L .EQ. 1 ) THEN
                     L0 = 3
                  ELSE
                     L0 = L - 1
                  ENDIF

C                 LES 2 SOMMETS DE L'ARETE L DE NF2 A TRAITER
                  NS1 = NOSOTR( L0 )
                  NS2 = NOSOTR( L  )

                  IF( L .NE. NC ) THEN

C                    LES 2 SOMMETS DE L'ARETE L DE NF2 A TRAITER
C                    NS1-NS2 EST ELLE UNE ARETE DU CF?
                     CALL ARARCF( NS1, NS2, NBCF, N1ARCF, NOARCF,
     %                            NFLCF, NAVANT, NA1 )
                     IF( NA1 .NE. 0 ) THEN
C                       PAS DE TRAVERSEE DE L'ARETE DE L'AUTRE COTE DU CF
ccc              print *,'fapeetst: arete',ns1,ns2,' du CF. face non empilee'
                        GOTO 156
                     ELSE
C                       RECHERCHE DANS ETOILE1 DE LA FACE DIFFERENTE
C                       DE NF2 OPPOSEE A L'ARETE L
                        CALL ARETOI( NS1, NS2, NOTETR, N1FEOC, NFETOI,
     %                               NF1, NC1 )
                        IF( NF1 .LE. 0 ) THEN
C                          LA FACE A DEJA ETE TRAITEE ET MISE DANS ETOILE2
ccc              print *,'fapeetst: arete',ns1,ns2,' NON dans l''etoile'
                           GOTO 156
                        ENDIF

C                       CETTE FACE OPPOSEE PAR L'ARETE L EST EMPILEE
                        IF( LHPILE .GE. MXPILE ) THEN
                           print *,'fapeetst: AUGMENTER MXPILE=',MXPILE
                           GOTO 156
                        ENDIF
                        LHPILE = LHPILE + 1
C                       LA FACE NFETOI
                        LAPILE( 1, LHPILE ) = NF1
C                       LE NUMERO DE L'ARETE 1:1-2, 2:2-3 3:3-1 DE TRAVERSEE
                        LAPILE( 2, LHPILE ) = NC1
ccc               print *,'fapeetst: arete',ns1,ns2,' empilee LHPILE=',LHPILE

                     ENDIF

                  ELSE
ccc               print *,'fapeetst: arete',ns1,ns2,' depilee'

                  ENDIF

 156           ENDDO

C              ALLER EN HAUT DE PILE
               GOTO 153

            ENDIF

C           AJOUT DES NBTRCF TRIANGLES DU CF POUR FERMER LA DEMI-ETOILE 2
C           -------------------------------------------------------------
            CALL AJTRCF( NS1ACF, NS2ACF, PTXYZD,
     %                   NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                   NOTETR, N2FEOC, N1FEVI, MXFETO, NFETOI, IERR )

C           VERIFICATION FINALE DE LA DEMI-ETOILE N2FEOC: A supprimer apres tests
C           TOUTE ARETE DE L'ETOILE APPARTIENT ELLE A 2n FACES?
            CALL VETAFET( NOTETR, N2FEOC, NFETOI, NBFETO, NBARPB )
            IF( NBARPB .GT. 0 ) THEN
               PRINT*,'fapeetst: FACE PERDUE',NFP,NBARPB,
     %             ' ARETES NON DANS 2n FACES'
C              TRACE DE LA DEMI-ETOILE 2 EN VERSION 2 DE NFETOI
               TRACTE = .TRUE.
               CALL TRFETO2( PTXYZD, NOTETR, N2FEOC, NFETOI,
     %                       0,      NOTRCF, LEFACO, NO0FAR )
            ENDIF

            CALL SOUSETO( PTXYZD, NOTETR, N2FEOC, NFETOI,
     %                    MXSSET, NBSSET, N1SSET, IERR )
            IF( IERR .NE. 0 ) THEN
               PRINT*,'fapeetst: sortie souseto avec IERR=',IERR
               GOTO 9900
            ENDIF

C           AJOUT DES NBTRCF TRIANGLES DU CF POUR FERMER LA DEMI-ETOILE 1
C           -------------------------------------------------------------
            CALL AJTRCF( NS2ACF, NS1ACF, PTXYZD,
     %                   NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                   NOTETR, N1FEOC, N1FEVI, MXFETO, NFETOI, IERR )
C           VERIFICATION FINALE DE LA DEMI-ETOILE N1FEOC: A supprimer apres tests
C           TOUTE ARETE DE L'ETOILE APPARTIENT ELLE A 2n FACES?
            CALL VETAFET( NOTETR, N1FEOC, NFETOI, NBFETO, NBARPB )
            IF( NBARPB .GT. 0 ) THEN
               PRINT*,'fapeetst: FACE PERDUE',NFP,NBARPB,
     %             ' ARETES NON DANS 2n FACES'
C              TRACE DE LA DEMI-ETOILE 1 EN VERSION 2 DE NFETOI
               TRACTE = .TRUE.
               CALL TRFETO2( PTXYZD, NOTETR, N2FEOC, NFETOI,
     %                       0,      NOTRCF, LEFACO, NO0FAR )
            ENDIF

            CALL SOUSETO( PTXYZD, NOTETR, N1FEOC, NFETOI,
     %                    MXSSET, NBSSET, N1SSET, IERR )
            IF( IERR .NE. 0 ) THEN
               PRINT*,'fapeetst: SORTIE souseto avec IERR=',IERR
               GOTO 9900
            ENDIF

C        2.7) CALCUL DU VOLUME DES TETRAEDRES DU TABLEAU NOTECF
C             C-A-D DES 2 DEMI-ETOILES OU DE L'ETOILE ENCOCHEE
C             C-A-D DE L'ETOILE ACTUELLE DES TETRAEDRES NOTECF
C             COMPRESSION DES TETRAEDRES SUPPRIMES DES 2 DEMI
C             PROTECTION DES TETRAEDRES INITIAUX DANS NOTECFPR
C        -----------------------------------------------------
C        LE VOLUME DES 2 DEMI-ETOILES AVANT LEUR TETRAEDRISATION
         VOLET0 = 0D0
         QUAMIN = 2.0
         QUAMOY = 0.0

C        N1TETS SERA MIS A JOUR AVEC LES NOUVEAUX TETRAEDRES
ccc         print*,'fapeetst: LISTE des TETRAEDRES A RE-TETRAEDRISER'

         NBTECFPR = 0
         NBT = 0
         DO K = 1, NBTECF

C           NUMERO NOTETR DU K-EME TETRAEDRE DE L'ETOILE
            NTE = NOTECF( K )

C           PROTECTION EN CAS DE PROBLEME LORS DE LA TETRAEDRISATION
C           DE L'ETOILE POUR REVENIR EN ARRIERE
            IF( NTE .GT. 0 ) THEN
               NBT = NBT + 1
               NOTECF  ( NBT ) = NTE
               NOTECFPR( NBT ) = NTE
C              SON VOLUME ET SA QUALITE
               CALL QUATETD( PTXYZD(1,NOTETR(1,NTE)),
     %                       PTXYZD(1,NOTETR(2,NTE)),
     %                       PTXYZD(1,NOTETR(3,NTE)),
     %                       PTXYZD(1,NOTETR(4,NTE)),
     %                       ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
ccc            IF( VOLUTE .LE. 0D0 ) THEN
ccc               print*,'fapeetst: TETRAEDRE INITIAL',NTE,
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
         NBTECF   = NBT
         NBTECFPR = NBT
         VMOYEN = VOLET0 / NBT
         QUAMOY = QUAMOY / NBT
ccc         print *,'fapeetst: FACE PERDUE',NFP,' NBTET0=',NBTECFPR,
ccc     %           ' VOLET0=',VOLET0,' VMOYEN0=',VMOYEN,
ccc     %           ' QualMOY0=',QUAMOY,' QualMIN0=',QUAMIN

C        TRACE DES NBTECF TETRAEDRES A RE-TETRAEDRISER
         KTITRE = '      TETRAEDRES INITIAUX A RE-TETRAEDRISER'
         WRITE(KTITRE(1:5),'(I5)') NBTECFPR
         CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECFPR, NOTECFPR, NOTETR )


C        2.8) TETRAEDRISATION DES FACES DES NBSSET SOUS-ETOILES
C             SANS OU AVEC DES POINTS INTERNES AJOUTES
C        ------------------------------------------------------
C        VOLUME APRES TETRAEDRISATION DE L'ETOILE
         NBREPR = 0
         VOLET1 = 0D0
         CALL TETRETO( QUAMINEX, QTEAME, KTITRE,
     %                 MXSOMM, NBSOMM0,NBSOMM, PTXYZD,
     %                 NPSOFR, VMOYEN,
     %                 MXTETR, N1TEVI, NOTETR, NUDTETR, N1TETS,
     %                 NBTRCF, NOTRCF, INFACO, MXFACO,  LEFACO, NO0FAR,
     %                 N1FEVI, MXFETO, NFETOI,
     %                 MXSSET, NBSSET, N1SSET, NBTECFPR,
     %                 MXTECF, NBTECR, NOTECF, VOLET1, IERR )

C        RE CALCUL DU VOLUME DES TETRAEDRES CREES pour VERIFICATION
         VOLET1 = 0D0
         QUAMIN = 2.0
         QUAMOY = 0.0
ccc         print*,'fapeetst: LISTE des TETRAEDRES DE LA RE-TETRAEDRISATION'
         VOLET1 = 0D0
         NBT = 0
         DO K = 1, NBTECR

C           NUMERO NOTETR DU K-EME TETRAEDRE DE L'ETOILE
            NTE = NOTECF( K )

C           PROTECTION EN CAS DE PROBLEME LORS DE LA TETRAEDRISATION
C           DE L'ETOILE POUR REVENIR EN ARRIERE
            IF( NTE .GT. 0 ) THEN
               NBT = NBT + 1
               NOTECF( NBT ) = NTE
C              SON VOLUME ET SA QUALITE
               CALL QUATETD( PTXYZD(1,NOTETR(1,NTE)),
     %                       PTXYZD(1,NOTETR(2,NTE)),
     %                       PTXYZD(1,NOTETR(3,NTE)),
     %                       PTXYZD(1,NOTETR(4,NTE)),
     %                       ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
ccc            IF( VOLUTE .LE. 0D0 ) THEN
ccc               print*,'fapeetst: TETRAEDRE FINAL  ',NTE,
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

C        VOLUME MOYEN DES TETRAEDRES A RE-TETRAEDRISER
         NBTECR = NBT
         VMOYEN = VOLET1 / NBT
         QUAMOY = QUAMOY / NBT
ccc         print *,'fapeetst: FACE PERDUE',NFP,' NBTECR=',NBTECR,
ccc     %           ' VOLET1=',VOLET1,' VMOYEN1=',VMOYEN,
ccc     %           ' QualMOY1=',QUAMOY,' QualMIN1=',QUAMIN

         IF( IERR .NE. 0 ) THEN
            GOTO 9900
         ENDIF

C        CONTROLE DU NOMBRE NBTECR DE TETRAEDRES CREES DANS L'ETOILE
         IF( NBTECR .GE. MXTECF ) THEN
            print*,'???????????????????????????????????????????????????'
           print*,'PB fapeetst: TROP DE TETRAEDRES CREES     ',NBTECR
           print*,'PB fapeetst: NOMBRE DE TETRAEDRES INITIAUX',NBTECFPR
            print*,'BOUCLE INFINIE DE CREATION???'
            print*,'???????????????????????????????????????????????????'
            KTITRE='PB:       TETRAEDRES EN BOUCLE INFINIE DE CREATION?'
            GOTO 9900
         ENDIF


C        2.9) CONTROLE DU VOLUME DE L'ETOILE AVANT et APRES
C        --------------------------------------------------
         print *,'fapeetst: NFP=',NFP,' VOLET0=',VOLET0,
     %         ' pour',NBTECFPR,' TETRAEDRES INITIAUX'
         print *,'fapeetst: NFP=',NFP,' VOLET1=',VOLET1,
     %         ' pour',NBTECR,' TETRAEDRES FINAUX'

         IF( ABS( VOLET1 - VOLET0 ) .GT. 1D-5 * VOLET0 ) THEN

          print*,'?????????????????????????????????????????????????????'
       print*,'PROBLEME fapeetst: VOLUMES TROP DIFFERENTS VOLET0=',
     %              VOLET0,'   VOLET1=',VOLET1
       print*,'PROBLEME fapeetst: VOLUMES TROP DIFFERENTS VOLET1=',
     %              VOLET1     
         print*,'?????????????????????????????????????????????????????'
          KTITRE='PB:       TETRAEDRES DE VOLUME DIFFERENT DES INITIAUX 
     %V0=               V1=           '
            WRITE(KTITRE(58:71),'(G14.8)') VOLET0
            WRITE(KTITRE(76:89),'(G14.8)') VOLET1
            GOTO 9900

         ENDIF


C        2.10) COMPLETION DES CHAINAGES DES TETRAEDRES CREES DANS
C              LES 2 DEMI-ETOILES OU L'UNIQUE ETOILE ENCOCHEE
C        --------------------------------------------------------
         CALL MJOPTE( NBTECR, NOTECF, N1TETS, NOTETR, MXTETR,
     %                N1TEVI, PTXYZD, NBFANR )
         IF( NBFANR .LE. 0 ) GOTO 240

C        DESSOUS PEUT ETRE A SUPPRIMER
         LHPILE = NBTECR
 202     IF( LHPILE .GT. 0 ) THEN

C           LE HAUT DE LA PILE DES TETRAEDRES CREES DANS L'ETOILE
            NTE = ABS( NOTECF( LHPILE ) )

C           LE TETRAEDRE EST DEPILE
            LHPILE = LHPILE - 1

            IF( NTE .LE. 0 ) GOTO 202

            DO 220 J=1,4

C              FACE J DU TETRAEDRE NTE et NUMERO DE SES 3 SOMMETS
               DO K=1,3
                  NOSOTR(K) = NOTETR( NOSOFA(K,J), NTE )
               ENDDO
C              TRI CROISSANT DES 3 SOMMETS DE LA FACE J DE NTE
               CALL TRI3NO( NOSOTR, NOSOTR )

C              LA FACE J DE NTE EST ELLE UNE FACE PERDUE DANS LEFACO?
               DO 205 M = 1, NBTRCF0

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
C              NOSOFA/1,3,2, 2,3,4, 3,1,4, 4,1,2/
 208           IF( NOTETR(J+4,NTE) .GT. 0 ) GOTO 220

C              RECHERCHE DE CETTE FACE NOSOTR PARMI LES FACES
C              DES TETRAEDRES CREES DANS L'ETOILE
               DO 210 I=LHPILE,1,-1

                  NTE1 = ABS( NOTECF( I ) )
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

               PRINT *, 'fapeetst: FACE ',NOSOTR, 'du TETRAEDRE',
     %                  (NOTETR(M,NTE),M=1,8),' SANS TETRAEDRE OPPOSE'

 220        ENDDO

C           RETOUR EN HAUT DE PILE
            GOTO 202
         ENDIF

C        VERIFICATION DE LA COMPLETION DES TETRAEDRES OPPOSES
 240     CALL MJOPTE( NBTECR, NOTECF, N1TETS, NOTETR, NUDTETR,
     %                N1TEVI, PTXYZD, NBFANR )
         PRINT*, 'fapeetst: sortie mjopte avec NBFANR=',NBFANR
ccc         CALL VEOPTE( NBTECR, NOTECF, NOTETR, PTXYZD, NBFANR )
         IF( NBFANR .GT. 0 ) THEN
            IERR = NBFANR
            GOTO 9900
         ENDIF

C        DESTRUCTION DES NBTECFPR TETRAEDRES NOTECFPR DE PROTECTION
         DO K=NBTECFPR,1,-1
            NTE = ABS( NOTECFPR( K ) )
            NOTETR(1,NTE) = 0
            NOTETR(5,NTE) = N1TEVI
            N1TEVI = NTE
         ENDDO


C        2.11) SUPPRESSION DANS NOFAPE DES FACES PERDUES RETROUVEES OU A EVITER
C        ----------------------------------------------------------------------
         DO 250 K=1,NBTRCF0

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


C        2.12) QUALITE et TRACE DES TETRAEDRES REMPLISSANT L'ETOILE OU ACTUELS
C              MISE A JOUR DU NO DE TETRAEDRE DES FACES DE LEFACO
C        ---------------------------------------------------------------------
         IF( IERR .EQ. 0 .AND. NBTECR .GT. 0 ) THEN

            QUAMOY = 0.0
            QUAMIN = 2.0
            DO K = 1, NBTECR

C              NUMERO NOTETR DU TETRAEDRE CREE
               NTE = ABS( NOTECF( K ) )

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
ccc109900       FORMAT('FIN fapeetst: FACE PERDUE',I5,' TETRAEDRE CREE',I4,
ccc     %             ' NOTETR(',I7,')=',4I7,3X,4I7 )

            ENDDO
            QUAMOY = QUAMOY / NBTECR

            print *,'fapeetst: FACE PERDUE',NFP,
     %              ' avec NBTRCF0=',NBTRCF0,
     %              ' a DONNE',NBTECR,' tetra de QUALITE MOYENNE',
     %              QUAMOY,' MINIMALE',QUAMIN,' avec',NBSOMM-NBSOMM0,
     %              ' SOMMETS AJOUTES pour',
     %              NBTECFPR,' tetra initiaux'

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
                     PRINT *,'fapeetst: FACE',NTR,
     %                '     RETROUVEE dans LEFACO=',
     %                    (LEFACO(L,NTR),L=1,11)
                     GOTO 965

                  ENDIF

               ENDIF

C              LA FACE NTR DE LEFACO N'EST PAS RETROUVEE
               LEFACO(11,NTR) = 0
               PRINT *,'fapeetst: Face',NTR,
     %                 ' NON RETROUVEE dans LEFACO=',
     %                 (LEFACO(L,NTR),L=1,11)

 965        ENDDO

C           TRACE DES NBTECR TETRAEDRES CREES POUR CETTE FACE PERDUE NFP
            KTITRE='      TETRAEDRES CREES dans ETOILE FINALE QUALITE MO
     %Y=        QUALITE MIN=      '
            WRITE(KTITRE(1:5),    '(I5)') NBTECR
            WRITE(KTITRE(55:61),'(F7.5)') QUAMOY
            WRITE(KTITRE(75:81),'(F7.5)') QUAMIN
            CALL TRFETO6( KTITRE,  PTXYZD,
     %                    NBTRCF0, NOTRCF, LEFACO,
     %                    NBTECR,  NOTECF, NOTETR )

         ENDIF

ccc         IF( N1TEVI .LE. 0 ) THEN
ccc            CALL MJN1TEVI( 'fapeetst', MXTETR, NOTETR, N1TEVI, NUDTETR )
ccc            IF( N1TEVI .LE. 0 ) THEN
ccc               PRINT*,'fapeetst: SATURATION DU TABLEAU NOTETR. MXTETR=',
ccc     %                 MXTETR
ccc               IERR = 2
ccc               GOTO 9900
ccc            ENDIF
ccc         ENDIF


C        PASSAGE A LA FACE LEFACO PERDUE SUIVANTE
C        ----------------------------------------
         GOTO 1

C        ===========================================================
C        AFFICHAGE DU PROBLEME RENCONTRE & ABANDON DE LA FACE PERDUE
C        ===========================================================
 9900    print*,'fapeetst: FINALEMENT ABANDON de la FACE PERDUE',NFP,
     %  ' avec NBTRCF=',NBTRCF0,' et',
     %   NBSOMM-NBSOMM0,' SOMMETS AJOUTES  IERR=',IERR
         NBABAN = NBABAN + 1

C        TRACE DES TETRAEDRES INITIAUX et CREES AVANT ABANDON
         TRACTE0 = TRACTE
         TRACTE = .TRUE.
         KTITRE='fapeetst:           TETRAEDRES INITIAUX suite ABANDON'
         WRITE(KTITRE(10:16),'(I7)') NBTECFPR
         CALL TRFETO8( KTITRE,   PTXYZD,
     %                 NBTRCF,   NOTRCF,  LEFACO, NO0FAR,
     %                 NBTECFPR, NOTECFPR,
     %                 NBTECR,   NOTECF,  NOTETR )
C        ABANDON DU TRAITEMENT DE LA FACE PERDUE
         TRACTE = TRACTE0

         IF( NBTECFPR .GT. 0 ) THEN

C           DETRUIRE  LES NBTECR   TETRAEDRES NOTECF,
C           RESTAURER LES NBTECFPR TETRAEDRES NOTECFPR
C           METTRE A JOUR N1TEVI, NOTETR, N1TETS
            CALL RSTETO( PTXYZD, NBTECFPR, NOTECFPR, NBTECR,  NOTECF,
     %                   MXTETR, N1TEVI,   NOTETR,   NUDTETR, N1TETS,
     %                   NBFANR )

C           MISE A JOUR DE LEFACO POUR LES TETRAEDRES A DETRUIRE
            DO K=1,NBTECR
               NTE = NOTECF( K )
C              SUPPRIMER EVENTUELLEMENT LES 4 FACES DU TETRAEDRE NTE
C              DANS LE TABLEAU LEFACO
               CALL SUTELEFA( NTE, NOTETR, INFACO, MXFACO, LEFACO )
            ENDDO

C           MISE A JOUR DE LEFACO POUR LES TETRAEDRES A RESTAURER
            DO K=1,NBTECFPR
               NTE = NOTECFPR( K )
C              AJOUTER EVENTUELLEMENT LES 4 FACES DU TETRAEDRE NTE
C              DANS LE TABLEAU LEFACO
               CALL AJTELEFA( NTE, NOTETR, INFACO, MXFACO, LEFACO )
            ENDDO

C           PLUS AUCUN TETRAEDRE CREE APRES leur DESTRUCTION
            NBTECR = 0
            IF( NBFANR .GT. 0 ) THEN
               PRINT*, 'fapeetst: PB RESTAURATION TETRA INITIAUX avec mj
     %opte => NBFANR=',NBFANR,' ?'
            ENDIF

         ENDIF

C        NOMBRE DE TETRAEDRES CREES
         NBTECR = 0
         NBTECF = 0

         DO 950 K = 1, NBTRCF0

C           LA FACE PERDUE ET ABANDONNEE EST EFFACEE DANS NOFAPE
            NTR = ABS( NOTRCF( K ) )

C           LE SIGNE + LUI EST RENDU
            NOTRCF( K ) = NTR

            DO J = 1, NBFAPE
               IF( NTR .EQ. ABS( NOFAPE(J) ) ) THEN
                  NOFAPE( J ) = -NTR
                  GOTO 950
               ENDIF
            ENDDO

 950     ENDDO

C        RECUPERATION DES POINTS AJOUTES POUR RETROUVER CETTE FACE PERDUE
         NBSOMM = NBSOMM0
         GOTO 1

      ENDIF

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
      PRINT*,'fapeetst: SORTIE avec NBSOMM=',NBSOMM,' et ',NBABAN,
     %' ABANDONS de FACES PERDUES sur',NBFAPE,
     %' FACES INITIALEMENT PERDUES ++++++++++++++++++++++++++++++++++++'
      PRINT*

      TRACTE = .FALSE.
      RETURN
      END
