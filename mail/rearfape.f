      SUBROUTINE REARFAPE( MXSOMM, PTXYZD, N1TETS,
     %                     MXTETR, NOTETR, N1TEVI, NUDTETR,
     %                     NBTRCF, NOTRCF, INFACO, MXFACO, LEFACO,
     %                     IVOLTE, NVOLTE, MXTE1S, NOTE1S,
     %                     MXPILE, NUPILE, MXFETO, NFETOI,
     %                     NBTARRE, IERR,  NAFACO, NFLPER )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ESSAI DE RECUPERATION D'ARETES PERDUES DES NBTRCF FACES
C -----    PERDUES NOTRCF DE LEFACO par
C          2Tetraedres avec 1Face Commune -> 1Arete commune a 3Tetraedres
C          ou RE TETRAEDRISATION DES TETRAEDRES INTERSECTES PAR L'ARETE
C          ou AJOUT du POINT MILIEU D'UNE ARETE PERDUE de la TETRAEDRISATION

C ENTREES:
C --------
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS PERMIS POUR LA TETRAEDRISATION
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES DANS LE TABLEAU NOTETR

C NBTRCF : NOMBRE DE FACES LEFACO PERDUES AYANT AU MOINS UNE ARETE COMMUNE
C          ET BORDE d'UN CF D'ARETES DANS LA TETRAEDRISATION
C NOTRCF : NO LEFACO DES NBTRCF FACES PERDUES D'ARETES A RETROUVER DANS
C          UN TETRAEDRE ACTUEL
C INFACO : =1 LE TABLEAU LEFACO DOIT ETRE PRESENT
C MXFACO : NOMBRE MAXIMAL DE FACES DECLARABLES DANS LEFACO
C LEFACO : LES 3 SOMMETS, 2 MATERIAUX, 3 FACES VOISINES ET CHAINAGE
C          DES FACES TRIANGULAIRES DU CONTOUR ET INTERFACES

C IVOLTE : 0 PAS       DE TABLEAU NVOLTE A L'APPEL
C          1 EXISTENCE DU TABLEAU NVOLTE A L'APPEL
C NVOLTE : NUMERO DU VOLUME (1 A NBVOPA) DE CHAQUE TETRAEDRE
C          -1 SI VOLUME INCONNU
C MXTE1S : MAX DE MOTS DU TABLEAU NOTE1S
C MXPILE : MAX DE MOTS DU TABLEAU NUPILE
C MXFETO : NOMBRE MAXIMAL DE FACES DU TABLEAU NFETOI

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
C NOTE1S : TABLEAU DES FACES CONTENANT UN SOMMET
C NUPILE : TABLEAU DU NO PTXYZD DES SOMMETS D'ARETE POUR RECUPERER
C          DES ARETES DANS LA TETRAEDRISATION
C          puis NUMERO NOTETR DES TETRAEDRES INTERSECTES PAR UNE ARETE

C SORTIES:
C --------
C NBSOMM : NUMERO DU DERNIER SOMMET AJOUTE
C NBTARRE: NOMBRE D'ARETES DES NBTRCF FACES PERDUES RECUPEREES
C          DANS LA TETRAEDRISATION
C IERR   : 0 SI PAS D'ERREUR
C          1 SI RETOUR EN ARRIERE A PROGRAMMER
C          2 SI SATURATION DU TABLEAU NOTETR
C          3 SI SATURATION DU TABLEAU NUPILE
C         17 SI DEMANDE d'AJOUT du MILIEU de l'ARETE NAFACO de LEFACO(NFLPER)
C NAFACO : NUMERO DE L'ARETE DE LA FACE NFLPER A AJOUTER SI IERR=11
C NFLPER : NUMERO DE LA FACE LEFACO D'ARETE NAFACO DE MILIEU A AJOUTER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint Pierre du Perray          Fevrier 2020
C....................................................................012
      PARAMETER        (MXXYZP=96, MXAR2E=192)
      include"./incl/gsmenu.inc"
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0
      CHARACTER*96      KTITRE
      DOUBLE PRECISION  PTXYZD(1:4,1:MXSOMM),
     %                  CBTR(3), XYZP(3,MXXYZP),
     %                  V, VOLTET, VOLET0, VOLET1
      INTEGER           NOTRCF(NBTRCF), NOTETR(8,MXTETR),
     %                  N1TETS(MXSOMM), LEFACO(11,0:MXFACO),
     %                  NVOLTE(*), NOTE1S(MXTE1S), NUPILE(MXPILE),
     %                  NFETOI(5,MXFETO), NOAR2E(3,MXAR2E)
      INTEGER           NTNOUV(3), NOSOTR(3),
     %                  NOSOFATE(3,4)
C     NO DES SOMMETS DES 4 FACES POUR QUE LES SOMMETS SOIENT VUS
C     DE L'INTERIEUR DU TETRAEDRE DANS LE SENS DIRECT I.E.
C     LEUR VECTEUR NORMAL POINTE VERS L'INTERIEUR DU TETRAEDRE
      DATA              NOSOFATE/ 1,2,3, 2,4,3, 3,4,1, 4,2,1 /

      LORBITE = 1
      TRACTE0 = TRACTE
      NBTARRE = 0

C     BOUCLE SUR LES NBTRCF FACES LEFACO PERDUES
C     ==========================================
      DO 1000 NFP = 1, NBTRCF

C     NUMERO LEFACO DE LA FACE PERDUE NFP
      NFLPER = NOTRCF( NFP )

ccc      PRINT*,'rearfape: DEBUT ......... FACE PERDUE',NFLPER,' St:',
ccc     %       (LEFACO(K,NFLPER),K=1,3),'  ..............................'

C     BOUCLE SUR LES 3 ARETES DE LA FACE PERDUE NFLPER DE LEFACO
C     POUR DETECTER LES ARETES NON DANS LA TETRAEDRISATION
C     ----------------------------------------------------------
      NBXYZP  = 0
      NBARTE0 = 0
      NO1ERA  = 0

C     RECHERCHE D'UNE ARETE DE LA FACE PERDUE NFLPER NON ARETE
C     DE LA TETRAEDRISATION
      DO I1=1,3

C        L'ARETE I1 DE LA FACE PERDUE NFLPER DE SOMMETS NS1 NS2
         IF( I1 .EQ. 3 ) THEN
            I2 = 1
         ELSE
            I2 = I1+1
         ENDIF
         NS1 = LEFACO( I1, NFLPER )
         NS2 = LEFACO( I2, NFLPER )

C        L'ARETE NS1-NS2 EST ELLE DANS LA TETRAEDRISATION?
         CALL TETR1A( NS1,    NS2,    N1TETS, NOTETR,
     %                NBTE1A, MXTE1S, NOTE1S, IERR )

         IF( NBTE1A .GT. 0 ) THEN
C           OUI: L'ARETE NS1-NS2 EST DEJA DANS LA TETRAEDRISATION
            NBARTE0 = NBARTE0 + 1
ccc            PRINT*,'rearfape: l''ARETE INITIALE',NS1,NS2,
ccc     %             ' EST dans',NBTE1A,' TETRAEDRES'
         ELSE
            IF( NO1ERA .EQ. 0 ) THEN
C              NO DE LA PREMIERE ARETE HORS TETRAEDRISATION
               NO1ERA = I1
            ENDIF
         ENDIF

      ENDDO
ccc      PRINT*,'rearfape: la FACE',NFLPER,' a',NBARTE0,
ccc     %       ' ARETES INITIALES dans la TETRAEDRISATION'

      IF( NO1ERA .EQ. 0 ) THEN
C        LES 3 ARETES DE LA FACE NFLPER SONT DANS LA TETRAEDRISATION
         GOTO 1000
      ENDIF


C     BOUCLE SUR LES 3 ARETES DE LA FACE PERDUE NFLPER DE LEFACO
C     ----------------------------------------------------------
      DO 100 NAFACO = NO1ERA, 3

C        L'ARETE NAFACO DE LA FACE PERDUE NFLPER A POUR SOMMETS NS1 NS2
         NBPENS12 = 0
         NBCHNS12 = 0
         IF( NAFACO .EQ. 3 ) THEN
            I2 = 1
         ELSE
            I2 = NAFACO+1
         ENDIF
         NS1 = LEFACO( NAFACO, NFLPER )
         NS2 = LEFACO(     I2, NFLPER )
         NS3 = 0

C        PREMIERE METHODE POUR RECUPERER L'ARETE NAFACO DE LA FACE PERDUE NFLPER
C        =======================================================================
 8       METHOD = 1
ccc         print*
ccc         PRINT*,'rearfape: Methode 1 sur l''ARETE NS1=',NS1,' NS2=',NS2,
ccc     %          ' Nombre PERMUTATIONS NS1<->NS2 =',NBCHNS12

C        INITIALISATION DE LA PILE DES SOMMETS NS3 POUR FAIRE
C        APPARAITRE L'ARETE NS1-NS3 DANS LA TETRAEDRISATION
         LHPILE = 1
         NUPILE( LHPILE ) = NS2

C        L'ARETE NS1-NS2 EST ELLE DANS LA TETRAEDRISATION?
C        -------------------------------------------------
 10      CALL TETR1A( NS1,    NS2,    N1TETS, NOTETR,
     %                NBTE1A, MXTE1S, NOTE1S, IERR )

         IF( NBTE1A .GT. 0 ) THEN
C           OUI: L'ARETE NS1-NS2 EST RETROUVEE DANS LA TETRAEDRISATION
ccc            PRINT*,'rearfape: FIN de PILE: ARETE NS1=',NS1,'- NS2=',NS2,
ccc     %             ' EST DANS',NBTE1A,' TETRAEDRES'

C           TRACE DES NBTE1A TETRAEDRES D'ARETE NS1-NS2
          KTITRE='rearfape:        tetraedres d''ARETE       -         '
            WRITE(KTITRE(11:16),'(I6)') NBTE1A
            WRITE(KTITRE(38:43),'(I6)') NS1
            WRITE(KTITRE(46:51),'(I6)') NS2
            CALL SANSDBL(  KTITRE, NC )
            CALL TRATE2AR( KTITRE(1:NC), PTXYZD, NBTE1A, NOTE1S, NOTETR,
     %                     NS1, NS2, 0, 0 )
            GOTO 99
         ENDIF

         IF( LHPILE .LE. 0 ) THEN
            PRINT*,'rearfape: NUPILE VIDE ???  ARETE NS1=',NS1,
     %             ' NS2=',NS2,' NS3=',NS3
            GOTO 52
         ENDIF

C        LE SOMMET FINAL NS3 DE L'ARETE NS1-NS3 A FAIRE APPARAITRE
C        DANS LA TETRAEDRISATION EST EN HAUT DE PILE
C        LA METHODE 2Tetra->3Tetra EST UTILISEE RECURSIVEMENT...
C        ---------------------------------------------------------
         NS3 = NUPILE( LHPILE )

C        NON: L'ARETE NS1-NS2 N'EST PAS DANS LA TETRAEDRISATION
C        RECHERCHE DES NBTNS1 TETRAEDRES DE SOMMET NS1
         CALL TETR1S( NS1,    N1TETS, NOTETR,
     %                NBTNS1, MXTE1S, NOTE1S, IERR )
         IF( NBTNS1 .LE. 0 ) THEN
            PRINT*,'rearfape: Bizarre AUCUN TETRAEDRE de SOMMET',NS1
            GOTO 40
         ENDIF

C        PARCOURS DES TETRAEDRES DE SOMMET NS1 POUR TROUVER
C        LA FACE PERFOREE PAR L'ARETE NS1-NS3
C        --------------------------------------------------
         DO 30 M = 1, NBTNS1

C           LE TETRAEDRE M DE SOMMET NS1
            NTE = NOTE1S( M )
            IF( NOTETR(1,NTE) .LE. 0 ) GOTO 30

C           RECHERCHE DU TETRAEDRE NTOP OPPOSE A LA FACE OPPOSEE
C           AU SOMMET NS1 DE NTE: RECHERCHE DU NO DE SOMMET NS1 DANS NTE
            DO K=1,4
               IF( NS1 .EQ. NOTETR(K,NTE) ) GOTO 20
            ENDDO
            GOTO 30

C           LE NO DE LA FACE OPPOSEE AU SOMMET K EST K MOD 4 + 1
 20         IF( K .NE. 4 ) THEN
               NFTE=K+1
            ELSE
               NFTE=1
            ENDIF

C           LA FACE NFTE DU TETRAEDRE NTE EST ELLE PERFOREE
C           PAR L'ARETE NS1-NS3?
C           NUMERO DES 3 SOMMETS DE LA FACE NFTE
            NOSOTR(1) = NOTETR( NOSOFATE(1,NFTE), NTE )
            NOSOTR(2) = NOTETR( NOSOFATE(2,NFTE), NTE )
            NOSOTR(3) = NOTETR( NOSOFATE(3,NFTE), NTE )

C           LE TETRAEDRE NTE A T IL POUR ARETE NS1-NS3?
            DO I=1,3
               IF( NOSOTR(I) .EQ. NS3 ) THEN
C                 OUI: LE TETRAEDRE NTE A POUR ARETE NS1-NS3
C                 L'ARETE NS1-NS2 EST INVERSEE ET REDEPART
                  NBCHNS12 = NBCHNS12 + 1
                  IF( NBCHNS12 .GE. 4 ) THEN
C                    LE PROCESS BOUCLE. PASSAGE A LA METHODE SUIVANTE
                     GOTO 52
                  ENDIF
                  K   = NS1
                  NS1 = NS2
                  NS2 = K
                  GOTO 8
               ENDIF
            ENDDO

C           INTERSECTION OU NON DE NS1-NS3 AVEC LA FACE NFTE DE NTE?
            CALL VOLU2T3T( NS1, NS3, NOSOTR, PTXYZD, NONOUI )
C           NONOUI: -1 NS1 EST UN SOMMET DU TRIANGLE NOSOTR
C                   -2 NS3 EST UN SOMMET DU TRIANGLE NOSOTR
C                    0 SI NS1-NS3 N'INTERSECTE PAS LE TRIANGLE NOSOTR
C                    1 NS1-NS3 INTERSECTE LE TRIANGLE FERME (ARETES COMPRISES)
C                      car LE VOLUME 2T = VOLUME 3T en VALEUR ABSOLUE
C            verifier que NONOUI et LINTER SONT COHERENTS...

C           CETTE FACE NOSOTR EST ELLE INTERSECTEE PAR L'ARETE NS1-NS3?
            CALL INARTR( PTXYZD(1,NS1), PTXYZD(1,NS3),
     %                   PTXYZD( 1, NOSOTR(1) ),
     %                   PTXYZD( 1, NOSOTR(2) ),
     %                   PTXYZD( 1, NOSOTR(3) ),
     %                   LINTER, XYZP(1,1), CBTR )
C           LINTER : -2 SI S1=S3
C           -1 SI S1-S3 PARALLELE AU PLAN DU TRIANGLE NOSOTR
C            0 SI S1-S3 N'INTERSECTE PAS LE TRIANGLE ENTRE CES 3 SOMMETS
C            1 SI S1-S3   INTERSECTE     LE TRIANGLE ENTRE CES 3 SOMMETS
C              XYZP: 3 COORDONNEES DU POINT D'INTERSECTION SI LINTER=1
C              CBTR: 3 COORDONNEES BARYCENTRIQUES DE XYZP DANS LE TRIANGLE
            IF( LINTER .LT. 0 ) THEN
C              POUR LE TEST DE VERIFICATION QUI SUIT... A SUPPRIMER ENSUITE
               LINTER = 0
            ENDIF
            IF( NONOUI .NE. LINTER ) THEN
       print*,'rearfape: RESULTATS DIFFERENTS entre volu2t3t et inartr!'
            print*,'rearfape: a voir! NONOUI=',NONOUI,' LINTER=',LINTER,
     %             ' NS1=',NS1,' NS2=',NS2,' NOSOTR=',NOSOTR,
     %             ' CBTR=',CBTR
            ENDIF

            IF( LINTER .LE. 0 ) THEN
C              PAS D'INTERSECTION DE NS1-NS3 AVEC LE TRIANGLE NOSOTR
C              PASSAGE AU TETRAEDRE DE SOMMET NS1 SUIVANT
               GOTO 30
            ENDIF

C           LINTER=1: LA FACE NOSOTR EST INTERSECTEE PAR L'ARETE NS1-NS3
C                     AU POINT XYZ(3)
C           NTOP LE TETRAEDRE OPPOSE A NS1 ET LE NUMERO NFOP DE LA FACE
C           QUI EST LA FACE NFTE DU TETRAEDRE NTE PERCEE PAR L'ARETE NS1-NS3
            CALL NOFAOP( NFTE, NTE, NOTETR, NFOP, NTOP )
C           LE SOMMET OPPOSE NSOP EST NFOP - 1
            IF( NFOP .NE. 1 ) THEN
               NSOP = NFOP - 1
            ELSE
               NSOP = 4
            ENDIF

C           LE NO NS4 DU SOMMET DU TETRAEDRE NTOP OPPOSE A SA FACE NFOP
            NS4 = NOTETR( NSOP, NTOP )

C           TRACE DES NBTNS1 TETRAEDRES DE SOMMET NS1, NTOP ET LES 2 ARETES
            KTITRE='rearfape:        tetraedres ARETES NS1=       NS2=    
     %     NS3=       '
            WRITE(KTITRE(11:16),'(I6)') NBTNS1
            WRITE(KTITRE(40:45),'(I6)') NS1
            WRITE(KTITRE(51:56),'(I6)') NS2
            WRITE(KTITRE(62:67),'(I6)') NS3
            CALL SANSDBL(   KTITRE, NC )
            CALL TRATES2AR( KTITRE(1:NC), PTXYZD, NTE, NTOP,
     %                      NBTNS1, NOTE1S, NOTETR, NS1, NS2, NS1, NS3,
     %                      0,      XYZP )

            IF( NS4 .EQ. NS3 ) THEN

C              ICI L'ARETE PERDUE NS1 NS3 EST LA "DIAGONALE"
C              DES 2 TETRAEDRES NTE et NTOP: 2T -> 3T?
C              ---------------------------------------------
C              SI L'ARETE NS1-NS3 PERFORE LA FACE COMMUNE AUX
C              2 TETRAEDRES ALORS ECHANGER 2 TETRAEDRES => 3 TETRAEDRES
C              D'ARETE COMMUNE NS1-NS3 QUI APPARTIENT A LA TETRAEDRISATION
               CALL CH2T3T( INFACO, MXFACO, LEFACO, 0,
     %                      IVOLTE, NVOLTE, NTE,    NFTE,   NTOP,
     %                      PTXYZD, N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                      NTNOUV, IERR )
C              IERR=1 SI PAS DE TETRAEDRE DE L'AUTRE COTE DE LA FACE
C                   2 SI SATURATION DU TABLEAU NOTETR
C                   3 SI LE VOLUME de NTE+NTOP EST DIFFERENT DU VOLUME
C                        DES 3 AUTRES TETRAEDRES
C                   5 SI LA FACE EST DANS LEFACO ET POUR LA CONSERVER
C                        PAS D'ECHANGE EFFECTUE
               IF( IERR .EQ. 2 ) THEN
C                 SATURATION DU TABLEAU NOTETR DES TETRAEDRES
                  GOTO 9900
               ENDIF

               IF( IERR .EQ. 0 ) THEN
C                 2T=>3T A ETE CORRECTEMENT EXECUTE
C                 L'ARETE NS1-NS3 A ETE CREEE DANS LA TETRAEDRISATION
C                 ---------------------------------------------------
ccc                  PRINT*,'rearfape: METHODE=',METHOD,' ARETE NS1=',NS1,
ccc     %                   '- NS3=',NS3,' VUE dans la TETRAEDRISATION'
C                 LE NO NOTETR DES 3 TETRAEDRES D'ARETE NS1-NS3
C                 EST STOCKE DANS NTNOUV(1:3)

C                 LE SOMMET NS3 EST DEPILE
                  LHPILE = LHPILE - 1
                  GOTO 10
               ENDIF

               IERR = 0
               GOTO 30

            ENDIF

C           NS3 N'EST PAS LE SOMMET OPPOSE NS4 A LA FACE COMMUNE NFTE DE NTE
C           IL FAUT FAIRE APPARAITRE L'ARETE NS1-NS4 DANS LA TETRAEDRISATION
            IF( LHPILE .GE. MXPILE ) THEN
               PRINT*,'rearfape: TABLEAU NUPILE SATURE MXPILE=',MXPILE
               IERR = 3
               GOTO 9999
            ENDIF
            LHPILE = LHPILE + 1
            NUPILE( LHPILE ) = NS4
            GOTO 10

 30      ENDDO

C        PROBLEME: L'ARETE NS1-NS3 EST ELLE DANS LA TETRAEDRISATION?
         CALL TETR1A( NS1,    NS3,    N1TETS, NOTETR,
     %                NBTE1A, MXTE1S, NOTE1S, IER )

 40      IF( NBTE1A .LE. 0 ) THEN
C           OUI: L'ARETE NS1-NS3 NON RETROUVEE DANS LA TETRAEDRISATION
ccc            PRINT*,'rearfape: ARETE PERDUE NS1=',NS1,' NS2=',NS2,
ccc     %             ' PAS de FACE DE TETRAEDRE PERCEE PAR L''ARETE NS1=',
ccc     %              NS1,' NS3=',NS3,'! ...'

C           TRACE DES NBTNS1 TETRAEDRES DE SOMMET NS1
            KTITRE='rearfape: PB PAS de FACE PERCEE par l''ARETE        
     %-          '
            WRITE(KTITRE(45:50),'(I6)') NS1
            WRITE(KTITRE(55:60),'(I6)') NS3
            CALL SANSDBL(  KTITRE, NC )
            CALL TRATE2AR( KTITRE(1:NC), PTXYZD, NBTNS1, NOTE1S, NOTETR,
     %                     NS1, NS2, NS1, NS3 )
         ENDIF


C        SECONDE METHODE POUR RECUPERER L'ARETE NAFACO DE LA FACE PERDUE NFLPER
C        ======================================================================
C        RECHERCHE DE TETRAEDRE EN TETRAEDRES INTERSECTES PAR L'ARETE NS1-NS2
C        AVEC ESSAI POUR 2 TETRAEDRES SUIVANTS DE 2T -> 3T
 52      METHOD = 2
         NBPASS55 = 0
ccc         print*
ccc         PRINT*,'rearfape: Methode 2 sur l''ARETE NS1=',NS1,' NS2=',NS2,
ccc     %          ' Nombre CHANGEMENTS NS1<->NS2 =',NBPENS12

C        NOMBRE DE PASSAGES AU DEBUT DU CALCUL DES INTERSECTIONS DE NS1-NS2
 55      NBPASS55 = NBPASS55 + 1
C        RECHERCHE DES NBTNS1 TETRAEDRES DE SOMMET NS1
         CALL TETR1S( NS1,    N1TETS, NOTETR,
     %                NBTNS1, MXTE1S, NOTE1S, IERR )
         IF( NBTNS1 .LE. 0 ) THEN
            PRINT*,'rearfape: Bizarre AUCUN TETRAEDRE de SOMMET',NS1
            GOTO 100
         ENDIF

C        PARCOURS DES TETRAEDRES DE SOMMET NS1 POUR EN TROUVER
C        UN DE FACE INTERSECTEE PAR L'ARETE NS1-NS2
C        -----------------------------------------------------
         NBXYZP = 0
         NB2T3T = 0
         DO 75 M = 1, NBTNS1

C           LE TETRAEDRE M DE SOMMET NS1
            NTE = NOTE1S( M )
            IF( NOTETR(1,NTE) .LE. 0 ) GOTO 75

C           RECHERCHE DU TETRAEDRE NTOP OPPOSE A LA FACE OPPOSEE
C           AU SOMMET NS1 DE NTE: RECHERCHE DU NO DE SOMMET NS1 DANS NTE
            DO K=1,4
               IF( NS1 .EQ. NOTETR(K,NTE) ) GOTO 60
            ENDDO
            GOTO 75

C           LE NO DE LA FACE OPPOSEE AU SOMMET K EST K MOD 4 + 1
 60         IF( K .NE. 4 ) THEN
               NFTE=K+1
            ELSE
               NFTE=1
            ENDIF

C           LA FACE NFTE DU TETRAEDRE NTE EST ELLE PERFOREE
C           PAR L'ARETE NS1-NS2?
C           NUMERO DES 3 SOMMETS DE LA FACE NFTE
            NOSOTR(1) = NOTETR( NOSOFATE(1,NFTE), NTE )
            NOSOTR(2) = NOTETR( NOSOFATE(2,NFTE), NTE )
            NOSOTR(3) = NOTETR( NOSOFATE(3,NFTE), NTE )

C           LE TETRAEDRE NTE A T IL POUR ARETE NS1-NS2?
            DO I=1,3
               IF( NOSOTR(I) .EQ. NS2 ) THEN
C                 OUI: LE TETRAEDRE NTE A POUR ARETE NS1-NS2
                  GOTO 99
               ENDIF
            ENDDO

C           INTERSECTION OU NON DE NS1-NS2 AVEC LA FACE NFTE DE NTE?
            CALL VOLU2T3T( NS1, NS2, NOSOTR, PTXYZD, NONOUI )
C           NONOUI: -1 NS1 EST UN SOMMET DU TRIANGLE NOSOTR
C                   -2 NS2 EST UN SOMMET DU TRIANGLE NOSOTR
C                    0 SI NS1-NS2 N'INTERSECTE PAS LE TRIANGLE NOSOTR
C                    1 NS1-NS2 INTERSECTE LE TRIANGLE FERME (ARETES COMPRISES)
C                      car LE VOLUME 2T = VOLUME 3T en VALEUR ABSOLUE
C            verifier que NONOUI et LINTER SONT COHERENTS...

C           CETTE FACE NOSOTR EST ELLE INTERSECTEE PAR L'ARETE NS1-NS2?
            CALL INARTR( PTXYZD(1,NS1), PTXYZD(1,NS2),
     %                   PTXYZD( 1, NOSOTR(1) ),
     %                   PTXYZD( 1, NOSOTR(2) ),
     %                   PTXYZD( 1, NOSOTR(3) ),
     %                   LINTER, XYZP(1,1), CBTR )
C           LINTER : -2 SI S1=S2
C           -1 SI S1-S2 PARALLELE AU PLAN DU TRIANGLE NOSOTR
C            0 SI S1-S2 N'INTERSECTE PAS LE TRIANGLE ENTRE CES 3 SOMMETS
C            1 SI S1-S2   INTERSECTE     LE TRIANGLE ENTRE CES 3 SOMMETS
C              XYZP: 3 COORDONNEES DU POINT D'INTERSECTION SI LINTER=1
C              CBTR: 3 COORDONNEES BARYCENTRIQUES DE XYZP1 DANS LE TRIANGLE
            IF( LINTER .LT. 0 ) THEN
C              POUR LE TEST DE VERIFICATION QUI SUIT... A SUPPRIMER ENSUITE
               LINTER = 0
            ENDIF
            IF( NONOUI .NE. LINTER ) THEN
              print*,'rearfape: verif NONOUI=',NONOUI,' LINTER=',LINTER,
     %               ' NS1=',NS1,' NS2=',NS2,' NOSOTR=',NOSOTR
            ENDIF

            IF( LINTER .LE. 0 ) THEN
C              PAS D'INTERSECTION DE NS1-NS2 AVEC LE TRIANGLE NOSOTR
C              PASSAGE AU TETRAEDRE DE SOMMET NS1 SUIVANT
               GOTO 75
            ENDIF

C           LINTER=1: LA FACE NFE DU TETRAEDRE NTE EST INTERSECTEE
C           PAR L'ARETE NS1-NS2 AU POINT XYZP(1)
C           ------------------------------------------------------
            NBXYZP = 1
            NUPILE( NBXYZP ) = NTE
ccc            print*
ccc            print*,'rearfape: arete',ns1,ns2,' 1-er Pt INTERSECTION',
ccc     %             (xyzp(kk,NBXYZP),kk=1,3)
ccc           print*,'rearfape: tetraedre',nte,':',(notetr(kk,nte),kk=1,8)
ccc            print*,'rearfape: face',nfte,':',nosotr,' cbtr=',cbtr

C           NTOP LE TETRAEDRE OPPOSE A NS1 ET LE NUMERO NFOP DE LA FACE
C           QUI EST LA FACE NFTE DU TETRAEDRE NTE PERCEE PAR L'ARETE NS1-NS2
            CALL NOFAOP( NFTE, NTE, NOTETR, NFOP, NTOP )
            IF( NTOP .EQ. 0 .OR. NOTETR(1,NTOP) .EQ. 0 ) THEN
C              PROBLEME LE TETRAEDRE NTOP EST DESACTIVE...
               PRINT*,'rearfape: PROBLEME: le TETRAEDRE(',NTOP,'):',
     %                 (NOTETR(kk,NTOP),kk=1,8),' EST DESACTIVE...'
               GOTO 100
            ENDIF

            NTOP0 = NTE
C           LE SOMMET OPPOSE NSOP EST NFOP - 1
            IF( NFOP .NE. 1 ) THEN
               NSOP = NFOP - 1
            ELSE
               NSOP = 4
            ENDIF
C           LE NO NS3 DU SOMMET DU TETRAEDRE NTOP OPPOSE A SA FACE NFOP
            NS3 = NOTETR( NSOP, NTOP )

C           TRACE DES NBTNS1 TETRAEDRES DE SOMMET NS1, NTOP ET L'ARETE
            IF( NBXYZP .GE. 32 ) THEN
               tracte = .true.
            ENDIF
            KTITRE='rearfape:        tetraedres de SOMMET NS1=      . AR
     %ETE NS1-NS2=       Nb XYZP=       '
            WRITE(KTITRE(11:16),'(I6)') NBTNS1
            WRITE(KTITRE(43:48),'(I6)') NS1
            WRITE(KTITRE(65:70),'(I6)') NS2
            WRITE(KTITRE(80:85),'(I6)') NBXYZP
            CALL SANSDBL(   KTITRE, NC )
            CALL TRATES2AR( KTITRE(1:NC), PTXYZD, NTE, NTOP,
     %                      NBTNS1, NOTE1S, NOTETR, NS1, NS2, 0, 0,
     %                      NBXYZP, XYZP )

C           RECHERCHE DE L'AUTRE POINT d'INTERSECTION de l'ARETE NS1-NS2
C           AVEC UNE AUTRE FACE DU TETRAEDRE NTOP
C           ------------------------------------------------------------
 65         DO 70 NFOP2 = 1, 4

               IF( NFOP2 .EQ. NFOP ) GOTO 70
C              NUMERO DES 3 SOMMETS DE LA FACE NFOP2
               NOSOTR(1) = NOTETR( NOSOFATE(1,NFOP2), NTOP )
               NOSOTR(2) = NOTETR( NOSOFATE(2,NFOP2), NTOP )
               NOSOTR(3) = NOTETR( NOSOFATE(3,NFOP2), NTOP )

C              CETTE FACE NOSOTR EST ELLE INTERSECTEE PAR L'ARETE NS1-NS2?
               CALL INARTR( PTXYZD(1,NS1), PTXYZD(1,NS2),
     %                      PTXYZD( 1, NOSOTR(1) ),
     %                      PTXYZD( 1, NOSOTR(2) ),
     %                      PTXYZD( 1, NOSOTR(3) ),
     %                      LINTER, XYZP(1,NBXYZP+1), CBTR )
C              LINTER : -2 SI S1=S2
C              -1 SI S1-S2 PARALLELE AU PLAN DU TRIANGLE NOSOTR
C               0 SI S1-S2 N'INTERSECTE PAS LE TRIANGLE ENTRE CES 3 SOMMETS
C               1 SI S1-S2   INTERSECTE     LE TRIANGLE ENTRE CES 3 SOMMETS
C                    XYZP: 3 COORDONNEES DU POINT D'INTERSECTION SI LINTER=1
C                    CBTR: 3 COORDONNEES BARYCENTRIQUES DE XYZP DANS
C                           LE TRIANGLE

               IF( LINTER .LE. 0 ) THEN
C                 PAS D'INTERSECTION DE NS1-NS2 AVEC LA FACE NFOP2 DE NTOP
C                 PASSAGE A LA FACE SUIVANTE DE NTOP
                  GOTO 70
               ENDIF

C              LINTER=1: LA FACE NFOP2 DE NTOP EST INTERSECTEE PAR
C                        L'ARETE NS1-NS2 AU POINT XYZP(NBXYZP)
ccc               print*
ccc           print*,'rearfape: arete',ns1,ns2,' Pt INTERSECTION',NBXYZP+1,
ccc     %              (xyzp(kk,NBXYZP+1),kk=1,3)
ccc         print*,'rearfape: tetraedre',ntop,':',(notetr(kk,ntop),kk=1,8)
ccc            print*,'rearfape: face',nfop2,':',nosotr,' cbtr=',cbtr

C              NTOPOP LE TETRAEDRE OPPOSE A LA FACE NFOP2 DU TETRAEDRE NTOP
C              ET LE NUMERO NFOPOP DE LA FACE
               CALL NOFAOP( NFOP2, NTOP, NOTETR, NFOPOP, NTOPOP )
               IF( NTOPOP .EQ. NTOP0 ) THEN
C                 RETROUVE LE TETRAEDRE PRECEDANT NTOP
                  GOTO 70
               ENDIF

               IF( NBXYZP .GE. MXXYZP ) THEN
C                 SATURATION DES TABLEAUX XYZP
                  GOTO 94
               ENDIF
C              UN POINT D'INTERSECTION DE PLUS
               NBXYZP = NBXYZP + 1
C              NO NOTETR DU TETRAEDRE DE SECONDE FACE INTERSECTEE PAR NS1-NS2
               NUPILE( NBXYZP ) = NTOP

C              LE SOMMET OPPOSE NSOPOP EST NFOPOP - 1
               IF( NFOPOP .NE. 1 ) THEN
                  NSOPOP = NFOPOP - 1
               ELSE
                  NSOPOP = 4
               ENDIF

C              LE NO NS4 DU SOMMET DU TETRAEDRE NTOPOP OPPOSE A SA FACE NFOPOP
               NS4 = NOTETR( NSOPOP, NTOPOP )

C              TRACE DES NBXYZP TETRAEDRES INTERSECTION AVEC L'ARETE NS1-NS2
               KTITRE='rearfape:        tetraedres INTERSECTION avec l''
     %ARETE NS1=        NS2=       '
               WRITE(KTITRE(11:16),'(I6)') NBXYZP
               WRITE(KTITRE(60:65),'(I6)') NS1
               WRITE(KTITRE(72:77),'(I6)') NS2
               CALL SANSDBL(   KTITRE, NC )
               CALL TRATES2AR( KTITRE(1:NC), PTXYZD, NTOP, NTOPOP,
     %                         NBXYZP, NUPILE, NOTETR, NS1,NS2, NS3,NS4,
     %                         NBXYZP, XYZP )

C              ESSAI D'ECHANGER les 2 TETRAEDRES NTOP NTOPOP => EN LES
C              3 TETRAEDRES D'ARETE COMMUNE NS3-NS4
               CALL CH2T3T( INFACO, MXFACO, LEFACO, 0,
     %                      IVOLTE, NVOLTE, NTOP,   NFOP2,  NTOPOP,
     %                      PTXYZD, N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                      NTNOUV, IERR )
C              IERR=1 SI PAS DE TETRAEDRE DE L'AUTRE COTE DE LA FACE NFOP2
C                   2 SI SATURATION DU TABLEAU NOTETR
C                   3 SI LE VOLUME de NTE+NTOP EST DIFFERENT DU VOLUME
C                        DES 3 AUTRES TETRAEDRES
C                   5 SI LA FACE EST DANS LEFACO ET POUR LA CONSERVER
C                        PAS D'ECHANGE EFFECTUE

               IF( IERR .EQ. 2 ) THEN
C                 SATURATION DU TABLEAU NOTETR DES TETRAEDRES
                  GOTO 9900
               ENDIF

               IF( IERR .EQ. 0 ) THEN
C                 NTOP+NTOPOP 2T => 3T d'ARETE COMMUNE NS3-NS4 A ETE EFFECTUEE
C                 L'ARETE NS3-NS4 A ETE CREEE DANS LA TETRAEDRISATION
                  NB2T3T = NB2T3T + 1
ccc                  PRINT*,'rearfape: METHODE=',METHOD,' ARETE NS3=',NS3,
ccc     %                   '- NS4=',NS4,' MISE dans TETRAEDRES par 2T->3T'
C                 LES TETRAEDRES NTOP NTOPOP SONT DEVENUS NTNOUV(1:3)
C                 RETOUR AU SOMMET NS1
                  GOTO 55
               ENDIF

               IF( NS4 .EQ. NS2 ) THEN

C                 L'ARETE NS1-NS2 EST ELLE DANS LA TETRAEDRISATION?
                  CALL TETR1A( NS1,    NS2,    N1TETS, NOTETR,
     %                         NBTE1A, MXTE1S, NOTE1S, IER )
                  IF( NBTE1A .GT. 0 ) THEN
C                    OUI: ARETE NS1-NS2 RETROUVEE DANS LA TETRAEDRISATION
                     GOTO 99
                  ELSE
C                    NON: UN CHANGEMENT 2T->3T A T IL ETE REALISE?
                     IF( NB2T3T .EQ. 0 ) THEN
C                       NON: PAS DE CHANGEMENT 2T->3T
                        IF( NBPENS12 .EQ. 0 ) THEN
C                          PERMUTATION de NS1-NS2 et REDEMARRAGE
                           NBPENS12 = NBPENS12 + 1
                           K   = NS1
                           NS1 = NS2
                           NS2 = K
                           GOTO 52
                        ELSE
C                          UN TETRAEDRE FINAL D'INTERSECTION MENANT A NS2
                           IF( NBXYZP .GE. MXXYZP ) THEN
C                             SATURATION DES TABLEAUX XYZP
                              GOTO 94
                           ENDIF
                           NBXYZP = NBXYZP + 1
C                          NO NOTETR DU TETRAEDRE DE SOMMET NS2
                           NUPILE( NBXYZP ) = NTOPOP
C                          PASSAGE A LA 3-EME METHODE POUR RECUPERER
C                          L'ARETE NS1-NS2
                           GOTO 80
                        ENDIF
                     ENDIF

C                    OUI: UN CHANGEMENT 2T->3T A ETE REALISE
ccc                     PRINT*,'rearfape: ARETE NS1=',NS1,'-',NS2,
ccc     %                      ' Nb de PASSAGES en NS1',NBPASS55
C                    REDEMARRAGE A PARTIR du SOMMET NS1
                     GOTO 55
                  ENDIF

               ENDIF

C              POURSUITE DES INTERSECTIONS DE L'ARETE NS1-NS2
               NTOP0 = NTOP
               NTOP  = NTOPOP
               NFOP  = NFOPOP
               GOTO 65

 70         ENDDO

ccc            PRINT*,'rearfate: PAS de SECOND POINT d''INTERSECTION du TET
ccc     %RAEDRE',NTOP,(NOTETR(kk,NTOP),kk=1,8)

 75      ENDDO


C        TROISIEME METHODE POUR RECUPERER L'ARETE NAFACO DE LA FACE PERDUE NFLPER
C        A PARTIR DES TETRAEDRES INTERSECTES PAR L'ARETE NS1-NS2
C        ========================================================================
 80      METHOD = 3
         NBTECR = 0
ccc         print*
ccc         PRINT*,'rearfape: Methode 3 sur l''ARETE NS1=',NS1,' NS2=',NS2,
ccc     %          ' Nombre de TETRAEDRES INTERSECTES=',NBXYZP
ccc         DO K=1,NBXYZP
ccc            NTE = NUPILE(K)
ccc          PRINT*,'rearfape: TETRAEDRE(',NTE,'):',(NOTETR(kk,NTE),kk=1,8)
ccc         ENDDO

C        TRACE DES NBXYZP TETRAEDRES INTERSECTION AVEC L'ARETE NS1-NS2
         KTITRE='rearfape:        tetraedres INTERSECTION avec L''ARETE 
     %NS1=         NS2=       '
         WRITE(KTITRE(11:16),'(I6)') NBXYZP
         WRITE(KTITRE(59:64),'(I6)') NS1
         WRITE(KTITRE(73:78),'(I6)') NS2
         CALL SANSDBL( KTITRE, NC )
         CALL TRATES2AR( KTITRE(1:NC), PTXYZD, NUPILE(1), NUPILE(NBXYZP)
     %                  ,NBXYZP, NUPILE, NOTETR, NS1, NS2, 0, 0,
     %                   NBXYZP-1, XYZP )

C        VOLUME DES NBXYZP TETRAEDRES INTERSECTES PAR L'ARETE NS1-NS2
         VOLET0 = 0D0
         DO K = 1, NBXYZP
            NTE = NUPILE( K )
            V = VOLTET( PTXYZD( 1, NOTETR(1,NTE) ),
     %                  PTXYZD( 1, NOTETR(2,NTE) ),
     %                  PTXYZD( 1, NOTETR(3,NTE) ),
     %                  PTXYZD( 1, NOTETR(4,NTE) ) )
            VOLET0 = VOLET0 + V
         ENDDO

C        CONSTRUCTION DE L'ETOILE DES FACES SIMPLES DES NBXYZP TETRAEDRES
C        EN VERSION 1
         CALL CRFETOI1( NBXYZP, NUPILE, NOTETR,
     %                  MXFETO, N1FEOC, N1FEVI, NFETOI )
C        N1FEOC: POINTEUR SUR LA PREMIERE FACE DE L'ETOILE
C                CHAINAGE SUIVANT DANS NFETOI(5,*)
C        N1FEVI: POINTEUR SUR LA PREMIERE FACE VIDE DE L'ETOILE
C                CHAINAGE SUIVANT DANS NFETOI(5,*)
C                =-1 SI SATURATION DU TABLEAU NFETOI
C        NFETOI: VERSION 1
C          1: NUMERO DU TETRAEDRE DANS NOTETR DE LA FACE DE L'ETOILE
C          2: NUMERO LOCAL AU TETRAEDRE DE LA FACE DE L'ETOILE
C             UN SIGNE NEGATIF INDIQUE UN TRAITEMENT EFFECTUE
C          3: NON UTILISE ICI
C          4: NON UTILISE ICI
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES

         CALL NFETOI12( NOTETR, N1FEOC, NFETOI )
C        EN SORTIE NFETOI EN VERSION 2
C          1: NUMERO NOTETR DU TETRAEDRE EXTERIEUR OPPOSE A LA FACE
C          2: NUMERO PTXYZD DU 1-ER SOMMET DE LA FACE
C          3: NUMERO PTXYZD DU 2-ME SOMMET DE LA FACE
C          4: NUMERO PTXYZD DU 3-ME SOMMET DE LA FACE
C             S1S2xS1S3 EST LA NORMALE DIRIGEE VERS L'INTERIEUR DE
C             L'ETOILE
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES

C        PARCOURS DES FACES SIMPLES DE L'ETOILE DE L'ARETE NS1-NS2
C        POUR VERIFIER QUE LA FACE
C        SOIT A UN SOMMET NS1 ou NS2
C        SOIT LE VOLUME FACE-NS1>0 et LE VOLUME FACE-NS2>0
C        ---------------------------------------------------------
C        PREMIERE FACE SIMPLE DE L'ETOILE DANS NFETOI
         NF = N1FEOC

 82      IF( NF .GT. 0 ) THEN

C           LA FACE NF A T ELLE POUR SOMMET NS1?
            DO K=1,3
               IF( NFETOI( 1+K, NF ) .EQ. NS1 ) THEN
                  GOTO 83
               ENDIF
            ENDDO

C           LA FACE NF A T ELLE POUR SOMMET NS2?
            DO K=1,3
               IF( NFETOI( 1+K, NF ) .EQ. NS2 ) THEN
                  GOTO 83
               ENDIF
            ENDDO

C           VOLUME DU TETRAEDRE FACE+NS1
            V = VOLTET( PTXYZD( 1, NFETOI(2,NF) ),
     %                  PTXYZD( 1, NFETOI(3,NF) ),
     %                  PTXYZD( 1, NFETOI(4,NF) ),
     %                  PTXYZD( 1, NS1 ) )
            IF( V .LT. 0D0 ) THEN
               GOTO 94
            ENDIF

C           VOLUME DU TETRAEDRE FACE+NS2
            V = VOLTET( PTXYZD( 1, NFETOI(2,NF) ),
     %                  PTXYZD( 1, NFETOI(3,NF) ),
     %                  PTXYZD( 1, NFETOI(4,NF) ),
     %                  PTXYZD( 1, NS2 ) )
            IF( V .LT. 0D0 ) THEN
               GOTO 94
            ENDIF

C           LA FACE SIMPLE SUIVANTE DE L'ETOILE
 83         NF = NFETOI( 5, NF )
            GOTO 82

         ENDIF

C        TRACE DES FACES DE L'ETOILE DE L'ARETE NS1-NS2
ccc         tracte = .true.
         CALL TRFETO2( PTXYZD, NOTETR, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR )

C        LES FACES SE PRESENTENT BIEN PAR RAPPORT A L'ARETE NS1-NS2
C        CONSTRUCTION DES ARETES DOUBLES DE L'ETOILE ET DE SOMMETS
C        DIFFERENTS DE NS1 et NS2
C       -----------------------------------------------------------
         NBAR2E = 0

C        PREMIERE FACE SIMPLE DE L'ETOILE DANS NFETOI
         NF = N1FEOC

 84      IF( NF .GT. 0 ) THEN

C           PARCOURS DES 3 ARETES DE LA FACE
            DO 85 NA = 1, 3

C              LE NO DES 2 SOMMETS DE L'ARETE
               NS3 = NFETOI( 1+NA, NF )
               IF( NA .EQ. 3 ) THEN
                  N = 1
               ELSE
                  N = NA + 1
               ENDIF
               NS4 = NFETOI( 1+N, NF )

               IF( NS3 .EQ. NS1 .OR. NS3 .EQ. NS2 ) GOTO 85
               IF( NS4 .EQ. NS1 .OR. NS4 .EQ. NS2 ) GOTO 85

C              L'ARETE NS3<NS4 EST ELLE DEJA STOCKEE?
               IF( NS3 .GT. NS4 ) THEN
                  K   = NS3
                  NS3 = NS4
                  NS4 = K
               ENDIF
               DO K = 1, NBAR2E
                  IF( NS3.EQ.NOAR2E(1,K) .AND. NS4.EQ.NOAR2E(2,K) ) THEN
                     NOAR2E(3,K) = NOAR2E(3,K) + 1
                     GOTO 85
                  ENDIF
               ENDDO

C              L'ARETE NS3-NS4 NON RETROUVEE EST AJOUTEE
               NBAR2E = NBAR2E + 1
               NOAR2E( 1, NBAR2E ) = NS3
               NOAR2E( 2, NBAR2E ) = NS4
               NOAR2E( 3, NBAR2E ) = 1

 85         ENDDO
           
C           LA FACE SIMPLE SUIVANTE DE L'ETOILE
            NF = NFETOI( 5, NF )
            GOTO 84

         ENDIF

C        AFFICHAGE POUR VERIFICATION DU TABLEAU NOAR2E
ccc         print*,'rearfape: Methode 3 No des sommets des',NBAR2E,
ccc     %          ' ARETES DOUBLES des FACES de l''ETOILE'
ccc         DO K= 1, NBAR2E
ccc            print*,'rearfape: No ARETE 2',K,': NS',(NOAR2E(kk,K),kk=1,3)
ccc         ENDDO

C        CONSTRUCTION DES TETRAEDRES D'ARETE NS1-NS2 + ARETE 2 de l'ETOILE
C        -----------------------------------------------------------------
         NBTECR = 0
         VOLET1 = 0D0
         DO NA = 1, NBAR2E

            IF( N1TEVI .LE. 0 ) THEN
C              SATURATION DU TABLEAU NOTETR DES TETRAEDRES
               GOTO 9900
            ENDIF
C           LE NOUVEAU TETRAEDRE
            NTE = N1TEVI
C           LE NOUVEAU PREMIER TETRAEDRE VIDE
            N1TEVI = NOTETR( 5, N1TEVI )

C           MISE A JOUR DU DERNIER TETRAEDRE OCCUPE DANS NOTETR
            NUDTETR = MAX( NUDTETR, NTE )

C           UN TETRAEDRE CREE D'ARETE NS1-NS2 DE PLUS
            NBTECR = NBTECR + 1
            NOTE1S( NBTECR ) = NTE

C           LES 2 SOMMETS DE L'ARETE DOUBLE NA DE L'ETOILE
            NS3 = NOAR2E( 1, NA )
            NS4 = NOAR2E( 2, NA )

C           VOLUME DU TETRAEDRE NS1234
            V = VOLTET( PTXYZD( 1, NS1 ), PTXYZD( 1, NS2 ),
     %                  PTXYZD( 1, NS3 ), PTXYZD( 1, NS4 ) )

            IF( V .LT. 0D0 ) THEN
C              PERMUTATION DES SOMMETS NS3 NS4
               K   = NS3
               NS3 = NS4
               NS4 = K
            ENDIF

            VOLET1 = VOLET1 + ABS( V )

C           LES 4 SOMMETS DU TETRAEDRE NTE D'ARETES NS1-NS2 et NS3-NS4
            NOTETR( 1, NTE ) = NS1
            NOTETR( 2, NTE ) = NS2
            NOTETR( 3, NTE ) = NS3
            NOTETR( 4, NTE ) = NS4

C           LES TETRAEDRES OPPOSES AUX 4 FACES SONT INCONNUS
            NOTETR( 5, NTE ) = -1
            NOTETR( 6, NTE ) = -1
            NOTETR( 7, NTE ) = -1
            NOTETR( 8, NTE ) = -1

C           MISE A JOUR D'UN TETRAEDRE DU SOMMET N1TETS
            DO K=1,4
               N1TETS( NOTETR(K,NTE) ) = NTE
            ENDDO

         ENDDO

C        TRACE DES NBTECR TETRAEDRES CREES d'ARETE NS1-NS2
ccc         tracte = .true.
         KTITRE='rearfape:        TETRAEDRES NOUVEAUX d''ARETE NS1=     
     %   NS2=        '
         WRITE(KTITRE(11:16),'(I6)') NBTECR
         WRITE(KTITRE(51:56),'(I6)') NS1
         WRITE(KTITRE(63:68),'(I6)') NS2
         CALL SANSDBL( KTITRE, NC )
         CALL TRATES2AR( KTITRE(1:NC), PTXYZD, 0, 0,
     %                   NBTECR, NOTE1S, NOTETR, NS1, NS2, 0, 0,
     %                   NBXYZP-1, XYZP )


C        COMPARAISON DU VOLUME DES TETRAEDRES DE L'ETOILE DE L'ARETE NS1-NS2
         IF( ABS(VOLET1-VOLET0) .GT. 1D-5 * VOLET0 ) THEN
            PRINT*,'rearfape: TROP GRANDE DIFFERENCE des VOLUMES V0=',
     %              VOLET0,' V1=',VOLET1
            GOTO 92
         ENDIF


C        MISE A JOUR DES TETRAEDRES OPPOSES DES NBTECR TETRAEDRES CREES
C        PAR AJOUT A NBTEOP DES TETRAEDRES OPPOSES AUX FACES DE L'ETOILE
         NBTEOP = NBTECR
C        PREMIERE FACE SIMPLE DE L'ETOILE DANS NFETOI
         NF = N1FEOC
 87      IF( NF .GT. 0 ) THEN

C           LE TETRAEDRE OPPOSE A LA FACE NF
            NTOP = NFETOI( 1, NF )
            IF( NTOP .GT. 0 ) THEN
               NBTEOP = NBTEOP + 1
               NOTE1S( NBTEOP ) = NTOP
            ENDIF
           
C           LA FACE SIMPLE SUIVANTE DE L'ETOILE
            NF = NFETOI( 5, NF )
            GOTO 87

         ENDIF
C        SUPPRESSION DES NO de TETRAEDRES DOUBLONS dans NOTE1S
         CALL UNITABL( NOTE1S, NBTEOP )

         CALL MJOPTE( NBTEOP, NOTE1S, N1TETS, NOTETR, NUDTETR,
     %                N1TEVI, PTXYZD, NBFANR )
C        NBFANR : NOMBRE DE FACES DE TETRAEDRE OPPOSE NON RETROUVE
         IF( NBFANR .GT. 0 ) THEN
            PRINT*,'rearfape:',NBFANR,' FACES OPPOSES DES TETRAEDRES CRE
     %ES NON RETROUVEES'
            GOTO 92
         ENDIF


C        LES NBTECR TETRAEDRES RECOUVRENT CORRECTEMENT L'ETOILE INITIALE
C        D'ARETE NS1-NS2
C        ---------------------------------------------------------------
         DO K = 1, NBTECR
            NTE = NOTE1S( K )
C           AJOUTER EVENTUELLEMENT LES 4 FACES DU TETRAEDRE NTE
C           DANS LE TABLEAU LEFACO
            CALL AJTELEFA( NTE, NOTETR, INFACO, MXFACO, LEFACO )
         ENDDO

C        DESTRUCTION DES NBXYZP TETRAEDRES INTERSECTES PAR L'ARETE NS1-NS2
C        ILS REDEVIENNENT VIDES DANS NOTETR
         DO K = 1, NBXYZP
            NTE = NUPILE( K )
            DO L=1,8
               NOTETR( L, NTE ) = 0
            ENDDO
            NOTETR( 5, NTE ) = N1TEVI
C           NTE REDEVIENT LE PREMIER TETRAEDRE VIDE DANS NOTETR
            N1TEVI = NTE
         ENDDO
         GOTO 99


C        RESTAURATION DES TETRAEDRES INITIAUX DE L'ETOILE D'ARETE NS1-NS2
C        ----------------------------------------------------------------
 92      PRINT*,'rearfape: Methode 3 ECHOUEE. REGENERATION DES',NBXYZP,
     %          ' TETRAEDRES INITIAUX'
C        SUPPRESSION DES NBTECR CREES DANS NOTETR
         DO K = 1, NBTECR
            NTE = NOTE1S( K )
            DO L=1,8
               NOTETR( L, NTE ) = 0
            ENDDO
            NOTETR( 5, NTE ) = N1TEVI
C           NTE REDEVIENT LE PREMIER TETRAEDRE VIDE DANS NOTETR
            N1TEVI = NTE
         ENDDO

C        MISE A JOUR DE N1TETS DES SOMMETS DES TETRAEDRES INITIAUX
         DO K = 1, NBXYZP
            NTE = NUPILE( K )
C           MISE A JOUR D'UN TETRAEDRE DU SOMMET N1TETS
            DO N=1,4
               N1TETS( NOTETR(N,NTE) ) = NTE
            ENDDO
         ENDDO

C        MISE A JOUR DES TETRAEDRES OPPOSES DES NBXYZP TETRAEDRES INITIAUX
         CALL MJOPTE( NBXYZP, NUPILE, N1TETS, NOTETR, NUDTETR,
     %                N1TEVI, PTXYZD, NBFANR )
C        NBFANR : NOMBRE DE FACES DE TETRAEDRE OPPOSE NON RETROUVE
         IF( NBFANR .GT. 0 ) THEN
            PRINT*,'rearfape:',NBFANR,
     %           ' FACES OPPOSES DES TETRAEDRES INITIAUX NON REGENEREES'
         ENDIF


C        METHODE 4 POUR RECUPERER L'ARETE NAFACO DE LA FACE
C        PERDUE NFLPER A PARTIR DES TETRAEDRES INTERSECTES PAR L'ARETE NS1-NS2
C        =====================================================================
 94      METHOD = 4
ccc         tracte = .true.
         print*
         PRINT*,'rearfape: Methode 4  AJOUT du MILIEU de l''ARETE NS1=',
     %           NS1,' NS2=',NS2,' ARETE',NAFACO,
     %          ' de la FACE PERDUE LEFACO',NFLPER
C        TRACE DES NBXYZP TETRAEDRES INTERSECTION AVEC L'ARETE NS1-NS2
ccc         tracte = .true.
         KTITRE='rearfape:        tetraedres INTERSECTION avec L''ARETE 
     %NS1=         NS2=       . Ajout du MILIEU'
         WRITE(KTITRE(11:16),'(I6)') NBXYZP
         WRITE(KTITRE(59:64),'(I6)') NS1
         WRITE(KTITRE(73:78),'(I6)') NS2
         CALL SANSDBL( KTITRE, NC )
         CALL TRATES2AR( KTITRE(1:NC), PTXYZD, 0, 0,
     %                   NBXYZP, NUPILE, NOTETR, NS1, NS2, 0, 0,
     %                   NBXYZP-1, XYZP )
C        RETOUR VERS LA TETRAEDRISATION DE CE MILIEU D'ARETE
         IERR = 17
         GOTO 9000


C        POUR FINIR, L'ARETE NS1-NS2 EST ELLE DANS LA TETRAEDRISATION?
C        -------------------------------------------------------------
 99      CALL TETR1A( NS1,    NS2,    N1TETS, NOTETR,
     %                NBTE1A, MXTE1S, NOTE1S, IERR )
         IF( NBTE1A .LE. 0 ) THEN

C           NON: L'ARETE NS1-NS2 NON RETROUVEE DANS LA TETRAEDRISATION
C           ..........................................................
            PRINT*,'rearfape: FIN METHODE=',METHOD,' ARETE NS1=',NS1,
     %         ' NS2=',NS2,' N''EST PAS une ARETE de la TETRAEDRISATION'
C           RECHERCHE DES NBTNS1 TETRAEDRES DE SOMMET NS1
            CALL TETR1S( NS1,    N1TETS, NOTETR,
     %                   NBTNS1, MXTE1S, NOTE1S, IERR )
C           RECHERCHE DES NBTNS1 TETRAEDRES DE SOMMET NS2
            CALL TETR1S( NS2,    N1TETS, NOTETR,
     %                   NBTNS2, MXTE1S, NOTE1S(NBTNS1+1), IERR )
            NBTNS12 = NBTNS1 + NBTNS2
            CALL UNITABL( NOTE1S, NBTNS12 )

C           TRACE DES NBTNS1 TETRAEDRES DE SOMMET NS1 et NS2
        KTITRE='rearfape: ARETE        -           DANS AUCUN TETRAEDRE'
            WRITE(KTITRE(17:22),'(I6)') NS1
            WRITE(KTITRE(26:31),'(I6)') NS2
            CALL SANSDBL(  KTITRE, NC )
            CALL TRATES2AR( KTITRE(1:NC), PTXYZD, 0, 0,
     %                      NBTNS12, NOTE1S, NOTETR, NS1, NS2, 0, 0,
     %                      NBXYZP-1, XYZP )

         ELSE

C           OUI: L'ARETE NS1-NS2 EST RETROUVEE DANS LA TETRAEDRISATION
C           ..........................................................
            PRINT*,'rearfape: Fin METHODE=',METHOD,' ARETE NS1=',NS1,
     %             '- NS2=',NS2,' RETROUVEE dans',NBTE1A,' TETRAEDRES'
            print*

C           TRACE DES NBTE1A TETRAEDRES D'ARETE NS1-NS2
          KTITRE='rearfape:        tetraedres d''ARETE       -         '
            WRITE(KTITRE(11:16),'(I6)') NBTE1A
            WRITE(KTITRE(38:43),'(I6)') NS1
            WRITE(KTITRE(46:51),'(I6)') NS2
            CALL SANSDBL(  KTITRE, NC )
            CALL TRATES2AR( KTITRE(1:NC), PTXYZD, 0, 0, NBTE1A, NOTE1S,
     %                      NOTETR, NS1, NS2, 0, 0, NBXYZP-1, XYZP )

         ENDIF

 100  ENDDO


C     CALCUL DU NOMBRE FINAL D'ARETES DE LA FACE NFLPER
C     DANS LA TETRAEDRISATION
C     -------------------------------------------------
      PRINT*,'rearfape: la FACE PERDUE LEFACO',NFLPER,' St:',
     %       (LEFACO(I,NFLPER),I=1,3),' IERR=',IERR
      NBARTE1 = 0
      DO I = 1, 3
C        L'ARETE I DE LA FACE PERDUE NFLPER A POUR SOMMETS NS1 NS2
         IF( I .EQ. 3 ) THEN
            I2 = 1
         ELSE
            I2 = I+1
         ENDIF
         NS1 = LEFACO( I , NFLPER )
         NS2 = LEFACO( I2, NFLPER )
C        L'ARETE NS1-NS2 EST ELLE DANS LA TETRAEDRISATION?
         CALL TETR1A( NS1,    NS2,    N1TETS, NOTETR,
     %                NBTE1A, MXTE1S, NOTE1S, IER )
         IF( NBTE1A .GT. 0 ) THEN
C           OUI: L'ARETE NS1-NS2 EST DANS LA TETRAEDRISATION
            NBARTE1 = NBARTE1 + 1
            PRINT*,'rearfape: l''ARETE FINALE',NS1,NS2,
     %             ' EST dans',NBTE1A,' TETRAEDRES'
         ENDIF
      ENDDO

      PRINT*,'rearfape:',NBARTE0,
     %       ' ARETES INITIALES dans la TETRAEDRISATION'
      PRINT*,'rearfape:',NBARTE1,
     %       ' ARETES FINALES   dans la TETRAEDRISATION'
      NBARRE = NBARTE1 - NBARTE0
      IF( NBARRE .GT. 0 ) THEN
         PRINT*,'rearfape:',NBARRE,
     %          ' ARETES RECUPEREES de la FACE PERDUE LEFACO',NFLPER
      ENDIF
      PRINT*,'FIN______________________________________________________'
      PRINT*

C     NOMBRE TOTAL D'ARETES PERDUES RECUPEREES DE NOTRCF
      NBTARRE = NBTARRE + NBARRE

C     PASSAGE A LA FACE LEFACO PERDUE SUIVANTE de NOTRCF
 1000 ENDDO

 9000 PRINT*,'rearfape: AU FINAL RECUPERATION de',NBTARRE,' ARETES des',
     %        NBTRCF,' FACES PERDUES pour la FACE',NFLPER
      GOTO 9999


 9900 IERR = 2
      print*,'rearfape: SATURATION DU TABLEAU NOTETR MXTETR=',MXTETR


 9999 TRACTE = TRACTE0
      RETURN
      END
