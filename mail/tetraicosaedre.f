      SUBROUTINE TETRAICOSAEDRE( PTXYZD, MXTETR, N1TEVI, NOTETR,
     %                           N1TETS, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CREATION DES 20 PREMIERS TETRAEDRES RECOUVRANT L'ICOSAEDRE
C -----    ENGLOBANT ET INITIALISATION DES TETRAEDRES VIDES
C          LE CENTRE DE L'ICOSAEDRE EST LE 13-EME POINT DE PTXYZD
C
C ENTREE :
C --------
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE ALENTOUR
C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES DECLARABLES DANS NOTETR
C
C SORTIES:
C --------
C N1TEVI : NUMERO DU PREMIER TETRAEDRE VIDE DANS LE CHAINAGE NOTETR(5,*)
C NOTETR : SOMMET1   < SOMMET2    SOMMET3     SOMMET4
C          ORDRE CLASSIQUE 1 2 3 DIRECT EN BAS ET 4 EN HAUT
C          LE SOMMET 1 EST LE PLUS PETIT
C          LE VOLUME EST POSITIF => UN ORDRE
C
C          TETRAEDRE1  TETREDRE2  TETRAEDRE3  TETREDRE4  VOISINS
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C N1TETS : N1TETS(I) NUMERO NOTETR D'UN TETRAEDRE DE SOMMET I PTXYZD
C IERR   : 0 SI PAS D'ERREUR
C          1 TOUS LES POINTS SONT COPLANAIRES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 2006
C2345X7..............................................................012
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  PTXYZD(1:4,1:*)
      INTEGER           NOTETR(8,*),
     &                  N1TETS(*)
      DOUBLE PRECISION  VOLUTE, CENTRE(4)

ccc      INTEGER           NOSOFATE(3,4)
ccc      DATA              NOSOFATE/ 1,2,3, 2,4,3, 3,4,1, 4,2,1 /
C
C     COORDONNEES DES 12 SOMMETS DE L'ICOSAEDRE ET DU BARYCENTRE
      DO 1 I=1,13
         PRINT *,'ICOSAEDRE XYZ', I,'=',(PTXYZD(J,I),J=1,3)
 1    CONTINUE
C     LE POLE NORD EST LE SOMMET
C     LES 5 POINTS DU PLAN // AU DESSUS  DE L'EQUATEUR SONT 2 3 4 5 6
C     LES 5 POINTS DU PLAN // AU DESSOUS DE L'EQUATEUR SONT 7 8 9 10 11
C     LE POLE SUD EST LE SOMMET 12
C     LE CENTRE DE L'ICOSAEDRE EST LE SOMMET 13
C
C     LES TETRAEDRES 1 A 5 SUPERIEURS
      NS4 = 13
      N   = 0
      DO 10 I=1,5
C
C        LE NO DES 4 SOMMETS DU TETRAEDRE I
         NS1 = 1
         NS2 = I+2
         IF( NS2 .EQ. 7 ) NS2 = 2
         NS3 = I+1
         N = N + 1
         NOTETR(1,N) = NS1
         NOTETR(2,N) = NS2
         NOTETR(3,N) = NS3
         NOTETR(4,N) = NS4
C
C        CALCUL DU CENTRE DE LA BOULE CIRCONSCRITE ET DU CARRE DU RAYON
         CALL CERABOCI( PTXYZD(1,NS1), PTXYZD(1,NS2),
     %                  PTXYZD(1,NS3), PTXYZD(1,NS4),
     %                  CENTRE, VOLUTE, IERR )
         IF( IERR .NE. 0 ) RETURN

 10   CONTINUE
C
C     LES TETRAEDRES 6 A 10 D'ARETES SUR L'EQUATEUR SUPERIEUR
      DO 20 I=1,5
C
C        LE NO DES 4 SOMMETS DU TETRAEDRE
         NS1 = I+6
         NS2 = I+1
         NS3 = I+2
         IF( NS3 .EQ. 7 ) NS3 = 2
         N = N + 1
         NOTETR(1,N) = NS1
         NOTETR(2,N) = NS2
         NOTETR(3,N) = NS3
         NOTETR(4,N) = NS4
C
C        CALCUL DU CENTRE DE LA BOULE CIRCONSCRITE ET DU CARRE DU RAYON
         CALL CERABOCI( PTXYZD(1,NS1), PTXYZD(1,NS2),
     %                  PTXYZD(1,NS3), PTXYZD(1,NS4),
     %                  CENTRE, VOLUTE, IERR )
         IF( IERR .NE. 0 ) RETURN

 20   CONTINUE
C
C     LES TETRAEDRES 11 A 15 D'ARETES SUR L'EQUATEUR INFERIEUR
      DO 30 I=1,5
C
C        LE NO DES 4 SOMMETS DU TETRAEDRE
         NS1 = I+2
         IF( NS1 .EQ. 7 ) NS1 = 2
         NS3 = I+6
         NS2 = NS3 + 1
         IF( NS2 .EQ. 12 ) NS2 = 7
         N = N + 1
         NOTETR(1,N) = NS1
         NOTETR(2,N) = NS2
         NOTETR(3,N) = NS3
         NOTETR(4,N) = NS4
C
C        CALCUL DU CENTRE DE LA BOULE CIRCONSCRITE ET DU CARRE DU RAYON
         CALL CERABOCI( PTXYZD(1,NS1), PTXYZD(1,NS2),
     %                  PTXYZD(1,NS3), PTXYZD(1,NS4),
     %                  CENTRE, VOLUTE, IERR )
         IF( IERR .NE. 0 ) RETURN
 30   CONTINUE
C
C     LES TETRAEDRES 16 A 20 INFERIEURS
      DO 40 I=1,5
C
C        LE NO DES 4 SOMMETS DU TETRAEDRE
         NS1 = 12
         NS2 = I+6
         NS3 = NS2 + 1
         IF( NS3 .EQ. 12 ) NS3 = 7
         N = N + 1
         NOTETR(1,N) = NS1
         NOTETR(2,N) = NS2
         NOTETR(3,N) = NS3
         NOTETR(4,N) = NS4
C
C        CALCUL DU CENTRE DE LA BOULE CIRCONSCRITE ET DU CARRE DU RAYON
         CALL CERABOCI( PTXYZD(1,NS1), PTXYZD(1,NS2),
     %                  PTXYZD(1,NS3), PTXYZD(1,NS4),
     %                  CENTRE, VOLUTE, IERR )
         IF( IERR .NE. 0 ) RETURN

 40   CONTINUE
C
CCCC     MISE A JOUR DU NO DE TETRAEDRE DE L'AUTRE COTE DES FACES
CCCC     ========================================================
CCC      DO 45 N=1,20
CCC         DO 42 K=5,8
CCC            NOTETR(K,N) = 0
CCC 42      CONTINUE
CCC 45   CONTINUE
CCCC
CCC      DO 90 N=1,20
CCC         DO 80 K=1,4
CCC            NS1 = NOTETR( NOSOFATE(1,K), N )
CCC            NS2 = NOTETR( NOSOFATE(2,K), N )
CCC            NS3 = NOTETR( NOSOFATE(3,K), N )
CCCC           CETTE FACE DE SOMMETS NS1-NS2-NS3 EST ELLE DANS UN AUTRE TETRAEDR
CCC            DO 70 NN=N+1,20
CCC               DO 60 KK=1,4
CCC                  NNS1 = NOTETR( NOSOFATE(1,KK), NN )
CCC                  NNS2 = NOTETR( NOSOFATE(2,KK), NN )
CCC                  NNS3 = NOTETR( NOSOFATE(3,KK), NN )
CCC                  IF( NS1. EQ. NNS1 .OR.
CCC     %                NS1. EQ. NNS2 .OR.
CCC     %                NS1. EQ. NNS3 ) THEN
CCC                     IF( NS2. EQ. NNS1 .OR.
CCC     %                   NS2. EQ. NNS2 .OR.
CCC     %                   NS2. EQ. NNS3 ) THEN
CCC                        IF( NS3. EQ. NNS1 .OR.
CCC     %                      NS3. EQ. NNS2 .OR.
CCC     %                      NS3. EQ. NNS3 ) THEN
CCC                           NOTETR(4+K ,N ) = NN
CCC                           NOTETR(4+KK,NN) = N
CCC                        ENDIF
CCC                     ENDIF
CCC                  ENDIF
CCC 60            CONTINUE
CCC 70         CONTINUE
CCC 80      CONTINUE
CCC 90   CONTINUE
CCCC
CCC      WRITE(IMPRIM,*)
CCC10090 FORMAT('      NOTETR(',I1,',',I2,') = ',I2)
CCC10091 FORMAT('C')
CCC      DO 96 N=1,20
CCC         DO 94 K=5,8
CCC            WRITE(IMPRIM,10090) K,N,NOTETR(K,N)
CCC 94      CONTINUE
CCC         WRITE(IMPRIM,10091)
CCC 96   CONTINUE
C
C     RESULTAT DES LIGNES CCC QUI PRECEDENT...
C     ========================================
      NOTETR(5, 1) =  0
      NOTETR(6, 1) =  6
      NOTETR(7, 1) =  5
      NOTETR(8, 1) =  2
C
      NOTETR(5, 2) =  0
      NOTETR(6, 2) =  7
      NOTETR(7, 2) =  1
      NOTETR(8, 2) =  3
C
      NOTETR(5, 3) =  0
      NOTETR(6, 3) =  8
      NOTETR(7, 3) =  2
      NOTETR(8, 3) =  4
C
      NOTETR(5, 4) =  0
      NOTETR(6, 4) =  9
      NOTETR(7, 4) =  3
      NOTETR(8, 4) =  5
C
      NOTETR(5, 5) =  0
      NOTETR(6, 5) = 10
      NOTETR(7, 5) =  4
      NOTETR(8, 5) =  1
C
      NOTETR(5, 6) =  0
      NOTETR(6, 6) =  1
      NOTETR(7, 6) = 11
      NOTETR(8, 6) = 15
C
      NOTETR(5, 7) =  0
      NOTETR(6, 7) =  2
      NOTETR(7, 7) = 12
      NOTETR(8, 7) = 11
C
      NOTETR(5, 8) =  0
      NOTETR(6, 8) =  3
      NOTETR(7, 8) = 13
      NOTETR(8, 8) = 12
C
      NOTETR(5, 9) =  0
      NOTETR(6, 9) =  4
      NOTETR(7, 9) = 14
      NOTETR(8, 9) = 13
C
      NOTETR(5,10) =  0
      NOTETR(6,10) =  5
      NOTETR(7,10) = 15
      NOTETR(8,10) = 14
C
      NOTETR(5,11) =  0
      NOTETR(6,11) = 16
      NOTETR(7,11) =  6
      NOTETR(8,11) =  7
C
      NOTETR(5,12) =  0
      NOTETR(6,12) = 17
      NOTETR(7,12) =  7
      NOTETR(8,12) =  8
C
      NOTETR(5,13) =  0
      NOTETR(6,13) = 18
      NOTETR(7,13) =  8
      NOTETR(8,13) =  9
C
      NOTETR(5,14) =  0
      NOTETR(6,14) = 19
      NOTETR(7,14) =  9
      NOTETR(8,14) = 10
C
      NOTETR(5,15) =  0
      NOTETR(6,15) = 20
      NOTETR(7,15) = 10
      NOTETR(8,15) =  6
C
      NOTETR(5,16) =  0
      NOTETR(6,16) = 11
      NOTETR(7,16) = 17
      NOTETR(8,16) = 20
C
      NOTETR(5,17) =  0
      NOTETR(6,17) = 12
      NOTETR(7,17) = 18
      NOTETR(8,17) = 16
C
      NOTETR(5,18) =  0
      NOTETR(6,18) = 13
      NOTETR(7,18) = 19
      NOTETR(8,18) = 17
C
      NOTETR(5,19) =  0
      NOTETR(6,19) = 14
      NOTETR(7,19) = 20
      NOTETR(8,19) = 18
C
      NOTETR(5,20) =  0
      NOTETR(6,20) = 15
      NOTETR(7,20) = 16
      NOTETR(8,20) = 19
C
C     DEBUT DES CHAINAGES DES TETRAEDRES OCCUPES ET VIDES
      N1TEVI = 21
C
C     CHAINAGE DES TETRAEDRES VIDES SUIVANTS
      DO 100 I=N1TEVI , MXTETR
C        LE NUMERO DU PREMIER SOMMET EST MARQUE TETRAEDRE VIDE
         NOTETR(1,I) = 0
C        LE TETRAEDRE VIDE SUIVANT
         NOTETR(5,I) = I + 1
 100  CONTINUE
C     LA FIN DU CHAINAGE DES TETRAEDRES SUIVANTS VIDES
      NOTETR(5,MXTETR) = 0
C
C     UN NUMERO DE TETRAEDRE POUR CHACUN DES 13 SOMMETS DE L'ICOSAEDRE
      N1TETS( 1) = 1
      N1TETS( 2) = 6
      N1TETS( 3) = 2
      N1TETS( 4) = 3
      N1TETS( 5) = 4
      N1TETS( 6) = 5
      N1TETS( 7) = 16
      N1TETS( 8) = 17
      N1TETS( 9) = 18
      N1TETS(10) = 19
      N1TETS(11) = 15
      N1TETS(12) = 20
      N1TETS(13) = 18
C     LE TETRAEDRE 18 CONTIENT LE MIN DE L'HEXAEDRE MIN-MAX
C
      RETURN
      END
