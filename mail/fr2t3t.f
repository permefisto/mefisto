      SUBROUTINE FR2T3T( PTXYZD, N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                   MXFACO, LEFACO, MXTE1S, MNTE1S,
     %                   NBFAPE, NOFAPE, NBARRE, MXTEFA, NOTEFA,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECUPERATION DES ARETES PERDUES PAR 2T -> 3T MULTIPLES
C -----
C
C ENTREES:
C --------
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C MXFACO : NOMBRE MAXIMAL DE FACES DECLARABLES DANS LEFACO
C LEFACO : FACE=TRIANGLE DE LA PEAU OU DES INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1, VOLUME2 DE LA FACE
C          NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C
C          ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C          => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C          LEFACO(9,*)  -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C          HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C          LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C          NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C          SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C          NF  = LEFACO( 9, NF )  ...
C          LEFACO(10,*) PREMIERE FACE DANS LE HACHAGE
C          LEFACO(11,.) = NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE, 0 SINON
CCCC       LEFACO(12,.) = NO FACEOC DE 1 A NBFACES D'OC
C
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
C NBFAPE : NOMBRE DE FACES PERDUES DE LEFACO DU CONTOUR AVANT ET APRES
C          NON PRESENTES DANS LA TETRAEDRISATION
C NOFAPE : NUMERO DANS LEFACO DE LA FACE PERDUE
C MXTE1S : NOMBRE MAXIMAL DE MOTS DU TABLEAU MNTE1S
C MNTE1S : TABLEAU DES TETRAEDRES D'UN SOMMET
C NOTEFA : TABLEAU AUXILIAIRE DE MXTEFA ENTIERS
C
C SORTIE :
C --------
C NBARRE : NOMBRE D'ARETES RETROUVEES
C IERR   : 0 SI PAS D'ERREUR
C          1 SI RETOUR EN ARRIERE A PROGRAMMER
C          2 SI SATURATION DU TABLEAU NOTETR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC       JUIN 2002
C....................................................................012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  PTXYZD(1:4,1:*)
      INTEGER           NOTETR(8,*),
     %                  N1TETS(*),
     %                  LEFACO(11,0:MXFACO),
     %                  NOFAPE(1:NBFAPE),
     %                  MNTE1S(MXTE1S),
     %                  NOTEFA(MXTEFA),
     %                  NVOLTE(1)
      INTEGER           NTNO23(3),
     %                  NOSOTR(3),
     %                  NOSOFATE(3,4)
C     NO DES SOMMETS DES FACES POUR QUE VU DU SOMMET MANQUANT
C     LES SOMMETS SOIENT VUS DANS LE SENS DIRECT
      DATA              NOSOFATE/ 1,2,3, 2,4,3, 3,4,1, 4,2,1 /
C
C     LA BOUCLE SUR LES FACES PERDUES
C     ===============================
CCC      MNTE1S = 0
CCC      MXTE1S = 0
      NBARRE = 0
      DO 1000 NFP=1,NBFAPE
C
C        LE NUMERO LEFACO DE LA FACE PERDUE
         NFPE = NOFAPE(NFP)
         IF( NFPE .LE. 0 ) GOTO 1000
C
C        BOUCLE SUR LES 3 ARETES DE LA FACE PERDUE
C        -----------------------------------------
         DO 900 I=1,3
            IF( I .EQ. 3 ) THEN
               I1 = 1
            ELSE
               I1 = I+1
            ENDIF
C           LES NUMEROS DES 2 SOMMETS DE L'ARETE I DE LA FACE PERDUE NFPE
            NS1 = LEFACO(I ,NFPE)
            NS2 = LEFACO(I1,NFPE)
C
C           CETTE ARETE I EST ELLE DANS LA TETRAEDRISATION?
            CALL TETR1A( NS1,    NS2,    N1TETS, NOTETR,
     %                   NBTE1A, MXTEFA, NOTEFA, IER )
            IF( NBTE1A .LE. 0 ) THEN
C
C              L'ARETE N'EST PAS DANS LA TETRAEDRISATION
C              RECHERCHE DES NBF1 TETRAEDRES DE SOMMET NS1
               CALL TETR1S( NS1,  N1TETS, NOTETR,
     %                      NBF1, MXTE1S, MNTE1S, IERR )
               IF( IERR .NE. 0 ) GOTO 1000
C
C              RECHERCHE DU TETRAEDRE OPPOSE AU SOMMET NS1
               MN = 0
               DO 90 JJ=1,NBF1
                  MN = MN + 1
C                 LE TETRAEDRE DE SOMMET NS1
                  NT = MNTE1S(MN)
                  IF( NOTETR(1,NT) .LE. 0 ) GOTO 90
C
C                 TETRAEDRE NON TRAITE: RECHERCHE DU SOMMET NS1 DANS NT
                  DO 20 K=1,4
                     IF( NS1 .EQ. NOTETR(K,NT) ) GOTO 30
 20               CONTINUE
                  GOTO 90
C
C                 LA FACE OPPOSEE AU SOMMET K EST K MOD 4 + 1
 30               IF( K .NE. 4 ) THEN
                     NF=K+1
                  ELSE
                     NF=1
                  ENDIF
C
C                 LE TETRAEDRE OPPOSE A NS1 ET LE NUMERO DE LA FACE
                  CALL NOFAOP( NF, NT, NOTETR, NFOP, NTOP )
C                 LE SOMMET OPPOSE NSOP EST NFOP - 1
                  IF( NFOP .NE. 1 ) THEN
                     NSOP = NFOP - 1
                  ELSE
                     NSOP = 4
                  ENDIF
                  IF( NOTETR(NSOP,NTOP) .NE. NS2 ) GOTO 90
C
C                 ICI L'ARETE PERDUE NS1 NS2 EST LA "DIAGONALE" DES 2 TETRAEDRES
C                 ==============================================================
C                 L'ARETE NS1-NS2 PERFORE-T-ELLE LA FACE COMMUNE AUX 2 TETRAEDRE
C                 ALORS 2 TETRAEDRES => 3 TETRAEDRES D'ARETE COMMUNE NS1-NS2
                  CALL CH2T3T( 1,      MXFACO, LEFACO, 0, 0, NVOLTE,
     %                         NT,     NF,     NTOP,
     %                         PTXYZD, N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                         NTNO23, IERR )
C                 IERR = 1 SI PAS DE TETRAEDRE DE L'AUTRE COTE DE LA FACE
C                        2 SI SATURATION DU TABLEAU NOTETR
C                        3 SI LE VOLUME DE NT+NTOP EST DIFFERENT DU VOLUME
C                             DES 3 AUTRES TETRAEDRES
                  IF( IERR .EQ. 2 ) THEN
C                    SATURATION DES TETRAEDRES
                     GOTO 9000
                  ELSE IF( IERR .NE. 0 ) THEN
                     IERR = 0
                     GOTO 90
                  ELSE
C                    2T=>3T A ETE CORRECTEMENT EXECUTE
                     NBARRE = NBARRE + 1
                     GOTO 800
                  ENDIF
C
 90            CONTINUE
               GOTO 800
C
C              CHANGEMENT DES TETRAEDRES POUR FAIRE APPARAITRE L'ARETE NS1-NS2
C              ===============================================================
               NT = NT0
               NF = NF0
C
C              L'ARETE NS1-NS2 PERFORE LA FACE COMMUNE AUX 2 TETRAEDRES
C              2 TETRAEDRES => 3 TETRAEDRES D'ARETE COMMUNE NS1-NS2
               CALL CH2T3T( 1,      MXFACO, LEFACO, 0, 0, NVOLTE,
     %                      NT,     NF,     NTOP,
     %                      PTXYZD, N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                      NTNO23, IERR )
C              IERR = 1 SI PAS DE TETRAEDRE DE L'AUTRE COTE DE LA FACE
C                     2 SI SATURATION DU TABLEAU NOTETR
C                     3 SI LE VOLUME DE NT+NTOP EST DIFFERENT DU VOLUME
C                          DES 3 AUTRES TETRAEDRES
               IF( IERR .EQ. 2 ) THEN
C                 SATURATION DES TETRAEDRES
                  GOTO 9000
               ELSE IF( IERR .NE. 0 ) THEN
                  IF( NT .EQ. NT0 ) THEN
                     IERR = 0
                     GOTO 800
                  ELSE
C                    PROBLEME: 2T=>3T A DEJA ETE FAIT ET TILT AVANT NS2
                     WRITE(IMPRIM,*) 'fr2t3t: RETOUR EN ARRIERE A FAIRE'
                     IERR = 1
                     GOTO 9000
                  ENDIF
               ENDIF
C
C              ICI: 2T=>3T A ETE CORRECTEMENT EXECUTE
C              RECHERCHE PARMI LES 3 TETRAEDRES GENERES DU TETRAEDRE
C              DONT UNE FACE EST PERFOREE PAR L'ARETE NS1-NS2
               NF0 = 0
  120          DO 150 J=1,3
C                 LE NOUVEAU TETRAEDRE J
                  NT = NTNO23(J)
                  DO 140 NF=1,4
C                    NUMERO DES 3 SOMMETS DE LA FACE NF DE NT
                     NOSOTR(1) = NOTETR( NOSOFATE(1,NF), NT )
                     NOSOTR(2) = NOTETR( NOSOFATE(2,NF), NT )
                     NOSOTR(3) = NOTETR( NOSOFATE(3,NF), NT )
                     CALL VOLU2T3T( NS1, NS2, NOSOTR, PTXYZD, NONOUI )
                     IF( NONOUI .NE. 1 ) GOTO 140
C
C                    LE TETRAEDRE OPPOSE A LA FACE NF DU TETRAEDRE NT
                     CALL NOFAOP( NF, NT, NOTETR, NFOP, NTOP )
C                    LE SOMMET OPPOSE NSOP EST NFOP - 1
                     IF( NFOP .NE. 1 ) THEN
                        NSOP = NFOP - 1
                     ELSE
                        NSOP = 4
                     ENDIF
                     IF( NOTETR(NSOP,NTOP) .EQ. NS2 ) NF0 = 1
C
C                    L'ARETE NS1-NS2 PERFORE LA FACE NF DE NT
C                    2 TETRAEDRES => 3 TETRAEDRES
                     CALL CH2T3T( 1, MXFACO, LEFACO, 0, 0, NVOLTE,
     %                          NT,     NF,     NTOP,
     %                          PTXYZD, N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                          NTNO23, IERR )
                     IF( IERR .EQ. 2 ) THEN
C                       SATURATION DES TETRAEDRES
                        GOTO 9000
                     ELSE IF( IERR .NE. 0 ) THEN
                        IF( NT .EQ. NT0 ) THEN
                           IERR = 0
                           GOTO 800
                        ELSE
C                          PROBLEME: 2T=>3T A DEJA ETE FAIT ET TILT AVANT NS2
                           WRITE(IMPRIM,*)
     %                     'fr2t3t: RETOUR EN ARRIERE A PROGRAMMER'
                           IERR = 1
                           GOTO 9000
                        ENDIF
                     ENDIF
C                    NTOP EST IL LE DERNIER TETRAEDRE?
                     IF( NF0 .NE. 0 ) THEN
C                       OUI: SORTIE POUR CETTE ARETE ENTIEREMENT RETROUVEE
                        NBARRE = NBARRE + 1
                        GOTO 800
                     ELSE
C                       NON: PASSAGE AU TETRAEDRE SUIVANT
                        GOTO 120
                     ENDIF
 140              CONTINUE
 150           CONTINUE
               WRITE(IMPRIM,*)'fr2t3t: ANOMALIE PAS DE TETRAEDRE PERCE'
 800           CONTINUE
            ENDIF
 900     CONTINUE
 1000 CONTINUE

 9000 PRINT*,'fr2t3t  :',NBARRE,' ARETES de FACES PERDUES RETROUVEES'

C     TRACE DES FACES PERDUES
      CALL TRFAPE( NBFAPE, NOFAPE, MXFACO, LEFACO, PTXYZD )

      RETURN
      END
