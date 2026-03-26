      SUBROUTINE BA2T3T2T( PTXYZD, MXFACO, LEFACO, IVOLTE, NVOLTE,
     %                     NOTETR, N1TEVI, NUDTETR, N1TETS,
     %                     MXPILE, LHPILE, LAPILE,
     %                     NBEC2T3T, NBEC3T2T )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AMELIORER LA QUALITE DES LHPILE FACES DE TETRAEDRES
C -----    PAR ECHANGE 2T->3T ou 3T->2T des TETRAEDRES

C ENTREES:
C --------
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
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

C IVOLTE : 0 PAS       DE TABLEAU NVOLTE A L'APPEL
C          1 EXISTENCE DU TABLEAU NVOLTE A L'APPEL
C NVOLTE : NUMERO DU VOLUME (1 A NBVOPA) DE CHAQUE TETRAEDRE
C          -1 SI VOLUME INCONNU

C MXPILE : NOMBRE MAXIMAL DE FACES DE TETRAEDRES DU TABLEAU LAPILE
C LHPILE : NOMBRE DE FACES DE TETRAEDRES A TRAITER DANS LAPILE

C MODIFIES:
C ---------
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NOTETR
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE
C N1TETS : NUMERO DANS NOTETR D'UN TETRAEDRE DE CHAQUE SOMMET
C LAPILE : LAPILE(1,*) NUMERO DANS NOTETR DES NBFATE TETRAEDRES A TRAITER
C          LAPILE(2,*) NUMERO LOCAL 1 A 4 DE LA FACE A ESSAYER D'ECHANGER

C SORTIES:
C --------
C NBEC2T3T: NOMBRE D'ECHANGES DES TETRAEDRES 2T->3T
C NBEC3T2T: NOMBRE D'ECHANGES DES TETRAEDRES 3T->2T
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR  : ALAIN PERRONNET             St PIERRE du PERRAY Fevrier 2018
C23456...............................................................012
      INTEGER           NOTETR(8,*), N1TETS(*), LAPILE(2,MXPILE),
     %                  NTEOLD(3), NTENEW(3), NVOLTE(*)
      DOUBLE PRECISION  PTXYZD(4,*)
C     NO (DE 1 A 9) DES 3 ARETES DES 4 FACES DU TETRAEDRE
      INTEGER           NOARFATE(3,4)
      DATA              NOARFATE / 1,2,3, 2,5,6, 3,6,4, 5,1,4 /

C     ESSAI D'AMELIORER LES NBFATE TETRAEDRES PAR LEURS FACES
C     =======================================================
      NBEC2T3T = 0
      NBEC3T2T = 0

 10   IF( LHPILE .GT. 0 ) THEN

C        RECUPERATION DU TETRAEDRE ET DE SA FACE
         NTE = LAPILE( 1, LHPILE )
         NFE = LAPILE( 2, LHPILE )

         LHPILE = LHPILE - 1

         IF( NTE .LE. 0 .OR. NFE .LE. 0 ) GOTO 10
         IF( NOTETR(1,NTE) .LE. 0 ) GOTO 10

C        LA FACE NF DE NTE EST ELLE DANS LEFACO?
         CALL NULEFT( NFE, NTE, NOTETR, MXFACO, LEFACO, NFLE )
         IF( NFLE .GT. 0 ) THEN
C           OUI: LA FACE EST DANS LEFACO
C           CONSERVATION DE CETTE FACE LEFACO EN NE PASSANT PAS AU DELA
            GOTO 10
         ENDIF

C        TENTATIVE DE CHANGER 2T->3T SI MIN QUALITE EST MAXIMISE
C        -------------------------------------------------------
         CALL CH2T3T( 1,      MXFACO, LEFACO, 1,      IVOLTE, NVOLTE,
     %                NTE,    NFE,    NTOP,
     %                PTXYZD, N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                NTENEW, IER )

         IF( IER .EQ. 0 ) THEN

C           REUSSITE: 2T->3T
            NBEC2T3T = NBEC2T3T + 1

ccc         PRINT*,'ba2t3t2t: NTE=',NTE,' NFE=',NFE,
ccc                ' 2T=>3T AMELIORE l''ETOILE',
ccc     %           NTENEW,' Qualites=',Q(1),Q(2),' =>',Q(3),Q(4),Q(5)

C           NTENEW=NO NOTETR DES 3 TETRAEDRES QUI REMPLACENT NT et NTOP
C           LES 2 FACES EXTERIEURES 1 et 4 DE CES 3 TETRAEDRES SONT EMPILEES
C           CF ch2t3t.f
            IF( LHPILE+6 .GT. MXPILE ) GOTO 9990
            DO K = 1, 3

C              LE NOUVEAU TETRAEDRE K
               NT = NTENEW( K )

C              RECHERCHE DU TETRAEDRE OPPOSE NTOP ET DE SA FACE NFOP
C              A LA FACE 1 DU TETRAEDRE NT
               CALL NOFAOP( 1, NT, NOTETR,  NFOP, NTOP )
               IF( NFOP .GT. 0 ) THEN
                  LHPILE = LHPILE + 1
                  LAPILE(1,LHPILE) = NTOP
                  LAPILE(2,LHPILE) = NFOP
               ENDIF

C              RECHERCHE DU TETRAEDRE OPPOSE NTOP ET DE SA FACE NFOP
C              A LA FACE 4 DU TETRAEDRE NT
               CALL NOFAOP( 4, NT, NOTETR,  NFOP, NTOP )
               IF( NFOP .GT. 0 ) THEN
                  LHPILE = LHPILE + 1
                  LAPILE(1,LHPILE) = NTOP
                  LAPILE(2,LHPILE) = NFOP
               ENDIF

            ENDDO

C           RETOUR EN HAUT DE PILE
            GOTO 10

         ENDIF


C        TENTATIVE DE CHANGER 3T->2T SI MIN QUALITE EST MAXIMISE
C        -------------------------------------------------------
C        PARCOURS DES 3 ARETES DE LA FACE NFE DU TETRAEDRE NTE
         DO N = 1, 3

C           NO DE L'ARETE N DE LA FACE NFE DE NTE
            NAE = NOARFATE( N, NFE )
            CALL CH3T2T( 1,   MXFACO, LEFACO, 1, IVOLTE, NVOLTE,
     %                   NTE, NAE,  PTXYZD,
     %                   N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                   NTEOLD, NTENEW, IER )

            IF( IER .EQ. 0 ) THEN

C              REUSSITE: 3T->2T
               NBEC3T2T = NBEC3T2T + 1

ccc            PRINT*,'ba2t3t2t: 3T->2T AMELIORE l''ETOILE',
ccc     %              NTENEW(1), NTENEW(2),
ccc     %             ' Qualites=',Q(1),Q(2),Q(3),' ->',Q(4),Q(5)

C              NTENEW: NUMERO NOTETR DES 2 TETRAEDRES QUI REMPLACENT NTEOLD(1:3)
C              LES 3 FACES EXTERIEURES (2 3 4) DE CES 2 TETRAEDRES SONT EMPILEES
C              CF ch3t2t.f
               IF( LHPILE+6 .GT. MXPILE ) GOTO 9990
               DO K = 1, 2

C                 LE NOUVEAU TETRAEDRE K
                  NT = NTENEW( K )

C                 RECHERCHE DU TETRAEDRE OPPOSE NTOP ET DE SA FACE NFOP
C                 A LA FACE 2 DU TETRAEDRE NT
                  CALL NOFAOP( 2, NT, NOTETR,  NFOP, NTOP )
                  IF( NFOP .GT. 0 ) THEN
                     LHPILE = LHPILE + 1
                     LAPILE(1,LHPILE) = NTOP
                     LAPILE(2,LHPILE) = NFOP
                  ENDIF

C                 RECHERCHE DU TETRAEDRE OPPOSE NTOP ET DE SA FACE NFOP
C                 A LA FACE 3 DU TETRAEDRE NT
                  CALL NOFAOP( 3, NT, NOTETR,  NFOP, NTOP )
                  IF( NFOP .GT. 0 ) THEN
                     LHPILE = LHPILE + 1
                     LAPILE(1,LHPILE) = NTOP
                     LAPILE(2,LHPILE) = NFOP
                  ENDIF

C                 RECHERCHE DU TETRAEDRE OPPOSE NTOP ET DE SA FACE NFOP
C                 A LA FACE 4 DU TETRAEDRE NT
                  CALL NOFAOP( 4, NT, NOTETR,  NFOP, NTOP )
                  IF( NFOP .GT. 0 ) THEN
                     LHPILE = LHPILE + 1
                     LAPILE(1,LHPILE) = NTOP
                     LAPILE(2,LHPILE) = NFOP
                  ENDIF

               ENDDO

C              RETOUR EN HAUT DE PILE
               GOTO 10

            ENDIF

         ENDDO

C        RETOUR EN HAUT DE PILE
         GOTO 10

      ENDIF

      RETURN
 

 9990 PRINT*,'ba2t3t2t: SATURATION LAPILE MXPILE=',MXPILE,' A AUGMENTER'
      RETURN
      END
