      SUBROUTINE VD2F2F( COANPL, PTXYZD, N1TETS, NOTETR,
     %                   MXFACO, LEFACO, N1FASC,
     %                   MXTE1S, NOTE1S, NB2F2F )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ECHANGE DU QUADRILATERE COPLANAIRE FORME PAR 2 TRIANGLES
C -----    PERDUS AVEC LES 2 AUTRES TRIANGLES SI L'AUTRE DIAGONALE
C          SE TROUVE DANS LA TETRAEDRISATION

C ENTREES:
C --------
C COANPL : SEUIL DU COSINUS DE L'ANGLE FORME PAR LES NORMALES AUX
C          2 FACES ET AU DESSUS DUQUEL LES FACES
C          SONT CONSIDEREES COPLANAIRES
C          ( 0.9848   => 10   DEGRES )
C          ( 0.99     => 8.11 DEGRES )
C          ( 0.9962   => 5    DEGRES )
C          ( 0.99756  => 4    DEGRES )
C          ( 0.99863  => 3    DEGRES )
C          ( 0.999    => 2.56 DEGRES )
C          ( 0.9999   => 0.8  DEGRES )
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C N1TETS : N1TETS(I) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET I
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C MXFACO : NOMBRE MAXIMAL DE FACES DECLARABLES DANS LEFACO
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
CCCC       LEFACO(12,*) NO DE FACE OC
C N1FASC : N1FASC(NS)=NUMERO (DANS LEFACO) D'UNE FACE DE SOMMET NS

C MXTE1S : MAX DE MOTS DU TABLEAU NOTE1S
C NOTE1S : TABLEAU AUXILIAIRE

C SORTIE :
C --------
C NB2F2F : NOMBRE D'ECHANGES  2Faces Perdues -> 2Faces Retrouvees
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC       JUIN 2002
C MODIFS : ALAIN PERRONNET  Saint PIERRE du PERRAY         Decembre 2019
C....................................................................012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  PTXYZD(1:4,1:*)
      INTEGER           NOTETR(8,*),
     %                  N1TETS(*),
     %                  LEFACO(11,0:MXFACO),
     %                  N1FASC(*),
     %                  NOTE1S(MXTE1S)

C     NOMBRE D'ECHANGES  2Faces Perdues -> 2Faces Retrouvees
      NB2F2F = 0

C     ====================================================================
C     BOUCLE SUR LES FACES LEFACO PERDUES I.E. NON DANS LA TETRAEDRISATION
C     ====================================================================
      DO 10 NFLPER = 1, MXFACO

C        LA FACE EST ELLE VIDE?
         IF( LEFACO( 1, NFLPER ) .EQ. 0 ) GOTO 10

C        LA FACE NFLPER DE LEFACO EST ELLE VRAIMENT PERDUE?
C        NTE EST LE NO NOTETR DU TETRAEDRE DE FACE NFLPER DANS LEFACO
         NTE = LEFACO( 11, NFLPER )
         IF( NTE .GT. 0 ) THEN

C           LE TETRAEDRE NTE A T IL VRAIMENT POUR FACE LA FACE NFLPER?
            CALL NO1F1T( LEFACO(1,NFLPER), NOTETR(1,NTE), NF )
            IF( NF .GT. 0 ) THEN

C              LA FACE NFLPER APPARTIENT AU TETRAEDRE NTE donc
C              CETTE FACE N'EST PAS PERDUE
               GOTO 10

            ELSE

C              NTE N'A PAS LA FACE NFLPER DE LEFACO
C              -> NFLPER SANS TETRAEDRE EST CONSIDEREE PERDUE
               LEFACO( 11, NFLPER ) = 0

            ENDIF

         ENDIF

C        TENTATIVE D'ECHANGE  2Faces Perdues -> 2Faces Retrouvees
         CALL VD2F2FFP( COANPL, PTXYZD, N1TETS, NOTETR,
     %                  NFLPER, MXFACO, LEFACO, N1FASC,
     %                  MXTE1S, NOTE1S, NONOUI )

C        NOMBRE D'ECHANGES  2Faces Perdues -> 2Faces Retrouvees
         NB2F2F = NB2F2F + NONOUI

 10   ENDDO

      IF( LANGAG .EQ. 0 ) THEN
        PRINT*,'vd2f2f:',NB2F2F,' ECHANGES 2F PERDUES ->2F RETROUVEES'
      ELSE
        PRINT*,'vd2f2f:',NB2F2F,' times 2 LOST FACES -> 2 REFOUND FACES'
      ENDIF

      RETURN
      END
