      SUBROUTINE ETOIL1TET( NPt,    POINT, PTXYZD, NUDTETR, NOTETR,
     %                      NOTET1, COBARY )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TROUVER PARMI LES TETRAEDRES NUDTETR A 1 UN TETRAEDRE NOTET1
C -----    CONTENANT LE POINT NPt (POINT NE PEUT ETRE UN SOMMET)
C
C ENTREES:
C --------
C NPt    : NUMERO PTXYZD DU POINT A TETRAEDRISER
C POINT  : 3 COORDONNEES DU POINT A TETRAEDRISER
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C
C SORTIES:
C --------
C NOTET1 : >0 NUMERO NOTETR DU TETRAEDRE CONTENANT LE POINT NPt
C          =0 SI PAS DE TETRAEDRE CONTENANT LE POINT NPt
C COBARY : LES 4 COORDONNEES BARYCENTRIQUES DE POINT DANS NOTET1
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC et St Pierre du Perray    Aout 2014
C23456...............................................................012
      INTEGER           NOTETR(8,NUDTETR)
      DOUBLE PRECISION  PTXYZD(1:4,1:*), POINT(3), COBARY(4),
     %                  CBSOM, CBSOMIN, V

      PRINT*
      PRINT*,'etoil1tet: Bizarre NPt=',NPt,' XYZ=',POINT,
     %       ' est RECHERCHE DANS la TOTALITE des',NUDTETR,' TETRAEDRES'

C     RECHERCHE D'UN TETRAEDRE NOTET1 CONTENANT LE POINT
C     ==================================================
      CBSOMIN = 1D100
      NTEMIN  = 0

      DO NOTET1 = NUDTETR, 1, -1

         IF( NOTETR(1,NOTET1) .GT. 0 ) THEN

            CALL COBATE( POINT,PTXYZD,NOTETR(1,NOTET1), V,COBARY, IERR )
            IF( IERR .EQ. 0 ) THEN

               CBSOM = ABS( COBARY(1) ) + ABS( COBARY(2) )
     %               + ABS( COBARY(3) ) + ABS( COBARY(4) )

ccc               IF( CBSOM .LE. 1.001D0 ) THEN
         IF( COBARY(1) .GE. -0.001D0 .AND. COBARY(1) .LE. 1.001D0 .AND.
     %       COBARY(2) .GE. -0.001D0 .AND. COBARY(2) .LE. 1.001D0 .AND.
     %       COBARY(3) .GE. -0.001D0 .AND. COBARY(3) .LE. 1.001D0 .AND.
     %       COBARY(4) .GE. -0.001D0 .AND. COBARY(4) .LE. 1.001D0 ) THEN
                  PRINT*,'etoil1tet: NPt=',NPt,
     %                   ' EST RETROUVE dans le TETRAEDRE=',NOTET1,
     %                   ' avec CoorBARY=',COBARY,' CBSom=',CBSOM
                  GOTO 9000
               ENDIF

C              STOKAGE DU PLUS PROCHE DE 1D0
               IF( CBSOM .LT. CBSOMIN ) THEN
                  CBSOMIN = CBSOM
                  NTEMIN  = NOTET1
               ENDIF

            ENDIF

         ENDIF

      ENDDO

      IF( NTEMIN .LE. 0 .OR. CBSOMIN .GT. 1.001D0 ) THEN
C        PAS DE TETRAEDRE CONTENANT LE POINT NPt
         PRINT*,'etoil1tet: NPt=', NPt,
     %          ' EXTERNE AUX', NUDTETR,' TETRAEDRES'
         NOTET1 = 0
      ELSE
C        NTEMIN LE TETRAEDRE LE PLUS PROCHE DU POINT NPt
         CALL COBATE( POINT, PTXYZD, NOTETR(1,NTEMIN), V, COBARY, IERR )
         PRINT*,'etoil1tet: NPt=', NPt,
     %          ' PROCHE du TETRAEDRE',NTEMIN,' de VOLUME=',V,
     %          ' COBARY=',COBARY,' CBSOMin=',CBSOMIN
         NOTET1 = NTEMIN
      ENDIF
      PRINT*

 9000 RETURN
      END
