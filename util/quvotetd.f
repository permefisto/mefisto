      SUBROUTINE QUVOTETD( QUAMINEX, PTXYZD, NBTEQV, NOTEQV, NOTETR,
     %                     QUAMIN,   QUAMOY, NBTEQM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA QUALITE MINIMALE DES NBTEQV TETRAEDRES DE NOTETR
C -----

C ENTREES:
C --------
C QUAMINEX:QUALITE AU DESSOUS DE LAQUELLE UN TETRAEDRE N'EST PAS CREE ou
C          QUALITE DES TETRAEDRES AU DESSOUS DE LAQUELLE UNE
C          AMELIORATION DE LA QUALITE DES TETRAEDRES EST DEMANDEE
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NBTEQV : NOMBRE DE TETRAEDRES DE QUALITE A CALCULER
C NOTEQV : NUMERO DANS NOTETR DES NBTEQV TETRAEDRES
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412

C SORTIES :
C ---------
C QUAMIN : QUALITE MINIMALE DES NBTEQV TETRAEDRES
C QUAMOY : QUALITE MOYENNE  DES NBTEQV TETRAEDRES
C NBTEQM : NOMBRE DE TETRAEDRES DE QUALITE MEDIOCRE AU DESSOUS DE QUAMINEX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2017
C2345X7..............................................................012
      DOUBLE PRECISION  PTXYZD(4,*),
     %                  ARMIN, ARMAX, SURFTR(4), VOLUTE
      INTEGER           NOTEQV(NBTEQV), NOTETR(8,*)

      QUAMIN = 2.
      QUAMOY = 0.
      NBTEQM = 0
      NBTE   = 0
      DO N=1,NBTEQV

C        LE TETRAEDRE DE NOTETR
         NTE = NOTEQV( N )

         IF( NTE .GT. 0 .AND. NOTETR(1,NTE) .GT. 0 ) THEN

            CALL QUATETD( PTXYZD(1,NOTETR(1,NTE)),
     %                    PTXYZD(1,NOTETR(2,NTE)),
     %                    PTXYZD(1,NOTETR(3,NTE)),
     %                    PTXYZD(1,NOTETR(4,NTE)),
     %                    ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )

            IF( QUALTE .LE. QUAMINEX ) THEN
               NBTEQM = NBTEQM + 1
ccc               PRINT *,'quvotetd: TETRAEDRE(',NTE,')=',
ccc     %                 (NOTETR(I,NTE),I=1,8),' V=',VOLUTE,' Q=',QUALTE
            ENDIF

            IF( QUALTE .GE. 0 ) THEN
               NBTE   = NBTE + 1
               QUAMOY = QUAMOY + QUALTE
            ENDIF

            IF( QUALTE .LT. QUAMIN ) THEN
               QUAMIN = QUALTE
            ENDIF

         ENDIF

      ENDDO

      IF( NBTE .GT. 0 ) THEN
         QUAMOY = QUAMOY / NBTE
      ENDIF

ccc      PRINT*,'quvotetd: Qualite MIN=',QUAMIN,' MOYENNE=',QUAMOY,
ccc     %       ' des',NBTEQV,' TETRAEDRES dont',NBTEQM,
ccc     %       ' de QUALITE MEDIOCRE <=',QUAMINEX

      RETURN
      END
