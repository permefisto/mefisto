      SUBROUTINE VOLTETRA( PTXYZD, MXTETR, NOTETR, NBTETR, NUDTETR,
     %                     VOLTOTAL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:    CALCUL DU VOLUME TOTAL DES NBTETR TETRAEDRES DE NOTETR
C ----

C ENTREES:
C --------
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C MXTETR : NOMBRE MAXIMUM DE TETRAEDRES UTISABLES DANS LE TABLEAU NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412

C SORTIES:
C --------
C NBTETR :  NOMBRE DE TETRAEDRES DE NOTETR
C NUDTETR:  NUMERO DU DERNIER TETRAEDRE DE NOTETR
C VOLTOTAL: VOLUME TOTAL DES NBTETR TETRAEDRES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET                St PIERRE du PERRAY  Juin 2018
C2345X7..............................................................012
      INTEGER           NOTETR(8,MXTETR)
      DOUBLE PRECISION  PTXYZD(4,*), VOLTOTAL, VOLTET, V

      NBTETR = 0
      VOLTOTAL = 0D0
      DO NTE = 1, MXTETR

         IF( NTE .GT. 0 .AND. NOTETR(1,NTE ) .GT. 0 ) THEN

C           VOLUME DU TETRAEDRE NTE
            V = VOLTET( PTXYZD( 1, NOTETR(1,NTE) ),
     %                  PTXYZD( 1, NOTETR(2,NTE) ),
     %                  PTXYZD( 1, NOTETR(3,NTE) ),
     %                  PTXYZD( 1, NOTETR(4,NTE) ) )

C           VOLUME CUMULE
            VOLTOTAL = VOLTOTAL + V

            NBTETR  = NBTETR + 1
            NUDTETR = NTE

         ENDIF

      ENDDO

      RETURN
      END
