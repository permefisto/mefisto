      SUBROUTINE AJTEOPTD( QUAMINEX, PTXYZD, NBGRTE, NOGRTE, NOTETR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  AJOUTER AUX NBGRTE TETRAEDRES LES TETRAEDRES OPPOSES AUX
C -----  4 FACES DES TETRAEDRES DEGENERES

C ENTREES:
C --------
C QUAMINEX: QUALITE AU DESSOUS DE LAQUELLE UN TETRAEDRE N'EST PAS CREE
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NBGRTE : NOMBRE INITIAL DU GROUPE DE TETRAEDRES
C NOGRTE : NUMERO DANS NOTETR DES NBGRTE TETRAEDRES
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412

C SORTIES :
C ---------
C NBGRTE : NOMBRE FINAL DE TETRAEDRES DU GROUPE DE TETRAEDRES
C NOGRTE : NUMERO DANS NOTETR DES NBGRTE TETRAEDRES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2017
C2345X7..............................................................012
      DOUBLE PRECISION  PTXYZD(4,*),
     %                  ARMIN, ARMAX, SURFTR(4), VOLUTE, V
      INTEGER           NOGRTE(NBGRTE), NOTETR(8,*)

      NBGRTE0 = NBGRTE
      DO N=1,NBGRTE

C        LE TETRAEDRE DE NOTETR
         NTE = NOGRTE( N )

         IF( NTE .GT. 0 ) THEN

            CALL QUATETD( PTXYZD(1,NOTETR(1,NTE)),
     %                    PTXYZD(1,NOTETR(2,NTE)),
     %                    PTXYZD(1,NOTETR(3,NTE)),
     %                    PTXYZD(1,NOTETR(4,NTE)),
     %                    ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )

            IF( QUALTE .LE. QUAMINEX ) THEN

               PRINT*,'ajteoptd: AJOUT des tetraedres opposes au TETRAED
     %RE(',NTE,')=',(NOTETR(I,NTE),I=1,8),' V=',VOLUTE,' Q=',QUALTE

               DO 10 K=1,4

C                 TETRAEDRE DE NOTETR OPPOSE A LA FACE K DE NTE
                  NTEOP = NOTETR( 4+K, NTE )
                  IF( NTEOP .GT. 0 ) THEN

                     CALL QUATETD( PTXYZD(1,NOTETR(1,NTEOP)),
     %                             PTXYZD(1,NOTETR(2,NTEOP)),
     %                             PTXYZD(1,NOTETR(3,NTEOP)),
     %                             PTXYZD(1,NOTETR(4,NTEOP)),
     %                             ARMIN, ARMAX, SURFTR, V, Q )

                     IF( Q .LT. 0.01 ) GOTO 10

C                    NTEOP EST DE BONNE QUALITE
                     DO L=1,NBGRTE
                        IF( NTEOP .EQ. NOGRTE(L) ) GOTO 10
                     ENDDO

C                    NTEOP N'EST PAS DANS NOGRTE. IL EST AJOUTE
                     NBGRTE = NBGRTE + 1
                     NOGRTE( NBGRTE ) = NTEOP

                  ENDIF
 10            ENDDO

            ENDIF

         ENDIF

      ENDDO

      PRINT*,'ajteoptd: AJOUT de',NBGRTE-NBGRTE0,' TETRAEDRES OPPOSES'
      RETURN
      END
