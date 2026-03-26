      SUBROUTINE RETETROP( NTE, MXTETR, NOTETR, NBTEOP )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     RETROUVER NTE COMME TETRAEDRE OPPOSE A COMBIEN DE TETRAEDRES
C ----

C ENTREES:
C --------
C NTE    : NUMERO NOTETR DU TETRAEDRE A RETROUVER OPPOSE
C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES DU TABLEAU NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C          NOTETR(5,N) = TETRAEDRE VIDE SUIVANT SI NOTETR(1,NT)=0

C SORTIE :
C --------
C NBTEOP : =0 SI NTE N'EST PAS UN TETRAEDRE OPPOSE
C          >0 NOMBRE DE TETRAEDRES DONT LE TETRAEDRE NTE EST L'OPPOSE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY          Juillet 2020
C2345X7..............................................................012
      INTEGER NOTETR(8,MXTETR)

      NBTEOP = 0

      PRINT*
      PRINT*,'retetrop: le TETRAEDRE',NTE,':',
     %       (NOTETR(K,NTE),K=1,8),' EST OPPOSE'

      DO NT = 1, MXTETR

         IF( NOTETR(1,NT) .GT. 0 ) THEN

            DO NFTE=1,4

               NTEOP = NOTETR( 4+NFTE, NT )
               IF( NTEOP .EQ. NTE ) THEN
                  NBTEOP = NBTEOP + 1
                  PRINT*,'retetrop: au TETRAEDRE',NT,':',
     %                   (NOTETR(K,NT),K=1,8)
               ENDIF

            ENDDO
         ENDIF
      ENDDO

      PRINT*,'retetrop: le TETRAEDRE',NTE,' EST OPPOSE a',NBTEOP,
     %       ' TETRAEDRES'

      RETURN
      END


