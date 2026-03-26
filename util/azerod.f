      SUBROUTINE AZEROD( NB, DA )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INITIALISATION A ZERO D UN TABLEAU DE NB
C ----- VARIABLES REELLES EN DOUBLE PRECISION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C MODIFS : PERRONNET ALAIN LJLL UPMC & St Pierre du Perray    AVRIL 2013
C23456---------------------------------------------------------------012
      INTEGER           NB, I
      DOUBLE PRECISION  DA(1:NB)

C        EXECUTION AVEC UN SEUL THREAD
         DO I = 1, NB
            DA(I) = 0.D0
         ENDDO

      RETURN
      END
