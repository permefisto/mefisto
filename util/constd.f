      SUBROUTINE CONSTD( DCTE, L, DA )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INITIALISATION A DCTE D UN TABLEAU DE L
C ----- VARIABLES REELLES EN DOUBLE PRECISION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN Laboratoire J-L. LIONS UPMC PARIS  Avril 2007
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  DCTE, DA(1:L)
C
      DO 1 I = 1 , L
         DA(I) = DCTE
    1 CONTINUE
      END
