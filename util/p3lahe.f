       SUBROUTINE P3LAHE( X0, X1, X2, X3, DX0, DX1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER LA DERIVEE AU POINT 0 ET - LA DERIVEE AU POINT 1
C -----   D'UN POLYNOME P DE DEGRE 3 LAGRANGE
C         TEL QUE P(0)=X0, P(1/3)=X1, P(2/3)=X2, P(1)=X3
C
C ENTREES :
C ---------
C X0, X1, X2, X3 : LES 4 POINTS DE PASSAGE DU POLYNOME
C                  TELS QUE P(0)=X0, P(1/3)=X1, P(2/3)=X2, P(1)=X3
C SORTIES :
C ---------
C DX0    : +P'(0)
C DX1    : -P'(1)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1996
C2345X7..............................................................012
      DX0 = (18 * X1 - 11 * X0 -  9 * X2 +  2 * X3 ) * 0.5
      DX1 = ( 2 * X0 -  9 * X1 + 18 * X2 - 11 * X3 ) * 0.5
      END
