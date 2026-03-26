      SUBROUTINE M33INV( A, DET, AM1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE INVERSE D'UNE MATRICE 3 x 3
C -----
C          ATTENTION: LA MATRICE A DOIT ETRE DIFFERENTE DE AM1 !
C ENTREES:
C --------
C A      : LA MATRICE INITIALE 3 x 3
C
C SORTIES:
C --------
C DET    : LE DETERMINANT DE LA MATRICE A
C          0D0 SI MATRICE NON INVERSIBLE ET AM1 N'EST PAS CALCULEE
C AM1    : LA MATRICE INVERSE  3 x 3 de A
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     OCTOBRE 1994
C2345X7..............................................................012
      DOUBLE PRECISION  A(3,3), AM1(3,3), DET, D, DETM33
C
C     LE DETERMINANT DE A
      DET = DETM33( A(1,1), A(1,2), A(1,3),
     %              A(2,1), A(2,2), A(2,3),
     %              A(3,1), A(3,2), A(3,3) )
C
      IF( ABS(DET) .LT. 1D-40 )  THEN
C        MATRICE CONSIDEREE COMME NON INVERSIBLE
         DET = 0D0
         RETURN
      ENDIF
      D = 1D0 / DET
C
C     LES 9 COEFFICIENTS DE LA MATRICE INVERSE
      AM1(1,1) = ( A(2,2) * A(3,3) - A(3,2) * A(2,3) ) * D
      AM1(2,1) = ( A(2,3) * A(3,1) - A(3,3) * A(2,1) ) * D
      AM1(3,1) = ( A(2,1) * A(3,2) - A(3,1) * A(2,2) ) * D
C
      AM1(1,2) = ( A(1,3) * A(3,2) - A(1,2) * A(3,3) ) * D
      AM1(2,2) = ( A(1,1) * A(3,3) - A(1,3) * A(3,1) ) * D
      AM1(3,2) = ( A(1,2) * A(3,1) - A(1,1) * A(3,2) ) * D
C
      AM1(1,3) = ( A(1,2) * A(2,3) - A(2,2) * A(1,3) ) * D
      AM1(2,3) = ( A(1,3) * A(2,1) - A(2,3) * A(1,1) ) * D
      AM1(3,3) = ( A(1,1) * A(2,2) - A(2,1) * A(1,2) ) * D
C
      RETURN
      END
