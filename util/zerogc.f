      SUBROUTINE ZEROGC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  MODIFIER LA PRECISION POUR LA CONVERGENCE DU GRADIENT CONJUGUE
C -----  OU LE TEST DE DEFINIE POSITIVITE DE LA FACTORISATION DE CHOLESKY
C
C        EPZERO EST DANS LE COMMON / EPSSSS /
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        JUIN 1999
C2345X7..............................................................012
      COMMON / EPSSSS / EPZERO, EPSXYZ
C
C     LECTURE DU MOT CLE A TRAITER
C     ----------------------------
      CALL INVITE( 92 )
      NCVALS = 5
      R      = EPZERO
      CALL LIRRSP( NCVALS , R )
      IF( R .GE. 0 ) THEN
          EPZERO = R
      ENDIF
C
      RETURN
      END
