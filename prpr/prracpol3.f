      PROGRAM PRRACPOL3
      DOUBLE PRECISION  RAC1, RAC2, RAC3

C     P = X**3
      CALL RACPOL3( 1D0, 0D0, 0D0, 0D0, 1D5,
     %              NBRAC, RAC1, RAC2, RAC3 )

C     P= 3 (X-1) (X-2) (X-3)
      CALL RACPOL3( 3D0, -18D0, 33D0, -18D0, -1D5,
     %              NBRAC, RAC1, RAC2, RAC3 )

C     P= (X-7) (X**2 + X + 1)
      CALL RACPOL3( 1D0, -6D0, -6D0, -7D0, 1D5,
     %              NBRAC, RAC1, RAC2, RAC3 )

      STOP
      END


      SUBROUTINE RACPOL3( A, B, C, D, RAC0, NBRAC, RAC1, RAC2, RAC3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DES 1 OU 3 RACINES REELLES DU POLYNOME DU 3-EME DEGRE
C ----    P(X) = A X**3 + B X**2 + C X + D
C         PAR LA METHODE DE NEWTON

C ENTREES:
C --------
C A, B, C, D : LES 4 COEFFICIENTS DU POLYNOME DU 3-EME DEGRE
C RAC0       : ESTIMATION du CENTRE DE L'INTERVALLE ENCADRANT RAC1

C SORTIES:
C --------
C NBRAC  : NOMBRE DE RACINES REELLES DU POLYNOME
C RAC1, RAC2, RAC3: NBRAC RACINES REELLES DU POLYNOME DU 3-EME DEGRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN LJLL UPMC Paris & Veulettes          Mai 2014
C234567..............................................................012
      DOUBLE PRECISION  A, B, C, D, RAC0, RAC1, RAC2, RAC3,
     %                  X, X1, P, P1, BETA, GAMA
      INTEGER           NBRAC


      IF( A .EQ. 0D0 ) GOTO 100

      X = RAC0

 5    NBITER = 0
 10   IF( NBITER .LE. 1000 ) THEN
         NBITER = NBITER + 1
         X1 = X - ( ( (   A*X +   B ) * X + C ) * X + D )
     %            / ( ( 3*A*X + 2*B ) * X + C )
         PRINT *,'RACPOL3: ITER=',NBITER,' X=', X,' X1=',X1

         IF( ABS(X1-X) .GT. ABS(X1)*1D-10 ) THEN
C           NON CONVERGENCE
            X = X1
            GOTO 10
         ELSE
C           CONVERGENCE
            GOTO 50
         ENDIF
      ENDIF

C     PAS DE CONVERGENCE => DICHOTOMIE
      X  = -1D3
      X1 =  1D3

 20   P  = ( ( A*X  + B ) * X  + C ) * X  + D
      P1 = ( ( A*X1 + B ) * X1 + C ) * X1 + D

      IF( P * P1 .GT. 0D0 ) THEN
         X  = X  * 1D3
         X1 = X1 * 1D3
         PRINT *,'RACPOL3: ITER=',NBITER,' INTERVALLE X=', X,' X1=',X1,
     %                                   '  P=',P,' P1=',P1
         GOTO 20
      ENDIF

      IF( P .GE. P1 ) THEN
         X = X1
      ENDIF
      GOTO 5

C     CALCUL DES 2 AUTRES RACINES
C     A X**3 + B X**2 + C X + D = A (X-X1) * ( X**2 + (X1+B/A)*X - D/(A*X1) )
 50   RAC1 = X1

      IF( X1 .EQ. 0D0 ) THEN
C        X1 EST SOLUTION TRIPLE
         RAC2 = RAC1
         RAC3 = RAC1
         NBARC = 3
         GOTO 9000
      ENDIF

      BETA = X1 + B / A
      GAMA = - D / ( A * X1 )

C     RACINES DU TRINOME DU SECOND DEGRE EN FACTEUR
      P = BETA * BETA - 4D0 * GAMA
      IF( P .GE. 0D0 ) THEN
         IF( BETA .GE. 0 ) THEN
            RAC2 = ( -BETA -SQRT( P ) ) / 2D0
         ELSE
            RAC2 = ( -BETA + SQRT( P ) ) / 2D0
         ENDIF
         RAC3 = - RAC2 - BETA
         NBRAC = 3
      ELSE
C        UNE SEULE RACINE MISE A 1D100 DES AUTRES
         RAC2 = 1D100
         RAC3 = 1D100
         NBRAC = 1
      ENDIF
      GOTO 9000

C     A=0 => TRINOME DU SECOND DEGRE   B*X**2 + C*X + D=0
 100  P = C * C - 4D0 * B * D
      IF( P .GE. 0D0 ) THEN
         IF( C .GE. 0 ) THEN
            RAC1 = -C - SQRT( P ) / (2D0*B)
         ELSE
            RAC1 = -C + SQRT( P ) / (2D0*B)
         ENDIF
         RAC2 = RAC1 - C / B
         NBRAC = 2
      ELSE
         NBRAC = 0
      ENDIF

C     AFFICHAGE FINAL DES RACINES
 9000 PRINT *,'RACPOL3: NBRAC=',NBRAC,' RAC1=',RAC1,
     %        ' RAC2=',RAC2,' RAC3=',RAC3
      PRINT *
      RETURN
      END
