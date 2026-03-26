      SUBROUTINE RACPOL3( A, B, C, D, RAC0,
     %                    NBRAC, RAC1, RAC2, RAC3 )
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
C                   CLASSEES PAR VALEURS ABSOLUES CROISSANTES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN LJLL UPMC Paris & Veulettes          Mai 2014
C234567..............................................................012
      DOUBLE PRECISION  A, B, C, D, RAC0, RAC1, RAC2, RAC3,
     %                  X, X1, XM, P, P1, PM, BETA, GAMA
      INTEGER           NBRAC

      NBITER = 0

      IF( A .EQ. 0D0 ) THEN
C        A=0  => TRINOME DU SECOND DEGRE   B*X**2 + C*X + D=0
         IF( B .EQ. 0D0 ) THEN
C           A=0 B=0 => BINOME DU PREMIER DEGRE  C*X + D=0

            IF( C .EQ. 0D0 ) THEN
C              A=0 B=0 C=0 =>  D=0
               IF( D .EQ. 0 ) THEN
C                 TOUT REEL EST RACINE
                  NBRAC = -1
                  RAC1 = 0D0
               ELSE
C                 PAS DE RACINE
                  NBRAC = 0
                  RAC1  = 0D0
               ENDIF
            ELSE
C              A=0 B=0 => BINOME DU PREMIER DEGRE  C*X + D=0
               RAC1  = -D / C
               NBRAC = 1
            ENDIF
         ELSE
C           A=0  => TRINOME DU SECOND DEGRE   B*X**2 + C*X + D=0
            P = C * C - 4D0 * B * D
            IF( P .GE. 0D0 ) THEN
               IF( C .GE. 0 ) THEN
                  RAC1 = -C - SQRT( P ) / (2D0*B)
               ELSE
                  RAC1 = -C + SQRT( P ) / (2D0*B)
               ENDIF
               RAC2  = RAC1 - C / B
               NBRAC = 2
            ELSE
               NBRAC = 0
               RAC1  = 0D0
            ENDIF
         ENDIF

         GOTO 9000
      ENDIF

C     P(X) = A X**3 + B X**2 + C X + D   COMPLET
      IF( D .EQ. 0D0 ) THEN
         X1 = 0D0
         GOTO 50
      ENDIF

C     VALEUR INITIALE POUR LE SCHEMA DE NEWTON
      X = RAC0

      NBITER = 0

 10   IF( NBITER .LE. 1000 ) THEN
         NBITER = NBITER + 1
         X1 = X - ( ( (   A*X +   B ) * X + C ) * X + D )
     %            / ( ( 3*A*X + 2*B ) * X + C )
ccc         PRINT *,'RACPOL3: ITER=',NBITER,' X=', X,' X1=',X1

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
C     RECHERCHE D'UN INTERVALLE DE SIGNE OPPOSE DU POLYNOME
      X1 = 1D3
      X  = -X1

 20   P  = ( ( A*X  + B ) * X  + C ) * X  + D
      P1 = ( ( A*X1 + B ) * X1 + C ) * X1 + D

      IF( P * P1 .GT. 0D0 ) THEN
         X1 = X1 * 1D3
         X  =-X1
ccc         PRINT *,'RACPOL3: ITER=',NBITER,' INTERVALLE X=', X,' X1=',X1,
ccc     %                                   '  P=',P,' P1=',P1
         GOTO 20
      ENDIF

C     DICHOTOMIE DE L'INTERVALLE
 30   XM = ( X + X1 ) / 2
      PM = ( ( A*XM + B ) * XM + C ) * XM + D
      IF( P * PM .LE. 0D0 ) THEN
         X1 = XM
         P1 = PM
      ELSE
         X = XM
         P = PM
      ENDIF
      IF( ABS(X1-X) .GT. 1D-6*ABS(X1) ) GOTO 30

C     CALCUL DES 2 AUTRES RACINES A PARTIR DE X1
C     A X**3 + B X**2 + C X + D = A * (X-X1) * ( X**2 + BETA * X + GAMA )
 50   RAC1 = X1
      IF( X1 .EQ. 0D0 ) THEN
C        A X**3 + B X**2 + C X (+D=0) = A * X * ( X**2 + B/A * X + C/A )
         BETA = B / A
         GAMA = C / A
      ELSE
C        A X**3 + B X**2 + C X + D = A (X-X1) * ( X**2 + (X1+B/A)*X - D/(A*X1) )
         BETA = X1 + B / A
         GAMA = - D / ( A * X1 )
      ENDIF

C     RACINES DU TRINOME DU SECOND DEGRE EN FACTEUR DE (X-X1)
      P = BETA * BETA - 4D0 * GAMA
      IF( P .GE. 0D0 ) THEN
         IF( BETA .GE. 0 ) THEN
            RAC2 = ( -BETA -SQRT( P ) ) / 2D0
         ELSE
            RAC2 = ( -BETA + SQRT( P ) ) / 2D0
         ENDIF
         RAC3  = - RAC2 - BETA
         NBRAC = 3
      ELSE
C        UNE SEULE RACINE MISE A -1D89 DE RAC2 RAC3
         RAC2 = -1D89
         RAC3 = -1D89
         NBRAC = 1
      ENDIF

      IF( NBRAC .EQ. 3 ) THEN
C        RANGEMENT PAR VALEURS ABSOLUES CROISSANTES

         IF( ABS(RAC1) .GT. ABS(RAC2) ) THEN
            X    = RAC1
            RAC1 = RAC2
            RAC2 = X
         ENDIF

         IF( ABS(RAC1) .GT. ABS(RAC3) ) THEN
            X    = RAC1
            RAC1 = RAC3
            RAC3 = X
         ENDIF

         IF( ABS(RAC2) .GT. ABS(RAC3) ) THEN
            X    = RAC2
            RAC2 = RAC3
            RAC3 = X
         ENDIF
      ENDIF

cccC     AFFICHAGE FINAL SI 3 RACINES REELLES
ccc 9000 IF( NBRAC .EQ. 3 ) THEN
ccc         PRINT 19000, RAC0, ITER, NBRAC, RAC1, RAC2, RAC3
ccc19000    FORMAT('RACPOL3: RAC0=',G15.7,' ITER=',I4,'  NBRAC=',I1,
ccc     %          ' RAC1=',G15.7,' RAC2=',G15.7,' RAC3=',G15.7 )
ccc      ENDIF

 9000 RETURN
      END
