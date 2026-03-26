      SUBROUTINE FQ1INV( X, Y, S, XC, YC, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DES 2 COORDONNEES (XC,YC) DANS LE CARRE (0,1)
C -----   IMAGE PAR F:CARRE UNITE-->QUADRANGLE APPARTENANT A Q1**2
C         PAR UNE RESOLUTION DIRECTE DUE A NICOLAS THENAULT
C
C ENTREES:
C --------
C X,Y   : COORDONNEES DU POINT IMAGE DANS LE QUADRANGLE DE SOMMETS S
C S     : LES 2 COORDONNEES DES 4 SOMMETS DU QUADRANGLE
C
C SORTIES:
C --------
C XC,YC : COORDONNEES DANS LE CARRE DONT L'IMAGE PAR F VAUT (X,Y)
C IERR  : 0 SI CALCUL SANS ERREUR, 1 SI QUADRANGLE DEGENERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: THENAULT TULENEW  ANALYSE NUMERIQUE PARIS        JANVIER 1998
C MODIFS : PERRONNET ALAIN   ANALYSE NUMERIQUE PARIS        JANVIER 1998
C MODIFS : PERRONNET ALAIN   ANALYSE NUMERIQUE PARIS        MARS    2006
C234567..............................................................012
c     zero a ajuster eventuellement
      double precision zero
      parameter       (zero=1d-15)
C
      REAL             S(1:2,1:4), DIST(2)
      DOUBLE PRECISION A,B,C,D,ALPHA,BETA,GAMMA,DELTA,X0,Y0,T(2),U,V,W
      INTRINSIC        REAL
C
      A = S(1,1)
      B = S(1,2) - S(1,1)
      C = S(1,4) - S(1,1)
      D = S(1,1) - S(1,2) + S(1,3) - S(1,4)
C
      ALPHA = S(2,1)
      BETA  = S(2,2) - S(2,1)
      GAMMA = S(2,4) - S(2,1)
      DELTA = S(2,1) - S(2,2) + S(2,3) - S(2,4)
C
      U = BETA  * C - B * GAMMA
      IF( ABS(U) .LE. ZERO ) THEN
C        QUADRANGLE DEGENERE
         IERR = 1
         RETURN
      ENDIF
      V = DELTA * C - D * GAMMA
      W = B * DELTA - BETA * D
C
      X0 = C * (Y-ALPHA) - GAMMA * (X-A)
      Y0 = B * (Y-ALPHA) - BETA  * (X-A)
C
      A = V  * W
      B = U  * U - W * X0 - V * Y0
      C = X0 * Y0
C
      IF( ABS(A) .GT. ZERO ) THEN
C
         DELTA = SQRT( B*B-4*A*C )
         IF( B .GE. 0.0 ) THEN
            T(2) = -B - DELTA
         ELSE
            T(2) = -B + DELTA
         ENDIF
C        LA RACINE DE PLUS GRANDE VALEUR ABSOLUE
C       (ELLE DONNE LE PLUS SOUVENT LE POINT EXTERIEUR AU CARRE UNITE
C        DONC A TESTER EN SECOND POUR REDUIRE LES CALCULS)
         T(2) = T(2) / ( 2 * A )
C        CALCUL DE LA SECONDE RACINE A PARTIR DE LA SOMME => PLUS STABLE
         T(1) = - B/A - T(2)
C
         DO 10 I=1,2
C
C           LA SOLUTION I DONNE T ELLE UN POINT INTERNE AU CARRE UNITE?
            XC = REAL( ( X0 - V * T(I) ) / U )
            YC = REAL( ( W * T(I) - Y0 ) / U )
            IF( 0.0 .LE. XC .AND. XC .LE. 1.0 ) THEN
               IF( 0.0 .LE. YC .AND. YC .LE. 1.0 ) GOTO 9000
            ENDIF
C
C           LE POINT (XC,YC) N'EST PAS DANS LE CARRE UNITE
C           CELA PEUT ETRE DU AUX ERREURS D'ARRONDI
C           => CHOIX PAR LE MINIMUM DE LA DISTANCE AUX BORDS DU CARRE
            DIST(I) = MAX( 0.0, -XC, XC-1.0, -YC, YC-1.0 )
C
 10      CONTINUE
C
         IF( DIST(1) .GT. DIST(2) ) THEN
C           F(XC,YC) POUR LA RACINE 2 EST PLUS PROCHE DE X,Y
C           XC YC SONT DEJA CALCULES
            GOTO 9000
         ENDIF
C
      ELSE IF ( ABS(B) .GT. ZERO ) THEN
         T(1) = - C / B
      ELSE
         T(1) = 0
      ENDIF
C
C     LES 2 COORDONNEES DU POINT DANS LE CARRE UNITE
      XC = REAL( ( X0 - V * T(1) ) / U )
      YC = REAL( ( W * T(1) - Y0 ) / U )
C
 9000 IERR = 0
      RETURN
      END
