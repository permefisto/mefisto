      SUBROUTINE INTARTR0( S1, S2, P1, P2, P3, LINTER, PTI, CBTR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 3 COORDONNEES DU POINT PTI, INTERSECTION DU
C -----    SEGMENT S1-S2 ET DU TRIANGLE DEFINI PAR SES SOMMETS P1, P2, P3
C          DANS LE CAS OU PTI EST INTERNE AU TRIANGLE P1P2P3

C          PAS DE RECHERCHE DE LA POSITION DE PTI SUR LA DROITE S1 S2

C ENTREES:
C --------
C S1,S2   : LES 2 POINTS QUI DEFINISSENT L'ARETE
C P1,P2,P3: LES 3 POINTS QUI DEFINISSENT LE TRIANGLE

C SORTIES:
C --------
C LINTER : -2 SI S1=S2  PAS DE CALCUL DE PTI
C          -1 SI S1-S2 PARALLELE AU PLAN DU TRIANGLE  PAS DE CALCUL DE PTI
C           0 SI S1-S2 N'INTERSECTE PAS LE TRIANGLE ENTRE CES 3 SOMMETS
C           1 SI S1-S2   INTERSECTE     LE TRIANGLE ENTRE CES 3 SOMMETS
C             C-A-D  PTI EST INTERNE AU TRIANGLE OU SUR UNE DE SES 3 ARETES
C PTI    : LES 3 COORDONNEES DU POINT D'INTERSECTION SI LINTER=1
C CBTR   : LES 3 COORDONNEES BARYCENTRIQUES DE PTI DANS LE TRIANGLE P1P2P3
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY   JUILLET 2009
C MODIFS: ALAIN PERRONNET LJLL UPMC & VEULETTES SUR MER   SEPTEMBRE 2014
C MODIFS: ALAIN PERRONNET Saint Pierre du Perray           Novembre 2019
C MODIFS: ALAIN PERRONNET Saint Pierre du Perray            Octobre 2020
C23456...............................................................012
      IMPLICIT  NONE
      DOUBLE PRECISION  S1(3),S2(3), P1(3),P2(3),P3(3), PTI(3), CBTR(3)
      DOUBLE PRECISION  A, B, C, D, E, COSA, F, S, R
      INTEGER           LINTER

C     LE CARRE DE LA DISTANCE S1-S2
      E = (S2(1)-S1(1))**2 + (S2(2)-S1(2))**2 + (S2(3)-S1(3))**2
      IF( E .LE. 0D0 ) THEN
         LINTER = -2
         GOTO 9999
      ENDIF

C     LE VECTEUR NORMAL AU PLAN DES 3 POINTS
      A = ( P2(2) - P1(2) ) * ( P3(3) - P1(3) )
     %  - ( P2(3) - P1(3) ) * ( P3(2) - P1(2) )
      B = ( P2(3) - P1(3) ) * ( P3(1) - P1(1) )
     %  - ( P2(1) - P1(1) ) * ( P3(3) - P1(3) )
      C = ( P2(1) - P1(1) ) * ( P3(2) - P1(2) )
     %  - ( P2(2) - P1(2) ) * ( P3(1) - P1(1) )
      D = A * A + B * B + C * C

C     LE COSINUS DE L'ANGLE ENTRE LE VECTEUR NORMAL ET S1-S2
      F = (S2(1)-S1(1))*A + (S2(2)-S1(2))*B + (S2(3)-S1(3))*C

C     COSINUS( S1-S2 , NORMALE AU PLAN DU TRIANGLE )
      COSA = F / SQRT(E*D)

      IF( ABS(COSA) .LT. 1D-4 ) THEN
C        COSA=1D-5 => A=89.999427 degres
C        COSA=1D-4 => A=89.994270 degres
C        COSA=1D-3 => A=89.942704 degres
C        S1-S2 PARALLELE AU PLAN
         LINTER = -1
         GOTO 9999
      ENDIF

C     LE POINT D'INTERSECTION DE L'ARETE S1-S2 AVEC LE PLAN
      S = ( (P1(1)-S1(1))*A + (P1(2)-S1(2))*B + (P1(3)-S1(3))*C ) / F
      PTI(1) = S1(1) + S * ( S2(1) - S1(1) )
      PTI(2) = S1(2) + S * ( S2(2) - S1(2) )
      PTI(3) = S1(3) + S * ( S2(3) - S1(3) )

C     CE POINT PTI EST IL INTERNE AU TRIANGLE?
C     CALCUL DES 3 COORDONNEES BARYCENTRIQUES DE PTI

C     TRIANGLE PTI-P2 PTI-P3
      A = ( P2(2) - PTI(2) ) * ( P3(3) - PTI(3) )
     %  - ( P2(3) - PTI(3) ) * ( P3(2) - PTI(2) )
      B = ( P2(3) - PTI(3) ) * ( P3(1) - PTI(1) )
     %  - ( P2(1) - PTI(1) ) * ( P3(3) - PTI(3) )
      C = ( P2(1) - PTI(1) ) * ( P3(2) - PTI(2) )
     %  - ( P2(2) - PTI(2) ) * ( P3(1) - PTI(1) )
      R = ( A * A + B * B + C * C ) / D
      IF( R .GT. 1.007D0 ) GOTO 9900
      CBTR(1) = SQRT( R )

C     TRIANGLE PTI-P3 PTI-P1
      A = ( P3(2) - PTI(2) ) * ( P1(3) - PTI(3) )
     %  - ( P3(3) - PTI(3) ) * ( P1(2) - PTI(2) )
      B = ( P3(3) - PTI(3) ) * ( P1(1) - PTI(1) )
     %  - ( P3(1) - PTI(1) ) * ( P1(3) - PTI(3) )
      C = ( P3(1) - PTI(1) ) * ( P1(2) - PTI(2) )
     %  - ( P3(2) - PTI(2) ) * ( P1(1) - PTI(1) )
      R = ( A * A + B * B + C * C ) / D
      IF( R .GT. 1.007D0 ) GOTO 9900
      CBTR(2) = SQRT( R )

C     TRIANGLE PTI-P1 PTI-P2
      A = ( P1(2) - PTI(2) ) * ( P2(3) - PTI(3) )
     %  - ( P1(3) - PTI(3) ) * ( P2(2) - PTI(2) )
      B = ( P1(3) - PTI(3) ) * ( P2(1) - PTI(1) )
     %  - ( P1(1) - PTI(1) ) * ( P2(3) - PTI(3) )
      C = ( P1(1) - PTI(1) ) * ( P2(2) - PTI(2) )
     %  - ( P1(2) - PTI(2) ) * ( P2(1) - PTI(1) )
      R = ( A * A + B * B + C * C ) / D
      IF( R .GT. 1.007D0 ) GOTO 9900
      CBTR(3) = SQRT( R )

      IF( CBTR(1) .GE. -0.001D0 .AND. CBTR(1) .LE. 1.001D0 .AND.
     %    CBTR(2) .GE. -0.001D0 .AND. CBTR(2) .LE. 1.001D0 .AND.
     %    CBTR(3) .GE. -0.001D0 .AND. CBTR(3) .LE. 1.001D0 ) THEN
C        LE POINT EST DANS LE TRIANGLE OU SUR UNE DE SES 3 ARETES
         LINTER = 1
         GOTO 9999
      ENDIF

C     La DROITE S1-S2 N'INTERSECTE PAS LE TRIANGLE ENTRE CES 3 SOMMETS
 9900 LINTER = 0

 9999 RETURN
      END
