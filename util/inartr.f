      SUBROUTINE INARTR( S1, S2, P1, P2, P3, LINTER, PTINT, CBT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 3 COORDONNEES DU POINT PTINT, INTERSECTION DU
C -----    SEGMENT S1-S2 ET DU TRIANGLE DEFINI PAR SES SOMMETS P1, P2, P3
C
C ENTREES:
C --------
C S1,S2   : LES 2 POINTS QUI DEFINISSENT L'ARETE
C P1,P2,P3: LES 3 POINTS QUI DEFINISSENT LE TRIANGLE
C
C SORTIES:
C --------
C LINTER : -3 SI S1=S2 PAS DE CALCUL DE PTINT
C          -2 SI S1-S2 EST DANS  LE PLAN DU TRIANGLE
C          -1 SI S1-S2 PARALLELE AU PLAN DU TRIANGLE SANS ETRE DEDANS
C           0 SI S1-S2 N'INTERSECTE PAS LE TRIANGLE ENTRE CES 3 SOMMETS
C                      ET ENTRE S1-S2
C           1 SI S1-S2   INTERSECTE     LE TRIANGLE ENTRE CES 3 SOMMETS
C                      ET ENTRE S1-S2
C           PTINT EST DANS LE TRIANGLE OU SUR UNE DE SES 3 ARETES
C PTINT  : LES 3 COORDONNEES DU POINT D'INTERSECTION SI LINTER=1
C CBT    : LES 3 COORDONNEES BARYCENTRIQUES DE PTINT DANS LE TRIANGLE P1P2P3
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY   JUILLET 2009
C MODIFS: ALAIN PERRONNET LJLL UPMC & VEULETTES SUR MER   SEPTEMBRE 2014
C MODIFS: ALAIN PERRONNET St Pierre du Perray               Octobre 2020
C23456...............................................................012
      DOUBLE PRECISION  S1(3), S2(3), P1(3), P2(3), P3(3), PTINT(3),
     %                  CBA(3), CBT(3), CBT123,
     %                  A, B, C, D, E, COSA, F, S

      LINTER = 0

C     LE CARRE DE LA DISTANCE S1-S2
      E = (S2(1)-S1(1))**2 + (S2(2)-S1(2))**2 + (S2(3)-S1(3))**2
      IF( E .LE. 0D0 ) THEN
         LINTER = -3
         GOTO 9999
      ENDIF

C     LE VECTEUR NORMAL AU PLAN DES 3 POINTS
      A = ( P2(2) - P1(2) ) * ( P3(3) - P1(3) )
     %  - ( P2(3) - P1(3) ) * ( P3(2) - P1(2) )
      B = ( P2(3) - P1(3) ) * ( P3(1) - P1(1) )
     %  - ( P2(1) - P1(1) ) * ( P3(3) - P1(3) )
      C = ( P2(1) - P1(1) ) * ( P3(2) - P1(2) )
     %  - ( P2(2) - P1(2) ) * ( P3(1) - P1(1) )

C     LE CARRE DE LA NORME DU VECTEUR NORMAL
      D = A * A + B * B + C * C

C     LE COSINUS DE L'ANGLE ENTRE LE VECTEUR NORMAL ET S1-S2
      F = (S2(1)-S1(1))*A + (S2(2)-S1(2))*B + (S2(3)-S1(3))*C

C     COSINUS( S1-S2 , NORMALE AU PLAN DU TRIANGLE )
      COSA = F / SQRT(E*D)

      IF( ABS(COSA) .LT. 1D-4 ) THEN

C        S1-S2 EST PARALLELE AU PLAN P1 P2 P3
C        COSA=1D-5 => A=89.999427 degres
C        COSA=1D-4 => A=89.994270 degres
C        COSA=1D-3 => A=89.942704 degres
         LINTER = -1

C        S1-S2 EST IL DANS LE PLAN P1 P2 P3?
C        I.E. S1 EST IL SUR LE PLAN?
         S = 0D0
         GOTO 10

      ENDIF

C     LE POINT D'INTERSECTION DE L'ARETE S1-S2 AVEC LE PLAN
      S = ( (P1(1)-S1(1))*A + (P1(2)-S1(2))*B + (P1(3)-S1(3))*C ) / F

 10   PTINT(1) = S1(1) + S * ( S2(1) - S1(1) )
      PTINT(2) = S1(2) + S * ( S2(2) - S1(2) )
      PTINT(3) = S1(3) + S * ( S2(3) - S1(3) )
C
C     CE POINT PTINT EST IL INTERNE AU TRIANGLE?
C     CALCUL DES 3 COORDONNEES BARYCENTRIQUES DE PTINT

C     TRIANGLE PTINT-P2 PTINT-P3
      A = ( P2(2) - PTINT(2) ) * ( P3(3) - PTINT(3) )
     %  - ( P2(3) - PTINT(3) ) * ( P3(2) - PTINT(2) )
      B = ( P2(3) - PTINT(3) ) * ( P3(1) - PTINT(1) )
     %  - ( P2(1) - PTINT(1) ) * ( P3(3) - PTINT(3) )
      C = ( P2(1) - PTINT(1) ) * ( P3(2) - PTINT(2) )
     %  - ( P2(2) - PTINT(2) ) * ( P3(1) - PTINT(1) )
      CBT(1) = SQRT( ( A * A + B * B + C * C ) / D )

C     TRIANGLE PTINT-P3 PTINT-P1
      A = ( P3(2) - PTINT(2) ) * ( P1(3) - PTINT(3) )
     %  - ( P3(3) - PTINT(3) ) * ( P1(2) - PTINT(2) )
      B = ( P3(3) - PTINT(3) ) * ( P1(1) - PTINT(1) )
     %  - ( P3(1) - PTINT(1) ) * ( P1(3) - PTINT(3) )
      C = ( P3(1) - PTINT(1) ) * ( P1(2) - PTINT(2) )
     %  - ( P3(2) - PTINT(2) ) * ( P1(1) - PTINT(1) )
      CBT(2) = SQRT( ( A * A + B * B + C * C ) / D )

C     TRIANGLE PTINT-P1 PTINT-P2
      A = ( P1(2) - PTINT(2) ) * ( P2(3) - PTINT(3) )
     %  - ( P1(3) - PTINT(3) ) * ( P2(2) - PTINT(2) )
      B = ( P1(3) - PTINT(3) ) * ( P2(1) - PTINT(1) )
     %  - ( P1(1) - PTINT(1) ) * ( P2(3) - PTINT(3) )
      C = ( P1(1) - PTINT(1) ) * ( P2(2) - PTINT(2) )
     %  - ( P1(2) - PTINT(2) ) * ( P2(1) - PTINT(1) )
      CBT(3) = SQRT( ( A * A + B * B + C * C ) / D )

      IF( LINTER .EQ. -1 ) THEN

C        S1-S2 EST PARALLELE AU PLAN P1 P2 P3
C        S1 EST IL UN POINT DU PLAN P1P2P3?
C        CALCUL DE LA DISTANCE DE S1 AU PLAN
         CALL DIPTPL( S1, P1, P2, P3, CBT123 )

         IF( CBT123 .LE. 0.00001D0 * SQRT(D) ) THEN
C           S1 EST DANS LE PLAN -> S1-S2 EST DANS LE PLAN P1P2P3
            LINTER = -2
         ENDIF
         GOTO 9999

      ENDIF

C     COORDONNEE BARYCENTRIQUE DE PTINT SUR L'ARETE S1-S2
C     AJOUT 29/09/2014 DU TEST SUR PTINT INTERIEUR A S1-S2
      N = 0
      DO K=1,3
         IF( S2(K) .NE. S1(K) ) THEN
            CBA(K) = ( PTINT(K) - S1(K) ) / ( S2(K) - S1(K) )
            N = N + 1
         ELSE
            CBA(K) = 0D0
         ENDIF
      ENDDO
      CBA(1) = ( CBA(1) + CBA(2) + CBA(3) ) / N

      CBT123 = ABS(CBT(1)) + ABS(CBT(2)) + ABS(CBT(3))

      IF( CBT123 .LE. 1.00001D0  .AND.
     %    0D0 .LE. CBA(1) .AND. CBA(1) .LE. 1D0 ) THEN
C        LE POINT EST DANS LE TRIANGLE OU SUR UNE DE SES 3 ARETES
C        ENTRE LES SOMMETS S1-S2
         LINTER = 1
      ELSE
         LINTER = 0
      ENDIF

 9999 RETURN
      END
