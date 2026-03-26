      SUBROUTINE LON6AR( P1, P2, P3, P4, LONARE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DE LA LONGUEUR DES 6 ARETES DU TETRAEDRE DE
C -----   SOMMETS P1 P2 P3 P4

C ENTREES:
C --------
C P1,P2,P3,P4 : LES 3 COORDONNEES DES 4 SOMMETS DU TETRAEDRE

C SORTIE :
C --------
C LONARE : LONGUEUR DES 6 ARETES DU TETRAEDRE DE SOMMETS P1P2P3P4
C          DANS L'ORDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS UPMC     JANVIER 1992
C2345X7..............................................................012
      REAL           P1(3),P2(3),P3(3),P4(3),LONARE(6)

      DISTAN(X,Y,Z)= SQRT( X*X + Y*Y + Z*Z )

C     LONGUEUR DE L'ARETE 12
      LONARE(1) = DISTAN( P2(1)-P1(1), P2(2)-P1(2), P2(3)-P1(3) )

C     LONGUEUR DE L'ARETE 23
      LONARE(2) = DISTAN( P3(1)-P2(1), P3(2)-P2(2), P3(3)-P2(3) )

C     LONGUEUR DE L'ARETE 31
      LONARE(3) = DISTAN( P3(1)-P1(1), P3(2)-P1(2), P3(3)-P1(3) )

C     LONGUEUR DE L'ARETE 41
      LONARE(4) = DISTAN( P4(1)-P1(1), P4(2)-P1(2), P4(3)-P1(3) )

C     LONGUEUR DE L'ARETE 42
      LONARE(5) = DISTAN( P4(1)-P2(1), P4(2)-P2(2), P4(3)-P2(3) )

C     LONGUEUR DE L'ARETE 43
      LONARE(6) = DISTAN( P4(1)-P3(1), P4(2)-P3(2), P4(3)-P3(3) )

      RETURN
      END
