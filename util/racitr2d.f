      SUBROUTINE RACITR2D( P1, P2, P3, R, S )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCULER LE RAYON R DU CERCLE CIRCONSCRIT DU TRIANGLE PLAN
C -----     R = ABC/4S  ET LA SURFACE S DU TRIANGLE
C
C ENTREES :
C ---------
C P1,P2,P3: LES 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C
C SORTIES :
C ---------
C R       : RAYON DU CERCLE CIRCONSCRIT AU TRIANGLE
C           PLUS GRAND REEL MACHINE SI LE TRIANGLE EST DEGENERE (SURFACE=0)
C S       : LA SURFACE DU TRIANGLE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1992
C2345X7..............................................................012
      REAL   P1(2),P2(2),P3(2)
C
C     LA LONGUEUR DES 3 COTES
      A = SQRT( (P2(1)-P1(1))**2 + (P2(2)-P1(2))**2 )
      B = SQRT( (P3(1)-P2(1))**2 + (P3(2)-P2(2))**2 )
      C = SQRT( (P1(1)-P3(1))**2 + (P1(2)-P3(2))**2 )
C
C     LA SURFACE DU TRIANGLE 2D
      S = (P2(1)-P1(1)) * (P3(2)-P1(2)) - (P2(2)-P1(2)) * (P3(1)-P1(1))
      S = ABS( S ) / 2
C
      IF( S .EQ. 0 ) THEN
C
C        TRIANGLE DEGENERE
         R = RINFO( 'GRAND' )
C
      ELSE
C
C        RAYON DU CERCLE CIRCONSCRIT AU TRIANGLE  R = ABC/4S
         R = A * B * C / ( 4.0 * ABS(S) )
C
      ENDIF
C
      RETURN
      END
