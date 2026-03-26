      SUBROUTINE QUADCXR( P1, P2, P3, P4, NONOUI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :      LE QUADRANGLE P1 P2 P3 P4 EST IL CONVEXE ?
C -----      LES 2 DECOUPAGES EN 2 TRIANGLES SONT ILS POSSIBLES?

C ENTREES:
C --------
C P1,P2,P3,P4 : LES 3 COORDONNEES DES 4 SOMMETS DU QUADRANGLE

C SORTIE :
C --------
C NONOUI : 1 LE QUADRANGLE EST CONVEXE
C          0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS UPMC     FEVRIER 1993
C2345X7..............................................................012
      REAL               P1(3),  P2(3),  P3(3),  P4(3)
      DOUBLE PRECISION  DP1(3), DP2(3), DP3(3), DP4(3)

C     PASSAGE REEL -> DOUBLE PRECISION DES XYZ DES 4 SOMMETS
      DO K = 1, 3
         DP1(K) = P1(K)
         DP2(K) = P2(K)
         DP3(K) = P3(K)
         DP4(K) = P4(K)
      ENDDO

      CALL QUADCXD( DP1, DP2, DP3, DP4, NONOUI )

      RETURN
      END
