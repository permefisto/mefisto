      SUBROUTINE CBPTTRR( P1, P2, P3, PT, CBPT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCULER LES 3 COORDONNEES BARYCENTRIQUES DU POINT PT
C -----     DANS LE TRIANGLE DE SOMMETS P1 P2 P3 DE R3
C           VERSION AVEC DES REELS SIMPLE PRECISION

C ENTREES :
C ---------
C P1,P2,P3: LES 3 SOMMETS DU TRIANGLE
C PT      : LE POINT DE COORDONNEES BARYCENTRIQUES A CALCULER
C
C SORTIE :
C --------
C CBPT   : LES 3 COORDONNEES BARYCENTRIQUE DU POINT PT DANS LE TRIANGLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain  Saint Pierre du Perray             Mars 2020
C2345X7..............................................................012
      REAL               P1(3),  P2(3),  P3(3),  PT(3)
      DOUBLE PRECISION  DP1(3), DP2(3), DP3(3), DPT(3), CBPT(3)

C     PASSAGE REAL -> DOUBLE PRECISION
      DO K=1,3
         DP1( K ) = DBLE( P1( K ) )
         DP2( K ) = DBLE( P2( K ) )
         DP3( K ) = DBLE( P3( K ) )
         DPT( K ) = DBLE( PT( K ) )
      ENDDO

C     CALCULS EN DOUBLE PRECISION
      CALL CBPTTR( DP1, DP2, DP3, DPT, CBPT )

      RETURN
      END
