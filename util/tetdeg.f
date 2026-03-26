       LOGICAL FUNCTION TETDEG( P1 , P2 , P3 , P4 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LE TETRAEDRE P1 P2 P3 P4 EST IL DEGENERE ?
C -----
C
C ENTREES:
C --------
C P1 P2 P3 P4 : LES 3 COORDONNEES DES 4 SOMMETS
C
C SORTIE :
C --------
C TETDEG : .TRUE. OU .FALSE. SELON DEGENERESCENCE OU NON DU TETRAEDRE
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS OCTOBRE 1987
C...........................................................................
      INTEGER           P1(3),P2(3),P3(3),P4(3)
      DOUBLE PRECISION  PP(9),D,V,VOLTET
C
C     FORMATION DU SYSTEME P(I)-P(1)   I=2,3,4
      I = 1
      DO 10 J=1,3
C
C              RANGEMENT              VALEURS
C             ( 1  4  7 )     ( X2-X1  Y2-Y1  Z2-Z1 )
C        PP = ( 2  5  8 )   = ( X3-X1  Y3-Y1  Z3-Z1 )
C             ( 3  6  9 )     ( X4-X1  Y4-Y1  Z4-Z1 )
C
         D           = DBLE( P1(J) )
         PP( I     ) = DBLE( P2(J) ) - D
         PP( I + 1 ) = DBLE( P3(J) ) - D
         PP( I + 2 ) = DBLE( P4(J) ) - D
         I           = I + 3
 10   CONTINUE
C
C     CALCUL DU DETERMINANT DE LA MATRICE OU 6 * VOLUME DU TETRAEDRE
      VOLTET = PP(1) * ( PP(5) * PP(9) - PP(6) * PP(8) )
     %       + PP(2) * ( PP(6) * PP(7) - PP(4) * PP(9) )
     %       + PP(3) * ( PP(4) * PP(8) - PP(5) * PP(7) )
C
C     LE TETRAEDRE EST IL DEGENERE ?
C     V= ( P1P2 ** 2 ) * ( P1P3 ** 2 ) * ( P1P4 ** 2 )
      V = ( PP(1) ** 2 + PP(4) ** 2 + PP(7) ** 2 )
     %  * ( PP(2) ** 2 + PP(5) ** 2 + PP(8) ** 2 )
     %  * ( PP(3) ** 2 + PP(6) ** 2 + PP(9) ** 2 )
C
      TETDEG =  VOLTET * VOLTET  .LE.  1D-9 * V
      END
