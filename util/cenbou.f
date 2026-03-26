      SUBROUTINE CENBOU( NP1 , NP2 , NP3 , NP4 , LECENT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES 3 COORDONNEES DU CENTRE DE LA BOULE CIRCONSCRITE
C -----    AU TETRAEDRE NP1 NP2 NP3 NP4 ET DU CARRE DE SON RAYON
C
C ENTREES:
C --------
C NP1 NP2 NP3 NP4 : LES 3 COORDONNEES DES 4 SOMMETS
C
C SORTIE :
C --------
C LECENT : 3 COORDONNEES DU CENTRE ET RAYON
C          SI LE TETRAEDRE EST DEGENERE ALORS LECENT(1 2 3 4 )=0. EN SORTIE
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS JANVIER 1987
C...............................................................................
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NP1(3),NP2(3),NP3(3),NP4(3),LECENT(4)
C
      DOUBLE PRECISION  PP(9),BB(3),D,D1,D2,D3,V
C
C     FORMATION DU SYSTEME NP1 PI . LECENT = NP1 PI . (NP1+PI)/2
      BB(1) = 0D0
      BB(2) = 0D0
      BB(3) = 0D0
      I     = 1
C
      DO 10 J=1,3
C
C              RANGEMENT              VALEURS
C             ( 1  4  7 )     ( X2-X1  Y2-Y1  Z2-Z1 )
C        PP = ( 2  5  8 )   = ( X3-X1  Y3-Y1  Z3-Z1 )
C             ( 3  6  9 )     ( X4-X1  Y4-Y1  Z4-Z1 )
C
         D           = DBLE( NP1(J) )
         PP( I     ) = DBLE( NP2(J) ) - D
         PP( I + 1 ) = DBLE( NP3(J) ) - D
         PP( I + 2 ) = DBLE( NP4(J) ) - D
C
C        LE SECOND MEMBRE
         BB( 1 ) = BB( 1 ) + PP( I ) * ( NP2( J ) + D )
         BB( 2 ) = BB( 2 ) + PP(I+1) * ( NP3( J ) + D )
         BB( 3 ) = BB( 3 ) + PP(I+2) * ( NP4( J ) + D )
         I       = I + 3
 10   CONTINUE
C
C     LA MATRICE EST MULTIPLIEE PAR 2 POUR EVITER UNE DIVISION PAR 2
      DO 15 I = 1 , 9
         PP( I ) = PP( I ) + PP( I )
 15   CONTINUE
C
C     CALCUL DU DETERMINANT DE LA MATRICE OU 6 * VOLUME DU TETRAEDRE
      D = PP(1) * ( PP(5) * PP(9) - PP(6) * PP(8) )
     %  + PP(2) * ( PP(6) * PP(7) - PP(4) * PP(9) )
     %  + PP(3) * ( PP(4) * PP(8) - PP(5) * PP(7) )
C
C     LE TETRAEDRE EST IL DEGENERE ?
C     V= ( NP1NP2 ** 2 ) * ( NP1NP3 ** 2 ) * ( NP1NP4 ** 2 )
C
      V = 1D0
      DO 20 I=1,9,3
         V = V * ( PP(I) ** 2 + PP(I+1) ** 2 + PP(I+2) ** 2 )
 20   CONTINUE
      IF( D * D .LE. 1D-9 * V ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR CENBOU : TETRAEDRE DEGENERE'
         CALL LEREUR
C
         WRITE(IMPRIM,10020) (NP1(I),NP2(I),NP3(I),NP4(I),I=1,3)
10020    FORMAT(' X=',4I14/' Y=',4I14/' Z=',4I14/)
         LECENT(1) = 0
         LECENT(2) = 0
         LECENT(3) = 0
         LECENT(4) = 0
         RETURN
      ENDIF
C
C     LES 3 COORDONNEES DU CENTRE DE LA BOULE CIRCONSCRITE AU TETRAEDRE
C     LE PROBLEME DE LA PRECISION IMPOSE LES DOUBLE PRECISION
C     UNE / AU LIEU DE UNE / PLUS UNE *  => MOINS D'ERREUR D'ARRONDIS
      D1 = ( BB(1) * ( PP(5) * PP(9) - PP(6) * PP(8) )
     %     + BB(2) * ( PP(6) * PP(7) - PP(4) * PP(9) )
     %     + BB(3) * ( PP(4) * PP(8) - PP(5) * PP(7) ) ) / D
C
      D2 = ( BB(1) * ( PP(3) * PP(8) - PP(2) * PP(9) )
     %     + BB(2) * ( PP(1) * PP(9) - PP(3) * PP(7) )
     %     + BB(3) * ( PP(2) * PP(7) - PP(1) * PP(8) ) ) / D
C
      D3 = ( BB(1) * ( PP(2) * PP(6) - PP(3) * PP(5) )
     %     + BB(2) * ( PP(3) * PP(4) - PP(1) * PP(6) )
     %     + BB(3) * ( PP(1) * PP(5) - PP(2) * PP(4) ) ) / D
C
C     LES COORDONNEES DU CENTRE EN ENTIERS
      LECENT(1) = NINT( D1 )
      LECENT(2) = NINT( D2 )
      LECENT(3) = NINT( D3 )
C
C     LE CARRE DU RAYON DE LA BOULE CIRCONSCRITE
      D = ( NP1(1) - D1 ) ** 2
     %  + ( NP1(2) - D2 ) ** 2
     %  + ( NP1(3) - D3 ) ** 2
C
C     LE RAYON TRADUIT EN ENTIER
      D         = SQRT( D )
      LECENT(4) = NINT( D )
C
C     LE RAYON PAR RAPPORT A NP2 NP3 NP4
CCC      D  = ( NP2(1) - D1 ) ** 2
CCC     %   + ( NP2(2) - D2 ) ** 2
CCC     %   + ( NP2(3) - D3 ) ** 2
CCC      N = NINT( SQRT( D ) )
CCC      WRITE(IMPRIM,19001) 1,LECENT(4),2,N
CCC      D  = ( NP3(1) - D1 ) ** 2
CCC     %   + ( NP3(2) - D2 ) ** 2
CCC     %   + ( NP3(3) - D3 ) ** 2
CCC      N = NINT( SQRT( D ) )
CCC      WRITE(IMPRIM,19001) 3,N
CCC      D  = ( NP4(1) - D1 ) ** 2
CCC     %   + ( NP4(2) - D2 ) ** 2
CCC     %   + ( NP4(3) - D3 ) ** 2
CCC      N = NINT( SQRT( D ) )
CCC      WRITE(IMPRIM,19001) 4,N
CCC19001 FORMAT( ' RAYON',I1,'=',I14 )
C
C     AFFICHAGE DE LA PRECISION
CCC      CALL INTBOU( LECENT , NP2 , NP1 , D1 )
CCC      CALL INTBOU( LECENT , NP3 , NP1 , D2 )
CCC      CALL INTBOU( LECENT , NP4 , NP1 , D3 )
CCC      D = LECENT(4)
CCC      D = D * D
CCC      WRITE(IMPRIM,19000) D1/D , D2/D , D3/D
CCC19000 FORMAT(' PRECISION DANS CENBOU :',3G15.6 )
      END
