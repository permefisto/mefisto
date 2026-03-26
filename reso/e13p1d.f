      SUBROUTINE E13P1D( XYZ, XYZBAR, DELTA, DFM1, DP )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LA VALEUR des XYZ du BARYCENTRE, DELTA, DFM1, DP
C ----- DU TETRAEDRE COURANT DEFINI PAR LES 3 COORDONNEES DE SES 4 SOMMETS

C ENTREE :
C --------
C XYZ    : 3 COORDONNEES DES 4 SOMMETS DU TETRAEDRE COURANT

C SORTIES:
C --------
C XYZBAR : 3 COORDONNEES XYZ DU BARYCENTRE DU TETRAEDRE DE SOMMETS XYZ
C DELTA  : DETERMINANT DE LA MATRICE DF JACOBIENNE
C DFM1   : LA MATRICE JACOBIENNE INVERSE  3 x 3
C DP     : GRADIENT (CONSTANT) DES 4 POLYNOMES DE BASE P1 SUR LE
C          TETRAEDRE A MULTIPLIER PAR LES 4 VALEURS AUX 4 SOMMETS
C          POUR OBTENIR LE GRADIENT SUR LE TETRAEDRE P1
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C23456---------------------------------------------------------------012
      REAL              XYZ(4,3)
      DOUBLE PRECISION  XYZBAR(3), DELTA, DFM1(3,3), DP(3,4), DF(3,3)

C     LE BARYCENTRE DU TETRAEDRE
      XYZBAR(1) = ( XYZ(1,1) + XYZ(2,1) + XYZ(3,1) + XYZ(4,1) ) * 0.25D0
      XYZBAR(2) = ( XYZ(1,2) + XYZ(2,2) + XYZ(3,2) + XYZ(4,2) ) * 0.25D0
      XYZBAR(3) = ( XYZ(1,3) + XYZ(2,3) + XYZ(3,3) + XYZ(4,3) ) * 0.25D0

C     LA MATRICE JACOBIENNE
      DF(1,1) = XYZ(2,1) - XYZ(1,1)
      DF(1,2) = XYZ(2,2) - XYZ(1,2)
      DF(1,3) = XYZ(2,3) - XYZ(1,3)

      DF(2,1) = XYZ(3,1) - XYZ(1,1)
      DF(2,2) = XYZ(3,2) - XYZ(1,2)
      DF(2,3) = XYZ(3,3) - XYZ(1,3)

      DF(3,1) = XYZ(4,1) - XYZ(1,1)
      DF(3,2) = XYZ(4,2) - XYZ(1,2)
      DF(3,3) = XYZ(4,3) - XYZ(1,3)

C     [DF]-1
      CALL M33INV( DF, DELTA, DFM1 )

C     LE GRADIENT  DP(i,j) = d Pj / dxi
C     ou Pj = j-eme COORDONNEE BARYCENTRIQUE sur le tetraedre de base
      DP(1,1) = - DFM1(1,1) - DFM1(1,2) - DFM1(1,3)
      DP(1,2) =   DFM1(1,1)
      DP(1,3) =   DFM1(1,2)
      DP(1,4) =   DFM1(1,3)

      DP(2,1) = - DFM1(2,1) - DFM1(2,2) - DFM1(2,3)
      DP(2,2) =   DFM1(2,1)
      DP(2,3) =   DFM1(2,2)
      DP(2,4) =   DFM1(2,3)

      DP(3,1) = - DFM1(3,1) - DFM1(3,2) - DFM1(3,3)
      DP(3,2) =   DFM1(3,1)
      DP(3,3) =   DFM1(3,2)
      DP(3,4) =   DFM1(3,3)

      RETURN
      END
