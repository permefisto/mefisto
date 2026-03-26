      SUBROUTINE CBPTTR( P1, P2, P3, PT, CBPT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCULER LES 3 COORDONNEES BARYCENTRIQUES DU POINT PT
C -----     DANS LE TRIANGLE DE SOMMETS P1 P2 P3 DE R3
C
C ENTREES :
C ---------
C P1,P2,P3: LES 3 SOMMETS DU TRIANGLE
C PT      : LE POINT DE COORDONNEES BARYCENTRIQUES A CALCULER
C
C SORTIE :
C --------
C CBPT   : LES 3 COORDONNEES BARYCENTRIQUE DU POINT PT DANS LE TRIANGLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain LJLL UPMC& St Pierre du Perray SEPTEMBRE 2011
C2345X7..............................................................012
      DOUBLE PRECISION  P1(3), P2(3), P3(3), PT(3), CBPT(3)
      DOUBLE PRECISION  A, B, C, XN, YN, ZN, DN
C
C     LE VECTEUR PRODUIT VECTORIEL S1S2 x S1S3
      XN = ( P2(2) - P1(2) ) * ( P3(3) - P1(3) )
     %   - ( P2(3) - P1(3) ) * ( P3(2) - P1(2) )
      YN = ( P2(3) - P1(3) ) * ( P3(1) - P1(1) )
     %   - ( P2(1) - P1(1) ) * ( P3(3) - P1(3) )
      ZN = ( P2(1) - P1(1) ) * ( P3(2) - P1(2) )
     %   - ( P2(2) - P1(2) ) * ( P3(1) - P1(1) )
C     LE CARRE DU MODULE DU VECTEUR NORMAL ou 2 FOIS LA SURFACE
      DN = XN * XN + YN * YN + ZN * ZN
C
C     LE VECTEUR PRODUIT VECTORIEL PTS2 x PTS3 du TRIANGLE PT-P2 PT-P3
      A = ( P2(2) - PT(2) ) * ( P3(3) - PT(3) )
     %  - ( P2(3) - PT(3) ) * ( P3(2) - PT(2) )
      B = ( P2(3) - PT(3) ) * ( P3(1) - PT(1) )
     %  - ( P2(1) - PT(1) ) * ( P3(3) - PT(3) )
      C = ( P2(1) - PT(1) ) * ( P3(2) - PT(2) )
     %  - ( P2(2) - PT(2) ) * ( P3(1) - PT(1) )
      CBPT(1) = SQRT( ( A * A + B * B + C * C ) / DN )
C     SIGNE DU PRODUIT SCALAIRE (XN,YN,ZN) . (A,B,C)
      IF( XN*A + YN*B + ZN*C .LT. 0D0 ) CBPT(1) = - CBPT(1)
C
C     LE VECTEUR PRODUIT VECTORIEL PTS3 x PTS1 du TRIANGLE PT-P3 PT-P1
      A = ( P3(2) - PT(2) ) * ( P1(3) - PT(3) )
     %  - ( P3(3) - PT(3) ) * ( P1(2) - PT(2) )
      B = ( P3(3) - PT(3) ) * ( P1(1) - PT(1) )
     %  - ( P3(1) - PT(1) ) * ( P1(3) - PT(3) )
      C = ( P3(1) - PT(1) ) * ( P1(2) - PT(2) )
     %  - ( P3(2) - PT(2) ) * ( P1(1) - PT(1) )
      CBPT(2) = SQRT( ( A * A + B * B + C * C ) / DN )
C     SIGNE DU PRODUIT SCALAIRE (XN,YN,ZN) . (A,B,C)
      IF( XN*A + YN*B + ZN*C .LT. 0D0 ) CBPT(2) = - CBPT(2)
C
C     PRODUIT VECTORIEL PTS1 x PTS2 du TRIANGLE PT-P1 PT-P2
      A = ( P1(2) - PT(2) ) * ( P2(3) - PT(3) )
     %  - ( P1(3) - PT(3) ) * ( P2(2) - PT(2) )
      B = ( P1(3) - PT(3) ) * ( P2(1) - PT(1) )
     %  - ( P1(1) - PT(1) ) * ( P2(3) - PT(3) )
      C = ( P1(1) - PT(1) ) * ( P2(2) - PT(2) )
     %  - ( P1(2) - PT(2) ) * ( P2(1) - PT(1) )
      CBPT(3) = SQRT( ( A * A + B * B + C * C ) / DN )
C     SIGNE DU PRODUIT SCALAIRE (XN,YN,ZN) . (A,B,C)
      IF( XN*A + YN*B + ZN*C .LT. 0D0 ) CBPT(3) = - CBPT(3)
C
      RETURN
      END
