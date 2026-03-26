      REAL FUNCTION TR2MED( P1, P2, P3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER LA MEDIANE DU SOMMET P1 DU TRIANGLE P1 P2 P3
C -----
C
C ENTREES :
C ---------
C P1,P2,P3 : LES 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C               SENS DIRECT POUR UNE SURFACE >0
C SORTIES :
C ---------
C TR2MED  : MEDIANE ENTRE LE SOMMET P1 ET LE COTE P2 P3
C           0 SI P1=P2 ou P1=P3
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1992
C2345X7..............................................................012
      REAL   P1(2),P2(2),P3(2)
C
C     LA MEDIANE P1 => P2+P3 / 2
      TR2MED = SQRT( ( (P2(1)+P3(1))*0.5 - P1(1) )**2
     %       +       ( (P2(2)+P3(2))*0.5 - P1(2) )**2 )
      END
