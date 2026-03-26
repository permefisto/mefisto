      SUBROUTINE AXOYZ( X,Y,Z,  YAX,ZAX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LES COORDONNEES YAX ZAX DANS LE PLAN D'AXONOMETRIE
C -----
C
C ENTREES :
C ---------
C X,Y,Z  : LES 3 COORDONNEES DU POINT
C
C SORTIES :
C ---------
C YAX,ZAX : LES COORDONNEES DANS LE PLAN D'AXONOMETRIE DU POINT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS        MAI 1990
C23456---------------------------------------------------------------012
      include"./incl/cm4vue.inc"
      INTRINSIC   REAL
C
      YAX = REAL( PTAXE3(1,2) * (X-PV(1)) + PTAXE3(2,2) * (Y-PV(2))
     %          + PTAXE3(3,2) * (Z-PV(3)) )
      ZAX = REAL( PTAXE3(1,3) * (X-PV(1)) + PTAXE3(2,3) * (Y-PV(2))
     %          + PTAXE3(3,3) * (Z-PV(3)) )
C
      RETURN
      END
