      SUBROUTINE CH3D3R( XYZ1 , D2D3 , XY , XYZ )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : EFFECTUER LE CHANGEMENT DE BASE DEFINI PAR LA MATRICE D2D3
C -----
C  (X       )    (          ) ( XX )    ( XX ) ( T         )( X      )
C  (Y  - P1 ) =  (  D2D3    ) ( YY ) ;  ( YY )=(   D2D3    )( Y - P1 )
C  (Z       )    (          ) ( ZZ )    ( ZZ ) (           )( Z      )
C
C     X  , Y  , Z  COORDONNEES DANS LE REPERE 3D
C     XX , YY , ZZ COORDONNEES DANS LE REPERE 2D DU PLAN P1 P2 P3
C
C ENTREES :
C ---------
C XYZ1   : LE POINT P1 ORIGINE DU SYSTEME D'AXES DU PLAN
C D2D3   : LA MATRICE DE PASSAGE DEFINIE CI DESSUS
C XY     : LES 3 COORDONNEES DU POINT DANS LE SYSTEME D'AXES NOUVEAU
C          ( XX , YY , ZZ )
C
C SORTIE :
C --------
C XYZ    : LES 3 COORDONNEES DANS LE SYSTEME D'AXES 3D ANCIEN( X Y Z )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  SEPTEMBRE 1990
C........................................................................
      DOUBLE PRECISION  D2D3(3,3),D
      REAL              XYZ1(3),XYZ(3),XY(3)
      INTRINSIC         REAL
C
C     MULTIPLICATION PAR LA MATRICE D2D3 AVEC ZZ = 0.
      DO 20 I=1,3
         D = 0.D0
         DO 10 J=1,3
            D = D + D2D3( I , J ) * XY( J )
 10      CONTINUE
         XYZ( I ) = REAL( D + XYZ1( I ) )
 20   CONTINUE
C
      RETURN
      END
