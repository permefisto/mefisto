      SUBROUTINE CV3D3D( D2D3 , XYZ , XY3 )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : EFFECTUER LE CHANGEMENT DE BASE DEFINI PAR LA MATRICE D2D3
C ----- POUR UN VECTEUR (ET NON PAS UN POINT cf CH3D3D POUR CELA!)
C  ( X )    (        ) ( XX )    ( XX ) ( T         )( X )
C  ( Y ) =  (  D2D3  ) ( YY ) ;  ( YY )=(   D2D3    )( Y )
C  ( Z )    (        ) ( ZZ )    ( ZZ ) (           )( Z )
C
C     X  , Y  , Z  COORDONNEES DANS LE REPERE 3D
C     XX , YY , ZZ COORDONNEES DANS LE REPERE 3D DU PLAN P1 P2 P3 ET
C                                                DIRECTION ORTHOGONALE
C
C ENTREES :
C ---------
C XYZ1   : LE POINT P1 ORIGINE DU SYSTEME D'AXES DU PLAN
C D2D3   : LA MATRICE DE PASSAGE DEFINIE CI DESSUS
C XYZ    : LES 3 COORDONNEES DU VECTEUR DANS LE SYSTEME D'AXES 3D ANCIEN
C          ( X , Y , Z )
C
C SORTIE :
C --------
C XY3    : LES 3 COORDONNEES DANS LE SYSTEME D'AXES 3D NOUVEAU (XX YY ZZ)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS          MAI 1990
C........................................................................
      DOUBLE PRECISION D2D3(3,3),D
      REAL             XYZ(3), XY3(3)
C
C     MULTIPLICATION PAR LA TRANSPOSEE DE D2D3
      DO 20 I=1,3
         D = 0.D0
         DO 10 J=1,3
            D = D + D2D3(J,I) * XYZ(J)
 10      CONTINUE
         XY3( I ) = REAL( D )
 20   CONTINUE
C
      RETURN
      END
