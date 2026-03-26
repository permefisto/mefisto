      SUBROUTINE CR3D3D( XYZ0, D2D3, XYZN, XYZA )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  EFFECTUER LE CHANGEMENT DE BASE DEFINI PAR LA MATRICE D2D3
C -----
C  (X       )    (          ) ( XX )    ( XX ) ( T         )( X      )
C  (Y  - P0 ) =  (  D2D3    ) ( YY ) ;  ( YY )=(   D2D3    )( Y - P0 )
C  (Z       )    (          ) ( ZZ )    ( ZZ ) (           )( Z      )
C
C     X  , Y  , Z  COORDONNEES DANS LE REPERE 3D ANCIEN
C     XX , YY , ZZ COORDONNEES DANS LE REPERE 3D NOUVEAU
C
C ENTREES :
C ---------
C XYZ0   : LE POINT P1 ORIGINE DU SYSTEME D'AXES DU PLAN
C D2D3   : LA MATRICE DE PASSAGE DEFINIE CI DESSUS
C XYZN   : LES 3 COORDONNEES DU POINT DANS LE SYSTEME D'AXES 3D NOUVEAU
C          ( XX , YY , ZZ )
C
C SORTIE :
C --------
C XYZA   : LES 3 COORDONNEES DANS LE SYSTEME D'AXES 3D ANCIEN( X Y Z )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire JL LIONS UPMC Paris Decembre 2006
C....................................................................012
      DOUBLE PRECISION D2D3(3,3),D
      REAL             XYZ0(3), XYZA(3), XYZN(3)
      INTRINSIC        REAL
C
C     MULTIPLICATION PAR LA MATRICE D2D3
      DO 20 I=1,3
         D = 0.D0
         DO 10 J=1,3
            D = D + D2D3( I , J ) * XYZN( J )
 10      CONTINUE
         XYZA( I ) = REAL( D + XYZ0( I ) )
 20   CONTINUE
C
      RETURN
      END
