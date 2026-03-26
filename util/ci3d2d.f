      SUBROUTINE CI3D2D( XYZ1 , D2D3 , XYZ , XY )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : EFFECTUER LE CHANGEMENT DE BASE DEFINI PAR LA MATRICE D2D3
C -----
C  (X       )    (          ) ( XX )    ( XX ) ( T         )( X      )
C  (Y  - P1 ) =  (  D2D3    ) ( YY ) ;  ( YY )=(   D2D3    )( Y - P1 )
C  (Z       )    (          ) ( ZZ )    ( ZZ ) (           )( Z      )
C
C          AVEC DANS LES APPLICATIONS ZZ = 0.
C     X  , Y  , Z  COORDONNEES DANS LE REPERE 3D
C     XX , YY , ZZ COORDONNEES DANS LE REPERE 2D DU PLAN P1 P2 P3
C
C ENTREES :
C ---------
C XYZ1   : LE POINT P1 ORIGINE DU SYSTEME D'AXES DU PLAN
C D2D3   : LA MATRICE DE PASSAGE DEFINIE CI DESSUS
C XYZ    : LES 3 COORDONNEES DU POINT DANS LE SYSTEME D'AXES 3D ANCIEN
C          ( X , Y , Z )
C
C SORTIE :
C --------
C XY     : LES 2 COORDONNEES DANS LE SYSTEME D'AXES 2D NOUVEAU(XX YY 0.)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  JANVIER 1986
C........................................................................
      DOUBLE PRECISION D2D3(3,3),D,XYZ1(3),XYZ(3),XY(2)
C
C     MULTIPLICATION PAR LA TRANSPOSEE DE D2D3  ZZ = 0.
      DO 20 I=1,2
         D = 0.D0
         DO 10 J=1,3
            D = D + D2D3(J,I) * ( XYZ(J) - XYZ1(J) )
 10      CONTINUE
         XY( I ) = D
 20   CONTINUE
      END
