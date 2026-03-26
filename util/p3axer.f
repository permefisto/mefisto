      SUBROUTINE P3AXER( POINT, AXE, IERR  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : A PARTIR DE 3 POINTS FORMER LES AXES ORTHONORMES TELS QUE
C ----- AXE( . , 1 ) = P1P2 / NORME(  P1P2  )
C       AXE( . , 3 ) = AXE( . , 1 ) PRODUIT VECTORIEL P1P3 / NORME(P1P3)
C       AXE( . , 2 ) = AXE( . , 3 ) PRODUIT VECTORIEL AXE( . , 1 )
C
C PARAMETRE DONNEE :
C ------------------
C POINT  : POINT( I , J )=I-EME COORDONNEE DU J-EME POINT
C
C PARAMETRE RESULTAT :
C --------------------
C AXE   : AXE( I , J )=I-EME COORDONNE DU J-EME AXE ORTHONORME
C IERR  : 1 SI UN DES 3 AXES A POUR NORME ZERO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS         MAI 1990
C ......................................................................
      DOUBLE PRECISION AXE(3,*)
      REAL             POINT(3,*)
C
C     AXE( . , 1 )
C     ============
      DO 10 I=1 , 3
           AXE( I , 1 ) = POINT( I , 2 ) - POINT( I , 1 )
           AXE( I , 2 ) = POINT( I , 3 ) - POINT( I , 1 )
   10 CONTINUE
      CALL NORME1( 3, AXE( 1 , 1 ), IERR )
      IF( IERR .NE. 0 ) RETURN
      CALL NORME1( 3, AXE( 1 , 2 ), IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     AXE( . , 3 )
C     ============
      CALL PROVEC( AXE( 1 , 1 ) , AXE( 1 , 2 ) , AXE( 1 , 3 ) )
      CALL NORME1( 3, AXE( 1 , 3 ), IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     AXE( . , 2 )
C     ============
      CALL PROVEC( AXE( 1 , 3 ) , AXE( 1 , 1 ) , AXE( 1 , 2 ) )
      CALL NORME1( 3, AXE( 1 , 2 ), IERR )
      END
