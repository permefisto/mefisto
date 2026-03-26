      SUBROUTINE PRPTPLD( PT, PT1, PT2, PT3, PTPROJ, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE POINT PROJECTION PTPROJ DU POINT PT
C -----    SUR LE PLAN DEFINI PAR LES 3 POINTS PT1, PT2, PT3

C ENTREES:
C --------
C PT     : X Y Z DU POINT
C PT1, PT2, PT3 : LES 3 POINTS DU PLAN OU 3 SOMMETS D'UN TRIANGLE

C SORTIES:
C --------
C PTPROJ : LES 3 COORDONNEES DU POINT DE PROJECTION
C IERR   : =0 SI POINT PROJETE CORRECTEMENT CALCULE
C          >0 SINON  ( 3 POINTS ALIGNES )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1993
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Janvier 2017
C2345X7..............................................................012
      DOUBLE PRECISION  PT(3),PT1(3),PT2(3),PT3(3),PTPROJ(3),VN(3)
      DOUBLE PRECISION  A(3,4), PROSCD

C     LE VECTEUR NORMAL AU PLAN NON NORMALISE A 1
      CALL VECNOR3D( PT1, PT2, PT3, VN )

      IF( PROSCD( VN, VN, 3 ) .LE. 0D0 ) THEN
         IERR = 1
         GOTO 9999
      ENDIF

C     FORMATION DE LA MATRICE ET SECOND MEMBRE DU SYSTEME LINEAIRE
C     PT  PTPROJ . PT1 PT2 = 0
C     PT  PTPROJ . PT1 PT3 = 0
C     PT1 PTPROJ . (NORMALE AU PLAN PT1 PT2 PT3 ) = 0
      A(1,1) = PT2(1) - PT1(1)
      A(1,2) = PT2(2) - PT1(2)
      A(1,3) = PT2(3) - PT1(3)

      A(2,1) = PT3(1) - PT1(1)
      A(2,2) = PT3(2) - PT1(2)
      A(2,3) = PT3(3) - PT1(3)

      A(3,1) = VN(1)
      A(3,2) = VN(2)
      A(3,3) = VN(3)

      A(1,4) = PT(1)  * A(1,1) + PT(2)  * A(1,2) + PT(3)  * A(1,3)
      A(2,4) = PT(1)  * A(2,1) + PT(2)  * A(2,2) + PT(3)  * A(2,3)
      A(3,4) = PT1(1) * A(3,1) + PT1(2) * A(3,2) + PT1(3) * A(3,3)

C     RESOLUTION DU SYSTEME LINEAIRE
      CALL GAUSPT( 3, 1, A, IERR )
      IF( IERR .NE. 0 ) GOTO 9999

C     LE POINT PROJETE
      PTPROJ(1) = A(1,4)
      PTPROJ(2) = A(2,4)
      PTPROJ(3) = A(3,4)

 9999 RETURN
      END
