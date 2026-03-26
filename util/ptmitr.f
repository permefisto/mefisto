      SUBROUTINE PTMITR( NBTRIA, NOSOTR, XYZSOM, P, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES COORDONNEES DU POINT P MINIMISANT LA SOMME
C -----    DES CARRES DES DISTANCES AUX NBTRIA TRIANGLES
C
C ENTREES:
C --------
C NBTRIA : NOMBRE DE TRIANGLES
C NOSOTR : NOSOTR(I,J) NUMERO XYZSOM DU I-EME SOMMET DU TRIANGLE J
C XYZSOM : COORDONNEES X Y Z DES SOMMETS DES TRIANGLES
C
C SORTIES:
C --------
C P      : 3 COORDONNEES DU POINT MINIMISANT
C IERR   : 0 PAS D ERREUR , 1 SI 3 SOMMETS D'UN TRIANGLE SONT ALIGNES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JANVIER 1992
C....................................................................012
      INTEGER           NOSOTR(3,NBTRIA)
      REAL              XYZSOM(3,*),P(3)
      REAL              COORPO(3,3),A(4)
      DOUBLE PRECISION  MAT(3,4),DEN,A1,A2,A3,A4
      INTRINSIC         REAL
C
C     MISE A ZERO DE LA MATRICE
      CALL AZEROD( 12, MAT )
C
      DO 50 J=1,NBTRIA
C
C        COEFFICIENTS DE L'EQUATION DU PLAN DU TRIANGLE
         DO 20 I=1,3
            DO 10 K=1,3
               COORPO(K,I) = XYZSOM(K,NOSOTR(I,J))
 10         CONTINUE
 20      CONTINUE
         CALL EQPLAN( COORPO, A, IERR )
         IF( IERR .GT. 0 ) RETURN
C
C        CONTRIBUTION DU TRIANGLE AU SYSTEME
         DEN = ( A(1)**2 + A(2)**2 + A(3)**2 )
         IF( DEN .LE. 0D0 ) THEN
            IERR = 1
            RETURN
         ENDIF
         DEN = 1.0 / DEN
         A1  = A(1) * DEN
         A2  = A(2) * DEN
         A3  = A(3) * DEN
         A4  =-A(4) * DEN
C
         MAT(1,1) = MAT(1,1) + A1 * A(1)
         MAT(1,2) = MAT(1,2) + A2 * A(1)
         MAT(1,3) = MAT(1,3) + A3 * A(1)
         MAT(1,4) = MAT(1,4) + A4 * A(1)
C
         MAT(2,1) = MAT(2,1) + A1 * A(2)
         MAT(2,2) = MAT(2,2) + A2 * A(2)
         MAT(2,3) = MAT(2,3) + A3 * A(2)
         MAT(2,4) = MAT(2,4) + A4 * A(2)
C
         MAT(3,1) = MAT(3,1) + A1 * A(3)
         MAT(3,2) = MAT(3,2) + A2 * A(3)
         MAT(3,3) = MAT(3,3) + A3 * A(3)
         MAT(3,4) = MAT(3,4) + A4 * A(3)
 50   CONTINUE
C
C     RESOLUTION DU SYSTEME
      CALL GAUSPT( 3, 1, MAT, IERR )
      IF( IERR .NE. 0 ) RETURN
C
      P(1) = REAL( MAT(1,4) )
      P(2) = REAL( MAT(2,4) )
      P(3) = REAL( MAT(3,4) )
      IERR = 0
C
      RETURN
      END
