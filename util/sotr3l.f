      SUBROUTINE SOTR3L( RLONGC, X, Y, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES 2 COORDONNEES DU SOMMET SUPERIEUR DU TRIANGLE
C -----    DEFINI PAR LA LONGUEUR DE SES 3 COTES
C
C ENTREES:
C --------
C RLONGC : LONGUEUR DES 3 ARETES DES 3 COTES DU TRIANGLE
C          RLONGC(1) >= RLONGC(2) ET RLONGC(1) >= RLONGC(3)
C
C SORTIES:
C --------
C X, Y   : COORDONNEES DU SOMMET SUPERIEUR DU TRIANGLE
C          DE SOMMET 1 L'ORIGINE, DE SOMMET 2 (RLONGC(1),0)
C IERR   : 0 SI LE TRIANGLE EST CONSTRUCTIBLE
C          1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS       SEPTEMBRE 1993
C234567..............................................................012
      REAL  RLONGC( 3 )
C
      X = (RLONGC(1)**2+RLONGC(3)**2-RLONGC(2)**2) / (2.0*RLONGC(1))
      Y = RLONGC(3)**2 - X**2
      IF( Y .LE. 0 ) THEN
C        IMPOSSIBILITE DE CONSTRUIRE LE TRIANGLE
         IERR = 1
      ELSE
         Y = SQRT( RLONGC(3)**2 - X**2 )
         IERR = 0
      ENDIF
      END
