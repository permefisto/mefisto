      FUNCTION FAECTR( P )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DU FACTEUR D'ECRASEMENT DU TRIANGLE DE SOMMETS
C ----- DE COORDONNEES P(1 A 2 , 1 A 3 )
C
C ENTREE :
C --------
C P      : P(I,J) I-EME COORDONNEE DU J-EME SOMMET
C
C SORTIE :
C --------
C FAECTR : MIN( HAUTEUR I / LONGUEUR COTE I )
C          I=1,2,3
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS NOVEMBRE 1986
C...............................................................................
      REAL  P(2,3)
C
C     BOUCLE SUR LES 3 ARETES DU TRIANGLE
C     ===================================
      FAECTR = 1.
      DO 10 I=1,3
C
C        LE SOMMET AVANT I
         I0 = I - 1
         IF( I0 .LE. 0 ) I0 = 3
C
C        LE SOMMET APRES I
         I2 = I + 1
         IF( I2 .GT. 3 ) I2 = 1
C
         X10 = P(1,I)  - P(1,I0)
         Y10 = P(2,I)  - P(2,I0)
         X20 = P(1,I2) - P(1,I0)
         Y20 = P(2,I2) - P(2,I0)
C
C        ECRASE : HAUTEUR / LONGUEUR DU COTE I
C               : ABS( P(I-1) P(I)  PRODUIT VECTORIEL P(I-1) P(I+1) )
C                  /    CARRE DE LA LONGUEUR DE P(I-1) P(I+1)
C
         ECRASE = ABS( X10 * Y20 - Y10 * X20 )
     %            /  ( X20 * X20 + Y20 * Y20 )
C
C        CALCUL DU MINIMUM SUR LES 3 COTES
         IF( ECRASE .LT. FAECTR ) THEN
            FAECTR = ECRASE
         ENDIF
 10   CONTINUE
      END
