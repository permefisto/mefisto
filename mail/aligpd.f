      SUBROUTINE ALIGPD( N0, N1, NOANCP, PXYD, EPSX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ALIGNER EN Y LES POINTS ENTRE NOANCP(N0) ET NOANCP(N1)
C -----    DU TABLEAU PXYD
C
C ENTREES:
C --------
C N0     : NUMERO DU PREMIER POINT
C N1     : NUMERO DU DERNIER POINT
C NOANCP : NUMERO DANS PXYD DE CHAQUE POINT
C EPSX   : PRECISION EN X POUR REPERER LE ZERO
C
C ENTREES ET SORTIES :
C --------------------
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JANVIER 1994
C....................................................................012
      DOUBLE PRECISION  PXYD(3,*), EPSX, D
      INTEGER           NOANCP(*)
C
      IF( N1 .GT. N0+1 ) THEN
C
C        AJUSTEMENT DE Y DES POINTS ALIGNES  (N0+1 A N1-1)
         NN0 = NOANCP(N0)
         NN1 = NOANCP(N1)
         D   = PXYD(1,NN1) - PXYD(1,NN0)
         IF( ABS(D) .GT. EPSX ) THEN
C           LE COEFFICIENT DIRECTEUR DE LA DROITE N0-N1
            D = (PXYD(2,NN1)-PXYD(2,NN0)) / D
            DO 10 K=N0+1,N1-1
C              L'ORDONNEE DU POINT K EST MODIFIEE
               NN = NOANCP(K)
               PXYD(2,NN) = PXYD(2,NN0) + D * (PXYD(1,NN)-PXYD(1,NN0))
 10         CONTINUE
         ENDIF
      ENDIF
      END
