      SUBROUTINE NLSEKeMeRe( NBLIGN, CoefKe, Ke, CoefMe, Me, CoefRe, Re,
     %                       KeMeRe )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : NLSE: CONSTRUIRE LA MATRICE KeMe A PARTIR DES MATRICES Ke Me Re
C -----       Re est LA MATRICE DE ROTATION

C               [-CoefKe Ke,             CoefMe Me + CoefRe Re]
C    [KeMere] =
C               [ CoefMe Me + CoefRe Re, CoefKe Ke            ]

C ENTREES:
C --------
C NBLIGN : NOMBRE DE LIGNES et COLONNES DES MATRICES Ke Me
C CoefKe : Coefficient Multiplicatif de la MATRICE Ke
C Ke     : Matrice symetrique pleine rangee par lignes
C          KeMe(2*NBLIGN*(2*NBLIGN+1)/2)

C CoefMe : Coefficient Multiplicatif de la MATRICE Me
C Me     : Matrice symetrique pleine rangee par lignes
C          Me(NBLIGN*(NBLIGN+1)/2)

C CoefRe : Coefficient Multiplicatif de la MATRICE 
C Re     : Matrice pleine rangee par lignes
C          Re(NBLIGN,NBLIGN)
C
C SORTIE :
C --------
C KeMeRe : Matrice(NBLIGN*NBLIGN) VUE CI DESSUS stockee SYMETRIQUE PLEINE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Novembre 2013
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      DOUBLE PRECISION  CoefKe, Ke(1:*),
     %                  CoefMe, Me(1:*),
     %                  CoefRe, Re(1:NBLIGN,1:NBLIGN),
     %                  KeMeRe(1:*), S
      INTEGER           NBLIGN, FINB11, K21, K22, K11, I, J
C
C     NUMERO DU DERNIER COEFFICIENT DU BLOC 1 1
      FINB11 = NBLIGN * ( NBLIGN + 1 ) / 2
      K21    = FINB11
      K22    = FINB11
      K11    = 0
C
      DO I = 1, NBLIGN
C
C        LE COEFFICIENT AVANT LE PREMIER DE LA LIGNE I DU BLOC 2x2
         K22 = K22 + NBLIGN
C
         DO J = 1, I
C
C           POSITION DANS LE BLOC 1 1 DU COEFFICIENT de KeMeRe
            K11 = K11 + 1
C
C           LE BLOC 1 1: -CoefKe Ke
            S = CoefKe * Ke(K11)
            KeMeRe( K11 ) = - S
C           LE BLOC 2 2: CoefKe Ke
            K22 = K22 + 1
            KeMeRe( K22 ) = S
C
C           LE BLOC 2 1 Partie inferieure:  CoefMe Me+ CoefRe * Re(I,J)
            S = CoefMe * Me( K11 ) + CoefRe * Re(I,J)
            KeMeRe( K21+J ) = S
C
C           LE BLOC 2 1 Partie superieure:  CoefMe Me+ CoefRe * Re(J,I)
            S = CoefMe * Me( K11 ) + CoefRe * Re(J,I)
            KeMeRe( (NBLIGN+J)*(NBLIGN+J-1)/2+I ) = S
C
         ENDDO
C
C        LA LIGNE SUIVANTE DU BLOC 2 1
         K21 = K21 + NBLIGN + I
C
      ENDDO
C
      RETURN
      END
