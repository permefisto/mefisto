      SUBROUTINE NLSEKeMe( NBLIGN, CoefKe, Ke, CoefMe, Me, KeMe )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    NLSE: CONSTRUIRE LA MATRICE KeMe A PARTIR DES MATRICES
C -----          Ke et Me EN LES MULTIPLIANT PAR LEURS COEFFICIENTS
C
C             [-CoefKe Ke, CoefMe Me ]
C    [KeMe] =
C             [ CoefMe Me, CoefKe Ke ]
C
C ENTREES:
C --------
C NBLIGN : NOMBRE DE LIGNES et COLONNES DES MATRICES Ke Me
C CoefKe : Coefficient Multiplicatif de la MATRICE Ke
C Ke     : Matrice symetrique pleine rangee par lignes
C          KeMe(2*NBLIGN*(2*NBLIGN+1)/2)
C CoefMe : Coefficient Multiplicatif de la MATRICE Me
C Me     : Matrice symetrique pleine rangee par lignes
C          Me(NBLIGN*(NBLIGN+1)/2)
C
C SORTIE :
C --------
C KeMe   : Matrice VUE CI DESSUS stockee SYMETRIQUE PLEINE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      DOUBLE PRECISION  CoefKe, Ke(1:*), CoefMe, Me(1:*), KeMe(1:*), S
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
C           POSITION DANS LE BLOC 1 1 DU COEFFICIENT de KeMe
            K11 = K11 + 1
C
C           LE BLOC 1 1: -CoefKe Ke
            S = CoefKe * Ke(K11)
            KeMe( K11 ) = - S
C           LE BLOC 2 2: CoefKe Ke
            K22 = K22 + 1
            KeMe( K22 ) = S
C
C           LE BLOC 2 1 Partie inferieure:  CoefMe Me
            S = CoefMe * Me( K11 )
            KeMe( K21+J ) = S
C           LE BLOC 2 1 Partie superieure:  CoefMe Me
            KeMe( (NBLIGN+J)*(NBLIGN+J-1)/2+I ) = S
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
