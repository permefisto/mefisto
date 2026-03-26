      SUBROUTINE NIVOEF6CUB( NA6CUB, NOEF, NIVO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : NO DES 6 INDICES DIFFERENCES FINIES DU 6-CUBE DE NUMERO NOEF
C -----
C    CALCUL DU NO (DE 0 a NA6CUB-1) DE L'EF SELON LE NIVEAU DES COORDONNEES
C    NOEF - 1 = NIVO(6) * NX**5 + NIVO(5) * NX**4 + ... + NIVO(1) * NX**0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC LJLL                        OCTOBRE 2005
C2345X7..............................................................012
      INTEGER NIVO(6)
C
      NA   = NA6CUB ** 6
      N    = NOEF - 1
      NOSO = 0
      DO 10 I=6,1,-1
         NA      = NA / NA6CUB
         NIVO(I) = N  / NA
         N       = N  - NIVO(I) * NA
 10   CONTINUE
      RETURN
      END
