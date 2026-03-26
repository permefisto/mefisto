      INTEGER FUNCTION N1ST6C( I, J, K, L, M, N )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : NO ELEMENTAIRE DU SOMMET(I,J,K,L,M,N) DANS LE 6-CUBE
C -----
C CALCUL DU NO DE L'EF NEF SELON LES NIVO
C NEF - 1 = NIVO(6) * NA**5 + NIVO(5) * NA**4 + ... + NIVO(1) * NA**0
C           NIVO(.) de 0 a NA-1
C NOSO = NO DU PREMIER SOMMET DU 6-CUBE NUELEM  NS=NA+1
C NOSO -1 = NIVO(6) * NS**5 + NIVO(5) * NS**4 + ... + NIVO(1)* NS**0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC LJLL                        OCTOBRE 2005
C2345X7..............................................................012
      N1ST6C = ((((N * 2 + M) * 2 + L) * 2 + K) * 2 + J) * 2 + I + 1
      RETURN
      END
