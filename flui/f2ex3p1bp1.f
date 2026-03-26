      SUBROUTINE F2EX3P1BP1( NONOEF, DELTA, Rho, VITRA, BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TRIANGLE BREZZI FORTIN
C -----    POUR LA VITESSE TRANSPORTEE
C
C ENTREES:
C --------
C NONOEF : NUMERO GLOBAL DES 4 NOEUDS DU TRIANGLE (3 SOMMETS + BARYCENTRE)
C DELTA  : JACOBIEN DU TRIANGLE
C Rho    : DENSITE VOLUMIQUE DE MASSE
C VITRA  : DL DES 2 COMPOSANTES DE LA VITESSE TRANSPORTEE
C
C SORTIES:
C --------
C BE     : VECTEUR ELEMENTAIRE DE LA VITESSE TRANSPORTEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Decembre 2010
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  DELTA, Rho, VITRA(2,*), BE(4,2), S
      INTEGER           NONOEF(4)
C
C     COEFFICIENT DES SOMMETS AVEC LE POIDS D'INTEGRATION NUMERIQUE
      S = Rho * DELTA / 24.D0
      DO K=1,2
         DO I=1,4
            BE(I,K) = S * VITRA( K, NONOEF(I) )
         ENDDO
C        CORRECTION POUR LE POIDS D'INTEGRATION NUMERIQUE
         BE(4,K) = BE(4,K) * 9.D0
      ENDDO
C
      RETURN
      END
