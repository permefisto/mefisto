      SUBROUTINE F3EX3P1BP1( NONOEF, DELTA, Rho, VITRA, BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TETRAEDRE BREZZI FORTIN
C -----    POUR LA VITESSE TRANSPORTEE
C
C ENTREES:
C --------
C NONOEF : NUMERO GLOBAL DES 5 NOEUDS DU TETRAEDRE
C DELTA  : JACOBIEN DU TETRAEDRE
C Rho    : DENSITE VOLUMIQUE DE MASSE
C VITRA  : DL DES 3 COMPOSANTES DE LA VITESSE TRANSPORTEE
C
C SORTIES:
C --------
C BE     : VECTEUR ELEMENTAIRE DE LA VITESSE TRANSPORTEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Decembre 2010
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  DELTA, Rho, VITRA(3,*), BE(5,3), S
      INTEGER           NONOEF(5)
C
C     COEFFICIENT DES SOMMETS AVEC LE POIDS D'INTEGRATION NUMERIQUE
      S = Rho * DELTA / 120.D0
      DO K=1,3
         DO I=1,5
            BE(I,K) = S * VITRA( K, NONOEF(I) )
         ENDDO
C        CORRECTION POUR LE POIDS D'INTEGRATION NUMERIQUE
         BE(5,K) = BE(5,K) * 16.D0
      ENDDO
C
      RETURN
      END
