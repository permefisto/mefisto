      SUBROUTINE ST81C6C( I, J, K, L, M, N,
     %                    IND1, VIND1, IND2, VIND2, IND3, VIND3,
     %                    NU3C, NS3C6C )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LE NO DES 8 SOMMETS DU 3-CUBE NU3C TEL QUE
C L'INDICE IND1 (I ou J ou K ou L ou M ou N) A LA VALEUR VIND1 (0 ou 1)
C L'INDICE IND2 (I ou J ou K ou L ou M ou N) A LA VALEUR VIND2 (0 ou 1)
C L'INDICE IND3 (I ou J ou K ou L ou M ou N) A LA VALEUR VIND3 (0 ou 1)
C
C SORTIES:
C --------
C NU3C   : NUMERO DU 3-CUBE (de 1 a 160) DANS LE 6-CUBE
C NS3C6C : NS3C6C(P,NU3C) NO DU P-EME SOMMET DU 3-CUBE NU3C D'UN 6-CUBE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: : ALAIN PERRONNET  LJLL UPMC PARIS                OCTOBRE 2005
C2345X7..............................................................012
      INTEGER  I, J, K, L, M, N,  IND1, VIND1, IND2, VIND2, IND3, VIND3
      INTEGER  NS3C6C(8,160)
C
C     3-CUBE <=> IND1=VIND1 IND2=VIND2 IND3=VIND3 A FIXER
      NS6C = 0
      NS3C = 0
      DO 16 N=0,1
         DO 15 M=0,1
            DO 14 L=0,1
               DO 13 K=0,1
                  DO 12 J=0,1
                     DO 11 I=0,1
C                       LE NO DE SOMMET DU 6-CUBE
                        NS6C = NS6C + 1
C
                        IF( IND1.EQ.VIND1  .AND.
     %                      IND2.EQ.VIND2  .AND.
     %                      IND3.EQ.VIND3 ) THEN
C                           LE SOMMET EST SUR LE 3-CUBE
                            NS3C = NS3C + 1
                            NS3C6C( NS3C, NU3C ) = NS6C
                        ENDIF
C
 11                  CONTINUE
 12               CONTINUE
 13            CONTINUE
 14         CONTINUE
 15      CONTINUE
 16   CONTINUE
C
      RETURN
      END
