      SUBROUTINE BSPLS1( DEGREX, LUX, LTX, LRX, UX,
     %                   TX, RX, NORXTX )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCULER LES NOEUDS TX ET RX  A PARTIR DES NOEUDS UX
C -----  ET LE TABLEAU POINTEUR SUR LE DERNIER TX DE CHAQUE RX
C        MORALEMENT LE NOEUD UX(J) EST AU MILIEU DU SUPPORT
C                   DE B(J,DEGREX)
C
C ENTREES:
C --------
C DEGREX : DEGRE  EN X DES POLYNOMES DE LA B-SPLINE
C LUX    : NOMBRE-1 DE NOEUDS D'INTERPOLATION EN X
C LTX    : NOMBRE-1 DE NOEUDS EN X POUR LES BX(J,M)
C LRX    : NOMBRE-1 DE NOEUDS TX IDENTIFIES EN X
C UX     : LES VALEURS DES NOEUDS D'INTERPOLATION EN X
C
C SORTIES:
C --------
C TX     : LES VALEURS DES NOEUDS DE BX(J,M) EN X
C RX     : LES VALEURS DES NOEUDS IDENTIFIES DE TX
C NORXTX : LE DERNIER TX DE CHAQUE RX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       JUIN  1990
C2345X7..............................................................012
      INTEGER  DEGREX, NORXTX(0:LRX)
      REAL     UX(0:LUX),
     %         RX(0:LRX),
     %         TX(0:LTX+DEGREX)

C     LE PREMIER NOEUD EN X DE MULTIPLICITE DEGREX
C     => PASSAGE PAR LE POINT
      DO J=0,DEGREX
         TX(J) = UX(0)
      ENDDO
      RX(0)     = UX(0)
      NORXTX(0) = DEGREX
      NOI       = DEGREX
C
C     LES NOEUDS INTERMEDIAIRES DE MULTIPLICITE 1
C     => CLASSE DE CONTINUITE = DEGREX-1
      DO I = 1, LRX - 1
C        CONTINUITE MAXIMALE DEGREX-1
         SS = 0
         DO J = I, I + DEGREX - 1
            SS = SS + UX(J)
         ENDDO
         NOI       = NORXTX(I-1) + 1
         NORXTX(I) = NOI
         TX(NOI)   = SS / DEGREX
         RX(I)     = TX(NOI)
      ENDDO
C
C     LE DERNIER NOEUD DE MULTIPLICITE DEGREX
C     => PASSAGE PAR LE POINT
      DO J=1,DEGREX+1
         TX(NOI+J) = UX(LUX-1)
      ENDDO
      NORXTX(LRX) = NOI + DEGREX + 1
      RX(LRX)     = UX(LUX-1)

      RETURN
      END
