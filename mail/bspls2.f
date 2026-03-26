       SUBROUTINE BSPLS2( UXI , DEGREX , LTX , LRX ,
     %                    TX , RX , NORXTX , BX , J )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL des BX ( -DEGREX:0 , DEGREX ) (UXI)
C -----  POUR UNE LIGNE B-SPLINE
C
C ENTREES:
C --------
C UXI    : VALEUR OU LES BX(J,M) SONT CALCULEES
C DEGREX : DEGRE  EN X DES POLYNOMES DE LA B-SPLINE
C LTX    : NOMBRE-1 DE NOEUDS EN X POUR LES BX(J,M)
C LRX    : NOMBRE-1 DE NOEUDS TX IDENTIFIES EN X
C TX     : LES VALEURS DES NOEUDS DE BX(J,M) EN X
C RX     : LES VALEURS DES NOEUDS IDENTIFIES DE TX
C NORXTX : LE DERNIER TX DE CHAQUE RX
C
C SORTIES:
C --------
C BX     : BX(-DEGREX:1)(UXI)
C J      : L'INTERVALLE [RX(J),RX(J+1)[ CONTIENT UXI
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       JUIN  1990
C2345X7..............................................................012
      INTEGER  NORXTX(0:LRX),DEGREX
      REAL     RX(0:LRX),
     %         TX(0:LTX+DEGREX),
     %         BX(-DEGREX:1)
C
C     QUEL EST L'INTERVALLE [RX(J),RX(J+1)[ CONTENANT UXI ?
      DO 10 J=0,LRX-1
         IF( UXI .LT. RX(J+1) ) GOTO 30
 10   CONTINUE
      J = LRX
C
C     LE NUMERO DU DERNIER NOEUD TX AU POINT RX(J)
 30   NOI = NORXTX( J )
C
C     CALCUL DES BX(J,M()(UXI)
C     ========================
      DO 40 JX=-DEGREX,1
         BX(JX) = 0
 40   CONTINUE
      BX(0) = 1.0
C
      DO 60 M=1,DEGREX
         DO 50 JX=-M,0
C           LE VRAI INDICE DE B
            JB = NOI + JX
C           EVALUXATION DES FRACTIONS DE LA FORMULE
            IF( TX(JB+M) .NE. TX(JB) ) THEN
               U1 = ( UXI - TX(JB) ) / ( TX(JB+M) - TX(JB) )
            ELSE
               U1 = 0
            ENDIF
            JBM1 = JB + M + 1
            IF( TX(JBM1) .NE. TX(JB+1) ) THEN
               U2 = ( TX(JBM1) - UXI )
     %            / ( TX(JBM1) - TX(JB+1) )
            ELSE
               U2 = 0
            ENDIF
            BX(JX) = U1 * BX(JX) + U2 * BX(JX+1)
 50      CONTINUE
 60   CONTINUE
      END
