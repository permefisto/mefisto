      SUBROUTINE DDP2TR( P2, DDP2 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES DERIVEES SECONDES DU POLYNOME P2 DE DEGRE 2
C -----    VALABLE SEULEMENT POUR LES POLYNOMES DE DEGRE 2 SUR LE TRIANGLE
C ENTREES:
C --------
C P2     : POLYNOME P2(0:2,0:2,6)
C
C SORTIE :
C --------
C DDP2   : TABLEAU DDP2(i,j,k)=ddP2k/dxi dxj
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M University at QATAR     Mars 2012
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  P2( 0:2, 0:2, 6 ), DDP2(2,2,6)
      DOUBLE PRECISION  DP(48), DDP(48)
C
      DO K=1,6
C
         print *
         DO J = 1,2
C
C           DERIVATION DP = dP2k/dxj
            CALL PN2DDE( J, 2+1, P2(0,0,K), DP )
C
            DO I = 1,2
C
C              DERIVATION DP = dP2k/dxj => d2P2k/dxi dxj CONSTANTE
               CALL PN2DDE( I, 2+1, DP, DDP )
C
C              CONSTANTE
               DDP2(I,J,K) = DDP(1)
               PRINT 10001,I,J,K,DDP2(I,J,K)
C
            ENDDO
C
         ENDDO
      ENDDO
C
      PRINT *
      PRINT 10000
10000 FORMAT('DDP2(i,j,k)=DDP2k /dxi dxj SUR LE TRIANGLE UNITE' )
C
      DO K=1,6
C
         DO J = 1,2
            DO I = 1, 2
               IF( ABS(DDP2(I,J,K)) .LT. 1D-14 ) DDP2(I,J,K)=0D0
            ENDDO
            PRINT 10001,(I,J,K,DDP2(I,J,K),I=1,2)
         ENDDO
C
      ENDDO
C
10001 FORMAT( 2('DDP2(',I1,',',I1,',',I1,')=',D25.17,'   ') )
C
      PRINT 10002,(((DDP2(I,J,K),I=1,2),J=1,2),K=1,6)
10002 FORMAT('     %',D25.17,',',D25.17,',')
C
      RETURN
      END
