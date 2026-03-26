      SUBROUTINE PN2PNPNDPN( N, NBPN, PN,  PNPNDPN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL DES INTEGRALES PN DPN DPN dx dy SUR LE TRIANGLE UNITE
C -----
C
C ENTREES:
C --------
C N      : DEGRE  DU POLYNOME  PN
C NBPN   : NOMBRE DE POLYNOMES PN
C PN     : POLYNOME PN( 0:N, 0:N, NBPN )
C
C SORTIE :
C --------
C PNPNDPN: TABLEAU DES INTEGRALES PNi PNj dPNl/dxk dx dy
C          SUR LE TRIANGLE RECTANGLE UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M University at QATAR     Mars 2012
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  PN( 0:N, 0:N, NBPN ),
     %                  PNPNDPN( NBPN, NBPN, 2, NBPN )
      DOUBLE PRECISION  DP(100), P(500), Q(500)
C
      DO L = 1, NBPN
C
         DO K = 1, 2
C
C           DERIVATION DP = dPNl/dxk
            CALL PN2DDE( K, N+1, PN(0,0,L), DP )
C
            DO J = 1, NBPN
C
C              PRODUIT PNj dPNl/dxk
               CALL PN2DPR( 0,N+1,PN(0,0,J), 0,N+1,DP, N+N+1,P )
C
               DO I = 1, NBPN
C
C                 PRODUIT PNi PNj dPNl/dxk
                  CALL PN2DPR( 0,N+1,PN(0,0,I), 0,N+N+1,P, N+N+N+1,Q )
C
C                 INTEGRALE DE PNi dPNj/dxk SUR LE TRIANGLE UNITE
                  CALL PN2DIN( 0,N+N+N+1,Q, PNPNDPN(I,J,K,L) )
C
               ENDDO
C
            ENDDO
C
         ENDDO
C
      ENDDO
C
      PRINT *
      PRINT 10000,N,N,N
10000 FORMAT('INTEGRALE sur le TRIANGLE UNITE de P',I1,' P',I1,' DP',
     %I1,' dx dy')
C
      DO L = 1, NBPN
         DO K=1,2
            DO J = 1, NBPN
               DO I = 1, NBPN
                  IF( ABS(PNPNDPN(I,J,K,L)) .LT. 1D-14 )
     %                    PNPNDPN(I,J,K,L)=0D0
               ENDDO
               PRINT 10001,(I,J,K,L,PNPNDPN(I,J,K,L),I=1,NBPN)
            ENDDO
         ENDDO
      ENDDO
10001 FORMAT(4('PNPNDPN(',I2,',',I2,',',I1,',',I2,')=',D25.17,'   '))
C
      PRINT 10002,((((PNPNDPN(I,J,K,L),I=1,NBPN),J=1,NBPN),K=1,2),
     %               L=1,NBPN)
10002 FORMAT('     %',D25.17,',',D25.17,',')
C
      RETURN
      END
