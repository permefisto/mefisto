      SUBROUTINE PN2PMDPNDPN( M, NBPM, PM,  N, NBPN, PN,  PMDPNDPN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES INTEGRALES PM DPN DPN dx dy SUR LE TRIANGLE UNITE
C -----
C
C ENTREES:
C --------
C M      : DEGRE  DU POLYNOME  PM
C NBPM   : NOMBRE DE POLYNOMES PM
C PM     : POLYNOME PM(0:M,0:M,NBPM)
C
C N      : DEGRE  DU POLYNOME  PN
C NBPM   : NOMBRE DE POLYNOMES PM
C PN     : POLYNOME PM(0:N,0:N,NBPN)
C
C SORTIE :
C --------
C PMDPNDPN : TABLEAU DES INTEGRALES PMi dPNj/dxk dPNjj/dxl  dx dy
C            SUR LE TRIANGLE RECTANGLE UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M University at QATAR  Fevrier 2012
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  PM( 0:M, 0:M, NBPM ), PN( 0:N, 0:N, NBPN ),
     %                  PMDPNDPN( NBPM, 2,NBPN, 2,NBPN )
      DOUBLE PRECISION  DP(200), P(1000), DQ(200), Q(1000)
C
      DO J = 1, NBPN
         DO K = 1, 2
C           DERIVATION DP = dPNj/dxk
            CALL PN2DDE( K, N+1, PN(0,0,J), DP )
C
            DO I = 1, NBPM
C              PRODUIT PMi dPNj/dxk
               CALL PN2DPR( 0,M+1,PM(0,0,I), 0,N+1,DP, M+N+1,P )
C
               DO JJ = 1, NBPN
                  DO L=1,2
C                    DERIVATION DQ = dPNjj/dxl
                     CALL PN2DDE( L, N+1, PN(0,0,JJ), DQ )
C
C                    PRODUIT PMi dPNj/dxk dPNjj/dxl
                     CALL PN2DPR( 0,M+N+1,P, 0,N+1,DQ,  M+N+1+N,Q )
C
C                    INTEGRALE DE PMi dPNj/dxk SUR LE TRIANGLE UNITE
                     CALL PN2DIN( 0,M+N+1+N,Q, PMDPNDPN(I,K,J,L,JJ) )
                  ENDDO
               ENDDO
C
            ENDDO
C
         ENDDO
      ENDDO
C
      PRINT *
      PRINT 10000,M,N,N
10000 FORMAT('INTEGRALE sur TRIANGLE UNITE de P',I1,' DP',I1,' DP',I1,
     %' dx dy')
C
      DO JJ = 1, NBPN
         DO L=1,2
            DO J = 1, NBPN
               DO K = 1, 2
                  DO I = 1, NBPM
                     IF( ABS(PMDPNDPN(I,K,J,L,JJ)) .LT. 1D-14 )
     %                       PMDPNDPN(I,K,J,L,JJ)=0D0
                  ENDDO
                  PRINT 10001,(I,K,J,L,JJ,PMDPNDPN(I,K,J,L,JJ),I=1,NBPM)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
10001 FORMAT(3('PNDPMDPM(',I2,',',I1,',',I2,',',I1,',',I2,')=',D25.17,
     %'   '))
C
      PRINT 10002,(((((PMDPNDPN(I,K,J,L,JJ),I=1,NBPM),K=1,2),J=1,NBPN),
     %                                                L=1,2),JJ=1,NBPN)
10002 FORMAT('     %',D25.17,',',D25.17,',')
C
      RETURN
      END
