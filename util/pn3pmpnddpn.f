      SUBROUTINE PN3PMPNDDPN( M, NBPM, PM,  N, NBPN, PN, PMPNDDPN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES INTEGRALES PM PN DDPN dx dy SUR LE TRIANGLE UNITE
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
C PMPNDDPN : TABLEAU DES INTEGRALES PMi PNj ddPNjj/dxkdxl  dx dy
C            SUR LE TRIANGLE RECTANGLE UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M University at QATAR  Fevrier 2012
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  PM( 0:M, 0:M, 0:M, NBPM ),
     %                  PN( 0:N, 0:N, 0:N, NBPN ),
     %                  PMPNDDPN( NBPM, NBPN, 3,3, NBPN )
      DOUBLE PRECISION  DP(300), P(1000), DDP(300), Q(1000)
C
      DO J = 1, NBPN
C
         DO I = 1, NBPM
C           PRODUIT PMi PNj
            CALL PN3DPR( 0,M+1,PM(0,0,0,I), 0,N+1,PN(0,0,0,J), M+N+1,P )
C
            DO JJ = 1, NBPN
               DO L=1,3
C                 DERIVATION DP = dPNjj/dxl
                  CALL PN3DDE( L, N+1, PN(0,0,0,JJ), DP )
C
                  DO K=1,3
C                    DERIVATION DDP = ddPNjj/dxkdxl
                     CALL PN3DDE( K, N+1, DP, DDP )
C
C                    PRODUIT PMi PNj ddPNjj/dxldxk
                     CALL PN3DPR( 0,M+N+1,P, 0,N+1,DDP,  M+N+1+N,Q )
C
C                    INTEGRALE DE PMi PNj ddPNjj/dxldxk SUR LE TRIANGLE UNITE
                     CALL PN3DIN( 0,M+N+1+N,Q, PMPNDDPN(I,J,K,L,JJ) )
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
         DO L=1,3
            DO K = 1, 3
               DO J = 1, NBPN
                  DO I = 1, NBPM
                     IF( ABS(PMPNDDPN(I,J,K,L,JJ)) .LT. 1D-14 )
     %                       PMPNDDPN(I,J,K,L,JJ)=0D0
                  ENDDO
                  PRINT 10001,(I,J,K,L,JJ,PMPNDDPN(I,J,K,L,JJ),I=1,NBPM)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
10001 FORMAT(3('PMPNDDPN(',I2,',',I2,',',I1,',',I1,',',I2,')=',D25.17,
     %'   '))
C
      PRINT 10002,(((((PMPNDDPN(I,J,K,L,JJ),I=1,NBPM),J=1,NBPN),
     %                                      K=1,3),L=1,3),JJ=1,NBPN)
10002 FORMAT('     %',D25.17,',',D25.17,',')
C
      RETURN
      END
