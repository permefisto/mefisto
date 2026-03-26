      SUBROUTINE PN2DPN( N, NBPN, PN, IntDPN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECHERCHE DES INTEGRALES DPN dx dy SUR LE TRIANGLE UNITE
C -----
C ENTREES:
C --------
C N      : DEGRE  DU POLYNOME  PN
C NBPN   : NOMBRE DE POLYNOMES PN
C PN     : POLYNOME PN(0:N,0:N, NBPN)
C
C SORTIE :
C --------
C IntDPN : TABLEAU DES INTEGRALES dPNj/dxk dx dy SUR LE TRIANGLE UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Decembre 2012
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  PN( 0:N, 0:N, NBPN ), IntDPN( 2, NBPN )
      DOUBLE PRECISION  DP(125)
C
      DO J = 1, NBPN
         DO K = 1, 2
C
C           DERIVATION DP = dPNj/dxk
            CALL PN2DDE( K, N+1, PN(0,0,J), DP )
C
C           INTEGRALE DE PNi dPNj/dxk SUR LE TRIANGLE UNITE
            CALL PN2DIN( 0, N+1, DP, IntDPN(K,J) )
C
         ENDDO
C
      ENDDO
C
      PRINT *
      PRINT 10000,N
10000 FORMAT('PN2DPN: INTEGRALE sur TRIANGLE UNITE DP',I1,' dx dy')
      DO J = 1, NBPN
         DO K = 1, 2
            IF( ABS(INTDPN(K,J)) .LT. 1D-14 ) IntDPN(K,J)=0D0
            PRINT 10001,K,J,IntDPN(K,J)
         ENDDO
      ENDDO
10001 FORMAT(2('PN2DPN(',I1,',',I2,')=',D25.17,'   '))
C
      PRINT 10002,((IntDPN(K,J),K=1,2),J=1,NBPN)
10002 FORMAT('     %',D25.17,',',D25.17,',')
C
      RETURN
      END
