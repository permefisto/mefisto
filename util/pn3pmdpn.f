      SUBROUTINE PN3PMDPN( M, NBPM, PM, N, NBPN, PN, PMDPN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   RECHERCHE DES INTEGRALES PM DPN dx dy dz SUR LE TETRAEDRE UNITE
C -----
C
C ENTREES:
C --------
C M      : DEGRE  DU POLYNOME  PM
C NBPM   : NOMBRE DE POLYNOMES PM
C PM     : POLYNOME PM(0:M,0:M,NBPM)
C
C N      : DEGRE  DU POLYNOME  PN
C NBPN   : NOMBRE DE POLYNOMES PN
C PN     : POLYNOME PM(0:N,0:N,NBPN)
C
C SORTIE :
C --------
C PMDPN : TABLEAU DES INTEGRALES PMi dPNj/dxk dx dy SUR LE TETRAEDRE UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY    Avril 2011
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  PM( 0:M, 0:M, 0:M, NBPM ),
     %                  PN( 0:N, 0:N, 0:N, NBPN ),
     %                  PMDPN( NBPM, 3, NBPN )
      DOUBLE PRECISION  DP(500), P(1000)
C
      DO J = 1, NBPN
         DO K = 1, 3
C
C           DERIVATION DP = dPNj/dxk
            CALL PN3DDE( K, N+1, PN(0,0,0,J), DP )
C
            DO I = 1, NBPM
C
C              PRODUIT PMi dPNj/dxk
               CALL PN3DPR( 0, M+1, PM(0,0,0,I),  0, N+1, DP,  M+N+1,P )
C
C              INTEGRALE DE PMi dPNj/dxk SUR LE TETRAEDRE UNITE
               CALL PN3DIN( 0, M+N+1, P, PMDPN(I,K,J) )
C
            ENDDO
C
         ENDDO
      ENDDO
C
      PRINT *
      PRINT 10000,M,N
10000 FORMAT('INTEGRALE sur TETRAEDRE UNITE P',I1,' DP',I1,' dx dy dz')
      DO J = 1, NBPN
         DO K = 1, 3
            DO I = 1, NBPM
               IF( ABS(PMDPN(I,K,J)) .LT. 1D-14 ) PMDPN(I,K,J)=0D0
            ENDDO
            PRINT 10001,(I,K,J,PMDPN(I,K,J),I=1,NBPM)
         ENDDO
      ENDDO
10001 FORMAT(3('PNDPM(',I2,',',I1,',',I2,')=',D25.17,'   '))
C
      PRINT 10002,(((PMDPN(I,K,J),I=1,NBPM),K=1,3),J=1,NBPN)
10002 FORMAT('     %',D25.17,',',D25.17,',')
C
      RETURN
      END
