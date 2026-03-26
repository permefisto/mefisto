      SUBROUTINE PN3DPMDPM( M, NBPM, PM, tDPMDPM )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECHERCHE DES INTEGRALES DPM DPM dx dy dz SUR LE TETRAEDRE UNITE
C -----
C
C ENTREES:
C --------
C M      : DEGRE  DU POLYNOME  PM
C NBPM   : NOMBRE DE POLYNOMES PM
C PM     : POLYNOME PM(0:M,0:M,NBPM)
C
C SORTIE :
C --------
C TDPMDPM : TABLEAU DES INTEGRALES dPMi/dxk dPMj/dxl dx dy dz
C           SUR LE TETRAEDRE UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY    Avril 2011
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  PM( 0:M, 0:M, 0:M, NBPM ),
     %                  tDPMDPM( 3, NBPM, 3, NBPM )
      DOUBLE PRECISION  DPi(100), DPj(100), P(1000)
C
      DO J = 1, NBPM
         DO L = 1, 3
C
C           DERIVATION DP = dPMj/dxl
            CALL PN3DDE( L, M+1, PM(0,0,0,J), DPj )
C
            DO I = 1, NBPM
               DO K = 1, 3
C
C                 DERIVATION DP = dPMi/dxk
                  CALL PN3DDE( K, M+1, PM(0,0,0,I), DPi )
C
C                 PRODUIT dPMi/dxk dPMj/dxl
                  CALL PN3DPR( 0, M+1, DPi, 0, M+1, DPJ, M+M+1, P )
C
C                 INTEGRALE DE PMi dPMj/dxk SUR LE TETRAEDRE UNITE
                  CALL PN3DIN( 0, M+M+1, P, tDPMDPM(K,I,L,J) )
C
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
      PRINT *
      PRINT 10000,M,M
10000 FORMAT('INTEGRALE sur TETRAEDRE UNITE DP',I1,' DP',I1,
     %       ' dx dy dz')
      DO J = 1, NBPM
         DO L = 1, 3
            DO I = 1, NBPM
               DO K=1,3
                 IF(ABS(tDPMDPM(K,I,L,J)).LT.1D-14) tDPMDPM(K,I,L,J)=0D0
               ENDDO
               PRINT 10001,(K,I,L,J,tDPMDPM(K,I,L,J),K=1,3)
            ENDDO
         ENDDO
      ENDDO
10001 FORMAT(3('DPMDPM(',I1,',',I2,',',I1,',',I2,')=',D25.17,'   '))
C
      PRINT 10002,((((tDPMDPM(K,I,L,J),K=1,3),I=1,NBPM),L=1,3),J=1,NBPM)
10002 FORMAT('     %',D25.17,',',D25.17,',')
C
      RETURN
      END
