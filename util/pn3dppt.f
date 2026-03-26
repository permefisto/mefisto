      SUBROUTINE PN3DPPT( M, NBPM, PM, NBPOINT, XYPOINT, DPS )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA DERIVEE DES POLYNOMES EN DES POINTS
C -----
C
C ENTREES:
C --------
C M      : DEGRE  DU POLYNOME  PM
C NBPM   : NOMBRE DE POLYNOMES PM
C PM     : POLYNOME PM(0:M,0:M,NBPM)
C NBPOINT: NOMBRE DE POINTS DE R3
C XYPOINT: X,Y DES NBPOINT POINTS
C
C SORTIE :
C --------
C DPS    : DPj/Dxl ( XYPOINT(K) ) VALEUR DE LA DERIVEE PREMIERE
C          DANS LA DIRECTION l DU POLYNOME j AU POINT k
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR  Fevrier 2012
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  PM( 0:M, 0:M, 0:M, NBPM ), XYPOINT(3, NBPOINT),
     %                  DPS( 3, NBPM, NBPOINT )
      DOUBLE PRECISION  DPj( 100 )
C
C     CALCUL DES DPS(L,J,K)=dPj/dxl(Pt k)
      DO J = 1, NBPM
C
         DO L = 1, 3
C
C           DERIVATION DP = dPMj/dxl
            CALL PN3DDE( L, M+1, PM(0,0,0,J), DPj )
C
            DO K = 1, NBPOINT
C
C              DPj/Dxl ( XYPOINT(K) )
               CALL PN3DVA( M+1, DPj,
     %                      XYPOINT(1,K), XYPOINT(2,K), XYPOINT(3,K),
     %                      DPS(L,J,K) )
C
            ENDDO
C
         ENDDO
C
      ENDDO
C
C     AFFICHAGE DES DPS(L,J,K)=dPj/dxl(Pt k)
      PRINT *
      PRINT 10000,M,NBPOINT
10000 FORMAT('VALEUR de DP',I1,' aux ',I2,' POINTS de R3' )
      DO K=1,NBPOINT
         DO J = 1, NBPM
            DO L = 1, 3
               IF(ABS(DPS(L,J,K)).LT.1D-14) DPS(L,J,K)=0D0
            ENDDO
            PRINT 10001,(L,J,K,DPS(L,J,K),L=1,3)
         ENDDO
      ENDDO
10001 FORMAT(3('DPS(',I1,',',I2,',',I1,')=',D25.17,'   '))
C
      PRINT 10002,((( DPS(L,J,K),L=1,3),J=1,NBPM),K=1,NBPOINT)
10002 FORMAT('     %',D25.17,',',D25.17,',' )
C
      RETURN
      END
