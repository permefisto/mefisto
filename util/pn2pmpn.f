      SUBROUTINE PN2PMPN( M, NBPM, PM, N, NBPN, PN, tPMPN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECHERCHE DES INTEGRALES PM PN dx dy SUR LE TRIANGLE UNITE
C -----
C
C ENTREES:
C --------
C M      : DEGRE  DU POLYNOME  PM
C NBPM   : NOMBRE DE POLYNOMES PM
C PM     : POLYNOME PM(0:M,0:M,NBPM)
C
C N      : DEGRE DU POLYNOME PN
C NBPN   : NOMBRE DE POLYNOMES PN
C PN     : POLYNOME PM(0:N,0:N,NBPN)
C
C SORTIE :
C --------
C tPMPN  : TABLEAU DES INTEGRALES PMi PNj dx dy SUR LE TRIANGLE UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY    Avril 2011
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  PM( 0:M, 0:M, NBPM ), PN( 0:N, 0:N, NBPN ),
     %                  tPMPN( NBPM, NBPN )
      DOUBLE PRECISION  P(100)
C
      DO I = 1, NBPM
         DO J = 1, NBPN
C
C           PRODUIT PMi PNJ
            CALL PN2DPR( 0, M+1, PM(0,0,I),  0, N+1, PN(0,0,J),
     %                   M+N+1, P )
C
C           INTEGRALE DE PMi PNj SUR LE TRIANGLE UNITE
            CALL PN2DIN( 0, M+N+1, P, tPMPN(I,J) )
C
         ENDDO
      ENDDO
C
      PRINT *
      PRINT 10000,M,N
10000 FORMAT('INTEGRALE sur TRIANGLE UNITE de P',I1,' P',I1,' dx dy')
      DO J = 1, NBPN
         DO I = 1, NBPM
            IF( ABS(tPMPN(I,J)) .LT. 1D-14 ) tPMPN(I,J) = 0D0
         ENDDO
         PRINT 10001,(I,J,tPMPN(I,J),I=1,NBPM)
      ENDDO
10001 FORMAT(3('PNDPM(',I2,',',I2,')=',D25.17,'   '))
C
      PRINT 10002,((tPMPN(I,J),I=1,NBPM),J=1,NBPN)
10002 FORMAT('     %',D25.17,',',D25.17,',')
C
      RETURN
      END
