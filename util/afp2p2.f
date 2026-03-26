      SUBROUTINE AFP2P2( NDEG, NBP, tPP )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  AFFICHER LE TABLEAU tPP(NBP,NBP)
C -----
C
C ENTREES:
C --------
C NDEG   : DEGRE DES POLYNOMES
C NBP    : NOMBRE DE POLYNOMES
C TPP    : INTEGRALES PNdeg PNdeg  dX   sur EF UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY    Avril 2011
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  TPP(NBP,NBP)
C
      PRINT *
      PRINT 10000,NDEG, NDEG
10000 FORMAT('INTEGRALE sur EF UNITE de P',I1,' P',I1,' dX')
      DO I = 1, NBP
         PRINT 10001,(I,J,tPP(I,J),J=1,NBP)
      ENDDO
10001 FORMAT(4('tPP(',I2,',',I2,')=',D25.17,'   '))
C
      RETURN
      END
