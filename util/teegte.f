      SUBROUTINE TEEGTE( NOSOTE1, NOSOTE2, NONOUI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: LES 4 SOMMETS DE NOSOTE1 SONT ILS LES 4 SOMMETS DE NOSOTE2?
C ----

C ENTREES:
C --------
C NOSOTE1 : 4 NUMEROS DE SOMMETS DU TETRAEDRE 1
C NOSOTE2 : 4 NUMEROS DE SOMMETS DU TETRAEDRE 2

C SORTIE :
C --------
C NONOUI : =1 SI LES 4 SOMMETS DE NOSOTE1 SONT CEUX DE NOSOTE2
C          =0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray          Novembre 2018
C2345X7..............................................................012
      INTEGER  NOSOTE1(4), NOSOTE2(4)

      DO 10 I1 = 1, 4
         N1 = NOSOTE1( I1 )
         DO I2 = 1, 4
            N2 = NOSOTE2( I2 )
            IF( N1 .EQ. N2 ) GOTO 10
         ENDDO
         GOTO 9000
 10   ENDDO

C     LES 4 SOMMETS DE NOSOTE1 SONT CEUX DE NOSOTE2
      NONOUI = 1
      GOTO 9999

C     LES 4 SOMMETS DE NOSOTE1 NE SONT PAS LES 4 SOMMETS DE NOSOTE2
 9000 NONOUI = 0

 9999 RETURN
      END
