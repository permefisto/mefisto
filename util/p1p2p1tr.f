      SUBROUTINE P1P2P1TR( P1, P2, P1P2P1 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES INTEGRALES P1 P2 P1 dx dy SUR LE TRIANGLE UNITE
C -----
C
C ENTREES:
C --------
C P1     : POLYNOME P1(0:1,0:1,3)
C P2     : POLYNOME P2(0:2,0:2,6)
C
C SORTIE :
C --------
C P1P2P1 : TABLEAU DES INTEGRALES P1i P2j P1k dx dy
C          SUR LE TRIANGLE RECTANGLE UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M University at QATAR     Mars 2012
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  P1(0:1,0:1,3), P2(0:2,0:2,6), P1P2P1(3,6,3)
      DOUBLE PRECISION  P(1000), Q(1000)
C
      DO K = 1, 3
C
         DO J = 1, 6
C
C           PRODUIT P2j P1k
            CALL PN2DPR( 0,2+1,P2(0,0,J), 0,1+1,P1(0,0,K), 2+1+1,P )
C
            DO I = 1, 3
C
C              PRODUIT P1i P2j P1k
               CALL PN2DPR( 0,1+1,P1(0,0,I), 0,2+1+1,P, 1+1+2+1,Q )
C
C              INTEGRALE P1i P2j P1k dx dy
               CALL PN2DIN( 0,1+1+2+1,Q, P1P2P1(I,J,K) )
               IF( ABS(P1P2P1(I,J,K)) .LT. 1D-14 ) P1P2P1(I,J,K)=0D0
C
            ENDDO
C
         ENDDO
C
      ENDDO
C
      PRINT *
      PRINT 10000
10000 FORMAT('INTEGRALE sur TRIANGLE UNITE de P1 P2 P1 dx dy')
C
C
      DO K = 1, 3
         DO J = 1, 6
            PRINT 10001,(I,J,K,P1P2P1(I,J,K),I=1,3)
         ENDDO
      ENDDO
10001 FORMAT( 3('P1P2P1(',I1,',',I1,',',I1,')=',D25.17,'   ') )
C
      PRINT 10002,(((P1P2P1(I,J,K),I=1,3),J=1,6),K=1,3)
10002 FORMAT('     %',D25.17,',',D25.17,',')
C
      RETURN
      END
