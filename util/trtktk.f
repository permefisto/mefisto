      SUBROUTINE TRTKTK( M1 , M2 , L )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRANSFERT DU TABLEAU M1 DE L MOTS DANS LE TABLEAU M2
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS    NOVEMBRE 1983
C2345X7..............................................................012
      INTEGER  M1(1),M2(1)

      DO I = 1,L
         M2(I) = M1(I)
      ENDDO

      RETURN
      END
