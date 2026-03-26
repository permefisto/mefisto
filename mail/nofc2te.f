      SUBROUTINE  NOFC2TE( NT1, NT2, NBSOTE, NOSOTE, NF1, NF2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LE NUMERO LOCAL DE LA FACE COMMUNE AUX 2 TETRAEDRES
C -----    NT1 NT2 DE SOMMETS DEFINIS DANS NOSOTE
C
C ENTREES:
C --------
C NT1    : NUMERO DU PREMIER TETRAEDRE
C NT2    : NUMERO DU SECOND  TETRAEDRE
C NOSOTE : NUMERO DES 4 SOMMETS DES TETRAEDRES
C
C SORTIES:
C --------
C NF1 : NUMERO DE 1 A 4 DE LA FACE DU TETRAEDRE NT1 SI COMMUNE
C NF2 : NUMERO DE 1 A 4 DE LA FACE DU TETRAEDRE NT2 SI COMMUNE
C      FACE 1:1,2,3   2:2,4,3   3:3,4,1   4:4,2,1
C      0 SI FACE NON RETROUVEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JANVIER 1992
C....................................................................012
      INTEGER      NOSOTE(1:NBSOTE,*),NS(3)

C     PARCOURS DES 4 FACES DE NT1
      DO 10 NF1=1,4

C        LES 3 SOMMETS DE LA FACE NF1 DU TETRAEDRE NT1
         IF( NF1 .EQ. 4 ) THEN
            I1 = 1
         ELSE
            I1 = NF1 + 1
         ENDIF
         IF( I1 .EQ. 4 ) THEN
            I2 = 1
         ELSE
            I2 = I1 + 1
         ENDIF
         NS(1) = NOSOTE( NF1, NT1 )
         NS(2) = NOSOTE( I1 , NT1 )
         NS(3) = NOSOTE( I2 , NT1 )
C
C        TRI CROISSANT DES 3 SOMMETS
         CALL TRI3NO( NS, NS )
C
C        NS EST ELLE LA FACE NF2 DU TETRAEDRE NT2 ?
         CALL NO1F1T( NS, NOSOTE(1,NT2), NF2 )
         IF( NF2 .GT. 0 ) RETURN
 10   CONTINUE
C
C     PAS DE FACE COMMUNE
      NF1 = 0
      END
