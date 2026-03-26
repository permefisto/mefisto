      SUBROUTINE PN3DIN(NPOUQ,N1,P, VINTEG)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INTEGRATION DU POLYNOME PRODUIT P SUR L ELEMENT DE REFERENCE
C -----
C
C PARAMETRES D ENTREE :
C ---------------------
C NPOUQ  : 0 SI INTEGRATION SUR LE TETRAEDRE RECTANGLE UNITE
C        : 1 SI INTEGRATION SUR L  HEXAEDRE  RECTANGLE UNITE
C        : 2 SI INTEGRATION SUR LE PENTAEDRE RECTANGLE UNITE
C N1     : DEGRE+1 DU POLYNOME P
C P      : TABLEAU A 3 INDICES CONTENANT LES COEFFICIENTS DU POLYNOME P
C          P(I,J,K)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C
C PARAMETRE-RESULTAT :
C --------------------
C VINTEG : LA VALEUR DE L INTEGRALE DE P SUR L ELEMENT DE REFERENCE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET LAN189 PARIS ET IRIA  OCTOBRE 1979
C ......................................................................
      DOUBLE PRECISION P(N1,N1,N1), V,VINTEG, FAI1,FAJ1,FAK1,FAIJ,FAIJK
C
C     AIGUILLAGE SELON L ELEMENT DE REFERENCE
C
      VINTEG = 0.D0
      GOTO( 1000 , 2000 , 1000 ), NPOUQ
C
C     LE TETRAEDRE ET LE PENTAEDRE
C
 1000 FAI1 = 1.D0
C
C     FAI1 = FACTORIELLE(I-1)
C     FAIJ = FACTORIELLE(I+J)
C     FAIJK= FACTORIELLE(I+J+K)
C
      DO 1001 I=1,N1
         FAI1 = FAI1 * MAX0(1,I-1)
         FAJ1 = 1.D0
         FAIJ = FAI1 * I
         DO 1002 J=1,N1
            FAJ1 = FAJ1 * MAX0(1,J-1)
            V    = FAI1 * FAJ1
            FAIJ = FAIJ * (I+J)
            IF( NPOUQ .EQ. 2 ) GOTO 1100
C
C           LE TETRAEDRE
C
            FAK1 = 1.D0
            FAIJK= FAIJ
            DO 1003 K=1,N1
               FAK1   = FAK1  * MAX0(1,K-1)
               FAIJK  = FAIJK * (I+J+K)
               VINTEG = VINTEG + P(I,J,K) * V * FAK1 / FAIJK
 1003       CONTINUE
            GOTO 1002
C
C           LE PENTAEDRE
C
 1100       DO 1103 K=1,N1
               VINTEG = VINTEG + P(I,J,K) * V / (FAIJ * K)
 1103       CONTINUE
 1002    CONTINUE
 1001 CONTINUE
      GOTO 9000
C
C     L HEXAEDRE
C
 2000 DO 2001 I=1,N1
         DO 2002 J=1,N1
            DO 2003 K=1,N1
               VINTEG = VINTEG + P(I,J,K) / (I * J * K)
 2003       CONTINUE
 2002    CONTINUE
 2001 CONTINUE
C
 9000 RETURN
      END
