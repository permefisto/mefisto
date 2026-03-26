      SUBROUTINE PN2DIN( NPOUQ, N1, P,  VINTEG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INTEGRATION DU POLYNOME PRODUIT P SUR L ELEMENT DE REFERENCE
C -----
C
C PARAMETRES D ENTREE :
C ---------------------
C NPOUQ  : 0 SI INTEGRATION SUR LE TRIANGLE RECTANGLE UNITE
C        : 1 SI INTEGRATION SUR LE CARRE UNITE
C N1     : DEGRE+1 DU POLYNOME P
C P      : TABLEAU A 2 INDICES CONTENANT LES COEFFICIENTS DU POLYNOME P
C          P(I,J)=COEFFICIENT DE X**(I-1) Y**(J-1)
C
C PARAMETRE-RESULTAT :
C --------------------
C VINTEG : LA VALEUR DE L INTEGRALE DE P SUR L ELEMENT DE REFERENCE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET LAN189 PARIS ET IRIA  OCTOBRE 1979
C ......................................................................
      DOUBLE PRECISION P(N1,N1),VINTEG,FAI1,FAJ1,FAIJ
C
C     AIGUILLAGE SELON L ELEMENT DE REFERENCE
      VINTEG = 0.D0

      IF( NPOUQ .GT. 0 ) GOTO 20

C     LE TRIANGLE
C     ------------
      FAI1 = 1.D0

C     FAI1 = FACTORIELLE(I-1)
C     FAIJ = FACTORIELLE(I+J)

      DO I=1,N1
         FAI1 = FAI1 * MAX0(1,I-1)
         FAJ1 = 1.D0
         FAIJ = FAI1 * I
         DO J=1,N1
            FAJ1 = FAJ1 * MAX0(1,J-1)
            FAIJ = FAIJ * (I+J)
            VINTEG = VINTEG + P(I,J) * FAI1 * FAJ1 / FAIJ
         ENDDO
      ENDDO
      GOTO 9000

C     LE CARRE
C     --------
 20   DO I=1,N1
         DO J=1,N1
            VINTEG = VINTEG + P(I,J) / (I * J)
         ENDDO
      ENDDO

 9000 RETURN
      END
