      SUBROUTINE PN2DDE( NODERI, N1, P, DP )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECHERCHE DES COEFFICIENTS DU POLYNOME DERIVE DP DE P
C ----- ATTENTION : A L'APPEL DP DOIT ETRE DIFFERENT DE P
C
C PARAMETRES D ENTREE :
C ---------------------
C NODERI : 1 DP EST LE POLYNOME DERIVE PAR RAPPORT A X DU POLYNOME P
C          2 DP EST LE POLYNOME DERIVE PAR RAPPORT A Y DU POLYNOME P
C N1     : DEGRE+1 DU POLYNOME P EN CHACUNE DES VARIABLES
C P      : TABLEAU A 2 INDICES CONTENANT LES COEFFICIENTS DU POLYNOME P
C          P(I,J)=COEFFICIENT DE X**(I-1) Y**(J-1)
C
C PARAMETRE-RESULTAT :
C --------------------
C DP     : LE POLYNOME DERIVE DE P PAR RAPPORT A X(NODERI)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1979
C ......................................................................
      include"./incl/gsmenu.inc"
      DOUBLE PRECISION  P(N1,N1), DP(N1,N1)
      COMMON / UNITES / LECTEU , IMPRIM , NUNITE(30)

      IF(N1 .LE. 1) GOTO 1000

      IF( NODERI .GE. 2 ) GOTO 200

C     DERIVATION SELON X
      DO J=1,N1
         DO I=2,N1
            I1 = I - 1
            DP(I1,J) = I1 * P(I,J)
         ENDDO
         DP(N1,J) = 0.D0
      ENDDO
      GOTO 900

C     DERIVATION SELON Y
  200 DO I=1,N1
         DO J=2,N1
            I1 = J - 1
            DP(I,I1) = I1 * P(I,J)
         ENDDO
         DP(I,N1) = 0.D0
      ENDDO

  900 RETURN

C     ERREURS
 1000 NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') N1
      KERR(1) ='PN2DDE:P DE DIMENSION INCORRECTE N1='//KERR(MXLGER)(1:4)
      CALL LEREUR

      STOP
      END
