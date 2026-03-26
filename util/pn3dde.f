      SUBROUTINE PN3DDE( NODERI, N1, P, DP )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECHERCHE DES COEFFICIENTS DU POLYNOME DERIVE DP / NODERI DE P
C -----
C
C PARAMETRES D ENTREE :
C ---------------------
C NODERI : 1 DP EST LE POLYNOME DERIVE PAR RAPPORT A X DU POLYNOME P
C          2 DP EST LE POLYNOME DERIVE PAR RAPPORT A Y DU POLYNOME P
C          3 DP EST LE POLYNOME DERIVE PAR RAPPORT A Z DU POLYNOME P
C N1     : DEGRE+1 DU POLYNOME P EN CHACUNE DES VARIABLES
C P      : TABLEAU A 3 INDICES CONTENANT LES COEFFICIENTS DU POLYNOME P
C          P(I,J,K)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C
C PARAMETRE-RESULTAT :
C --------------------
C DP     : LE POLYNOME DERIVE DE P PAR RAPPORT A X(NODERI)
C
C ATTENTION : DP DOIT ETRE DIFFERENT DE P
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET LAN189 PARIS ET IRIA  OCTOBRE 1979
C ......................................................................
      include"./incl/gsmenu.inc"
      DOUBLE PRECISION P(N1,N1,N1), DP(N1,N1,N1)
C
      IF(N1 .LE. 1) GOTO 1000
      GOTO( 100 , 200 , 300 ), NODERI
C
C     DERIVATION SELON X
C
 100  DO 10 K=1,N1
         DO 20 J=1,N1
            DO 30 I=2,N1
               M1 = I - 1
               DP(M1,J,K) = M1 * P(I,J,K)
 30         CONTINUE
            DP(N1,J,K) = 0.D0
 20      CONTINUE
 10   CONTINUE
      GOTO 900
C
C     DERIVATION SELON Y
C
 200  DO 11 K=1,N1
         DO 21 I=1,N1
            DO 31 J=2,N1
               M1 = J - 1
               DP(I,M1,K) = M1 * P(I,J,K)
 31         CONTINUE
            DP(I,N1,K) = 0.D0
 21      CONTINUE
 11   CONTINUE
      GOTO 900
C
C     DERIVATION SELON Z
C
 300  DO 12 J=1,N1
         DO 22 I=1,N1
            DO 32 K=2,N1
               M1 = K - 1
               DP(I,J,M1) = M1 * P(I,J,K)
 32         CONTINUE
            DP(I,J,N1) = 0.D0
 22      CONTINUE
 12   CONTINUE
      GOTO 900
C
C     ERREURS
 1000 NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') N1
      KERR(1) = 'DEGRE INCORRECT=' // KERR(MXLGER)(1:4)
      CALL LEREUR
C
 900  RETURN
      END
