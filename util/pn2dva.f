      SUBROUTINE PN2DVA( N1, P, X, Y, VAPXY )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : VALEUR DU POLYNOME PRODUIT P AU POINT X , Y
C -----
C
C PARAMETRES D ENTREE :
C ---------------------
C N1     : DEGRE+1 DU POLYNOME P
C P      : TABLEAU A 2 INDICES CONTENANT LES COEFFICIENTS DU POLYNOME P
C          P(I,J)=COEFFICIENT DE X**(I-1) Y**(J-1)
C X,Y    : COORDONNEES DU POINT OU P DOIT ETRE CALCULE
C
C PARAMETRE-RESULTAT :
C --------------------
C VAPXY  : LA VALEUR DE P AU POINT X , Y ,Z
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1979
C23456---------------------------------------------------------------012
      DOUBLE PRECISION P(N1,N1), VAPXY, X, Y, XI1, YI1
C
      VAPXY = 0.D0
      XI1 = 1.D0
      DO 1 I=1,N1
         YI1 = 1.D0
         DO 2 J=1,N1
            VAPXY = VAPXY + P(I,J) * XI1 * YI1
            YI1 = YI1 * Y
    2    CONTINUE
         XI1 = XI1 * X
    1 CONTINUE
C
      IF (DABS(VAPXY) .LT. 1.D-14) VAPXY = 0.D0
      RETURN
      END
