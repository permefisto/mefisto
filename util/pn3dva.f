      SUBROUTINE PN3DVA( N1, P, X, Y, Z,  VAPXYZ )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : VALEUR DU POLYNOME P AU POINT X, Y, Z
C -----
C
C PARAMETRES D ENTREE :
C ---------------------
C N1     : DEGRE+1 DU POLYNOME P
C P      : TABLEAU A 3 INDICES CONTENANT LES COEFFICIENTS DU POLYNOME P
C          P(I,J,K)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C X,Y,Z  : COORDONNEES DU POINT OU P DOIT ETRE CALCULE
C
C PARAMETRE-RESULTAT :
C --------------------
C VAPXYZ : LA VALEUR DE P AU POINT X, Y, Z
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1995
C ......................................................................
      DOUBLE PRECISION P(N1,N1,N1),VAPXYZ,X,Y,Z,XI1,YI1,ZI1,XYI1
C
      VAPXYZ = 0.D0
      XI1 = 1.D0
      DO 1 I=1,N1
         YI1 = 1.D0
         DO 2 J=1,N1
            XYI1 = XI1 * YI1
            ZI1  = 1.D0
            DO 3 K=1,N1
               VAPXYZ = VAPXYZ + P(I,J,K) * XYI1 * ZI1
               ZI1    = ZI1    * Z
 3          CONTINUE
            YI1 = YI1 * Y
 2       CONTINUE
         XI1 = XI1 * X
 1    CONTINUE
C
      IF (DABS(VAPXYZ) .LT. 1.D-14) VAPXYZ = 0.D0
C
      RETURN
      END
