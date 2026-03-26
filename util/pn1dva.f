C/UPDATE ADD NAME=PN1DVA,SSI=82121623
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                   SP PN1DVA
C                   ---------
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : VALEUR DU POLYNOME PRODUIT P AU POINT X
C -----
C
C PARAMETRES D ENTREE :
C ---------------------
C N1     : DEGRE+1 DU POLYNOME P
C P      : TABLEAU A 1 INDICE CONTENANT LES COEFFICIENTS DU POLYNOME P
C          P( I )=COEFFICIENT DE X**(I-1)
C X      : COORDONNEE DU POINT OU P DOIT ETRE CALCULE
C
C PARAMETRE-RESULTAT :
C --------------------
C VAPX   : LA VALEUR DE P AU POINT X
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET LAN189 PARIS ET INRIA FEVRIER 1981
C ......................................................................
      SUBROUTINE PN1DVA(N1,P,X, VAPX)
C ......................................................................
      DOUBLE PRECISION P(N1),VAPX,X,DABS
C
      VAPX = P( N1 )
      IF( N1 .LE. 1 ) GOTO 100
      N    = N1 - 1
            DO 10 I=1,N
            VAPX = VAPX * X + P( N1 - I )
   10       CONTINUE
C
      IF (DABS(VAPX) .LT. 1.D-14) VAPX = 0.D0
  100 RETURN
      END
