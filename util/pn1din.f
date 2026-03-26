C/UPDATE ADD NAME=PN1DIN,SSI=82121623
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                   SP PN1DIN
C                   ---------
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INTEGRATION DU POLYNOME PRODUIT P SUR L ELEMENT DE REFERENCE
C -----
C
C PARAMETRES D ENTREE :
C ---------------------
C N1     : DEGRE+1 DU POLYNOME P
C P      : TABLEAU A 1 INDICE CONTENANT LES COEFFICIENTS DU POLYNOME P
C          P( I )=COEFFICIENT DE X**(I-1)
C
C PARAMETRE-RESULTAT :
C --------------------
C VINTEG : LA VALEUR DE L INTEGRALE DE P SUR L ELEMENT DE REFERENCE(0,1)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET LAN189 PARIS ET INRIA FEVRIER 1981
C ......................................................................
      SUBROUTINE PN1DIN (N1,P, VINTEG)
C ......................................................................
      DOUBLE PRECISION P(N1),VINTEG
C
      VINTEG = 0.D0
           DO 1 I=1,N1
           VINTEG = VINTEG + P( I ) / I
    1      CONTINUE
      RETURN
      END
