C/UPDATE ADD NAME=PN1DPR,SSI=82121623
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                   SP PN1DPR
C                   ---------
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECHERCHE DES COEFFICIENTS DU POLYNOME PRODUIT P3 DE P1 PAR P2
C -----
C
C PARAMETRES D ENTREE :
C ---------------------
C          POUR I=1 OU 2
C NI     : DEGRE+1 DU POLYNOME PI AU SENS CI DESSUS
C PI     : TABLEAU A 1 INDICE CONTENANT LES COEFFICIENTS DU POLYNOME PI
C          PI( J )=COEFFICIENT DE X**(J-1)
C N3     : DIMENSION DU POLYNOME PRODUIT P3
C
C PARAMETRE-RESULTAT :
C --------------------
C P3     : LE POLYNOME PRODUIT DE P1 PAR P2
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET LAN189 PARIS ET INRIA  FEVRIER 1981
C ......................................................................
      SUBROUTINE PN1DPR(N1,P1,N2,P2,N3,P3)
C ......................................................................
      DOUBLE PRECISION P1(N1),P2(N2),P3(N3)
C
C     L INITIALISATION A ZERO DE P3
C
            DO 1 I=1,N3
            P3( I ) = 0D0
    1       CONTINUE
C
C     LE PRODUIT
C
            DO 20 I=1,N1
            IP = I - 1
                 DO 10 J=1,N2
                 IPP = IP + J
                 P3( IPP ) = P3( IPP ) + P1( I ) * P2( J )
   10            CONTINUE
   20       CONTINUE
C
      RETURN
      END
