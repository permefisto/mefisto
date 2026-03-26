      SUBROUTINE CHSGND(N,TABLE)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                          SP CHSGND
C                          ---------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: CHANGER LE SIGNE DU VECTEUR TABLE
C ----
C
C PARAMETRES D ENTREE:
C --------------------
C N      : NOMBRE DE LIGNES DU VECTEUR
C TABLE  : LE VECTEUR DONNE
C
C PARAMETRE MODIFIE  :
C --------------------
C TABLE  : LE VECTEUR RESULTAT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR: P. JOLY  ANALYSE NUMERIQUE PARIS      MAI 1990
C.......................................................................
      DOUBLE PRECISION TABLE(N)
C
      DO 1 I = 1 , N
         TABLE(I) = - TABLE(I)
1     CONTINUE
      RETURN
      END
