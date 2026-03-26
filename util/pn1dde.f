      SUBROUTINE PN1DDE ( N1,P,  DP )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECHERCHE DES COEFFICIENTS DU POLYNOME DERIVE DP DE P
C -----
C
C PARAMETRES D ENTREE :
C ---------------------
C N1     : DEGRE+1 DU POLYNOME P A UNE VARIABLE
C P      : TABLEAU A 1 INDICE CONTENANT LES COEFFICIENTS DU POLYNOME P
C          P( I )=COEFFICIENT DE X**(I-1)
C
C PARAMETRE-RESULTAT :
C --------------------
C DP     : LE POLYNOME DERIVE DE P PAR RAPPORT A X
C
C ATTENTION : DP DOIT ETRE DIFFERENT DE P
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET LAN189 PARIS ET INRIA   FEVRIER 1981
C ......................................................................
      include"./incl/gsmenu.inc"
      DOUBLE PRECISION P(N1),DP(N1)
      COMMON / UNITES / LECTEU , IMPRIM , NUNITE(30)
C
      IF( N1 .GE. 1 ) THEN
C
C        DERIVATION SELON X
C        ==================
         DO 10 I=2,N1
            I1 = I - 1
            DP( I1 ) = I1 * P( I )
   10    CONTINUE
         DP( N1 ) = 0D0
      ELSE
C
C        ERREURS
C        =======
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') N1
         KERR(1) ='PN1DDE:P DE DIMENSION N1='//KERR(MXLGER)(1:4)
         CALL LEREUR
      ENDIF
      RETURN
      END
