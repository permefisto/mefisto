      REAL FUNCTION VALP3H( U, F0, F1, T1, T2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULE L'INTERPOLATION HERMITE P3 CROISEE POUR LA VALEUR DU
C -----    PARAMETRE U DE [0,1] ET LES DEGRES DE LIBERTE
C          LES 2 POINTS F0,F1 ET LES DERIVEES T1,T2 EN CES 2 POINTS
C          C'EST A DIRE
C                     VALP3H(0) = F0 = F(0)
C                     VALP3H(1) = F1 = F(1)
C                     d VALP3H / du (0)(1-0) = T1 = DF(0) (1-0)
C                     d VALP3H / du (1)(0-1) = T2 = DF(1) (0-1) = -DF(1)
C          ATTENTION: d VALP3H / du (1) = -T2 ...!
C
C ENTREES:
C --------
C U      : VALEUR DU PARAMETRE (COMPRISE ENTRE 0 ET 1)
C F0,F1  : F(0)
C F1     : F(1)
C T1,T2  : DF(0) (1-0)
C T2     : DF(1) (0-1)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS        MAI 1996
C2345X7..............................................................012
      U2 = U * U
      U3 = U * U2
      V  = 3*U2 - 2*U3
      VALP3H = V * (F1 - F0) + F0 + (U3-U2) * (T1 - T2) + T1 * (U - U2)
      END
