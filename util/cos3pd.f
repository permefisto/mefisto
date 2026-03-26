      DOUBLE PRECISION FUNCTION COS3PD( P1 , P2 , P3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE COSINUS DE L'ANGLE DEFINI PAR 3 POINTS DE R**3
C -----
C ENTREES:
C --------
C P1 P2 P3 : LES 3 COORDONNEES DES 3 POINTS
C
C SORTIE :
C --------
C COS3PD : COSINUS DE L'ANGLE( P1-P2 , P1-P3 )
C          -10.D0 SI P1 EST CONFONDU AVEC P2 OU P3
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS OCTOBRE 1987
C..............................................................................
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  P1(3),P2(3),P3(3)
      DOUBLE PRECISION  PS,D12,D13,X12,X13
C
C     COSINUS = P1P2  PRODUIT SCALAIRE   P1P3
C               -----------------------------
C               LONGUEUR P1P2 * LONGUEUR P1P3
      PS  = 0.D0
      D12 = 0.D0
      D13 = 0.D0
C
      DO 10 I=1,3
         X12 = P2(I) - P1(I)
         X13 = P3(I) - P1(I)
         PS  = PS  + X12 * X13
         D12 = D12 + X12 ** 2
         D13 = D13 + X13 ** 2
 10   CONTINUE
C
      COS3PD = D12 * D13
      IF( COS3PD .LE. 0.D0 ) THEN
C        ERREUR  P1 CONFONDU AVEC P2 OU P3
         NBLGRC(NRERR) = 1
         KERR(1) ='ERREUR COS3PD: 2 POINTS CONFONDUS SUR 3'
         CALL LEREUR
         WRITE(IMPRIM,*) ' P1=',P1,' P2=',P2,' P3=',P3
         COS3PD = -10.D0
      ELSE
         COS3PD = PS / ( SQRT( COS3PD ) )
      ENDIF
      END
