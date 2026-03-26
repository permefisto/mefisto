      SUBROUTINE CENCER( XY1, XY2, XY3, CETRIA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DES COORDONNEES DU CENTRE DU CERCLE CIRCONSCRIT
C ----- DU TRIANGLE DEFINI PAR SES 3 SOMMETS DE COORDONNEES
C       XY1 XY2 XY3 AINSI QUE LE CARRE DU RAYON DE CE CERCLE
C
C ENTREES :
C ---------
C XY1 XY2 XY3 : LES 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C
C SORTIE :
C --------
C CETRIA : CETRIA(1)=ABCISSE  DU CENTRE
C          CETRIA(2)=ORDONNEE DU CENTRE
C          CETRIA(3)=CARRE DU RAYON,  1E28 SI TRIANGLE DEGENERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  AVRIL 1986
C...............................................................................
      PARAMETER        (EPSURF=5E-5)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      DOUBLE PRECISION  X1, Y1, X21, Y21, X31, Y31, X31Y21, X21Y31, 
     %                  AIRE2, XA21, YA21, XC, YC, DBLE
      REAL              XY1(2), XY2(2), XY3(2), CETRIA(3)
      INTRINSIC         REAL
C
C     LE CALCUL DE 2 FOIS L'AIRE DU TRIANGLE
C     ATTENTION L'ORDRE DES 3 SOMMETS EST DIRECT OU NON
      X1  = DBLE( XY1(1) )
      Y1  = DBLE( XY1(2) )
      X21 = DBLE( XY2(1) ) - X1
      Y21 = DBLE( XY2(2) ) - Y1
      X31 = DBLE( XY3(1) ) - X1
      Y31 = DBLE( XY3(2) ) - Y1
      X31Y21 = X31 * Y21
      X21Y31 = X21 * Y31
      AIRE2  = X21Y31 - X31Y21
C
C     RECHERCHE D'UN TEST RELATIF PEU COUTEUX
C     POUR REPERER LA DEGENERESCENCE DU TRIANGLE
      XA21 = ABS(X21)
      YA21 = ABS(Y21)
      XC   = ABS(X31)
      YC   = ABS(Y31)
      XA21 = XA21 + XC
      YA21 = YA21 + YC
      IF( ABS(AIRE2) .LE. EPSURF*XA21*YA21 ) THEN
C        TRIANGLE DE QUALITE TROP FAIBLE
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR CENCER:TRIANGLE DEGENERE'
         CALL LEREUR
         WRITE(IMPRIM,10000)  XY1,XY2,XY3,AIRE2
10000 FORMAT( 3(' X=',G15.6,' Y=',G15.6)/' 2*AIRE=',G15.6)
         CETRIA(1) = 0
         CETRIA(2) = 0
         CETRIA(3) = 1E28
         RETURN
      ENDIF
C
C     LES 2 COORDONNEES DU CENTRE
      XC = ( X21Y31 * ( X1 + XY2(1) ) - X31Y21 * ( X1 + XY3(1) )
     %     - Y21 * Y31 * ( XY3(2) - XY2(2) )  ) * 0.5D0 / AIRE2
      YC = Y21 * Y21
C
      IF( YC .GE. 1D-4 * ( X21 * X21 + YC ) .OR. Y31 .EQ. 0.D0 ) THEN
         YC = X21 / Y21 * ( ( X1 + XY2(1) ) * 0.5D0 - XC )
     %        + ( Y1 + XY2(2) ) * 0.5D0
      ELSE
         YC = X31 / Y31 * ( ( X1 + XY3(1) ) * 0.5D0 - XC )
     %        + ( Y1 + XY3(2) ) * 0.5D0
      ENDIF
C
      CETRIA(1) = REAL( XC )
      CETRIA(2) = REAL( YC )
C
C     LE CARRE DU RAYON
      CETRIA(3) =  REAL( (X1-XC) ** 2 + (Y1-YC) ** 2 )
C
      RETURN
      END
