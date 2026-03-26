      SUBROUTINE CENCED( XY1, XY2, XY3, CETRIA, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DES COORDONNEES DU CENTRE DU CERCLE CIRCONSCRIT
C ----- DU TRIANGLE DEFINI PAR SES 3 SOMMETS DE COORDONNEES
C       XY1 XY2 XY3 AINSI QUE LE CARRE DU RAYON DE CE CERCLE
C
C ENTREES :
C ---------
C XY1 XY2 XY3 : LES 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C IERR   : <0  => PAS D'AFFICHAGE SI TRIANGLE DEGENERE
C          >=0 =>       AFFICHAGE SI TRIANGLE DEGENERE
C
C SORTIE :
C --------
C CETRIA : CETRIA(1)=ABCISSE  DU CENTRE
C          CETRIA(2)=ORDONNEE DU CENTRE
C          CETRIA(3)=CARRE DU RAYON   1D28 SI TRIANGLE DEGENERE
C IERR   : 0 SI TRIANGLE NON DEGENERE
C          1 SI TRIANGLE DEGENERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS        JUIN 1995
C2345X7..............................................................012
      DOUBLE PRECISION  EPSURF
      PARAMETER        (EPSURF=1D-7)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  X1,Y1,X21,Y21,X31,Y31,
     %                  AIRE2,XC,YC,ROT,
     %                  XY1(2),XY2(2),XY3(2),CETRIA(3)
C
C     LE CALCUL DE 2 FOIS L'AIRE DU TRIANGLE
C     ATTENTION L'ORDRE DES 3 SOMMETS EST DIRECT OU NON
      X1  = XY1(1)
      X21 = XY2(1) - X1
      X31 = XY3(1) - X1
C
      Y1  = XY1(2)
      Y21 = XY2(2) - Y1
      Y31 = XY3(2) - Y1
C
      AIRE2  = X21 * Y31 - X31 * Y21
C
C     RECHERCHE D'UN TEST RELATIF PEU COUTEUX
C     POUR REPERER LA DEGENERESCENCE DU TRIANGLE
      IF( ABS(AIRE2) .LE.
     %    EPSURF*(ABS(X21)+ABS(X31))*(ABS(Y21)+ABS(Y31)) ) THEN
C        TRIANGLE DE QUALITE TROP FAIBLE
         IF( IERR .GE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'ERREUR CENCED: TRIANGLE DEGENERE'
            CALL LEREUR
            WRITE(IMPRIM,10000)  XY1,XY2,XY3,AIRE2
         ENDIF
10000 FORMAT( 3(' X=',G24.16,' Y=',G24.16/),' AIRE*2=',G24.16)
         CETRIA(1) = 0D0
         CETRIA(2) = 0D0
         CETRIA(3) = 1D28
         IERR = 1
         RETURN
      ENDIF
C
C     LES 2 COORDONNEES DU CENTRE INTERSECTION DES 2 MEDIATRICES
C     X = (X1+X2)/2 + LAMBDA * (Y2-Y1)
C     Y = (Y1+Y2)/2 - LAMBDA * (X2-X1)
C     X = (X1+X3)/2 + ROT    * (Y3-Y1)
C     Y = (Y1+Y3)/2 - ROT    * (X3-X1)
C     ==========================================================
      ROT = ((XY2(1)-XY3(1))*X21 + (XY2(2)-XY3(2))*Y21) / (2 * AIRE2)
C
      XC = ( X1 + XY3(1) ) * 0.5D0 + ROT * Y31
      YC = ( Y1 + XY3(2) ) * 0.5D0 - ROT * X31
C
      CETRIA(1) = XC
      CETRIA(2) = YC
C
C     LE CARRE DU RAYON
      CETRIA(3) = (X1-XC) ** 2 + (Y1-YC) ** 2
C
C     PAS D'ERREUR RENCONTREE
      IERR = 0

      RETURN
      END
