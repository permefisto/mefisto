      SUBROUTINE XYZTRT( U, V, XYZ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    T : CARRE => QUADRANGLE TRANSFINI  T(U,V)=(X,Y,Z)
C -----    CALCULER T(U,V)=(X,Y,Z)
C          ICI LES ARETES DES 4 COTES DU QUADRANGLE SONT SUPPOSEES
C          DROITES (P1)
C
C ATTENTION: LES VARIABLES DU COMMON /S09S01/ ...
C            DOIVENT ETRE INITIALISEES!    CF LE SP QUNSAL
C
C ENTREES:
C --------
C U,V    : LES 2 COORDONNEES DU POINT DU CARRE UNITE
C
C SORTIES:
C --------
C XYZ    : LES 3 COORDONNEES DU POINT IMAGE SUR LE QUADRANGLE TRANSFINI
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS         FEVRIER 1997
C234567..............................................................012
      include"./incl/a___xyzsommet.inc"
C     ATTENTION: LE COMMON SUIVANT DOIT ETRE INITIALISE!
C                CF LE SP QUNSAL
      COMMON /S09S01/  MNCOCU, NBSOCQ(4), NUCOTQ(4), MNSOCQ(4),
     %                 XYZ4ST(3,4), XYZL(3,4), XY(4), ST(2,4), UNMOT
      REAL             XYZ(3)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (RMCN(1),MCN(1))
C
C     ( X, Y, Z ) = T( U, V )
C     PASSAGE DES COORDONNEES DU POINT DU CARRE UNITE
C     AUX COORDONNEES DU POINT DU QUADRANGLE TRANSFINI DANS R**3
      XY(1) = U
      XY(2) = V
C
      XY(3) = 1.0 - XY(1)
      XY(4) = 1.0 - XY(2)
C
C     ADRESSE DU DEBUT DE LA VALEUR DU PARAMETRE SUR CHACUN DES 4 COTES
      MNC = MNCOCU
      N   = 0
      DO 40 I=1,4
C
C        LE NUMERO EFFECTIF DU I-EME COTE DANS NUCOTQ
         NC = ABS( NUCOTQ(I) )
C
C        RECHERCHE DE L'INTERVALLE I0 CONTENANT XY(I) SUR LE COTE I
         CALL INTCOT( XY(I), NBSOCQ(NC)-1, RMCN(MNC), I0 )
C
C        CALCUL DE L'ADRESSE DES XYZ DES 2 POINTS EXTREMITES DE L'ARETE
C        SUR LE QUADRANGLE TRANSFINI CORRESPONDANT A l'INTERVALLE I0 DU COTE I
         IF( NUCOTQ(I) .GT. 0 ) THEN
C           SENS DU CONTOUR FERME
            MNL0 = MNSOCQ(NC) + WYZSOM + 3 * I0
            MNL1 = MNL0 + 3
         ELSE
C           SENS INVERSE DU CONTOUR FERME
            MNL0 = MNSOCQ(NC) + WYZSOM + 3 * (NBSOCQ(NC) - 1 - I0)
            MNL1 = MNL0 - 3
         ENDIF
C
C        ADRESSE DU PARAMETRE DES 2 EXTREMITES DE L'INTERVALLE I0
C        DU PARAMETRE DU COTE I
         MNC0 = MNCOCU + N + I0
         MNC1 = MNC0 + 1
C
C        LES 3 COORDONNEES DU POINT SUR LE COTE I DU QUADRANGLE COURBE
         DO 30 K=1,3
            XYZL(K,I) = (  (XY(I)-RMCN(MNC0)) * RMCN(MNL1+K-1)
     %                   + (RMCN(MNC1)-XY(I)) * RMCN(MNL0+K-1) )
     %                  /  ( RMCN(MNC1) - RMCN(MNC0) )
 30      CONTINUE
C
C        PASSAGE AU COTE SUIVANT
         MNC = MNC + NBSOCQ(NC)
         N   = N   + NBSOCQ(NC)
 40   CONTINUE
C     LES 3 COORDONNEES (X,Y,Z) DU POINT SUR LE QUADRANGLE TRANSFINI
      DO 50 K=1,3
         XYZ(K) = XY(3) * ( XYZL(K,4) - XY(2)*XYZ4ST(K,4) )
     %          + XY(1) * ( XYZL(K,2) - XY(4)*XYZ4ST(K,2) )
     %          + XY(4) * ( XYZL(K,1) - XY(3)*XYZ4ST(K,1) )
     %          + XY(2) * ( XYZL(K,3) - XY(1)*XYZ4ST(K,3) )
 50   CONTINUE
      END
