      SUBROUTINE XYZTTT( U, V, XYZ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    T : TRIANGLE RECTANGLE UNITE => TRIANGLE ALGEBRIQUE
C -----    CALCULER T(U,V)=(X,Y,Z)
C          ICI LES ARETES DES 3 COTES DU TRIANGLE SONT SUPPOSEES DROITES (P1)
C          SEUL LE SP TRNSAL PEUT EXECUTER CE SP CAR DE NOMBREUX
C          PARAMETRES SONT CACHES DANS LE COMMON /S09S01/
C
C ENTREES:
C --------
C U,V    : LES 2 COORDONNEES DU POINT DU TRIANGLE RECTANGLE UNITE
C
C SORTIES:
C --------
C XYZ    : LES 3 COORDONNEES DU POINT IMAGE SUR LE TRIANGLE ALGEBRIQUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS            MARS 1998
C234567..............................................................012
      include"./incl/a___xyzsommet.inc"
C     ATTENTION: LE COMMON SUIVANT DOIT ETRE INITIALISE!
C                CF LE SP TRNSAL
      COMMON /S09S01/  MNCOCU, NBSOCT(4), NUCOTE(4), MNSOCT(4),
     %                 XYZ3ST(3,4), XYZL(3,4), CB(4), ST(2,4), UNMOT
      REAL             XYZ(3), XYZCOT(3,3,2)
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (RMCN(1),MCN(1))
C
C     ( X, Y, Z ) = T( U, V )
C     PASSAGE DES COORDONNEES DU POINT DU TRIANGLE RECTANGLE UNITE
C     AUX COORDONNEES DU POINT DU TRIANGLE ALGEBRIQUE DANS R**3
C     LES 3 COORDONNEES BARYCENTRIQUES
      CB(2) = U
      CB(3) = V
      CB(1) = 1.0 - CB(2) - CB(3)
C
      DO 20 K=1,3
         IF( ABS(CB(K)-1.0) .LE. 5.E-5 ) THEN
C           C'EST LE SOMMET K
            XYZ(1) = XYZ3ST(1,K)
            XYZ(2) = XYZ3ST(2,K)
            XYZ(3) = XYZ3ST(3,K)
            RETURN
         ENDIF
 20   CONTINUE
C
C     ADRESSE DU DEBUT DE LA VALEUR DU PARAMETRE SUR CHACUN DES 3 COTES
      MNC = MNCOCU
      N   = 0
      DO 40 I=1,3
C
C        LE NUMERO EFFECTIF DU I-EME COTE DANS NUCOTE
         NC = ABS( NUCOTE(I) )
C
C        RECHERCHE DE L'INTERVALLE I0 CONTENANT CB(IP) SUR LE COTE I
         IF( I .EQ. 3 ) THEN
            IP = 1
         ELSE
            IP = I + 1
         ENDIF
C
         DO 35 L=1,2
            IF( L .EQ. 1 ) THEN
               COOBA = CB(IP)
            ELSE
               COOBA = 1.0 - CB(I)
            ENDIF
C
C           RECHERCHE DE L'INTERVALLE I0 CONTENANT COOBA SUR LE COTE I
            CALL INTCOT( COOBA, NBSOCT(NC)-1, RMCN(MNC), I0 )
C
            IF( NUCOTE(I) .GT. 0 ) THEN
C              SENS DU CONTOUR FERME
               MNL0 = MNSOCT(NC) + WYZSOM + 3 * I0
               MNL1 = MNL0 + 3
            ELSE
C              SENS INVERSE DU CONTOUR FERME
               MNL0 = MNSOCT(NC) + WYZSOM + 3 * (NBSOCT(NC) - 1 - I0)
               MNL1 = MNL0 - 3
            ENDIF
C
C           ADRESSE DU PARAMETRE DES 2 EXTREMITES DE L'INTERVALLE I0
C           DU PARAMETRE DU COTE I
            MNC0 = MNCOCU + N + I0
            MNC1 = MNC0 + 1
C
C           LES 3 COORDONNEES DU POINT SUR LE COTE I DU TRIANGLE ALGEBRIQUE
            DO 30 K=1,3
               XYZCOT(K,I,L) = ( (COOBA-RMCN(MNC0)) * RMCN(MNL1+K-1)
     %                         + (RMCN(MNC1)-COOBA) * RMCN(MNL0+K-1) )
     %                       / ( RMCN(MNC1) - RMCN(MNC0) )
 30         CONTINUE
C
 35      CONTINUE
C
C        PASSAGE AU COTE SUIVANT
         MNC = MNC + NBSOCT(NC)
         N   = N   + NBSOCT(NC)
 40   CONTINUE
C
C     LES 3 COORDONNEES DU SOMMET SUR LE TRIANGLE COURBE
      DO 50 K=1,3
         XYZ(K) = CB(1)*( XYZCOT(K,1,1) + XYZCOT(K,3,2) - XYZ3ST(K,1) )
     %          + CB(2)*( XYZCOT(K,2,1) + XYZCOT(K,1,2) - XYZ3ST(K,2) )
     %          + CB(3)*( XYZCOT(K,3,1) + XYZCOT(K,2,2) - XYZ3ST(K,3) )
 50   CONTINUE
      END
