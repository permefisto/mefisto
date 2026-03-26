      SUBROUTINE BSPARA ( NUTYSU , LUX , LUY , POINTS ,
     %                    NBAXQB , NBAYQB , RAGXQB , RAGYQB , UX , UY )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER DES TABLEAUX DE PARAMETRES POUR UNE SURFACE
C -----    B-SPLINE DEFINIE PAR DES POINTS OU DES LIGNES B-SPLINE.
C
C ENTREES:
C --------
C NUTYSU : =3 CORRESPOND A UNE SURFACE DEFINIE PAR DES POINTS.
C          =4 CORRESPOND A UNE SURFACE DEFINIE PAR DES LIGNES B-SPLINE.
C LUX    : NOMBRE DE NOEUDS D'INTERPOLATION EN X - 1
C LUY    : NOMBRE DE NOEUDS D'INTERPOLATION EN Y - 1
C POINTS : LES 3 COORDONNEES DES POINTS D'INTERPOLATION (NUTYSU=3).
C NBAXQB : NOMBRE D'ARETES EN X
C NBAYQB : NOMBRE D'ARETES EN Y
C RAGXQB : RAISON DE LA PROGRESSION GEOMETRIQUE D'ESPACEMENT EN X
C RAGYQB : RAISON DE LA PROGRESSION GEOMETRIQUE D'ESPACEMENT EN Y
C
C SORTIES:
C --------
C UX     : TABLEAU DES PARAMETRES EN X
C UY     : TABLEAU DES PARAMETRES EN Y
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS : SINOQUET FRENOD   ANALYSE NUMERIQUE UPMC  PARIS FEVRIER 1991
C2345X7.................................................................
      IMPLICIT     INTEGER(W)
      include"./incl/a_ligne__bspline.inc"
      include"./incl/pp.inc"
      COMMON       MCN(MOTMCN)
      REAL         RMCN(1)
      EQUIVALENCE (MCN(1),RMCN(1))
      REAL         POINTS(1:LUX,1:LUY,1:3),
     %             UX(0:LUX),
     %             UY(0:LUY)
C
      IF ( NUTYSU .EQ. 3 ) THEN
C
C         SURFACE DEFINIE PAR DES POINTS
C         ------------------------------
          COEF  = 0.
          UX(0) = 0
          IF( RAGXQB .EQ. 1 .AND. MOD(NBAXQB,LUX-1) .EQ. 0 ) THEN
C
C            UX A DES VALEURS REGULIERES POUR ASSURER LA CONTINUITE
             DO 10 IX=1,LUX
                UX( IX ) = IX
 10          CONTINUE
C
          ELSE
C
C            CALCUL DE LA LONGUEUR MOYENNE ENTRE LES POINTS SUCCESSIFS
             DO 1 IX=2,LUX
                COEF = 0.
                DO 2 IY=1,LUY
                   A = ( POINTS(IX,IY,1) - POINTS(IX-1,IY,1) )**2
     %               + ( POINTS(IX,IY,2) - POINTS(IX-1,IY,2) )**2
     %               + ( POINTS(IX,IY,3) - POINTS(IX-1,IY,3) )**2
C                  CALCUL DE LA MOYENNE (COEF)
                   COEF = COEF + SQRT( A )
    2           CONTINUE
                UX(IX-1) = UX(IX-2) + COEF / LUY
    1        CONTINUE
C            VALEUR A INITIALISER POUR EVITER UNE ERREUR DANS UN CALCUL
             UX(LUX) = UX(LUX-1) + COEF
          ENDIF
C
          UY(0) = 0
          IF( RAGYQB .EQ. 1 .AND. MOD(NBAYQB,LUY-1) .EQ. 0 ) THEN
C
C            UY A DES VALEURS REGULIERES POUR ASSURER LA CONTINUITE
             DO 11 IY=1,LUY
                UY( IY ) = IY
 11          CONTINUE
C
          ELSE
C
C            CALCUL DE LA LONGUEUR MOYENNE ENTRE LES POINTS SUCCESSIFS
             DO 3 IY=2,LUY
                COEF = 0.
                DO 4 IX=1,LUX
                   A = ( POINTS(IX,IY,1) - POINTS(IX,IY-1,1) )**2
     %               + ( POINTS(IX,IY,2) - POINTS(IX,IY-1,2) )**2
     %               + ( POINTS(IX,IY,3) - POINTS(IX,IY-1,3) )**2
C                  CALCUL DE LA MOYENNE (COEF)
                   COEF = COEF + SQRT( A )
    4           CONTINUE
                UY(IY-1) = UY(IY-2) + COEF / LUX
    3        CONTINUE
C            VALEUR A INITIALISER POUR EVITER UNE ERREUR DANS UN CALCUL
             UY(LUY) = UY(LUY-1) + COEF
          ENDIF
C
      ELSE IF( NUTYSU .EQ. 4 ) THEN
C
C         SURFACE DEFINIE PAR DES LIGNES B-SPLINE
C         ---------------------------------------
C         CALCUL DES COEFFICIENTS EN X
          DO 5 IX=0,LUX
             UX(IX) = IX
 5        CONTINUE
C         CALCUL DES COEFFICIENTS EN Y
          DO 7 IY=0,LUY
             UY(IY) = IY
 7        CONTINUE
C
      ENDIF
      END
