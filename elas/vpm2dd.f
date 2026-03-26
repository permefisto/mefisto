      SUBROUTINE VPM2DD( AS, VP1,VP2, XV1,YV1 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DES VALEURS ET VECTEURS PROPRES D-UNE MATRICE 2X2
C ----- SYMETRIQUE EN DOUBLE PRECISION
C
C ENTREES:
C --------
C AS        : LES COEFFICIENTS DE LA MATRICE
C                        (  AS(1)  ,  AS(2)  )
C             AS    =    (                   )
C                        (  AS(2)  ,  AS(3)  )
C
C SORTIES:
C --------
C VP1,VP2   : LES VALEURS PROPRES
C XV1,YV1   : COMPOSANTES DU VECTEUR PROPRE ASSOCIE A LA PLUS GRANDE
C             VALEUR PROPRE EN MODULE ET DE NORME 1
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN  PERRONNET  ANALYSE NUMERIQUE PARIS 6  SEPTEMBRE 1984
C ....................................................................
      include"./incl/gsmenu.inc"
      DOUBLE PRECISION  EPS, EPSN
      PARAMETER        (EPS=1D-12, EPSN=1D-20)
      DOUBLE PRECISION  AS(3),VP1,VP2,XV1,YV1
      DOUBLE PRECISION  A(3),DELTA,BB,CC,DNORM
C
C     NORMALISATION DES COEFFICIENTS DE LA MATRICE
      DNORM = 0D0
      DO 1 I=1,3
         DNORM = MAX( DNORM, ABS( AS(I) ) )
 1    CONTINUE
C
C     SI LES COEFFICIENTS DE AS SONT TOUS TRES FAIBLES => VP NULLES
C     =============================================================
 3    IF( DNORM .LE. EPSN ) THEN
         VP1 = 0.D0
         VP2 = 0.D0
         XV1 = 1.D0
         YV1 = 0.D0
         RETURN
      ENDIF
C
C     MISE A L'ECHELLE
      A(1) = AS(1) / DNORM
      A(2) = AS(2) / DNORM
      A(3) = AS(3) / DNORM
C
C     LES RACINES DU TRINOME AA * X * X + BB * X + CC = 0 SONT CALCULEES
C     ==================================================================
C     AA    = 1.D0
      BB    = -(A(1) + A(3))
      CC    = A(1) * A(3) - A(2) * A(2)
      DELTA = BB * BB - 4.D0 * CC
C
      IF( DELTA .LT. 0.D0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = ' ERREUR: VALEUR PROPRE COMPLEXE'
         CALL LEREUR
         DNORM = 0D0
         GOTO 3
      ENDIF
C
      DELTA = SQRT( DELTA )
      IF( ABS(BB) .LE. EPS ) THEN
C
C        RACINE DOUBLE
         VP1 =  DELTA * 0.5D0
         VP2 =  VP1
         GOTO 20
C
      ENDIF
C
      IF( BB .LT. 0.D0 ) THEN
C
C        LA VALEUR DE PLUS GRANDE NORME (VALEUR POSITIVE)
         VP1 = (-BB + DELTA) * 0.5D0
C
      ELSE
C
C        LA VALEUR DE PLUS GRANDE NORME (VALEUR NEGATIVE)
         VP1 = (-BB - DELTA) * 0.5D0
C
      ENDIF
C
C     LA SECONDE VALEUR PROPRE PAR LE PRODUIT DES RACINES DU TRINOME
      VP2 = CC / VP1
C
C     CALCUL DU PREMIER VECTEUR PROPRE
C     ================================
 20   BB  = A(1) - VP1
C     REMISE A L'ECHELLE INVERSE
      VP1 = VP1 * DNORM
      VP2 = VP2 * DNORM
C
      IF( ABS( ABS(BB) - ABS(A(2)) ) .LT. EPS ) THEN
C
C        LES 2 COEFFICIENTS DE BB XV1 + A(2) YV1 = 0 SONT EGAUX
         IF( ABS(VP2) .LE. ABS(VP1) ) THEN
            XV1 = 1.D0
            YV1 = 0.D0
         ELSE
            XV1 = 0.D0
            YV1 = 1.D0
         ENDIF
         RETURN
C
      ELSE
C
         IF( ABS(BB) .LT. ABS(A(2)) ) THEN
C
C           DIVISION PAR LE PLUS GRAND DES 2 COEFFICIENTS
            XV1 = 1.D0
            YV1 = -BB / A(2)
C
         ELSE
C
C           DIVISION PAR LE PLUS GRAND DES 2 COEFFICIENTS
            XV1 = -A(2) / BB
            YV1 = 1.D0
C
         ENDIF
      ENDIF
C
C     NORMALISATION A 1 DU VECTEUR PROPRE
      DNORM = SQRT( XV1 ** 2 + YV1 ** 2 )
      XV1   = XV1 / DNORM
      YV1   = YV1 / DNORM
      RETURN
      END
