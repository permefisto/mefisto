      SUBROUTINE CONTOU( XE , YE , N )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ORIENTER LE CONTOUR DE L'ENVELOPPE CONVEXE DE N POINTS
C -----
C       ** ATTENTION : DANS CETTE VERSION N EST LIMITE A 10 ! ***
C
C ENTREES:
C --------
C XE, YE : COORDONNEES DES N POINTS EN ENTREE
C
C SORTIES:
C --------
C XE, YE : COORDONNEES DES N POINTS EN SORTIE
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C XS, YS, COXY
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY  ANALYSE NUMERIQUE PARIS UPMC       JUIN 1994
C....................................................................012
      PARAMETER ( NMAX = 10 )
      COMMON / EPSSSS / EPZERO,EPSXYZ
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DIMENSION         XE(N),YE(N),XS(NMAX),YS(NMAX),COXY(NMAX)
C
C     TEST DE CONFORMITE
C     ==================
      IF ( N.GT.NMAX ) THEN
         WRITE(IMPRIM,1000) N,NMAX
1000     FORMAT(//,'SP CONTOU : NOMBRE DE POINTS',I7,
     &   ' PLUS GRAND QUE ',I7,' =>  PAS D''ORIENTATION DU CONTOUR',//)
         RETURN
      ENDIF
C
C     COPIE DES COORDONNEES
C     =====================
      DO 50 I=1,N
         XS(I)=XE(I)
         YS(I)=YE(I)
 50   CONTINUE
C
C     RECHERCHE DU POINT N1 DANS XE, YE D'ORDONNEE MINIMALE
C     =====================================================
      Y1 = YE(1)
      N1 = 1
      DO 10 I=2,N
         Y = YE(I)
         IF( Y .GT. Y1 + EPSXYZ ) GOTO 10
         IF( Y .LT. Y1 - EPSXYZ ) THEN
            Y1 = Y
            N1 = I
         ELSE
C           A ORDONNEES EGALES,
C           ON CHOISIT L'ABSCISSE LA PLUS PETITE
            IF( XE(I) .LT. XE(N1) ) THEN
               N1 = I
            ENDIF
         ENDIF
 10   CONTINUE
C                                    --> ->
C     CALCUL DE -COSINUS DES ANGLES (N1M,OX)
C     ======================================
      X1 = XE(N1)
      Y1 = YE(N1)
      DO 20 I=1,N
         IF( I .NE. N1 ) THEN
            X = X1 - XE(I)
            Y = YE(I) - Y1
            DENO = X * X + Y * Y
            IF ( DENO .LT. EPZERO ) THEN
               COXY(I) = 0
            ELSE
               COXY(I) = X / SQRT( DENO )
               IF( ABS(COXY(I)) .LT. EPZERO ) THEN
                  COXY(I) = 0
               ENDIF
            ENDIF
         ELSE
            COXY(I) = -1.01
         ENDIF
 20   CONTINUE
C
C     TRI DES COSINUS SELON LEUR VALEUR CROISSANTE
C     ============================================
      DO 30 K=N,2,-1
         INDIC=0
         DO 31 I=1,K-1
            IP1=I+1
            IF(COXY(I).GT.COXY(IP1)+EPSXYZ) THEN
C              ECHANGE
               INDIC=1
               COXYI=COXY(I)
               COXY(I)=COXY(IP1)
               COXY(IP1)=COXYI
               XSI=XS(I)
               XS(I)=XS(IP1)
               XS(IP1)=XSI
               YSI=YS(I)
               YS(I)=YS(IP1)
               YS(IP1)=YSI
            END IF
 31       CONTINUE
         IF(INDIC.EQ.0) GOTO 40
 30   CONTINUE
C
C     REAGENCEMENT DES POINTS EN
C     TOURNANT DANS LE SENS DIRECT
C     ============================
 40   DO 60 I=1,N
         XE(I)=XS(I)
         YE(I)=YS(I)
 60   CONTINUE
C
      RETURN
      END
