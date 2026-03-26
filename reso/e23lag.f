      SUBROUTINE E23LAG( NBPOLY, NBPOF, NOPOF, POLYF, DPOLYF, X,
     &                   GL, DGL, DELTA)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES COORDONNEES D'UN POINT D'UNE FACE
C -----    ET DU JACOBIEN DE L APPLICATION G : FACE -> FACE DE R ** 3
C
C ENTREES:
C --------
C NBPOLY : NOMBRE DE POINTS DE L'EF TRIDIMENSIONNEL COMPLET
C NBPOF  : NOMBRE DE POINTS DEFINISSANT L APPLICATION G DE LA FACE
C NOPOF  : NO DANS L ELEMENT DES NBPOF POINTS DE LA FACE
C POLYF  : POLYF(I)=VALEUR DU I-EME POLYNOME EN CE POINT
C DPOLYF : DPOLYF(I,J)=VALEUR DE LA I-EME DERIVEE DU J-EME
C                      POLYNOME EN CE POINT
C X      : COORDONNEES DES POINTS DE L ELEMENT  (NBPOLY,3)
C
C SORTIES:
C --------
C GL     : GL(3) LES 3 COORDONNEES DU POINT SUR LA FACE COURANTE
C DGL    : DGL(2,3) LES 2 DERIVEES DE G EN CE POINT DE LA FACE ACTUELLE
C          DGL(I,J) = D(GL)(J) / DXI
C DELTA  : JACOBIEN DE L APPLICATION G EN CE POINT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
      REAL             X(NBPOLY,3)
      INTEGER          NOPOF(NBPOF)
      DOUBLE PRECISION POLYF(NBPOF),
     &                 DPOLYF(2,NBPOF),
     &                 GL(3),
     &                 DGL(2,3),
     &                 DELTA, PPP, DP1, DP2, XX
C
C     MISE A ZERO DES TABLEAUX GL ET DGL
C     ---------------------------------
      DO 10 I=1,3
         GL ( I ) = 0D0
         DGL(1,I) = 0D0
         DGL(2,I) = 0D0
   10 CONTINUE
C
C     CALCUL DE G : R **2 -> R ** 3 ET DE SES DERIVEES AUX POINTS
C     D INTEGRATION DE LA FACE
C     ----------------------------------------------------------
      DO 30 J=1,NBPOF
         K   = NOPOF(J)
         PPP = POLYF(J)
         DP1 = DPOLYF(1,J)
         DP2 = DPOLYF(2,J)
         DO 20 I=1,3
            XX    = X(K,I)
            GL(I) = GL(I) + PPP * XX
            DGL(1,I) = DGL(1,I) + DP1 * XX
            DGL(2,I) = DGL(2,I) + DP2 * XX
   20    CONTINUE
   30 CONTINUE
C
C     CALCUL DU JACOBIEN DELTA DE G
C     -----------------------------
      CALL JAR2R3( DGL , DELTA )
      END
