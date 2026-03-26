      SUBROUTINE E22LAG( NBPOLY, NBPOLA, NOPOAR, POLYA, DPOLYA, X,
     &                   GL, DGL, DELTA )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DES COORDONNEES D UN POINT D UNE ARETE
C -----   ET DU JACOBIEN DE L APPLICATION G : SEGMENT -> ARETE DE R ** 2
C
C ENTREE :
C --------
C NBPOLY : NBRE DE POINTS DE L ELEMENT BIDIMENSIONNEL COMPLET
C NBPOLA : NBRE DE POINTS DEFINISSANT L APPLICATION G DE L ARETE
C NOPOAR : NO DANS L ELEMENT DES NBPOLA POINTS DE L ARETE
C POLYA  : POLYA(I)=VALEUR DU I-EME POLYNOME EN CE POINT
C DPOLYA : DPOLYA(I,J)=VALEUR DE LA I-EME DERIVEE DU J-EME
C                      POLYNOME EN CE POINT
C X      : COORDONNEES DES POINTS DE L ELEMENT  (NBPOLY,2)
C
C SORTIE :
C --------
C GL     : LES 2 COORDONNEES DU POINT SUR L ARETE COURANTE
C DGL    : LA DERIVEE DE G EN CE POINT DE L ARETE ACTUELLE
C          DGL(J) = D(GL)(J)  J=1 OU 2
C DELTA  : JACOBIEN DE L APPLICATION G EN CE POINT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1995
C23456---------------------------------------------------------------012
      REAL             X(NBPOLY,2)
      INTEGER          NOPOAR(1)
      DOUBLE PRECISION POLYA(NBPOLA),DPOLYA(NBPOLA),GL(2),DGL(2),
     &                 DELTA,PPP,DPP,XX,DSQRT
C
C     MISE A ZERO DES TABLEAUX GL ET DGL
C     ----------------------------------
      DO 10 I=1,2
         GL ( I ) = 0D0
         DGL( I ) = 0D0
   10 CONTINUE
C
C     CALCUL DE G : (0,1) -> R ** 2 ET DE SA DERIVEE AUX POINTS
C     D INTEGRATION DE L ARETE
C     ----------------------------------------------------------
      DO 30 J=1,NBPOLA
         K   = NOPOAR(J)
         PPP = POLYA(J)
         DPP = DPOLYA(J)
         DO 20 I=1,2
            XX     = X(K,I)
            GL (I) = GL (I) + PPP * XX
            DGL(I) = DGL(I) + DPP * XX
   20    CONTINUE
   30 CONTINUE
C
C     CALCUL DU JACOBIEN DELTA DE G
C     -----------------------------
      DELTA = DSQRT ( DGL(1) ** 2 + DGL(2) ** 2 )

      RETURN
      END
