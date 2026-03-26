      SUBROUTINE SPH1S8( N,      RAYOSP, X, Y, Z,
     %                   XYZSOM, NOEUDS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT    : GENERER LE MAILLAGE DE LA SURFACE D'UN 1/8EME DE SPHERE
C --------
C ENTREES:
C --------
C N      : NOMBRE D'ARETES SUR CHAQUE COTE DE TRIANGLE
C RAYOSP : RAYON DE LA SPHERE
C X      : COORDONNEE X DU CENTRE DE LA SPHERE
C Y      : COORDONNEE Y DU CENTRE DE LA SPHERE
C Z      : COORDONNEE Z DU CENTRE DE LA SPHERE
C
C SORTIES:
C --------
C XYZSOM : TABLEAU DES COORDONNES DES NOEUDS
C NOEUDS : TABLEAU D'AFFECTATION DES SOMMETS AUX ELEMENTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS : A. CUBIER ET F. DEVUYST DEA A.N. UPMC PARIS    DECEMBRE 1989
C MODIFS  : A. PERRONNET  ANALYSE NUMERIQUE UPMC  PARIS        JUIN 1996
C234567--------------------------------------------------------------012
      DIMENSION XYZSOM(3, *)
      DIMENSION NOEUDS(4, *)
      INTEGER N, A, B, C
      REAL X, Y, Z, RAYOSP
C
      A = 1
      B = N*(N+1)/2 +1
      C = (N+1)*(N+2)/2
C
      XYZSOM(1,A) = 0.
      XYZSOM(2,A) = 0.
      XYZSOM(3,A) = 1.
      XYZSOM(1,B) = 1.
      XYZSOM(2,B) = 0.
      XYZSOM(3,B) = 0.
      XYZSOM(1,C) = 0.
      XYZSOM(2,C) = 1.
      XYZSOM(3,C) = 0.
C
      DO 10 I=1,N
         NUM1 = I*(I+1)/2+1
         NUM2 = (I+1)*(I+2)/2
         RISN = FLOAT(I)/FLOAT(N)
         DO 20 K=1,3
         XYZSOM(K,NUM1) = ( 1.-RISN )*XYZSOM(K,A)
     %                    + RISN*XYZSOM(K,B)
  20     CONTINUE
         DO 40 K=1,3
         XYZSOM(K,NUM2) = ( 1.-RISN )*XYZSOM(K,A)
     %                    + RISN*XYZSOM(K,C)
  40     CONTINUE
         DO 60 J=1,I-1
         RJSI = FLOAT(J)/FLOAT(I)
         DO 70 K=1,3
         XYZSOM(K,NUM1+J) = ( 1.-RJSI )*XYZSOM(K,NUM1)
     %                      + RJSI*XYZSOM(K,NUM2)
  70     CONTINUE
  60     CONTINUE
  10  CONTINUE
C
C     PROJECTION DU NOEUDS SUR LA SPHERE
C     ==================================
      DO 50 I=A,C
         RNORME = SQRT( XYZSOM(1,I)**2+XYZSOM(2,I)**2+XYZSOM(3,I)**2)
         DO 80 K=1,3
            XYZSOM(K,I) = RAYOSP*XYZSOM(K,I)/RNORME
  80     CONTINUE
  50  CONTINUE
C
C     TRANSLATION APPLIQUEE AUX NOEUDS
C     ================================
      DO 90 I=A,C
         XYZSOM(1,I) = XYZSOM(1,I) + X
         XYZSOM(2,I) = XYZSOM(2,I) + Y
         XYZSOM(3,I) = XYZSOM(3,I) + Z
  90  CONTINUE
C
C     CREATION ET NUMEROTATION DES ELEMENTS
C     =====================================
      NBELEM = 0
      DO 110 I=0,N-1
         NUM1=I*(I+1)/2+1
         NUM2=(I+1)*(I+2)/2+1
         DO 120 J=0,I
            NBELEM = NBELEM + 1
            NOEUDS(1,NBELEM) = NUM1 + J
            NOEUDS(2,NBELEM) = NUM2 + J
            NOEUDS(3,NBELEM) = NUM2 + J + 1
            NOEUDS(4,NBELEM) = 0
            IF ( J .NE. I ) THEN
               NBELEM = NBELEM + 1
               NOEUDS(1,NBELEM) = NUM2 + J + 1
               NOEUDS(2,NBELEM) = NUM1 + J + 1
               NOEUDS(3,NBELEM) = NUM1 + J
               NOEUDS(4,NBELEM) = 0
            ENDIF
 120     CONTINUE
 110  CONTINUE

      RETURN
      END
