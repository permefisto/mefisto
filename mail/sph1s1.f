      SUBROUTINE SPH1S1( N,      RAYOSP,  X,     Y,    Z,
     %                   T,      T0,     TBG,    TBD,
     %                   NBSOMT, NBTRIA, XYZSOM, NOSOMT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT    : GENERER LA TRIANGULATION REGULIERE DE LA SURFACE D'UNE SPHERE
C --------
C ENTREES:
C --------
C N      : NOMBRE D'ARETES SUR CHAQUE COTE DE TRIANGLE
C RAYOSP : RAYON DE LA SPHERE
C X      : COORDONNEE X DU CENTRE DE LA SPHERE
C Y      : COORDONNEE Y DU CENTRE DE LA SPHERE
C Z      : COORDONNEE Z DU CENTRE DE LA SPHERE

C TABLEAUX AUXILIAIRES:
C ---------------------
C T,T0   : TABLEAUX( (N+1) * (2*N+1) + 1 ) D'ENTIERS
C TBG,TFD: TABLEAUX( 3 * N + 2 )           D'ENTIERS

C SORTIES:
C --------
C NBSOMT : NOMBRE DE SOMMETS   GENERES
C NBTRIA : NOMBRE DE TRIANGLES GENERES ( = 20 * N**2 )
C XYZSOM : TABLEAU DES COORDONNEES DES SOMMETS
C NOSOMT : TABLEAU DU NUMERO DES 3 SOMMETS DES TRIANGLES (0 EN POSITION 4)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS : A. CUBIER ET F. DEVUYST DEA A.N. UPMC PARIS    DECEMBRE 1989
C MODIFS  : A. PERRONNET  ANALYSE NUMERIQUE  UPMC PARIS        JUIN 1996
C234567--------------------------------------------------------------012
      REAL              XYZSOM(3,NBSOMT)
      INTEGER           NOSOMT(4,NBTRIA)
      INTEGER           N, A, B, C, D, E, F
      INTEGER           T(0:*), T0(0:*), TBG(0:*), TBD(0:*)
      REAL              X, Y, Z, RAYOSP
      DOUBLE PRECISION  PI5, PHI, CPHI, SPHI, CPS5, SPS5, CDPS5, SDPS5

      NBNOLO = (N+1)*(2*N+1)

C     PHI    = PI * 26.5650511771D0 / 180D0
      PHI    = ATAN( 0.5D0 )
      CPHI   = COS (PHI)
      SPHI   = SIN (PHI)

      PI5    = ATAN( 1D0 ) * 4D0 / 5D0
      CPS5   = COS( PI5 )
      SPS5   = SIN( PI5 )
      CDPS5  = COS( 2D0 * PI5 )
      SDPS5  = SIN( 2D0 * PI5 )

C     NUMERO DES POINTS DE L'ICOSAEDRE
C     ================================
      A = 1
      C = (N+1)*(N+2)/2
      B = C - N
      E = B + (N+1)*(N+2) - N - 2
      D = E - N
      F = D + C - 1

C     COORDONNEES DES POINTS DE L'ICOSAEDRE
C     =====================================
      XYZSOM(1,A) = 0.
      XYZSOM(2,A) = 0.
      XYZSOM(3,A) = 1.
      XYZSOM(1,B) = REAL(  CPS5 * CPHI )
      XYZSOM(2,B) = REAL(- SPS5 * CPHI )
      XYZSOM(3,B) = REAL(  SPHI )
      XYZSOM(1,C) = REAL(  CPS5 * CPHI )
      XYZSOM(2,C) = REAL(  SPS5 * CPHI )
      XYZSOM(3,C) = REAL(  SPHI )
      XYZSOM(1,D) = REAL(  CDPS5 * CPHI )
      XYZSOM(2,D) = REAL(- SDPS5 * CPHI )
      XYZSOM(3,D) = REAL(- SPHI )
      XYZSOM(1,E) = REAL(  CPHI )
      XYZSOM(2,E) = 0.
      XYZSOM(3,E) = REAL(- SPHI )
      XYZSOM(1,F) = 0.
      XYZSOM(2,F) = 0.
      XYZSOM(3,F) = -1.

C     CONSTRUCTION DES TABLEAUX TBG ET TBD
C     ====================================
      TBD(0) = 0
      DO I=1,N
         TBD(I) = TBD(I-1)+I
      ENDDO
      DO I=1,N+1
         TBD(N+I) = TBD(N+I-1)+N+1
      ENDDO
      DO I=1,N
         J=N-I+1
         TBD(2*N+1+I) = TBD(2*N+I) + J
      ENDDO

      TBG(0) = 1
      DO I=1,N+1
         TBG(I) = TBG(I-1)+I-1
      ENDDO
      DO I=1,N
         TBG(N+1+I) = TBG(N+I)+N+1
      ENDDO
      DO I=1,N
         J=N-I+2
         TBG(2*N+1+I) = TBG(2*N+I) +J
      ENDDO

      IC = NBNOLO + 1
      NBTRIA = 0

C     BOUCLE SUR LES BLOCS DE L'ICOSAEDRE
C     ===================================
      DO 500 M=1,5
      IF (M.EQ.1) THEN
         DO I=1, NBNOLO
            T(I) = I
         ENDDO
      ELSE
        ICTBLP = 1
        ICTBLO = 1
        IC0 = IC
        DO I=1,NBNOLO
           T0(I) = T(I)
        ENDDO
        DO I=1,NBNOLO
           IF (I.EQ.TBG(ICTBLO)) THEN
              T(I) = T(TBD(ICTBLO))
              ICTBLO = ICTBLO + 1
           ELSE IF ( I.EQ.TBD(ICTBLP) .AND. M.EQ.5 ) THEN
              T(I) = TBG(ICTBLP)
              ICTBLP = ICTBLP + 1
           ELSE
              T(I) = IC
              IC = IC + 1
           ENDIF
           IF (ICTBLO.EQ.2) ICTBLP=2
           ENDDO
      ENDIF

C     CALCUL DES COORDONNEES DES NOSOMT
C     =================================
      IF (M .EQ. 1 ) THEN

      DO I=1,N
         NUM1 = I*(I+1)/2+1
         NUM2 = (I+1)*(I+2)/2
         DO K=1,3
            RISN = FLOAT(I)/FLOAT(N)
            XYZSOM (K,NUM1) = (1.-RISN)*XYZSOM(K,A)+RISN*XYZSOM(K,B)
            XYZSOM (K,NUM2) = (1.-RISN)*XYZSOM(K,A)+RISN*XYZSOM(K,C)
         ENDDO
         DO J=1,I-1
            RJSI = FLOAT(J)/FLOAT(I)
            DO K=1,3
               XYZSOM (K,NUM1+J) = (1.-RJSI)*XYZSOM(K,NUM1)
     %                            +RJSI*XYZSOM(K,NUM2)
            ENDDO
         ENDDO
      ENDDO

      DO I=1,N
         NUM1 = B + I*(N+1)
         NUM2 = C + I*(N+1)
         NUM3 = NUM1 + I
         DO K=1,3
            RISN = FLOAT(I)/FLOAT(N)
            XYZSOM (K,NUM1) = (1.-RISN)*XYZSOM(K,B)+RISN*XYZSOM(K,D)
            XYZSOM (K,NUM2) = (1.-RISN)*XYZSOM(K,C)+RISN*XYZSOM(K,E)
            XYZSOM (K,NUM3) = (1.-RISN)*XYZSOM(K,B)+RISN*XYZSOM(K,E)
         ENDDO
         DO J=1,I-1
            RJSI = FLOAT(J)/FLOAT(I)
            DO K=1,3
               XYZSOM (K,NUM1+J) = (1.-RJSI)*XYZSOM(K,NUM1)
     %                            +RJSI*XYZSOM(K,NUM3)
            ENDDO
         ENDDO
         DO J=N-I-1,1,-1
            RJSI = FLOAT(J)/FLOAT(N-I)
            DO K=1,3
               XYZSOM (K,NUM3+J) = (1.-RJSI)*XYZSOM(K,NUM3)
     %                            + RJSI*XYZSOM(K,NUM2)
            ENDDO
         ENDDO

      ENDDO

      NUM1 = D + N + 1
      NUM2 = E + N
      DO I=N-1,1,-1
         DO K=1,3
            RISN = FLOAT(I)/FLOAT(N)
            XYZSOM (K,NUM1) = (1.-RISN)*XYZSOM(K,F)+RISN*XYZSOM(K,D)
            XYZSOM (K,NUM2) = (1.-RISN)*XYZSOM(K,F)+RISN*XYZSOM(K,E)
         ENDDO
         DO J=1,I-1
            RJSI = FLOAT(J)/FLOAT(I)
            DO K=1,3
               XYZSOM (K,NUM1+J) = (1.-RJSI)*XYZSOM(K,NUM1)
     %                            +RJSI*XYZSOM(K,NUM2)
            ENDDO
         ENDDO
         NUM1 = NUM1 + I + 1
         NUM2 = NUM2 + I
      ENDDO

C     PROJECTION SUR LA SPHERE DE RAYON R
C     ===================================
      DO I=1,F
         RNORME = SQRT( XYZSOM(1,I)**2+XYZSOM(2,I)**2+XYZSOM(3,I)**2 )
         XYZSOM(1,I) = XYZSOM(1,I) * RAYOSP / RNORME
         XYZSOM(2,I) = XYZSOM(2,I) * RAYOSP / RNORME
         XYZSOM(3,I) = XYZSOM(3,I) * RAYOSP / RNORME
      ENDDO

C FIN DU CAS M=1

      ELSE
C
      DO I=1,NBNOLO
         X1 = XYZSOM(1,T0(I))
         Y1 = XYZSOM(2,T0(I))
         Z1 = XYZSOM(3,T0(I))
         IF (T(I).GE.IC0) THEN
            XYZSOM(1,T(I)) = REAL( X1*CDPS5 - Y1*SDPS5 )
            XYZSOM(2,T(I)) = REAL( X1*SDPS5 + Y1*CDPS5 )
            XYZSOM(3,T(I)) = Z1
         ENDIF
      ENDDO

      ENDIF

C     CREATION DU NUMERO DES SOMMETS DES TRIANGLES DU MAILLAGE
C     ========================================================
      DO I=0,N-1
         NUM1 = I*(I+1)/2 +1
         NUM2 = (I+1)*(I+2)/2+1
         DO J=0,I
            NBTRIA = NBTRIA + 1
            NOSOMT(1,NBTRIA) = T(NUM1+J)
            NOSOMT(2,NBTRIA) = T(NUM2+J)
            NOSOMT(3,NBTRIA) = T(NUM2+J+1)
            NOSOMT(4,NBTRIA) = 0
            IF ( J.NE.I ) THEN
               NBTRIA = NBTRIA + 1
               NOSOMT(1,NBTRIA) = T(NUM2+J+1)
               NOSOMT(2,NBTRIA) = T(NUM1+J+1)
               NOSOMT(3,NBTRIA) = T(NUM1+J)
               NOSOMT(4,NBTRIA) = 0
            ENDIF
         ENDDO
      ENDDO
C
      DO I=0,N-1
         NUM1 = B + I*(N+1)
         NUM2 = B + (I+1)*(N+1)
         DO J=0,N-1
            NBTRIA = NBTRIA + 1
            NOSOMT(1,NBTRIA) = T(NUM1+J)
            NOSOMT(2,NBTRIA) = T(NUM2+J)
            NOSOMT(3,NBTRIA) = T(NUM2+J+1)
            NOSOMT(4,NBTRIA) = 0
            NBTRIA = NBTRIA + 1
            NOSOMT(1,NBTRIA) = T(NUM2+J+1)
            NOSOMT(2,NBTRIA) = T(NUM1+J+1)
            NOSOMT(3,NBTRIA) = T(NUM1+J)
            NOSOMT(4,NBTRIA) = 0
         ENDDO
      ENDDO

      INCRE1 = 0
      INCRE2 = N+1
      INC2   = N+1
      DO I=N,1,-1
         NUM1 = D + INCRE1
         NUM2 = D + INCRE2
         DO J=0,I-1
            NBTRIA = NBTRIA + 1
            NOSOMT(1,NBTRIA) = T(NUM1+J)
            NOSOMT(2,NBTRIA) = T(NUM2+J)
            NOSOMT(3,NBTRIA) = T(NUM1+J+1)
            NOSOMT(4,NBTRIA) = 0
            IF ( J.NE.(I-1) ) THEN
               NBTRIA = NBTRIA + 1
               NOSOMT(1,NBTRIA) = T(NUM1+J+1)
               NOSOMT(2,NBTRIA) = T(NUM2+J)
               NOSOMT(3,NBTRIA) = T(NUM2+J+1)
               NOSOMT(4,NBTRIA) = 0
            ENDIF
         ENDDO
         INCRE1 = INCRE2
         INC2   = INC2 - 1
         INCRE2 = INCRE2 + INC2
      ENDDO

C     FIN DE LA BOUCLE SUR LES BLOCS
C     ==============================
 500  ENDDO

C     TRANSLATION DE LA SPHERE AU POINT (X,Y,Z)
C     =========================================
      NBSOMT = IC - 1
      DO I=1,NBSOMT
         XYZSOM(1,I) = XYZSOM (1,I) + X
         XYZSOM(2,I) = XYZSOM (2,I) + Y
         XYZSOM(3,I) = XYZSOM (3,I) + Z
      ENDDO

      RETURN
      END
