      SUBROUTINE JACOBI(A,N,NP,D,V,NROT)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DES VALEURS ET VECTEURS PROPRES D-UNE MATRICE NxN
C ----- SYMETRIQUE EN DOUBLE PRECISION PAR LA METHODE DE JACOBI
C       NUMERICAL RECIPES PAGE 356-348 Press Flannery Teukosky Vetterling
C
C ENTREES:
C --------
C A      : LA MATRICE SYMETRIQUE STOCKEE NPxNP
C N      : NOMBRE DE LIGNES ET COLONNES ET DE VECTEURS PROPRES A CALCULER
C NP     : DECLARATION A(NP,NP) ET DES AUTRES TABLEAUX
C
C SORTIES:
C --------
C D      : D(N) LES N VALEURS PROPRES
C V      : V(N,K) LE K-EME VECTEUR PROPRE ORTHONORMALISE
C NROT   : NOMBRE DE ROTATIONS DE LA METHODE DE JACOBI
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS:  NUMERICAL RECIPES PAGE 356-348 Press Flannery Teukosky Vetterling
C23456...............................................................12
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NMAX=100)
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)
C
      DO 12 IP=1,N
        DO 11 IQ=1,N
          V(IP,IQ)=0D0
11      CONTINUE
        V(IP,IP)=1.D0
12    CONTINUE
      DO 13 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.D0
13    CONTINUE
      NROT=0
      DO 24 I=1,50
        SM=0.D0
        DO 15 IP=1,N-1
          DO 14 IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
14        CONTINUE
15      CONTINUE
        IF(SM.EQ.0.D0)RETURN
        IF(I.LT.4)THEN
          TRESH=0.2D0*SM/N**2
        ELSE
          TRESH=0.D0
        ENDIF
        DO 22 IP=1,N-1
          DO 21 IQ=IP+1,N
            G=100.D0*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP)))
     *         .AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ)=0.
            ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(ABS(H)+G.EQ.ABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5D0*H/A(IP,IQ)
                T=1.D0/(ABS(THETA)+SQRT(1.D0+THETA**2))
                IF(THETA.LT.0.)T=-T
              ENDIF
              C=1.D0/SQRT(1+T**2)
              S=T*C
              TAU=S/(1.D0+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=0.D0
              DO 16 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
16            CONTINUE
              DO 17 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
17            CONTINUE
              DO 18 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
18            CONTINUE
              DO 19 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
19            CONTINUE
              NROT=NROT+1
            ENDIF
21        CONTINUE
22      CONTINUE
        DO 23 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.D0
23      CONTINUE
24    CONTINUE
      WRITE(IMPRIM,*) 'SP JACOBI:50 iterations should never happen'
      RETURN
      END
