      SUBROUTINE PN6DVA( N1, P, X, Y, Z, U, V, W,  VAP6 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : VALEUR DU POLYNOME PRODUIT P AU POINT X, Y, Z, U, V, W
C -----
C
C ENTREES:
C --------
C N1     : DEGRE+1 DU POLYNOME P OU NOMBRE DE SES COEFFICIENTS
C P      : TABLEAU A 6 INDICES CONTENANT LES COEFFICIENTS DU POLYNOME P
C          P(I,J,K,L,M,N)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C                                        U**(L-1) V**(M-1) W**(N-1)
C X,Y,Z,U,V,W : 6 COORDONNEES DU POINT OU P DOIT ETRE CALCULE
C
C SORTIES:
C --------
C VAP6   : LA VALEUR DE P AU POINT X, Y, Z, U, V, W
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University           Juillet 2005
C ......................................................................
      DOUBLE PRECISION P(N1,N1,N1,N1,N1,N1),
     %                 VAP6, X,Y,Z,U,V,W,
     %                 XP, YP, ZP, UP, VP, WP,
     %                 VWP, UVWP, ZUVWP, YZUVWP
C
      VAP6 = 0.D0
      WP   = 1.D0
      DO 1 N=1,N1
         VP = 1.D0
         DO 2 M=1,N1
            VWP = VP * WP
            UP  = 1.D0
            DO 3 L=1,N1
               UVWP = UP * VWP
               ZP = 1.D0
               DO 4 K=1,N1
                  ZUVWP = ZP * UVWP
                  YP = 1.D0
                  DO 5 J=1,N1
                     YZUVWP = YP * ZUVWP
                     XP  = 1.D0
                     DO 6 I=1,N1
                        VAP6 = VAP6 + P(I,J,K,L,M,N) * XP * YZUVWP
                        XP  = XP * X
 6                   CONTINUE
                     YP = YP * Y
 5                CONTINUE
                  ZP = ZP * Z
 4             CONTINUE
               UP = UP * U
 3          CONTINUE
            VP = VP * V
 2       CONTINUE
         WP = WP * W
 1    CONTINUE
C
      IF (DABS(VAP6) .LT. 1.D-14) VAP6 = 0.D0
      RETURN
      END
