      SUBROUTINE MAGCVL( NCODSA, NTDL, LPLIGN, LPCOLO, AG, ALPHA, X, P )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU PRODUIT: P = P + ALPHA * A * U
C -----    ALPHA REEL DOUBLE, A MATRICE MORSE, X P VECTEURS GLOBAUX
C
C ENTREES:
C --------
C NCODSA : CODE DE STOCKAGE
C NTDL   : ORDRE DE LA MATRICE AG
C LPLIGN : LES POINTEURS SUR LES LIGNES DE LA MATRICE MORSE AG
C LPCOLO : LES NUMEROS DES COLONNES DES COEFFICIENTS DE LA MATRICE MORSE
C AG     : LES COEFFICIENTS DE LA MATRICE MORSE AG
C ALPHA  : REEL DOUBLE PRECISION
C X      : LE VECTEUR A MULTIPLIER
C
C SORTIES:
C --------
C P      : LE VECTEUR FINAL P + ALPHA * A * U
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY  ANALYSE NUMERIQUE UPMC PARIS        JANVIER 1999
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  AG(1:*), X(NTDL), P(NTDL), S
      INTEGER           LPLIGN(NTDL+1), LPCOLO(1:*)
C
      IF ( NCODSA .GT. 0 ) THEN
C
C        ------   MATRICE MORSE SYMETRIQUE   ------
         J1 = 1
         DO 2 I=1,NTDL
            J2=LPLIGN(I+1) - 1
            IF(J1.GT.J2) GOTO 4
            DO 3 J=J1,J2
               NOCOL = LPCOLO(J)
               P(I)     = P(I)     + ALPHA*X(NOCOL) * AG(J)
               P(NOCOL) = P(NOCOL) + ALPHA*X(I)     * AG(J)
 3          CONTINUE
C           ---   LE COEFFICIENT DIAGONAL   ---
 4          J2 = J2 + 1
            P(I) = P(I) + ALPHA*X(I) * AG(J2)
            J1 = J2 + 1
 2       CONTINUE
C
      ELSE
C
C        ------   MATRICE MORSE NON SYMETRIQUE   ------
         J1 = 1
         DO 101 I=1,NTDL
            J2 = LPLIGN(I+1)
            S = 0.D0
            DO 102 J=J1,J2
               NOCOL = LPCOLO(J)
               S     = S + ALPHA*X(NOCOL) * AG(J)
102         CONTINUE
            P(I) = S
            J1 = J2 + 1
101      CONTINUE
      ENDIF
C
      RETURN
      END
