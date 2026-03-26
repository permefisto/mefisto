      SUBROUTINE MAGCVE ( NUINIT, CTE,    NTDL,
     %                    NCODSA, LPLIGN, LPCOLO, AG,  X,
     %                    Y )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:    Y=CTE*AG*X ou Y=Y+CTE*AG*X avec
C ----    Y,X VECTEURS (NTDL)  AG MATRICE MORSE     SYMETRIQUE
C                              OU MATRICE DIAGONALE
C                              OU MATRICE MORSE NON SYMETRIQUE
C
C ENTREES:
C --------
C NUINIT : 0 => Y=  CTE*AG*X
C          1 => Y=Y+CTE*AG*X
C CTE    : LA CONSTANTE DE MULTIPLICATION
C NTDL   : NOMBRE DE LIGNES DE LA MATRICE ET DU VECTEUR
C NCODSA : 0 SI LA MATRICE A EST DIAGONALE
C          1 SI LA MATRICE EST SYMETRIQUE PROFIL
C         -1 SI LA MATRICE EST NON-SYMETRIQUE PROFIL
C LPLIGN : LES POINTEURS SUR LES LIGNES DE LA MATRICE MORSE AG
C LPCOLO : LES NUMEROS DES COLONNES DES COEFFICIENTS DE LA MATRICE MORSE
C AG     : LES COEFFICIENTS DE LA MATRICE MORSE AG
C X      : VECTEUR DE NTDL COMPOSANTES
C
C SORTIE ou MODIFIE:
C ------------------
C Y      : VECTEUR DE NTDL COMPOSANTES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1999
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  CTE, AG(1:*), X(NTDL), Y(NTDL), S, AUX
      INTEGER           LPLIGN(0:NTDL), LPCOLO(1:*)
C
      IF( CTE .EQ. 1D0 ) THEN
C
         IF( NCODSA .EQ. 0 ) THEN
C
C           A  DIAGONALE
C           ============
            IF( NUINIT .EQ. 0 ) THEN
               DO 2 I=1,NTDL
                  Y(I) = AG(I) * X(I)
 2             CONTINUE
            ELSE
               DO 4 I=1,NTDL
                  Y(I) = Y(I) + AG(I) * X(I)
 4             CONTINUE
            ENDIF
            RETURN
         ENDIF
C
C        INITIALISATION A 0 OU NON DU TABLEAU P
         IF( NUINIT .EQ. 0 ) THEN
            DO 5 I=1,NTDL
               Y(I) = 0.D0
 5          CONTINUE
         ENDIF
C
         IF ( NCODSA .GT. 0 ) THEN
C
C           MATRICE MORSE SYMETRIQUE
C           ========================
            J1 = 1
            DO 12 I=1,NTDL
               J2 = LPLIGN(I) - 1
               DO 10 J=J1,J2
                  NOCOL    = LPCOLO(J)
                  Y(I)     = Y(I)     + X(NOCOL) * AG(J)
                  Y(NOCOL) = Y(NOCOL) + X(I)     * AG(J)
 10            CONTINUE
C              LE COEFFICIENT DIAGONAL
               J2   = J2 + 1
               Y(I) = Y(I) + X(I) * AG(J2)
               J1   = J2 + 1
 12         CONTINUE
C
         ELSE
C
C           MATRICE MORSE NON SYMETRIQUE
C           ============================
            J1 = 1
            DO 21 I=1,NTDL
               J2 = LPLIGN(I)
               S = 0.D0
               DO 22 J=J1,J2
                  S = S + X(LPCOLO(J)) * AG(J)
 22            CONTINUE
               Y(I) = S
               J1   = J2 + 1
 21         CONTINUE
         ENDIF
C
      ELSE
C
C        LA CTE N'EST PAS 1D0 => PRODUIT
C        ===============================
         IF( NCODSA .EQ. 0 ) THEN
C
C           A  DIAGONALE
C           ============
            IF( NUINIT .EQ. 0 ) THEN
               DO 32 I=1,NTDL
                  Y(I) = CTE * AG(I) * X(I)
 32            CONTINUE
            ELSE
               DO 34 I=1,NTDL
                  Y(I) = Y(I) + CTE * AG(I) * X(I)
 34            CONTINUE
            ENDIF
            RETURN
         ENDIF
C
C        INITIALISATION A 0 OU NON DU TABLEAU P
         IF( NUINIT .EQ. 0 ) THEN
            DO 35 I=1,NTDL
               Y(I) = 0.D0
 35         CONTINUE
         ENDIF
C
         IF ( NCODSA .GT. 0 ) THEN
C
C           MATRICE MORSE SYMETRIQUE
C           ========================
            J1 = 1
            DO 42 I=1,NTDL
               AUX = CTE * X(I)
               J2  = LPLIGN(I) - 1
               DO 40 J=J1,J2
                  NOCOL    = LPCOLO(J)
                  Y(I)     = Y(I)     + CTE * X(NOCOL) * AG(J)
                  Y(NOCOL) = Y(NOCOL) + AUX * AG(J)
 40            CONTINUE
C              LE COEFFICIENT DIAGONAL
               J2   = J2 + 1
               Y(I) = Y(I) + AUX * AG(J2)
               J1   = J2 + 1
 42         CONTINUE
C
         ELSE
C
C           MATRICE MORSE NON SYMETRIQUE
C           ============================
            J1 = 1
            DO 60 I=1,NTDL
               J2 = LPLIGN(I)
               S  = 0.D0
               DO 50 J=J1,J2
                  S = S + CTE * X(LPCOLO(J)) * AG(J)
 50            CONTINUE
               Y(I) = S
               J1   = J2 + 1
 60         CONTINUE
         ENDIF
      ENDIF
C
      RETURN
      END
