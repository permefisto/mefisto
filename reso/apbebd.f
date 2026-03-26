      SUBROUTINE APBEBD( NCODSA, MU, A, X, NDSM, NTDL, Y )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:    Y=A*X  Y,X TABLEAUX (NDSM,NTDL)  A MATRICE PROFIL SYMETRIQUE
C ----                                       OU DIAGONALE
C                                            OU NON SYMETRIQUE
C
C ATTENTION: SP DIFFERENT DE MAPRVE CAR X ET Y SONT LES TABLEAUX TRANSPOSES
C
C ENTREES:
C --------
C NCODSA : 0 SI LA MATRICE A EST DIAGONALE
C          1 SI LA MATRICE EST SYMETRIQUE PROFIL
C         -1 SI LA MATRICE EST NON-SYMETRIQUE PROFIL
C MU     : POINTEUR SUR CHAQUE COEFFICIENT DIAGONAL
C A      : MATRICE PROFIL D ORDRE NTDL ( CF TMS PROFIL )
C B      : TABLEAU (NDSM,NTDL) A MULTIPLIER PAR A
C
C SORTIES:
C --------
C Y      : TABLEAU (NDSM,NTDL)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    SEPTEMBRE 1977
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  A(*),X(NDSM,NTDL),Y(NDSM,NTDL)
      INTEGER           MU(1:NTDL+1)
C
C     AIGUILLAGE SELON NCODSA
C     =======================
      IF( NCODSA .EQ. 0 ) THEN
C
C        A  DIAGONALE
C        ------------
         DO 5 I=1,NTDL
            DO 2 NDS=1,NDSM
               Y(NDS,I) = A(I) * X(NDS,I)
   2        CONTINUE
   5     CONTINUE
         RETURN
      ENDIF
C
C     INITIALISATION A 0 DU TABLEAU Y
C     -------------------------------
      DO 11 J=1,NTDL
         DO 12 I=1,NDSM
            Y(I,J) = 0D0
   12    CONTINUE
   11 CONTINUE
C
      IF( NCODSA .GT. 0 ) THEN
C
C        A SYMETRIQUE . MU PAR D.L
C        -------------------------
         JA = 0
         DO 28 J=1,NTDL
            JMI = MU(J+1) - MU(J)
            I   = J - JMI
            DO 23 JD=1,JMI-1
               JA = JA + 1
               I  = I  + 1
               DO 22 NDS=1,NDSM
                  Y(NDS,J) = Y(NDS,J) + A(JA) * X(NDS,I)
                  Y(NDS,I) = Y(NDS,I) + A(JA) * X(NDS,J)
   22          CONTINUE
   23       CONTINUE
C
            JA = JA + 1
            DO 25 NDS=1,NDSM
               Y(NDS,J) = Y(NDS,J) + A(JA) * X(NDS,J)
   25       CONTINUE
   28    CONTINUE
C
      ELSE
C
C        A NON-SYMETRIQUE . MU PAR D.L
C        -----------------------------
         DO 31 I=1,NTDL
            MUI = MU(I)
            IH  = (MU(I+1)-MUI) / 2
            IF( IH .EQ. 0 ) GOTO 34
            MUJ = MUI + IH
            ILI = I   - IH
            DO 32 J=1,IH
               IMUI = MUI + J
               JMUJ = MUJ + J
               DO 33 NDS=1,NDSM
                  Y(NDS,I  )=Y(NDS,I  )+A(IMUI)*X(NDS,ILI)
                  Y(NDS,ILI)=Y(NDS,ILI)+A(JMUJ)*X(NDS,I  )
   33          CONTINUE
               ILI = ILI + 1
   32       CONTINUE
C
   34       DO 35 NDS=1,NDSM
               Y(NDS,I) = Y(NDS,I) + A(MU(I+1)) * X(NDS,I)
   35       CONTINUE
   31    CONTINUE
      ENDIF
C
      RETURN
      END
