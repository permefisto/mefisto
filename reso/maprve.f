      SUBROUTINE MAPRVE( NUINIT, CTE, NTDL,
     %                   NCODSA, MU,  A,    X,
     %                   Y )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:    Y=CTE*A*X ou Y=Y+CTE*A*X avec
C ----    Y,X VECTEURS (NTDL)  A MATRICE PROFIL     SYMETRIQUE
C                             OU MATRICE DIAGONALE
C                             OU MATRICE PROFIL NON SYMETRIQUE
C
C ENTREES:
C --------
C NUINIT : 0 => Y=  CTE*A*X
C          1 => Y=Y+CTE*A*X
C CTE    : REEL DOUBLE PRECISION
C NTDL   : ORDRE DE LA MATRICE ET DES VECTEURS
C
C NCODSA : 0 SI LA MATRICE A EST DIAGONALE
C          1 SI LA MATRICE EST SYMETRIQUE PROFIL
C         -1 SI LA MATRICE EST NON-SYMETRIQUE PROFIL
C MU     : POINTEUR SUR CHAQUE COEFFICIENT DIAGONAL DE LA MATRICE A
C A      : MATRICE PROFIL D ORDRE NTDL
C X      : VECTEUR DE NTDL COMPOSANTES
C
C SORTIE ou MODIFIE:
C ------------------
C Y      : VECTEUR DE NTDL COMPOSANTES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1999
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  CTE, A(1:*), X(NTDL), Y(NTDL), AUX
      INTEGER           MU(1:NTDL+1)
C
      IF( CTE .EQ. 1D0 ) THEN
C
         IF( NCODSA .EQ. 0 ) THEN
C
C           A  DIAGONALE
C           ------------
            IF( NUINIT .EQ. 0 ) THEN
               DO 5 I=1,NTDL
                  Y(I) = A(I) * X(I)
    5          CONTINUE
            ELSE
               DO 8 I=1,NTDL
                  Y(I) = Y(I) + A(I) * X(I)
    8          CONTINUE
            ENDIF
            RETURN
         ENDIF
C
C        INITIALISATION A 0 DU TABLEAU Y
C        -------------------------------
         IF( NUINIT .EQ. 0 ) THEN
            DO 10 I=1,NTDL
               Y(I) = 0D0
   10       CONTINUE
         ENDIF
C
         IF( NCODSA .GT. 0 ) THEN
C
C           A SYMETRIQUE . MU PAR D.L
C           -------------------------
            JA = 0
            DO 30 J=1,NTDL
               JMI = MU(J+1) - MU(J)
               I   = J - JMI
               DO 20 JD=1,JMI-1
                  JA = JA + 1
                  I  = I  + 1
                  Y(J) = Y(J) + A(JA) * X(I)
                  Y(I) = Y(I) + A(JA) * X(J)
   20          CONTINUE
               JA   = JA + 1
               Y(J) = Y(J) + A(JA) * X(J)
   30       CONTINUE
C
         ELSE
C
C           A NON-SYMETRIQUE . MU PAR D.L
C           -----------------------------
            DO 50 I=1,NTDL
               MUI = MU(I)
               IH  = ( MU(I+1) - MUI ) / 2
               IF( IH .GT. 0 ) THEN
                  MUJ = MUI + IH
                  II  = I   - IH
                  DO 40 J=1,IH
                     IMUI  = MUI + J
                     JMUJ  = MUJ + J
                     Y(I ) = Y(I ) + A(IMUI) * X(II)
                     Y(II) = Y(II) + A(JMUJ) * X(I )
                     II    = II + 1
   40             CONTINUE
               ENDIF
               Y(I) = Y(I) + A(MU(I+1)) * X(I)
   50       CONTINUE
         ENDIF
C
      ELSE
C
C        LA CTE N'EST PAS 1D0 => PRODUIT
C        ===============================
         IF( NCODSA .EQ. 0 ) THEN
C
C           A  DIAGONALE
C           ------------
            IF( NUINIT .EQ. 0 ) THEN
               DO 105 I=1,NTDL
                  Y(I) = CTE * A(I) * X(I)
  105          CONTINUE
            ELSE
               DO 108 I=1,NTDL
                  Y(I) = Y(I) + CTE * A(I) * X(I)
  108          CONTINUE
            ENDIF
            RETURN
         ENDIF
C
C        INITIALISATION A 0 DU TABLEAU Y
C        -------------------------------
         IF( NUINIT .EQ. 0 ) THEN
            DO 110 I=1,NTDL
               Y(I) = 0D0
  110       CONTINUE
         ENDIF
C
         IF( NCODSA .GT. 0 ) THEN
C
C           A SYMETRIQUE . MU PAR D.L
C           -------------------------
            JA = 0
            DO 130 J=1,NTDL
               JMI = MU(J+1) - MU(J)
               I   = J - JMI
               AUX = CTE * X(J)
               DO 120 JD=1,JMI-1
                  JA = JA + 1
                  I  = I  + 1
                  Y(J) = Y(J) + A(JA) * CTE * X(I)
                  Y(I) = Y(I) + A(JA) * AUX
  120          CONTINUE
               JA = JA + 1
               Y(J) = Y(J) + A(JA) * AUX
  130       CONTINUE
C
         ELSE
C
C           A NON-SYMETRIQUE . MU PAR D.L
C           -----------------------------
            DO 150 I=1,NTDL
               AUX = CTE * X(I)
               MUI = MU(I)
               IH  = ( MU(I+1) - MUI ) / 2
               IF( IH .GT. 0 ) THEN
                  MUJ = MUI + IH
                  II  = I   - IH
                  DO 140 J=1,IH
                     IMUI  = MUI + J
                     JMUJ  = MUJ + J
                     Y(I ) = Y(I ) + A(IMUI) * CTE * X(II)
                     Y(II) = Y(II) + A(JMUJ) * AUX
                     II    = II + 1
  140             CONTINUE
               ENDIF
               Y(I) = Y(I) + A(MU(I+1)) * AUX
  150       CONTINUE
         ENDIF
      ENDIF
C
      RETURN
      END
