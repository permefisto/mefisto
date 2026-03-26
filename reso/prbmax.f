      SUBROUTINE PRBMAX(NCODSA,NDSM,NTDL,MU,A,X,Y)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: CALCULER LE RESIDU Y = Y - A * X DU SYSTEME LINEAIRE
C ----
C      Y , X TABLEAUX (NDSM,NTDL)  A MATRICE PROFIL SYMETRIQUE
C                                  OU DIAGONALE
C                                  OU NON SYMETRIQUE
C PARAMETRES D ENTREE:
C --------------------
C NCODSA : 0 SI LA MATRICE A EST DIAGONALE
C          1 SI LA MATRICE EST SYMETRIQUE PROFIL
C         -1 SI LA MATRICE EST NON-SYMETRIQUE PROFIL
C NDSM   : NOMBRE DE CAS OU SECONDS MEMBRES
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE OU ORDRE DE LA MATRICE
C MU     : POINTEUR SUR CHAQUE COEFFICIENT DIAGONAL
C A      : MATRICE PROFIL  D ORDRE NTDL( CF S.D MUA)
C X      : TABLEAU(NTDL,NDSM) A MULTIPLIER PAR A
C
C PARAMETRE MODIFIE  :
C --------------------
C Y      : TABLEAU(NTDL,NDSM)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR: A.PERRONNET ANALYSE NUMERIQUE PARIS  JUILLET 1984
C.......................................................................
      DOUBLE PRECISION A(*), X(NTDL,NDSM), Y(NTDL,NDSM)
      INTEGER          MU(*)
C
C     AIGUILLAGE SELON NCODSA
C     -----------------------
C
      IF(NCODSA.EQ.0) GOTO 40
C
      IF(NCODSA.LT.0) GOTO 30
C
C     A SYMETRIQUE . MU PAR D.L
C     ------------------------
C
      JA=0
           DO 21 J=1,NTDL
           JMI=MU(J+1)-MU(J)
           IF(JMI.LE.1) GOTO 22
           I=J-JMI
           JMI=JMI-1
                DO 23 JD=1,JMI
                JA=JA+1
                I=I+1
                      DO 24 NDS=1,NDSM
                      Y(J,NDS)=Y(J,NDS)-A(JA)*X(I,NDS)
                      Y(I,NDS)=Y(I,NDS)-A(JA)*X(J,NDS)
   24                 CONTINUE
   23           CONTINUE
   22      JA=JA+1
                 DO 25 NDS=1,NDSM
                 Y(J,NDS)=Y(J,NDS)-A(JA)*X(J,NDS)
   25            CONTINUE
   21      CONTINUE
       RETURN
C
C      A NON-SYMETRIQUE . MU PAR D.L
C      -----------------------------
C
   30       DO 31 I=1,NTDL
            MUI=MU(I)
            IH=(MU(I+1)-MUI)/2
            IF(IH.EQ.0) GOTO 34
            MUJ=MUI+IH
            ILI=I-IH
                 DO 32 J=1,IH
                 IMUI=MUI+J
                 JMUJ=MUJ+J
                      DO 33 NDS=1,NDSM
                      Y(I,NDS)=Y(I,NDS)-A(IMUI)*X(ILI,NDS)
                      Y(ILI,NDS)=Y(ILI,NDS)-A(JMUJ)*X(I,NDS)
   33                 CONTINUE
                 ILI=ILI+1
   32            CONTINUE
C
   34            DO 35 NDS=1,NDSM
                 Y(I,NDS)=Y(I,NDS)-A(MU(I+1))*X(I,NDS)
   35            CONTINUE
   31       CONTINUE
      RETURN
C
C      A  DIAGONALE
C      ------------
C
   40       DO 41 I=1,NTDL
                  DO 42 NDS=1,NDSM
                  Y(I,NDS)=Y(I,NDS)-A(I)*X(I,NDS)
   42             CONTINUE
   41       CONTINUE
C
      RETURN
      END
