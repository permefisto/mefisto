      SUBROUTINE INCHGC ( NTDL,
     %                    LPLIGN, LPCOLO, AG,
     %                    LPLIGC, LPCOLC, AGC, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FACTORISATION INCOMPLETE PAR NIVEAU DE CHOLESKY AGC ~ L * TL
C -----

C ENTREES:
C --------
C NTDL   : ORDRE DE LA MATRICE AG (ET AGC)
C LPLIGN : POINTEUR SUR LE COEFFICIENT DIAGONAL DE CHAQUE LIGNE DE AG
C LPCOLO : NO COLONNE DE CHAQUE COEFFICIENT STOCKE DE AG
C AG     : MATRICE INITIALE NON FACTORISEE
C LPLIGC : POINTEUR SUR LE COEFFICIENT DIAGONAL DE CHAQUE LIGNE DE AGC
C LPCOLC : NO COLONNE DE CHAQUE COEFFICIENT STOCKE DE AGC

C SORTIES:
C --------
C AGC    : MATRICE FACTORISEE INCOMPLETE <=> MATRICE DE PRECONDITIONNEMENT
C IERR   : =0 SI PAS D'ERREUR DETECTEE
C          =1 SI LA FACTORISATION EST INSTABLE PAR AU MOINS UN
C             PIVOT<=0 LE PIVOT QUI PRECEDE EST PRIS A SA PLACE
C          =2 SI LA FACTORISATION EST INSTABLE PAR
C             UN NOMBRE DE PIVOTS INCORRECTS DEPASSE (>NBPIMA=10)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY    ANALYSE NUMERIQUE UPMC PARIS      JUILLET 1989
C23456---------------------------------------------------------------012
      PARAMETER ( NBPIMA = 10 )
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      COMMON / EPSSSS / EPZERO, EPSXYZ
      DOUBLE PRECISION  AG(*), AGC(*),  S, PIVOT, SQRT
      INTEGER           LPLIGN(NTDL+1), LPCOLO(*)
      INTEGER           LPLIGC(NTDL+1), LPCOLC(*)

 1000 FORMAT(' ATTENTION DANS INCHGC: LE PIVOT',I8,' VAUT ',
     %G13.6,' =<',G13.6/' LE PIVOT QUI PRECEDE EST PRIS A SA PLACE')
 1500 FORMAT(' ERREUR INCHGC : LIGNE ',I8,
     %' L''INDICE DE COLONNE',I8,' N''EST PAS RETROUVE')
 2000 FORMAT(/,' INCHGC: FACTORISATION INSTABLE ABANDONNEE')

C     INITIALISATION
C     --------------
      PIVOT = SQRT( AG(LPLIGN(2)) )
      IF( PIVOT .LT. EPZERO ) THEN
         PIVOT = 1D0
      ENDIF

C     MISE A ZERO DE LA MATRICE DE PRECONDITIONNEMENT
      DO K=1,LPLIGC(NTDL+1)
         AGC(K) = 0.D0
      ENDDO

      DO I=1,NTDL
         IC1=LPLIGC(I)
         IC2=LPLIGC(I+1)
         IC0=IC1
         DO 101 KT=LPLIGN(I)+1,LPLIGN(I+1)-1
            IC1=IC0
            DO 102 KC=IC1+1,IC2-1
               IF(LPCOLO(KT).EQ.LPCOLC(KC)) THEN
                  AGC(KC)=AG(KT)
                  IC0=KC
                  GOTO 101
               ENDIF
 102        ENDDO
            WRITE (IMPRIM,1500) I,LPCOLO(KT)
 101     ENDDO
         AGC(IC2)=AG( LPLIGN(I+1) )
      ENDDO

      NBPIVO = 0

C     BOUCLE SUR LES LIGNES
C     ---------------------
      DO 100 I=1,NTDL
         S  = 0.D0
         I1 = LPLIGC(I)
         I2 = LPLIGC(I+1)
         K1 = I1 + 1
         IF(I2 .LT. K1 ) GOTO 100
         IF(I2 .EQ. K1 ) GOTO 2

C        IL EXISTE DES COEFFICIENTS I1+1 a I2-1 NON DIAGONAUX SUR LA LIGNE I
C        -------------------------------------------------------------------
         K2 = I2 - 1

C        LA BOUCLE SUR CES COEFFICIENTS
C        ------------------------------
         DO 4 K = K1 , K2
           S = 0.D0
C          J NO DE LA COLONNE A TRAITER
           J  = LPCOLC(K)
           J1 = LPLIGC(J)
           J2 = LPLIGC(J+1)
           IF(J2 .LE. J1+1) GOTO 5

C          IL EXISTE DES COEFFICIENTS DANS LA LIGNE J
C          ------------------------------------------
C          IC  NO DU COEFFICIENT COURANT DE LA LIGNE I
C          JC  NO DU COEFFICIENT COURANT DE LA LIGNE J
C          IJC NO DE LA COLONNE DU COEFFICIENT IC
C          JJC NO DE LA COLONNE DU COEFFICIENT JC
           IC  = I1
           JC  = J1

 11        IC  = IC + 1
           IJC = LPCOLC(IC)
           JC  = JC + 1
           JJC = LPCOLC(JC)

 10        IF(IJC .GT. JJC) GOTO 9
           IF(IJC .EQ. JJC) GOTO 8

           IF(JJC .GE. J) GOTO 5
           IC  = IC + 1
           IJC = LPCOLC(IC)
           GOTO 10

 9         IF(IJC .GE. J) GOTO 5
           JC  = JC + 1
           JJC = LPCOLC(JC)
           GOTO 10

C          2 COEFFICIENTS ONT MEME NUMERO DE COLONNE
C          -----------------------------------------
 8         IF(JJC .GE. J) GOTO 5
           S = S + AGC(IC) * AGC(JC)
           GOTO 11

C           LE COEFFICIENT L(I,J)
C           ---------------------
 5          AGC(K) = (AGC(K) - S) / AGC(J2)
 4       ENDDO

C        LE COEFFICIENT L(I,I)
C        ---------------------
         S = 0.D0
         DO K = K1 , K2
            S = S + AGC(K) ** 2
         ENDDO

 2       S = AGC(I2) - S
         IF ( S .LE. 0D0 ) GOTO 17
         IF ( DABS(S) .LE. EPZERO*AGC(I2) )  GOTO 17

C        PIVOT CORRECT
C        -------------
         AGC(I2) = SQRT( S )
         PIVOT   = AGC(I2)
         GOTO 100

C        PIVOT INCORRECT => CHANGEMENT DE NIVEAU
C        ---------------------------------------
 17      NBPIVO = NBPIVO + 1
         IF( NBPIVO .LE. NBPIMA ) THEN
            WRITE (IMPRIM,1000) I,S,EPZERO
            AGC(I2) = PIVOT
            IERR = 1
         ELSE
C           NOMBRE DE PIVOTS INCORRECTS DEPASSE
            WRITE (IMPRIM,2000)
            IERR = 2
           RETURN
         END IF
 100  ENDDO

      RETURN
      END
