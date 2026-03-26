      SUBROUTINE CHOLPR( NTDL, NCODSA, LPDIAG, A0,
     %                   A,    IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: FACTORISER UNE MATRICE SYMETRIQUE DEFINIE POSITIVE SOUS LA
C ---- FORME  A = L * TL DITE DE CHOLESKY
C      A ET L SONT  TRIANGULAIRES INFERIEURES STOCKEES SOUS FORME PROFIL
C      C-A-D PAR LIGNES DU 1-ER COEFFICIENT NON NUL AU COEFFICIENT
C      DIAGONAL. LE PROFIL EST LE MEME POUR A ET L
C      A0 ET A=L PEUVENT ETRE CONFONDUS EN ENTREE
C
C      VERSION REELLE DOUBLE PRECISION
C      VERSION MODIFIEE POUR PRENDRE EN COMPTE LE CAS A SINGULIERE
C
C ENTREES:
C --------
C NTDL   : ORDRE DE LA MATRICE A
C NCODSA : 0     SI MATRICE DIAGONALE
C          NON 0 SI MATRICE SYMETRIQUE NON DIAGONALE
C LPDIAG : SI NCODSA>0 ALORS POINTEUR SUR LE COEFFICIENT DIAGONAL DE D.L
C          LPDIAG(0)=0, LPDIAG(I)=ADRESSE DANS A DU I-EME COEF DIAG
C                                 DIAGONAL SI A NON DIAGONALE
C A0     : MATRICE A FACTORISER SOUS FORME L * TL AVEC STOCKAGE PROFIL
C
C SORTIES:
C --------
C A      : MATRICE FACTORISEE L POUR TOUTE VALEUR DE NCODSA
C IERR   : CODE D'ERREUR 0 SI PAS D'ERREUR, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MAI  1989
C MODIFIE: PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS        MARS 1990
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      COMMON / EPSSSS / EPZERO, EPSXYZ
      DOUBLE PRECISION  SA, A0(*), A(*), PIVOT
      INTEGER           LPDIAG(0:NTDL)
      INTRINSIC         SQRT, ABS
C
10010 FORMAT(' ATTENTION - CHOLPR :',I8,'EME PIVOT=',G13.6,' < 0.',
     %/,' LE PIVOT PRECEDENT',G13.6,' EST PRIS A LA PLACE')
10020 FORMAT(/,' ATTENTION - CHOLPR : FACTORISATION INSTABLE')
C
C     SEUIL POUR DETECTER LES PIVOTS INCORRECTS
C     =========================================
      IERR = 0
      IF( A0(1) .LE. 0D0 ) THEN
         IERR = 3
         RETURN
      ENDIF
C
      PIVOT = SQRT( ABS(A0(1)) )
      IF( PIVOT .LT. EPZERO ) THEN
         PIVOT = 1.D0
      ENDIF
C
      IF( NCODSA .EQ. 0 ) GOTO 1000
C
C     MATRICE SYMETRIQUE NON DIAGONALE
C     ================================
      NBPIMA = 10
      NBPIVO = 0
C     I INDIQUE LES COLONNES IP
C     J INDIQUE LES LIGNES   JP
      DO 50 JP = 1, NTDL
         JP1  = JP - 1
         JA0  = LPDIAG(JP1)
         NBCJ = LPDIAG(JP ) - JA0 - 1
         JPMI = JP - NBCJ
C
C        CALCUL DES COEFFICIENTS NON DIAGONAUX
         DO 10 IP = JPMI, JP1
C
C           COEFFICIENT NON DIAGONAL L(JP,IP)
            IA0  = LPDIAG(IP-1)
            NBCI = LPDIAG(IP  ) - IA0 - 1
            IPMI = IP - NBCI
            IF(IPMI.GT.JPMI) THEN
               IA = IA0
               JA = JA0 + IPMI - JPMI
               KP = IPMI
            ELSE
               IA = IA0 + JPMI - IPMI
               JA = JA0
               KP = JPMI
            ENDIF
            NBC = IP - KP
            SA  = 0.D0
            DO 5 KD = 1, NBC
               SA = SA + A(IA+KD) * A(JA+KD)
  5         CONTINUE
            NBC  = NBC + 1
            JA   = JA  + NBC
            A(JA) = ( A0(JA) - SA ) / A(IA+NBC)
ccc         WRITE(IMPRIM,*) 'L(',JP,',',IP,')=',A(JA)
 10      CONTINUE
C
C        COEFFICIENT DIAGONAL L(JP,JP)
         SA = 0.D0
         DO 20 KD = 1, NBCJ
            SA = SA + A(JA0+KD)**2
  20     CONTINUE
         JA = JA0 + NBCJ + 1
         SA = A0(JA) - SA
         IF( SA .LE. 0.D0 ) THEN
C           COEFFICIENT DIAGONAL NEGATIF
            NBPIVO = NBPIVO + 1
            IF( NBPIVO .LE. NBPIMA ) THEN
               WRITE (IMPRIM,10010) JP, SA, PIVOT
               IERR = 1
               SA = PIVOT
            ELSE
               WRITE (IMPRIM,10020)
               IERR = 2
              RETURN
            END IF
         ENDIF
C
         IF( SA - EPZERO * A0(JA) .LE. 0.D0 ) THEN
C           ON EST AU DESSOUS DU SEUIL DE PRECISION
            IERR = 1
            SA   = PIVOT
         ELSE
C           SAUVEGARDE DU PIVOT ACTUEL
            PIVOT = SA
         ENDIF
         A(JA) = SQRT( SA )
ccc
ccc      WRITE(IMPRIM,*) 'L(',JP,',',JP,')=',A(JA)
ccc      WRITE(IMPRIM,*)
C
 50   CONTINUE
      RETURN
C
C     MATRICE DIAGONALE
C     =================
 1000 DO 1001 IP = 1, NTDL
C
         IF( A0(IP) .LE. 0.D0 ) THEN
C
C           COEFFICIENT DIAGONAL INCORRECT
            IERR = 1
            WRITE(IMPRIM,10010) IP, A0(IP), PIVOT
            A(IP) = PIVOT
C
         ELSE
C
C           COEFFICIENT DIAGONAL CORRECT
            A(IP) = SQRT( A0(IP) )
            PIVOT = A(IP)
C
         ENDIF
 1001 CONTINUE
C
      RETURN
      END
