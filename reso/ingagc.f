      SUBROUTINE INGAGC( NTDL,   LPLIGN, LPCOLO, AG,
     %                   LPLIGC, LPDILU, LPCOLC, AGC,
     %                   NIVEAU, CX, IND, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :         FACTORISATION INCOMPLETE PAR NIVEAU DE GAUSS AG # L * U
C -----
C
C PARAMETRES D ENTREE :
C ---------------------
C NTDL          : ORDRE DE LA MATRICE
C LPLIGN,LPCOLO : POINTEURS ASSOCIES A LA MATRICE MORSE AG
C AG            : MATRICE DU SYSTEME LINEAIRE
C NIVEAU        : NIVEAU DE LA FACTORISATION INCOMPLETE
C CX,IND        : TABLEAUX DE TRAVAIL
C
C PARAMETRE RESULTAT :
C ------------------
C AGC   : MATRICE FACTORISEE
C LPLIGC,LPDILU,LPCOLC : POINTEURS ASSOCIES A LA MATRICE MORSE AGC
C IERR          : CODE D'ERREUR ( = 1 SI LA FACTORISATION EST INSTABLE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Pascal JOLY LABORATOIRE D'ANALYSE NUMERIQUE PARIS 6  MAI 1989
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  AG, AGC, CX, PIVOT
      DIMENSION         LPLIGN(NTDL+1), LPCOLO(*), AG(*)
      DIMENSION         LPLIGC(NTDL+1), LPCOLC(*), AGC(*)
      DIMENSION         LPDILU(NTDL), IND(NTDL), CX(NTDL)
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      COMMON / EPSSSS / EPZERO, EPSXYZ
C
      WRITE (IMPRIM,1000) NIVEAU
C
C     INITIALISATION
C     --------------
      EPS   = EPZERO
      SIGMA = 1.D0
      PIVOT = DSQRT( AG(LPLIGN(2)) ) * SIGMA
      IF( PIVOT .LT. EPS ) THEN
         PIVOT = SIGMA
      ENDIF
      DO 200 K=1,LPLIGC(NTDL+1)
         AGC(K) = 0.D0
 200  CONTINUE
      NBPIVO = 0
C
C     BOUCLE SUR LES LIGNES
C     ---------------------
      DO 1 I=1,NTDL
         KC1 = LPLIGC(I) + 1
         KC2 = LPLIGC(I+1)
C
C        INITIALISATION DE IND ET CX
         DO 2 K=KC1,KC2
            J=LPCOLC(K)
            IND(J)=I
            CX(J)=0.D0
2        CONTINUE
C
C        COPIE DE LA LIGNE I DE AG DANS CX
         K1 = LPLIGN(I) + 1
         K2 = LPLIGN(I+1)
         DO 3 K=K1,K2
            J=LPCOLO(K)
            CX(J)=AG(K)
3        CONTINUE
C
C        STABILISATION PAR AMPLIFICATION DE LA DIAGONALE
         CX(I)=CX(I)*SIGMA
C
C        FACTORISATION INCOMPLETE :
C        BOUCLE D'ELIMINATION SUR LA PARTIE
C        TRIANGULAIRE INFERIEURE DE AGC
         DO 4 K=KC1,LPDILU(I)
            J=LPCOLC(K)
            CX(J)=CX(J)/AGC(LPLIGC(J+1))
            DO 5 L=LPDILU(J)+1,LPLIGC(J+1)-1
               JJ=LPCOLC(L)
               IF(IND(JJ).EQ.I) CX(JJ)=CX(JJ)-CX(J)*AGC(L)
 5          CONTINUE
 4       CONTINUE
C
C        COPIE DE CX DANS LA LIGNE I DE AGC
         DO 6 K=KC1,KC2-1
            J=LPCOLC(K)
            AGC(K)=CX(J)
6        CONTINUE
C
C        LE COEFFICIENT DIAGONAL
         IF(DABS(CX(I)).LT.EPSXYZ) THEN
              NBPIVO=NBPIVO+1
              CX(I)=PIVOT
         ELSE
              PIVOT=CX(I)
         END IF
         AGC(LPLIGC(I+1))=CX(I)
C
1     CONTINUE
C
C     CONTROLE DES PIVOTS
C     -------------------
      WRITE (IMPRIM,3000) LPLIGN(NTDL+1),LPLIGC(NTDL+1)
      IF(NBPIVO.GT.0) THEN
         IERR = 1
         WRITE (IMPRIM,2000) NBPIVO,NTDL
      ENDIF
C
      RETURN
1000  FORMAT(/,1X,131(1H*),
     S       /,' FACTORISATION INCOMPLETE DE A DE NIVEAU',I4/,
     S       1X,43(1H-)/)
2000  FORMAT(I6,' COEFFICIENTS DIAGONAUX SUR',I6,' SONT NULS')
3000  FORMAT(' NOMBRE DE COEFFICIENTS NON NULS DE AG ',I6/
     S       ' NOMBRE DE COEFFICIENTS NON NULS DE AGC',I6)
      END
