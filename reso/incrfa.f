      SUBROUTINE INCRFA( NTDL, LPLIGA, LPCOLA, AG, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FACTORISATION INCOMPLETE PAR NIVEAU DE CROUT AG # L * D * TL
C -----    SUR ELLE-MEME MATRICE MORSE

C ENTREES:
C --------
C NTDL   : ORDRE DE LA MATRICE AG
C LPLIGA : POINTEUR SUR LE COEFFICIENT DIAGONAL DE CHAQUE LIGNE DE AG
C LPCOLA : NO COLONNE DE CHAQUE COEFFICIENT STOCKE DE AG

C MODIFIE:
C --------
C AG     : MATRICE AVANT ET APRES FACTORISATION

C SORTIE :
C --------
C IERR   : =0 SI PAS D'ERREUR DETECTEE
C          =1 SI LA FACTORISATION EST INSTABLE PAR AU MOINS UN
C             ABS(PIVOT)<=eps LE PIVOT QUI PRECEDE EST PRIS A SA PLACE
C          =2 SI LA FACTORISATION EST INSTABLE PAR
C             UN NOMBRE DE PIVOTS INCORRECTS DEPASSE (>NBPIMA=10)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS     JUILLET 1989
C MODIFS : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris OCTOBRE 2007
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  MARS    2009
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      COMMON / EPSSSS / EPZERO, EPSXYZ
      INTEGER           LPLIGA(NTDL+1), LPCOLA(*)
      DOUBLE PRECISION  AG(*)
      DOUBLE PRECISION  S, PIVOT

 1000 FORMAT(' ATTENTION DANS INCRFA: LE PIVOT',I8,' VAUT ',
     %G13.6,' =<eps=',G13.6/' LE PIVOT=',G13.6,
     %' QUI PRECEDE EST PRIS A SA PLACE')
 2000 FORMAT(/,' INCRFA: FACTORISATION INSTABLE ABANDONNEE')

      print *
      print *,'MATRICE AG debut INCRFA   #############################'
      DO i=1,15
         print 10010,(i,LPCOLA(m),ag(m),m=LPLIGA(i)+1, LPLIGA(i+1))
      enddo
10010 format(5('  ag(',i3,',',i3,')=',D15.6))

C     FACTORISATION INCOMPLETE DE CROUT
C     =================================
      NBPIMA = 10
      NBPIVO = 0
C     LE PREMIER COEFFICIENT DIAGONAL DE AG
      PIVOT = AG( LPLIGA(2) )
      IF( ABS(PIVOT) .EQ. 0D0 ) THEN
         PIVOT = 1D0
      ENDIF

C     BOUCLE SUR LES LIGNES
C     ---------------------
      DO 1 I=1,NTDL
         S  = 0.D0

C        PREMIER ET DERNIER COEFFICIENT DE LA LIGNE I
         I1 = LPLIGA(I)
         I2 = LPLIGA(I+1)
         IF( I2 .EQ. I1+1 ) GOTO 2

C        IL EXISTE DES COEFFICIENTS NON DIAGONAUX SUR LA LIGNE I
C        -------------------------------------------------------
         K1 = I1 + 1
         K2 = I2 - 1

C        LA BOUCLE SUR CES COEFFICIENTS
C        ------------------------------
         DO 4 K = K1, K2
           S = 0.D0
C          J NO DE LA COLONNE A TRAITER
           J = LPCOLA(K)

C          PREMIER ET DERNIER COEFFICIENT DE LA LIGNE J
           J1 = LPLIGA(J)
           J2 = LPLIGA(J+1)
           IF( J2 .EQ. J1+1 ) GOTO 5

C          IL EXISTE DES COEFFICIENTS DANS LA LIGNE J
C          ------------------------------------------
C          IC  NO DU COEFFICIENT COURANT DE LA LIGNE I
C          JC  NO DU COEFFICIENT COURANT DE LA LIGNE J
C          IJC NO DE LA COLONNE DU COEFFICIENT IC
C          JJC NO DE LA COLONNE DU COEFFICIENT JC
           IC  = I1
           JC  = J1

 11        IC  = IC + 1
           IJC = LPCOLA(IC)

           JC  = JC + 1
           JJC = LPCOLA(JC)

 10        IF( IJC .GT. JJC ) GOTO 9
           IF( IJC .EQ. JJC ) GOTO 8

           IF( JJC .GE. J ) GOTO 5
           IC  = IC + 1
           IJC = LPCOLA(IC)
           GOTO 10

 9         IF( IJC .GE. J ) GOTO 5
           JC  = JC + 1
           JJC = LPCOLA(JC)
           GOTO 10

C          2 COEFFICIENTS ONT MEME NUMERO DE COLONNE
C          -----------------------------------------
 8         IF( JJC .GE. J ) GOTO 5
           S = S + AG(IC) * AG(LPLIGA(JJC+1)) * AG(JC)
           GOTO 11

C           LE COEFFICIENT L(I,J)
C           ---------------------
 5          AG(K) = (AG(K) - S) / AG(J2)

 4       ENDDO

C        LE COEFFICIENT L(I,I)
C        ---------------------
         S = 0.D0
         DO 12 K = K1 , K2
            JJC = LPCOLA(K)
            S = S + AG( LPLIGA(JJC+1) ) * AG(K) ** 2
 12      ENDDO

 2       S = AG(I2) - S
         IF ( DABS(S) .LE. 0D0 )   GOTO 17
         IF ( DABS(S) .LE. EPZERO*AG(I2) )  GOTO 17

C        PIVOT CORRECT
C        -------------
         AG(I2) = S
         PIVOT  = S
         GOTO 1

C        PIVOT INCORRECT => CHANGEMENT DE NIVEAU A DEMANDER
C        --------------------------------------------------
 17      NBPIVO = NBPIVO + 1
         IF( NBPIVO .LE. NBPIMA ) THEN
            WRITE (IMPRIM,1000) I, S, EPZERO, PIVOT
            AG(I2) = PIVOT
            IERR = 1
         ELSE
C           NOMBRE DE PIVOTS INCORRECTS DEPASSE
            WRITE (IMPRIM,2000)
            IERR = 2
           RETURN
         ENDIF
 1    ENDDO

      print *
      print *,'MATRICE AG # L D tL fin INCRFA ####################'
      do i=1,15
         print 10010,(i,LPCOLA(m),ag(m),m=LPLIGA(i)+1, LPLIGA(i+1))
      enddo

      RETURN
      END
