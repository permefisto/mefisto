      SUBROUTINE DRCR1D( NTDL, NCODSA, MUDL, A, NDSM, B0, NIVEAU,
     +                   B,    IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: DESCENTE ET/OU REMONTEE D UN SYSTEME LINEAIRE FACTORISE SOUS
C ---- LA FORME  A * B = L * D * TL * B = B0 DITE DE CROUT
C      L EST UNE MATRICE TRIANGULAIRE INFERIEURE A DIAGONALE UNITE
C      D EST UNE MATRICE DIAGONALE INVERSIBLE
C      STOCKEES SOUS FORME PROFIL DANS A
C      CRMC1D DOIT ETRE EXECUTE AUPARAVANT
C
C ENTREES:
C --------
C NTDL   : ORDRE DE LA MATRICE A
C NCODSA : =-1 MATRICE PROFIL NON SYMETRIQUE
C          = 0 MATRICE DIAGONALE
C          = 1 MATRICE PROFIL     SYMETRIQUE
C MUDL   : SI NCODSA NON NUL ALORS
C          MUDL EST LE TABLEAU DE POINTEURS SUR CHAQUE COEFFICIENT DIAGONAL DE A
C          MUDL(0) = 0, MUDL(I) = ADRESSE DANS A DU I-EME COEFFICIENT
C                                 DIAGONAL  SI A EST NON DIAGONALE
C A      : MATRICE L, D RESULTATS DE LA FACTORISATION DE CROUT
C          ISSUE DE L EXECUTION DE CRMC1D
C NDSM   : NOMBRE DE SECONDS MEMBRES
C B0     : TABLEAU B0(NDSM, NTDL) DES NDSM SECONDS MEMBRES
C
C NIVEAU : 0  TL*B=B0      REMONTEES SEULEMENT VALABLE SI B=B0 SINON ERREUR
C          1  L *B=B0      DESCENTES SEULEMENT
C          2  L*D*X=B0     DESCENTES ET PRODUITS DIAGONAUX
C          3  L*D*TL*B=B0  RESOLUTIONS TOTALES
C          4  D*TL*B=B0    PRODUITS DIAGONAUX ET REMONTEES
C                          SEULEMENT VALABLE SI B0 = B EN ENTREE
C
C SORTIES:
C --------
C B      : TABLEAU B(NDSM, NTDL) DES SOLUTIONS
C          SI NIVEAU =0 OU 4 B0 ET B DOIVENT ETRE IDENTIQUES A L APPEL
C IERR   : 0 SI PAS D'ERREUR DETECTEE, 1 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET UPMC ANALYSE NUMERIQUE PARIS        AOUT 2006
C2345X7..............................................................012
      INTEGER           MUDL(0:NTDL)
      DOUBLE PRECISION  A(*), B(NDSM,NTDL), B0(NDSM,NTDL)
      DOUBLE PRECISION  S
c
ccc      print *
ccc      do 1 i=1,ntdl
ccc         print *,'DRCR1D: ligne',i,' ADiag=',A(MUDL(i)),'  b=',b0(1,i)
ccc 1    continue
C
      IF( NCODSA .EQ. 0 ) THEN
C
C        MATRICE DIAGONALE
C        =================
         CALL AXEBDD( NTDL, NDSM, A, B0,  B, IERR )
         RETURN
C
      ENDIF
C
C     MATRICE NON DIAGONALE:
C     ======================
      IF( NIVEAU .EQ. 0 ) GOTO 300
      IF( NIVEAU .EQ. 4 ) GOTO 200
C
C     LES DESCENTES
C     -------------
      NCDIAG0 = 0
      DO 10 JP = 1,NTDL
C        NO DU COEFFICIENT DIAGONAL
         NCDIAG = MUDL(JP)
C        NOMBRE DE COEFFICIENTS DE LA LIGNE JP
         NBCL = NCDIAG - NCDIAG0
C        NO-1 1-ER COLONNE
         JPMI = JP - NBCL
C
         DO 8 NC = 1,NDSM
            S = 0.D0
            DO 5 K = 1,NBCL-1
               S = S + A(NCDIAG0+K) * B(NC,JPMI+K)
    5       CONTINUE
            B(NC,JP) = B0(NC,JP) - S
    8    CONTINUE
C
C        PASSAGE A LA LIGNE SUIVANTE DE LA MATRICE
         NCDIAG0 = NCDIAG
   10 CONTINUE
C
      IF( NIVEAU .EQ. 1 ) GOTO 900
C
C     D*Z = Y
C     -------
  200 DO 25 JP = 1,NTDL
         S  = 1.D0 / A( MUDL(JP) )
         DO 20 NC = 1, NDSM
            B(NC,JP) = B(NC,JP) * S
   20    CONTINUE
   25 CONTINUE
C
      IF( NIVEAU .EQ. 2 ) GOTO 900
C
C     LES REMONTEES
C     -------------
 300  DO 50 JP = NTDL,2,-1
C        NO COEFFICIENT DE LA DIAGONALE
         NCDIAG = MUDL(JP)
C        NOMBRE-1 DE COEFFICIENTS DE LA LIGNE JP
         NBCL = NCDIAG - MUDL(JP-1) - 1
         IF( NBCL .GT. 0 ) THEN
            DO 40 K = 1, NBCL
               JPK = JP - K
               DO 35 NC = 1,NDSM
                  B(NC,JPK) = B(NC,JPK) - A(NCDIAG-K) * B(NC,JP)
   35          CONTINUE
   40       CONTINUE
         ENDIF
   50 CONTINUE
C
  900 RETURN
ccc      print *,'sortie DRCR1D: solution X'
ccc      do 901 i=1,ntdl
ccc         print *,'X',i,'=',b(1,i)
ccc 901  continue
C
      END
