      SUBROUTINE DRCH1D( NTDL, NCODSA, LPDIAG, A, NDSM, B0, NIVEAU,
     %                   B   )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: DESCENTE ET/OU REMONTEE D UN SYSTEME LINEAIRE FACTORISE SOUS
C ---- LA FORME  A * B = L * TL * B = B0  DITE DE CHOLESKY
C      L EST UNE MATRICE TRIANGULAIRE INFERIEURE RANGEE SOUS FORME
C      PROFIL . CHMC1D DOIT ETRE EXECUTE AVANT.
C      VERSION REELLE DOUBLE PRECISION
C
C ENTREES:
C --------
C NTDL   : ORDRE DE LA MATRICE A
C NCODSA : 0 SI MATRICE DIAGONALE, >0 SI MATRICE SYMETRIQUE NON DIAGONALE
C LPDIAG : SI NCODSA>0 ALORS POINTEUR SUR CHAQUE COEFFICIENT DIAGONAL DE A
C          LPDIAG(1) = 0, LPDIAG(I+1) = ADRESSE DANS A DU I-EME COEFFICIENT
C                                       DIAGONAL SI A NON DIAGONALE
C A      : MATRICE L RESULTAT DE LA FACTORISATION DE CHOLESKY ISSUE DE CHOLPR
C NDSM   : NOMBRE DE SECONDS MEMBRES DANS B0 (ET B)
C B0     : TABLEAU B0(NDSM,NTDL) DES NDSM SECONDS MEMBRES
C NIVEAU : =0  TL*B=B0   REMONTEES SEULEMENT VALABLE SI B=B0 SINON ERREUR
C          =1  L *B=B0   DESCENTES SEULEMENT
C          >1  L*TL*B=B0 DESCENTES ET REMONTEES
C
C SORTIES:
C --------
C B      : TABLEAU B(NTDL,NDSM) DES SOLUTIONS
C          SI NIVEAU > 1 B0 ET B PEUVENT ETRE IDENTIQUES A L APPEL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1989
C23456---------------------------------------------------------------012
      DOUBLE PRECISION S,A,B,B0
      DIMENSION        A(1),B(NTDL,NDSM),B0(NTDL,NDSM),LPDIAG(NTDL+1)
C
      IF( NCODSA .EQ. 0 ) GOTO 1000
C
C     MATRICE NON DIAGONALE
C     =====================
C
C     LES DESCENTES
C     -------------
      IF(NIVEAU.EQ.0) GOTO 29
      JA=1
      DO 1 JP=1,NTDL
         IA   = LPDIAG(JP+1) - LPDIAG(JP) - 1
         IB   = JP - IA
         JA0  = JA
         DO 21 NC=1,NDSM
            S  = 0.D0
            JA = JA0
            DO 3 IAU = IB , JP-1
               S  = S + A(JA) * B(IAU,NC)
               JA = JA+1
    3       CONTINUE
            B(JP,NC) = ( B0(JP,NC) - S ) / A(JA)
   21    CONTINUE
         JA=JA+1
    1 CONTINUE
      IF(NIVEAU.EQ.1) GOTO 100
C
C     LES REMONTEES
C     -------------
   29 IA = LPDIAG(NTDL+1) + 1
      DO 6 JP=NTDL,1,-1
         IAU = LPDIAG(JP+1)-LPDIAG(JP)-1
         IA  = IA - 1
         S   = 1.D0 / A(IA)
         DO 11 NC=1,NDSM
            B(JP,NC) = B(JP,NC) * S
   11    CONTINUE
         DO 9 JA=1,IAU
            IA  = IA-1
            S   = A(IA)
            DO 10 NC=1,NDSM
               B(JP-JA,NC) = B(JP-JA,NC) - B(JP,NC) * S
   10       CONTINUE
    9    CONTINUE
    6 CONTINUE
C
  100 RETURN
C
C     MATRICE DIAGONALE
C     =================
C
C     LES DESCENTES
C     -------------
 1000 IF(NIVEAU.EQ.0) GOTO 1100
      DO 1001 I=1,NTDL
         DO 1002 NC=1,NDSM
            B(I,NC)=B0(I,NC)/A(I)
 1002    CONTINUE
 1001 CONTINUE
C
C     LES REMONTEES
C     -------------
 1100 IF(NIVEAU.EQ.1) RETURN
      DO 1101 I=1,NTDL
         DO 1102 NC=1,NDSM
            B(I,NC)=B(I,NC)/A(I)
 1102    CONTINUE
 1101 CONTINUE
C
      RETURN
      END
