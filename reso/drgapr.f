      SUBROUTINE DRGAPR( LPDIAG, A, B0, NTDL, NIVEAU,  B )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DESCENTE ET/OU REMONTEE D'UN SYSTEME DONT LA MATRICE EST
C -----    FACTORISEE SOUS LA FORME A = L * U  DITE DE GAUSS
C          L EST UNE MATRICE TRIANGULAIRE INFERIEURE A DIAGONALE UNITE
C          U EST UNE MATRICE TRIANGULAIRE SUPERIEURE
C          L ET U SONT RANGEES DANS A SOUS FORME PROFIL
C          A(I,J) NON NUL => A(J,I) NON NUL  (PROFIL SYMETRIQUE)
C          A NON SYMETRIQUE LPDIAG POINTE SUR LES COEFFICIENTS DIAGONAUX
C          LE SP FAGAPR A DU ETRE EXECUTE AUPARAVANT
C
C ENTREES:
C --------
C LPDIAG : LPDIAG(0)=0, LPDIAG(I)= ADRESSE A DU I-EME COEF DIAGONAL
C A      : MATRICE FACTORISEE A=L*U AVEC LE RANGEMENT
C          1 3    11   |    L PARTIE TRIANGULAIRE INFERIEURE
C          2 4  6 12   |    U PARTIE TRIANGULAIRE SUPERIEURE+DIAGONALE
C            5  7 13   V
C          8 9 10 14
C            --->
C B0     : LE SECOND MEMBRE B0(NTDL)
C NTDL   : L'ORDRE DE LA MATRICE A=LU
C NIVEAU : 0  U*B=B0    REMONTEE SEULE EXACTE SEULEMENT SI B = B0
C          1  L*B=B0    DESCENTE SEULE
C          2  L*U*B=B0  DESCENTE et REMONTEE
C
C SORTIE :
C --------
C B      :  TABLEAU DES RESULTATS  B(NTDL)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LAN 189 PARIS ET IRIA            OCTOBRE 1977
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY    AVRIL 2010
C.......................................................................
      INTEGER           LPDIAG(0:NTDL)
      DOUBLE PRECISION  A(1:*), B0(1:*), B(1:*), S
C
C     DESCENTE ET/OU REMONTEE?
      IF( NIVEAU .EQ. 0 ) GOTO 10
C
      print *
      print *,'AVANT la DESCENTE de LU x=b:'
      print 10000,(i,b0(i),i=1,20)
10000 FORMAT(5(i4,':',D25.17))
C
C     LA DESCENTE
C     -----------
      LPDIAGI = 0
      DO I = 1, NTDL
         LPDIAGI1 = LPDIAG(I)
         LH = ( LPDIAGI1 - LPDIAGI ) / 2
         IB = I -LH -1
         S  = 0D0
         DO K = 1, LH
            S = S + A(LPDIAGI+K) * B(IB+K)
         ENDDO
         B(I) = B0(I) - S
         LPDIAGI = LPDIAGI1
      ENDDO
C
C     REMONTEE OU PAS?
      IF( NIVEAU .EQ. 1 ) RETURN
C
      print *
      print *,'APRES la DESCENTE de LU x=b:'
      print 10000,(i,b(i),i=1,20)
C
C     LA REMONTEE
C     -----------
 10   LPDIAGI1 = LPDIAG(NTDL)
      DO I = NTDL, 1, -1
         LPDIAGI = LPDIAG(I-1)
         LH = ( LPDIAGI1 - LPDIAGI ) / 2
C
C        LA SOLUTION I
         B(I) = B(I) / A(LPDIAGI1)
C
C        CONTRIBUTION DE LA SOLUTION B(I) AUX LIGNES DE B QUI PRECEDENT
         DO K=1,LH
            B(I-K) = B(I-K) - B(I) * A(LPDIAGI1-K)
         ENDDO
         LPDIAGI1 = LPDIAGI
      ENDDO
C
      print *
      print *,'APRES la REMONTEE de LU x=b:'
      print 10000,(i,b(i),i=1,20)
      print *
C
      RETURN
      END
