      SUBROUTINE FAGAPR( LPDIAG, A0, NTDL, EPS, NENTRE,
     %                   A, NRETOU )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     FACTORISATION DE GAUSS A=L.U POUR UNE MATRICE INVERSIBLE
C ----     NON SYMETRIQUE MAIS A PROFIL SYMETRIQUE
C          A(I,J) NON NUL => A(J,I) NON NUL
C          L EST UNE MATRICE TRIANGULAIRE INFERIEURE A DIAGONALE UNITE
C          U EST UNE MATRICE TRIANGULAIRE SUPERIEURE
C          L ET U SONT RANGEES DANS A SOUS FORME PROFIL
C
C ENTREES:
C --------
C LPDIAG : LPDIAG(0)=0, LPDIAG(I)=ADRESSE DU I-EME COEFFICIENT DIAGONAL
C A0     : MATRICE AVANT FACTORISATION AVEC LE RANGEMENT
C          1 3    11   |    L PARTIE TRIANGULAIRE INFERIEURE
C          2 4  6 12   |    U PARTIE TRIANGULAIRE SUPERIEURE+DIAGONALE
C            5  7 13   V
C          8 9 10 14
C            --->
C
C NTDL   : ORDRE DES MATRICES A0 ET A
C EPS    : PRECISION DE LA FACTORISATION
C NENTRE : =0 RETOUR AU PROGRAMME DES LA DETECTION D'UN |PIVOT| < EPS
C          =1 LES CALCULS SE POURSUIVENT DANS TOUS LES CAS
C
C SORTIES:
C --------
C A      : MATRICE L ET  U  APRES FACTORISATION
C NRETOU : =0 SI AUCUN PIVOT (A0(IDIAGONAL)-SA)/A0(IDIAGONAL)<EPS
C          >0 NOMBRE DE PIVOTS SOUS LE SEUIL DE PRECISION EPS
C
C REMARQUE:
C ---------
C LES MATRICES A0 ET A PEUVENT ETRE CONFONDUES A L APPEL C-A-D
C QUE LA FACTORISATION PEUT SE FAIRE SUR LA MATRICE ELLE-MEME
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray  Avril 2010
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      DOUBLE PRECISION  A0(*), A(*)
      INTEGER           LPDIAG(0:*)
      DOUBLE PRECISION  AA, SL, SU, SDU, DAU, DABS
      COMMON /UNITES/   LECTEU, IMPRIM, NUNITE(30)
C
C
C     AFFICHAGE DU PROFIL DE LA MATRICE A AVANT FACTORISATION
C     -------------------------------------------------------
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)'Factorisation LU avec ntdl=',ntdl,
     %   ' Code Entree=', nentre,' PRECISION EPS=',eps
         WRITE(IMPRIM,*)'La matrice A AVANT factorisation'
      ELSE
         WRITE(IMPRIM,*)'LU Factorization with',ntdl,
     %   ' dof  In Code=', nentre,' PRECISION EPS=',eps
         WRITE(IMPRIM,*)'The A matrice BEFORE factorization'
      ENDIF
      CALL AFPROF( 20, LPDIAG, A )
C
      NRETOU  = 0
      LPDIAGI = 0
      DO I=1,NTDL
C
C        TRAITEMENT DE LA I-EME LIGNE DE L ET COLONNE DE U   (A=L*U)
         LPDIAGI1 = LPDIAG(I)
         IH = ( LPDIAGI1 - LPDIAGI ) / 2
C        IH = NOMBRE DE COEFFICIENTS DE LA I-EME LIGNE (DIAGONALE EXCLUE)
         SDU = 0.D0
         IF( IH .GT. 0 ) THEN
C
C           IL EXISTE DES COEFFICIENTS NON DIAGONAUX
            IMI = I - IH
C           IMI = N0 DE LA 1-ERE COLONNE NON NULLE DE LA I-EME LIGNE
            IAL0 = LPDIAGI
            IAU0 = LPDIAGI + IH
            IMI1 = IMI-1
            LPDIAGJ = LPDIAG(IMI1)
C
C           BOUCLE SUR LES COLONNES  DE LA I-EME  LIGNE DE L
C           BOUCLE SUR LES LIGNES    DE LA I-EME  COLONNE DE U
            DO JJ = 1, IH
               J = IMI1 + JJ
               LPDIAGJ1 = LPDIAG(J)
               JH  = ( LPDIAGJ1 - LPDIAGJ ) / 2
               JMI = J - JH
               JAL0 = LPDIAGJ
               JAU0 = LPDIAGJ + JH
C              JAL0 POINTE EN TETE DE LA J-EME LIGNE
C              JAU0 POINTE EN TETE DE LA J-EME COLONNE
               IS = IMI - JMI
               IF (IS .LE. 0 ) THEN
                  MA  = JMI
                  IAL = IAL0-IS
                  JAL = JAL0
                  IAU = IAU0-IS
                  JAU = JAU0
               ELSE
                  MA  = IMI
                  IAL = IAL0
                  JAL = JAL0+IS
                  IAU = IAU0
                  JAU = JAU0+IS
               ENDIF
C
               SL = 0.D0
               SU = 0.D0
               K1 = J - MA
               DO K = 1, K1
C
C                 L(I,J)
                  SL = SL + A(IAL+K) * A(JAU+K)
C
C                 U(J,I)
                  SU = SU + A(JAL+K) * A(IAU+K)
C
               ENDDO
C
               IAL = IAL +K1 +1
               IAU = IAU +K1 +1
               A(IAL) = (A0(IAL) - SL) / A(LPDIAGJ1)
               A(IAU) =  A0(IAU) - SU
C
               LPDIAGJ = LPDIAGJ1
            ENDDO
C
C           TRAITEMENT DU I-EME COEFFICIENT DIAGONAL
            DO K = 1,IH
               SDU = SDU + A(IAL0+K) * A(IAU0+K)
            ENDDO
C
         ENDIF
C
         AA  = A0( LPDIAGI1 )
         DAU = AA - SDU
C
C        TEST SUR LA PRECISION DU PIVOT
         IF( DABS(DAU) .LE. EPS*DABS(AA) ) THEN
C
C           COEFFICIENT AU DESSOUS DE LA PRECISION => PB?
            IF( LANGAG .EQ. 0 ) THEN
               WRITE (IMPRIM,10000) I, EPS, DAU, AA
            ELSE
               WRITE (IMPRIM,20000) I, EPS, DAU, AA
            ENDIF
            NRETOU = NRETOU + 1
            IF( NENTRE .EQ. 0 ) RETURN
C
         ENDIF
C
C        U DIAGONAL
         A(LPDIAGI1) = DAU
C
         LPDIAGI = LPDIAGI1
      ENDDO
C
10000 FORMAT(120('%')/'FACTORISATION A VERIFIER: COEFFICIENT DIAGONAL',
     %I10,' EST EN RELATIF AU DESSOUS DE LA PRECISION'/
     %'EPS=',G15.7,' DiagLU=',G15.7,' DiagA0= ',G15.7/120('%'))
C
20000 FORMAT(120('%')/'VERIFY the FACTORIZATION: DIAGONAL COEFFICIENT',
     %I10,' IS UNDER the PRECISION'/
     %'EPS=',G15.7,' DiagLU=',G15.7,' DiagA0= ',G15.7/120('%'))
C
C     AFFICHAGE DU PROFIL DE LA MATRICE L U (BLOC (20,20))
C     -------------------------------------
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
       WRITE(IMPRIM,*)'La matrice LU APRES factorisation  Code Retour=',
     %                 NRETOU
      ELSE
         WRITE(IMPRIM,*)'The LU matrice AFTER factorization  Out Code=',
     %                   NRETOU
      ENDIF
      CALL AFPROF( 20, LPDIAG, A )
C
      RETURN
      END
