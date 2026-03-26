      SUBROUTINE CRMC1D( MUDL, A0, NTDL, EPS, NENTRE,
     %                   A,    NRETOU )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: FACTORISER UNE MATRICE SYMETRIQUE INVERSIBLE SOUS LA
C ---- FORME  A = L * D * TL  DITE DE CROUT OU GAUSS SYMETRIQUE
C      L EST TRIANGULAIRE INFERIEURE A DIAGONALE UNITE
C      D EST DIAGONALE. A0 PUIS D,L SONT STOCKEES SOUS FORME PROFIL
C      C-A-D PAR LIGNES DU 1-ER COEFFICIENT NON NUL AU COEFFICIENT DIAGONAL
C      A0 ET A=(L,D) EN SORTIE,PEUVENT ETRE CONFONDUES EN ENTREE
C
C ENTREES:
C --------
C MUDL   : POINTEUR SUR LE COEFFICIENT DIAGONAL DE D.L
C          MUDL( 1 ) = 0
C          MUDL( 2 ) = NTDL SI LA MATRICE EST DIAGONALE (PAS DE FACTORISATION A)
C          MUDL(I+1) = ADRESSE DANS A DU I-EME COEFFICIENT DIAGONAL
C                      SI LA MATRICE A0 N'EST PAS DIAGONALE
C A0     : MATRICE A FACTORISER L * D * TL STOCKAGE PROFIL
C NTDL   : ORDRE DE LA MATRICE A
C EPS    : SEUIL AU DESSOUS DUQUEL LA FACTORISATION EST INCORRECTE
C NENTRE : =0 RETOUR AU PROGRAMME APPELANT SI ABS(PIVOT)<EPS
C          =1 LES CALCULS SE POURSUIVENT SAUF SI PIVOT=0
C
C SORTIES:
C --------
C A      : MATRICE FACTORISEE L ET D SI A0 EST NON DIAGONALE
C          MATRICE EGALE A A0 SI DIAGONALE > 0 ( MUDL(2)=NTDL )
C NRETOU : 0 SI AUCUN PIVOT<EPS
C          1 SI AU MOINS UN PIVOT<EPS
C          ON ENTEND PAR PIVOT (A0(IDIAGONAL)-SA)/A0(IDIAGONAL)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET      ANALYSE NUMERIQUE UPMC PARIS    AOUT 1998
C23456---------------------------------------------------------------012
      INTEGER           MUDL(1+NTDL), NTDL, NENTRE, NRETOU
      DOUBLE PRECISION  A0(1:*), A(1:*)
      REAL              EPS
      DOUBLE PRECISION  SA
C
      NRETOU = 0
C
      IF( MUDL(2) .EQ. NTDL ) THEN
C
C        MATRICE DIAGONALE
C        =================
         DO IP=1,NTDL
            A(IP) = A0(IP)
            IF( A0(IP) .LE. EPS ) THEN
               NRETOU = 1
               IF( NENTRE .NE. 1 ) RETURN
            ENDIF
         ENDDO
C
      ELSE
C
C        MATRICE NON DIAGONALE
C        =====================
10000    FORMAT(' crmc1d: [L] [D] t[L] line',I10)
         DO JP=1,NTDL
            IF( MOD( JP, 100 000 ) .EQ. 0 ) PRINT 10000, JP
            SA   = 0.D0
            JP1  = JP - 1
            JA0  = MUDL(JP)
            IAU2 = MUDL(JP+1) - JA0 - 1
            JPMI = JP - IAU2
C
C           CALCUL DES COEFFICIENTS NON DIAGONAUX
            DO IP=JPMI,JP1
               IA0  = MUDL(IP)
               IAU3 = MUDL(IP+1) - IA0 - 1
               IPMI = IP-IAU3
               IF( IPMI .GE. JPMI ) THEN
                  IA = IA0
                  JA = JA0 + IPMI - JPMI
                  KP = IPMI
               ELSE
                  IA = IA0 + JPMI - IPMI
                  JA = JA0
                  KP = JPMI
               ENDIF
               IAU = IP - KP
               SA  = 0.D0
               DO KD=1,IAU
                  SA = SA + A(IA+KD) * A(MUDL(KP+KD)) * A(JA+KD)
               ENDDO
               IAU4 = IAU + 1
               JA   = JA  + IAU4
               A(JA) = ( A0(JA) - SA ) / A(IA+IAU4)
            ENDDO
C
C           LE COEFFICIENT DIAGONAL
            SA = 0.D0
            DO KD=1,IAU2
               SA = SA + A(JA0+KD)**2 * A(MUDL(JPMI+KD))
            ENDDO
C
            JA = JA0 + IAU2 + 1
            SA = A0(JA) - SA
C
            IF( ABS(SA) .LE. EPS*ABS(A0(JA)) ) THEN
C              ON EST AU DESSOUS DU SEUIL DE PRECISION
               NRETOU=1
               IF( NENTRE .NE. 1 ) RETURN
            ENDIF
C
C           VALEUR DU COEFFICIENT DIAGONAL JP DE LA MATRICE D
            A(JA) = SA
         ENDDO
      ENDIF
C
      PRINT 19999,NTDL
19999 FORMAT(' crmc1d: [L] [D] t[L] line',I10,' ENDED')

      RETURN
      END
