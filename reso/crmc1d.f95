      SUBROUTINE CRMC1D( MUDL, A0, NTDL, EPS, NENTRE, &
                         A,    NRETOU )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT: FACTORISER UNE MATRICE SYMETRIQUE INVERSIBLE SOUS LA
! ---- FORME  A = L * D * TL  DITE DE CROUT OU GAUSS SYMETRIQUE
!      L EST TRIANGULAIRE INFERIEURE A DIAGONALE UNITE
!      D EST DIAGONALE. A0 PUIS D,L SONT STOCKEES SOUS FORME PROFIL
!      C-A-D PAR LIGNES DU 1-ER COEFFICIENT NON NUL AU COEFFICIENT DIAGONAL
!      A0 ET A=(L,D) EN SORTIE,PEUVENT ETRE CONFONDUES EN ENTREE

! ENTREES:
! --------
! MUDL   : POINTEUR SUR LE COEFFICIENT DIAGONAL DE D.L
!          MUDL( 1 ) = 0
!          MUDL( 2 ) = NTDL SI LA MATRICE EST DIAGONALE (PAS DE FACTORISATION A)
!          MUDL(I+1) = ADRESSE DANS A DU I-EME COEFFICIENT DIAGONAL
!                      SI LA MATRICE A0 N'EST PAS DIAGONALE
! A0     : MATRICE A FACTORISER L * D * TL STOCKAGE PROFIL
! NTDL   : ORDRE DE LA MATRICE A
! EPS    : SEUIL AU DESSOUS DUQUEL LA FACTORISATION EST INCORRECTE
! NENTRE : =0 RETOUR AU PROGRAMME APPELANT SI ABS(PIVOT)<EPS
!          =1 LES CALCULS SE POURSUIVENT SAUF SI PIVOT=0

! SORTIES:
! --------
! A      : MATRICE FACTORISEE L ET D SI A0 EST NON DIAGONALE
!          MATRICE EGALE A A0 SI DIAGONALE > 0 ( MUDL(2)=NTDL )
! NRETOU : 0 SI AUCUN PIVOT<EPS
!          1 SI AU MOINS UN PIVOT<EPS
!          ON ENTEND PAR PIVOT (A0(IDIAGONAL)-SA)/A0(IDIAGONAL)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR: ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         AOUT 1998
! MODIFS: ALAIN PERRONNET Saint Pierre du Perray              AVRIL 2023
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
   include"./incl/threads.inc"
      INTEGER           MUDL(1:*), NTDL, NENTRE, NRETOU
      DOUBLE PRECISION  A0(1:*), A(1:*)
      REAL              EPS
      DOUBLE PRECISION  SA

      NRETOU = 0

      IF( MUDL(2) .EQ. NTDL ) THEN

!        MATRICE DIAGONALE
!        =================
         DO IP=1,NTDL
            A(IP) = A0(IP)
            IF( A0(IP) .LE. EPS ) THEN
               NRETOU = 1
               IF( NENTRE .NE. 1 ) RETURN
            ENDIF
         ENDDO

      ELSE

!        MATRICE NON DIAGONALE
!        =====================
10000    FORMAT(' crmc1d.f95: start [L] [D] t[L] of lines',I12)
10001    FORMAT(' crmc1d.f95:       [L] [D] t[L] at line ',I12)
         PRINT 10000,NTDL
         DO JP=1,NTDL
            IF( MOD( JP, 100000 ) .EQ. 0 ) PRINT 10001, JP
            SA   = 0.D0
            JP1  = JP - 1
            JA0  = MUDL(JP)
            IAU2 = MUDL(JP+1) - JA0 - 1
            JPMI = JP - IAU2

!           CALCUL DES COEFFICIENTS NON DIAGONAUX
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
!///////////////////////////////////////////////////////////////////////
!EXECUTION AVEC NBTHREADS
!$OMP PARALLEL PRIVATE(KD) SHARED(MUDL,A,IAU,IA,KP,JA)
!$OMP DO REDUCTION(+:SA)
               DO KD=1,IAU
                  SA = SA + A(IA+KD) * A(MUDL(KP+KD)) * A(JA+KD)
               ENDDO
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////
               IAU4 = IAU + 1
               JA   = JA  + IAU4
               A(JA) = ( A0(JA) - SA ) / A(IA+IAU4)
            ENDDO

!           LE COEFFICIENT DIAGONAL
            SA = 0.D0
!///////////////////////////////////////////////////////////////////////
!EXECUTION AVEC NBTHREADS
!$OMP PARALLEL PRIVATE(KD) SHARED(MUDL,A,IAU2,JA0,JPMI)
!$OMP DO REDUCTION(+:SA)
            DO KD=1,IAU2
               SA = SA + A(JA0+KD)**2 * A(MUDL(JPMI+KD))
            ENDDO
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////
            JA = JA0 + IAU2 + 1
            SA = A0(JA) - SA

            IF( ABS(SA) .LE. EPS*ABS(A0(JA)) ) THEN
!              ON EST AU DESSOUS DU SEUIL DE PRECISION
               NRETOU=1
               IF( NENTRE .NE. 1 ) RETURN
            ENDIF

!           VALEUR DU COEFFICIENT DIAGONAL JP DE LA MATRICE D
            A(JA) = SA
         ENDDO
      ENDIF

      PRINT 19999,NTDL
19999 FORMAT(' crmc1d.f95: End   [L] [D] t[L] of lines',I12)

      RETURN
      END
