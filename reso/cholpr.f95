SUBROUTINE CHOLPR( NTDL, NCODSA, LPDIAG, A0, &
                   A,    IERR )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT: FACTORISER UNE MATRICE SYMETRIQUE DEFINIE POSITIVE SOUS LA
! ---- FORME  A = L * TL DITE DE CHOLESKY
!      A ET L SONT  TRIANGULAIRES INFERIEURES STOCKEES SOUS FORME PROFIL
!      C-A-D PAR LIGNES DU 1-ER COEFFICIENT NON NUL AU COEFFICIENT
!      DIAGONAL. LE PROFIL EST LE MEME POUR A ET L
!      A0 ET A=L PEUVENT ETRE CONFONDUS EN ENTREE

!      VERSION REELLE DOUBLE PRECISION
!      VERSION MODIFIEE POUR PRENDRE EN COMPTE LE CAS A SINGULIERE

! ENTREES:
! --------
! NTDL   : ORDRE DE LA MATRICE A
! NCODSA : 0     SI MATRICE DIAGONALE
!          NON 0 SI MATRICE SYMETRIQUE NON DIAGONALE
! LPDIAG : SI NCODSA>0 ALORS POINTEUR SUR LE COEFFICIENT DIAGONAL DE D.L
!          LPDIAG(0)=0, LPDIAG(I)=ADRESSE DANS A DU I-EME COEF DIAG
!                                 DIAGONAL SI A NON DIAGONALE
! A0     : MATRICE A FACTORISER SOUS FORME L * TL AVEC STOCKAGE PROFIL

! SORTIES:
! --------
! A      : MATRICE FACTORISEE L POUR TOUTE VALEUR DE NCODSA
! IERR   : CODE D'ERREUR 0 SI PAS D'ERREUR, NON NUL SINON
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MAI  1989
! MODIFS : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS        MARS 1990
! MODIFS : PERRONNET ALAIN LJLL UPMC & St Pierre du Perray    AVRIL 2013
! 3456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      COMMON / EPSSSS / EPZERO, EPSXYZ
      DOUBLE PRECISION  SA, A0(*), A(*), PIVOT, PROSCD
      INTEGER           LPDIAG(0:NTDL)
      INTRINSIC         SQRT, ABS

10010 FORMAT(' ATTENTION - CHOLPR :',I8,'EME PIVOT=',G13.6,' < 0.', &
      /,' LE PIVOT PRECEDENT',G13.6,' EST PRIS A LA PLACE')
10020 FORMAT(/,' ATTENTION - CHOLPR : FACTORISATION INSTABLE')

!     SEUIL POUR DETECTER LES PIVOTS INCORRECTS
!     =========================================
      IERR = 0
      IF( A0(1) .LE. 0D0 ) THEN
         IERR = 3
         RETURN
      ENDIF

      PIVOT = SQRT( ABS(A0(1)) )
      IF( PIVOT .LT. EPZERO ) THEN
         PIVOT = 1.D0
      ENDIF

      IF( NCODSA .EQ. 0 ) GOTO 1000

!     MATRICE SYMETRIQUE NON DIAGONALE
!     ================================
      NBPIMA = 10
      NBPIVO = 0
!     I INDIQUE LES COLONNES IP
!     J INDIQUE LES LIGNES   JP
      DO JP = 1, NTDL
         JP1  = JP - 1
         JA0  = LPDIAG(JP1)
         NBCJ = LPDIAG(JP ) - JA0 - 1
         JPMI = JP - NBCJ

!        CALCUL DES COEFFICIENTS NON DIAGONAUX
         DO IP = JPMI, JP1

!           COEFFICIENT NON DIAGONAL L(JP,IP)
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

            SA  = PROSCD( A(IA+1), A(JA+1), NBC )
!!!            SA = 0D0
!!!            DO KD = 1, NBC
!!!               SA = SA + A(IA+KD) * A(JA+KD)
!!!            ENDDO

            NBC  = NBC + 1
            JA   = JA  + NBC
            A(JA) = ( A0(JA) - SA ) / A(IA+NBC)
!!!         WRITE(IMPRIM,*) 'L(',JP,',',IP,')=',A(JA)

         ENDDO

!        COEFFICIENT DIAGONAL L(JP,JP)
         SA = PROSCD( A(JA0+1), A(JA0+1), NBCJ )
!!!         SA = 0D0
!!!         DO KD = 1, NBCJ
!!!            SA = SA + A(JA0+KD)**2
!!!         ENDDO

         JA = JA0 + NBCJ + 1
         SA = A0(JA) - SA
         IF( SA .LE. 0.D0 ) THEN
!           COEFFICIENT DIAGONAL NEGATIF
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

         IF( SA - EPZERO * A0(JA) .LE. 0.D0 ) THEN
!           ON EST AU DESSOUS DU SEUIL DE PRECISION
            IERR = 1
            SA   = PIVOT
         ELSE
!           SAUVEGARDE DU PIVOT ACTUEL
            PIVOT = SA
         ENDIF
         A(JA) = SQRT( SA )

!!!      WRITE(IMPRIM,*) 'L(',JP,',',JP,')=',A(JA)
!!!      WRITE(IMPRIM,*)

      ENDDO
      RETURN

!     CHOLESKY SUR LA MATRICE DIAGONALE A
!     ===================================
 1000 DO IP = 1, NTDL

         IF( A0(IP) .LE. 0.D0 ) THEN

!           COEFFICIENT DIAGONAL INCORRECT
            IERR = 1
            WRITE(IMPRIM,10010) IP, A0(IP), PIVOT
            A(IP) = PIVOT

         ELSE

!           COEFFICIENT DIAGONAL CORRECT
            A(IP) = SQRT( A0(IP) )
            PIVOT = A(IP)

         ENDIF
      ENDDO

      RETURN
END SUBROUTINE CHOLPR
