      SUBROUTINE ASMEPC( NBDL,   NODL, &
                         NCODAE, AES,    AE, &
                         NCODAG, LPDIAG, AG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT :    ASSEMBLER UNE MATRICE ELEMENTAIRE EN LA MATRICE GLOBALE DE
! -----    PROFIL SYMETRIQUE OU NON TOUTE EN MEMOIRE CENTRALE
!          (DIFFERENCE AVEC LE SP ASABPC PAS DE VECTEUR ELEMENTAIRE)

! ENTREES:
! --------
! NBDL   : NOMBRE DE       DEGRES DE LIBERTE DE L ELEMENT FINI
! NODL   : NUMERO DES NBDL DEGRES DE LIBERTE DE L ELEMENT FINI
! NCODAE : CODE DE STOCKAGE DE LA MATRICE ELEMENTAIRE
!           0 : MATRICE DIAGONALE
!          -1 : MATRICE NON SYMETRIQUE
!           1 : MATRICE SYMETRIQUE
! AES    : MATRICE ELEMENTAIRE SYMETRIQUE ( NBDL * (NBDL+1) / 2 )
!          OU DIAGONALE (NBDL)
! AE     : MATRICE ELEMENTAIRE PLEINE ( NBDL, NBDL )
! NCODAG : CODE DE STOCKAGE DE LA MATRICE GLOBALE PROFIL
!           0 : MATRICE DIAGONALE
!          -1 : MATRICE NON SYMETRIQUE
!           1 : MATRICE SYMETRIQUE
! LPDIAG : TABLEAU POINTEUR SUR LES COEFFICIENTS DIAGONAUX

! MODIFIE:
! --------
! AG     : MATRICE PROFIL TOUTE EN MEMOIRE CENTRALE
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : A. PERRONNET LABORATOIRE D'ANALYSE NUMERIQUE PARIS   MAI 1989
! ......................................................................
!$    use OMP_LIB
      INTEGER           NODL(NBDL), LPDIAG(0:*)
      DOUBLE PRECISION  AES(1:*), AE(NBDL,NBDL), AG(1:*)

      IF( NCODAE .EQ. 0 .AND. NCODAG .EQ. 0 ) THEN

!        LA MATRICE ELEMENTAIRE ET GLOBALE SONT DIAGONALES
!        =================================================
         DO I=1,NBDL

!           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IEG = NODL( I )
!           ASSEMBLAGE DE AE( IEG )
!$OMP ATOMIC
            AG( IEG ) = AG( IEG ) + AES( I )

         ENDDO
         GOTO 9999

      ELSE IF( NCODAE .GT. 0 .AND. NCODAG .GT. 0 ) THEN

!        LA MATRICE ELEMENTAIRE ET GLOBALE SONT SYMETRIQUES
!        ==================================================
         L = 0
         DO I=1,NBDL

!           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IEG = NODL( I )

!           ASSEMBLAGE DE AES DANS AG
!           LA BOUCLE SUR LES COLONNES DE LA MATRICE ELEMENTAIRE
            DO J=1,I
               JEG = NODL( J )
               IF( JEG .LT. IEG ) THEN
!                 JEG <  IEG
                  KEG = LPDIAG( IEG ) - IEG + JEG
               ELSE
!                 JEG >= IEG
                  KEG = LPDIAG( JEG ) - JEG + IEG
               ENDIF
!              ASSEMBLAGE DE AE( I , J )
               L = L + 1
!$OMP ATOMIC
               AG( KEG ) = AG( KEG ) + AES( L )
            ENDDO

         ENDDO
         GOTO 9999

      ELSE IF( NCODAE .LT. 0 .AND. NCODAG .LT. 0 ) THEN

!        LA MATRICES ELEMENTAIRE ET GLOBALE SONT NON SYMETRIQUES
!        =======================================================
         DO I=1,NBDL

!           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IEG = NODL( I )

!           LA BOUCLE SUR LES COLONNES DE LA MATRICE ELEMENTAIRE
            DO J=1,NBDL

!              LE NO GLOBAL DU J-EME DEGRE DE LIBERTE LOCAL
               JEG = NODL( J )

!              ADRESSE DE AE(I,J) DANS LA MATRICE GLOBALE DE AG(IEG,JEG)
               IF( JEG .LT. IEG ) THEN
!                 JEG <  IEG
                  KEG = (LPDIAG(IEG)+LPDIAG(IEG-1)+1)/2 - IEG + JEG
               ELSE
!                 JEG >= IEG
                  KEG = LPDIAG(JEG) - JEG + IEG
               ENDIF

!              ASSEMBLAGE DE AE(I,J)
!$OMP ATOMIC
               AG( KEG ) = AG( KEG ) + AE( I , J )

            ENDDO

         ENDDO
         GOTO 9999

      ELSE IF( NCODAE .EQ. 0 .AND. NCODAG .NE. 0 ) THEN

!        LA MATRICE ELEMENTAIRE EST DIAGONALE MAIS
!        LA MATRICE GLOBALE EST SYMETRIQUE OU NON SYMETRIQUE
!        ===================================================
         DO I=1,NBDL

!           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IEG = NODL( I )
!           SA POSITION SUR LA DIAGONALE DE LA MATRICE GLOBALE PROFIL
            KEG = LPDIAG( IEG )
!           ASSEMBLAGE DE AES( IEG )
!$OMP ATOMIC
            AG( KEG ) = AG( KEG ) + AES( I )

         ENDDO

      ENDIF

9999  RETURN
      END
