      SUBROUTINE ASMEPC( NBDL,   NODL,
     %                   NCODAE, AES,    AE,
     %                   NCODAG, LPDIAG, AG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ASSEMBLER UNE MATRICE ELEMENTAIRE EN LA MATRICE GLOBALE DE
C -----    PROFIL SYMETRIQUE OU NON TOUTE EN MEMOIRE CENTRALE
C          (DIFFERENCE AVEC LE SP ASABPC PAS DE VECTEUR ELEMENTAIRE)

C ENTREES:
C --------
C NBDL   : NOMBRE DE       DEGRES DE LIBERTE DE L ELEMENT FINI
C NODL   : NUMERO DES NBDL DEGRES DE LIBERTE DE L ELEMENT FINI
C NCODAE : CODE DE STOCKAGE DE LA MATRICE ELEMENTAIRE
C           0 : MATRICE DIAGONALE
C          -1 : MATRICE NON SYMETRIQUE
C           1 : MATRICE SYMETRIQUE
C AES    : MATRICE ELEMENTAIRE SYMETRIQUE ( NBDL * (NBDL+1) / 2 )
C          OU DIAGONALE (NBDL)
C AE     : MATRICE ELEMENTAIRE PLEINE ( NBDL, NBDL )
C NCODAG : CODE DE STOCKAGE DE LA MATRICE GLOBALE PROFIL
C           0 : MATRICE DIAGONALE
C          -1 : MATRICE NON SYMETRIQUE
C           1 : MATRICE SYMETRIQUE
C LPDIAG : TABLEAU POINTEUR SUR LES COEFFICIENTS DIAGONAUX

C MODIFIE:
C --------
C AG     : MATRICE PROFIL TOUTE EN MEMOIRE CENTRALE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A. PERRONNET LABORATOIRE D'ANALYSE NUMERIQUE PARIS   MAI 1989
C ......................................................................
      INTEGER           NODL(NBDL), LPDIAG(0:*)
      DOUBLE PRECISION  AES(1:*), AE(NBDL,NBDL), AG(1:*)

      IF( NCODAE .EQ. 0 .AND. NCODAG .EQ. 0 ) THEN

C        LES MATRICES ELEMENTAIRE ET GLOBALE SONT DIAGONALES
C        ===================================================
         DO I=1,NBDL
C           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IEG = NODL( I )
C           ASSEMBLAGE DE AE( IEG )
            AG( IEG ) = AG( IEG ) + AES( I )
         ENDDO
         GOTO 9999

      ELSE IF( NCODAE .GT. 0 .AND. NCODAG .GT. 0 ) THEN

C        LES MATRICES ELEMENTAIRE ET GLOBALE SONT SYMETRIQUES
C        ====================================================
         L = 0
         DO I=1,NBDL
C           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IEG = NODL( I )

C           ASSEMBLAGE DE AES DANS AG
C           LA BOUCLE SUR LES COLONNES DE LA MATRICE ELEMENTAIRE
            DO J=1,I
               JEG = NODL( J )
               IF( JEG .LT. IEG ) THEN
C                 JEG <  IEG
                  KEG = LPDIAG( IEG ) - IEG + JEG
               ELSE
C                 JEG >= IEG
                  KEG = LPDIAG( JEG ) - JEG + IEG
               ENDIF
C              ASSEMBLAGE DE AE( I , J )
               L         = L + 1
               AG( KEG ) = AG( KEG ) + AES( L )
            ENDDO
         ENDDO
         GOTO 9999

      ELSE IF( NCODAE .LT. 0 .AND. NCODAG .LT. 0 ) THEN

C        LES MATRICES ELEMENTAIRE ET GLOBALE SONT NON SYMETRIQUES
C        ========================================================
         DO I=1,NBDL

C           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IEG = NODL( I )

C           LA BOUCLE SUR LES COLONNES DE LA MATRICE ELEMENTAIRE
            DO J=1,NBDL

C              LE NO GLOBAL DU J-EME DEGRE DE LIBERTE LOCAL
               JEG = NODL( J )

C              ADRESSE DE AE(I,J) DANS LA MATRICE GLOBALE DE AG(IEG,JEG)
               IF( JEG .LT. IEG ) THEN
C                 JEG <  IEG
                  KEG = (LPDIAG(IEG)+LPDIAG(IEG-1)+1)/2 - IEG + JEG
               ELSE
C                 JEG >= IEG
                  KEG = LPDIAG(JEG) - JEG + IEG
               ENDIF

C              ASSEMBLAGE DE AE(I,J)
               AG( KEG ) = AG( KEG ) + AE( I , J )
            ENDDO
         ENDDO
         GOTO 9999

      ELSE IF( NCODAE .EQ. 0 .AND. NCODAG .NE. 0 ) THEN

C        LA MATRICE ELEMENTAIRE EST DIAGONALE MAIS
C        LA MATRICE GLOBALE EST SYMETRIQUE OU NON SYMETRIQUE
C        ===================================================
         DO I=1,NBDL
C           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IEG = NODL( I )
C           SA POSITION SUR LA DIAGONALE DE LA MATRICE GLOBALE PROFIL
            KEG = LPDIAG( IEG )
C           ASSEMBLAGE DE AES( IEG )
            AG( KEG ) = AG( KEG ) + AES( I )
         ENDDO

      ENDIF

 9999 RETURN
      END
