      SUBROUTINE ASMEGC( NBDLEF, NODLEF, NCODAE, AES, AE, &
                         NCODAG, LPLIGN, LPCOLO, AG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT :   ASSEMBLER LA MATRICE ELEMENTAIRE EN LA MATRICE GLOBALE
! -----   MORSE SYMETRIQUE OU NON TOUTE EN MEMOIRE CENTRALE
!         (DIFFERENCE AVEC LE SP ASABGC PAS DE VECTEUR ELEMENTAIRE)

! ENTREES:
! --------
! NBDLEF : NOMBRE DE   DEGRES DE LIBERTE DE L ELEMENT FINI
! NODLEF : NO DES NBDLEF DEGRES DE LIBERTE DE L ELEMENT FINI
! NCODAE : CODE DE STOCKAGE DE LA MATRICE ELEMENTAIRE
!           0 : MATRICE DIAGONALE
!          -1 : MATRICE NON SYMETRIQUE
!           1 : MATRICE SYMETRIQUE
! AES    : MATRICE ELEMENTAIRE SYMETRIQUE (NBDLEF * (NBDLEF+1) / 2 )
! AE     : MATRICE ELEMENTAIRE  (NBDLEF,NBDLEF)
! NCODAG : CODE DE STOCKAGE DE LA MATRICE GLOBALE MORSE
!           0 : MATRICE DIAGONALE
!          -1 : MATRICE NON SYMETRIQUE
!           1 : MATRICE SYMETRIQUE
! LPLIGN : TABLEAU POINTEUR SUR LES COEFFICIENTS DIAGONAUX
! LPCOLO : TABLEAU DES NUMEROS DE COLONNES

! SORTIES:
! --------
! AG     : MATRICE MORSE TOUTE EN MEMOIRE CENTRALE
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR :   A. PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS   Mai 1989
! MODIF OMP: A. PERRONNET LJLL UPMC & St Pierre du Perray  Decembre 2012
! MODIFS:    A. PERRONNET           & St Pierre du Perray       Mai 2023
! ......................................................................
!$ use OMP_LIB
      INTEGER           NODLEF(1:NBDLEF), LPLIGN(0:*), LPCOLO(1:*)
      DOUBLE PRECISION  AES(1:*), AE(1:NBDLEF,1:NBDLEF), AG(1:*), S

      IF( NCODAE .EQ. 0 .AND. NCODAG .EQ. 0 ) THEN

!        LES MATRICES ELEMENTAIRE ET GLOBALE SONT DIAGONALES
!        ===================================================
         DO I=1,NBDLEF

!           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IGL = NODLEF( I )
!           ASSEMBLAGE DE AE( IGL )
!$OMP ATOMIC
            AG( IGL ) = AG( IGL ) + AES( I )

         ENDDO
         GOTO 9999

      ELSE IF( NCODAE .GT. 0 .AND. NCODAG .GT. 0 ) THEN

!        LES MATRICES ELEMENTAIRE ET GLOBALE SONT SYMETRIQUES
!        ====================================================
         L = 0
         DO 40 I=1,NBDLEF

!           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IGL = NODLEF( I )

!           ASSEMBLAGE DE AES DANS AG
!           LA BOUCLE SUR LES COLONNES DE LA MATRICE ELEMENTAIRE
            DO 30 J=1,I

!              LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
               JGL = NODLEF( J )

!              ADRESSES DANS LA MATRICE MORSE
               IF( JGL .LE. IGL ) THEN
                  NC1  = LPLIGN(IGL-1) + 1
                  NC2  = LPLIGN(IGL)
                  NCGL = JGL
               ELSE
                  NC1  = LPLIGN(JGL-1) + 1
                  NC2  = LPLIGN(JGL)
                  NCGL = IGL
               ENDIF

!              RECHERCHE DU NUMERO DE COLONNE NCGL PAR DICHOTOMIE
!              ENTRE NC1 et NC2
               NCI=NC1
               NCS=NC2

!              NCM NO COLONNE AU MILIEU
 25            NCM = (NCI + NCS) / 2

               IF( LPCOLO(NCM) .EQ. NCGL ) GOTO 29

               IF( LPCOLO(NCM) .GT. NCGL ) THEN
                  NCS=NCM
               ELSE
                  NCI=NCM
               ENDIF

               IF( NCI+1 .NE. NCS ) GOTO 25

!              NCS = NCI+1  A TRAITER A CAUSE DE LA DIVISION ENTIERE de NCM
               IF( LPCOLO(NCS) .EQ. NCGL ) THEN
                  NCM = NCS
                  GOTO 29
               ENDIF

               IF( LPCOLO(NCI) .EQ. NCGL ) THEN
                   NCM = NCI
                   GOTO 29
               ENDIF

             print*,'asmegc.f95: no ligne          =',IGL
             print*,'asmegc.f95: no colonne        =',JGL
             print*,'asmegc.f95: no colonne INCONNU=',NCGL
             print*,'asmegc.f95: lpcolo voisins=',(LPCOLO(k),k=NC1,NC2)

!              RECHERCHE BESTIALE DU NUMERO DE COLONNE NCGL
               DO NCM=NC1,NC2
                  IF( LPCOLO(NCM) .EQ. NCGL ) GOTO 29
               ENDDO

!              ASSEMBLAGE DE AES( I, J ) dans AG( IGL, JGL )
 29            L = L + 1
               S = AES( L )
!$OMP ATOMIC
               AG( NCM ) = AG( NCM ) + S

 30         ENDDO

 40      ENDDO
         GOTO 9999

      ELSE IF( NCODAE .LT. 0 .AND. NCODAG .LT. 0 ) THEN

!        LES MATRICES ELEMENTAIRE ET GLOBALE SONT NON SYMETRIQUES
!        ========================================================
         DO I=1,NBDLEF

!           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IGL = NODLEF( I )

!           ADRESSES DANS LA MATRICE MORSE
            NC1 = LPLIGN(IGL-1) + 1
            NC2 = LPLIGN(IGL)

!           LA BOUCLE SUR LES COLONNES DE LA MATRICE ELEMENTAIRE
            DO 130 J=1,NBDLEF

!              LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
               JGL = NODLEF( J )

!              RECHERCHE DU NUMERO DE COLONNE JGL PAR DICHOTOMIE
               NCI=NC1
               NCS=NC2

!              RECHERCHE DU NUMERO DE COLONNE NCGL PAR DICHOTOMIE
!              ENTRE NC1 et NC2
               NCI=NC1
               NCS=NC2

!              NCM NO COLONNE AU MILIEU
 125           NCM = (NCI + NCS) / 2

               IF( LPCOLO(NCM) .EQ. JGL ) GOTO 129

               IF( LPCOLO(NCM) .GT. JGL ) THEN
                  NCS=NCM
               ELSE
                  NCI=NCM
               ENDIF

               IF( NCI+1 .NE. NCS ) GOTO 125

!              NCS = NCI+1  A TRAITER A CAUSE DE LA DIVISION ENTIERE de NCM
               IF( LPCOLO(NCS) .EQ. JGL ) THEN
                  NCM = NCS
                  GOTO 129
               ENDIF

               IF( LPCOLO(NCI) .EQ. JGL ) THEN
                   NCM = NCI
                   GOTO 129
               ENDIF

             print*,'asmegc.f95: NO ligne          =',IGL
             print*,'asmegc.f95: NO colonne INCONNU=',JGL
             print*,'asmegc.f95: lpcolo voisins    =',(LPCOLO(k),k=NC1,NC2)

 129           S = AE( I , J )
!              ASSEMBLAGE DE AE(I,J) DANS AG(IGL,JGL)
!$OMP ATOMIC
               AG( NCM ) = AG( NCM ) + S

 130        ENDDO

         ENDDO
         GOTO 9999

      ELSE IF( NCODAE .EQ. 0 .AND. NCODAG .NE. 0 ) THEN

!        LA MATRICE ELEMENTAIRE EST DIAGONALE MAIS
!        LA MATRICE GLOBALE EST SYMETRIQUE OU NON SYMETRIQUE
!        ===================================================
         DO I=1,NBDLEF

!           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IGL = NODLEF( I )
!           SA POSITION SUR LA DIAGONALE DE LA MATRICE GLOBALE MORSE
            IGL = LPLIGN( IGL )
!           ASSEMBLAGE DE AES( IGL )
!$OMP ATOMIC
            AG( IGL ) = AG( IGL ) + AES( I )

         ENDDO

      ENDIF

 9999 RETURN
      END
