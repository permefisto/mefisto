SUBROUTINE AGVITTH( Rho,    DtMhu,  XYZNOE, &
                    NUTYEL, NBNOEF, NBELEM, NUNOEF, &
                    NORESO, NBNOE,  NBCVG,  LP2LIGN, LP2COLO, VG )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CONSTRUIRE ET ASSEMBLER LA MATRICE GLOBALE DU PROBLEME
! ----     VG=( Rho - dt Mhu Laplacien ) = Integrale Rho P2 P2
!                                        + dtMhu grad P2 grad P2 dX
!          avec INTEGRATION EXACTE SUR L'ELEMENT FINI
!
! ENTREES:
! --------
! Rho    : COEFFICIENT DE LA MATRICE DE MASSE      P2  P2
! DtMhu  : COEFFICIENT DE LA MATRICE DE VISCOSITE DP2 DP2
! XYZNOE : XYZNOE(3,NBNOE)  3 COORDONNEES DES NOEUDS DU MAILLAGE
! NUTYEL : NUMERO DU TYPE D'EF ( 15 TRIANGLE, 20 TETRAEDRE TAYLOR-HOOD )
! NBNOEF : NOMBRE DE NOEUDS D'UN EF DE CE TYPE
!          ( 6 POUR LE TRIANGLE, 10 POUR LE TETRAEDRE )
! NBELEM : NOMBRE D'EF DU MAILLAGE
! NUNOEF : NUNOEF(NBELEM,NBNOEF) NO DES 6 NOEUDS DES NBELEM EF

! NORESO : CODE RESOLUTION DES SYSTEMES LINEAIRES AVEC LA MATRICE VG
!          =1 FACTORISATION COMPLETE DE CROUT SUR MATRICE PROFIL VG
!          =2 GRADIENT CONJUGUE AVEC UN STOCKAGE MORSE DE VG
! NBNOE  : NOMBRE DE NOEUDS DU MAILLAGE (SOMMETS ET MILIEUX DES ARETES)
! NBCVG  : NOMBRE DE COEFFICIENTS DE VG
! LP2LIGN: POINTEURS SUR LES COEFFICIENTS DIAGONAUX DE VG
! LP2COLO: NUMERO DE COLONNE DES COEFFICIENTS DE LA MATRICE MORSE VG
!
! SORTIE :
! --------
! VG     : MATRICE GLOBALE ASSEMBLEE A PARTIR DES MATRICES ELEMENTAIRES
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Octobre 2012
! MODIFS : ALAIN PERRONNET             St Pierre du Perray    Avril 2023
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
      IMPLICIT NONE
      include"./incl/threads.inc"
      include"./incl/donflu.inc95"
      include"./incl/a___npef.inc95"

      INTEGER           NORESO, NBNOE, NBCVG, LP2LIGN(0:NBNOE), &
                        LP2COLO(NBCVG), NUTYEL, NBNOEF, NBELEM, &
                        NUNOEF(NBELEM,NBNOEF)
      DOUBLE PRECISION  VG(NBCVG), Rho, DtMhu
      REAL              XYZNOE(3,NBNOE)

      INTEGER           NUELEM, NONOEF(10), &
                        I, J, L, IEG, JEG, KEG
      DOUBLE PRECISION  AE(55), S

!     MISE A ZERO DE LA MATRICE VG
      CALL AZEROD( NBCVG, VG )

!     BOUCLE SUR LES EF de TAYLOR-HOOD 2d ou 3d
!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( AE, S ) &
!$OMP PRIVATE( NUELEM, NONOEF, I, J, L, IEG, JEG, KEG )
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS )

      DO NUELEM = 1, NBELEM

!        NO DES NOEUDS P2 DE L'ELEMENT FINI NUELEM
         DO I = 1, NBNOEF
            NONOEF(I) = NUNOEF(NUELEM,I)
         ENDDO

         IF( NUTYEL .EQ. 15 ) THEN

!           TRIANGLE TAYLOR-HOOD en 2d
!           CALCUL DE LA MATRICE DE VISCOSITE Ae DE LA VITESSE P2
            CALL FAE2P2P2( Rho, DtMhu, NONOEF, NBNOE, XYZNOE, AE )

         ELSE

!           TETRAEDRE TAYLOR-HOOD en 3d
!           CALCUL DE LA MATRICE DE VISCOSITE Ae DE LA VITESSE P2
            CALL FAE3P2P2( Rho, DtMhu, NONOEF, NBNOE, XYZNOE, AE )

         ENDIF

!!!        print*
!!!        print*,'agvitth.95: EF',NUELEM,' NONOEF=',(NONOEF(L),L=1,NBNOEF)
!!!        print*,'agvitth.95: AE=',(AE(L),L=1,21)

!        ASSEMBLAGE DE LA MATRICE ELEMENTAIRE DANS LA MATRICE GLOBALE
         IF( NORESO .EQ. 1 ) THEN

!           VG MATRICE SYMETRIQUE PROFIL
!!!            CALL ASMEPC( NBNOEF, NONOEF, 1, AE, AE, 1, &
!!!                         LP2LIGN, VG )
            L = 0
            DO I = 1, NBNOEF
!              LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
               IEG = NONOEF( I )
!              ASSEMBLAGE DE AE DANS VG
!              LA BOUCLE SUR LES COLONNES DE LA MATRICE ELEMENTAIRE
               DO J = 1, I
                  JEG = NONOEF( J )
                  IF( JEG .LT. IEG ) THEN
!                    JEG <  IEG
                     KEG = LP2LIGN( IEG ) - IEG + JEG
                  ELSE
!                    JEG >= IEG
                     KEG = LP2LIGN( JEG ) - JEG + IEG
                  ENDIF
!                 ASSEMBLAGE DE AE( I , J ) DANS VG( KEG )
                  L = L + 1
                  S = AE( L )
!$OMP ATOMIC
                  VG( KEG ) = VG( KEG ) + S
               ENDDO
            ENDDO

         ELSE

!           VG MATRICE SYMETRIQUE MORSE CONDENSEE
!           LES MATRICES ELEMENTAIRE ET GLOBALE SONT ICI SYMETRIQUES
            CALL ASMEGC( NBNOEF, NONOEF, 1, AE, AE, 1, &
                         LP2LIGN, LP2COLO, VG )

         ENDIF
!!!         call affvect( 'agvitth.f95: MATRICE VITESSE VG', 5,    VG )

      ENDDO     ! FIN DE LA BOUCLE DES EF POUR ASSEMBLER VG
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!!!       call affvect( 'agvitth.f95: MATRICE VITESSE VG', 5,    VG )
!!!       call afl1ve(  'agvitth.f95: MATRICE VITESSE VG', NBCVG, VG )

      RETURN
END SUBROUTINE AGVITTH
