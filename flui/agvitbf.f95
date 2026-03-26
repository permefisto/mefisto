SUBROUTINE AGVITBF( DELTAT, Rho,    Mhu,    XYZNOE, &
                    NUTYEL, NBNOEF, NBELEM, NUNOEF, NOBARY, &
                    NORESO, NBSOM,  NBCVG,  LPLIGN, LPCOLO, &
                    VG,     AGGAUSS )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     CONSTRUIRE ET ASSEMBLER LA MATRICE GLOBALE DU PROBLEME
! ----     VG=( Rho - dt Mhu Laplacien ) = Integrale Rho P1B P1B
!                                        + deltat Mhu grad P1B grad P1B dX
!          avec INTEGRATION EXACTE SUR L'ELEMENT FINI
!
! ENTREES:
! --------
! DELTAT : PAS DE TEMPS DU SCHEMA
! Rho    : COEFFICIENT DE LA MATRICE DE MASSE      P1B  P1B
! Mhu    : COEFFICIENT DE LA MATRICE DE VISCOSITE DP1B DP1B
! XYZNOE : XYZNOE(3,NBSOM)  3 COORDONNEES DES SOMMETS DU MAILLAGE
! NUTYEL : NUMERO DU TYPE D'EF ( 13 TRIANGLE, 19 TETRAEDRE BREZZI-FORTIN )
! NBNOEF : NOMBRE DE SOMMETS D'UN EF DE CE TYPE
!          ( 3 POUR LE TRIANGLE, 4 POUR LE TETRAEDRE)
! NBELEM : NOMBRE D'EF DU MAILLAGE
! NUNOEF : NUNOEF(NBELEM,NBNOEF) NO DES SOMMETS DES NBELEM EF
! NOBARY : NUMERO DU DU DL BARYCENTRE = NBNOEF + 1 = NBSOMMET 1 EF + 1
!
! NORESO : CODE RESOLUTION DES SYSTEMES LINEAIRES AVEC LA MATRICE VG
!          1 FACTORISATION COMPLETE DE CHOLESKY ET MATRICE PROFIL VG
!          2 GRADIENT CONJUGUE AVEC UN STOCKAGE MORSE DE VG
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
! NBCVG  : NOMBRE DE COEFFICIENTS DE VG
!          (SI NORESO=1   NBCVG=0
!          (SI NORESO=2   NBCVG=NOMBRE DE COLONNES STOCKEES)
! LPLIGN : POINTEURS SUR LES COEFFICIENTS DIAGONAUX DE VG
! LPCOLO : NUMERO DE COLONNE DES COEFFICIENTS DE LA MATRICE MORSE VG
!
! SORTIE :
! --------
! VG     : MATRICE GLOBALE ASSEMBLEE
! AGGAUSS: COEFFICIENTS DE LA MATRICE ELEMENTAIRE ELIMINES PAR GAUSS
!          POUR LES NBELEM ELEMENTS FINIS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Octobre 2012
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
      IMPLICIT NONE
      include"./incl/threads.inc"
      INTEGER           NORESO, NBSOM, NBCVG, LPLIGN(0:NBSOM), &
                        LPCOLO(NBCVG), NUTYEL, NBNOEF, NBELEM, &
                        NUNOEF(NBELEM,NBNOEF)
      DOUBLE PRECISION  VG(NBCVG), AGGAUSS(NOBARY,NBELEM), &
                        DELTAT, Rho, Mhu, DtMhu
      REAL              XYZNOE(3,NBSOM)

      INTEGER           M, NUELEM, NONOEF(5), NOBARY, &
                        I, J, L, IEG, JEG, KEG
      DOUBLE PRECISION  AE(15), S

!     MISE A ZERO DE LA MATRICE VG
      CALL AZEROD( NBCVG, VG )

!     Coefficient de dt Mhu Integrale Grad P1B Grad P1B dx
      DtMhu = DELTAT * Mhu

!     BOUCLE SUR LES EF de BREZZI-FORTIN
!//////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( NUELEM, M, NONOEF, AE, S ) &
!$OMP PRIVATE( I, J, L, IEG, JEG, KEG )
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS )

!     LE STOCKAGE DES LIGNES ELIMINEES PAR GAUSS AUX BARYCENTRES
      DO NUELEM = 1, NBELEM

!        NO DES NOEUDS DE L'ELEMENT FINI NUELEM
         DO M = 1, NBNOEF
            NONOEF(M) = NUNOEF(NUELEM,M)
         ENDDO
         NONOEF(NOBARY) = NBSOM + NUELEM

!        CALCUL DE LA MATRICE avec INTEGRATION EXACTE SUR L'EF de
!        Integrale ( Rho P1B P1B + dt Mhu Grad P1B Grad P1B ) dx
!        ELIMINATION DE GAUSS DU DL VITESSE DU BARYCENTRE (1 COMPOSANTE)
!        SUR LA MATRICE ELEMENTAIRE
         IF( NUTYEL .EQ. 13 ) THEN
!           TRIANGLE BF 2D
            CALL AEBFV2( Rho, DtMhu, NONOEF, NBSOM, XYZNOE, &
                         AE,  AGGAUSS(1,NUELEM) )
         ELSE
!           TETRAEDRE BF 3D
            CALL AEBFV3( Rho, DtMhu, NONOEF, NBSOM, XYZNOE, &
                         AE,  AGGAUSS(1,NUELEM) )
         ENDIF

!        ASSEMBLAGE DE LA MATRICE ELEMENTAIRE DANS LA MATRICE GLOBALE
         IF( NORESO .EQ. 1 ) THEN

!           MATRICE SYMETRIQUE PROFIL
!!!            CALL ASMEPC( NBNOEF, NONOEF, 1, AE, AE, 1, LPLIGN, VG )
            L = 0
            DO I=1,NBNOEF
!              LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
               IEG = NONOEF( I )
!              ASSEMBLAGE DE AE DANS VG
!              LA BOUCLE SUR LES COLONNES DE LA MATRICE ELEMENTAIRE
               DO J=1,I
                  JEG = NONOEF( J )
                  IF( JEG .LT. IEG ) THEN
!                    JEG <  IEG
                     KEG = LPLIGN( IEG ) - IEG + JEG
                  ELSE
!                    JEG >= IEG
                  KEG = LPLIGN( JEG ) - JEG + IEG
                  ENDIF
!                 ASSEMBLAGE DE AE( I , J ) DANS VG( KEG )
                  L = L + 1
                  S = AE( L )
!$OMP ATOMIC
                  VG( KEG ) = VG( KEG ) + S
               ENDDO
            ENDDO

         ELSE

!           MATRICE SYMETRIQUE MORSE
!           LES MATRICES ELEMENTAIRE ET GLOBALE SONT ICI SYMETRIQUES
            CALL ASMEGC( NBNOEF, NONOEF, 1, AE, AE, 1, &
                         LPLIGN, LPCOLO, VG )
         ENDIF

      ENDDO     ! FIN DE LA BOUCLE DES EF POUR ASSEMBLER VG
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!!!      call affvect( 'agvitbf.f95: MATRICE VITESSE VG COMPOSANTE', 20, VG )
!!!      call affvect( 'agvitbf.f95: MATRICE VITESSE AG GAUSS', 20, AGGAUSS )

      RETURN
END SUBROUTINE AGVITBF
