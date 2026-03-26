      SUBROUTINE CL2VED( NBVAR, ALFA1, B01, ALFA2, B02, B )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT :    B <= ALFA1 * B01 + ALFA2 * B02
! -----    B , B01 , B02 VECTEURS DE NBVAR VARIABLES DOUBLE PRECISION
!          ALFA1 , ALFA2 REELS DOUBLE PRECISION
!          VARIANTE OpenMP f95
!
! ENTREES:
! --------
! NBVAR  : NOMBRE DE VARIABLES DES VECTEURS B, B01, B02
! ALFA1  : PREMIER REEL DOUBLE PRECISION DE LA COMBINAISON LINEAIRE
! B01    : PREMIER VECTEUR REEL DOUBLE PRECISION DE LA COMBINAISON LINEAIRE
! ALFA2  : SECOND  REEL DOUBLE PRECISION DE LA COMBINAISON LINEAIRE
! B02    : SECOND  VECTEUR REEL DOUBLE PRECISION DE LA COMBINAISON LINEAIRE
!
! SORTIE :
! --------
! B      : VECTEUR RESULTAT DE LA COMBINAISON LINEAIRE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : PERRONNET ALAIN LJLL UPMC & St Pierre du Perray  Octobre 2012
!2345X7..............................................................012
!$ USE OMP_LIB
IMPLICIT  NONE
DOUBLE PRECISION  B01(NBVAR), B02(NBVAR), B(NBVAR), ALFA1, ALFA2
INTEGER           NBVAR, K

!///////////////////////////////////////////////////////////////////////
!!$OMP PARALLEL 
!!$OMP MASTER
!      print *,'CL2VEDomp: Nombre THREADS=',OMP_GET_NUM_THREADS()
!!$OMP END MASTER
!!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!     TESTS SUIVANT LES VALEURS DE ALFA1 ET ALFA2
!     POUR REDUIRE LE NOMBRE DES OPERATIONS
!     =============================================
      IF( ALFA1 .EQ. 0.D0 ) THEN

         IF( ALFA2 .EQ. 1.D0 ) THEN
!           ALFA1=0.D0 ET ALFA2=1.D0
            DO K = 1, NBVAR
               B(K) = B02(K)
            ENDDO

         ELSE

!           ALFA1=0.D0 ET ALFA2 QUELCONQUE
!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DO PRIVATE(K) SHARED(NBVAR,B,ALFA2,B02)
            DO K = 1, NBVAR
               B(K) = ALFA2 * B02(K)
            ENDDO
!$OMP END PARALLEL DO
!///////////////////////////////////////////////////////////////////////

         ENDIF
!
      ELSE IF( ALFA2 .EQ. 0.D0 ) THEN

         IF( ALFA1 .EQ. 1.D0 ) THEN
!           ALFA1=1.D0 ET ALFA2=0.D0
            DO K = 1, NBVAR
               B(K) = B01(K)
            ENDDO

         ELSE

!           ALFA1 QUELCONQUE ET ALFA2=0.D0
!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DO PRIVATE(K) SHARED(NBVAR,B,ALFA1,B01)
            DO K = 1, NBVAR
               B(K) = ALFA1 * B01(K)
            ENDDO
!$OMP END PARALLEL DO
!///////////////////////////////////////////////////////////////////////

         ENDIF

      ELSE IF( ALFA1 .EQ. 1.D0 ) THEN

!        ALFA1=1.D0
         IF( ALFA2 .EQ. 1.D0 ) THEN

!           ALFA1=1.D0,ALFA2=1.D0
!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DO PRIVATE(K) SHARED(NBVAR,B,B01,B02)
            DO K = 1, NBVAR
               B(K) = B01(K) + B02(K)
            ENDDO
!$OMP END PARALLEL DO
!///////////////////////////////////////////////////////////////////////

         ELSE

!           ALFA1=1.D0 ET ALFA2<>1.D0
!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DO PRIVATE(K) SHARED(NBVAR,B,B01,ALFA2,B02)
            DO K = 1, NBVAR
               B(K) = B01(K) + ALFA2 * B02(K)
            ENDDO
!$OMP END PARALLEL DO
!///////////////////////////////////////////////////////////////////////

         ENDIF

      ELSE

!        ALFA1<>1D0
         IF( ALFA2 .EQ. 1.D0 ) THEN

!           ALFA1<>1.D0 ET ALFA2=1.D0
!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DO PRIVATE(K) SHARED(NBVAR,B,ALFA1,B01,B02)
            DO K = 1, NBVAR
               B(K) = ALFA1 * B01(K) + B02(K)
            ENDDO
!$OMP END PARALLEL DO
!///////////////////////////////////////////////////////////////////////
!
         ELSE

!           ALFA1<>1.D0 ET ALFA2<>1.D0
!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DO PRIVATE(K) SHARED(NBVAR,B,ALFA1,B01,ALFA2,B02)
            DO K = 1, NBVAR
               B(K) = ALFA1 * B01(K) + ALFA2 * B02(K)
            ENDDO
!$OMP END PARALLEL DO
!///////////////////////////////////////////////////////////////////////
!
         ENDIF

      ENDIF

      RETURN
END SUBROUTINE CL2VED
