      SUBROUTINE SOM2VED( NBVAR, B01, B02, B )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT :    B <= B01 + B02
! -----    B , B01 , B02 VECTEURS DE NBVAR VARIABLES DOUBLE PRECISION
!          VARIANTE OpenMP f95
!
! ENTREES:
! --------
! NBVAR  : NOMBRE DE VARIABLES DES VECTEURS B, B01, B02
! B01    : PREMIER VECTEUR REEL DOUBLE PRECISION DE LA SOMME
! B02    : SECOND  VECTEUR REEL DOUBLE PRECISION DE LA SOMME
!
! SORTIE :
! --------
! B      : VECTEUR SOMME DES 2 VECTEURS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : PERRONNET ALAIN LJLL UPMC & St Pierre du Perray Novembre 2012
!2345X7..............................................................012
!$ USE OMP_LIB
   IMPLICIT  NONE
   include"./incl/threads.inc"
   DOUBLE PRECISION  B01(NBVAR), B02(NBVAR), B(NBVAR)
   INTEGER           NBVAR, K

      IF( NBVAR .LE. MININDSTHR ) THEN

         !EXECUTION AVEC UN SEUL THREAD
         DO K = 1, NBVAR
            B(K) = B01(K) + B02(K)
         ENDDO

      ELSE

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL PRIVATE(K) SHARED(NBVAR,B,B01,B02,NBTHREADS)
!$OMP DO SCHEDULE( STATIC, NBVAR/NBTHREADS )
         DO K = 1, NBVAR
            B(K) = B01(K) + B02(K)
         ENDDO
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////
      ENDIF

      RETURN
END SUBROUTINE SOM2VED
