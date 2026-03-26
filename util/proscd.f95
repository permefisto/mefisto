DOUBLE PRECISION FUNCTION PROSCD( V1, V2, N )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT : PRODUIT SCALAIRE DE 2 VECTEURS V1 V2 REELS DOUBLE PRECISION
! ----- DE N COMPOSANTES AVEC NBTHREAD THREADS
!       V1 et V2 SONT DANS LA MEMOIRE PARTAGEE A TOUS LES THREADS
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEURS: Philippe PARNAUDEAU  LJLL UPMC Paris             Octobre 2012
!          Alain    PERRONNET   LJLL UPMC Paris & St Pierre du Perray
!.......................................................................
!$ USE OMP_LIB
   IMPLICIT  NONE
   include"./incl/threads.inc"
   INTEGER           N, K
   DOUBLE PRECISION  V1(N), V2(N), PS, S

      PS = 0.d0

      IF( N .LE. 0 ) GOTO 9999

      IF( N .LE. MININDSTHR .OR. N .LE. NBTHREADS*MININD1THR ) THEN

         !EXECUTION AVEC UN SEUL THREAD
         DO K = 1, N
            PS = PS + V1(K) * V2(K)
         ENDDO

      ELSE

!///////////////////////////////////////////////////////////////////////
!EXECUTION AVEC NBTHREADS
!$OMP PARALLEL PRIVATE(K) SHARED(N,V1,V2)
!$OMP DO REDUCTION(+:PS)

         DO K = 1, N
            PS = PS + V1(K) * V2(K)
         ENDDO

!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

      ENDIF

 9999 PROSCD = PS
!     print *,'PROSCD: (V1,V2)=',PROSCD,' pour',N,' composantes'

      RETURN
END FUNCTION PROSCD
