      SUBROUTINE AZEROI( NB, NTAB )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INITIALISATION A ZERO D UN TABLEAU NTAB DE NB VARIABLES ENTIERES
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C MODIFS : PERRONNET ALAIN LJLL UPMC & St Pierre du Perray    AVRIL 2013
C23456---------------------------------------------------------------012
C$    USE OMP_LIB
      include"./incl/threads.inc"
      INTEGER  NB, NTAB(1:NB), I

      IF( NB .LE. MININDSTHR .OR. NB .GT. 2**30 ) THEN

C        EXECUTION AVEC UN SEUL THREAD
         DO I = 1 , NB
            NTAB( I ) = 0
         ENDDO

      ELSE

C///////////////////////////////////////////////////////////////////////
C$OMP PARALLEL PRIVATE( I )
C$OMP DO SCHEDULE( STATIC, NB/NBTHREADS )
         DO I = 1, NB
            NTAB( I ) = 0
         ENDDO
C$OMP END DO
C$OMP END PARALLEL
C///////////////////////////////////////////////////////////////////////

      ENDIF

      RETURN
      END
