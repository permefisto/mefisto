      SUBROUTINE TRTATD( D1, D2, NBV )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRANSFERT DU TABLEAU D1 DE NBV VARIABLES REELLES DOUBLE
C -----  PRECISION DANS LE TABLEAU D2
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS    NOVEMBRE 1983
C2345X7..............................................................012
C$    USE OMP_LIB
      include"./incl/threads.inc"
      INTEGER           NBV, I
      DOUBLE PRECISION  D1(1:NBV), D2(1:NBV)

      IF( NBV .LE. MININDSTHR ) THEN

         DO I = 1, NBV
            D2(I) = D1(I)
         ENDDO

      ELSE

C///////////////////////////////////////////////////////////////////////
C$OMP PARALLEL PRIVATE( I )
C$OMP DO SCHEDULE( STATIC, NBV/NBTHREADS )
         DO I = 1, NBV
            D2(I) = D1(I)
         ENDDO
C$OMP END DO
C$OMP END PARALLEL
C///////////////////////////////////////////////////////////////////////

      ENDIF

      RETURN
      END
