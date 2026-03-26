      SUBROUTINE MULTPROC( NBDOF, NBCOEF, ROW, COLUMN, A, B,
     &                     X,     IERR )
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    THE MULTI-PROCESSORS METHOD to SOLVE A x = b
C -----    Version muette pour permettre l'edition de liens
C
C INPUT  :
C --------
C NBDOF  : NUMBER OF ROWS (and COLUMNS) of the MATRIX A
C NBCOEF : TOTAL NUMBER OF THE STORED COEFFICIENTS OF THE MATRIX A
C
C ROW    : NUMBER OF LINE   I OF EVERY COEFFICIENT OF A
C COLUMN : NUMBER OF COLUMN J OF EVERY COEFFICIENT OF A
C A      : NON ZERO COEFFICIENTS AIJ OF THE MATRIX A
C B      : SECOND MEMBER VECTOR of A x = b
C
C OUTPUT :
C --------
C X      : SOLUTION VECTOR of A x = B
C IERR   : ERROR CODE : 0 IF NO ERROR, >0 ELSE
C23456---------------------------------------------------------------012
C AUTHOR : Alain PERRONNET LJLL UPMC & St Pierre du PERRAY Decembre 2010
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/langue.inc"
      INTEGER           LECTEU, IMPRIM, NUNITE
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           IERR,   NBDOF,  NBCOEF
      INTEGER           ROW(1:NBCOEF)
      INTEGER           COLUMN(1:NBCOEF)
      REAL(8)           A(1:NBCOEF), B(1:NBDOF), X(1:NBDOF)
C
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'SP MULTPROC: VERSION NON PROGRAMMEE'
         WRITE(IMPRIM,*) 'CHOISIR UNE AUTRE METHODE DE RESOLUTION'
         WRITE(IMPRIM,*) 'SP UTILE POUR PERMETTRE L''EDITION DE LIENS'
      ELSE
         WRITE(IMPRIM,*) 'SP MULTPROC: NOT PROGRAMMED METHOD'
         WRITE(IMPRIM,*) 'CHOOSE AN OTHER SOLVER OF LINEAR SYSTEM'
         WRITE(IMPRIM,*) 'SUBROUTINE USEFUL TO PERMIT THE BINDING'
      ENDIF
C
      WRITE(IMPRIM,*) ROW(1), COLUMN(1), A(1), B(1), X(1)
C
      IERR = 1
      RETURN
      END
