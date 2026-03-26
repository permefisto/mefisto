      SUBROUTINE AXEBDD( NTDL, NDSM, A, B, X, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE X TEL QUE A * X = B AVEC
C -----    A MATRICE DIAGONALE
C          B(NDSM,NTDL) LES SECONDS MEMBRES
C          X(NDSM,NTDL) LES NDSM SOLUTIONS
C
C ENTREES:
C --------
C NTDL   : ORDRE DE LA MATRICE A
C NDSM   : NOMBRE DE SECONDS MEMBRES
C A      : MATRICE DIAGONALE A(NTDL)
C B      : NDSM SECONDS MEMBRES
C
C SORTIES:
C -------
C X     : NDSM SOLUTIONS DE A * X = B
C IERR  : 0 SI PAS D'ERREUR DETECTEE, 1 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        AOUT 1978
C ......................................................................
      PARAMETER        (EPS = 1.E-06)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  A(NTDL),B(NDSM,NTDL),X(NDSM,NTDL)
C
    3 FORMAT(' ERREUR AXEBDD : LE',I12,'-EME COEFF DIAGONAL DE ',
     +       'A EST NUL (=',G15.7,').A*X=B NON INVERSIBLE'/
     +        1X,130('%'))
C
      IERR = 0
      DO 10 I=1,NTDL
C
         IF( ABS( A(I) ) .LE. EPS ) THEN
C
            WRITE (IMPRIM,3) I,A(I)
            NBLGRC(NRERR) = 2
            KERR(1) = 'ERREUR AXEBDD: UN COEFFICIENT DIAGONAL NUL'
            KERR(2) = 'DIVISION INTERDITE'
            CALL LEREUR
            IERR = 1
            RETURN
C
         ELSE
            DO 2 N=1,NDSM
               X(N,I)=B(N,I) / A(I)
    2       CONTINUE
C
         ENDIF
C
   10 CONTINUE
      RETURN
      END
