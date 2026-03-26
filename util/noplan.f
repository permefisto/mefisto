      SUBROUTINE NOPLAN( NPPLAN , COEFNO , IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCULER LES COEFFICIENTS DE LA DIRECTION NORMALE AU PLAN
C -----  DEFINI PAR LES 3 POINTS DE NUMERO NPPLAN
C
C ENTREE :
C --------
C NPPLAN : NUMERO DES 3 POINTS DEFINISSANT LE PLAN
C
C SORTIES :
C ---------
C COEFNO : LES 3 COEFFICIENTS DE LA DIRECTION NORMALE AU PLAN
C          LA NORME EUCLIDIENNE DE COEFNO VAUT 1.
C IERR   : 0 SI PAS D'ERREUR ; 1 SINON
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1988
C....................................................................
      include"./incl/pp.inc"
      COMMON  RMCN(MOTMCN)
      INTEGER NPPLAN(3)
      REAL    XYZ(3,3),COEFNO(3),COEPLA(4)
C
C     LES COORDONNEES DES 3 POINTS DU PLAN
      DO 10 I=1,3
         CALL XYZPOI( NPPLAN(I) , MN , IERR )
         IF( IERR .NE. 0 ) RETURN
         XYZ(1,I) = RMCN( MN )
         XYZ(2,I) = RMCN( MN + 1 )
         XYZ(3,I) = RMCN( MN + 2 )
 10   CONTINUE
C
C     L'EQUATION DU PLAN
      CALL EQPLAN( XYZ , COEPLA , IERR )
      IF( IERR .GT. 0 ) RETURN
C
C     NORMALISATION A 1 DE LA DIRECTION NORMALE
      S = 0.
      DO 20 I=1,3
         S = S + COEPLA(I) ** 2
 20   CONTINUE
      S = 1. / SQRT(S)
      DO 30 I=1,3
         COEFNO(I) = COEPLA(I) * S
 30   CONTINUE
      END
