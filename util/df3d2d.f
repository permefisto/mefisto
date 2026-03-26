      SUBROUTINE DF3D2D( P1 , P2 , P3 , D2D3 , IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : GENERER LA MATRICE DE PASSAGE DU SYSTEME D'AXE ORIGINAL
C ----- AU SYSTEME D'AXES P1P2 , ORTHOGONAL A P1P2 DANS LE PLAN P1 P2 P3
C       ORTHOGONAL AU PLAN P1 P2 P3
C
C  (X      )    (         ) ( XX )    ( XX ) ( T       )( X      )
C  (Y - P1 ) =  (   D2D3  ) ( YY ) ;  ( YY )=(   D2D3  )( Y - P1 )
C  (Z      )    (         ) ( ZZ )    ( ZZ ) (         )( Z      )
C
C          AVEC DANS LES APPLICATIONS ZZ = 0.
C     X  , Y  , Z  COORDONNEES DANS LE REPERE 3D
C     XX , YY , ZZ COORDONNEES DANS LE REPERE 2D DU PLAN P1 P2 P3
C
C ENTREES :
C ---------
C P1 P2 P3 : 3 POINTS DEFINISSANT LE PLAN 2D
C
C SORTIE :
C --------
C D2D3   : MATRICE DE PASSAGE D'UN SYSTEME A L'AUTRE
C IERR   : 0 SI PAS D'ERREUR ET LES 3 POINTS NON COLINEAIRES
C          1 SI LES 3 POINTS SONT ALIGNES  (D2D3 NON CALCULE)
C         -1 SI LES 3 POINTS SONT PROCHES D'ETRE ALIGNES
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS      JANVIER 1986
C........................................................................
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/gsmenu.inc"
      DOUBLE PRECISION D2D3(3,3),POINT(3,3)
      REAL             P1(3),P2(3),P3(3),RPOINT(3,3),RD2D3(3,3)
      EQUIVALENCE      (POINT,RPOINT)
C
C     TRANSFERT DES 3 POINTS DANS LE TABLEAU RPOINT
      DO 10 I=1,3
         RPOINT(I,1) = P1(I)
         RPOINT(I,2) = P2(I)
         RPOINT(I,3) = P3(I)
 10   CONTINUE
C
C     EQUATION DU PLAN FORME PAR LES 3 POINTS
      CALL EQPLAN( RPOINT , RD2D3 , IERR )
      IF( IERR .GT. 0 ) THEN
C        LES 3 POINTS SONT ALIGNES
CCC         NBLGRC(NRERR) = 1
CCC         KERR(1) = 'DF3D2D: 3 POINTS ALIGNES'
CCC         CALL LEREUR
CCC         WRITE(IMPRIM,10000) ((RPOINT(I,J),I=1,3),J=1,3)
CCC10000 FORMAT(' LES 3 POINTS SONT ALIGNES' /
CCC     %3(' X=',G15.6,' Y=',G15.6,' Z=',G15.6/) )
         RETURN
      ENDIF
C
C     TRANSFORMATION EN DOUBLE PRECISION DES POINTS
      DO 20 I=1,3
         POINT(I,1) = P1(I)
         POINT(I,2) = P2(I)
         POINT(I,3) = P3(I)
 20   CONTINUE
C
C     LES AXES ORTHONORMAUX P1P2 ,...
      CALL P3AXE3( POINT, D2D3, IERR )
      END
