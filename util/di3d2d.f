      SUBROUTINE DI3D2D( P1 , P2 , P3 , D2D3 , IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : GENERER LA MATRICE DE PASSAGE DU SYSTEME D'AXE ORIGINAL
C ----- AU SYSTEME D'AXES P1P2 , ORTHOGONAL A P1P2 DANS LE PLAN P1 P2 P3
C       ORTHOGONAL AU PLAN P1 P2 P3
C
C  (X       )    (          ) ( XX )    ( XX ) ( T         )( X      )
C  (Y  - P1 ) =  (  D2D3    ) ( YY ) ;  ( YY )=(   D2D3    )( Y - P1 )
C  (Z       )    (          ) ( ZZ )    ( ZZ ) (           )( Z      )
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
C          1 SI LES 3 POINTS SONT ALIGNES   (D2D3 NON CALCULE)
C         -1 SI LES 3 POINTS SONT PROCHES D'ETRE ALIGNES
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  JANVIER 1986
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/gsmenu.inc"
      DOUBLE PRECISION  D2D3(3,3),DPOINT(3,3)
      DOUBLE PRECISION  P1(3),P2(3),P3(3)
C
C     TRANSFERT DES 3 POINTS DANS LE TABLEAU DPOINT
      DO 10 I=1,3
         DPOINT(I,1) = P1(I)
         DPOINT(I,2) = P2(I)
         DPOINT(I,3) = P3(I)
 10   CONTINUE
C
C     EQUATION DU PLAN FORME PAR LES 3 POINTS
      CALL EQPLAD( DPOINT, D2D3, IERR )
      IF( IERR .GT. 0 ) THEN
C        LES 3 POINTS SONT ALIGNES OU PROCHES DE L'ETRE
CCC         NBLGRC(NRERR) = 1
CCC         KERR(1) = 'SP DI3D2D: 3 POINTS ALIGNES'
CCC         CALL LEREUR
CCC         WRITE(IMPRIM,10000) ((DPOINT(I,J),I=1,3),J=1,3)
CCC10000 FORMAT(3(' X=',G17.9,' Y=',G17.9,' Z=',G17.9/) )
         RETURN
      ENDIF
C
C     LES AXES ORTHONORMAUX P1P2 ,...
      CALL P3AXE3( DPOINT, D2D3, IERR )
      END
