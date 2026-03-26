      SUBROUTINE REELAS( NYOBJT,NUOBJT,NDIM,XPI,YPI,ZPI,MNYOUN,  ELAS)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT : REMPLIR LE TABLEAU ELAS COEFFICIENTS DU TENSEUR SYMETRIQUE
C ----- D'ELASTICITE SOUS FORME D'UN VECTEUR EN UN POINT D'INTEGRATION
C       NUMERIQUE
C
C ENTREES :
C ---------
C NYOBJT : NUMERO DU TYPE DE L'OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C NDIM   : 2  OU 3 DIMENSION DE L'ESPACE
C          23 SI TRAITEMENT AXISYMETRIQUE
C
C XPI,YPI,ZPI : LES 3 COORDONNEES DU POINT D'INTEGRATION
C MNYOUN : ADRESSE MCN DU TABLEAU 'YOUNG'
C
C SORTIE :
C --------
C ELAS : LE TENSEUR SYMETRIQUE DE L'ELASTICITE RANGE SOUS LA FORME
C        1: E11  2: E21  3:E22  4:E31, ...
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS  JUILLET 1989
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/donela.inc"
      include"./incl/a___young.inc"
      include"./incl/ctemps.inc"
      DOUBLE PRECISION XPI,YPI,ZPI,ELAS(1:*),PXYZ(6)
      DOUBLE PRECISION YOUNG,POISON,D1,D2,D3
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))
C
C     LE TYPE DES DONNEES DE L'ELASTICITE
      LTYOUN = MCN( MNYOUN + WTYOUN )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTYOUN .EQ. 10 ) THEN
C
C        MATERIAU ANISOTROPE CONSTANT
C        ============================
         DO 10 I=1,MCN(MNYOUN+WVELAS)
            ELAS(I) = RMCN( MNYOUN + WEELAS + I - 1 )
 10      CONTINUE
         RETURN
C
      ELSE IF( LTYOUN .GE. 1 .AND. LTYOUN .LE. 4 ) THEN
C
C        MATERIAU HOMOGENE ISOTROPE
C        ==========================
         YOUNG  = RMCN( MNYOUN + WYOUNG )
         POISON = RMCN( MNYOUN + WOISON )
C
      ELSE
C
C        FONCTIONS UTILISATEURS
C        ======================
         PXYZ(1) = TEMPS
         PXYZ(2) = XPI
         PXYZ(3) = YPI
         PXYZ(4) = ZPI
         PXYZ(5) = NYOBJT
         PXYZ(6) = NUOBJT
         CALL FONVAL( MCN(MNYOUN+WFYOUN), 6, PXYZ,  NCODEV, YOUNG )
         CALL FONVAL( MCN(MNYOUN+WFPOIS), 6, PXYZ,  NCODEV, POISON )
      ENDIF
C
C     TRANSFORMATION DE YOUNG ET POISSON EN TENSEUR D'ELASTICITE
C     ==========================================================
C     LE CAS AXISYMETRIQUE
      IF( NDIM .EQ. 23 ) LTYOUN = 2
      LTYOUN = ABS( LTYOUN )
C
C     LES COEFFICIENTS NON NULS DE LA MATRICE D'ELASTICITE
      D3       = YOUNG / ( 1.D0 + POISON )
      D1       = D3    / ( 1.D0 - 2.D0 * POISON )
      D2       = D1 * POISON
      D1       = D1 * ( 1.D0 - POISON )
      D3       = D3 * 0.5D0
C
      IF( LTYOUN .EQ. 1 ) THEN
C
C        MATERIAU 3D HOMOGENE ISOTROPE
C        LES 3 PREMIERS COEFFICIENTS DIAGONAUX
         ELAS( 1) = D1
         ELAS( 3) = D1
         ELAS( 6) = D1
C        LES 3 COEFFICIENTS NON DIAGONAUX
         ELAS( 2) = D2
         ELAS( 4) = D2
         ELAS( 5) = D2
C        LES 3 DERNIERS COEFFICIENTS DIAGONAUX
         ELAS(10) = D3
         ELAS(15) = D3
         ELAS(21) = D3
C        LES COEFFICIENTS NULS
         ELAS( 7) = 0D0
         ELAS( 8) = 0D0
         ELAS( 9) = 0D0
         ELAS(11) = 0D0
         ELAS(12) = 0D0
         ELAS(13) = 0D0
         ELAS(14) = 0D0
         ELAS(16) = 0D0
         ELAS(17) = 0D0
         ELAS(18) = 0D0
         ELAS(19) = 0D0
         ELAS(20) = 0D0
C
      ELSE IF( LTYOUN .EQ. 2 ) THEN
C
C        MATERIAU AXISYMETRIQUE HOMOGENE ISOTROPE
         ELAS( 1) = D1
         ELAS( 3) = D1
         ELAS( 6) = D1
         ELAS( 2) = D2
         ELAS( 4) = D2
         ELAS( 5) = D2
         ELAS( 7) = 0.D0
         ELAS( 8) = 0.D0
         ELAS( 9) = 0.D0
         ELAS(10) = D3
C
      ELSE IF( LTYOUN .EQ. 3 ) THEN
C
C        MATERIAU 2D CONTRAINTES PLANES
         D1  = YOUNG / ( 1D0 - POISON * POISON )
         ELAS( 1) = D1
         ELAS( 3) = D1
         ELAS( 2) = D1 * POISON
         ELAS( 4) = 0.D0
         ELAS( 5) = 0.D0
         ELAS( 6) = D3
C
      ELSE IF( LTYOUN .EQ. 4 ) THEN
C
C        MATERIAU 2D DEFORMATIONS PLANES
         ELAS( 1) = D1
         ELAS( 2) = D2
         ELAS( 3) = D1
         ELAS( 4) = 0.D0
         ELAS( 5) = 0.D0
         ELAS( 6) = D3
      ENDIF
      END
