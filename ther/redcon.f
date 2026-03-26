      SUBROUTINE REDCON( NYOBJT, NUOBJT, NBCOOR, XYZP, MNDCON,
     %                   DCOND)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LE TABLEAU DCOND COEFFICIENTS DU GRADIENT DU TENSEUR
C -----    DE CONDUCTIVITE SOUS FORME D'UN VECTEUR
C          EN UN POINT D'INTEGRATION NUMERIQUE A L'INSTANT TEMPS
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE D'OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C XYZP   : LES NBCOOR COORDONNEES DU POINT D'INTEGRATION
C MNDCON : ADRESSE MCN DU TABLEAU 'DCONDUCTIVITE'
C
C SORTIES:
C --------
C DCOND : GRADIENT DU TENSEUR SYMETRIQUE DE CONDUCTIVITE RANGE SOUS LA FORME
C         1: C11  2: C21  3:C22  4:C31, ...
C         LE SECOND INDICE POUR DX, DY, DZ
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1995
C....6...............................................................012
      include"./incl/donthe.inc"
      include"./incl/a___dconductivite.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      DOUBLE PRECISION XYZP(NBCOOR),DCOND(1:6,1:3),DCONDU(3),PXYZ(11)
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))
C
C     LE TYPE DES DONNEES DE LA CONDUCTIVITE
      LTDCON = MCN( MNDCON + WTDCON )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTDCON .EQ. 10 ) THEN
C
C        MATERIAU ANISOTROPE CONSTANT
C        ============================
         DO 10 I=1,MCN(MNDCON+WVDCON)
            DO 5 J=1,3
               DCOND(I,J) = RMCN( MNDCON + WEDCON + J + 3 * I - 4 )
 5          CONTINUE
 10      CONTINUE
         RETURN
C
      ELSE IF( LTDCON .GE. 1 .AND. LTDCON .LE. 3 ) THEN
C
C        MATERIAU HOMOGENE ISOTROPE
C        ==========================
         DO 15 J=1,3
            DCONDU(J) = RMCN( MNDCON + WCONDU + J - 1 )
 15      CONTINUE
C
      ELSE IF( LTDCON .EQ. -1 ) THEN
C
C        FONCTION UTILISATEUR
C        ====================
         PXYZ(1) = TEMPS
         N = 1
         DO 22 I=1,NBCOOR
            N = N + 1
            PXYZ(N) = XYZP(I)
 22      CONTINUE
         PXYZ(N+1) = NYOBJT
         PXYZ(N+2) = NUOBJT
C        TEMPERATURE CALCULEE ET STOCKEE DANS LE COMMON cthet.inc
         N = N + 3
         PXYZ(N) = TEMPEL
         DO 25 J=1,3
C           LA DERIVEE PAR RAPPORT A XJ
            PXYZ(7) = J
            CALL FONVAL( MCN(MNDCON+WFDCON), N, PXYZ,
     %                   NCODEV, DCONDU(J) )
25       CONTINUE
C
      ELSE
C
C        ERREUR DE DONNEES
         DCONDU(1) = 0.D0
         DCONDU(2) = 0.D0
         DCONDU(3) = 0.D0
      ENDIF
C
C     CREATION DU TENSEUR DE DCONDUCTIVITE ISOTROPE HOMOGENE
C     =====================================================
C     MATERIAU 3D 2D HOMOGENE ISOTROPE
C     MATERIAU AXISYMETRIQUE HOMOGENE ISOTROPE
      DO 35 J=1,3
C        LES 3 PREMIERS COEFFICIENTS DIAGONAUX
         DCOND( 1, J ) = DCONDU(J)
         DCOND( 3, J ) = DCONDU(J)
         DCOND( 6, J ) = DCONDU(J)
C        LES 3 COEFFICIENTS NON DIAGONAUX
         DCOND( 2, J ) = 0.D0
         DCOND( 4, J ) = 0.D0
         DCOND( 5, J ) = 0.D0
 35   CONTINUE
C
      RETURN
      END
