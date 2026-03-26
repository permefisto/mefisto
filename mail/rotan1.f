      SUBROUTINE ROTAN1( PT1, PT2, NBSOM, SOMMET,
     %                   NBSOAX, NOSOAX, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCUL ET LOCALISATION DES SOMMETS SUR L'AXE PT1-PT2
C -----
C ENTREES :
C ---------
C PT1 PT2 : LES 3 COORDONNEES DES 2 POINTS DE DEFINITION DE L'AXE
C NBSOM   : LE NOMBRE DE SOMMETS A COMPARER
C SOMMET  : LES 3 COORDONNEES DES NBSOM SOMMETS
C
C SORTIES :
C ---------
C NBSOAX  : LE NOMBRE DE SOMMETS SUR L'AXE
C NOSOAX  : LE NUMERO DU SOMMET S'IL N'EST PAS SUR L'AXE PT1-PT2
C           -LE NUMERO DU SOMMET S'IL EST SUR L'AXE PT1-PT2
C IERR    : 1 SI PT1=PT2 AXE MAL DEFINI
C           0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS  SEPTEMBRE 1990
C23456---------------------------------------------------------------012
      COMMON / EPSSSS / EPZERO, EPSXYZ
      REAL              PT1(3), PT2(3), SOMMET(3,NBSOM)
      INTEGER           NOSOAX(NBSOM)
C
C     LES 2 POINTS PT1-PT2 SONT ILS IDENTIQUES
      CALL  XYZIDE( PT1 , PT2 , IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     LA LONGUEUR PT1-PT2
      D12 = DIST2P( PT1, PT2 ) * EPSXYZ
C
      NBSOAX = 0
      DO 10 N=1,NBSOM
C
C        LE SOMMET N EST-IL SUR L'AXE PT1-PT2 ?
         D3 = DISTPD( SOMMET(1,N), PT1, PT2 )
         IF( D3 .GT. D12 ) THEN
C
C           SOMMET NON SUR L'AXE
            NOSOAX( N ) = N
C
         ELSE
C
C           SOMMET SUR L'AXE
            NBSOAX = NBSOAX + 1
            NOSOAX( N ) = -N
         ENDIF
 10   CONTINUE
      END
