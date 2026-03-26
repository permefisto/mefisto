      SUBROUTINE RECPRE( NYOBJT, NUOBJT, XPI, YPI, ZPI, MNCPRE,
     %                   COPRESS )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    RETROUVER LA VALEUR DU COEFFICIENT DU GRADIENT DE LA PRESSION
C -----    DU FLUIDE EN UN POINT (XPI, YPI, ZPI) A L'INSTANT TEMPS
C          DANS NAVIER-STOKES STANDARD CE COEFFICIENT VAUT 1
C          CE COEFFICIENT PERMET DE TESTER NAVIER-STOKES EN 2D
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE D'OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C
C XPI,YPI,ZPI: LES 3 COORDONNEES DU POINT D'INTEGRATION
C MNCPRE : ADRESSE MCN DU TABLEAU 'COEFPRESSION'
C
C SORTIES:
C --------
C COPRESS: VALEUR du COEFFICIENT DE LA PRESSION AU POINT XPI, YPI, ZPI
C          SI ERREUR DE DONNEES 1D0 EST RETOURNE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du Perray Decembre 2009
C....6...............................................................012
      include"./incl/donflu.inc"
      include"./incl/a___coefpression.inc"
      include"./incl/ctemps.inc"
C
      COMMON /UNITES/ LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION XPI,YPI,ZPI,PXYZ(6)
      DOUBLE PRECISION COPRESS
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))
C
C     LE TYPE DES DONNEES DE COEFFICIENT DE LA PRESSION
      LTCPRE = MCN( MNCPRE + WTCPRE )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTCPRE .EQ. 1 ) THEN
C
C        MATERIAU HOMOGENE ISOTROPE
C        ==========================
         COPRESS = RMCN( MNCPRE + WPRESS )
C
      ELSE IF( LTCPRE .EQ. -1 ) THEN
C
C        FONCTION UTILISATEUR
C        ====================
         PXYZ(1) = TEMPS
         PXYZ(2) = XPI
         PXYZ(3) = YPI
         PXYZ(4) = ZPI
         PXYZ(5) = NYOBJT
         PXYZ(6) = NUOBJT
         CALL FONVAL( MCN(MNCPRE+WFCPRE), 6, PXYZ,
     %                NCODEV, COPRESS )
C        NCODEV: 0 COPRESS N'EST PAS INITIALISE EN SORTIE
C                1 COPRESS   EST     INITIALISE EN SORTIE
         IF( NCODEV .NE. 1 ) COPRESS = 1.D0
C
      ELSE
C
C        ERREUR DE DONNEES
         COPRESS = 1.D0
C
      ENDIF
C
      RETURN
      END
