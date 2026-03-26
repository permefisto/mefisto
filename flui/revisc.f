      SUBROUTINE REVISC( NYOBJT, NUOBJT, XPI, YPI, ZPI, MNVISC,
     &                   VISCOS )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    RETROUVER LA VALEUR DE LA VISCOSITE DU FLUIDE
C -----    EN UN POINT (XPI, YPI, ZPI) A L'INSTANT TEMPS
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE D'OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C
C XPI,YPI,ZPI : LES 3 COORDONNEES DU POINT D'INTEGRATION
C MNVISC : ADRESSE MCN DU TABLEAU 'VISCOSITE'
C
C SORTIES:
C --------
C VISCOS : VALEUR DE LA VISCOSITE AU POINT XPI, YPI, ZPI
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1998
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/donflu.inc"
      include"./incl/a___viscosite.inc"
      include"./incl/ctemps.inc"
C
      COMMON /UNITES/ LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION XPI,YPI,ZPI,PXYZ(6)
      DOUBLE PRECISION VISCOS
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))
C
C     LE TYPE DES DONNEES DE LA VISCOSITE
      LTVISC = MCN( MNVISC + WTVISC )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTVISC .EQ. 1 ) THEN
C
C        MATERIAU HOMOGENE ISOTROPE
C        ==========================
         VISCOS = RMCN( MNVISC + WISCOS )
C
      ELSE IF( LTVISC .EQ. -1 ) THEN
C
C        FONCTION UTILISATEUR
C        ====================
         PXYZ(1) = TEMPS
         PXYZ(2) = XPI
         PXYZ(3) = YPI
         PXYZ(4) = ZPI
         PXYZ(5) = NYOBJT
         PXYZ(6) = NUOBJT
         CALL FONVAL( MCN(MNVISC+WFVISC), 6, PXYZ,
     %                NCODEV, VISCOS )
C        NCODEV: 0 VISCOS N'EST PAS INITIALISEE EN SORTIE
C                1 VISCOS   EST     INITIALISEE EN SORTIE
C
      ELSE
C
C        ERREUR DE DONNEES
         VISCOS = 0.D0
C
      ENDIF
      END
