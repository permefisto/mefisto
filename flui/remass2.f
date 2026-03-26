      SUBROUTINE REMASS2( NYOBJT, NUOBJT, XPI, YPI, ZPI, MNMASS,
     &                    DMASSE )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    RETROUVER LA VALEUR DE LA DENSITE DE MASSE DU FLUIDE
C -----    EN UN POINT (XPI, YPI, ZPI) A L'INSTANT TEMPS
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE D'OBJET (1:POINT, 2:LIGNE, 3:SURFACE, 4:VOLUME )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C
C XPI,YPI,ZPI : LES 3 COORDONNEES DU POINT D'INTEGRATION
C MNMASS : ADRESSE MCN DU TABLEAU 'MASSE'
C
C SORTIES:
C --------
C DMASSE : VALEUR DE LA DENSITE DE MASSE AU POINT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : OLIVIER CIONI ANALYSE NUMERIQUE UPMC PARIS    JANVIER 2001
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/donflu.inc"
      include"./incl/a___masse.inc"
      include"./incl/ctemps.inc"
C
      DOUBLE PRECISION XPI,YPI,ZPI,PXYZ(6)
      DOUBLE PRECISION DMASSE
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))
C
C     LE TYPE DES DONNEES DE LA MASSE
      LTDMAS = MCN( MNMASS + WTDMAS )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTDMAS .EQ. 1 ) THEN
C
C        MATERIAU HOMOGENE ISOTROPE
C        ==========================
         DMASSE = RMCN( MNMASS + WMASSE )
C
      ELSE IF( LTDMAS .EQ. -1 ) THEN
C
C        FONCTION UTILISATEUR
C        ====================
         PXYZ(1) = TEMPS
         PXYZ(2) = XPI
         PXYZ(3) = YPI
         PXYZ(4) = ZPI
         PXYZ(5) = NYOBJT
         PXYZ(6) = NUOBJT
         CALL FONVAL( MCN(MNMASS+WFDMAS), 6, PXYZ,
     %                NCODEV, DMASSE )
C
      ELSE
C
C        ERREUR DE DONNEES
         DMASSE = 0.D0
C
      ENDIF
      END
