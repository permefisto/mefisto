      SUBROUTINE RECBOU( NYOBJT, NUOBJT, XPI, YPI, ZPI, MNCBOU,
     %                   CoBOUSSI )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT : RETROUVER LA VALEUR DU COEFFICIENT DU COEFFICIENT DE PROPORTIONNALITE
C ----- DE LA DIFFERENCE DE TEMPERATURE A LA VARIATION DE LA
C       MASSE VOLUMIQUE dans L'APPROXIMATION de BOUSSINESQ
C       Rho0-Rho = Rho CoBOUSSI (T-T0)
C       DU FLUIDE EN UN POINT (XPI, YPI, ZPI) A L'INSTANT TEMPS

C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE D'OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C
C XPI,YPI,ZPI: LES 3 COORDONNEES DU POINT DE CALCUL
C MNCBOU : ADRESSE MCN DU TABLEAU 'COBOUSSI'
C
C SORTIES:
C --------
C CoBOUSSI: VALEUR du COEFFICIENT DU COEFFICIENT DE PROPORTIONNALITE
C           DE LA DIFFERENCE DE TEMPERATURE A LA VARIATION DE LA
C           MASSE VOLUMIQUE dans L'APPROXIMATION de BOUSSINESQ
C           Rho0-Rho=Rho CoBOUSSI (T-T0)
C           AU POINT XPI, YPI, ZPI
C           SI ERREUR DE DONNEES CoBOUSSI=0D0 EST RETOURNE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du Perray            Avril 2022
C....6...............................................................012
      include"./incl/donflu.inc"
      include"./incl/a___coefboussinesq.inc"
      include"./incl/ctemps.inc"
C
      COMMON /UNITES/ LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION XPI,YPI,ZPI,PXYZ(6)
      DOUBLE PRECISION CoBOUSSI
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))
C
C     LE TYPE DES DONNEES DE COEFFICIENT DE BOUSSINESQ
      LTCBOU = MCN( MNCBOU + WTCBOU )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTCBOU .EQ. 1 ) THEN
C
C        MATERIAU HOMOGENE ISOTROPE
C        ==========================
         COBOUSSI = RMCN( MNCBOU + WBOUSS )
C
      ELSE IF( LTCBOU .EQ. -1 ) THEN
C
C        FONCTION UTILISATEUR
C        ====================
         PXYZ(1) = TEMPS
         PXYZ(2) = XPI
         PXYZ(3) = YPI
         PXYZ(4) = ZPI
         PXYZ(5) = NYOBJT
         PXYZ(6) = NUOBJT
         CALL FONVAL( MCN(MNCBOU+WFCBOU), 6, PXYZ,
     %                NCODEV, CoBOUSSI )
C        NCODEV: 0 COBOUSSI N'EST PAS INITIALISE EN SORTIE
C                1 COBOUSSI   EST     INITIALISE EN SORTIE
C
      ELSE
C
C        ERREUR DE DONNEES
         CoBOUSSI = 0.D0
C
      ENDIF
C
      RETURN
      END
