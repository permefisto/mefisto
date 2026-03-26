      SUBROUTINE RECOIN( NYOBJT,NUOBJT,NBCOMP,XPI,YPI,ZPI,MNCOIN,
     %                   CONINI )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LE TABLEAU CONINI COEFFICIENTS DU TENSEUR SYMETRIQUE
C -----    DES CONTRAINTES INITIALES EN UN POINT D'INTEGRATION NUMERIQUE
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE DE L'OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C NBCOMP : 3 EN 2D , 4 EN AXISYMETRIQUE , 6 EN 3D
C XPI,YPI,ZPI : LES 3 COORDONNEES DU POINT D'INTEGRATION
C MNCOIN : ADRESSE MCN DU TABLEAU 'CONTRINIT'
C
C SORTIE :
C --------
C CONINI : LE TENSEUR SYMETRIQUE DES CONTRAINTES INITIALES DANS L'ORDRE
C          EN 2D:            1:SXX 2:SYY 3:SXY
C          EN AXISYMETRIQUE: 1:SZZ 2:SRR 3:STETA 4:SRZ
C          EN 3D:            1:SXX 2:SYY 3:SZZ 4:SXY 5:SYZ 6:SZX
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1990
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/donela.inc"
      include"./incl/a___contrinit.inc"
      include"./incl/ctemps.inc"
C
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))
C
      DOUBLE PRECISION XPI,YPI,ZPI,CONINI(1:NBCOMP),PXYZ(7)
C
C     LE TYPE DES DONNEES
C     ===================
      LTCOIN = MCN( MNCOIN + WTCOIN )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTCOIN .EQ. 1 ) THEN
C
C        CONSTANTE
C        =========
         MN = MNCOIN + WONINI
         DO 10 J=1,NBCOMP
            CONINI(J) = RMCN( MN )
            MN = MN + 1
 10      CONTINUE
C
      ELSE
C
C        FONCTION UTILISATEUR
C        ====================
         PXYZ(1) = TEMPS
         PXYZ(2) = XPI
         PXYZ(3) = YPI
         PXYZ(4) = ZPI
         PXYZ(5) = NYOBJT
         PXYZ(6) = NUOBJT
         DO 30 J=1,NBCOMP
            PXYZ(7) = J
            CALL FONVAL( MCN(MNCOIN+WFCOIN), 7, PXYZ,
     %                   NCODEV, CONINI(J) )
 30      CONTINUE
      ENDIF
      END
