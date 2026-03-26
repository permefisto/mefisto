      SUBROUTINE REECHA( NYOBJT, NUOBJT, NBCOOR, XYZP, MNECHA,  ECHANG )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LA VALEUR DU COEFFICIENT D'ECHANGE DE L'OBJET
C -----    DE TYPE NYOBJT ET DE NUMERO NUOBJT AU POINT (XPI,YPI,ZPI)
C          EN UN POINT D'INTEGRATION NUMERIQUE A L'INSTANT TEMPS
C
C ENTREES :
C ---------
C NYOBJT : NUMERO DU TYPE D'OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C XYZP   : LES NBCOOR COORDONNEES DU POINT D'INTEGRATION
C MNECHA : ADRESSE MCN DU TABLEAU 'ECHANGE'
C
C SORTIE :
C --------
C ECHANG : DENSITE DE ECHANGE
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1990
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/donthe.inc"
      include"./incl/a___echange.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      DOUBLE PRECISION  XYZP(NBCOOR), ECHANG
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  XYZ(10)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
C     LE TYPE DES DONNEES DE L'ECHANGE
      LTECHA = MCN( MNECHA + WTECHA )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTECHA .EQ. 1 ) THEN
C
C        MATERIAU HOMOGENE ISOTROPE
C        ==========================
         ECHANG = RMCN( MNECHA + WCHANG )
C
      ELSE
C
C        FONCTION UTILISATEUR
C        ====================
         XYZ(1) = TEMPS
         N = 1
         DO 5 I=1,NBCOOR
            N = N + 1
            XYZ(N) = XYZP(I)
 5       CONTINUE
         XYZ(N+1) = NYOBJT
         XYZ(N+2) = NUOBJT
C        TEMPERATURE CALCULEE ET STOCKEE DANS LE COMMON cthet.inc
         N = N + 3
         XYZ(N) = TEMPEL
         CALL FONVAL( MCN(MNECHA+WFECHA), N, XYZ,
     %                NCODEV, ECHANG )
C
      ENDIF
C
      RETURN
      END
