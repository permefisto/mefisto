      SUBROUTINE REPREIN( NYOBJT, NUOBJT, XPI, YPI, ZPI, MNPRIN,
     %                    PRESSION )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LA VARIABLE PRESSION AVEC LA VALEUR DE LA PRESSION
C -----    LUE DANS PREFLUIN
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C XPI,YPI,ZPI : LES 3 COORDONNEES DU POINT D'INTEGRATION NUMERIQUE
C MNPRIN : ADRESSE MCN DU TABLEAU 'PREFLUIN'
C
C SORTIE :
C --------
C PRESSION: LA VALEUR DE LA PRESSION A BLOQUER
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris     Mai 2007
C23456---------------------------------------------------------------012
      IMPLICIT INTEGER(W)
      include"./incl/donflu.inc"
      include"./incl/a___prefluin.inc"
      include"./incl/ctemps.inc"
      DOUBLE PRECISION XPI,YPI,ZPI,XYZ(6)
      DOUBLE PRECISION PRESSION
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))
C
C     REMPLISSAGE SELON LE TYPE DE CONDITION INITIALE DE LA PRESSION
      IF( MCN( MNPRIN + WTPRE0 ) .EQ. 1 ) THEN
C
C        VALEURS CONSTANTES POUR CE FLUIDE
C        =================================
         PRESSION = RMCN( MNPRIN + WAPRE0 )
C
      ELSE
C
C        FONCTION UTILISATEUR
C        ====================
         XYZ(1) = TEMPS
         XYZ(2) = XPI
         XYZ(3) = YPI
         XYZ(4) = ZPI
         XYZ(5) = NYOBJT
         XYZ(6) = NUOBJT
         CALL FONVAL( MCN(MNPRIN+WFPRE0), 6, XYZ,
     %                NCODEV, PRESSION )
      ENDIF
      END
