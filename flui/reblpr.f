      SUBROUTINE REBLPR( NYOBJT,NUOBJT,XPI,YPI,ZPI,MNBLPR,
     %                   BLPRES )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LA VARIABLE BLPRES AVEC LA VALEUR DE LA PRESSION
C -----    A BLOQUER
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C NDSM   : NOMBRE DE CAS DE CHARCES
C XPI,YPI,ZPI : LES 3 COORDONNEES DU POINT D'INTEGRATION NUMERIQUE
C MNBLPR : ADRESSE MCN DU TABLEAU 'BLPRESSION'
C
C SORTIE :
C --------
C BLPRES : LA VALEUR DE LA PRESSION A BLOQUER
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : BENHAMADOUCHE SOFIANE                            JANVIER 1998
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/donflu.inc"
      include"./incl/a___blpression.inc"
      include"./incl/ctemps.inc"
      DOUBLE PRECISION XPI,YPI,ZPI,XYZ(6)
      DOUBLE PRECISION BLPRES
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))
C
C     REMPLISSAGE SELON LE TYPE DU BLOCAGE DE LA PRESSION
      IF( MCN( MNBLPR + WTBLPR ) .EQ. 1 ) THEN
C
C        VALEURS CONSTANTES POUR CE FLUIDE
C        =================================
         BLPRES = RMCN( MNBLPR + WLPRES )
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
         CALL FONVAL( MCN(MNBLPR+WFBLPR), 6, XYZ,
     %                NCODEV, BLPRES )
      ENDIF
      END
