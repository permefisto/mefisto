        SUBROUTINE LISTLO( MNSOLI, RLONG )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA SOMME DES LONGUEURS DES ARETES D'UNE
C -----    LIGNE STRUCTUREE SUPPOSEE EXISTANTE ET SANS ERREUR
C
C ENTREES:
C --------
C MNSOLI : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA LIGNE
C          CF '~/TD/D/A___XYZSOMMET'
C
C SORTIES:
C --------
C RLONG  : SOMME DES LONGUEURS DES ARETES ENTRE LES POINTS SUCCESSIFS
C          D'UNE LIGNE STRUCTUREE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS       SEPTEMBRE 1993
C234567..............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/a___xyzsommet.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
C     LA BOUCLE SUR LES SOMMETS
      RLONG = 0.0
      MN0   = MNSOLI + WYZSOM
      DO 10 I = 2, MCN( MNSOLI + WNBSOM )
         MN1   = MN0 + 3
         RLONG = RLONG + DIST2P( RMCN(MN0), RMCN(MN1) )
         MN0   = MN1
 10   CONTINUE
      END
