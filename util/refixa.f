      SUBROUTINE REFIXA( NYOBJT, NUOBJT, XPI,YPI,ZPI, MNFIXA,
     %                   NBCOFI, FIXAT )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LE TABLEAU FIXAT DES FIXATIONS EXERCEES SUR L'OBJET
C -----
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C XPI,YPI,ZPI : LES 3 COORDONNEES DU POINT D'INTEGRATION NUMERIQUE
C MNFIXA : ADRESSE MCN DU TABLEAU 'FIXATION'
C
C SORTIES:
C --------
C NBCOFI : LE NOMBRE DE COMPOSANTES DU DEPLACEMENT A FIXER
C FIXAT  : LE TABLEAU (NBCOFI) DES FIXATIONS EXERCEES SUR L'OBJET
C
C ATTENTION: SEULES NBCOFI COMPOSANTES DES DEPLACEMENTS A FIXER
C ---------- PARMI LES NDIM POSSIBLES SONT INITIALISEES!
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1990
C....6...............................................................012
      include"./incl/donela.inc"
      include"./incl/a___fixation.inc"
      include"./incl/ctemps.inc"
      DOUBLE PRECISION  XPI,YPI,ZPI,XYZ(7)
      DOUBLE PRECISION  FIXAT(NBCOFI)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
C     LE NOMBRE DE COMPOSANTES DE FIXATIONS A DEFINIR
      NBCOFI = MCN( MNFIXA + WBCOFI )
C
C     REMPLISSAGE SELON LE TYPE DES DONNEES DES FIXATIONS
      MN = MNFIXA + WUCOFI + NBCOFI
      IF( MCN( MNFIXA + WTFIXA ) .EQ. 1 ) THEN
C
C        VALEURS CONSTANTES POUR CET OBJET
C        =================================
         DO 10 J=1,NBCOFI
C           LA VALEUR DE FIXATION EST EMPILEE
            FIXAT( J ) = RMCN( MN - 1 + J )
 10      CONTINUE
C
      ELSE
C
C        FONCTIONS UTILISATEURS
C        ======================
         XYZ(1) = TEMPS
         XYZ(2) = XPI
         XYZ(3) = YPI
         XYZ(4) = ZPI
         XYZ(5) = NYOBJT
         XYZ(6) = NUOBJT
         DO 20 J=1,NBCOFI
C           LE NUMERO DE LA COMPOSANTE DU DEPLACEMENT A FIXER
            XYZ(7) = MCN( MNFIXA + WUCOFI + J - 1 )
            CALL FONVAL( MCN(MN), 7, XYZ,
     %                   NCODEV, FIXAT(J) )
 20      CONTINUE
C
      ENDIF
C
      RETURN
      END
