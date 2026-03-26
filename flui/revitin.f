      SUBROUTINE REVITIN ( NYOBJT, NUOBJT, XPI, YPI, ZPI, MNVIIN,
     %                     VITEIN )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LE TABLEAU VITEIN DES VITESSES DU FLUIDE
C -----    au point XPI,YPI,ZPI
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C XPI,YPI,ZPI : LES 3 COORDONNEES DU POINT D'INTEGRATION NUMERIQUE
C MNVIIN : ADRESSE MCN DU TABLEAU 'VITFLUIN'
C
C SORTIE :
C --------
C VITEIN : LE TABLEAU (NBCOVI) DES VITESSES INITIALES
C
C ATTENTION: SEULES NBCOVI COMPOSANTES DES VITESSES PARMI
C ---------- LES NDIM POSSIBLES SONT INITIALISEES!
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris     Mai 2007
C23456---------------------------------------------------------------012
      include"./incl/a___vitfluin.inc"
      include"./incl/ctemps.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      DOUBLE PRECISION XPI, YPI, ZPI, XYZ(7)
      DOUBLE PRECISION VITEIN(3)
C
C     LE NOMBRE DE COMPOSANTES DE VITESSES A DEFINIR
      NBCOVI = MCN( MNVIIN + WBCOVI )
C
C     REMPLISSAGE SELON LE TYPE DES DONNEES DES VITESSES
      IF( MCN( MNVIIN + WTVIF0 ) .EQ. 1 ) THEN
C
C        VALEURS CONSTANTES POUR CET OBJET
C        =================================
         MN = MNVIIN + WUCOVI - 1
         DO 10 J=1,NBCOVI
C           NUMERO DE LA COMPOSANTE DE LA VITESSE
            N = MCN( MN + J )
C           LA VALEUR DE LA VITESSE INITIALE
            VITEIN( N ) = RMCN( MN + NBCOVI + J )
 10      CONTINUE
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
         MN = MNVIIN + WUCOVI
C        NO DE LA FONCTION VITESSE INITIALE DU FLUIDE
         NOFONC = MCN( MN + NBCOVI )
         DO 20 J=1,NBCOVI
C           NUMERO DE LA COMPOSANTE DE LA VITESSE
            N  = MCN( MN - 1 + J )
            XYZ(7) = N
            CALL FONVAL( NOFONC, 7, XYZ,
     %                   NCODEV, VITEIN(N) )
 20      CONTINUE
      ENDIF
C
      RETURN
      END
