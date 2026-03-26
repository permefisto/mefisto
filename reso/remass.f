      SUBROUTINE REMASS( NYOBJT,NUOBJT, NBCOOR, XYZPI, MNMASS,  DMASSE )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LA VALEUR DE LA DENSITE DE MASSE AU POINT XPI,YPI,ZPI
C -----    A L'INSTANT TEMPS
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE DE L'OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C NBCOOR : NOMBRE DE COORDONNEES DES POINTS D'INTEGRATION 3 ou 6
C          EN 2D FOURNIR NBCOOR=3 ET LA VALEUR XYZPI(3)=0
C XYZPI  : LES NBCOOR COORDONNEES DU POINT DE CALCUL DE LA DENSITE DE MASSE
C MNMASS : ADRESSE MCN DU TABLEAU 'MASSE'
C          PAR DEFAUT DE DONNEE (MNMASS=0) => DMASSE=1D0
C
C SORTIE :
C --------
C DMASSE : DENSITE DE MASSE
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JUILLET 1989
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/donela.inc"
      include"./incl/a___masse.inc"
      include"./incl/ctemps.inc"
C
      DOUBLE PRECISION  XYZPI(NBCOOR), DMASSE
      DOUBLE PRECISION  XYZ(9)
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
C     PAR DEFAUT DE DONNEE MNMASS = 0 => DMASSE=1D0
      IF( MNMASS .LE. 0 ) THEN
         DMASSE = 1D0
         RETURN
      ENDIF
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
C        FONCTION UTILISATEUR(temps, x,y,z, no_type_object, materiau)
C        ============================================================
         XYZ(1) = TEMPS
         N = 1
         DO 5 I=1,NBCOOR
            N = N + 1
            XYZ( N ) = XYZPI(I)
 5       CONTINUE
         XYZ(N+1) = NYOBJT
         XYZ(N+2) = NUOBJT
         N = N + 2
         CALL FONVAL( MCN(MNMASS+WFDMAS), N, XYZ,
     %                NCODEV, DMASSE )
C
      ELSE
C
C        ERREUR DE DONNEES
         DMASSE = 0.D0
C
      ENDIF
      END
