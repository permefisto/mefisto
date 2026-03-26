      SUBROUTINE RECATHMA( NYOBJT, NUOBJT, NBCOOR, XYZPI, MNCHMA,  Cp )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    CALCULER la CHALEUR MASSIQUE ou la CHALEUR SPECIFIQUE
C -----    ou renommee CAPACITE THERMIQUE MASSIQUE en J/Kelvin/kg
C          (Notee Cp dans les FORMULES)
C          EN UN POINT XYZPI A L'INSTANT TEMPS

C ENTREES :
C ---------
C NYOBJT : NUMERO DU TYPE D'OBJET 1:POINT, 2:LIGNE, 3:SURFACE, 4:VOLUME
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C
C NBCOOR : NOMBRE DE COORDONNEES DES POINTS D'INTEGRATION 3 ou 6
C          EN 1D FOURNIR NBCOOR=3 ET LES VALEURS XYZPI(2:3)=0
C          EN 2D FOURNIR NBCOOR=3 ET LA VALEUR XYZPI(3)=0
C XYZPI  : LES NBCOOR COORDONNEES DU POINT DE CALCUL DE LA CHALEUR MASSIQUE
C MNCHMA : ADRESSE MCN DU TABLEAU 'CHALEUR MASSIQUE'

C SORTIE :
C --------
C Cp     : LA CHALEUR MASSIQUE (ou CAPACITE THERMIQUE MASSIQUE)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1999
C....6...............................................................012
      include"./incl/a___chaleurmassique.inc"
      include"./incl/donthe.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      DOUBLE PRECISION XYZPI(NBCOOR), Cp, PXYZ(10)
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))

C     LE TYPE DES DONNEES DE LA CHALEUR MASSIQUE
C     ==========================================
      IF( MNCHMA .LE. 0 ) THEN
         Cp = 0D0
         GOTO 9999
      ENDIF
C
C     REMPLISSAGE SELON LE TYPE DE DONNEE DE LA CHALEUR MASSIQUE
      LTCHMA = MCN( MNCHMA + WTCHMA )
      IF( LTCHMA .EQ. 1 ) THEN

C        MATERIAU HOMOGENE ISOTROPE
C        ==========================
         Cp = RMCN( MNCHMA + WHAMAS )

      ELSE IF( LTCHMA .EQ. -1 ) THEN

C        FONCTION UTILISATEUR
C        F(temps, x,y,z, no_type_object, materiau, temperature)
C        ======================================================
         PXYZ(1) = TEMPS
         N = 1
         DO 5 I=1,NBCOOR
            N = N + 1
            PXYZ(N) = XYZPI(I)
 5       CONTINUE
         PXYZ(N+1) = NYOBJT
         PXYZ(N+2) = NUOBJT
C        CAS LINEAIRE OU NON LINEAIRE
         N = N + 3
         PXYZ(N) = TEMPEL

C        CALCUL DE LA VALEUR DE LA CHALEUR MASSIQUE
         CALL FONVAL( MCN(MNCHMA+WFCHMA), N, PXYZ,
     %                NCODEV, Cp )

      ELSE

C        ERREUR SUR LA DONNEE DE LA CHALEUR MASSIQUE
         Cp = 0.D0

      ENDIF

 9999 RETURN
      END
