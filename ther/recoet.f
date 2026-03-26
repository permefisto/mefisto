      SUBROUTINE RECOET( NYOBJT, NUOBJT, NBCOOR, XYZPI,
     &                   MNCOET, COEFTE )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    CALCULER LE COEFFICIENT DE LA TEMPERATURE
C ---      EN UN POINT A L'INSTANT TEMPS
C
C ENTREES :
C ---------
C NYOBJT : NUMERO DU TYPE D'OBJET 1:POINT, 2:LIGNE, 3:SURFACE, 4:VOLUME
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C
C NBCOOR : NOMBRE DE COORDONNEES DES POINTS D'INTEGRATION 3 ou 6
C          EN 2D FOURNIR NBCOOR=3 ET LA VALEUR XYZPI(3)=0D0
C XYZPI  : LES NBCOOR COORDONNEES DU POINT DE CALCUL DU COEFFICIENT
C MNCOET : ADRESSE MCN DU TABLEAU a___coeftemperature
C
C SORTIE :
C --------
C COEFTE : LA VALEUR DU COEFFICIENT DE LA TEMPERATURE
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        AOUT 1998
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/donthe.inc"
      include"./incl/a___coeftemperature.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      DOUBLE PRECISION  COEFTE,PXYZ(10),XYZPI(NBCOOR)
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))
C
C     LE TYPE DES DONNEES DU COEFFICIENT DE LA TEMPERATURE
      LTCOET = MCN( MNCOET + WTCOET )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTCOET .EQ. 1 ) THEN
C
C        MATERIAU HOMOGENE ISOTROPE
C        ==========================
         COEFTE = RMCN( MNCOET + WOEFTE )
         RETURN
C
      ELSE IF( LTCOET .EQ. -1 ) THEN
C
C        FONCTION UTILISATEUR(t,x,y,z,ntyobj,nuobj,Temperature)
C        ====================
         PXYZ(1) = TEMPS
         N = 1
         DO 5 I=1,NBCOOR
            N = N+ 1
            PXYZ( N ) = XYZPI(I)
 5       CONTINUE
         PXYZ(N+1) = NYOBJT
         PXYZ(N+2) = NUOBJT
C        CAS LINEAIRE OU NON LINEAIRE
         N = N + 3
         PXYZ(N) = TEMPEL
C
C        CALCUL DE LA VALEUR DU COEFFICIENT DE LA TEMPERATURE
         CALL FONVAL( MCN(MNCOET+WFCOET), N, PXYZ,
     %                NCODEV, COEFTE )
C
      ELSE
C
C        ERREUR DE DONNEES
         COEFTE = 0.D0

      ENDIF
C
      RETURN
      END
