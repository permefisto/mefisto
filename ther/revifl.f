      SUBROUTINE REVIFL( NYOBJT, NUOBJT, NDIM, NBCOOR, XYZP, MNVIFL,
     %                   VITEFL )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LA VALEUR DE LA VITESSE DU FLUIDE
C -----    DANS LE CAS OU ELLE EST CONSTANTE
C                      OU DEFINIE PAR UNE FONCTION UTILISATEUR
C          AU NOEUD DE COORDONNEES XNO,YNO,ZNO A L'INSTANT TEMPS POUR LE CAS 1
C
C ENTREES :
C ---------
C NYOBJT : NUMERO DU TYPE D'OBJET (1:POINT, 2:LIGNE, ..., 5:OBJET )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C NDIM   : 1 ou 2 ou 3 DIMENSION DE L'ESPACE  ( 2 SI PB AXISYMETRIQUE )
C NBCOOR : NOMBRE DE COORDONNEES DU NOEUD ( 3 ou 6 )
C XYZP   : LES NBCOOR COORDONNEES DU NOEUD
C MNVIFL : ADRESSE MCN DU TABLEAU 'VITESSEFLUIDE'
C
C SORTIE :
C --------
C VITEFL : TABLEAU (NDIM) DES NDIM COMPOSANTES DE LA VITESSE
C          AU NOEUD XNO, YNO, ZNO A L'INSTANT TEMPS
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1998
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/donela.inc"
      include"./incl/a___vitessefluide.inc"
      include"./incl/ctemps.inc"
      include"./incl/gsmenu.inc"
      include"./incl/cnonlin.inc"
      include"./incl/cthet.inc"
C
      DOUBLE PRECISION  XYZP(NBCOOR), VITEFL(1:NDIM)
      DOUBLE PRECISION  XYZ(11)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))

      IF( MNVIFL .LE. 0 ) THEN
C        PAS DE DONNEE DU TMS a___vitessefluide
         MNVIFL = 0
         RETURN
      ENDIF

C     TYPE DES DONNEES DE LA VITESSE DU FLUIDE
      LTVIFL = MCN( MNVIFL + WTVIFL )

C ======================================================================
C ATTENTION: LE CAS LTVIFL N'EST PAS TRAITE ICI
C            VOIR ts2lag.f ou ts3lag.f POUR LE TRAITER
C ======================================================================

C     REMPLISSAGE SELON LE TYPE
      IF( LTVIFL .EQ. 1 ) THEN
C
C        VITESSE CONSTANTE
C        =================
         IF( MCN( MNVIFL + WBVIFL ) .NE. NDIM ) THEN
            WRITE(KERR(4)(1:1),'(I1)') NDIM
            NBLGRC(NRERR) = 3
            KERR(1) ='ERREUR: NOMBRE COMPOSANTES VITESSE DU FLUIDE'
            KERR(2) ='NON EGAL A ' // KERR(4)(1:1)
            KERR(3) ='VITESSE NULLE IMPOSEE'
            CALL LEREUR
            DO I=1,NDIM
               VITEFL(I) = 0D0
            ENDDO
         ELSE
            DO I=1,NDIM
               VITEFL(I) = RMCN( MNVIFL + WAVIFL - 1 + I )
            ENDDO
         ENDIF
C
      ELSE IF( LTVIFL .EQ. -1 ) THEN
C
C        FONCTION UTILISATEUR
C        ====================
         XYZ(1) = TEMPS
         N = 1
         DO I=1,NBCOOR
            N = N + 1
            XYZ(N) = XYZP(I)
         ENDDO
         XYZ(N+1) = NYOBJT
         XYZ(N+2) = NUOBJT
C        TEMPERATURE CALCULEE ET STOCKEE DANS LE COMMON cthet.inc
         N = N + 3
         XYZ(N) = TEMPEL
C
C        LE NO DE LA DIRECTION
         N = N + 1
C
         DO I=1,NDIM
C           LE NUMERO DE LA COMPOSANTE DE LA VITESSE
            XYZ(N) = I
            CALL FONVAL( MCN(MNVIFL+WFVIFL), N, XYZ,
     %                   NCODEV, VITEFL(I) )
         ENDDO
C
      ENDIF
      END
