      SUBROUTINE REVIIN( NYOBJT, NUOBJT, NDIM, XNO,YNO,ZNO, MNVIT0,
     %                   VITES0 )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LA VALEUR DU GRADIENT DU DEPLACEMENT
C -----    C-A-D LA VITESSE A L'INSTANT TEMPS
C         (EN FAIT C'EST LE PLUS SOUVENT A L'INSTANT INITIAL TEMPS=0!)
C          DANS LE CAS OU IL EST CONSTANT EN TOUS LES NOEUDS
C          OU DEFINI PAR UNE FONCTION UTILISATEUR EN UN NOEUD
C          DE COORDONNEES XNO,YNO,ZNO A L'INSTANT TEMPS ET POUR LE CAS 1
C
C ENTREES :
C ---------
C NYOBJT : NUMERO DU TYPE D'OBJET (1:POINT, 2:LIGNE, ..., 5:OBJET )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C NDIM   : 1 OU 2 OU 3 DIMENSION DE L'ESPACE ( 2 SI PB AXISYMETRIQUE )
C XNO,YNO,ZNO : LES 3 COORDONNEES DU NOEUD
C MNVIT0 : ADRESSE MCN DU TABLEAU 'VITESSEINIT'
C
C SORTIE :
C --------
C VITES0 : TABLEAU (NDIM) DES NDIM COMPOSANTES DE LA VITESSE
C          AU NOEUD XNO, YNO, ZNO A L'INSTANT TEMPS
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1998
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/donela.inc"
      include"./incl/a___vitesseinit.inc"
      include"./incl/ctemps.inc"
      include"./incl/gsmenu.inc"
      include"./incl/cnonlin.inc"
C
      DOUBLE PRECISION  XNO, YNO, ZNO, VITES0(1:NDIM)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  XYZ(7)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
C     TYPE DES DONNEES DE LA VITESSE OU GRADIENT DU DEPLACEMENT AU NOEUD
      LTVIT0 = MCN( MNVIT0 + WTVIT0 )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTVIT0 .EQ. 1 ) THEN
C
C        VITESSE CONSTANTE
C        =================
         IF( MCN( MNVIT0 + WBVIT0 ) .NE. NDIM ) THEN
            WRITE(KERR(4)(1:1),'(I1)') NDIM
            NBLGRC(NRERR) = 3
            KERR(1) ='ERREUR: NOMBRE DE COMPOSANTES VITESSE INITIALE'
            KERR(2) ='NON EGAL A ' // KERR(4)(1:1)
            KERR(3) ='VITESSE NULLE IMPOSEE'
            CALL LEREUR
            DO 10 I=1,NDIM
               VITES0(I) = 0D0
 10         CONTINUE
         ELSE
            DO 20 I=1,NDIM
               VITES0(I) = RMCN( MNVIT0 + WAVIT0 - 1 + I )
 20         CONTINUE
         ENDIF
C
      ELSE IF( LTVIT0 .EQ. -1 ) THEN
C
C        FONCTION UTILISATEUR
C        ====================
         XYZ(1) = TEMPS
         XYZ(2) = XNO
         XYZ(3) = YNO
         XYZ(4) = ZNO
         XYZ(5) = NYOBJT
         XYZ(6) = NUOBJT
         DO 30 I=1,NDIM
C           LE NUMERO DE LA COMPOSANTE DE LA VITESSE
            XYZ(7) = I
            CALL FONVAL( MCN(MNVIT0+WFVIT0), 7, XYZ,
     %                   NCODEV, VITES0(I) )
 30      CONTINUE
C
      ENDIF
      END
