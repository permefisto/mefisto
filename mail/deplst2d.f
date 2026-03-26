       SUBROUTINE DEPLST2D( MNXYZS,  NUMST0, X0, Y0, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DEPLACER UN SOMMET DANS UN MAILLAGE 2D DEFINI PAR XYZSOMMET
C -----
C ENTREES:
C --------
C MNXYZS : ADRESSE DU TABLEAU XYZSOMMET DE LA SURFACE
C
C SORTIES:
C --------
C NUMST0 : NUMERO DU SOMMET DEPLACE
C X0, Y0 : COORDONNEES INITIALES DU SOMMET NUMST0
C IERR   : 0 SI PAS D'ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Labo J-L. LIONS  UPMC  PARIS   SEPTEMBRE 2007
C2345X7..............................................................012
       include"./incl/a___xyzsommet.inc"
      include"./incl/pp.inc"
       COMMON         MCN(MOTMCN)
       REAL          RMCN(1)
       EQUIVALENCE  (RMCN(1),MCN(1))
       REAL          COORD(1:2)
C
       IERR = 0
C
C      CLIC D'UN SOMMET SUR L'ECRAN
C      ============================
       CALL SAIPTC( NOTYEV,  NX, NY, NOCHAR )
       IF( NOTYEV .LE. 0 ) THEN
C          -1 => CARACTERE TAPE AU CLAVIER AVEC NOCHAR
C           0 => ABANDON
          IERR = 1
          RETURN
       ENDIF
C
C      COORDONNEES OBJET DU SOMMET CLIQUE
C      ==================================
       X = XOB2PX( NX )
       Y = YOB2PX( NY )
C
C      RECHERCHE DU SOMMET LE PLUS PROCHE NUMST0
C      =========================================
       MN = MNXYZS + WYZSOM
       COORD(1) = RMCN(MN)
       COORD(2) = RMCN(MN+1)
       DISTANCE = SQRT( (X-COORD(1))**2 + (Y-COORD(2))**2 )
       NUMST0  = 1
       DISTMP = DISTANCE
C      ON BOUCLE SUR LES SOMMETS DU MAILLAGE
       DO 10 I=2,MCN(MNXYZS+WNBSOM)
          MN = MN + 3
          COORD(1) = RMCN(MN)
          COORD(2) = RMCN(MN+1)
          DISTMP = SQRT( (X-COORD(1))**2+(Y-COORD(2))**2 )
          IF( DISTMP .LT. DISTANCE ) THEN
             DISTANCE = DISTMP
             NUMST0 = I
          ENDIF
 10    CONTINUE
C
C      CLIC D'UNE POSITION NOUVELLE DU SOMMET SUR L'ECRAN
C      ==================================================
       CALL SAIPTC( NOTYEV,  NX, NY, NOCHAR )
       IF( NOTYEV .LE. 0 ) THEN
C          -1 => CARACTERE TAPE AU CLAVIER AVEC NOCHAR
C           0 => ABANDON
          IERR = 1
          NUMST0 = 0
          RETURN
       ENDIF
C
C      ANCIENNES COORDONNEES OBJET DU POINT CLIQUE
C      ===========================================
       MN = MNXYZS + WYZSOM + 3*NUMST0 - 3
       X0 = RMCN(MN  )
       Y0 = RMCN(MN+1)
C
C      NOUVELLES COORDONNEES OBJET DU POINT CLIQUE
C      ===========================================
       RMCN(MN  ) = XOB2PX( NX )
       RMCN(MN+1) = YOB2PX( NY )
       RMCN(MN+2) = 0
C
       RETURN
       END
