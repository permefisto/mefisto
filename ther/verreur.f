      SUBROUTINE VERREUR(NBCOOI, NOFOTI, MODECO, MNPOGE,
     %                   NDSM,   NTDL,   TEMPER, TIMES,
     %                   ERRMIN, NOEMIN, NCAMIN, ERRMAX, NOEMAX, NCAMAX,
     %                   MOERRE, MNERRE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER L'ERREUR |ERREUR EXACTE(Noeud)-ERREUR CALCULEE(Noeud)|
C -----    A PARTIR DE LA DONNEE DE LA FONCTION UTILISATEUR
C          ou TEMPERATURE_EXACTE(t,x,y,z)
C          ou TEMPERATURE_EXACTE(t,x,y,z,u,v,w)
C          ou DEPLACEMENT_EXACT(t,x,y,z,nc)
C          ou PARTIE_PARTIE_REELLE_EXACTE(t,x,y,z)
C          ou PARTIE_IMAGINAIRE_EXACTE(t,x,y,z)
C
C ENTREES:
C --------
C NBCOOI : NOMBRE DE COORDONNEES D'UN NOEUD DU MAILLAGE INITIAL (3 ou 6)
C NOFOTI : NUMERO DE LA FONCTION TEMPERATURE_EXACTE(t,x,y,z)
C          DANS LE LEXIQUE DES FONCTIONS DE L'UTILISATEUR
C NDIM   : DIMENSION DE L'ESPACE DE L'OBJET (2 OU 3 OU 6)
C KNOMOB : NOM DE L'OBJET
C MODECO : MODE DE TRACE DES VECTEURS
C          =1 CE SONT DES TEMPERATURES
C          =2 CE SONT DES VECTEURS PROPRES
C          =3 CE SONT DES PRESSIONS P1 A PARTIR D'UNE INTERPOLATION
C             SOIT P2 SOIT P1+BULLE P3
C          =4 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C          =8 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C             DU MODULE D'UNE ONDE COMPLEXE
C            (CALCUL DU MODULE DE L'ERREUR COMPLEXE A FAIRE)
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS FINIS DANS CETTE TOPOLOGIE
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE
C MNELEM : ADRESSE MCN DES TABLEAUX ELEMENTS FINIS DE CETTE TOPOLOGIE
C MNPOGE : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C MNNOEU : ADRESSE MCN DU TABLEAU NOEUDS D'INTERPOLATION DU MAILLAGE
C NDPGST : CODE TRAITEMENT DES SOMMETS POINTS NOEUDS DU MAILLAGE
C          0 : NOEUDS=POINTS=SOMMETS
C          1 : NOEUDS=POINTS#SOMMETS
C          2 : NOEUDS#POINTS=SOMMETS
C          3 : NOEUDS#POINTS#SOMMETS
C NDSM   : NOMBRE DE CAS OU SECONDS MEMBRES DU SYSTEME LINEAIRE
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE
C TEMPER : TEMPERATURE DES NDSM CAS AVEC NTDL NOEUDS
C TIMES  : TEMPS DU CALCUL DES NDSM VECTEURS
C
C SORTIES:
C --------
C ERRMIN : ERREUR MINIMALE EN UN NOEUD
C NOEMIN : NUMERO DU NOEUD OU L'ERREUR EST MINIMALE
C NCAMIN : NUMERO DU CAS DU MINIMUM
C ERRMAX : ERREUR MAXIMALE EN UN NOEUD
C NOEMAX : NUMERO DU NOEUD OU L'ERREUR EST MAXIMALE
C NCAMAX : NUMERO DU CAS DU MAXIMUM
C MOERRE : NOMBRE DE MOTS DU TABLEAU DES ERREURS
C MNERRE : ADRESSE MCN DU TMC ERREUR AVEC LA STRUCTURE DE TMS
C          DE 'VECTEUR"...' SANS LES TEMPS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Juillet 2011
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
ccc      include"./incl/inteel.inc"
      include"./incl/ctemps.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (DMCN(1),RMCN(1),MCN(1))
      DOUBLE PRECISION  TEMPER(NTDL,NDSM)
      REAL              TIMES(NDSM)
      DOUBLE PRECISION  D, TEXACT, IEXACT, DPARAF(8)
      INTEGER           NOFOWE(2)
C
      IERR    = 0
      MNERRE  = 0
      NOPROJ  = 0
C
C     NOMBRE DE MOTS D'UN REEL DOUBLE PRECISION
      MOREE2 = MOTVAR(6)
C
C     NOMBRE DE COORDONNEES DES NOEUDS
      NBCOOR = MCN(MNPOGE+WBCOOP)
C
C     CALCUL DU TABLEAU DES ERREURS A L'ADRESSE MNERRE
C     CALCUL DU NOMBRE DE PARAMETRES DE LA FONCTION EXACTE
C     TEMPERATURE_EXACTE(t,x,y,z) ou DEPLACEMENT_EXACT(t,x,y,z,nc)
C     ou PARTIE_PARTIE_REELLE_EXACTE(t,x,y,z)
C     ou PARTIE_IMAGINAIRE_EXACTE(t,x,y,z)
      NOFOWE(1) = NOFOPREX()
      NOFOWE(2) = NOFOPIEX()
C
      IF( NOFOTEEX() .EQ. NOFOTI ) THEN
         IF( NBCOOI .LE. 3 ) THEN
C           TEMPERATURE_EXACTE(t,x,y,z)
            NBARGS = 4
         ELSE
C           TEMPERATURE_EXACTE(t,u,v,w,x,y,z)
            NBARGS = 7
         ENDIF
      ELSE IF( NOFODEEX() .EQ. NOFOTI ) THEN
C        DEPLACEMENT_EXACT(t,x,y,z,nc)
         NBARGS = 5
      ELSE IF( NOFOWE(1) .EQ. NOFOTI ) THEN
C        PARTIE_PARTIE_REELLE_EXACTE(t,x,y,z)
         NBARGS = 4
      ELSE IF( NOFOWE(2) .EQ. NOFOTI ) THEN
C        PARTIE_IMAGINAIRE_EXACTE(t,x,y,z)
         NBARGS = 4
      ELSE
C        ERREUR PAS DE TRACE
         IERR = 1
         RETURN
      ENDIF
C
C     CREATION DU TMC ERREUR AVEC LA STRUCTURE DE TMS DE 'VECTEUR"...'
      MOERRE = WECTEU + (NTDL*MOREE2)*NDSM
      CALL TNMCDC( 'ENTIER', MOERRE, MNERRE )
C deftms ~>>>VECTEUR 'TABLEAU de VECTEURS'
C variable muet DATE  'date de creation du tableau'         reel2  ;
C variable muet NUTD  'numero du tableau descripteur'       entier ;
C variable NBCOVE 'Nombre de composantes de chaque VECTEUR' entier ;
C variable NBVECT 'Nombre de VECTEURS'                      entier ;
C variable NBCPIN 'Nombre d''informations complementaires'  entier ;
C tableau  VECTEU(1..NBCOVE,1..NBVECT) 'Composantes des VECTEURS' reel2 ;
C tableau  CPINFO(1..NBCPIN) 'Complements,temps,...' reel ;
C fintms ;
      MCN( MNERRE + WBCOVE ) = NTDL
      MCN( MNERRE + WBVECT ) = NDSM
      MCN( MNERRE + WBCPIN ) = NDSM
C
C     OPTIONS DU TRACE DE L'ERREUR D'UN NCAS
      ERRMIN  = 1E28
      ERRMAX  =-1E28
C
C     BOUCLE SUR LES VECTEURS SOLUTIONS
C     =================================
      MNE     =(MNERRE + WECTEU - 1)/MOREE2
      DO 20 NCAS = 1, NDSM
C
C        LE TEMPS ACTUEL
         TEMPS = TIMES(NCAS)
C
         MNS = MNPOGE + WYZPOI - 1
C
         DO 15 NS=1,MCN(MNPOGE+WNBPOI)
C
C           CALCUL DES PARAMETRES DE LA SOLUTION EXACTE EN CE NOEUD
            DPARAF(1) = TEMPS
            DO K=1,NBCOOR
               DPARAF(1+K) = RMCN( MNS + K )
            ENDDO
C           NUMERO DE LA SEULE COMPOSANTE 1 ICI
            DPARAF(2+NBCOOR) = 1D0
C
            IF( MODECO .NE. 8 ) THEN
C              1 SEULE FONCTION UTILISATEUR A CALCULER
C                 TEMPERATURE_EXACTE(TEMPS,X,Y,Z)
C              ou TEMPERATURE_EXACTE(TEMPS,X,Y,Z,U,V,W)
C              ou DEPLACEMENT_EXACT (TEMPS,X,Y,Z,NOCOMP)
C              ou PARTIE_REELLE_EXACTE(TEMPS,X,Y,Z)
C              ou PARTIE_IMAGINAIRE_EXACTE(TEMPS,X,Y,Z)
               CALL FONVAL( NOFOTI, NBARGS, DPARAF, NCODEV, TEXACT )
            ELSE
C              2 FONCTIONS UTILISATEUR A CALCULER ET MODULE A CALCULER
C              FONCTION PARTIE_REELLE_EXACTE(t,x,y,z)
               CALL FONVAL( NOFOWE(1), NBARGS, DPARAF, NCODEV, TEXACT )
C              FONCTION PARTIE_IMAGINAIRE_EXACTE(t,x,y,z)
               CALL FONVAL( NOFOWE(2), NBARGS, DPARAF, NCODEV, IEXACT )
C              LE MODULE
               TEXACT = SQRT( TEXACT**2 + IEXACT**2 )
            ENDIF
            IF( NCODEV .GT. 0 ) THEN
C
C              D L'ERREUR ABSOLUE AU NOEUD NS
               D = ABS( TEXACT - TEMPER(NS,NCAS) )
               DMCN( MNE + NS ) = D
C
               IF( D .LT. ERRMIN ) THEN
                  NOEMIN = NS
                  ERRMIN = REAL( D )
                  NCAMIN = NCAS
               ENDIF
               IF( D .GT. ERRMAX ) THEN
                  NOEMAX = NS
                  ERRMAX = REAL( D )
                  NCAMAX = NCAS
               ENDIF
            ENDIF
C
C           COORDONNEES DU POINT SUIVANT
            MNS = MNS + NBCOOR
 15      CONTINUE
C
C        LE PROCHAIN VECTEUR
         MNE = MNE + NTDL
C
 20   CONTINUE
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNERRE) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNERRE + MOTVAR(6) ) = NONMTD( '~>>>VECTEUR' )
C
C     AFFICHAGE DE L'ERREUR MIN ET MAX
      MNS = MNPOGE + WYZPOI - NBCOOR - 1
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10020) NCAMIN, ERRMIN, NOEMIN,
     %                      (RMCN(MNS+NBCOOR*NOEMIN+K),K=1,NBCOOR)
         WRITE(IMPRIM,10021) NCAMAX, ERRMAX, NOEMAX,
     %                      (RMCN(MNS+NBCOOR*NOEMAX+K),K=1,NBCOOR)
      ELSE
         WRITE(IMPRIM,20020) NCAMIN, ERRMIN, NOEMIN,
     %                      (RMCN(MNS+NBCOOR*NOEMIN+K),K=1,NBCOOR)
         WRITE(IMPRIM,20021) NCAMAX, ERRMAX, NOEMAX,
     %                      (RMCN(MNS+NBCOOR*NOEMAX+K),K=1,NBCOOR)
      ENDIF
10020 FORMAT(/' VECTEUR',I8,
     %        ' ERREUR MINIMALE=',G15.7,' au NOEUD',I12,' COOR ',6G15.7)
10021 FORMAT( ' VECTEUR',I8,
     %        ' ERREUR MAXIMALE=',G15.7,' au NOEUD',I12,' COOR ',6G15.7)
20020 FORMAT(/' VECTOR',I8,
     %        ' MINIMUM ERROR=',G15.7,' at NODE',I12,' COOR ',6G15.7)
20021 FORMAT( ' VECTOR',I8,
     %        ' MAXIMUM ERROR=',G15.7,' at NODE',I12,' COOR ',6G15.7)
C
      RETURN
      END
