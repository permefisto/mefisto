      SUBROUTINE NEKTON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INTERFACER MEFISTO AVEC NEKTON POUR UNE QUADRANGULATION 2D
C ----- D'UN OBJET DONT L'INTERPOLATION LAGRANGE DEGRE 1 A ETE FAITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1993
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a_transfo__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___texte.inc"
      include"./incl/msvaau.inc"
      include"./incl/ponoel.inc"
      include"./incl/nbcamo.inc"
      include"./incl/lu.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      COMMON / UNITES /  LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      INTEGER            NOOBVC,NOOBSF(6),NOOBLA(12),NOOBPS(8)
      REAL               RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      CHARACTER*4        NOMELE(2)
      CHARACTER*10       NMTYOB,KNM
      CHARACTER*24       KNOMOB,KNOM
      CHARACTER*5        KFRMLG
      CHARACTER*(NBCALI) KLIG
      CHARACTER*(NBCAMO) CHARX
      LOGICAL            LEXIST,LOPEN
      CHARACTER*1        CHAR
      REAL               CENTRE(2)
C
C     NOM DE L'OBJET A TRAITER
 10   CALL INVITE( 45 )
      IERR   = 0
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNOMOB )
      IF( NCVALS .EQ. -1 ) RETURN
C
C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
      CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNOMOB
         KERR(2) ='ERREUR: OBJET INCONNU'
         CALL LEREUR
         CALL LXIM( NTOBJE )
         GOTO 10
      ENDIF
C
C     TRACE DE L'OBJET
      NBLGRC(NRHIST) = 1
      KHIST(1) = 'OBJET: ' // KNOMOB
      CALL LHISTO
      IF( INTERA .GE. 1 ) THEN
C        TRACE DE L'OBJET SI CONSOLE GRAPHIQUE INTERACTIVE
         CALL T1OBJE( KNOMOB )
      ENDIF
      WRITE(IMPRIM,10001) KNOMOB
10001 FORMAT(/' NOM OBJET POUR NEKTON: ',A)
C
C     RECHERCHE DU TABLEAU DEFINITION
      CALL LXTSOU( NTLXOB , 'DEFINITION' , NTDFOB , MNDFOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNOMOB
         KERR(2) ='ERREUR: OBJET INCONNU'
         CALL LEREUR
         GOTO 10
      ENDIF
C
C     ADRESSAGE DES ADRESSES DES TABLEAUX NPEF"xxxx DE CET OBJET
      MNPIFT = 0
      MNXX   = 0
      MNYY   = 0
      MNZZ   = 0
      MNCCUR = 0
      MNCURV = 0
      MNXYZ  = 0
      MNCALG = 0
      MNCOLG = 0
      MNELEM = 0
      MXTYEL = 3
      CALL TNMCDC( 'ENTIER' , 2*MXTYEL , MNELEM )
      MNTELE = MNELEM + MXTYEL
C
C     RECHERCHE DES TABLEAUX XYZSOMMET XYZNOEUD XYZPOINT ASSOCIES A L'OBJET
      CALL NDPGEL( NTLXOB , NTTOPO , MNTOPO ,
     %             NTPOGE , MNPOGE , NTNOEU , MNNOEU ,
     %             NBTYEL , MCN(MNTELE) , MCN(MNELEM) , IERR )
      IF( IERR .NE. 0 ) GOTO 9990
      IF( NBTYEL .GT. 1 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR: TYPE D''EF DIFFERENT DE QUADRANGLE'
         CALL LEREUR
         GOTO 9990
      ENDIF
C
C     LE NOMBRE DE POINTS DU MAILLAGE DE L'OBJET
      NBPOI = MCN( MNPOGE + WNBPOI )
C     NDIM LA DIMENSION 2 OU 3 DE L'ESPACE DES COORDONNEES
      CALL DIMCOO( NBPOI , MCN(MNPOGE+WYZPOI) , NDIM )
C     LE NOMBRE DE NOEUDS DU MAILLAGE DE L'OBJET
      NBNOEU = MCN(MNNOEU+WNBNOE)
C
      WRITE(IMPRIM,10210) NDIM,NBNOEU
10210 FORMAT(' DIMENSION 2 OU 3 DE L''ESPACE',T32,'=',I6/
     %' NOMBRE DE NOEUDS'                    ,T32,'=',I6)
      IF( NDIM .NE. 2 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'DIMENSION DE L''ESPACE NON REDUITE A 2'
         CALL LEREUR
         GOTO 9990
      ENDIF
      IF( NBNOEU .LE. 4 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'NOMBRE INCORRECT DE NOEUDS'
         CALL LEREUR
         GOTO 9990
      ENDIF
C
C     LE TYPE DE MAILLAGE ET LES OBJETS INTERNES ET AUX LIMITES
      NDPGST = MCN( MNTOPO + WDPGST )
      NBOBIN = MCN( MNTOPO + WBOBIN )
      NBOBCL = MCN( MNTOPO + WBOBCL )
C
C     L'ADRESSE MCN DU DEBUT DU TABLEAU NUOBIN
      MNOBIN = MNTOPO + WMTYEL + NBTYEL
C
C     L'ADRESSE MCN DU DEBUT DU TABLEAU NUOBCL
      MNOBCL = MNOBIN + MOTVAR(13) * NBOBIN
C
C     BOUCLE SUR LES OBJETS "INTERNES"
C     ================================
      MN     = MNOBIN - 2
      DO 200 I=1,NBOBIN
C        LE TYPE DE L'OBJET
         MN   = MN + 2
         NYOB = MCN( MN )
         IF( NYOB .NE. 3 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = KNOMOB
            KERR(2) ='ERREUR: OBJET NON UNE SURFACE'
            CALL LEREUR
            IERR = IERR + 1
            GOTO 200
         ENDIF
         NUOB = MCN( MN + 1 )
C        OUVERTURE DE L'OBJET
         CALL LXNLOU( NTMN(NYOB) , NUOB , NTOB , MNOB )
         IF( NTOB .LE. 0 ) THEN
            KNM = NMTYOB( NYOB )
            CALL NMOBNU( KNM , NUOB , KNOM )
            NBLGRC(NRERR) = 2
            KERR(1) = KNOM
            KERR(2) = 'ERREUR: OBJET INCONNU'
            CALL LEREUR
            IERR = IERR + 1
            GOTO 200
         ENDIF
 200  CONTINUE
      IF( IERR .NE. 0 ) GOTO 9990
C
C     BOUCLE SUR LES OBJETS AUX LIMITES DE L'OBJET
C     RECHERCHE DU NUMERO MIN ET MAX DES LIGNES
C     ============================================
      NLMIN = 1 000 000
      NLMAX = -1
      MN = MNOBCL - 2
      DO 210 I=1,NBOBCL
C        LE TYPE DE L'OBJET
         MN   = MN + 2
         NYOB = MCN( MN )
         IF( NYOB .NE. 2 ) GOTO 210
         NUOB = MCN( MN + 1 )
C        OUVERTURE DE L'OBJET
         CALL LXNLOU( NTMN(NYOB) , NUOB , NTOB , MNOB )
         IF( NTOB .LE. 0 ) THEN
            KNM = NMTYOB( NYOB )
            CALL NMOBNU( KNM , NUOB , KNOM )
            NBLGRC(NRERR) = 1
            KERR(1) = 'ERREUR :'//KNM//' INCONNU ' // KNOM
            CALL LEREUR
            IERR = IERR + 1
            GOTO 210
         ENDIF
         NLMIN = MIN( NLMIN, NUOB )
         NLMAX = MAX( NLMAX, NUOB )
 210  CONTINUE
      IF( NLMAX .LT. NLMIN ) THEN
         NLMAX = 0
         NLMIN = 0
      ENDIF
C
C     LES TABLEAUX CODE ET CARACTERISTIQUES NEKTON DES LIGNES
C     =======================================================
      NBLG = NLMAX - NLMIN + 1
C     LE CODE DES COURBES
      CALL TNMCDC( 'ENTIER', NBLG, MNCOLG )
C     PAR DEFAUT LE CODE NEKTON DE LA LIGNE EST ' '
      DO 215 MN = MNCOLG, MNCOLG-1+NBLG
         MCN( MN ) = ICHAR( ' ' )
 215  CONTINUE
C     LES CARACTERISTIQUES DES COURBES
      CALL TNMCDC( 'REEL', 6*NBLG, MNCALG )
      CALL AZEROR( 6*NBLG, RMCN(MNCALG) )
C     LA PILE DES NUMEROS DE FONCTIONS
      CALL TNMCDC( 'ENTIER', 3*NBLG, MNPIFT )
      CALL AZEROI( 3*NBLG, MCN(MNPIFT) )
      NBPIFT = 0
C
C     CREATION DES CODES ET CARACTERISTIQUES NEKTON DES LIGNES
C     ========================================================
      MNOB = MNOBCL - 2
      DO 280 I=1,NBOBCL
C        LE TYPE DE L'OBJET
         MNOB = MNOB + 2
         NYOB = MCN( MNOB )
         IF( NYOB .NE. 2 ) GOTO 280
         NUOB = MCN( MNOB + 1 )
         IF( NUOB .LE. 0 ) GOTO 280
C        OUVERTURE DE L'OBJET
         CALL LXNLOU( NTMN(NYOB) , NUOB , NTOB , MN )
         KNM = NMTYOB( NYOB )
         CALL NMOBNU( KNM , NUOB , KNOM )
C        OUVERTURE DU TABLEAU DEFINITION
         CALL LXTSOU( NTOB , 'DEFINITION' , NTDF, MNDF )
         IF( NTDF .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = KNOM
            KERR(2) ='ERREUR:'//KNM//' NON DEFINIE '
            CALL LEREUR
            IERR = IERR + 1
            GOTO 280
         ELSE
C           LES CARACTERISTIQUES DE LA LIGNE
            NUTYLI = MCN( MNDF + WUTYLI )
C
C           CF ~/td/d/a_ligne__definition
CCC   deftms ~>LIGNE>>DEFINITION 'la definition d''une ligne'
CCC   variable NTYTRL 'transformation (I pour IDENTITE)' ^~>TRANSFO ;
CCC   variable NUTYLI 'numero du type de la ligne'    entier
CCC   ( 0  : 'liste de points definis par 3 coordonnees'
CCC   , 1  : 'liste de points nommes par l''utilisateur'
CCC   , 2  : 'segment de droite'
CCC   , 3  : 'arc de cercle'
CCC   , 4  : 'cercle'
CCC   , 10 : 'donnee des tableaux XYZSOMMET et NSEF'
CCC   , 11 : 'B-spline ouverte uniforme d''interpolation'
CCC   , 12 : 'B-spline ouverte uniforme polynomiale'
CCC   , 13 : 'B-spline ouverte non uniforme polynomiale'
CCC   , 14 : 'B-spline ouverte non uniforme rationnelle'
CCC   , 21 : 'B-spline fermee  uniforme d''interpolation'
CCC   , 22 : 'B-spline fermee  uniforme polynomiale'
CCC   , 23 : 'B-spline fermee  non uniforme polynomiale'
CCC   , 24 : 'B-spline fermee  non uniforme rationnelle'
CCC   , 28 : 'ligne avec conge'
CCC   , 35 : 'extraction d''aretes d''une ligne'
CCC   , 36 : 'extraction d''une sous ligne sur critere'
CCC   , 37 : 'extraction d''une ligne d''une surface'
CCC   , 40 : 'restructuration d''une ligne non structuree'
CCC   , 41 : 'destructuration d''une ligne structuree'
CCC   , 45 : 'intersection 2 cylindres axes Z et (X ou Y)'
CCC   , 50 : 'transformation d''une ligne existante'
CCC   , 51 : 'union de plusieurs lignes'
CCC   , 81 : 'renommer une ligne'
CCC   , 82 : 'tuer cette ligne' );
C
C           LE CODE NEKTON DE LA LIGNE ( PAR DEFAUT => LIGNE DROITE )
C           ---------------------------------------------------------
            MN = MNCOLG + NUOB - NLMIN
            IF( NUTYLI.EQ.3 .OR. NUTYLI.EQ.4 .OR. NUTYLI.EQ.8 ) THEN
C
C              ARC DE CERCLE OU CERCLE
CCC  3: variable NBARLI 'nombre d''aretes de la ligne' entier ;
CCC     variable RAIGEO 'raison geometrique espacant les sommets' reel ;
CCC     variable NUPTIN 'nom point initial de la ligne' ^~>POINT ;
CCC     variable NUPTFI 'nom point final   de la ligne' ^~>POINT ;
CCC     variable NUPTEN 'nom point sur l''arc et entre les 2 extremites' ^~>POINT ;

CCC  4: variable NBARLI 'nombre d''aretes de la ligne' entier ;
CCC     variable RAIGEO 'raison geometrique espacant les sommets' reel ;
CCC     variable NUPTIN 'nom point initial de l''arc' ^~>POINT ;
CCC     variable NUPTFI 'nom point final   de l''arc' ^~>POINT ;
CCC     variable NUPTAR 'nom point du plan et plus proche du centre' ^~>POINT ;
CCC     variable RAYARC 'rayon du cercle passant par les 2 points' reel ;

CCC  8: variable NBARLI 'nombre d''aretes de la ligne' entier ;  
CCC     variable NUTYCI 'numero du type du cercle' entier
CCC  ( 1: 'cercle defini par 3 points sur le cercle'
CCC  , 2: 'cercle defini par le centre, un point du cercle, un point du plan'
CCC  , 3: 'cercle defini par le centre, le rayon, le plan a X ou Y ou Z=Cte' ) ;
CCC  cas NUTYCI 
CCC   1: variable NUPT1C 'nom 1-er  point sur le cercle' ^~>POINT ;
CCC      variable NUPT2C 'nom 2-eme point sur le cercle' ^~>POINT ;
CCC      variable NUPT3C 'nom 3-eme point sur le cercle' ^~>POINT ;
CCC   2: variable NUPTCE 'nom du point centre du cercle' ^~>POINT ;
CCC      variable NUPTSC 'nom du point sur le cercle'    ^~>POINT ;
CCC      variable NUPDPL 'nom du point pour definir le plan du cercle' ^~>POINT ;
CCC   3: variable NUPTCE 'nom du point centre du cercle' ^~>POINT ;
CCC      variable RAYDCI 'rayon du cercle' reel ;
CCC      variable NUPLCT 'numero du plan X ou Y ou Z=Cte' entier
CCC    ( 1: 'cercle dans le plan YZ avec X=centre(1)=Cte'
CCC    , 2: 'cercle dans le plan ZX avec Y=centre(2)=Cte'
CCC    , 3: 'cercle dans le plan XY avec Z=centre(3)=Cte' ) ;
CCC  fincas ; 

               IF( NUTYLI .EQ. 3 ) THEN

C                 ARC DE CERCLE
C                 LES 2 POINTS EXTREMITES ET LE POINT INTERMEDIAIRE
                  CALL NPTXYZ( MCN( MNDF + WUPTIN ), MNXYZ1, N )
                  IERR = IERR + N
                  CALL NPTXYZ( MCN( MNDF + WUPTEN ), MNXYZ2, N )
                  IERR = IERR + N
                  CALL NPTXYZ( MCN( MNDF + WUPTFI ), MNXYZ3, N )
                  IERR = IERR + N
                  IF( IERR .GT. 0 ) GOTO 9990
C                 CALCUL DU RAYON ET DU CENTRE DU CERCLE PASSANT PAR CES 3 POINTS
                  CALL RACE3P( RMCN(MNXYZ1), RMCN(MNXYZ2), RMCN(MNXYZ3),
     %                         RAYON, CENTRE, IERR )

               ELSE IF( NUTYLI .EQ. 4 ) THEN

C                 ARC DE CERCLE
C                 LES 2 POINTS EXTREMITES ET LE POINT DU PLAN DU CERCLE
                  CALL NPTXYZ( MCN( MNDF + WUPTIN ), MNXYZ1, N )
                  IERR = IERR + N
                  CALL NPTXYZ( MCN( MNDF + WUPTAR ), MNXYZ2, N )
                  IERR = IERR + N
                  CALL NPTXYZ( MCN( MNDF + WUPTFI ), MNXYZ3, N )
                  IERR = IERR + N
                  IF( IERR .GT. 0 ) GOTO 9990
C                 CALCUL DU RAYON ET DU CENTRE DU CERCLE PASSANT PAR CES 3 POINTS
                  CALL RACE3P( RMCN(MNXYZ1), RMCN(MNXYZ2), RMCN(MNXYZ3),
     %                         RAYON, CENTRE, IERR )

               ELSE IF( NUTYLI .EQ. 8 ) THEN

C                 CERCLE
                  NUTYCI = MCN( MNDF + WUTYCI )
                  IF( NUTYCI .EQ. 1 ) THEN

C                    CERCLE DEFINI PAR 3 POINTS SUR LE CERCLE
                     NUPT1C = MCN( MNDF + WUPT1C )
                     NUPT2C = MCN( MNDF + WUPT2C )
                     NUPT3C = MCN( MNDF + WUPT3C )
                     CALL NPTXYZ( MCN( MNDF + WUPT1C ), MNXYZ1, N )
                     IERR = IERR + N
                     CALL NPTXYZ( MCN( MNDF + WUPT2C ), MNXYZ2, N )
                     IERR = IERR + N
                     CALL NPTXYZ( MCN( MNDF + WUPT3C ), MNXYZ3, N )
                     IERR = IERR + N
                     IF( IERR .GT. 0 ) GOTO 9990
C                    CALCUL DU RAYON ET DU CENTRE DU CERCLE PASSANT PAR CES 3 PO
                     CALL RACE3P(RMCN(MNXYZ1),RMCN(MNXYZ2),RMCN(MNXYZ3),
     %                           RAYON, CENTRE, IERR )
                     IF( IERR .GT. 0 ) GOTO 9990

                  ELSE IF( NUTYCI .EQ. 2 ) THEN

C                    CERCLE DEFINI PAR SON CENTRE ET UN POINT DU CERCLE
                     CALL NPTXYZ( MCN( MNDF + WUPTCE ), MNXYZ1, N )
                     IERR   = IERR + N
                     CALL NPTXYZ( MCN( MNDF + WUPTSC ), MNXYZ2, N )
                     IERR   = IERR + N
                     IF( IERR .GT. 0 ) GOTO 9990
                     CENTRE(1) = RMCN(MNXYZ1)
                     CENTRE(2) = RMCN(MNXYZ1+1)
                     RAYON     = DIST2P( RMCN(MNXYZ1), RMCN(MNXYZ2) )

                  ELSE

C                    CERCLE DEFINI PAR LE CENTRE, LE RAYON, LE PLAN A X OU Y OU Z=CTE
                     CALL NPTXYZ( MCN( MNDF + WUPTCE ), MNXYZ1, N )
                     CENTRE(1) = RMCN(MNXYZ1)
                     CENTRE(2) = RMCN(MNXYZ1+1)
                     RAYON = RMCN( MNDF + WAYDCI )

                  ENDIF

               ENDIF
C
C              ARC DE CERCLE OU CERCLE
               MCN( MN ) = ICHAR( 'C' )
C              CARACTERISTIQUE NEKTON => SON RAYON ET SON CENTRE(1:2)
               MNC = MNCALG + 6 * ( NUOB - NLMIN )
               RMCN( MNC   ) = RAYON
               RMCN( MNC+1 ) = CENTRE(1)
               RMCN( MNC+2 ) = CENTRE(2)
               CALL AZEROR( 3, RMCN(MNC+3) )
C
            ELSE IF( NUTYLI .GE. 11 .AND. NUTYLI .LE. 24 ) THEN
C
C              B-SPLINE
               NBLGRC(NRERR) = 2
               KERR(1) = KNOM
               KERR(2) ='ERREUR: SPLINE NON INTEGRE A NEKTON'
               CALL LEREUR
CCC               MCN( MN ) = ICHAR( 'S' )
CCCC              CARACTERISTIQUE NEKTON => ???
CCC               MNC = MNCALG + 6 * ( NUOB - NLMIN )
CCC               CALL AZEROR( 6, RMCN(MNC) )
            ENDIF
C
C           LA TRANSFORMATION DE LA LIGNE
            NTYTRL = MCN( MNDF + WTYTRL )
            IF( NTYTRL .GT. 1 ) THEN
C              LIGNE SUBISSANT LA TRANSFORMATION NTYTRL
C              CARACTERISTIQUE NEKTON => 3 NUMEROS DES FONCTIONS
               CALL LXNLOU( NTTRAN, NTYTRL, NTP, MNP )
               IF( NTP .LE. 0 ) THEN
                  WRITE(KERR(MXLGER)(1:10),'(I10)') NTYTRL
                  NBLGRC(NRERR) = 1
                  KERR(1)='ERREUR: TRANSFORMATION INCONNUE NUMERO '
     %                  // KERR(MXLGER)(1:10)
                  CALL LEREUR
                  IERR = 5
                  GOTO 9990
               ENDIF
CCC    deftms ~>TRANSFO>>DEFINITION
CCC      'la definition d''une transformation mathematique'
CCC    variable NUTYTR 'numero du type d''une transformation'    entier
CCC         ( 1  : 'identite de R3 -> R3'
CCC         , 2  : 'transformation de R3 -> R3 definie par 3 fonctions'
CCC         , 3  : 'isometrie de R3 -> R3'
CCC         , 4  : 'projection dans R3' ) ;
CCC    cas NUTYTR
CCC    2 : variable NUFONX 'nom fonction X=F1(x,y,z)' ^~>FONCTION ;
CCC        variable NUFONY 'nom fonction Y=F2(x,y,z)' ^~>FONCTION ;
CCC        variable NUFONZ 'nom fonction Z=F3(x,y,z)' ^~>FONCTION ;
C
               CALL LXTSOU( NTP, 'DEFINITION', NTXYZ, MNXYZ )
               IF( NTXYZ .LE. 0 ) THEN
                  CALL NMOBNU( 'TRANSFO' , NTYTRL , KERR(3) )
                  LL = NUDCNB(KERR(3))
                  NBLGRC(NRERR) = 2
                  KERR(1) = 'TRANSFORMATION : ' // KERR(3)(1:LL)
                  KERR(2) = 'ERREUR: DEFINITION INCONNUE'
                  CALL LEREUR
                  IERR = 5
                  GOTO 9990
               ENDIF
C              LE TYPE DE LA TRANSFORMATION
               NUTYTR = MCN(MNXYZ+WUTYTR)
               IF( NUTYTR .NE. 2 ) THEN
                  CALL NMOBNU( 'TRANSFO' , NTYTRL , KERR(3) )
                  LL = NUDCNB(KERR(3))
                  NBLGRC(NRERR) = 2
                  KERR(1) = 'TRANSFORMATION : ' // KERR(3)(1:LL)
                  KERR(2) = 'ERREUR: TYPE DIFFERENT DE 2 INTERDIT'
                  CALL LEREUR
                  IERR = 5
                  GOTO 9990
               ENDIF
C              LE NUMERO DES FONCTIONS EST EMPILE
               DO 250 N=0,2
                  NF = MCN( MNXYZ + WUFONX + N )
                  DO 240 J=1,NBPIFT
                     IF( NF .EQ. MCN(MNPIFT-1+J) ) GOTO 250
 240              CONTINUE
                  MCN( MNPIFT + NBPIFT ) = NF
                  NBPIFT = NBPIFT + 1
 250           CONTINUE
C
C              LIGNE SUBISSANT LA TRANSFORMATION NTYTRL COMPATIBLE NEKTON
               MCN( MN ) = ICHAR( 'F' )
C              LE NUMERO DES 3 FONCTIONS DANS LE LEXIQUE DES FONCTIONS
               MNC = MNCALG + 6 * ( NUOB - NLMIN )
               RMCN( MNC   ) = MCN( MNXYZ + WUFONX )
               RMCN( MNC+1 ) = MCN( MNXYZ + WUFONY )
               RMCN( MNC+2 ) = MCN( MNXYZ + WUFONZ )
               CALL AZEROR( 3, RMCN(MNC+3) )
            ENDIF
         ENDIF
 280  CONTINUE
      IF( IERR .NE. 0 ) GOTO 9990
C
C     CALCUL DU NOMBRE D'EF QUADRANGLES
      NBEF = 0
      DO 300 NOTYEL = 1 , NBTYEL
C        L'ADRESSE DU TABLEAU NPEF
         MNELE = MCN( MNELEM - 1 + NOTYEL )
C        LE NOMBRE D'ELEMENTS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
         NBEF  = NBEF + NBELEM
 300  CONTINUE
      WRITE(IMPRIM,10300) NBEF
10300 FORMAT(' NOMBRE DE QUADRANGLES',T32,'=',I6)
C
C     DECLARATION DES TABLEAUX POUR UN EF
C     ===================================
      CALL TNMCDC( 'REEL', 8*3, MNXYZ )
      CALL AZEROR( 8*3, RMCN(MNXYZ) )
      MNX = MNXYZ - 1
      MNY = MNX + 8
      MNZ = MNY + 8
C
C     DECLARATION DES TABLEAUX NEKTON
C     ===============================
      CALL TNMCDC( 'REEL', NBEF*8, MNXX )
      CALL TNMCDC( 'REEL', NBEF*8, MNYY )
      CALL TNMCDC( 'REEL', NBEF*8, MNZZ )
C     LE CODAGE DES COURBES
      CALL TNMCDC( 'ENTIER', 8*NBEF, MNCCUR )
      DO 350 MN=MNCCUR,MNCCUR+8*NBEF-1
C        PAR DEFAUT LA COURBE EST UNE DROITE
         MCN(MN) = ICHAR( ' ' )
 350  CONTINUE
C     LES CARACTERISTIQUES DES COURBES
      CALL TNMCDC( 'REEL', 6*8*NBEF, MNCURV )
      CALL AZEROR( 6*8*NBEF, RMCN(MNCURV) )
C
C     LA BOUCLE DE REMPLISSAGE DUES TABLEAUX 'NEKTON' EF PAR EF
C     =========================================================
      NUEF = 0
      DO 500 NOTYEL = 1 , NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF
         MNELE = MCN( MNELEM - 1 + NOTYEL )
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        LE NOMBRE D'ELEMENTS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES ELEMENTS
         MNNDEL = MNELE + WUNDEL
         MNPGEL = MNNDEL
         IF( NDPGST .GE. 2 ) THEN
            MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
         ENDIF
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL , NOMELE )
         CALL ELTYCA( NUTYEL )
         CALL AZEROI( 12, NOOBLA )
         CALL AZEROI(  6, NOOBSF )
C
         IF( NBPOE .NE. 4 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'EF DIFFERENT D''UN QUADRANGLE'
            CALL LEREUR
            GOTO 9990
         ENDIF
C
C        LES ADRESSAGES POUR RECUPERER LES INFORMATIONS
C        SELON LE TYPE DE L'ELEMENT
         GOTO(401,401,401,401,400,400,400,400,400,400,
     &        400,400,401,400,401,401,400,401,400,400,
     &        400,400,400,400)                       , NUTYEL
C
C        ERREUR
 400     NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR: ELEMENT '// NOMELE(1)
     &             // NOMELE(2) //' NON PROGRAMME'
         CALL LEREUR
         GOTO 9990
C
C        LA BOUCLE SUR LES EF DE CE TYPE NUTYEL
 401     DO 495 NUELEM = 1 , NBELEM
C
C           LE NUMERO DE L'EF
            NUEF = NUEF + 1
C
C           LES NOEUDS DE L'ELEMENT
C           -----------------------
CCC            CALL EFNOEU( MNELE , NUELEM , NBNDEL , MCN(IANOEF) )
C
C           LES POINTS GEOMETRIQUES DE L'ELEMENT
C           ------------------------------------
CCC            CALL EFPOGE( MNELE , NUELEM , NBPGEF , MCN(IAPOEF) )
C
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
C           ----------------------------------------
            CALL EFPLSV( MNELE  , NUELEM ,
     %                   NOVCEL , NOSFEL , NOLAEL , NOPSEL ,
     %                   NOOBVC , NOOBSF , NOOBLA , NOOBPS )
C
C           LES COORDONNEES DES POINTS DE L'ELEMENT
C           ---------------------------------------
            DO 440 I=1,NBPOE
C              LE NUMERO DU I-EME POINT DE L'ELEMENT NUELEM
               N   = MCN( MNPGEL-1 + NUELEM + NBELEM*(I-1) )
               MNC = MNPOGE + WYZPOI + (N-1) * 3
               RMCN( MNX + I ) = RMCN( MNC     )
               RMCN( MNY + I ) = RMCN( MNC + 1 )
               RMCN( MNZ + I ) = RMCN( MNC + 2 )
C              LE NO DE LA LIGNE DU COTE I DU QUADRANGLE EST DANS NOOBLA
 440        CONTINUE
C
C           REMPLISSAGE DES TABLEAUX NEKTON
C           X PUIS Y PUIS Z DES SOMMETS DES QUADRANGLES
            DO 460 I=1,8
C              L'ABSCISSE, L'ORDONNEE ET LA COTE Z  ( PREVUE POUR LE 3D )
               RMCN( MNXX - 1 + NUEF + (I-1)*NBEF ) = RMCN( MNX + I )
               RMCN( MNYY - 1 + NUEF + (I-1)*NBEF ) = RMCN( MNY + I )
               RMCN( MNZZ - 1 + NUEF + (I-1)*NBEF ) = RMCN( MNZ + I )
C              LE CODE DE COURBE ( NO LIGNE = NOOBLA(I) )
               IF( NOOBLA(I) .GT. 0 ) THEN
                  MN = MNCOLG + NOOBLA(I) - NLMIN
                  MCN( MNCCUR - 1 + I + 8 * (NUEF-1) ) = MCN( MN )
C                 LES CARACTERISTIQUES DE LA COURBE
                  MN =  MNCALG + 6 * ( NOOBLA(I) - NLMIN )
                  MNC = MNCURV + 6 * ( I - 1 ) + 48 * ( NUEF - 1 )
                  DO 450 K=0,5
                     RMCN(MNC+K) = RMCN(MN+K)
 450              CONTINUE
               ENDIF
 460        CONTINUE
 495     CONTINUE
 500  CONTINUE
C
C     DECLARATION OUVERTURE DU FICHIER NEKTON
C     =======================================
      INQUIRE( FILE='NEKTON', EXIST=LEXIST, OPENED=LOPEN )
      IF( LEXIST ) THEN
C        LE FICHIER 'NEKTON' EXISTE
         IF( .NOT. LOPEN ) THEN
C           OUVERTURE DU FICHIER 'NEKTON'
            CALL TRUNIT( NFNEKT )
            OPEN( FILE='NEKTON', UNIT=NFNEKT, STATUS='OLD' )
         ENDIF
C        DESTRUCTION DU FICHIER
         CLOSE( NFNEKT, STATUS='DELETE' )
      ENDIF
C
C     CREATION DU FICHIER 'NEKTON'
      CALL TRUNIT( NFNEKT )
      OPEN( UNIT=NFNEKT, ERR=9995, STATUS='NEW',
     %      FILE='NEKTON' , ACCESS='SEQUENTIAL' , FORM='FORMATTED' )
C     LE NOMBRE D'EF
      WRITE(NFNEKT,'(I10)') NBEF
      WRITE(NFNEKT,'(5E15.7)') (RMCN(MN),MN=MNXX,MNXX+8*NBEF-1)
      WRITE(NFNEKT,'(5E15.7)') (RMCN(MN),MN=MNYY,MNYY+8*NBEF-1)
      WRITE(NFNEKT,'(5E15.7)') (RMCN(MN),MN=MNZZ,MNZZ+8*NBEF-1)
C     LES CODES DES COURBES
      WRITE(NFNEKT, '(50(A1,1X))' )
     %                     (CHAR(MCN(MN)),MN=MNCCUR,MNCCUR+8*NBEF-1)
C     LES CARACTERISTIQUES DES COURBES
      WRITE(NFNEKT,'(5E15.7)') (RMCN(MN),MN=MNCURV,MNCURV+48*NBEF-1)
C
C     LE TEXTE DES FONCTIONS
      WRITE(KFRMLG(1:5),10700) NBCALI
10700 FORMAT('(A',I2,')')
      DO 750 I=1,NBPIFT
C        LE NUMERO DE LA FONCTION DANS LE LEXIQUE DES FONCTIONS
         NF = MCN(MNPIFT - 1 + I)
         WRITE(IMPRIM,10701) NF
10701    FORMAT(/' FONCTION NUMERO ',I4)
         CALL LXNLOU( NTFONC, NF, NTOB , MNOB )
         CALL NMOBNU( 'FONCTION', NF, KNOM )
C        OUVERTURE DU TABLEAU FONCTION>NOM>TEXTE
         CALL LXTSOU( NTOB , 'TEXTE' , NTDF, MNDF )
         IF( NTDF .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = KNOM
            KERR(2) ='ERREUR: FONCTION DE TEXTE NON DEFINI'
            CALL LEREUR
            IERR = IERR + 1
         ENDIF
         NL = MCN(MNDF+WBLITX)
         WRITE(NFNEKT,'(2I10)') NF,NL
C        LE NOMBRE DE MOTS POUR STOCKER UNE LIGNE DE NBCALI (KLG) CARACTERES
         MOT = MCN(MNDF+WBCLTX)
         MN  = MNDF + WECATX
         DO 710 J=1,NL
            NC1 = 0
            MN0 = MN
            DO 705 K=1,MOT
C              ICI: NBCAMO CARACTERES PAR MOT DE MCN
               NC2 = NC1 + NBCAMO
               KLIG(NC1+1:NC2) = CHARX( MCN(MN) )
               NC1 = NC2
               MN  = MN + 1
 705        CONTINUE
            WRITE(NFNEKT,KFRMLG) KLIG
            WRITE(IMPRIM, '(A)') KLIG
 710     CONTINUE
 750  CONTINUE
      CLOSE( NFNEKT )
C
C     DESTRUCTION DES TABLEAUX INUTILES
C     ---------------------------------
 9990 IF( MNELEM .GT. 0 ) CALL TNMCDS( 'ENTIER' , 2*MXTYEL , MNELEM )
      IF( MNCOLG .GT. 0 ) CALL TNMCDS( 'ENTIER', NBLG, MNCOLG )
      IF( MNCALG .GT. 0 ) CALL TNMCDS( 'REEL', 6*NBLG, MNCALG )
      IF( MNXYZ  .GT. 0 ) CALL TNMCDS( 'REEL', 8*3, MNXYZ )
      IF( MNXX   .GT. 0 ) CALL TNMCDS( 'REEL', NBEF*8, MNXX )
      IF( MNYY   .GT. 0 ) CALL TNMCDS( 'REEL', NBEF*8, MNYY )
      IF( MNZZ   .GT. 0 ) CALL TNMCDS( 'REEL', NBEF*8, MNZZ )
      IF( MNCCUR .GT. 0 ) CALL TNMCDS( 'ENTIER', 8*NBEF, MNCCUR )
      IF( MNCURV .GT. 0 ) CALL TNMCDS( 'REEL', 6*8*NBEF, MNCURV )
      IF( MNPIFT .GT. 0 ) CALL TNMCDS( 'REEL', 3*NBLG, MNPIFT )
      RETURN
C
C     PROBLEME SUR LE FICHIER
 9995 WRITE(IMPRIM,*) 'FICHIER NEKTON NON OUVRABLE'
      GOTO 9990
      END
