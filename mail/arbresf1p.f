      SUBROUTINE ARBRESF1P( NBZCER, NUPCZC, XYZCZC,RAYZC,NBARLC,
     %                      NULGZC, MNXYZC0,
     %                      NOPIEC, NBZCP, NUZCPI, NCTRIZ,
     %                      NOPTIM, MXPTIM, NBPTIM, XYZDPTIM,
     %                      NUSFTRPI, NOTYPE, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   GENERER UNE TRIANGULATION DE LA SURFACE D'UNE PIECE D'UN ARBRE
C -----   DEFINIE PAR LA LISTE DE SES Z-CERCLES (dans un PLAN Z=Cte)
C         ATTENTION: UN Z-CERCLE PEUT ETRE REDUIT A UN POINT SON CENTRE
C                    DANS CE CAS LE NOMBRE DE SOMMETS DOIT ETRE EGAL A 0
C         LE TYPE D'UNE PIECE PEUT ETRE
C         1: CONE SOMMET BAS  (defini par 1 Pt=ZC + 1 ZC)
C         2: CONE SOMMET HAUT (defini par 1 ZC + 1 Pt=ZC)
C         3: TRONC DE CONE    (defini par 2 ZC)
C         4: n Z-CERCLES RELIES a UN Z-CERCLE  (defini par n + 1 ZC)
C         5: UN Z-CERCLE RELIE  a n  Z-CERCLES (defini par 1 + n ZC)

C ENTREES:
C --------
C NOFOTI : NO DANS LE LX DES FONCTIONS DE LA FONCTION TAILLE_IDEALE(X,Y,Z)
C NBZCER : NOMBRE DE Z-CERCLES DE L'ARBRE
C NUPCZC : NUMERO DANS LE LX POINTS DES NBZCER CENTRES DES CERCLES
C XYZCZC : XYZ DU CENTRE DES NBZCER Z-CERCLES
C RAYZC  : RAYON DES NBZCER Z-CERCLES
C NBARLC : >2 NOMBRE D'ARETES DE CHACUN DES NBZCER CERCLES
C          =0 INDIQUE UN CERCLE REDUIT A SON POINT CENTRE
C NULGZC : NUMERO DANS LE LEXIQUE LIGNES DES NBZCER CERCLES
C MNXYZC0: ADRESSE MCN DU TABLEAU DE L'ADRESSE MNXYZS du TMS XYZSOMMET
C          DES NBZCER LIGNES Z-CERCLES

C NOPIEC : NUMERO DE LA PIECE DE L'ARBRE
C NBZCP  : NOMBRE DE Z-CERCLES DE LA PIECE NOPIEC
C NUZCPI : NUMERO 1 a NBZCP du Z-CERCLE de chaque PIECE

C NOPTIM : =0 PAS D'AJOUT DE POINTS A IMPOSER DANS UNE FUTURE TETRAEDRISATION
C          =1       AJOUT DE POINTS A IMPOSER DANS UNE FUTURE TETRAEDRISATION
C MXPTIM  : MAXIMUM DE POINTS A TETRAEDRISER STOCKABLES DANS XYZDPTIM
C XYZDPTIM: X Y Z Distance Souhaitee aux sommets voisins DES POINTS
C           INITIAUX A TETRAEDRISER ENSUITE

C AUXILIAIRES:
C ------------
C NCTRIZ : NUMERO DE 1 A NBZCER DES Z-CERCLES DE LA PIECE NOPIEC
C          A TRIER SELON Z CROISSANT DU CENTRE DE LEUR CERCLE

C SORTIES :
C ---------
C NBPTIM  : NOMBRE DE POINTS A TETRAEDRISER ENSUITE
C XYZDPTIM: X Y Z Distance Souhaitee aux sommets voisins DES POINTS
C           A TETRAEDRISER ENSUITE
C NUSFTRPI: NUMERO DANS LE LEXIQUE SURFACES DE LA SURFACE TRIANGULEE
C           DE LA PIECE NOPIEC DE L'ARBRE
C NOTYPE  : LE TYPE DE LA PIECE est NOTYPE
C           =1: CONE SOMMET BAS
C           =2: CONE SOMMET HAUT
C           =3: TRONC DE CONE
C           =4: n Z-CERCLES RELIES a UN Z-CERCLE
C           =5: UN Z-CERCLE RELIE  a n  Z-CERCLES
C            puis TRIANGULATION DE SA SURFACE LIMITEE PAR NBZCP Z-CERCLES
C IERR    : =0 SI PAS D'ERREUR DETECTEE, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  Saint PIERRE du PERRAY             Juin 2019
C23456...............................................................012
      PARAMETER      (MXZCER=256)
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_point__definition.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/darete.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"

      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))

      REAL             RAYZC(NBZCER),  XYZCZC(3,NBZCER),
     %                 XYZDPTIM(4,MXPTIM)
      INTEGER          NUPCZC(NBZCER), NBARLC(NBZCER), NULGZC(NBZCER)
      INTEGER          NUZCPI(NBZCP), NCTRIZ(NBZCP)

      INTEGER          NULGZCP(MXZCER)
      DOUBLE PRECISION XYZCZCP(3,MXZCER), RAYZCP(MXZCER), TAIDAR(MXZCER)
      REAL             AMPLIRC(0:MXZCER)

      CHARACTER*8      KNC
      CHARACTER*24     KNMLGZC, KNMPISF, KNM, KNMLGZCP, KNMTRIA

      DOUBLE PRECISION PI, PIS3, PIS6, DECANGL, XYZBAR(4), RAYONK,
     %                 R, R1, R2, RR2, RAPP, RAPP1, RCPMIN, TAR,
     %                 D, DMR, DMRMIN, D4, DCPMIN, FR, FR2,
     %                 ZH, ZHTE, ZHSCE, ANGLE, ANGLED, COSANG, SINANG,
     %                 XYZD(4), XYZDR(4), XYZ(3), XYZSCE(3),
     %                 XYZ1(3), XYZ2(3), XYZ4(3), XYZ5(3),
     %                 P1PLAN(3), P2PLAN(3), P3PLAN(3)

      IERR    = 0
      Pi      = ATAN( 1D0 ) * 4D0
      PIS3    = Pi / 3D0
      PIS6    = Pi / 6D0
      MNXYZCP = 0
      MNFRST  = 0
      NBPTIM  = 0
      NOTYPE  = 0

      print*
      IF( NBZCP .LT. 2 ) THEN
         PRINT*,'arbresf1p:',NBZCP,
     %          '<2 NOMBRE INCORRECT DE Z-CERCLES de la PIECE',NOPIEC
         IERR = 2
         GOTO 9900
      ENDIF

      PRINT*,'arbresf1p: TRIANGULATION de la SURFACE de la PIECE',NOPIEC
     %,' a partir de ses',NBZCP,' Z-CERCLES'
      DO N=1,NBZCP
C        NOM DU CENTRE DU N-ieme Z-CERCLE de la PIECE NOPIEC
         NOZCPI =  NUZCPI( N )
         CALL NMOBNU( 'POINT',  NUPCZC( NOZCPI ), KNM )
         PRINT*,'arbresf1p: CENTRE du Z-CERCLE: ',KNM,
     %          ' de RAYON=',RAYZC( NOZCPI ),' de XYZ=',
     %          (XYZCZC(K,NOZCPI),K=1,3)
      ENDDO

C     RECHERCHE DU TYPE NOTYPE DE LA PIECE
C     1: CONE SOMMET BAS
C     2: CONE SOMMET HAUT
C     3: TRONC DE CONE
C     4: n Z-CERCLES RELIES a UN Z-CERCLE
C     5: UN Z-CERCLE RELIE  a n  Z-CERCLES
C     ------------------------------------

C     TABLEAU DU NO 1 A NBZCER DES NBZCP Z-CERCLES DE LA PIECE
      DO N = 1, NBZCP
         NCTRIZ( N ) = NUZCPI( N )
      ENDDO

C     RANGEMENT CROISSANT SELON LA COTE Z DES Z-CERCLES
      DO N = 1, NBZCP-1
         N1ZC = NCTRIZ( N )
         ZMIN = XYZCZC( 3, N1ZC )
         NMIN = N
         DO N2 = N+1, NBZCP
            N2ZC = NCTRIZ( N2 )
            Z2   = XYZCZC( 3, N2ZC )
            IF( Z2 .LT. ZMIN ) THEN
               ZMIN = Z2
               NMIN = N2
            ENDIF
         ENDDO
         IF( NMIN .GT. N ) THEN
C           PERMUTATION N1ZC et N2ZC
            NCTRIZ( N    ) = NCTRIZ( NMIN )
            NCTRIZ( NMIN ) = N1ZC
         ENDIF
      ENDDO

C     RECHERCHE DU TYPE DE LA PIECE NOPIEC SELON LE NOMBRE DE Z-CERCLES
      IF( NBZCP .EQ. 2 ) THEN

         IF( NBARLC( NCTRIZ(1) ) .LE. 0 ) THEN

            IF( NBARLC( NCTRIZ(2) ) .LE. 0 ) THEN
C              ERREUR: JONCTION ENTRE 2 CERCLES REDUITS A LEUR CENTRE
               PRINT *,'arbresf1p: PIECE',NOPIEC,
     %        ' Z-CERCLE',NCTRIZ(1),' avec',NBARLC(NCTRIZ(1)),' ARETES',
     %        ' Z-CERCLE',NCTRIZ(2),' avec',NBARLC(NCTRIZ(2)),' ARETES',
     %        ' donc REDUITS A LEUR CENTRE'
               PRINT *,'=> PAS de CREATION d''UNE SURFACE d''UN SEGMENT'
               IERR = 3
               GOTO 9900
            ENDIF

C           PIECE DE TYPE 1: CONE SOMMET BAS
            NOTYPE = 1
            GOTO 11

         ENDIF

         IF( NBARLC( NCTRIZ( 2 ) ) .LE. 0 ) THEN
C           PIECE DE TYPE 2: CONE SOMMET HAUT
            NOTYPE = 2
            GOTO 12
         ENDIF

         IF( NBARLC( NCTRIZ( 1 ) ) .GT. 0  .AND.
     %       NBARLC( NCTRIZ( 2 ) ) .GT. 0 ) THEN
C            PIECE DE TYPE 3: TRONC DE CONE
             NOTYPE = 3
             GOTO 20
         ENDIF

         PRINT*,'arbresf1p: PIECE',NOPIEC,' COMPOSEE de',NBZCP,
     %          ' Z-CERCLES INCORRECTS'
         IERR = 4
         GOTO 9900

      ENDIF

C     LA PIECE NOPIEC A PLUS DE 2 Z-CERCLES
      GOTO 40


C     PIECE DE TYPE 1: CONE SOMMET BAS
C     ================================
 11   NCERCL = NCTRIZ( 2 )
      NCSTCO = NCTRIZ( 1 )
      GOTO 13

C     PIECE DE TYPE 2: CONE SOMMET HAUT
C     =================================
 12   NCERCL = NCTRIZ( 1 )
      NCSTCO = NCTRIZ( 2 )
         
C     CONE: LE Z-CERCLE NCSTCO EST REDUIT A SON CENTRE
C     JONCTION DU CENTRE NCSTCO AUX 2 EXTREMITES DES ARETES DU CERCLE NCERCL
C     ----------------------------------------------------------------------
C     CONSTRUCTION DU LEXIQUE de la SURFACE du CONE
C     NOM DE LA SURFACE DE CE CONE NOPIEC
 13   WRITE( KNC(1:8), '(I8)' ) NOPIEC
C     RETRAIT DES BLANCS
      CALL SANSBL( KNC, NBCAR )
      KNMPISF = 'CONE_' // KNC(1:NBCAR) // '_AD '

C     CONSTRUCTION DU TMS LEXIQUE de la SURFACE KNMPISF
C     SI CETTE SURFACE EXISTE, ELLE EST DETRUITE
      CALL LXLXOU( NTSURF, KNMPISF, NT1PISF, MN1PISF )
      IF(MN1PISF.GT.0) CALL LXTSDS( NTSURF, KNMPISF )
      CALL LXLXDC( NTSURF, KNMPISF, 24, 8 )
      CALL LXLXOU( NTSURF, KNMPISF, NT1PISF, MN1PISF )

C     NUMERO DE LA SURFACE DANS LE LEXIQUE DES SURFACES
      CALL NUOBNM( 'SURFACE', KNMPISF, NUSFTRPI )

C     CONSTRUCTION DU TMS a_surface__definition du CONE
C     LA SURFACE A UN TMS DEFINITION DE SURFACE DE TYPE 32
      CALL LXTNDC( NT1PISF, 'DEFINITION', 'MOTS', WULGJN+1 )
      CALL LXTSOU( NT1PISF, 'DEFINITION', NT1PISFD, MN1PISFD )

C     TRANSFORMATION (I pour IDENTITE)
      MCN( MN1PISFD + WTYTRS ) = 1
C     TYPE DE LA SURFACE 32: CONE  LIAISON 1SOMMET-> 1LIGNE
      MCN( MN1PISFD + WUTYSU ) = 32

C     NBARPL = nombre d''aretes point=CENTRE HAUT->ligne CERCLE NCERCL
C     DISTANCE ENTRE LES CENTRES DES 2 CERCLES
      D = ( XYZCZC(1,NCERCL) - XYZCZC(1,NCSTCO) ) **2
     %  + ( XYZCZC(2,NCERCL) - XYZCZC(2,NCSTCO) ) **2
     %  + ( XYZCZC(3,NCERCL) - XYZCZC(3,NCSTCO) ) **2
      D = SQRT( D )

C     LONGUEUR DE L'ARETE SUR LE CERCLE NCERCL
      R = 2 * Pi * RAYZC( NCERCL ) / NBARLC( NCERCL )

C     LONGUEUR DE L'ARETE SUR LA HAUTEUR DU CONE
ccc      TAR = MIN( DARETE, (DARETE+R)/2, 3*R )
      TAR = MIN( DARETE, (DARETE+R)/2, 2*R )

C     NOMBRE D'ARETES SUR LE SEGMENT JOIGNANT LES 2 CENTRES
      NBARPL = INT( D / TAR ) + 1
      MCN( MN1PISFD + WBARPL ) = NBARPL

C     NUPTJN = nom du point=CENTRE HAUT a joindre aux aretes de NCERCL
      MCN( MN1PISFD + WUPTJN ) = NUPCZC( NCSTCO )
C     NULGJN = nom de la ligne a joindre au point
      MCN( MN1PISFD + WULGJN ) = NULGZC( NCERCL )

C     LA DATE DE CREATION
      CALL ECDATE( MCN(MN1PISFD) )

C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN(MN1PISFD+MOTVAR(6))=NONMTD('~>SURFACE>>DEFINITION')

C     CONSTRUCTION DU MAILLAGE DE LA SURFACE CONE
      CALL SUEX32( NT1PISF, MCN(MN1PISFD),
     %             NTNSEFS, MNNSEFS, NTXYZS, MNXYZS, IERR )

      print*,'arbresf1p: TRIANGULATION de la SURFACE',NOPIEC,
     %       ' de NOM: ',KNMPISF,
     %       ' avec',MCN(MNNSEFS+WBEFOB),' TRIANGLES',
     %        MCN(MNXYZS+WNBSOM),' SOMMETS'

C     EVENTUELLE CREATION DES POINTS A IMPOSER DANS LA TETRAEDRISATION
C     LE LONG DE L'AXE JOIGNANT LE CENTRE DES 2 CERCLES BAS et HAUT
C     EN FONCTION de ARETGR
      IF( XYZCZC( 3,NCTRIZ(1) ) .GT. XYZCZC( 3,NCTRIZ(2) ) ) THEN
C        RACINES
         NCTRON = NCTRIZ( NBZCP )
         NDBRAN = 0
      ELSE
C        BRANCHES
         NCTRON = NCTRIZ( 1 )
         NDBRAN = 1
      ENDIF
      NBBRAN = NBZCP - 1

      GOTO 70


C     PIECE DE TYPE 3: TRONC DEFINI PAR 2 Z-CERCLES
C     =============================================
 20   NCBAS = NCTRIZ( 1 )
      NCHAU = NCTRIZ( 2 )

C     LE Z-CERCLE BAS COMPORTE LE PLUS GRAND NOMBRE D'ARETES DES 2
      IF( NBARLC( NCHAU ) .GT. NBARLC( NCBAS ) ) THEN
C        PERMUTATION DES 2 Z-CERCLES
         K     = NCHAU
         NCHAU = NCBAS
         NCBAS = K
      ENDIF

C     DIFFERENCE DES NOMBRES D'ARETES DES 2 Z-CERCLES
      NBARDI = NBARLC( NCBAS ) - NBARLC( NCHAU )

C     DIFFERENCE DE HAUTEUR EN Z DES 2 Z-CERCLES
      ZH = ABS( XYZCZC( 3, NCHAU ) - XYZCZC( 3, NCBAS ) )

C     LONGUEUR DE L'ARETE SUR LE CERCLE NCHAU
      R1 = 2 * Pi * RAYZC( NCHAU ) / NBARLC( NCHAU )

C     LONGUEUR DE L'ARETE SUR LE CERCLE NCBAS
      R2 = 2 * Pi * RAYZC( NCBAS ) / NBARLC( NCBAS )

C     LONGUEUR DE L'ARETE MOYENNE SUR LES CERCLES
      R = ( R1 + R2 ) / 2

C     LONGUEUR DE L'ARETE SUR LA HAUTEUR ZH
ccc      TAR = MIN( DARETE, (DARETE+R)/2, 3*R )
      TAR = MIN( DARETE, (DARETE+R)/2, 2*R )

C     NOMBRE D'ARETES DARETE DANS CETTE HAUTEUR ZH
      NBARZH = INT( ZH / TAR ) + 1

      IF( NBARDI .LE. 1 .AND. NBARZH .LE. 1 ) THEN
C        LES 2 Z-CERCLES ONT MEME NOMBRE D'ARETES ou +-1 et
C        LE NOMBRE D'ARETES EN Z EST <=1
C        PAS D'AJOUT DE Z-CERCLES INTERMEDIAIRES
         NBLGJD = 2
C        NUMERO LX DES 2 Z-CERCLES
         NULGZCP( 1 ) = NULGZC( NCBAS )
         NULGZCP( 2 ) = NULGZC( NCHAU )
         GOTO 30
      ENDIF

C     AJOUT DE NBLGJD-1 Z-CERCLES INTERMEDIAIRES ENTRE NCBAS ET NCHAU
C     NBLGJD nombre de lignes d''aretes a joindre
      NBLGJD = MAX( NBARDI, NBARZH ) + 1

C     NUMERO LX DES 2 Z-CERCLES EXTREMES
      NULGZCP( 1 )      = NULGZC( NCBAS )
      NULGZCP( NBLGJD ) = NULGZC( NCHAU )

C     NOMBRE D'ARETES DU Z-CERCLE NCBAS
      NZCARE = NBARLC( NCBAS )

      DO NC = 2, NBLGJD-1

C           CONSTRUCTION DU Z-CERCLE NC INTERMEDIAIRE
C           -----------------------------------------

C           CONSTRUCTION DU POINT CENTRE DU Z-CERCLE NC
C           ...........................................
            WRITE( KNC(1:8), '(I8)' ) NC
C           RETRAIT DES CARACTERES BLANCS
            CALL SANSBL( KNC, NBCAR )
            KNM = 'PCEZC_' // KNC(1:NBCAR) // '_AD '

C           CONSTRUCTION DU TMS LEXIQUE DU POINT CENTRE
C           SI CE POINT EXISTE, IL EST DETRUIT
            CALL LXLXOU( NTPOIN, KNM, NTCEZCP, MNCEZCP )
            IF( MNCEZCP .GT. 0 ) CALL LXTSDS( NTPOIN, KNM )
            CALL LXLXDC( NTPOIN, KNM, 24, 8 )
            CALL LXLXOU( NTPOIN, KNM, NTCEZCP, MNCEZCP )

C           NUMERO DU POINT CENTRE DU Z-CERCLE DANS LE LX DES POINTS
            CALL NUOBNM( 'POINT', KNM, NUPTCEZC )

C           CONSTRUCTION DU TMS XYZSOMMET DU POINT
            CALL LXTNDC( NTCEZCP, 'XYZSOMMET', 'ENTIER', WYZSOM+3 )
            CALL LXTSOU( NTCEZCP, 'XYZSOMMET',  NTXYZCEZC, MNXYZCEZC )
C           LE NOMBRE DE COORDONNEES PAR SOMMET
            MCN( MNXYZCEZC + WBCOOR ) = 3
C           LE NOMBRE DE SOMMETS
            MCN( MNXYZCEZC + WNBSOM) = 1
C           LE NOMBRE DE TANGENTES
            MCN( MNXYZCEZC + WNBTGS) = 0
C           LES 3 COORDONNEES DU CENTRE DU Z-CERCLE INTERMEDIAIRE PAR
C           INTERPOLATION LINEAIRE ENTRE LES 2 CENTRES DES Z-CERCLES
            RAP  = REAL( NC-1 ) / REAL( NBLGJD-1 )
            RAP1 = 1.0 - RAP
            DO K = 1, 3
               RMCN( MNXYZCEZC+WYZSOM-1+K ) = RAP1 * XYZCZC( K, NCBAS )
     %                                      + RAP  * XYZCZC( K, NCHAU )
            ENDDO
C           LA DATE DE CREATION
            CALL ECDATE( MCN(MNXYZCEZC) )
C           LE NUMERO DU TABLEAU DESCRIPTEUR
            MCN( MNXYZCEZC + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )

C           CONSTRUCTION DE LA LIGNE DU Z-CERCLE NC
C           .......................................
            WRITE( KNC(1:8), '(I8)' ) NC
C           RETRAIT DES CARACTERES BLANCS
            CALL SANSBL( KNC, NBCAR )
            KNMLGZCP = 'TRCLZC_' // KNC(1:NBCAR) // '_AD '

C           CONSTRUCTION DU TMS LEXIQUE de la LIGNE KNMLGZCP
C           SI CETTE LIGNE EXISTE, ELLE EST DETRUITE
            CALL LXLXOU( NTLIGN, KNMLGZCP, NTLGCP, MNLGCP )
            IF( MNLGCP .GT. 0 ) CALL LXTSDS( NTLIGN, KNMLGZCP )
            CALL LXLXDC( NTLIGN, KNMLGZCP, 24, 8 )
            CALL LXLXOU( NTLIGN, KNMLGZCP, NTLGCP, MNLGCP )

C           NUMERO DE LA LIGNE Z-CERCLE DANS LE LX DES LIGNES
            CALL NUOBNM( 'LIGNE', KNMLGZCP, NULGZCP(NC) )

C           CONSTRUCTION DU TMS a_ligne__definition du Z-CERCLE NC
C           LA LIGNE A UN TMS DEFINITION DE LIGNE DE TYPE 8:
C           CERCLE de R3 de TYPE NUTYCI=3 DEFINI PAR
C           LE CENTRE, RAYON, PLAN X ou Y ou Z = Cte
            CALL LXTNDC( NTLGCP, 'DEFINITION', 'MOTS', WUPLCT+1  )
            CALL LXTSOU( NTLGCP, 'DEFINITION', NT1CDE, MN1CDE )

C           TRANSFORMATION (I pour IDENTITE)
            MCN( MN1CDE + WTYTRL ) = 1
C           TYPE DE LA LIGNE: 8: CERCLE DE R3
            MCN( MN1CDE + WUTYLI ) = 8
C           NOMBRE D'ARETES DU CERCLE NC
            IF( NZCARE .GT. NBARLC(NCHAU) ) THEN
C              PERTE D'UNE ARETE
               NZCARE = NZCARE - 1
C           ELSE
C              NZCARE EST INCHANGE
            ENDIF
            MCN( MN1CDE + WBARLI ) = NZCARE
C           NUTYCI numero du type du cercle
            MCN( MN1CDE + WUTYCI ) = 3
C           NUPTCE numero du point centre du cercle nc
            MCN( MN1CDE + WUPTCE ) = NUPTCEZC
C           RAYDCI rayon du cercle nc
            RMCN( MN1CDE + WAYDCI ) = RAP1 * RAYZC( NCBAS )
     %                              + RAP  * RAYZC( NCHAU )
C           NUPLCT numero du plan XY a Z=Cte=CENTRE(3)
            MCN( MN1CDE + WUPLCT ) = 3
C           LA DATE DE CREATION
            CALL ECDATE( MCN(MN1CDE) )
C           LE NUMERO DU TABLEAU DESCRIPTEUR
            MCN( MN1CDE + MOTVAR(6) ) = NONMTD('~>LIGNE>>DEFINITION')

C           CONSTRUCTION DU MAILLAGE EN ARETES DU CERCLE NC
            CALL LIEX08( NTLGCP,  MCN(MN1CDE), RMCN(MN1CDE),
     %                   NTNSEFC, MNNSEFC, NTXYZC, MNXYZC, IERR )

            IF( IERR .NE. 0 ) THEN
               IERR = IERR + 10
               PRINT*,'arbresf1p: ERREUR dans le MAILLAGE du CERCLE',NC,
     %                 KNMLGZC
               IERR = 7
               GOTO 9900
            ENDIF

C           SUPPRESSION DES TANGENTES AU CERCLE
            MCN( MNXYZC  + WNBTGS ) = 0
            MCN( MNNSEFC + WBTGEF ) = 0
            MCN( MNNSEFC + WBEFAP ) = 0
            MCN( MNNSEFC + WBEFTG ) = 0

            IF( LANGAG .EQ. 0 ) THEN
               print*,'arbresf1p: ARETISATION du CERCLE intermediaire',
     %                 NC,' de NOM: ',KNMLGZCP,
     %                ' avec',MCN(MNNSEFC+WBEFOB),' ARETES et',
     %                 MCN(MNXYZC+WNBSOM),' SOMMETS est la LIGNE',
     %                 NULGZCP(NC)
            ELSE
               print*,'arbresf1p: EDGES of ADDED CIRCLE',NC,
     %                ' of NAME: ',KNMLGZCP,
     %                ' with',MCN(MNNSEFC+WBEFOB),' EDGES and',
     %                 MCN(MNXYZC+WNBSOM),' VERTICES is the LINE',
     %                 NULGZCP(NC)
            ENDIF

C           SUPPRESSION DU POINT CENTRE DU Z-CERCLE NC INTERMEDIAIRE
            IF( MNCEZCP .GT. 0 ) CALL LXTSDS( NTPOIN, KNM )

      ENDDO

 
C     LES CERCLES NCHAU et NCBAS ONT AU MOINS 3 SOMMETS
C     TRIANGULATION-QUADRANGULATION D'UNE SURFACE
C     PAR JONCTION DES EXTREMITES DES ARETES DE 2 LIGNES
C     --------------------------------------------------

C     CONSTRUCTION DU LEXIQUE de la SURFACE TRONC de CONE
C     NOM DE LA SURFACE DE CE TRONC DE CONE NOPIEC
 30   WRITE( KNC(1:8), '(I8)' ) NOPIEC
C     RETRAIT DES BLANCS
      CALL SANSBL( KNC, NBCAR )
      KNMPISF = 'TCS_' // KNC(1:NBCAR) // '_AD '

C     CONSTRUCTION DU TMS LEXIQUE de la SURFACE KNMPISF
C     SI CETTE SURFACE EXISTE, ELLE EST DETRUITE
      CALL LXLXOU( NTSURF, KNMPISF, NT1PISF, MN1PISF )
      IF( MN1PISF .GT. 0 ) CALL LXTSDS( NTSURF,KNMPISF )
      CALL LXLXDC( NTSURF, KNMPISF, 24, 8 )
      CALL LXLXOU( NTSURF, KNMPISF, NT1PISF, MN1PISF )

C     NUMERO DE LA SURFACE DANS LE LEXIQUE DES SURFACES
      CALL NUOBNM( 'SURFACE', KNMPISF, NUSFTRPI )

C     CONSTRUCTION DU TMS a_surface__definition du TRONC DE CONE
C     LA SURFACE A UN TMS DEFINITION DE SURFACE DE TYPE 33
      CALL LXTNDC( NT1PISF, 'DEFINITION', 'MOTS', WULGJN+NBLGJD )
      CALL LXTSOU( NT1PISF, 'DEFINITION', NT1PISFD, MN1PISFD )

C     TRANSFORMATION (I pour IDENTITE)
      MCN( MN1PISFD + WTYTRS ) = 1
C     TYPE DE LA SURFACE 33: TRONC DE CONE 
      MCN( MN1PISFD + WUTYSU ) = 33
C     NBLGJD nombre de lignes d''aretes a joindre
      MCN( MN1PISFD + WBLGJD ) = NBLGJD
C     NULGJD(1..NBLGJD) 'numero du nom des lignes a joindre
      DO NC = 1, NBLGJD
         MCN( MN1PISFD + WULGJD -1 + NC ) = NULGZCP( NC )
      ENDDO
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MN1PISFD) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN(MN1PISFD+MOTVAR(6))= NONMTD('~>SURFACE>>DEFINITION')

C     CONSTRUCTION DU MAILLAGE DE LA SURFACE TRONC DE CONE
      CALL SUEX33( NT1PISF, MCN(MN1PISFD),
     %             NTNSEFS, MNNSEFS, NTXYZS, MNXYZS, IERR )

      print*,'arbresf1p: TRIANGULATION de la SURFACE',NOPIEC,
     %       ' de NOM: ',KNMPISF,
     %       ' avec',MCN(MNNSEFS+WBEFOB),' TRIANGLES-QUADRANGLES',
     %         MCN(MNXYZS+WNBSOM),' SOMMETS'

C     DESTRUCTION DES EVENTUELLES LIGNES INTERMEDIAIRES CREEES
      DO NC = 2, NBLGJD-1
C        NUMERO LX DE LA LIGNE
         NULX = MCN( MN1PISFD + WULGJD -1 + NC )
         IF( NULX .GT. 0 ) THEN
C           NOM DE LA LIGNE CERCLE PROJETE DANS LE LX LIGNES
            CALL NMOBNU( 'LIGNE', NULX, KNM )
            CALL LXLXOU( NTLIGN, KNM, NT1LZC, MN1LZC )
            IF( MN1LZC .GT. 0 ) CALL LXTSDS( NTLIGN, KNM )
         ENDIF
      ENDDO

C     EVENTUELLE CREATION DES POINTS A IMPOSER DANS LA TETRAEDRISATION
C     LE LONG DE L'AXE JOIGNANT LE CENTRE DES 2 CERCLES BAS et HAUT
C     EN FONCTION de ARETGR
      IF( XYZCZC( 3,NCTRIZ(1) ) .GT. XYZCZC( 3,NCTRIZ(2) ) ) THEN
C        RACINES
         NCTRON = NCTRIZ( NBZCP )
         NDBRAN = 0
      ELSE
C        BRANCHES
         NCTRON = NCTRIZ( 1 )
         NDBRAN = 1
      ENDIF
      NBBRAN = NBZCP - 1

      GOTO 70


C     LA PIECE EST DE TYPE 4 ou 5
C     LA PIECE NOPIEC A PLUS DE 2 Z-CERCLES => RACINES ou BRANCHES
C     ============================================================
C     RECHERCHE DU Z-CERCLE DE COTE DIFFERENTE DES AUTRES Z-CERCLES
C     C'EST LA PREMIERE OU LA DERNIERE SUITE AU TRI CROISSANT DES Z
 40   IF( XYZCZC( 3,NCTRIZ(1) ) .EQ. XYZCZC( 3,NCTRIZ(2) ) ) THEN
C        RACINES
         NCTRON = NCTRIZ( NBZCP )
         NDBRAN = 0
         NOTYPE = 5
      ELSE
C        BRANCHES
         NCTRON = NCTRIZ( 1 )
         NDBRAN = 1
         NOTYPE = 4
      ENDIF
      NBBRAN = NBZCP - 1

C     LE CERCLE NCTRON A NBBRAN Z-CERCLES BRANCHES
C     ============================================
C     LE CERCLE NCTRON EST IL REDUIT A SON CENTRE?
      IF( NBARLC( NCTRON ) .LE. 0 ) THEN
C        ERREUR: JONCTION ENTRE 1 CENTRE ET PLUSIEURS CERCLES
         PRINT *,'arbresf1p: PIECE',NOPIEC,' le Z-CERCLE',NCTRON,
     %           ' REDUIT a son CENTRE EST a JOINDRE a',
     %            NBBRAN,' CERCLES BRANCHES'
         PRINT *,'SITUATION INCORRECTE'
         IERR = 5
         GOTO 9900
      ENDIF

C     RAYON INITIAL DU CERCLE NCTRON
      RAYZC0 = RAYZC( NCTRON )

C     AMPLIFICATION DU RAYON DU CERCLE NCTRON
      AMPLIRC(0) = 1.

C     SOMME DES RAYONS DES CERCLES DE LA PIECE
      SRAY = RAYZC( NCTRON )
      DO NOBRA = 1, NBBRAN

C        NZCBRA DE 1 A NBZCER EST UN Z-CERCLE BRANCHE OPPOSE a NCTRON
         NZCBRA = NCTRIZ( NDBRAN + NOBRA )

C        SOMME DES RAYONS DES CERCLES
         SRAY = SRAY + RAYZC( NZCBRA )

C        AMPLIFICATION DU RAYON DU Z-CERCLE BRANCHE PROJETE NOBRA
         AMPLIRC( NOBRA ) = 1.

      ENDDO

C     PROJECTION-CONTRACTION DES NBBRAN CERCLES BRANCHES DANS LE
C     CERCLE NCTRON SANS INTERSECTION DE CES CERCLES PROJETES
C     TRIANGULATION suex09 du CERCLE NCTRON MOINS LES CERCLES BRAN PROJETES
C     PROJECTION DES SOMMETS DES CERCLES PROJETES SUR LES CERCLES BRANCHES
C     PROJECTION SUR LE TRONC DE CONE DES SOMMETS INTERNES A LA
C     TRIANGULATION et NON MODIFICATION DES SOMMETS SUR LE CERCLE NCTRON
C     ---------------------------------------------------------------------
C     CONSTRUCTION DU SOMMET SCE INFERIEUR DU CONE ENGLOBANT
C     POUR LA PROJECTION A PARTIR DES XYZ DU CENTRE DU CERCLE NCTRON
C     ZH HAUTEUR EN Z ENTRE LE PLAN NCTRON ET LE SOMMET DU CONE ENGLOBANT
ccc      ZH = RAYZC( NCTRON ) / 2
      ZH = SRAY / 3
      IF( NOTYPE .EQ. 4 ) THEN
         ZHSCE = -ZH
      ELSE
         ZHSCE =  ZH
      ENDIF

C     XYZ DU POINT FOCAL SCE
C     ----------------------
      XYZSCE( 1 ) = XYZCZC( 1, NCTRON )
      XYZSCE( 2 ) = XYZCZC( 2, NCTRON )
      XYZSCE( 3 ) = XYZCZC( 3, NCTRON ) + ZHSCE

C     LES 3 POINTS DEFINISSANT LE PLAN DU CERCLE NCTRON
C     -------------------------------------------------
      P1PLAN( 1 ) = XYZCZC( 1, NCTRON )
      P1PLAN( 2 ) = XYZCZC( 2, NCTRON )
      P1PLAN( 3 ) = XYZCZC( 3, NCTRON )

      P2PLAN( 1 ) =  P1PLAN( 1 ) + RAYZC0
      P2PLAN( 2 ) =  P1PLAN( 2 )
      P2PLAN( 3 ) =  P1PLAN( 3 )

      P3PLAN( 1 ) =  P1PLAN( 1 )
      P3PLAN( 2 ) =  P1PLAN( 2 ) + RAYZC0
      P3PLAN( 3 ) =  P1PLAN( 3 )

C     TABLEAU DE PROTECTION DU TMS XYZSOMMET DU CERCLE NCTRON
      CALL TNMCDC( 'ENTIER', NBBRAN, MNXYZCP )

C     XYZCZCP DU POINT INTERSECTION DU PLAN DU CERCLE NCTRON AVEC
C     LA DROITE CENTRE DU CERCLE BRANCHE-POINT FOCAL XYZSCE
C     -----------------------------------------------------------
      DO NOBRA = 1, NBBRAN

C           NZCBRA EST UN BRAN DE NCTRON
            NZCBRA = NCTRIZ( NDBRAN + NOBRA )

C           PROJECTION DU CENTRE DU CERCLE NZCBRA SUR LE PLAN DU CERCLE NCTRON
            XYZ1( 1 ) = XYZCZC( 1, NZCBRA )
            XYZ1( 2 ) = XYZCZC( 2, NZCBRA )
            XYZ1( 3 ) = XYZCZC( 3, NZCBRA )

C           POINT INTERSECTION DU PLAN Z DU CERCLE NCTRON AVEC
C           LA DROITE SOMMET INFERIEUR XYZSCE-XYZ1
            CALL INDRPL( XYZSCE, XYZ1, P1PLAN, P2PLAN, P3PLAN,
     %                   XYZCZCP(1,NOBRA), NOCODE )
            IF( NOCODE .NE. 0 ) THEN
               PRINT*,'arbresf1p: PROJECTION du CENTRE du CERCLE,',
     %                 NZCBRA,' NON CALCULABLE'
               IERR = 9
               GOTO 9900
            ENDIF

C           RAYON DU CERCLE PROJETE NZCBRA SUR NCTRON
C           NUMERO DE LA LIGNE DANS LX LIGNES
            NULXLF = NULGZC( NZCBRA )

            IF( NULXLF .GT. 0 ) THEN
C              DISTANCE CENTRE DU CERCLE NZCBRA Pt XYZSCE
               D = ( XYZCZC(1,NZCBRA) - XYZSCE(1) ) **2
     %           + ( XYZCZC(2,NZCBRA) - XYZSCE(2) ) **2
     %           + ( XYZCZC(3,NZCBRA) - XYZSCE(3) ) **2
               D = SQRT( D )

C              DISTANCE CENTRE DU CERCLE PROJETE NOBRA Pt SCE
               D4 = ( XYZCZCP(1,NOBRA) - XYZSCE(1) ) **2
     %            + ( XYZCZCP(2,NOBRA) - XYZSCE(2) ) **2
     %            + ( XYZCZCP(3,NOBRA) - XYZSCE(3) ) **2
               D4 = SQRT( D4 )

               RAYZCP( NOBRA ) = RAYZC( NZCBRA ) * D4 / D
            ELSE

C              CERCLE PROJETE REDUIT AU CENTRE PROJETE
               RAYZCP( NOBRA ) = 0D0

            ENDIF

      ENDDO

C     VERIFICATION: LE CENTRE DU CERCLE PROJETE EST IL
C                   DANS LE CERCLE NCTRON?
C     ------------------------------------------------
 22   DO NOBRA = 1, NBBRAN

C           NZCBRA NUMERO DU Z-CERCLE NON NCTRON
            NZCBRA = NCTRIZ( NDBRAN + NOBRA )

            D = ( XYZCZCP(1,NOBRA) - XYZCZC(1,NCTRON) ) **2
     %        + ( XYZCZCP(2,NOBRA) - XYZCZC(2,NCTRON) ) **2
     %        + ( XYZCZCP(3,NOBRA) - XYZCZC(3,NCTRON) ) **2
            D = SQRT( D )

            IF( D .GE. RAYZC(NCTRON) + RAYZCP(NOBRA) ) THEN

C              LE CERCLE PROJETE NOBRA EST EXTERNE AU CERCLE NCTRON
               PRINT*
               PRINT*,'arbresf1p: PIECE',NOPIEC,': CERCLE',NCTRON,
     %                ' le CERCLE PROJETE',NOBRA,
     %                ' lui est EXTERNE D=',D,
     %                ' Rayons=',RAYZC(NCTRON),RAYZCP(NOBRA),
     %       ' => Dilatation du CERCLE ou Contraction du CERCLE PROJETE'
               PRINT*,'arbresf1p: CERCLE NON PROJETE',NCTRON,' CENTRE',
     %                (XYZCZC(K,NCTRON),K=1,3)
               PRINT*,'arbresf1p: CERCLE     PROJETE',NOBRA,' CENTRE',
     %                (XYZCZCP(K,NOBRA),K=1,3)

C              AMPLIFICATION DU RAYON DU CERCLE NCTRON
               R = ( D*1.25D0 + 2 * RAYZCP(NOBRA) ) / RAYZC(NCTRON)
               AMPLIRC(0) = REAL( AMPLIRC(0) * R )

C              RAYON AMPLIFIE DU CERCLE NCTRON
               RAYZC( NCTRON ) = REAL( RAYZC( NCTRON ) * R )

            ELSE IF( D .GT. RAYZC(NCTRON)-RAYZCP(NOBRA) ) THEN

C              LE CERCLE PROJETE NOBRA N'EST PAS INTERNE AU CERCLE NCTRON
C              AMPLIFICATION DU RAYON DU CERCLE NCTRON
C              OU CONTRACTION DU CERCLE PROJETE NOBRA
               PRINT*,'arbresf1p: PIECE',NOPIEC,': CERCLE',NCTRON,
     %                ' le CERCLE PROJETE',NOBRA,
     %                ' l''INTERSECTE D=',D,
     %                ' Rayons=',RAYZC(NCTRON),RAYZCP(NOBRA),
     %       ' => DILATATION du CERCLE ou CONTRACTION du CERCLE PROJETE'
                PRINT*,'arbresf1p: CERCLE NON PROJETE',NCTRON,' CENTRE',
     %                          (XYZCZC(K,NCTRON),K=1,3)
                 PRINT*,'arbresf1p: CERCLE     PROJETE',NOBRA,' CENTRE',
     %                          (XYZCZCP(K,NOBRA),K=1,3)

C              AMPLIFICATION DU RAYON DU CERCLE NCTRON
               R = ( D*1.25D0 + 2 * RAYZCP(NOBRA) ) / RAYZC(NCTRON)
               AMPLIRC(0) = REAL( AMPLIRC(0) * R )

C              RAYON AMPLIFIE DU CERCLE NCTRON
               RAYZC( NCTRON ) = REAL( RAYZC( NCTRON ) * R )

            ENDIF

      ENDDO


C     VERIFICATION DE NON INTERSECTION DES CERCLES PROJETES ENTRE EUX
C     ---------------------------------------------------------------
      DO NOBRA1 = 1, NBBRAN-1

            DO 26 NOBRA2 = NOBRA1+1, NBBRAN

               IF( RAYZCP(NOBRA1) .LE. 0  .AND.
     %             RAYZCP(NOBRA2) .LE. 0 ) GOTO 26

C              AU MOINS UN DES RAYONS DES 2 CERCLES EST NON NUL
               D = ( XYZCZCP(1,NOBRA2) - XYZCZCP(1,NOBRA1) ) **2
     %           + ( XYZCZCP(2,NOBRA2) - XYZCZCP(2,NOBRA1) ) **2
     %           + ( XYZCZCP(3,NOBRA2) - XYZCZCP(3,NOBRA1) ) **2
               D = SQRT( D )

               IF( D .LT. RAYZCP(NOBRA1) + RAYZCP(NOBRA2) ) THEN

C                 LES 2 CERCLES S'INTERSECTENT
                  IF( RAYZCP(NOBRA1) .GT. 0 .AND.
     %                RAYZCP(NOBRA2) .GT. 0 ) THEN

C                    CONTRACTION DU RAYON NON NUL DES 2 CERCLES
                     PRINT*
                  PRINT*,'arbresf1p: PIECE',NOPIEC,': CERCLES PROJETES',
     %                       NOBRA1,NOBRA2,
     %                      ' en INTERSECTION D=',D,
     %                      ' Rayons=',RAYZCP(NOBRA1),RAYZCP(NOBRA2),
     %                      ' => CONTRACTION des 2 RAYONS'
                     PRINT*,'arbresf1p: CERCLE PROJETE',NOBRA1,' CENTRE'
     %                     ,(XYZCZCP(K,NOBRA1),K=1,3)
                     PRINT*,'arbresf1p: CERCLE PROJETE',NOBRA2,' CENTRE'
     %                     ,(XYZCZCP(K,NOBRA2),K=1,3)

C                    CONTRACTION DU RAYON DU CERCLE PROJETE NOBRA1
ccc                     D = 0.9D0 * D
                     D = 0.8D0 * D
                     IF( RAYZCP(NOBRA1) .GT. 0 ) THEN
                        R = ( D-RAYZCP(NOBRA2) ) / RAYZCP(NOBRA1)
                        AMPLIRC(NOBRA1) = REAL( AMPLIRC(NOBRA1) * R )
                        RAYZCP( NOBRA1 ) = RAYZCP( NOBRA1 ) * R
                     ENDIF

C                    CONTRACTION DU RAYON DU CERCLE PROJETE NOBRA2
                     IF( RAYZCP(NOBRA2) .GT. 0 ) THEN
                        R = ( D-RAYZCP(NOBRA1) ) / RAYZCP(NOBRA2)
                        AMPLIRC(NOBRA2) = REAL( AMPLIRC(NOBRA2) * R )
                        RAYZCP( NOBRA2 ) = RAYZCP( NOBRA2 ) * R
                     ENDIF

                     PRINT*,'arbresf1p: CONTRACTION DU RAYON DES CERCLES
     % PROJETES',NOBRA1,NOBRA2
                     GOTO 22

                  ELSE

C                    UN DES 2 RAYONS EST NUL
C                    IL FAUT SORTIR LE CENTRE DU CERCLE PROJETE DE RAYON NUL
C                    DE L'AUTRE CERCLE DE RAYON NON NUL
                     IF( RAYZCP(NOBRA1) .LE. 0 ) THEN
                        NC0 = NOBRA1
                        NC1 = NOBRA2
                     ELSE
                        NC0 = NOBRA2
                        NC1 = NOBRA1
                     ENDIF

C                    CONTRACTION DU RAYON DU CERCLE NC1
C                    POUR EN SORTIR LE CENTRE DU CERCLE NC1
ccc                     R = 0.8D0 * D / RAYZCP(NC1)
                     R = 0.75D0 * D / RAYZCP(NC1)
                     AMPLIRC(NC1) = REAL( AMPLIRC(NC1) * R )
                     RAYZCP(NC1) = 0.75D0 * D
ccc                     RAYZCP(NC1) = 0.8D0 * D
                     PRINT*,'arbresf1p: SORTIE du CENTRE PROJETE',NC0,
     %                      ' du CERCLE PROJETE',NC1,
     %                      ' PAR CONTRACTION',R,' de son RAYON'
                     GOTO 22

                  ENDIF

               ENDIF

 26         ENDDO

      ENDDO

C     ICI LES CERCLES PROJETES SONT INTERNES AU CERCLE NCTRON
C     et NE S'INTERSECTENT PAS
C     -------------------------------------------------------
      PRINT*
      PRINT*,'arbresf1p: PIECE',NOPIEC
      PRINT*,'arbresf1p: CERCLE',NCTRON,' ENGLOBANT   de    CENTRE',
     %          (XYZCZC(K,NCTRON),K=1,3),' RAYON        =',RAYZC(NCTRON)

      DO NOBRA = 1, NBBRAN
C        NZCBRA NUMERO DU Z-CERCLE NON PROJETE NON NCTRON
         NZCBRA = NCTRIZ( NDBRAN + NOBRA )
         PRINT*,'arbresf1p: CERCLE',NZCBRA,' CENTRE DU CERCLE PROJETE',
     %         (XYZCZCP(K,NOBRA),K=1,3),' RAYON PROJETE=',RAYZCP(NOBRA)
      ENDDO

C     CONSTRUCTION DES SOMMETS DES NBBRAN CERCLES PROJETES
C     SUR LE CERCLE NCTRON C-A-D INTERSECTION DU PLAN P1PLAN-P2PLAN-P3PLAN
C     DU CERCLE NCTRON ET DE LA DROITE POINT DU CERCLE NOBRA-POINT FOCAL
C     --------------------------------------------------------------------
C     NOMBRE DE Z-CERCLES REDUIT A UN SOMMET D'UN CONE
      NBZCSTCO = 0
      NBLBRA   = 0
      DO 32 NOBRA = 1, NBBRAN

C           NZCBRA DE 1 A NBZCER EST UN Z-CERCLE BRANCHE DE NCTRON
            NZCBRA = NCTRIZ( NDBRAN + NOBRA )

C           NUMERO DE LA LIGNE DANS LX LIGNES
            NULXLF = NULGZC( NZCBRA )

            IF( NULXLF .LE. 0 ) THEN
C              LIGNE REDUITE AU SOMMET d'un CONE
               NBZCSTCO = NBZCSTCO + 1
               GOTO 32
            ENDIF

C           OUVERTURE DU LX DE CETTE LIGNE BRAN NON SOMMET d'un CONE
            NBLBRA = NBLBRA + 1
            CALL LXNLOU( NTLIGN, NULXLF, NT1LZC, MN1LZC )

C           RESTAURATION DES TABLEAUX XYZSOMMET ET NSEF
            CALL LXTSOU( NT1LZC, 'XYZSOMMET', NTSOLI, MNSOLI )
            IF( NTSOLI .LE. 0 ) THEN
               PRINT*,'arbresf1p: LIGNE',NULXLF,' SANS TMS XYZSOMMET'
               IERR = 6
               GOTO 9900
            ENDIF
            CALL LXTSOU( NT1LZC, 'NSEF', NTARLI, MNARLI )
            IF( NTARLI .LE. 0 ) THEN
               PRINT*,'arbresf1p: LIGNE',NULXLF,' SANS TMS NSEF'
               IERR = 7
               GOTO 9900
            ENDIF

C           PROJECTION DES SOMMETS DES ARETES DE LA LIGNE CERCLE BRAN
C           LE NOMBRE DE SOMMETS DE LA LIGNE BRAN
            NBSOLI = MCN( MNSOLI + WNBSOM )
C           LE NOMBRE D'ARETES DE LA LIGNE BRAN
            NBARLI = MCN( MNARLI + WBEFOB )
               
            MNXYZLP = 0
            CALL TNMCDC( 'REEL', 3*NBSOLI, MNXYZLP )

            MNXS1 = MNSOLI + WYZSOM -1
            MNXS2 = MNXYZLP -1
            DO I = 1, NBSOLI

C              PROJECTION DU SOMMET I DU CERCLE SUR LE PLAN DU CERCLE NCTRON
               DO K=1,3
                  XYZ1( K ) = RMCN( MNXS1 + K )
               ENDDO

C              POINT INTERSECTION DU PLAN Z DU CERCLE NCTRON AVEC
C              LA DROITE SOMMET INFERIEUR XYZSCE-XYZ1
               CALL INDRPL( XYZSCE, XYZ1, P1PLAN, P2PLAN, P3PLAN,
     %                      XYZ2, NOCODE )
               IF( NOCODE .NE. 0 ) THEN
                  PRINT*,'arbresf1p: PROJECTION du POINT du CERCLE,',
     %                    NZCBRA,' NON CALCULABLE'
                  IERR = 8
                  GOTO 9900
               ENDIF

C              PRISE EN COMPTE DE L'AMPLIFICATION(ou CONTRACTION)
C              DU CERCLE PROJETE
               DO K=1,3
                  RMCN( MNXS2 + K ) = REAL( XYZCZCP(K,NOBRA)
     %                 + ( XYZ2(K)-XYZCZCP(K,NOBRA) ) * AMPLIRC(NOBRA) )
               ENDDO

               MNXS1 = MNXS1 + 3
               MNXS2 = MNXS2 + 3

            ENDDO

C           KNMLGZC NOM DE LA LIGNE MOMENTANEE DU CERCLE NZCBRA PROJETE
            WRITE( KNC(1:8), '(I8)' ) NZCBRA
C           RETRAIT DES CARACTERES BLANCS
            CALL SANSBL( KNC, NBCAR )
            KNMLGZCP = 'LCPR_' // KNC(1:NBCAR) // '_AD '

C           CONSTRUCTION DU TMS LEXIQUE de la LIGNE KNMLGZCP
C           SI CETTE LIGNE EXISTE, ELLE EST DETRUITE
            CALL LXLXOU( NTLIGN, KNMLGZCP, NTLGCP, MNLGCP )
            IF( MNLGCP .GT. 0 ) CALL LXTSDS( NTLIGN, KNMLGZCP )
            CALL LXLXDC( NTLIGN, KNMLGZCP, 24, 8 )
            CALL LXLXOU( NTLIGN, KNMLGZCP, NTLGCP, MNLGCP )

C           NUMERO DE LA LIGNE CERCLE PROJETE DANS LE LX DES LIGNES
            CALL NUOBNM( 'LIGNE', KNMLGZCP, NULGZCP(NOBRA) )

C           CONSTRUCTION DU TMS a_ligne__definition du CERCLE PROJETE
C           LA LIGNE EST DEFINIE PAR SES TABLEAUX XYZSOMMET ET NSEF OPTION 10
            CALL LXTNDC( NTLGCP, 'DEFINITION', 'ENTIER', WUTSSL+1 )
            CALL LXTSOU( NTLGCP, 'DEFINITION',  NTDFLP , MNDFLP )
C           LA TRANSFORMATION
            MCN( MNDFLP + WTYTRL ) = 1
C           LE TYPE DE LA LIGNE
            MCN( MNDFLP + WUTYLI ) = 10

C           CONSTRUCTION DU TMS XYZSOMMET DE LA LIGNE PROJETEE
            CALL LXTNDC(NTLGCP,'XYZSOMMET', 'ENTIER',WYZSOM+3*NBSOLI)
            CALL LXTSOU(NTLGCP,'XYZSOMMET',  NTSOLP, MNSOLP )
C           SAUVEGARDE DE L'ADRESSE DU TMS XYZSOMMET
            MCN( MNXYZCP-1+NOBRA) = MNSOLP
C           LE NOMBRE DE COORDONNEES PAR SOMMET
            MCN( MNSOLP + WBCOOR ) = 3
C           LE NOMBRE DE SOMMETS
            MCN( MNSOLP + WNBSOM) = NBSOLI
C           LE NOMBRE DE TANGENTES
            MCN( MNSOLP + WNBTGS) = 0
C           COPIE DES XYZ DES POINTS PROJETES
            CALL TRTATA(RMCN(MNXYZLP), RMCN(MNSOLP+WYZSOM), 3*NBSOLI)
C           LA DATE DE CREATION
            CALL ECDATE( MCN(MNSOLP) )
C           LE NUMERO DU TABLEAU DESCRIPTEUR
            MCN( MNSOLP + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )

C           CONSTRUCTION DU TMS NSEF DE LA LIGNE PROJETEE
            CALL LXTNDC( NTLGCP, 'NSEF', 'ENTIER', WUSOEF+2*NBARLI )
            CALL LXTSOU( NTLGCP, 'NSEF',  NTARLP, MNARLP )
C           TYPE DE L'OBJET : LIGNE CERCLE = LIGNE CERCLE PROJETEE
            CALL TRTATA( RMCN(MNARLI), RMCN(MNARLP), WUSOEF+2*NBARLI )
C           SUPPRESSION DES TANGENTES POUR LA LIGNE PROJETEE
            MCN( MNARLP + WBTGEF ) = 0
            MCN( MNARLP + WBEFAP ) = 0
            MCN( MNARLP + WBEFTG ) = 0
C           LA DATE DE CREATION DU TMS NSEF DE LA LIGNE PROJETEE
            CALL ECDATE( MCN(MNARLP) )
C           LE NUMERO DU TABLEAU DESCRIPTEUR
            MCN( MNARLP + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )

C           DEFINITION EST COMPLETE
C           TABLEAU XYZSOMMET DES SOMMETS DU TMS DEFINITION
            MCN( MNDFLP + WUTSOL ) = NTSOLP
C           TABLEAU NSEF NO DES SOMMETS DES ARETES
            MCN( MNDFLP + WUTSSL ) = NTARLP
C           LA DATE DE CREATION DU TMS DEFINITION DU CERCLE PROJETE
C           A FAIRE APRES LA CONSTRUCTION DES TMS NSEF XYZSOMMET
C           A CAUSE DE LA MISE A JOUR AUTOMATIQUE
            CALL ECDATE( MCN(MNDFLP) )
C           LE NUMERO DU TABLEAU DESCRIPTEUR
            MCN( MNDFLP + MOTVAR(6) ) = NONMTD('~>LIGNE>>DEFINITION')

            print*,'arbresf1p: ARETISATION du Z-CERCLE PROJETE',NZCBRA,
     %             ' de NOM: ',KNMLGZCP,' avec',MCN(MNARLP+WBEFOB),
     %             ' ARETES', MCN(MNSOLP+WNBSOM),
     %             ' SOMMETS est la LIGNE',NULGZCP(NOBRA)

            CALL TNMCDS( 'REEL', 3*NBSOLI, MNXYZLP )

 32   ENDDO


C     AMPLIFICATION DU RAYON DU Z-CERCLE NCTRON ET DE SES SOMMETS
C     -----------------------------------------------------------
      print*,'arbresf1p: AMPLIFICATION=',AMPLIRC(0),
     %       ' du RAYON',NCTRON,' =',RAYZC0,' =>',RAYZC(NCTRON)

C     OUVERTURE DU LX DE CE Z-CERCLE NCTRON
      NULXLF = NULGZC( NCTRON )
      CALL LXNLOU( NTLIGN, NULXLF, NT1LZC, MN1LZC )

C     OUVERTURE DU TABLEAU XYZSOMMET DU Z-CERCLE NCTRON
      CALL LXTSOU( NT1LZC, 'XYZSOMMET', NTSOLI, MNSOLI )
      IF( NTSOLI .LE. 0 ) THEN
         PRINT*,'arbresf1p: LIGNE',NULXLF,' SANS TMS XYZSOMMET'
         IERR = 6
         GOTO 9900
      ENDIF

C     LE NOMBRE DE SOMMETS DE LA LIGNE
      NBSOLI = MCN( MNSOLI + WNBSOM )

C     SAUVEGARDE DES XYZ DES NBSOLI SOMMETS DU TMS XYZSOMMET
      MNXYZLP = 0
      CALL TNMCDC( 'REEL', 3*NBSOLI, MNXYZLP )
      MNXS1 = MNSOLI + WYZSOM
      CALL TRTATA( RMCN(MNXS1), RMCN(MNXYZLP), 3*NBSOLI )

C     AMPLIFICATION DES XYZ DES SOMMETS DE LA LIGNE CERCLE NCTRON
      MNXS1 = MNSOLI + WYZSOM - 1
      DO I = 1, NBSOLI
         DO K=1,3
            RMCN(MNXS1+K) = XYZCZC(K,NCTRON)
     %                  +( RMCN(MNXS1+K)-XYZCZC(K,NCTRON) ) * AMPLIRC(0)
         ENDDO
         MNXS1 = MNXS1 + 3
      ENDDO


C     CONSTRUCTION DES NBZCSTCO  POINTS SOMMET PROJETE DES CONES
C     POUR DEVENIR DES SOMMETS IMPOSES DANS LA TRIANGULATION
C     DU CERCLE NCTRON
C     --------------------------------------------------------
      IF( NBZCSTCO .GT. 0 ) THEN

            NBZCSTCO = 0
            DO NOBRA = 1, NBBRAN

C              NZCBRA DE 1 A NBZCER EST UN Z-CERCLE BRANCHE DE NCTRON
               NZCBRA = NCTRIZ( NDBRAN + NOBRA )

C              NUMERO DE LA LIGNE DANS LX LIGNES
               NULXLF = NULGZC( NZCBRA )

               IF( NULXLF .LE. 0 ) THEN

C                 LIGNE REDUITE AU POINT-SOMMET d'un CONE
C                 CONSTRUCTION DU POINT CENTRE PROJETE SUR NCTRON
                  NBZCSTCO = NBZCSTCO + 1
                  WRITE( KNC(1:8), '(I8)' ) NBZCSTCO
C                 RETRAIT DES CARACTERES BLANCS
                  CALL SANSBL( KNC, NBCAR )
                  KNM = 'PCPR_' // KNC(1:NBCAR) // '_AD '

C                 CONSTRUCTION DU TMS LEXIQUE DU POINT CENTRE
C                 SI CE POINT EXISTE, IL EST DETRUIT
                  CALL LXLXOU( NTPOIN, KNM, NTCEZCP, MNCEZCP )
                  IF( MNCEZCP .GT. 0 ) CALL LXTSDS( NTPOIN, KNM )
                  CALL LXLXDC( NTPOIN, KNM, 24, 8 )
                  CALL LXLXOU( NTPOIN, KNM, NTCEZCP, MNCEZCP )

C                 NUMERO DU POINT CENTRE DU Z-CERCLE PROJETE
C                 DANS LE LX DES POINTS
                  CALL NUOBNM( 'POINT', KNM, NPt )
                  NULGZCP(NOBRA) = -NPt

C                 CONSTRUCTION DU TMS XYZSOMMET DU POINT
                  CALL LXTNDC(NTCEZCP, 'XYZSOMMET', 'ENTIER', WYZSOM+3 )
                  CALL LXTSOU(NTCEZCP, 'XYZSOMMET', NTXYZCEZC,MNXYZCEZC)
C                 LE NOMBRE DE COORDONNEES PAR SOMMET
                  MCN( MNXYZCEZC + WBCOOR ) = 3
C                 LE NOMBRE DE SOMMETS
                  MCN( MNXYZCEZC + WNBSOM) = 1
C                 LE NOMBRE DE TANGENTES
                  MCN( MNXYZCEZC + WNBTGS) = 0

C                 LES 3 COORDONNEES DU CENTRE DU Z-CERCLE PROJETE
                  DO K = 1, 3
                     RMCN(MNXYZCEZC+WYZSOM-1+K) = REAL(XYZCZCP(K,NOBRA))
                  ENDDO

C                 TAILLE IDEALE DES ARETES SOUHAITEES ISSUES DE CE POINT
                  TAIDAR(NBZCSTCO) = MIN( DARETE/8, RAYZC(NCTRON)/12 )

                  PRINT*,'arbresf1p: CERCLE PROJETE=CENTRE=POINT=',KNM,
     %                (XYZCZCP(K,NOBRA),K=1,3),' RAYON=',RAYZCP(NOBRA),
     %                ' ARETE AUTOUR=',TAIDAR(NBZCSTCO)
 
C                 LA DATE DE CREATION
                  CALL ECDATE( MCN(MNXYZCEZC) )
C                 LE NUMERO DU TABLEAU DESCRIPTEUR
                  MCN( MNXYZCEZC + MOTVAR(6) ) = NONMTD('~>>>XYZSOMMET')

               ENDIF

            ENDDO
      ENDIF


C     CONSTRUCTION DES POINTS A MI-DISTANCE ENTRE 2 CERCLES PROJETES
C     IMPOSES ENSUITE COMME POINTS AVEC UNE DISTANCE SOUHAITEE PETITE
C     POUR OBTENIR PLUS DE TRIANGLES ENTRE LES CERCLES 2 A 2
C     LE TOUT DANS LE  CERCLE NCTRON
C     ---------------------------------------------------------------
      NBZCSTCO0 = NBZCSTCO
      DO NOBRA1 = 1, NBBRAN-1
            DO 34 NOBRA2 = NOBRA1+1, NBBRAN

C              DISTANCE ENTRE LES 2 CENTRES PROJETES
               D = ( XYZCZCP(1,NOBRA2) - XYZCZCP(1,NOBRA1) ) **2
     %           + ( XYZCZCP(2,NOBRA2) - XYZCZCP(2,NOBRA1) ) **2
     %           + ( XYZCZCP(3,NOBRA2) - XYZCZCP(3,NOBRA1) ) **2
               D = SQRT( D )

C              DISTANCE ENTRE LES 2 CERCLES
               R = D - ( RAYZCP(NOBRA1) + RAYZCP(NOBRA2) )
               IF( R .GT. DARETE*1.25D0 ) GOTO 34

C              CONSTRUCTION DU POINT A MI-DISTANCE ENTRE LES 2 Z-CERCLES
               NBZCSTCO = NBZCSTCO + 1

               WRITE( KNC(1:8), '(I8)' ) NBZCSTCO
C              RETRAIT DES CARACTERES BLANCS
               CALL SANSBL( KNC, NBCAR )
               KNM = 'PMID_' // KNC(1:NBCAR) // '_AD '

C              CONSTRUCTION DU TMS LEXIQUE DU POINT CENTRE
C              SI CE POINT EXISTE, IL EST DETRUIT
               CALL LXLXOU( NTPOIN, KNM, NTCEZCP, MNCEZCP )
               IF( MNCEZCP .GT. 0 ) CALL LXTSDS( NTPOIN, KNM )
               CALL LXLXDC( NTPOIN, KNM, 24, 8 )
               CALL LXLXOU( NTPOIN, KNM, NTCEZCP, MNCEZCP )

C              NUMERO DU POINT MI-DISTANCE DANS LE LX DES POINTS
               CALL NUOBNM( 'POINT', KNM, NPt )
               NULGZCP( NBBRAN + NBZCSTCO ) = -NPt

C              CONSTRUCTION DU TMS XYZSOMMET DU POINT
               CALL LXTNDC(NTCEZCP, 'XYZSOMMET', 'ENTIER', WYZSOM+3 )
               CALL LXTSOU(NTCEZCP, 'XYZSOMMET', NTXYZCEZC,MNXYZCEZC)
C              LE NOMBRE DE COORDONNEES PAR SOMMET
               MCN( MNXYZCEZC + WBCOOR ) = 3
C              LE NOMBRE DE SOMMETS
               MCN( MNXYZCEZC + WNBSOM) = 1
C              LE NOMBRE DE TANGENTES
               MCN( MNXYZCEZC + WNBTGS) = 0

C              LE POINT SUR LA DROITE DES CENTRES DES 2 CERCLES PROJETES
C              ET SUR LE CERCLE NOBRA1 ET NOBRA2
               DO K=1,3
                  XYZ1(K) =  XYZCZCP(K,NOBRA1) + RAYZCP(NOBRA1) / D
     %                    * (XYZCZCP(K,NOBRA2) - XYZCZCP(K,NOBRA1))
                  XYZ2(K) =  XYZCZCP(K,NOBRA2) + RAYZCP(NOBRA2) / D
     %                    * (XYZCZCP(K,NOBRA1) - XYZCZCP(K,NOBRA2))
               ENDDO
C              LES 3 COORDONNEES DU POINT A MI-DISTANCE DES 2 CERCLES
               DO K = 1, 3
                  RMCN(MNXYZCEZC+WYZSOM-1+K)=REAL((XYZ1(K)+XYZ2(K))/2)
               ENDDO

C              TAILLE IDEALE DES ARETES SOUHAITEES ISSUES DU POINT
C              A MI-DISTANCE DES 2 CERCLES
               TAIDAR(NBZCSTCO) = MIN( R/8, DARETE/8, RAYZC(NCTRON)/12 )

               PRINT*,'arbresf1p: POINT MI-DISTANCE 2 CERCLES=',KNM,
     %                (RMCN(MNXYZCEZC+WYZSOM-1+K),K=1,3),
     %                ' ARETE AUTOUR=',TAIDAR(NBZCSTCO)
               PRINT*,'arbresf1p: CERCLE de CENTRE',
     %                (XYZCZCP(K,NOBRA1),K=1,3),' RAYON=',RAYZCP(NOBRA1)
               PRINT*,'arbresf1p: CERCLE de CENTRE',
     %                (XYZCZCP(K,NOBRA2),K=1,3),' RAYON=',RAYZCP(NOBRA2)

C              LA DATE DE CREATION
               CALL ECDATE( MCN(MNXYZCEZC) )
C              LE NUMERO DU TABLEAU DESCRIPTEUR
               MCN( MNXYZCEZC + MOTVAR(6) ) = NONMTD('~>>>XYZSOMMET')

 34         ENDDO
      ENDDO

      print*,'arbresf1p: au total',NBZCSTCO,
     %       ' SOMMETS IMPOSES dans la TRIANGULATION du CERCLE',NCTRON,
     %       ' CENTRE',(XYZCZC(K,NCTRON),K=1,3),
     %       ' de RAYON AMPLIFIE=',RAYZC(NCTRON)

C     CONSTRUCTION DU MAILLAGE DE LA SURFACE PLANE
C     DU CERCLE NCTRON MOINS LES CERCLES BRAN PROJETES
C     ------------------------------------------------
C     CONSTRUCTION DU LEXIQUE de la SURFACE TRONC de CONES
C     NOM DE LA SURFACE DE CE TRONC DE CONE NOPIEC
      WRITE( KNC(1:8), '(I8)' ) NOPIEC
C     RETRAIT DES BLANCS
      CALL SANSBL( KNC, NBCAR )
      KNMPISF = 'S1CNC_' // KNC(1:NBCAR) // '_AD '

C     CONSTRUCTION DU TMS LEXIQUE de la SURFACE KNMPISF
C     SI CETTE SURFACE EXISTE, ELLE EST DETRUITE
      CALL LXLXOU( NTSURF, KNMPISF, NT1PISF, MN1PISF )
      IF( MN1PISF .GT. 0 ) CALL LXTSDS( NTSURF, KNMPISF )
      CALL LXLXDC( NTSURF, KNMPISF, 24, 8 )
      CALL LXLXOU( NTSURF, KNMPISF, NT1PISF, MN1PISF )

C     NUMERO DE LA SURFACE DANS LE LEXIQUE DES SURFACES
      CALL NUOBNM( 'SURFACE', KNMPISF, NUSFTRPI )

C     CONSTRUCTION DU TMS a_surface__definition du TRONC DE CONE
C     LA SURFACE A UN TMS DEFINITION DE SURFACE DE TYPE 9
      CALL LXTNDC( NT1PISF, 'DEFINITION', 'MOTS',
     %             WULFTR+1+NBLBRA+NBZCSTCO*2 )
      CALL LXTSOU( NT1PISF, 'DEFINITION', NT1PISFD, MN1PISFD )

C     TRANSFORMATION (I pour IDENTITE)
      MCN( MN1PISFD + WTYTRS ) = 1
C     TYPE DE LA SURFACE 9:
C     TRIANGULATION1 DE LIGNES FERMEES PAR TRIANGLES EQUILATERAUX
      MCN( MN1PISFD + WUTYSU ) = 9
C     ARETMX taille max des aretes des triangles equilateraux
C     ARETGR = DARETE EST LA VALEUR PAR DEFAUT DES ARETES DU MAILLAGE
C     REDUCTION POUR AVOIR ASSEZ DE SOMMETS PRES DES CERCLES BRANCHES
      RMCN( MN1PISFD + WRETMX ) = REAL( DARETE/4 )
C     NBLFTR 'nombre de lignes fermees contour de la surface
      MCN( MN1PISFD + WBLFTR ) = 1 + NBLBRA
C     NBPTIT 'nombre de points internes futurs sommets'
      MCN( MN1PISFD + WBPTIT ) = NBZCSTCO

C     NULFTR(1..NBLFTR) 'nom des lignes fermees'
      MCN( MN1PISFD + WULFTR ) = NULGZC( NCTRON )
      NBLBRA = 0
      DO NOBRA = 1, NBBRAN
C           NZCBRA DE 1 A NBZCER EST UN Z-CERCLE BRANCHE DE NCTRON
            NZCBRA = NCTRIZ( NDBRAN + NOBRA )
C           NUMERO DE LA LIGNE DANS LX LIGNES
            NULXLF = NULGZC( NZCBRA )
            IF( NULXLF .GT. 0 ) THEN
C              LIGNE NON REDUITE AU SOMMET d'un CONE
               NBLBRA = NBLBRA + 1
               MCN( MN1PISFD + WULFTR + NBLBRA ) = NULGZCP( NOBRA )
            ENDIF
      ENDDO

C     NUPTIT(1..NBZCSTCO) 'nom des points internes futurs sommets'
      NB = 0
      DO NOBRA = 1, NBBRAN
C           NZCBRA DE 1 A NBZCER EST UN Z-CERCLE BRANCHE DE NCTRON
            NZCBRA = NCTRIZ( NDBRAN + NOBRA )
C           NUMERO DE LA LIGNE DANS LX LIGNES
            NULXLF = NULGZC( NZCBRA )
            IF( NULXLF .LE. 0 ) THEN
C              LIGNE REDUITE AU POINT-SOMMET-CENTRE PROJETE d'un CONE
               NB = NB + 1
               MCN( MN1PISFD +WULFTR +NBLBRA +NB ) = -NULGZCP(NOBRA)
            ENDIF
      ENDDO
      DO N = NB+1, NBZCSTCO
C           POINT A MI-DISTANCE ENTRE LES 2 CERCLES
            MCN( MN1PISFD +WULFTR +NBLBRA +N ) = -NULGZCP(NBBRAN+N)
      ENDDO

C     DISOSV(1..NBPTIT) 'wished distance to the near vertices'
C     TAILLE IDEALE DES ARETES SOUHAITEES ISSUES DE CE POINT
      NB = 0
      DO NOBRA = 1, NBBRAN
C           NZCBRA DE 1 A NBZCER EST UN Z-CERCLE BRANCHE DE NCTRON
            NZCBRA = NCTRIZ( NDBRAN + NOBRA )
C           NUMERO DE LA LIGNE DANS LX LIGNES
            NULXLF = NULGZC( NZCBRA )
            IF( NULXLF .LT. 0 ) THEN
C              LIGNE REDUITE AU POINT-SOMMET-CENTRE PROJETE d'un CONE
               NB = NB + 1
               RMCN( MN1PISFD +WULFTR +NBLBRA +NBZCSTCO +NB ) =
     %               REAL( TAIDAR(NB) )
            ENDIF
      ENDDO
      DO N = NB+1, NBZCSTCO
C        POINT A MI-DISTANCE ENTRE LES 2 CERCLES
         RMCN(MN1PISFD +WULFTR +NBLBRA +NBZCSTCO +N)= REAL(TAIDAR(N))
      ENDDO

C     LA DATE DE CREATION
      CALL ECDATE( MCN(MN1PISFD) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN(MN1PISFD+MOTVAR(6))= NONMTD('~>SURFACE>>DEFINITION')

C     CONSTRUCTION DE LA TRIANGULATION DE LA SURFACE PLANE
C     1Cercle de Projection NCTRON - nCercles BRAN Projetes
C     -----------------------------------------------------
      CALL SUEX09( 9, NT1PISF, MCN(MN1PISFD), RMCN(MN1PISFD),
     %             NTNSEFS, MNNSEFS, NTXYZS, MNXYZS, IERR )

      IF( IERR .NE. 0 ) THEN
         CALL NMOBNU( 'POINT', NUPCZC(NCTRON), KNM )
         PRINT*,'arbresf1p: PIECE',NOBRA,' Surface ',KNMPISF,
     %          ' NON TRIANGULABLE dans le PLAN du Z-CERCLE ',
     %            NCTRON,' de CENTRE ',KNM
         IERR = IERR + 100
         GOTO 9900
      ENDIF

C     TRACE DE LA TRIANGULATION DANS LE PLAN DU CERCLE NCTRON
C     MISE A JOUR DES COORDONNEES EXTREMES A CELLE DE LA SURFACE
      NOTYVI = 0
      INIEXT = 0
      CALL MAJEXT( MNXYZS )
C     LA VISEE DE L'OBJET SELON LES PARAMETRES DU COMMON / TRVARI /
C     LA VISEE 3D A PARTIR DE COOEXT
      CALL VISEE0
C     LONGITUDE ET LATITUDE DE LA DIRECTION DE VISEE
      CALL LONLAT( 0.0, 84.0 )
      LORBITE = 1
      IAVFAC  = 1
C     TRACE DE LA QUALITE DES EF
      LCRITR  = 1
C     LA PALETTE DES QUALITES
      CALL PALCDE( 12 )
      CALL TRAFAC( KNMPISF, NUSFTRPI, MNNSEFS, MNXYZS )

C     NOMBRE DE SOMMETS DE LA TRIANGULATION DANS LE PLAN DE NCTRON
      NBSTSF = MCN( MNXYZS + WNBSOM )

C     TABLEAU DU NUMERO DU CERCLE DU SOMMET ou ZERO DE CHAQUE SOMMET
C     DE LA TRIANGULATION
      MNFRST = 0
      CALL TNMCDC( 'ENTIER', NBSTSF, MNFRST )
      CALL AZEROI( NBSTSF, MCN(MNFRST) )
      MNFRST1 = MNFRST - 1

C     NOMBRE DE TRIANGLES DE LA TRIANGULATION
      NBTRIA = MCN( MNNSEFS + WBEFOB )

      print*,'arbresf1p: TRIANGULATION de la SURFACE',NOPIEC,
     %       ' de NOM: ',KNMPISF,
     %       ' avec',NBTRIA,' TRIANGLES',NBSTSF,' SOMMETS'


C     RELEVEMENT DES SOMMETS DES TRIANGLES ENTRE LES CERCLES PROJETES
C     POUR PASSER EN Z DU CERCLE HAUT AUX NBBRAN CERCLES BRAN
C     ---------------------------------------------------------------

C     IDENTIFICATION DES SOMMETS DE LA TRIANGULATION
C     COMME SOMMETS DES CERCLES PROJETES BRAN PUIS
C     DEVENANT LES SOMMETS DES CERCLES BRAN
      DO 38 NOBRA = 1, NBBRAN

C           NZCBRA EST UN NO DE Z-CERCLE BRANCHE DE NCTRON NON PROJETE
            NZCBRA = NCTRIZ( NDBRAN + NOBRA )

C           NUMERO DE LA LIGNE Z-CERCLE DANS LE LEXIQUE DES LIGNES
            NULXLF = NULGZC( NZCBRA )

            IF( NULXLF .LE. 0 ) GOTO 35

C           LIGNE CERCLE PROJETE NON REDUIT A UN POINT SOMMET de CONE
C           ---------------------------------------------------------
C           RECUPERATION DU TMS XYZSOMMET DU CERCLE NON PROJETE BRANCHE NZCBRA
            MNSOL0 = MCN( MNXYZC0 -1 + NZCBRA )

C           RECUPERATION DU TMS XYZSOMMET DU CERCLE PROJETE BRANCHE NOBRA
            MNSOLP = MCN( MNXYZCP -1 + NOBRA )

C           LE NOMBRE DE SOMMETS DE LA LIGNE CERCLE BRANCHE PROJETE
            NBSOLI = MCN( MNSOLP + WNBSOM )

C           PARCOURS DES SOMMETS DU CERCLE PROJETE NOBRA
            MNSL = MNSOLP + WYZSOM -1
            DO NSL = 1, NBSOLI

C              LES 3 COORDONNEES DU SOMMET NSL DU CERCLE PROJETE NOBRA
               XNSL = RMCN( MNSL + 1 )
               YNSL = RMCN( MNSL + 2 )
               ZNSL = RMCN( MNSL + 3 )

C              RECHERCHE DU SOMMET NSL DU CERCLE PROJETE PARMI LES
C              NBSTSF SOMMETS DE LA TRIANGULATION
               DIPRO = 1E38
               NSPRO = 0
               MNST  = MNXYZS + WYZSOM -1

               DO NS = 1, NBSTSF
                  DIST = ( RMCN( MNST + 1 ) - XNSL ) **2
     %                 + ( RMCN( MNST + 2 ) - YNSL ) **2
     %                 + ( RMCN( MNST + 3 ) - ZNSL ) **2
                  IF( DIST .LT. DIPRO ) THEN
                     DIPRO = DIST
                     NSPRO = NS
                  ENDIF
                  MNST = MNST + 3
               ENDDO

C              DISTANCE MINIMALE NSL-SOMMET DE LA TRIANGULATION
               DIPRO = SQRT( DIPRO )

               IF( DIPRO .GT. 1E-3 * RAYZC(NZCBRA) ) THEN
C                 NSL SOMMET NON IDENTIFIE DANS LA TRIANGULATION
                  print*
                  print*,'arbresf1p: Probleme SOMMET',NSL,
     %                   ' DU CERCLE',NZCBRA,
     %                  ' NON IDENTIFIE a un SOMMET de la TRIANGULATION'
                  print*
                  GOTO 36
               ENDIF

C              NSL SOMMET IDENTIFIE DANS LA TRIANGULATION
C              SES XYZ DU CERCLE BRAN NON PROJETE
C              SONT IMPOSES AU SOMMET NSPRO DE LA TRIANGULATION
               MNSPRO = MNXYZS + WYZSOM + 3*NSPRO -4
               MNSL0  = MNSOL0 + WYZSOM + 3*NSL   -4

               IF( DIPRO .NE. 0. ) THEN
                  print*
                  print*,'arbresf1p: CERCLE PROJETE',NZCBRA,' de',NBSOLI
     %                  ,' SOMMETS'
                  print*,'arbresf1p: sommet NSL  =',NSL,
     %                   ' AVANT xyz=',XNSL,YNSL,ZNSL,
     %             ' a la DISTANCE MIN=',DIPRO,' du sommet NSPRO=',NSPRO
                  print*,'arbresf1p: sommet nspro=',nspro,
     %                   ' AVANT xyz=',(rmcn( mnspro + k ),k=1,3)
               ENDIF

C              IDENTIFICATION-PROJECTION DU SOMMET DU CERCLE PROJETE
C              SUR LE CERCLE INITIAL
               DO K=1,3
                  RMCN( MNSPRO + K ) = RMCN( MNSL0 + K )
               ENDDO

               IF( DIPRO .NE. 0D0 ) THEN
                  print*,'arbresf1p: sommet nspro=',nspro,
     %                   ' APRES xyz=',(rmcn( mnspro + k ),k=1,3)
               ENDIF

C              NUMERO DU CERCLE DU SOMMET NSPRO
               MCN( MNFRST1 + NSPRO ) = NZCBRA

C              PASSAGE AU SOMMET SUIVANT DE LA LIGNE NZCBRA
 36            MNSL = MNSL + 3

            ENDDO

            GOTO 38


C           LIGNE CERCLE PROJETE REDUIT A UN POINT SOMMET de CONE
C           -----------------------------------------------------
C           LES 3 COORDONNEES DU SOMMET-CERCLE PROJETE NOBRA
 35         XNSL = REAL( XYZCZCP( 1, NOBRA ) )
            YNSL = REAL( XYZCZCP( 2, NOBRA ) )
            ZNSL = REAL( XYZCZCP( 3, NOBRA ) )

C           RECHERCHE DU CENTRE PROJETE CERCLE PARMI LES
C           NBSTSF SOMMETS DE LA TRIANGULATION
            DIPRO = 1E28
            NSPRO = 0
            MNST  = MNXYZS + WYZSOM -1

            DO NS = 1, NBSTSF
               DIST = ( RMCN( MNST + 1 ) - XNSL ) **2
     %              + ( RMCN( MNST + 2 ) - YNSL ) **2
     %              + ( RMCN( MNST + 3 ) - ZNSL ) **2
               IF( DIST .LT. DIPRO ) THEN
                  DIPRO = DIST
                  NSPRO = NS
               ENDIF
               MNST = MNST + 3
            ENDDO

C           DISTANCE MINIMALE NSL-SOMMET DE LA TRIANGULATION
            DIPRO = SQRT( DIPRO )

            IF( DIPRO .GT. 1E-3 * RAYZC0 ) THEN
C              POINT-SOMMET NON IDENTIFIE DANS LA TRIANGULATION
               print*
               print*,'arbresf1p: Probleme POINT-CERCLE',NZCBRA,
     %                ' CENTRE',XNSL,YNSL,ZNSL,
     %                ' NON IDENTIFIE a un SOMMET de la TRIANGULATION'
               GOTO 38
            ENDIF
    
C           NUMERO DU CERCLE DU SOMMET NSPRO
            MCN( MNFRST1 + NSPRO ) = NZCBRA

C           RETOUR AUX XYZ DU SOMMET DU CONE
            MNST = MNXYZS + WYZSOM + 3 * NSPRO - 4
            DO K=1,3
               RMCN( MNST + K ) = XYZCZC( K, NZCBRA )
            ENDDO

            print*,'arbresf1p: sommet de CONE nspro=',nspro,
     %             ' APRES xyz=',(rmcn( mnst + k ),k=1,3),' IMPOSE'

 38   ENDDO


C     RELEVEMENT DES SOMMETS DE LA TRIANGULATION NON SUR LES CERCLES
C     TOUT SOMMET D'UN CERCLE PROJETE REDONNE LE POINT DU CERCLE BRAN
C     TOUT SOMMET DU CERCLE NCTRON EST INVARIANT
C     TOUT AUTRE SOMMET EST PROJETE SUR LE TRONC DE CONE NCTRON-NCBRANCHES
C     --------------------------------------------------------------------
      NBST0 = 0
      MNS = MNXYZS + WYZSOM -1
      DO NS = 1, NBSTSF

         IF( MCN( MNFRST1 + NS ) .GT. 0 ) GOTO 50

C        XYZ DU SOMMET NS NON SUR UN Z-CERCLE EN DOUBLE PRECISION
         XYZ( 1 ) = RMCN( MNS+1 )
         XYZ( 2 ) = RMCN( MNS+2 )
         XYZ( 3 ) = RMCN( MNS+3 )

C        LE SOMMET NS EST IL SUR LE CERCLE NCTRON?
         D = ( XYZCZC(1,NCTRON) - XYZ(1) ) **2
     %     + ( XYZCZC(2,NCTRON) - XYZ(2) ) **2
     %     + ( XYZCZC(3,NCTRON) - XYZ(3) ) **2
         D = D / ( RAYZC( NCTRON )**2 )

         IF( D .GT. 0.999D0 .AND. D .LT. 1.001D0 ) THEN

C           OUI: RETOUR AUX COORDONNEES DU SOMMET NS DU Z-CERCLE
C                NCTRON AVANT AMPLIFICATION DE SON RAYON
C           ....................................................
            NBST0 = NBST0 + 1
            DO K=1,3
               RMCN(MNS+K) = XYZCZC(K,NCTRON)
     %                + ( RMCN(MNS+K)-XYZCZC(K,NCTRON) ) / AMPLIRC(0)
            ENDDO

C           NUMERO APPARTENANCE AU CERCLE NCTRON
            MCN( MNFRST1 + NS ) = NCTRON
            GOTO 50

         ENDIF

C        NON: RECHERCHE DU CERCLE BRAN PROJETE NCPMIN
C             LE PLUS PROCHE DU SOMMET NS
C        ............................................
         NCPMIN = 0
         DMRMIN = 1D100
         DO NOBRA = 1, NBBRAN

C              DISTANCE**2 DU SOMMET NS AU CENTRE DU CERCLE PROJETE NOBRA
               D = SQRT( ( XYZCZCP(1,NOBRA) - XYZ(1) ) **2
     %                 + ( XYZCZCP(2,NOBRA) - XYZ(2) ) **2
     %                 + ( XYZCZCP(3,NOBRA) - XYZ(3) ) **2 )

C              DISTANCE - RAYON DU CERCLE PROJETE BRAN
               DMR = D - RAYZCP( NOBRA )
               IF( DMR .LT. DMRMIN ) THEN
                  DMRMIN = D
                  NCPMIN = NOBRA
               ENDIF

         ENDDO

C        LE SOMMET NS EST PROJETE SUR LE TRONC DE CONE DEFINI
C        PAR LE CERCLE NCTRON et LE CERCLE NZCBRAM NON PROJETE
         NZCBRAM = NCTRIZ( NDBRAN + NCPMIN )

         DCPMIN = ( XYZCZCP(1,NCPMIN) - XYZ(1) ) **2
     %          + ( XYZCZCP(2,NCPMIN) - XYZ(2) ) **2
     %          + ( XYZCZCP(3,NCPMIN) - XYZ(3) ) **2
         DCPMIN = SQRT( DCPMIN )

C        CALCUL DE L'ANGLE EN RADIANS NS-CENTRE DU CERCLE
C        PROJETE NCPMIN DANS LE PLAN XY DU Z-CERCLE NCTRON
         XYZ1( 1 ) = XYZCZCP( 1, NCPMIN ) + RAYZC( NCTRON )
         XYZ1( 2 ) = XYZCZCP( 2, NCPMIN )
         XYZ1( 3 ) = XYZCZCP( 3, NCPMIN )

         ANGLE  = ANGLED( XYZCZCP(1,NCPMIN), XYZ1, XYZ )
         COSANG = COS( ANGLE )
         SINANG = SIN( ANGLE )

C        LE POINT INTERSECTION DU CERCLE NCTRON ET DU SEGMENT
C        SOMMET NS - CENTRE DU CERCLE PROJETE NCPMIN EST OBTENU
C        PAR DES ITERATIONS DE DICHOTOMIE POUR CALCULER LE R
C        TEL QUE LA FONCTION F(R) SOIT NULLE
C        F(R)=( XYZCZCP(1,NCPMIN) + R * COSANG - XYZCZC(1,NCTRON) )**2
C            +( XYZCZCP(2,NCPMIN) + R * SINANG - XYZCZC(2,NCTRON) )**2
C            -  RAYZC(NCTRON)**2

C        CHOIX INITIAL DE R1 et R2 TELS QUE: F(R1)<0 et F(R2)>0
C        LE POINT SUR LE CERCLE PROJETE NCPMIN D'ANGLE ANGLE
C        EST INTERNE AU CERCLE NCTRON => F(R1)<0
C        POINT AU DELA DU CERCLE NCTRON
         RR2 = RAYZC( NCTRON ) ** 2

         R2  = 2 * RAYZC( NCTRON )
         FR2 = ( XYZCZCP(1,NCPMIN) +R2 * COSANG -XYZCZC(1,NCTRON) ) **2
     %       + ( XYZCZCP(2,NCPMIN) +R2 * SINANG -XYZCZC(2,NCTRON) ) **2
     %       -  RR2

C        RAYON DU CERCLE PROJETE LE PLUS PROCHE DE NS
         RCPMIN = RAYZCP( NCPMIN )
         R1     = RCPMIN

C        LES ITERATIONS POUR TROUVER F(R)*F(R2)<0
         ITERDICO = 0

 39      ITERDICO = ITERDICO + 1
         FR = ( XYZCZCP(1,NCPMIN) +R1 * COSANG -XYZCZC(1,NCTRON) ) **2
     %      + ( XYZCZCP(2,NCPMIN) +R1 * SINANG -XYZCZC(2,NCTRON) ) **2
     %      -  RR2

         IF( FR * FR2 .GT. 0 ) THEN

               IF( ITERDICO .GT. 16 ) THEN
                  print*,'arbresf1p: NON OBTENU R1=',R1,' et R2=',R2,
     %                   ' TELS QUE F(R1)*F(R2)<0 APRES',
     %                    ITERDICO,' ITERATIONS. PIECE=',NOPIEC
                  print*,'arbresf1p: CERCLE ENGLOBANT',NCTRON,
     %                   ' CENTRE',(XYZCZC(K,NCTRON),K=1,3),
     %                   ' RAYON',RAYZC(NCTRON)
                  print*,'arbresf1p: CERCLE PROJETE',NZCBRAM,
     %                   ' CENTRE',(XYZCZCP(K,NCPMIN),K=1,3),
     %                   ' RAYON',RAYZCP(NCPMIN)
                  GOTO 50
               ENDIF

C              R1 MAL INITIALISE
               IF( R1 .EQ. 0D0 ) THEN
                  R1 = - RAYZC( NCTRON ) / 32
               ELSE
                  R1 = -2 * R1
               ENDIF
               GOTO 39

         ENDIF

C        LES ITERATIONS DE DICHOTOMIE POUR OBTENIR F(R)=0
         ITERDICO = 0

 41      ITERDICO = ITERDICO + 1
            R = ( R1 + R2 ) / 2

            FR =( XYZCZCP(1,NCPMIN) + R * COSANG -XYZCZC(1,NCTRON) ) **2
     %         +( XYZCZCP(2,NCPMIN) + R * SINANG -XYZCZC(2,NCTRON) ) **2
     %         -  RR2

            IF( FR .EQ. 0D0 ) GOTO 42

            IF( FR * FR2 .GT. 0 ) THEN
               R2  = R
               FR2 = FR
            ELSE
               R1 = R
            ENDIF

            IF( ITERDICO .GT. 32 ) THEN
               print*,'arbresf1p: NON CONVERGENCE APRES',
     %                 ITERDICO,' ITERATIONS de DICHOTOMIE'
              print*,'arbresf1p: R1=',R1,' R2=',R2,' FR=',FR,' FR2=',FR2
               print*,'arbresf1p: SUITE AVEC R=',(R1+R2)/2
            ELSE
               IF( ( R2-R1 .GT. 1D-5 * RAYZC(NCTRON) ) ) THEN
C                 UNE ITERATION DE PLUS DE DICHOTOMIE
                  GOTO 41
               ENDIF
            ENDIF

C        CONVERGENCE ASSUREE: XYZ4 SUR LE CERCLE NCTRON
C        ET LA DROITE XYZ-CENTRE DU CERCLE PROJETE
         R = ( R1 + R2 ) / 2

 42      XYZ4( 1 ) = XYZCZCP( 1, NCPMIN ) + R * COSANG
         XYZ4( 2 ) = XYZCZCP( 2, NCPMIN ) + R * SINANG
         XYZ4( 3 ) = XYZCZCP( 3, NCPMIN )

C        DISTANCE XYZ4 ET CENTRE DU CERCLE PROJETE NCPMIN 
         D4 = ( XYZCZCP(1,NCPMIN) - XYZ4(1) ) **2
     %      + ( XYZCZCP(2,NCPMIN) - XYZ4(2) ) **2
     %      + ( XYZCZCP(3,NCPMIN) - XYZ4(3) ) **2
         D4 = SQRT( D4 )

         IF( RCPMIN .EQ. 0D0 ) THEN
C           XYZ5 EST LE SOMMET DU CONE
C           LE SOMMET NS EST SUR L'ARETE XYZ4-SOMET DU CONE
            XYZ5( 1 ) = XYZCZC( 1, NZCBRAM )
            XYZ5( 2 ) = XYZCZC( 2, NZCBRAM )
            XYZ5( 3 ) = XYZCZC( 3, NZCBRAM )
         ELSE
C           XYZ5 LE POINT INTERSECTION DU CERCLE NZCBRAM ET DU SEGMENT
C           CENTRE DU CERCLE NZCBRAM - POINT DU CERCLE NZCBRAM ET DE MEME ANGLE
            XYZ5( 1 ) = XYZCZC(1,NZCBRAM) + RAYZC(NZCBRAM) * COSANG
            XYZ5( 2 ) = XYZCZC(2,NZCBRAM) + RAYZC(NZCBRAM) * SINANG
            XYZ5( 3 ) = XYZCZC(3,NZCBRAM)
         ENDIF

C        LE SOMMET NS EST POSITIONNE SUR LA DROITE XYZ4-XYZ5
C        C-A-D SUR LE TRONC DE CONE CERCLE NCTRON- CERCLE NZCBRAM
C        AVEC LE RAPPORT
         RAPP = ( DCPMIN - RCPMIN ) / ( D4 - RCPMIN )

cccC        SI LE CERCLE NZCBRAM EST REDUIT A SON CENTRE, RAPP EST
cccC        AUGMENTE POUR RAPPROCHER LE SOMMET DE XYZ4
ccc         IF( RCPMIN .EQ. 0D0 ) THEN
cccC           EMPLOI D'UNE PARABOLE AU DESSUS DE LA DROITE (0,0)-(1,1)
cccC           F(X) = X ( a X + b )
cccC           F(0) = 0 et F(1) = 1 = a + b   => b = 1-a
cccC           F(m) = m ( a m + 1 - a ) = a m ( m - 1 ) + m
cccC           a = ( F(m) - m ) / (m ( m - 1 ) )
cccC           CHOIX m=1/4 F(m)=3/4 => a =( 3/4 - 1/4 ) / ( 1/4 ( 1/4 - 1 ) )
cccC                                        = -8 / 3
ccc            RAPP = RAPP * ( -8D0/3D0 * RAPP - 5D0/3D0 ) 
ccc         ENDIF

         IF( RAPP .LE. 1D0 ) THEN

C           LE NOUVEAU SOMMET NS NON SUR LES 2 CERCLES NCTRON NZCBRAM
            RAPP1 = 1D0 - RAPP
            DO K=1,3
               XYZ( K ) = RAPP1 * XYZ5( K ) + RAPP * XYZ4( K )
               RMCN( MNS+K ) = REAL( XYZ( K ) )
            ENDDO

         ELSE

C           LE SOMMET NS EST CHOISI SUR L'ARETE XYZ5-XYZ
            print*,'arbresf1p: RAPP=',RAPP,' >1 POURQUOI?...'
            print*,'arbresf1p: RCPMIN=',RCPMIN,' DCPMIN=',DCPMIN,
     %                        ' D4=',D4
            RAPP  = D4 / DCPMIN
            RAPP1 = 1D0 - RAPP
            DO K=1,3
               XYZ( K ) = RAPP1 * XYZ5( K ) + RAPP * XYZ4( K )
               RMCN( MNS+K ) = REAL( XYZ( K ) )
            ENDDO

         ENDIF

C        PASSAGE AU SOMMET NS SUIVANT
 50      MNS = MNS + 3

      ENDDO


C     VERIFICATION: TOUS LES SOMMETS DU Z-CERCLE NCTRON SONT ILS RETROUVES?
C     ---------------------------------------------------------------------
      NBST1 = 0
      DO NS = 1, NBSTSF
         IF( MCN( MNFRST1 + NS ) .EQ. NCTRON ) THEN
C           NS EST SUR LE CERCLE NCTRON
            NBST1 = NBST1 + 1
         ENDIF
      ENDDO

      IF( NBST0 .NE. NBARLC(NCTRON) .OR. NBST0 .NE. NBST1 ) THEN
         PRINT*
         CALL NMOBNU( 'POINT', NUPCZC(NCTRON), KNM )
         PRINT*,'arbresf1p:',NBARLC(NCTRON),' SOMMETS du CERCLE',
     %           NCTRON,' de POINT CENTRE ',KNM
         PRINT*,'arbresf1p:',NBST0,' RETROUVES DANS LA TRIANGULATION de 
     %la PIECE',NOPIEC
         PRINT*,'arbresf1p:',NBST1,' RETROUVES SUR LE CERCLE',NCTRON
         PRINT*,'arbresf1p: PROBLEME CES 3 NOMBRES DEVRAIENT ETRE EGAUX'
         PRINT*
      ENDIF

C     TRACE DE LA TRIANGULATION FINALE
C     MISE A JOUR DES COORDONNEES EXTREMES A CELLE DE LA SURFACE
      NOTYVI = 0
      INIEXT = 0
      CALL MAJEXT( MNXYZS )
C     LA VISEE DE L'OBJET SELON LES PARAMETRES DU COMMON / TRVARI /
      CALL VISEE1
      LORBITE = 1
      IAVFAC  = 1
C     TRACE DES EF
      LCRITR  = 0
C     LA PALETTE DES GRIS
      CALL PALCDE( 10 )
      CALL TRAFAC( KNMPISF, NUSFTRPI, MNNSEFS, MNXYZS )


C     BARYCENTRAGE DES SOMMETS INTERNES NON SUR UN CERCLE
C     ---------------------------------------------------
      DO NBITERBA = 1, 10
         CALL BASTMAEF( NBSTSF, RMCN(MNXYZS+WYZSOM), MCN(MNFRST),
     %                  4, NBTRIA, MCN(MNNSEFS+WUSOEF) )
      ENDDO

C     TRACE DE LA TRIANGULATION FINALE
C     MISE A JOUR DES COORDONNEES EXTREMES A CELLE DE LA SURFACE
      NOTYVI = 0
      INIEXT = 0
      CALL MAJEXT( MNXYZS )
C     LA VISEE DE L'OBJET SELON LES PARAMETRES DU COMMON / TRVARI /
      CALL VISEE1
      LORBITE = 1
      IAVFAC  = 1
C     TRACE DES EF
      LCRITR  = 0
C     LA PALETTE DES GRIS
      CALL PALCDE( 10 )
      CALL TRAFAC( KNMPISF, NUSFTRPI, MNNSEFS, MNXYZS )


C     RETOUR AU RAYON DU Z-CERCLE NCTRON INITIAL AVANT AMPLIFICATION
C     --------------------------------------------------------------
      RAYZC( NCTRON ) = RAYZC0

C     RESTAURATION DU TABLEAU XYZSOMMET DU Z-CERCLE NCTRON
C     LE NOMBRE DE SOMMETS DE LA LIGNE
      NBSOLI = MCN( MNSOLI + WNBSOM )
      MNXS1  = MNSOLI + WYZSOM
      CALL TRTATA( RMCN(MNXYZLP), RMCN(MNXS1), 3*NBSOLI )
      CALL TNMCDS( 'REEL', 3*NBSOLI, MNXYZLP )

C     SUPPRESSION DES TABLEAUX AUXILIAIRES DEVENUS INUTILES
C     -----------------------------------------------------
C     SUPPRESSION DES ADRESSES DES TMS XYZSOMMET DES CERCLES PROJETES
      IF( MNXYZCP .GT. 0 ) CALL TNMCDS( 'ENTIER', NBBRAN, MNXYZCP )

C     SUPPRESSION DU NUMERO DU CERCLE DES SOMMETS
      CALL TNMCDS( 'ENTIER', NBSTSF, MNFRST )

C     DESTRUCTION DES LIGNES CERCLES PROJETES
C     DESTRUCTION DES POINTS SOMMETS CONES PROJETES
      DO NOBRA = 1, NBBRAN

C        NOM DE LA LIGNE OU POINT
         NULX = NULGZCP( NOBRA )

         IF( NULX .GT. 0 ) THEN

C           NOM DE LA LIGNE CERCLE PROJETE DANS LE LX LIGNES
            CALL NMOBNU( 'LIGNE', NULX, KNM )
            CALL LXLXOU( NTLIGN, KNM, NT1LZC, MN1LZC )
            IF( MN1LZC .GT. 0 ) CALL LXTSDS( NTLIGN, KNM )

         ELSE IF( NULX .LT. 0 ) THEN

C           NOM DU POINT SOMMET CONE PROJETE
            CALL NMOBNU( 'POINT', -NULX, KNM )
            CALL LXLXOU( NTPOIN, KNM, NT1LZC, MN1LZC )
            IF( MN1LZC .GT. 0 ) CALL LXTSDS( NTPOIN, KNM )

         ENDIF

      ENDDO

      DO N = NBZCSTCO0+1, NBZCSTCO
C        NOM DU POINT A MI-DISTANCE ENTRE LES 2 CERCLES
         CALL NMOBNU( 'POINT', -NULGZCP(NBBRAN+N), KNM )
         CALL LXLXOU( NTPOIN, KNM, NT1LZC, MN1LZC )
         IF( MN1LZC .GT. 0 ) CALL LXTSDS( NTPOIN, KNM )
      ENDDO


C     LE MAILLAGE EST IL UNE QUAD-TRIANGULATION?
C     C-A-D CONTIENT IL DES QUADRANGLES?
C     ------------------------------------------
C     NOMBRE DE TRIANGLES-QUADRANGLES
 70   NBTRQU = MCN( MNNSEFS + WBEFOB )

C     RECHERCHE DU NOMBRE DE QUADRANGLES
      NBTRIA = 0
      NBQUAD = 0
      MNEF = MNNSEFS + WUSOEF +3
      DO N = 1, NBTRQU
         IF( MCN(MNEF) .GT. 0 ) THEN
            NBQUAD = NBQUAD + 1
         ELSE
            NBTRIA = NBTRIA + 1
         ENDIF
         MNEF = MNEF + 4
      ENDDO

      IF( NBQUAD .GT. 0 ) THEN

C        TRIANGULATION D'UNE TRIA-QUADRANGULATION
C        ----------------------------------------
         PRINT*,'arbresf1p: PIECE',NOPIEC,' NOTYPE=',NOTYPE,
     %           NBQUAD,' QUADRANGLES et',
     %           NBTRIA,' TRIANGLES => TRIANGULATION 1Q->2T'

C        NUMERO DE LA SURFACE DANS LE LEXIQUE DES SURFACES
         CALL NUOBNM( 'SURFACE', KNMPISF, NUSFTRPI )

C        NOM DE LA SURFACE AVEC SEULEMENT DES TRIANGLES
         CALL SANSBL(  KNMPISF, NBCAR )
         KNMTRIA = KNMPISF(1:NBCAR) // '_TRIA '

C        CONSTRUCTION DU TMS LEXIQUE de la SURFACE KNMTRIA
C        SI CETTE SURFACE EXISTE, ELLE EST DETRUITE
         CALL LXLXOU( NTSURF, KNMTRIA, NTTRIA, MNTRIA )
         IF(MNTRIA.GT.0) CALL LXTSDS( NTSURF, KNMTRIA )
         CALL LXLXDC( NTSURF, KNMTRIA, 24, 8 )
         CALL LXLXOU( NTSURF, KNMTRIA, NTTRIA, MNTRIA )

C        CONSTRUCTION DU TMS a_surface__definition de la PIECE TRIANGULEE
C        LA SURFACE A UN TMS DEFINITION DE SURFACE DE TYPE 31
         CALL LXTNDC(NTTRIA, 'DEFINITION', 'MOTS', WDECOU+1 )
         CALL LXTSOU(NTTRIA, 'DEFINITION', NTTRIAD, MNTRIAD)

C        TRANSFORMATION (I pour IDENTITE)
         MCN( MNTRIAD + WTYTRS ) = 1
C        TYPE DE LA SURFACE 32: TRIANGULATION D'UNE QUADRANGULATION
         MCN( MNTRIAD + WUTYSU ) = 31

C        NUSUQU nom de la surface quadrangulee a trianguler
         MCN( MNTRIAD + WUSUQU ) = NUSFTRPI

C        NDECOU 'code decoupage des quadrangles'
C        ( 1 : 'creation d''une ARETE S1 S3'
C        , 2 : 'creation d''une ARETE S2 S4' 
C        , 3 : 'COUPE du PLUS GRAND ANGLE du quadrangle'
C        , 4 : 'MAX du MINIMUM des qualites des triangles' )
         MCN( MNTRIAD + WDECOU ) = 4

C        LA DATE DE CREATION
         CALL ECDATE( MCN(MNTRIAD) )

C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN(MNTRIAD+MOTVAR(6))=NONMTD('~>SURFACE>>DEFINITION')

C        GENERER LA TRIANGULATION D'UNE QUADRANGULATION
         CALL SUEX31( NTTRIA, MCN( MNTRIAD ),
     %                NTNSEFS, MNNSEFS, NTXYZS, MNXYZS, IERR )

C        DESTRUCTION DE LA SURFACE AVEC QUADRANGLES
         CALL LXLXOU( NTSURF, KNMPISF, NTTR, MNTR )
         IF(MNTR.GT.0) CALL LXTSDS( NTSURF, KNMPISF )

C        NUMERO DE LA SURFACE TRIANGULEE DANS LE LEXIQUE DES SURFACES
         CALL NUOBNM( 'SURFACE', KNMTRIA, NUSFTRPI )

         NBTRQU = MCN( MNNSEFS + WBEFOB )
         PRINT*,'arbresf1p: PIECE',NOPIEC,' NOTYPE=',NOTYPE,
     %           NBTRQU,' TRIANGLES'

      ENDIF


C     EVENTUELLE CREATION DES POINTS A IMPOSER DANS LA TETRAEDRISATION
C     LE LONG DE SEGMENTS D'EXTREMITES LES CENTRES DES CERCLES ET
C     EN FONCTION de ARETGR
C     ================================================================
      IF( NOPTIM .EQ. 0 ) GOTO 9999

C     CALCUL DU BARYCENTRE XYZBAR DE LA TRIANGULATION
      DO K=1,3
         XYZBAR( K ) = 0
      ENDDO

C     NBTRIA Nombre de TRIANGLES DE LA SURFACE TRIANGULEE
      NBTRIA = MCN( MNNSEFS + WBEFOB )
      MNTR   = MNNSEFS + WUSOEF - 1
      DO N=1,NBTRIA
C        LE TRIANGLE N
         DO I=1,3
C           LE NO DE SOMMET I DU TRIANGLE N
            NS  = MCN( MNTR + I )
            MNS = MNXYZS + WYZSOM + 3*NS - 4
            DO K=1,3
               XYZBAR( K ) = XYZBAR( K ) + RMCN( MNS + K )
            ENDDO
         ENDDO
         MNTR = MNTR + 4
      ENDDO
      DO K=1,3
         XYZBAR( K ) = XYZBAR( K ) / (3*NBTRIA)
      ENDDO

C     CALCUL DE TAILLE_IDEALE AU BARYCENTRE XYZBAR
      CALL TAILIDEA( NOFOTI, XYZBAR, NCODEV, XYZBAR(4) )

C     HAUTEUR D'UN TETRAEDRE EQUILATERAL D'ARETE XYZBAR(4)
      ZHTE = XYZBAR(4) * SQRT( 2D0 / 3D0 )

C     LES POINTS DANS DES PLANS Z=Cte LE LONG DES NBBRAN BRANCHES
C     -----------------------------------------------------------
      DO NOBRA = 1, NBBRAN

C        NZCBRA NUMERO DE 1 A NBZCER EST UN Z-CERCLE BRANCHE
         NZCBRA = NCTRIZ( NDBRAN + NOBRA )

C        DISTANCE ENTRE LES CENTRES DES 2 CERCLES NCTRON et NZCBRA
         D = ( XYZCZC(1,NZCBRA) - XYZCZC(1,NCTRON) ) **2
     %     + ( XYZCZC(2,NZCBRA) - XYZCZC(2,NCTRON) ) **2
     %     + ( XYZCZC(3,NZCBRA) - XYZCZC(3,NCTRON) ) **2
         D = SQRT( D )

C        NOMBRE DE FOIS DE L'ARETE SUR CE SEGMENT JOIGNANT LES 2 CENTRES
C        SOIT LE NOMBRE DE COUCHES A Z=CTE ENTRE LES 2 CENTRES
         NBCOUCHZ = INT( D / ZHTE ) + 1

C        LONGUEUR DE L'ARETE SUR LE CERCLE NZCBRA
         R = 2 * Pi * RAYZC( NZCBRA ) / NBARLC( NZCBRA )

C        LONGUEUR DE L'ARETE SUR L'ARETE DES 2 CENTRES
         TAR = MIN( DARETE, (DARETE+R)/2, 2*R )

C        NOMBRE D'ARETES DARETE DANS CETTE HAUTEUR ZH
         NBARZH = INT( D / TAR ) + 1

C        NOMBRE DE FOIS DE L'ARETE SUR CE SEGMENT JOIGNANT LES 2 CENTRES
C        SOIT LE NOMBRE DE COUCHES A Z=CTE ENTRE LES 2 CENTRES
         NBCOUCHZ = MAX( NBCOUCHZ, NBARZH )

         IF( NBCOUCHZ .GT. 1 ) THEN
            DO 80 NOCOUCH = 1, NBCOUCHZ-1

C              LE POINT NOCOUCH A IMPOSER SUR LE SEGMENT JOIGNANT LES 2 CENTRES
C              ----------------------------------------------------------------
               R1 = DBLE( NOCOUCH ) / DBLE( NBCOUCHZ )
               DO N=1,3
                  XYZD( N ) = XYZCZC( N, NCTRON ) + R1
     %                      * ( XYZCZC(N,NZCBRA) - XYZCZC(N,NCTRON) )
               ENDDO
               XYZD( 4 ) = TAR

C              CALCUL DE TAILLE_IDEALE( XYZD )
               IF( NOFOTI .GT. 0 ) THEN
C                 LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
                  CALL FONVAL( NOFOTI, 3, XYZD, NCODEV, XYZD(4) )
                  IF( NCODEV .LE. 0 ) THEN
                     XYZD( 4 ) = DARETE
                  ENDIF
               ENDIF
C              INTERDICTION D'UNE VALEUR NEGATIVE OU NULLE
C              EN DISTANCE SOUHAITEE
               IF( XYZD(4) .LE. 0D0 ) THEN
                  XYZD(4) = DARETE
               ENDIF

C              EXISTE T IL UN POINT IMPOSE TROP PROCHE DE
C              CELUI A IMPOSER?
               DO NPAJ = 1, NBPTIM
                  R = ( XYZDPTIM( 1, NPAJ ) - XYZD( 1 ) ) **2
     %              + ( XYZDPTIM( 2, NPAJ ) - XYZD( 2 ) ) **2
     %              + ( XYZDPTIM( 3, NPAJ ) - XYZD( 3 ) ) **2
                  R = SQRT( R )
                  IF( R .LT. 0.7D0 * XYZD(4) ) THEN
C                    ABANDON DU POINT XYZD
                     GOTO 80
                  ENDIF
               ENDDO

C              AJOUT DE XYZD AUX NBPTIM POINTS IMPOSES
               IF( NBPTIM .GT. MXPTIM ) THEN
                  GOTO 9100
               ENDIF
               NBPTIM  = NBPTIM + 1
               NBPTIM1 = NBPTIM
               DO N=1,4
                  XYZDPTIM( N, NBPTIM1 ) = REAL( XYZD( N ) )
               ENDDO


               IF( NBPTIM .GT. 0 ) GOTO 80
C              CI-DESSOUS INACTIF DURANT LA MISE AU POINT


C              AJOUT DES 6 SOMMETS D'HEXAGONES DE CENTRE CE POINT
C              --------------------------------------------------
C              LE RAYON INTERPOLE DANS LE PLAN Z ENTRE LES 2 CERCLES
               RAYONI = ( NOCOUCH * RAYZC(NCTRON)
     %                + (NBCOUCHZ-NOCOUCH) * RAYZC(NZCBRA) ) / NBCOUCHZ

               IF( RAYONI .LT. 1.26 * XYZDPTIM(4,NBPTIM1) ) THEN
C                 ABANDON DU POINT XYZD
                  GOTO 80
               ENDIF

C              NBPSR NOMBRE DE POINTS A IMPOSER DANS LE RAYON INTERPOLE
C              ET DANS LE PLAN Z DU POINT CENTRE SUR LE SEGMENT
               NBPSR = INT( RAYONI / XYZDPTIM(4,NBPTIM1) ) - 1
               IF( NBPSR .LE. 0 ) GOTO 80

C              AJOUT DES 6 SOMMETS DES NBPSR HEXAGONES
               IF( NBPTIM+6*NBPSR .GT. MXPTIM ) THEN
                  GOTO 9100
               ENDIF

C              POUR UNE ROTATION ET DECALAGE DE COUCHE EN COUCHE
               NOCOUCHMOD2 = MOD( NOCOUCH-1, 2 )

               DO 77 KK = 1, NBPSR

C                 LE POINT KK A IMPOSER
                  RAYONK  = ( DBLE(KK) - NOCOUCHMOD2 / 2D0 ) * XYZD(4)
                  DECANGL = NOCOUCHMOD2 * PIS6

                  DO M=1,6

C                    LE SOMMET M DE L'HEXAGONE
                     ANGLE = (M-1) * PIS3 + DECANGL
                     XYZDR( 1 ) = XYZD( 1 ) + RAYONK * COS( ANGLE )
                     XYZDR( 2 ) = XYZD( 2 ) + RAYONK * SIN( ANGLE )
                     XYZDR( 3 ) = XYZD( 3 )
                     XYZDR( 4 ) = XYZD( 4 )

C                    CALCUL DE TAILLE_IDEALE
                     IF( NOFOTI .GT. 0 ) THEN
C                       LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
                        CALL FONVAL( NOFOTI, 3, XYZDR, NCODEV, XYZDR(4))
                        IF( NCODEV .LE. 0 ) THEN
                           XYZDR( 4 ) = DARETE
                        ENDIF
                     ENDIF
C                    INTERDICTION D'UNE VALEUR NEGATIVE OU NULLE
C                    EN DISTANCE SOUHAITEE
                     IF( XYZDR(4) .LE. 0D0 ) THEN
                        XYZDR(4) = DARETE
                     ENDIF

cccC                    EXISTE T IL UN POINT IMPOSE TROP PROCHE
cccC                    DE CELUI A IMPOSER?
ccc                     DO NPAJ = 1, NBPTIM
ccc                        R = ( XYZDPTIM( 1, NPAJ ) - XYZDR( 1 ) ) **2
ccc     %                    + ( XYZDPTIM( 2, NPAJ ) - XYZDR( 2 ) ) **2
ccc     %                    + ( XYZDPTIM( 3, NPAJ ) - XYZDR( 3 ) ) **2
ccc                        R = SQRT( R )
ccc                        IF( R .LT. 0.7D0 * XYZD(4) ) THEN
cccC                          ABANDON DU POINT XYZD
ccc                           GOTO 80
ccc                        ENDIF
ccc                     ENDDO

C                    AJOUT EFFECTIF DU POINT A LA LISTE
                     NBPTIM  = NBPTIM + 1
                     DO N=1,4
                        XYZDPTIM( N, NBPTIM )= REAL( XYZDR( N ) )
                     ENDDO

                  ENDDO

 77            ENDDO


C              AJOUT DES POINTS INTERIEURS AU TRIANGLE LIMITE PAR
C              2 RAYONS DE L'HEXAGONE MAXIMAL
C              --------------------------------------------------
               IF( NBPTIM + 6*NBPSR*(NBPSR-1)/2 .GT. MXPTIM ) THEN
                  GOTO 9100
               ENDIF

C              A PARTIR DU SECOND HEXAGONE AJOUTE
               DO NHEX = 2, NBPSR

                  DO MTR=1,6

C                    TRIANGLE MTR DE L'HEXAGONE MAXIMAL
C                    NO DU POINT AJOUTE SUR LE RAYON1
                     NBP1 = NBPTIM1 + (NHEX-1) * 6 + MTR
C                    NO DU POINT AJOUTE SUR LE RAYON2
                     IF( MTR .NE. 6 ) THEN
                        NBP2 = NBP1 + 1
                     ELSE
                        NBP2 = NBPTIM1 + (NHEX-1) * 6 + 1
                     ENDIF

                     DO 78 NPSAR = 1, NHEX-1

                        R = DBLE( NPSAR ) / DBLE( NHEX )
                        DO N=1,4
                           XYZDR( N ) = XYZDPTIM( N, NBP1 ) 
     %                         + R * (XYZDPTIM(N,NBP2)-XYZDPTIM(N,NBP1))
                        ENDDO

C                       CALCUL DE TAILLE_IDEALE
                        IF( NOFOTI .GT. 0 ) THEN
C                          LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
                           CALL FONVAL( NOFOTI, 3, XYZDR,
     %                                  NCODEV, XYZDR(4))
                           IF( NCODEV .LE. 0 ) THEN
                              XYZDR( 4 ) = DARETE
                           ENDIF
                        ENDIF
C                       INTERDICTION D'UNE VALEUR NEGATIVE
C                       OU NULLE EN DISTANCE SOUHAITEE
                        IF( XYZDR(4) .LE. 0D0 ) THEN
                           XYZDR(4) = DARETE
                        ENDIF

C                       EXISTE T IL UN POINT IMPOSE TROP PROCHE
C                       DE CELUI A IMPOSER?
                        DO NPAJ = 1, NBPTIM

                           R = ( XYZDPTIM( 1, NPAJ ) - XYZDR( 1 ) ) **2
     %                       + ( XYZDPTIM( 2, NPAJ ) - XYZDR( 2 ) ) **2
     %                       + ( XYZDPTIM( 3, NPAJ ) - XYZDR( 3 ) ) **2
                           R = SQRT( R )

                           IF( R .LT. 0.7D0 * XYZD(4) ) THEN
C                             ABANDON DU POINT XYZD
                              GOTO 78
                           ENDIF

                        ENDDO

C                       AJOUT EFFECTIF DU POINT A LA LISTE
                        NBPTIM = NBPTIM + 1
                        DO N=1,4
                           XYZDPTIM( N, NBPTIM )= REAL( XYZDR( N ) )
                        ENDDO

 78                  ENDDO

                  ENDDO

               ENDDO

 80         ENDDO
         ENDIF

      ENDDO

 900  PRINT*
      PRINT*,'arbresf1p: PIECE',NOPIEC,' NOTYPE=',NOTYPE,' :',
     %        NBPTIM,' POINTS IMPOSES AJOUTES'
      DO K=1, NBPTIM
         print*,'arbresf1p: Pt IMPOSE a TETRAEDRISER',K,
     %          ': XYZD=',(XYZDPTIM(N,K),N=1,4)
      ENDDO
      GOTO 9999


 9100 print*,'arbresf1p: SATURATION du TABLEAU XYZDPTIM'
      print*,'arbresf1p: AUGMENTER MXPTIM=',MXPTIM
      GOTO 900


 9900 PRINT*,'arbref1p: PIECE',NOPIEC,' NOTYPE=',NOTYPE,
     %' SORTIE avec IERR=',IERR,' => PAS de TRIANGULATION'

 9999 RETURN
      END

