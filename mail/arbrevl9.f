      SUBROUTINE ARBREVL9( NBZCER, NUPCZC, XYZCZC, RAYZC,NBARLC,
     %                     NBPIEC, NTZCPI, NBZCPI, NUZCPI,
     %                     NULGZC, NUSFZC, NCTRIZ,
     %                     NOPTIM, MXPTIM, NBPTIM, XYZDPTIM,
     %                     NUVLTEPI, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER UNE TETRAEDRISATION DU VOLUME D'UN ARBRE
C -----    DEFINI PAR NBPIEC PIECES (TRONCONS) DEFINIES PAR LA LISTE
C          DE SES Z-CERCLES
C          ATTENTION: UN Z-CERCLE PEUT ETRE REDUIT A UN POINT SON CENTRE
C                     DANS CE CAS LE NOMBRE DE SES ARETES DOIT ETRE EGAL A 0
C          LE TYPE D'UNE PIECE PEUT ETRE
C          1: CONE SOMMET BAS  (defini par 1 Pt=ZC + 1 ZC)
C          2: CONE SOMMET HAUT (defini par 1 ZC + 1 Pt=ZC)
C          3: TRONC DE CONE    (defini par 2 ZC)
C          4: n Z-CERCLES RELIES a UN Z-CERCLE  (defini par n + 1 ZC)
C          5: UN Z-CERCLE RELIE  a n  Z-CERCLES (defini par 1 + n ZC)

C ENTREES:
C --------
C NOFOTI : NUMERO DANS LE LX FONCTIONS DE LA FONCTION TAILLE_IDEALE(X,Y,Z)
C NBZCER : NOMBRE DE Z-CERCLES DE L'ARBRE
C NUPCZC : NUMERO LX POINTS DES NBZCER CENTRES DES CERCLES
C XYZCZC : XYZ DU CENTRE DES NBZCER Z-CERCLES
C RAYZC  : RAYON DES NBZCER Z-CERCLES
C NBARLC : >2 NOMBRE D'ARETES DE CHACUNE DES NBZCER CERCLES
C          =0 INDIQUE UN CERCLE REDUIT A SON POINT CENTRE

C NBPIEC : nombre de PIECES ou TRONCONS de l'ARBRE
C NTZCPI : nombre total de Z-CERCLES decrivant les PIECES de l'ARBRE
C NBZCPI : nombre de Z-CERCLES de chaque PIECE de l'ARBRE
C NUZCPI : NUMERO 1 a NBZCER du Z-CERCLE de chaque PIECE
C MXPTIM : MAXIMUM DE POINTS A TETRAEDRISER STOCKABLES DANS XYZDPTIM

C AUXILIAIRES:
C ------------
C NMLGZC : NOM DES NBZCER CERCLES LIGNES ARETISEES
C NULGZC : NUMERO DANS LE LEXIQUE LIGNES   DES NBZCER CERCLES
C NUSFZC : NUMERO DANS LE LEXIQUE SURFACES DES NBZCER CERCLES
C NCTRIZ : NUMERO DE 1 A NBZCER DES Z-CERCLES D'UNE PIECE DE L'ARBRE
C          A TRIER SELON Z CROISSANT DE LEUR CERCLE

C SORTIES:
C --------
C NOPTIM   : =1 AJOUT DE POINTS A IMPOSER DANS UNE FUTURE TETRAEDRISATION
C NBPTIM   : NOMBRE DE POINTS A TETRAEDRISER ENSUITE INITIAUX
C XYZDPTIM : X Y Z Distance Souhaitee aux sommets voisins DES POINTS
C            INITIAUX A TETRAEDRISER ENSUITE
C NUVLTEPI : NO DANS LE LEXIQUE VOLUMES DES NBPIEC VOLUMES TETRAEDRISES
C IERR     : =0 SI PAS D'ERREUR DETECTEE, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  Saint PIERRE du PERRAY             Juin 2019
C23456...............................................................012
      PARAMETER      (MXZCER=256, QTEAME=0.08)
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_point__definition.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/darete.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(MOTMCN)
      EQUIVALENCE       (MCN(1),RMCN(1))

      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE,TRACTE0

      REAL              RAYZC(NBZCER), XYZCZC(3,NBZCER),
     %                  XYZDPTIM(4,MXPTIM)
      INTEGER           NUPCZC(NBZCER), NBARLC(NBZCER)

      INTEGER           NBZCPI(NBPIEC), NUZCPI(NTZCPI)

      INTEGER           NULGZC(NBZCER), NUSFZC(NBZCER), NCTRIZ(NBZCER),
     %                  NUSFTRPI, NUVLTEPI(NBPIEC)

      CHARACTER*8       KNC
      CHARACTER*24      KNMLGZC, KNMSFZC, KNM, KNMSFPI, KNMSFPA,
     %                  KNMPTIM, KNMV8PI, KNMV9PI

      IERR    = 0
      TRACTE0 = TRACTE
      MNXYZCP = 0
      MNFRST  = 0

C     NOMBRE DE POINTS IMPOSES A TETRAEDRISER ENSUITE
      NBPTIM = 0

C     ANGL2P = angle en degres de coplanearite entre 2 faces adjacentes
ccc      ANGL2P = 45.
ccc      ANGL2P = 31.
      ANGL2P = 31.

      PRINT*
      PRINT 10000
      PRINT*,'arbrevl9: TETRAEDRISATION du VOLUME d''un ARBRE a partir d
     %e',NBPIEC,' PIECES et',NBZCER,' CERCLES d''un plan Z=Cte'
      PRINT 10000

C     CONSTRUCTION du LX de la LIGNE et du TMS DEFINITION et XYZSOMMET
C     DES ARETES DES NBZCER Z-CERCLES DECLARES
C     ================================================================
      CALL TNMCDC( 'ENTIER', NBZCER, MNXYZC0 )
      DO NC = 1, NBZCER

C        CERCLE NC A TRAITER ou REDUIT AU POINT CENTRE?
         IF( NBARLC( NC ) .EQ. 0 ) THEN

            CALL NMOBNU( 'POINT', NUPCZC(NC), KNM )
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'arbrevl9: ARETISATION du CERCLE',NC,
     %                ' REDUIT au POINT CENTRE: ',KNM,
     %                ' NON FAITE'
            ELSE
               PRINT*,'arbrevl9: EDGES of CIRCLE',NC,
     %                ' REDUCED to one POINT CENTRE: ',KNM
            ENDIF

C           NUMERO DE LA LIGNE DANS LE LX DE LA LIGNE NC
C           DEVIENT LE NUMERO DU POINT CENTRE SOMMET D'UN CONE
            NULGZC( NC ) = -NC

         ELSE IF( NBARLC( NC ) .GE. 3 ) THEN

C           AU MOINS 3 SOMMETS DEMANDES SUR LE CERCLE NC
            IF( RAYZC( NC ) .LE. 0. ) THEN
               PRINT*,'arbrevl9: Z-CERCLE',NC,' avec', NBARLC(NC),
     %                ' ARETES et un RAYON',RAYZC(NC),' INCOHERENT'
               IERR = 10
               GOTO 900
            ENDIF

C           KNMLGZC NOM DE LA LIGNE MOMENTANEE DU CERCLE NC
            WRITE( KNC(1:8), '(I8)' ) NC
C           RETRAIT DES CARACTERES BLANCS
            CALL SANSBL( KNC, NBCAR )
            KNMLGZC = 'LZC_' // KNC(1:NBCAR) // '_AD '

C           CONSTRUCTION DU TMS LEXIQUE de la LIGNE KNMLGZC
C           SI CETTE LIGNE EXISTE, ELLE EST DETRUITE
            CALL LXLXOU( NTLIGN, KNMLGZC, NT1LZC, MN1LZC )
            IF( MN1LZC .GT. 0 ) CALL LXTSDS( NTLIGN, KNMLGZC )
            CALL LXLXDC( NTLIGN, KNMLGZC, 24, 8 )
            CALL LXLXOU( NTLIGN, KNMLGZC, NT1LZC, MN1LZC )

C           NUMERO DE LA LIGNE DANS LE LX DE LA LIGNE NC
            CALL NUOBNM( 'LIGNE', KNMLGZC, NULGZC(NC) )

C           CONSTRUCTION DU TMS a_ligne__definition du CERCLE NC
C           LA LIGNE A UN TMS DEFINITION DE LIGNE DE TYPE 8:
C           CERCLE de R3 de TYPE NUTYCI=3 DEFINI PAR
C           LE CENTRE, RAYON, PLAN X ou Y ou Z = Cte
            CALL LXTNDC( NT1LZC, 'DEFINITION', 'MOTS', WUPLCT+1  )
            CALL LXTSOU( NT1LZC, 'DEFINITION', NT1CDE, MN1CDE )

C           TRANSFORMATION (I pour IDENTITE)
            MCN( MN1CDE + WTYTRL ) = 1
C           TYPE DE LA LIGNE: 8: CERCLE DE R3
            MCN( MN1CDE + WUTYLI ) = 8
C           NOMBRE D'ARETES DU CERCLE NC
            MCN( MN1CDE + WBARLI ) = NBARLC( NC )
C           NUTYCI numero du type du cercle
            MCN( MN1CDE + WUTYCI ) = 3
C           NUPTCE numero du point centre du cercle nc
            MCN( MN1CDE + WUPTCE ) = NUPCZC( NC )
C           RAYDCI rayon du cercle nc
            RMCN( MN1CDE + WAYDCI ) = RAYZC( NC )
C           NUPLCT numero du plan XY a Z=Cte=CENTRE(3)
            MCN( MN1CDE + WUPLCT ) = 3
C           LA DATE DE CREATION
            CALL ECDATE( MCN(MN1CDE) )
C           LE NUMERO DU TABLEAU DESCRIPTEUR
            MCN( MN1CDE + MOTVAR(6) ) = NONMTD('~>LIGNE>>DEFINITION')

C           CONSTRUCTION DU MAILLAGE EN ARETES DU CERCLE NC
            CALL LIEX08( NT1LZC,  MCN(MN1CDE), RMCN(MN1CDE),
     %                   NTNSEFC, MNNSEFC, NTXYZC, MNXYZC, IERR )

C           SAUVEGARDE DE L'ADRESSE DU TMS XYZSOMMET DU CERCLE NC
            MCN( MNXYZC0 -1 + NC ) = MNXYZC

            IF( IERR .NE. 0 ) THEN
               IERR = IERR + 10
               PRINT*,'arbrevl9: ERREUR dans le MAILLAGE du CERCLE',NC
     %                  ,KNMLGZC
               GOTO 900
            ENDIF

C           SUPPRESSION DES TANGENTES AU CERCLE
            MCN( MNXYZC  + WNBTGS ) = 0
            MCN( MNNSEFC + WBTGEF ) = 0
            MCN( MNNSEFC + WBEFAP ) = 0
            MCN( MNNSEFC + WBEFTG ) = 0

            IF( LANGAG .EQ. 0 ) THEN
               print*,'arbrevl9: ARETISATION du CERCLE',NC,' de NOM: ',
     %                 KNMLGZC,' avec',MCN(MNNSEFC+WBEFOB),' ARETES',
     %                 MCN(MNXYZC+WNBSOM),' SOMMETS est la LIGNE',
     %                 NULGZC(NC)
            ELSE
               print*,'arbrevl9: EDGES of CIRCLE',NC,
     %                ' of NAME: ',KNMLGZC,
     %                ' with',MCN(MNNSEFC+WBEFOB),' EDGES and',
     %                 MCN(MNXYZC+WNBSOM),' VERTICES is the LINE',
     %                 NULGZC(NC)
            ENDIF

         ELSE IF( NBARLC( NC ) .GT. 0 ) THEN

            PRINT *,' arbrevl9: CERCLE', NC,' avec',NBARLC(NC),
     %              ' SOMMETS DEMANDES au LIEU de =0 ou >2'
            IERR = 1
            GOTO 900

ccc      ELSE
cccC        CERCLE REDUIT A UN CENTRE=POINT=SOMMET

         ENDIF

      ENDDO


C     CONSTRUCTION du LX de la SURFACE et des TMS DEFINITION XYZSOMMET
C     NSEF de la TRIANGULATION des NBZCER Z-CERCLES DECLARES
C     ================================================================
      NBSOMMZC = 0
      NBTRIAZC = 0
      DO NC = 1, NBZCER

C        CERCLE NC REDUIT AU POINT CENTRE?
         IF( NBARLC( NC ) .LE. 0 ) THEN

C           OUI: CERCLE NC REDUIT AU POINT CENTRE  PAS DE TRIANGULATION
            CALL NMOBNU( 'POINT', NUPCZC(NC), KNM )
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'arbrevl9: TRIANGULATION du CERCLE',NC,
     %                ' REDUITE au POINT CENTRE: ',KNM
            ELSE
               PRINT*,'arbrevl9: TRIANGULATION of CIRCLE',NC,
     %                ' REDUCED to one POINT CENTRE: ',KNM
            ENDIF

C           NUMERO DE LA SURFACE DANS LE LX DE LA SURFACE NC
C           DEVIENT LE NUMERO DU POINT CENTRE SOMMET D'UN CONE
            NUSFZC( NC ) = -NC

         ELSE

C           KNMSFZC NOM DE LA SURFACE MOMENTANEE DU CERCLE NC
            WRITE( KNC(1:8), '(I8)' ) NC
C           RETRAIT DES CARACTERES BLANCS
            CALL SANSBL( KNC, NBCAR )
            KNMSFZC = 'SZC_' // KNC(1:NBCAR) // '_AD '

C           CONSTRUCTION DU TMS LEXIQUE de la SURFACE KNMSFZC
C           SI CETTE SURFACE EXISTE, ELLE EST DETRUITE
            CALL LXLXOU( NTSURF, KNMSFZC, NT1SZC, MN1SZC )
            IF( MN1SZC .GT. 0 ) CALL LXTSDS( NTSURF, KNMSFZC )
            CALL LXLXDC( NTSURF, KNMSFZC, 24, 8 )
            CALL LXLXOU( NTSURF, KNMSFZC, NT1SZC, MN1SZC )

C           NUMERO DE LA SURFACE DANS LE LX DE LA SURFACE NC
            CALL NUOBNM( 'SURFACE', KNMSFZC, NUSFZC(NC) )

C           CONSTRUCTION DU TMS a_surface__definition du Z-CERCLE
C           LA SURFACE A UN TMS DEFINITION DE SURFACE DE TYPE 9
            CALL LXTNDC( NT1SZC, 'DEFINITION', 'MOTS', WULFTR+1 )
            CALL LXTSOU( NT1SZC, 'DEFINITION', NT1SZCD, MN1SZCD )

C           TRANSFORMATION (I pour IDENTITE)
            MCN( MN1SZCD + WTYTRS ) = 1
C           TYPE DE LA SURFACE 9:
C           TRIANGULATION1 DE LIGNES FERMEES PAR TRIANGLES EQUILATERAUX
            MCN( MN1SZCD + WUTYSU ) = 9
C           ARETMX taille max des aretes des triangles equilateraux
C           ARETGR = DARETE EST LA VALEUR PAR DEFAUT DES ARETES DU MAILLAGE
C           REDUCTION POUR AVOIR ASSEZ DE SOMMETS DANS LE Z-CERCLE
ccc            RMCN( MN1SZCD + WRETMX ) = REAL( DARETE/4 )
            RMCN( MN1SZCD + WRETMX ) = REAL( DARETE/3 )
C           NBLFTR 'nombre de lignes fermees contour de la surface
            MCN( MN1SZCD + WBLFTR ) = 1
C           NBPTIT 'nombre de points internes futurs sommets'
            MCN( MN1SZCD + WBPTIT ) = 0
C           NULFTR(1..NBLFTR) 'nom des lignes fermees'
            MCN( MN1SZCD + WULFTR ) = NULGZC(NC)

C           LA DATE DE CREATION
            CALL ECDATE( MCN(MN1SZCD) )
C           LE NUMERO DU TABLEAU DESCRIPTEUR
            MCN( MN1SZCD + MOTVAR(6) ) = NONMTD('~>SURFACE>>DEFINITION')

C           CONSTRUCTION DU MAILLAGE EN ARETES DU CERCLE NC
            CALL SUEX09( 9, NT1SZC,  MCN(MN1SZCD), RMCN(MN1SZCD),
     %                   NTNSEFC, MNNSEFC, NTXYZC, MNXYZC, IERR )

            IF( IERR .NE. 0 ) THEN
               IERR = IERR + 10
               PRINT*,'arbrevl9: ERREUR dans la TRIANGULATION du CERCLE'
     %               ,NC,' de NOM ',KNMSFZC
               GOTO 900
            ENDIF

C           SUPPRESSION EVENTUELLE DES TANGENTES AU CERCLE
            MCN( MNXYZC  + WNBTGS ) = 0
            MCN( MNNSEFC + WBTGEF ) = 0
            MCN( MNNSEFC + WBEFAP ) = 0
            MCN( MNNSEFC + WBEFTG ) = 0

            NBST     = MCN(MNXYZC+WNBSOM)
            NBSOMMZC = NBSOMMZC + NBST
            NBTRIA   = MCN(MNNSEFC+WBEFOB)
            NBTRIAZC = NBTRIAZC + NBTRIA

            print*,'arbrevl9: La SURFACE du Z-CERCLE',NC,
     %             ' de NOM de LIGNE ',KNMLGZC,
     %             ' est TRIANGULEE avec',NBST,' SOMMETS',
     %             ' et',NBTRIA,' TRIANGLES. Son NOM est ',KNMSFZC

         ENDIF

      ENDDO

      PRINT*
      PRINT*,'arbrevl9: Les',NBZCER,' Z-CERCLES sont TRIANGULES avec',
     %NBSOMMZC,' sommets et',NBTRIAZC,' triangles'


C     CREATION DE LA TRIANGULATION DE LA SURFACE TOTALE DE CHAQUE PIECE
C     =================================================================
C     NO DU DERNIER Z-CERCLE AVANT CEUX DE LA PIECE A TRAITER
      NDPIEC = 0

C     CREATION EVENTUELLE DE POINTS INTERNES IMPOSES A LA TETRAEDRISATION
      NOPTIM = 1
10000 FORMAT(166('@'))

      DO NOPIEC = 1, NBPIEC

         PRINT*
         PRINT 10000
         PRINT*,'arbrevl9: Tetraedrisation (voex09) de la PIECE',NOPIEC
         PRINT 10000

C        NOMBRE DE Z-CERCLES DE LA PIECE NOPIEC
         NBZCP = NBZCPI( NOPIEC )

         IF( NBZCP .LT. 2 ) THEN
            PRINT*,'arbrevl9:',NBZCP,'<2 NOMBRE INCORRECT DE Z-CERCLES d
     %e la PIECE',NOPIEC
            IERR = 2
            GOTO 900
         ENDIF

C        RECHERCHE DU TYPE NOTYPE DE LA PIECE
C        1: CONE SOMMET BAS
C        2: CONE SOMMET HAUT
C        3: TRONC DE CONE
C        4: n Z-CERCLES RELIES a UN Z-CERCLE
C        5: UN Z-CERCLE RELIE  a n  Z-CERCLES
C        ET TRIANGULATION DE SA SURFACE LIMITEE PAR NBZCP Z-CERCLES

         CALL ARBRESF1P( NBZCER, NUPCZC, XYZCZC, RAYZC, NBARLC,
     %                   NULGZC, MNXYZC0,
     %                   NOPIEC, NBZCP,  NUZCPI(NDPIEC+1), NCTRIZ,
     %                   NOPTIM, MXPTIM, NBPTIM, XYZDPTIM,
     %                   NUSFTRPI, NOTYPE, IERR )
         IF( IERR .NE. 0 ) GOTO 900


C        CONSTRUCTION DE LA TRIANGULATION DE LA PIECE NOPIEC PAR UNION
C        DE SES SURFACES LATERALE ET NBZCP SECTIONS Z-CERCLES
C        -------------------------------------------------------------
C        NBZCP  NOMBRE DE Z-CERCLES DE LA PIECE NOPIEC
C        NDPIEC NO DU DERNIER Z-CERCLE AVANT CEUX DE LA PIECE A TRAITER
C        NCTRIZ TABLEAU DU NO 1 A NBZCER DES NBZCP Z-CERCLES DE LA PIECE

C        KNMSFPI NOM DE LA SURFACE UNION DE LA PIECE NOPIEC
         WRITE( KNC(1:8), '(I8)' ) NOPIEC
C        RETRAIT DES CARACTERES BLANCS
         CALL SANSBL( KNC, NBCAR )
         KNMSFPI = 'PIS_' // KNC(1:NBCAR) // '_SV89 '

C        CONSTRUCTION DU TMS LEXIQUE de la SURFACE KNMSFPI
C        SI CETTE SURFACE EXISTE, ELLE EST DETRUITE
         CALL LXLXOU( NTSURF, KNMSFPI, NTSFPI, MNSFPI )
         IF( MNSFPI .GT. 0 ) CALL LXTSDS( NTSURF, KNMSFPI )
         CALL LXLXDC( NTSURF, KNMSFPI, 24, 8 )
         CALL LXLXOU( NTSURF, KNMSFPI, NTSFPI, MNSFPI )

C        NUMERO DE LA SURFACE UNION DE LA PIECE DANS LE LX DES SURFACES
         CALL NUOBNM( 'SURFACE', KNMSFPI, NUSFPI )

C        CONSTRUCTION DU TMS a_surface__definition
C        LA SURFACE A UN TMS DEFINITION DE SURFACE DE TYPE 51
         NBVRZC = 0
         DO N = 1, NBZCP
            NOZC = NCTRIZ( N )
            IF( NBARLC( NOZC ) .GT. 0 ) THEN
               NBVRZC = NBVRZC + 1
               NCTRIZ( NBVRZC ) = NOZC
            ENDIF
         ENDDO

C        NOMBRE DE VRAIES SURFACES (NON SOMMET DE CONE)
         NBSUUN = NBVRZC + 1 

         CALL LXTNDC( NTSFPI, 'DEFINITION', 'MOTS', WUSUUN+NBSUUN )
         CALL LXTSOU( NTSFPI, 'DEFINITION', NTSFPID, MNSFPID )
C        TRANSFORMATION (I pour IDENTITE)
         MCN( MNSFPID + WTYTRS ) = 1
C        TYPE DE LA SURFACE 51: UNION DE PLUSIEURS SURFACES EN G0-CONTINUITE
         MCN( MNSFPID + WUTYSU ) = 51
C        NBSUUN 'nombre de surfaces de l'union'
         MCN( MNSFPID + WBSUUN ) = NBSUUN
C        NUSUUN(1..NBSUUN) 'nom des surfaces'
         DO N = 1, NBVRZC
C           NUMERO DE SURFACE TRIANGULEE DU Z-CERCLE DE LA PIECE NOPIEC
            MCN( MNSFPID + WUSUUN -1 + N ) = NUSFZC( NCTRIZ( N ) )
         ENDDO
C        NUMERO SURFACE TRIANGULEE LATERALE DE LA PIECE NOPIEC
         MCN( MNSFPID + WUSUUN -1 + NBSUUN ) = NUSFTRPI
C        LA DATE DE CREATION
         CALL ECDATE( MCN(MNSFPID) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNSFPID + MOTVAR(6) ) = NONMTD('~>SURFACE>>DEFINITION')

C        CONSTRUCTION DU MAILLAGE DE LA SURFACE TOTALE DE LA PIECE NOPIEC
C        PAR UNE UNION de SES NBZCP+1 SURFACES
         CALL UNPLSV( 51,   3, NUSFPI, NBSUUN, MCN(MNSFPID+WUSUUN),
     %                NTNSEFS, MNNSEFS, NTXYZSS, MNXYZSS,
     %                NTUNIOS, MNUNIOS, IERR )

         IF( IERR .EQ. 0 ) THEN
C           IMPRESSION DE LA QUALITE DE LA TRIANGULATION
            CALL IMPQUA( 3, KNMSFPI, MNNSEFS, MNXYZSS, NBEFMQ, QUAMIN,
     %                   SURVOLEF )
         ENDIF

C        TRACE DU MAILLAGE DE LA SURFACE NUSFPI UNION DES SURFACES
C        MISE A JOUR DES COORDONNEES EXTREMES A CELLE DE LA SURFACE
         NOTYVI = 0
         INIEXT = 0
         CALL MAJEXT( MNXYZSS )
         CALL EFFACE
C        LA VISEE DE L'OBJET SELON LES PARAMETRES DU COMMON / TRVARI /
C        LA VISEE 3D A PARTIR DE COOEXT
         CALL VISEE0
C        LONGITUDE ET LATITUDE DE LA DIRECTION DE VISEE
         CALL LONLAT( 0.0, 20.0 )
         LORBITE = 1
         IAVFAC  = 1
C        TRACE DE LA QUALITE DES EF
         LCRITR  = 1
C        LA PALETTE DES QUALITES
         CALL PALCDE( 12 )
         PREDUF = 5.
         IF( TRACTE ) THEN
            CALL TRAFAC( KNMSFPI, NUSFPI, MNNSEFS, MNXYZSS )
         ENDIF


C        CONSTRUCTION DE LA TRIANGULATION AMELIOREE DE LA SURFACE
C        SI SON TYPE EST 4 ou 5 (CONE et TRONC de CONE NON AMELIORE)
C        -----------------------------------------------------------
         IF( NOTYPE .LT. 4 ) THEN
C           PAS D'AMELIORATION DE LA TRIANGULATION DU CONE ou TRONC de CONE
            KNMSFPA = KNMSFPI 
            NUSFPA  = NUSFPI
            GOTO 50
         ENDIF

C        SURFACE DE TYPE NOTYPE
C        4: n Z-CERCLES RELIES a UN Z-CERCLE
C        5: UN Z-CERCLE RELIE  a n  Z-CERCLES

C        KNMSFPA NOM DE LA SURFACE UNION DE LA PIECE NOPIEC
         WRITE( KNC(1:8), '(I8)' ) NOPIEC
C        RETRAIT DES CARACTERES BLANCS
         CALL SANSBL( KNC, NBCAR )
         KNMSFPA = 'PIS_' // KNC(1:NBCAR) // '_S9A '

C        CONSTRUCTION DU TMS LEXIQUE de la SURFACE KNMSFPA
C        SI CETTE SURFACE EXISTE, ELLE EST DETRUITE
         CALL LXLXOU( NTSURF, KNMSFPA, NTSFPA, MNSFPA )
         IF( MNSFPA .GT. 0 ) CALL LXTSDS( NTSURF, KNMSFPA )
         CALL LXLXDC( NTSURF, KNMSFPA, 24, 8 )
         CALL LXLXOU( NTSURF, KNMSFPA, NTSFPA, MNSFPA )

C        NUMERO DE LA SURFACE UNION DE LA PIECE DANS LE LX DES SURFACES
         CALL NUOBNM( 'SURFACE', KNMSFPA, NUSFPA )

         CALL LXTNDC( NTSFPA, 'DEFINITION', 'MOTS', WNGL2P+1 )
         CALL LXTSOU( NTSFPA, 'DEFINITION', NTSFPAD, MNSFPAD )
C        TRANSFORMATION (I pour IDENTITE)
         MCN( MNSFPAD + WTYTRS ) = 1
C        TYPE DE LA SURFACE 30: AMELIORATION DE LA QUALITE D'UNE TRIANGULATION 3D
         MCN( MNSFPAD + WUTYSU ) = 30
C        NUSUQU nom de la triangulation a ameliorer ^~>SURFACE ;
         MCN( MNSFPAD + WUSUQU ) = NUSFPI
C        QUALMN 'qualite au dessous de laquelle ameliorer' reel;
         RMCN( MNSFPAD + WUALMN ) = QTEAME

C        ANGL2P 'degres du petit angle de coplanearite' reel;
         RMCN( MNSFPAD + WNGL2P ) = ANGL2P

C        LA DATE DE CREATION
         CALL ECDATE( MCN(MNSFPAD) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNSFPAD + MOTVAR(6) ) = NONMTD('~>SURFACE>>DEFINITION')

         CALL SUEX30( NUSFPA,  NTSFPA,  MCN(MNSFPAD), RMCN(MNSFPAD),
     %                NTNSEFS, MNNSEFS, NTXYZSS, MNXYZSS, IERR )

         IF( IERR .EQ. 0 ) THEN
C           IMPRESSION DE LA QUALITE DE LA TRIANGULATION
            CALL IMPQUA( 3, KNMSFPA, MNNSEFS, MNXYZSS, NBEFMQ, QUAMIN,
     %                   SURVOLEF )
         ENDIF

ccc         IF( MCN( MNNSEFS + WUTFMA ) .NE. 1 ) THEN
cccC           SURFACE TRIANGULEE NON FERMEE
         TRACTE = .TRUE.
ccc         ENDIF
         IF( TRACTE ) THEN
C           TRACE DE LA TRIANGULATION DE LA SURFACE AMELIOREE
            CALL TRAFAC( KNMSFPA, NUSFPA, MNNSEFS, MNXYZSS )
         ENDIF
         TRACTE = TRACTE0


C        CETTE SURFACE NUSFPA LIMITE UN VOLUME DE TYPE 8
C        -----------------------------------------------
C        KNMV8PI NOM DU VOLUME DE LA PIECE NOPIEC de TYPE 8
 50      WRITE( KNC(1:8), '(I8)' ) NOPIEC
C        RETRAIT DES CARACTERES BLANCS
         CALL SANSBL( KNC, NBCAR )
         KNMV8PI = 'PIV_' // KNC(1:NBCAR) // '_V8 '
C        CONSTRUCTION DU TMS LEXIQUE du volume KNMV8PI
C        SI CE VOLUME EXISTE, IL EST DETRUIT
         CALL LXLXOU( NTVOLU, KNMV8PI, NTV8PI, MNV8PI )
         IF( MNV8PI .GT. 0 ) CALL LXTSDS( NTVOLU, KNMV8PI )
         CALL LXLXDC( NTVOLU, KNMV8PI, 24, 8 )
         CALL LXLXOU( NTVOLU, KNMV8PI, NTV8PI, MNV8PI )
C        NUMERO DU VOLUME TETRAEDRISE DANS LE LX DES VOLUMES
         CALL NUOBNM( 'VOLUME', KNMV8PI, NUV8PI )

C        CONSTRUCTION DU TMS a_volume__definition de la PIECE NOPIEC V8
C        LA VOLUME A UN TMS DEFINITION DE VOLUME DE TYPE 8
         CALL LXTNDC( NTV8PI, 'DEFINITION', 'MOTS', WUSF1V+1 )
         CALL LXTSOU( NTV8PI, 'DEFINITION', NTV8PID, MNV8PID )
C        TRANSFORMATION (I pour IDENTITE)
         MCN( MNV8PID + WTYTRV ) = 1
C        TYPE DU VOLUME 8: VOLUME FERME POUR TETRAEDRISATION
         MCN( MNV8PID + WUTYVO ) = 8
C        NBSF1V NUMBER of closed surfaces of the volume
         MCN( MNV8PID + WBSF1V ) = 1
C        NUSF1V(1..NBSF1V) 'NAMES of the closed surfaces' ^~>SURFACE ;
         MCN( MNV8PID + WUSF1V ) = NUSFPA
C        LA DATE DE CREATION
         CALL ECDATE( MCN(MNV8PID) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNV8PID + MOTVAR(6) ) = NONMTD('~>VOLUME>>DEFINITION')


C        TETRAEDRISATION V9 DE LA PIECE NOPIEC A PARTIR DE LA TRIANGULATION
C        DE SA SURFACE KNMSFPA NUSFPA FRONTIERE DU VOLUME V8
C        ------------------------------------------------------------------
C        KNMV9PI NOM DU VOLUME DE LA PIECE NOPIEC de TYPE 9
         WRITE( KNC(1:8), '(I8)' ) NOPIEC
C        RETRAIT DES CARACTERES BLANCS
         CALL SANSBL( KNC, NBCAR )
         KNMV9PI = 'PIV_' // KNC(1:NBCAR) // '_V9 '
C        CONSTRUCTION DU TMS LEXIQUE du volume KNMV9PI
C        SI CE VOLUME EXISTE, IL EST DETRUIT
         CALL LXLXOU( NTVOLU, KNMV9PI, NTV9PI, MNV9PI )
         IF( MNV9PI .GT. 0 ) CALL LXTSDS( NTVOLU, KNMV9PI )
         CALL LXLXDC( NTVOLU, KNMV9PI, 24, 8 )
         CALL LXLXOU( NTVOLU, KNMV9PI, NTV9PI, MNV9PI )

C        NUMERO DU VOLUME TETRAEDRISE DANS LE LX DES VOLUMES
         CALL NUOBNM( 'VOLUME', KNMV9PI, NUV9PI )
         NUVLTEPI( NOPIEC ) = NUV9PI

         IF( MCN( MNNSEFS + WUTFMA ) .NE. 1 ) THEN
C           SURFACE TRIANGULEE NON FERMEE DU VOLUME -> PAS DE TETRAEDRISATION
            PRINT*
            PRINT*,'arbrevl9:',KNMSFPA,' est une SURFACE TRIANGULEE NON 
     %FERMEE'
            PRINT*,'arbrevl9:',KNMV9PI,' son VOLUME N''EST PAS TETRAEDR
     %ISE'
            PRINT*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            PRINT*
            GOTO 800
         ENDIF

C        CONSTRUCTION DU TMS a_volume__definition de la PIECE NOPIEC V9
C        LA VOLUME A UN TMS DEFINITION DE VOLUME DE TYPE 9
         CALL LXTNDC( NTV9PI, 'DEFINITION', 'MOTS', WUVOPA+1+2*NBPTIM )
         CALL LXTSOU( NTV9PI, 'DEFINITION', NTV9PID, MNV9PID )
C        TRANSFORMATION (I pour IDENTITE)
         MCN( MNV9PID + WTYTRV ) = 1
C        TYPE DU VOLUME 9: TETRAEDRISATION1 DE VOLUMES D'OPTION 8
         MCN( MNV9PID + WUTYVO ) = 9
C        MXSOVO nombre maximal de sommets de la tetraedrisation
         MCN( MNV9PID + WXSOVO ) = 100000
C        ARETGR longueur de l''arete de la grille reguliere
         RMCN( MNV9PID + WRETGR ) = REAL( DARETE )
C        ANEN2P 'max (degrees) angle of coplanearity
         RMCN( MNV9PID + WNEN2P ) = ANGL2P
C        NBPTIV nombre de points internes fournis par l''utilisateur
         MCN( MNV9PID + WBPTIV ) = NBPTIM
C        NBVOPA nombre de volumes de type 8 de la partition
         MCN( MNV9PID + WBVOPA ) = 1
C        NUVOPA(1..NBVOPA) nom du volume 8 du volume v9 ^~>VOLUME
         MCN( MNV9PID + WUVOPA ) = NUV8PI

C        tableau NUPTIV(1..NBPTIV) et DISSOV(1..NBPTIV)
         DO K=1,NBPTIM

C           CONSTRUCTION DU POINT IMPOSE XYZDPTIM(K) DANS LA TETRAEDRISATION
C           CONSTRUCTION DU LEXIQUE du POINT de NOM KNMPTIM
            WRITE( KNC(1:8), '(I8)' ) K
C           RETRAIT DES BLANCS
            CALL SANSBL( KNC, NBCAR )
            KNMPTIM = 'PTIM_' // KNC(1:NBCAR) // '_AD '

C           CONSTRUCTION DU TMS LEXIQUE du POINT KNMPTIM
C           SI CE POINT EXISTE, IL EST DETRUIT
            CALL LXLXOU( NTPOIN, KNMPTIM, NT1PTIM, MN1PTIM )
            IF(MN1PTIM.GT.0) CALL LXTSDS( NTPOIN, KNMPTIM )
            CALL LXLXDC( NTPOIN, KNMPTIM, 24, 8 )
            CALL LXLXOU( NTPOIN, KNMPTIM, NT1PTIM, MN1PTIM )

C           NUMERO DU POINT K DANS LE LEXIQUE DES POINTS
            CALL NUOBNM( 'POINT', KNMPTIM, NUPTIM )

C           CONSTRUCTION de son TMS a_point__definition
C           LE POINT A UN TMS DEFINITION DE POINT DE TYPE 1
            CALL LXTNDC( NT1PTIM, 'DEFINITION', 'MOTS', WOORPO+2 )
            CALL LXTSOU( NT1PTIM, 'DEFINITION', NT1PTIMD, MN1PTIMD )
C           TRANSFORMATION (I pour IDENTITE)
            MCN( MN1PTIMD + WTYTRP ) = 1
C           TYPE DU POINT 1:  XYZ FOURNIS
            MCN( MN1PTIMD + WUTYPO ) = 1
C           XYZ
            RMCN( MN1PTIMD + WOORPO    ) = XYZDPTIM( 1, K )
            RMCN( MN1PTIMD + WOORPO +1 ) = XYZDPTIM( 2, K )
            RMCN( MN1PTIMD + WOORPO +2 ) = XYZDPTIM( 3, K )
C           LA DATE DE CREATION
            CALL ECDATE( MCN(MN1PTIMD) )
C           LE NUMERO DU TABLEAU DESCRIPTEUR
            MCN( MN1PTIMD + MOTVAR(6) ) = NONMTD('~>POINT>>DEFINITION')

C           CONSTRUCTION DU TMX XYZSOMMET du POINT NUPTIM
            CALL POEX01( NT1PTIM,  MCN(MN1PTIMD),
     %                   NTSOPTIM, MNSOPTIM, IERR )

C           tableau NUPTIV(1..NBPTIV)
C          'noms des points internes au volume partition' ^~>POINT ;
            MCN( MNV9PID + WUVOPA + K ) = NUPTIM

C           tableau DISSOV(1..NBPTIV)
C          'distance souhaitee aux sommets voisins' reel ;
            RMCN( MNV9PID + WUVOPA + NBPTIM + K ) = XYZDPTIM( 4, K )

         ENDDO

C        LA DATE DE CREATION
         CALL ECDATE( MCN(MNV9PID) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNV9PID + MOTVAR(6) ) = NONMTD('~>VOLUME>>DEFINITION')

C        TETRAEDRISATION TYPE 9 DE VOLUMES D'OPTION 8
         CALL VOEX09( NUV9PI,  NTV9PI, MCN(MNV9PID), RMCN(MNV9PID),
     %                NTNSEFV, MNNSEFV, NTXYZSV, MNXYZSV, IERR )

C        TRACE DES EF DU VOLUME NUV9PI SI PAS D'ERREUR
         IF( IERR .EQ. 0 ) THEN

C           AFFICHER L'HISTOGRAMME DES QUALITES DES EF
C           LA QUALITE MOYENNE, MINIMALE ET L'ECART TYPE A 1
            CALL IMPQUA( 4, KNMV9PI, MNNSEFV, MNXYZSV, NBEFMQ, QUAMIN,
     %                   SURVOLEF )

C           MISE A JOUR DES COORDONNEES EXTREMES A CELLE DE LA SURFACE
            NOTYVI = 0
            INIEXT = 0
            CALL MAJEXT( MNXYZSV )
C           LA VISEE DE L'OBJET SELON LES PARAMETRES DU COMMON / TRVARI /
C           LA VISEE 3D A PARTIR DE COOEXT
            CALL VISEE0
C           LONGITUDE ET LATITUDE DE LA DIRECTION DE VISEE
            CALL LONLAT( 0., 22. )
            LORBITE = 1
            IAVFAC  = 1
C           TRACE DE LA QUALITE DES EF
            LCRITR  = 1
C           LA PALETTE DES QUALITES
            CALL PALCDE( 12 )
            PREDUF = 12.
            tracte =.true.
            IF( TRACTE ) THEN
               CALL TRACUB( KNMV9PI, NUV9PI, MNNSEFV, MNXYZSV )
            ENDIF
            tracte = tracte0

         ELSE

C           ERREUR RENCONTREE -> PAS DE TETRAEDRISATION

         PRINT*,'///////////////////////////////////////////////////////
     %/////////////////////////////////////////////////////////////////'
            IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'arbrevl9: le VOLUME ',KNMV9PI,' N''EST PAS TETRAEDRISE'
            ELSE
         PRINT*,'arbrevl9: the VOLUME ',KNMV9PI,' IS NOT TETRAHEDRIZED'
            ENDIF
         PRINT*,'///////////////////////////////////////////////////////
     %/////////////////////////////////////////////////////////////////'
            PRINT*

         ENDIF

C        DESTRUCTION DES NBPTIM POINTS IMPOSES DANS LA TETRAEDRISATION
         DO K = 1, NBPTIM
C          'nom du point internes au volume partition' ^~>POINT ;
            NUPTIM = MCN( MNV9PID + WUVOPA + K )
C           NOM DU POINT K DANS LE LEXIQUE DES POINTS
            CALL NMOBNU( 'POINT', NUPTIM, KNMPTIM )
C           DESTRUCTION DU LEXIQUE DU POINT K IMPOSE
            CALL LXTSDS( NTPOIN, KNMPTIM )
         ENDDO

C        PASSAGE A LA PIECE SUIVANTE
 800     NDPIEC = NDPIEC + NBZCP 

      ENDDO


C     SUPPRESSION DES ADRESSES DES TMS XYZSOMMET DES CERCLES
C     ------------------------------------------------------
 900  IF( MNXYZC0.GT.0 ) CALL TNMCDS( 'ENTIER', NBZCER, MNXYZC0 )

C     DESTRUCTION DES LIGNES et SURFACES des NBZCER Z-CERCLES
C     -------------------------------------------------------
      DO NC = 1, NBZCER

C        NOM DE LA LIGNE DU CERCLE NC
         NULX = NULGZC( NC )
         IF( NULX .GT. 0 ) THEN
C           NOM DE LA LIGNE CERCLE DANS LE LX LIGNES
            CALL NMOBNU( 'LIGNE', NULX, KNM )
            CALL LXLXOU( NTLIGN, KNM, NT1LZC, MN1LZC )
            IF( MN1LZC .GT. 0 ) CALL LXTSDS( NTLIGN, KNM )
         ENDIF

C        NOM DE LA SURFACE DU CERCLE NC
         NULX = NUSFZC( NC )
         IF( NULX .GT. 0 ) THEN
C           NOM DE LA SURFACE CERCLE DANS LE LX SURFACES
            CALL NMOBNU( 'SURFACE', NULX, KNM )
            CALL LXLXOU( NTSURF, KNM, NT1SZC, MN1SZC )
            IF( MN1SZC .GT. 0 ) CALL LXTSDS( NTSURF, KNM )
         ENDIF

      ENDDO

      TRACTE = TRACTE0
      RETURN
      END
