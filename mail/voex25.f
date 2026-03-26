      SUBROUTINE VOEX25( NUARBR, LADEFI, RADEFI,
     %                   NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :      GENERER UNE TETRAEDRISATION D'UN ARBRE ET DE SES RACINES
C -----      A PARTIR DE TRONCS DE CONE DEFINIS PAR DES CERCLES

C ENTREES:
C --------
C NUARBR   : NUMERO DU VOLUME ARBRE DANS LE LEXIQUE DES VOLUMES
C LADEFI   : TABLEAU ENTIER DE DEFINITION DU VOLUME ARBRE
C RADEFI   : TABLEAU REEL   DE DEFINITION DU VOLUME ARBRE
C            CF '~/td/d/a_volume__definition'

C SORTIES:
C --------
C NTNSEF : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES
C MNNSEF : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES
C            CF '~/td/d/a___nsef'
C NTXYZS : NUMERO      DU TMS 'XYZSOMMET' DU VOLUME
C MNXYZS : ADRESSE MCN DU TMS 'XYZSOMMET' DU VOLUME
C            CF '~/td/d/a___xyzsommet'
C IERR     : = 0 SI PAS D'ERREUR
C            <>0 SI ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  Saint PIERRE du PERRAY             Juin 2019
C23456...............................................................012
      PARAMETER      (MXZCER=256, MXPTIM=8192)
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/darete.inc"

      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER          LADEFI(0:*)
      REAL             RADEFI(0:*)

C     DECLARATION MAXIMALE DES TABLEAUX NECESSAIRES
      REAL             RAYZCV(MXZCER), XYZCZC(3,MXZCER),
     %                 XYZDPTIM(4,MXPTIM)
      INTEGER          NUPCZC(MXZCER), NBARLC(MXZCER)
      INTEGER          NULGZC(MXZCER), NUSFZC(MXZCER), NCTRIZ(MXZCER), 
     %                 NUVLTEPI(MXZCER)
      CHARACTER*24     KNM
      DOUBLE PRECISION DEUXPI

      DEUXPI = ATAN( 1D0 ) * 8D0

C     RECUPERATION DES DONNEES DU TMS a_volume__definition
C     -----------------------------------------------------
      IERR     = 0
      MNSTCO   = 0
      MNNUZCPI = 0
      MNLADEFI = 0

C     VOLUME d''UN ARBRE et SES RACINES a PARTIR de CERCLES: TYPE 25
ccc   NUTYVO = LADEFI( WUTYVO )

C     NBPIEC 'nombre de PIECES ou TRONCONS de l''ARBRE'
      NBPIEC = LADEFI( WBPIEV )

C     NTZCPI 'nombre total de Z-CERCLES decrivant les PIECES'
      NTZCPI = LADEFI( WTZCPV )

C     NBZCER 'nombre de lignes Z-CERCLES'
      NBZCER = LADEFI( WBZCER )

      IF( NBZCER .LT. 4 ) THEN
         PRINT *,'voex25: ARBRE avec NBZCER=',NBZCER,
     %           '<4 un NOMBRE de CERCLES INSUFFISANT'
         IERR = 8
         GOTO 9999
      ENDIF

      IF( NBZCER .GT. MXZCER ) THEN
         PRINT *,'voex25: AUGMENTER MXZCER=',MXZCER,
     %           ' AU DELA de NBZCER=',NBZCER
         IERR = 9
         GOTO 9999
      ENDIF

C     EXISTENCE OU NON DE LA FONCTION 'TAILLE_IDEALE'
C     -----------------------------------------------
      NOFOTI = NOFOTIEL()
      IF( LANGAG .EQ. 0 ) THEN
         IF( NOFOTI .GT. 0 ) THEN
C           LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
            WRITE(IMPRIM,10000) 'AVEC'
         ELSE
C           LA FONCTION TAILLE_IDEALE(X,Y,Z) N'EXISTE PAS
            WRITE(IMPRIM,10000) 'SANS'
            WRITE(IMPRIM,10001) DARETE
         ENDIF
      ELSE
         IF( NOFOTI .GT. 0 ) THEN
C           LA FONCTION EDGE_LENGTH(X,Y,Z) EXISTE
            WRITE(IMPRIM,20000) '  WITH '
         ELSE
C           LA FONCTION EDGE_LENGTH(X,Y,Z) N'EXISTE PAS
            WRITE(IMPRIM,20000) 'WITHOUT'
            WRITE(IMPRIM,20001) DARETE
         ENDIF
      ENDIF
10000 FORMAT('voex25: TETRAEDRISATION ',A,
     %       ' FONCTION TAILLE_IDEALE(X,Y,Z) (ou EDGE_LENGTH(X,Y,Z))' )
20000 FORMAT('voex25: TETRAHEDRISATION ',A,
     %       ' FUNCTION EDGE_LENGTH(X,Y,Z) (or TAILLE_IDEALE(X,Y,Z))' )
10001 FORMAT(' LONGUEUR de l''ARETE SOUHAITEE =',G14.6)
20001 FORMAT(' WISHED EDGE LENGTH =',G14.6)


C     LES TABLEAUX COPIES DE LADEFI SONT CREES POUR POUVOIR ETRE MODIFIES
c     tableau  NUPCZV(1..NBZCER) 'nom des centres des Z-cercles' ^~>POINT ;
      MN = WUPZCV -1
      DO NC = 1, NBZCER

C        NUMERO DU POINT-CENTRE DANS LE LEXIQUE DES POINTS
         NUPCZC( NC ) = LADEFI( MN + NC )

C        RECUPERATION DES 3 COORDONNEES DU POINT CENTRE DU CERCLE NC
         CALL LXNLOU( NTPOIN, NUPCZC( NC ), NTLXZC, MNP )
         IF( NTLXZC .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'POINT INCONNU'
            CALL LEREUR
            IERR = 6
            GOTO 9999
         ENDIF

C        RESTAURATION DU TMS XYZSOMMET DU POINT NUPCZC( NC )
         CALL LXTSOU( NTLXZC, 'XYZSOMMET', NTXYZC, MNXYZC )
         IF( NTXYZC .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'POINT SANS SOMMET'
            CALL LEREUR
            IERR = 7
            GOTO 9999
         ENDIF

C        LES 3 COORDONNEES DU POINT CENTRE DU Z-CERCLE NC
         XYZCZC( 1, NC ) = RMCN( MNXYZC + WYZSOM )
         XYZCZC( 2, NC ) = RMCN( MNXYZC + WYZSOM + 1 )
         XYZCZC( 3, NC ) = RMCN( MNXYZC + WYZSOM + 2 )

      ENDDO
 
c     tableau  RAYZCV(1..NBZCER)  'rayon des cercles' reel ;
      MN = MN + NBZCER
      DO NC = 1, NBZCER

C        RAYON DU CERCLE NC
         RAYZCV( NC ) = RADEFI( MN + NC )

C        CALCUL DU NOMBRE D'ARETES DU CERCLE NC EN FONCTION DE
C        LA LONGUEUR DE L'ARETE PAR DEFAUT CONNUE DARETE
         IF( RAYZCV( NC ) .EQ. 0D0 ) THEN
C           CERCLE REDUIT A SON POINT CENTRE: PAS D'ARETE
            NB = 0
         ELSE
            NB = NINT( DEUXPI * RAYZCV( NC ) / DARETE )
            NB = MAX( 7, NB )
         ENDIF
         NBARLC( NC ) = NB

      ENDDO

C     tableau NBZCPI(1..NBPIEC) 'nombre de Z-CERCLES de chaque PIECE'
      MNNBZCPI = MN + NBZCER + 1

C     tableau PZCVPI(1..NTZCPI) 'nom du CENTRE des Z-CERCLES de chaque PIECE'
      MNPZCVPI = MNNBZCPI + NBPIEC

C     tableau numero de 1 a NBZCER de chaque Z-CERCLE des PIECES de l'ARBRE
      CALL TNMCDC( 'ENTIER', NTZCPI, MNNUZCPI )
      MN = MNNUZCPI
      DO 10 K = 1, NTZCPI

C        NUMERO DU K-EME POINT CENTRE DU Z-CERCLE DANS LE LEXIQUE POINT
         NUZC = LADEFI( MNPZCVPI-1+K )

         DO NC = 1, NBZCER
C           NUMERO DU POINT CENTRE DANS LE LEXIQUE POINT
            IF( NUZC .EQ. NUPCZC( NC ) ) THEN
C               NUZC EST LE CENTRE DU Z-CERCLE NC
               MCN( MN ) = NC
               MN = MN + 1
               GOTO 10
            ENDIF
         ENDDO

C        CE POINT NUZC N'EST PAS UN CENTRE DE Z-CERCLES
         CALL NMOBNU( 'POINT', NUZC, KNM )
         PRINT*,'voex25: LE POINT ',KNM,
     %          ' N''EST PAS UN POINT CENTRE DES Z-CERCLES'
         IERR = 10

 10   ENDDO

      IF( IERR .NE. 0 ) GOTO 9999

C     CONSTRUCTION DU VOLUME DE L'ARBRE et de SES RACINES
C     GENERATION DE LA TETRAEDRISATION DES NBPIEC PIECES
C     ---------------------------------------------------
      NOPTIM = 1
      CALL ARBREVL9( NBZCER, NUPCZC, XYZCZC, RAYZCV, NBARLC,
     %               NBPIEC, NTZCPI, LADEFI(MNNBZCPI), MCN(MNNUZCPI),
     %               NULGZC, NUSFZC, NCTRIZ,
     %               NOPTIM, MXPTIM, NBPTIM, XYZDPTIM,
     %               NUVLTEPI, IERR )
      IF( IERR .NE. 0 ) GOTO 9999


C     UNION DES NBPIEC VOLUMES TRONCS DE CONE CREEES ENTRE CERCLES
C     ------------------------------------------------------------
      IF( NBPIEC .GE. 2 ) THEN

C        CONSTRUCTION DU MAILLAGE DU VOLUME UNION des TRONCS DE CONE
         CALL UNPLSV( 51,  4, NUARBR, NBPIEC, NUVLTEPI,
     %                NTNSEF, MNNSEF, NTXYZS, MNXYZS,
     %                NTUNIO, MNUNIO, IERR )

C        DESTRUCTION DU VOLUME DES NBPIEC PIECES
         DO N = 1, NBPIEC

C           NUMERO DU VOLUME
            NULX = NUVLTEPI( N )
            IF( NULX .GT. 0 ) THEN
C              NOM DU VOLUME DANS LE LEXIQUE DES VOLUMES
               CALL NMOBNU( 'VOLUME', NULX, KNM )
               CALL LXLXOU( NTVOLU, KNM, NTSF, MNSF )
ccc               IF( MNSF .GT. 0 ) CALL LXTSDS( NTVOLU, KNM )
            ENDIF

         ENDDO

      ENDIF


C     1 VOLUME MULTI-MATERIAUX, UNION DE PLUSIEURS VOLUMES,
C     DEVIENT UN SEUL VOLUME MONO-MATERIAU
C     IL SUFFIT DE DETRUIRE les tableaux a_volume__materiaux et union
C     ATTENTION: Cette OPERATION N'EST PAS REVERSIBLE!
C     ---------------------------------------------------------------
      CALL VOEX54( NUARBR, NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )


C     MISE A JOUR de la LONGUEUR des ARETES des TETRAEDRES SELON
C     SOIT la LONGUEUR>0 des ARETES par DEFAUT (code 0 du Menu DEBUT)
C     SOIT la FONCTION UTILISATEUR TAILLE_IDEALE(x,y,z) DONNEE ou 
C               the USER''s Function EDGE_LENGTH(x,y,z) GIVEN
C     ---------------------------------------------------------------
C     DECLARATION DU TABLEAU LADEFi=RADEFI de DEFINITION du VOLUME
C     POUR L'APPEL de VOEX26
      MOLADEFI = WYZDIM + 4 * NBPTIM
      CALL TNMCDC( 'REEL', MOLADEFI, MNLADEFI )

C     variable NUVOIN no du volume de longueur d'arete a ameliorer
      MCN( MNLADEFI + WUVOIN ) = NUARBR

C     variable NBPTIM nombre points internes imposes par l''utilisateur
      MCN( MNLADEFI + WBPTIM ) = NBPTIM

      IF( NBPTIM .GT. 0 ) THEN
C        tableau XYZDIM(1..4,1..NBPTIM) XYZ et distance souhaitee aux voisins
C        RADEFI(WYZDIM:WYZDIM+NBPTIM-1)
         CALL TRTATA( XYZDPTIM, MCN( MNLADEFI + WYZDIM ), MOLADEFI )
      ENDIF

C     variable NTYTRV 'transformation (I pour IDENTITE)' ^~>TRANSFO ;
      MCN( MNLADEFI + WTYTRV ) = 1

C     variable NUTYVO 'numero du type du volume' entier 
      MCN( MNLADEFI + WUTYVO ) = 26

C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNLADEFI) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNXYZS + MOTVAR(6) ) = NONMTD( '~>VOLUME>>DEFINITION' )

      CALL VOEX26( NUARBR, MCN( MNLADEFI), RMCN( MNLADEFI),
     %             NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )


C     FIN: RECUPERATION DE LA MEMOIRE DES TABLEAUX DECLARES DANS MCN
 9999 IF( MNNUZCPI .GT. 0 ) CALL TNMCDS( 'ENTIER', NTZCPI, MNNUZCPI )
      IF( MNLADEFI .GT. 0 ) CALL TNMCDS( 'REEL', MOLADEFI, MNLADEFI )

      RETURN
      END
