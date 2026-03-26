      SUBROUTINE SUEX25( NUARBR, LADEFI, RADEFI,
     %                   NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER UNE TRIANGULATION D'UN ARBRE ET DE SES RACINES
C -----    A PARTIR DE TRONCS DE CONE DEFINIS PAR DES CERCLES
C          APPARTENANT A DES PLANS Z=Cte

C ENTREES:
C --------
C NUARBR : NUMERO DE LA SURFACE ARBRE DANS LE LEXIQUE DES SURFACES
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA SURFACE ARBRE
C RADEFI : TABLEAU REEL   DE DEFINITION DE LA SURFACE ARBRE
C          CF '~/td/d/a_surface__definition'

C SORTIES:
C --------
C NTNSEF : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES
C MNNSEF : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES
C            CF '~/td/d/a___nsef'
C NTXYZS : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNXYZS : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C            CF '~/td/d/a___xyzsommet'
C IERR   : = 0 SI PAS D'ERREUR
C          <>0 SI ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  Saint PIERRE du PERRAY             Juin 2019
C23456...............................................................012
      PARAMETER      (MXZCER=256)
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/darete.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL           RMCN(1)
      EQUIVALENCE    (MCN(1),RMCN(1))
      INTEGER        LADEFI(0:*)
      REAL           RADEFI(0:*)

C     DECLARATION MAXIMALE DES TABLEAUX NECESSAIRES
      REAL           RAYC(MXZCER),   XYZCZC(3,MXZCER)
      INTEGER        NUPTZC(MXZCER), NBARLC(MXZCER)
      INTEGER        NULGZC(MXZCER), NCTRIZ(MXZCER), NUSFTRPI(MXZCER)
      CHARACTER*24   KNM
      DOUBLE PRECISION DEUXPI

      DEUXPI = ATAN( 1D0 ) * 8D0

C     RECUPERATION DES DONNEES DU TMS a_surface__definition
C     -----------------------------------------------------
      IERR     = 0
      MNNUZCPI = 0

C     SURFACE d''UN ARBRE et SES RACINES a PARTIR de CERCLES: TYPE 25
ccc   NUTYSU = LADEFI( WUTYSU )

C     NBPIEC 'nombre de PIECES ou TRONCONS de l''ARBRE'
      NBPIEC = LADEFI( WBPIEC )

C     NTZCPI 'nombre total de Z-CERCLES decrivant les PIECES'
      NTZCPI = LADEFI( WTZCPI )

C     NBZCER 'nombre de lignes Z-CERCLES'
      NBZCER = LADEFI( WBLGZC )

      IF( NBZCER .LT. 4 ) THEN
         PRINT *,'suex25: ARBRE avec NBZCER=',NBZCER,
     %           '<4 un NOMBRE de CERCLES INSUFFISANT'
         IERR = 8
         GOTO 9999
      ENDIF

      IF( NBZCER .GT. MXZCER ) THEN
         PRINT *,'suex25: AUGMENTER MXZCER=',MXZCER,
     %           ' AU DELA de NBZCER=',NBZCER
         IERR = 9
         GOTO 9999
      ENDIF

C     LES TABLEAUX COPIES DE LADEFI SONT CREES POUR POUVOIR ETRE MODIFIES
c     tableau  NUPTZC(1..NBZCER) 'nom des centres des cercles' ^~>POINT ;
      MN = WUPTZC -1
      DO NC = 1, NBZCER

C        NUMERO DU POINT-CENTRE DANS LE LEXIQUE DES POINTS
         NUPTZC( NC ) = LADEFI( MN + NC )

C        RECUPERATION DES 3 COORDONNEES DU POINT CENTRE DU CERCLE NC
         CALL LXNLOU( NTPOIN, NUPTZC( NC ), NTLXZC, MNP )
         IF( NTLXZC .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'POINT INCONNU'
            CALL LEREUR
            IERR = 6
            GOTO 9999
         ENDIF

C        RESTAURATION DU TMS XYZSOMMET DU POINT NUPTZC( NC )
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
 
c     tableau  RAYZC(1..NBZCER)  'rayon des cercles' reel ;
      MN = MN + NBZCER
      DO NC = 1, NBZCER

C        RAYON DU CERCLE NC
         RAYC( NC ) = RADEFI( MN + NC )

C        CALCUL DU NOMBRE D'ARETES DU CERCLE NC EN FONCTION DE
C        LA LONGUEUR DE L'ARETE PAR DEFAUT CONNUE DARETE
         IF( RAYC( NC ) .EQ. 0D0 ) THEN
C           CERCLE REDUIT A SON POINT CENTRE: PAS D'ARETE
            NB = 0
         ELSE
            NB = NINT( DEUXPI * RAYC( NC ) / DARETE )
            NB = MAX( 5, NB )
         ENDIF
         NBARLC( NC ) = NB

      ENDDO

C     tableau NBZCPI(1..NBPIEC) 'nombre de Z-CERCLES de chaque PIECE'
      MNNBZCPI = MN + NBZCER + 1

C     tableau PTZCPI(1..NTZCPI) 'nom du CENTRE des Z-CERCLES de chaque PIECE'
      MNPTZCPI = MNNBZCPI + NBPIEC

C     tableau numero de 1 a NBZCER de chaque Z-CERCLE des PIECES de l'ARBRE
      CALL TNMCDC( 'ENTIER', NTZCPI, MNNUZCPI )
      MN = MNNUZCPI
      DO 10 K = 1, NTZCPI

C        NUMERO DU K-EME POINT CENTRE DU Z-CERCLE DANS LE LEXIQUE POINT
         NUZC = LADEFI( MNPTZCPI-1+K )

         DO NC = 1, NBZCER
C           NUMERO DU POINT CENTRE DANS LE LEXIQUE POINT
            IF( NUZC .EQ. NUPTZC( NC ) ) THEN
C               NUZC EST LE CENTRE DU Z-CERCLE NC
               MCN( MN ) = NC
               MN = MN + 1
               GOTO 10
            ENDIF
         ENDDO

C        CE POINT NUZC N'EST PAS UN CENTRE DE Z-CERCLES
         CALL NMOBNU( 'POINT', NUZC, KNM )
         PRINT*,'suex25: LE POINT ',KNM,
     %          ' N''EST PAS UN POINT CENTRE DES Z-CERCLES'
         IERR = 10

 10   ENDDO

      IF( IERR .NE. 0 ) GOTO 9999

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
10000 FORMAT('suex25: TRIANGULATION ',A,
     %       ' FONCTION TAILLE_IDEALE(X,Y,Z) (ou EDGE_LENGTH(X,Y,Z))' )
20000 FORMAT('suex25: TRIANGULATION ',A,
     %       ' FUNCTION EDGE_LENGTH(X,Y,Z) (or TAILLE_IDEALE(X,Y,Z))' )
10001 FORMAT(' LONGUEUR de l''ARETE SOUHAITEE =',G14.6)
20001 FORMAT(' WISHED EDGE LENGTH =',G14.6)


C     CONSTRUCTION DE LA SURFACE DE L'ARBRE et de SES RACINES
C     -------------------------------------------------------
      CALL ARBRESF( NBZCER, NUPTZC, XYZCZC, RAYC, NBARLC,
     %              NBPIEC, NTZCPI, LADEFI(MNNBZCPI), MCN(MNNUZCPI),
     %              NULGZC, NCTRIZ, NUSFTRPI, IERR )
      IF( IERR .NE. 0 ) GOTO 9999


C     UNION DES NBPIEC SURFACES TRONCS DE CONE CREEES ENTRE CERCLES
C     -------------------------------------------------------------
      IF( NBPIEC .GE. 2 ) THEN

C        CONSTRUCTION DU MAILLAGE DE LA SURFACE UNION des TRONCS DE CONE
         CALL UNPLSV( 51,  3, NUARBR, NBPIEC, NUSFTRPI,
     %                NTNSEF, MNNSEF, NTXYZS, MNXYZS,
     %                NTUNIO, MNUNIO, IERR )

C        DESTRUCTION DE LA SURFACE DES NBPIEC PIECES
         DO N = 1, NBPIEC

C           NUMERO DE LA SURFACE
            NULX = NUSFTRPI( N )

            IF( NULX .GT. 0 ) THEN
C              NOM DE LA SURFACE DANS LE LEXIQUE DES SURFACES
               CALL NMOBNU( 'SURFACE', NULX, KNM )
               CALL LXLXOU( NTSURF, KNM, NTSF, MNSF )
               IF( MNSF .GT. 0 ) CALL LXTSDS( NTSURF, KNM )
            ENDIF

         ENDDO

C        SUR LA TRIANGULATION DE LA SURFACE TOTALE DE L'ARBRE
C        BARYCENTRAGE DES SOMMETS HORMIS LES SOMMETS DES CONES
C        C-A-D LES CENTRES DES Z-CERCLES DE RAYON NUL
C        -----------------------------------------------------
         NBSOM = MCN( MNXYZS + WNBSOM )
         CALL TNMCDC( 'ENTIER', NBSOM+NBZCER, MNSTCO )
         CALL AZEROI( NBSOM+NBZCER, MCN(MNSTCO) )

C        LISTAGE DES CENTRES DES Z-CERCLES DE NOMBRE D'ARETES NUL
         MNSTCC = MNSTCO + NBSOM
         NBSTCO = 0
         DO NC = 1, NBZCER
            IF( NBARLC( NC ) .LE. 0 ) THEN
               MCN( MNSTCC + NBSTCO ) = NC
               NBSTCO = NBSTCO + 1
            ENDIF
         ENDDO

         DO 20 N = 1, NBSTCO

C           NUMERO DE 1 A NBZCER DU POINT CENTRE DES Z-CERCLES
C           DE NOMBRE D'ARETES NUL
            NSTC = MCN( MNSTCC -1 + N )

            MNS = MNXYZS + WYZSOM
            DO NS = 1, NBSOM
C              IDENTIFICATION DES 2 SOMMETS?
               CALL XYZIDE( XYZCZC(1,NSTC), RMCN(MNS), IDENTQ )
               IF( IDENTQ .NE. 0 ) THEN
C                 OUI: SOMMETS IDENTIFIES
                  MCN( MNSTCO - 1 + NS ) = NSTC
                  GOTO 20
               ENDIF
               MNS = MNS + 3
            ENDDO

 20      ENDDO

C        BARYCENTRAGES DES SOMMETS DE LA TRIANGULATION
C        EXCEPTES LES SOMMETS DE CONES
         DO NBITERBA = 1, 5
            CALL BASTMAEF( NBSOM, RMCN( MNXYZS+WYZSOM), MCN( MNSTCO ),
     %                     4, MCN(MNNSEF+WBEFOB), MCN(MNNSEF+WUSOEF) )
         ENDDO

         CALL TNMCDS( 'ENTIER', NBSOM+NBZCER, MNSTCO )

      ELSE

         PRINT *,'suex25: ARBRE avec NBPIEC=',NBPIEC,
     %           '<2 TRIANGULATIONS. un NOMBRE INSUFFISANT'
         IERR = 12
         GOTO 9999

      ENDIF

 9999 IF( MNNUZCPI.GT.0 ) CALL TNMCDS( 'ENTIER', NTZCPI, MNNUZCPI )

      RETURN
      END
