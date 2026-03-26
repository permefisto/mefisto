      SUBROUTINE ARBRESF( NBZCER, NUPCZC, XYZCZC, RAYZC, NBARLC,
     %                    NBPIEC, NTZCPI, NBZCPI, NUZCPI,
     %                    NULGZC, NCTRIZ, NUSFTRPI, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER UNE TRIANGULATION DE LA SURFACE D'UN ARBRE
C -----    DEFINI PAR NBPIEC PIECES (TRONCONS) DEFINIES PAR LA LISTE
C          DE SES Z-CERCLES
C          ATTENTION: UN Z-CERCLE PEUT ETRE REDUIT A UN POINT SON CENTRE
C                     DANS CE CAS LE NOMBRE DE SOMMETS DOIT ETRE EGAL A 0
C          LE TYPE D'UNE PIECE PEUT ETRE
C          1: CONE SOMMET BAS  (defini par 1 Pt=ZC + 1 ZC)
C          2: CONE SOMMET HAUT (defini par 1 ZC + 1 Pt=ZC)
C          3: TRONC DE CONE    (defini par 2 ZC)
C          4: n Z-CERCLES RELIES a UN Z-CERCLE  (defini par n + 1 ZC)
C          5: UN Z-CERCLE RELIE  a n  Z-CERCLES (defini par 1 + n ZC)

C ENTREES:
C --------
C NBZCER : NOMBRE DE Z-CERCLES DE L'ARBRE
C NUPCZC : NUMERO LX POINTS DES NBZCER CENTRES DES CERCLES
C XYZCZC : XYZ DU CENTRE DES NBZCER Z-CERCLES
C RAYZC  : RAYON DES NBZCER Z-CERCLES
C NBARLC : >2 NOMBRE D'ARETES DE CHACUNE DES NBZCER CERCLES
C          =0 INDIQUE UN CERCLE REDUIT A SON POINT CENTRE

C NBPIEC : nombre de PIECES ou TRONCONS de l'ARBRE
C NTZCPI : nombre total de Z-CERCLES nommees pou DEFINIR
C          les NBPIEC PIECES de l'ARBRE
C NBZCPI : nombre de Z-CERCLES de chaque PIECE de l'ARBRE
C NUZCPI : NUMERO 1 a NBZCER du Z-CERCLE de chaque PIECE

C AUXILIAIRES:
C ------------
C KNMLGZC: NOM DES NBZCER CERCLES LIGNES ARETISEES
C  NULGZC: NUMERO DANS LE LEXIQUE LIGNES DES NBC CERCLES
C NCTRIZ : NUMERO DE 1 A NBZCER DES Z-CERCLES D'UNE PIECE DE L'ARBRE
C          A TRIER SELON Z CROISSANT DE LEUR CERCLE

C SORTIES:
C --------
C NUSFTRPI: NUMERO DANS LE LEXIQUE SURFACES DES NBPIEC SURFACES TRIANGULEES
C IERR    : =0 SI PAS D'ERREUR DETECTEE, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  Saint PIERRE du PERRAY             Juin 2019
C23456...............................................................012
      PARAMETER      (MXZCER=256, MXPTIM=1)
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_point__definition.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"

      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))

      REAL             RAYZC(NBZCER),  XYZCZC(3,NBZCER),
     %                 XYZDPTIM(4,MXPTIM)
      INTEGER          NUPCZC(NBZCER), NBARLC(NBZCER)
      INTEGER          NBZCPI(NBPIEC), NUZCPI(NTZCPI)
      INTEGER          NULGZC(NBZCER), NCTRIZ(NBZCER),
     %                 NUSFTRPI(NBPIEC)

      CHARACTER*8      KNC
      CHARACTER*24     KNMLGZC, KNM

      IERR    = 0
      MNXYZCP = 0
      MNFRST  = 0
      NBPTIM  = 0

      print*
      PRINT*,'arbresf: TRIANGULATION de la SURFACE d''un ARBRE a partir 
     %de',NBPIEC,' PIECES et',NBZCER,' Z-CERCLES'

C     NOMBRE DE SURFACES JOIGNANT AU MOINS 2 CERCLES
      NBPISF = 0

C     CONSTRUCTION DES XYZ DES SOMMETS DES NBZCER Z-CERCLES DECLARES
C     =============================================================
      CALL TNMCDC( 'ENTIER', NBZCER, MNXYZC0 )
      DO NC = 1, NBZCER

C        CERCLE NC A TRAITER ou REDUIT AU POINT CENTRE?
         IF( NBARLC( NC ) .EQ. 0 ) THEN

            CALL NMOBNU( 'POINT', NUPCZC(NC), KNM )
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'arbresf: ARETISATION du Z-CERCLE',NC,
     %                ' REDUIT au POINT CENTRE: ',KNM,
     %                ' NON FAITE'
            ELSE
               PRINT*,'arbresf: EDGES of Z-CIRCLE',NC,
     %                ' REDUCED to one POINT CENTRE: ',KNM
            ENDIF

C           NUMERO DE LA LIGNE DANS LE LX DE LA LIGNE NC
C           DEVIENT LE NUMERO DU POINT CENTRE SOMMET D'UN CONE
            NULGZC( NC ) = -NC

         ELSE IF( NBARLC( NC ) .GE. 3 ) THEN

C           AU MOINS 3 SOMMETS DEMANDES SUR LE CERCLE NC
            IF( RAYZC(NC) .LE. 0. ) THEN
               PRINT*,'arbresf: Z-CERCLE',NC,' avec', NBARLC(NC),
     %                ' ARETES et un RAYON',RAYZC(NC),' INCOHERENT'
               IERR = 10
               GOTO 9990
            ENDIF

C           KNMLGZC NOM DE LA LIGNE MOMENTANEE DU CERCLE NC
            WRITE( KNC(1:8), '(I8)' ) NC
C              RETRAIT DES CARACTERES BLANCS
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
               PRINT*,'arbresf: ERREUR dans le MAILLAGE du CERCLE',NC
     %               ,KNMLGZC
               GOTO 9990
            ENDIF

C           SUPPRESSION DES TANGENTES AU CERCLE
            MCN( MNXYZC  + WNBTGS ) = 0
            MCN( MNNSEFC + WBTGEF ) = 0
            MCN( MNNSEFC + WBEFAP ) = 0
            MCN( MNNSEFC + WBEFTG ) = 0

            IF( LANGAG .EQ. 0 ) THEN
               print*,'arbresf: ARETISATION du Z-CERCLE',NC,' de NOM: ',
     %                 KNMLGZC,' avec',MCN(MNNSEFC+WBEFOB),' ARETES',
     %                 MCN(MNXYZC+WNBSOM),' SOMMETS est la LIGNE',
     %                 NULGZC(NC)
            ELSE
               PRINT*,'arbresf: EDGES of Z-CIRCLE',NC, ' of NAME: ',
     %                 KNMLGZC,' with',MCN(MNNSEFC+WBEFOB),' EDGES',
     %                 MCN(MNXYZC+WNBSOM),' VERTICES is the LINE',
     %                 NULGZC(NC)
            ENDIF

         ELSE IF( NBARLC( NC ) .GT. 0 ) THEN

            PRINT *,' arbresf: CERCLE', NC,' avec',NBARLC(NC),
     %              ' SOMMETS DEMANDES au LIEU de =0 ou >2'
            IERR = 1
            GOTO 9990

ccc         ELSE

cccC           CERCLE REDUIT A UN CENTRE=POINT=SOMMET

         ENDIF

      ENDDO
      PRINT*


C     CREATION DE LA TRIANGULATION DE LA SURFACE DES NBPIEC PIECES
C     ============================================================
C     NOMBRE DE PIECES TRAITEES
      NBPISF = 0

C     NO DU DERNIER Z-CERCLE AVANT CEUX DE LA PIECE A TRAITER
      NDPIEC = 0

      DO NOPIEC = 1, NBPIEC

C        NOMBRE DE Z-CERCLES DE LA PIECE NOPIEC
         NBZCP = NBZCPI( NOPIEC )

         IF( NBZCP .LT. 2 ) THEN
            PRINT*,'arbresf:',NBZCP,
     %             '<2 NOMBRE INCORRECT DE Z-CERCLES de la PIECE',NOPIEC
            IERR = 2
            GOTO 9990
         ENDIF

C        RECHERCHE DU TYPE NOTYPE DE LA PIECE
C        1: CONE SOMMET BAS
C        2: CONE SOMMET HAUT
C        3: TRONC DE CONE
C        4: n Z-CERCLES RELIES a UN Z-CERCLE
C        5: UN Z-CERCLE RELIE  a n  Z-CERCLES
C        ET TRIANGULATION DE SA SURFACE LIMITEE PAR NBZCP Z-CERCLES

C        PAS DE CREATION DE POINTS INTERNES A IMPOSER
         NOPTIM = 0
         CALL ARBRESF1P( NBZCER, NUPCZC, XYZCZC, RAYZC, NBARLC,
     %                   NULGZC, MNXYZC0,
     %                   NOPIEC, NBZCP,  NUZCPI(NDPIEC+1), NCTRIZ,
     %                   NOPTIM, MXPTIM, NBPTIM, XYZDPTIM,
     %                   NUSFTRPI(NOPIEC), NOTYPE, IERR )
         IF( IERR .NE. 0 ) GOTO 9990

C        PASSAGE A LA PIECE SUIVANTE
         NDPIEC = NDPIEC + NBZCP 

      ENDDO

C     SUPPRESSION DES ADRESSES DES TMS XYZSOMMET DES CERCLES
C     ======================================================
 9990 CALL TNMCDS( 'ENTIER', NBZCER, MNXYZC0 )

C     DESTRUCTION DES LIGNES Z-CERCLES
      DO NC = 1, NBZCER

C        NOM DE LA LIGNE DU CERCLE NC
         NULX = NULGZC( NC )
         IF( NULX .GT. 0 ) THEN
C           NOM DE LA LIGNE CERCLE DANS LE LX LIGNES
            CALL NMOBNU( 'LIGNE', NULX, KNM )
            CALL LXLXOU( NTLIGN, KNM, NT1LZC, MN1LZC )
            IF( MN1LZC .GT. 0 ) CALL LXTSDS( NTLIGN, KNM )
         ENDIF

      ENDDO

      RETURN
      END

