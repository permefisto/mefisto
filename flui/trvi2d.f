      SUBROUTINE TRVI2D( KNOMOB, NBNOEU, XYZNOE, MNNPEF,
     %                   NCAS0,  NCAS1,  vitx, vity,
     %                   VITMIN, VITMOY, VITMAX,
     %                   TIMES,  CMVITE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LA VITESSE DES NOEUDS D'UN FLUIDE 2D AVEC DES FLECHES
C -----

C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET FLUIDE
C NBNOEU : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
C XYZNOE : 3 COORDONNEES DES NOEUDS VITESSE DU FLUIDE
C MNNPEF : ADRESSE DE L'ADRESSE DU TABLEAU NPEF"TYPE EF (TRIA 2P1D ou 2P2C)

C NCAS0:NCAS1 : LES CAS TRAITES
C VITX   : VITESSE EN X AUX NOEUDS des CAS NCAS0 a NCAS1
C VITY   : VITESSE EN Y AUX NOEUDS des CAS NCAS0 a NCAS1
C VITMIN : NORME MINIMALE DE LA VITESSE
C VITMOY : NORME MOYENNE  DE LA VITESSE
C VITMAX : NORME MAXIMALE DE LA VITESSE

C TIMES  : LES TEMPS DES NBVPFILE VECTEURS VITESSE
C CMVITE : NOMBRE DE CM par UNITE de MODULE de la VITESSE

C SORTIES:
C --------
C VITMIN : NORME MINIMALE DE LA VITESSE
C VITMOY : NORME MOYENNE  DE LA VITESSE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    Novembre 2000
C MODIFS : ALAIN PERRONNET Saint PIERRE du PERRAY              Mars 2021
C MODIFS : ALAIN PERRONNET Saint PIERRE du PERRAY           Fevrier 2022
C MODIFS : ALAIN PERRONNET Saint PIERRE du PERRAY              Mars 2023
C23456---------------------------------------------------------------012
      PARAMETER        (LIGCON=0, LIGTIR=1 )
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ponoel.inc"
      include"./incl/a___npef.inc"
      include"./incl/ctemps.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)

      CHARACTER*(*)     KNOMOB
      CHARACTER*24      NOMFGIF
      INTEGER           NONOEF(6)
      REAL              XYZNOE(3,NBNOEU), TIMES(NCAS0:NCAS1)

      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ), dimension(NCAS0:NCAS1) :: vitx, vity
      DOUBLE PRECISION  VITNORM, VITMIN, VITMOY, VITMAX,
     %                           VITMI,  VITMO,  VITMA

      INTRINSIC         REAL, INT, SQRT

C     LA PALETTE 11: ARC EN CIEL
      CALL PALCDE(11)

C     LE NOMBRE-1 DE COULEURS DISPONIBLES
      NBCOUL = NDCOUL - N1COUL

C     EPAISSEUR DU TRAIT DES FLECHES
      NEPFLE = 0

C     NOM DU FICHIER VIDEO  SELON MODECO // 'Arow'
      CALL VIDEONM( 13, 'Arow', NOMFGIF )

C     OPTIONS DE LA VISEE POUR VOIR TOUTES LES NPAFLE FLECHES
C     =======================================================
      IF( NPAFLE .LE. 0 ) NPAFLE = 1

C     INITIALISATION DE TRANSLATION ZOOM
      IF( LORBITE .NE. 0 ) THEN
         CALL ZOOM2D0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 9000
      ENDIF

      VITMIN = 0D0
      VITMOY = 0D0

C     ======================================================================
C     BOUCLE SUR LES CAS (DIFFERENTS TEMPS) DES VITESSES A TRACER
C     ======================================================================
 20   DO 100 NCAS = NCAS0, NCAS1

         CALL EFFACEMEMPX

C        TEMPS DE CALCUL DU VECTEUR NCAS
         TEMPS = TIMES( NCAS )

C        TRACE DES AXES 2D
C        -----------------
         CALL TRAXE2

C        TRACE DES ARETES DES EF
C        -----------------------
         CALL XVEPAISSEUR( 1 )
         CALL XVTYPETRAIT( NTLAFR )

C        MNELE : ADRESSE DU TABLEAU NPEF"TYPE EF (TRIA 2P1D ou 2P2C)
         MNELE  = MCN( MNNPEF )

C        NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELFI = MCN( MNELE + WBELEM )

C        LE NUMERO DU TYPE DE L'ELEMENT FINI (TRIA 2P1D ou 2P2C)
         NUTYEL = MCN( MNELE + WUTYEL )

C        NUTYEF NO DU TYPE DU TRIANGLE ELEMENT FINI FLUIDE
C        NUTYEF = 1 SI TYPE BREZZI-FORTIN avec TRIA 2P1D
C               = 2 SI TYPE TAYLOR-HOOD   avec TRIA 2P2C

         DO NUELEM=1,NBELFI
C           LE NUMERO DES NOEUDS DE L'EF NUELEM
            CALL EFNOEU( MNELE, NUELEM, NBNOEF, NONOEF )
            K = 3
            DO I=1,3
               NS1 = NONOEF(K)
               NS2 = NONOEF(I)
               CALL TRAIT2D( NCOUAF, XYZNOE(1,NS1) , XYZNOE(2,NS1),
     %                               XYZNOE(1,NS2) , XYZNOE(2,NS2) )
               K = I
            ENDDO
         ENDDO

C        ----------------------------------------
C        TRACE DES VITESSES SOUS FORME DE FLECHES
C        ----------------------------------------
         CALL XVEPAISSEUR( NEPFLE )
         CALL XVTYPETRAIT( LIGCON )

         VITMO = 0D0
         VITMI = 1D100
         VITMA =-1D100

         NBF = 0
         DO NS1=1,NBNOEU

            IF( MOD(NS1,NPAFLE) .EQ. 0 ) THEN

C              LES 2 COORDONNEES DU NOEUD
               X = XYZNOE(1,NS1)
               Y = XYZNOE(2,NS1)

C              NORME DE LA VITESSE
               VITNORM = SQRT( vitx(NCAS)%dptab(NS1) **2
     %                       + vity(NCAS)%dptab(NS1) **2 )

               VITMO = VITMO + VITNORM
               IF( VITNORM .LT. VITMI ) VITMI = VITNORM
               IF( VITNORM .GT. VITMA ) VITMA = VITNORM

C              LONGUEUR CM DE LA FLECHE NBF DE LA VITESSE
               NBF = NBF + 1
               COXF  = REAL( vitx(NCAS)%dptab(NS1) * CMVITE )
               COYF  = REAL( vity(NCAS)%dptab(NS1) * CMVITE )
               CONCM = REAL( VITNORM * CMVITE )

C              COULEUR DU BOIS DE LA FLECHE
               NOCOUL = INT( VITNORM / VITMAX * NBCOUL + N1COUL )

C              LE TRACE DE LA FLECHE DE LA VITESSE
               CALL T2FLEC( NOCOUL, X, Y, CONCM, COXF, COYF )

            ENDIF

         ENDDO

C        VITESSE MOYENNE DU NCAS
         VITMO  = VITMO / NBF

         VITMOY = VITMOY + VITMO
         IF( VITMI .LT. VITMIN ) VITMIN = VITMI

C        EFFACEMENT DE LA LEGENDE SUR POSTSCRIPT
         IF ( LASOPS .NE. 0 ) THEN
            IF ( LASOPS .EQ. 1 ) THEN
               LASOPS = -11
            ELSE
               IF ( LASOPS .EQ. 2 ) THEN
                  LASOPS = -12
               ELSE
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'TRVI2D: MAUVAISE VALEUR de LASOPS'
                     KERR(2) = 'ARRET du TRACE POSTSCRIPT'
                  ELSE
                     KERR(1) = 'TRVI2D: BAD VALUE of LASOPS'
                     KERR(2) = 'STOP of the POSTSCRIPT DRAWING'
                  ENDIF
                  CALL LEREUR
                  LASOPS = 0
               ENDIF
            ENDIF
            CALL XVPOSTSCRIPT(LASOPS)
            LASOPS = - LASOPS
            CALL XVPOSTSCRIPT(LASOPS)
         ENDIF

C        LE TRACE DE LA LEGENDE DES VITESSES et DU TITRE FINAL
C        =====================================================
C        RETOUR AUX PARAMETRES INITIAUX
         CALL XVEPAISSEUR( 1 )
         CALL XVTYPETRAIT( LIGCON )

C        CONSTRUCTION DU TITRE ET TRACE
         CALL LEGVIT( KNOMOB, NCAS, VITMO, VITMI, VITMA, CMVITE )

C        MISE SUR FICHIER NomfgifBoImage.xwd puis NomfgifNoImage.jpg
C        DE LA PIXMAP de la FENETRE X11 ACTUELLE
         CALL VIDEO1( NOMFGIF, NCAS )
C
C        ATTENDRE POUR LIRE LE TRACE
         CALL ATTENDSEC( TEMP2TRAC )
C
C     FIN DE LA BOUCLE SUR LES CAS
 100  CONTINUE


C     CONSTRUCTION FINALE DU FICHIER.gif
      CALL VIDEOFIN( NOMFGIF )

      VITMOY = VITMOY / ( NCAS1 - NCAS0 + 1 )
C
C     RETOUR POUR UNE NOUVELLE VISEE
C     ------------------------------
      IF( LORBITE .NE. 0 ) THEN
         IF( NCAS0 .EQ. NCAS1 ) THEN
C           ZOOM  BOUTON ENFONCE et DEPLACE
            CALL ZOOM2D1( NOTYEV )
         ELSE
C           ZOOM  BOUTON ENFONCE et DEPLACE et RELACHE
            CALL ZOOM2D3( NOTYEV )
         ENDIF
         IF( NOTYEV .EQ. 0 ) GOTO 9000
         GOTO 20
      ELSE
         CALL CLICSO
      ENDIF
C
C     RETOUR AU TRACE NORMAL POUR POSTSCRIPT
 9000 IF( LASOPS.NE.0 ) THEN
        LASOPS = LASOPS -10
        CALL XVPOSTSCRIPT(LASOPS)
      ENDIF
C
      RETURN
      END
