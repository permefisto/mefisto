      SUBROUTINE TTTSUPA2D( KNOMOB, NTLXOB, MNNPEF, MNXYZP,
     %                      NBNOVI, NCAS0,  NCAS1,  NBPART, MNPART,
     %                      vitx,   vity,   VITMAX, TIMES )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACER 2D DU SUIVI DE PARTICULES DANS UN FLUIDE 2D MAILLE EN
C -----  TRIANGLES DURANT L'INTERVALLE DE TEMPS CONNAISSANT LA VITESSE
C        (TAYLOR-HOOD et BREZZI-FORTIN SONT DES EF NON ISOPARAMETRIQUES)

C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET
C NTLXOB : NO TMS DU LEXIQUE DE L'OBJET KNOMOB
C MNNPEF : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE
C MNXYZP : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
C NCAS0  : NUMERO DU PREMIER TEMPS PARMI LES NBVPFILE TEMPS ET VECTEUR
C NCAS1  : NUMERO DU DERNIER TEMPS PARMI LES NBVPFILE TEMPS >NCAS0
C NBPART : NOMBRE DE PARTICULES A SUIVRE
C MNPART : RMCN ADRESSE DES NBPART XYZ + VITESSE XYZ INITIALE
C          + RAYON + TEMPS DE DEPART de chaque BOULE-PARTICULE
C vitx   : TABLEAU(ncas0:ncas1) de pointeurs sur le tableau de la vitesseX(NBNOVI)
C vity   : TABLEAU(ncas0:ncas1) de pointeurs sur le tableau de la vitesseY(NBNOVI)
C VITMAX : MODULE MAXIMAL DE LA VITESSE EN UN NOEUD
C TIMES  : INSTANT DU CALCUL DES NBVPFILE VECTEURS VITESSE+PRESSION
C          STOCKE SUR FICHIERS DANS LE REPERTOIRE DU PROJET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  OCTOBRE 2010
C23456---------------------------------------------------------------012
      PARAMETER      ( NBPAST=1000, MXSEGM=100*NBPAST, LIGCON=0 )
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___arete.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ctemps.inc"
      include"./incl/xyzext.inc"
      include"./incl/nctyef.inc"

      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     KNOMOB
      CHARACTER*24      NOMFGIF

      REAL              TIMES(NCAS0:NCAS1)

      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(NCAS0:NCAS1) :: vitx, vity
      DOUBLE PRECISION  VITMAX
      DOUBLE PRECISION  PROSCD, CB1,CB2,CB3, FBASE(6), VITEF(6)
      DOUBLE PRECISION  XP0, YP0, XP1, YP1, VXYP0(2), V, VPMAX, VTPMAX
      DOUBLE PRECISION  DELTAT, PARCOUR
      INTRINSIC         SQRT, REAL, INT

      IF( NCAS0 .GE. NCAS1 .OR. NBPART .LE. 0 ) RETURN
      NBSOMT = 0
      VMIN   = 0
      VMAX   = REAL( VITMAX )
      VPMAX  = 0D0
C     VITESSE MAXIMALE des PARTICULES: VALEUR par DEFAUT
C     MISE A JOUR DES LE SECOND CALCUL
      VTPMAX = 9.81D0 * 2D0

C     NOM DU FICHIER VIDEO // 'path'
      CALL VIDEONM( 5, 'path', NOMFGIF )

C     NBCOOR = NOMBRE DE COORDONNEES DES POINTS 3D ou 6D
      NBCOOR = MCN(MNXYZP+WBCOOP )
      IF( NBCOOR .NE. 3 ) GOTO 9900

C     NBPOI   NOMBRE DE POINTS DU MAILLAGE DE L'OBJET
      NBPOI = MCN(MNXYZP+WNBPOI)

C     RECUPERATION DU TMS DES ARETES DE L'OBJET 2D
C     ============================================
      CALL LXTSOU( NTLXOB, 'ARETE', NTARET, MNARET )

      IF( NTARET .LE. 0 ) THEN
C        LE TABLEAU N'EXISTE PAS => IL EST CREE
C        CALCUL PAR HACHAGE DES ARETES DE L'OBJET A PARTIR DE TOPO+NPEF"...
C        NECESSAIRE POUR CONNAITRE LES EF ADJACENTS PAR ARETE
         CALL HAC2AF( KNOMOB, 3, NTARET, MNARET, IERR )
         IF( IERR .NE. 0 ) GOTO 9900
      ENDIF
      IF( NTARET .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'OBJET SANS TABLEAU DES ARETES'
         ELSE
            KERR(1) = 'OBJECT WITHOUT THE ARRAY OF EDGES'
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF

C     LE NOMBRE D'ENTIERS PAR ARETE
      MOARET = MCN( MNARET + WOARET )
C     LA MAJORATION DU NOMBRE D'ARETES
      MXARET = MCN( MNARET + WXARET )
C     LE NOMBRE D'ARETES FRONTALIERES
      NBARFB = MCN( MNARET + WBARFB )
C     LE NOMBRE D'ARETES INTERFACES
      NBARIN = MCN( MNARET + WBARIN )
C     LE NUMERO MINIMAL DE LIGNE DE L'OBJET
      NUMILF = MCN( MNARET + WUMILF )
C     LE NUMERO MAXIMAL DE LIGNE DE L'OBJET
      NUMXLF = MCN( MNARET + WUMXLF )
C     LE NUMERO DE LA PREMIERE ARETE FRONTALIERE NON SUR LIGNES DE L'OBJET
      L1ARFB = MCN( MNARET + W1ARFB )
C     LE NUMERO DE LA PREMIERE ARETE INTERFACE NON SUR LIGNES DE L'OBJET
      L1ARIN = MCN( MNARET + W1ARIN )
C     ADRESSE MCN DU 1-ER MOT DU TABLEAU LARETE
      MNLARE = MNARET + W1LGFR + NUMXLF - NUMILF + 1
C     LARETE : TABLEAU DES ARETES DU MAILLAGE
C     LARETE(1,I)= NO DU 1-ER  SOMMET DE L'ARETE (0 SI PAS D'ARETE)
C     LARETE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C     LARETE(3,I)= CHAINAGE HACHAGE SUR ARETE SUIVANTE
C     LARETE(4,I)= NUMERO DU 1-ER TRIANGLE CONTENANT CETTE ARETE
C                  0 SI PAS DE 1-ER  TRIANGLE
C     LARETE(5,I)= NUMERO DU 2-EME TRIANGLE CONTENANT CETTE ARETE
C                  0 SI PAS DE 2-EME TRIANGLE
C     LARETE(6,I)= NUMERO DANS LARETE DE L'ARETE SUIVANTE
C              SOIT DANS LE CHAINAGE D'UNE LIGNE J ENTRE NUMILF ET NUMX
C              SOIT DANS LE CHAINAGE DES ARETES FRONTALIERES
C              0 SI C'EST LA DERNIERE

C     L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF DE CE TYPE D'EF
      MNELE = MCN( MNNPEF )

C     LE NUMERO DU TYPE DES ELEMENTS FINIS
      NUTYEL = MCN( MNELE + WUTYEL )

C     LE NOMBRE DE TELS ELEMENTS FINIS
      NBELEM = MCN( MNELE + WBELEM )

C     LES CARACTERISTIQUES DE L'ELEMENT FINI
      CALL ELTYCA( NUTYEL )

C     NO D'INTERPOLATION DES COMPOSANTES DE LA VITESSE
C     NOMBRE DE NOEUDS VITESSE DE L'EF
      IF( NUTYEL .EQ. 13 ) THEN
C        TRIANGLE BREZZI-FORTIN
         NOINTE = 2003
         NBNOE  = 4
         NBNOED = 3
         NBSOMT = NBNOVI - NBELEM
      ELSE IF( NUTYEL .EQ. 15 ) THEN
C        TRIANGLE TAYLOR-HOOD
         NOINTE = 2002
         NBNOE  = 6
         NBNOED = NBNOE
      ELSE
C        EF NON TRAITE ICI
         IERR = 1
         GOTO 9900
      ENDIF

C     ECRAN EFFACE
      CALL EFFACEMEMPX

C     NOMBRE DE COULEURS DISPONIBLES
      NBCOUL = NDCOUL - N1COUL + 1

C     LA PALETTE 11: ARC EN CIEL
      CALL PALCDE(11)

C     COULEUR PAR DEFAUT DES ARETES DES FACES FRONTIERE
      NCOAFR = NCGRIC

C     OPTIONS DE LA VISEE POUR VOIR LES PARTICULES
C     ============================================
 20   CALL VISE2D( NMTCL )
      IF( NMTCL .LT. 0 ) GOTO 9900

C     INITIALISATION DE TRANSLATION ZOOM
      NBORBIT = 0
      IF( LORBITE .NE. 0 ) THEN
         CALL ZOOM2D0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 20
      ENDIF

C     TRACE DES AXES 2D
C     -----------------
 30   CALL TRAXE2

C     TRACE DE TOUTES LES ARETES DU TABLEAU LARETE
C     --------------------------------------------
      CALL TRARM2D( NCOUAF, MOARET, MXARET, MCN(MNLARE),
     %              RMCN(MNXYZP+WYZPOI) )

C     TRACE DES ARETES CHAINEES DANS LE TABLEAU LARETE
C     ------------------------------------------------
C     TRACE DES ARETES FRONTALIERES NON SUR LES LIGNES DE L'OBJET
      MNXYZ = MNXYZP + WYZPOI

      CALL TRARCH2D(NCOAFR, MOARET, L1ARFB, MCN(MNLARE), RMCN(MNXYZ))

C     TRACE DES ARETES INTERFACES NON SUR   LES LIGNES DE L'OBJET
      CALL TRARCH2D(NCOAFR, MOARET, L1ARIN, MCN(MNLARE), RMCN(MNXYZ))

C     TRACE DES ARETES DES LIGNES FRONTALIERES OU INTERFACES
C     DEFINIES DANS L'OBJET PAR L'UTILISATEUR
      MN1LGF = MNARET + W1LGFR - NUMILF
      DO I = NUMILF, NUMXLF
C        TETE DE CHAINAGE DES ARETES DE LA LIGNE I
         L1AR = MCN(MN1LGF+I)
         CALL TRARCH2D(NCOAFR, MOARET, L1AR, MCN(MNLARE),RMCN(MNXYZ))
      ENDDO

C     PAS DE TEMPS POUR INTEGRER LA VITESSE ET RECOUVRIR L'INTERVALLE DE TEMPS
C     ------------------------------------------------------------------------
      DELTAT = ( TIMES(NCAS1) - TIMES(NCAS0) ) / NBPAST
      PRINT *
      IF( LANGAG .EQ. 0 ) THEN
         PRINT *, 'Du TEMPS',TIMES(NCAS0),' au ',TIMES(NCAS1),
     %            ' AVEC',NBPAST,' SOUS-PAS DE TEMPS',DELTAT
      ELSE
         PRINT *, 'From TIME',TIMES(NCAS0),' to ',TIMES(NCAS1),
     %            ' with',NBPAST,' TIME SUB-STEPS',DELTAT
      ENDIF

C     =======================================================================
C     BOUCLE SUR LES PARTICULES K
C     =======================================================================

C     LES TRAITS DU PARCOURS DE LA PARTICULE SONT TRACES AVEC 3 EPAISSEURS
      CALL XVEPAISSEUR( 3 )

      DO 200 K = 1, NBPART

         IF( LANGAG .EQ. 0 ) THEN
            PRINT *, 'tttsupa2d: Debut du parcours de la particule',K
         ELSE
            PRINT *, 'tttsupa2d: the particle',K,' starts'
         ENDIF

C        INITIALISATIONS POUR LA PARTICULE K
         NBSEGM  = 0
         PARCOUR = 0D0

C        LES 2 COORDONNEES INITIALES DE LA PARTICULE K
         MNPA = MNPART + 8*K -8
         XP0  = RMCN(MNPA  )
         YP0  = RMCN(MNPA+1)

C        LES 2 COMPOSANTES DE LA VITESSE INITIALE DE LA PARTICULE K
         XV0PART = RMCN(MNPA+3)
         YV0PART = RMCN(MNPA+4)

C        TEMPS du DEPART DE LA PARTICULE K
         TDEPART = RMCN(MNPA+7)
         IF( TDEPART .GE.  TIMES(NCAS1) ) THEN
C           TEMPS DU DEPART AU DELA DES TEMPS EN MEMOIRE
C           RECUL AU 1ER TEMPS DE L'ETUDE
            TDEPART = TIMES( NCAS0 )
         ENDIF

C        TRACE DE LA PARTICULE K AU POINT DE DEPART ET AFFICHAGE
         XV0 = REAL( XP0 )
         YV0 = REAL( YP0 )
         CALL SYMBOLE2D( NCNOIR, XV0, YV0, 'D' )

C        RECHERCHE EXHAUSTIVE DU TRIANGLE NEF CONTENANT LE POINT XYP0
C        ------------------------------------------------------------
         CALL RETRIAXY0( XP0, YP0, MNNPEF, MNXYZP,
     %                   NEF, CB1, CB2, CB3, IERR )
         IF( IERR .GT. 0 ) GOTO 200

C        LA BOUCLE SUR LES PAS DE TEMPS POUR LA PARTICULE K
C        ==================================================
C        RECHERCHE DE L'INTERVALLE DE TEMPS CONTENANT TDEPART
         DO NCAS = NCAS0, NCAS1-1
            IF(TDEPART.GE.TIMES(NCAS).AND.TDEPART.LT.TIMES(NCAS+1)) THEN
               GOTO 5
            ENDIF
         ENDDO

C        TDEPART N'EST PAS DANS L'INTERVALLE : NCAS0 EST SON TEMPS DE DEPART
         TDEPART = TIMES( NCAS0 )
         NCAS = NCAS0

C        LE TEMPS INITIAL TDEPART EST ENTRE TIMES(NCAS) et TIMES(NCAS+1)
 5       TEMPS0 = TIMES( NCAS   )
         TEMPS1 = TIMES( NCAS+1 )
         DIFTEMP= TEMPS1 - TEMPS0
         TEMPS  = TDEPART

         DO 100 NUPAST = 1, NBPAST

C        TEMPS AU PAS DE TEMPS NUPAST
         TEMPS = REAL( TEMPS + DELTAT )
         IF( TEMPS .GT. TEMPS1 ) THEN

C           FIN DU PARCOURS PAR ATTEINTE DU TEMPS MAXIMUM?
            NCAS = NCAS + 1
            IF( NCAS .EQ. NCAS1 ) GOTO 110

C           NON
            TEMPS0 = TEMPS1
            TEMPS1 = TIMES( NCAS+1 )
C           LE TEMPS EST ENTRE LES TEMPS STOCKES NCAS ET NCAS+1
            DIFTEMP= TEMPS1 - TEMPS0

         ENDIF

C        ICI L'EF NEF CONTIENT LE POINT XYP0:
C        RECHERCHE DU POINT XYP1 = XYP0 + VITESSE * DELTAT
C        -------------------------------------------------
C        CALCUL DE LA VITESSE AU POINT XYP0 SELON L'INTERPOLATION
C        LA VALEUR DES NBNOE FONCTIONS DE BASE EN (XD,YD,0D0)
         CALL INTERP( NOINTE, CB2, CB3, 0D0, NBNOE, FBASE )

C        LA VALEUR DES 2 COMPOSANTES DE LA VITESSE EN CE POINT XYP0
C        INTERPOLATION LINEAIRE ENTRE LES TEMPS NCAS EST NCAS+1
         COEF0 = ( TEMPS1 - TEMPS  ) / DIFTEMP
         COEF1 = ( TEMPS  - TEMPS0 ) / DIFTEMP
         DO J=1,2
            MNN = MNELE + WUNDEL -1 + NEF
            DO L=1,NBNOE

C              LE NUMERO DU NOEUD L DU TRIANGLE NEF
               NOE = MCN(MNN)
               MNN = MNN + NBELEM
C              CAS DU TRIANGLE BREZZI-FORTIN
               IF( NOINTE .EQ. 2003 .AND. L .EQ. 4 ) NOE=NBSOMT+NEF

C              INTERPOLATION LINEAIRE ENTRE LES 2 TEMPS STOCKES
               IF( J .EQ. 1 ) THEN
                  VITEF(L) = COEF0 * vitx( NCAS   )%dptab( NOE )
     %                     + COEF1 * vitx( NCAS+1 )%dptab( NOE )
               ELSE
                  VITEF(L) = COEF0 * vity( NCAS   )%dptab( NOE )
     %                     + COEF1 * vity( NCAS+1 )%dptab( NOE )
               ENDIF
            ENDDO
            VXYP0(J) = PROSCD( FBASE, VITEF, NBNOE )
         ENDDO
         V = SQRT( VXYP0(1)**2 + VXYP0(2)**2 )

C        NORME DE LA VITESSE MAXIMALE DES PARTICULES SUR L'INTERVALLE DE TEMPS
         VPMAX  = MAX( VPMAX, V )

C        LE POINT XP1 YP1 OU ARRIVE LA PARTICULE APRES DELTAT
C        ----------------------------------------------------
         XP1 = XP0 + DELTAT * ( VXYP0(1) + XV0PART )
         YP1 = YP0 + DELTAT * ( VXYP0(2) + YV0PART )

C        RECHERCHE DE L'EF CONTENANT LE POINT XP1 YP1
C        --------------------------------------------
         CALL RETRIAXY1( XP0, YP0, XP1, YP1, MNNPEF, MNXYZP,
     %                   MOARET, MXARET, MNLARE,
     %                   NEF, CB1, CB2, CB3, IERR )
         IF( IERR .EQ. 0 ) THEN

            IF( NEF .GT. 0 ) THEN

C              LE POINT XYP1 EST DANS LE TRIANGLE NEF OU SUR UNE ARETE
C              -------------------------------------------------------
C              TRACE DU SEGMENT XP0-YP0  XP1-YP1
               XV0 = REAL( XP0 )
               YV0 = REAL( YP0 )
               XV1 = REAL( XP1 )
               YV1 = REAL( YP1 )

C              TRACE DU SEGMENT  XP0-YP0 XP1-YP1 AVEC UNE COULEUR
C              SELON LA VITESSE V DU FLUIDE EN XP0-YP0
               NCOUL = INT( (V-VMIN) / (VMAX-VMIN) * NBCOUL + N1COUL )

C              SELON LA VITESSE DE LA PARTICULE
               NCOUL=INT( N1COUL + V/VTPMAX * NBCOUL )

               IF( N1COUL .LE. NCOUL .AND. NCOUL .LE. NDCOUL ) THEN
                  CALL TRAIT2D( NCOUL, XV0, YV0, XV1, YV1 )
               ENDIF
C
               NBSEGM  = NBSEGM + 1
               PARCOUR = PARCOUR + SQRT( (XP1-XP0)**2 + (YP1-YP0)**2 )

               IF( NBSEGM .GT. MXSEGM ) THEN
                  IF( LANGAG .EQ. 0 ) THEN
                     PRINT *, 'PARTICULE',K,' au TEMPS',TEMPS,
     %                        ' MAXIMUM SEGMENTS ATTEINT =',NBSEGM
                  ELSE
                     PRINT *, 'PARTICLE',K,' at TIME',TEMPS,
     %                        'DRAWN SEGMENTS MAXIMUM ATTAINED=',NBSEGM
                  ENDIF
C                 PASSAGE A LA PARTICULE SUIVANTE
                  GOTO 110
               ENDIF
C
C              PASSAGE AU DELTAT SUIVANT
               XP0 = XP1
               YP0 = YP1
               GOTO 100

            ELSE

C              PAS DE TRIANGLE ADJACENT => LA PARTICULE SORT DU MAILLAGE
C              ---------------------------------------------------------
               IF( LANGAG .EQ. 0 ) THEN
                 PRINT *, 'PARTICULE',K,' au TEMPS',TEMPS,
     %                    ' SORT du MAILLAGE au POINT',
     %                     XP0,YP0,' avec une VITESSE',V
               ELSE
                 PRINT *, 'PARTICLE',K,' at TIME',TEMPS,
     %                    ' EXITS of MESH at POINT',
     %                     XP0,YP0,' with a VELOCITY',V
               ENDIF

C              TRACE DE LA PARTICULE AU POINT D'EXIT
               XV0 = REAL( XP0 )
               YV0 = REAL( YP0 )
               CALL SYMBOLE2D( NCNOIR, XV0, YV0, 'E' )

C              PASSAGE A LA PARTICULE SUIVANTE
               GOTO 110

            ENDIF

         ELSE

C           PROBLEME POUR TROUVER NEF CONTENANT XYP1
C           PASSAGE A LA PARTICULE SUIVANTE
            GOTO 110

         ENDIF

C        FIN DE LA BOUCLE SUR LES PAS DE TEMPS
 100     CONTINUE

C        BILAN DU PARCOURS DE LA PARTICULE K
 110     IF( NBORBIT .GT. 0 ) GOTO 200
         IF( LANGAG  .EQ. 0 ) THEN
            PRINT *, 'PARTICULE',K,
     %               ': X0=',RMCN(MNPA  ),' Y0=', RMCN(MNPA+1),
     %               ' VX0=',RMCN(MNPA+3),' VY0=',RMCN(MNPA+4),
     %               '   ',NBSEGM,' SEGMENTS TRACES',
     %               '   PARCOURS',PARCOUR,' VITESSE FINALE',V
         ELSE
            PRINT *, 'PARTICLE',K,
     %               ': X0=',RMCN(MNPA  ),' Y0=', RMCN(MNPA+1),
     %               ' VX0=',RMCN(MNPA+3),' VY0=',RMCN(MNPA+4),
     %               '   ',NBSEGM,' DRAWN SEGMENTS',
     %               '   TRAVEL',PARCOUR,' LAST VELOCITY',V
         ENDIF

C        FIN DU SUIVI DE LA PARTICULE K
 200  CONTINUE

      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'NORME de la VITESSE MAXIMALE d''UNE PARTICULE=',VPMAX
      ELSE
         PRINT*,'The PARTICLE MAXIMUM NORM of VELOCITY=',VPMAX
      ENDIF
      VTPMAX = VPMAX

C     LE TRACE DU TITRE FINAL
      VMIN = 0.
      VMAX = REAL( VTPMAX )
      CALL LEGZONTH( KNOMOB, NCAS0, 0, 5, VMIN, VMAX, VMIN, VMAX )

C     MISE SUR FICHIER NomfgifBoImage.xwd puis NomfgifNoImage.jpg
C     DE LA PIXMAP de la FENETRE X11 ACTUELLE
      CALL VIDEO1( NOMFGIF, NCAS0 )

C     ATTENDRE POUR LIRE LE TRACE
      CALL ATTENDSEC( TEMP2TRAC )

C     CONSTRUIRE le FICHIER VIDEO Nomfic.gif A PARTIR DES FICHIERS
C     CONSTRUITS de NOMS NomfgifNoImag.jpg
      CALL VIDEOFIN( NOMFGIF )

C     RETOUR POUR UNE NOUVELLE VISEE
      NBORBIT = NBORBIT + 1
      IF( LORBITE .NE. 0 ) THEN
C        ZOOM  BOUTON ENFONCE et DEPLACE
         CALL ZOOM2D1( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 20
         GOTO 30
      ELSE
C        POUR LIRE LE TRACE AVANT D'AFFICHER UN MENU
         CALL CLICSO
      ENDIF
      GOTO 20

C     SORTIE DU TRACE DU SUIVI DES PARTICULES EN 2D
C     =============================================
 9900 RETURN
      END
