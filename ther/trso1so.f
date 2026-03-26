      SUBROUTINE TRSO1SO(INTERP, NOPROJ, MODECO, KNOMOB, NTLXOB, MNDFOB,
     %                   NBCOMP, NCAS0,  NCAS1,
     %                   NTYP,   SOLUTION, dptemp,
     %                   SOLMIN, SOLMAX, TIMES  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     TRACER DE LA SOLUTION SCALAIRE NCAS0 A NCAS1 EN CHAQUE NOEUD
C -----    (TEMPERATURE ou PRESSION ou ...) SUR UNE OU PLUSIEURS
C           SURFACES 3D NOMMEES DE L'OBJET 3D
C
C ENTREES:
C --------
C INTERP : NO D'INTERPOLATION A PRENDRE EN COMPTE
C          1 POLYNOME DE LAGRANGE DE DEGRE 1
C          2 POLYNOME DE LAGRANGE DE DEGRE 2
C          3 POLYNOME DE BREZZI-FORTIN
C NOPROJ : TYPE DE PROJECTION 0 CI-DESSOUS FIXE LA COORDONNEE A ZERO
C          -1 PAS DE PROJECTION TRAITEMENT en XYZ NORMAL
C           1 : 'X Y Z 0 0 0'
C           2 : 'X Y 0 U 0 0'
C           3 : 'X 0 0 U V 0'
C           4 : '0 0 0 U V W'
C MODECO : MODE DE TRACE DES VECTEURS DU TMS D'ADRESSE MNDEPL
C          = 1 CE SONT DES ISO-SOLUTIONS
C          = 2 CE SONT DES MODES PROPRES
C          = 3 CE SONT DES PRESSIONS P1 A PARTIR D'UNE INTERPOLATION
C              SOIT en 2D P2 SOIT P1+BULLE P3 OU en 3D P1 OU en 3D P2
C          = 4 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C          = 5 CE SONT DES SOLUTIONS         AUX NOEUDS D'UN MAILLAGE
C          = 6 CE SONT DES NORMES DE VITESSE AUX NOEUDS D'UN MAILLAGE
C          = 7 FONCTION COURANT  DANS UN FLUIDE
C          = 8 MODULE            D'UNE ONDE COMPLEXE NLSE
C          = 9 PARTIE REELLE     D'UNE ONDE COMPLEXE NLSE
C          =10 PARTIE IMAGINAIRE D'UNE ONDE COMPLEXE NLSE
C          =11 ERREURS DU MODULE DE LA VITESSE     D'UN FLUIDE
C          =12 ERREURS DU MODULE DE LA PRESSION P1 D'UN FLUIDE
C
C KNOMOB : NOM DE L'OBJET
C NTLXOB : NO TMS DU LEXIQUE DE L'OBJET
C MNDFOB : ADRESSE MCN DU TABLEAU DE DEFINITION DE L'OBJET
C
C NBCOMP : NOMBRE DE COMPOSANTE D'UN VECTEUR SOLUTION
C NCAS0:NCAS1 : NOMBRE DE VECTEURS SOLUTIONS
C NTYP    : =0 EMPLOI de SOLUTION
C           =1 EMPLOI de dptemp
C SOLUTION: VECTEUR"SOLUTION(NBCOMP,NCAS0:NCAS1)
C dptemp  : TABLEAU DES NCAS0:NCAS1 TABLEAU(NBCOMP) SOLUTIONS

C NCAS0  : NUMERO DU PREMIER CAS A TRAITER PARMI LES NCAS0:NCAS1 VECTEURS
C NCAS1  : NUMERO DU DERNIER CAS A TRAITER PARMI LES NCAS0:NCAS1 VECTEURS
C SOLMIN : SOLUTION MINIMALE DES NCAS0:NCAS1 CAS
C SOLMAX : SOLUTION MAXIMALE DES NCAS0:NCAS1 CAS
C TIMES  : TEMPS DU CALCUL   DES NCAS0:NCAS1 VECTEURS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY SEPTEMBRE 2010
C23456---------------------------------------------------------------012
      PARAMETER   ( MXTYEL=7 )
      PARAMETER   ( LIGCON=0 )
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___face.inc"
      include"./incl/a___aretefr.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ctemps.inc"
      include"./incl/donthe.inc"
      include"./incl/xvfontes.inc"
C
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      CHARACTER*(*)     KNOMOB
      CHARACTER*120     KNOM
      CHARACTER*24      NOMFGIF
      INTEGER           NUMIOB(4),  NUMAOB(4), MNDOEL(4),  MXDOEL(4)
      INTEGER           NONOEF(27), NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NONOFK(8)
      REAL              TIMES(NCAS0:NCAS1)
      DOUBLE PRECISION  SOLUTION(NBCOMP,NCAS0:NCAS1)
      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(NCAS0:NCAS1) :: dptemp

C     NOM DU FICHIER VIDEO  SELON MODECO // 'isos'
      CALL VIDEONM( MODECO, 'isos', NOMFGIF )
C
      MOREE2 = MOTVAR(6)
      SOLMIN0 = SOLMIN
      SOLMAX0 = SOLMAX
      MNBARY = 0
      MNNUFA = 0
      MNFACT = 0
C
C     LECTURE DU NOMBRE ET DES NOMS DES SURFACES DE TRACE DE LA SOLUTION
C     ------------------------------------------------------------------
      NBNOSU = 0
      MNNOSU = 0
      CALL LINMSURF( MNDFOB,  NBNOSU, MNNOSU )
      IF( NBNOSU .LE. 0 ) GOTO 9900
C
C     RECHERCHE DES TABLEAUX XYZSOMMET XYZNOEUD XYZPOINT NPEF" DE l'OBJET
C     FIND THE TMS XYZSOMMET XYZNOEUD XYZPOINT NPEF"xxxx  of the OBJECT
C     Cf $MEFISTO/td/da/a___xyznoeud   a___npef
C     ===================================================================
      CALL MIMAOB(      1, NTLXOB, MXDOTH,
     %             NTTOPO, MNTOPO,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %             MXTYEL, NBTYEL, MTNPEF, MNNPEF,
     %             NUMIOB, NUMAOB,
     %             NDPGST, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             MXDOEL, MNDOEL, IERR )
      IF( IERR .NE. 0 ) GOTO 9900
C     ICI NUMIOB CONTIENT LES 4 NUMEROS MINIMA DES PLSV
C         NUMAOB          LES 4 NUMEROS MAXIMA DES PLSV
C     MNDOEL LES 4 ADRESSES MCN DES TABLEAUX MIN MAX DES PLSV
C     MNDOEL THE 4 ADRESSES MCN of THE MIN MAX ARRAY MCN ADDRESSES OF PLSV
C     TABLEAUX DECRIVANT LA THERMIQUE DE L'OBJET COMPLET
C     SI EF BREZZI-FORTIN: XYZP EST REDUIT AUX SOMMETS
C                          XYZN CONTIENT LES SOMMETS ET BARYCENTRES
C
C     NOMBRE (3 ou 6) DES COORDONNEES DES NOEUDS
      NBCOOR = MCN( MNXYZN + WBCOON )
C
C     NOMBRE DE NOEUDS DU MAILLAGE DE L'OBJET
C     BREZZI-FORTIN => SOMMETS + BARYCENTRES
      NBNOEU = MCN( MNXYZN + WNBNOE )
C
C     DECLARATION DU TABLEAU DU NO SURFACE + NO SOMMETS DES FACES NOMMEES
      IF( INTERP .EQ. 1 ) THEN
C        4 SOMMETS POUR FACE P1 et Q1
         MOFACT = 4
      ELSE IF( INTERP .EQ. 2 ) THEN
C        4 SOMMETS + 4 MILIEUX POUR FACE P2 et Q2
         MOFACT = 8
      ELSE IF( INTERP .EQ. 3 ) THEN
C        3 SOMMETS + 0 POUR FACE P1 CAR BULLE=0 SUR LA FACE
         MOFACT = 4
      ENDIF
      MXFACT = NBNOEU
      CALL TNMCDC( 'ENTIER', MOFACT*MXFACT, MNFACT )
C
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
C     ----------------------------------------
      NBFACT = 0
      MNF    = MNFACT -1
      DO 40 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"TYPE EF
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM =  MCN( MNELE + WBELEM )
C
C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES EF DE CE TYPE
         MNNDEL = MNELE + WUNDEL
         MNPGEL = MNNDEL
         IF( NDPGST .GE. 2 ) THEN
C           POINTS DIFFERENTS DES NOEUDS
            MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
         ENDIF
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL, NOMELE )
         CALL ELTYCA( NUTYEL )
C
C        LA BOUCLE SUR LES EF DE CE TYPE NUTYEL
C        --------------------------------------
         DO 30 NUELEM = 1, NBELEM
C
C           LES NOEUDS DE L'ELEMENT FINI NUELEM
            CALL EFNOEU( MNELE, NUELEM, NBNDEL, NONOEF )
C
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
            CALL EFPLSV( MNELE , NUELEM,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                   NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C           PARCOURS DES FACES DE L'EF NUELEM
            DO NF = 1, NFACE
C
C              LA FACE NF EST ELLE SUR UNE DES SURFACES NOMMEES
               DO K = 1, NBNOSU
C
C                 NO DE LA SURFACE NOMMEE K
                  NOSURF = MCN(MNNOSU-1+K)
                  IF( NOSURF .EQ. NOOBSF(NF) ) THEN
C
C                    OUI: LA FACE EST SUR LA SURFACE NOSURF => AJOUT
                     IF( NBFACT .GE. MXFACT ) THEN
C                       TABLEAU DE TAILLE INSUFFISANTE
                        CALL TNMCAU( 'ENTIER',  MOFACT*MXFACT,
     %                         2*MOFACT*MXFACT, MOFACT*NBFACT, MNFACT )
                        MXFACT = 2 * MXFACT
                        MNF    = MNFACT + MOFACT*NBFACT -1
                     ENDIF
C
C                    UNE FACE DE PLUS A TRAITER
                     IF( INTERP .EQ. 1 .OR. INTERP .EQ. 3 ) THEN
C
C                       INTERPOLATION DEGRE 1 AVEC LES SEULS SOMMETS
C                       LES SOMMETS DE LA FACE NF DE L'EF
                        NB = NBSOFA(NF)
                        DO N = 1, NB
                           MCN( MNF + N ) = NONOEF( NOSOFA(N,NF) )
                        ENDDO
                        IF( NB .EQ. 3 ) MCN( MNF + 4 ) = 0
C
                     ELSE
C
C                       INTERPOLATION DEGRE 2 AVEC LES NOEUDS
C                       RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE NF
                        CALL ELNOFA( NUTYEL, NF, NB, NONOFK )
C
C                       LES NB NOEUDS DE LA FACE NF DE L'EF
                        DO N = 1, NB
                           MCN( MNF + N ) = NONOEF( NONOFK(N) )
                        ENDDO
                        IF( NB .EQ. 6 ) THEN
C                          TRIANGLE P2
                           MCN( MNF + 7 ) = 0
                           MCN( MNF + 8 ) = 0
                        ENDIF
C
                     ENDIF
C
C                    UNE FACE DE PLUS
                     NBFACT = NBFACT + 1
                     MNF    = MNF + MOFACT
C
                  ENDIF
               ENDDO
            ENDDO
 30      CONTINUE
 40   CONTINUE
C
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'NOMBRE DE FACES DES SURFACES NOMMEES=',NBFACT
         WRITE(IMPRIM,*) 'VECTEUR SOLUTIONS TRAITES',NCAS0,' A ', NCAS1
      ELSE
         WRITE(IMPRIM,*) 'NUMBER of FACES of NAMED SURFACES=',NBFACT
         WRITE(IMPRIM,*) 'TREATED SOLUTION VECTORS ',NCAS0,' TO ',NCAS1
      ENDIF
C
C     TABLEAU DES BARYCENTRES DES FACES DES SURFACES
      IF( MNBARY .GT. 0 ) CALL TNMCDS( 'REEL',   NBFACT, MNBARY )
      IF( MNNUFA .GT. 0 ) CALL TNMCDS( 'ENTIER', NBFACT, MNNUFA )
      CALL TNMCDC( 'REEL',   NBFACT, MNBARY )
      CALL TNMCDC( 'ENTIER', NBFACT, MNNUFA )
C
C     OPTIONS DE TRACE DES FACES DES SURFACES
C     ---------------------------------------
 95   CALL LIMTCL( 'solsurf3', NMTCL )
      IF( NMTCL .LT. 0 ) GOTO 9900
      IF( NMTCL .EQ. 5 ) THEN
C
C        COULEUR des ARETES de la FRONTIERE
C        ..................................
         CALL LIMTCL( 'couleur0' , I )
         IF( I .EQ. -1 ) GOTO 100
         IF( I .EQ. -2 ) THEN
C           COULEUR INVISIBLE A NE PAS TRACER
            NCOAFR = -2
         ELSE IF( I .EQ. 0 ) THEN
C           LA COULEUR NOIRE
            NCOAFR = 0
         ELSE
            NCOAFR = N1COEL + I
         ENDIF
         GOTO 95
C
      ELSE IF( NMTCL .EQ. 6 ) THEN
C
C        TYPE du TRAIT des ARETES de la FRONTIERE
C        ........................................
         CALL LIMTCL( 'typtrait' , I )
         IF( I .EQ. -1 ) GOTO 100
         NTLAFR = I
         GOTO 95
C
      ENDIF
C
ccc   IF( NMTCL .EQ. 90 ) GOTO 100
C
C     ------------------------------------------------------------
C     OPTIONS DE VISEE POUR VOIR LA SOLUTION SUR LES SURFACES
C     ------------------------------------------------------------
 100  CALL VISE3D( NMTCL )
      IF( NMTCL .LT. 0 ) GOTO 95
C
C     INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
      IF( LORBITE .NE. 0 ) THEN
         CALL ORBITE0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 100
      ENDIF
C
C     FORMATION DES TABLEAUX NUMERO DES FACES ET DISTANCE AXONOMETRIQUE
C     POUR LES FACES DES SURFACES NOMMEES DE L'OBJET
      NBP = 0
 120  NBP = NBP + 1
      CALL TRBFSFN( NBFACT, MOFACT, MXFACT, MCN(MNFACT),
     %              RMCN(MNXYZN+WYZNOE),
     %              MCN(MNNUFA), RMCN(MNBARY) )
C
C     CALCUL DU MIN ET MAX DE LA SOLUTION SUR LES SURFACES NOMMEES
      IF( NBP .EQ. 1 ) THEN
         CALL TRSO1MIMX( NBFACT, MCN(MNNUFA),
     %                   MOFACT, MXFACT, MCN(MNFACT),
     %                   NBCOOR, RMCN(MNXYZN+WYZNOE),
     %                   NBCOMP, NCAS0, NCAS1,
     %                   NTYP,   SOLUTION, dptemp,
     %                   SOLMIN, SOLMAX )
C        POUR PRENDRE EN COMPTE LES NOUVEAUX MIN MAX DES COORDONNEES
         NOTYVI = 0
C        TRACE DEMANDE
         CALL VISEE1
      ENDIF
C
C     LE TRI PAR TAS DE CETTE DISTANCE
C     LA FACE LA PLUS PROCHE EST LA PREMIERE TRACEE
      CALL TRITRP( NBFACT, RMCN(MNBARY), MCN(MNNUFA) )
C
C     TRACE DES DIFFERENTS CAS CALCULES
      DO NCAS = NCAS0, NCAS1
C
C        L'ECRAN PIXELS EST EFFACE (PAS DE SCINTILLEMENT)
         CALL EFFACEMEMPX
C
C        TRACE DES AXES 3D
         CALL TRAXE3
C
C        TRACE DES ARETES FRONTALIERES ET DE LA SOLUTION SUR
C        LES FACES FRONTALIERES
         CALL TRSO1SO3( NBFACT, MCN(MNNUFA),
     %                  MOFACT, MXFACT, MCN(MNFACT),
     %                  NBCOOR, RMCN(MNXYZN+WYZNOE),
     %                  NBCOMP, NCAS0, NCAS1,
     %                  NTYP, SOLUTION, dptemp, NCAS,
     %                  SOLMIN, SOLMAX )
C
C        EFFACEMENT DE LA LEGENDE DU TRACE POSTSCRIPT
         CALL EFLEGPOS
C
C        LE TRACE DU TITRE FINAL
         KNOM = 'OBJET: ' // KNOMOB
         I    = NUDCNB( KNOM )
         CALL XVCOULEUR( NCNOIR )
         CALL XVTEXTE( KNOM(1:I), I, 50, 30 )
C
C        LE TRACE DE LA LEGENDE : COULEURS => VALEURS
         NBCOUL = NDCOUL - N1COUL
         NCPAS  = NBCOUL / 10
         TPAS   = (SOLMAX-SOLMIN) / 10
         T      = SOLMIN
C
C        TRACE DE 11 VALEURS
         NCOUL = N1COUL
         NX    = LAPXFE - 160
         NY    = LHPXFE - 30
         DO I=0,10
            CALL XVCOULEUR( NCOUL )
            CALL XVRECTANGLE( NX, NY, 30, 10 )
            WRITE( KNOM(1:10), '(G10.3)' ) T
            CALL XVTEXTE( KNOM(1:10), 10, NX+40, NY+10 )
            NCOUL = NCOUL + NCPAS
            T     = T  + TPAS
            NY    = NY - 15
         ENDDO
C
C        TRACE DU MIN ET MAX DE LA SOLUTION
         KNOM =' '
         WRITE(KNOM,10000) NCAS0, NCAS1, SOLMIN, SOLMAX
10000    FORMAT('From CASE',I5,' to',I5,' MINIMUM=',G14.6,'  MAXIMUM=',
     %           G14.6)
         I=NUDCNB( KNOM )
C        SAUVEGARDE DU NUMERO DE LA FONTE DE CARACTERES ACTUELS
         NOFONT0 = NOFONT
C        CHANGEMENT DE POLICE DE CARACTERES POUR UNE DE 24 PIXELS DE HAUT
         LHPXCA = 20
         CALL CHOIXFONTE( LHPXCA )
         CALL XVCOULEUR( NCJAUN )
         CALL XVTEXTE( KNOM(1:I), I, 50, 100 )
C        RESTAURATION DU NUMERO DE LA FONTE DE CARACTERES ACTUELS
         CALL CHARGEFONTE( NOFONT0 )
C
C        RETOUR AU TRACE NORMAL POUR POSTSCRIPT
         IF ( LASOPS .NE. 0 ) THEN
           LASOPS = LASOPS - 10
           CALL XVPOSTSCRIPT(LASOPS)
         ENDIF
C
C        RETOUR AUX PARAMETRES INITIAUX
         CALL XVEPAISSEUR( 1 )
         CALL XVTYPETRAIT( LIGCON )
C
C        DEFINITION DU TITRE ET FIN DU TRACE
         TEMPS = TIMES( NCAS )
         CALL LETITR( NOPROJ, MODECO, NCAS, TEMPS, KNOM )
C
C        MISE SUR FICHIER NomfgifBoImage.xwd puis NomfgifNoImage.jpg
C        DE LA PIXMAP de la FENETRE X11 ACTUELLE
         CALL VIDEO1( NOMFGIF, NCAS )
C
C        ATTENDRE POUR LIRE LE TRACE
         CALL ATTENDSEC( TEMP2TRAC )
C
C        FIN DE LA BOUCLE SUR LES CAS
      ENDDO
C
C     CONSTRUIRE le FICHIER VIDEO Nomfic.gif A PARTIR DES FICHIERS
C     CONSTRUITS de NOMS NomfgifNoImag.jpg
C     ------------------------------------------------------------
      CALL VIDEOFIN( NOMFGIF )
C
C     RETOUR POUR UNE NOUVELLE VISEE
C     ------------------------------
      IF( LORBITE .NE. 0 ) THEN
         IF( NCAS0 .EQ. NCAS1 ) THEN
C           ORBITE BOUTON ENFONCE et DEPLACE
            CALL ORBITE1( NOTYEV )
         ELSE
C           ORBITE BOUTON ENFONCE et DEPLACE et RELACHE
            CALL ORBITE3( NOTYEV )
         ENDIF
         IF( NOTYEV .EQ. 0 ) GOTO 100
         GOTO 120
      ELSE
C        POUR LIRE LE TRACE AVANT D'AFFICHER UN MENU
         CALL CLICSO
      ENDIF
      GOTO 100
C
C     SORTIE DU TRACE DES ARETES ET FACES
C     ===================================
 9900 SOLMIN = SOLMIN0
      SOLMAX = SOLMAX0
      IF( MNFACT .GT. 0 ) CALL TNMCDS( 'ENTIER', MOFACT*MXFACT, MNFACT )
      IF( MNBARY .GT. 0 ) CALL TNMCDS( 'REEL',   NBFACT, MNBARY )
      IF( MNNUFA .GT. 0 ) CALL TNMCDS( 'ENTIER', NBFACT, MNNUFA )
      IF( MNNOSU .GT. 0 ) CALL TNMCDS( 'ENTIER', NBNOSU, MNNOSU )
      RETURN
      END
