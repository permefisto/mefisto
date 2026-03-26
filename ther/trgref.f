      SUBROUTINE TRGREF( NOPROJ, KNOMOB, NTLXOB,
     %                   NBTYEL, MNELEM, MNPOGE, NDPGST,
     %                   NCAS0,  NCAS1,  NDIM,   NMTCL,  MODECO, NOPT,
     %                   GRADMX, GRADCM, CMFLEC, CMPGRA, TIMES )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DES GRADIENTS DES EF D'UN OBJET 1D 2D 3D ou 6D
C -----
C ENTREES:
C --------
C NOPROJ : SI OBJET EN 6D
C          TYPE DE PROJECTION 0 CI-DESSOUS FIXE LES COORDONNEES A ZERO
C          -1 PAS DE PROJECTION TRAITEMENT en XYZ NORMAL
C           1 : 'X Y Z 0 0 0'
C           2 : 'X Y 0 U 0 0'
C           3 : 'X 0 0 U V 0'
C           4 : '0 0 0 U V W'
C KNOMOB : NOM DE L'OBJET
C NTLXOB : NO LEXIQUE DE L'OBJET
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS FINIS DANS CETTE TOPOLOGIE
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE
C MNELEM : ADRESSE MCN DES TABLEAUX ELEMENTS FINIS DE CETTE TOPOLOGIE
C MNPOGE : ADRESSE MCN DU TABLEAU XYZPOINT
C NDPGST : CODE TRAITEMENT DES SOMMETS POINTS NOEUDS DU MAILLAGE
C          0 : NOEUDS=POINTS=SOMMETS
C          1 : NOEUDS=POINTS#SOMMETS
C          2 : NOEUDS#POINTS=SOMMETS
C          3 : NOEUDS#POINTS#SOMMETS
C NCAS0  : NUMERO DU PREMIER JEU DE SOLUTION A TRACER
C NCAS1  : NUMERO DU DERNIER JEU DE SOLUTION A TRACER
C NDIM   : DIMENSION DE L'ESPACE DE TRAVAIL 1 2 3 ou 6
C NMTCL  : NO DE LA DERNIERE OPTION CHOISIE DANS LE MENU
C MODECO : MODE DE TRACE DES VECTEURS DU TMS D'ADRESSE MNDEPL
C          =1 CE SONT TEMPERATURES
C          =2 CE SONT DES MODES PROPRES
C NOPT   : NO D'OPTION DE CALCUL DE L'UNITE DE GRADIENT
C GRADMX : MAXIMUM EN VALEUR ABSOLUE DES GRADIENTS
C GRADCM : GRADIENT DE 1 CM
C CMFLEC : NOMBRE DE CM DE LA PLUS GRANDE FLECHE DES GRADIENTS
C CMPGRA : NOMBRE DE CM POUR L'UNITE DE GRADIENT
C TIMES  : TEMPS DU CALCUL DES NBVECT VECTEURS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire JL LIONS UPMC Paris Decembre 2006
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2010
C23456...............................................................012
      PARAMETER    ( LIGCON=0, LIGTIR=1 )
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___dtemperature.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/a___aretefr.inc"
      include"./incl/ctemps.inc"
      include"./incl/xyzext.inc"
      include"./incl/traaxe.inc"
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (RMCN(1),MCN(1))
      CHARACTER*(*)     KNOMOB
      CHARACTER*90      KNOM
      CHARACTER*4       NOMELE(2)
      REAL              TIMES(NCAS0:NCAS1)
C
      MOREE2 = MOTVAR(6)
C
C     EXECUTION DU TRACE DES GRADIENTS DE LA TEMPERATURE
C     ==================================================
      IF( NOPT .EQ. 1 .AND. CMFLEC .LE. 0. ) GOTO 9999
      IF( NOPT .EQ. 2 .AND. GRADCM .LE. 0. ) GOTO 9999
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10310) NCAS0, NCAS1
      ELSE
         WRITE(IMPRIM,20310) NCAS0, NCAS1
      ENDIF
10310 FORMAT(' JEUX DE GRADIENTS TRACES =',I6,' A ',I6/)
20310 FORMAT(' NUMBER of DRAWN GRADIENT CASES =',I6,' to ',I6/)
C
C     MISE A JOUR DE L'ECHELLE
C     CMPGRA : NOMBRE DE CM POUR L'UNITE DE GRADIENT
      IF( NOPT .EQ. 1 ) THEN
         CMPGRA = CMFLEC / GRADMX
      ELSE IF( NOPT .EQ. 2 ) THEN
         CMPGRA = 1. / GRADCM
      ENDIF
C
C     LA PREPARATION DU TRACE 3D
C     --------------------------
      IF( NDIM .EQ. 3 ) THEN
C
C        CREATION OU REDECOUVERTE DU TMS  OBJET>>>FACE
         CALL HACHOB( KNOMOB, 4, NTFAOB, MNFAOB, IERR )
C
C        CREATION OU REDECOUVERTE DU TMS  OBJET>>>ARETEFR
C        DES ARETES DES FACES FRONTALIERES DE L'OBJET
         CALL HACHAF( KNOMOB, 0, NTFAOB, MNFAOB,
     %                NTAFOB, MNAFOB, I )
C
C        LE NOMBRE D'ENTIERS PAR ARETE FRONTALIERE
         MOARFR = MCN( MNAFOB + WOARFR )
C        LE NUMERO DANS LAREFR DE LA PREMIERE ARETE FRONTALIERE
         L1ARFR = MCN( MNAFOB + W1ARFR )
C
      ENDIF
C
C     LA VISEE INITIALE
C     -----------------
 320  CALL VISEE( NMTCL )
      IF( NMTCL .LT. 0 ) GOTO 9000
C
      IF( NDIM  .LE. 2 ) THEN
C        INITIALISATION DU ZOOM DEPLACEMENT
         CALL ZOOM2D0( NOTYEV )
      ELSE
C        INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
         CALL ORBITE0( NOTYEV )
      ENDIF
      IF( NOTYEV .EQ. 0 ) GOTO 320
C
C     BOUCLE SUR LES DIFFERENTS CAS A TRACER
C     ======================================
 340  DO NCAS = NCAS0, NCAS1
C
C        TEMPS DE CALCUL DU VECTEUR NCAS
         TEMPS = TIMES( NCAS )
C
C        L'ECRAN PIXELS EST EFFACE (PAS DE SCINTILLEMENT)
         CALL EFFACEMEMPX

C        TRACE EFFECTIF DES AXES
         NETAXE = 0
C
         IF( NDIM .EQ. 2 ) THEN
C
C           TRACE DES AXES
            CALL TRAXE2
C
         ELSE
C
C           TRACE DES AXES
            CALL TRAXE3
C
C           LE TRACE DES ARETES FRONTALIERES EN LIGNE NON EPAISSIE
            CALL XVEPAISSEUR( 1 )
C           LIGNE TIRETEE OU NON
            CALL XVTYPETRAIT( NTLAFR )
            CALL TRARFR( NCOUAF, MOARFR, L1ARFR, MCN(MNAFOB+WAREFR),
     %                   MCN(MNPOGE+WYZPOI) )
         ENDIF
C
C        BOUCLE SUR LES DIFFERENTS TYPES D'ELEMENTS FINIS DU MAILLAGE
         DO 350 I = 0, NBTYEL-1
C
C           L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF"TYPE_EF
            MNELE = MCN( MNELEM + I )
C
C           LE NOM DU TABLEAU GRADIENTS ASSOCIE
            NUTYEL = MCN( MNELE + WUTYEL )
C
C           LES CARACTERISTIQUES DE L'ELEMENT FINI
            CALL ELNUNM( NUTYEL, NOMELE )
            CALL ELTYCA( NUTYEL )
C
C           LE NOMBRE DE TELS ELEMENTS
            NBELEM = MCN(MNELE + WBELEM )
C
C           L'ADRESSE MCN DU TABLEAU 'NPEF' POUR CE TYPE D'EF
            MNPGEL = MNELE + WUNDEL
            IF( NDPGST .GE. 2 ) THEN
               MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
            ENDIF
C
            KNOM = 'DTEMPERATURE"' // NOMELE(2)
C           OUVERTURE DU TABLEAU
            CALL LXTSOU( NTLXOB, KNOM, NTGRAD, MNGRAD )
            IF( NTGRAD .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'OBJET SANS GRADIENT de TEMPERATURE'
               ELSE
                  KERR(1) = 'OBJECT WITHOUT the TEMPERATURE GRADIENT'
               ENDIF
               CALL LEREUR
               GOTO 350
            ENDIF
C
            IF( NDIM .LE. 2 ) THEN
C
C              TRACE DES ARETES DES EF 2D DE CE TYPE
C              -------------------------------------
               CALL XVTYPETRAIT( NTLAFR )
               CALL XVEPAISSEUR( 1 )
               CALL TRAREF( NDIM, MCN(MNPOGE+WNBPOI),MCN(MNPOGE+WYZPOI),
     %                      NBELEM, MCN(MNPGEL) )
C
            ENDIF
C
C           TRACE EFFECTIF DES FLECHES DU GRADIENT DE TEMPERATURE
C           -----------------------------------------------------
C           LIGNES EPAISSIES
            CALL XVEPAISSEUR( 2 )
            CALL XVTYPETRAIT( LIGCON )
C
C           LES DIFFERENTS TABLEAUX ET VARIABLES DE 'DTEMPERATURE'
C           AFIN DE BENEFICIER DES INDICES FORTRAN EN CLAIR
            NBELEM = MCN(MNGRAD+WBELFD)
            NDIMES = MCN(MNGRAD+WDIMED)
            NBJECA = MCN(MNGRAD+WBJECD)
            NBPIEF = MCN(MNGRAD+WBPEFD)
            MNCOPG = MNGRAD+WTEMPE+MOREE2*NBELEM*NDIMES*NBPIEF*NBJECA
            CALL TRGR123( NBELEM, NBPIEF, NDIMES, NBJECA,
     %                    MCN(MNCOPG), MCN(MNGRAD+WTEMPE),
     %                    NCAS, CMPGRA )
 350     CONTINUE
C
C        RETOUR AUX PARAMETRES INITIAUX
         CALL XVEPAISSEUR( 1 )
C
C        LE TRACE DU TITRE SELON MODECO
C        ------------------------------
         KNOM = 'OBJET: ' // KNOMOB
         I    = NUDCNB( KNOM )
         CALL XVCOULEUR( NCNOIR )
         CALL XVTEXTE( KNOM(1:I), I, 50, 30 )
C
         IF( LANGAG .EQ. 0 ) THEN
            KNOM = '1CM de FLECHE='
            WRITE( KNOM(15:27), '(G13.5)' ) 1.0 / CMPGRA
            KNOM(28:47) = ' UNITES du GRADIENT '
         ELSE
            KNOM = '1CM of ARROW='
            WRITE( KNOM(15:27), '(G13.5)' ) 1.0 / CMPGRA
            KNOM(28:48) = ' UNITIES of GRADIENT '
         ENDIF
         I = NUDCNB( KNOM )
         CALL XVTEXTE( KNOM(1:I), I, 50, 70 )
C
         IF( MODECO .EQ. 1 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               KNOM='CAS      Le GRADIENT de TEMPERATURE au TEMPS  '
            ELSE
               KNOM='CASE     The GRADIENT of TEMPERATURE at TIME  '
            ENDIF
            WRITE( KNOM(5:8),   '(I4)'    ) NCAS
            WRITE( KNOM(46:60), '(G15.7)' ) TEMPS
         ELSE
            IF( LANGAG .EQ. 0 ) THEN
         KNOM='VALEUR PROPRE               : Le GRADIENT de TEMPERATURE'
            ELSE
           KNOM='EIGENVALUE                  : The TEMPERATURE GRADIENT'
            ENDIF
            WRITE( KNOM(15:29), '(G15.7)' ) TEMPS
         ENDIF
C
C        AJOUT DE TEXTE SI OBJET 6D
         CALL PROJ6C( NOPROJ, KNOM )
C
C        TRACE DANS LA FENETRE
         CALL TRFINS( KNOM )
C
C        ATTENDRE POUR LIRE LE TRACE
         CALL ATTENDSEC( TEMP2TRAC )
C
C        FIN DE LA BOUCLE SUR LES CAS  ------------------------
      ENDDO
C
C     RETOUR POUR UNE NOUVELLE VISEE?
C     ------------------------------
      IF( LORBITE .EQ. 0 ) THEN
C        UN CLIC POUR DONNER LE TEMPS DE VOIR LE TRACE SANS MENU
         CALL CLICSO
         GOTO 320
      ENDIF
C
      IF( NDIM .GE. 3 ) THEN
         IF( NCAS0 .EQ. NCAS1 ) THEN
C           ORBITE BOUTON ENFONCE et DEPLACE
            CALL ORBITE1( NOTYEV )
         ELSE
C           ORBITE BOUTON 1 ou 2 ou 3 ENFONCE et DEPLACE et RELACHE
            CALL ORBITE3( NOTYEV )
         ENDIF
      ELSE
         IF( NCAS0 .EQ. NCAS1 ) THEN
C           ZOOM  BOUTON 1 ou 3 ENFONCE et DEPLACE
            CALL ZOOM2D1( NOTYEV )
         ELSE
C           ZOOM  BOUTON ENFONCE et DEPLACE et RELACHE
            CALL ZOOM2D3( NOTYEV )
         ENDIF
      ENDIF
      IF( NOTYEV .EQ. 0 ) GOTO 320
      GOTO 340
C
C     SORTIE SANS ERREUR
C     ==================
 9000 RETURN
C
C     ERREUR
C     ======
 9999 WRITE(KERR(MXLGER-3)(1:13),'(I13)')   NOPT
      WRITE(KERR(MXLGER-2)(1:13),'(G13.5)') CMFLEC
      WRITE(KERR(MXLGER-1)(1:13),'(G13.5)') GRADCM
      WRITE(KERR(MXLGER)(1:13),'(G13.5)')   CMPGRA
      NBLGRC(NRERR) = 5
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'ERREUR :  OPTION '
     &           // KERR(MXLGER-3)(1:13)
         KERR(2) = 'FLECHE MAXIMALE en CM '
     &           // KERR(MXLGER-2)(1:13)
         KERR(3) = 'GRADIENT d''UN CM de TRACE '
     &           // KERR(MXLGER-1)(1:13)
      ELSE
         KERR(1) = 'ERROR :  OPTION '
     &           // KERR(MXLGER-3)(1:13)
         KERR(2) = 'MAXIMUM ARROW in CM '
     &           // KERR(MXLGER-2)(1:13)
         KERR(3) = 'GRADIENT of 1 CM of DRAWING'
     &           // KERR(MXLGER-1)(1:13)
      ENDIF
      KERR(4) = 'CM / GRADIENT '
     &        // KERR(MXLGER)(1:13)
      CALL LEREUR
C
      RETURN
      END
