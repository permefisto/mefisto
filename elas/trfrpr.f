      SUBROUTINE TRFRPR( KNOMOB, MODECO, NDIM,
     %                   NBTYEL, MNTOPO, MNNPEF, MNXYZN,
     %                   NBVECT, MNDEPL, MNTIME, NCAS,  AMPLID )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER L'OBJET DEFORME PAR SES DEPLACEMENTS AMPLIFIES
C -----    EN MOUVEMENT SELON LA FREQUENCE PROPRE DE L'OBJET
C
C ENTREES :
C ---------
C KNOMOB : NOM DE L'OBJET
C MODECO : MODE DE TRACE DES VECTEURS DU TMS D'ADRESSE MNDEPL
C          =1 CE SONT DEPLACEMENTS
C          =2 ou 3 CE SONT DES MODES PROPRES
C NDIM   : DIMENSION DE L'ESPACE =2 OU 3
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS FINIS DANS CETTE TOPOLOGIE
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE
C MNNPEF : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE
C MNXYZN : ADRESSE MCN DU TABLEAU XYZNOEUD DE L'OBJET
C NBVECT : NOMBRE DE VECTEURS DEPLACEMENTS
C MNDEPL : ADRESSE MCN DU TABLEAU DES VECTEURS DEPLACEMENT
C MNTIME : ADRESSE MCN DU TABLEAU DES TEMPS OU LES DEPLACEMENTS
C          ONT ETE CALCULES
C          0 SI PAS DE STOCKAGE
C
C MODIFIES EVENTUELLEMENT :
C -------------------------
C NCAS   : NUMERO DU CAS A TRAITER
C AMPLID : FACTEUR D'AMPLIFICATION DES DEPLACEMENTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Paris & VEULETTES sur MER AVRIL 2009
C23456---------------------------------------------------------------012
      PARAMETER        (LIGCON=0, NBVUES=5)
      IMPLICIT INTEGER (W)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
      include"./incl/ctemps.inc"
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      REAL             RMCN(1)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      DOUBLE PRECISION  DEXMAX, DECMAX, PI
      CHARACTER*(*)     KNOMOB
      CHARACTER*160     KNOM
      CHARACTER*4       NOMELE(2)
      INTRINSIC         REAL
C
C     SI PAS DE FREQUECE PROPRE => RETOUR
      IF( MODECO .LT. 2 ) RETURN
C
C     NOMBRE DE COMPOSANTES DU VECTEUR DEPLACEMENT
      NTDL   = MCN(MNDEPL+WBCOVE)
C     NOMBRE DE VECTEURS DEPLACEMENT
      NBVECT = MCN(MNDEPL+WBVECT)
C     NOMBRE DE NOEUDS DU MAILLAGE DE L'OBJET
      NBNOEU = MCN(MNXYZN+WNBNOE)
C     DIMENSION DE L'ESPACE DE L'OBJET (2 OU 3)
      NDIMLI = NTDL / NBNOEU
C     LE CAS TRAITE DOIT ETRE DANS LES LIMITES ACCEPTABLES
      NCAS   = MIN( NCAS, NBVECT )
      NCAS   = MAX( NCAS, 1 )
C
C     PI = 3.14159...
      PI = ATAN( 1D0 ) * 4D0
C
C     LA PALETTE PAR DEFAUT
      CALL PALCDE( 10 )
C
C     LECTURE DES DONNEES DU TRACE DES DEPLACEMENTS
C     =============================================
 100  CALL LIMTCL( 'tracdepl', NMTCL )
      IF( NMTCL .LE.  0 ) GOTO 9000
      IF( NMTCL .EQ. 50 ) GOTO 190
      IF( NMTCL .EQ. 90 ) GOTO 200
C
      GOTO( 110, 120, 100, 100, 150, 160, 170, 100, 140 ), NMTCL
C
C     NUMERO DU CAS A VISUALISER
 110  NCVALS = 4
      CALL INVITE( 84 )
      CALL LIRENT( NCVALS, NCAS )
      IF( NCVALS .EQ. -1 ) GOTO 100
C     PROTECTION DU NUMERO DE CAS A TRACER
      NCAS = MAX( 1, ABS(NCAS) )
      NCAS = MIN( NBVECT, NCAS )
C     LES TEMPS ONT ILS ETE STOCKES?
      IF( MNTIME .GT. 0 ) THEN
C        OUI: LE TEMPS INITIAL EST CELUI DU VECTEUR"DEPLACT
         TEMPS = RMCN( MNTIME + NCAS )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'TRACE au TEMPS ',TEMPS
         ELSE
            WRITE(IMPRIM,*) 'DRAWING at TIME ',TEMPS
         ENDIF
      ELSE
C        TEMPS INITIAL SUPPOSE NUL
         TEMPS = 0
      ENDIF
      GOTO 100
C
C     AMPLIFICATION DES DEPLACEMENTS
 120  NCVALS = 0
      CALL INVITE( 2 )
      CALL LIRRSP( NCVALS, AMPLID )
      IF( NCVALS .EQ. -1 ) GOTO 100
C     AMPLID PEUT ETRE NEGATIF POUR LE TRACE DES MODES PROPRES
      IF( AMPLID .EQ. 0. ) AMPLID = 1.
      GOTO 100
C
C     AFFICHAGE DES DEPLACEMENTS DU CAS NCAS
 140  CALL AFDEPL( NBNOEU, NCAS,   MNXYZN,
     %             NDIMLI, NBNOEU, NBVECT, MCN(MNDEPL+WECTEU),
     %             DECMAX, NOFOTI, DEXMAX )
      GOTO 100
C
C     COULEUR des ARETES du MAILLAGE
 150  CALL LIMTCL( 'couleur0' , I )
      IF( I .EQ. -1 ) THEN
         GOTO 9000
      ELSE IF( I .EQ. -2 ) THEN
         NCOUAF = -2
      ELSE IF( I .EQ. 0 ) THEN
C        COULEUR NOIRE
         NCOUAF = 0
      ELSE
C        COULEUR RESERVEE
         NCOUAF = N1COEL + I
      ENDIF
      GOTO 100
C
C     TYPE du TRAIT des ARETES du MAILLAGE
 160  CALL LIMTCL( 'typtrait' , I )
      IF( I .EQ. -1 ) GOTO 100
      NTLAFR = I
      GOTO 100
C
C     COULEUR des ARETES DEFORMEES
 170  CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 9000
      ELSE IF( I .EQ. 0 ) THEN
C        COULEUR NOIRE
         NCOUAD = 0
      ELSE
         NCOUAD = N1COEL + I
      ENDIF
      GOTO 100
C
C     EFFACER LE TRACE ACTUEL
 190  CALL EFFACE
C     PLUS D'ITEMS VISIBLES
      CALL ITEMS0
      CALL TRAXES
      GOTO 100
C
C     EXECUTION DU TRACE DE LA PIECE DEFORMEE
C     LE CODE DE TRAITEMENT DE L'INTERPOLATION DE L'OBJET
C     NDPGST : CODE TRAITEMENT DE L ELEMENT
C                 0 : NOEUDS=POINTS=SOMMETS
C                 1 : NOEUDS=POINTS#SOMMETS
C                 2 : NOEUDS#POINTS=SOMMETS
C                 3 : NOEUDS#POINTS#SOMMETS
 200  NDPGST = MCN( MNTOPO + WDPGST )
      IF( NDPGST .GE. 2 ) THEN
          NBLGRC(NRERR) = 2
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1) = 'TRACE AVEC NOEUDS # POINTS'
     %              // ' A L''AIDE DE FONCTIONS'
             KERR(2) = 'ERREUR : OPTION NON PROGRAMMEE'
          ELSE
             KERR(1) = 'DRAWING with NODES # POINTS FROM FUNCTIONS'
             KERR(2) = 'ERROR: OPTION NOT PROGRAMMED'
          ENDIF
          CALL LEREUR
          RETURN
      ENDIF
C
      IF( NDIM .LE. 2 ) THEN
C
C        ********
C        OBJET 2D
C        ********
C
C        OPTIONS DE LA VISEE POUR VOIR L'OBJET ET SES DEPLACEMENTS
 250     CALL VISE2D( NMTCL )
         IF( NMTCL .LT. 0 ) GOTO 100
C
C        INITIALISATION DE TRANSLATION ZOOM
         IF( LORBITE .NE. 0 ) THEN
            CALL ZOOM2D0( NOTYEV )
            IF( NOTYEV .EQ. 0 ) GOTO 250
         ENDIF
C
 275     IF( LORBITE .NE. 0 ) THEN
            CALL ZOOM2D1( NOTYEV )
            IF( NOTYEV .EQ. 0 ) GOTO 250
         ENDIF
C
C        TRACE DES NBVUES DE LA VIBRATION EN FREQUENCE PROPRE
         DO 350 NOVUE = 0, NBVUES
C
C           TRACE DES AXES 2D
            CALL TRAXE2
C
C           FACTEUR DE DEFORMATION POUR LA FREQUENCE PROPRE
            AMPLDE = REAL( AMPLID * COS( ( PI * NOVUE ) / NBVUES ) )
C
C           TRACE DES EF DEFORMES D'UN OBJET 2D
C           BOUCLE SUR LES DIFFERENTS TYPES D'ELEMENTS FINIS DU MAILLAGE
            DO 300 I = 0, NBTYEL-1
C
C              L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF"NOMEF
               MNELE = MCN( MNNPEF + I )
C              LE NUMERO DU TYPE D'EF
               NUTYEL = MCN( MNELE + WUTYEL )
C              LES CARACTERISTIQUES DE L'ELEMENT FINI
               CALL ELNUNM( NUTYEL, NOMELE )
C              LES CARACTERISTIQUES DE L'ELEMENT FINI
               CALL ELTYCA( NUTYEL )
C
C              TRACE EFFECTIF DES ELEMENTS FINIS AVEC ET SANS DEPLACEMENTS
               NBELEM = MCN(MNELE + WBELEM )
C              CE SP POUR BENEFICIER DES INDICES EN CLAIR
               CALL TRFRP2( AMPLDE, NCAS, NDIMLI,
     %                      NBELEM, MCN(MNELE+WBNDEL),MCN(MNELE+WUNDEL),
     %                      MCN(MNXYZN+WYZNOE),
     %                      MCN(MNDEPL+WBCOVE), MCN(MNDEPL+WBVECT),
     %                      MCN(MNDEPL+WECTEU) )
 300        CONTINUE
C
            IF( LANGAG .EQ. 0 ) THEN
               KNOM = 'AMPLIFICATION des DEPLACEMENTS= '
            ELSE
               KNOM = 'AMPLIFICATION of DISPLACEMENTS= '
            ENDIF
            WRITE( KNOM(32:44), '(G13.5)' ) AMPLID
            I = NUDCNB( KNOM )
            CALL XVTEXTE( KNOM(1:I), I, 50, 60 )
C
            IF( LANGAG .EQ. 0 ) THEN
               KNOM = 'OBJET DEFORME CAS '
            ELSE
               KNOM = 'DEFORMED OBJECT CASE '
            ENDIF
            I = NUDCNB( KNOM )
            WRITE( KNOM(I+1:I+4), '(I4)' ) NCAS
C
            I = NUDCNB( KNOM )
            IF( LANGAG .EQ. 0 ) THEN
               KNOM(I+1:I+12) = '  FREQUENCE '
            ELSE
               KNOM(I+1:I+12) = '  FREQUENCY '
            ENDIF
            I = NUDCNB( KNOM )
            WRITE( KNOM(I+1:I+14), '(G14.6)' ) TEMPS
            I = NUDCNB( KNOM )
            KNOM(I+1:I+4) = ' Hz '
C
C           TRACE FORCE DU TITRE
            IAVTIT = 1
            CALL TRFINS( KNOM )
C
 350     CONTINUE
C
C        RETOUR AU TRACE NORMAL POUR POSTSCRIPT
         IF ( LASOPS.NE.0 ) THEN
            LASOPS = LASOPS -10
            CALL XVPOSTSCRIPT(LASOPS)
         ENDIF
C
         IF( LORBITE .NE. 0 ) GOTO 275
         CALL CLICSO
         GOTO 250
C
      ELSE
C
C        ********
C        OBJET 3D
C        ********
C
C        TRACE DES ARETES FRONTALIERES NON DEFORMEES ET DES
C        FACES FRONTALIERES DEFORMEES (DEPLACEMENTS AMPLIFIES) D'UN OBJET 3D
C        ===================================================================
         CALL TRFRP3( PI, NBVUES, AMPLID, KNOMOB, MNXYZN,
     %                NCAS, NTDL, NBVECT, MCN(MNDEPL+WECTEU) )
C
      ENDIF
      GOTO 100
C
9000  RETURN
      END
