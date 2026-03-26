      SUBROUTINE TRERTH( NCAS,   NDIM,   NTLXOB,
     %                   NBTYEL, MNELEM, MNPOGE, NDPGST )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ESTIMATEURS D'ERREUR DU CAS NCAS EN REPRESENTANT
C -----    LE SAUT DES FLUX NORMAUX DE CHALEUR
C          AUX POINTS D'INTEGRATION NUMERIQUE DES FACES (OU ARETES)
C          LA NORME DU RESIDU F-DIV(K GRAD TEMPERATURE) DANS L2
C
C ENTREES :
C ---------
C NCAS   : NUMERO DU JEU DE TEMPERATURE A TRACER
C NDIM   : DIMENSION DE L'ESPACE DE L'OBJET ( 2 OU 3 )
C KNOMOB : NOM DE L'OBJET
C NTLXOB : NUMERO DU TMS LEXIQUE DE L'OBJET A TRAITER
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS DANS CETTE TOPOLOGIE
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE
C MNELEM : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE
C MNPOGE : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C NDPGST : CODE TRAITEMENT DES SOMMETS POINTS NOEUDS DU MAILLAGE
C          0 : NOEUDS=POINTS=SOMMETS
C          1 : NOEUDS=POINTS#SOMMETS
C          2 : NOEUDS#POINTS=SOMMETS
C          3 : NOEUDS#POINTS#SOMMETS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1995
C23456---------------------------------------------------------------012
      IMPLICIT INTEGER (W)
      PARAMETER     ( LIGCON=0, LIGTIR=1 )
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__erreurth.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___arete.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ctemps.inc"
      include"./incl/xyzext.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      DOUBLE PRECISION DMCN(1)
      EQUIVALENCE     (DMCN(1),MCN(1))
      CHARACTER*160     KNOM
      CHARACTER*4       NOMELE(2)
      REAL              HEXSEC(6,2)
C
      IF( NDIM .GT. 2 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
         KERR(1)='ESTIMATEUR ERREUR EN 2D et PB STATIONNAIRE SEULEMENT'
         ELSE
       KERR(1)='ERROR ESTIMATOR ONLY AVAILABLE in 2D and STEADY PROBLEM'
         ENDIF
         NBLGRC(NRERR) = 1
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
C
C     EN 2D : RECUPERATION DU TMS DES ARETES DE L'OBJET 2D
C     ====================================================
      CALL LXTSOU( NTLXOB, 'ARETE', NTARET, MNARET )
      IF( NTARET .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'OBJET SANS TABLEAU DES ARETES'
         ELSE
            KERR(1) = 'OBJECT WITHOUT THE ARRAY OF EDGES'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
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
C     LE NUMERO DE LA PREMIERE ARETE FRONTALIERE
      L1ARFB = MCN( MNARET + W1ARFB )
C     LE NUMERO DE LA PREMIERE ARETE INTERFACE
      L1ARIN = MCN( MNARET + W1ARIN )
C     ADRESSE MCN DU 1-ER MOT DU TABLEAU LARETE
      MNLARE = MNARET + W1LGFR + NUMXLF - NUMILF + 1
C
C     LE TMS OBJET>>ERREURTH
C     ======================
      CALL LXTSOU( NTLXOB, 'ERREURTH', NTERTH, MNERTH )
      IF( NTERTH .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'PAS D''ESTIMATEUR D''ERREUR POUR CET OBJET'
            KERR(2) = 'RESOUDRE LE PB THERMIQUE AVANT'
         ELSE
            KERR(1) = 'ERROR ESTIMATOR NOT COMPUTED'
            KERR(2) = 'SOLVE THERMAL PROBLEM BEFORE'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
      NBPTAF = MCN( MNERTH + WBPTAF )
      NDSM   = MCN( MNERTH + WBCAAF )
      IF( NCAS .GT. NDSM ) THEN
         WRITE(KERR(3)(1:6), '(I6)') NCAS
         WRITE(KERR(3)(7:12),'(I6)') NDSM
         NBLGRC(NRERR) = 4
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ESTIMATEUR d''ERREUR NON CALCULE pour le CAS'
     %              // KERR(3)(1:6)
            KERR(2) = 'NUMERO du DERNIER CAS CALCULE' // KERR(3)(7:12)
            KERR(3) = 'DERNIER CAS CALCULE est A IMPOSER si'
            KERR(4) = 'le PROBLEME est STATIONNAIRE SEULEMENT'
         ELSE
            KERR(1) = 'ERROR ESTIMATOR NOT COMPUTED for the CASE'
     %              // KERR(3)(1:6)
            KERR(2) = 'LAST COMPUTED CASE NUMBER' // KERR(3)(7:12)
            KERR(3) = 'IMPOSE THIS LAST COMPUTED CASE NUMBER if'
            KERR(4) = 'the PROBLEM is STEADY ONLY'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
      NBTTEF = MCN( MNERTH + WBTTEF )
C     LA REPARTITION INTERNE DES SOUS-TABLEAUX ET LEUR ADRESSE MCN
      MNSFLU = MNERTH + WFLUAF
      MNERAF = MNSFLU + MOREE2 * 2 * NBPTAF * MXARET * NDSM
      MNEREF = MNERAF + MOREE2 * MXARET * NDSM
      MNH1EF = MNEREF + MOREE2 * NBTTEF * NDSM
      MNEETH = MNH1EF + MOREE2 * NBTTEF * NDSM
      MNH1TE = MNEETH + MOREE2 * NDSM
      MNEEH1 = MNH1TE + MOREE2 * NDSM
      MNXYZC = MNEEH1 + MOREE2 * NDSM
C
C     QUELQUES INITIALISATIONS
C     ========================
      IERR   = 0
      CMPFLU = 0.0
      CMFLEC = 2.5
C
C     RECHERCHE DU MIN ET MAX DES COORDONNEES DES POINTS DE L'OBJET
C     =============================================================
      CALL MXXYZT( MOARET, MXARET, MCN(MNLARE),
     %             NBPTAF, MCN(MNXYZC), HEXSEC )
C
C     RECHERCHE DU MAXIMUM EN VALEUR ABSOLUE DES ESTIMATEURS D'ERREUR
C     ===============================================================
      CALL MXSAFL( NCAS,   NDSM,
     %             MOARET, MXARET, MCN(MNLARE),
     %             NBPTAF, MCN(MNSFLU),
     %             SFLUMX )
      IF( SFLUMX .LE. 0 ) SFLUMX = 1.0
      CMPFLU = 1.0 / SFLUMX
C
C     LECTURE DES DONNEES POUR DEFINIR LA TAILLE DES FLECHES
C     ======================================================
 100  CALL LIMTCL( 'tracerrt', NMTCL )
      IF( NMTCL .LE.  0 ) GOTO 9000
C
      GOTO( 110, 120, 100, 100, 150, 160, 170 ), NMTCL
C
C     NOMBRE DE CM POUR TRACER LA FLECHE MAXIMALE
C     -------------------------------------------
 110  NOPT   = 1
      NCVALS = 0
      CALL INVITE( 66 )
      CALL LIRRSP( NCVALS, CMFLEC )
      IF( NCVALS .EQ. -1 ) GOTO 100
C     CMFLEC : NOMBRE DE CM DU DE LA FLECHE MAXIMALE
      IF( CMFLEC .LT. 0. ) CMFLEC = -CMFLEC
      GOTO 300
C
C     1CM  VAUT EN FLUX
C     -----------------
 120  NOPT   = 2
      NCVALS = 0
      CALL INVITE( 102 )
      CALL LIRRSP(NCVALS,SFLUMX)
      IF( NCVALS .EQ. -1 ) GOTO 100
      IF( SFLUMX .EQ. 0. ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SAUT FLUX NORMAL NUL A CORRIGER >0'
         ELSE
            KERR(1) = 'NULL JUMP of NORMAL FLUX. To be >0'
         ENDIF
         CALL LEREUR
         GOTO 100
      ENDIF
      SFLUMX = ABS( SFLUMX )
      GOTO 300
C
C     COULEUR des ARETES du MAILLAGE
C     ------------------------------
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
C     ------------------------------------
 160  CALL LIMTCL( 'typtrait' , I )
      IF( I .EQ. -1 ) GOTO 100
      NTLAFR = I
      GOTO 100
C
C     COULEUR des FLECHES
C     -------------------
 170  CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 9000
      ELSE IF( I .EQ. 0 ) THEN
C        COULEUR NOIRE
         NCOUFL = 0
      ELSE
         NCOUFL = N1COEL + I
      ENDIF
      GOTO 100
C
C     EXECUTION DU TRACE DES FLUX NORMAUX DE TEMPERATURE
C     ==================================================
 300  IF( NOPT .EQ. 1 .AND. CMFLEC .LE. 0. ) GOTO 9999
      IF( NOPT .EQ. 2 .AND. SFLUMX .LE. 0. ) GOTO 9999
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10310) NCAS
      ELSE
         WRITE(IMPRIM,20310) NCAS
      ENDIF
10310 FORMAT(' JEU DE TEMPERATURE TRACE =',I5/)
20310 FORMAT(' TEMPERATURE CASE',I5,' IS DRAWN'/)
C
C     MISE A JOUR DE L'ECHELLE
C     CMPFLU : NOMBRE DE CM POUR L'UNITE DE FLUX
      IF( NOPT .EQ. 1 ) THEN
         CMPFLU = CMFLEC / SFLUMX
      ELSE IF( NOPT .EQ. 2 ) THEN
         CMPFLU = 1. / SFLUMX
      ENDIF
C
C     LA PREPARATION DU TRACE
 320  CALL VISE2D( NMTCL )
      IF( NMTCL   .LT. 0 ) GOTO 100
      IF( LORBITE .EQ. 0 ) GOTO 350
C
C     INITIALISATION DU ZOOM DEPLACEMENT
      CALL ZOOM2D0( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 320
      GOTO 350
C
C     ZOOM OU TRANSLATION ACTIFS
 330  CALL ZOOM2D1( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 320
C
C     BOUCLE SUR LES DIFFERENTS TYPES D'ELEMENTS FINIS DU MAILLAGE
C     ============================================================
 350  NUEF = 0
      DO 400 I = 0, NBTYEL-1
C
C        PARAMETRES PAR DEFAUT
         CALL XVEPAISSEUR( 1 )
         CALL XVTYPETRAIT( LIGCON )
C
C        L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF"TYPE_EF
         MNELE = MCN( MNELEM + I )
C
C        LE NOM DU TABLEAU FLUX ASSOCIE
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL, NOMELE )
         CALL ELTYCA( NUTYEL )
C
C        LE NOMBRE DE TELS ELEMENTS
         NBELEM = MCN(MNELE + WBELEM )
C
C        L'ADRESSE MCN DU TABLEAU 'NPEF' POUR CE TYPE D'EF
         MNPGEL = MNELE + WUNDEL
         IF( NDPGST .GE. 2 ) THEN
            MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
         ENDIF
C
C        LE TRACE DE L'ESTIMATEUR D'ERREUR DES EF DE CE TYPE
C        ---------------------------------------------------
         CALL TREEEF( NBELEM, MCN(MNPOGE+WYZPOI), MCN(MNPGEL),
     %                NCAS,   NDSM,
     %                NBTTEF, MCN(MNEREF), MCN(MNH1EF),
     %                MCN(MNH1TE), MCN(MNEEH1), NUEF )
CCCC
CCCC        LE TRACE DES ARETES DES EF 2D DE CE TYPE
CCCC        ----------------------------------------
CCC         CALL XVTYPETRAIT( NTLAFR )
CCC         CALL XVEPAISSEUR( 1 )
CCC         CALL TRAREF( NDIM, MCN(MNPOGE+WNBPOI), MCN(MNPOGE+WYZPOI),
CCC     %                NBELEM, MCN(MNPGEL) )
C
 400  CONTINUE
C
C     TRACE DE l'ESTIMATEUR D'ERREUR EN 2D
C     ------------------------------------
C     LIGNES EPAISSIES
      CALL XVEPAISSEUR( 2 )
      CALL XVTYPETRAIT( LIGCON )
C
C     TRACE DE L'ESTIMATEUR D'ERREUR EN 2D
C     ------------------------------------
      CALL TREEAR( NCAS,   NDSM, MOARET, MXARET, MCN(MNLARE),
     %             NBPTAF, MCN(MNSFLU), MCN(MNXYZC), CMPFLU )
C
C     RETOUR AUX PARAMETRES INITIAUX
      CALL XVEPAISSEUR( 1 )
C
C     LE TRACE DU TITRE
C     -----------------
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = '1CM de FLECHE='
      ELSE
         KNOM = '1CM of ARROW ='
      ENDIF
      WRITE( KNOM(15:27), '(G13.5)' ) 1.0 / CMPFLU
      IF( LANGAG .EQ. 0 ) THEN
         KNOM(28:72) = 'UNITES du SAUT de FLUX NORMAL de TEMPERATURE'
      ELSE
         KNOM(28:72) = 'UNITIES of JUMP of TEMPERATURE NORMAL FLUX'
      ENDIF
      I = NUDCNB( KNOM )
      CALL XVTEXTE( KNOM(1:I), I, 50, 70 )
C
C     TRACE DE LA LEGENDE DE L'ESTIMATEUR D'ERREUR
      NX     = LAPXFE - 300
      NY     = 350
      CALL XVCOULEUR( NCNOIR )
C
C     ESTIMATEUR D'ERREUR
      R = REAL( DMCN((MNEETH-1)/2+NCAS) )
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'ESTIMATEUR ERREUR  = '
      ELSE
         KNOM = 'ESTIMATOR of ERROR = '
      ENDIF
      WRITE(KNOM(22:33),'(G12.3)') R
      CALL XVTEXTE( KNOM(1:33), 33, NX, NY )
C
C     NORME H1 DE LA TEMPERATURE
      NY   = NY + 20
      R    = REAL( DMCN((MNH1TE-1)/2+NCAS) )
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'NORMEH1 TEMPERATURE= '
      ELSE
         KNOM = 'NORM H1 TEMPERATURE= '
      ENDIF
      WRITE(KNOM(22:33),'(G12.3)') R
      CALL XVTEXTE( KNOM(1:33), 33, NX, NY )
C
C     RAPPORT ESTIMATEUR ERREUR / NORME H1 DE LA TEMPERATURE
      NY   = NY + 20
      IF( ABS(DMCN((MNH1TE-1)/2+NCAS)) .LE. 1D-28 ) THEN
C        /0 EVITEE
         R = 0
      ELSE
         R = REAL( DMCN((MNEETH-1)/2+NCAS) / DMCN((MNH1TE-1)/2+NCAS) )
      ENDIF
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'ERREUR / NORME H1  ='
      ELSE
         KNOM = 'ERROR  / NORME H1  ='
      ENDIF
      WRITE(KNOM(22:33),'(G12.3)') R
      CALL XVTEXTE( KNOM(1:33), 33, NX, NY )
C
C     ESTIMATEUR ERREUR/NORME H1
      NY   = NY + 30
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'COULEURS  ERREUR/NORME H1'
      ELSE
         KNOM = 'COLORS   ERROR / NORM H1'
      ENDIF
      CALL XVTEXTE( KNOM(1:26), 26, NX, NY )
C
C     LE NOMBRE DE COULEURS DISPONIBLES DANS LA PALETTE
      NBCOUL = NDCOUL - N1COUL
      H      = NBCOUL / 10
      R      = REAL( DMCN((MNEEH1-1)/2+NCAS) )
      NX     = LAPXFE - 200
      NY     = LHPXFE - 30
C
      DO 600 I=0,10
C        LA COULEUR DU RAPPORT ESTIMATEUR / NORME H1 SUR LES EF
         NCOUL = N1COUL + NINT( H * I )
         CALL XVCOULEUR( NCOUL )
         CALL XVRECTANGLE( NX, NY, 30, 10 )
         WRITE( KNOM(1:4), '(I4)' ) I
         KNOM(5:7) = ' : '
         WRITE( KNOM(8:17), '(G10.3)' ) I * R / 10.
         CALL XVTEXTE( KNOM(1:17), 17, NX+30, NY+10 )
         NY = NY - 15
 600  CONTINUE
C
C     FIN DU TRACE
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'CAS      ESTIMATEUR D''ERREUR EN TEMPERATURE '
      ELSE
         KNOM = 'CASE     ESTIMATOR of TEMPERATURE ERROR      '
      ENDIF
      WRITE( KNOM(5:8),   '(I4)'    ) NCAS
      WRITE( KNOM(45:58), '(G14.6)' ) TEMPS
      CALL TRFINS( KNOM )
C
C     RETOUR A LA DEFINITION D'UNE NOUVELLE VUE
      IF( LORBITE .NE. 0 ) GOTO 330
      CALL CLICSO
      GOTO 320
C
C     SORTIE SANS ERREUR
C     ==================
 9000 RETURN
C
C     ERREUR
C     ======
 9999 WRITE(KERR(MXLGER-3)(1:13),'(I13)')   NOPT
      WRITE(KERR(MXLGER-2)(1:13),'(G13.5)') CMFLEC
      WRITE(KERR(MXLGER-1)(1:13),'(G13.5)') SFLUMX
      WRITE(KERR(MXLGER)(1:13),'(G13.5)')   CMPFLU
      NBLGRC(NRERR) = 5
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'ERREUR :  OPTION '
     &           // KERR(MXLGER-3)(1:13)
         KERR(2) = 'FLECHE MAXIMALE EN CM '
     &           // KERR(MXLGER-2)(1:13)
         KERR(3) = 'FLUX D''UN CM DE TRACE '
     &           // KERR(MXLGER-1)(1:13)
      ELSE
         KERR(1) = 'ERROR :  OPTION '
     &           // KERR(MXLGER-3)(1:13)
         KERR(2) = 'CM of MAX ARROW '
     &           // KERR(MXLGER-2)(1:13)
         KERR(3) = 'FLUX of 1 CM of DRAWING '
     &           // KERR(MXLGER-1)(1:13)
      ENDIF
      KERR(4) = 'CM / FLUX'
     &        // KERR(MXLGER)(1:13)
      CALL LEREUR
      GOTO 100
      END
