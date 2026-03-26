      SUBROUTINE ELADON( NUMIOB, NBOBIN, MNOBIN, NBOBCL, MNOBCL, MNDOEL,
     %                   IEMASS, IEYOUN, IEDILA, IECOED, IECOIN, IEFOIN,
     %                   IEFIXA, IEFOCL, IEFOPO,
     %                   IEDEIN, IEVIIN, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LES ADRESSES MCN DES DONNEES DE L'ELASTICITE
C -----    DES   SV "OBJETS INTERNES"    DE L'OBJET
C          DES PLS  "OBJETS AUX LIMITES" DE L'OBJET
C
C ENTREES:
C --------
C NUMIOB : NUMERO MINIMAL DU PLSV DANS LA DEFINITION DE L'OBJET
C NBOBIN : NOMBRE DE VOLUMES en 3D, SURFACES en 2D DE L'OBJET
C MNOBIN : ADRESSE MCN DU DEBUT DU TABLEAU NUOBIN (PARTIE DE TOPOLOGIE)
C NBOBCL : NOMBRE DE PLS en 3D, PL en 2D DE L'OBJET
C MNOBCL : ADRESSE MCN DU DEBUT DU TABLEAU NUOBCL (PARTIE DE TOPOLOGIE)
C MNDOEL : LES 4 ADRESSES MCN DES TABLEAUX DES ADRESSES DES
C          TABLEAUX DECRIVANT LES DONNEES ELASTICITE DE L'OBJET COMPLET
C
C SORTIES :
C ---------
C IEMASS : NOMBRE DE TMS MASSE            DES SV DE L'OBJET RETROUVES
C IEYOUN : NOMBRE DE TMS YOUNG (+POISSON) DES SV DE L'OBJET RETROUVES
C IEDILA : NOMBRE DE TMS DILATATION       DES SV DE L'OBJET RETROUVES
C IECOED : NOMBRE DE TMS COEFDEPLACEMENT  DES SV DE L'OBJET RETROUVES
C IECOIN : NOMBRE DE TMS CONTRINIT        DES SV DE L'OBJET RETROUVES
C IEFOIN : NOMBRE DE TMS FORCE "INTERNE"  DES SV DE L'OBJET RETROUVES
C
C IEFIXA : NOMBRE DE TMS FIXATION            DES PLS DE L'OBJET RETROUVES
C IEFOCL : NOMBRE DE TMS FORCE "AUX LIMITES" DES PLS DE L'OBJET RETROUVES
C IEFOPO : NOMBRE DE TMS FORCE               DES P   DE L'OBJET RETROUVES
C
C IEDEIN : NOMBRE DE TMS DEPLACTINIT         DES PLSV DE L'OBJET RETROUVES
C IEVIIN : NOMBRE DE TMS VITESSEINIT         DES PLSV DE L'OBJET RETROUVES
C
C IERR   : 0 SI PAS D'ERREUR, >0 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS       MARS 1999
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/donela.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           NUMIOB(4), MNDOEL(4)
C
      CHARACTER*80      KNOM,NOMTS
      CHARACTER*10      KNM
C
C     ============================
C     BOUCLE SUR LES SV "INTERNES"
C     ============================
      IERR   = 0
C     NOMBRE DE TMS RETROUVES SUR LES PLSV DE L'OBJET
      IEMASS = 0
      IEYOUN = 0
      IEDILA = 0
      IECOED = 0
      IECOIN = 0
      IEFOIN = 0
      IEFIXA = 0
      IEFOCL = 0
      IEFOPO = 0
      IEDEIN = 0
      IEVIIN = 0
C
C     POINTEUR SUR LE DEBUT DE OBIN
      MN = MNOBIN - 2
C
      DO 50 I=1,NBOBIN
C
C        LE TYPE DU SV INTERNE I
         MN   = MN + 2
         NYOB = MCN( MN )
         NUOB = MCN( MN + 1 )
C
C        LE NOM KNM DU SV
C        LE NOM DU SV DE TYPE KNM ET DE NO NUOB
C        LE PREFIXE NOMTS "PLSV>NOM>"
C        LE NUMERO DU CARACTERE DERRIERE LE DERNIER > DE NOMTS
         CALL NMPLSV( NYOB, NUOB,  KNM, KNOM, NOMTS, L )
C
         IF( NYOB .LE. 2 ) THEN
C           POINT OU LIGNE OBJET INTERNE => SA DIMENSION EST <=2 => ERREUR
C           RECUPERATION DANS KNOM DU NOM DU PLSV DE TYPE KNM ET DE NO NUOB
            NBLGRC(NRERR) = 2
            KERR(1) = KNM // ' : ' // KNOM
            IF( LANGAG .EQ. 0 ) THEN
               KERR(2) = 'NE PEUT ETRE INTERNE'
            ELSE
               KERR(2) = 'CAN''T BE INTERNAL'
            ENDIF
            CALL LEREUR
            IERR = IERR + 1
            GOTO 50
         ENDIF
C
C        OUVERTURE DU SV INTERNE I
         CALL LXNLOU( NTMN(NYOB), NUOB, NTOB, MNOB )
         IF( NTOB .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: ' // KNOM
               KERR(2) = KNM // ' INTERNE INCONNU'
            ELSE
               KERR(1) = 'ERROR: ' // KNOM
               KERR(2) = KNM // ' INTERNAL UNKNOWN'
            ENDIF
            CALL LEREUR
            IERR = IERR + 1
            GOTO 50
         ENDIF
C
C        L'ADRESSE DU DEBUT DES TABLEAUX DES DONNEES DU SV
C        -------------------------------------------------
         MN1 = MNDOEL( NYOB )
         MN1 = MN1 + MXDOEL * ( NUOB - NUMIOB(NYOB) ) - 1
C
C        OUVERTURE DU TABLEAU MASSE
C        --------------------------
         CALL LXTSOU( NTOB, 'MASSE', NT, MCN(MN1+LPMASS) )
         IF( NT .GT. 0 ) THEN
             IEMASS = IEMASS + 1
             NOMTS(L:80)='MASSE'
             CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU YOUNG
C        --------------------------
         CALL LXTSOU( NTOB, 'YOUNG', NT, MCN(MN1+LPYOUN) )
         IF( NT .GT. 0 ) THEN
             IEYOUN = IEYOUN + 1
             NOMTS(L:80)='YOUNG '
             CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU DILATATION
C        -------------------------------
         CALL LXTSOU( NTOB, 'DILATATION', NT, MCN(MN1+LPDILA) )
         IF( NT .GT. 0 ) THEN
             IEDILA = IEDILA + 1
             NOMTS(L:80)='DILATATION'
             CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU FORCE INTERNE
C        ----------------------------------
         CALL LXTSOU( NTOB, 'FORCE', NT, MCN(MN1+LPFORC) )
         IF( NT .GT. 0 ) THEN
             IEFOIN = IEFOIN + 1
             NOMTS(L:80)='FORCE '
             CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU COEFDEPLACEMENT
C        ------------------------------------
         CALL LXTSOU( NTOB, 'COEFDEPLACEMENT', NT, MCN(MN1+LPCOED) )
         IF( NT .GT. 0 ) THEN
             IECOED = IECOED + 1
             NOMTS(L:80)='COEFDEPLACEMENT'
             CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU CONTRINIT
C        ------------------------------
         CALL LXTSOU( NTOB, 'CONTRINIT', NT, MCN(MN1+LPCOIN) )
         IF( NT .GT. 0 ) THEN
             IECOIN = IECOIN + 1
             NOMTS(L:80)='CONTRINIT'
             CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU DEPLACEMENT INITIAL
C        ----------------------------------------
         CALL LXTSOU( NTOB, 'DEPLACTINIT', NT, MCN(MN1+LPDEIN) )
         IF( NT .GT. 0 ) THEN
             IEDEIN = IEDEIN + 1
             NOMTS(L:80)='DEPLACTINIT'
             CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU VITESSE INITIALE
C        -------------------------------------
         CALL LXTSOU( NTOB, 'VITESSEINIT', NT, MCN(MN1+LPVIIN) )
         IF( NT .GT. 0 ) THEN
             IEVIIN = IEVIIN + 1
             NOMTS(L:80)='VITESSEINIT'
             CALL AFTSTD( NOMTS )
         ENDIF
C
 50   CONTINUE
C
C     =========================================
C     BOUCLE SUR LES PLS AUX LIMITES DE L'OBJET
C     =========================================
C     POINTEUR SUR LE DEBUT DE OBCL
      MN = MNOBCL - 2
C
      DO 90 I=1,NBOBCL
C
C        LE TYPE DE PLS DU PLS I
         MN   = MN + 2
         NYOB = MCN( MN )
         NUOB = MCN( MN + 1 )
C
C        LE NOM KNM DU PLS
C        LE NOM DU PLS DE TYPE KNM ET DE NO NUOB
C        LE PREFIXE "PLSV>NOM>"
         CALL NMPLSV( NYOB, NUOB,  KNM, KNOM, NOMTS, L )
C
C        OUVERTURE DU PLS I
         CALL LXNLOU( NTMN(NYOB), NUOB, NTOB, MNOB )
         IF( NTOB .LE. 0 ) THEN
            NBLGRC(NRERR) = 3
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: AUX LIMITES '
               KERR(2) = KNOM
               KERR(3) = KNM // ' INCONNU'
            ELSE
               KERR(1) = 'ERROR: ON BOUNDARY '
               KERR(2) = KNOM
               KERR(3) = KNM // ' UNKNOWN'
            ENDIF
            CALL LEREUR
            IERR = IERR + 1
            GOTO 90
         ENDIF
C
C        L'ADRESSE DU DEBUT DES DONNEES DU PLS I
         MN1 = MNDOEL( NYOB )
         MN1 = MN1 + MXDOEL * ( NUOB - NUMIOB(NYOB) ) - 1
C
C        OUVERTURE DU TABLEAU FORCE FRONTALIERE
C        --------------------------------------
         CALL LXTSOU( NTOB, 'FORCE', NT, MCN(MN1+LPFORC) )
         IF( NT .GT. 0 ) THEN
C           NOMBRE DE FORCES AUX PLS DE L'OBJET
            IEFOCL = IEFOCL + 1
            IF( NYOB .EQ. 1 ) THEN
C              NOMBRE DE FORCES PONCTUELLES DE L'OBJET
               IEFOPO = IEFOPO + 1
            ENDIF
            NOMTS(L:80)='FORCE '
            CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU FIXATION
C        -----------------------------
         CALL LXTSOU( NTOB, 'FIXATION', NT, MCN(MN1+LPFIXA) )
         IF( NT .GT. 0 ) THEN
             IEFIXA = IEFIXA + 1
             NOMTS(L:80)='FIXATION '
             CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU DEPLACEMENT INITIAL
C        ----------------------------------------
         CALL LXTSOU( NTOB, 'DEPLACTINIT', NT, MCN(MN1+LPDEIN) )
         IF( NT .GT. 0 ) THEN
             IEDEIN = IEDEIN + 1
             NOMTS(L:80)='DEPLACTINIT'
             CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU VITESSE INITIALE
C        -------------------------------------
         CALL LXTSOU( NTOB, 'VITESSEINIT', NT, MCN(MN1+LPVIIN) )
         IF( NT .GT. 0 ) THEN
             IEVIIN = IEVIIN + 1
             NOMTS(L:80)='VITESSEINIT'
             CALL AFTSTD( NOMTS )
         ENDIF
C
 90   CONTINUE
      RETURN
      END
