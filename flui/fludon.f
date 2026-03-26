      SUBROUTINE FLUDON( NUMIOB, NBOBIN, MNOBIN, NBOBCL, MNOBCL, MNDOEL,
     %                   IEMASS, IEVISC, IECPRE, IEFOIN, IEVTIN, IEPRIN,
     %                   IEFOCL, IEFOPO, IEBLVI, IEBLPR, IEVIAN, IECBOU,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LES ADRESSES MCN DES DONNEES DU FLUIDE
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
C          TABLEAUX DECRIVANT LES DONNEES FLUIDE DE L'OBJET COMPLET
C
C SORTIES :
C ---------
C IEMASS : NOMBRE DE TMS MASSE           DES SV DE L'OBJET RETROUVES
C IEVISC : NOMBRE DE TMS VISCOSITE       DES SV DE L'OBJET RETROUVES
C IECPRE : NOMBRE DE TMS COEFPRESSION    DES SV DE L'OBJET RETROUVES
C IEFOIN : NOMBRE DE TMS FORCE "INTERNE" DES SV DE L'OBJET RETROUVES
C
C IEVTIN : NOMBRE DE TMS VITESSE  INITIALE DES SV DE L'OBJET RETROUVES
C IEPRIN : NOMBRE DE TMS PRESSION INITIALE DES SV DE L'OBJET RETROUVES
C
C IEFOCL : NOMBRE DE TMS FORCE "AUX LIMITES" DES PLS DE L'OBJET RETROUVES
C IEFOPO : NOMBRE DE TMS FORCE               DES P   DE L'OBJET RETROUVES
C
C IEBLVI : NOMBRE DE TMS BLVITESSE   DES PLS DE L'OBJET TROUVE
C IEBLPR : NOMBRE DE TMS BLPRESSION  DES PLS DE L'OBJET TROUVE
C
C IEVIAN : NOMBRE DE TMS VITESSE ANGULAIRE DES SV DE L'OBJET RETROUVES
C IECBOU : NOMBRE DE TMS COEFBOUSSINESQ DES SV DE L'OBJET RETROUVES
C
C
C IERR   : 0 SI PAS D'ERREUR, >0 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS : SOFIANE BENHAMADOUCHE  VINCENT BOYER            JANVIER 1999
C MODIFS  : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   NOVEMBRE 2000
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/donflu.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (RMCN(1),MCN(1))
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
      IEVISC = 0
      IECPRE = 0
      IEFOIN = 0
      IEFOCL = 0
      IEFOPO = 0
      IEBLVI = 0
      IEBLPR = 0
      IEVTIN = 0
      IEPRIN = 0
      IEVIAN = 0
      IECBOU = 0
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
               KERR(2) = KNM // ' UNKNOWN INTERNAL'
            ENDIF
            CALL LEREUR
            IERR = IERR + 1
            GOTO 50
         ENDIF
C
C        L'ADRESSE DU DEBUT DES TABLEAUX DES DONNEES DU SV
C        -------------------------------------------------
         MN1 = MNDOEL( NYOB )
         MN1 = MN1 + MXDOFL * ( NUOB - NUMIOB(NYOB) ) - 1
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
C        OUVERTURE DU TABLEAU VISCOSITE
C        ------------------------------
         CALL LXTSOU( NTOB, 'VISCOSITE', NT, MCN(MN1+LPVISC) )
         IF( NT .GT. 0 ) THEN
             IEVISC = IEVISC + 1
             NOMTS(L:80)='VISCOSITE'
             CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU COEFPRESSION
C        ---------------------------------
         CALL LXTSOU( NTOB, 'COEFPRESSION', NT, MCN(MN1+LPCPRE) )
         IF( NT .GT. 0 ) THEN
             IECPRE = IECPRE + 1
             NOMTS(L:80)='COEFPRESSION'
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
C        OUVERTURE DU TABLEAU VITESSE INITIALE
C        -------------------------------------
         CALL LXTSOU( NTOB, 'VITFLUIN', NT, MCN(MN1+LPVITF) )
         IF( NT .GT. 0 ) THEN
             IEVTIN = IEVTIN + 1
             NOMTS(L:80)='VITFLUIN'
            CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU PRESSION INITIALE
C        --------------------------------------
         CALL LXTSOU( NTOB, 'PREFLUIN', NT, MCN(MN1+LPPREF) )
         IF( NT .GT. 0 ) THEN
             IEPRIN = IEPRIN + 1
             NOMTS(L:80)='PREFLUIN'
             CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU VITESSE ANGULAIRE
C        --------------------------------------
         CALL LXTSOU( NTOB, 'VITESSEANGULAIRE', NT, MCN(MN1+LPVIAN) )
         IF( NT .GT. 0 ) THEN
             IEVIAN = IEVIAN + 1
             NOMTS(L:80)='VITESSEANGULAIRE'
             CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU coefboussines
C        -------------------------------
         CALL LXTSOU( NTOB, 'COEFBOUSSINESQ',
     %                NT, MCN(MN1+LPCBOU) )
         IF( NT .GT. 0 ) THEN
            IECBOU = IECBOU + 1
            NOMTS(L:80)='COEFBOUSSINESQ'
            CALL AFTSTD( NOMTS )
         ENDIF

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
               KERR(1) = 'ERREUR: AUX LIMITES'
               KERR(2) = KNOM
               KERR(3) = KNM // ' INCONNU'
            ELSE
               KERR(1) = 'ERROR: At BOUNDARY'
               KERR(2) = KNOM
               KERR(3) = KNM // ' is UNKNOWN'
            ENDIF
            CALL LEREUR
            IERR = IERR + 1
            GOTO 90
         ENDIF
C
C        L'ADRESSE DU DEBUT DES DONNEES DU PLS I
         MN1 = MNDOEL( NYOB )
         MN1 = MN1 + MXDOFL * ( NUOB - NUMIOB(NYOB) ) - 1
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
C        OUVERTURE DU TABLEAU CONDITION LIMITE SUR LES VITESSES
C        ------------------------------------------------------
         CALL LXTSOU( NTOB, 'BLVITESSE', NT, MCN(MN1+LPBLVI) )
         IF( NT .GT. 0 ) THEN
             IEBLVI = IEBLVI + 1
             NOMTS(L:80)='BLVITESSE'
             CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU CONDITIONS LIMITES SUR LA PRESSION
C        -------------------------------------------------------
         CALL LXTSOU( NTOB, 'BLPRESSION', NT, MCN(MN1+LPBLPR) )
         IF( NT .GT. 0 ) THEN
             IEBLPR = IEBLPR + 1
             NOMTS(L:80)='BLPRESSION'
             CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU VITESSE INITIALE
C        -------------------------------------
         CALL LXTSOU( NTOB, 'VITFLUIN', NT, MCN(MN1+LPVITF) )
         IF( NT .GT. 0 ) THEN
             IEVTIN = IEVTIN + 1
             NOMTS(L:80)='VITFLUIN'
             CALL AFTSTD( NOMTS )
         ENDIF
C
C        OUVERTURE DU TABLEAU PRESSION INITIALE
C        --------------------------------------
         CALL LXTSOU( NTOB, 'PREFLUIN', NT, MCN(MN1+LPPREF) )
         IF( NT .GT. 0 ) THEN
             IEPRIN = IEPRIN + 1
             NOMTS(L:80)='PREFLUIN'
             CALL AFTSTD( NOMTS )
         ENDIF
C
 90   CONTINUE
      RETURN
      END
