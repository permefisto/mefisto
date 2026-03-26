      SUBROUTINE AFVIPRER( KNOMOB, KNOMFIC, NDIM, NUTYEL, MNXYZN, NODDL,
     %                     NBNOVI, NBNOPR, NTDLVP, NTDLTE,
     %                     NCAS0,  NCAS1,  NOEUDV1, NOEUDV2,
     %                     TIMES,  VXYZPN, TEMPER, NOFOVI, NOFOPR,
     %                     VITMIN, VITMAX, VITMOY, PREMIN,PREMAX,PREMOY)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LES VITESSES ET LES PRESSIONS CALCULEES
C -----    ET SI FONCTION VITESSE_EXACTE(t,x,y,z,nc) est DONNEE
C                L'ERREUR SUR LE MODULE DE LA VITESSE AUX NOEUDS
C          ET SI FONCTION PRESSION_EXACTE(t,x,y,z) est DONNEE
C                L'ERREUR SUR LA PRESSION AUX SOMMETS DES EF
C          POUR TOUS LES PAS DE TEMPS STOKES de NCAS0 a NCAS1
C
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET CONTENANT LE FLUIDE
C KNOMFIC: NOM DU FICHIER SUPPORT D'UN VECTEUR VITESSE+PRESSION
C NDIM   : DIMENSION 2 OU 3 DE L'ESPACE DE L'OBJET
C NUTYEL : 13 TRIANGLE  BREZZI FORTIN 2D
C          15 TRIANGLE  TAYLOR-HOOD   2D
C          19 TETRAEDRE BREZZI FORTIN 3D
C          20 TETRAEDRE TAYLOR-HOOD   3D
C MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD' DU MAILLAGE DE L'OBJET
C NODDL  : TABLEAU DU NUMERO DU DERNIER D.L. DE CHAQUE NOEUD
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
C          SOMMETS + BARYCENTRES        des TETRAEDRES POUR BREZZI-FORTIN
C          SOMMETS + MILIEUX DES ARETES des TETRAEDRES POUR TAYLOR-HOOD
C NBNOPR : NOMBRE DE NOEUDS SUPPORT DE LA PRESSION
C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTES EN VITESSES PRESSIONS

C NCAS0  : PREMIER VECTEUR VITESSE-PRESSION A TRAITER
C NCAS1  : DERNIER VECTEUR VITESSE-PRESSION A TRAITER
C NOEUDV1: NUMERO DU PREMIER NOEUD D'ERREUR A AFFICHER
C NOEUDV2: NUMERO DU DERNIER NOEUD D'ERREUR A AFFICHER

C TIMES  : LES NCAS0:NCAS1 TEMPS DES CALCULS DE LA VITESSE-PRESSION
C VXYZPN : UNE CARTE NCAS DES VITESSES-PRESSIONS AUX NOEUDS DU MAILLAGE

C NOFOVI : NUMERO DE LA FONCTION VITESSE_EXACTE(t,x,y,z,nc)
C          =0 SI ELLE N'EST PAS DONNEE PAR L'UTILISATEUR
C NOFOPR : NUMERO DE LA FONCTION PRESSION_EXACTE(t,x,y,z)
C          =0 SI ELLE N'EST PAS DONNEE PAR L'UTILISATEUR

C SORTIES:
C --------
C VITMIN : MINIMUM DE LA NORME DE LA VITESSE CALCULEE EN UN NOEUD
C VITMAX : MAXIMUM DE LA NORME DE LA VITESSE CALCULEE EN UN NOEUD
C VITMOY : MOYENNE DE LA NORME DE LA VITESSE CALCULEE EN TOUS LES NOEUDS
C PREMIN : MINIMUM DE LA PRESSION CALCULEE EN UN NOEUD DU MAILLAGE
C PREMAX : MAXIMUM DE LA PRESSION CALCULEE EN UN NOEUD DU MAILLAGE
C PREMOY : MOYENNE DE LA PRESSION CALCULEE EN TOUS LES SOMMETS DU MAILLAGE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR  Fevrier 2012
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray              Aout 2020
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray               Mai 2021
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/ctemps.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      CHARACTER*(*)     KNOMOB
      CHARACTER*80      KNOMFIC
      INTEGER           NODDL(0:NBNOVI)
      REAL              TIMES(NCAS0:NCAS1)
      DOUBLE PRECISION  VXYZPN(NTDLVP), TEMPER(*),
     %                  VITNOR, VITMIN, VITMAX, VITMOY,
     %                  PREMIN, PREMAX, PREMOY

C     CALCUL DES ERREURS AVEC UNE SOLUTION EXACTE TAYLOR-GREEN VERTEX
      DOUBLE PRECISION  XP, YP, ZP, DPARAF(5),
     %                  XVelocity, YVelocity, ZVelocity, PresExac,
     %                  VEXNORM, VEXAMX,   VERROR,   VERRELA,  VERRLMX,
     %                  VERRMX,  VEXASOM,  VERRSOM,
     %                  PRERRELA,PRERRLMX, PRERRABS, PRERRAMX, PRECALC,
     %                  PREXSOM, PRERRSOM, PREXMX,   PREXMI,
     %                  PCERVIT, PCERVITMOY, PCERPRE,  PCERPREMOY,
     %                  PRCAMI,  PRCAMX
      INTRINSIC         SQRT

C     RECHERCHE DE LA FONCTION UTILISATEUR DONNANT LA
C     VITESSE_EXACTE(t,x,y,z,nocomp) ou EXACT_VELOCITY(t,x,y,z,nocomp)
      NOFOVI = NOFOVITE()

C     RECHERCHE DE LA FONCTION UTILISATEUR DONNANT LA
C     PRESSION_EXACTE(t,x,y,z) ou EXACT_PRESSURE(t,x,y,z)
      NOFOPR = NOFOPRES()

C     ERREUR RELATIVE MOYENNE EN TEMPS SUR LA VITESSE
      PCERVITMOY = 0D0
C     ERREUR RELATIVE MOYENNE EN TEMPS SUR LA PRESSION
      PCERPREMOY = 0D0
      PRERRELA   = 0D0

C     BOUCLE SUR LES TEMPS DES CALCULS
C     ================================
      MNXYZ = MNXYZN + WYZNOE - 3

      DO 1000 NCAS = NCAS0, NCAS1

C     LE TEMPS DU CALCUL
      TEMPS = TIMES(NCAS)

C     TRANSFERT DE LA VITESSE-PRESSION AU TEMPS NCAS PAR NOEUDS
C     FICHIER DU REPERTOIRE PROJET DANS LE TABLEAU MNVXYZP0
      CALL LIFIVIPRTE( KNOMOB,  TEMPS,  NCAS,   NAVSTO, NBPASDT,
     %                 NDIM,    NBNOVI, NBNOPR, NBNOTE,
     %                 NTDLVP,  VXYZPN, NTDLTE, TEMPER,
     %                 KNOMFIC, IERR )
C
C     L'ADRESSE MCN DES COORDONNEES DES NOEUDS VITESSE
C     ICI LE TMS XYZNOEUD A ETE ENRICHI DES BARYCENTRES DES EF
C     DANS LES CAS DES EF DE BREZZI-FORTIN
      MN = MNXYZ

      IF( NOFOVI .LE. 0 .OR. NOFOPR .LE. 0 ) THEN

C     =============================================================
C     AFFICHAGE DES COMPOSANTES CALCULEES DE LA VITESSE ET PRESSION
C     NON CONNAISSANCE DE LA VITESSE EXACTE OU PRESSION EXACTE
C     =============================================================

      IF( LANGAG .EQ. 0 ) THEN

C        AFFICHAGE EN FRANCAIS
         WRITE(IMPRIM,10019) TEMPS, NOEUDV1, NOEUDV2, NBNOVI, NTDLVP
10019    FORMAT(//'Au temps ',G15.6,' la VITESSE et PRESSION des noeuds'
     %   ,I9,' a',I9,' parmi',I9,' NOEUDS et',I10,' DL:'/140('-'))
C
         DO 20 I = NOEUDV1, NOEUDV2
            MN  = MN  + 3
            NDL = NODDL(I-1)
            ND  = NODDL(I) - NDL
            IF( NUTYEL .EQ. 19 .OR.  NUTYEL .EQ. 20 ) THEN

C              3D BF ou TH : VITESSE"PRESSION
               IF( ND .GT. NDIM ) THEN
                  WRITE(IMPRIM,10023) I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND)
               ELSE
                  WRITE(IMPRIM,10033) I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND)
               ENDIF
c
            ELSE
C
C              2D: NOMBRE DE DL EN LE NOEUD I
               IF( ND .GT. NDIM ) THEN
                  WRITE(IMPRIM,10022) I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND)
               ELSE
                  WRITE(IMPRIM,10032) I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND)
               ENDIF
            ENDIF
 20      ENDDO

10022 FORMAT('Noeud ',I9,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VITESSE VX=',G15.7,' VY=',G15.7,' PRESSION=',G15.7)

10032 FORMAT('Noeud ',I9,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VITESSE VX=',G15.7,' VY=',G15.7)

10023 FORMAT('Noeud ',I9,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VITESSE VX=',G15.7,' VY=',G15.7,' VZ=',G15.7,
     %'  PRESSION=',G15.7)

10033 FORMAT('Noeud ',I9,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VITESSE VX=',G15.7,' VY=',G15.7, ' VZ=',G15.7 )
C
      ELSE
C
C        ENGLISH PRINTING
         WRITE(IMPRIM,20019) TEMPS, NOEUDV1, NOEUDV2, NBNOVI, NTDLVP
20019    FORMAT(//'At time ',G15.6,': VELOCITY and PRESSURE at Nodes',
     %     I9,' to',I9,' among',I9,' Nodes and ',I10,' DoF:'/140('-'))
C
         DO 30 I = NOEUDV1, NOEUDV2
            MN  = MN  + 3
            NDL = NODDL(I-1)
            ND  = NODDL(I) - NDL
C
C           AFFICHAGE DES COMPOSANTES DE LA VITESSE ET
C           EVENTUELLEMENT DE LA PRESSION
            IF( NDIM .EQ. 3 ) THEN
C
C              ELEMENT FINI 3D
               IF( ND .GT. NDIM ) THEN
C                 NOEUD AVEC PRESSION
                  WRITE(IMPRIM,20023) I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND)
               ELSE
C                 NOEUD SANS PRESSION
                  WRITE(IMPRIM,20033) I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND)
               ENDIF
c
            ELSE
C
C              ELEMENT FINI 2D
               IF( ND .GT. NDIM ) THEN
C                 NOEUD AVEC PRESSION
                  WRITE(IMPRIM,20022) I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND)
               ELSE
C                 NOEUD SANS PRESSION
                  WRITE(IMPRIM,20032) I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND)
               ENDIF
            ENDIF
 30      ENDDO

      ENDIF

20022 FORMAT('Node',I10,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VELOCITY VX=',G15.7,' VY=',G15.7,'  PRESSURE=',G15.7)
C
20032 FORMAT('Node',I10,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VELOCITY VX=',G15.7,' VY=',G15.7)
C
20023 FORMAT('Node',I10,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VELOCITY VX=',G15.7,' VY=',G15.7,' VZ=',G15.7,
     %'  PRESSURE=',G15.7)
C
20033 FORMAT('Node',I10,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VELOCITY VX=',G15.7,' VY=',G15.7,' VZ=',G15.7)

      ELSE

C     ==================================================================
C     AFFICHAGE DES ERREURS SUR LA VITESSE ET LA PRESSION
C     ICI SONT SUPPOSEES CONNUES ET APPELABLES LES FONCTIONS UTILISATEUR
C     VITESSE_EXACTE(t,x,y,z,nocomp) ou EXACT_VELOCITY(t,x,y,z,nocomp)
C     PRESSION_EXACTE(t,x,y,z) ou EXACT_PRESSURE(t,x,y,z)
C     ==================================================================
      IF( LANGAG .EQ. 0 ) THEN
C        AFFICHAGE EN FRANCAIS
         WRITE(IMPRIM,10110) TEMPS, NOEUDV1, NOEUDV2, NBNOVI, NTDLVP
10110    FORMAT(//'Au temps ',G15.6,' ERREURS sur la VITESSE et PRESSION
     % aux Noeuds ',I9,' a ',I9,' parmi',I9,
     %' Noeuds et ',I10,' D.L.):'/140('-'))
      ELSE
C        ENGLISH PRINTING
         WRITE(IMPRIM,20110) TEMPS, NOEUDV1, NOEUDV2, NBNOVI, NTDLVP
20110    FORMAT(//'At time ',G15.6,' ERRORS on VELOCITY and PRESSURE at
     %Nodes',I9,' to ',I9,' among',I9,' Nodes and ',I10,' DoF:'
     %    /140('-'))
      ENDIF

C     INITIALISATIONS POUR LE CALCUL DES ERREURS
      VEXASOM = 0D0
      VERRSOM = 0D0
      VERRMX  =-1D100
      VEXAMX  = 0D0
      VERRLMX =-1D100
      IVERMX  = 0

      PREXSOM  = 0D0
      PRCAMI   = 1D100
      PRCAMX   =-1D100
      PREXMI   = 1D100
      PREXMX   =-1D100
      PRERRLMX =-1D100
      PRERRAMX = 0D0
      PRERRSOM = 0D0

C     NOMBRE DE NOEUDS VITESSE+PRESSION IMPRIMES
      NBVPIMP = 0

      DO 120 I=1,NBNOVI

C        ADRESSE DES COORDONNEES DU NOEUD I
         MN = MN + 3
         XP = RMCN(MN)
         YP = RMCN(MN+1)
         ZP = RMCN(MN+2)

C        NOMBRE DE DL AU NOEUD I
         NDL = NODDL(I-1)
         ND  = NODDL(I) - NDL

         DPARAF(1) = TEMPS
         DPARAF(2) = XP
         DPARAF(3) = YP
         DPARAF(4) = ZP

         IF( ND .GT. NDIM ) THEN

C           PRESSION_EXACTE(t,x,y,z,nocomp) ou
C           EXACT_VELOCITY(t,x,y,z,nocomp) AU NOEUD ET AU TEMPS NCAS
C           LA PRESSION_EXACTE(temps,xp,yp,zp) AU NOEUD I
C           ---------------------------------------------
            CALL FONVAL( NOFOPR, 4, DPARAF, NCODEV, PresExac )

C           CALCUL DES ERREURS LOCALES SUR LA PRESSION AU NOEUD I
            PREXMI  = MIN( PREXMI, PresExac )
            PREXMX  = MAX( PREXMX, PresExac )
            PREXSOM = PREXSOM + ABS( PresExac )
C
            PRECALC = VXYZPN(NODDL(I))
            PRCAMI  = MIN( PRCAMI, PRECALC )
            PRCAMX  = MAX( PRCAMX, PRECALC )
C
            PRERRABS = ABS( PRECALC - PresExac )
            PRERRSOM = PRERRSOM + PRERRABS
            PRERRAMX = MAX( PRERRAMX, PRERRABS )

C           MAXIMUM DE L'ERREUR RELATIVE EN UN NOEUD DU MAILLAGE
            PRERRELA = PRERRAMX / (PREXMX-PREXMI)
C
         ENDIF
C
C        LA VITESSE EXACTE(temps,xp,yp,zp,NoComposante)
C        ----------------------------------------------
C        NUMERO DE LA COMPOSANTE TRAITEE DE LA VITESSE
         DPARAF(5) = 1
C        COMPOSANTE X DE LA VITESSE EXACTE AU NOEUD I AU TEMPS NCAS
         CALL FONVAL( NOFOVI, 5, DPARAF, NCODEV, XVelocity )

         DPARAF(5) = 2
C        COMPOSANTE Y DE LA VITESSE EXACTE AU NOEUD I AU TEMPS NCAS
         CALL FONVAL( NOFOVI, 5, DPARAF, NCODEV, YVelocity )

         IF( NDIM .GE. 3 ) THEN
            DPARAF(5) = 3
C           COMPOSANTE Z DE LA VITESSE EXACTE AU NOEUD I AU TEMPS NCAS
            CALL FONVAL( NOFOVI, 5, DPARAF, NCODEV, ZVelocity )
         ELSE
            ZVelocity = 0D0
         ENDIF

C        CALCUL DES ERREURS LOCALES SUR LA VITESSE
         IF( NDIM .GE. 3 ) THEN
            VEXNORM= SQRT(  XVelocity**2 + YVelocity**2 + ZVelocity**2 )
            VERROR = SQRT( (XVelocity-VXYZPN(NDL+1))**2
     %                   + (YVelocity-VXYZPN(NDL+2))**2
     %                   + (ZVelocity-VXYZPN(NDL+3))**2 )
         ELSE
            VEXNORM = SQRT(  XVelocity**2 + YVelocity**2 )
            VERROR  = SQRT( (XVelocity-VXYZPN(NDL+1))**2
     %                    + (YVelocity-VXYZPN(NDL+2))**2 )
         ENDIF
C        SOMME DES NORMES DE LA VITESSE EXACTE AUX NOEUDS
         VEXASOM = VEXASOM + VEXNORM
C        MAXIMUM DES NORMES DE LA VITESSE EXACTE AUX NOEUDS
         VEXAMX  = MAX( VEXAMX, VEXNORM )
C        SOMME   DE L'ERREUR SUR LA NORME DE LA VITESSE AUX NOEUDS
         VERRSOM = VERRSOM + VERROR
C        MAXIMUM DE L'ERREUR SUR LA NORME DE LA VITESSE AUX NOEUDS
         IF( VERROR .GT. VERRMX ) THEN
            VERRMX = VERROR
            IVERMX = I
         ENDIF
C
         IF( VEXNORM .LT. 1D-6) THEN
            VERRELA = 0D0
         ELSE
            VERRELA = VERROR / VEXNORM
         ENDIF
C        MAXIMUM DE L'ERREUR RELATIVE EN UN NOEUD DU MAILLAGE
         VERRLMX  = MAX( VERRLMX, VERRELA )

C        LIMITATION DU NOMBRE DE LIGNES AFFICHEES
C        A AU MOINS 10 AVEC PRESSION C-A-D AUX SOMMETS
C        ---------------------------------------------
         IF( I .GT. NBNOPR/2 .and. I .LE. NBNOPR/2+50 .and.
     %       NBVPIMP .LE. 10 ) THEN

         IF( NDIM .EQ. 3 ) THEN

C           ELEMENT FINI 3D
            IF( ND .GT. NDIM ) THEN
C              NOEUD AVEC PRESSION
               NBVPIMP = NBVPIMP + 1
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10023) I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND)
                  WRITE(IMPRIM,10324)    (RMCN(MN+K),K=0,2),
     %                  XVelocity, YVelocity, ZVelocity,
     %                  PresExac, VERRELA*100, PRERRELA*100
               ELSE
                  WRITE(IMPRIM,20023) I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND)
                  WRITE(IMPRIM,20324)    (RMCN(MN+K),K=0,2),
     %                  XVelocity, YVelocity, ZVelocity,
     %                  PresExac, VERRELA*100, PRERRELA*100
               ENDIF

            ELSE
C              NOEUD SANS PRESSION
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10033) I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND)
                  WRITE(IMPRIM,10334)    (RMCN(MN+K),K=0,2),
     %                  XVelocity, YVelocity, ZVelocity, VERRELA*100
               ELSE
                  WRITE(IMPRIM,20033) I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND)
                  WRITE(IMPRIM,20334)    (RMCN(MN+K),K=0,2),
     %                  XVelocity, YVelocity, ZVelocity, VERRELA*100
               ENDIF

            ENDIF

         ELSE

C           ELEMENT FINI 2D
            IF( ND .GT. NDIM ) THEN
C              NOEUD AVEC PRESSION
               NBVPIMP = NBVPIMP + 1
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10022) I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND)
                  WRITE(IMPRIM,10124)    (RMCN(MN+K),K=0,2),
     %                  XVelocity, YVelocity,
     %                  PresExac, VERRELA*100, PRERRELA*100
               ELSE
                  WRITE(IMPRIM,20022) I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND)
                  WRITE(IMPRIM,20124)    (RMCN(MN+K),K=0,2),
     %                  XVelocity, YVelocity,
     %                  PresExac, VERRELA*100, PRERRELA*100
               ENDIF

            ELSE
C              NOEUD SANS PRESSION
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10032) I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND)
                  WRITE(IMPRIM,10134)    (RMCN(MN+K),K=0,2),
     %                  XVelocity, YVelocity, VERRELA*100
               ELSE
                  WRITE(IMPRIM,20032) I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND)
                  WRITE(IMPRIM,20134)    (RMCN(MN+K),K=0,2),
     %                  XVelocity, YVelocity, VERRELA*100
               ENDIF

            ENDIF
         ENDIF
         ENDIF

 120  ENDDO

C     AFFICHAGES FINAUX DES ERREURS
C     -----------------------------
      PCERVIT = VERRSOM / VEXASOM * 100
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10501) VEXAMX
10501    FORMAT('Max|VITESSE|Exacte(Noeud)=',G14.6)
         WRITE(IMPRIM,10502) VERRMX, VERRMX/VEXAMX*100, PCERVIT
10502 FORMAT('MaxErreur|Vitesse|(Noeud)=',g14.6,
     % '  MaxErreur|Vitesse|(Noeud)/Max|Vitesse|(Noeud)=',g14.6,'%',
     % '  SomErreur|Vitesse|(Noeud)/Som|Vitesse|(Noeud)=',G14.6,'%')
      ELSE
         WRITE(IMPRIM,20501) VEXAMX
20501    FORMAT('Max|Velocity|Exact(Node)=',G14.6)
         WRITE(IMPRIM,20502) VERRMX, VERRMX/VEXAMX*100, PCERVIT
20502 FORMAT('Max|Velocity|Error(Node)=',G14.6,
     % '  Max|Velocity|Error(Node)/Max|Velocity|(Node)=',G14.6,'%',
     % '  Sum|Velocity|Error(Node)/Sum|Velocity|(Node)=',G14.6,'%')
      ENDIF

C     ADRESSE DES COORDONNEES DU NOEUD IVERMX
      MN = MNXYZ + 3*IVERMX
      XP = RMCN(MN)
      YP = RMCN(MN+1)
      ZP = RMCN(MN+2)

      DPARAF(1) = TEMPS
      DPARAF(2) = XP
      DPARAF(3) = YP
      DPARAF(4) = ZP

C     LA VITESSE EXACTE(temps,xp,yp,zp,NoComposante)
C     NUMERO DE LA COMPOSANTE TRAITEE DE LA VITESSE
      DPARAF(5) = 1
C     COMPOSANTE X DE LA VITESSE EXACTE AU NOEUD I AU TEMPS NCAS
      CALL FONVAL( NOFOVI, 5, DPARAF, NCODEV, XVelocity )

      DPARAF(5) = 2
C     COMPOSANTE Y DE LA VITESSE EXACTE AU NOEUD I AU TEMPS NCAS
      CALL FONVAL( NOFOVI, 5, DPARAF, NCODEV, YVelocity )

      IF( NDIM .GE. 3 ) THEN
         DPARAF(5) = 3
C        COMPOSANTE Z DE LA VITESSE EXACTE AU NOEUD I AU TEMPS NCAS
         CALL FONVAL( NOFOVI, 5, DPARAF, NCODEV, ZVelocity )
      ELSE
         ZVelocity = 0D0
      ENDIF

      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10499) VEXAMX, VERRMX, IVERMX, XP,YP,ZP,
     %                       XVelocity, YVelocity, ZVelocity,
     %                      (VXYZPN(NODDL(IVERMX-1)+K),K=1,NDIM)
      ELSE
         WRITE(IMPRIM,20499) VEXAMX, VERRMX, IVERMX, XP,YP,ZP,
     %                       XVelocity, YVelocity, ZVelocity,
     %                      (VXYZPN(NODDL(IVERMX-1)+K),K=1,NDIM)
      ENDIF

10499 FORMAT('Max|VITESSE|(Noeuds)=',G14.6,
     %'  Max|VITESSE|ERREUR(Noeuds)=',G14.6,
     %' au noeud',I9,' X=',G14.6,' Y=',G14.6,' Z=',G14.6/
     %'Vx Exact=',G14.6,' Vy Exact=',G14.6,' Vz Exact=',G14.6/
     %'Vx Calcu=',G14.6,' Vy Calcu=',G14.6,' Vz Calcu=',G14.6)
20499 FORMAT('Max|VELOCITY|(Nodes)=',G14.6,
     %' Max|VELOCITY|ERROR(Nodes)=',G14.6,
     %' at node',I9,' X=',G14.6,' Y=',G14.6,' Z=',G14.6/
     %'Vx Exact=',G14.6,' Vy Exact=',G14.6,' Vz Exact=',G14.6/
     %'Vx Compu=',G14.6,' Vy Compu=',G14.6,' Vz Compu=',G14.6)
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10505) PREXMI, PREXMX, PREXMX-PREXMI
         WRITE(IMPRIM,10515) PRCAMI, PRCAMX, PRCAMX-PRCAMI
      ELSE
         WRITE(IMPRIM,20505) PREXMI, PREXMX, PREXMX-PREXMI
         WRITE(IMPRIM,20515) PRCAMI, PRCAMX, PRCAMX-PRCAMI
      ENDIF

10505 FORMAT(/'Min Exact     PRESSION=',G14.6,
     %      '  Max Exact     PRESSION=',G14.6,
     %      '  Max-Min Exact PRESSION=',G14.6)
10515 FORMAT( 'Min Calcul    PRESSION=',G14.6,
     %      '  Max Calcul    PRESSION=',G14.6,
     %      '  Max-Min CalculPRESSION=',G14.6)
20505 FORMAT(/'Min Exact     PRESSURE=',G14.6,
     %      '  Max Exact     PRESSURE=',G14.6,
     %      '  Max-Min Exact PRESSURE=',G14.6)
20515 FORMAT( 'Min Computed  PRESSURE=',G14.6,
     %      '  Max Computed  PRESSURE=',G14.6,
     %      '  Max-Min Compu PRESSURE=',G14.6)

      PCERPRE = PRERRSOM / PREXSOM * 100
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10503) PRERRAMX, PRERRAMX/(PREXMX-PREXMI)*100,
     %                       PCERPRE
      ELSE
         WRITE(IMPRIM,20503) PRERRAMX, PRERRAMX/(PREXMX-PREXMI)*100,
     %                       PCERPRE
      ENDIF
10503 FORMAT( 'Max ERREUR    PRESSION=',G14.6,
     %      '  MaxPressionERREUR/Max-MinPression=',G14.6,'%',
     %       ' SumPressionERREUR/SumPression=',G14.6,'%')
20503 FORMAT( 'Max ERROR     PRESSURE=',G14.6,
     %      '  MaxPressureERROR/Max-MinPressure=',G14.6,'%',
     %       ' SumPressureERROR/SumPressure=',G14.6,'%')

10124 FORMAT('Exacte SOLUTION: X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VITESSE VX=',G15.7,' VY=',G15.7,' PRESSION=',G15.7,
     %'  V-ERREUR=',G12.4,'%   P-ERREUR=',G12.4,'%'/)

10324 FORMAT('Exacte SOLUTION: X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VITESSE VX=',G15.7,' VY=',G15.7,' VZ=',G15.7,
     %'  PRESSION=',G15.7,
     %'  V-ERREUR=',G12.4,'%   P-ERREUR=',G12.4,'%'/)

20124 FORMAT('Exact SOLUTION: X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VELOCITY VX=',G15.7,' VY=',G15.7,'  PRESSURE=',G15.7,
     %'  V-ERROR=',G12.4,'%   P-ERROR=',G12.4,'%'/)

20324 FORMAT('Exact SOLUTION: X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VELOCITY VX=',G15.7,' VY=',G15.7,' VZ=',G15.7,
     %'  PRESSURE=',G15.7,' V-ERROR=',G12.4,'%   P-ERROR=',G12.4,'%'/)

10134 FORMAT('Exacte SOLUTION: X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VITESSE VX=',G15.7,' VY=',G15.7,25X,
     %'  V-ERREUR=',G12.4,'%'/)

10334 FORMAT('Exacte SOLUTION: X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VITESSE VX=',G15.7,' VY=',G15.7,' VZ=',G15.7,25X,
     %'  V-ERREUR=',G12.4,'%'/)

20134 FORMAT('Exact SOLUTION: X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VELOCITY VX=',G15.7,' VY=',G15.7,26X,
     %'  V-ERROR=',G12.4,'%'/)

20334 FORMAT('Exact SOLUTION: X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VELOCITY VX=',G15.7,' VY=',G15.7,' VZ=',G15.7,25X,
     %'  V-ERROR=',G12.4,'%'/)

C     ERREUR RELATIVE MOYENNE EN TEMPS SUR LA VITESSE
      PCERVITMOY = PCERVITMOY + PCERVIT
C     ERREUR RELATIVE MOYENNE EN TEMPS SUR LA PRESSION
      PCERPREMOY = PCERPREMOY + PCERPRE

      ENDIF

C     ============================================
C     Min Max de ||Vitesse|| et PRESSION CALCULEES
C     ============================================
      VITMOY =  0D0
      VITMIN =  1D100
      VITMAX = -1D100
      PREMIN =  1D100
      PREMAX = -1D100
      PREMOY = 0D0
      NBST   = 0

      DO 90 I=1,NBNOVI
         NDL = NODDL(I-1)

C        NORME DE LA VITESSE CALCULEE AU NOEUD I
         IF( NDIM .EQ. 2 ) THEN
            VITNOR = SQRT( VXYZPN(NDL+1)**2
     %                   + VXYZPN(NDL+2)**2 )
         ELSE
            VITNOR = SQRT( VXYZPN(NDL+1)**2
     %                   + VXYZPN(NDL+2)**2
     %                   + VXYZPN(NDL+3)**2 )
         ENDIF
C        NORME MOYENNE DE LA VITESSE
         VITMOY = VITMOY + VITNOR
C        NORME MIN et MAX DE LA VITESSE CALCULEE
         IF( VITNOR .LT. VITMIN ) VITMIN = VITNOR
         IF( VITNOR .GT. VITMAX ) THEN
            VITMAX = VITNOR
            NOVMAX = I
         ENDIF

C        PRESSION CALCULEE AU NOEUD I
         IF( NODDL(I)-NDL .GT. NDIM ) THEN
            NBST    = NBST + 1
            PRECALC = VXYZPN(NODDL(I))
            PREMOY  = PREMOY + PRECALC
C           PRECALC CALCULEE MIN et MAX EN UN NOEUD
            IF( PRECALC .LT. PREMIN ) PREMIN=PRECALC
            IF( PRECALC .GT. PREMAX ) PREMAX=PRECALC
         ENDIF
 90   ENDDO

C     NORME MOYENNE DE LA VITESSE
      VITMOY = VITMOY / NBNOVI

C     PRESSION MOYENNE
      PREMOY = PREMOY / NBST

C     COORDONNEES DU NOEUD DE VITESSE MAXIMALE
      MN = MNXYZ + 3*NOVMAX
      XP = RMCN(MN)
      YP = RMCN(MN+1)
      ZP = RMCN(MN+2)

      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10090)
     %         TEMPS, VITMOY, VITMAX, NOVMAX, XP, YP, ZP,
     %         TEMPS, PREMOY, PREMIN, PREMAX, PREMAX-PREMIN
      ELSE
         WRITE(IMPRIM,20090)
     %         TEMPS, VITMOY, VITMAX, NOVMAX, XP, YP, ZP,
     %         TEMPS, PREMOY, PREMIN, PREMAX, PREMAX-PREMIN
      ENDIF

10090 FORMAT(/'Au Temps',  G13.5,' |VITESSE|Moyenne=',G14.6,
     %' |VITESSE|Max=',    G14.6,' au noeud',I9,' XYZ=',3G14.6/
     %        'Au Temps',  G13.5,' PRESSION Moyenne=',G14.6,
     %' PRESSION Min=',    G14.6,' PRESSION Max=',    G14.6,
     %' PRESSION Max-Min=',G14.6/140('-')/)

20090 FORMAT(/'At Time', G13.5,' |VELOCITY|Mean=',G14.6,
     %' |VELOCITY|Max=', G14.6,' at node',I9,' XYZ=',3G14.6/
     %        'At Time', G13.5,
     %' PRESSURE  Mean=',G14.6,' PRESSURE  Min=',G14.6,
     %' PRESSURE  Max=', G14.6,' PRESSURE Max-Min=',G14.6/140('-')/)

 1000 ENDDO


      IF( NOFOVI .GT. 0 .AND. NOFOPR .GT. 0 ) THEN
C        ERREUR RELATIVE MOYENNE EN TEMPS SUR LA VITESSE
         PCERVITMOY = PCERVITMOY / ( NCAS1-NCAS0+1)
C        ERREUR RELATIVE MOYENNE EN TEMPS SUR LA PRESSION
         PCERPREMOY = PCERPREMOY / ( NCAS1-NCAS0+1)
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10999) TIMES(NCAS0), TIMES(NCAS1),
     %                          PCERVITMOY, PCERPREMOY
         ELSE
            WRITE(IMPRIM,20999) TIMES(NCAS0), TIMES(NCAS1),
     %                          PCERVITMOY, PCERPREMOY
         ENDIF
      ENDIF

10999 FORMAT(142('+')/'Au total: du TEMPS ',G13.6,' a ',G13.6,
     %' ERREUR MOYENNE RELATIVE de la VITESSE=',G14.6,
     %'%   ERREUR MOYENNE RELATIVE de la PRESSION=',G14.6,'%'/142('+')/)

20999 FORMAT(142('+')/'Finally: From TIME ',G13.6,' to ',G13.6,
     %' VELOCITY RELATIVE MEAN ERROR=',G14.6,
     %'%    PRESSURE RELATIVE MEAN ERROR=',G14.6,'%'/142('+')/)

      RETURN
      END
