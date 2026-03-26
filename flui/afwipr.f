      SUBROUTINE AFWIPR( NUTYEL, NDIM,   NBNOVI, MNXYZN, NODDL,
     %                   NTDLVP, NCAS,   NBVPAF, WITPRE,
     %                   VITMIN, VITMAX, VITMOY,
     %                   PREMIN, PREMAX, PREMOY )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LES VITESSES WIT ET LES PRESSIONS CALCULEES PRE
C -----    ET LES PRESSIONS CALCULEES DU CAS NCAS PARMI LES CAS CALCULES
C
C ENTREES:
C --------
C TEMPS  : TEMPS DU VECTEUR VITESSE+PRESSION dans ./incl/ctemps.inc
C NUTYEL : 13 TRIANGLE  BREZZI FORTIN 2D
C          15 TRIANGLE  TAYLOR-HOOD   2D
C          19 TETRAEDRE BREZZI FORTIN 3D
C          20 TETRAEDRE TAYLOR-HOOD   3D
C NDIM   : DIMENSION 2 OU 3 DE L'ESPACE DE L'OBJET
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
C          SOMMETS + BARYCENTRES        des TETRAEDRES POUR BREZZI-FORTIN
C          SOMMETS + MILIEUX DES ARETES des TETRAEDRES POUR TAYLOR-HOOD
C MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD' DU MAILLAGE DE L'OBJET
C NODDL  : TABLEAU DU NUMERO DU DERNIER D.L. DE CHAQUE NOEUD
C NTDLVP   : NOMBRE TOTAL DE DEGRES DE LIBERTES EN VITESSES PRESSIONS
C NCAS   : NUMERO DE LA CARTE DES VITESSES_PRESSIONS A AFFICHER
C NBVPAF : NOMBRE DE NOEUDS DE VITESSE-PRESSION A AFFICHER
C WITPRE : LES CARTES DES VITESSES-PRESSIONS AUX NOEUDS DU MAILLAGE
C
C SORTIES:
C --------
C VITMIN : MINIMUM DE LA NORME DE LA VITESSE CALCULEE EN UN NOEUD
C VITMAX : MAXIMUM DE LA NORME DE LA VITESSE CALCULEE EN UN NOEUD
C VITMOY : MOYENNE DE LA NORME DE LA VITESSE CALCULEE EN TOUS LES NOEUDS
C PREMIN : MINIMUM DE LA PRESSION CALCULEE EN UN NOEUD DU MAILLAGE
C PREMAX : MAXIMUM DE LA PRESSION CALCULEE EN UN NOEUD DU MAILLAGE
C PREMOY : MOYENNE DE LA PRESSION CALCULEE EN TOUS LES SOMMETS DU MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Juillet 2011
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
C
      INTEGER           NODDL(0:NBNOVI)
      DOUBLE PRECISION  WITPRE(NTDLVP,1:*)
      DOUBLE PRECISION  VITNOR, VITMIN, VITMAX, VITMOY,
     %                  PREMIN, PREMAX, PREMOY
C
C     CALCUL DES ERREURS AVEC UNE SOLUTION EXACTE TAYLOR-GREEN VERTEX
      DOUBLE PRECISION  XP, YP, ZP, PARAMF(5),
     %                  XWelocity, YWelocity, ZWelocity, Pressure,
     %                  VEXNORM, VEXAMX,   VERROR,  VERRELA, VERRLMX,
     %                  VERRMX,  VEXASOM,  VERRSOM,
     %                  PRERRELA,PRERRLMX, PRERRABS, PRERRAMX, PRECALC,
     %                  PREXSOM, PRERRSOM, PREXMX,   PREXMI,
     %                  PRCAMI,  PRCAMX
      INTRINSIC         SQRT
C
C     RECHERCHE DES FONCTIONS UTILISATEUR DONNANT LA
C     WITESSE_EXACTE(t,x,y,z,nocomp) ou EXACT_WELOCITY(t,x,y,z,nocomp)
      NOFOVI = NOFOWITE()
C
C     RECHERCHE DE LA FONCTION UTILISATEUR DONNANT LA
C     PRESSION_EXACTE(t,x,y,z) ou EXACT_PRESSURE(t,x,y,z)
      NOFOPR = NOFOPRES()
      PRERRELA = 0D0
C
C     L'ADRESSE MCN DES COORDONNEES DES NOEUDS VITESSE
C     ICI LE TMS XYZNOEUD A ETE ENRICHI DES BARYCENTRES DES EF
C     DANS LES CAS DES EF DE BREZZI-FORTIN
C     MN = MNXYZN + WYZNOE - 3
C
      IF( NOFOVI .LE. 0 .OR. NOFOPR .LE. 0 ) THEN
C
C     =============================================================
C     AFFICHAGE DES COMPOSANTES CALCULEES DE LA WITESSE ET PRESSION
C     NON CONNAISSANCE DE LA VITESSE EXACTE OU PRESSION EXACTE
C     =============================================================
ccc      NOEUD1 = NBNOVI/2+1
ccc      NOEUD2 = NOEUD1 + NBVPAF
ccc      NOEUD2 = MIN( NOEUD2, NBNOVI )
ccc
ccc      NBVP   = MIN( NBNOVI, NBVPAF )
ccc      NOEUD1 = NBNOVI - NBVP + 1
ccc      NOEUD2 = NBNOVI
      NOEUD1 = NBNOVI/2
      NOEUD2 = NOEUD1 + 10
      NOEUD2 = MIN( NOEUD2, NBNOVI )
C
      IF( LANGAG .EQ. 0 ) THEN
C
C        AFFICHAGE EN FRANCAIS
         WRITE(IMPRIM,10019) TEMPS,NOEUD1,NOEUD2,NBNOVI,NTDLVP
10019    FORMAT(/'Au Temps ',G15.6,': les WITESSES et PRESSIONS des',
     %    'NOEUDS',I8,' a ',I8,'/',I9,' NOEUDS (Au total',I9,' DL):')
C
         DO 20 I = NOEUD1, NOEUD2
            MN  = MNXYZN + WYZNOE - 3 + 3*I
            NDL = NODDL(I-1)
            ND  = NODDL(I) - NDL
            IF( NUTYEL .EQ. 19 .OR.  NUTYEL .EQ. 20 ) THEN
C
C              3D BF ou TH : WITESSE"PRESSION
               IF( ND .GT. NDIM ) THEN
                  WRITE(IMPRIM,10023) I, (RMCN(MN+K),K=0,2),
     %                      (WITPRE(NDL+K,NCAS),K=1,ND)
               ELSE
                  WRITE(IMPRIM,10033) I, (RMCN(MN+K),K=0,2),
     %                      (WITPRE(NDL+K,NCAS),K=1,ND)
               ENDIF
c
            ELSE
C
C              2D: NOMBRE DE DL EN LE NOEUD I
               IF( ND .GT. NDIM ) THEN
                  WRITE(IMPRIM,10022) I, (RMCN(MN+K),K=0,2),
     %                 (WITPRE(NDL+K,NCAS),K=1,ND)
               ELSE
                  WRITE(IMPRIM,10032) I, (RMCN(MN+K),K=0,2),
     %                 (WITPRE(NDL+K,NCAS),K=1,ND)
               ENDIF
            ENDIF
 20      CONTINUE
C
10022 FORMAT('Noeud',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WITESSE VX=',G15.7,' WY=',G15.7,' PRESSION=',G15.7)
C
10032 FORMAT('Noeud',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WITESSE VX=',G15.7,' WY=',G15.7)
C
10023 FORMAT('Noeud',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WITESSE VX=',G15.7,' WY=',G15.7,' WZ=',G15.7,
     %'  PRESSION=',G15.7)
C
10033 FORMAT('Noeud',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WITESSE VX=',G15.7,' WY=',G15.7, ' WZ=',G15.7 )
C
      ELSE
C
C        ENGLISH PRINTING
         WRITE(IMPRIM,20019) TEMPS,NOEUD1,NOEUD2,NBNOVI,NTDLVP
20019 FORMAT(/'At Time ',G15.6,': the VELOCITIES and PRESSURES of nodes'
     %,I8,' to ',I8,'/',I8,' NODES (Total:',I9,' DoF):')
C
         DO 30 I = NOEUD1, NOEUD2
            MN  = MNXYZN + WYZNOE - 3 + 3*I
            NDL = NODDL(I-1)
            ND  = NODDL(I) - NDL
C
C           AFFICHAGE DES COMPOSANTES DE LA WITESSE ET
C           EVENTUELLEMENT DE LA PRESSION
            IF( NDIM .EQ. 3 ) THEN
C
C              ELEMENT FINI 3D
               IF( ND .GT. NDIM ) THEN
C                 NOEUD AVEC PRESSION
                  WRITE(IMPRIM,20023) I, (RMCN(MN+K),K=0,2),
     %                      (WITPRE(NDL+K,NCAS),K=1,ND)
               ELSE
C                 NOEUD SANS PRESSION
                  WRITE(IMPRIM,20033) I, (RMCN(MN+K),K=0,2),
     %                      (WITPRE(NDL+K,NCAS),K=1,ND)
               ENDIF
c
            ELSE
C
C              ELEMENT FINI 2D
               IF( ND .GT. NDIM ) THEN
C                 NOEUD AVEC PRESSION
                  WRITE(IMPRIM,20022) I, (RMCN(MN+K),K=0,2),
     %                 (WITPRE(NDL+K,NCAS),K=1,ND)
               ELSE
C                 NOEUD SANS PRESSION
                  WRITE(IMPRIM,20032) I, (RMCN(MN+K),K=0,2),
     %                 (WITPRE(NDL+K,NCAS),K=1,ND)
               ENDIF
            ENDIF
 30      CONTINUE
C
      ENDIF
C
20022 FORMAT('Node',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WELOCITY VX=',G15.7,' WY=',G15.7,'  PRESSURE=',G15.7)
C
20032 FORMAT('Node',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WELOCITY VX=',G15.7,' WY=',G15.7)
C
20023 FORMAT('Node',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WELOCITY VX=',G15.7,' WY=',G15.7,' WZ=',G15.7,
     %'  PRESSURE=',G15.7)
C
20033 FORMAT('Node',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WELOCITY VX=',G15.7,' WY=',G15.7,' WZ=',G15.7)
C
      ELSE
C
C     ==================================================================
C     AFFICHAGE DES ERREURS SUR LA WITESSE ET LA PRESSION
C     ICI SONT SUPPOSEES CONNUES ET APPELABLES LES FONCTIONS UTILISATEUR
C     VITESSE_EXACTE(t,x,y,z,nocomp) ou EXACT_VELOCITY(t,x,y,z,nocomp)
C     PRESSION_EXACTE(t,x,y,z) ou EXACT_PRESSURE(t,x,y,z)
C     ==================================================================
      NBVP = MIN( NBNOVI, NBVPAF )
      IF( LANGAG .EQ. 0 ) THEN
C        AFFICHAGE EN FRANCAIS
         WRITE(IMPRIM,10110) TEMPS,NBVP,NBNOVI,NTDLVP
10110    FORMAT(/'Au TEMPS ',G15.6,': les ERREURS sur les WITESSES et PR
     %ESSIONS de ',I8,'/',I8,' NOEUDS (Au total',I8,' DL):'/108(1H-))
      ELSE
C        ENGLISH PRINTING
         WRITE(IMPRIM,20110) TEMPS,NBVP,NBNOVI,NTDLVP
20110    FORMAT(/'At TIME ',G15.6,': the ERRORS on VELOCITIES and PRESSU
     %RES of ',I8,'/',I8,' NODES (Total:',I8,' DoF):')
      ENDIF
C
C     INITIALISATIONS POUR LE CALCUL DES ERREURS
      VEXASOM = 0D0
      VERRSOM = 0D0
      VERRMX  =-1D100
      VEXAMX  = 0D0
      VERRLMX =-1D100
      IVERMX  = 0
C
      PREXSOM  = 0D0
      PRCAMI   = 1D100
      PRCAMX   =-1D100
      PREXMI   = 1D100
      PREXMX   =-1D100
      PRERRLMX =-1D100
      PRERRAMX = 0D0
      PRERRSOM = 0D0
      IPMAX    = 0
C
      DO 120 I=1,NBNOVI
C
C        ADRESSE DES COORDONNEES DU NOEUD I
         MN = MNXYZN + WYZNOE - 3 + 3*I
         XP = RMCN(MN)
         YP = RMCN(MN+1)
         ZP = RMCN(MN+2)
C
C        NOMBRE DE DL AU NOEUD I
         NDL = NODDL(I-1)
         ND  = NODDL(I) - NDL
C
         IF( ND .GT. NDIM ) THEN
C
C           CALCUL DE LA PRESSION_EXACTE(temps,xp,yp,zp) AU NOEUD I
            PARAMF(1) = TEMPS
            PARAMF(2) = XP
            PARAMF(3) = YP
            PARAMF(4) = ZP
            CALL FONVAL( NOFOPR, 4, PARAMF, NCODEV, Pressure )
            IF( NCODEV .LE. 0 ) RETURN
C
C           CALCUL DES ERREURS LOCALES SUR LA PRESSION AU NOEUD I
            PREXMI  = MIN( PREXMI, Pressure )
            PREXMX  = MAX( PREXMX, Pressure )
            PREXSOM = PREXSOM + ABS( Pressure )
C
            PRECALC = WITPRE(NODDL(I),NCAS)
            PRCAMI  = MIN( PRCAMI, PRECALC )
            PRCAMX  = MAX( PRCAMX, PRECALC )
C
            PRERRABS = ABS( PRECALC - Pressure )
            PRERRSOM = PRERRSOM + PRERRABS
            PRERRAMX = MAX( PRERRAMX, PRERRABS )
C
            PRERRELA = ABS( PRERRABS/Pressure )
            IF( ABS(Pressure) .GT. 1D-6 .AND.
     %         PRERRELA .GT. PRERRLMX ) THEN
               PRERRLMX = PRERRELA
               IPMAX = I
            ENDIF
         ENDIF
C
C        CALCUL DE LA VITESSE EXACTE(temps,xp,yp,zp,NoComposante)
         PARAMF(1) = TEMPS
         PARAMF(2) = XP
         PARAMF(3) = YP
         PARAMF(4) = ZP
         PARAMF(5) = 1
         CALL FONVAL( NOFOVI, 5, PARAMF, NCODEV, XWelocity )
         IF( NCODEV .LE. 0 ) RETURN
C
         PARAMF(5) = 2
         CALL FONVAL( NOFOVI, 5, PARAMF, NCODEV, YWelocity )
         IF( NCODEV .LE. 0 ) RETURN
C
         IF( NDIM .GE. 3 ) THEN
            PARAMF(5) = 3
            CALL FONVAL( NOFOVI, 5, PARAMF, NCODEV, ZWelocity )
            IF( NCODEV .LE. 0 ) RETURN
         ELSE
            ZWelocity = 0D0
         ENDIF
C
C        CALCUL DES ERREURS LOCALES SUR LA VITESSE
         IF( NDIM .GE. 3 ) THEN
            VEXNORM= SQRT(  XWelocity**2 + YWelocity**2 + ZWelocity**2 )
            VERROR = SQRT( (XWelocity-WITPRE(NDL+1,NCAS))**2
     %                   + (YWelocity-WITPRE(NDL+2,NCAS))**2
     %                   + (ZWelocity-WITPRE(NDL+3,NCAS))**2 )
         ELSE
            VEXNORM = SQRT(  XWelocity**2 + YWelocity**2 )
            VERROR  = SQRT( (XWelocity-WITPRE(NDL+1,NCAS))**2
     %                    + (YWelocity-WITPRE(NDL+2,NCAS))**2 )
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
         VERRLMX = MAX( VERRLMX, VERRELA )
C
C        LIMITATION DU NOMBRE DE LIGNES AFFICHEES
C        ----------------------------------------
ccc         IF( I .GE. NBNOVI/2 .and. I .LE. NBNOVI/2+NBVP ) THEN
         IF( I .GE. NBNOVI-15 .and. I .LE. NBNOVI ) THEN
C
         IF( NDIM .EQ. 3 ) THEN
C
C           ELEMENT FINI 3D
            IF( ND .GT. NDIM ) THEN
C              NOEUD AVEC PRESSION
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10023) I, (RMCN(MN+K),K=0,2),
     %                      (WITPRE(NDL+K,NCAS),K=1,ND)
                  WRITE(IMPRIM,10324)    (RMCN(MN+K),K=0,2),
     %                  XWelocity, YWelocity, ZWelocity,
     %                  Pressure, VERRELA*100, PRERRELA*100
               ELSE
                  WRITE(IMPRIM,20023) I, (RMCN(MN+K),K=0,2),
     %                      (WITPRE(NDL+K,NCAS),K=1,ND)
                  WRITE(IMPRIM,20324)    (RMCN(MN+K),K=0,2),
     %                  XWelocity, YWelocity, ZWelocity,
     %                  Pressure, VERRELA*100, PRERRELA*100
               ENDIF
            ELSE
C              NOEUD SANS PRESSION
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10033) I, (RMCN(MN+K),K=0,2),
     %                      (WITPRE(NDL+K,NCAS),K=1,ND)
                  WRITE(IMPRIM,10334)    (RMCN(MN+K),K=0,2),
     %                  XWelocity, YWelocity, ZWelocity, VERRELA*100
               ELSE
                  WRITE(IMPRIM,20033) I, (RMCN(MN+K),K=0,2),
     %                      (WITPRE(NDL+K,NCAS),K=1,ND)
                  WRITE(IMPRIM,20334)    (RMCN(MN+K),K=0,2),
     %                  XWelocity, YWelocity, ZWelocity, VERRELA*100
               ENDIF
            ENDIF
c
         ELSE
C
C           ELEMENT FINI 2D
            IF( ND .GT. NDIM ) THEN
C              NOEUD AVEC PRESSION
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10022) I, (RMCN(MN+K),K=0,2),
     %                 (WITPRE(NDL+K,NCAS),K=1,ND)
                  WRITE(IMPRIM,10124)    (RMCN(MN+K),K=0,2),
     %                  XWelocity, YWelocity,
     %                  Pressure, VERRELA*100, PRERRELA*100
               ELSE
                  WRITE(IMPRIM,20022) I, (RMCN(MN+K),K=0,2),
     %                 (WITPRE(NDL+K,NCAS),K=1,ND)
                  WRITE(IMPRIM,20124)    (RMCN(MN+K),K=0,2),
     %                  XWelocity, YWelocity,
     %                  Pressure, VERRELA*100, PRERRELA*100
               ENDIF
            ELSE
C              NOEUD SANS PRESSION
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10032) I, (RMCN(MN+K),K=0,2),
     %                 (WITPRE(NDL+K,NCAS),K=1,ND)
                  WRITE(IMPRIM,10134)    (RMCN(MN+K),K=0,2),
     %                  XWelocity, YWelocity, VERRELA*100
               ELSE
                  WRITE(IMPRIM,20032) I, (RMCN(MN+K),K=0,2),
     %                 (WITPRE(NDL+K,NCAS),K=1,ND)
                  WRITE(IMPRIM,20134)    (RMCN(MN+K),K=0,2),
     %                  XWelocity, YWelocity, VERRELA*100
               ENDIF
            ENDIF
         ENDIF
         ENDIF
C
 120  CONTINUE
C
      WRITE(IMPRIM,20501) VEXAMX
20501 FORMAT('MAX NODE EXACT WELOCITY NORM =',G14.6)
      WRITE(IMPRIM,20502) VERRMX, VERRMX/VEXAMX*100, VERRLMX,
     %                    VERRSOM/VEXASOM*100
20502 FORMAT('MAX ABSOLUTE   WELOCITY ERROR=',G14.6,
     %      ' MAX RELATIVE WELOCITY ERROR=',G14.6,'%',
     %     '  MAX NODE RELATIVE WELOCITY ERROR=',G14.6,'%',
     %     '  MEAN RELATIVE WELOCITY ERROR=',G14.6,'%')
C
C     ADRESSE DES COORDONNEES DU NOEUD IVERMX
      MN = MNXYZN + WYZNOE - 3 + 3*IVERMX
      XP = RMCN(MN)
      YP = RMCN(MN+1)
      ZP = RMCN(MN+2)
C
C     CALCUL DE LA VITESSE EXACTE(temps,xp,yp,zp,NoComposante)
      PARAMF(1) = TEMPS
      PARAMF(2) = XP
      PARAMF(3) = YP
      PARAMF(4) = ZP
      PARAMF(5) = 1
      CALL FONVAL( NOFOVI, 5, PARAMF, NCODEV, XWelocity )
      IF( NCODEV .LE. 0 ) RETURN
C
      PARAMF(5) = 2
      CALL FONVAL( NOFOVI, 5, PARAMF, NCODEV, YWelocity )
      IF( NCODEV .LE. 0 ) RETURN
C
      IF( NDIM .GE. 3 ) THEN
         PARAMF(5) = 3
         CALL FONVAL( NOFOVI, 5, PARAMF, NCODEV, ZWelocity )
         IF( NCODEV .LE. 0 ) RETURN
      ELSE
         ZWelocity = 0D0
      ENDIF
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10498) VEXAMX, VERRMX, IVERMX, XP,YP,ZP,
     %                       XWelocity, YWelocity, ZWelocity,
     %                      (WITPRE(NODDL(IVERMX-1)+K,NCAS),K=1,NDIM)
      ELSE
         WRITE(IMPRIM,20498) VEXAMX, VERRMX, IVERMX, XP,YP,ZP,
     %                       XWelocity, YWelocity, ZWelocity,
     %                      (WITPRE(NODDL(IVERMX-1)+K,NCAS),K=1,NDIM)
      ENDIF
C
10498 FORMAT('MAX NORME WITESSE AUX NOEUDS=',G14.6,
     %'  MAX ERREUR ABSOLUE WITESSE EN UN NOEUD=',G14.6,
     %' AU NOEUD',i7,' X=',G14.6,' Y=',G14.6,' Z=',G14.6/
     %'Wx Exact=',G14.6,' Wy Exact=',G14.6,' Wz Exact=',G14.6/
     %'Wx Calcu=',G14.6,' Wy Calcu=',G14.6,' Wz Calcu=',G14.6)
20498 FORMAT('MAX MESH NODE  WELOCITY NORM =',G14.6,
     %' MAX ABSOLUTE NODE WELOCITY ERROR=',G14.6,
     %' at NODE',i7,' X=',G14.6,' Y=',G14.6,' Z=',G14.6/
     %'Wx Exact=',G14.6,' Wy Exact=',G14.6,' Wz Exact=',G14.6/
     %'Wx Compu=',G14.6,' Wy Compu=',G14.6,' Wz Compu=',G14.6)
C
      WRITE(IMPRIM,20505) PREXMI, PREXMX, PREXMX-PREXMI
20505 FORMAT(/'MIN EXACT     PRESSURE=',G14.6,
     %      '  MAX EXACT     PRESSURE=',G14.6,
     %      '  MAX-MIN EXACT PRESSURE=',G14.6)
      WRITE(IMPRIM,20515) PRCAMI, PRCAMX, PRCAMX-PRCAMI
20515 FORMAT( 'MIN COMPUTED  PRESSURE=',G14.6,
     %      '  MAX COMPUTED  PRESSURE=',G14.6,
     %      '  MAX-MIN COMPU PRESSURE=',G14.6)
      WRITE(IMPRIM,20503) PRERRAMX, PRERRAMX/(PREXMX-PREXMI)*100,
     %            PRERRSOM/PREXSOM*100
20503 FORMAT( 'MAX ABS PRESSURE ERROR=',G14.6,
     %      '  REL MAX PRESSURE ERROR=',G14.6,'%',
     %       ' REL MEANPRESSURE ERROR=',G14.6,'%')
C
10124 FORMAT('EXACT SOL : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WITESSE VX=',G15.7,' WY=',G15.7,'  PRESSION=',G15.7,
     %'  V-ERREUR=',G12.4,'%   P-ERREUR=',G12.4,'%'/)
C
10324 FORMAT('EXACT SOL : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WITESSE VX=',G15.7,' WY=',G15.7,' WZ=',G15.7,
     %'  PRESSION=',G15.7,
     %'  V-ERREUR=',G12.4,'%   P-ERREUR=',G12.4,'%'/)
C
20124 FORMAT('EXACT SOL : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WELOCITY VX=',G15.7,' WY=',G15.7,'  PRESSURE=',G15.7,
     %'  V-ERROR=',G12.4,'%   P-ERROR=',G12.4,'%'/)
C
20324 FORMAT('EXACT SOL : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WELOCITY VX=',G15.7,' WY=',G15.7,' WZ=',G15.7,
     %'  PRESSURE=',G15.7,' V-ERROR=',G12.4,'%   P-ERROR=',G12.4,'%'/)
C
10134 FORMAT('EXACT SOL : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WITESSE VX=',G15.7,' WY=',G15.7,25X,
     %'  V-ERREUR=',G12.4,'%'/)
C
10334 FORMAT('EXACT SOL : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WITESSE VX=',G15.7,' WY=',G15.7,' WZ=',G15.7,25X,
     %'  V-ERREUR=',G12.4,'%'/)
C
20134 FORMAT('EXACT SOL : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WELOCITY VX=',G15.7,' WY=',G15.7,26X,
     %'  V-ERROR=',G12.4,'%'/)
C
20334 FORMAT('EXACT SOL : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  WELOCITY VX=',G15.7,' WY=',G15.7,' WZ=',G15.7,25X,
     %'  V-ERROR=',G12.4,'%'/)
C
C     CALCUL DE LA PRESSION_EXACTE(temps,xp,yp,zp) AU NOEUD IPMAX
      IF( IPMAX .GT. 0 ) THEN
C
C        ADRESSE DES COORDONNEES DU NOEUD IPMAX
         MN = MNXYZN + WYZNOE - 3 + 3*IPMAX
         XP = RMCN(MN)
         YP = RMCN(MN+1)
         ZP = RMCN(MN+2)
C
C        PRESSION EXACTE AU NOEUD IPMAX
         PARAMF(1) = TEMPS
         PARAMF(2) = XP
         PARAMF(3) = YP
         PARAMF(4) = ZP
         CALL FONVAL( NOFOPR, 4, PARAMF, NCODEV, Pressure )
C
C        PRESSION CALCULEE AU NOEUD IPMAX
         PRECALC = WITPRE(NODDL(IPMAX),NCAS)
C
C        CALCUL DES ERREURS LOCALES SUR LA PRESSION AU NOEUD IPMAX
         PRERRABS = ABS( PRECALC - Pressure )
         PRERRLMX = ABS( PRERRABS / Pressure )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10698) Pressure, PRECALC, PRERRLMX*100,
     %                          IPMAX, XP,YP,ZP
         ELSE
            WRITE(IMPRIM,20698) Pressure, PRECALC, PRERRLMX*100,
     %                          IPMAX, XP,YP,ZP
         ENDIF
C
      ENDIF
C
10698 FORMAT('PRESSION EXACTE=',G14.6,' PRESSION CALCULEE=',G14.6,
     %' MAX ERREUR RELATIVE=',G12.4,'%',
     %' au NOEUD',I8,' X=',G14.6,' Y=',G14.6,' Z=',G14.6)
20698 FORMAT('EXACT PRESSURE=',G14.6,' COMPUTED PRESSURE=',G14.6,
     %' RELATIVE ERROR MAX=',G12.4,'%',
     %' at NODE',I8,' X=',G14.6,' Y=',G14.6,' Z=',G14.6)
C
      ENDIF
C
C     ==========================================
C     MIN MAX de |Witesse| et PRESSION CALCULEES
C     ==========================================
      VITMOY =  0D0
      VITMIN =  1D100
      VITMAX = -1D100
      PREMIN =  1D100
      PREMAX = -1D100
      PREMOY = 0D0
      NBST   = 0
C
      DO 90 I=1,NBNOVI
         NDL = NODDL(I-1)
C
C        NORME DE LA WITESSE CALCULEE AU NOEUD I
         IF( NDIM .EQ. 2 ) THEN
            VITNOR = SQRT( WITPRE(NDL+1,NCAS)**2
     %                   + WITPRE(NDL+2,NCAS)**2 )
         ELSE
            VITNOR = SQRT( WITPRE(NDL+1,NCAS)**2
     %                   + WITPRE(NDL+2,NCAS)**2
     %                   + WITPRE(NDL+3,NCAS)**2 )
         ENDIF
C        NORME MOYENNE DE LA WITESSE
         VITMOY = VITMOY + VITNOR
C        NORME MIN et MAX DE LA WITESSE CALCULEE
         IF( VITNOR .LT. VITMIN ) VITMIN = VITNOR
         IF( VITNOR .GT. VITMAX ) THEN
            VITMAX = VITNOR
            NOVMAX = I
         ENDIF
C
C        PRESSION CALCULEE AU NOEUD I
         IF( NODDL(I)-NDL .GT. NDIM ) THEN
            NBST    = NBST + 1
            PRECALC = WITPRE(NODDL(I),NCAS)
            PREMOY  = PREMOY + PRECALC
C           PRECALC CALCULEE MIN et MAX EN UN NOEUD
            IF( PRECALC .LT. PREMIN ) PREMIN=PRECALC
            IF( PRECALC .GT. PREMAX ) PREMAX=PRECALC
         ENDIF
 90   CONTINUE
C
C     NORME MOYENNE DE LA WITESSE
      VITMOY = VITMOY / NBNOVI
C
C     PRESSION MOYENNE
      PREMOY = PREMOY / NBST
C
C     COORDONNEES DU NOEUD DE WITESSE MAXIMALE
      MN = MNXYZN + WYZNOE - 3 + 3*NOVMAX
      XP = RMCN(MN)
      YP = RMCN(MN+1)
      ZP = RMCN(MN+2)
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10090)
     %         TEMPS, VITMOY, VITMAX, NOVMAX, XP, YP, ZP,
     %         TEMPS, PREMOY, PREMIN, PREMAX, PREMAX-PREMIN
      ELSE
         WRITE(IMPRIM,20090)
     %         TEMPS, VITMOY, VITMAX, NOVMAX, XP, YP, ZP,
     %         TEMPS, PREMOY, PREMIN, PREMAX, PREMAX-PREMIN
      ENDIF
C
10090 FORMAT(/'Au Temps',  G13.5,
     %' |WITESSE|Moyenne=',G14.6,
     %' |WITESSE|Maxi=',   G14.6,' au noeud',I9,' XYZ=',3G14.6/
     %        'Au Temps',  G13.5,
     %' PRESSION Moyenne=',G14.6,' PRESSION Mini=',G14.6,
     %' PRESSION Maxi=',G14.6,' PRESSION Max-Min=',G14.6 )
C
20090 FORMAT(/'At Time', G13.5,
     %' |WELOCITY|Mean=',G14.6,
     %' |WELOCITY|Maxi=',G14.6,' at node',I9,' XYZ=',3G14.6/
     %        'At Time', G13.5,
     %' PRESSURE  Mean=',G14.6,' PRESSURE  Mini=',G14.6,
     %' PRESSURE  Maxi=',G14.6,' PRESSURE Max-Min=',G14.6 )
C
      RETURN
      END
