      SUBROUTINE SUEX27( NTLXSU, LADEFI, RADEFI,
     %                   NTCUSU, MNCUSU, NTSOCU, MNSOCU, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     GENERER LE MAILLAGE DE LA SURFACE D'UN SECTEUR ANGULAIRE
C -----     D'UN TROC DE CONE OU CYLINDRE
C
C ENTREES :
C --------
C NTLXSU  : NUMERO DU TABLEAU TS DU LEXIQUE DES SURFACES
C LADEFI  : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
C RADEFI  : TABLEAU REEL   DE DEFINITION DE LA SURFACE
C           CF '~td/d/a_surface__definition'
C
C SORTIES :
C ---------
C NTCUSU  : NUMERO      DU TMS 'NSEF' DES NUMEROS DES ELEMENTS FINIS
C MNCUSU  : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES ELEMENTS FINIS
C           CF '~td/d/a___nsef'
C NTSOCU  : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOCU  : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C           CF '~td/d/a___xyzsommet'
C IERR    : 0 SI PAS D'ERREUR
C         > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: BATOUCHE_EMMEL MARS DEA D'ANALYSE NUMERIQUE      Fevrier  1991
C MODIFS: PERRONNET Alain LJLL UPMC & St Pierre du Perray  Novembre 2011
C2345X7..............................................................012
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_surface__definition.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
C
      IERR = 0
C
C     PARAMETRES DU MAILLAGE
C     ======================
      NBAHCY = LADEFI(WBAHCS)
      NBSECY = LADEFI(WBACCS)
      RAYOCS = RADEFI(WAYHAU)
      RAYOCI = RADEFI(WAYBAS)
      NUPBCY = LADEFI(WTSABA)
      NUPHCY = LADEFI(WTSAHA)
      ANGREC = RADEFI(WNGLCS)
      RAGHCY = RADEFI(WGHACS)
      RAGACY = RADEFI(WGSACS)
C
C     VERIFICATION DES PARAMETRES
C     ===========================
      IF ( ANGREC.LT. 0. )  THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:5),'(I5)') ANGREC
         KERR(1) = 'ANGLE DE RECOUVREMENT INCORRECT ='
     %              // KERR(MXLGER)(1:5)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF (NBAHCY .LT. 2 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:5),'(I5)') NBAHCY
         KERR(1) =  'NOMBRE D ARETES SUR H INCORRECT ='
     %              // KERR(MXLGER)(1:5)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF( ((NBSECY .LT. 3 ).AND.(ANGREC .GT. 180.0))
     %         .OR. (NBSECY .LE. 1) ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:5),'(I5)') NBSECY
         KERR(1) =  'NOMBRE DE SECTEURS INCORRECT ='
     %              // KERR(MXLGER)(1:5)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF( RAYOCS .LE. 0.0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15),'(G15.7)') RAYOCS
         KERR(1) =  'RAYON SUPERIEUR INCORRECT ='
     %              // KERR(MXLGER)(1:15)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF( RAYOCI .LE. 0.0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15),'(G15.7)') RAYOCI
         KERR(1) =  'RAYON INFERIEUR INCORRECT ='
     %              // KERR(MXLGER)(1:15)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF(RAGHCY.LE.0) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15),'(G15.7)') RAGHCY
         KERR(1) =  'RAGHCY INCORRECT ='
     %              // KERR(MXLGER)(1:15)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF(RAGACY.LE.0) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:15),'(G15.7)') RAGACY
         KERR(1) =  'RAGACY INCORRECT ='
     %              // KERR(MXLGER)(1:15)
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     LIMITATION A 360 DEGRES DE L'ANGLE DU SECTEUR ANGULAIRE
      ANGREC = MIN( ANGREC, 360.0 )
C
C     RECUPERATION DES 3 COORDONNEES DU CENTRE INFERIEUR DU CYLINDRE
C     ==============================================================
      CALL LXNLOU( NTPOIN, NUPBCY, NTLXSC, MN )
      CALL LXTSOU( NTLXSC, 'XYZSOMMET', NTSOSC, MNSOSC )
      MNCB = MNSOSC + WYZSOM
C
C     RECUPERATION DES 3 COORDONNEES DU CENTRE SUPERIEUR DU CYLINDRE
C     ==============================================================
      CALL LXNLOU( NTPOIN, NUPHCY, NTLXSC, MN )
      CALL LXTSOU( NTLXSC, 'XYZSOMMET', NTSOSC, MNSOSC )
      MNCH = MNSOSC + WYZSOM
C
      CALL XYZIDE( RMCN(MNCB), RMCN(MNCH), IDENTQ )
      IF( IDENTQ .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'CENTRES INFERIEUR et SUPERIEUR CONFONDUS'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     GENERATION DES XYZ DES SOMMETS
C     ==============================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET'
      NBSOM = (NBAHCY+1) * (NBSECY+1)
      IF( ANGREC .EQ. 360.0 ) THEN
          NBSOM = NBSOM-(NBAHCY+1)
      ENDIF
      CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'MOTS', WYZSOM+3*NBSOM )
      CALL LXTSOU( NTLXSU, 'XYZSOMMET', NTSOCU, MNSOCU )
C
C     GENERATION DES NUMEROS DES NOEUDS DES EF
C     ========================================
C     CONSTRUCTION DU TABLEAU 'NSEF'
      NBCUSU= NBAHCY*NBSECY
      CALL LXTNDC(NTLXSU, 'NSEF', 'ENTIER', WUSOEF+4*NBCUSU )
      CALL LXTSOU(NTLXSU, 'NSEF', NTCUSU,   MNCUSU )
C
C     GENERATION DES SOMMETS ET ELEMENTS DU CYLINDRE OU CONE
C     ======================================================
      CALL CYLCONE( NBSECY, NBAHCY, RMCN(MNCB), RMCN(MNCH), RAYOCI,
     %              RAYOCS, RAGHCY, RAGACY, ANGREC,
     %              RMCN(MNSOCU+WYZSOM), MCN(MNCUSU+WUSOEF), IERR )
      IF( IERR .NE. 0 ) THEN
C        PAS DE TMS XYZSOMMET ET NSEF
         NTCUSU=0
         MNCUSU=0
         NTSOCU=0
         MNSOCU=0
         GOTO 9999
      ENDIF
C
C     MISE A JOUR DES TABLEAUX
C     ========================
C     MISE A JOUR DU TABLEAU 'XYZSOMMET' DE LA SURFACE
C     NBSOM 'nombre de sommets'
      MCN( MNSOCU + WNBSOM ) = NBSOM
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOCU) )
C     LE NOMBRE DE TANGENTES
      MCN( MNSOCU + WNBTGS ) = 0
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSOCU + WBCOOR ) = 3
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOCU + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     MISE A JOUR DU TABLEAU 'NSEF' DE LA SURFACE
C     TYPE DE L'OBJET : SURFACE
      MCN( MNCUSU + WUTYOB ) = 3
C     LE TYPE NON-FERME DE FERMETURE DU MAILLAGE
      MCN( MNCUSU + WUTFMA ) = 0
C     PAS DE TANGENTES STOCKEES
      MCN( MNCUSU + WBTGEF ) = 0
      MCN( MNCUSU + WBEFAP ) = 0
      MCN( MNCUSU + WBEFTG ) = 0
C     NUMERO DU TYPE DE MAILLAGE : NON STRUCTURE
      MCN( MNCUSU + WUTYMA ) = 0
C     NBSOEF 'nombre de sommets par sous-objets'
      MCN( MNCUSU + WBSOEF ) = 4
C     NBEFOB 'nombre de sous-objets de l''objet'
      MCN( MNCUSU + WBEFOB ) = NBCUSU
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNCUSU) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNCUSU + MOTVAR(6) ) = NONMTD ( '~>>>NSEF' )
C
 9999 RETURN
      END
