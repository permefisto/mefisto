      SUBROUTINE VOEX14( NTLXVO, LADEFI, RADEFI,
     &                   NTCUVO, MNCUVO, NTSOCU, MNSOCU, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     GENERER LE MAILLAGE D'UN VOLUME PAR EXTRUSION D'UNE SURFACE
C ----      C'EST A DIRE PAR COUCHES OU TRANCHES
C
C ENTREES :
C --------
C NTLXVO  : NUMERO DU TABLEAU TS DU LEXIQUE DU VOLUME
C LADEFI  : TABLEAU ENTIER DE DEFINITION DU VOLUME
C RADEFI  : TABLEAU REEL   DE DEFINITION DU VOLUME
C           CF '~TD/D/A_VOLUME__DEFINITION'
C
C SORTIES :
C ---------
C NTCUVO  : NUMERO      DU TMS 'NSEF' DES NUMEROS DES ELEMENTS
C MNCUVO  : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES ELEMENTS
C           CF '~td/d/a___nsef'
C NTSOCU  : NUMERO      DU TMS 'XYZSOMMET' DU VOLUME
C MNSOCU  : ADRESSE MCN DU TMS 'XYZSOMMET' DU VOLUME
C           CF '~td/d/a___xyzsommet'
C IERR    : 0 SI PAS D'ERREUR
C         > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1989
C MODIFS : PERRONNET ALAIN Laboratoire J-L. LIONS UPMC PARIS   Mars 2007
C23456+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
      IMPLICIT INTEGER (W)
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
      INTEGER           NOSOEL(64),NTLIGS(3)
      INTEGER           MNSOMS(3)
C
      IERR   = 0
      MNPLTR = 0
      MNVECT = 0
      MNZ    = 0
C
C     LE TYPE D'EXTRUSION
      NUTYCS = LADEFI( WUTYCS )
C
C     LE NUMERO DE LA SURFACE D'EXTRUSION
      NUSUCS = LADEFI( WUSUCS )
C     LE LEXIQUE DE LA SURFACE
      CALL LXNLOU( NTSURF, NUSUCS, NTLXSU, MNLXSU )
      IF( NTLXSU .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SURFACE INCONNUE'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C     LES TMS NSEF ET SOMMETS
      CALL LXTSOU( NTLXSU, 'NSEF', NTFASU, MNFASU )
      CALL LXTSOU( NTLXSU, 'XYZSOMMET'   , NTSOSU, MNSOSU )
      IF( NTFASU .LE. 0 .OR. NTSOSU .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SURFACE NON MAILLEE'
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C     LE NOMBRE DE SOMMETS DE LA SURFACE
      NBSOMS = MCN( MNSOSU+WNBSOM )
C
C     LES PARAMETRES DES NO SOMMET DE LA SURFACE
      CALL NSEFPA( MCN(MNFASU),
     %             NUTYMS, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBFASU,
     %             NX    , NY    , NZ    ,
     %             IERR   )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SURFACE DE MAILLAGE INCORRECT'
         CALL LEREUR
         IERR = 7
         RETURN
      ENDIF
C
C     RECHERCHE DES TRIEDRES DE CHAQUE SURFACE SELON NUTYCS
      IF( NUTYCS .EQ. 0 ) THEN
C
C        CHAQUE SURFACE EST DEFINIE PAR UN TRIEDRE A PARTIR DE 3 POINTS
C        ==============================================================
C        LE NOMBRE DE TRANCHES
         NBTRCS = LADEFI( WBTRCS )
         IF( NBTRCS .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'NOMBRE INCORRECT DE TRANCHES'
            CALL LEREUR
            IERR = 4
            RETURN
         ENDIF
C
C        LE TABLEAU DES ADRESSES DES 3 COORDONNEES DES 3 * NBTRCS POINTS
C        DEFINISSANT LE PLAN DES NBTRCS SURFACES FINALES DES TRANCHES
         CALL TNMCDC( 'ENTIER', 3*NBTRCS+3, MNPLTR )
         DO 10 I=1,3*NBTRCS+3
            NP = LADEFI(WPTRCS-1+I)
            CALL LXNLOU( NTPOIN, NP, NT, MN )
            IF( NT .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:10),'(I10)') NP
               KERR(1) ='POINT INCONNU' // KERR(MXLGER)(1:10)
               CALL LEREUR
               IERR = 5
               GOTO 10
            ENDIF
            CALL LXTSOU( NT, 'XYZSOMMET', NTS, MNS )
            IF( NTS .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:10),'(I10)') NP
               KERR(1) =' POINT SANS COORDONNEES' // KERR(MXLGER)(1:10)
               CALL LEREUR
               IERR = 6
               GOTO 10
            ENDIF
            MCN( MNPLTR - 1 + I ) = MNS + WYZSOM
 10      CONTINUE
C
      ELSE IF( NUTYCS .EQ. 1 ) THEN
C
C        CHAQUE SURFACE EST DEFINIE PAR LES POINTS DE 3 LIGNES
C        =====================================================
C        RECHERCHE DES 3 LIGNES
         NBTRCS = 0
         DO 20 I=0,2
C           LE NUMERO DE LA LIGNE
            NL = LADEFI(WLTRCS+I)
            IF( NL .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') NL
               KERR(1) =  'LIGNE INCORRECTE'//KERR(MXLGER)(1:4)
               CALL LEREUR
               IERR = 9
               GOTO 20
            ENDIF
C           OUVERTURE DE LA LIGNE
            CALL LXNLOU( NTLIGN, NL, NT, MN )
            IF( NT .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') NL
               KERR(1) = 'LIGNE INCONNUE' // KERR(MXLGER)(1:10)
               CALL LEREUR
               IERR = 10
               GOTO 20
            ENDIF
            NTLIGS(I+1) = NT
C           OUVERTURE DU TABLEAU SOMMETS
            CALL LXTSOU( NT, 'XYZSOMMET', NTS, MNS )
            IF( NTS .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') NL
               KERR(1) =  'LIGNE'//KERR(MXLGER)(1:4)//' SANS SOMMETS'
               CALL LEREUR
               IERR = 6
               GOTO 20
            ENDIF
C           OUVERTURE DU TABLEAU NSEF
            CALL LXTSOU( NT, 'NSEF', NTA, MNA )
            IF( NTA .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') NL
               KERR(1) =  'LIGNE' // KERR(MXLGER)(1:4)
     %                    //' SANS NSEF'
               CALL LEREUR
               IERR = 6
               GOTO 20
            ENDIF
C
C           LES PARAMETRES DES NO SOMMET DE LA LIGNE
            CALL NSEFPA( MCN(MNA),
     %                   NUTYML, NBSOLL, NBSOSL, NBTG,
     %                   LDAP,   LDNG,   LDTG,   NBARLI,
     %                   NXL  , NYL  , NZL  ,
     %                   IERR  )
            IF( IERR .NE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') NL
               KERR(1) = 'MAILLAGE INCORRECT DE LA LIGNE'
     %                 // KERR(MXLGER)(1:4)
               CALL LEREUR
               IERR = 7
               RETURN
            ENDIF
            IF( NBTRCS .LE. 0 ) THEN
               NBTRCS = NBARLI
            ELSE
               IF( NBTRCS .NE. NBARLI ) THEN
                  NBLGRC(NRERR) = 3
                  WRITE(KERR(MXLGER)(1:4),'(I4)') NBARLI
                  WRITE(KERR(MXLGER)(11:14),'(I4)') NBTRCS
                  KERR(1) = 'LES 3 LIGNES DOIVENT AVOIR'
                  KERR(2) = 'LE MEME NOMBRE D''ARETE'
                  KERR(3) = 'ICI '//KERR(MXLGER)(1:4) //' ET '
     %                     //KERR(MXLGER)(11:14)
                  CALL LEREUR
                  IERR = 11
                  GOTO 20
               ENDIF
            ENDIF
 20      CONTINUE
C
         IF( NBTRCS .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(MXLGER)(1:4),'(I4)') NBTRCS
            KERR(1) = 'NOMBRE INCORRECT DE TRANCHES='
     %                // KERR(MXLGER)(1:4)
            CALL LEREUR
            IERR = 4
         ENDIF
         IF( IERR .NE. 0 ) GOTO 9990
C
C        LE TABLEAU DES ADRESSES DES 3 COORDONNEES DES 3 * NBTRCS POINTS
C        DEFINISSANT LE PLAN DES NBTRCS SURFACES FINALES DES TRANCHES
         CALL TNMCDC( 'ENTIER', 3*NBTRCS+3, MNPLTR )
         DO 50 J=1,3
C           LA LIGNE J
            NT = NTLIGS(J)
C           OUVERTURE DU TABLEAU XYZSOMMET
            CALL LXTSOU( NT, 'XYZSOMMET', NTS, MNS )
C           OUVERTURE DU TABLEAU NSEF
            CALL LXTSOU( NT, 'NSEF', NTA, MNA )
C
C           TENTATIVE DE STRUCTURATION DE LA LIGNE
            CALL LIGSTR( NT, NTA, MNA, NTS, MNS, IERR )
            IF( IERR .NE. 0 .AND. IERR .NE. 2 ) GOTO 50
C
C           LE NOMBRE DE SOMMETS DE LA LIGNE
            NBSOML = MCN( MNS+WNBSOM )
C
C           LES PARAMETRES DES NO SOMMET DE LA LIGNE
            CALL NSEFPA( MCN(MNA),
     %                   NUTYML, NBSOLL, NBSOSL, NBTG,
     %                   LDAP,   LDNG,   LDTG,   NBARLI,
     %                   NXL  , NYL  , NZL  ,
     %                   IER  )
            MNSOMS(J) = MNS + WYZSOM
C
C           L'ADRESSE MCN DES SOMMETS DE LA LIGNE
            MN   = MNPLTR - 1 + J
            MNSO = MNSOMS(J)
            DO 30 I=1,NBARLI
               MCN( MN ) = MNSO
               MN   = MN   + 3
               MNSO = MNSO + 3
 30         CONTINUE
C
            IF( IERR .EQ. 2 ) THEN
C              LE DERNIER SOMMET EST LE PREMIER
               MCN( MN ) = MNSOMS(J)
               IERR = 0
            ELSE
C              LE DERNIER SOMMET
               MCN( MN ) = MNSO
            ENDIF
 50      CONTINUE
C
      ELSE IF( NUTYCS .GE. 2 .AND. NUTYCS .LE. 5 ) THEN
C
C        COTE MAXIMALE OU COTE DE CHAQUE TRANCHE
C        VECTEUR TRANSLATION CONSTANT OU VARIABLE PAR TRANCHE
C        LE NOMBRE DE TRANCHES
         NBTRCS = LADEFI(WBTRCS)
C        DECLARATION D'UN VECTEUR 3*NBTRCS + NBTRCS+1
         I = 4*NBTRCS+1
         CALL TNMCDC( 'REEL', I, MNVECT )
         MNZ = MNVECT + 3 * NBTRCS
         CALL AZEROR( I, RMCN(MNVECT) )
C
C        REMPLISSAGE DU TABLEAU VECTOR SELON NUTYCS
         IF( NUTYCS .EQ. 2 ) THEN
C           LA COTE MAXIMALE
            RMCN( MNVECT - 1 + 3 * NBTRCS ) = RADEFI( WMAXCS )
         ELSE IF( NUTYCS .EQ. 3 ) THEN
C           LA COTE DE CHAQUE SURFACE DES TRANCHES
            DO 52 J=1,NBTRCS
               RMCN(MNVECT-1+3*J) = RADEFI(WCTRCS-1+J)
 52         CONTINUE
         ELSE IF( NUTYCS .EQ. 4 ) THEN
C           LE VECTEUR CONSTANT DE TRANSLATION
            RMCN(MNVECT  ) = RADEFI(WRCOCS)
            RMCN(MNVECT+1) = RADEFI(WRCOCS+1)
            RMCN(MNVECT+2) = RADEFI(WRCOCS+2)
         ELSE IF( NUTYCS .EQ. 5 ) THEN
C           LE VECTEUR DE TRANSLATION DE CHAQUE SURFACE
            CALL TRTATA( RADEFI(WRVACS), RMCN(MNVECT), 3*NBTRCS )
         ENDIF
C
      ELSE IF( NUTYCS .EQ. 6 ) THEN
C
C        LE NOMBRE DE TRANCHES
         NBTRCS = LADEFI(WBTRCS)
C
      ELSE
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYCS
         KERR(1) = 'OPTION INCONNUE' // KERR(MXLGER)(1:4)
         CALL LEREUR
         IERR = 13
      ENDIF
      IF( IERR .NE. 0 ) GOTO 9990
C
C     LE VOLUME EST IL UN TORE C-A-D
C     ------------------------
C     A PRIORI CE N'EST PAS UN TORE
      LETORE = 0
C     LA DERNIERE SURFACE SE TROUVE ETRE AUSSI LA PREMIERE?
      IF( MNPLTR .GT. 0 ) THEN
         MN = MNPLTR + 3 * NBTRCS
         IF( MCN(MNPLTR)   .EQ. MCN(MN)   .AND.
     %       MCN(MNPLTR+1) .EQ. MCN(MN+1) .AND.
     %       MCN(MNPLTR+2) .EQ. MCN(MN+2) ) THEN
C           IL S'AGIT D'UN TORE
            LETORE = 1
         ENDIF
      ENDIF
C
C     GENERATION DES NSEF DU VOLUME
C     =============================
C     CONSTRUCTION DU TABLEAU 'NSEF' DE CE VOLUME
      NBCUVO = NBTRCS * NBFASU
      CALL LXTNDC(NTLXVO,'NSEF','ENTIER',WUSOEF+8*NBCUVO)
      CALL LXTSOU(NTLXVO,'NSEF',NTCUVO,MNCUVO)
C
C     MISE A JOUR DU TABLEAU 'NSEF' DE CE VOLUME
C     TYPE DE L'OBJET : VOLUME
      MCN( MNCUVO + WUTYOB ) = 4
C
C     LE TYPE INCONNU DE FERMETURE DU MAILLAGE
      MCN( MNCUVO + WUTFMA ) = -1
C
C     PAS DE TANGENTES STOCKEES
      MCN( MNCUVO + WBTGEF ) = 0
      MCN( MNCUVO + WBEFAP ) = 0
      MCN( MNCUVO + WBEFTG ) = 0
C
C     NUMERO DU TYPE DE MAILLAGE : NON STRUCTURE
      MCN( MNCUVO + WUTYMA ) = 0
C
C     NBSOEF 'NOMBRE DE SOMMETS PAR NSEF'
      MCN( MNCUVO + WBSOEF ) = 8
C
C     NBCUVO  LE NOMBRE DE CUBES DU VOLUME
      MCN( MNCUVO + WBEFOB ) = NBCUVO
C
C     LA BOUCLE SUR LES TRANCHES DU VOLUME
      LEDEC0 = 0
      LEDEC1 = NBSOMS
      MNC    = MNCUVO + WUSOEF - 1
      DO 300 NBT = 1, NBTRCS
C
C        LA BOUCLE SUR LES NSEF DU MAILLAGE DE LA SURFACE
         DO 100 N=1,NBFASU
C           LE NUMERO DES NBSOEF SOMMETS DU SOUSOBJET N
            CALL NSEFNS( N     , NUTYMS, NBSOEF, NBTGEF,
     %                   LDAPEF, LDNGEF, LDTGEF,
     %                   MNFASU, NX, NY, NZ,
     %                   NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
            IF( NCOGEL .NE. 3 .AND. NCOGEL .NE. 4 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') NCOGEL
               KERR(1) = 'SURFACE AVEC UN ELEMENT FINI DE CODE'
     %                  //KERR(MXLGER)(1:4)
               CALL LEREUR
               IERR = 8
               CALL LXTSDS( NTLXVO, 'NSEF' )
               GOTO 9990
            ENDIF
C           LE NOMBRE DE SOMMETS DE CET ELEMENT
            NBSO = NCOGEL
            DO 60 J=1,NBSO
               MCN( MNC + J        ) = NOSOEL(J) + LEDEC0
               MCN( MNC + J + NBSO ) = NOSOEL(J) + LEDEC1
 60         CONTINUE
            DO 70 J = NBSO+NBSO+1, 8
               MCN( MNC + J ) = 0
 70         CONTINUE
C           LE CUBE SUIVANT
            MNC = MNC + 8
 100     CONTINUE
C
C        LA SURFACE SUIVANTE
         LEDEC0 = LEDEC1
         IF( NBT .EQ. NBTRCS-1 .AND. LETORE .EQ. 1 ) THEN
            LEDEC1 = 0
         ELSE
            LEDEC1 = LEDEC1 + NBSOMS
         ENDIF
 300  CONTINUE
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNCUVO) )
C
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNCUVO + MOTVAR(6) ) = NONMTD ( '~>>>NSEF' )
C
C     GENERATION DES SOMMETS DU VOLUME
C     ================================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CE VOLUME
C     LE NOMBRE DE SOMMETS SELON TORE OU NON
      IF( LETORE .EQ. 0 ) THEN
C        NON TORE
         NBSOM = (NBTRCS+1) * NBSOMS
      ELSE
C        TORE
         NBSOM =  NBTRCS * NBSOMS
      ENDIF
      CALL LXTNDC( NTLXVO, 'XYZSOMMET', 'MOTS', WYZSOM+3*NBSOM )
      CALL LXTSOU( NTLXVO, 'XYZSOMMET', NTSOCU, MNSOCU )
C
C     CALCUL DES 3 COORDONNEES DES NBSOM  SOMMETS DU VOLUME
C     CALCUL DES 3 COMPOSANTES DES NBTGSV TANGENTES DU VOLUME(NUTYCS=0,1)
      NBTGSU = 0
      CALL VOE114( NUTYCS, LETORE, NBTRCS, MCN(MNPLTR),
     %             LADEFI(WONXCS),
     %             MCN(MNVECT), MCN( MNZ ),
     %             NBSOMS, RMCN(MNSOSU+WYZSOM), RMCN(MNSOCU+WYZSOM),
     %             NBTGSU, RMCN(MNSOSU+WYZSOM+3*NBSOMS),
     %             RMCN(MNSOCU+WYZSOM+3*NBSOM),
     %             IERR   )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1)='PB DE GENERATION DES SOMMETS DU VOLUME'
         CALL LEREUR
         CALL LXTSDS( NTLXVO, 'XYZSOMMET' )
         GOTO 9990
      ENDIF
C
C     NBSOM 'NOMBRE DE SOMMETS'
      MCN( MNSOCU + WNBSOM ) = NBSOM
C
C     NOMBRE DE TANGENTES
      MCN( MNSOCU + WNBTGS ) = 0
C
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSOCU + WBCOOR ) = 3
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOCU) )
C
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOCU + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     VERIFICATION DU VOLUME POSITIF DES PENTAEDRES ET HEXAEDRES
C     PERMUTATIONS DE SOMMETS POUR LE RENDRE POSITIF SINON
      CALL VOLPLUS( MNSOCU, MNCUVO )
C
C     DESTRUCTION DES TABLEAUX DEVENUS INUTILES
 9990 IF( MNPLTR .GT. 0 ) CALL TNMCDS( 'ENTIER', 3*NBTRCS+3, MNPLTR )
      IF( MNVECT .GT. 0 ) CALL TNMCDS( 'REEL'  , 4*NBTRCS+1, MNVECT )
C
      RETURN
      END
