      SUBROUTINE SUEX14( NTLXSU, LADEFI, RADEFI,
     &                   NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     GENERER LE MAILLAGE D'UNE SURFACE PAR EXTRUSION D'UNE LIGNE
C ----      STRUCTUREE C'EST A DIRE PAR COUCHES OU TRANCHES
C
C ENTREES :
C --------
C NTLXSU  : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE
C LADEFI  : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
C RADEFI  : TABLEAU REEL   DE DEFINITION DE LA SURFACE
C           CF '~/TD/D/A_SURFACE__DEFINITION'
C
C SORTIES :
C ---------
C NTFASU  : NUMERO      DU TMS 'NSEF' DES NUMEROS DES ELEMENTS
C MNFASU  : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES ELEMENTS
C           CF '~/TD/D/A___NSEF'
C NTSOFA  : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOFA  : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C           CF '~/TD/D/A___XYZSOMMET'
C IERR    : 0 SI PAS D'ERREUR
C         > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR  : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   NOVEMBRE 1996
C23456+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
      IMPLICIT INTEGER (W)
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
      INTEGER           NOSOEL(4),NTLIGS(3)
      INTEGER           MNSOMS(3)
C
      IERR   = 0
      MNPLTR = 0
      MNVECT = 0
      MNZ    = 0
      NBTGSU = 0
C
C     LE TYPE D'EXTRUSION
      NUTYCL = LADEFI( WUTYCL )
C
C     LE NUMERO DE LA LIGNE D'EXTRUSION
      NULICL = LADEFI( WULICL )
C     LE LEXIQUE DE LA LIGNE
      CALL LXNLOU( NTLIGN, NULICL, NTLXLI, MNLXLI )
      IF( NTLXLI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE A DEPLACER INCONNUE'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C     LES TMS NSEF ET SOMMETS DE LA LIGNE D'EXTRUSION
      CALL LXTSOU( NTLXLI, 'NSEF', NTARLI, MNARLI )
      CALL LXTSOU( NTLXLI, 'XYZSOMMET', NTSOLI, MNSOLI )
      IF( NTARLI .LE. 0 .OR. NTSOLI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE A DEPLACER NON MAILLEE'
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C     LE NOMBRE DE SOMMETS DE LA LIGNE A DEPLACER
      NBSOML = MCN( MNSOLI + WNBSOM )
C     LE NOMBRE DE TANGENTES DE LA LIGNE
      NBTGSL = MCN( MNSOLI + WNBTGS )
C
C     LES PARAMETRES DES NO SOMMET DE LA LIGNE A DEPLACER
      CALL NSEFPA( MCN(MNARLI),
     %             NUTYML, NBSOEL, NBSOEF, NBTGEL,
     %             LDAPEF, LDNGEF, LDTGEF, NBARLI,
     %             NX   , NY   , NZ   ,
     %             IERR   )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE A DEPLACER DE MAILLAGE INCORRECT'
         CALL LEREUR
         IERR = 7
         RETURN
      ENDIF
C     LIGNE FERMEE
      NUTFML = MCN( MNARLI + WUTFMA )
C
C     RECHERCHE DES TRIEDRES DE CHAQUE LIGNE SELON NUTYCL
      IF( NUTYCL .EQ. 0 ) THEN
C
C        CHAQUE LIGNE EST DEFINIE PAR UN TRIEDRE A PARTIR DE 3 POINTS
C        ==============================================================
C        LE NOMBRE DE TRANCHES
         NBTRCL = LADEFI( WBTRCL )
         IF( NBTRCL .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'NOMBRE INCORRECT DE TRANCHES'
            CALL LEREUR
            IERR = 4
            RETURN
         ENDIF
C
C        LE TABLEAU DES ADRESSES DES 3 COORDONNEES DES 3 * NBTRCL POINTS
C        DEFINISSANT LE PLAN DES NBTRCL LIGNES FINALES DES TRANCHES
         CALL TNMCDC( 'ENTIER', 3*NBTRCL+3, MNPLTR )
         DO 10 I=1,3*NBTRCL+3
            NP = LADEFI(WPTRCL-1+I)
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
               KERR(1) ='POINT SANS COORDONNEES' // KERR(MXLGER)(1:10)
               CALL LEREUR
               IERR = 6
               GOTO 10
            ENDIF
            MCN( MNPLTR - 1 + I ) = MNS + WYZSOM
 10      CONTINUE
C
      ELSE IF( NUTYCL .EQ. 1 ) THEN
C
C        CHAQUE LIGNE EST DEFINIE PAR LES POINTS DE 3 LIGNES
C        ===================================================
C        RECHERCHE DES 3 LIGNES
         NBTRCL = 0
         DO 20 I=0,2
C           LE NUMERO DE LA LIGNE
            NL = LADEFI(WLTRCL+I)
            IF( NL .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') I+1
               KERR(1) =  'LIGNE GUIDE INCORRECTE '//KERR(MXLGER)(1:4)
               CALL LEREUR
               IERR = 9
               GOTO 20
            ENDIF
C           OUVERTURE DE LA LIGNE
            CALL LXNLOU( NTLIGN, NL, NT, MN )
            IF( NT .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') I+1
               KERR(1) = 'LIGNE GUIDE INCONNUE ' // KERR(MXLGER)(1:10)
               CALL LEREUR
               IERR = 10
               GOTO 20
            ENDIF
            NTLIGS(I+1) = NT
C           OUVERTURE DU TABLEAU SOMMETS
            CALL LXTSOU( NT, 'XYZSOMMET', NTS, MNS )
            IF( NTS .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') I+1
               KERR(1) =  'LIGNE GUIDE '//KERR(MXLGER)(1:4)//
     %                    ' SANS SOMMETS'
               CALL LEREUR
               IERR = 6
               GOTO 20
            ENDIF
C           OUVERTURE DU TABLEAU NSEF
            CALL LXTSOU( NT, 'NSEF', NTA, MNA )
            IF( NTA .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') I+1
               KERR(1) =  'LIGNE GUIDE ' // KERR(MXLGER)(1:4)
     %                    //' SANS NSEF'
               CALL LEREUR
               IERR = 6
               GOTO 20
            ENDIF
C
C           LES PARAMETRES DES NO SOMMET DE LA LIGNE
            CALL NSEFPA( MCN(MNA),
     %                   NUTYL, NBSOLL, NBSOSL, NBTG,
     %                   LDAP,  LDNG, LDTG, NBARSL,
     %                   NXL  , NYL  , NZL  ,
     %                   IERR  )
            IF( IERR .NE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') I+1
               KERR(1) = 'MAILLAGE INCORRECT DE LA LIGNE GUIDE '
     %                 // KERR(MXLGER)(1:4)
               CALL LEREUR
               IERR = 7
               RETURN
            ENDIF
            IF( NBTRCL .LE. 0 ) THEN
               NBTRCL = NBARSL
            ELSE
               IF( NBTRCL .NE. NBARSL ) THEN
                  NBLGRC(NRERR) = 3
                  WRITE(KERR(MXLGER)(1:4),'(I4)') NBARSL
                  WRITE(KERR(MXLGER)(11:14),'(I4)') NBTRCL
                  KERR(1) = 'LES 3 LIGNES DOIVENT ASUIR'
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
         IF( NBTRCL .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(MXLGER)(1:4),'(I4)') NBTRCL
            KERR(1) = 'NOMBRE INCORRECT DE TRANCHES='
     %                // KERR(MXLGER)(1:4)
            CALL LEREUR
            IERR = 4
         ENDIF
         IF( IERR .NE. 0 ) GOTO 9990
C
C        LE TABLEAU DES ADRESSES DES 3 COORDONNEES DES 3 * NBTRCL POINTS
C        DEFINISSANT LE PLAN DES NBTRCL LIGNES FINALES DES TRANCHES
         CALL TNMCDC( 'ENTIER', 3*NBTRCL+3, MNPLTR )
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
            IF( IERR .NE. 0 .AND. IERR .NE. 2 ) GOTO 9990
C
C           LES PARAMETRES DES NO SOMMET DE LA LIGNE
            CALL NSEFPA( MCN(MNA),
     %                   NUTYL, NBSOLL, NBSOSL, NBTG,
     %                   LDAP,  LDNG,   LDTG,   NBARSL,
     %                   NXL,   NYL,    NZL,
     %                   IER   )
C
C           BOUCLE SUR LES ARETES DE LA LIGNE
            MNSOMS(J) = MNS + WYZSOM
C
C           L'ADRESSE MCN DES SOMMETS DE LA LIGNE
            MN   = MNPLTR - 1 + J
            MNSO = MNSOMS(J)
            DO 30 I=1,NBARSL
               MCN( MN ) = MNSO
               MN   = MN   + 3
               MNSO = MNSO + 3
 30         CONTINUE
            IF( IERR .EQ. 2 ) THEN
C              LE DERNIER SOMMET EST LE PREMIER
               MCN( MN ) = MNSOMS(J)
               IERR = 0
            ELSE
C              LE DERNIER SOMMET
               MCN( MN ) = MNSO
            ENDIF
 50      CONTINUE
         IF( IERR .NE. 0 ) GOTO 9990
C
      ELSE IF( NUTYCL .GE. 2 .AND. NUTYCL .LE. 5 ) THEN
C
C        COTE MAXIMALE OU COTE DE CHAQUE TRANCHE
C        VECTEUR TRANSLATION CONSTANT OU VARIABLE PAR TRANCHE
C        LE NOMBRE DE TRANCHES
         NBTRCL = LADEFI(WBTRCL)
C        DECLARATION D'UN VECTEUR 3*NBTRCL + NBTRCL+1
         I = 4*NBTRCL+1
         CALL TNMCDC( 'REEL', I, MNVECT )
         MNZ = MNVECT + 3 * NBTRCL
         CALL AZEROR( I, RMCN(MNVECT) )
C
C        REMPLISSAGE DU TABLEAU VECTOR SELON NUTYCL
         IF( NUTYCL .EQ. 2 ) THEN
C           LA COTE MAXIMALE
            RMCN( MNVECT - 1 + 3 * NBTRCL ) = RADEFI( WMAXCL )
         ELSE IF( NUTYCL .EQ. 3 ) THEN
C           LA COTE DE CHAQUE LIGNE DES TRANCHES
            DO 52 J=1,NBTRCL
               RMCN(MNVECT-1+3*J) = RADEFI(WCTRCL-1+J)
 52         CONTINUE
         ELSE IF( NUTYCL .EQ. 4 ) THEN
C           LE VECTEUR CONSTANT DE TRANSLATION
            RMCN(MNVECT  ) = RADEFI(WRCOCL)
            RMCN(MNVECT+1) = RADEFI(WRCOCL+1)
            RMCN(MNVECT+2) = RADEFI(WRCOCL+2)
         ELSE IF( NUTYCL .EQ. 5 ) THEN
C           LE VECTEUR DE TRANSLATION DE CHAQUE LIGNE
            CALL TRTATA( RADEFI(WRVACL), RMCN(MNVECT), 3*NBTRCL )
         ENDIF
C
      ELSE IF( NUTYCL .EQ. 6 ) THEN
C
C        LE NOMBRE DE TRANCHES
         NBTRCL = LADEFI(WBTRCL)
C
      ELSE
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYCL
         KERR(1) = 'OPTION INCONNUE' // KERR(MXLGER)(1:4)
         CALL LEREUR
         IERR = 13
      ENDIF
      IF( IERR .NE. 0 ) GOTO 9990
C
C     LE SURFACE EST ELLE UN TORE C-A-D
C     ---------------------------------
C     LA DERNIERE LIGNE SE TROUVE ETRE AUSSI LA PREMIERE?
      MN = MNPLTR + 3 * NBTRCL
      IF( MCN(MNPLTR) .EQ. MCN(MN) .AND. MCN(MNPLTR+1) .EQ. MCN(MN+1)
     %    .AND. MCN(MNPLTR+2) .EQ. MCN(MN+2) ) THEN
C        IL S'AGIT D'UN TORE
         LETORE = 1
      ELSE
C        CE N'EST PAS UN TORE
         LETORE = 0
      ENDIF
C
C     GENERATION DES NSEF DE LA SURFACE
C     ========================================
C     CONSTRUCTION DU TABLEAU 'NSEF' DE CE SURFACE
      NBFASU = NBTRCL * NBARLI
      NBTL   = 0
      IF( NBTGSL .GT. 0 .AND. NUTYCL .NE. 6 ) THEN
         NBTGEF = 8
         NBEFAP = NBFASU
         NBEFTG = NBFASU
         IF( NUTYCL .LE. 1 ) THEN
C           LES TANGENTES CHANGENT POUR CHAQUE TRANCHE
            NBTGSU = NBTGSL * (1+NBTRCL-LETORE)
            NBTL   = NBTGSL
         ELSE IF( NUTYCL .LE. 5 ) THEN
C           TRANSLATION SEULEMENT OU TANGENTE INCALCULABLE
            NBTGSU = NBTGSL
         ENDIF
      ELSE
C        TANGENTE INCALCULABLE
         NBTGSU = 0
         NBEFAP = 0
         NBEFTG = 0
         NBTGEF = 0
      ENDIF
C
C     LE NOMBRE DE MOTS DU TMS 'NSEF'
      N = WUSOEF + 4*NBFASU + NBEFAP + NBEFTG * ( 1+NBTGEF )
      CALL LXTNDC(NTLXSU,'NSEF','ENTIER',N)
      CALL LXTSOU(NTLXSU,'NSEF',NTFASU,MNFASU)
      MNEFAP = MNFASU + WUSOEF + 4*NBFASU
      MNNTGS = MNEFAP + NBEFAP + NBEFTG
C
C     MISE A JOUR DU TABLEAU 'NSEF' DE CE SURFACE
C     TYPE DE L'OBJET : SURFACE
      MCN( MNFASU + WUTYOB ) = 3
C     LE TYPE INCONNU DE FERMETURE DU MAILLAGE
      MCN( MNFASU + WUTFMA ) = -1
C     NOMBRE DE SOMMETS PAR EF
      MCN( MNFASU + WBSOEF ) = 4
C     LE NOMBRE DE TANGENTES STOCKEES PAR EF : SURFACE C1 OU C0
      MCN( MNFASU + WBTGEF ) = NBTGEF
C     NBFASU  LE NOMBRE DE FACES DE LA SURFACE
      MCN( MNFASU + WBEFOB ) = NBFASU
C     LE NOMBRE D'EF A TG
      MCN( MNFASU + WBEFTG ) = NBEFTG
C     NOMBRE DES EF AVEC POINTEUR SUR EF A TG
      MCN( MNFASU + WBEFAP ) = NBEFAP
C     NUMERO DU TYPE DE MAILLAGE : NON STRUCTURE
      MCN( MNFASU + WUTYMA ) = 0
C
C     LA BOUCLE SUR LES TRANCHES DE LA SURFACE
C     ----------------------------------------
      LEDEC0 = 0
      LEDEC1 = NBSOML
      NBTGS0 = 0
      NBTGS1 = NBTL
      MNC    = MNFASU + WUSOEF - 1
      MNT    = MNNTGS - 1
      DO 300 NBT = 1, NBTRCL
C
C        SI TORE ET DERNIERE TRANCHE => DERNIERE LIGNE=PREMIERE LIGNE
         IF( LETORE .EQ. 1 .AND. NBT .EQ. NBTRCL ) THEN
            LEDEC1 = 0
            NBTGS1 = 0
         ENDIF
C
C        LA BOUCLE SUR LES ARETES DU MAILLAGE DE LA LIGNE
         DO 100 N=1,NBARLI
C           LE NUMERO DES NBSOEF SOMMETS DE L'ARETE N
            CALL NSEFNS( N     , NUTYML, NBSOEF, NBTGEL,
     %                   LDAPEF, LDNGEF, LDTGEF,
     %                   MNARLI, NX, NY, NZ,
     %                   NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
            IF( NCOGEL .NE. 2 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') NCOGEL
               KERR(1) = 'LIGNE AVEC UN ELEMENT FINI DE CODE'
     %                 // KERR(MXLGER)(1:4)
               CALL LEREUR
               IERR = 8
               CALL LXTSDS( NTLXSU, 'NSEF' )
               GOTO 9990
            ENDIF
C           LE NUMERO DES 4 SOMMETS DE LA FACE
            MCN( MNC + 1 ) = NOSOEL(1) + LEDEC0
            MCN( MNC + 2 ) = NOSOEL(2) + LEDEC0
            MCN( MNC + 3 ) = NOSOEL(2) + LEDEC1
            MCN( MNC + 4 ) = NOSOEL(1) + LEDEC1
C
C           LES EVENTUELLES 8 TANGENTES DE LA FACE
            IF( NBTGSL .GT. 0 .AND. NUTYCL .NE. 6 ) THEN
C              LA TANGENTE A LA COURBE AU SOMMET 1
               IF( NOSOEL(3) .GT. 0 ) THEN
                  MCN( MNT + 1 ) = NOSOEL(3) + NBTGS0
                  MCN( MNT + 8 ) = NOSOEL(3) + NBTGS1
               ELSE IF( NOSOEL(3) .LT. 0 ) THEN
                  MCN( MNT + 1 ) = NOSOEL(3) - NBTGS0
                  MCN( MNT + 8 ) = NOSOEL(3) - NBTGS1
               ELSE
                  MCN( MNT + 1 ) = 0
                  MCN( MNT + 8 ) = 0
               ENDIF
               MCN( MNT + 2 ) = 0
               MCN( MNT + 3 ) = 0
C              LA TANGENTE A LA COURBE AU SOMMET 2
               IF( NOSOEL(4) .GT. 0 ) THEN
                  MCN( MNT + 4 ) = NOSOEL(4) + NBTGS0
                  MCN( MNT + 5 ) = NOSOEL(4) + NBTGS1
               ELSE IF( NOSOEL(4) .LT. 0 ) THEN
                  MCN( MNT + 4 ) = NOSOEL(4) - NBTGS0
                  MCN( MNT + 5 ) = NOSOEL(4) - NBTGS1
               ELSE
                  MCN( MNT + 4 ) = 0
                  MCN( MNT + 5 ) = 0
               ENDIF
               MCN( MNT + 6 ) = 0
               MCN( MNT + 7 ) = 0
               MNT = MNT + 8
            ENDIF
C
C           LA FACE SUIVANTE
            MNC = MNC + 4
 100     CONTINUE
C
C        LA LIGNE SUIVANTE
         LEDEC0 = LEDEC1
         LEDEC1 = LEDEC1 + NBSOML
         NBTGS0 = NBTGS1
         NBTGS1 = NBTGS1 + NBTL
 300  CONTINUE
C
      IF( NBTGSU .GT. 0 ) THEN
C
C        LE POINTEUR SUR LES EF A TG
         DO 310 N=1,NBEFAP
            MCN(MNEFAP-1+N) = N
 310     CONTINUE
C
C        LE CODE GEOMETRIQUE 0
         NBT = MNEFAP-1+NBEFAP
         DO 320 N=1,NBEFTG
            MCN(NBT+N) = 0
 320     CONTINUE
C
      ENDIF
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNFASU) )
C
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFASU + MOTVAR(6) ) = NONMTD ( '~>>>NSEF' )
C
C     GENERATION DES SOMMETS ET DES TANGENTES DE LA SURFACE
C     ============================================ =========
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CE SURFACE
C     LE NOMBRE DE SOMMETS SELON TORE OU NON
      IF( LETORE .EQ. 0 ) THEN
C        NON TORE
         NBSOM = (NBTRCL+1) * NBSOML
      ELSE
C        TORE
         NBSOM =  NBTRCL * NBSOML
      ENDIF
      N = WYZSOM + 3 * ( NBSOM + NBTGSU )
      CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'MOTS', N )
      CALL LXTSOU( NTLXSU, 'XYZSOMMET', NTSOFA, MNSOFA )
C
C     CALCUL DES 3 COORDONNEES DES NBSOM SOMMETS DE LA SURFACE
      CALL VOE114( NUTYCL, LETORE, NBTRCL, MCN(MNPLTR),
     %             LADEFI(WONXCL),
     %             MCN(MNVECT), MCN( MNZ ),
     %             NBSOML, RMCN(MNSOLI+WYZSOM), RMCN(MNSOFA+WYZSOM),
     %             NBTGSL, RMCN(MNSOLI+WYZSOM+3*NBSOML),
     %             RMCN(MNSOFA+WYZSOM+3*NBSOM),
     %             IERR   )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1)='PB DANS LA GENERATION DES SOMMETS DE LA SURFACE'
         CALL LEREUR
         CALL LXTSDS( NTLXSU, 'XYZSOMMET' )
         GOTO 9990
      ENDIF
C
C     CALCUL DES 3 COORDONNEES DES TANGENTES DE LA SURFACE
      IF( NUTYCL .GE. 1 .AND. NUTYCL .LE. 5 ) THEN
C        TYPE 2 3 4 5 => TRANSLATION SIMPLE LES TGS SONT CELLES DE LA LIGNE
         CALL TRTATA( RMCN(MNSOLI+WYZSOM+3*NBSOML),
     %                RMCN(MNSOFA+WYZSOM+3*NBSOM),
     %                3 * NBTGSL )
      ENDIF
C
C     NBSOM 'NOMBRE DE SOMMETS'
      MCN( MNSOFA + WNBSOM ) = NBSOM
C     LE NOMBRE DE TANGENTES STOCKEES DE CETTE SURFACE
      MCN( MNSOFA + WNBTGS ) = NBTGSU
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSOFA + WBCOOR ) = 3
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOFA) )
C
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOFA + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
      IF( LETORE .EQ. 0 .AND. NUTFML .EQ. 0 .AND.
     %    NBSOML .EQ. NBARLI+1 ) THEN
C
C        SI LE MAILLAGE EST UN QUADRANGLE ALORS IL EST RESTRUCTURE
C        =========================================================
         NBEFOB = MCN( MNFASU + WBEFOB )
         NBSOEF = MCN( MNFASU + WBSOEF )
         NBEFAP = MCN( MNFASU + WBEFAP )
         NBEFTG = MCN( MNFASU + WBEFTG )
         NBTGEF = MCN( MNFASU + WBTGEF )
C
C        TRANSLATION DE 4*NBEFOB AU DELA DU NUMERO DES SOMMETS
         I = WUSOEF + NBSOEF * NBEFOB - WBARYQ - 1
         N = MNFASU + WUSOEF + NBSOEF * NBEFOB
         DO 200 MN = N, N+NBEFAP+NBEFTG*(1+NBTGEF)-1
            MCN(MN-I) = MCN(MN)
200      CONTINUE
C
C        NUMERO DU TYPE DU MAILLAGE : QUADRANGLE STRUCTURE
         MCN ( MNFASU + WUTYMA ) = 4
C        VARIABLE NBARXQ : NOMBRE DE SEGMENTS SUIVANT X
         MCN ( MNFASU + WBARXQ ) = NBARLI
C        VARIABLE NBARYQ : NOMBRE DE SEGMENTS SUIVANT Y
         MCN ( MNFASU + WBARYQ ) = NBTRCL
C
C        REDUCTION DU TMS 'NSEF'
         CALL TAMSRA( NTFASU, WBARYQ+1+NBEFAP+NBEFTG*(1+NBTGEF) )
      ENDIF
C
C     DESTRUCTION DES TABLEAUX DEVENUS INUTILES
 9990 IF( MNPLTR .GT. 0 ) CALL TNMCDS( 'ENTIER', 3*NBTRCL+3, MNPLTR )
      IF( MNVECT .GT. 0 ) CALL TNMCDS( 'REEL'  , 4*NBTRCL+1, MNVECT )
      END
