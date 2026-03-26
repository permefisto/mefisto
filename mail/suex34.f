        SUBROUTINE SUEX34( NTLXSU, LADEFI,
     %                     NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE D'UNE QUADRANGULATION STRUCTUREE
C -----    EXTRAITE D'UN HEXAEDRE STRUCTURE AVEC OU SANS TANGENTES
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE A CREER
C LADEFI : TABLEAU DE DEFINITION DE LA SURFACE PARTITIONNEE
C          CF '~/TD/D/A_SURFACE__DEFINITION'
C
C SORTIES:
C --------
C NTFASU : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES
C          CF '~/TD/D/A___NSEF'
C NTSOFA : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOFA : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF '~/TD/D/A___XYZSOMMET'
C IERR   : =0 SI PAS D'ERREUR
C          >0 EN CAS D'ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : CHRISTOPHE DOURSAT ANALYSE NUMERIQUE PARIS     SEPTEMBRE 1989
C AUTEUR : ALAIN PERRONNET    ANALYSE NUMERIQUE PARIS     NOVEMBRE  1996
C.......................................................................
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL           RMCN(1)
      EQUIVALENCE   (MCN(1),RMCN(1))
      INTEGER        LADEFI(0:*)
      INTEGER        NBTGFA(6), NOTGFA(8,6)
C
      IERR = 0
C
C     RESTAURATION DU MAILLAGE DU VOLUME INITIAL
C     ==========================================
      CALL LXNLOU( NTVOLU, LADEFI(WUVOST), NTVOST, MN )
      IF( NTVOST .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SUEX34: VOLUME INCONNU'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     RESTAURATION DES TABLEAUX SOMMETS ET NSEF
      CALL LXTSOU( NTVOST, 'XYZSOMMET', NTSOVO, MNSOVO )
      IF( NTSOVO .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SUEX34: VOLUME SANS SOMMETS'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
      CALL LXTSOU( NTVOST, 'NSEF', NTSSVO, MNSSVO )
      IF( NTSSVO .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SUEX34: VOLUME SANS NSEF'
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C     LE NOMBRE DE SOMMETS DU VOLUME
      NBSOMV = MCN( MNSOVO + WNBSOM )
C     LE NOMBRE DE TANGENTES
      NBTGSV = MCN( MNSOVO + WNBTGS )
C
C     L'ADRESSE DU DEBUT DES COORDONNEES DES SOMMETS
      MNCOVO = MNSOVO + WYZSOM
C     TYPE DE L'OBJET : VOLUME
      NTYPVO = MCN(MNSSVO + WUTYOB )
      IF( NTYPVO .NE. 4 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SUEX34: TYPE DIFFERENT DE VOLUME'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     NUMERO DU TYPE DU MAILLAGE : HEXAEDRE STRUCTURE
      NTYPMA = MCN ( MNSSVO + WUTYMA )
      IF (NTYPMA.NE.7) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SUEX34: VOLUME NON STRUCTURE'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     VARIABLE NBARXH : NOMBRE DE SEGMENTS SUIVANT X
      NBARXH = MCN ( MNSSVO + WBARXH )
      NBS1   = NBARXH + 1
C     VARIABLE NBARYH : NOMBRE DE SEGMENTS SUIVANT Y
      NBARYH = MCN ( MNSSVO + WBARYH )
      NBS2   = NBARYH + 1
C     VARIABLE NBARZH : NOMBRE DE SEGMENTS SUIVANT Z
      NBARZH = MCN ( MNSSVO + WBARZH )
      NBS3   = NBARZH + 1
C
C     LES TANGENTES
      NBEFOB = MCN( MNSSVO + WBEFOB )
      NBEFAP = MCN( MNSSVO + WBEFAP )
      NBEFTG = MCN( MNSSVO + WBEFTG )
      NBTGEF = MCN( MNSSVO + WBTGEF )
      MNEFAP = MNSSVO + WBARZH + 8*NBEFOB
      MNCGEF = MNEFAP + NBEFAP
      MNTGEF = MNCGEF + NBEFTG
C
C     CALCUL ET VERIFICATION DES BORNES DE LA SURFACE EXTRAITE
C     ========================================================
      IMIN = LADEFI(WI2MIN)
      IMAX = LADEFI(WI2MAX)
      JMIN = LADEFI(WJ2MIN)
      JMAX = LADEFI(WJ2MAX)
      KMIN = LADEFI(WK2MIN)
      KMAX = LADEFI(WK2MAX)
C
      IF ((IMIN.LT.1).OR.(JMIN.LT.1).OR.(KMIN.LT.1)) GOTO 9990
      IF ((IMAX.GT.NBS1).OR.(JMAX.GT.NBS2).OR.(KMAX.GT.NBS3)) GOTO 9990
      IDIF = IMAX-IMIN
      JDIF = JMAX-JMIN
      KDIF = KMAX-KMIN
      IF ((IDIF.LT.0).OR.(JDIF.LT.0).OR.(KDIF.LT.0)) GOTO 9990
      PROD = IDIF*JDIF*KDIF
      IF (PROD.NE.0) THEN
         NBLGRC(NRERR) = 1
         KERR(1)='SUEX34:MIN doit etre <= MAX'
         CALL LEREUR
         IERR = 1
         GOTO 9990
      ELSE
        PROD = IDIF*JDIF+IDIF*KDIF+JDIF*KDIF
        IF (PROD.EQ.0) THEN
           NBLGRC(NRERR) = 1
           KERR(1)=
     %'SUEX34:MIN MAX DOIVENT ETRE EGAUX POUR AU PLUS 1 des I ou J ou K'
           CALL LEREUR
           IERR = 1
           GOTO 9990
        ENDIF
      ENDIF
C
C     NOMBRE DE SOMMETS DES FACES DE LA SURFACE
      NBSOFA = (IDIF+1)*(JDIF+1)*(KDIF+1)
C
C     NOMBRE MAXIMAL DE TANGENTES DES FACES DE LA SURFACE EXTRAITE
      NBFASU = 1
      IF( IDIF .GT. 0 ) NBFASU = IDIF
      IF( JDIF .GT. 0 ) NBFASU = NBFASU * JDIF
      IF( KDIF .GT. 0 ) NBFASU = NBFASU * KDIF
      IF( NBTGSV .LE. 0 ) THEN
C        PAS DE TANGENTES POUR LE VOLUME ET LA SURFACE
         NBTGSF = 0
         NBTGEF = 0
         NBEFTG = 0
      ELSE
C        8 TANGENTES MAXIMUM PAR FACE
         NBTGSF = 8 * NBFASU
         NBTGEF = 8
      ENDIF
C
C     GENERATION DES SOMMETS DE LA SURFACE CREEE: TMS XYZSOMMET
C     =========================================================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CETTE SURFACE
      CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'MOTS',
     %             WYZSOM+3*(NBSOFA+NBTGSF) )
      CALL LXTSOU( NTLXSU, 'XYZSOMMET',  NTSOFA, MNSOFA )
C
C     L ADRESSE DU DEBUT DES COORDONNES DES SOMMETS
      MNCOSU = MNSOFA + WYZSOM
C
C     REMPLISSAGE DES COORDONNEES DES SOMMETS DE LA SURFACE
      DO 40 K=KMIN,KMAX
        DO 30 J=JMIN,JMAX
          DO 20 I=IMIN,IMAX
            NEULVO = (K-1)*(NBS1*NBS2)+(J-1)*NBS1+I
            NEULSU = (K-KMIN)*(IDIF+1)*(JDIF+1)+(J-JMIN)*(IDIF+1)
     S               +I-IMIN+1
            IADVO = MNCOVO-1+3*(NEULVO-1)
            IADSU = MNCOSU-1+3*(NEULSU-1)
            DO 10 N=1,3
              RMCN(IADSU+N) = RMCN(IADVO+N)
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
   40 CONTINUE
C
C     LES TANGENTES
      NBTGS1 = 0
      IF( NBTGSV .GT. 0 ) THEN
C
C        LE NUMERO DES 8 TANGENTES DES 6 FACES DE L'HEXAEDRE
         CALL TGFACU( 7, NBTGFA, NOTGFA )
C
C        RESERVATION DU TABLEAU DES 8 NUMEROS  DES TANGENTES DES FACES
         CALL TNMCDC( 'ENTIER', NBFASU * 9, MNEFA )
C        MNEFA ADRESSE AVANT LES POINTEURS SUR LES EF A TG
C        MNEFT ADRESSE AVANT LES 8 NUMEROS DE TANGENTES DES FACES A TG
         MNEFT  = MNEFA + NBFASU
         NBEF   = 0
         NBEFTG = 0
C
C        LES 8 NUMEROS DES TANGENTES A EXTRAIRE
         IF( KDIF .EQ. 0 ) THEN
C           KMIN=KMAX => FACES 1 OU 4 DES HEXAEDRES A EXTRAIRE
            IF( KMIN .EQ. 1 ) THEN
C              LE NUMERO DANS L'HEXAEDRE DES 8 TANGENTES DE LA FACE 1
               NUFACE = 1
            ELSE
C              LE NUMERO DANS L'HEXAEDRE DES 8 TANGENTES DE LA FACE 4
               NUFACE = 4
            ENDIF
            I1 = IMIN
            I2 = IMAX
            J1 = JMIN
            J2 = JMAX
            K  = 3
C
         ELSE IF( JDIF .EQ. 0 ) THEN
C           JMIN=JMAX => FACES 2 OU 5 DES HEXAEDRES A EXTRAIRE
            IF( JMIN .EQ. 1 ) THEN
C              LE NUMERO DANS L'HEXAEDRE DES 8 TANGENTES DE LA FACE 2
               NUFACE = 2
            ELSE
C              LE NUMERO DANS L'HEXAEDRE DES 8 TANGENTES DE LA FACE 5
               NUFACE = 5
            ENDIF
            I1 = IMIN
            I2 = IMAX
            J1 = KMIN
            J2 = KMAX
            K  = 2
C
         ELSE IF( IDIF .EQ. 0 ) THEN
C           IMIN=IMAX => FACES 3 OU 6 DES HEXAEDRES A EXTRAIRE
            IF( IMIN .EQ. 1 ) THEN
C              LE NUMERO DANS L'HEXAEDRE DES 8 TANGENTES DE LA FACE 3
               NUFACE = 3
            ELSE
C              LE NUMERO DANS L'HEXAEDRE DES 8 TANGENTES DE LA FACE 6
               NUFACE = 6
            ENDIF
            I1 = JMIN
            I2 = JMAX
            J1 = KMIN
            J2 = KMAX
            K  = 1
C
C           PARCOURS DES FACES EXTRAITES
            DO 80 J=J1,J2
               DO 70 I=I1,I2
C                 LE NUMERO DE L'HEXAEDRE
                  GOTO( 51, 52, 53 ),K
C                 IMIN=IMAX
 51               NH = IMIN + NBARXH * ( (I-1) + NBARYH * (J-1) )
                  GOTO 60
C                 JMIN=JMAX
 52               NH = I + NBARXH * ( (JMIN-1) + NBARYH * (J-1) )
                  GOTO 60
C                 KMIN=KMAX
 53               NH = I + NBARXH * ( (J-1) + NBARYH * (KMIN-1) )
C
C                 LE NUMERO DES TANGENTES
 60               NUEFTG = MCN(MNEFAP+NH)
                  NBEF   = NBEF  + 1
                  IF( NUEFTG .LE. 0 ) THEN
C                    HEXAEDRE SANS TG
                     MCN( MNEFA + NBEF ) = 0
                     CALL AZEROI( 8, MCN(MNEFT+1) )
                  ELSE
C                    HEXAEDRE AVEC TG
                     MN    = MNTGEF + 24 * NUEFTG - 24
                     MNEFT = MNEFT  + 8
                     N     = 0
                     DO 65 L=1,8
                        NUTG = MCN( MN + NOTGFA(L,NUFACE) )
                        MCN(MNEFT+L) = NUTG
                        IF( NUTG .NE. 0 ) N = N + 1
 65                  CONTINUE
                     IF( N .EQ. 0 ) THEN
C                       FACE SANS TG
                        MNEFT = MNEFT - 8
                        MCN( MNEFA + NBEF ) = 0
                     ELSE
C                       FACE AVEC TG
                        NBEFTG = NBEFTG + 1
                        MCN( MNEFA + NBEF ) = NBEFTG
                     ENDIF
                  ENDIF
 70            CONTINUE
 80         CONTINUE
         ENDIF
C
C        RECENSEMENT DES TANGENTES DES EF EXTRAITS
C        =========================================
C        LE NOUVEAU NUMERO DES TANGENTES RECENSEES
         CALL TNMCDC( 'ENTIER', NBTGSF+1, MNNEWT )
         MN     = MNEFA + NBFASU
         DO 100 I=MN, MN+8*NBEFTG-1
C           LE NUMERO DE LA TANGENTE EST RECENSE
            N = ABS( MCN(I) )
            MCN(MNNEWT+N) = N
 100     CONTINUE
         DO 110 I=1,NBTGSF
C           LE NUMERO ANCIEN DE LA TANGENTE
            N = MCN(MNNEWT+I)
            IF( N .NE. 0 ) THEN
C              LE NUMERO NOUVEAU DE LA TANGENTE
               NBTGS1 = NBTGS1 + 1
               MCN(MNNEWT+I) = NBTGS1
            ENDIF
 110     CONTINUE
         IF( NBTGS1 .EQ. 0 ) THEN
C           SURFACE SANS TG
            NBEFAP = 0
            NBEFTG = 0
            NBTGEF = 0
            GOTO 150
         ENDIF
C
C        LES 3 COMPOSANTES DES TANGENTES DU MAILLAGE
         MN0 = MNSOVO + WYZSOM - 3 + 3 * NBSOMV
         MN1 = MNSOFA + WYZSOM + 3 * NBSOFA
         DO 120 I=1,NBTGSF
            N = MCN(MNNEWT+I)
            IF( N .NE. 0 ) THEN
               MN  = MN0 + 3 * I
               RMCN(MN1  ) = RMCN(MN  )
               RMCN(MN1+1) = RMCN(MN+1)
               RMCN(MN1+2) = RMCN(MN+2)
               MN1 = MN1 + 3
            ENDIF
 120     CONTINUE
C
C        +- LE NUMERO DES NBTGEF TG
         MN0 = MNEFA + NBFASU
         MN3 = MN0
         DO 140 L=0,NBTGEF-1
C           +- L'ANCIEN NUMERO DE LA TG
            N = MCN(MN0+L)
            IF( N .LT. 0 ) THEN
               LESIGN = -1
               N      = -N
            ELSE
               LESIGN = 1
            ENDIF
C           +- LE NOUVEAU NUMERO DE LA TG
            MCN(MN3) = LESIGN * MCN(MNNEWT+N)
            MN3 = MN3 + 1
 140     CONTINUE
C
C        DESTRUCTION DU TABLEAU DEVENU INUTILE
 150     CALL TNMCDS( 'ENTIER', NBTGSF+1, MNNEWT )
      ENDIF
C
C     LE NOMBRE DE SOMMETS
      MCN( MNSOFA + WNBSOM ) = NBSOFA
C     LE NOMBRE DE TANGENTES DE LA SURFACE
      MCN( MNSOFA + WNBTGS ) = NBTGS1
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSOFA + WBCOOR ) = 3
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOFA) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOFA + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
      CALL TAMSRA( NTSOFA, WYZSOM + 3*(NBSOFA+NBTGS1) )
C
C     GENERATION DES NUMEROS DES SOMMETS DES EF: TMS NSEF
C     ===================================================
C     CONSTRUCTION DU TABLEAU 'NSEF' DE CETTE SURFACE
      CALL LXTNDC( NTLXSU, 'NSEF', 'ENTIER',
     %             1+WBARYQ+NBEFAP+NBEFTG*(1+NBTGEF) )
      CALL LXTSOU( NTLXSU, 'NSEF',  NTFASU,  MNFASU )
C     TYPE DE L'OBJET : SURFACE
      MCN( MNFASU + WUTYOB ) = 3
C     NUMERO DU TYPE DU MAILLAGE: QUADRANGLE STRUCTURE
      MCN( MNFASU + WUTYMA ) = 4
      IF( IDIF .EQ. 0 ) THEN
        IDIF = JDIF
        JDIF = KDIF
      ELSE
        IF( JDIF .EQ. 0 ) THEN
          JDIF = KDIF
        ENDIF
      ENDIF
C     LE TYPE INCONNU DE FERMETURE DU MAILLAGE
      MCN( MNFASU + WUTFMA ) = -1
C     LE NOMBRE DE SOMMETS PAR FACE
      MCN( MNFASU + WBSOEF ) = 4
C     LE NOMBRE DE TANGENTES STOCKEES PAR EF
      MCN( MNFASU + WBTGEF ) = NBTGEF
      MCN( MNFASU + WBEFAP ) = NBEFAP
      MCN( MNFASU + WBEFTG ) = NBEFTG
C     LE NOMBRE D'EF DE LA SURFACE
      MCN( MNFASU + WBEFOB ) = IDIF * JDIF
C     VARIABLE NBARXQ : NOMBRE DE SEGMENTS SUIVANT X
      MCN( MNFASU + WBARXQ ) = IDIF
C     VARIABLE NBARYQ : NOMBRE DE SEGMENTS SUIVANT Y
      MCN( MNFASU + WBARYQ ) = JDIF
C
      IF( NBEFTG .GT. 0 ) THEN
C        LES NUMEROS DES EF A TG
         MN = MNFASU + WBARYQ + 1
         CALL TRTATA( MCN(MNEFA), MCN(MN), NBFASU )
C        LE CODE GEOMETRIQUE
         DO 170 I=1,NBEFTG
            MCN(MN+I) = 0
 170     CONTINUE
C        LES NUMEROS DES TANGENTES DES FACES A TG
         CALL TRTATA( MCN(MNEFA+NBFASU), MCN(MN+NBEFTG), 8*NBEFTG )
C        DESTRUCTION DU TABLEAU DEVENU INUTILE
         CALL TNMCDS( 'ENTIER', NBFASU * 9, MNEFA )
      ENDIF
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNFASU) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFASU + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
      RETURN
C
C     ERREUR DANS LES DONNEES
 9990 NBLGRC(NRERR) = 4
      KERR(1) = 'SUEX34: DEFINITION INCORRECTE DES BORNES DES SOMMETS'
      WRITE(KERR(MXLGER)(1:10),'(I10)') NBS1
      KERR(2) = 'I COMPRIS ENTRE 1 ET ' // KERR(MXLGER)(1:10)
      WRITE(KERR(MXLGER)(1:10),'(I10)') NBS2
      KERR(3) = 'J COMPRIS ENTRE 1 ET ' // KERR(MXLGER)(1:10)
      WRITE(KERR(MXLGER)(1:10),'(I10)') NBS3
      KERR(4) = 'K COMPRIS ENTRE 1 ET ' // KERR(MXLGER)(1:10)
      CALL LEREUR
      IERR = 1
      RETURN
      END
