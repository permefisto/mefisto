      SUBROUTINE SUEX13( NTLXSU, LADEFI,
     &                   NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     GENERER LA QUADRANGULATION D'UNE LIGNE FERMEE AYANT UN
C ----      NOMBRE PAIR D'ARETES ET DE SOMMET LE BARYCENTRE
C          (LA LIGNE FERMEE DOIT ETRE CONVEXE, CONNEXE SANS INTERSECTION)
C
C ENTREES :
C --------
C NTLXSU  : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE
C LADEFI  : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
C           CF '~/TD/D/A_SURFACE__DEFINITION'
C
C SORTIES :
C ---------
C NTFASU  : NUMERO      DU TMS 'NSEF' DES NUMEROS DES ELEMENTS FINIS
C MNFASU  : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES ELEMENTS FINIS
C           CF '~/TD/D/A___NSEF'
C NTSOFA  : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOFA  : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C           CF '~/TD/D/A___XYZSOMMET'
C IERR    : 0 SI PAS D'ERREUR
C         > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR  : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    JANVIER 2006
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
      REAL              XYZB(3)
      INTEGER           NOSOEL(8)
C
      IERR   = 0
C
C     LE NUMERO DE LA LIGNE CONVEXE FERMEE
      NULICE = LADEFI( WULICE )
C     LE LEXIQUE DE LA LIGNE
      CALL LXNLOU( NTLIGN, NULICE, NTLXLI, MNLXLI )
      IF( NTLXLI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE CONVEXEINCONNUE'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C     LES TMS NSEF ET SOMMETS DE LA LIGNE CONVEXE FERMEE
      CALL LXTSOU( NTLXLI, 'NSEF', NTARLI, MNARLI )
      CALL LXTSOU( NTLXLI, 'XYZSOMMET', NTSOLI, MNSOLI )
      IF( NTARLI .LE. 0 .OR. NTSOLI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE CONVEXE FERMEE NON MAILLEE'
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C     LE NOMBRE DE SOMMETS DE LA LIGNE A DEPLACER
      NBSOML = MCN( MNSOLI + WNBSOM )
      IF( MOD(NBSOML,2) .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'NOMBRE IMPAIR ARETES DE LA LIGNE CONVEXE'
         CALL LEREUR
         IERR = 4
         RETURN
      ENDIF
C
C     LE NOMBRE DE TANGENTES DE LA LIGNE
      NBTGSL = MCN( MNSOLI + WNBTGS )
C
C     LES PARAMETRES DES NO SOMMET DE LA LIGNE A DEPLACER
      CALL NSEFPA( MCN(MNARLI),
     %             NUTYML, NBSOEL, NBSOEF, NBTGEL,
     %             LDAPEF, LDNGEF, LDTGEF, NBARLI,
     %             NX, NY, NZ,
     %             IERR   )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE CONVEXE FERMEE: MAILLAGE INCORRECT'
         CALL LEREUR
         IERR = 5
         RETURN
      ENDIF
C
C     LIGNE FERMEE?
      IF( MCN( MNARLI + WUTFMA ) .NE. 1 .OR.
     %    NBSOML .NE. NBARLI ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE CONVEXE FERMEE NON FERMEE'
         CALL LEREUR
         IERR = 6
         RETURN
      ENDIF
C
C     GENERATION DES SOMMETS ET DES TANGENTES DE LA SURFACE
C     =====================================================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CE SURFACE
C     LE NOMBRE DE SOMMETS
      NBSOM  = 1 + NBSOML
      NBCOOR = 3
      N = WYZSOM + NBCOOR * ( NBSOM + NBTGSL )
      CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'MOTS', N )
      CALL LXTSOU( NTLXSU, 'XYZSOMMET', NTSOFA, MNSOFA )
C
C     CALCUL DES COORDONNEES DU BARYCENTRE DES SOMMETS DU CONVEXE FERMEE
      DO 2 K=1,NBCOOR
         XYZB(K) = 0
 2    CONTINUE
      MNC = MNSOLI + WYZSOM
      DO 10 N=1,NBSOML
         DO 5 K=1,NBCOOR
            XYZB(K) = XYZB(K) + RMCN(MNC)
            MNC = MNC + 1
 5       CONTINUE
 10   CONTINUE
      DO 20 K=1,NBCOOR
         XYZB(K) = XYZB(K) / NBSOML
 20   CONTINUE
C
C     COPIE DES 3 COORDONNEES DES SOMMETS DE LA LIGNE CONVEXE FERMEE
      CALL TRTATA( RMCN(MNSOLI+WYZSOM),
     %             RMCN(MNSOFA+WYZSOM),
     %             NBCOOR * NBSOML )
      MNC = MNSOFA + WYZSOM + NBCOOR * NBSOML
      RMCN(MNC  ) = XYZB(1)
      RMCN(MNC+1) = XYZB(2)
      RMCN(MNC+2) = XYZB(3)
C
C     COPIE DES 3 COMPOSANTES DES TGS DE LA LIGNE CONVEXE FERMEE
      MNC = MNC + NBCOOR
      CALL TRTATA( RMCN(MNSOLI+WYZSOM+NBCOOR*NBSOML),
     %             RMCN(MNC),
     %             NBCOOR * NBTGSL )
C
C     NBSOM 'NOMBRE DE SOMMETS'
      MCN( MNSOFA + WNBSOM ) = NBSOM
C     NBCOOR 'NOMBRE DE COORDONNEES
C     MCN( MNSOFA + WBCOOR ) = NBCOOR
C     LE NOMBRE DE TANGENTES STOCKEES DE CETTE SURFACE
      MCN( MNSOFA + WNBTGS ) = NBTGSL
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOFA) )
C
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOFA + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     GENERATION DES NSEF DE LA SURFACE
C     =================================
C     CONSTRUCTION DU TABLEAU 'NSEF' DE CE SURFACE
      NBFASU = NBSOML / 2
      IF( NBTGSL .GT. 0 ) THEN
         NBTGEF = 8
         NBTGSU = NBTGSL
         NBEFAP = NBFASU
         NBEFTG = NBFASU
      ELSE
C        TANGENTE INCALCULABLE
         NBTGEF = 0
         NBTGSU = 0
         NBEFAP = 0
         NBEFTG = 0
      ENDIF
C
C     LE NOMBRE DE MOTS DU TMS 'NSEF'
      N = WUSOEF + 4*NBFASU + NBEFAP + NBEFTG * ( 1+NBTGEF )
      CALL LXTNDC( NTLXSU, 'NSEF', 'ENTIER', N )
      CALL LXTSOU( NTLXSU, 'NSEF', NTFASU, MNFASU )
      MNEFAP = MNFASU + WUSOEF + 4*NBFASU
      MNNTGS = MNEFAP + NBEFAP + NBEFTG
C
C     LA BOUCLE SUR LES COUPLES D'ARETES DU MAILLAGE DE LA LIGNE
      MNC = MNFASU + WUSOEF - 1
      MNT = MNNTGS - 1
      NBETG = 0
      DO 100 N=1,NBFASU
C
C           LE NUMERO DES NBSOEF SOMMETS DE L'ARETE 2*N-1
            CALL NSEFNS( 2*N-1 , NUTYML, NBSOEF, NBTGEL,
     %                   LDAPEF, LDNGEF, LDTGEF,
     %                   MNARLI, NX, NY, NZ,
     %                   NCOGEL, NUGEEF, NUEFTG, NOSOEL(1), IERR )
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
C
C           LE NUMERO DES NBSOEF SOMMETS DE L'ARETE 2*N
            CALL NSEFNS( 2*N,    NUTYML, NBSOEF, NBTGEL,
     %                   LDAPEF, LDNGEF, LDTGEF,
     %                   MNARLI, NX, NY, NZ,
     %                   NCOGEL, NUGEEF, NUEFTG, NOSOEL(5), IERR )
C           LE NUMERO DES 4 SOMMETS DE LA FACE
            MCN( MNC + 1 ) = NOSOEL(1)
            MCN( MNC + 2 ) = NOSOEL(2)
            MCN( MNC + 3 ) = NOSOEL(6)
            MCN( MNC + 4 ) = NBSOM
C
C           LES EVENTUELLES 8 TANGENTES DE LA FACE
            IF( NBTGSL .GT. 0 ) THEN
C
C              IL EXISTE DES EF A TG
               IF( NOSOEL(3) .NE. 0 .OR.
     %             NOSOEL(4) .NE. 0 .OR.
     %             NOSOEL(7) .NE. 0 .OR.
     %             NOSOEL(8) .NE. 0 ) THEN
C                 NOUVEAU EF A TG
                  NBETG = NBETG + 1
                  MCN(MNEFAP-1+N) = NBETG
C                 LE CODE GEOMETRIQUE 0  => 'C1 degre 3'
                  MCN(MNEFAP-1+NBEFAP+N) = 0
               ELSE
C                 EF NON A TG
                  MCN(MNEFAP-1+N) = 0
               ENDIF
               DO 30 K=1,8
                  MCN( MNT + K ) = 0
 30            CONTINUE
C              LA TANGENTE A LA COURBE AU SOMMET 1
               IF( NOSOEL(3) .NE. 0 ) THEN
                  MCN( MNT + 1 ) = NOSOEL(3)
               ENDIF
C              LA TANGENTE A LA COURBE AU SOMMET 2
               IF( NOSOEL(4) .NE. 0 ) THEN
                  MCN( MNT + 4 ) = NOSOEL(4)
               ENDIF
C              LA TANGENTE A LA COURBE AU SOMMET 2
               IF( NOSOEL(7) .NE. 0 ) THEN
                  MCN( MNT + 3 ) = NOSOEL(7)
               ENDIF
C              LA TANGENTE A LA COURBE AU SOMMET 3
               IF( NOSOEL(8) .NE. 0 ) THEN
                  MCN( MNT + 6 ) = NOSOEL(8)
               ENDIF
               MNT = MNT + 8
            ENDIF
C
C           LA FACE SUIVANTE
            MNC = MNC + 4
 100  CONTINUE
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
      MCN( MNFASU + WBEFTG ) = NBETG
C     NOMBRE DES EF AVEC POINTEUR SUR EF A TG
      MCN( MNFASU + WBEFAP ) = NBEFAP
C     NUMERO DU TYPE DE MAILLAGE : NON STRUCTURE
      MCN( MNFASU + WUTYMA ) = 0
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNFASU) )
C
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFASU + MOTVAR(6) ) = NONMTD ( '~>>>NSEF' )
C
 9990 RETURN
      END
