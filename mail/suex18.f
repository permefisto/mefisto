        SUBROUTINE SUEX18( NTLXSU, LADEFI, RADEFI,
     %                     NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE DE LA SURFACE D'UN 8EME OU DE LA SPHERE
C -----
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SPHERE
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
C          CF $MEFISTO/td/d/a_surface__definition
C RADEFI : TABLEAU REEL DE DEFINITION DE LA SURFACE
C          CF $MEFISTO/td/d/a_surface__definition
C
C SORTIES:
C --------
C NTNSEF : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES
C MNNSEF : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES
C          CF $MEFISTO/td/d/a___nsef
C NTXYZS : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNXYZS : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF $MEFISTO/td/d/a___xyzsommet
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS   : A. CUBIER ET F. DEVUYST DEA A.N. UPMC PARIS  DECEMBRE 1989
C MODIF TGS : A. PERRONNET  ANALYSE NUMERIQUE UPMC  PARIS      JUIN 1996
C234567--------------------------------------------------------------012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
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
C
      IERR = 0
C
      NCAS = LADEFI(WUTYSU)
C
C     PARAMETRES DU MAILLAGE
C     ======================
      NBTRSP = LADEFI(WBTRSP)
      RAYOSP = RADEFI(WAYOSP)
      NUPCSP = LADEFI(WUPCSP)
C
C     VERIFICATION DES PARAMETRES
C     ===========================
      IF ( NBTRSP.LE.0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:10),'(I10)') NBTRSP
         KERR(1) = 'NOMBRE APPROXIMATIF DE TRIANGLES INCORRECT ='
     %             // KERR(MXLGER)(1:10)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF ( RAYOSP.LE.0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:14),'(G14.6)') RAYOSP
         KERR(1) = 'RAYON DE LA SPHERE INCORRECT =' //
     %              KERR(MXLGER)(1:14)
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     RECUPERATION DES TROIS COORDONNEES DU CENTRE DE LA SPHERE
C     =========================================================
      CALL LXNLOU( NTPOIN,  NUPCSP  , NTLXCS, MN )
      CALL LXTSOU( NTLXCS, 'XYZSOMMET', NTCECS, MNCECS )
      MN = MNCECS + WYZSOM
      X = RMCN (MN     )
      Y = RMCN (MN + 1 )
      Z = RMCN (MN + 2 )
C
C     CALCUL NOMBRE EXACT DE TRIANGLES ET DE SOMMETS
C     ==============================================
      IF ( NCAS .EQ. 18 ) THEN
         NBTRSP = NBTRSP / 20
      ENDIF
C     CALCUL DU N EFFECTIF TEL QUE NBTRSP ~ N**2
      N = 1
  10  IF ( N**2 .LE. NBTRSP ) THEN
         N=N+1
         GOTO 10
      ENDIF
C
      IF ( ABS(NBTRSP-N**2).GT.ABS(NBTRSP-(N-1)**2) ) N = N - 1
      NBTRIA = N**2
C
      IF ( NCAS .EQ. 17 ) THEN
         NBSOMT = (N+1)*(N+2)/2
      ELSE
         NBSOMT = 10 * N**2 + 2
         NBTRIA = NBTRIA*20
      ENDIF
C
C     LE NOMBRE DE TANGENTES (IL EST POSSIBLE DE LE DIVISER PAR 2)
C     ======================
      NBTGS = 6 * NBTRIA
C
C     GENERATION DES SOMMETS
C     ======================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CETTE SURFACE
      CALL LXTNDC ( NTLXSU, 'XYZSOMMET', 'MOTS',
     %              WYZSOM+3*(NBSOMT+NBTGS) )
      CALL LXTSOU ( NTLXSU, 'XYZSOMMET', NTXYZS, MNXYZS )
C
C     GENERATION DU TABLEAU NSEF
C     ==========================
C     CONSTRUCTION DU TABLEAU 'NSEF' DE CETTE SURFACE
      CALL LXTNDC( NTLXSU, 'NSEF', 'MOTS', WUSOEF+14*NBTRIA )
      CALL LXTSOU( NTLXSU, 'NSEF',  NTNSEF, MNNSEF )
C
C     GENERATION DES SOMMETS ET TRIANGLES DE LA SURFACE
C     =================================================
      IF( NCAS .EQ. 17 ) THEN
C
C        1/8 DE SPHERE
         CALL SPH1S8( N, RAYOSP, X, Y, Z,
     %                MCN(MNXYZS+WYZSOM),
     %                MCN(MNNSEF+WUSOEF) )
C
      ELSE IF ( NCAS .EQ. 18 ) THEN
C
C        LA SPHERE COMPLETE
         NBNOLO = (N+1) * (2*N+1) + 1
         NBNOTB = 3 * N + 2
         CALL TNMCDC( 'ENTIER', NBNOLO, MNTATL )
         CALL TNMCDC( 'ENTIER', NBNOLO, MNTATO )
         CALL TNMCDC( 'ENTIER', NBNOTB, MNTATG )
         CALL TNMCDC( 'ENTIER', NBNOTB, MNTATD )
         CALL SPH1S1( N, RAYOSP, X, Y, Z,
     %                MCN(MNTATL), MCN(MNTATO),
     %                MCN(MNTATG), MCN(MNTATD),
     %                NBSOMT, NBTRIA,
     %                MCN(MNXYZS+WYZSOM),
     %                MCN(MNNSEF+WUSOEF) )
         CALL TNMCDS( 'ENTIER', NBNOLO, MNTATL )
         CALL TNMCDS( 'ENTIER', NBNOLO, MNTATO )
         CALL TNMCDS( 'ENTIER', NBNOTB, MNTATG )
         CALL TNMCDS( 'ENTIER', NBNOTB, MNTATD )
C
      ENDIF
C
C     MISE A JOUR DU TABLEAU XYZSOMMET DE CETTE SURFACE
C     =================================================
C     NBSOM 'NOMBRE DE SOMMETS'
      MCN( MNXYZS + WNBSOM ) = NBSOMT
C
C     LE NOMBRE DE TANGENTES DU MAILLAGE
      MCN( MNXYZS + WNBTGS ) = NBTGS
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNXYZS + WBCOOR ) = 3
C
C     CALCUL DES XYZ DES 8 TANGENTES ( + 2 FOIS 0 ) DES NBTRIA TRIANGLES
C     ------------------------------------------------------------------
      CALL TGSPHE( RAYOSP, X, Y, Z, NBTRIA, MCN(MNNSEF+WUSOEF),
     %             MCN(MNXYZS+WYZSOM), MCN(MNXYZS+WYZSOM+3*NBSOMT) )
C
C     LE NOM DU TABLEAU DESCRIPTEUR
      MCN( MNXYZS + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNXYZS) )
C
C     MISE A JOUR DU TABLEAU NSEF DE CETTE SURFACE
C     ============================================
C     TYPE DE L'OBJET : SURFACE
      MCN ( MNNSEF + WUTYOB ) = 3
C
C     LE TYPE DE FERMETURE DU MAILLAGE
      IF( NCAS .EQ. 17 ) THEN
C        1/8 SPHERE => SURFACE NON FERMEE
         MCN( MNNSEF + WUTFMA ) = 0
      ELSE
C        1 SPHERE => SURFACE FERMEE
         MCN( MNNSEF + WUTFMA ) = 1
      ENDIF
C
C     8 TANGENTES PAR EF
      MCN( MNNSEF + WBTGEF ) = 8
C     TOUS LES TRIANGLES POINTENT SUR UN EF A TG
      MCN( MNNSEF + WBEFAP ) = NBTRIA
C     TOUS LES TRIANGLES SONT DES EF A TG
      MCN( MNNSEF + WBEFTG ) = NBTRIA
C
C     NUMERO DU TYPE DE MAILLAGE : NON STRUCTURE
      MCN ( MNNSEF + WUTYMA ) = 0
C
C     NBSOEF 'NOMBRE DE SOMMETS PAR NSEF'
      MCN ( MNNSEF + WBSOEF ) = 4
C
C     NBEFOB 'NOMBRE D'EF DE L'OBJET'
      MCN ( MNNSEF + WBEFOB ) = NBTRIA
C
C     LE POINTEUR SUR LES EF A TG
      MN = MNNSEF + WUSOEF + 4*NBTRIA - 1
      DO 30 I=1,NBTRIA
         MCN(MN+I) = I
 30   CONTINUE
      MN = MN + NBTRIA
C
C     LE CODE GEOMETRIQUE DE L'EF: ICI SPHERE => 11
      DO 40 I=1,NBTRIA
         MCN(MN+I) = 11
 40   CONTINUE
      MN = MN + NBTRIA
C
C     LE NUMERO DES 8 TANGENTES DES NBTRIA TRIANGLES
      K = 0
      DO 60 I=1,NBTRIA
         DO 50 J=1,6
            K = K + 1
            MCN(MN+J) = K
 50      CONTINUE
         MCN(MN+7) = 0
         MCN(MN+8) = 0
         MN = MN + 8
 60   CONTINUE
C
C     LA DATE DE CREATION
      CALL ECDATE ( MCN( MNNSEF ) )
C
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN ( MNNSEF + MOTVAR(6) ) = NONMTD ( '~>>>NSEF' )
C
C     REDUCTION DU NOMBRE DE TANGENTES
C     ================================
      CALL MOINTG( NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
CCCC
CCCC     PROJECTION DES TGS SUR LE PLAN A DISTANCE MINIMALE
CCCC     INUTILE CAR REDONNE LES MEMES TANGENTES => C'EST RASSURANT!
CCCC     ===========================================================
CCC      CALL SUC0G1( NTLXSU, 0,      I,
CCC     %             NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
      RETURN
      END
