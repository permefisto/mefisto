      SUBROUTINE NEFMEF2D( KNOMOB, NTLXOB, NBSURFACES,
     %                     NBSOMMET, MNXYZS,
     %                     NBEF, NBNOEF, MNNOST, NBST1EF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CONSTRUIRE LES TMS XYZSOMMET et NSEF DES PLSV DE L'OBJET
C ----- CREATION DES VOLUMES, SURFACES, LIGNES et POINTS MEFISTO
C       de L'OBJET NEF deja lu SUR LE FICHIER NomObjet.XyzNoEF
C
C ENTREES:
C --------
C KNOMOB    : NOM DE L'OBJET
C NTLXOB    : NO DU TMS DE L'OBJET
C NBSURFACES: DE L'OBJET
C NBSOMMET  : DES EF DES VOLUMES DE L'OBJET
C MNXYZS    : POINTEUR MCN SUR XYZ1,...XYZnbsommet DES SOMMETS DES EF
C NBEF      : NOMBRE D'EF DES VOLUMES DE L'OBJET
C NBNOEF    : NOMBRE DE NUMEROS PAR EF
C MNNOST    : POINTEUR MCN SUR LES NOS DES EF
C NBST1EF   : NOMBRE DE SOMMETS POUR 1 EF (3 POUR UN TRIANGLE)
C
C SORTIES:
C --------
C          PLSV DES OBJETS SONT CREES DANS LEURS LEXIQUES RESPECTIFS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE LJLL UPMC PARIS   JUIN 2004
C2345X7..............................................................012
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
c
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a_point__definition.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/xyzext.inc"
C
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     KNOMOB
      CHARACTER*24      KNOM
      INTEGER           NGS(2)
C
C     NOMBRE D'EF DE CHACUNE DES NBSURFACES
      MNFAF1 = 0
      NEWOLD = 0
      MNBTEV = 0
      CALL TNMCDC( 'ENTIER', NBSURFACES, MNBTEV )
      CALL AZEROI( NBSURFACES, MCN(MNBTEV) )
      MNBTE1 = MNBTEV - 1
C
      MNNEF  = MNNOST - 1 - NBNOEF
      DO 10 N = 1, NBEF
         MNNEF = MNNEF + NBNOEF
C        NO DE SURFACE
         NSF = MCN(MNNEF+NBNOEF)
C        UN EF DE PLUS POUR CETTE SURFACE NSF
         MCN(MNBTE1+NSF) = MCN(MNBTE1+NSF) + 1
 10   CONTINUE
C
C     CONSTRUCTION DES NBSURFACES
      IF( NBSURFACES .LE. 0 ) GOTO 222

C     NO DU TMS DES SURFACES, NO DU TMS NSEF, ADRESSE MCN DU TMS NSEF, NB EF
      CALL TNMCDC( 'ENTIER', NBSURFACES, NTSFOBS )
      CALL TNMCDC( 'ENTIER', NBSURFACES, NTNSEFS )
      CALL TNMCDC( 'ENTIER', NBSURFACES, MNNSEFS )
      CALL TNMCDC( 'ENTIER', NBSURFACES, MNNOLXS )
C
      NBSOEFV = 4
      DO 50 NSF=1,NBSURFACES
C
         N = NBCHIF( NSF )
         KERR(1) = '(I )'
         WRITE( KERR(1)(3:3), '(I1)' ) N
         KNOM = 'SURF   '
         WRITE( KNOM(5:4+N), KERR(1)(1:4) ) NSF
C
         CALL LXLXDC( NTSURF, KNOM, 24, 8 )
         CALL LXLXOU( NTSURF, KNOM, NTSFOB, MNSFOB )
         MCN(NTSFOBS-1+NSF) = NTSFOB
C
C        RECHERCHE DU NUMERO DU NOM KNOM DANS LE LEXIQUE DES SURFACES
         CALL LXNMNO( NTSURF, KNOM, MCN(MNNOLXS-1+NSF), MNLX )
C
C        LE TMS DEFINITION
         CALL LXTNDC( NTSFOB, 'DEFINITION', 'MOTS', WUTSSV+1 )
         CALL LXTSOU( NTSFOB, 'DEFINITION', NTDFOB, MNDFOB )
C
C        LE TMS NSEF
         NBTRIA = MCN(MNBTE1+NSF)
         CALL LXTNDC( NTSFOB,'NSEF','MOTS', WUSOEF+NBSOEFV*NBTRIA )
         CALL LXTSOU( NTSFOB,'NSEF',
     %                MCN(NTNSEFS-1+NSF), MCN(MNNSEFS-1+NSF) )
C
C        SURFACE DEFINIE PAR SES TMS XYZSOMMET NSEF
         MCN(MNDFOB+WTYTRV) = 1
         MCN(MNDFOB+WUTYVO) = 10
C         MCN(MNDFOB+WUTSOV) = calcule plus tard
         MCN(MNDFOB+WUTSSV) = MCN(NTNSEFS-1+NSF)
C
C        LA DATE DU TMS DEFINITION DE LA SURFACE
         CALL ECDATE( MCN(MNDFOB) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNDFOB + MOTVAR(6) ) = NONMTD( '~>SURFACE>>DEFINITION' )
 50   CONTINUE
C
C     LE NUMERO DES SOMMETS DES TETRAEDRES ET TRIANGLES DANS NSEF
      CALL AZEROI( NBSURFACES, MCN(MNBTEV) )
C     HACHAGE DES ARETES DES TRIANGLES DES FACES
C     L1ARET : NOMBRE DE MOTS PAR ARETE DU TABLEAU NARET
C     L2ARET : NOMBRE DE ARETES DU TABLEAU NARET
C     MNARET : ADRESSE DANS M DU TABLEAU NARET DES ARETES DU MAILLAGE
C              EN SORTIE NARET(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C                        NARET(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C                        NARET(3,I)= CHAINAGE HACHAGE SUR ARETE SUIVANTE
C                        NARET(4,I)= NUMERO 1-ERE FACE CONTENANT CETTE ARETE
C                        NARET(5,I)= NUMERO 2-EME FACE CONTENANT CETTE ARETE
      L1ARET = 5
      L2ARET = 3 * NBEF / 2 + 1024
      L1     = L1ARET * L2ARET
      MNARET = 0
      CALL TNMCDC( 'ENTIER', L1, MNARET )
C
C     LE TABLEAU DES ARETES EST INITIALISE A ZERO
      CALL AZEROI( L1, MCN(MNARET) )
C
C     LA 1-ERE ARETE LIBRE EST LA DERNIERE DU TABLEAU
      LIBREF = L2ARET
C
      MNNEF = MNNOST - 1 - NBNOEF
      DO 80 N = 1, NBEF
         MNNEF = MNNEF + NBNOEF
C        NO DE SURFACE
         NSF = MCN(MNNEF+NBNOEF)
C
C        LE TRIANGLE: SES 3 SOMMETS ET 0
         NEF = MCN(MNBTEV-1+NSF)
         MN  = MCN(MNNSEFS-1+NSF) + WUSOEF + NEF * NBSOEFV - 1
         DO 62 I=1,NBST1EF
            MCN( MN + I ) = MCN( MNNEF + I )
 62      CONTINUE
         DO 63 I=NBST1EF+1,NBSOEFV
            MCN( MN + I ) = 0
 63      CONTINUE
         MCN(MNBTEV-1+NSF) = MCN(MNBTEV-1+NSF) + 1
C
C        LES ARETES DE L'EF
         NSA1 = MCN( MN + NBST1EF )
         DO 70 J=1,3
C           LE NO DES 2 SOMMETS DE L'ARETE
            NSA2 = MCN( MN + J )
C           NS1<NS2 DE L'ARETE ET NO DE LA FACE
            IF( NSA1 .LT. NSA2 ) THEN
               NGS(1) = NSA1
               NGS(2) = NSA2
            ELSE
               NGS(1) = NSA2
               NGS(2) = NSA1
            ENDIF
            CALL HACHAG( 2,NGS,L1ARET,L2ARET,MCN(MNARET),3,
     &                   LIBREF,NOAR )
C           L'ARETE RETROUVEE OU AJOUTEE
            MNA = MNARET + L1ARET * ( ABS(NOAR) - 1 )
C           STOCKAGE DU NO EF N DANS NARET(4,...;NOAR)
            L1 = 2
 67         L1 = L1 + 1
            IF( L1 .LT. L1ARET ) THEN
               IF( MCN(MNA + L1) .NE. 0 ) GOTO 67
C              NO SURFACE LIBRE RECOIT LE NO DE SURFACE DU TRIANGLE
               MCN( MNA + L1 ) = NSF
            ELSE
               WRITE(IMPRIM,10066) (NGS(L1),L1=1,2),
     %                             (MCN(MNA+L1),L1=3,4)
10066          FORMAT(' ARETE DE SOMMETS ',2I9,
     %                ' DANS PLUS DES 2 FACES',2I6 )
            ENDIF
            NSA1 = NSA2
 70      CONTINUE
C
 80   CONTINUE
C
C     RENUMEROTATION CONTINUE DES SOMMETS DES SURFACES
C     NO ANCIEN DES SOMMETS
      CALL TNMCDC( 'ENTIER', 1+NBSOMMET, NEWOLD )
C
      DO 100 NSF=1,NBSURFACES
C
C        MISE A JOUR DU TMS NSEF D'ADRESSE MCN
         MNSEF = MCN(MNNSEFS-1+NSF)
C
C        variable NUTYOB        'numero de type de l''objet'    entier
C        ( 1 : 'point' , 2 : 'ligne' , 3 : 'surface' , 4 : 'volume' )  ;
         MCN( MNSEF + WUTYOB ) = 3
C
C        variable NUTFMA 'Ligne ou Surface fermee ou non-ferme (P et V sont inco
C        ( -1 : 'inconnu' , 0 : 'non-ferme' , 1 : 'ligne ou surface fermee' ) ;
         MCN( MNSEF + WUTFMA ) = -1
C
C        variable NBSOEF 'nombre de sommets par EF' entier
C        ( 1 : 'sommet' , 2 : 'arete' , 4 : 'face' , 8 : 'cube' ) ;
         MCN( MNSEF + WBSOEF ) = NBSOEFV
C
C        variable NBTGEF 'nombre de tangentes par EF' entier
C        ( 0 : '0 tg par EF' , 1 : '1 tg par sommet' , 2 : '2 tg par arete' ,
C          8 : '8 tg par face' ,  24 : '24 tg par cube' )  ;
         MCN( MNSEF + WBTGEF ) = 0
C
C        variable NBEFOB 'nombre des EF du PLSV'  entier ;
         NBEFOB = MCN(MNBTEV-1+NSF)
         MCN( MNSEF + WBEFOB ) = NBEFOB
         MCN( MNSEF + WBEFTG ) = 0
         MCN( MNSEF + WBEFAP ) = 0
C
C        variable NUTYMA 'numero de type du maillage' entier
C        ( 0 : 'non structure' , ... )
         MCN( MNSEF + WUTYMA ) = 0
C
         DO 82 N=0,NBSOMMET
            MCN(NEWOLD+N) = 0
 82      CONTINUE
C
C        MARQUAGE DES SOMMETS UTILES DE LA SURFACE NSF
         MN = MNSEF + WUSOEF - 1
         DO 86 N=1,NBEFOB
            DO 84 I=1,NBST1EF
               NS = MCN( MN + I )
               MCN(NEWOLD+NS) = 1
 84         CONTINUE
            MN = MN + NBSOEFV
 86      CONTINUE
C
C        NOUVEAU NUMERO DES SOMMETS
         MN  = MNSEF + WUSOEF - 1
         NBS = 0
         DO 90 NS=1,NBSOMMET
            IF( MCN(NEWOLD+NS) .GT. 0 ) THEN
               NBS = NBS + 1
               MCN(NEWOLD+NS) = NBS
            ENDIF
 90      CONTINUE
C
C        CREATION DU TABLEAU XYZSOMMET DE LA SURFACE NSF
         NTSFOB = MCN(NTSFOBS-1+NSF)
C
C        LE TMS XYZSOMMET DE LA SURFACE NSF
         CALL LXTNDC( NTSFOB, 'XYZSOMMET', 'MOTS', WYZSOM+NBS*3 )
         CALL LXTSOU( NTSFOB, 'XYZSOMMET', NTXYZ, MNXYZ )
         MCN( MNXYZ + WBCOOR ) = 3
         MCN( MNXYZ + WNBSOM ) = NBS
         MCN( MNXYZ + WNBTGS ) = 0
         MN0 = MNXYZS - 3
         MN1 = MNXYZ + WYZSOM - 3
         DO 92 NS=1,NBSOMMET
            N = MCN(NEWOLD+NS)
            IF( N .GT. 0 ) THEN
               M0 = MN0 + 3*NS
               M1 = MN1 + 3*N
               RMCN( M1     ) = RMCN( M0     )
               RMCN( M1 + 1 ) = RMCN( M0 + 1 )
               RMCN( M1 + 2 ) = RMCN( M0 + 2 )
            ENDIF
 92      CONTINUE
C        LA DATE DU TMS XYZSOMMET DE LA SURFACE NSF
         CALL ECDATE( MCN(MNXYZ) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNXYZ + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C        MISE A JOUR DU NO DE SOMMET DANS NSEF DE LA SURFACE NSF
         MN = MNSEF + WUSOEF - 1
         DO 96 N=1,NBEFOB
            DO 94 I=1,4
               NS  = MCN( MN + I )
               NEW = MCN(NEWOLD+NS)
               MCN( MN + I ) = NEW
 94         CONTINUE
            MN = MN + NBSOEFV
 96      CONTINUE
C
C        LA DATE DU TMS NSEF DE LA SURFACE NSF
         CALL ECDATE( MCN(MNSEF) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNSEF + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C        ON COMPLETE LE NO DE TMS XYZSOMMET ET NSEF DANS LA DEFINITION DE LA SUR
         CALL LXTSOU( NTSFOB, 'DEFINITION', NTDFOBS, MNDFOBS )
         MCN(MNDFOBS+WUTSOV) = NTXYZ
C        MCN(MNDFOBS+WUTSSV) = MCN(NTNSEFS-1+NSF)  deja fait
C
 100  CONTINUE
C
C     COMPTAGE DES LIGNES=ARETES APPARTENANT A AU MOINS 2 SURFACES
      CALL TNMCDC( 'ENTIER', (NBSURFACES+1)**2, MNFAFA )
      CALL AZEROI( (NBSURFACES+1)**2, MCN(MNFAFA) )
      MNFAF1 = MNFAFA - 1
      MNA    = MNARET - L1ARET
      DO 220 N=1,L2ARET
         MNA  = MNA + L1ARET
         NSA1 = MCN( MNA )
         IF( NSA1 .EQ. 0 ) GOTO 220
         NF1  = MCN( MNA+3 )
         NF2  = MCN( MNA+4 )
         IF( NF1 .EQ. NF2 ) THEN
C           ARETE DOUBLE D'UNE MEME FACE => ELIMINEE
            MCN( MNA ) = 0
         ELSE
C           ARETE DOUBLE SUR 2 FACES DE LA LIGNE NF1<NF2
            IF( NF1 .GT. NF2 ) THEN
               NF  = NF1
               NF1 = NF2
               NF2 = NF
               MCN( MNA+3 ) = NF1
               MCN( MNA+4 ) = NF2
            ENDIF
C           LA LIGNE NF1-NF2 EST MARQUEE EXISTANTE
            MN = MNFAF1 + NF1 + NBSURFACES*NF2
            MCN( MN ) = MCN( MN ) + 1
         ENDIF
 220  CONTINUE
C
C     NOMBRE DE LIGNES
 222  NBLIGNES = 0
      DO 240 NF2=1,NBSURFACES
         DO 230 NF1=0,NF2-1
            MN = MNFAF1 + NF1 + NBSURFACES*NF2
            IF( MCN(MN) .GT. 0 ) THEN
               NBLIGNES = NBLIGNES + 1
            ENDIF
 230     CONTINUE
 240  CONTINUE
C
C     CONSTRUCTION DES NBLIGNES LIGNES
C     --------------------------------
      IF( NBLIGNES .LE. 0 ) GOTO 301

C     NO DU TMS DES LIGNES, NO DU TMS NSEF, ADRESSE MCN DU TMS NSEF, NB EF
      CALL TNMCDC( 'ENTIER', NBLIGNES, NTLIOBL )
      CALL TNMCDC( 'ENTIER', NBLIGNES, NTNSEFL )
      CALL TNMCDC( 'ENTIER', NBLIGNES, MNNSEFL )
      CALL TNMCDC( 'ENTIER', NBLIGNES, MNNOLXL )
C
      NBSOEFL = 2
      NL = 0
      DO 260 NF2=1,NBSURFACES
         DO 250 NF1=0,NF2-1
            MN = MNFAF1 + NF1 + NBSURFACES * NF2
            NBARET = MCN( MN )
            IF( NBARET .GT. 0 ) THEN
               NL = NL + 1
C              ATTENTION: LE NOMBRE D'ARETES DEVIENT LE NO DE LA LIGNE
               MCN( MN ) = NL
C
               N = NBCHIF( NL )
               KERR(1) = '(I )'
               WRITE( KERR(1)(3:3), '(I1)' ) N
               KNOM = 'L     '
               WRITE( KNOM(2:1+N), KERR(1)(1:4) ) NL
               M = N + 2
C
               IF( NF1 .NE. 0 ) THEN
                  KNOM(M:M+1) = '_S'
                  M = M + 2
                  N = NBCHIF( NF1 )
                  KERR(1) = '(I )'
                  WRITE( KERR(1)(3:3), '(I1)' ) N
                  WRITE( KNOM(M:M-1+N), KERR(1)(1:4) ) NF1
                  M = M + N
               ENDIF
C
               IF( NF2 .NE. 0 ) THEN
                  KNOM(M:M+1) = '_S'
                  M = M + 2
                  N = NBCHIF( NF2 )
                  KERR(1) = '(I )'
                  WRITE( KERR(1)(3:3), '(I1)' ) N
                  WRITE( KNOM(M:M-1+N), KERR(1)(1:4) ) NF2
               ENDIF
C
C              LEXIQUE DE LA LIGNE
               CALL LXLXDC( NTLIGN, KNOM, 24, 8 )
               CALL LXLXOU( NTLIGN, KNOM, NTLGOB, MNLGOB )
               MCN(NTLIOBL-1+NL) = NTLGOB
C
C              RECHERCHE DU NUMERO DU NOM KNOM DANS LE LEXIQUE DES LIGNES
               CALL LXNMNO( NTLIGN, KNOM, MCN(MNNOLXL-1+NL), MNLX )
C
C              LE TMS DEFINITION DE LA LIGNE
               CALL LXTNDC( NTLGOB, 'DEFINITION', 'MOTS', WUTSSL+1 )
               CALL LXTSOU( NTLGOB, 'DEFINITION', NTDFOB, MNDFOB )
C
C              LE TMS NSEF DE LA LIGNE NF1-NF2
               CALL LXTNDC( NTLGOB,'NSEF','MOTS',WUSOEF+NBSOEFL*NBARET )
               CALL LXTSOU( NTLGOB,'NSEF', NTSEF, MNSEF )
               MCN(NTNSEFL-1+NL) = NTSEF
               MCN(MNNSEFL-1+NL) = MNSEF
C
C              LIGNE DEFINIE PAR SES TMS XYZSOMMET NSEF
               MCN(MNDFOB+WTYTRL) = 1
               MCN(MNDFOB+WUTYLI) = 10
C               MCN(MNDFOB+WUTSOL) = calcule plus tard
               MCN(MNDFOB+WUTSSL) = MCN(NTNSEFL-1+NL)
C
C              LA DATE DU TMS DEFINITION DE LA LIGNE
               CALL ECDATE( MCN(MNDFOB) )
C              LE NUMERO DU TABLEAU DESCRIPTEUR
               MCN( MNDFOB + MOTVAR(6) ) = NONMTD('~>LIGNE>>DEFINITION')
            ENDIF
 250     CONTINUE
 260  CONTINUE
C
C     CONSTRUCTION DANS NSEF DES ARETES DES LIGNES
      CALL TNMCDC( 'ENTIER', NBLIGNES, MNBARL )
      CALL AZEROI( NBLIGNES, MCN(MNBARL) )
      MNA = MNARET - L1ARET
      DO 280 N=1,L2ARET
         MNA  = MNA + L1ARET
         NSA1 = MCN(MNA)
         IF( NSA1 .GT. 0 ) THEN
C           ARETE APPARTENANT A AU MOINS 2 FACES
            NSA2 = MCN( MNA + 1 )
            NF1  = MCN(MNA+3)
            NF2  = MCN(MNA+4)
C           NUMERO DE LA LIGNE
            MN  = MNFAF1 + NF1 + NBSURFACES*NF2
            NL  = MCN( MN )
C           NO ACTUEL DE LA DERNIERE ARETE DE CETTE LIGNE
            NEF = MCN(MNBARL-1+NL)
            MN  = MCN(MNNSEFL-1+NL) + WUSOEF + NEF * NBSOEFL
C           STOCKAGE DANS NSEF DE LA LIGNE
            MCN( MN     ) = NSA1
            MCN( MN + 1 ) = NSA2
C           UNE ARETE DE PLUS POUR CETTE LIGNE
            MCN(MNBARL-1+NL) = NEF + 1
         ENDIF
 280  CONTINUE
C
C     RENUMEROTATION CONTINUE DES SOMMETS DES LIGNES
      CALL TNMCDC( 'ENTIER', NBLIGNES, MNLXYZ )
      DO 300 NL=1,NBLIGNES
C
C        MISE A JOUR DU TMS NSEF D'ADRESSE MCN DE LA LIGNE NL
         MNSEF = MCN(MNNSEFL-1+NL)
C        variable NUTYOB  'numero de type de l''objet' 2 : 'ligne'
         MCN( MNSEF + WUTYOB ) = 2
C
C        variable NUTFMA 'Ligne ou Surface fermee ou non-ferme (P et V sont inco
C        ( -1 : 'inconnu' , 0 : 'non-ferme' , 1 : 'ligne ou surface fermee' ) ;
         MCN( MNSEF + WUTFMA ) = -1
C
C        variable NBSOEF 'nombre de sommets par EF' 2 : 'arete'
         MCN( MNSEF + WBSOEF ) = NBSOEFL
C
C        variable NBTGEF 'nombre de tangentes par EF' entier
C        ( 0 : '0 tg par EF' , 1 : '1 tg par sommet' , 2 : '2 tg par arete' ,
C          8 : '8 tg par face' ,  24 : '24 tg par cube' )  ;
         MCN( MNSEF + WBTGEF ) = 0
C
C        variable NBEFOB 'nombre des EF du PLSV'  entier ;
         NBEFOB = MCN(MNBARL-1+NL)
         MCN( MNSEF + WBEFOB ) = NBEFOB
         MCN( MNSEF + WBEFTG ) = 0
         MCN( MNSEF + WBEFAP ) = 0
C
C        variable NUTYMA numero de type du maillage 0 : 'non structure'
         MCN( MNSEF + WUTYMA ) = 0
C
         DO 282 N=0,NBSOMMET
            MCN(NEWOLD+N) = 0
 282     CONTINUE
C
C        MARQUAGE DES SOMMETS UTILES DE LA LIGNE NL
         MN = MNSEF + WUSOEF - 1
         DO 286 N=1,NBEFOB
            DO 284 I=1,NBSOEFL
               NS = MCN( MN + I )
               MCN(NEWOLD+NS) = 1
 284        CONTINUE
            MN = MN + NBSOEFL
 286     CONTINUE
C
C        NOUVEAU NUMERO DES SOMMETS
         MN  = MNSEF + WUSOEF - 1
         NBS = 0
         DO 290 NS=1,NBSOMMET
            IF( MCN(NEWOLD+NS) .GT. 0 ) THEN
               NBS = NBS + 1
               MCN(NEWOLD+NS) = NBS
            ENDIF
 290     CONTINUE
C
C        CREATION DU TABLEAU XYZSOMMET DE LA LIGNE NL
         NTLGOB = MCN(NTLIOBL-1+NL)
C
C        LE TMS XYZSOMMET
         CALL LXTNDC( NTLGOB, 'XYZSOMMET', 'MOTS', WYZSOM+NBS*3 )
         CALL LXTSOU( NTLGOB, 'XYZSOMMET', NTXYZ, MNXYZ )
         MCN(MNLXYZ-1+NL) = MNXYZ
         MCN( MNXYZ + WBCOOR ) = 3
         MCN( MNXYZ + WNBSOM ) = NBS
         MCN( MNXYZ + WNBTGS ) = 0
         MN0 = MNXYZS - 3
         MN1 = MNXYZ + WYZSOM - 3
         DO 292 NS=1,NBSOMMET
            N = MCN(NEWOLD+NS)
            IF( N .GT. 0 ) THEN
               M0 = MN0 + 3*NS
               M1 = MN1 + 3*N
               RMCN( M1     ) = RMCN( M0     )
               RMCN( M1 + 1 ) = RMCN( M0 + 1 )
               RMCN( M1 + 2 ) = RMCN( M0 + 2 )
            ENDIF
 292     CONTINUE
C        LA DATE DU TMS XYZSOMMET DE LA LIGNE NL
         CALL ECDATE( MCN(MNXYZ) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNXYZ + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C        MISE A JOUR DU NO DE SOMMET DANS NSEF DE LA LIGNE NL
         MN = MNSEF + WUSOEF - 1
         DO 296 N=1,NBEFOB
            DO 294 I=1,NBSOEFL
               NS  = MCN( MN + I )
               NEW = MCN(NEWOLD+NS)
               MCN( MN + I ) = NEW
 294        CONTINUE
            MN = MN + NBSOEFL
 296     CONTINUE
C
C        LA DATE DU TMS NSEF DE LA LIGNE NL
         CALL ECDATE( MCN(MNSEF) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNSEF + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C        ON COMPLETE LE NO DE TMS XYZSOMMET ET NSEF DANS LA DEFINITION DE LA LIG
         CALL LXTSOU( NTLGOB, 'DEFINITION', NTDFOBS, MNDFOBS )
         MCN(MNDFOBS+WUTSOL) = NTXYZ
C        MCN(MNDFOBS+WUTSSL) = MCN(NTNSEFL-1+NL)  deja fait
C
 300  CONTINUE
C
C     CONSTRUCTION DES POINTS=SOMMETS EXTREMITES D'AU MOINS 2 LIGNES
 301  NBPOINTS = 0
      MXPOINTS = 6 * NBLIGNES
      IF( MXPOINTS .LE. 0 ) GOTO 400
C
C     NO DU TMS DES LIGNES, NO DU TMS NSEF, ADRESSE MCN DU TMS NSEF, NB EF
      CALL TNMCDC( 'REEL', 3*MXPOINTS, MNSTL )
      CALL TNMCDC( 'ENTIER', MXPOINTS, MNBLGST )
      CALL AZEROI( MXPOINTS, MCN(MNBLGST) )
      CALL TNMCDC( 'ENTIER', 2*MXPOINTS, MNNOLGP )
      CALL AZEROI( 2*MXPOINTS, MCN(MNNOLGP) )
C
C     CONSTRUCTION DES XYZ DES SOMMETS EXTREMITES DES LIGNES AVEC IDENTIFICATION
      NBPOINTS = 0
      DO 350 NL=1,NBLIGNES
C
C        RECHERCHE DES EXTREMITES SIMPLES DE LA LIGNE NL
C        QUI PEUVENT ETRE DIFFERENTES DU NO 1 ET NBSOM!
         MNXYZ  = MCN( MNLXYZ-1+NL )
         NBSOM  = MCN( MNXYZ + WNBSOM )
         MNXYZ1 = 0
         MNNBAS = 0
         CALL TNMCDC( 'ENTIER', NBSOM, MNNBAS )
         CALL AZEROI( NBSOM, MCN(MNNBAS) )
C
         MNSEF = MCN( MNNSEFL-1+NL )
         NBEFL = MCN( MNSEF + WBEFOB )
         MN    = MNSEF + WUSOEF
         DO 305 J=1,NBEFL
C           LE NO DES 2 SOMMETS DE L'ARETE J DE LA LIGNE NL
            NSA1 = MCN( MN   )
            NSA2 = MCN( MN+1 )
            MN   = MN + 2
C           NOMBRE D'ARETES AUXQUELLES LE SOMMET APPARTIENT
            MCN( MNNBAS-1+NSA1 ) = MCN( MNNBAS-1+NSA1 ) + 1
            MCN( MNNBAS-1+NSA2 ) = MCN( MNNBAS-1+NSA2 ) + 1
 305     CONTINUE
C
         NBS1 = 0
         DO 308 J=1,NBSOM
            IF( MCN( MNNBAS-1+J ) .EQ. 1 ) NBS1 = NBS1 + 1
 308     CONTINUE
         IF( NBS1 .LE. 0 ) GOTO 345
C
C        TABLEAU DES 3 COORDONNEES DES SOMMETS SIMPLES
         CALL TNMCDC( 'REEL', 3*NBS1, MNXYZS1 )
         MNXYZ0 = MNXYZS1
         DO 310 J=1,NBSOM
            IF( MCN( MNNBAS-1+J ) .EQ. 1 ) THEN
               MN1 = MCN( MNLXYZ-1+NL ) + WYZSOM - 3 + 3 * J
               RMCN( MNXYZ0   ) = RMCN( MN1   )
               RMCN( MNXYZ0+1 ) = RMCN( MN1+1 )
               RMCN( MNXYZ0+2 ) = RMCN( MN1+2 )
               MNXYZ0 = MNXYZ0 + 3
            ENDIF
 310     CONTINUE
C
C        IDENTIFICATION DES SOMMETS SIMPLES UNIQUES
         MNXYZ0 = MNXYZS1 - 3
         DO 340 J=1,NBS1
            MNXYZ0 = MNXYZ0 + 3
            MN     = MNSTL
            DO 330 K=1,NBPOINTS
C              IDENTIFICATION DES 2 SOMMETS?
               CALL XYZIDE( RMCN(MNXYZ0), RMCN(MN), IDENT )
               IF( IDENT .NE. 0 ) THEN
C                 SOMMET DEJA IDENTIFIE. COMPTEUR+1
                  MCN( MNBLGST-1+K ) = MCN( MNBLGST-1+K ) + 1
                  MCN( MNNOLGP-1+K*2 ) = NL
                  GOTO 340
               ENDIF
               MN = MN + 3
 330        CONTINUE
C
C           NOUVEAU SOMMET A AJOUTER AU TABLEAU MNSTL
            RMCN( MN   ) = RMCN( MNXYZ0   )
            RMCN( MN+1 ) = RMCN( MNXYZ0+1 )
            RMCN( MN+2 ) = RMCN( MNXYZ0+2 )
            MCN( MNBLGST+  NBPOINTS ) = 1
            MCN( MNNOLGP+2*NBPOINTS ) = NL
            NBPOINTS = NBPOINTS + 1
 340     CONTINUE
         CALL TNMCDS( 'REEL', 3*NBS1, MNXYZS1 )
C
 345     CALL TNMCDS( 'ENTIER', NBSOM, MNNBAS )
 350  CONTINUE
C
C     ICI NBPOINTS EST LE NOMBRE DE XYZ EXTREMITES DES LIGNES IDENTIFIES
      IF( NBPOINTS .GT. 0 ) THEN
      CALL TNMCDC( 'ENTIER', NBPOINTS, MNNOLXP )
C
      NBP = 0
      DO 370 K=1,NBPOINTS
CCC         IF( MCN( MNBLGST-1+K ) .GE. 2 ) THEN
C           CE SOMMET PRODUIT UN POINT => TMS DEFINITION ET XYZSOMMET
            NBP = NBP + 1
C
C           NOM DU POINT
            N = NBCHIF( NBP )
            KERR(1) = '(I )'
            WRITE( KERR(1)(3:3), '(I1)' ) N
            KNOM = 'P     '
            WRITE( KNOM(2:1+N), KERR(1)(1:4) ) NBP
            M = N + 2
C
            NF1 = MCN( MNNOLGP-2+2*K )
            IF( NF1 .NE. 0 ) THEN
               KNOM(M:M+1) = '_L'
               M = M + 2
               N = NBCHIF( NF1 )
               KERR(1) = '(I )'
               WRITE( KERR(1)(3:3), '(I1)' ) N
               WRITE( KNOM(M:M-1+N), KERR(1)(1:4) ) NF1
               M = M + N
            ENDIF
C
            NF2 = MCN( MNNOLGP-1+2*K )
            IF( NF2 .NE. 0 ) THEN
               KNOM(M:M+1) = '_L'
               M = M + 2
               N = NBCHIF( NF2 )
               KERR(1) = '(I )'
               WRITE( KERR(1)(3:3), '(I1)' ) N
               WRITE( KNOM(M:M-1+N), KERR(1)(1:4) ) NF2
            ENDIF
C
C           LEXIQUE DU POINT
            CALL LXLXDC( NTPOIN, KNOM, 24, 8 )
            CALL LXLXOU( NTPOIN, KNOM, NTPTOB, MNPTOB )
C
C           RECHERCHE DU NUMERO DU NOM KNOM DANS LE LEXIQUE DES POINTS
            CALL LXNMNO( NTPOIN, KNOM, MCN(MNNOLXP-1+NBP), MNLX )
C
C           LE TMS DEFINITION DU POINT
            CALL LXTNDC( NTPTOB, 'DEFINITION', 'MOTS', WOORPO+3 )
            CALL LXTSOU( NTPTOB, 'DEFINITION', NTDFPT, MNDFPT )
            MCN( MNDFPT + WUTYPO ) = 1
            MN0 = MNSTL - 3 + 3 * K
            RMCN( MNDFPT + WOORPO     ) = RMCN( MN0 )
            RMCN( MNDFPT + WOORPO + 1 ) = RMCN( MN0+1 )
            RMCN( MNDFPT + WOORPO + 2 ) = RMCN( MN0+2 )
C           LA DATE DU TMS XYZSOMMET DE LA LIGNE NL
            CALL ECDATE( MCN(MNDFPT) )
C           LE NUMERO DU TABLEAU DESCRIPTEUR
            MCN( MNDFPT + MOTVAR(6) ) = NONMTD( '~>POINT>>DEFINITION' )
C
C           LE TMS XYZSOMMET DU POINT K
            CALL LXTNDC( NTPTOB, 'XYZSOMMET', 'MOTS', WYZSOM+3 )
            CALL LXTSOU( NTPTOB, 'XYZSOMMET', NTXYZ, MNXYZ )
            MCN( MNXYZ + WBCOOR ) = 3
            MCN( MNXYZ + WNBSOM ) = 1
            MCN( MNXYZ + WNBTGS ) = 0
            MN0 = MNSTL - 3 + 3 * K
            MN1 = MNXYZ + WYZSOM
            RMCN(MN1  ) = RMCN( MN0   )
            RMCN(MN1+1) = RMCN( MN0+1 )
            RMCN(MN1+2) = RMCN( MN0+2 )
C           LA DATE DU TMS XYZSOMMET DE LA LIGNE NL
            CALL ECDATE( MCN(MNXYZ) )
C           LE NUMERO DU TABLEAU DESCRIPTEUR
            MCN( MNXYZ + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
CCC         ENDIF
 370  CONTINUE
      ENDIF
C
      CALL TNMCDS( 'ENTIER', 2*MXPOINTS, MNNOLGP )
      CALL TNMCDS( 'REEL',   3*MXPOINTS, MNSTL   )
      CALL TNMCDS( 'ENTIER',   MXPOINTS, MNBLGST )
      CALL TNMCDS( 'ENTIER',   NBLIGNES, MNLXYZ )
C
C
C     CONSTRUCTION DU TABLEAU DEFINITION DE L'OBJET FINAL POUR MEFISTO
 400  L1 = MOTVAR( NUMTYP('TYPEOBJET') )
      NBDOBJ = NBSURFACES + NBLIGNES + NBPOINTS
      NBOBPR = NBDOBJ
      CALL LXTNDC( NTLXOB, 'DEFINITION', 'MOTS', WTYOBJ+L1*NBDOBJ )
      CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )
C
C     mode de traitement NDOUNO: classique 2D ou 3D
      MCN( MNDFOB+WDOUNO ) = 0
C     nombre des Points Lignes Surfaces Volumes Objets de l'objet
      MCN( MNDFOB+WBDOBJ ) = NBDOBJ
C
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) 'OBJET: ', KNOMOB, '============================='
      WRITE(IMPRIM,*) 'NOMBRE DE SURFACES =',NBSURFACES
      MNOB = MNDFOB+WBDOBJ
      DO 530 N=1,NBSURFACES
C        LE NUMERO DU TYPE DE PLSV: SURFACE
         MCN( MNOB + 1 ) = 3
C        LE NUMERO DU PLSV DANS SON LEXIQUE
         MCN( MNOB + 2 ) = MCN(MNNOLXS-1+N)
         MNOB = MNOB + 2
 530  CONTINUE
C
      WRITE(IMPRIM,*) 'NOMBRE DE LIGNES   =',NBLIGNES
      DO 550 N=1,NBLIGNES
C        LE NUMERO DU TYPE DE PLSV: LIGNE
         MCN( MNOB + 1 ) = 2
C        LE NUMERO DU PLSV DANS SON LEXIQUE
         MCN( MNOB + 2 ) = MCN(MNNOLXL-1+N)
         MNOB = MNOB + 2
 550  CONTINUE
C
      WRITE(IMPRIM,*) 'NOMBRE DE POINTS   =',NBPOINTS
      DO 560 N=1,NBPOINTS
C        LE NUMERO DU TYPE DE PLSV: POINT
         MCN( MNOB + 1 ) = 1
C        LE NUMERO DU PLSV DANS SON LEXIQUE
         MCN( MNOB + 2 ) = MCN(MNNOLXP-1+N)
         MNOB = MNOB + 2
 560  CONTINUE
C
C     LA DATE
      CALL ECDATE( MCN(MNDFOB) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNDFOB + MOTVAR(6) ) = NONMTD( '~>OBJET>>DEFINITION' )
C
C     DESTRUCTIONS DE TABLEAUX MC
      CALL TNMCDS( 'ENTIER', L1ARET*L2ARET, MNARET  )
      IF( NBSURFACES .GT. 0 ) THEN
         CALL TNMCDS( 'ENTIER', NBSURFACES,    MNBTEV  )
         CALL TNMCDS( 'ENTIER', NBSURFACES,    MNNOLXS )
         CALL TNMCDS( 'ENTIER', NBSURFACES,    NTSFOBS )
         CALL TNMCDS( 'ENTIER', NBSURFACES,    NTNSEFS )
         CALL TNMCDS( 'ENTIER', NBSURFACES,    MNNSEFS )
         CALL TNMCDS( 'ENTIER',(NBSURFACES+1)**2, MNFAFA )
      ENDIF
      IF( NEWOLD .GT. 0 ) CALL TNMCDS( 'ENTIER', 1+NBSOMMET, NEWOLD )
      IF( NBLIGNES .GT. 0 ) THEN
         CALL TNMCDS( 'ENTIER', NBLIGNES,      NTLIOBL )
         CALL TNMCDS( 'ENTIER', NBLIGNES,      NTNSEFL )
         CALL TNMCDS( 'ENTIER', NBLIGNES,      MNNSEFL )
         CALL TNMCDS( 'ENTIER', NBLIGNES,      MNBARL  )
         CALL TNMCDS( 'ENTIER', NBLIGNES,      MNNOLXL )
      ENDIF
      IF( NBPOINTS.GT.0 ) THEN
         CALL TNMCDS( 'ENTIER', NBPOINTS, MNNOLXP )
      ENDIF

      RETURN
      END
