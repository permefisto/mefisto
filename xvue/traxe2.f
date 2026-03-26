      SUBROUTINE TRAXE2
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACER DES 2 AXES SELON LES COORDONNEES XOBMIN XOBMAX et
C ----- YOBMIN YOBMAX DU COMMON TRVARI          VERSION xvue
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS       AOUT 1997
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/traaxe.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      CHARACTER*16      KFORM
      CHARACTER*24      KVAL0, KVAL
C
      IF( INTERA .LE. 0 ) RETURN
C
      IF( NTRAXE .LE. 0 ) GOTO 9900
ccc      IF( NETAXE .GT. 0 ) GOTO 9900
C
C     L'ECRAN A POUR COORDONNEES OBJETS XOBMIN XOBMAX YOBMIN YOBMAX
C     LE TRACE DE X ET Y EN LIGNES EPAISSIES 2 FOIS
      CALL XVEPAISSEUR( 3 )
      IF( NDCOUL .EQ. N1COUL+1 ) THEN
C        TRACE BLANC DES AXES
C        NOIR ET BLANC
         NCX = NCBLAN
         NCY = NCBLAN
      ELSE
C        TRACE DES AXES
         NCX = NCROUG
         NCY = NCVERT
      ENDIF
C
C     ECART EN % DEFINIT LA MARGE ENTRE LES AXES ET LE BORD ECRAN
      ECARTY = (YOBMAX - YOBMIN) * 0.006
      ECARTX = (XOBMAX - XOBMIN) * 0.006
C
C     L'AXE DES X
C     -----------
      CALL TRAIT2D( NCX, XOBMAX, YOBMIN+ECARTY*3,
     %                   XOBMIN, YOBMIN+ECARTY*3 )
C     LE CARACTERE X
      CALL TEXTE2D( NCNOIR, XOBMAX-ECARTX*6, YOBMIN+ECARTY*7, '-X->' )
C
C     L'AXE DES Y
C     -----------
      CALL TRAIT2D( NCY, XOBMIN+ECARTX, YOBMIN,
     %                   XOBMIN+ECARTX, YOBMAX )
C     LE CARACTERE Y
      CALL TEXTE2D( NCNOIR, XOBMIN+ECARTX*2,
     %                      YOBMAX-ECARTY*6, 'Y' )
C
C     LE TRACE DES POINTS SUR L'AXE X
C     -------------------------------
      CALL GRADUA( 1.0, XOBMIN, XOBMAX, PREMIER, PAS, NBPAS, NBDFOR )
C     LE FORMAT D'ECRITURE DE CHAQUE VALEUR
      WRITE( KFORM,10010 ) NBDFOR+8,NBDFOR
10010 FORMAT('(G',I2,'.',I1,') ')
ccc      print *,'traxe2 X: KFORM=',KFORM,'  NBDFOR=',NBDFOR
C     L'ORDONNEE DE L'AXE DES X
      Y = YOBMIN + ECARTY*3
C
      DO 10 I=0,NBPAS
C        L'ABSCISSE DU NOMBRE A TRACER
         X = PREMIER + I * PAS
C        TRACE D'UN TRAIT VERTICAL
         CALL TRAIT2D( NCX, X, Y-ECARTY, X, Y+ECARTY )
C        LA VALEUR A TRACER EN EVITANT LES 1E-7 ...
         IF( ABS(X) .LT. 1E-3 * ABS(PAS) ) X=0
         WRITE( KVAL0, KFORM ) X
         CALL REELCA( KVAL0, KVAL )
ccc      print *,'traxe2 Y: KFORM=',KFORM,'  X=',X
         CALL SANSBL( KVAL, NBC )
         CALL TEXTE2D( NCGRIS, X, Y+ECARTY*0.5, KVAL(1:NBC) )
 10   CONTINUE
C
C     LE TRACE DES POINTS SUR L'AXE Y
C     -------------------------------
      CALL GRADUA( 1.0, YOBMIN, YOBMAX, PREMIER, PAS, NBPAS, NBDFOR )
C     LE FORMAT D'ECRITURE DE CHAQUE VALEUR
      WRITE( KFORM,10010 ) NBDFOR+8,NBDFOR
ccc      print *,'traxe2 Y: KFORM=',KFORM,'  NBDFOR=',NBDFOR
C     L'ABSCISSE DU TRACE DE L'AXE DES Y
      X = XOBMIN + ECARTX
C
      DO 20 I=0,NBPAS
C        L'ORDONNEE DU NOMBRE A TRACER
         Y = PREMIER + I * PAS
C        TRACE D'UN TRAIT HORIZONTAL
         CALL TRAIT2D( NCY, X-ECARTX, Y, X+ECARTX, Y )
C        LA VALEUR A TRACER EN EVITANT LES 1E-7 ...
         IF( ABS(Y) .LT. 1E-3 * ABS(PAS) ) Y=0
ccc         print *,'traxe2 Y: KFORM=',KFORM,'  Y=',Y
         WRITE( KVAL0, KFORM ) Y
         CALL REELCA( KVAL0, KVAL )
ccc         print *,'traxe2 Y: KFORM=',KFORM,'  Y=',Y,' KVAL=',KVAL
         CALL SANSBL( KVAL, NBC )
         CALL TEXTE2D( NCGRIS, XOBMIN+ECARTX*0.8, Y+ECARTY*0.5,
     %                 KVAL(1:NBC) )
 20   CONTINUE
C
CCCC     TRACE SUR LA FENETRE
CCC 9900 CALL MEMPXFENETRE   ici CCC EVITE LE SCINTILLEMENT
C
 9900 CALL XVEPAISSEUR( 0 )
ccc      NETAXE = 2  ne sert pas ailleurs...
      RETURN
      END
