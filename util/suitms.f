      SUBROUTINE SUITMS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SUIVI DES FICHIERS DE LA MEMOIRE SECONDAIRE
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1990
C2345X7..............................................................012
      include"./incl/pp.inc"
      include"./incl/ppmck.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/msvaau.inc"
      include"./incl/trvari.inc"
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
C
      REAL              RMCN(MOTMCN)
      COMMON             MCN(MOTMCN)
      EQUIVALENCE       (MCN(1),RMCN(1))
C
      COMMON / UNITES /  LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / MSSFTA /  MSSF(28),NTADAM
C
      CHARACTER*160      KNOM, KNOMFI
      CHARACTER*48       FICHNM
      CHARACTER*9        TYPNUM
C
C     LECTURE DU MOT CLE A TRAITER
C     ----------------------------
 10   CALL LIMTCL( 'suivitms', NMTCL )
      IF( NMTCL .LE. 0 ) GOTO 9000
      GOTO( 100, 200, 300, 400, 500, 600, 700), NMTCL
C
C     TABLEAU TMS_A_AFFICHER
C     ======================
 100  CALL INVITE( 57 )
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNOM )
      IF( NCVALS .EQ. -1 ) GOTO 9000
      CALL AFTSTD( KNOM )
      GOTO 10
C
C     ENTREE_TMS
C     ==========
 200  CALL INVITE( 58 )
      NCVALS = 0
      CALL LIRCAR(NCVALS, KNOM )
      IF( NCVALS .EQ. -1 ) GOTO 9000
C     DEDUCTION DU TABLEAU DESCRIPTEUR DE CE TABLEAU NOM DE TMS
      CALL NUTDTS( KNOM, NUTD )
      IF( NUTD .LE. 0 ) THEN
         KERR(1) = 'PAS DE TABLEAU DESCRIPTEUR DE'
         N = NUDCNB( KNOM )
         KERR(2) = KNOM(1:N)
         NBLGRC(NRERR) = 2
         CALL LEREUR
         GOTO 10
      ENDIF
C     ENTREE DU TMS
      CALL MOTSTD( DICOTD(NUTD), KNOM, NRETOU )
      GOTO 10
C
C     AFFICHER LA VALEUR D'UNE VARIABLE DE TABLEAU TMS
C     ================================================
 300  CALL INVITE( 64 )
      NCVALS = 0
      CALL LIRCAR(NCVALS, KNOM )
      IF( NCVALS .EQ. -1 ) GOTO 9000
      CALL VATSTD( KNOM, NOTYPE, MNTMS, LDTMS, N )
      IF( MNTMS .LE. 0 ) THEN
         N = NUDCNB( KNOM )
         KERR(1)  = ' VARIABLE INEXISTANTE:' // KNOM(1:N)
         NBLGRC(NRERR) = 1
         CALL LEREUR
         GOTO 10
      ENDIF
C     AFFICHAGE DE LA VARIABLE
      CALL AFVATS( KNOM, NOTYPE, MNTMS, LDTMS )
      GOTO 10
C
C     IMPOSER LA VALEUR D'UNE VARIABLE D'UN TMS
C     =========================================
 400  CALL INVITE( 65 )
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNOM )
      IF( NCVALS .EQ. -1 ) GOTO 9000
      CALL VATSTD( KNOM, NOTYPE, MNTMS, LDTMS, N )
      IF( MNTMS .LE. 0 ) THEN
         N = NUDCNB( KNOM )
         KERR(1) = 'VARIABLE INEXISTANTE:' // KNOM(1:N)
         NBLGRC(NRERR) = 1
         CALL LEREUR
         GOTO 10
      ENDIF
C     AFFICHAGE DE LA VALEUR ACTUELLE DE LA VARIABLE
      CALL AFVATS( KNOM, NOTYPE, MNTMS, LDTMS )
C     AFFECTATION DE LA NOUVELLE VALEUR DE LA VARIABLE
      CALL INVITD( 'NOUVELLE VALEUR ' // TYPNUM(NOTYPE) )
      NCVALS = 0
      IF( NOTYPE .EQ. 3 .OR. NOTYPE .EQ. 4 .OR. NOTYPE .EQ. 11) THEN
         CALL LIRENT( NCVALS, MCN(MNTMS+LDTMS) )
      ELSE IF( NOTYPE .EQ. 5 ) THEN
         CALL LIRRSP( NCVALS, RMCN(MNTMS+LDTMS) )
      ELSE IF( NOTYPE .EQ. 6 ) THEN
         CALL LIRRDP( NCVALS, RMCN(MNTMS+LDTMS) )
      ELSE IF( NOTYPE .EQ. 12 ) THEN
         CALL LIRRSP( NCVALS, RMCN(MNTMS+LDTMS) )
         CALL LIRRSP( NCVALS, RMCN(MNTMS+LDTMS+1) )
         CALL LIRRSP( NCVALS, RMCN(MNTMS+LDTMS+2) )
      ELSE
         WRITE(KERR(MXLGER)(1:4),'(I4)') NOTYPE
         KERR(1) = 'TYPE NON TRAITE' // KERR(MXLGER)(1:4)
         NBLGRC(NRERR) = 1
         CALL LEREUR
         GOTO 10
      ENDIF
      IF( NCVALS .EQ. -1 ) GOTO 9000
C     LA DATE DU TABLEAU EST CHANGEE
      CALL ECDATE( MCN(MNTMS) )
      GOTO 10
C
C     CHANGER UN NOM DANS UN LEXIQUE
C     ==============================
 500  CALL INVITE( 50 )
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNOM )
      IF( NCVALS .LE. 0 ) GOTO 10
C     OUVERTURE DU LEXIQUE
      CALL LXXXOU( NTADAM, KNOM, NTLX, MNLX )
      IF( NTLX .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         N = NUDCNB( KNOM )
         KERR(1) = KNOM(1:N)
         KERR(2) = 'LEXIQUE INCONNU'
         CALL LEREUR
         GOTO 500
      ENDIF
C     L'ANCIEN NOM
      CALL INVITE( 4 )
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNOM )
      IF( NCVALS .LE. 0 ) GOTO 10
C     CE NOM EXISTE-T-IL ?
      CALL LXNMNO( NTLX, KNOM, N, MN )
      IF( N .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         N = NUDCNB( KNOM )
         KERR(1) = KNOM(1:N)
         KERR(2) = 'ANCIEN NOM INCONNU'
         CALL LEREUR
         GOTO 500
      ENDIF
C     LE NOUVEAU NOM
      CALL INVITE( 81 )
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNOMFI )
      IF( NCVALS .LE. 0 ) GOTO 10
C     LE CHANGEMENT DE NOM
      CALL LXNMNM( NTLX, KNOM, KNOMFI )
      GOTO 10
C
C     TRANSFERT DE TMS VERS UN AUTRE ORDINATEUR
C     =========================================
C     LE NUMERO D'UNITE LOGIQUE
 600  CALL TRUNIT( NFFICA )
      KNOMFI = 'CRAY'
      OPEN( UNIT=NFFICA, STATUS='UNKNOWN',
     %      FILE=KNOMFI, ACCESS='SEQUENTIAL',
     %      FORM='FORMATTED', IOSTAT=I )
      IF( I .NE. 0 ) THEN
         KERR(1) ='FICHIER POUR LE CRAY NON OUVRABLE'
         N = NUDCNB( KNOMFI )
         KERR(2) = KNOMFI(1:N)
         NBLGRC(NRERR) = 2
         CALL LEREUR
         GOTO 10
      ENDIF
C
C     POSITION EN TETE DE FICHIER
      REWIND( UNIT=NFFICA, IOSTAT=I )
      IF( I .NE. 0 ) THEN
         N = NUDCNB( KNOMFI )
         KERR(1) = KNOMFI(1:N)
         KERR(2) = 'FICHIER NON REMBOBINABLE'
         NBLGRC(NRERR) = 1
         CALL LEREUR
         GOTO 10
      ENDIF
C
 660  CALL INVITE( 59 )
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNOM )
      IF( NCVALS .EQ. -1 ) THEN
C        FERMETURE DU FICHIER CRAY
         CLOSE( UNIT=NFFICA, STATUS='KEEP' )
C
CCCC        VERIFICATION DU CONTENU DU FICHIER CRAY
CCC         OPEN( UNIT=NFFICA, STATUS='OLD',
CCC     %         FILE=KNOMFI, ACCESS='SEQUENTIAL',
CCC     %         FORM='FORMATTED', IOSTAT=I )
CCC         IF( I .NE. 0 ) THEN
CCC            KERR(1) = 'FICHIER ' // KNOMFI //
CCC     %                ' ECRIT POUR LE CRAY NON OUVRABLE'
CCC            NRETOU = 1
CCC            NBLGRC(NRERR) = 1
CCC            CALL LEREUR
CCC            GOTO 690
CCC         ENDIF
CCCC
CCC 670     READ( UNIT=NFFICA, FMT='(A)', IOSTAT=I ) KFICA
CCC         IF( I .EQ. 0 ) THEN
CCC            KERR(1) = KFICA
CCC            NBLGRC(NRERR) = 1
CCC            CALL LEREUR
CCC            GOTO 670
CCC         ENDIF
CCC         CLOSE( UNIT=NFFICA, STATUS='KEEP' )
         GOTO 10
      ENDIF
C
C     CONVERSION DES VALEURS NUMERIQUES DU TMS KNOM
C     ECRITURE DES CARACTERES SUR LE FICHIER FICA
      LIRECR = 0
      LCFICA = 0
      CALL NCTSTD( LIRECR, KNOM, NFFICA, NRETOU )
      IF( NRETOU .NE. 0 ) THEN
         KERR(1) = 'FICHIER CRAY INCORRECT ET DETRUIT'
         NBLGRC(NRERR) = 1
         CALL LEREUR
         CLOSE( UNIT=NFFICA, STATUS='DELETE' )
         GOTO 10
      ENDIF
      GOTO 660
C
C     RECUPERATION DES TMS VENUS D'UN AUTRE ORDINATEUR
C     ================================================
 700  CALL TRUNIT( NFFICA )
C     LE NOM DU FICHIER : NOM CONDENSE DU PROJET SUFFIXE PAR 'CIBLE'
      KNOMFI = FICHNM(9)
C     LES CARACTERES AU DELA DE . SONT ECRASES PAR 'CIBLE'
      I = INDEX( KNOMFI, '.' )
      KNOMFI(I+1:I+6) = 'CIBLE'
C
      OPEN( UNIT=NFFICA, STATUS='OLD',
     %      FILE=KNOMFI, ACCESS='SEQUENTIAL',
     %      FORM='FORMATTED', IOSTAT=I )
      IF( I .NE. 0 ) THEN
         KERR(1) ='FICHIER VENU DU CRAY NON OUVRABLE'
         N = NUDCNB( KNOMFI )
         KERR(2) = KNOMFI(1:N)
         NBLGRC(NRERR) = 2
         CALL LEREUR
         GOTO 10
      ENDIF
C
      CALL CNTSTD( NFFICA, NRETOU )
      IF( NRETOU .EQ. 1001 ) THEN
         NRETOU = 0
      ELSE IF( NRETOU .NE. 0 ) THEN
         KERR(1) = 'PROBLEME DANS LE FICHIER CIBLE'
         NBLGRC(NRERR) = 1
         CALL LEREUR
      ENDIF
      GOTO 10
C
 9000 RETURN
      END
